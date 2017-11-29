/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:

    Arno Mayrhofer (CFDEMresearch GmbH, Linz)
    Christoph Kloss (DCS Computing GmbH, Linz)

    Copyright 2016-     CFDEMresearch GmbH, Linz
    Copyright 2014-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#include <mpi.h>
#include <string.h>
#include <stdio.h>
#include "fix_contact_property_atom.h"
#include "fix_property_atom.h"
#include "pair_gran.h"
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "force.h"
#include "pair.h"
#include "update.h"
#include "timer.h"
#include "modify.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixContactPropertyAtom::FixContactPropertyAtom(LAMMPS *lmp, int narg, char **arg) :
  FixContactHistory(lmp, narg, arg),
  fix_nneighs_full_(0),
  build_neighlist_(true),
  reset_each_ts_(true)
{
    bool hasargs = true;
    while(iarg_ < narg && hasargs)
    {
        hasargs = false;
        if (strcmp(arg[iarg_],"reset") == 0) {
            if (narg < iarg_+2)
                error->fix_error(FLERR,this,"not enough arguments for keyword 'reset'");
            iarg_++;
            if (strcmp(arg[iarg_], "no") == 0)
                reset_each_ts_ = false;
            else if (strcmp(arg[iarg_], "yes") != 0)
                error->fix_error(FLERR,this,"keyword 'reset' requires yes or no as value");
            iarg_++;
            hasargs = true;
        } else if(strcmp(style,"contactproperty/atom") == 0) {
            char *errmsg = new char[strlen(arg[iarg_])+50];
            sprintf(errmsg,"unknown keyword or wrong keyword order: %s", arg[iarg_]);
            error->fix_error(FLERR,this,errmsg);
            delete []errmsg;
        }
    }
    // avoid reseting npartners after the initial compute in setup
    just_created = false;

    // TODO throw error if newton is set
}

/* ---------------------------------------------------------------------- */

FixContactPropertyAtom::~FixContactPropertyAtom()
{
}

/* ---------------------------------------------------------------------- */

void FixContactPropertyAtom::post_create()
{
    // do not do this for wall
    if(strcmp(style,"contactproperty/atom"))
        return;

    fix_nneighs_full_ = static_cast<FixPropertyAtom*>(modify->find_fix_id("nneighs_full"));
    if (!fix_nneighs_full_)
    {
        const char **fixarg = new const char*[9];
        fixarg[0]="nneighs_full";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="nneighs_full";
        fixarg[4]="scalar";
        fixarg[5]="no"; // restart
        fixarg[6]="yes"; // fw comm
        fixarg[7]="no"; //rev comm
        fixarg[8]="0.0";
        modify->add_fix(9,const_cast<char**>(fixarg));
        fix_nneighs_full_ = static_cast<FixPropertyAtom*>(modify->find_fix_id("nneighs_full"));
        delete []fixarg;
    }
}

/* ---------------------------------------------------------------------- */

int FixContactPropertyAtom::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= MIN_PRE_FORCE;
  mask |= PRE_NEIGHBOR;
  mask |= PRE_EXCHANGE;
  mask |= MIN_PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixContactPropertyAtom::init()
{
  FixContactHistory::init();

  if(!force->pair_match("gran", 0))
          error->fix_error(FLERR,this,"Please use a granular pair style for this fix");
  pair_gran_ = static_cast<PairGran*>(force->pair_match("gran", 0));

  comm_forward = 20*dnum_;
}

/* ---------------------------------------------------------------------- */

void FixContactPropertyAtom::setup_pre_exchange()
{
    // in case of run 100, run 100 is executed after first sequence of run 100
    
    //if(0 == neighbor->ago)
    //    pre_exchange();
}
/* ---------------------------------------------------------------------- */

void FixContactPropertyAtom::min_setup_pre_exchange()
{
    
    if(0 == neighbor->ago)
        pre_exchange();
}

/* ---------------------------------------------------------------------- */

void FixContactPropertyAtom::pre_exchange()
{
   // set maxtouch = max # of partners of any owned atom
   // bump up comm->maxexchange_fix if necessary

   int nlocal = atom->nlocal;

   maxtouch_ = 0;
   for (int i = 0; i < nlocal; i++) maxtouch_ = MAX(maxtouch_,npartner_[i]);
   comm->maxexchange_fix = MAX(comm->maxexchange_fix,(dnum_+1)*maxtouch_+1);
}

/* ----------------------------------------------------------------------
   for wall:
   need to execute here since neighlist is refreshed in setup_pre_force(int foo)
   if this is not done, data will get out-of-sync

   need not care about overwriting hist values since always have pointer
   to most current data
------------------------------------------------------------------------- */

void FixContactPropertyAtom::setup_pre_force(int dummy)
{
   pre_neighbor();
   pre_force(0);
}

/* ---------------------------------------------------------------------- */

void FixContactPropertyAtom::pre_neighbor()
{
    build_neighlist_ = true;
    
}

/* ----------------------------------------------------------------------
   allocate storage
------------------------------------------------------------------------- */

void FixContactPropertyAtom::pre_force(int dummy)
{
    if (reset_each_ts_)
        clear();
}

/* ---------------------------------------------------------------------- */

void FixContactPropertyAtom::clear()
{
    int nall = atom->nlocal+atom->nghost;
    //double **x = atom->x;

    double *nneighs_full = fix_nneighs_full_->vector_atom;

    // reset number of partners every time-step
    vectorZeroizeN(npartner_,nall);

    // other stuff to do only upon neigh list rebuild
    if(!build_neighlist_)
        return;
    build_neighlist_ = false;

    int nneighs_next;

    int inum = pair_gran_->list->inum;
    int * ilist = pair_gran_->list->ilist;
    int * numneigh = pair_gran_->list->numneigh;
    int ** firstneigh = pair_gran_->list->firstneigh;

    ipage_->reset();
    dpage_->reset();

    vectorZeroizeN(nneighs_full,nall);

    // tally # of neighs to full nneighs

    for (int ii = 0; ii < inum; ii++) {
      const int i = ilist[ii];
      int * const jlist = firstneigh[i];
      const int jnum = numneigh[i];
      for (int jj = 0; jj < jnum; jj++) {
        const int j = jlist[jj] & NEIGHMASK;

        nneighs_full[i] += 1.;
        nneighs_full[j] += 1.;
      }
    }

    fix_nneighs_full_->do_forward_comm();

    // allocate mem for owned and ghost
    for (int i = 0; i < nall; i++)
    {
        nneighs_next = static_cast<int>(nneighs_full[i]);

        npartner_[i] = 0;

        partner_[i] = ipage_->get(nneighs_next);
        contacthistory_[i] = dpage_->get(nneighs_next*dnum_);
        
        if (0 == partner_[i] || 0 == contacthistory_[i])
            error->one(FLERR,"Contact history overflow, boost neigh_modify one");
        vectorInitializeN(partner_[i],nneighs_next,-1);
        vectorZeroizeN(contacthistory_[i],nneighs_next*dnum_);

   }
}

/* ---------------------------------------------------------------------- */

void FixContactPropertyAtom::min_setup_pre_force(int dummy)
{
  pre_force(0);
}

/* ---------------------------------------------------------------------- */

void FixContactPropertyAtom::min_pre_force(int dummy)
{
  pre_force(0);
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixContactPropertyAtom::grow_arrays(int nmax)
{
  FixContactHistory::grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixContactPropertyAtom::copy_arrays(int i, int j, int delflag)
{
  // just copy pointers for partner and shearpartner
  // b/c can't overwrite chunk allocation inside ipage,dpage
  // incoming atoms in unpack_exchange just grab new chunks
  // so are orphaning chunks for migrating atoms
  // OK, b/c will reset ipage,dpage on next reneighboring

  FixContactHistory::copy_arrays(i,j,delflag);
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixContactPropertyAtom::unpack_exchange(int nlocal, double *buf)
{
  // allocate new chunks from ipage,dpage for incoming values
  
  int m = 0;

  npartner_[nlocal] = ubuf(buf[m++]).i;
  maxtouch_ = MAX(maxtouch_,npartner_[nlocal]);
  partner_[nlocal] = ipage_->get(npartner_[nlocal]);
  contacthistory_[nlocal] = dpage_->get(dnum_*npartner_[nlocal]);

  if (partner_[nlocal] == NULL || contacthistory_[nlocal] == NULL)
      error->one(FLERR,"Contact history overflow, boost neigh_modify one");

  for (int n = 0; n < npartner_[nlocal]; n++) {
    partner_[nlocal][n] = ubuf(buf[m++]).i;
    
    for (int d = 0; d < dnum_; d++) {
      contacthistory_[nlocal][n*dnum_+d] = buf[m++];
      
    }
  }

  /*
  for (int n = npartner_[nlocal]; n < nneighs; n++) {
    partner_[nlocal][n] = -1;
    for (int d = 0; d < dnum_; d++) {
      contacthistory_[nlocal][n*dnum_+d] = 0.;
    }
  }
  */

  return m;
}

/* ----------------------------------------------------------------------
   forward and backward comm to be used by other fixes as needed
------------------------------------------------------------------------- */

void FixContactPropertyAtom::do_forward_comm()
{
    timer->stamp();
    comm->forward_comm_variable_fix(this);
    timer->stamp(TIME_COMM);
}

/* ---------------------------------------------------------------------- */

int FixContactPropertyAtom::pack_comm(int n, int *list, double *buf,
                             int pbc_flag, int *pbc)
{
    int i,j;

    //we dont need to account for pbc here
    int m = 0;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = ubuf(npartner_[j]).d;
      for (int np = 0; np < npartner_[j]; np++) {
        buf[m++] = ubuf(partner_[j][np]).d;
        for (int d = 0; d < dnum_; d++) {
          buf[m++] = contacthistory_[j][np*dnum_+d];
        }
      }
    }
    return m;
}

/* ---------------------------------------------------------------------- */

void FixContactPropertyAtom::unpack_comm(int n, int first, double *buf)
{
   
      int i,m,last;
      m = 0;
      last = first + n;
      for (i = first; i < last; i++) {
           npartner_[i] = ubuf(buf[m++]).i;
           for (int np = 0; np < npartner_[i]; np++) {
               partner_[i][np] = ubuf(buf[m++]).i;
               for (int d = 0; d < dnum_; d++) {
                  contacthistory_[i][np*dnum_+d] = buf[m++];
               }
           }
      }
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixContactPropertyAtom::unpack_restart(int nlocal, int nth)
{
  // ipage = NULL if being called from granular pair style init()

  if (ipage_ == NULL) allocate_pages();

  // skip to Nth set of extra values

  double **extra = atom->extra;

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  // allocate new chunks from ipage,dpage for incoming values
  
  int d;

  npartner_[nlocal] = ubuf(extra[nlocal][m++]).i;
  maxtouch_ = MAX(maxtouch_,npartner_[nlocal]);
  partner_[nlocal] = ipage_->get(npartner_[nlocal]);
  contacthistory_[nlocal] = dpage_->get(npartner_[nlocal]*dnum_);

  if (partner_[nlocal] == NULL || contacthistory_[nlocal] == NULL)
      error->one(FLERR,"Contact history overflow, boost neigh_modify one");

  for (int n = 0; n < npartner_[nlocal]; n++) {
    partner_[nlocal][n] = ubuf(extra[nlocal][m++]).i;
    for (d = 0; d < dnum_; d++) {
      contacthistory_[nlocal][n*dnum_+d] = extra[nlocal][m++];
      
    }
  }
}

/* ----------------------------------------------------------------------
   pack state of Fix into one write
------------------------------------------------------------------------- */

void FixContactPropertyAtom::write_restart(FILE *fp)
{
    FixContactHistory::write_restart(fp);
}
