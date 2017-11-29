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
    This file is from LAMMPS, but has been modified. Copyright for
    modification:

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz

    Copyright of original file:
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

#include <mpi.h>
#include <string.h>
#include <stdio.h>
#include <cmath>
#include <algorithm>
#include "fix_contact_history.h"
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "force.h"
#include "pair_gran.h"
#include "update.h"
#include "modify.h"
#include "memory.h"
#include "math_extra_liggghts.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixContactHistory::FixContactHistory(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  dnum_(0),
  variablename_(0),
  newtonflag_(0),
  history_id_(0),
  index_decide_noncontacting_(-1),
  npartner_(0),
  partner_(0),
  contacthistory_(0),
  maxtouch_(0),
  pair_gran_(0),
  computeflag_(0),
  pgsize_(0),
  oneatom_(0),
  ipage_(0),
  dpage_(0)
{
  restart_global = 1;
  restart_peratom = 1;
  create_attribute = 1;

  // perform initial allocation of atom-based arrays
  // register with atom class

  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);

  // initialize npartner to 0 so neighbor list creation is OK the 1st time

  std::fill_n(npartner_, atom->nmax, 0);

  //=====================
  // parse args
  //=====================

  if(narg < 4)
   error->fix_error(FLERR,this,"not enough parameters");

  iarg_ = 3;

  if(!MathExtraLiggghts::is_int(arg[iarg_]))
  {
    int n = strlen(arg[iarg_]) + 1;
    variablename_ = new char[n];
    strcpy(variablename_,arg[iarg_]);
    iarg_++;
    
  }
  else
  {
    int n = strlen("contacthistory") + 1;
    variablename_ = new char[n];
    strcpy(variablename_,"contacthistory");
  }

  // read dnum
  dnum_ = atoi(arg[iarg_++]);
  
  if(dnum_ < 0)
    error->fix_error(FLERR,this,"dnum must be >=0");

  // do not proceed for derived classes
  
  if(!strstr(style,"property") && strcmp(style,"contacthistory"))
    return;

  // parse args
  if(narg < 6)
    error->fix_error(FLERR,this,"not enough parameters");

  // read newtonflag
  if(narg-iarg_ < 2*dnum_)
     error->fix_error(FLERR,this,"not enough parameters: need to specify an id and a newtonflag for each dnum");

  newtonflag_ = new int[dnum_];
  history_id_ = (char**) memory->smalloc((dnum_)*sizeof(char*),"FixContactHistory:history_id");

  for(int i = 0 ; i < dnum_; i++)
  {
    
    history_id_[i] = new char[strlen(arg[iarg_])+1];
    strcpy(history_id_[i],arg[iarg_++]);
    newtonflag_[i] = atoi(arg[iarg_++]);
    if(newtonflag_[i] != 0 && newtonflag_[i] != 1)
        error->fix_error(FLERR,this,"newtonflag must be either 0 or 1");

  }
}

/* ---------------------------------------------------------------------- */

FixContactHistory::~FixContactHistory()
{
  // unregister this fix so atom class doesn't invoke it any more

  atom->delete_callback(id,0);
  atom->delete_callback(id,1);

  // delete locally stored arrays

  memory->destroy(npartner_);
  memory->sfree(partner_);
  memory->sfree(contacthistory_);
  if(ipage_) delete [] ipage_;
  if(dpage_) delete [] dpage_;

  if(variablename_) delete [] variablename_;
  if(newtonflag_) delete [] newtonflag_;

  if(history_id_)
  {
      for(int i = 0; i < dnum_; i++)
          delete [] (history_id_[i]);
      memory->sfree(history_id_);
  }
}

/* ---------------------------------------------------------------------- */

int FixContactHistory::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  mask |= MIN_PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixContactHistory::init()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,
               "Contact history requires atoms have IDs");

  // if (-1 == index_decide_noncontacting && 1. < neighbor->contactHistoryDistanceFactor)
  //    error->fix_error(FLERR,this,"have to call set_index_decide_noncontacting() function");

  if(0 == strcmp(style,"contacthistory"))
  {
      
      pair_gran_ = static_cast<PairGran*>(force->pair_match("gran", 1));

      if(pair_gran_ && dnum_!=static_cast<PairGran*>(pair_gran_)->dnum())
      {
          printf("WARNING: FixContactHistory: Found a PairGran 'gran', but dnum (=%d) and pair_gran->dnum (=%d) is NOT consistent! \n", 
		         dnum_, static_cast<PairGran*>(pair_gran_)->dnum());

          if(!pair_gran_ || dnum_ != static_cast<PairGran*>(pair_gran_)->dnum_pair())
              pair_gran_ = static_cast<PairGran*>(force->pair_match("gran_bubble", 1));
          if(!pair_gran_ || dnum_ != static_cast<PairGran*>(pair_gran_)->dnum_pair())
              pair_gran_ = static_cast<PairGran*>(force->pair_match("bubble", 1));

          if(!pair_gran_)
              error->fix_error(FLERR,this,"Please use a pair style 'gran', 'gran_bubble' or 'bubble' for fix contacthistory, or check your settings for dnum (they may be inconsistent)");

          if(dnum_ != static_cast<PairGran*>(pair_gran_)->dnum_pair())
          {
              printf("FixContactHistory: PairGran 'gran' dnum (=%d) and pair_gran->dnum (=%d) is NOT consistent! \n",
	             dnum_,  static_cast<PairGran*>(pair_gran_)->dnum());
              error->fix_error(FLERR,this,"internal error");
          }
      }
      int dim;
      computeflag_ = (int *) pair_gran_->extract("computeflag",dim);
  }

  allocate_pages();
}

/* ----------------------------------------------------------------------
  create pages if first time or if neighbor pgsize/oneatom has changed
  note that latter could cause shear history info to be discarded
------------------------------------------------------------------------- */

void FixContactHistory::allocate_pages()
{
  int create = 0;
  if (ipage_ == NULL) create = 1;
  if (pgsize_ != neighbor->pgsize) create = 1;
  if (oneatom_ != neighbor->oneatom) create = 1;

  if (create) {
    delete [] ipage_;
    delete [] dpage_;

    pgsize_ = neighbor->pgsize;
    oneatom_ = neighbor->oneatom;
    int nmypage = comm->nthreads;
    ipage_ = new MyPage<int>[nmypage];
    dpage_ = new MyPage<double>[nmypage];
    for (int i = 0; i < nmypage; i++) {
      ipage_[i].init(oneatom_,pgsize_);
      dpage_[i].init(oneatom_*std::max(1,dnum_),pgsize_);
    }
  }
}

/* ----------------------------------------------------------------------
   called by setup of run or minimize
   called by write_restart or write_data as input script command
   only invoke pre_exchange() if neigh list stores more current history info
     than npartner/partner arrays in this fix
   that will only be case if pair->compute() has been invoked since
     update of npartner/npartner
   this logic avoids 2 problems:
     run 100; write_restart; run 100
       setup_pre_exchange is called twice (by write_restart and 2nd run setup)
       w/out a neighbor list being created in between
     read_restart; run 100
       setup_pre_exchange called by run setup whacks restart shear history info
------------------------------------------------------------------------- */

void FixContactHistory::setup_pre_exchange()
{
  if (*computeflag_)
  {
      
      pre_exchange();
  }
  // computeflag is re-set in pre_exchange()
}

/* ----------------------------------------------------------------------
   copy contacthistory partner info from neighbor lists to atom arrays
   so can be migrated or stored with atoms
------------------------------------------------------------------------- */

void FixContactHistory::pre_exchange()
{
  int i,j,ii,jj,m,n,inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *contact_flag,**first_contact_flag;
  double *hist,*allhist,**firsthist;

  // re-set computeflag: most current info is now in atom arrays
  *computeflag_ = 0;

  // nlocal may include atoms added since last neigh build

  int nmax = atom->nmax;

  // zero npartner for all current atoms
  // clear 2 page data structures

  std::fill_n(npartner_, nmax, 0);

  ipage_->reset();
  dpage_->reset();

  // 1st loop over neighbor list
  // calculate npartner for each owned atom
  // nlocal_neigh = nlocal when neigh list was built, may be smaller than nlocal

  int *tag = atom->tag;
  NeighList *list = pair_gran_->list;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  first_contact_flag = list->listgranhistory->firstneigh;
  firsthist = list->listgranhistory->firstdouble;

  int nlocal_neigh = 0;
  if (inum) nlocal_neigh = ilist[inum-1] + 1;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    
    contact_flag = first_contact_flag[i];

    for (jj = 0; jj < jnum; jj++) {
      if (contact_flag[jj]) {
        
        npartner_[i]++;
        j = jlist[jj];
        j &= NEIGHMASK;
        if (j < nlocal_neigh) npartner_[j]++;
      }
    }
  }

  // get page chunks to store atom IDs and shear history for my atoms
  
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    
    n = npartner_[i];
    partner_[i] = ipage_->get(n);
    contacthistory_[i] = dpage_->get(dnum_*n);
    
    if (partner_[i] == NULL || contacthistory_[i] == NULL)
      error->one(FLERR,"Contact history overflow, boost neigh_modify one");
  }

  // 2nd loop over neighbor list
  // store atom IDs and shear history for my atoms
  // re-zero npartner to use as counter for all my atoms

  std::fill_n(npartner_, nmax, 0);

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    allhist = firsthist[i];
    jnum = numneigh[i];
    contact_flag = first_contact_flag[i];

    for (jj = 0; jj < jnum; jj++) {
      if (contact_flag[jj]) {
        hist = &allhist[dnum_*jj];
        j = jlist[jj];
        j &= NEIGHMASK;
        m = npartner_[i];
        partner_[i][m] = tag[j];
        
        for (int d = 0; d < dnum_; d++) {
          contacthistory_[i][m*dnum_+d] = hist[d];
        }
        
        npartner_[i]++;
        if (j < nlocal_neigh) {
          
          m = npartner_[j];

          partner_[j][m] = tag[i];
          for (int d = 0; d < dnum_; d++) {
            if(newtonflag_[d])
              contacthistory_[j][m*dnum_+d] = -hist[d];
            else
              contacthistory_[j][m*dnum_+d] =  hist[d];
          }
          
          npartner_[j]++;
        }
      }
    }
  }

  // set maxtouch = max # of partners of any owned atom
  // bump up comm->maxexchange_fix if necessary
  maxtouch_ = 0;
  int nlocal = atom->nlocal;
  if(nlocal > 0) maxtouch_ = *std::max_element(npartner_, npartner_+nlocal);

  comm->maxexchange_fix = MAX(comm->maxexchange_fix,(dnum_+1)*maxtouch_+1);
}

/* ---------------------------------------------------------------------- */

void FixContactHistory::min_setup_pre_exchange()
{
  if (*computeflag_)
    pre_exchange();
  // computeflag is re-set in pre_exchange()
}

/* ---------------------------------------------------------------------- */

void FixContactHistory::min_pre_exchange()
{
  pre_exchange();
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixContactHistory::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax * sizeof(int);
  bytes += nmax * sizeof(int *);
  bytes += nmax * sizeof(double *);

  int nmypage = comm->nthreads;
  for (int i = 0; i < nmypage; i++) {
    bytes += ipage_[i].size();
    bytes += dpage_[i].size();
  }

  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixContactHistory::grow_arrays(int nmax)
{
  memory->grow(npartner_,nmax,"contact_history:npartner");
  partner_ = (int **) memory->srealloc(partner_,nmax*sizeof(int *),
                                      "contact_history:partner");
  typedef double (*sptype);
  contacthistory_ = (sptype *)
    memory->srealloc(contacthistory_,nmax*sizeof(sptype),
                     "contact_history:shearpartner");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixContactHistory::copy_arrays(int i, int j, int delflag)
{
  // just copy pointers for partner and shearpartner
  // b/c can't overwrite chunk allocation inside ipage,dpage
  // incoming atoms in unpack_exchange just grab new chunks
  // so are orphaning chunks for migrating atoms
  // OK, b/c will reset ipage,dpage on next reneighboring

  npartner_[j] = npartner_[i];
  partner_[j] = partner_[i];
  contacthistory_[j] = contacthistory_[i];
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixContactHistory::set_arrays(int i)
{
  npartner_[i] = 0;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixContactHistory::pack_exchange(int i, double *buf)
{
  // NOTE: how do I know comm buf is big enough if extreme # of touching neighs
  // Comm::BUFEXTRA may need to be increased

  int m = 0;
  buf[m++] = ubuf(npartner_[i]).d;
  for (int n = 0; n < npartner_[i]; n++) {
    
    buf[m++] = ubuf(partner_[i][n]).d;
    for (int d = 0; d < dnum_; d++) {
      buf[m++] = contacthistory_[i][n*dnum_+d];
    }
  }
  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixContactHistory::unpack_exchange(int nlocal, double *buf)
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
  return m;
}

/* ----------------------------------------------------------------------
   pack  state of Fix into one write
------------------------------------------------------------------------- */

void FixContactHistory::write_restart(FILE *fp)
{
  int n = 0;
  double list[6];
  list[n++] = static_cast<double>(dnum_);
  list[n++] = static_cast<double>(maxtouch_);

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixContactHistory::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  int unpack_dnum = static_cast<int> (list[n++]);
  //int unpack_maxtouch = static_cast<int> (list[n++]);

  if(unpack_dnum != dnum_)
    error->fix_error(FLERR,this,"saved simulation state used different contact history model - can not restart");
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixContactHistory::pack_restart(int i, double *buf)
{
  int m = 0;
  
  buf[m++] = (dnum_+1)*npartner_[i] + 2;
  buf[m++] = ubuf(npartner_[i]).d;
  for (int n = 0; n < npartner_[i]; n++) {
    buf[m++] = ubuf(partner_[i][n]).d;
    for (int d = 0; d < dnum_; d++) {
      buf[m++] = contacthistory_[i][n*dnum_+d];
      
    }
  }
  return m;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixContactHistory::unpack_restart(int nlocal, int nth)
{
  // ipage = NULL if being called from granular pair style init()

  if (ipage_ == NULL) allocate_pages();

  // skip to Nth set of extra values

  double **extra = atom->extra;

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  // allocate new chunks from ipage,dpage for incoming values
  npartner_[nlocal] = ubuf(extra[nlocal][m++]).i;
  maxtouch_ = MAX(maxtouch_,npartner_[nlocal]);
  partner_[nlocal] = ipage_->get(npartner_[nlocal]);
  contacthistory_[nlocal] = dpage_->get(npartner_[nlocal]*dnum_);
  if (partner_[nlocal] == NULL || contacthistory_[nlocal] == NULL)
      error->one(FLERR,"Contact history overflow, boost neigh_modify one");

  for (int n = 0; n < npartner_[nlocal]; n++) {
    partner_[nlocal][n] = ubuf(extra[nlocal][m++]).i;
    
    for (int d = 0; d < dnum_; d++) {
      contacthistory_[nlocal][n*dnum_+d] = extra[nlocal][m++];
    }
  }
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixContactHistory::maxsize_restart()
{
  // maxtouch_all = max # of touching partners across all procs

  int maxtouch_all;
  MPI_Allreduce(&maxtouch_,&maxtouch_all,1,MPI_INT,MPI_MAX,world);
  return (dnum_+1)*maxtouch_all + 2;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixContactHistory::size_restart(int nlocal)
{
  return (dnum_+1)*npartner_[nlocal] + 2;
}
