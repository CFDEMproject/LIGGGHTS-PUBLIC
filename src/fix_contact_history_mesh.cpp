/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "string.h"
#include "stdio.h"
#include "fix_contact_history_mesh.h"
#include "atom.h"
#include "fix_mesh_surface.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "force.h"
#include "pair.h"
#include "update.h"
#include "modify.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixContactHistoryMesh::FixContactHistoryMesh(LAMMPS *lmp, int narg, char **arg) :
  FixContactHistory(lmp, narg, arg),
  ipage1_(0),
  dpage1_(0),
  ipage2_(0),
  dpage2_(0),
  delpage_(0),
  delflag_(0),
  mesh_(0),
  fix_neighlist_mesh_(0),
  fix_nneighs_(0),
  build_neighlist_(true)
{
  
  // parse args
  Fix *f = modify->find_fix_id(arg[iarg_++]);
  if(!f || strncmp(f->style,"mesh/surface",12) )
    error->fix_error(FLERR,this,"wrong ID for fix mesh/surface");
  mesh_ = (static_cast<FixMeshSurface*>(f))->triMesh();
  fix_neighlist_mesh_ = (static_cast<FixMeshSurface*>(f))->meshNeighlist();

  swap_ = new double[dnum_];

  // initial allocation of delflag
  delflag_ = (bool **) memory->srealloc(delflag_,atom->nmax*sizeof(bool *),
                                      "contact_history:delflag");
}

/* ---------------------------------------------------------------------- */

FixContactHistoryMesh::~FixContactHistoryMesh()
{
  // delete locally stored arrays

  if(ipage1_) delete [] ipage1_;
  if(dpage1_) delete [] dpage1_;
  if(ipage2_) delete [] ipage2_;
  if(dpage2_) delete [] dpage2_;
  if(delpage_) delete [] delpage_;

  ipage_ = 0;
  dpage_ = 0;

  delete [] swap_;

  if(delflag_) memory->sfree(delflag_);
}

/* ---------------------------------------------------------------------- */

int FixContactHistoryMesh::setmask()
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

void FixContactHistoryMesh::init()
{
  FixContactHistory::init();

  char *fix_nneighs_name = new char[strlen(mesh_->mesh_id())+1+14];
  sprintf(fix_nneighs_name,"n_neighs_mesh_%s",mesh_->mesh_id());

  fix_nneighs_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(fix_nneighs_name,"property/atom","scalar",0,0,this->style));

  delete [] fix_nneighs_name;
}

/* ----------------------------------------------------------------------
  create pages if first time or if neighbor pgsize/oneatom has changed
  note that latter could cause shear history info to be discarded
------------------------------------------------------------------------- */

void FixContactHistoryMesh::allocate_pages()
{
  if ((ipage_ == NULL) || (pgsize_ != neighbor->pgsize) || (oneatom_ != neighbor->oneatom)) {

    bool use_first = ipage_ == ipage2_;

    delete [] ipage1_;
    delete [] dpage1_;
    delete [] ipage2_;
    delete [] dpage2_;
    delete [] delpage_;

    pgsize_ = neighbor->pgsize;
    oneatom_ = neighbor->oneatom;
    int nmypage = comm->nthreads;
    ipage1_ = new MyPage<int>[nmypage];
    dpage1_ = new MyPage<double>[nmypage];
    ipage2_ = new MyPage<int>[nmypage];
    dpage2_ = new MyPage<double>[nmypage];
    delpage_ = new MyPage<bool>[nmypage];
    for (int i = 0; i < nmypage; i++) {
      ipage1_[i].init(oneatom_,pgsize_);
      dpage1_[i].init(oneatom_,pgsize_);
      ipage2_[i].init(oneatom_,pgsize_);
      dpage2_[i].init(oneatom_,pgsize_);
      delpage_[i].init(oneatom_,pgsize_);
    }

    if(use_first)
    {
        ipage_ = ipage1_;
        dpage_ = dpage1_;
    }
    else
    {
        ipage_ = ipage2_;
        dpage_ = dpage2_;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixContactHistoryMesh::setup_pre_exchange()
{
    pre_exchange();
}
/* ---------------------------------------------------------------------- */

void FixContactHistoryMesh::min_setup_pre_exchange()
{
  pre_exchange();
}

/* ----------------------------------------------------------------------
   sort contacthistory so need
------------------------------------------------------------------------- */

void FixContactHistoryMesh::pre_exchange()
{
    
    if(!recent_restart)
        sort_contacts();

   // set maxtouch = max # of partners of any owned atom
   // bump up comm->maxexchange_fix if necessary

   int nlocal = atom->nlocal;

   maxtouch_ = 0;
   for (int i = 0; i < nlocal; i++) maxtouch_ = MAX(maxtouch_,npartner_[i]);
   comm->maxexchange_fix = MAX(comm->maxexchange_fix,4*maxtouch_+1);
}

/* ----------------------------------------------------------------------
   need to execute here since neighlist is refreshed in setup_pre_force(int foo)
   if this is not done, data will get out-of-sync

   need not care about overwriting hist values since always have pointer
   to most current data
------------------------------------------------------------------------- */

void FixContactHistoryMesh::setup_pre_force(int dummy)
{
   pre_neighbor();
   pre_force(0);
}

/* ---------------------------------------------------------------------- */

void FixContactHistoryMesh::pre_neighbor()
{
    build_neighlist_ = true;
    
}

/* ----------------------------------------------------------------------
   swap 2 pages data structures where contacthistory is stored
   do this here because exchange and sort can rag the data structure
------------------------------------------------------------------------- */

void FixContactHistoryMesh::pre_force(int dummy)
{
    if(!build_neighlist_)
        return;
    build_neighlist_ = false;

    cleanUpContactJumps();

    int nneighs_next;
    int nlocal = atom->nlocal;
    int *partner_prev;
    double *contacthistory_prev;
    double *nn = fix_nneighs_->vector_atom;

    MyPage<int>    *ipage_prev = ipage_;
    MyPage<double> *dpage_prev = dpage_;
    MyPage<int>    *ipage_next = (ipage_ == ipage1_) ? ipage2_ : ipage1_;
    MyPage<double> *dpage_next = (dpage_ == dpage1_) ? dpage2_ : dpage1_;

    ipage_next->reset();
    dpage_next->reset();

    for (int i = 0; i < nlocal; i++)
    {
        nneighs_next = static_cast<int>(round(nn[i]));

        partner_prev = partner_[i];
        contacthistory_prev = contacthistory_[i];

        partner_[i] = ipage_next->get(nneighs_next);
        vectorInitializeN(partner_[i],nneighs_next,-1);
        contacthistory_[i] = dpage_next->get(nneighs_next*dnum_);
        vectorZeroizeN(contacthistory_[i],nneighs_next*dnum_);

        if(npartner_[i] > nneighs_next)
        {
            
            error->one(FLERR,"internal error");
        }

        for(int ipartner = 0; ipartner < npartner_[i]; ipartner++)
        {
            if(partner_prev[ipartner] < 0)
                error->one(FLERR,"internal error");

            partner_[i][ipartner] = partner_prev[ipartner];
            vectorCopyN(&(contacthistory_prev[ipartner*dnum_]),&(contacthistory_[i][ipartner*dnum_]),dnum_);
            
        }
   }

    ipage_ = ipage_next;
    dpage_ = dpage_next;
}

/* ---------------------------------------------------------------------- */

void FixContactHistoryMesh::sort_contacts()
{
    int nlocal = atom->nlocal;
    int nneighs, first_empty, last_filled;

    for(int i = 0; i < nlocal; i++)
    {
        nneighs = static_cast<int>(round(fix_nneighs_->vector_atom[i]));

        if(0 == nneighs)
            continue;

        do
        {
            first_empty = last_filled = -1;

            for(int j = 0; j < nneighs; j++)
            {
                if(-1 == first_empty && -1 == partner_[i][j])
                    first_empty = j;
                if(partner_[i][j] >= 0)
                    last_filled = j;
            }

            if(first_empty > -1 && last_filled > -1 && first_empty < last_filled)
                swap(i,first_empty,last_filled,true);
        }
        while(first_empty > -1 && last_filled > -1 && first_empty < last_filled);
    }

}

/* ----------------------------------------------------------------------
     mark all contacts for deletion
------------------------------------------------------------------------- */

void FixContactHistoryMesh::markAllContacts()
{
    int nneighs;
    int nlocal = atom->nlocal;
    double *nn = fix_nneighs_->vector_atom;
    delpage_->reset();

    for(int i = 0; i < nlocal; i++)
    {
      nneighs = static_cast<int>(round(nn[i]));

      delflag_[i] = delpage_->get(nneighs);
      for(int ii = 0; ii < nneighs; ii++)
          delflag_[i][ii] = true;
    }
}

/* ---------------------------------------------------------------------- */

void FixContactHistoryMesh::cleanUpContacts()
{
    int nlocal = atom->nlocal;
    int nneighs;

    for(int i = 0; i < nlocal; i++)
    {
        nneighs = static_cast<int>(round(fix_nneighs_->vector_atom[i]));

        for(int j = 0; j < nneighs; j++)
        {
            // delete values
            
            if(delflag_[i][j])
            {
                if(partner_[i][j] > -1)
                {
                    npartner_[i]--;
                    
                }
                partner_[i][j] = -1;
                vectorZeroizeN(&(contacthistory_[i][j*dnum_]),dnum_);
            }
        }
        
    }
}
/* ---------------------------------------------------------------------- */

void FixContactHistoryMesh::cleanUpContactJumps()
{
    int nlocal = atom->nlocal;
    int iTri;

    for(int i = 0; i < nlocal; i++)
    {
        int ipartner = 0;
        while (ipartner < npartner_[i])
        {
            
            if(partner_[i][ipartner] < 0)
                error->one(FLERR,"internal error");

            iTri = mesh_->map(partner_[i][ipartner]);
            
            if(iTri > -1 && !fix_neighlist_mesh_->contactInList(iTri,i))
            {
                
                partner_[i][ipartner] = -1;
                vectorZeroizeN(&(contacthistory_[i][ipartner*dnum_]),dnum_);
                swap(i,ipartner,npartner_[i]-1,false);
                npartner_[i]--;

            }
            else
                ipartner++;
        }
    }
}

/* ---------------------------------------------------------------------- */

void FixContactHistoryMesh::reset_history()
{
    int nlocal = atom->nlocal;
    int nneighs;

    for(int i = 0; i < nlocal; i++)
    {
        nneighs = static_cast<int>(round(fix_nneighs_->vector_atom[i]));

        // zeroize values
        for(int j = 0; j < nneighs; j++)
            vectorZeroizeN(&(contacthistory_[i][j*dnum_]),dnum_);
    }
}

/* ---------------------------------------------------------------------- */

void FixContactHistoryMesh::min_setup_pre_force(int dummy)
{
  if (*computeflag_) pre_force(0);
  *computeflag_ = 0;
}

/* ---------------------------------------------------------------------- */

void FixContactHistoryMesh::min_pre_force(int dummy)
{
  pre_force(0);
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixContactHistoryMesh::grow_arrays(int nmax)
{
  FixContactHistory::grow_arrays(nmax);
  delflag_ = (bool **) memory->srealloc(delflag_,nmax*sizeof(bool *),
                                      "contact_history:delflag");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixContactHistoryMesh::copy_arrays(int i, int j, int delflag)
{
  // just copy pointers for partner and shearpartner
  // b/c can't overwrite chunk allocation inside ipage,dpage
  // incoming atoms in unpack_exchange just grab new chunks
  // so are orphaning chunks for migrating atoms
  // OK, b/c will reset ipage,dpage on next reneighboring

  FixContactHistory::copy_arrays(i,j,delflag);
  delflag_[j] = delflag_[i];
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixContactHistoryMesh::unpack_exchange(int nlocal, double *buf)
{
  // allocate new chunks from ipage,dpage for incoming values
  
  int m = 0;
  int nneighs = static_cast<int>(round(fix_nneighs_->vector_atom[nlocal]));

  npartner_[nlocal] = static_cast<int> (buf[m++]);
  maxtouch_ = MAX(maxtouch_,npartner_[nlocal]);
  partner_[nlocal] = ipage_->get(nneighs);
  contacthistory_[nlocal] = dpage_->get(dnum_*nneighs);

  for (int n = 0; n < npartner_[nlocal]; n++) {
    partner_[nlocal][n] = static_cast<int> (buf[m++]);
    
    for (int d = 0; d < dnum_; d++) {
      contacthistory_[nlocal][n*dnum_+d] = buf[m++];
      
    }
  }

  for (int n = npartner_[nlocal]; n < nneighs; n++) {
    partner_[nlocal][n] = -1;
    for (int d = 0; d < dnum_; d++) {
      contacthistory_[nlocal][n*dnum_+d] = 0.;
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixContactHistoryMesh::unpack_restart(int nlocal, int nth)
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

  npartner_[nlocal] = static_cast<int> (extra[nlocal][m++]);
  maxtouch_ = MAX(maxtouch_,npartner_[nlocal]);
  partner_[nlocal] = ipage_->get(npartner_[nlocal]);
  contacthistory_[nlocal] = dpage_->get(npartner_[nlocal]*dnum_);

  for (int n = 0; n < npartner_[nlocal]; n++) {
    partner_[nlocal][n] = static_cast<int> (extra[nlocal][m++]);
    for (d = 0; d < dnum_; d++) {
      contacthistory_[nlocal][n*dnum_+d] = extra[nlocal][m++];
      
    }
  }
}

/* ----------------------------------------------------------------------
   pack state of Fix into one write
------------------------------------------------------------------------- */

void FixContactHistoryMesh::write_restart(FILE *fp)
{
    FixContactHistory::write_restart(fp);

    sort_contacts();
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixContactHistoryMesh::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax * sizeof(int);
  bytes += nmax * sizeof(int *);
  bytes += nmax * sizeof(double *);

  int nmypage = comm->nthreads;
  for (int i = 0; i < nmypage; i++) {
    bytes += ipage1_[i].size();
    bytes += dpage1_[i].size();
    bytes += ipage2_[i].size();
    bytes += dpage2_[i].size();
    bytes += delpage_[i].size();
  }

  return bytes;
}
