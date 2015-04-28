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

#include "neighbor.h"
#include "neigh_list.h"
#include "atom.h"
#include "group.h"
#include "update.h"
#include "fix_contact_history.h" 
#include "error.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   granular particles
   N^2 / 2 search for neighbor pairs with partial Newton's 3rd law
   shear history must be accounted for when a neighbor pair is added
   pair added to list if atoms i and j are both owned and i < j
   pair added if j is ghost (also stored by proc owning j)
------------------------------------------------------------------------- */

void Neighbor::granular_nsq_no_newton(NeighList *list)
{
  int i,j,m,n,nn=0,bitmask=0,d; 
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double radi,radsum,cutsq;
  int *neighptr,*contact_flag_ptr = NULL;
  double *contact_hist_ptr = NULL;

  NeighList *listgranhistory;
  int *npartner = NULL,**partner = NULL;
  double **contacthistory = NULL; 
  int **first_contact_flag;
  double **first_contact_hist;
  MyPage<int> *ipage_contact_flag = NULL;
  MyPage<double> *dpage_contact_hist = NULL;
  int dnum = 0; 

  double **x = atom->x;
  double *radius = atom->radius;
  int *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  if (includegroup) {
    nlocal = atom->nfirst;
    bitmask = group->bitmask[includegroup];
  }

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;

  FixContactHistory *fix_history = list->fix_history; 
  if (fix_history) {
    npartner = fix_history->npartner_; 
    partner = fix_history->partner_; 
    contacthistory = fix_history->contacthistory_; 
    listgranhistory = list->listgranhistory;
    first_contact_flag = listgranhistory->firstneigh;
    first_contact_hist = listgranhistory->firstdouble;
    ipage_contact_flag = listgranhistory->ipage;
    dpage_contact_hist = listgranhistory->dpage;
    dnum = listgranhistory->dnum; 
  }

  int inum = 0;
  ipage->reset();
  if (fix_history) {
    ipage_contact_flag->reset();
    dpage_contact_hist->reset();
  }

  for (i = 0; i < nlocal; i++) {
    n = 0;
    neighptr = ipage->vget();
    if (fix_history) {
      nn = 0;
      contact_flag_ptr  = ipage_contact_flag->vget();
      contact_hist_ptr = dpage_contact_hist->vget();
    }

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];

    // loop over remaining atoms, owned and ghost

    for (j = i+1; j < nall; j++) {
      if (includegroup && !(mask[j] & bitmask)) continue;
      if (exclude && exclusion(i,j,type[i],type[j],mask,molecule)) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radsum = (radi + radius[j]) * contactDistanceFactor; 
      cutsq = (radsum+skin) * (radsum+skin);

      if (rsq <= cutsq) {
        neighptr[n] = j;

        if (fix_history) {
          if (rsq < radsum*radsum)
          {
            for (m = 0; m < npartner[i]; m++)
              if (partner[i][m] == tag[j]) break;
            if (m < npartner[i]) {
              contact_flag_ptr[n] = 1;
              for (d = 0; d < dnum; d++) {  
                contact_hist_ptr[nn++] = contacthistory[i][m*dnum+d];
              }
            } else {
              contact_flag_ptr[n] = 0;
              for (d = 0; d < dnum; d++) {  
                contact_hist_ptr[nn++] = 0.0;
              }
            }
          } else {
            contact_flag_ptr[n] = 0;
            for (d = 0; d < dnum; d++) {  
              contact_hist_ptr[nn++] = 0.0;
            }
          }
        }

        n++;
      }
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
    if (fix_history) {
      first_contact_flag[i] = contact_flag_ptr;
      first_contact_hist[i] = contact_hist_ptr;
      ipage_contact_flag->vgot(n);
      dpage_contact_hist->vgot(nn);
    }
  }

  list->inum = inum;
}

/* ----------------------------------------------------------------------
   granular particles
   N^2 / 2 search for neighbor pairs with full Newton's 3rd law
   no shear history is allowed for this option
   pair added to list if atoms i and j are both owned and i < j
   if j is ghost only me or other proc adds pair
   decision based on itag,jtag tests
------------------------------------------------------------------------- */

void Neighbor::granular_nsq_newton(NeighList *list)
{
  int i,j,n,itag,jtag,bitmask=0;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double radi,radsum,cutsq;
  int *neighptr;

  double **x = atom->x;
  double *radius = atom->radius;
  int *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  if (includegroup) {
    nlocal = atom->nfirst;
    bitmask = group->bitmask[includegroup];
  }

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;

  int inum = 0;
  ipage->reset();

  for (i = 0; i < nlocal; i++) {

    n = 0;
    neighptr = ipage->vget();

    itag = tag[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];

    // loop over remaining atoms, owned and ghost

    for (j = i+1; j < nall; j++) {
      if (includegroup && !(mask[j] & bitmask)) continue;

      if (j >= nlocal) {
        jtag = tag[j];
        if (itag > jtag) {
          if ((itag+jtag) % 2 == 0) continue;
        } else if (itag < jtag) {
          if ((itag+jtag) % 2 == 1) continue;
        } else {
          if (x[j][2] < ztmp) continue;
          if (x[j][2] == ztmp) {
            if (x[j][1] < ytmp) continue;
            if (x[j][1] == ytmp && x[j][0] < xtmp) continue;
          }
        }
      }

      if (exclude && exclusion(i,j,type[i],type[j],mask,molecule)) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radsum = (radi + radius[j]) * contactDistanceFactor;
      cutsq = (radsum+skin) * (radsum+skin);

      if (rsq <= cutsq) neighptr[n++] = j;
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }

  list->inum = inum;
}

/* ----------------------------------------------------------------------
   granular particles
   binned neighbor list construction with partial Newton's 3rd law
   include neighbors of ghost atoms, but no "special neighbors" for ghosts
   shear history must be accounted for when a neighbor pair is added
   owned and ghost atoms check own bin and other bins in stencil
   pair stored once if i,j are both owned and i < j
   pair stored by me if j is ghost (also stored by proc owning j)
   pair stored once if i,j are both ghost and i < j
   no contact history for neighbors of ghosts
------------------------------------------------------------------------- */

void Neighbor::granular_bin_no_newton_ghost(NeighList *list)
{
  int i,j,k,m,n,nn=0,ibin,d;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int xbin,ybin,zbin,xbin2,ybin2,zbin2;
  double radi,radsum,cutsq;
  int *neighptr,*contact_flag_ptr = NULL;
  double *contact_hist_ptr = NULL;

  NeighList *listgranhistory;
  int *npartner = NULL,**partner = NULL;
  double **contacthistory = NULL;
  int **first_contact_flag = NULL;
  double **first_contact_hist = NULL;
  MyPage<int> *ipage_contact_flag = NULL;
  MyPage<double> *dpage_contact_hist = NULL;
  int dnum = 0; 

  // bin local & ghost atoms

  bin_atoms();

  // loop over each atom, storing neighbors

  double **x = atom->x;
  double *radius = atom->radius;
  int *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int nstencil = list->nstencil;
  int *stencil = list->stencil;
  int **stencilxyz = list->stencilxyz;
  MyPage<int> *ipage = list->ipage;

  FixContactHistory *fix_history = list->fix_history; 
  if (fix_history) {
    npartner = fix_history->npartner_; 
    partner = fix_history->partner_; 
    contacthistory = fix_history->contacthistory_; 
    listgranhistory = list->listgranhistory;
    first_contact_flag = listgranhistory->firstneigh;
    first_contact_hist = listgranhistory->firstdouble;
    ipage_contact_flag = listgranhistory->ipage;
    dpage_contact_hist = listgranhistory->dpage;
    dnum = listgranhistory->dnum; 
  }

  int inum = 0;
  ipage->reset();
        if (fix_history) {
    ipage_contact_flag->reset();
    dpage_contact_hist->reset();
    }

  for (i = 0; i < nall; i++) {
    n = 0;
    neighptr = ipage->vget();
    if (fix_history) {
      nn = 0;
      contact_flag_ptr = ipage_contact_flag->vget();
      contact_hist_ptr = dpage_contact_hist->vget();
    }

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];

    // loop over all atoms in surrounding bins in stencil including self
    // only store pair if i < j
    // stores own/own pairs only once
    // stores own/ghost pairs on both procs
    // stores ghost/ghost pairs only once

    if (i < nlocal) {
      ibin = coord2bin(x[i]);

      for (k = 0; k < nstencil; k++) {
        for (j = binhead[ibin+stencil[k]]; j >= 0; j = bins[j]) {
          if (j <= i) continue;
          if (exclude && exclusion(i,j,type[i],type[j],mask,molecule)) continue;

          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx*delx + dely*dely + delz*delz;
          radsum = (radi + radius[j]) * contactDistanceFactor; 
          cutsq = (radsum+skin) * (radsum+skin);
          
          if (rsq <= cutsq) {
            neighptr[n] = j;
            
            if (fix_history) {
              if (rsq < radsum*radsum)
              {
                for (m = 0; m < npartner[i]; m++)
                  if (partner[i][m] == tag[j]) break;
                if (m < npartner[i]) {
                  contact_flag_ptr[n] = 1;
                  for (d = 0; d < dnum; d++) { 
                    contact_hist_ptr[nn++] = contacthistory[i][m*dnum+d];
                  }
                } else {
                   contact_flag_ptr[n] = 0;
                   for (d = 0; d < dnum; d++) { 
                     contact_hist_ptr[nn++] = 0.0;
                   }
                }
              } else {
                contact_flag_ptr[n] = 0;
                for (d = 0; d < dnum; d++) { 
                    contact_hist_ptr[nn++] = 0.0;
                }
              }
            }

            n++;
          }
        }
      }
    } else {
      ibin = coord2bin(x[i],xbin,ybin,zbin);
      for (k = 0; k < nstencil; k++) {
        xbin2 = xbin + stencilxyz[k][0];
        ybin2 = ybin + stencilxyz[k][1];
        zbin2 = zbin + stencilxyz[k][2];
        if (xbin2 < 0 || xbin2 >= mbinx ||
            ybin2 < 0 || ybin2 >= mbiny ||
            zbin2 < 0 || zbin2 >= mbinz) continue;
        for (j = binhead[ibin+stencil[k]]; j >= 0; j = bins[j]) {
          if (j <= i) continue;

          if (exclude && exclusion(i,j,type[i],type[j],mask,molecule)) continue;

          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx*delx + dely*dely + delz*delz;
          radsum = (radi + radius[j]) * contactDistanceFactor; 
          cutsq = (radsum+skin) * (radsum+skin);

          if (rsq <= cutsq) neighptr[n++] = j;
        }
      }
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
    if (fix_history && i < nlocal) {
      first_contact_flag[i] = contact_flag_ptr;
      first_contact_hist[i] = contact_hist_ptr;
      ipage_contact_flag->vgot(n);
      dpage_contact_hist->vgot(nn);
    } else {
      first_contact_flag[i] = 0;
      first_contact_hist[i] = 0;
    }
  }

  list->inum = inum;
}

/* ----------------------------------------------------------------------
   granular particles
   binned neighbor list construction with partial Newton's 3rd law
   shear history must be accounted for when a neighbor pair is added
   each owned atom i checks own bin and surrounding bins in non-Newton stencil
   pair stored once if i,j are both owned and i < j
   pair stored by me if j is ghost (also stored by proc owning j)
------------------------------------------------------------------------- */

void Neighbor::granular_bin_no_newton(NeighList *list)
{
  int i,j,k,m,n,nn=0,ibin,d;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double radi,radsum,cutsq;
  int *neighptr,*contact_flag_ptr = NULL;
  double *contact_hist_ptr = NULL;

  NeighList *listgranhistory;
  int *npartner = NULL,**partner = NULL;
  double **contacthistory = NULL;
  int **first_contact_flag = NULL;
  double **first_contact_hist = NULL;
  MyPage<int> *ipage_contact_flag = NULL;
  MyPage<double> *dpage_contact_hist = NULL;
  int dnum = 0; 

  // bin local & ghost atoms

  bin_atoms();

  // loop over each atom, storing neighbors

  double **x = atom->x;
  double *radius = atom->radius;
  int *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  if (includegroup) nlocal = atom->nfirst;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int nstencil = list->nstencil;
  int *stencil = list->stencil;
  MyPage<int> *ipage = list->ipage;

  FixContactHistory *fix_history = list->fix_history; 
  if (fix_history) {
    npartner = fix_history->npartner_; 
    partner = fix_history->partner_; 
    contacthistory = fix_history->contacthistory_; 
    listgranhistory = list->listgranhistory;
    first_contact_flag = listgranhistory->firstneigh;
    first_contact_hist = listgranhistory->firstdouble;
    ipage_contact_flag = listgranhistory->ipage;
    dpage_contact_hist = listgranhistory->dpage;
    dnum = listgranhistory->dnum; 
  }

  int inum = 0;
  ipage->reset();
        if (fix_history) {
    ipage_contact_flag->reset();
    dpage_contact_hist->reset();
    }

  for (i = 0; i < nlocal; i++) {
    n = 0;
    neighptr = ipage->vget();
    if (fix_history) {
      nn = 0;
      contact_flag_ptr = ipage_contact_flag->vget();
      contact_hist_ptr = dpage_contact_hist->vget();
    }

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    ibin = coord2bin(x[i]);

    // loop over all atoms in surrounding bins in stencil including self
    // only store pair if i < j
    // stores own/own pairs only once
    // stores own/ghost pairs on both procs

    for (k = 0; k < nstencil; k++) {
      for (j = binhead[ibin+stencil[k]]; j >= 0; j = bins[j]) {
        if (j <= i) continue;
        if (exclude && exclusion(i,j,type[i],type[j],mask,molecule)) continue;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        radsum = (radi + radius[j]) * contactDistanceFactor; 
        cutsq = (radsum+skin) * (radsum+skin);
        
        if (rsq <= cutsq) {
          neighptr[n] = j;
          
          if (fix_history) {
            
            if (rsq < radsum*radsum)
            {
              for (m = 0; m < npartner[i]; m++)
                if (partner[i][m] == tag[j]) break;
              if (m < npartner[i]) {
                contact_flag_ptr[n] = 1;
                
                for (d = 0; d < dnum; d++) { 
                  contact_hist_ptr[nn++] = contacthistory[i][m*dnum+d];
                }
              } else {
                 
                 contact_flag_ptr[n] = 0;
                 for (d = 0; d < dnum; d++) { 
                   contact_hist_ptr[nn++] = 0.0;
                 }
              }
            }
            
            else
            {
              
              contact_flag_ptr[n] = 0;
              for (d = 0; d < dnum; d++) { 
                contact_hist_ptr[nn++] = 0.0;
              }
            }
          }

          n++;
        }
      }
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
    if (fix_history) {
      first_contact_flag[i] = contact_flag_ptr;
      first_contact_hist[i] = contact_hist_ptr;
      ipage_contact_flag->vgot(n);
      dpage_contact_hist->vgot(nn);
    }
  }

  list->inum = inum;
}

/* ----------------------------------------------------------------------
   granular particles
   binned neighbor list construction with full Newton's 3rd law
   no shear history is allowed for this option
   each owned atom i checks its own bin and other bins in Newton stencil
   every pair stored exactly once by some processor
------------------------------------------------------------------------- */

void Neighbor::granular_bin_newton(NeighList *list)
{
  int i,j,k,n,ibin;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double radi,radsum,cutsq;
  int *neighptr;

  // bin local & ghost atoms

  bin_atoms();

  // loop over each atom, storing neighbors

  double **x = atom->x;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  if (includegroup) nlocal = atom->nfirst;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int nstencil = list->nstencil;
  int *stencil = list->stencil;
  MyPage<int> *ipage = list->ipage;

  int inum = 0;
  ipage->reset();

  for (i = 0; i < nlocal; i++) {

    n = 0;
    neighptr = ipage->vget();

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];

    // loop over rest of atoms in i's bin, ghosts are at end of linked list
    // if j is owned atom, store it, since j is beyond i in linked list
    // if j is ghost, only store if j coords are "above and to the right" of i

    for (j = bins[i]; j >= 0; j = bins[j]) {
      if (j >= nlocal) {
        if (x[j][2] < ztmp) continue;
        if (x[j][2] == ztmp) {
          if (x[j][1] < ytmp) continue;
          if (x[j][1] == ytmp && x[j][0] < xtmp) continue;
        }
      }

      if (exclude && exclusion(i,j,type[i],type[j],mask,molecule)) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radsum = (radi + radius[j]) * contactDistanceFactor; 
      cutsq = (radsum+skin) * (radsum+skin);

      if (rsq <= cutsq) neighptr[n++] = j;
    }

    // loop over all atoms in other bins in stencil, store every pair

    ibin = coord2bin(x[i]);
    for (k = 0; k < nstencil; k++) {
      for (j = binhead[ibin+stencil[k]]; j >= 0; j = bins[j]) {
        if (exclude && exclusion(i,j,type[i],type[j],mask,molecule)) continue;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        radsum = radi + radius[j];
        cutsq = (radsum+skin) * (radsum+skin);

        if (rsq <= cutsq) neighptr[n++] = j;
      }
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }

  list->inum = inum;
}

/* ----------------------------------------------------------------------
   granular particles
   binned neighbor list construction with Newton's 3rd law for triclinic
   no shear history is allowed for this option
   each owned atom i checks its own bin and other bins in triclinic stencil
   every pair stored exactly once by some processor
------------------------------------------------------------------------- */

void Neighbor::granular_bin_newton_tri(NeighList *list)
{
  int i,j,k,n,ibin;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double radi,radsum,cutsq;
  int *neighptr;

  // bin local & ghost atoms

  bin_atoms();

  // loop over each atom, storing neighbors

  double **x = atom->x;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  if (includegroup) nlocal = atom->nfirst;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int nstencil = list->nstencil;
  int *stencil = list->stencil;
  MyPage<int> *ipage = list->ipage;

  int inum = 0;
  ipage->reset();

  for (i = 0; i < nlocal; i++) {

    n = 0;
    neighptr = ipage->vget();

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];

    // loop over all atoms in bins in stencil
    // pairs for atoms j "below" i are excluded
    // below = lower z or (equal z and lower y) or (equal zy and lower x)
    //         (equal zyx and j <= i)
    // latter excludes self-self interaction but allows superposed atoms

    ibin = coord2bin(x[i]);
    for (k = 0; k < nstencil; k++) {
      for (j = binhead[ibin+stencil[k]]; j >= 0; j = bins[j]) {
        if (x[j][2] < ztmp) continue;
        if (x[j][2] == ztmp) {
          if (x[j][1] < ytmp) continue;
          if (x[j][1] == ytmp) {
            if (x[j][0] < xtmp) continue;
            if (x[j][0] == xtmp && j <= i) continue;
          }
        }

        if (exclude && exclusion(i,j,type[i],type[j],mask,molecule)) continue;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        radsum = (radi + radius[j]) * contactDistanceFactor; 
        cutsq = (radsum+skin) * (radsum+skin);

        if (rsq <= cutsq) neighptr[n++] = j;
      }
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }

  list->inum = inum;
}
