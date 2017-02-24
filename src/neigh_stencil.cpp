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
    This file is from LAMMPS
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

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   routines to create a stencil = list of bin offsets
   stencil = bins whose closest corner to central bin is within cutoff
   sx,sy,sz = bin bounds = furthest the stencil could possibly extend
   3d creates xyz stencil, 2d creates xy stencil
   for half list with newton off:
     stencil is all surrounding bins including self
     regardless of triclinic
   for half list with newton on:
     stencil is bins to the "upper right" of central bin
     stencil does not include self
   for half list with triclinic:
     stencil is all bins in z-plane of self and above, but not below
     in 2d is all bins in y-plane of self and above, but not below
     stencil includes self
   for full list:
     stencil is all surrounding bins including self
     regardless of newton on/off or triclinic
   for multi:
     create one stencil for each atom type
     stencil is same bin bounds as newton on/off, triclinic, half/full
     cutoff is not cutneighmaxsq, but max cutoff for that atom type
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_bin_2d_no_newton(NeighList *list,
                                             int sx, int sy, int sz)
{
  int i,j;
  int *stencil = list->stencil;
  int nstencil = 0;

  for (j = -sy; j <= sy; j++)
    for (i = -sx; i <= sx; i++)
      if (bin_distance(i,j,0) < cutneighmaxsq)
        stencil[nstencil++] = j*mbinx + i;
  list->nstencil = nstencil;
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_ghost_bin_2d_no_newton(NeighList *list,
                                                   int sx, int sy, int sz)
{
  
  int *stencil = list->stencil;
  int **stencilxyz = list->stencilxyz;
  int nstencil = 0;

  for (int j = -sy; j <= sy; j++) 
    for (int i = -sx; i <= sx; i++) 
      if (bin_distance(i,j,0) < cutneighmaxsq) {
        stencilxyz[nstencil][0] = i;
        stencilxyz[nstencil][1] = j;
        stencilxyz[nstencil][2] = 0;
        stencil[nstencil++] = j*mbinx + i;
      }

  list->nstencil = nstencil;
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_bin_3d_no_newton(NeighList *list,
                                             int sx, int sy, int sz)
{
  int i,j,k;
  int *stencil = list->stencil;
  int nstencil = 0;

  for (k = -sz; k <= sz; k++)
    for (j = -sy; j <= sy; j++)
      for (i = -sx; i <= sx; i++)
        if (bin_distance(i,j,k) < cutneighmaxsq)
          stencil[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
  
  list->nstencil = nstencil;
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_ghost_bin_3d_no_newton(NeighList *list,
                                                   int sx, int sy, int sz)
{
  int i,j,k;
  int *stencil = list->stencil;
  int **stencilxyz = list->stencilxyz;
  int nstencil = 0;

  for (k = -sz; k <= sz; k++)
    for (j = -sy; j <= sy; j++)
      for (i = -sx; i <= sx; i++)
        if (bin_distance(i,j,k) < cutneighmaxsq) {
          stencilxyz[nstencil][0] = i;
          stencilxyz[nstencil][1] = j;
          stencilxyz[nstencil][2] = k;
          stencil[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
        }

  list->nstencil = nstencil;
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_bin_2d_newton(NeighList *list,
                                          int sx, int sy, int sz)
{
  int i,j;
  int *stencil = list->stencil;
  int nstencil = 0;

  for (j = 0; j <= sy; j++)
    for (i = -sx; i <= sx; i++)
      if (j > 0 || (j == 0 && i > 0))
        if (bin_distance(i,j,0) < cutneighmaxsq)
          stencil[nstencil++] = j*mbinx + i;

  list->nstencil = nstencil;
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_bin_3d_newton(NeighList *list,
                                          int sx, int sy, int sz)
{
  int i,j,k;
  int *stencil = list->stencil;
  int nstencil = 0;

  for (k = 0; k <= sz; k++)
    for (j = -sy; j <= sy; j++)
      for (i = -sx; i <= sx; i++)
        if (k > 0 || j > 0 || (j == 0 && i > 0))
          if (bin_distance(i,j,k) < cutneighmaxsq)
            stencil[nstencil++] = k*mbiny*mbinx + j*mbinx + i;

  list->nstencil = nstencil;
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_bin_2d_newton_tri(NeighList *list,
                                              int sx, int sy, int sz)
{
  int i,j;
  int *stencil = list->stencil;
  int nstencil = 0;

  for (j = 0; j <= sy; j++)
    for (i = -sx; i <= sx; i++)
      if (bin_distance(i,j,0) < cutneighmaxsq)
        stencil[nstencil++] = j*mbinx + i;

  list->nstencil = nstencil;
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_bin_3d_newton_tri(NeighList *list,
                                              int sx, int sy, int sz)
{
  int i,j,k;
  int *stencil = list->stencil;
  int nstencil = 0;

  for (k = 0; k <= sz; k++)
    for (j = -sy; j <= sy; j++)
      for (i = -sx; i <= sx; i++)
        if (bin_distance(i,j,k) < cutneighmaxsq)
          stencil[nstencil++] = k*mbiny*mbinx + j*mbinx + i;

  list->nstencil = nstencil;
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_multi_2d_no_newton(NeighList *list,
                                               int sx, int sy, int sz)
{
  int i,j,n;
  double rsq,typesq;
  int *s;
  double *distsq;

  int *nstencil_multi = list->nstencil_multi;
  int **stencil_multi = list->stencil_multi;
  double **distsq_multi = list->distsq_multi;

  int ntypes = atom->ntypes;
  for (int itype = 1; itype <= ntypes; itype++) {
    typesq = cuttypesq[itype];
    s = stencil_multi[itype];
    distsq = distsq_multi[itype];
    n = 0;
    for (j = -sy; j <= sy; j++)
      for (i = -sx; i <= sx; i++) {
        rsq = bin_distance(i,j,0);
        if (rsq < typesq) {
          distsq[n] = rsq;
          s[n++] = j*mbinx + i;
        }
      }
    nstencil_multi[itype] = n;
  }
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_multi_3d_no_newton(NeighList *list,
                                               int sx, int sy, int sz)
{
  int i,j,k,n;
  double rsq,typesq;
  int *s;
  double *distsq;

  int *nstencil_multi = list->nstencil_multi;
  int **stencil_multi = list->stencil_multi;
  double **distsq_multi = list->distsq_multi;

  int ntypes = atom->ntypes;
  for (int itype = 1; itype <= ntypes; itype++) {
    typesq = cuttypesq[itype];
    s = stencil_multi[itype];
    distsq = distsq_multi[itype];
    n = 0;
    for (k = -sz; k <= sz; k++)
      for (j = -sy; j <= sy; j++)
        for (i = -sx; i <= sx; i++) {
          rsq = bin_distance(i,j,k);
          if (rsq < typesq) {
            distsq[n] = rsq;
            s[n++] = k*mbiny*mbinx + j*mbinx + i;
          }
        }
    nstencil_multi[itype] = n;
  }
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_multi_2d_newton(NeighList *list,
                                            int sx, int sy, int sz)
{
  int i,j,n;
  double rsq,typesq;
  int *s;
  double *distsq;

  int *nstencil_multi = list->nstencil_multi;
  int **stencil_multi = list->stencil_multi;
  double **distsq_multi = list->distsq_multi;

  int ntypes = atom->ntypes;
  for (int itype = 1; itype <= ntypes; itype++) {
    typesq = cuttypesq[itype];
    s = stencil_multi[itype];
    distsq = distsq_multi[itype];
    n = 0;
    for (j = 0; j <= sy; j++)
      for (i = -sx; i <= sx; i++)
        if (j > 0 || (j == 0 && i > 0)) {
          rsq = bin_distance(i,j,0);
          if (rsq < typesq) {
            distsq[n] = rsq;
            s[n++] = j*mbinx + i;
          }
        }
    nstencil_multi[itype] = n;
  }
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_multi_3d_newton(NeighList *list,
                                            int sx, int sy, int sz)
{
  int i,j,k,n;
  double rsq,typesq;
  int *s;
  double *distsq;

  int *nstencil_multi = list->nstencil_multi;
  int **stencil_multi = list->stencil_multi;
  double **distsq_multi = list->distsq_multi;

  int ntypes = atom->ntypes;
  for (int itype = 1; itype <= ntypes; itype++) {
    typesq = cuttypesq[itype];
    s = stencil_multi[itype];
    distsq = distsq_multi[itype];
    n = 0;
    for (k = 0; k <= sz; k++)
      for (j = -sy; j <= sy; j++)
        for (i = -sx; i <= sx; i++)
          if (k > 0 || j > 0 || (j == 0 && i > 0)) {
            rsq = bin_distance(i,j,k);
            if (rsq < typesq) {
              distsq[n] = rsq;
              s[n++] = k*mbiny*mbinx + j*mbinx + i;
            }
          }
    nstencil_multi[itype] = n;
  }
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_multi_2d_newton_tri(NeighList *list,
                                                int sx, int sy, int sz)
{
  int i,j,n;
  double rsq,typesq;
  int *s;
  double *distsq;

  int *nstencil_multi = list->nstencil_multi;
  int **stencil_multi = list->stencil_multi;
  double **distsq_multi = list->distsq_multi;

  int ntypes = atom->ntypes;
  for (int itype = 1; itype <= ntypes; itype++) {
    typesq = cuttypesq[itype];
    s = stencil_multi[itype];
    distsq = distsq_multi[itype];
    n = 0;
    for (j = 0; j <= sy; j++)
      for (i = -sx; i <= sx; i++) {
        rsq = bin_distance(i,j,0);
        if (rsq < typesq) {
          distsq[n] = rsq;
          s[n++] = j*mbinx + i;
        }
      }
    nstencil_multi[itype] = n;
  }
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_multi_3d_newton_tri(NeighList *list,
                                                int sx, int sy, int sz)
{
  int i,j,k,n;
  double rsq,typesq;
  int *s;
  double *distsq;

  int *nstencil_multi = list->nstencil_multi;
  int **stencil_multi = list->stencil_multi;
  double **distsq_multi = list->distsq_multi;

  int ntypes = atom->ntypes;
  for (int itype = 1; itype <= ntypes; itype++) {
    typesq = cuttypesq[itype];
    s = stencil_multi[itype];
    distsq = distsq_multi[itype];
    n = 0;
    for (k = 0; k <= sz; k++)
      for (j = -sy; j <= sy; j++)
        for (i = -sx; i <= sx; i++) {
          rsq = bin_distance(i,j,k);
          if (rsq < typesq) {
            distsq[n] = rsq;
            s[n++] = k*mbiny*mbinx + j*mbinx + i;
          }
        }
    nstencil_multi[itype] = n;
  }
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_full_bin_2d(NeighList *list,
                                   int sx, int sy, int sz)
{
  int i,j;
  int *stencil = list->stencil;
  int nstencil = 0;

  for (j = -sy; j <= sy; j++)
    for (i = -sx; i <= sx; i++)
      if (bin_distance(i,j,0) < cutneighmaxsq)
        stencil[nstencil++] = j*mbinx + i;

  list->nstencil = nstencil;
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_full_ghost_bin_2d(NeighList *list,
                                         int sx, int sy, int sz)
{
  int i,j;
  int *stencil = list->stencil;
  int **stencilxyz = list->stencilxyz;
  int nstencil = 0;

  for (j = -sy; j <= sy; j++)
    for (i = -sx; i <= sx; i++)
      if (bin_distance(i,j,0) < cutneighmaxsq) {
        stencilxyz[nstencil][0] = i;
        stencilxyz[nstencil][1] = j;
        stencilxyz[nstencil][2] = 0;
        stencil[nstencil++] = j*mbinx + i;
      }

  list->nstencil = nstencil;
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_full_bin_3d(NeighList *list,
                                   int sx, int sy, int sz)
{
  int i,j,k;
  int *stencil = list->stencil;
  int nstencil = 0;

  for (k = -sz; k <= sz; k++)
    for (j = -sy; j <= sy; j++)
      for (i = -sx; i <= sx; i++)
        if (bin_distance(i,j,k) < cutneighmaxsq)
          stencil[nstencil++] = k*mbiny*mbinx + j*mbinx + i;

  list->nstencil = nstencil;
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_full_ghost_bin_3d(NeighList *list,
                                         int sx, int sy, int sz)
{
  int i,j,k;
  int *stencil = list->stencil;
  int **stencilxyz = list->stencilxyz;
  int nstencil = 0;

  for (k = -sz; k <= sz; k++)
    for (j = -sy; j <= sy; j++)
      for (i = -sx; i <= sx; i++)
        if (bin_distance(i,j,k) < cutneighmaxsq) {
          stencilxyz[nstencil][0] = i;
          stencilxyz[nstencil][1] = j;
          stencilxyz[nstencil][2] = k;
          stencil[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
        }

  list->nstencil = nstencil;
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_full_multi_2d(NeighList *list,
                                     int sx, int sy, int sz)
{
  int i,j,n;
  double rsq,typesq;
  int *s;
  double *distsq;

  int *nstencil_multi = list->nstencil_multi;
  int **stencil_multi = list->stencil_multi;
  double **distsq_multi = list->distsq_multi;

  int ntypes = atom->ntypes;
  for (int itype = 1; itype <= ntypes; itype++) {
    typesq = cuttypesq[itype];
    s = stencil_multi[itype];
    distsq = distsq_multi[itype];
    n = 0;
    for (j = -sy; j <= sy; j++)
      for (i = -sx; i <= sx; i++) {
        rsq = bin_distance(i,j,0);
        if (rsq < typesq) {
          distsq[n] = rsq;
          s[n++] = j*mbinx + i;
        }
      }
    nstencil_multi[itype] = n;
  }
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_full_multi_3d(NeighList *list,
                                     int sx, int sy, int sz)
{
  int i,j,k,n;
  double rsq,typesq;
  int *s;
  double *distsq;

  int *nstencil_multi = list->nstencil_multi;
  int **stencil_multi = list->stencil_multi;
  double **distsq_multi = list->distsq_multi;

  int ntypes = atom->ntypes;
  for (int itype = 1; itype <= ntypes; itype++) {
    typesq = cuttypesq[itype];
    s = stencil_multi[itype];
    distsq = distsq_multi[itype];
    n = 0;
    for (k = -sz; k <= sz; k++)
      for (j = -sy; j <= sy; j++)
        for (i = -sx; i <= sx; i++) {
          rsq = bin_distance(i,j,k);
          if (rsq < typesq) {
            distsq[n] = rsq;
            s[n++] = k*mbiny*mbinx + j*mbinx + i;
          }
        }
    nstencil_multi[itype] = n;
  }
}
