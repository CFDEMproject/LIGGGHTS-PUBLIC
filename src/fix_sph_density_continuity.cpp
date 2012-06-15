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

/* ----------------------------------------------------------------------
Contributing author for SPH:
Andreas Aigner (CD Lab Particulate Flow Modelling, JKU)
andreas.aigner@jku.at
------------------------------------------------------------------------- */

#include "math.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "fix_sph_density_continuity.h"
#include "update.h"
#include "respa.h"
#include "atom.h"
#include "force.h"
#include "modify.h"
#include "pair.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "sph_kernels.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSPHDensityContinuity::FixSPHDensityContinuity(LAMMPS *lmp, int narg, char **arg) :
  FixSPH(lmp, narg, arg)
{

}

/* ---------------------------------------------------------------------- */

FixSPHDensityContinuity::~FixSPHDensityContinuity()
{

}

/* ---------------------------------------------------------------------- */

int FixSPHDensityContinuity::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSPHDensityContinuity::init()
{
  int i;

  FixSPH::init();

/*  // check if I am first SPH fix
  // if I am not first, something must be mis-ordered

  int first = 10000, mine = -1;
  for(int i = 0; i < modify->nfix; i++)
  {
      if(modify->fix[i] == this) mine = i;
      if(strncmp("sph",modify->fix[i]->style,3) == 0 && i < first) first = i;
  }

  if (first < mine) error->all(FLERR,"There must be only one fix sph/density, and it must come before other sph fixes"); */

  dtdensity = update->dt;
}

/* ---------------------------------------------------------------------- */

void FixSPHDensityContinuity::post_integrate()
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,r,rinv,s,gradWmag;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  double **v = atom->v;
  double *density = atom->density;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int newton_pair = force->newton_pair;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // need updated ghost positions and self contributions

  comm->forward_comm();

  // loop over neighbors of my atoms

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      if (!(mask[j] & groupbit)) continue;
      jtype = type[j];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq >= cutsq[itype][jtype]) continue;

      // calculate distance and normalized distance

      r = sqrt(rsq);
      rinv = 1./r;
      s = r * hinv;

      // calculate value for magnitude of grad W

      gradWmag = SPH_KERNEL_NS::sph_kernel_der(kernel_id,s,h,hinv);

      // add contribution of neighbor
      // have a half neigh list, so do it for both if necessary

      density[i] += delx * rinv * dtdensity * mass[jtype] * (v[i][0]-v[j][0]) * gradWmag;
      density[i] += dely * rinv * dtdensity * mass[jtype] * (v[i][1]-v[j][1]) * gradWmag;
      density[i] += delz * rinv * dtdensity * mass[jtype] * (v[i][2]-v[j][2]) * gradWmag;

      if (newton_pair || j < nlocal) {
        density[j] += -delx * rinv * dtdensity * mass[itype] * (v[j][0]-v[i][0]) * gradWmag;
        density[j] += -dely * rinv * dtdensity * mass[itype] * (v[j][1]-v[i][1]) * gradWmag;
        density[j] += -delz * rinv * dtdensity * mass[itype] * (v[j][2]-v[i][2]) * gradWmag;
      }
    }
  }

  // density is now correct, send to ghosts

  comm->forward_comm();

}

/* ---------------------------------------------------------------------- */

void FixSPHDensityContinuity::post_integrate_respa(int ilevel, int flag)
{
  if (flag) return;             // only used by NPT,NPH

  dtdensity = step_respa[ilevel];

  if (ilevel == 0) post_integrate();

}

/* ---------------------------------------------------------------------- */

void FixSPHDensityContinuity::reset_dt()
{
  dtdensity = update->dt;
}
