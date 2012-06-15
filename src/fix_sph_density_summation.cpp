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
#include "fix_sph_density_summation.h"
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

FixSPHDensitySum::FixSPHDensitySum(LAMMPS *lmp, int narg, char **arg) :
  FixSPH(lmp, narg, arg)
{

}

/* ---------------------------------------------------------------------- */

FixSPHDensitySum::~FixSPHDensitySum()
{

}

/* ---------------------------------------------------------------------- */

int FixSPHDensitySum::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSPHDensitySum::init()
{
  FixSPH::init();

}

/* ---------------------------------------------------------------------- */

void FixSPHDensitySum::post_integrate()
{
  int i,j,m,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,r,s,W;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  double *density = atom->density;
  double *mass = atom->mass;
  int newton_pair = force->newton_pair;

  // reset and add density contribution of self

  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++) {

    // this gets a value for W at self, perform error check

    W = SPH_KERNEL_NS::sph_kernel(kernel_id,0.,h,hinv);
    if (W < 0.)
    {
        fprintf(screen,"s = 0, W = %f\n",W);
        error->one(FLERR,"Illegal kernel used, W < 0");
    }

    // add contribution of self

    density[i] = mass[type[i]] * W;
  }

/*
  for (i = 0; i < nlocal; i++)
    fprintf(screen,"ts %d, particle %d after self: density %f\n",update->ntimestep,i,density[i]);
*/

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
      s = r * hinv;

      // this sets a value for W
      // error on kernel existence already performed in FixSPH::FixSPH()

      W = SPH_KERNEL_NS::sph_kernel(kernel_id,s,h,hinv);
      if (W < 0.)
      {
          fprintf(screen,"s = %f, W = %f\n",s,W);
          error->one(FLERR,"Illegal kernel used, W < 0");
      }

      // add contribution of neighbor
      // have a half neigh list, so do it for both if necessary

      density[i] += mass[jtype] * W;

      if (newton_pair || j < nlocal)
        density[j] += mass[itype] * W;
    }
  }

  // density is now correct, send to ghosts

  comm->forward_comm();

/*
  for (i = 0; i < nlocal; i++)
    fprintf(screen,"ts %d, particle %d after neigh: density %f\n",update->ntimestep,i,density[i]);

  error->all(FLERR,"end");
*/
}
