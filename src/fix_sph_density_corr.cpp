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
#include "fix_sph_density_corr.h"
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
#include "fix_property_atom.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSPHDensityCorr::FixSPHDensityCorr(LAMMPS *lmp, int narg, char **arg) :
  FixSPH(lmp, narg, arg)
{
  int iarg = 3;
  if (narg < iarg+1) error->fix_error(FLERR,this,"Not enough arguments");

  if (strcmp(arg[iarg],"shepard") == 0) {
    if (iarg+2 > narg) error->fix_error(FLERR,this,"Not enough arguments");
    if (strcmp(arg[iarg+1],"every") == 0) {
      every = force->inumeric(arg[iarg+2]);
      if (every <= 0) error->fix_error(FLERR,this,"every <= 0 not allowed");
      corrStyle = CORR_SHEPARD;
      iarg += 2;
    } else error->fix_error(FLERR,this,"");
  } else if (strcmp(arg[iarg],"mls") == 0) {
    error->fix_error(FLERR,this,"MLS correction is not implemented until now.");
    corrStyle = CORR_MLS;
  } else error->fix_error(FLERR,this,"Unknown style for fix sph/density/corr. Valid styles are 'shepard' or 'mls'");

  quantity_name = new char[strlen("corrKernel")+1];
  strcpy(quantity_name,"corrKernel");

  fix_quantity = NULL;

  peratom_flag = 1;
  size_peratom_cols = 0;
  peratom_freq = 1;
  time_depend = 1;

  scalar_flag = 1;
  global_freq = 1;

  ago = 0;

}

/* ---------------------------------------------------------------------- */

FixSPHDensityCorr::~FixSPHDensityCorr()
{
    delete []quantity_name;
}

/* ---------------------------------------------------------------------- */

void FixSPHDensityCorr::pre_delete(bool unfixflag)
{
    //unregister property/atom fixes
    if (fix_quantity) modify->delete_fix(quantity_name);
}

/* ---------------------------------------------------------------------- */

int FixSPHDensityCorr::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSPHDensityCorr::updatePtrs()
{
  quantity = fix_quantity->vector_atom;

  vector_atom = quantity;
}

/* ---------------------------------------------------------------------- */

void FixSPHDensityCorr::post_create()
{
  char **fixarg;
  fixarg=new char*[9];
  for (int kk=0;kk<9;kk++) fixarg[kk]=new char[30];

  if (fix_quantity==NULL) {
    strcpy(fixarg[0],quantity_name);
    fixarg[1]="all";
    fixarg[2]="property/atom";
    strcpy(fixarg[3],quantity_name);
    fixarg[4]="scalar";
    fixarg[5]="yes";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.";
    modify->add_fix(9,fixarg);
    fix_quantity=static_cast<FixPropertyAtom*>(modify->find_fix_property(quantity_name,"property/atom","scalar",0,0,style));
  }

  delete []fixarg;

  updatePtrs();
}

/* ---------------------------------------------------------------------- */

void FixSPHDensityCorr::init()
{
  FixSPH::init();
}

/* ---------------------------------------------------------------------- */

void FixSPHDensityCorr::post_integrate()
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

  int nlocal = atom->nlocal;

  updatePtrs();

  ago++;
  if (ago % every == 0) {

  ago = 0;

  // kernel normalization

  for (i = 0; i < nlocal; i++)
  {
    if (mask[i] & groupbit) {

      // this gets a value for W at self, perform error check

      W = SPH_KERNEL_NS::sph_kernel(kernel_id,0.,h,hinv);
      if (W < 0.)
      {
        fprintf(screen,"s = 0, W = %f\n",W);
        error->one(FLERR,"Illegal kernel used, W < 0");
      }

      quantity[i] = W * mass[type[i]] / density[i];
    }
  }

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

      quantity[i] += W * mass[jtype] / density[j];

      if (newton_pair || j < nlocal)
        quantity[j] += W * mass[itype] / density[i];
    }
  }

  // reset and add density contribution of self

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

  // normalize density
  for (i = 0; i < nlocal; i++) {
      density[i] = density[i]/quantity[i];
  }

  // density is now correct, send to ghosts

  comm->forward_comm();

//  fprintf(screen,"ts %d, particles are reinitaialized. \n",update->ntimestep);

  }

}
