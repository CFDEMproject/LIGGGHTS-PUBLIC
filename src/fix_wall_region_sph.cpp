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
#include "stdlib.h"
#include "string.h"
#include "fix_wall_region.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "region.h"
#include "lattice.h"
#include "update.h"
#include "output.h"
#include "respa.h"
#include "error.h"

#include "fix_wall_region_sph.h"
#include "sph_kernels.h"
#include "force.h"
#include "pair.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixWallRegionSph::FixWallRegionSph(LAMMPS *lmp, int narg, char **arg) :
  FixSPH(lmp, narg, arg)
{
  if (narg != 6) error->all(FLERR,"Illegal fix wall/region/sph command");

  int iarg = 3;

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;

  // parse args

  iregion = domain->find_region(arg[iarg]);
  if (iregion == -1) error->all(FLERR,"Fix wall/region/sph region ID does not exist");

  r0 = force->numeric(arg[iarg+1]);
  D  = force->numeric(arg[iarg+2]);

  iarg += 3;

  //only for force!
  cutoff = r0;
  if (cutoff <= 0.0) error->all(FLERR,"Fix wall/region/sph cutoff <= 0.0");

  eflag = 0;
  ewall[0] = ewall[1] = ewall[2] = ewall[3] = 0.0;
}

/* ---------------------------------------------------------------------- */

int FixWallRegionSph::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallRegionSph::init()
{

  FixSPH::init();

  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixWallRegionSph::setup(int vflag)
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixWallRegionSph::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixWallRegionSph::post_force(int vflag)
{
  int i,m,n;
  double fx,fy,fz;

  eflag = 0;
  ewall[0] = ewall[1] = ewall[2] = ewall[3] = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  Region *region = domain->regions[iregion];
  int onflag = 0;

  double s,rinv,gradWmag;
  int itype;
  double *density = atom->density;
  double *mass = atom->mass;
  double *q = atom->q;

  // region->match() insures particle is in region or on surface, else error
  // if returned contact dist r = 0, is on surface, also an error

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      itype = type[i];
      cutoff = sqrt(cutsq[itype][itype]);

      if (!region->match(x[i][0],x[i][1],x[i][2])) {
        onflag = 1;
        //fprintf(screen,"Particle %d with the Coordinates x= %f, y= %f, z= %f. \n",i,x[i][0],x[i][1],x[i][2]); //TEST OUTPUT
        continue;
      }

      // number of wall contacts
      n = region->surface(x[i][0],x[i][1],x[i][2],cutoff);

      for (m = 0; m < n; m++) {

        // check wall distance
        if (region->contact[m].r <= 0.0) {
          onflag = 1;
          //fprintf(screen,"Particle %d with the Coordinates x= %f, y= %f, z= %f has zero distance. \n",i,x[i][0],x[i][1],x[i][2]); // TEST OUTPUT
          continue;
        }

        // fwall, due to self influence
        // a small perturbation to avoid stacking particles
        s = region->contact[m].r * hinv;
        rinv = 1./region->contact[m].r;

        // calculate value for magnitude of grad W
        gradWmag = SPH_KERNEL_NS::sph_kernel_der(kernel_id,s,h,hinv);

        // zero order approximation - properties at wall assumed to be equal to properties at particle (thuis factor 2)
        fwall = - rinv * mass[itype] * mass[itype] * (2.*q[i]/(density[i]*density[i])) * gradWmag;

        fx = fwall * region->contact[m].delx;
        fy = fwall * region->contact[m].dely;
        fz = fwall * region->contact[m].delz;
        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
        ewall[1] -= fx;
        ewall[2] -= fy;
        ewall[3] -= fz;

        // fwall, due to repulsive force
        repulsivsph(region->contact[m].r);

        ewall[0] += eng;
        fx = fwall * region->contact[m].delx;
        fy = fwall * region->contact[m].dely;
        fz = fwall * region->contact[m].delz;
        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
        ewall[1] -= fx;
        ewall[2] -= fy;
        ewall[3] -= fz;
      }
    }

    if (onflag) {
      error->one(FLERR,"Particle on or inside fix wall/region/sph surface \n");
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixWallRegionSph::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixWallRegionSph::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   energy of wall interaction
------------------------------------------------------------------------- */

double FixWallRegionSph::compute_scalar()
{
  // only sum across procs one time

  if (eflag == 0) {
    MPI_Allreduce(ewall,ewall_all,4,MPI_DOUBLE,MPI_SUM,world);
    eflag = 1;
  }
  return ewall_all[0];
}

/* ----------------------------------------------------------------------
   components of force on wall
------------------------------------------------------------------------- */

double FixWallRegionSph::compute_vector(int n)
{
  // only sum across procs one time

  if (eflag == 0) {
    MPI_Allreduce(ewall,ewall_all,4,MPI_DOUBLE,MPI_SUM,world);
    eflag = 1;
  }
  return ewall_all[n+1];
}

/* ----------------------------------------------------------------------
   repulsiv force for sph simulation
   compute eng and fwall = magnitude of wall force
------------------------------------------------------------------------- */

void FixWallRegionSph::repulsivsph(double r)
{
  double rinv,frac,frac2,frac4;

  if (r <= r0) {
    rinv = 1.0/r;
    frac = r0*rinv;
    frac2 = frac*frac;
    frac4 = frac2*frac2;
    fwall = D * (frac4 - frac2) * rinv;
    eng = 0;
  } else {
    fwall = 0;
    eng = 0;
  }
}
