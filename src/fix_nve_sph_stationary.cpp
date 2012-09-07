/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "stdio.h"
#include "string.h"
#include "fix_nve_sph_stationary.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVESphStationary::FixNVESphStationary(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if ((atom->e_flag != 1) || (atom->rho_flag != 1))
    error->all(FLERR,
        "fix sph command requires atom_style with both energy and density");

  if (narg < 3)
    error->all(FLERR,"Illegal fix nve command");

  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixNVESphStationary::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVESphStationary::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixNVESphStationary::initial_integrate(int vflag)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *rho = atom->rho;
  double *drho = atom->drho;
  double *e = atom->e;
  double *de = atom->de;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      // same like src/USER-SPH/fix_meso !!
      e[i] += dtf * de[i]; // half-step update of particle internal energy
      //rho[i] += dtf * drho[i]; // ... and density
      rho[i] += 2.0 * dtf * drho[i]; // ... and density
    }
}

/* ---------------------------------------------------------------------- */

void FixNVESphStationary::final_integrate()
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *e = atom->e;
  double *de = atom->de;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      e[i] += dtf * de[i]; // half-step update of particle internal energy
    }
}

/* ---------------------------------------------------------------------- */

void FixNVESphStationary::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
