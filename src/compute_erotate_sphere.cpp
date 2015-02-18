/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   This file was modified with respect to the release in LAMMPS
   Modifications are Copyright 2009-2012 JKU Linz
                     Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
   the producer of the LIGGGHTS® software and the CFDEM®coupling software
   See http://www.cfdem.com/terms-trademark-policy for details.

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "compute_erotate_sphere.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "modify.h" 
#include "fix_multisphere.h" 
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

#define INERTIA 0.4          // moment of inertia prefactor for sphere

/* ---------------------------------------------------------------------- */

ComputeERotateSphere::ComputeERotateSphere(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute erotate/sphere command");

  scalar_flag = 1;
  extscalar = 1;

  // error check

  if (!atom->sphere_flag)
    error->all(FLERR,"Compute erotate/sphere requires atom style sphere");

  fix_ms = 0; 
}

/* ---------------------------------------------------------------------- */

void ComputeERotateSphere::init()
{
  pfactor = 0.5 * force->mvv2e * INERTIA;

  fix_ms =  static_cast<FixMultisphere*>(modify->find_fix_style("multisphere",0)); 
}

/* ---------------------------------------------------------------------- */

double ComputeERotateSphere::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  double **omega = atom->omega;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // sum rotational energy for each particle
  // point particles will not contribute, due to radius = 0.0

  double erotate = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && (!fix_ms || fix_ms->belongs_to(i) < 0)) 
      erotate += (omega[i][0]*omega[i][0] + omega[i][1]*omega[i][1] +
                  omega[i][2]*omega[i][2]) * radius[i]*radius[i]*rmass[i];

  MPI_Allreduce(&erotate,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  scalar *= pfactor;
  return scalar;
}
