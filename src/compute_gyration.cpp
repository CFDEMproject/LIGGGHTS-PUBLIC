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

#include <cmath>
#include "compute_gyration.h"
#include "update.h"
#include "atom.h"
#include "group.h"
#include "domain.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeGyration::ComputeGyration(LAMMPS *lmp, int &iarg, int narg, char **arg) :
  Compute(lmp, iarg, narg, arg)
{
  if (narg != iarg) error->all(FLERR,"Illegal compute gyration command");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 0;

  vector = new double[6];
}

/* ---------------------------------------------------------------------- */

ComputeGyration::~ComputeGyration()
{
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeGyration::init()
{
  masstotal = group->mass(igroup);
}

/* ---------------------------------------------------------------------- */

double ComputeGyration::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  double xcm[3];
  group->xcm(igroup,masstotal,xcm);
  scalar = group->gyration(igroup,masstotal,xcm);
  return scalar;
}

/* ----------------------------------------------------------------------
   compute the radius-of-gyration tensor of group of atoms
   around center-of-mass cm
   must unwrap atoms to compute Rg tensor correctly
------------------------------------------------------------------------- */

void ComputeGyration::compute_vector()
{
  invoked_vector = update->ntimestep;

  double xcm[3];
  group->xcm(igroup,masstotal,xcm);

  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  tagint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  double dx,dy,dz,massone;
  double unwrap[3];

  double rg[6];
  rg[0] = rg[1] = rg[2] = rg[3] = rg[4] = rg[5] = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];

      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - xcm[0];
      dy = unwrap[1] - xcm[1];
      dz = unwrap[2] - xcm[2];

      rg[0] += dx*dx * massone;
      rg[1] += dy*dy * massone;
      rg[2] += dz*dz * massone;
      rg[3] += dx*dy * massone;
      rg[4] += dx*dz * massone;
      rg[5] += dy*dz * massone;
    }
  MPI_Allreduce(rg,vector,6,MPI_DOUBLE,MPI_SUM,world);

  if (masstotal == 0.0) return;
  for (int i = 0; i < 6; i++) vector[i] /= masstotal;
}
