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

#include "string.h"
#include "compute_erotate_sphere_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "fix_multisphere.h" 

using namespace LAMMPS_NS;

#define INERTIA 0.4          // moment of inertia prefactor for sphere

/* ---------------------------------------------------------------------- */

ComputeErotateSphereAtom::
ComputeErotateSphereAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3)
    error->all(FLERR,"Illegal compute erotate/sphere//atom command");

  peratom_flag = 1;
  size_peratom_cols = 0;

  // error check

  if (!atom->sphere_flag)
    error->all(FLERR,"Compute erotate/sphere/atom requires atom style sphere");

  nmax = 0;
  erot = NULL;

  fix_ms = 0; 
}

/* ---------------------------------------------------------------------- */

ComputeErotateSphereAtom::~ComputeErotateSphereAtom()
{
  memory->destroy(erot);
}

/* ---------------------------------------------------------------------- */

void ComputeErotateSphereAtom::init()
{
  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"erotate/sphere/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute erotate/sphere/atom");

  pfactor = 0.5 * force->mvv2e * INERTIA;

  fix_ms =  static_cast<FixMultisphere*>(modify->find_fix_style("multisphere",0)); 
}

/* ---------------------------------------------------------------------- */

void ComputeErotateSphereAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow erot array if necessary

  if (atom->nlocal > nmax) {
    memory->destroy(erot);
    nmax = atom->nmax;
    memory->create(erot,nmax,"erotate/sphere/atom:erot");
    vector_atom = erot;
  }

  // compute rotational kinetic energy for each atom in group
  // point particles will have erot = 0.0, due to radius = 0.0

  double **omega = atom->omega;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit && (!fix_ms || fix_ms->belongs_to(i) < 0)) 
    {
      erot[i] = (omega[i][0]*omega[i][0] + omega[i][1]*omega[i][1] +
                 omega[i][2]*omega[i][2]) * radius[i]*radius[i]*rmass[i];
      erot[i] *= pfactor;
    } else erot[i] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeErotateSphereAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
