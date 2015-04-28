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

#include "string.h"
#include "compute_ke_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "fix_multisphere.h" 
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeKEAtom::ComputeKEAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute ke/atom command");

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  ke = NULL;
  fix_ms = NULL; 
}

/* ---------------------------------------------------------------------- */

ComputeKEAtom::~ComputeKEAtom()
{
  memory->destroy(ke);
}

/* ---------------------------------------------------------------------- */

void ComputeKEAtom::init()
{
  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"ke/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute ke/atom");
  fix_ms =  static_cast<FixMultisphere*>(modify->find_fix_style("multisphere",0)); 
}

/* ---------------------------------------------------------------------- */

void ComputeKEAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow ke array if necessary

  if (atom->nlocal > nmax) {
    memory->destroy(ke);
    nmax = atom->nmax;
    memory->create(ke,nmax,"ke/atom:ke");
    vector_atom = ke;
  }

  // compute kinetic energy for each atom in group

  double mvv2e = force->mvv2e;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  if (rmass)
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit && (!fix_ms || fix_ms->belongs_to(i) < 0)) { 
        ke[i] = 0.5 * mvv2e * rmass[i] *
          (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
      } else ke[i] = 0.0;
    }

  else
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        ke[i] = 0.5 * mvv2e * mass[type[i]] *
          (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
      } else ke[i] = 0.0;
    }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeKEAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
