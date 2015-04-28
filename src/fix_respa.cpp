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

#include "stdlib.h"
#include "fix_respa.h"
#include "atom.h"
#include "memory.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixRespa::FixRespa(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  // nlevels = # of rRESPA levels

  nlevels = force->inumeric(FLERR,arg[3]);

  // perform initial allocation of atom-based arrays
  // register with Atom class

  f_level = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
}

/* ---------------------------------------------------------------------- */

FixRespa::~FixRespa()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);

  // delete locally stored arrays

  memory->destroy(f_level);
}

/* ---------------------------------------------------------------------- */

int FixRespa::setmask()
{
  return 0;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixRespa::memory_usage()
{
  double bytes = atom->nmax*nlevels*3 * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixRespa::grow_arrays(int nmax)
{
  memory->grow(f_level,nmax,nlevels,3,"fix_respa:f_level");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixRespa::copy_arrays(int i, int j, int delflag)
{
  for (int k = 0; k < nlevels; k++) {
    f_level[j][k][0] = f_level[i][k][0];
    f_level[j][k][1] = f_level[i][k][1];
    f_level[j][k][2] = f_level[i][k][2];
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixRespa::pack_exchange(int i, double *buf)
{
  int m = 0;
  for (int k = 0; k < nlevels; k++) {
    buf[m++] = f_level[i][k][0];
    buf[m++] = f_level[i][k][1];
    buf[m++] = f_level[i][k][2];
  }
  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixRespa::unpack_exchange(int nlocal, double *buf)
{
  int m = 0;
  for (int k = 0; k < nlevels; k++) {
    f_level[nlocal][k][0] = buf[m++];
    f_level[nlocal][k][1] = buf[m++];
    f_level[nlocal][k][2] = buf[m++];
  }
  return m;
}
