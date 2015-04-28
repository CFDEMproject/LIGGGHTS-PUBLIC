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
#include "fix_enforce2d.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "respa.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixEnforce2D::FixEnforce2D(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal fix enforce2d command");
}

/* ---------------------------------------------------------------------- */

int FixEnforce2D::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixEnforce2D::init()
{
  if (domain->dimension == 3)
    error->all(FLERR,"Cannot use fix enforce2d with 3d simulation");
}

/* ---------------------------------------------------------------------- */

void FixEnforce2D::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    int nlevels_respa = ((Respa *) update->integrate)->nlevels;
    for (int ilevel = 0; ilevel < nlevels_respa; ilevel++) {
      ((Respa *) update->integrate)->copy_flevel_f(ilevel);
      post_force_respa(vflag,ilevel,0);
      ((Respa *) update->integrate)->copy_f_flevel(ilevel);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixEnforce2D::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixEnforce2D::post_force(int vflag)
{
  double **v = atom->v;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      v[i][2] = 0.0;
      f[i][2] = 0.0;
    }

  // for systems with omega/angmom/torque, zero x and y components

  if (atom->omega_flag) {
    double **omega = atom->omega;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        omega[i][0] = 0.0;
        omega[i][1] = 0.0;
      }
  }

  if (atom->angmom_flag) {
    double **angmom = atom->angmom;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        angmom[i][0] = 0.0;
        angmom[i][1] = 0.0;
      }
  }

  if (atom->torque_flag) {
    double **torque = atom->torque;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        torque[i][0] = 0.0;
        torque[i][1] = 0.0;
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixEnforce2D::post_force_respa(int vflag, int ilevel, int iloop)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixEnforce2D::min_post_force(int vflag)
{
  post_force(vflag);
}
