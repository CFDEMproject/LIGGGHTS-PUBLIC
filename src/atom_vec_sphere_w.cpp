/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
   the producer of the LIGGGHTS® software and the CFDEM®coupling software
   See http://www.cfdem.com/terms-trademark-policy for details.

   LIGGGHTS® is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include "atom_vec_sphere.h"
#include "domain_wedge.h"

#ifndef DOMAIN_WEDGE_REAL_H
#define DOMAIN_WEDGE_REAL_H

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

int AtomVecSphere::pack_border_vel_wedge(int n, int *list, double *buf,
                                     int pbc_flag, int *pbc)
{
    return 0;
}

/* ---------------------------------------------------------------------- */

int AtomVecSphere::pack_comm_vel_wedge(int n, int *list, double *buf,
                                   int pbc_flag, int *pbc)
{
    return 0;
}

#endif
