/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   This file was modified with respect to the release in LAMMPS
   Modifications are Copyright 2009-2012 JKU Linz
                     Copyright 2012-     DCS Computing GmbH, Linz

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#ifdef INTEGRATE_CLASS

IntegrateStyle(verlet/implicit,VerletImplicit)

#else

#ifndef LMP_VERLET_IMPLICIT_H
#define LMP_VERLET_IMPLICIT_H

#include "verlet.h"

namespace LAMMPS_NS {

class VerletImplicit : public Verlet  {
 public:
  VerletImplicit(class LAMMPS *, int, char **);
  virtual ~VerletImplicit() {}

  virtual void run(int);

};

}

#endif
#endif

/* ERROR/WARNING messages:

W: No fixes defined, atoms won't move

If you are not using a fix like nve, nvt, npt then atom velocities and
coordinates will not be updated during timestepping.

*/
