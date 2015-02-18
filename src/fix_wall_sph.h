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

/* ----------------------------------------------------------------------
Contributing author for SPH:
Andreas Aigner (CD Lab Particulate Flow Modelling, JKU)
andreas.aigner@jku.at
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(wall/sph,FixWallSph)

#else

#ifndef LMP_FIX_WALL_SPH_H
#define LMP_FIX_WALL_SPH_H

#include "fix_sph.h"

namespace LAMMPS_NS {

class FixWallSph : public FixSph {
 public:
  FixWallSph(class LAMMPS *, int, char **);
  ~FixWallSph();
  int setmask();
  void init();
  void setup(int vflag);
  void post_force(int vflag);
  void post_force_respa(int vflag, int ilevel, int iloop);

 protected:
  int wallstyle;
  double lo,hi,cylradius;
  double r0,D;

};

}

#endif
#endif
