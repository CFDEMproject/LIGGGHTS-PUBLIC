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

FixStyle(sph/density/summation,FixSPHDensitySum)

#else

#ifndef LMP_FIX_SPH_DENSITY_SUMMATION_H
#define LMP_FIX_SPH_DENSITY_SUMMATION_H

#include "fix_sph.h"

namespace LAMMPS_NS {

class FixSPHDensitySum : public FixSph {
 public:
  FixSPHDensitySum(class LAMMPS *, int, char **);
  ~FixSPHDensitySum();
  virtual int setmask();
  virtual void init();
  virtual void post_integrate();

 private:
  template <int> void post_integrate_eval();

};

}

#endif
#endif
