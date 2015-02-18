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
-------------------------------------------------------------------------
Contributing author for SPH:
Andreas Eitzlmayr (Institute for Process and Particle Engineering, TU Graz)
andreas.eitzlmayr@tugraz.at
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(sph/pressure,FixSPHPressure)

#else

#ifndef LMP_FIX_SPH_PRESSURE_H
#define LMP_FIX_SPH_PRESSURE_H

#include "fix_sph.h"

namespace LAMMPS_NS {

  enum {PRESSURESTYLE_ABSOLUT,PRESSURESTYLE_TAIT,PRESSURESTYLE_RELATIV};

class FixSPHPressure : public FixSph {
 public:
  FixSPHPressure(class LAMMPS *, int, char **);
  ~FixSPHPressure();
  int setmask();
  void init();
  void pre_force(int);

  double return_rho0() {
    if (pressureStyle == PRESSURESTYLE_ABSOLUT) return 0;
    else return rho0;
  };

 private:
  int pressureStyle;
  double B,rho0,rho0inv,gamma,P0;
};

}

#endif
#endif
