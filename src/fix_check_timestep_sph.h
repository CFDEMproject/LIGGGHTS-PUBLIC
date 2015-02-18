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

FixStyle(check/timestep/sph,FixCheckTimestepSph)

#else

#ifndef LMP_FIX_CHECK_TIMESTEP_SPH_H
#define LMP_FIX_CHECK_TIMESTEP_SPH_H

#include "fix_sph.h"

namespace LAMMPS_NS {

class FixCheckTimestepSph : public FixSph {
 public:
  FixCheckTimestepSph(class LAMMPS *, int, char **);
  int setmask();
  void init();
  void end_of_step();
  double compute_vector(int);

 private:
  class FixPropertyGlobal* cs;
  void calc_courant_estims();
  template <int> void calc_courant_estims_eval();

//  double r_min;
  double vmax; //max relatice velocity
  double mumax; //max viscous term
  double courant_time;
  double fraction_courant,fraction_skin;
  double fraction_courant_lim;
  bool warnflag;
};

}

#endif
#endif
