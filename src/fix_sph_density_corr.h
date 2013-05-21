/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
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

FixStyle(sph/density/corr,FixSphDensityCorr)

#else

#ifndef LMP_FIX_SPH_DENSITY_CORR_H
#define LMP_FIX_SPH_DENSITY_CORR_H

#include "fix_sph.h"

namespace LAMMPS_NS {

  enum {CORR_SHEPARD,CORR_MLS};

class FixSphDensityCorr : public FixSph {
 public:
  FixSphDensityCorr(class LAMMPS *, int, char **);
  ~FixSphDensityCorr();
  void pre_delete(bool unfixflag);
  virtual int setmask();
  void updatePtrs();
  void post_create();
  virtual void init();
  virtual void pre_force(int);

 private:
  template <int> void pre_force_eval();

  class FixPropertyAtom* fix_quantity;
  char *quantity_name;

  double quantity_0;
  double *quantity;

  int corrStyle;
  int every;
  int ago;

};

}

#endif
#endif
