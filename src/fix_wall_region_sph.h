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

FixStyle(wall/region/sph,FixWallRegionSph)

#else

#ifndef LMP_FIX_WALL_REGION_SPH_H
#define LMP_FIX_WALL_REGION_SPH_H

#include "fix_sph.h"

namespace LAMMPS_NS {

class FixWallRegionSph : public FixSPH {
 public:
  FixWallRegionSph(class LAMMPS *, int, char **);
  ~FixWallRegionSph() {}
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();
  double compute_vector(int);

 private:
  int iregion;
  double cutoff;
  int eflag;
  double ewall[4],ewall_all[4];
  int nlevels_respa;
  double dt;

  double eng,fwall;

  double r0,D; // coefficient for repulsivsph
  void repulsivsph(double);
};

}

#endif
#endif
