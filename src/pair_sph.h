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

#ifdef PAIR_CLASS

PairStyle(sph,PairSPH)

#else

#ifndef LMP_PAIR_SPH
#define LMP_PAIR_SPH

#include "pair.h"

namespace LAMMPS_NS {

class PairSPH : public Pair {
 friend class FixSPH;
 public:
  PairSPH(class LAMMPS *);
  ~PairSPH();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  virtual void init_style();
  void init_list(int, class NeighList *);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  int sph_kernel_id(){return kernel_id;}
  double sph_kernel_h(){return h;}

 protected:
  double **cut_sq(){return cutsq;}

  double h,hinv,cut_global;
  double **cut;
  void allocate();
  int kernel_id;
  char *kernel_style;

  int artVisc_flag, tensCorr_flag;
  double alpha,beta,cAB,eta; //artifical viscosity
  double epsilon,WdeltaPinv; //tensile instability
};

}

#endif
#endif
