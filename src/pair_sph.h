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

#ifdef PAIR_CLASS

#else

#ifndef LMP_PAIR_SPH_H
#define LMP_PAIR_SPH_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSph : public Pair {
 public:

  friend class FixSPH;

  PairSph(class LAMMPS *);
  ~PairSph();

  /* INHERITED FROM Pair */

  virtual void compute(int, int) = 0;
  virtual void settings(int, char **) = 0;
  virtual void setKernelAndLength(int, char **);
  virtual void coeff(int, char **) = 0;
  virtual void init_style();
  virtual void init_substyle() = 0;
  virtual void init_list(int, class NeighList *);
  virtual double init_one(int, int);
  virtual void write_restart(FILE *){}
  virtual void read_restart(FILE *){}
  virtual void write_restart_settings(FILE *){}
  virtual void read_restart_settings(FILE *){}
  //virtual void reset_dt();

  /* PUBLIC ACCESS FUNCTIONS */

  int sph_kernel_id(){return kernel_id;}
  int returnPairStyle(){return pairStyle_; };
  double returnViscosity() {return viscosity_; };

 protected:

  void allocate();
  virtual void updatePtrs();
  virtual void updateRadius();
  //virtual double interpDist(double, double);
  inline double interpDist(double disti, double distj) {return 0.5*(disti+distj);}

  class FixPropertyAtom* fppaSl; //fix for smoothing length
  class FixPropertyGlobal* fppaSlType; //fix for per type smoothing length
  double *sl;         // per atom smoothing length
  double **slComType; // common smoothing length in case of mass_type=1
  double sl_0;

  int kernel_id;
  char *kernel_style;

  double *onerad;
  double *maxrad;

  int mass_type; // flag defined in atom_vec*

  int pairStyle_;
  double viscosity_;

  // storage for force part caused by pressure gradient (grad P / rho):
  class FixPropertyAtom* fix_fgradP_;
  double **fgradP_;
};

}

#endif
#endif
