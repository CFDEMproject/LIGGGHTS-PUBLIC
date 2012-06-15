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

#ifndef LMP_FIX_SPH
#define LMP_FIX_SPH

#include "fix.h"

namespace LAMMPS_NS {

class FixSPH : public Fix {
 public:
  FixSPH(class LAMMPS *, int, char **);
  ~FixSPH();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  virtual void post_integrate() {};
  virtual void post_integrate_respa(int, int);

 protected:
  int iarg;
  int kernel_id;
  double h,hinv; //kernel constant
  double **cutsq;
  class NeighList *list;
  int nlevels_respa;
};

}

#endif
