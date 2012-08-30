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

#ifdef PAIR_CLASS

PairStyle(gran/hooke/history/stiffness,PairGranHookeHistorySimple)

#else

#ifndef LMP_PAIR_GRAN_HOOKE_HISTORY_SIMPLE_H
#define LMP_PAIR_GRAN_HOOKE_HISTORY_SIMPLE_H

#include "pair_gran_hooke_history.h"

namespace LAMMPS_NS {

class PairGranHookeHistorySimple : public PairGranHookeHistory {

 friend class FixWallGranHookeHistorySimple;

 public:

  PairGranHookeHistorySimple(class LAMMPS *);
  ~PairGranHookeHistorySimple();

  virtual void settings(int, char **);
  virtual void init_granular();

 protected:
  virtual void allocate_properties(int);
  virtual void deriveContactModelParams(int &, int &,double &, double &, double &,double &, double &, double &, double &,double &);

  //stiffness and damp parameters
  class FixPropertyGlobal *k_n1,*k_t1,*gamma_n1,*gamma_t1;
  double **k_n,**k_t,**gamma_n,**gamma_t;
  int damp_massflag;
};

}

#endif
#endif
