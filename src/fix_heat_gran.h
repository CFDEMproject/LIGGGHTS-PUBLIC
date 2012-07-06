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
#ifdef FIX_CLASS

FixStyle(heat/gran,FixHeatGran)

#else

#ifndef LMP_FIX_HEATGRAN_H
#define LMP_FIX_HEATGRAN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixHeatGran : public Fix {
 public:

  FixHeatGran(class LAMMPS *, int, char **);
  ~FixHeatGran();
  int setmask();
  void post_create();
  void pre_delete(bool unfixflag);
  void init();
  void updatePtrs();
  void post_force(int);
  double compute_scalar();
  void cpl_evaluate(class ComputePairGranLocal *);
  void register_compute_pair_local(class ComputePairGranLocal *ptr);
  void unregister_compute_pair_local(class ComputePairGranLocal *ptr);

 private:

  template <int> void post_force_eval(int,int);

  class FixPropertyAtom* fix_temp;
  class FixPropertyAtom* fix_heatFlux;
  class FixPropertyAtom* fix_heatSource;
  class FixPropertyGlobal* fix_conductivity;
  class FixScalarTransportEquation *fix_ste;

  class ComputePairGranLocal *cpl;

  double T0;              
  double *Temp;           
  double *heatFlux;       
  double *heatSource;     
  double *conductivity;

  // for heat transfer area correction
  int area_correction_flag;
  double const* const* deltan_ratio;

  class PairGran *pair_gran;
  int history_flag;
};

}

#endif
#endif
