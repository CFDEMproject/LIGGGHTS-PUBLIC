/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com


   Copyright (C): 2014 DCS Computing GmbH (www.dcs-computing.com), Linz, Austria
                  2014 Graz University of Technology (ippt.tugraz.at), Graz, Austria

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   Parts of the code were developped in the frame of the NanoSim project funded
   by the European Commission through FP7 Grant agreement no. 604656.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(couple/pascal,FixPaScalCouple)

#else

#ifndef LMP_FIX_PASCAL_COUPLE_H
#define LMP_FIX_PASCAL_COUPLE_H

#include "fix.h"
namespace PASCAL_NS { class PaScal; }

namespace LAMMPS_NS {

class FixPaScalCouple : public Fix  {
 public:
  FixPaScalCouple(class LAMMPS *, int narg, char **arg);
  ~FixPaScalCouple();
  virtual void post_create();
  void      updatePtrs();

  int       setmask();
  void      init();
  void      setup(int);

  void      pre_exchange();
  void      end_of_step();

  int*      get_liggghts_map(int &length);

  void*     find_pull_property(const char *name, const char *type, int &len1, int &len2);

  void*     find_push_property(const char *name, const char *type, int &len1, int &len2);

 protected:

  class FixPropertyAtom* fix_shellTemperature_;
  class FixPropertyAtom* fix_shellHeatFlux_;

 private:

  // data transfer is handled by this class
  class CfdDatacouplingSimple *dc_;

  bool      verbose_;
  int       couple_at_least_every_;
  bool      couple_this_step_;
  bool      pascal_setup_;
  double    time_;
  int       iarg_;    

  // PaScal Object
  PASCAL_NS::PaScal *pasc_;
};

} //end namespace
#endif
#endif

