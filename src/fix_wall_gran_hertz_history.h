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
   Contributing authors for original version: Leo Silbert (SNL), Gary Grest (SNL)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(wall/gran/hertz/history,FixWallGranHertzHistory)

#else

#ifndef LMP_FIX_WALL_GRAN_HERTZ_HISTORY_H
#define LMP_FIX_WALL_GRAN_HERTZ_HISTORY_H

#include "fix.h"
#include "fix_wall_gran_hooke_history.h"

namespace LAMMPS_NS {

class FixWallGranHertzHistory : public FixWallGranHookeHistory {
 public:
  FixWallGranHertzHistory(class LAMMPS *, int, char **);
  void init_granular();

 protected:

  virtual void deriveContactModelParams(int ip, double deltan,double meff_wall, double &kn, double &kt, double &gamman, double &gammat, double &xmu,double &rmu,double &vnnr);
};

}

#endif
#endif
