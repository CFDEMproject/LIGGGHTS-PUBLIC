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

#ifdef FIX_CLASS

FixStyle(region/variable,FixRegionVariable)

#else

#ifndef LMP_FIX_REGION_VARIABLE_H
#define LMP_FIX_REGION_VARIABLE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRegionVariable : public Fix {
 public:
  FixRegionVariable(class LAMMPS *, int, char **);
  ~FixRegionVariable();

  int setmask();
  void init();

  void write_restart(FILE *);
  void restart(char *);

  class Region* region();

 protected:
  int iarg;

  int step_start;
  double dt;

  int n_regions;
  double *steps;
  double *times;
  class Region **regions;
};

}

#endif
#endif
