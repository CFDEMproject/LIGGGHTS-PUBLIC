/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2015-     DCS Computing GmbH, Linz

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

#ifndef LMP_PROPERTIES_H
#define LMP_PROPERTIES_H

#include "lammps.h"
#include "pointers.h"

namespace LAMMPS_NS {

class Properties: protected Pointers
{
 public:

  Properties(LAMMPS *lmp);
  ~Properties();

  int max_type();
  void* find_property(const char *name, const char *type, int &len1, int &len2);
  inline class MultisphereParallel *ms_data() { return ms_data_;}

 private:

  // multisphere
  class FixMultisphere *ms_;
  class MultisphereParallel *ms_data_;

  int mintype,maxtype;
}; //end class

}

#endif
