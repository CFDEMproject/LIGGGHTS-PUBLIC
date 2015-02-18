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

#ifndef LMP_CFD_REGIONMODEL_H
#define LMP_CFD_REGIONMODEL_H

#include "pointers.h"

namespace LAMMPS_NS {

class CfdRegionmodel : protected Pointers {
 public:
  CfdRegionmodel(class LAMMPS *lmp, int jarg, int narg, char **arg,class FixCfdCoupling* fc) : Pointers(lmp)
  {
    UNUSED(narg);
    UNUSED(jarg);
    UNUSED(arg);
      this->fc = fc;
  }
  ~CfdRegionmodel() {}

  int get_iarg() {return iarg;}
  bool liggghts_is_active;

  virtual void init() {};
  virtual void rm_update() {};

 protected:
  int iarg;
  class FixCfdCoupling *fc;

  void add_push_property(const char *name, const char *type)
  {
     fc->add_push_property(name,type);
  }

  void add_pull_property(const char *name, const char *type)
  {
     fc->add_pull_property(name,type);
  }
};

}

#endif
