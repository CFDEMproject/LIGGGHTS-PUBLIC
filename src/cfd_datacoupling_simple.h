/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2012-2014 DCS Computing GmbH, Linz
   Copyright 2012-2014 Graz University of Technology

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#ifndef LMP_CFD_DATACOUPLING_SIMPLE_H
#define LMP_CFD_DATACOUPLING_SIMPLE_H

#include "cfd_datacoupling.h"

namespace LAMMPS_NS {

class CfdDatacouplingSimple : public CfdDatacoupling {
 public:

  CfdDatacouplingSimple(class LAMMPS *lmp, int jarg, int narg, char **arg, void* fix)
    : CfdDatacoupling(lmp,jarg,narg,arg,NULL) {}

  virtual ~CfdDatacouplingSimple() {}
  void exchange() {};

  void* find_pull_property(const char *name, const char *type, int &len1, int &len2)
  { return CfdDatacoupling::find_pull_property(name,type,len1,len2); }

  void* find_push_property(const char *name, const char *type, int &len1, int &len2)
  { return CfdDatacoupling::find_push_property(name,type,len1,len2); }
};

}

#endif

