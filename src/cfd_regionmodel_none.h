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

#ifndef LMP_CFD_REGIONMODEL_NONE_H
#define LMP_CFD_REGIONMODEL_NONE_H

#include "cfd_regionmodel.h"

namespace LAMMPS_NS {

class CfdRegionmodelNone : public CfdRegionmodel {
 public:
  CfdRegionmodelNone(class LAMMPS *, int, int, char **,class FixCfdCoupling* fc);
  ~CfdRegionmodelNone();
  void init();

 protected:
  class FixPropertyAtom *inRegion;
  class FixPropertyGlobal *outRegion;

  double *inregion;
  double *outregion;
  int nout;

  int nlocal_last;

  virtual void special_settings();
  virtual void rm_update();
};

}

#endif
