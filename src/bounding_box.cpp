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

/* ----------------------------------------------------------------------
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Philippe Seil (JKU Linz)
------------------------------------------------------------------------- */

#include "bounding_box.h"

namespace LAMMPS_NS
{

  BoundingBox::BoundingBox()
  : xLo(0.), xHi(0.), yLo(0.), yHi(0.), zLo(0.), zHi(0.), initGiven(false)
  {}
  BoundingBox::BoundingBox(double xLo_, double xHi_, double yLo_, double yHi_, double zLo_, double zHi_)
  : xLo(xLo_), xHi(xHi_), yLo(yLo_), yHi(yHi_), zLo(zLo_), zHi(zHi_), initGiven(true)
  {}

  BoundingBox::~BoundingBox()
  {}

  void BoundingBox::reset()
  {
    xLo = 0.; xHi = 0.;
    yLo = 0.; yHi = 0.;
    zLo = 0.; zHi = 0.;
    initGiven = false;
  }

} /* namespace LAMMPS_NS */
