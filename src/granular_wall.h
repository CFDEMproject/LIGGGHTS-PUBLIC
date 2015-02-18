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
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */
#ifndef GRANULAR_WALL_FIX_H
#define GRANULAR_WALL_FIX_H

#include "utils.h"
#include "contact_interface.h"

namespace LAMMPS_NS {
  class FixWallGran;
}

namespace LIGGGHTS {
using namespace LAMMPS_NS;

namespace Walls {
  class IGranularWall {
  public:
    typedef FixWallGran ParentType;
    virtual ~IGranularWall();
    virtual void settings(int nargs, char ** args) = 0;
    virtual void init_granular() = 0;
    virtual void compute_force(FixWallGran * fwg, ContactModels::CollisionData & cdata, double * v_wall) = 0;
  };

  class Factory : public Utils::AbstractFactory<IGranularWall> {
    Factory() {}
  public:
    static Factory & instance();
  };

}

}

#endif // GRANULAR_WALL_FIX_H
