/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:

    Richard Berger (JKU Linz)

    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef GRANULAR_WALL_FIX_H
#define GRANULAR_WALL_FIX_H

#include "utils.h"
#include "tri_mesh.h"
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
    virtual void compute_force(FixWallGran * fwg, ContactModels::SurfacesIntersectData & sidata, double * v_wall,class TriMesh *mesh = 0,int iTri = 0) = 0;
  };

  class Factory : public Utils::AbstractFactory<IGranularWall> {
    Factory() {}
  public:
    static Factory & instance();
  };

}

}

#endif // GRANULAR_WALL_FIX_H
