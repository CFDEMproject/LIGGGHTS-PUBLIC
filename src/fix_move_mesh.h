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

    Christoph Kloss (DCS Computing GmbH, Linz, JKU Linz)
    Philippe Seil (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(move/mesh,FixMoveMesh)
FixStyle(move/mesh/gran,FixMoveMesh) // for backward compatibility

#else

#ifndef FIX_MOVE_MESH_H
#define FIX_MOVE_MESH_H

#include "fix.h"
#include "container.h"

namespace LAMMPS_NS
{

  class FixMoveMesh : public LAMMPS_NS::Fix
  {
      public:

        FixMoveMesh(LAMMPS *lmp, int narg, char **arg);
        virtual ~FixMoveMesh();

        void post_create();

        void setup(int vflag);
        void pre_delete(bool unfixflag);
        int setmask();

        void initial_integrate(int);
        void final_integrate();

        void write_restart(FILE *);
        void restart(char *);

        void add_reference_point(double *point);
        void get_reference_point(double *point);
        void reset_reference_point();

        class AbstractMesh * mesh()
        { return mesh_; }

     protected:

        class FixMesh* fixMesh()
        { return fix_mesh_; }

      private:

        class FixMesh *fix_mesh_;
        class MeshMover *move_;
        class AbstractMesh *mesh_;

        double time_;
        double time_since_setup_;

        double reference_point_[3];
  };
} /* namespace LAMMPS_NS */
#endif
#endif
