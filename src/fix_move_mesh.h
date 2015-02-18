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
