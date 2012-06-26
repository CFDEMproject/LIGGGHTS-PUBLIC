/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
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

        void setup(int vflag);
        void setup_pre_force(int vflag);
        void pre_delete(bool unfixflag);
        int setmask();

        void initial_integrate(int);
        void pre_force(int vflag);
        void final_integrate();
        void pre_neighbor();

        void write_restart(FILE *);
        void restart(char *);

     protected:
        class FixMesh* fixMesh()
        { return fix_mesh_; }

      private:

        bool decide_rebuild();
        void store_node_pos();

        class FixMesh *fix_mesh_;
        class MeshMover *move_;
        class AbstractMesh *mesh_;

        bool neighListFresh_;

        double time_, time_since_setup_;
        class MultiVectorContainer<double,3,3> oldNodes;
  };
} /* namespace LAMMPS_NS */
#endif
#endif
