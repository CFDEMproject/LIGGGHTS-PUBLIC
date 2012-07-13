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

#else

#ifndef LMP_FIX_MESH_H
#define LMP_FIX_MESH_H

#include "fix.h"

namespace LAMMPS_NS
{
  class FixMesh : public Fix
  {
      public:

        FixMesh(LAMMPS *lmp, int narg, char **arg);
        virtual ~FixMesh();

        virtual void post_create();
        virtual void pre_delete(bool unfixflag);

        virtual void init() {}

        virtual int setmask();
        void setup_pre_force(int);

        void write_restart(FILE *fp);
        void restart(char *buf);

        virtual void pre_exchange();
        virtual void pre_force(int);
        virtual void final_integrate();

        int min_type();
        int max_type();

        class AbstractMesh* mesh()
        { return mesh_; }

        virtual bool surfaceVel()
        { return false; }

      protected:

        void create_mesh(char *mesh_fname);
        void create_mesh_restart();

        int iarg_;

        class AbstractMesh *mesh_;

        int atom_type_mesh_;

      private:

        void initialSetup();

        // mesh manipulation upon creation
        void moveMesh(double const dx, double const dy, double const dz);
        void rotateMesh(double const axisX, double const axisY, double const axisZ, double const phi);
        void scaleMesh(double const factor);

        // flag if mesh is already setup
        bool setupFlag_;

        // decides if parallel operations needed for
        // mesh on this time-step
        bool pOpFlag_;
  };

} /* namespace LAMMPS_NS */

#endif
#endif
