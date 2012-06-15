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

FixStyle(mesh,FixMesh)

// for backward compatibility
FixStyle(mesh/gran,FixMesh)

#else

#ifndef LMP_FIX_MESH_H
#define LMP_FIX_MESH_H

#include "fix.h"
#include "tri_mesh.h"
#include "fix_contact_history.h"
#include "fix_neighlist_mesh.h"
#include "custom_value_tracker.h"

namespace LAMMPS_NS
{
  class FixMesh : public Fix
  {
      //friend class FixContactTracker;

      public:

        FixMesh(LAMMPS *lmp, int narg, char **arg);
        virtual ~FixMesh();

        virtual void post_create();
        virtual void pre_delete(bool unfixflag);

        virtual int setmask();
        void setup_pre_force(int);

        int min_type();
        int max_type();

        virtual void pre_exchange();
        virtual void pre_force(int);

        void createNeighList();
        void createContactHistory(int dnum);

        int atomTypeWall()
        { return atom_type_wall_;}

        class FixContactHistory* contactHistory()
        { return fix_contact_history_;}

        class FixNeighlistMesh* meshNeighlist()
        { return fix_mesh_neighlist_;}

        class TriMesh *mesh()
        { return mesh_; }

        bool surfaceVel()
        { return surfaceVel_; }

      protected:

        int iarg_;

        class TriMesh *mesh_;

        bool surfaceVel_;

        int atom_type_wall_;
        class FixContactHistory *fix_contact_history_;
        class FixNeighlistMesh *fix_mesh_neighlist_;

      private:

        void initialSetup();

        void moveMesh(double const dx, double const dy, double const dz);
        void rotateMesh(double const axisX, double const axisY, double const axisZ, double const phi);
        void scaleMesh(double const factor);

        void initVel(double *conv_vel);
        void initAngVel(double *origin,double *axis,double omega);

        // flag if mesh is already setup
        bool setupFlag_;

        // decides if parallel operations needed for
        // mesh on this time-step
        bool pOpFlag_;
  };

} /* namespace LAMMPS_NS */

#endif /* LMP_FIX_MESH_H */
#endif /* FIX_CLASS */
