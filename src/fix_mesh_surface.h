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

FixStyle(mesh/surface,FixMeshSurface)
FixStyle(mesh/surface/planar,FixMeshSurface)

#else

#ifndef LMP_FIX_SURFACE_MESH_H
#define LMP_FIX_SURFACE_MESH_H

#include "fix_mesh.h"
#include "tri_mesh.h"
#include "fix_contact_history_mesh.h"
#include "fix_contact_property_atom.h"
#include "fix_neighlist_mesh.h"
#include "custom_value_tracker.h"

namespace LAMMPS_NS
{
  class FixMeshSurface : public FixMesh
  {
      public:

        FixMeshSurface(LAMMPS *lmp, int narg, char **arg);
        virtual ~FixMeshSurface();

        virtual void post_create();
        virtual void pre_delete(bool unfixflag);

        virtual void init() {}
        virtual void setup(int vflag) {}

        virtual int setmask();
        virtual void setup_pre_force(int);

        virtual void pre_force(int);
        virtual void final_integrate();

        virtual void createWallNeighList(int igrp);
        virtual class FixNeighlistMesh* createOtherNeighList(int igrp,const char *nId);
        void createContactHistory(int dnum);
        void createMeshforceContact();

        void deleteWallNeighList();
        void deleteContactHistory();
        void deleteMeshforceContact();

        inline bool trackStress()
        {return stress_flag_;}

        inline int atomTypeWall()
        { return atom_type_mesh_;}

        inline class FixContactHistoryMesh* contactHistory()
        { return fix_contact_history_mesh_;}

        inline class FixNeighlistMesh* meshNeighlist()
        { return fix_mesh_neighlist_;}

        inline bool hasNeighList()
        { return fix_mesh_neighlist_?true:false; }

        inline class FixContactPropertyAtomWall* meshforceContact()
        { return fix_meshforce_contact_;}

        inline class TriMesh *triMesh()
        { return static_cast<TriMesh*>(mesh()); }

        bool surfaceVel()
        { return (velFlag_ || angVelFlag_); }

        void dumpAdd()
        { n_dump_active_++; }

        void dumpRemove()
        { n_dump_active_--; }

      protected:

        class FixContactHistoryMesh *fix_contact_history_mesh_;
        class FixNeighlistMesh *fix_mesh_neighlist_;
        class FixContactPropertyAtomWall *fix_meshforce_contact_;

        // flag for stressanalysis
        bool stress_flag_;

      private:

        void initVel();
        void setVel();
        void initAngVel();
        void setAngVel();

        // surface velocity
        bool velFlag_;
        double vSurf_[3];

        // rotational surf vel
        bool angVelFlag_;
        double origin_[3];
        double axis_[3];
        double omegaSurf_;

        // flag if dump mesh/* or dump stl is active
        int n_dump_active_;

        // mesh curvature
        double curvature_;
  };

} /* namespace LAMMPS_NS */

#endif /* LMP_FIX_MESH_H */
#endif /* FIX_CLASS */
