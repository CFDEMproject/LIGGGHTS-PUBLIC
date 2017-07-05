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

    Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
    Philippe Seil   (JKU Linz)
    Arno Mayrhofer  (CFDEMresearch GmbH, DCS Computing GmbH, Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
    Copyright 2016-     CFDEMresearch GmbH, Linz
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
#include "mesh_module.h"
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <cstddef>

namespace LAMMPS_NS
{
  class FixMeshSurface : public FixMesh
  {
      public:

        FixMeshSurface(LAMMPS *lmp, int narg, char **arg);
        virtual ~FixMeshSurface();

        virtual void post_create_pre_restart();
        virtual void post_create();
        virtual void pre_delete(bool unfixflag);

        virtual void init();
        virtual void setup(int vflag);

        virtual int setmask();
        virtual void setup_pre_force(int);

        virtual void pre_force(int);
        virtual void initial_integrate(int);
        virtual void final_integrate();
        virtual void end_of_step();

        virtual int modify_param(int narg, char **arg);

        virtual double compute_vector (int);

        virtual void createWallNeighList(int igrp);
        virtual class FixNeighlistMesh* createOtherNeighList(int igrp,const char *nId);
        void createContactHistory(int dnum);
        void createMeshforceContact();

        void createMeshforceContactStress();
        void createMulticontactData();

        void deleteWallNeighList();
        void deleteOtherNeighList(const char *nId);
        void deleteAllOtherNeighList();
        void deleteContactHistory();
        void deleteMeshforceContact();

        void deleteMeshforceContactStress();
        void deleteMeshMulticontactData();

        MeshModule* get_module(std::string name);
        void add_particle_contribution(int ip, double *frc, double *delta, int iTri, double *v_wall);

        bool trackStress();

        inline int atomTypeWall()
        { return atom_type_mesh_;}

        inline class FixContactHistoryMesh* contactHistory()
        { return fix_contact_history_mesh_;}

        inline class FixNeighlistMesh* meshNeighlist()
        { return fix_mesh_neighlist_;}

        inline bool hasNeighList()
        { return fix_mesh_neighlist_?true:false; }

        inline bool hasOtherNeighList()
        { return list_other_neighlist_.size() > 0; }

        inline class FixContactPropertyAtomWall* meshforceContact()
        { return fix_meshforce_contact_;}

        inline class FixContactPropertyAtomWall* meshforceContactStress()
        { return fix_meshforce_contact_stress_;}

        inline class FixContactPropertyAtomWall* meshMulticontactData()
        { return fix_mesh_multicontact_data_;}

        inline class TriMesh *triMesh()
        { return static_cast<TriMesh*>(mesh()); }

        bool surfaceVel()
        { return (velFlag_ || angVelFlag_); }

        void dumpAdd()
        { n_dump_active_++; }

        void dumpRemove()
        { n_dump_active_--; }

        inline int get_groupbit() const
        { return groupbit; }

        inline int get_nevery() const
        { return nevery; }

        inline bigint* get_next_reneighbor_ptr()
        { return &next_reneighbor; }

        void set_global_freq(const unsigned int global_freq_)
        { global_freq = global_freq_; }

        void set_extvector(const unsigned int extvector_)
        { extvector = extvector_; }

        void set_time_depend(const unsigned int time_depend_)
        { time_depend = time_depend_; }

        void set_nevery(const unsigned int nevery_)
        { nevery = nevery_; }

      protected:

        int getCreateMeshTriCount()
        { return extrusion_tri_count_; }

        double * getCreateMeshTriNode(const int i)
        { return &(extrusion_tri_nodes_[i*3]); }

        class FixContactHistoryMesh *fix_contact_history_mesh_;
        class FixNeighlistMesh *fix_mesh_neighlist_;
        class FixContactPropertyAtomWall *fix_meshforce_contact_;
        class FixContactPropertyAtomWall *fix_meshforce_contact_stress_;
        class FixContactPropertyAtomWall *fix_mesh_multicontact_data_;

        // list of neighbor lists created by createOtherNeighList
        std::map<std::string, FixNeighlistMesh*> list_other_neighlist_;

      private:

        void initVel();
        void setVel();
        void initAngVel();
        void setAngVel();

        // surface velocity
        bool velFlag_;
        double vSurf_[3];
        char *vSurfStrX_, *vSurfStrY_, *vSurfStrZ_;
        int vSurfVarX_, vSurfVarY_, vSurfVarZ_;
        int vSurfStyleX_, vSurfStyleY_, vSurfStyleZ_;

        // rotational surf vel
        bool angVelFlag_;
        double origin_[3];
        double axis_[3];
        double omegaSurf_;
        char *omegaStr_;
        int omegaVar_, omegaStyle_;

        // flag if dump mesh/* or dump stl is active
        int n_dump_active_;

        // mesh curvature
        double curvature_;
        bool curvature_tolerant_;

        // extrude mesh
        bool extrude_mesh_;
        double extrusion_length_;
        int extrusion_tri_count_;
        double *extrusion_tri_nodes_;
        bool extrusion_created_;

        // active mesh modules
        typedef MeshModule *(*MeshModuleCreator)(LAMMPS *lmp, int &iarg_, int narg, char **arg, FixMeshSurface *fix_mesh);
        std::map<std::string, MeshModule*> active_mesh_modules;
        std::vector<std::string> mesh_module_order;

        template <typename T> static MeshModule *meshmodule_creator(LAMMPS *lmp, int &iarg_, int narg, char **arg, FixMeshSurface *fix_mesh);
  };

} /* namespace LAMMPS_NS */

#endif /* LMP_FIX_MESH_H */
#endif /* FIX_CLASS */
