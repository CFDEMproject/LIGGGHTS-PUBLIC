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

    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Philippe Seil (JKU Linz)
    Tóth János (MATE, Gödöllő)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

#else

#ifndef LMP_FIX_MESH_H
#define LMP_FIX_MESH_H

#include "fix_base_liggghts.h"
#include "fix_move_mesh.h"
#include "abstract_mesh.h"
#include <list>
#include "input_mesh_tri.h"

namespace LAMMPS_NS
{
    class FixMesh : public FixBaseLiggghts
    {
      public:

        FixMesh(LAMMPS *lmp, int narg, char **arg);
        virtual ~FixMesh();

        virtual void post_create();
        virtual void pre_delete(bool unfixflag);

        virtual void init();
        virtual void setup(int vflag);

        virtual int setmask();
        void setup_pre_force(int);

        void write_restart(FILE *fp);
        void restart(char *buf);

        virtual void pre_exchange();
        virtual void pre_force(int);
        virtual void final_integrate();

        void box_extent(double &xlo,double &xhi,double &ylo,double &yhi,double &zlo,double &zhi);

        int min_type();
        int max_type();

        void resetNodePosOrig(const FixMoveMesh * const caller);
        void move(const double * const dx, const FixMoveMesh * const caller);
        void rotate(const double dphi, const double * const axis, const double * const center, const FixMoveMesh * const caller);

        class AbstractMesh* mesh()
        { return mesh_; }

        virtual bool surfaceVel()
        { return false; }

        inline bool trackPerElementTemp()
        { return trackPerElementTemp_; }

        bool manipulated()
        { return manipulated_; }

        bool verbose()
        { return verbose_; }

        void register_move(FixMoveMesh * toInsert)
        { fixMoveMeshes_.push_back(toInsert); }

        void unregister_move(const FixMoveMesh * const toDelete)
        {
            std::list<FixMoveMesh *>::iterator it;
            for (it = fixMoveMeshes_.begin(); it != fixMoveMeshes_.end(); it++)
            {
                if (*it == toDelete)
                {
                    fixMoveMeshes_.erase(it);
                    return;
                }
            }
        }

      protected:

        // mesh manipulation upon creation
        virtual void moveMesh(double const dx, double const dy, double const dz);

        virtual void rotateMesh(double const axisX, double const axisY, double const axisZ, double const phi);
        virtual void scaleMesh(double const factor);

        void create_mesh(const char *mesh_file_or_generator_type, bool is_fix);
        void create_mesh_restart(const char *mesh_file_or_generator_type);
        InputMeshTri::GeneratedType generator_type(const char * type);

        int iarg_;

        int atom_type_mesh_;

        double temperature_mesh_;
        double mass_temperature_;

        bool trackPerElementTemp_;

      private:

        void handle_exclusion_list();

        void initialSetup();

        // mesh object
        class AbstractMesh *mesh_;

        // flag if mesh is already setup
        bool setupFlag_;

        // decides if parallel operations needed for
        // mesh on this time-step
        bool pOpFlag_;

        bool manipulated_;

        // flags and params to be passed to the mesh
        bool verbose_,autoRemoveDuplicates_;

        // mesh precision
        double precision_;

        // ignore features smaller than this size
        double min_feature_length_;

        // mesh correction
        FILE *element_exclusion_list_;
        bool read_exclusion_list_;
        int *exclusion_list_;
        int size_exclusion_list_;

        class FixPropertyGlobal *fix_capacity_;

        std::list<FixMoveMesh *> fixMoveMeshes_;

        // this friend class needs access to the moveMesh function
        friend class MeshModuleStress6DOF;
        friend class MeshModuleStress6DOFexternal;

        void parse_and_consume_generator_args(int *start_arg, int narg, char **args);
        bool is_valid_generator_type(const char * type);

        InputMeshTri::GeneratorParameters gen_params_;
    };

} /* namespace LAMMPS_NS */

#endif
#endif
