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
    (if not contributing author is listed, this file has been contributed
    by the core developer)

    Arno Mayrhofer (DCS Computing GmbH)

    Copyright 2016-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifdef MESHMODULE_CLASS

MeshModuleStyle(liquidtransfer,MeshModuleLiquidTransfer)

#else

#ifndef LMP_MESH_MODULE_LIQUIDTRANSFER_H
#define LMP_MESH_MODULE_LIQUIDTRANSFER_H

#include "mesh_module.h"

namespace LAMMPS_NS
{
    class MeshModuleLiquidTransfer : public MeshModule
    {
      public:

        MeshModuleLiquidTransfer(LAMMPS *lmp, int &iarg_, int narg, char **arg, FixMeshSurface *fix_mesh);
        ~MeshModuleLiquidTransfer();

        void post_create_pre_restart();
        virtual void post_create();

        virtual void init();
        virtual void setup(int vflag);
        virtual int setmask();

        virtual void setup_pre_force(int vflag);

        virtual void pre_force(int vflag);
        virtual void final_integrate();

        void add_source_contribution(const int iTri, const double liquidTransferred);

        void add_liquid_flux(const int iTri, const double liquidTransferred, const bool limit, const double maxLiquid);

        double get_wall_thickness()
        { return wall_thickness_; }

      private:

        void update_mesh_liquid_content();

        inline double& liquid_content(int i)
        { return (*liquid_content_)(i); }

        inline double& liquid_flux(int i)
        { return (*liquid_flux_)(i); }

        // per-element Temperature and flux
        ScalarContainer<double> *liquid_content_;
        ScalarContainer<double> *liquid_flux_;

        // per-particle temp and flux
        class FixPropertyAtom* fix_temp_;
        class FixPropertyAtom* fix_flux_;

        // check if wall liquid volume is to be limited
        bool limit_liquid_content_;
        double max_liquid_content_;

        // wall thickness and density
        double wall_thickness_;
        double initial_liquid_content_;

    };

} /* namespace LAMMPS_NS */

#endif
#endif
