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

#include "pointers.h"
#include "tri_mesh.h"
#include "scalar_container.h"
#include "fix_property_atom.h"
#include "error.h"
#include "fix_mesh_surface.h"
#include <string>

#ifndef LMP_MESH_MODULE_H
#define LMP_MESH_MODULE_H

namespace LAMMPS_NS
{
    class FixMeshSurface;

    class MeshModule : protected Pointers
    {
    public:
        ~MeshModule();

        virtual void post_create_pre_restart() {}
        virtual void post_create() {}
        virtual void init() {}
        virtual void setup(int vflag) {}
        virtual void setup_pre_force(int vflag) {}
        virtual void pre_force(int vflag) {}
        virtual void initial_integrate(int vflag) {}
        virtual void final_integrate_pre_comm() {}
        virtual void final_integrate() {}
        virtual void end_of_step() {}
        virtual double compute_vector(int n) { return 0.0; }
        virtual void add_particle_contribution(int ip, double *frc, double *delta, int iTri, double *v_wall) {}
        virtual int modify_param(int narg, char **arg) { return 0; }

        virtual int  setmask()
        { return 0; }

        TriMesh* get_mesh()
        { return mesh; }

        virtual int get_num_vector_components() const
        { return 0; }

    protected:
        MeshModule(LAMMPS *lmp, int &iarg_, int narg, char **arg, FixMeshSurface *fix_mesh_);
        FixMeshSurface *fix_mesh;
        TriMesh *mesh;
    };
}

#endif
