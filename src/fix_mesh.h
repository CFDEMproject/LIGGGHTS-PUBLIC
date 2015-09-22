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

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
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

        virtual void init();
        virtual void setup(int vflag) {}

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

        class AbstractMesh* mesh()
        { return mesh_; }

        virtual bool surfaceVel()
        { return false; }

        bool manipulated()
        { return manipulated_; }

        bool verbose()
        { return verbose_; }

      protected:

        // mesh manipulation upon creation
        virtual void moveMesh(double const dx, double const dy, double const dz);
        virtual void rotateMesh(double const axisX, double const axisY, double const axisZ, double const phi);
        virtual void scaleMesh(double const factor);

        void create_mesh(char *mesh_fname);
        void create_mesh_restart();

        int iarg_;

        int atom_type_mesh_;

        double mass_temperature_;

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

        // mesh correction
        FILE *element_exclusion_list_;
        bool read_exclusion_list_;
        int *exclusion_list_;
        int size_exclusion_list_;

        class FixPropertyGlobal *fix_capacity_;
  };

} /* namespace LAMMPS_NS */

#endif
#endif
