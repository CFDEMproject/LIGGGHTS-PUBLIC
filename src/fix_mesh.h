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

      private:

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
  };

} /* namespace LAMMPS_NS */

#endif
#endif
