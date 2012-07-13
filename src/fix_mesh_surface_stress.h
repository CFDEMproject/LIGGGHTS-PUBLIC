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

FixStyle(mesh/surface/stress,FixMeshSurfaceStress)

#else

#ifndef LMP_FIX_SURFACE_MESH_STRESS_H
#define LMP_FIX_SURFACE_MESH_STRESS_H

#include "fix_mesh_surface.h"

namespace LAMMPS_NS
{
  class FixMeshSurfaceStress : public FixMeshSurface
  {
      public:

        FixMeshSurfaceStress(LAMMPS *lmp, int narg, char **arg);
        virtual ~FixMeshSurfaceStress();

        virtual void post_create();

        virtual void init();
        virtual int setmask();

        void pre_force(int vflag);
        void final_integrate();

        double compute_vector(int n);

        void add_particle_contribution(int ip, double *frc, double *delta,
                                       int iTri, double *v_wall);

        // inline access

        inline bool trackWear()
        { return wear_flag_; }

        inline void f_total(double *_f)
        { vectorCopy3D(f_total_,_f); }

        inline void torque_total(double *_t)
        { vectorCopy3D(torque_total_,_t); }

      protected:

        // inline access

        inline double* f(int i)
        { return (*f_)(i); }

        inline double& sigma_n(int i)
        { return (*sigma_n_)(i); }

        inline double& sigma_t(int i)
        { return (*sigma_t_)(i); }

        inline double& wear(int i)
        { return (*wear_)(i); }

        inline double& wear_step(int i)
        { return (*wear_step_)(i); }

      private:

        // inititalization fcts
        void initStress();
        void initWear();

        void calc_total_force();

        // STRESS

        // stress flag in FixMeshSurface
        // reference point, total force and total torque
        double p_ref_[3];
        double f_total_[3], torque_total_[3];

        // per-element force and torque
        VectorContainer<double,3> *f_;
        ScalarContainer<double> *sigma_n_;
        ScalarContainer<double> *sigma_t_;

        // WEAR

        // flag for wear model and Finnie constant
        int wear_flag_;
        double const* const* k_finnie_;
        ScalarContainer<double> *wear_;
        ScalarContainer<double> *wear_step_;
  };

} /* namespace LAMMPS_NS */

#endif /* LMP_FIX_MESH_H */
#endif /* FIX_CLASS */
