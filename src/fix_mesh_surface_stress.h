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

        virtual void post_create_pre_restart();
        virtual void post_create();

        virtual void init();
        virtual void setup(int vflag) {}
        virtual int setmask();

        virtual void setup_pre_force(int vflag);

        virtual void pre_force(int vflag);
        virtual void final_integrate();

        virtual double compute_vector(int n);

        virtual void add_particle_contribution(int ip, double *frc,
                            double *delta, int iTri, double *v_wall);

        void add_global_external_contribution(double *frc);
        void add_global_external_contribution(double *frc,double *trq);

        // inline access

        inline bool trackWear()
        { return wear_flag_; }

        inline bool trackStress()
        { return stress_flag_; }

        inline void f_total(double *_f)
        { vectorCopy3D(f_total_,_f); }

        inline double f_total(int i)
        { return f_total_[i]; }

        inline double f_total_mag()
        { return vectorMag3D(f_total_); }

        inline void torque_total(double *_t)
        { vectorCopy3D(torque_total_,_t); }

        inline double torque_total(int i)
        { return torque_total_[i]; }

      protected:

        // STRESS
        // total force and total torque
        double f_total_[3], torque_total_[3]; 

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

        inline void set_p_ref(double *_p_ref)
        { p_ref_.set(0,_p_ref); }

        inline double p_ref(int i)
        { return p_ref_(0)[i]; }

      private:

        // inititalization fcts
        void regStress();
        void regWear();
        void zeroizeStress();
        void zeroizeWear();

        void calc_total_force();
        void add_gravity();

        // STRESS

        // stress flag in FixMeshSurface
        // reference point, total force and total torque
        VectorContainer<double,3> &p_ref_;

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

#endif
#endif
