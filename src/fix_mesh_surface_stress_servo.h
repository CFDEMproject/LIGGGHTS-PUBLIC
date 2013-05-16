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
   Andreas Aigner (JKU Linz)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(mesh/surface/stress/servo,FixMeshSurfaceStressServo)

#else

#ifndef LMP_FIX_MESH_SURFACE_STRESS_SERVO_H
#define LMP_FIX_MESH_SURFACE_STRESS_SERVO_H

#include "fix.h"
#include "input.h"
#include "math.h"
#include "fix_mesh_surface_stress.h"

namespace LAMMPS_NS {

class FixMeshSurfaceStressServo : public FixMeshSurfaceStress {

    public:

      FixMeshSurfaceStressServo(class LAMMPS *, int, char **);
      virtual ~FixMeshSurfaceStressServo();

      virtual void post_create();

      void init();
      int setmask();

      virtual void setup_pre_force(int vflag);
      void initial_integrate(int vflag);
      void add_particle_contribution(int ip, double *frc,
                            double *delta, int iTri, double *v_wall);
      void final_integrate();

      void reset_dt();
      double compute_vector(int n);

    private:

      void init_defaults();
      void error_checks();

      void limit_vel();
      void update_mass();
      void set_v_node();
      void set_v_node_rotate();
      double getMaxRad();

      int dim_;
      double axis_[3],totalPhi_;

      // properties of mesh

      VectorContainer<double,3> &xcm_;
      VectorContainer<double,3> &vcm_;
      VectorContainer<double,3> &omegacm_;
      VectorContainer<double,3> &xcm_orig_;

      // servo settings and controller

      bool mode_flag_;
      double vel_max_,vel_min_,set_point_,set_point_inv_,ctrl_output_max_,ctrl_output_min_,ratio_;
      char *sp_str_;
      int sp_var_, sp_style_;
      double *control_output_,*process_value_;
      int pv_flag_;
      double err_, sum_err_, old_process_value_;
      double kp_,ki_,kd_;

      // timesteps and flags for integration

      double dtf_,dtv_;

      bool int_flag_;
      int modify_param(int, char **);

      // position and velocity for each node

      double*** nodes_;
      MultiVectorContainer<double,3,3> &v_;

      // signum function
      template <typename T> int sgn(T val) {
          return (T(0) < val) - (val < T(0));
      }

      // for area calculation
      class ModifiedAndrew *mod_andrew_;

}; //end class

}

#endif
#endif
