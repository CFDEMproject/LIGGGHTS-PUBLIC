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

    virtual void post_create_pre_restart();
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
    int modify_param(int, char **);
    void resetIntegrator() {sum_err_ = 0;}

    // properties of mesh

    VectorContainer<double,3> &xcm_;
    VectorContainer<double,3> &vcm_;
    VectorContainer<double,3> &omegacm_;
    VectorContainer<double,3> &xcm_orig_;

    // position and velocity for each node

    double*** nodes_;
    MultiVectorContainer<double,3,3> *v_;

    // servo settings and controller

    double axis_[3],totalPhi_;
    double *ctrl_op_,*pv_vec_;
    double vel_max_,vel_min_,ctrl_op_max_,ctrl_op_min_,ratio_;
    double sp_mag_,sp_mag_inv_;
    double pv_mag_,old_pv_mag_;
    double err_, sum_err_;
    double kp_,ki_,kd_;

    // variable set point
    int sp_var_, sp_style_;
    char *sp_str_;

    // flags
    bool int_flag_;
    bool mode_flag_;
    int ctrl_style_;

    // timesteps
    double dtf_,dtv_;

    // for area calculation
    class ModifiedAndrew *mod_andrew_;

    // signum function
    template <typename T> int sgn(T val) {
      return (T(0) < val) - (val < T(0));
    }

  }; //end class

}

#endif
#endif
