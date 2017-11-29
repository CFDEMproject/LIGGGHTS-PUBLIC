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

    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Andreas Aigner (DCS Computing GmbH, JKU Linz)
    Arno Mayrhofer (DCS Computing GmbH)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifdef MESHMODULE_CLASS

MeshModuleStyle(servo,MeshModuleStressServo)

#else

#ifndef LMP_MESH_MODULE_STRESS_SERVO_H
#define LMP_MESH_MODULE_STRESS_SERVO_H

#include "fix.h"
#include "input.h"
#include <cmath>
#include "mesh_module_stress.h"
#include "mesh_module.h"

namespace LAMMPS_NS
{

  class MeshModuleStressServo : public MeshModule
  {

  public:

    MeshModuleStressServo(LAMMPS *lmp, int &iarg_, int narg, char **arg, FixMeshSurface *fix_mesh);
    virtual ~MeshModuleStressServo();

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

    inline int get_num_vector_components() const
    { return 3; }

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

    MeshModuleStress *mm_stress;

  }; //end class

}

#endif
#endif
