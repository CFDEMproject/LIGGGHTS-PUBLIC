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

    Copyright 2012-         DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#include <stdlib.h>
#include <string.h>
#include <vector>
#include "mesh_module_stress_servo.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "variable.h"
#include "error.h"
#include "vector_liggghts.h"
#include "fix_property_global.h"
#include "modified_andrew.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using MODIFIED_ANDREW_AUX::Circle;

#define EPSILON 1.0e-7
#define BIG 1000000.

//#define MIN(a,b) ((a) < (b) ? (a) : (b))

// identifier for variable set point
// identifier for controlled process value
enum{NONE,CONSTANT,EQUAL,FORCE,TORQUE};

/* ---------------------------------------------------------------------- */

MeshModuleStressServo::MeshModuleStressServo(LAMMPS *lmp, int &iarg_, int narg, char **arg, FixMeshSurface *fix_mesh) :
    MeshModule(lmp, iarg_, narg, arg, fix_mesh),
    xcm_(*mesh->prop().addGlobalProperty< VectorContainer<double,3> > ("xcm","comm_none","frame_invariant","restart_yes")),
    vcm_(*mesh->prop().addGlobalProperty< VectorContainer<double,3> > ("vcm","comm_none","frame_invariant","restart_yes")),
    omegacm_(*mesh->prop().addGlobalProperty< VectorContainer<double,3> > ("torquecm","comm_none","frame_invariant","restart_yes")),
    xcm_orig_(*mesh->prop().addGlobalProperty< VectorContainer<double,3> > ("xcm_orig","comm_none","frame_invariant","restart_yes")),

    v_(0),

    totalPhi_( 0.),
    ctrl_op_( 0),
    pv_vec_( 0),

    vel_max_(  0.),
    vel_min_(  0.),
    ctrl_op_max_(0.),
    ctrl_op_min_(0.),
    ratio_(0.),

    sp_mag_(0.),
    sp_mag_inv_(0.),
    pv_mag_(0.),
    old_pv_mag_(0.),

    err_(               0.),
    sum_err_(       0.),
    kp_(             0.01),
    ki_(             0.),
    kd_(             0.),

    sp_var_(     -1),
    sp_style_( NONE),
    sp_str_(     NULL),

    int_flag_( true),
    mode_flag_(false),
    ctrl_style_(NONE),
    mod_andrew_(new ModifiedAndrew(lmp))
{
    mm_stress = static_cast<MeshModuleStress*>(fix_mesh->get_module("stress"));
    if (!mm_stress)
        error->one(FLERR,"Mesh module \"servo\" requires mesh modules \"stress\"");

    if(!mm_stress->trackStress())
        error->one(FLERR,"stress = 'on' required");

    if(fix_mesh->verbose() && fix_mesh->manipulated())
        error->warning(FLERR,"Mesh has been scaled, moved, or rotated.\n"
                                     "Please note that values for 'com', 'vel' refer to the scaled, moved, or rotated configuration");

    // override default from base
    fix_mesh->set_time_depend(1); 

    // set defaults

    init_defaults();

    // parse further args

    bool hasargs = true;
    while(iarg_ < narg && hasargs)
    {
        hasargs = false;
        if(strcmp(arg[iarg_],"com") == 0) {
            if (narg < iarg_+4) error->one(FLERR,"not enough arguments for 'com'");
            ++iarg_;
            double _com[3];
            _com[0] = force->numeric(FLERR,arg[iarg_++]);
            _com[1] = force->numeric(FLERR,arg[iarg_++]);
            _com[2] = force->numeric(FLERR,arg[iarg_++]);
            xcm_.add(_com);
            mm_stress->set_p_ref(xcm_(0));
            hasargs = true;
        } else if(strcmp(arg[iarg_],"ctrlPV") == 0) {
            if (narg < iarg_+2) error->one(FLERR,"not enough arguments for 'ctrlPV'");
            if           (strcmp(arg[iarg_+1],"force") == 0) ctrl_style_ = FORCE;
            else if (strcmp(arg[iarg_+1],"torque") == 0) ctrl_style_ = TORQUE;
            else error->one(FLERR,"only 'force', 'torque' are valid arguments for ctrlPV");
            iarg_ = iarg_ + 2;
            hasargs = true;
        } else if(strcmp(arg[iarg_],"vel_max") == 0) {
            if (narg < iarg_+2) error->one(FLERR,"not enough arguments for 'vel_max'");
            ++iarg_;
            vel_max_ = force->numeric(FLERR,arg[iarg_++]);
            if(vel_max_ <= 0.)
                error->one(FLERR,"vel_max > 0 required");
            hasargs = true;
        } else if(strcmp(arg[iarg_],"target_val") == 0) {
            if (narg < iarg_+2) error->one(FLERR,"not enough arguments for 'target_val'");
            ++iarg_;
            if (strstr(arg[iarg_],"v_") == arg[iarg_]) {
                const int n = strlen(&arg[iarg_][2]) + 1;
                sp_str_ = new char[n];
                strcpy(sp_str_,&arg[iarg_][2]);
                sp_style_ = EQUAL;
            } else {
                sp_mag_ = -force->numeric(FLERR,arg[iarg_]); // the resultant force/torque/shear acts in opposite direction --> negative value
                if (sp_mag_ == 0.) error->one(FLERR,"'target_val' (desired force/torque) has to be != 0.0");
                sp_mag_inv_ = 1./fabs(sp_mag_);
                sp_style_ = CONSTANT;
            }
            ++iarg_;
            hasargs = true;
        } else if(strcmp(arg[iarg_],"axis") == 0) {
            if (narg < iarg_+4) error->one(FLERR,"not enough arguments for 'axis'");
            axis_[0] = force->numeric(FLERR,arg[iarg_+1]);
            axis_[1] = force->numeric(FLERR,arg[iarg_+2]);
            axis_[2] = force->numeric(FLERR,arg[iarg_+3]);
            // normalize axis
            vectorNormalize3D(axis_);
            iarg_ = iarg_+4;
            hasargs = true;
        } else if(strcmp(arg[iarg_],"kp") == 0) {
            if (narg < iarg_+2) error->one(FLERR,"not enough arguments");
            kp_ = force->numeric(FLERR,arg[iarg_+1]);
            iarg_ = iarg_+2;
            hasargs = true;
        } else if(strcmp(arg[iarg_],"ki") == 0) {
            if (narg < iarg_+2) error->one(FLERR,"not enough arguments");
            ki_ = force->numeric(FLERR,arg[iarg_+1]);
            iarg_ = iarg_+2;
            hasargs = true;
        } else if(strcmp(arg[iarg_],"kd") == 0) {
            if (narg < iarg_+2) error->one(FLERR,"not enough arguments");
            kd_ = force->numeric(FLERR,arg[iarg_+1]);
            iarg_ = iarg_+2;
            hasargs = true;
        } else if(strcmp(arg[iarg_],"mode") == 0) {
            if (narg < iarg_+2) error->one(FLERR,"not enough arguments");
            ++iarg_;
            if (strcmp("auto",arg[iarg_]) == 0) {
                mode_flag_ = true;
            }  else error->one(FLERR,"mode supports only auto");
            ++iarg_;
            hasargs = true;
        } else if(strcmp(arg[iarg_],"ratio") == 0) {
            if (narg < iarg_+2) error->one(FLERR,"not enough arguments");
            ratio_ = force->numeric(FLERR,arg[iarg_+1]);
            iarg_ = iarg_+2;
            hasargs = true;
        }
    }

    // store original position
    xcm_orig_.add(xcm_(0));
}

/* ---------------------------------------------------------------------- */

MeshModuleStressServo::~MeshModuleStressServo()
{
    delete [] sp_str_;
    delete mod_andrew_;
}

/* ---------------------------------------------------------------------- */

void MeshModuleStressServo::post_create_pre_restart()
{
    
    mesh->prop().addElementProperty< MultiVectorContainer<double,3,3> > ("v","comm_exchange_borders","frame_invariant","restart_no");
}

/* ---------------------------------------------------------------------- */

void MeshModuleStressServo::post_create()
{
    
    if (ctrl_style_ == FORCE)
        mesh->registerMove(false,true,false);
    else if(ctrl_style_ == TORQUE)
        mesh->registerMove(false,false,true);
    else
        error->one(FLERR,"Bad registration of upcoming move.");

    mesh->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v")->setAll(0.);
}

/* ---------------------------------------------------------------------- */

void MeshModuleStressServo::init_defaults()
{
    double zerovec[3] = {0., 0., 0.};
    vcm_.set(0,zerovec);
    omegacm_.set(0,zerovec);

    vectorZeroize3D(axis_);
}

/* ---------------------------------------------------------------------- */

void MeshModuleStressServo::error_checks()
{
    
    if(ctrl_style_ == NONE)
        error->one(FLERR,"please define 'ctrlPV' for the mesh");
    if(!xcm_.size())
        error->one(FLERR,"please define 'com' for the mesh");
    if(sp_style_ == CONSTANT && sp_mag_ == 0.)
        error->one(FLERR,"please define 'set_point' for the mesh");
    if(vel_max_ == 0.)
        error->one(FLERR,"please define 'vel_max' for the mesh");
    if(mode_flag_) {
        if(ratio_ == 0.)
            error->one(FLERR,"please define 'ratio' for the mesh, since you use the auto mode");
    } else {
        if(kp_ < 0. || ki_ < 0. || kd_ < 0.)
            error->one(FLERR,"kp, ki, and kd >= 0 required.");
        if(kp_ == 0. && ki_ == 0. && kd_ == 0.)
            error->one(FLERR,"kp, ki, and kd are zero. Please set a valid configuration");
    }

    if(mesh->nMove() > 1)
        error->one(FLERR,"this fix does not allow superposition with moving mesh fixes");

    // check if servo-wall is also a granular wall
    if (!fix_mesh->hasNeighList())
        error->one(FLERR,"The servo-wall requires a contact model. Therefore, it has to be used for a fix wall/gran too.");

    // no respa integrator
    if (strcmp(update->integrate_style,"respa") == 0)
        error->one(FLERR,"not respa-compatible");

}

/* ---------------------------------------------------------------------- */

void MeshModuleStressServo::init()
{
    
    mm_stress->set_p_ref(xcm_(0));

    // do some error checks
    error_checks();

    // get timestep
    reset_dt();

    // update ptrs
    v_ = mesh->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v");

    // check variables
    if (sp_str_) {
        sp_var_ = input->variable->find(sp_str_);
        if (sp_var_ < 0)
            error->one(FLERR,"Variable name does not exist");
        if (!input->variable->equalstyle(sp_var_))
            error->one(FLERR,"Variable is invalid style");
    }

    // controller parameters
    double r_min,r;
    const int nlocal = atom->nlocal;
    const double rPaMax = getMaxRad();

    r = r_min = BIG;
    for (int i = 0; i < nlocal; ++i) {
        r = atom->radius[i];
        r_min = MIN(r_min,r);
    }
    MPI_Min_Scalar(r_min,world);
    vel_min_ = ratio_*r_min/dtv_;
    
    // set pointers for controller
    switch (ctrl_style_) {
        case FORCE:
            pv_vec_ = mm_stress->f_total_ptr();
            ctrl_op_ = vcm_(0);
            ctrl_op_max_ = vel_max_;
            ctrl_op_min_ = vel_min_;
            break;
        case TORQUE:
            pv_vec_ = mm_stress->torque_total_ptr();
            ctrl_op_ = omegacm_(0);

            // find maximum distance axis-node
            if (rPaMax == 0)
                error->one(FLERR,"All mesh nodes are located at the rotation axis.");

            // maximum angular velocity
            ctrl_op_max_ = vel_max_/rPaMax;
            ctrl_op_min_ = vel_min_/rPaMax;
            break;
        default:
            error->one(FLERR,"This may not happen!");
            break;
    }

    // check maximal velocity
    const double skin = neighbor->skin;
    if(vel_max_ >= skin/(2.*dtv_))
        error->one(FLERR,"vel_max < skin/(2.*dt) required");

    // compute global number of contacts
    fix_mesh->meshNeighlist()->enableTotalNumContacts(true);
}

/* ---------------------------------------------------------------------- */

int MeshModuleStressServo::setmask()
{
    int mask = 0;
    mask |= INITIAL_INTEGRATE;
    mask |= FINAL_INTEGRATE;
    return mask;
}

/* ---------------------------------------------------------------------- */

void MeshModuleStressServo::setup_pre_force(int vflag)
{
    // set xcm_orig_
    xcm_orig_.set(0,xcm_(0));
    totalPhi_ = 0;
}

/* ---------------------------------------------------------------------- */

void MeshModuleStressServo::initial_integrate(int vflag)
{
    UNUSED(vflag);

    double dX[3],dx[3];

    // only if the wall should move
    if (int_flag_) {

        switch (ctrl_style_) {
            case FORCE:

                // update xcm by full step

                dx[0] = dtv_ * vcm_(0)[0];
                dx[1] = dtv_ * vcm_(0)[1];
                dx[2] = dtv_ * vcm_(0)[2];
                vectorAdd3D(xcm_(0),dx,xcm_(0));
                vectorSubtract3D(xcm_(0),xcm_orig_(0),dX);

                mesh->move(dX,dx);

                // update reference point to COM
                
                mm_stress->set_p_ref(xcm_(0));
                
                break;

            case TORQUE:
                const double incrementalPhi = dtv_ * vectorDot3D(omegacm_(0),axis_);
                totalPhi_ += incrementalPhi;

                //rotate the mesh
                mesh->rotate(totalPhi_,incrementalPhi,axis_,xcm_(0));

                break;
        }

    }

    // for area calculation
    // set vector of touching particles to zero
    //mod_andrew_.clear();
    mod_andrew_->clear_contacts();

}

/* ---------------------------------------------------------------------- */

void MeshModuleStressServo::final_integrate()
{
    // calcualte area
    // double area = mod_andrew_->area();
    mod_andrew_->area();

    double dfdt;

    // only if the wall should move
    if (int_flag_) {

        // variable force, wrap with clear/add
        if (sp_style_ == EQUAL) {

            modify->clearstep_compute();

            sp_mag_ = -input->variable->compute_equal(sp_var_);
            if (sp_mag_ == 0.) error->one(FLERR,"Set point (desired force/torque/shear) has to be != 0.0");
            sp_mag_inv_ = 1./fabs(sp_mag_);
            
            modify->addstep_compute(update->ntimestep + 1);

        }

        // auto mode - p-controller with variable proportional gain
        if (mode_flag_) {

            pv_mag_ = vectorDot3D(pv_vec_,axis_);
            err_ = (sp_mag_ - pv_mag_) * sp_mag_inv_;

            // Hard coded piecewise controller
            double ctrl_kp;
            const double err_low = 0.9;
            const double err_high = 1.0;
            const double ctrl_scale = 0.1;

            int totNumContacts = fix_mesh->meshNeighlist()->getTotalNumContacts();
            
            if (totNumContacts == 0) {
                // cruise mode
                ctrl_kp = ctrl_op_max_;
            } else {
                if (fabs(err_) <= err_low) {
                    ctrl_kp = ctrl_scale*ctrl_op_min_;
                } else if(fabs(err_) >= err_high) {
                    ctrl_kp = ctrl_op_min_;
                } else { // linear interpolation
                    ctrl_kp = ctrl_scale*ctrl_op_min_ + ((1-ctrl_scale)*ctrl_op_min_) * (fabs(err_)-err_low)/(err_high-err_low);
                }
            }

            vectorScalarMult3D(axis_,-ctrl_kp*err_,ctrl_op_);

        } else {

            // simple PID-controller

            // calc error and sum of the errors
            pv_mag_ = vectorDot3D(pv_vec_,axis_);
            err_ = (sp_mag_ - pv_mag_);
            sum_err_ += err_*dtv_;
            // derivative term
            // force used instead of the error in order to avoid signal spikes in case of change of the set point
            // de()/dt = - dforce()/dt for constant set point
            dfdt = -( pv_mag_ - old_pv_mag_)/dtv_;

            // vel points opposite to force vector
            const double ctrl_op_mag = -ctrl_op_max_ * (err_ * kp_ + sum_err_ * ki_ + dfdt * kd_) * sp_mag_inv_;
            vectorScalarMult3D(axis_,ctrl_op_mag,ctrl_op_);
            // save process value for next timestep
            old_pv_mag_ = pv_mag_;

        }

        limit_vel();

        switch (ctrl_style_) {
            case FORCE:
                set_v_node();
                break;
            case TORQUE:
                set_v_node_rotate();
                break;
        }

    }

}

/* ---------------------------------------------------------------------- */

void MeshModuleStressServo::limit_vel()
{

    double maxOutput;
    const double vmag = vectorMag3D(ctrl_op_);

    // saturation of the velocity
    int totNumContacts = fix_mesh->meshNeighlist()->getTotalNumContacts();
    if (mode_flag_ && totNumContacts > 0) {
        maxOutput = ctrl_op_min_;
    } else {
        maxOutput = ctrl_op_max_;
    }

    // saturation of the velocity

    if(vmag > maxOutput && vmag != 0) {
        const double factor = maxOutput / vmag;

        // scale output vector
        vectorScalarMult3D(ctrl_op_,factor);

        // anti-windup of the integral part (only if ki_>0)
        if (ki_ > 0) {
            const double ctrl_out_mag = vectorDot3D(ctrl_op_,axis_);
            sum_err_ = (-sgn(ctrl_out_mag) * sp_mag_ -err_*kp_)/ki_; //inverted controller equation
        }

    }
}

/* ---------------------------------------------------------------------- */

void MeshModuleStressServo::set_v_node()
{
    const int nall = mesh->size();
    const int nnodes = mesh->numNodes();

    for(int i = 0; i < nall; ++i)
        for(int j = 0; j < nnodes; ++j)
            v_->set(i,j,vcm_(0));

}

/* ---------------------------------------------------------------------- */

void MeshModuleStressServo::set_v_node_rotate()
{
    double node[3],rPA[3],vRot[3];

    const int nall = mesh->size();
    const int nnodes = mesh->numNodes();

    for(int i = 0; i < nall; ++i)
    {
        for(int j = 0; j < nnodes; ++j)
        {
            mesh->node(i, j, node);
            vectorSubtract3D(node,xcm_(0),rPA);
            vectorCross3D(omegacm_(0),rPA,vRot);
            v_->set(i,j,vRot);
        }
    }
}

/* ---------------------------------------------------------------------- */

void MeshModuleStressServo::reset_dt()
{
    dtv_ = update->dt;
    dtf_ = 0.5 * update->dt * force->ftm2v;
}

/* ----------------------------------------------------------------------
    called during wall force calc

    detected contacts are registered to contribute to the area of the servo
------------------------------------------------------------------------- */

void MeshModuleStressServo::add_particle_contribution(int ip, double *frc,
                                                                                                                    double *delta, int iTri, double *v_wall)
{
    //double *x = atom->x[ip];
    //double r = atom->radius[ip];

    //Circle c = {x[(1+dim_)%3], x[(2+dim_)%3], r};
    //mod_andrew_->add_contact(c);
}

/* ----------------------------------------------------------------------
    get maximal distance rotation axis - mesh nodes
------------------------------------------------------------------------- */

double MeshModuleStressServo::getMaxRad()
{
    double node[3],rPA[3],vRot[3];

    double rPaMax = 0;
    const int nall = mesh->size();
    const int nnodes = mesh->numNodes();

    for(int i = 0; i < nall; ++i)
    {
        for(int j = 0; j < nnodes; ++j)
        {
            mesh->node(i, j, node);
            vectorSubtract3D(node,xcm_(0),rPA);
            vectorCross3D(axis_,rPA,vRot);
            rPaMax = MAX(rPaMax,vectorMag3D(vRot));
        }
    }

    // get max for all processors
    MPI_Max_Scalar(rPaMax,world);
    
    return rPaMax;
}

/* ----------------------------------------------------------------------
     return total force or torque component on body or xcm
------------------------------------------------------------------------- */

double MeshModuleStressServo::compute_vector(int n)
{
        return xcm_(0)[n];
}

/* ---------------------------------------------------------------------- */

int MeshModuleStressServo::modify_param(int narg, char **arg)
{
    std::string arg0(arg[0]);
    std::size_t slash = arg0.find_first_of('/');
    std::string command(arg0.substr(slash == std::string::npos ? 0 : slash+1));

    if (command.compare("integrate") == 0) {
        if (narg < 2) error->one(FLERR,"not enough arguments for fix_modify 'integrate'");

        if (strcmp(arg[1],"start") == 0) {
            int_flag_ = true;
        } else if (strcmp(arg[1],"stop") == 0) {
            int_flag_ = false;
        } else
            error->one(FLERR,"wrong argument for fix_modify 'integrate'");

        return 2;
    } else if (command.compare("vel_max") == 0) {
        if (narg < 2) error->one(FLERR,"not enough arguments for 'vel_max'");
        vel_max_ = force->numeric(FLERR,arg[1]);
        if(vel_max_ <= 0.)
            error->one(FLERR,"vel_max > 0 required");

        return 2;
    } else if (command.compare("target_val") == 0) {
        if (narg < 2) error->one(FLERR,"not enough arguments for fix_modify 'target_val'");
        if (strstr(arg[1],"v_") == arg[1]) {
            const int n = strlen(&arg[1][2]) + 1;
            sp_str_ = new char[n];
            strcpy(sp_str_,&arg[1][2]);

            // check variables
            if (sp_str_) {
                sp_var_ = input->variable->find(sp_str_);
                if (sp_var_ < 0)
                    error->one(FLERR,"Variable name does not exist");
                if (input->variable->equalstyle(sp_var_)) sp_style_ = EQUAL;
                else error->one(FLERR,"Variable is invalid style");
            }

        } else {
            sp_mag_ = -force->numeric(FLERR,arg[1]); // the resultant force/torque/shear acts in opposite direction --> negative value
            if (sp_mag_ == 0.) error->one(FLERR,"'target_val' (desired force/torque) has to be != 0.0");
            sp_mag_inv_ = 1./fabs(sp_mag_);
            sp_style_ = CONSTANT;
        }

        resetIntegrator();

        return 2;
    } else if (command.compare("ctrlParam") == 0) {
        if (narg < 4) error->one(FLERR,"not enough arguments for fix_modify 'ctrlParam'");

        kp_ = force->numeric(FLERR,arg[1]);
        ki_ = force->numeric(FLERR,arg[2]);
        kd_ = force->numeric(FLERR,arg[3]);

        resetIntegrator();

        return 4;
    } else if (command.compare("ratio") == 0) {
        if (narg < 2)
            error->one(FLERR,"not enough arguments for fix_modify 'ratio'");

        if (!mode_flag_)
            error->warning(FLERR,"Modifying 'ratio' makes sense only when mode='auto'");

        ratio_ = force->numeric(FLERR,arg[1]);

        return 2;
    }

    return 0;
}
