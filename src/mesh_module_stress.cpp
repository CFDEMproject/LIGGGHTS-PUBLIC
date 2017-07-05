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

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author for surface velocity rotation only
   Evan Smuts (U Cape Town)
------------------------------------------------------------------------- */

#include "mesh_module_stress.h"
#include <stdio.h>
#include <string.h>
#include "error.h"
#include "force.h"
#include "modify.h"
#include "comm.h"
#include "math_extra.h"
#include "fix_property_global.h"
#include "fix_gravity.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define EPSILON 0.0001

/* ---------------------------------------------------------------------- */

MeshModuleStress::MeshModuleStress(LAMMPS *lmp, int &iarg_, int narg, char **arg, FixMeshSurface *fix_mesh)
:   MeshModule(lmp, iarg_, narg, arg, fix_mesh),
    
    updatedStresses_(false),
    p_ref_(*mesh->prop().addGlobalProperty<VectorContainer<double,3> >("p_ref","comm_none","frame_general","restart_yes")),
    f_(0),
    sigma_n_(0),
    sigma_t_(0),
    wear_flag_(0),
    k_finnie_(0),
    wear_(0),
    wear_step_(0),
    wear_increment_(NULL),
    store_wear_increment_(false)
{
    vectorZeroize3D(f_total_);
    vectorZeroize3D(torque_total_);
    vectorZeroize3D(f_total_old_);
    vectorZeroize3D(torque_total_old_);

    double zerovec[3] = {0.,0.,0.};
    mesh->prop().setGlobalProperty<VectorContainer<double,3> >("p_ref",zerovec);

    // override default from base
    stress_flag_ = true;

    fix_mesh->set_global_freq(1);
    fix_mesh->set_extvector(1);

    // parse further args

    bool hasargs = true;
    while(iarg_ < narg && hasargs)
    {
        hasargs = false;
        if (strcmp(arg[iarg_],"reference_point") == 0) {
            if (narg < iarg_+4)
                error->one(FLERR,"not enough arguments");
            if(fix_mesh->manipulated())
                error->warning(FLERR,"Mesh for fix mesh/surface/stress has been scaled, moved, or rotated.\n"
                               "Please note that values for 'reference_point' refer to the scaled, moved, or rotated configuration");
            iarg_++;
            double _p_ref[3];
            _p_ref[0] = force->numeric(FLERR,arg[iarg_++]);
            _p_ref[1] = force->numeric(FLERR,arg[iarg_++]);
            _p_ref[2] = force->numeric(FLERR,arg[iarg_++]);
            mesh->prop().setGlobalProperty<VectorContainer<double,3> >("p_ref",_p_ref);
            hasargs = true;
        } else if(strcmp(arg[iarg_],"stress") == 0) {
            if (narg < iarg_+2)
                error->one(FLERR,"not enough arguments");
            iarg_++;
            if(strcmp(arg[iarg_],"on") == 0)
                stress_flag_ = true;
            else if(strcmp(arg[iarg_],"off") == 0)
                stress_flag_ = false;
            else
                error->one(FLERR,"expecting 'on' or 'off' as stress argument");
            iarg_++;
            hasargs = true;
        } else if(strcmp(arg[iarg_],"wear") == 0) {
            if (narg < iarg_+2)
                error->one(FLERR,"not enough arguments");
            iarg_++;
            if(strcmp(arg[iarg_],"finnie") == 0)
                wear_flag_ = 1;
            else if(strcmp(arg[iarg_],"off") == 0)
                wear_flag_ = 0;
            else
                error->one(FLERR,"expecting 'finnie' or 'off' as wear argument");
            iarg_++;
            hasargs = true;
        }
    }
}

/* ---------------------------------------------------------------------- */

MeshModuleStress::~MeshModuleStress()
{

}

/* ---------------------------------------------------------------------- */

void MeshModuleStress::post_create_pre_restart()
{
    
    //Np -->register properties and set values for non-restart properties here

    if(stress_flag_)
        regStress();

    if(wear_flag_)
        regWear();
}

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void MeshModuleStress::post_create()
{
    
    if(stress_flag_)
        zeroizeStress();

    if(wear_flag_)
        zeroizeWear();
}

/* ---------------------------------------------------------------------- */

void MeshModuleStress::regStress()
{
    mesh->prop().addElementProperty<VectorContainer<double,3> >("f","comm_reverse","frame_invariant","restart_no");
    mesh->prop().addElementProperty<ScalarContainer<double> >("sigma_n","comm_none","frame_invariant","restart_no");
    mesh->prop().addElementProperty<ScalarContainer<double> >("sigma_t","comm_none","frame_invariant","restart_no");
}
/* ---------------------------------------------------------------------- */

void MeshModuleStress::zeroizeStress()
{
    mesh->prop().getElementProperty<VectorContainer<double,3> >("f")->setAll(0.);
    mesh->prop().getElementProperty<ScalarContainer<double> >("sigma_n")->setAll(0.);
    mesh->prop().getElementProperty<ScalarContainer<double> >("sigma_t")->setAll(0.);
}

/* ---------------------------------------------------------------------- */

void MeshModuleStress::regWear()
{
    mesh->prop().addElementProperty<ScalarContainer<double> >("wear","comm_exchange_borders","frame_invariant","restart_yes");
    mesh->prop().getElementProperty<ScalarContainer<double> >("wear")->setAll(0.);
    mesh->prop().addElementProperty<ScalarContainer<double> >("wear_step","comm_reverse","frame_invariant","restart_no");
    if (store_wear_increment_)
        mesh->prop().addElementProperty<ScalarContainer<double> >("wear_increment","comm_reverse","frame_invariant","restart_no");

}

/* ---------------------------------------------------------------------- */

void MeshModuleStress::zeroizeWear()
{
    mesh->prop().getElementProperty<ScalarContainer<double> >("wear_step")->setAll(0.);
}

/* ---------------------------------------------------------------------- */

void MeshModuleStress::init()
{
    if(stress_flag_)
    {
        f_ = mesh->prop().getElementProperty<VectorContainer<double,3> >("f");
        sigma_n_ = mesh->prop().getElementProperty<ScalarContainer<double> >("sigma_n");
        sigma_t_ = mesh->prop().getElementProperty<ScalarContainer<double> >("sigma_t");
        if(!f_ || !sigma_n_ || !sigma_t_)
            error->one(FLERR,"Internal error");
    }

    if(wear_flag_)
    {
        k_finnie_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property("k_finnie","property/global","peratomtypepair",atom->ntypes,atom->ntypes,"meshmodule/stress"))->get_array();
        wear_ = mesh->prop().getElementProperty<ScalarContainer<double> >("wear");
        wear_step_ = mesh->prop().getElementProperty<ScalarContainer<double> >("wear_step");
        if (store_wear_increment_)
            wear_increment_ = mesh->prop().getElementProperty<ScalarContainer<double> >("wear_increment");
        if(!wear_ || ! wear_step_ || (store_wear_increment_ && !wear_increment_))
            error->one(FLERR,"Internal error");
    }
}

/* ---------------------------------------------------------------------- */

int MeshModuleStress::setmask()
{
    int mask = 0;
    mask |= PRE_FORCE;
    mask |= FINAL_INTEGRATE;
    return mask;
}

/* ---------------------------------------------------------------------- */

void MeshModuleStress::pre_force(int vflag)
{
    if(trackStress())
    {
        vectorCopy3D(f_total_,f_total_old_);
        vectorCopy3D(torque_total_,torque_total_old_);
        vectorZeroize3D(f_total_);
        vectorZeroize3D(torque_total_);
        updatedStresses_ = false;
    }

    if (wear_flag_ && store_wear_increment_)
        mesh->prop().getElementProperty<ScalarContainer<double> >("wear_increment")->setAll(0.);

}

/* ---------------------------------------------------------------------- */

void MeshModuleStress::final_integrate()
{
    calc_total_force();
}

/* ----------------------------------------------------------------------
   called during wall force calc
------------------------------------------------------------------------- */

void MeshModuleStress::add_particle_contribution(int ip,double *frc,
                                double *delta,int iTri,double *v_wall)
{
    double E,c[3],v_rel[3],v_rel_mag,cos_gamma,sin_gamma,sin_2gamma;
    double contactPoint[3]={},surfNorm[3], tmp[3], tmp2[3];

    // do not include if not in fix group
    if(!(atom->mask[ip] & fix_mesh->get_groupbit())) return;

    double *x = atom->x[ip];
    double *v = atom->v[ip];

    vectorNegate3D(frc);

    vectorAdd3D(x,delta,contactPoint);

    if(trackStress())
    {
        
        // add contribution to triangle force
        vectorAdd3D(f(iTri),frc,f(iTri));

        // add contribution to total body force and torque
        vectorAdd3D(f_total_,frc,f_total_);
        vectorSubtract3D(contactPoint,p_ref_(0),tmp);
        
        vectorCross3D(tmp,frc,tmp2); // tmp2 is torque contrib
        vectorAdd3D(torque_total_,tmp2,torque_total_);
    }

    // add wear if applicable
    if(trackWear())
    {
        if (store_wear_increment_)
            wear_increment(iTri) = 0.0;

        vectorSubtract3D(contactPoint,x,c);

        // calculate relative velocity
        vectorSubtract3D(v,v_wall,v_rel);

        if(vectorDot3D(c,v_rel) < 0.) return;
        v_rel_mag = vectorMag3D(v_rel);

        // get surface normal
        
        mesh->surfaceNorm(iTri,surfNorm);

        // return if no relative velocity
        if(0.0000001 > v_rel_mag)
            return;

        sin_gamma = MathExtraLiggghts::abs(vectorDot3D(v_rel,surfNorm)) / (v_rel_mag);
        cos_gamma = MathExtraLiggghts::abs(vectorCrossMag3D(v_rel,surfNorm)) / (v_rel_mag);

        if(cos_gamma > 1.) cos_gamma = 1.;
        if(sin_gamma > 1.) sin_gamma = 1.;

        if(cos_gamma < EPSILON || 3.*sin_gamma > cos_gamma)
        {
            E = 0.33333 * cos_gamma * cos_gamma;
            
        }
        else
        {
            sin_2gamma = 2. * sin_gamma * cos_gamma;
            E = sin_2gamma - 3. * sin_gamma * sin_gamma;
            
        }
        const int atom_type_mesh = fix_mesh->atomTypeWall();
        E *= 2.*k_finnie_[atom_type_mesh-1][atom->type[ip]-1] * v_rel_mag * vectorMag3D(frc);
        
        const double part_wear_increment = E*update->dt / mesh->areaElem(iTri);
        if (store_wear_increment_)
            wear_increment(iTri) = part_wear_increment;
        wear_step(iTri) += part_wear_increment;
    }
}

/* ----------------------------------------------------------------------
   add external force (such as gravity)
   called by all procs, only proc0 adds
   has to be called before final_integrate()
------------------------------------------------------------------------- */

void MeshModuleStress::add_global_external_contribution(double *frc)
{
    vectorAdd3D(f_total_,frc,f_total_);
}

void MeshModuleStress::add_global_external_contribution(double *frc, double *trq)
{
    vectorAdd3D(f_total_,frc,f_total_);
    vectorAdd3D(torque_total_,trq,torque_total_);
}

/* ----------------------------------------------------------------------
   allreduce total force on tri
------------------------------------------------------------------------- */

void MeshModuleStress::calc_total_force()
{
    double surfNorm[3], invSurfArea, temp[3];
    int nTri = mesh->size();

    // add wear from this step to total wear
    if(trackWear())
    {
        for(int i = 0; i < nTri; i++)
        {
            wear(i) += wear_step(i);
            wear_step(i) = 0.;
        }
    }

    // calculate normal and shear stress
    if(trackStress())
    {
        // total force and torque on mesh

        MPI_Sum_Vector(f_total_,3,world);
        MPI_Sum_Vector(torque_total_,3,world);
        
        updatedStresses_ = true; // f_total_ and torque_total_ are now up-to-date

        for(int i = 0; i < nTri; i++)
        {
            // get element surface norm and area
            mesh->surfaceNorm(i,surfNorm);
            invSurfArea = 1./mesh->areaElem(i);

            // calculate normal force
            sigma_n(i) = vectorDot3D(f(i),surfNorm);
            
            // calculate tangential force
            vectorScalarMult3D(surfNorm,sigma_n(i),temp);
            vectorSubtract3D(f(i),temp,temp);
            sigma_t(i) = vectorMag3D(temp);

            // make both positive
            // necessary since orientation of surfNorm not known
            sigma_n(i) = MathExtraLiggghts::abs(sigma_n(i));
            sigma_t(i) = MathExtraLiggghts::abs(sigma_t(i));

            // divide by area so have stress
            sigma_n(i) *= invSurfArea;
            sigma_t(i) *= invSurfArea;
            
        }
    }
}

/* ----------------------------------------------------------------------
   return total force or torque component on body
------------------------------------------------------------------------- */

double MeshModuleStress::compute_vector(int n)
{
    if (n < 3)
        return updatedStresses_ ? f_total_[n] : f_total_old_[n];
    else if (n < 6)
        return updatedStresses_ ? torque_total_[n-3] : torque_total_old_[n-3];
    else if (n < 9)
        return p_ref_(0)[n-6];
    return 0.0;
}
