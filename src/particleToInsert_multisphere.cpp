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

#include "particleToInsert_multisphere.h"
#include <cmath>
#include "error.h"
#include "vector_liggghts.h"
#include "atom.h"
#include "atom_vec.h"
#include "modify.h"
#include "domain.h"
#include "comm.h"
#include "fix.h"
#include "fix_multisphere.h"
#include <string.h>
#include "math_extra_liggghts.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ParticleToInsertMultisphere::ParticleToInsertMultisphere(LAMMPS* lmp,int ns) : ParticleToInsert(lmp,ns)
{
    memory->create(displace,nparticles,3,"displace");
    memory->create(volumeweight,nparticles,"displace");

    for(int i = 0; i < nparticles; i++)
       vectorZeroize3D(displace[i]);

    fflag[0] = fflag[1] = fflag[2] = true;
    tflag[0] = tflag[1] = tflag[2] = true;
}

/* ---------------------------------------------------------------------- */

ParticleToInsertMultisphere::~ParticleToInsertMultisphere()
{
    memory->destroy(displace);
    memory->destroy(volumeweight);
}

/* ---------------------------------------------------------------------- */

int ParticleToInsertMultisphere::set_x_v_omega(double *x, double *v, double *omega, double *quat)
{
    double disp_glob[3];

    vectorCopy3D(x,xcm_ins);
    vectorCopy4D(quat,quat_ins);
    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);

    //if(!isIdentityQuat4D(quat_ins)) error->warning(FLERR,"quaternion rotation untested in ParticleToInsertMultisphere");

    MathExtraLiggghts::vec_quat_rotate(ex_space,quat);
    MathExtraLiggghts::vec_quat_rotate(ey_space,quat);
    MathExtraLiggghts::vec_quat_rotate(ez_space,quat);

    for(int j = 0; j < nparticles; j++)
    {
        MathExtraLiggghts::local_coosys_to_cartesian(disp_glob,displace[j],ex_space,ey_space,ez_space);
        vectorAdd3D(x,disp_glob,x_ins[j]);
    }

    return nparticles;
}

/* ---------------------------------------------------------------------- */

int ParticleToInsertMultisphere::check_near_set_x_v_omega(double *x,double *v, double *omega, double *quat, RegionNeighborList<interpolate_no> & neighList)//, double **xnear, int &nnear)
{

    // check every sphere against all others in xnear
    // if no overlap add to xnear
    double disp_glob[3]; //del[3], rsq, radsum;
    double ex_space_try[3], ey_space_try[3], ez_space_try[3];

    // rotate if needed
    // calculate x_ins for this quaternion
    // do this in a "try" step since we do not know if we will succeed

    //if(!isIdentityQuat4D(quat_ins)) error->one(FLERR,"quaternion rotation untested in ParticleToInsertMultisphere");
    MathExtraLiggghts::vec_quat_rotate(ex_space,quat,ex_space_try);
    MathExtraLiggghts::vec_quat_rotate(ey_space,quat,ey_space_try);
    MathExtraLiggghts::vec_quat_rotate(ez_space,quat,ez_space_try);
    for(int j = 0; j < nparticles; j++)
    {
        MathExtraLiggghts::local_coosys_to_cartesian(disp_glob,displace[j],ex_space_try,ey_space_try,ez_space_try);
        vectorAdd3D(x,disp_glob,x_ins[j]);
    }

    // check for overlap with nearby particles
    for(int j = 0; j < nparticles; j++)
    {
        if(neighList.hasOverlap(x_ins[j], radius_ins[j])) {
            return 0;
        }
    }

    /*
    for(int i = 0; i < nnear; i++)
    {
        for(int j = 0; j < nparticles; j++)
        {
           vectorSubtract3D(x_ins[j],xnear[i],del);
           rsq = vectorMag3DSquared(del);
           radsum = radius_ins[j] + xnear[i][3];

           // no success in overlap
           if (rsq <= radsum*radsum) return 0;
        }
    }*/

    // no overlap with any other - success

    vectorCopy3D(x,xcm_ins);
    vectorCopy4D(quat,quat_ins);
    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);

    // set axis to the value we succeeded for, x_ins already up to date
    vectorCopy3D(ex_space_try,ex_space);
    vectorCopy3D(ey_space_try,ey_space);
    vectorCopy3D(ez_space_try,ez_space);

    // add to xnear
    for(int j = 0; j < nparticles; j++)
    {
        neighList.insert(x_ins[j], radius_ins[j]);
        /*
        vectorCopy3D(x_ins[j],xnear[nnear]);
        xnear[nnear][3] = radius_ins[j];
        nnear++;
        */
    }

    return nparticles;
}

/* ---------------------------------------------------------------------- */

int ParticleToInsertMultisphere::insert()
{
    int inserted = 0;
    int nfix = modify->nfix;
    Fix **fix = modify->fix;

    // perform the actual particle insertion
    // add particles, set coordinate and radius
    // set group mask to "all" plus fix groups
    
    for(int i = 0; i < nparticles; i++)
    {
        
        //if (domain->is_in_extended_subdomain(xcm_ins))
        //{
            inserted++;
            atom->avec->create_atom(atom_type,x_ins[i]);
            
            int m = atom->nlocal - 1;
            atom->mask[m] = 1 | groupbit;
            atom->radius[m] = radius_ins[i];
            atom->density[m] = density_ins;

            atom->rmass[m] = mass_ins; 

            vectorZeroize3D(atom->v[m]);
            vectorZeroize3D(atom->omega[m]);
            vectorZeroize3D(atom->f[m]);
            vectorZeroize3D(atom->torque[m]);

            for (int j = 0; j < nfix; j++)
            {
               if (fix[j]->create_attribute) fix[j]->set_arrays(m);
            }
        //}
    }

    // now rigid body insertion
    
    int nlocal = atom->nlocal;

    if(modify->n_fixes_style("multisphere") != 1) {
        printf("Number of fix multisphere used: %d\n", modify->n_fixes_style("multisphere"));
        error->one(FLERR,"Multi-sphere particle inserted: You have to use exactly one fix multisphere.");
    }

    FixMultisphere *fix_multisphere = static_cast<FixMultisphere*>(modify->find_fix_style("multisphere",0));

    fix_multisphere->data().add_body(nparticles,xcm_ins,xcm_to_xbound,r_bound_ins, v_ins, omega_ins, mass_ins,
                                density_ins,atom_type,type_ms,inertia,ex_space,ey_space,ez_space,displace,fflag,tflag);

    // set displace correctly, set body to -2
    
    int i = 0;
    for(int isphere = nlocal-nparticles; isphere < nlocal; isphere++)
    {
        
        fix_multisphere->set_body_displace(isphere,displace[i],-2,volumeweight[i]);
        i++;
    }

    return inserted;
}

/* ---------------------------------------------------------------------- */

void ParticleToInsertMultisphere::scale_pti(double r_scale)
{
    error->one(FLERR,"scale_pti not implemented in ParticleToInsertMultisphere");
}

/* ---------------------------------------------------------------------- */

void ParticleToInsertMultisphere::random_rotate(double rn1,double rn2, double rn3)
{
    
    if(nparticles==1)return;

    double *vert_before_rot;
    double vert_after_rot[3]={0.0,0.0,0.0};

    double phix=rn1*2.*M_PI;
    double phiy=rn2*2.*M_PI;
    double phiz=rn3*2.*M_PI;

    double cos_phix = cos(phix);
    double cos_phiy = cos(phiy);
    double cos_phiz = cos(phiz);
    double sin_phix = sin(phix);
    double sin_phiy = sin(phiy);
    double sin_phiz = sin(phiz);

    for(int i=0;i<3;i++)
    {
        if     (i==0) vert_before_rot=ex_space;
        else if(i==1) vert_before_rot=ey_space;
        else if(i==2) vert_before_rot=ez_space;

        vert_after_rot[0] = vert_before_rot[0]*cos_phiy*cos_phiz+vert_before_rot[1]*(cos_phiz*sin_phix*sin_phiy-cos_phix*sin_phiz)+vert_before_rot[2]*(cos_phix*cos_phiz*sin_phiy+sin_phix*sin_phiz);
        vert_after_rot[1] = vert_before_rot[0]*cos_phiy*sin_phiz+vert_before_rot[2]*(-cos_phiz*sin_phix+cos_phix*sin_phiy*sin_phiz)+vert_before_rot[1]*(cos_phix*cos_phiz+sin_phix*sin_phiy*sin_phiz);
        vert_after_rot[2] = vert_before_rot[2]*cos_phix*cos_phiy+vert_before_rot[1]*cos_phiy*sin_phix-vert_before_rot[0]*sin_phiy;

        if     (i==0) for(int j=0;j<3;j++) ex_space[j]=vert_after_rot[j];
        else if(i==1) for(int j=0;j<3;j++) ey_space[j]=vert_after_rot[j];
        else if(i==2) for(int j=0;j<3;j++) ez_space[j]=vert_after_rot[j];
    }

    for(int i=0;i<nparticles;i++)
    {
        x_ins[i][0] = xcm_ins[0] + ex_space[0]*displace[i][0] +   ey_space[0]*displace[i][1] +   ez_space[0]*displace[i][2];
        x_ins[i][1] = xcm_ins[1] + ex_space[1]*displace[i][0] +   ey_space[1]*displace[i][1] +   ez_space[1]*displace[i][2];
        x_ins[i][2] = xcm_ins[2] + ex_space[2]*displace[i][0] +   ey_space[2]*displace[i][1] +   ez_space[2]*displace[i][2];
    }

}
