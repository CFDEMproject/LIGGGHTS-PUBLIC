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

#include "particleToInsert.h"
#include <cmath>
#include "error.h"
#include "update.h"
#include "domain.h"
#include "atom.h"
#include "atom_vec.h"
#include "fix_property_atom.h"
#include "vector_liggghts.h"
#include "math_extra_liggghts.h"
#include "modify.h"

using namespace LAMMPS_NS;

ParticleToInsert::ParticleToInsert(LAMMPS* lmp, int ns, FixPropertyAtom * const _fix_release) :
    Pointers(lmp),
    id_ins(-1),
    fix_release(_fix_release),
    fix_property(NULL),
    n_fix_property(0),
    fix_property_nentry(NULL),
    fix_property_value(NULL),
    fix_template_(NULL)
{
    groupbit = 0;

    distorder = -1;

    nparticles = ns;

    memory->create(x_ins,nparticles,3,"x_ins");
    radius_ins = new double[nparticles];

    atom_type_vector = new int[nparticles];
    atom_type_vector_flag = false;

    if (fix_release && fix_release->get_nvalues() <= 14)
        error->all(FLERR, "Invalid fix_release, more than 14 entries required");
}

/* ---------------------------------------------------------------------- */

ParticleToInsert::~ParticleToInsert()
{
        memory->destroy(x_ins);
        delete []radius_ins;
        delete []atom_type_vector;
        if (fix_property_value)
        {
            for (int i = 0; i < n_fix_property; i++)
                delete [] fix_property_value[i];
            delete [] fix_property_value;
        }
        if (fix_property)
            delete [] fix_property;
}

/* ---------------------------------------------------------------------- */

int ParticleToInsert::insert()
{
    // perform the actual insertion
    // add particles, set coordinate and radius
    // set group mask to "all" plus fix groups

    int inserted = 0;
    int nfix = modify->nfix;
    Fix **fix = modify->fix;

    for(int i = 0; i < nparticles; i++)
    {
        
        //if (domain->is_in_extended_subdomain(x_ins[i]))
        //{
                
                inserted++;
                if(atom_type_vector_flag)
                    atom->avec->create_atom(atom_type_vector[i],x_ins[i]);
                else
                    atom->avec->create_atom(atom_type,x_ins[i]);
                int m = atom->nlocal - 1;
                atom->mask[m] = 1 | groupbit;
                vectorCopy3D(v_ins,atom->v[m]);
                vectorCopy3D(omega_ins,atom->omega[m]);
                atom->radius[m] = radius_ins[i];
                atom->density[m] = density_ins;
                
                atom->rmass[m] = (1==nparticles)? (mass_ins) : (4.18879020479/*4//3*pi*/*radius_ins[i]*radius_ins[i]*radius_ins[i]*density_ins);

                //pre_set_arrays() called via FixParticleDistribution
                for (int j = 0; j < nfix; j++)
                   if (fix[j]->create_attribute) fix[j]->set_arrays(m);

                // apply fix property setting coming from fix insert
                // this overrides the set_arrays call above
                if(fix_property)
                {
                    for (int j = 0; j < n_fix_property; j++)
                    {
                        if (fix_property_nentry[j] == 1)
                            fix_property[j]->vector_atom[m] = fix_property_value[j][0];
                        else
                        {
                            for (int k = 0; k < fix_property_nentry[j]; k++)
                                fix_property[j]->array_atom[m][k] = fix_property_value[j][k];
                        }
                    }
                }
                if (fix_template_)
                    fix_template_->vector_atom[m] = (double)distorder;
                if (fix_release)
                    fix_release->array_atom[m][14] = (double) id_ins;
        //}
    }
    
    return inserted;
}

/* ---------------------------------------------------------------------- */

/*
int ParticleToInsert::check_near_set_x_v_omega(double *x,double *v, double *omega, double *quat, double **xnear, int &nnear)
{
    if(nparticles > 1)
        return check_near_set_x_v_omega_ms(x,v, omega,quat,xnear,nnear);

    // check sphere against all others in xnear
    // if no overlap add to xnear
    double del[3], rsq, radsum;

    vectorCopy3D(x,x_ins[0]);

    for(int i = 0; i < nnear; i++)
    {
        vectorSubtract3D(x_ins[0],xnear[i],del);
        rsq = vectorMag3DSquared(del);
        
/*
        radsum = radius_ins[0] + xnear[i][3];

        // no success in overlap
        if (rsq <= radsum*radsum) return 0;
    }

    // no overlap with any other - success

    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);

    // add to xnear
    vectorCopy3D(x_ins[0],xnear[nnear]);
    xnear[nnear][3] = radius_ins[0];
    nnear++;

    return 1;
}*/

/* ---------------------------------------------------------------------- */

/*
int ParticleToInsert::check_near_set_x_v_omega_ms(double *x,double *v, double *omega, double *quat, double **xnear, int &nnear)
{
    // x is position where insertion should take place
    // v and omega are the velocity and omega for the newly inserted particles
    double rel[3],xins_j_try[3];
    double del[3], rsq, radsum;

    // check insertion position, take quat into account
    // relative position of spheres to each other already stored at this point
    // check sphere against all others in xnear
    for(int j = 0; j < nspheres; j++)
    {
        // take orientation into account; x_bound_ins is in the global coordinate system
        // calculate xins_j_try for every sphere and check if would work
        vectorSubtract3D(x_ins[j],x_bound_ins,rel);
        MathExtraLiggghts::vec_quat_rotate(rel,quat);
        vectorAdd3D(rel,x,xins_j_try);

        for(int i = 0; i < nnear; i++)
        {
           vectorSubtract3D(xins_j_try,xnear[i],del);
           rsq = vectorMag3DSquared(del);
           radsum = radius_ins[j] + xnear[i][3];

           // no success in overlap
           if (rsq <= radsum*radsum)
            return 0;
        }
    }

    // no overlap with any other - success
    // set x_ins, v_ins and omega_ins
    for(int j = 0; j < nspheres; j++)
    {
        vectorSubtract3D(x_ins[j],x_bound_ins,rel);
        MathExtraLiggghts::vec_quat_rotate(rel,quat);
        vectorAdd3D(rel,x,x_ins[j]);
    }
    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);

    // add to xnear for future checks
    for(int j = 0; j < nspheres; j++)
    {
        vectorCopy3D(x_ins[j],xnear[nnear]);
        xnear[nnear][3] = radius_ins[j];
        nnear++;
    }

    return nparticles;
}*/

/* ---------------------------------------------------------------------- */

int ParticleToInsert::check_near_set_x_v_omega(double *x,double *v, double *omega, double *quat, RegionNeighborList<interpolate_no> & neighList)
{
    if(nparticles > 1)
        return check_near_set_x_v_omega_ms(x,v, omega,quat,neighList);

    vectorCopy3D(x,x_ins[0]);

    if(neighList.hasOverlap(x_ins[0], radius_ins[0])) {
        return 0;
    }

    // no overlap with any other - success

    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);

    neighList.insert(x_ins[0], radius_ins[0]);

    return 1;
}

/* ---------------------------------------------------------------------- */

int ParticleToInsert::check_near_set_x_v_omega_ms(double *x,double *v, double *omega, double *quat, RegionNeighborList<interpolate_no> & neighList)
{
    // x is position where insertion should take place
    // v and omega are the velocity and omega for the newly inserted particles
    double rel[3],xins_j_try[3];
    //double del[3], rsq, radsum;

    // check insertion position, take quat into account
    // relative position of spheres to each other already stored at this point
    // check sphere against all others in xnear
    for(int j = 0; j < nparticles; j++)
    {
        // take orientation into account; x_bound_ins is in the global coordinate system
        // calculate xins_j_try for every sphere and check if would work
        vectorSubtract3D(x_ins[j],x_bound_ins,rel);
        MathExtraLiggghts::vec_quat_rotate(rel,quat);
        vectorAdd3D(rel,x,xins_j_try);

        if(neighList.hasOverlap(xins_j_try, radius_ins[j])) {
            return 0;
        }
    }

    // no overlap with any other - success
    // set x_ins, v_ins and omega_ins
    for(int j = 0; j < nparticles; j++)
    {
        vectorSubtract3D(x_ins[j],x_bound_ins,rel);
        MathExtraLiggghts::vec_quat_rotate(rel,quat);
        vectorAdd3D(rel,x,x_ins[j]);
    }
    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);

    // add to xnear for future checks
    for(int j = 0; j < nparticles; j++)
    {
        neighList.insert(x_ins[j], radius_ins[j]);
    }

    return nparticles;
}

/* ---------------------------------------------------------------------- */

int ParticleToInsert::set_x_v_omega(double *x, double *v, double *omega, double *quat)
{
    double rel[3];

    // x is position where insertion should take place
    // v and omega are the velocity and omega for the newly inserted particles

    // add insertion position
    // relative position of spheres to each other already stored at this point
    // also take quat into account
    for(int j = 0; j < nparticles; j++)
    {
        // if only one sphere, then x_bound = x_ins and there is
        // no relevant orientation
        if(1 == nparticles)
            vectorAdd3D(x_ins[j],x,x_ins[j]);
        // if > 1 sphere, take orientation into account
        // x_bound_ins is in the global coordinate system
        else
        {
            vectorSubtract3D(x_ins[j],x_bound_ins,rel);
            MathExtraLiggghts::vec_quat_rotate(rel,quat);
            vectorAdd3D(rel,x,x_ins[j]);
        }
    }

    // set velocity and omega
    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);

    return nparticles;
}

/* ---------------------------------------------------------------------- */

void ParticleToInsert::scale_pti(double r_scale)
{
    double r_scale3 = r_scale*r_scale*r_scale;

    for(int i = 0; i < nparticles; i++)
    {
        radius_ins[i] *= r_scale;
        vectorScalarMult3D(x_ins[i],r_scale);
    }

    volume_ins *= r_scale3;
    mass_ins *= r_scale3;

    r_bound_ins *= r_scale;

    vectorScalarMult3D(x_bound_ins,r_scale);
}
