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

    Alexander Podlozhnyuk (DCS Computing GmbH, Linz)

    Copyright 2015-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifdef SUPERQUADRIC_ACTIVE_FLAG

#include "math_extra_liggghts.h"
#include "particleToInsertSuperquadric.h"
#include "particleToInsert.h"
#include <cmath>
#include "error.h"
#include "update.h"
#include "domain.h"
#include "atom.h"
#include "atom_vec.h"
#include "fix.h"
#include "vector_liggghts.h"
#include "modify.h"
#include "math_extra_liggghts_nonspherical.h"

using namespace LAMMPS_NS;

ParticleToInsertSuperquadric::ParticleToInsertSuperquadric(LAMMPS* lmp,int ns) : ParticleToInsert(lmp, ns)
{
}

/* ---------------------------------------------------------------------- */

ParticleToInsertSuperquadric::~ParticleToInsertSuperquadric()
{
}

/* ---------------------------------------------------------------------- */

int ParticleToInsertSuperquadric::insert()
{
    // perform the actual insertion
    // add particles, set coordinate and radius
    // set group mask to "all" plus fix groups

    int inserted = 0;
    int nfix = modify->nfix;
    Fix **fix = modify->fix;

    for(int i = 0; i < nparticles; i++) {

                inserted++;
                if(atom_type_vector_flag)
                    atom->avec->create_atom(atom_type_vector[i],x_ins[i]);
                else
                    atom->avec->create_atom(atom_type,x_ins[i]);
                int m = atom->nlocal - 1;
                atom->mask[m] = 1 | groupbit;
                vectorCopy3D(v_ins,atom->v[m]);
                vectorCopy3D(omega_ins,atom->omega[m]);
                vectorCopy4D(quat_ins, atom->quaternion[m]);
                atom->radius[m] = radius_ins[0];
                atom->density[m] = density_ins;
                atom->rmass[m] = mass_ins;

                atom->volume[m] = volume_ins;
                atom->area[m] = area_ins;
//Superquadric bonus-------------------------------------------------
                vectorCopy3D(shape_ins, atom->shape[m]);
                vectorCopy3D(inertia_ins, atom->inertia[m]);

                atom->blockiness[m][0] = blockiness_ins[0];
                atom->blockiness[m][1] = blockiness_ins[1];
                MathExtraLiggghtsNonspherical::omega_to_angmom(atom->quaternion[m], atom->omega[m], atom->inertia[m],atom->angmom[m]);
//-------------------------------------------------------------------
                //pre_set_arrays() called via FixParticleDistribution
                for (int j = 0; j < nfix; j++)
                   if (fix[j]->create_attribute) fix[j]->set_arrays(m);
        //}
    }

    return inserted;
}

/* ---------------------------------------------------------------------- */

int ParticleToInsertSuperquadric::check_near_set_x_v_omega(double *x,double *v, double *omega, double *quat, RegionNeighborList<interpolate_no> & neighList)
{
    // check sphere against all others in xnear
    // if no overlap add to xnear

    if(nparticles > 1)
        error->one(FLERR,"check_near_set_x_v_omega not implemented yet for nparticles>1");

    if(neighList.hasOverlap_superquadric(x, radius_ins[0], quat, shape_ins, blockiness_ins)) {
      return 0;
    }

    // no overlap with any other - success

    vectorCopy3D(x,x_ins[0]);
    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);
    vectorCopy4D(quat, quat_ins);

    // add to xnear
    neighList.insert_superquadric(x_ins[0], radius_ins[0], quat_ins, shape_ins, blockiness_ins);

    return 1;
}

/* ---------------------------------------------------------------------- */

int ParticleToInsertSuperquadric::set_x_v_omega(double *x, double *v, double *omega, double *quat)
{
    // set velocity and omega
    vectorCopy3D(x,x_ins[0]);
    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);
    vectorCopy4D(quat, quat_ins);

    return nparticles;
}

/* ---------------------------------------------------------------------- */

void ParticleToInsertSuperquadric::scale_pti(double r_scale)
{
    double r_scale3 = r_scale*r_scale*r_scale;

    for(int i = 0; i < nparticles; i++) {
        radius_ins[i] *= r_scale;
        shape_ins[0] *= r_scale;
        shape_ins[1] *= r_scale;
        shape_ins[2] *= r_scale;
        vectorScalarMult3D(x_ins[i],r_scale);
    }

    volume_ins *= r_scale3;
    mass_ins *= r_scale3;

    r_bound_ins *= r_scale;
}

#endif
