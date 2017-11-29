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
    Stefan Radl, TU Graz
    Copyright 2016 - TU Graz
------------------------------------------------------------------------- */

#include <cmath>
#include <stdlib.h>
#include <string.h>
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "modify.h"
#include "force.h"
#include "memory.h"
#include "update.h"
#include "error.h"
#include "fix_mesh_surface.h"
#include "fix_neighlist_mesh.h"
#include "fix_multisphere.h"
#include "fix_property_atom.h"
#include "mpi_liggghts.h"
#include "math_extra_liggghts.h"
#include "fix_massflow_mesh_sieve.h"
#include "math_const.h"
#include "irregular.h"
#include "atom_vec_ellipsoid.h"
#include "random_park.h"

#include <algorithm>

using namespace LAMMPS_NS;
using namespace MathExtraLiggghts;
using namespace MathConst;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMassflowMeshSieve::FixMassflowMeshSieve(LAMMPS *lmp, int narg, char **arg) :
  FixMassflowMesh(lmp, narg, arg),
  sieveMultiSphereCanPass_(false),
  verbose_(false), //developer to hard code here
  sieveSize_(-1),
  sieveSpacing_(-1),
  sieveStiffness_(-1),
  sieveDamping_(0),
  fix_applyForce_(0),
  random_(0)
{

    // random number generator, seed is hardcoded
    random_ = new RanPark(lmp,"15485863");

    bool hasargs = true;

    while(iarg_ < narg && hasargs)
    {
        hasargs = false;

        if(strcmp(arg[iarg_],"sieveMultiSphereCanPass") == 0) {
            sieveMultiSphereCanPass_     = true;
            iarg_++;
            hasargs = true;
        }
        else if ( strcmp(arg[iarg_],"sieveSize") == 0) {
            if(narg < iarg_+2)
                error->fix_error(FLERR,this,"not enough arguments for 'sieveSize'");
            iarg_++;
            sieveSize_ = atof(arg[iarg_++]);
            hasargs = true;
        }
        else if ( strcmp(arg[iarg_],"sieveSpacing") == 0) {
            if(narg < iarg_+2)
                error->fix_error(FLERR,this,"not enough arguments for 'sieveSpacing'");
            iarg_++;
            sieveSpacing_ = atof(arg[iarg_++]);
            hasargs = true;
        }
        else if ( strcmp(arg[iarg_],"sieveStiffness") == 0) {
            if(narg < iarg_+2)
                error->fix_error(FLERR,this,"not enough arguments for 'sieveStiffness'");
            iarg_++;
            sieveStiffness_ = atof(arg[iarg_++]);
            hasargs = true;
        }
        else if ( strcmp(arg[iarg_],"sieveDamping") == 0) {
            if(narg < iarg_+2)
                error->fix_error(FLERR,this,"not enough arguments for 'sieveDamping'");
            iarg_++;
            sieveDamping_ = atof(arg[iarg_++]);
            hasargs = true;
        }
        else
          error->fix_error(FLERR,this,"unknown argument.");
    }

    if(sieveSize_<0 || sieveSize_> 1e5)
        error->fix_error(FLERR,this,"'sieveSize' is not set, larger than 1e5, or less than zero.");

    if(sieveSpacing_<0 || sieveSpacing_> 1e5)
        error->fix_error(FLERR,this,"'sieveSpacing' is not set, larger than 1e5, or less than zero.");

    if(sieveStiffness_<0 || sieveStiffness_> 1e10)
        error->fix_error(FLERR,this,"'sieveStiffness' is not set, larger than 1e10, or less than zero.");
}

/* ---------------------------------------------------------------------- */

FixMassflowMeshSieve::~FixMassflowMeshSieve()
{
    if (random_) delete random_;
}

/* ---------------------------------------------------------------------- */
void FixMassflowMeshSieve::post_create()
{

    FixMassflowMesh::post_create();

    // add per-particle counter for displacements
    const char * fixarg[9];
    sprintf(fixidApplyForce_,"massflowSieve_%s",id);
    fixarg[0]=fixidApplyForce_;
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]=fixidApplyForce_;
    fixarg[4]="scalar";
    fixarg[5]="yes";
    fixarg[6]="no";
    fixarg[7]="no";
    fixarg[8]="-1";
    modify->add_fix(9,const_cast<char**>(fixarg));

    fix_applyForce_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(fixidApplyForce_,"property/atom","scalar",0,0,style));

}

/* ----------------------------------------------------------------------
   initialize this fix
------------------------------------------------------------------------- */
void FixMassflowMeshSieve::init()
{
    FixMassflowMesh::init();
}

/* ---------------------------------------------------------------------- */
int FixMassflowMeshSieve::setmask()
{
    int mask = FixMassflowMesh::setmask();
    mask |= POST_FORCE;
    mask |= PRE_EXCHANGE; //ensure pre_exchange always happens
    return mask;
}

/* ----------------------------------------------------------------------
   perform additinal force computations
------------------------------------------------------------------------- */
void FixMassflowMeshSieve::post_force(int vflag)
{
    int nlocal = atom->nlocal;
    double **x = atom->x;
    double **v = atom->v;
    double **f = atom->f;

    double *radius = atom->radius;
    int    *mask = atom->mask;
    double *counter = fix_counter_->vector_atom;
    double *applyForce = fix_applyForce_->vector_atom;

    double delta[3]={}, bary[3];
    double deltan;

    TriMesh *mesh = fix_mesh_->triMesh();
    int nTriAll = mesh->sizeLocal() + mesh->sizeGhost();

    // loop owned and ghost triangles
    // count only if owned particle
    for(int iTri = 0; iTri < nTriAll; iTri++)
    {
        const std::vector<int> & neighborList = fix_neighlist_->get_contact_list(iTri);
        const int numneigh = neighborList.size();
        double iTriDouble = double(iTri);

        for(int iNeigh = 0; iNeigh < numneigh; iNeigh++)
        {
            const int iPart = neighborList[iNeigh];

            // skip ghost particles
            if(iPart >= nlocal)
                continue;

            // skip particles not in fix group
            if (!(mask[iPart] & groupbit))
                continue;

            // in case of once_ == true, ignore everything which has been already counted
            if(compDouble(counter[iPart],2.))
                continue;

            int barySign;
            deltan = fix_mesh_->triMesh()->resolveTriSphereContactBary(iPart,iTri,radius[iPart],x[iPart],delta,bary,barySign);

            if(deltan <= 0)
            {
                //Check if first contact, and throw coin to determine passage
                if(applyForce[iPart]<0.0) //have first contact
                {
                    double randNumber = random_->uniform();
                    if(!sieveMultiSphereCanPass_ && (fix_volumeweight_ms_ && fix_volumeweight_ms_->vector_atom[iPart] > 0)) //have an unwanted multisphere collision --> apply force in any case!
                        applyForce[iPart] = 1+iTriDouble;
                    else if( randNumber  < sievePassProbability(radius[iPart]) )
                    {
                        applyForce[iPart] = 0;
                        continue;
                    }
                    else
                        applyForce[iPart] = 1+iTriDouble;

//                    printf("FixMassflowMeshSieve:: particle %d has an overlap with deltan: %g, delta: %g %g %g, radius: %g, applyForce: %g, randNumber: %g \n",
//                       iPart, deltan, delta[0],delta[1],delta[2],radius[iPart], applyForce[iPart],randNumber);

                }
                else if(compDouble(applyForce[iPart],0)) continue; //suppress force calculation

                double contactNormal[3];
                double invRadius = 1/radius[iPart];
                contactNormal[0] = delta[0] * invRadius;
                contactNormal[1] = delta[1] * invRadius;
                contactNormal[2] = delta[2] * invRadius;

                double
                normalVel = v[iPart][0] * contactNormal[0]
                          + v[iPart][1] * contactNormal[1]
                          + v[iPart][2] * contactNormal[2];

                const double Fn_damping = sieveDamping_  * normalVel;
                const double Fn_contact =-sieveStiffness_* deltan; //deltan is NEGATIVE!
                double Fn = Fn_damping + Fn_contact;
                if(Fn<0.0) //avoid attractive forces during contact
                  Fn = 0.0;

                f[iPart][0] -= Fn*contactNormal[0];
                f[iPart][1] -= Fn*contactNormal[1];
                f[iPart][2] -= Fn*contactNormal[2];
            }
            else //particle has no contact with this triangle now - but check if it had previously
            {
                if( compDouble(applyForce[iPart],1+iTriDouble) )  //previously a force was applied with this triangle, now release it!
                    applyForce[iPart] = -1;
                //else --> have force with other triangle, or no interaction. Thus, no action for this triangle
            }
        }
    }

}
/* ----------------------------------------------------------------------
   perform some operations here if needed
   done before exchange, borders, reneighbor
   so that ghost atoms and neighbor lists will be correct
------------------------------------------------------------------------- */
void FixMassflowMeshSieve::pre_exchange()
{
    FixMassflowMesh::pre_exchange();
}

double  FixMassflowMeshSieve::sievePassProbability(double radius)
{
    double dPass = sieveSize_ - 2.0 * radius;
    if( dPass > 0.0  )
        return dPass * dPass * PIOVERFOUR / (sieveSpacing_ * sieveSpacing_);
    else
        return 0;
};

