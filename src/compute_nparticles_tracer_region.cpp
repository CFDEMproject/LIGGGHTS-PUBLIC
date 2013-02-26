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

#include "stdio.h"
#include "stdlib.h"
#include "compute_nparticles_tracer_region.h"
#include "update.h"
#include "region.h"
#include "group.h"
#include "atom.h"
#include "modify.h"
#include "domain.h"
#include "mpi_liggghts.h"
#include "math_extra_liggghts.h"
#include "fix_property_atom_tracer.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeNparticlesTracerRegion::ComputeNparticlesTracerRegion(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  image_dim_(-1),       // turned off
  image_no_(-1),
  reset_marker_(true),
  iregion_count_(-1),
  idregion_count_(0),
  fix_tracer_(0)
{
    int iarg = 3;

    // parse args

    bool hasargs = true;
    while(iarg < narg && hasargs)
    {
        hasargs = false;
        if(strcmp(arg[iarg],"periodic") == 0) {
            if(narg < iarg+3)
                error->compute_error(FLERR,this,"not enough arguments for 'periodic'");
            iarg++;

            if(strcmp(arg[iarg],"x") == 0)
                image_dim_ = 0;
            else if(strcmp(arg[iarg],"y") == 0)
                image_dim_ = 1;
            else if(strcmp(arg[iarg],"z") == 0)
                image_dim_ = 2;
            else
                error->compute_error(FLERR,this,"'x', 'y' or 'z' expected after 'periodic'");
            iarg++;

            if(strcmp(arg[iarg],"all") == 0)
                image_no_ = -1;
            else image_no_ = atoi(arg[iarg]);

            if(!domain->periodicity[image_dim_])
                error->compute_error(FLERR,this,"using 'periodic' in a dimension which is not periodic");
            iarg++;

            hasargs = true;
        } else if(strcmp(arg[iarg],"reset_marker") == 0) {
            if(narg < iarg+2)
                error->compute_error(FLERR,this,"not enough arguments for 'reset_marker'");
            if(strcmp(arg[iarg+1],"yes") == 0)
                reset_marker_ = true;
            else if(strcmp(arg[iarg+1],"no") == 0)
                reset_marker_ = false;
            else
                error->compute_error(FLERR,this,"expecing 'yes' or 'no' for 'reset_marker'");
            iarg += 2;
            hasargs = true;
        } else if(strcmp(arg[iarg],"tracer") == 0) {
            if(narg < iarg+2)
                error->compute_error(FLERR,this,"not enough arguments for 'tracer'");
            int n = strlen(arg[iarg+1]) + 1;
            fix_tracer_name_ = new char[n];
            strcpy(fix_tracer_name_,arg[iarg+1]);
            iarg += 2;
            hasargs = true;
        } else if(strcmp(arg[iarg],"region_count") == 0) {
            if(narg < iarg+2)
                error->compute_error(FLERR,this,"not enough arguments for 'region_count'");
            iregion_count_ = domain->find_region(arg[iarg+1]);
            if(iregion_count_ == -1)
                error->compute_error(FLERR,this,"Region ID does not exist");
            int n = strlen(arg[iarg+1]) + 1;
            idregion_count_ = new char[n];
            strcpy(idregion_count_,arg[iarg+1]);
            iarg += 2;
            hasargs = true;
        } else
            error->compute_error(FLERR,this,"unknown keyword");
    }

    if(!fix_tracer_name_)
        error->compute_error(FLERR,this,"have to define 'tracer'");
    if(iregion_count_ < 0)
        error->compute_error(FLERR,this,"have to define 'region_count'");

    scalar_flag = 1;
}

/* ---------------------------------------------------------------------- */

ComputeNparticlesTracerRegion::~ComputeNparticlesTracerRegion()
{
    delete []idregion_count_;
    delete []fix_tracer_name_;
}

/* ---------------------------------------------------------------------- */

void ComputeNparticlesTracerRegion::init()
{
    iregion_count_ = domain->find_region(idregion_count_);
    if (iregion_count_ == -1)
        error->compute_error(FLERR,this,"Region ID (region_count) does not exist");

    FixPropertyAtom *fix_ppa;
    
    fix_ppa = static_cast<FixPropertyAtom*>(modify->find_fix_property(fix_tracer_name_,"property/atom","scalar",0,0,style));

    fix_tracer_ = dynamic_cast<FixPropertyAtomTracer*>(fix_ppa);
    if(!fix_tracer_)
        error->compute_error(FLERR,this,"need a tracer fix of type fix property/atom/tracer");
}

/* ---------------------------------------------------------------------- */

double ComputeNparticlesTracerRegion::compute_scalar()
{
    invoked_scalar = update->ntimestep;

    if(image_dim_ == -1)
        scalar = compute_scalar_eval<false>();
    else
        scalar = compute_scalar_eval<true>();

    return scalar;
}

/* ---------------------------------------------------------------------- */

template<bool IMAGE>
double ComputeNparticlesTracerRegion::compute_scalar_eval()
{
    int ncount;
    int nlocal = atom->nlocal;
    int *image = atom->image;
    double **x = atom->x;
    double *marker = fix_tracer_->vector_atom;
    Region *region = domain->regions[iregion_count_];

    // count all particles in region, taking image flag into account
    ncount = 0;
    for(int i = 0; i < nlocal; i++)
    {
        
        if(IMAGE && image_dim_ == 0 && ( ((image[i] & 1023) - 512      ) != image_no_) )
            continue;
        if(IMAGE && image_dim_ == 1 && ( ((image[i] >> 10 & 1023) - 512) != image_no_) )
            continue;
        if(IMAGE && image_dim_ == 2 && ( ((image[i] >> 20) - 512) != image_no_) )
            continue;

        if(MathExtraLiggghts::compDouble(marker[i]  ,1.0,1e-5) &&
           region->match(x[i][0],x[i][1],x[i][2]))
        {
            ncount++;
            
            if(reset_marker_)
                marker[i] = 0.0;
        }
    }

    MPI_Sum_Scalar(ncount,world);
    return static_cast<double>(ncount);
}
