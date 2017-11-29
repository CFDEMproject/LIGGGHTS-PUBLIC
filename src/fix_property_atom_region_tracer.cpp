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
    (if no contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2013-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#include <cmath>
#include <stdlib.h>
#include <string.h>
#include "atom.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "group.h"
#include "region.h"
#include "domain.h"
#include "mpi_liggghts.h"
#include "math_extra_liggghts.h"
#include "fix_property_atom_region_tracer.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPropertyAtomRegionTracer::FixPropertyAtomRegionTracer(LAMMPS *lmp, int narg, char **arg,bool parse) :
  FixPropertyAtom(lmp, narg, arg, false),
  iarg_(3),
  check_every_(10)
{
    
    if (strcmp(style, "property/atom/timetracer") == 0)
        error->fix_error(FLERR, this, "Style of fix property/atom/timetracer is now property/atom/regiontracer/time");
    // do the derived class stuff

    // parse args

    bool hasargs = true;
    while(iarg_ < narg && hasargs)
    {
        hasargs = false;

        if(strcmp(arg[iarg_],"add_region") == 0) {
            if(narg < iarg_+2)
                error->fix_error(FLERR,this,"not enough arguments for 'add_region'");
            iarg_++;
            int n = strlen(arg[iarg_]) + 1;
            char *idreg = new char[n];
            strcpy(idreg,arg[iarg_]);
            int ireg = domain->find_region(arg[iarg_++]);
            if (ireg == -1)
                error->fix_error(FLERR,this,"Region ID does not exist");
            iregion_.push_back(ireg);
            idregion_.push_back(idreg);
            hasargs = true;
        }  else if(strcmp(arg[iarg_],"check_region_every") == 0) {
            if(narg < iarg_+2)
                error->fix_error(FLERR,this,"not enough arguments for 'check_region_every'");
            iarg_++;
            check_every_ = atoi(arg[iarg_]);
            if(check_every_ < 0)
                error->fix_error(FLERR,this,"check_region_every > 0 required");
            iarg_++;
            hasargs = true;
        } else if(strcmp(style,"property/atom/regiontracer/time") == 0 )
            error->fix_error(FLERR,this,"unknown keyword");
    }

    // do the base class stuff
    // allocate one extra per-particle for each region

    int n = strlen(id) + 1;
    char *tracer_name = new char[n];
    strcpy(tracer_name,id);
    const int n_reg = iregion_.size();
    if (sizeof(double) == 8 && n_reg > 53)
        error->fix_error(FLERR, this, "Only 53 regions are allowed in a fix property/atom/region/tracer (on systems with 64 bit doubles)");
    else if (sizeof(double) == 4 && n_reg > 24)
        error->fix_error(FLERR, this, "Only 24 regions are allowed in a fix property/atom/region/tracer (on systems with 32 bit doubles)");

    const int num_args = 9 + (n_reg > 0 ? n_reg+1 : 0);
    char **baseargs = new char*[num_args]; // VS does not support variable-length arrays --> new
    baseargs[0] = tracer_name; 
    baseargs[1] = (char *) "all";
    baseargs[2] = (char *) "property/atom/tracer";
    baseargs[3] = tracer_name;
    if(0 == n_reg)
        baseargs[4] = (char *) "scalar"; 
    else
        baseargs[4] = (char *) "vector"; 
    baseargs[5] = (char *) "yes";   
    baseargs[6] = (char *) "yes";   
    baseargs[7] = (char *) "no";    
    baseargs[8] = (char *) "0.";            // Property for simulation domain
    for(int i = 0; i < n_reg; i++)
        baseargs[9+i] = (char *) "0.";      // Property for region
    if (n_reg > 0)
        baseargs[9+n_reg] = (char *) "0.";  // Bit map for region
    parse_args(num_args,(char**)baseargs);

    delete []baseargs;

    // settings
    nevery = 1;

    vector_flag = 1;
    size_vector = 1+n_reg;
    global_freq = check_every_;
    extvector = 1;
}

/* ---------------------------------------------------------------------- */

FixPropertyAtomRegionTracer::~FixPropertyAtomRegionTracer()
{
    for(size_t i = 0; i < idregion_.size(); i++)
        delete [] idregion_[i];
}

/* ----------------------------------------------------------------------
   initialize this fix
------------------------------------------------------------------------- */

void FixPropertyAtomRegionTracer::init()
{
    iregion_.clear();
    for(size_t i = 0; i < idregion_.size(); i++)
    {
        int ireg = domain->find_region(idregion_[i]);
        if (ireg == -1)
            error->fix_error(FLERR,this,"Region ID does not exist");
        iregion_.push_back(ireg);
    }

    int n_ms = modify->n_fixes_style("multisphere");
    if(n_ms > 0)
        error->fix_error(FLERR,this,"may not be used together with fix multisphere");
}

/* ---------------------------------------------------------------------- */

int FixPropertyAtomRegionTracer::setmask()
{
    int mask = FixPropertyAtom::setmask();
    mask |= END_OF_STEP;
    return mask;
}

/* ----------------------------------------------------------------------
   add residence time(s)
------------------------------------------------------------------------- */

void FixPropertyAtomRegionTracer::end_of_step()
{
    const int nlocal = atom->nlocal;
    double **x = atom->x;
    const int * const mask = atom->mask;

    // add multiple of dt in case check_region_every > 1
    const double dt = update->dt;
    const int n_reg = iregion_.size();

    if (n_reg == 0)
    {
        double * const data = this->vector_atom;
        for(int i = 0; i < nlocal; i++)
        {
            if (mask[i] & groupbit)
                data[i] += trace_multiplier(i)*dt;
        }
    }
    else // fix has region list
    {
        double * const * const data = this->array_atom;

        for(int i = 0; i < nlocal; i++)
        {
            if (mask[i] & groupbit)
            {
                const double add = trace_multiplier(i)*dt;
                data[i][0] += add;

                long region_map = (long)data[i][n_reg+1];
                for(int ireg = 0; ireg < n_reg; ireg++)
                {
                    const long region_flag = 1<<ireg;
                    // check if we need to update the region list
                    if (update->ntimestep % check_every_ == 0)
                    {
                        Region * const region = domain->regions[iregion_[ireg]];
                        if(region->match(x[i][0],x[i][1],x[i][2]))
                            region_map |= region_flag;
                        else
                            region_map &= ~region_flag;
                    }
                    if (region_map & region_flag)
                        data[i][1+ireg] += add;
                }
                data[i][n_reg+1] = static_cast<double>(region_map);
            }
        }
    }
}

/* ----------------------------------------------------------------------
   return average of residence time, n = 0..nvalues-1
------------------------------------------------------------------------- */

double FixPropertyAtomRegionTracer::compute_vector(int n)
{
    const int nlocal = atom->nlocal;
    const int * const mask = atom->mask;

    double value = 0.;
    double norm = 0.;

    for (int i = 0; i < nlocal; i++)
    {
        if (mask[i] & groupbit)
        {
            const double incr = array_atom ? array_atom[i][n] : vector_atom[i];
            // sum all values
            value += incr;
            // increment the normalizing factor
            if (incr > 0.0)
                norm += 1.0;
        }
    }

    MPI_Sum_Scalar(value,world);
    MPI_Sum_Scalar(norm,world);

    if (norm > 0.5)
        return value/atom->natoms;
    else
        return 0.0;
}
