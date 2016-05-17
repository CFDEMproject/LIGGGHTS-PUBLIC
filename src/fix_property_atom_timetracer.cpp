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

#include "math.h"
#include "stdlib.h"
#include "string.h"
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
#include "fix_property_atom_timetracer.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPropertyAtomTimeTracer::FixPropertyAtomTimeTracer(LAMMPS *lmp, int narg, char **arg,bool parse) :
  FixPropertyAtom(lmp, narg, arg, false),
  iarg_(3),
  check_every_(10)
{
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
        } else if(strcmp(style,"property/atom/timetracer") == 0 )
            error->fix_error(FLERR,this,"unknown keyword");
    }

    // do the base class stuff
    // allocate one extra per-particle for each region

    int n = strlen(id) + 1;
    char *tracer_name = new char[n];
    strcpy(tracer_name,id);
    int n_reg = iregion_.size();
    const char *baseargs[9+n_reg];
    baseargs[0] = tracer_name; 
    baseargs[1] = "all";
    baseargs[2] = "property/atom/tracer";
    baseargs[3] = tracer_name;
    if(0 == n_reg)
        baseargs[4] = "scalar"; 
    else
        baseargs[4] = "vector"; 
    baseargs[5] = "yes";    
    baseargs[6] = "yes";    
    baseargs[7] = "no";    
    baseargs[8] = "0.";
    for(int i = 0; i < n_reg; i++)
        baseargs[9+i] = "0.";
    parse_args(9+n_reg,(char**)baseargs);

    // settings
    nevery = check_every_;

    vector_flag = 1;
    size_vector = 1+n_reg;
    global_freq = check_every_;
    extvector = 1;
}

/* ---------------------------------------------------------------------- */

FixPropertyAtomTimeTracer::~FixPropertyAtomTimeTracer()
{
    for(size_t i = 0; i < idregion_.size(); i++)
        delete [] idregion_[i];
}

/* ----------------------------------------------------------------------
   initialize this fix
------------------------------------------------------------------------- */

void FixPropertyAtomTimeTracer::init()
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

int FixPropertyAtomTimeTracer::setmask()
{
    int mask = FixPropertyAtom::setmask();
    mask |= END_OF_STEP;
    return mask;
}

/* ----------------------------------------------------------------------
   add residence time(s)
------------------------------------------------------------------------- */

void FixPropertyAtomTimeTracer::end_of_step()
{
    int nlocal = atom->nlocal;
    double **x = atom->x;

    // add multiple of dt in case check_region_every > 1
    double dt_this = static_cast<double>(check_every_)*update->dt;
    int n_reg = iregion_.size();
    bool has_regions = n_reg > 0;

    if(!has_regions)
    {
        double *time = this->vector_atom;
        for(int i = 0; i < nlocal; i++)
            time[i] += dt_this;
    }
    else // has_regions
    {
        double **time = this->array_atom;
        Region *region = 0;

        for(int i = 0; i < nlocal; i++)
        {
            time[i][0] += dt_this;
            for(int ireg = 0; ireg < n_reg; ireg++)
            {
                region = domain->regions[iregion_[ireg]];
                if(region->match(x[i][0],x[i][1],x[i][2]))
                    time[i][1+ireg] += dt_this;
            }
        }
    }
}

/* ----------------------------------------------------------------------
   return average of residence time, n = 0..nvalues-1
------------------------------------------------------------------------- */

double FixPropertyAtomTimeTracer::compute_vector(int n)
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  double value = 0.;

  for (int i = 0; i < nlocal; i++)
  {
      if (mask[i] & groupbit)
      {
          if (array_atom) value += array_atom[i][n];
          else value += vector_atom[i];
      }
  }

  MPI_Sum_Scalar(value,world);
  return value/atom->natoms;
}
