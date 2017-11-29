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
#include "neighbor.h"
#include "mpi_liggghts.h"
#include "math_extra_liggghts.h"
#include "fix_property_atom_tracer.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPropertyAtomTracer::FixPropertyAtomTracer(LAMMPS *lmp, int narg, char **arg,bool parse) :
  FixPropertyAtom(lmp, narg, arg, false),
  iarg_(3),
  marker_style_(MARKER_DIRAC),
  step_(-1),
  check_every_(10),
  first_mark_(true),
  iregion_(-1),
  idregion_(0),
  nmarked_last_(0),
  nmarked_(0)
{
    // do the base class stuff

    int n = strlen(id) + 1;
    tracer_name_ = new char[n];
    strcpy(tracer_name_,id);
    const char *baseargs[9];
    baseargs[0] = tracer_name_; 
    baseargs[1] = "all";
    baseargs[2] = "property/atom/tracer";
    baseargs[3] = tracer_name_;
    baseargs[4] = "scalar"; 
    baseargs[5] = "yes";    
    baseargs[6] = "yes";    
    baseargs[7] = "no";    
    baseargs[8] = "0.";
    parse_args(9,(char**)baseargs);

    // do the derived class stuff

    // parse args

    bool hasargs = true;
    while(iarg_ < narg && hasargs)
    {
        hasargs = false;

        if(strcmp(arg[iarg_],"region_mark") == 0) {
            if(narg < iarg_+2)
                error->fix_error(FLERR,this,"not enough arguments for 'region_mark'");
            iarg_++;
            n = strlen(arg[iarg_]) + 1;
            idregion_ = new char[n];
            strcpy(idregion_,arg[iarg_]);
            iregion_ = domain->find_region(arg[iarg_++]);
            if (iregion_ == -1)
                error->fix_error(FLERR,this,"Region ID does not exist");
            hasargs = true;
        } else if(strcmp(arg[iarg_],"mark_step") == 0) {
            if(narg < iarg_+2)
                error->fix_error(FLERR,this,"not enough arguments for 'mark_step'");
            iarg_++;
            step_ = atoi(arg[iarg_++]);
            if(step_ < 0)
                error->fix_error(FLERR,this,"mark_step > 0 required");
            
            if(step_ < update->ntimestep)
                first_mark_ = false;
            hasargs = true;
        } else if(strcmp(arg[iarg_],"marker_style") == 0) {
            if(narg < iarg_+2)
                error->fix_error(FLERR,this,"not enough arguments for 'marker_style'");
            iarg_++;
            if(strcmp(arg[iarg_],"heaviside") == 0)
                marker_style_ = MARKER_HEAVISIDE;
            else if(strcmp(arg[iarg_],"dirac") == 0)
                marker_style_ = MARKER_DIRAC;
            else if(strcmp(arg[iarg_],"none") == 0)
                marker_style_ = MARKER_NONE;
            else
                error->fix_error(FLERR,this,"expecting 'heaviside' or 'dirac' after keyword 'marker_style'");
            iarg_++;
            hasargs = true;
        } else if(strcmp(arg[iarg_],"check_mark_every") == 0) {
            if(narg < iarg_+2)
                error->fix_error(FLERR,this,"not enough arguments for 'check_mark_every'");
            iarg_++;
            check_every_ = atoi(arg[iarg_]);
            if(check_every_ < 0)
                error->fix_error(FLERR,this,"check_mark_every > 0 required");
            iarg_++;
            hasargs = true;
        } else if(strcmp(style,"property/atom/tracer") == 0 )
            error->fix_error(FLERR,this,"unknown keyword");
    }

    // error checks for this class

    if(strcmp(style,"property/atom/tracer") == 0)
    {
        if (iregion_ == -1)
            error->fix_error(FLERR,this,"expecting keyword 'region_mark'");
        if (step_ == -1)
            error->fix_error(FLERR,this,"expecting keyword 'mark_step'");
    }

    // settings
    nevery = check_every_;

    scalar_flag = 1;
    global_freq = 1;
}

/* ---------------------------------------------------------------------- */

FixPropertyAtomTracer::~FixPropertyAtomTracer()
{
    delete []tracer_name_;
    if(idregion_) delete []idregion_;
}

/* ----------------------------------------------------------------------
   initialize this fix
------------------------------------------------------------------------- */

void FixPropertyAtomTracer::init()
{
    iregion_ = domain->find_region(idregion_);
    if (iregion_ == -1)
        error->fix_error(FLERR,this,"Region ID does not exist");
}

/* ---------------------------------------------------------------------- */

int FixPropertyAtomTracer::setmask()
{
    int mask = FixPropertyAtom::setmask();
    mask |= END_OF_STEP;
    return mask;
}

/* ----------------------------------------------------------------------
   mark particles
------------------------------------------------------------------------- */

void FixPropertyAtomTracer::end_of_step()
{
    int ts = update->ntimestep;

    if(ts < step_ || marker_style_ == MARKER_NONE || (marker_style_ == MARKER_DIRAC && !first_mark_))
        return;

    int nlocal = atom->nlocal;
    int *mask = atom->mask;
    double **x = atom->x;
    double *marker = this->vector_atom;
    Region *region = domain->regions[iregion_];

    int nmarked_this = 0;

    for(int i = 0; i < nlocal; i++)
    {
        if (!(mask[i] & groupbit))
            continue;

        if (!MathExtraLiggghts::compDouble(marker[i],1.0,1e-5) && region->match(x[i][0],x[i][1],x[i][2]))
        {
            marker[i] = 1.0;
            nmarked_this++;
        }
    }

    MPI_Sum_Scalar(nmarked_this,world);
    nmarked_ += nmarked_this;

    first_mark_ = false;
}

/* ----------------------------------------------------------------------
   return # marked particles
------------------------------------------------------------------------- */

double FixPropertyAtomTracer::compute_scalar()
{
    int result = nmarked_ - nmarked_last_;
    nmarked_last_ = nmarked_;
    return static_cast<double>(result);
}
