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

#include <string.h>
#include <stdlib.h>
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include <cmath>
#include "comm.h"
#include "vector_liggghts.h"
#include "fix_cfd_coupling.h"
#include "fix.h"
#include "fix_multisphere.h"
#include "fix_property_global.h"
#include "fix_property_atom.h"
#include "fix_cfd_coupling.h"
#include "style_cfd_datacoupling.h"
#include "cfd_regionmodel.h"
#include "style_cfd_regionmodel.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCfdCoupling::FixCfdCoupling(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  iarg_ = 3;

  rm_ = NULL;
  dc_ = NULL;

  // end_of_step is executed each ts
  nevery = 1;

  ts_create_ = update->ntimestep;

  couple_this_ = 0;
  couple_nevery_ = 0;

  // parse args

  if (narg < 4)
    error->fix_error(FLERR,this,"");

  if(0 == strcmp(arg[iarg_],"every") || 0 == strcmp(arg[iarg_],"couple_every"))
  {
    iarg_++;

    if (narg < 6)
        error->fix_error(FLERR,this,"not enough arguments");

    couple_nevery_ = atoi(arg[iarg_++]);
    if(couple_nevery_ < 0)
        error->fix_error(FLERR,this,"'every' value must be >=0");
  }
  // construct data coupling submodel and parse its args

  if (0) return;
  #define CFD_DATACOUPLING_CLASS
  #define CfdDataCouplingStyle(key,Class) \
  else if (strcmp(arg[iarg_],#key) == 0) dc_ = new Class(lmp,iarg_+1,narg,arg,this);
  #include "style_cfd_datacoupling.h"
  #undef CFD_DATACOUPLING_CLASS
  else error->fix_error(FLERR,this,"Unknown data coupling style - expecting 'file' or 'MPI'");

  if(!dynamic_cast<CfdDatacouplingMPI*>(dc_) && 0 == couple_nevery_)
    error->fix_error(FLERR,this,"expecting keyword 'couple_every' ");
  else if(dynamic_cast<CfdDatacouplingMPI*>(dc_) && 0 != couple_nevery_)
  {
    if(comm->me == 0) error->message(FLERR,"couple_every as specified in LIGGGHTS is overriden by calling external program",1);
  }

  iarg_ = dc_->get_iarg();

  bool hasargs = true;
  while (iarg_ < narg && hasargs)
  {
      hasargs = false;
      if(strcmp(arg[iarg_],"regionmodel") == 0)
      {
          hasargs = true;
          iarg_++;
          if (0) return;
          #define CFD_REGIONMODEL_CLASS
          #define CfdRegionStyle(key,Class) \
          else if (strcmp(arg[iarg_],#key) == 0) rm_ = new Class(lmp,iarg_+1,narg,arg,this);
          #include "style_cfd_regionmodel.h"
          #undef CFD_REGIONMODEL_CLASS
          else error->fix_error(FLERR,this,"Unknown cfd regionmodel style");
          iarg_ = rm_->get_iarg();
      }
  }

}

/* ---------------------------------------------------------------------- */

FixCfdCoupling::~FixCfdCoupling()
{
    delete rm_;
    delete dc_;
}

/* ---------------------------------------------------------------------- */

void FixCfdCoupling::post_create()
{
    if(dc_) dc_->post_create();
}

/* ---------------------------------------------------------------------- */

int FixCfdCoupling::setmask()
{
    int mask = 0;
    mask |= END_OF_STEP;
    mask |= POST_FORCE_RESPA;
    mask |= MIN_POST_FORCE;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCoupling::init()
{
    if(0 == atom->map_style)
      error->fix_error(FLERR,this,"requires an 'atom_modify map' command to allocate an atom map");

    // make sure there is only one cfdcoupling fix
    if(modify->n_fixes_style_strict(style) != 1)
      error->fix_error(FLERR,this,"More than one fix of style couple/cfd is not allowed");

    if (strcmp(update->integrate_style,"respa") == 0)
      nlevels_respa = ((Respa *) update->integrate)->nlevels;

    if(rm_) rm_->init();

    if(dc_) dc_->init();
}

/* ---------------------------------------------------------------------- */

void FixCfdCoupling::setup(int vflag)
{
    if (strstr(update->integrate_style,"verlet"))
      post_force(vflag);
    else {
      ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
      post_force_respa(vflag,nlevels_respa-1,0);
      ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
    }

    if(update->ntimestep == 0) end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixCfdCoupling::end_of_step()
{
    // only do this if coupling desired
    if(couple_nevery_ == 0) return;

    int ts = update->ntimestep;

    // set flag if couple the next time-step
    // will not be active this ts, since in eos
    if((ts+1) % couple_nevery_ || ts_create_ == ts+1)
        couple_this_ = 0;
    else
        couple_this_ = 1;

    // only execute if pushing or pulling is desired
    // do not execute on step of creation
    if(ts % couple_nevery_ || ts_create_ == ts) return;

    if(screen && comm->me == 0)
        fprintf(screen,"CFD Coupling established at step %d\n",ts);
    if(logfile && comm->me == 0)
        fprintf(logfile,"CFD Coupling established at step %d\n",ts);

    // call region model
    if(rm_) rm_->rm_update();

    // call data exchange model to exchane data
    dc_->exchange();

    // check if datatransfer was successful
    // dc_->check_datatransfer();
}

/* ---------------------------------------------------------------------- */

void FixCfdCoupling::min_setup(int vflag)
{
    post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixCfdCoupling::post_force_respa(int vflag, int ilevel, int)
{
    if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixCfdCoupling::min_post_force(int vflag)
{
    post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixCfdCoupling::add_push_property(const char *name, const char *type)
{
    dc_->add_push_property(name,type);
}

/* ---------------------------------------------------------------------- */

void FixCfdCoupling::add_pull_property(const char *name, const char *type)
{
    dc_->add_pull_property(name,type);
}
