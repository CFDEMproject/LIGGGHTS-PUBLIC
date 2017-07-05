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

    Arno Mayrhofer (DCS Computing GmbH)

    Copyright 2016-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#include "mesh_module_liquidtransfer.h"
#include <stdio.h>
#include <string.h>
#include "error.h"
#include "force.h"
#include "modify.h"
#include "comm.h"
#include "pair_gran.h"
#include "math_extra.h"
#include "fix_property_global.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

MeshModuleLiquidTransfer::MeshModuleLiquidTransfer(LAMMPS *lmp, int &iarg_, int narg, char **arg, FixMeshSurface *fix_mesh)
: MeshModule(lmp, iarg_, narg, arg, fix_mesh),
  liquid_content_(0),
  liquid_flux_(0),
  limit_liquid_content_(false),
  max_liquid_content_(0.0),
  wall_thickness_(0.),
  initial_liquid_content_(0.)
{
    // parse further args

    bool hasargs = true;
    while(iarg_ < narg && hasargs)
    {
        hasargs = false;

        if(strcmp(arg[iarg_],"wall_thickness") == 0) {
          if (iarg_+2 > narg) error->one(FLERR,"not enough arguments for keyword 'wall_thickness'");
          wall_thickness_ = force->numeric(FLERR,arg[iarg_+1]);
          iarg_ += 2;
          hasargs = true;
        } else if(strcmp(arg[iarg_],"initial_liquid_content") == 0) {
          if (iarg_+2 > narg) error->one(FLERR,"not enough arguments for keyword 'initial_liquid_content'");
          initial_liquid_content_ = force->numeric(FLERR,arg[iarg_+1]);
          iarg_ += 2;
          hasargs = true;
        }
    }

    if(wall_thickness_ < 1e-12)
        error->one(FLERR,"have to define 'wall_thickness'");
    if(initial_liquid_content_ < 1e-12)
        error->one(FLERR,"have to define 'initial_liquid_content'");
}

/* ---------------------------------------------------------------------- */

MeshModuleLiquidTransfer::~MeshModuleLiquidTransfer()
{

}

/* ---------------------------------------------------------------------- */

void MeshModuleLiquidTransfer::post_create_pre_restart()
{
    
    mesh->prop().addElementProperty<ScalarContainer<double> >("LiquidContent","comm_forward","frame_invariant","restart_yes");
    mesh->prop().addElementProperty<ScalarContainer<double> >("LiquidFlux","comm_reverse","frame_invariant","restart_no");
}

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void MeshModuleLiquidTransfer::post_create()
{
    
    //Np --> set values for no-restart properties here

    mesh->prop().getElementProperty<ScalarContainer<double> >("LiquidContent")->setAll(initial_liquid_content_);
    mesh->prop().getElementProperty<ScalarContainer<double> >("LiquidFlux")->setAll(0.);
}

/* ---------------------------------------------------------------------- */

void MeshModuleLiquidTransfer::init()
{
    liquid_content_ = mesh->prop().getElementProperty<ScalarContainer<double> >("LiquidContent");
    liquid_flux_ = mesh->prop().getElementProperty<ScalarContainer<double> >("LiquidFlux");
    if(!liquid_content_ || ! liquid_flux_)
      error->one(FLERR,"Internal error");
}

/* ---------------------------------------------------------------------- */

int MeshModuleLiquidTransfer::setmask()
{
    int mask = 0;
    mask |= PRE_FORCE;
    mask |= FINAL_INTEGRATE;
    return mask;
}

/* ---------------------------------------------------------------------- */

void MeshModuleLiquidTransfer::setup_pre_force(int vflag)
{
    liquid_flux_->setAll(0.);

}

/* ---------------------------------------------------------------------- */

void MeshModuleLiquidTransfer::setup(int vflag)
{ }

/* ---------------------------------------------------------------------- */

void MeshModuleLiquidTransfer::pre_force(int vflag)
{
    liquid_flux_->setAll(0.);
}

/* ---------------------------------------------------------------------- */

void MeshModuleLiquidTransfer::final_integrate()
{
    
    update_mesh_liquid_content();
}

/* ----------------------------------------------------------------------
   called during wall force calc
------------------------------------------------------------------------- */

void MeshModuleLiquidTransfer::add_liquid_flux(const int iTri, const double liquidTransferred, const bool limit, const double maxLiquid)
{
    limit_liquid_content_ = limit; // TODO get from registry
    max_liquid_content_ = maxLiquid;
    liquid_flux(iTri) += liquidTransferred;
}

/* ----------------------------------------------------------------------
   called during post force
------------------------------------------------------------------------- */

void MeshModuleLiquidTransfer::add_source_contribution(const int iTri, const double liquidTransferred)
{
    // liquidTransferred input has units vol/s => need to transfer to vol%/s
    liquid_flux(iTri) += liquidTransferred/(mesh->areaElem(iTri)*wall_thickness_);
}

/* ----------------------------------------------------------------------
   update liquid content of mesh elements
------------------------------------------------------------------------- */

void MeshModuleLiquidTransfer::update_mesh_liquid_content()
{
    int nTri = mesh->sizeLocal();
    double dt = update->dt;

    for(int iTri = 0; iTri < nTri; iTri++)
    {
        double lc = liquid_content(iTri);
        lc += liquid_flux(iTri) * dt;
        lc = fmax(0.0, fmin(1.0, lc));
        if (limit_liquid_content_)
            lc = fmin(lc, max_liquid_content_);
        liquid_content(iTri) = lc;
    }
}
