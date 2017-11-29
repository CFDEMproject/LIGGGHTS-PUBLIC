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

    Christoph Kloss (DCS Computing GmbH, Linz)
    Stefan Radl (TU Graz)

    Copyright 2015-     DCS Computing GmbH, Linz
    Copyright 2015-     TU Graz
------------------------------------------------------------------------- */

#include <string.h>
#include <stdlib.h>
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "group.h"
#include "comm.h"
#include <cmath>
#include "vector_liggghts.h"
#include "fix_cfd_coupling_convection_impl.h"
#include "fix_property_atom.h"
#include "neighbor.h"
#include "fix_scalar_transport_equation.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCfdCouplingConvectiveImpl::FixCfdCouplingConvectiveImpl(LAMMPS *lmp, int narg, char **arg) :  Fix(lmp, narg, arg)
{
    integrateHeatEqn_ = false; //default: DO NOT integrate heat equation, since we expect external program to do so
    forceExplicit_    = false; //default: do not force explicit calculation, thus do implicit calculation

    fix_coupling = NULL;
    fix_convectiveFlux  = fix_heatFluid  = fix_heatTransCoeff = NULL;
    fix_heatFlux = NULL;
    fix_temperature = NULL;
    T0 = 0;

    int iarg = 3;

    if(narg >= iarg + 1) 
    {
        if(strcmp(arg[iarg++],"integrateHeatEqn") != 0) 
            error->all(FLERR,"Fix couple/cfd/convectiveImpl: Expecting keyword 'integrateHeatEqn'");
        else
        {
            if(narg < iarg + 1) error->all(FLERR,"Fix couple/cfd/convectiveImpl: Wrong number of arguments after integrateHeatEqn. Provide 'true' or 'false' ");

            if(strcmp(arg[iarg],"false") == 0)
                integrateHeatEqn_ = false;
            else if (strcmp(arg[iarg],"true") == 0)
            {
                iarg++;
                integrateHeatEqn_ = true;
                if(strcmp(arg[iarg++],"T0") != 0) error->all(FLERR,"Fix couple/cfd/convectiveImpl: Expecting keyword 'T0'");
                if(narg < iarg+1) error->all(FLERR,"Fix couple/cfd/convectiveImpl: please specify a value after 'T0' ");
                T0 = atof(arg[iarg++]);
            }
            else
                error->all(FLERR,"Fix couple/cfd/convectiveImpl: Wrong argument after integrateHeatEqn. provide 'true' or 'false' ");
        }
        if(narg > iarg) 
        {
          if(strcmp(arg[iarg++],"forceExplicit") == 0) 
          {
            if(narg < iarg + 1) error->all(FLERR,"Fix couple/cfd/convectiveImpl: Wrong number of arguments after forceExplicit. provide 'true' or 'false' ");

            if(strcmp(arg[iarg],"false") == 0)
                forceExplicit_ = false;
            else if (strcmp(arg[iarg],"true") == 0)
            {
                iarg++;
                forceExplicit_ = true;
            }
            else
                error->all(FLERR,"Fix couple/cfd/convectiveImpl: Wrong argument after forceExplicit. provide 'true' or 'false' ");
          }
        }
    }
}

/* ---------------------------------------------------------------------- */

FixCfdCouplingConvectiveImpl::~FixCfdCouplingConvectiveImpl()
{
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingConvectiveImpl::pre_delete(bool unfixflag)
{
    if(fix_heatFluid)       modify->delete_fix("fix_heatFluid");
    if(fix_heatTransCoeff)  modify->delete_fix("fix_heatTransCoeff");
    if(fix_convectiveFlux)  modify->delete_fix("convectiveHeatFlux");
}

/* ---------------------------------------------------------------------- */

int FixCfdCouplingConvectiveImpl::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingConvectiveImpl::post_create()
{
  //  register fluid temperature and transfer coefficient
  if(!fix_heatFluid)
  {
        const char* fixarg[11];
        fixarg[0]="heatFluid";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="heatFluid";
        fixarg[4]="scalar";
        fixarg[5]="no";
        fixarg[6]="yes";
        fixarg[7]="no";
        fixarg[8]="0.";
        fix_heatFluid = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  if(!fix_heatTransCoeff)
  {
        const char* fixarg[11];
        fixarg[0]="heatTransCoeff";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="heatTransCoeff";
        fixarg[4]="scalar";
        fixarg[5]="no";
        fixarg[6]="yes";
        fixarg[7]="no";
        fixarg[8]="0.";
        fix_heatTransCoeff = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  //  register convective flux
  if(!fix_convectiveFlux)
  {
        const char* fixarg[11];
        fixarg[0]="convectiveHeatFlux";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="convectiveHeatFlux";
        fixarg[4]="scalar";
        fixarg[5]="no";
        fixarg[6]="yes";
        fixarg[7]="no";
        fixarg[8]="0.";
        fix_convectiveFlux = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }
  if(!integrateHeatEqn_) return; //this is usefil in case external tool performs integration

  //  add heat transfer model if not yet active
  FixScalarTransportEquation *fix_ste = modify->find_fix_scalar_transport_equation("heattransfer");
  if(!fix_ste)
  {
        const char *newarg[15];
        newarg[0] = "ste_heattransfer";
        newarg[1] = group->names[igroup];
        newarg[2] = "transportequation/scalar";
        newarg[3] = "equation_id";
        newarg[4] = "heattransfer";
        newarg[5] = "quantity";
        newarg[6] = "Temp";
        newarg[7] = "default_value";
        char arg8[30];
        sprintf(arg8,"%f",T0);
        newarg[8]  = arg8;
        newarg[9]  = "flux_quantity";
        newarg[10] = "heatFlux";
        newarg[11] = "source_quantity";
        newarg[12] = "heatSource";
        newarg[13] = "capacity_quantity";
        newarg[14] = "thermalCapacity";
        modify->add_fix(15,const_cast<char**>(newarg));
  }
  fix_ste = modify->find_fix_scalar_transport_equation("heattransfer");

  if(!forceExplicit_)
    fix_ste->register_implicit_fixes((char *)"heatFluid", 0.0, (char *)"heatTransCoeff", 0);
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingConvectiveImpl::init()
{
    // make sure there is only one fix of this style
    if(modify->n_fixes_style(style) != 1)
      error->all(FLERR,"More than one fix of style couple/cfd/convectiveImpl is not allowed");

    // find coupling fix
    fix_coupling = static_cast<FixCfdCoupling*>(modify->find_fix_style_strict("couple/cfd",0));
    if(!fix_coupling)
      error->all(FLERR,"Fix couple/cfd/convectiveImpl needs a fix of type couple/cfd");

    //values to send to OF, this is the SURFACE TEMPERATURE!
    fix_coupling->add_push_property("Temp","scalar-atom");

    //values to come from OF
    fix_coupling->add_pull_property("heatFluid","scalar-atom");
    fix_coupling->add_pull_property("heatTransCoeff","scalar-atom");
    fix_coupling->add_pull_property("convectiveHeatFlux","scalar-atom");

    fix_heatFluid = static_cast<FixPropertyAtom*>(modify->find_fix_property("heatFluid","property/atom","scalar",0,0,style));

    fix_heatTransCoeff = static_cast<FixPropertyAtom*>(modify->find_fix_property("heatTransCoeff","property/atom","scalar",0,0,style));
    fix_convectiveFlux = static_cast<FixPropertyAtom*>(modify->find_fix_property("convectiveHeatFlux","property/atom","scalar",0,0,style));
    if(!integrateHeatEqn_) return; //only access heat flux if needed

    fix_heatFlux    = static_cast<FixPropertyAtom*>(modify->find_fix_property("heatFlux","property/atom","scalar",0,0,style));
    fix_temperature = static_cast<FixPropertyAtom*>(modify->find_fix_property("Temp","property/atom","scalar",0,0,style));
}

/* ---------------------------------------------------------------------- */
void FixCfdCouplingConvectiveImpl::post_force(int)
{
  // communicate convective flux to ghosts, there might be new data
  if(0 == neighbor->ago)
  {
        fix_heatFluid->do_forward_comm();
        fix_heatTransCoeff->do_forward_comm();
        fix_convectiveFlux->do_forward_comm();
  }

  if(!integrateHeatEqn_) return; //only integrate if needed
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *Temp           = fix_temperature->vector_atom;
  double *convectiveFlux = fix_convectiveFlux->vector_atom;

  double *heatFlux      = fix_heatFlux->vector_atom;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      heatFlux[i] += convectiveFlux[i];

  if(!forceExplicit_) return;  //only add implicit terms to heatFlux if necessary

  //Force explicit calculation by adding implicit part to (explicit) flux
  double *radius         = atom->radius;
  double *heatFluid      = fix_heatFluid->vector_atom;
  double *heatTransCoeff = fix_heatTransCoeff->vector_atom;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
    {
      double currRadius = radius[i];
      heatFlux[i] += heatTransCoeff[i]
                  *  currRadius * currRadius * 12.5663706144 
                  * (heatFluid[i] - Temp[i]); //12.5663706144=4*pi

    }
}
