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

    Copyright 2015-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#include <string.h>
#include <stdlib.h>
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "neighbor.h"
#include "memory.h"
#include "modify.h"
#include "group.h"
#include "comm.h"
#include <cmath>
#include "vector_liggghts.h"
#include "fix_cfd_coupling_convection_species.h"
#include "fix_property_atom.h"
#include "fix_property_global.h"
#include "properties.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCfdCouplingConvectionSpecies::FixCfdCouplingConvectionSpecies(LAMMPS *lmp, int narg, char **arg) :  Fix(lmp, narg, arg)
{
    fix_coupling = NULL;
    fix_speciesConcentration = fix_convectiveFlux = fix_totalFlux = NULL;
    fix_speciesFluid = fix_speciesTransCoeff = NULL;

    int iarg = 3;

    if(narg < iarg + 4) error->all(FLERR,"Fix couple/cfd/speciesConvection: Wrong number of arguments");
    if(strcmp(arg[iarg++],"speciesName") != 0) error->all(FLERR,"Fix couple/cfd/speciesConvection: Expecting keyword 'speciesName'");
        strcpy(speciesName_,arg[iarg++]);
        sprintf(sourceName_,        "%sSource",speciesName_);
        sprintf(convectiveFluxName_,"%sFlux",speciesName_);
        sprintf(capacityName_,      "%sCapacity",speciesName_);
        sprintf(steName_,           "%sSTE",speciesName_);
        sprintf(totalFluxName_,     "%sTotalFlux",speciesName_);
        sprintf(speciesFluidName_,      "%sFluid",speciesName_);
        sprintf(speciesTransCoeffName_, "%sTransCoeff",speciesName_);

    if(strcmp(arg[iarg++],"species0") != 0) error->all(FLERR,"Fix couple/cfd/speciesConvection: Expecting keyword 'species0'");
        species0 = atof(arg[iarg++]);

    if(species0 < 0.) error->all(FLERR,"Fix couple/cfd/speciesConvection: species0 must be >= 0");
}

/* ---------------------------------------------------------------------- */
FixCfdCouplingConvectionSpecies::~FixCfdCouplingConvectionSpecies()
{
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingConvectionSpecies::pre_delete(bool unfixflag)
{
    if(fix_convectiveFlux) modify->delete_fix(convectiveFluxName_);
}

/* ---------------------------------------------------------------------- */

int FixCfdCouplingConvectionSpecies::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingConvectionSpecies::post_create()
{

  //  register species concentration
  if(!fix_speciesConcentration)
  {
        const char* fixarg[9];
        fixarg[0]=speciesName_;
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]=speciesName_;
        fixarg[4]="scalar";
        fixarg[5]="no";
        fixarg[6]="yes";
        fixarg[7]="no";
        fixarg[8]="0.";
        fix_speciesConcentration = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  //  register convective flux
  if(!fix_convectiveFlux)
  {
        const char* fixarg[9];
        fixarg[0]=convectiveFluxName_;
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]=convectiveFluxName_;
        fixarg[4]="scalar";
        fixarg[5]="no";
        fixarg[6]="yes";
        fixarg[7]="no";
        fixarg[8]="0.";
        fix_convectiveFlux = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  //  add species transfer model if not yet active
  if(!fix_speciesFluid)
  {
        const char* fixarg[9];
        fixarg[0]=speciesFluidName_;
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]=speciesFluidName_;
        fixarg[4]="scalar";
        fixarg[5]="no";
        fixarg[6]="yes";
        fixarg[7]="no";
        fixarg[8]="0.";
        fix_speciesFluid = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  //  register transfer coefficient
  if(!fix_speciesTransCoeff)
  {
        const char* fixarg[9];
        fixarg[0]=speciesTransCoeffName_;
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]=speciesTransCoeffName_;
        fixarg[4]="scalar";
        fixarg[5]="no";
        fixarg[6]="yes";
        fixarg[7]="no";
        fixarg[8]="0.";
        fix_speciesTransCoeff = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }
  FixScalarTransportEquation *fix_ste = modify->find_fix_scalar_transport_equation("speciesTransfer");
  if(!fix_ste)
  {
        int max_type = atom->get_properties()->max_type();
        Fix* capacityQty = modify->find_fix_property(capacityName_,"property/global","peratomtype",max_type,0,style,false);

        if(capacityQty==NULL)
        {
            sprintf(capacityName_, "none");
            if(comm->me==0)
                  printf("WARNING: FixCfdCouplingConvectionSpecies cannot locate capacity quantity. Thus, will assume you are not using capacity in scalar transport equation is '%s'.", capacityName_);

        }

        const char *newarg[15];
        newarg[0] = steName_;
        newarg[1] = group->names[igroup];
        newarg[2] = "transportequation/scalar";
        newarg[3] = "equation_id";
        newarg[4] = steName_;
        newarg[5] = "quantity";
        newarg[6] = speciesName_;
        newarg[7] = "default_value";
        char arg8[30];
        sprintf(arg8,"%f",species0);
        newarg[8] = arg8;
        newarg[9]  = "flux_quantity";
        newarg[10] = totalFluxName_;
        newarg[11] = "source_quantity";
        newarg[12] = sourceName_;
        newarg[13] = "capacity_quantity";
        newarg[14] = capacityName_;
        modify->add_fix(15,const_cast<char**>(newarg));
  }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingConvectionSpecies::init()
{
    // make sure there is only one fix of this style
    //if(modify->n_fixes_style(style) != 1)
    //  error->fix_error(FLERR,this,"More than one fix of this style is not allowed");

    // find coupling fix
    fix_coupling = static_cast<FixCfdCoupling*>(modify->find_fix_style_strict("couple/cfd",0));
    if(!fix_coupling)
      error->fix_error(FLERR,this,"needs a fix of type couple/cfd");

    //values to send to OF
    fix_coupling->add_push_property(speciesName_,"scalar-atom");

    //values to come from OF
    fix_coupling->add_pull_property(convectiveFluxName_,"scalar-atom");
    fix_coupling->add_pull_property(sourceName_,"scalar-atom");

    fix_coupling->add_pull_property(speciesFluidName_,"scalar-atom");
    fix_coupling->add_pull_property(speciesTransCoeffName_,"scalar-atom");

    fix_speciesConcentration  = static_cast<FixPropertyAtom*>(modify->find_fix_property(speciesName_,"property/atom","scalar",0,0,style));
    fix_convectiveFlux        = static_cast<FixPropertyAtom*>(modify->find_fix_property(convectiveFluxName_,"property/atom","scalar",0,0,style));
    fix_totalFlux             = static_cast<FixPropertyAtom*>(modify->find_fix_property(totalFluxName_,"property/atom","scalar",0,0,style));
    fix_speciesFluid          = static_cast<FixPropertyAtom*>(modify->find_fix_property(speciesFluidName_,"property/atom","scalar",0,0,style));
    fix_speciesTransCoeff     = static_cast<FixPropertyAtom*>(modify->find_fix_property(speciesTransCoeffName_,"property/atom","scalar",0,0,style));

    if(!fix_speciesConcentration || !fix_convectiveFlux)
      error->fix_error(FLERR,this,"could not find concentration and flux fix");

}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingConvectionSpecies::post_force(int)
{

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // communicate convective flux to ghosts, there might be new data
  if(0 == neighbor->ago)
        fix_convectiveFlux->do_forward_comm();

  double *totalFlux      = fix_totalFlux->vector_atom;
  double *convectiveFlux = fix_convectiveFlux->vector_atom;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      totalFlux[i] += convectiveFlux[i];
}
