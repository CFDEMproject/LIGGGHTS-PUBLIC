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

#include "fix_heat_gran.h"

#include "atom.h"
#include "fix_property_atom.h"
#include "fix_scalar_transport_equation.h"
#include "force.h"
#include "group.h"
#include "math_extra.h"
#include "modify.h"
#include "pair_gran.h"
#include "stdlib.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixHeatGran::FixHeatGran(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg){

  if ((!atom->radius_flag)||(!atom->rmass_flag)) error->all(FLERR,"Fix heat/gran needs per particle radius and mass");

  if (narg < 5)
    error->fix_error(FLERR,this,"not enough arguments");

  int iarg = 3;

  if(strcmp(arg[iarg++],"initial_temperature"))
    error->fix_error(FLERR,this,"expecting keyword 'initial_temperature'");
  T0 = atof(arg[iarg++]);

  fix_temp = fix_heatFlux = fix_heatSource = NULL;
  fix_ste = NULL;
  fix_directionalHeatFlux = NULL;
  peratom_flag = 1;      
  size_peratom_cols = 0; 
  peratom_freq = 1;

  scalar_flag = 1; 
  global_freq = 1; 

  cpl = NULL;
}

/* ---------------------------------------------------------------------- */

void FixHeatGran::post_create()
{
  // register directional flux
  fix_directionalHeatFlux = static_cast<FixPropertyAtom*>(modify->find_fix_property("directionalHeatFlux","property/atom","vector",3,0,this->style,false));
  if(!fix_directionalHeatFlux)
  {
    const char* fixarg[11];
    fixarg[0]="directionalHeatFlux";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="directionalHeatFlux";
    fixarg[4]="vector";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.";
    fixarg[9]="0.";
    fixarg[10]="0.";
    fix_directionalHeatFlux = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
  }

  fix_ste = modify->find_fix_scalar_transport_equation("heattransfer");
  if(!fix_ste)
  {
    const char * newarg[15];
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
    newarg[8] = arg8;
    newarg[9] = "flux_quantity";
    newarg[10] = "heatFlux";
    newarg[11] = "source_quantity";
    newarg[12] = "heatSource";
    newarg[13] = "capacity_quantity";
    newarg[14] = "thermalCapacity";
    modify->add_fix(15,(char**)newarg);
  }
}

/* ---------------------------------------------------------------------- */

void FixHeatGran::updatePtrs(){

  Temp = fix_temp->vector_atom;
  vector_atom = Temp; 

  heatFlux = fix_heatFlux->vector_atom;
  heatSource = fix_heatSource->vector_atom;
  directionalHeatFlux = fix_directionalHeatFlux->array_atom;
}

/* ---------------------------------------------------------------------- */

void FixHeatGran::init()
{
  
  int n_ms = modify->n_fixes_style("multisphere");
  if(n_ms > 0)
    error->fix_error(FLERR,this,"may not be used together with fix multisphere");

  if (!atom->radius_flag || !atom->rmass_flag)
    error->fix_error(FLERR,this,"must use a granular atom style ");

    // check if a fix of this style already exists
  if(modify->n_fixes_style(style) > 1)
    error->fix_error(FLERR,this,"cannot have more than one fix of this style");

  if(!force->pair_match("gran", 0))
    error->fix_error(FLERR,this,"needs a granular pair style to be used");

  pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));
  history_flag = pair_gran->is_history();

  fix_ste = modify->find_fix_scalar_transport_equation("heattransfer");
  if(!fix_ste) error->fix_error(FLERR,this,"needs a fix transportequation/scalar to work with");

  fix_temp = static_cast<FixPropertyAtom*>(modify->find_fix_property("Temp","property/atom","scalar",0,0,style));
  fix_heatFlux = static_cast<FixPropertyAtom*>(modify->find_fix_property("heatFlux","property/atom","scalar",0,0,style));
  fix_heatSource = static_cast<FixPropertyAtom*>(modify->find_fix_property("heatSource","property/atom","scalar",0,0,style));
  fix_directionalHeatFlux = static_cast<FixPropertyAtom*>(modify->find_fix_property("directionalHeatFlux","property/atom","vector",0,0,style));

  updatePtrs();
}

/* ---------------------------------------------------------------------- */

int FixHeatGran::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixHeatGran::initial_integrate(int vflag)
{
  updatePtrs();

  //reset heat flux
  //sources are not reset
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
  {
     if (mask[i] & groupbit)
     {
        directionalHeatFlux[i][0] = 0.;
        directionalHeatFlux[i][1] = 0.;
        directionalHeatFlux[i][2] = 0.;
     }
  }

  //update ghosts
  fix_directionalHeatFlux->do_forward_comm();
}

/* ---------------------------------------------------------------------- */

double FixHeatGran::compute_scalar()
{
    return fix_ste->compute_scalar();
}

/* ---------------------------------------------------------------------- */

void FixHeatGran::cpl_evaluate(class ComputePairGranLocal * cpl){

  char *mystyle = style;
  char *emsg = new char[100];
  sprintf(emsg, "Fix %s does not implement cpl_evaluate().\n", mystyle);
  error->all(FLERR, emsg);

}

/* ---------------------------------------------------------------------- */

void FixHeatGran::register_compute_pair_local(class ComputePairGranLocal *ptr){

  char *mystyle = style;
  char *emsg = new char[100];
  sprintf(emsg, "Fix %s does not implement register_compute_pair_local().\n", mystyle);
  error->all(FLERR, emsg);

}

/* ---------------------------------------------------------------------- */

void FixHeatGran::unregister_compute_pair_local(class ComputePairGranLocal *ptr){

  char *mystyle = style;
  char *emsg = new char[100];
  sprintf(emsg, "Fix %s does not implement unregister_compute_pair_local().\n", mystyle);
  error->all(FLERR, emsg);

}
