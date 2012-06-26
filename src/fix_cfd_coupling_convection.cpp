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

#include "string.h"
#include "stdlib.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "group.h"
#include "comm.h"
#include "math.h"
#include "vector_liggghts.h"
#include "fix_cfd_coupling_convection.h"
#include "fix_property_atom.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCfdCouplingConvection::FixCfdCouplingConvection(LAMMPS *lmp, int narg, char **arg) :  Fix(lmp, narg, arg)
{
    fix_coupling = NULL;
    fix_convectiveFlux  = fix_heatFlux = NULL;

    int iarg = 3;

    if(narg < iarg + 2) error->all(FLERR,"Fix couple/cfd/convection: Wrong number of arguments");
    if(strcmp(arg[iarg++],"T0") != 0) error->all(FLERR,"Fix couple/cfd/convection: Expecting keyword 'T0'");
    T0 = atof(arg[iarg++]);

    if(T0 < 0.) error->all(FLERR,"Fix couple/cfd/convection: T0 must be >= 0");
}

/* ---------------------------------------------------------------------- */

FixCfdCouplingConvection::~FixCfdCouplingConvection()
{

}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingConvection::pre_delete(bool unfixflag)
{
    if(fix_convectiveFlux) modify->delete_fix("convectiveHeatFlux");
}

/* ---------------------------------------------------------------------- */

int FixCfdCouplingConvection::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingConvection::post_create()
{
  //  register convective flux
  if(!fix_convectiveFlux)
  {
        char* fixarg[11];
        fixarg[0]="convectiveHeatFlux";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="convectiveHeatFlux";
        fixarg[4]="scalar"; 
        fixarg[5]="no";    
        fixarg[6]="yes";    
        fixarg[7]="no";    
        fixarg[8]="0.";
        fix_convectiveFlux = modify->add_fix_property_atom(9,fixarg,style);
  }

  //  add heat transfer model if not yet active
  FixScalarTransportEquation *fix_ste = modify->find_fix_scalar_transport_equation("heattransfer");
  if(!fix_ste)
  {
        char **newarg = new char*[15];
        newarg[0] = (char *) "ste_heattransfer";
        newarg[1] = group->names[igroup];
        newarg[2] = (char *) "transportequation/scalar";
        newarg[3] = (char *) "equation_id";
        newarg[4] = (char *) "heattransfer";
        newarg[5] = (char *) "quantity";
        newarg[6] = (char *) "Temp";
        newarg[7] = (char *) "default_value";
        newarg[8] = new char[30];
        sprintf(newarg[8],"%f",T0);
        newarg[9] = (char *) "flux_quantity";
        newarg[10] = (char *) "heatFlux";
        newarg[11] = (char *) "source_quantity";
        newarg[12] = (char *) "heatSource";
        newarg[13] = (char *) "capacity_quantity";
        newarg[14] = (char *) "thermalCapacity";
        modify->add_fix(15,newarg);
        delete [] newarg[8];
        delete [] newarg;
  }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingConvection::init()
{
    // make sure there is only one fix of this style
    if(modify->n_fixes_style(style) != 1)
      error->all(FLERR,"More than one fix of style couple/cfd/dust/simple is not allowed");

    // find coupling fix
    fix_coupling = static_cast<FixCfdCoupling*>(modify->find_fix_style_strict("couple/cfd",0));
    if(!fix_coupling)
      error->all(FLERR,"Fix couple/cfd/convection needs a fix of type couple/cfd");

    //values to send to OF
    fix_coupling->add_push_property("Temp","scalar-atom");

    //values to come from OF
    fix_coupling->add_pull_property("convectiveHeatFlux","scalar-atom");

    // heat transfer added heatFlux, get reference to it
    fix_heatFlux = static_cast<FixPropertyAtom*>(modify->find_fix_property("heatFlux","property/atom","scalar",0,0,style));

    fix_convectiveFlux = static_cast<FixPropertyAtom*>(modify->find_fix_property("convectiveHeatFlux","property/atom","scalar",0,0,style));
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingConvection::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *heatFlux = fix_heatFlux->vector_atom;
  double *convectiveFlux = fix_convectiveFlux->vector_atom;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      heatFlux[i] += convectiveFlux[i];
}
