/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
   the producer of the LIGGGHTS® software and the CFDEM®coupling software
   See http://www.cfdem.com/terms-trademark-policy for details.

   LIGGGHTS® is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include "string.h"
#include "stdlib.h"
#include "atom.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "math.h"
#include "fix_property_atom.h"
#include "fix_property_global.h"
#include "fix_cfd_coupling.h"
#include "cfd_regionmodel_none.h"

#define DELTA 10000

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

CfdRegionmodelNone::CfdRegionmodelNone(LAMMPS *lmp, int jarg,int narg, char **arg,FixCfdCoupling *fc)  :
  CfdRegionmodel(lmp, jarg, narg, arg,fc)
{
    iarg = jarg;
    //do something else here

    inRegion = NULL;
    outRegion = NULL;
    inregion = NULL;
    outregion = NULL;
    nout = 0;
    nlocal_last = 0;
}

/* ---------------------------------------------------------------------- */

CfdRegionmodelNone::~CfdRegionmodelNone()
{
    if(inRegion) modify->delete_fix("inRegion");
    if(outRegion) modify->delete_fix("outRegion");
}

/* ---------------------------------------------------------------------- */

void CfdRegionmodelNone::init()
{
    const char* fixarg[9];
    if(!inRegion)
    {
        fixarg[0]="inRegion";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="inRegion";
        fixarg[4]="scalar"; 
        fixarg[5]="yes";    
        fixarg[6]="no";    
        fixarg[7]="no";    
        fixarg[8]="1.";
        modify->add_fix(9,const_cast<char**>(fixarg));
    }

    inRegion = static_cast<FixPropertyAtom*>(modify->find_fix_property("inRegion","property/atom","scalar",1,0,"cfd_regionmodel none"));

    if(!outRegion)
    {
        fixarg[0]="outRegion";
        fixarg[1]="all";
        fixarg[2]="property/global";
        fixarg[3]="outRegion";
        fixarg[4]="vector";
        fixarg[5]="0.";
        fixarg[6]="0.";
        modify->add_fix(7,const_cast<char**>(fixarg));
    }

    outRegion = static_cast<FixPropertyGlobal*>(modify->find_fix_property("outRegion","property/global","vector",2,0,"cfd_regionmodel none"));

    special_settings();
}

/* ---------------------------------------------------------------------- */

void CfdRegionmodelNone::special_settings()
{
  //values to be transfered to OF
  add_push_property("inRegion","scalar");
  add_push_property("outRegion","globalvector");

  //values to come from OF
  //add_pull_property("inRegion","scalar");
}

/* ---------------------------------------------------------------------- */

void CfdRegionmodelNone::rm_update()
{
   nout = nlocal_last;
   int nlocal = atom->nlocal;

   outRegion->grow(nout,0);

   inregion = inRegion->vector_atom;
   outregion = outRegion->values;

   for(int i = 0; i < nout; i++)
      outregion[i] = 1.;

   for(int i = 0; i < nlocal; i++)
      inregion[i] = 1.;

   nlocal_last = nlocal;
}
