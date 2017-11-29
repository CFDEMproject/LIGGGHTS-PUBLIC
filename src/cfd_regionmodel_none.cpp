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
#include "error.h"
#include "memory.h"
#include "modify.h"
#include <cmath>
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
{/*
   nout = nlocal_last;
   int nlocal = atom->nlocal;

   outRegion->grow(nout,0);

   inregion = inRegion->vector_atom;
   outregion = outRegion->values;

   for(int i = 0; i < nout; i++)
      outregion[i] = 1.;

   for(int i = 0; i < nlocal; i++)
      inregion[i] = 1.;

   nlocal_last = nlocal;*/
}
