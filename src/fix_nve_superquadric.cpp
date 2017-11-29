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
    Alexander Podlozhnyuk, DCS Computing GmbH, Linz

    Copyright 2015-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

//#ifdef SUPERQUADRIC_ACTIVE_FLAG
#ifdef SUPERQUADRIC_ACTIVE_FLAG
#include <cmath>
#include <stdio.h>
#include <string.h>
#include "fix_nve_superquadric.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "respa.h"
#include "force.h"
#include "error.h"
#include "domain.h"
#include "math_extra_liggghts_nonspherical.h"
#include "fix_property_atom.h"

/* ---------------------------------------------------------------------- */

FixNVESuperquadric::FixNVESuperquadric(LAMMPS *lmp, int narg, char **arg) :
  FixNVEAsphereBase(lmp, narg, arg)
{

}

/* ---------------------------------------------------------------------- */

void FixNVESuperquadric::init()
{
  FixNVEAsphereBase::init();

  // error checks

  if (!atom->superquadric_flag )
    error->all(FLERR,"Fix nve/superquadric requires atom style superquadric");
}
#endif
