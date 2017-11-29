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

    Stefan Radl (TU Graz)

    Copyright 2014-     TU Graz
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   Description
    Fix for to deform a box, but ensures that it is initialized AFTER
    the fix 'contacthistory' (if any). This is to avoid troubles with
    contacthistory when flipping the box.
------------------------------------------------------------------------- */

#include <string.h>
#include <stdlib.h>
#include <cmath>
#include "fix_deform_check.h"
#include "atom.h"
#include "update.h"
#include "comm.h"
#include "irregular.h"
#include "domain.h"
#include "lattice.h"
#include "force.h"
#include "modify.h"
#include "math_const.h"
#include "kspace.h"
#include "input.h"
#include "variable.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */
FixDeformCheck::FixDeformCheck(LAMMPS *lmp, int narg, char **arg)
:
 FixDeform(lmp, narg, arg),
 iContactHist_(-1),
 noProblem_(false)
{
  //Check if contacthistory is present before building this one
  iContactHist_ = modify->find_fix("contacthistory");

}

/* ---------------------------------------------------------------------- */

FixDeformCheck::~FixDeformCheck()
{
}

/* ---------------------------------------------------------------------- */

void FixDeformCheck::init()
{
  FixDeform::init();

  if (iContactHist_ >= 0)
    noProblem_ = true;  //existed before this fix was build, so okey
  else
  {
     int iF = modify->find_fix("contacthistory");
     if(iF<0)
        noProblem_ = true;  //fix does not exist at all, so okey
  }

  if(!noProblem_)
     error->one(FLERR,"Cannot find fix contacthistory. This fix must be present to run this fix.\n"
                      "Placing 'run 0' before this fix might help.  \n\n");

}

