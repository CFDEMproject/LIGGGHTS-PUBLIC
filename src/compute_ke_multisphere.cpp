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
    (if no contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2014-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#include "mpi.h"
#include "string.h"
#include "compute_ke_multisphere.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "fix_multisphere.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeKEMultisphere::ComputeKEMultisphere(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  fix_ms_(0)
{
  if (narg != 3) error->compute_error(FLERR,this,"");

  scalar_flag = 1;
  extscalar = 1;
}

/* ---------------------------------------------------------------------- */

ComputeKEMultisphere::~ComputeKEMultisphere()
{
}

/* ---------------------------------------------------------------------- */

void ComputeKEMultisphere::init()
{
  fix_ms_ =  static_cast<FixMultisphere*>(modify->find_fix_style("multisphere",0));
  if (!fix_ms_)
    error->compute_error(FLERR,this,"fix multisphere does not exist");
}

/* ---------------------------------------------------------------------- */

double ComputeKEMultisphere::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  scalar = fix_ms_->extract_ke();
  
  return scalar;
}
