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

#include "fix_diam_max.h"
#include "modify.h"
#include "math_extra_liggghts.h"
#include "fix_particledistribution_discrete.h"
#include <cmath>
#include <algorithm>

using namespace LAMMPS_NS;
using namespace MathExtraLiggghts;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixDiamMax::FixDiamMax(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  maxrbound_(0.)
{
  scalar_flag = 1;

  FixParticledistributionDiscrete *fpdd;
  int nfix = modify->n_fixes_style("particledistribution/discrete");
  for(int i = 0; i < nfix; i++)
  {
      fpdd = static_cast<FixParticledistributionDiscrete*>(modify->find_fix_style("particledistribution/discrete",0));
      maxrbound_ = std::max(maxrbound_,fpdd->max_r_bound());
  }
}

/* ---------------------------------------------------------------------- */

FixDiamMax::~FixDiamMax()
{

}

/* ----------------------------------------------------------------------*/

int FixDiamMax::setmask()
{
    int mask = 0;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixDiamMax::init()
{
  maxrbound_ = 0.;

  FixParticledistributionDiscrete *fpdd;
  int nfix = modify->n_fixes_style("particledistribution/discrete");
  for(int i = 0; i < nfix; i++)
  {
      fpdd = static_cast<FixParticledistributionDiscrete*>(modify->find_fix_style("particledistribution/discrete",0));
      maxrbound_ = std::max(maxrbound_,fpdd->max_r_bound());
  }
}

/* ---------------------------------------------------------------------- */

double FixDiamMax::compute_scalar()
{
    return 2.*maxrbound_;
}
