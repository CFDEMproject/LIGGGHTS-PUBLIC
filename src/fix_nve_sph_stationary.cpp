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

    Andreas Aigner, JKU Linz
------------------------------------------------------------------------- */

#include "stdio.h"
#include "string.h"
#include "fix_nve_sph_stationary.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVESphStationary::FixNVESphStationary(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if ((atom->e_flag != 1) || (atom->rho_flag != 1))
    error->all(FLERR,
        "fix sph command requires atom_style with both energy and density");

  if (narg < 3)
    error->all(FLERR,"Illegal fix nve command");

  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixNVESphStationary::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVESphStationary::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixNVESphStationary::initial_integrate(int vflag)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *rho = atom->rho;
  double *drho = atom->drho;
  double *e = atom->e;
  double *de = atom->de;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      // same like src/USER-SPH/fix_meso !!
      e[i] += dtf * de[i]; // half-step update of particle internal energy
      //rho[i] += dtf * drho[i]; // ... and density
      rho[i] += 2.0 * dtf * drho[i]; // ... and density
    }
}

/* ---------------------------------------------------------------------- */

void FixNVESphStationary::final_integrate()
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *e = atom->e;
  double *de = atom->de;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      e[i] += dtf * de[i]; // half-step update of particle internal energy
    }
}

/* ---------------------------------------------------------------------- */

void FixNVESphStationary::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
