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

#include "stdio.h"
#include "string.h"
#include "fix_nve_sph.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVESph::FixNVESph(LAMMPS *lmp, int narg, char **arg) :
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

int FixNVESph::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVESph::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  if (strstr(update->integrate_style,"respa"))
    step_respa = ((Respa *) update->integrate)->step;
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixNVESph::initial_integrate(int vflag)
{
  double dtfm;

  // update v and x of atoms in group

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double **vest = atom->vest;
  double *rho = atom->rho;
  double *drho = atom->drho;
  double *e = atom->e;
  double *de = atom->de;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];

        // same like src/USER-SPH/fix_meso !!
        e[i] += dtf * de[i]; // half-step update of particle internal energy
        //rho[i] += dtf * drho[i]; // ... and density
        rho[i] += 2.0 * dtf * drho[i]; // ... and density

        // euler integration for velocity
        vest[i][0] = v[i][0] + 2.0 * dtfm * f[i][0];
        vest[i][1] = v[i][1] + 2.0 * dtfm * f[i][1];
        vest[i][2] = v[i][2] + 2.0 * dtfm * f[i][2];

        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];

        e[i] += dtf * de[i]; // half-step update of particle internal energy
        //rho[i] += dtf * drho[i]; // ... and density
        rho[i] += 2.0 * dtf * drho[i]; // ... and density

        // euler integration for velocity
        vest[i][0] = v[i][0] + 2.0 * dtfm * f[i][0];
        vest[i][1] = v[i][1] + 2.0 * dtfm * f[i][1];
        vest[i][2] = v[i][2] + 2.0 * dtfm * f[i][2];

        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVESph::final_integrate()
{
  double dtfm;

  // update v of atoms in group

  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  //double *rho = atom->rho; 
  //double *drho = atom->drho;
  double *e = atom->e;
  double *de = atom->de;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];

        e[i] += dtf * de[i]; // half-step update of particle internal energy
        //rho[i] += dtf * drho[i]; // ... and density
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];

        e[i] += dtf * de[i]; // half-step update of particle internal energy
        //rho[i] += dtf * drho[i]; // ... and density
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVESph::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;

  // innermost level - NVE update of v and x
  // all other levels - NVE update of v

  if (ilevel == 0) initial_integrate(vflag);
  else final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVESph::final_integrate_respa(int ilevel, int iloop)
{
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVESph::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
