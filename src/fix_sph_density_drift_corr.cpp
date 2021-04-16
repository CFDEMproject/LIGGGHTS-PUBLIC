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

    Andreas Aigner (JKU Linz)
    Andreas Eitzlmayr (TU Graz)

    Copyright 2009-2012 JKU Linz
    Copyright 2013-     TU Graz
------------------------------------------------------------------------- */

#include <cmath>
#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include "fix_sph_density_drift_corr.h"
#include "update.h"
#include "respa.h"
#include "atom.h"
#include "force.h"
#include "modify.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "sph_kernels.h"
#include "fix_property_atom.h"
#include "timer.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSphDensityDriftCorr::FixSphDensityDriftCorr(LAMMPS *lmp, int narg, char **arg) :
  FixSph(lmp, narg, arg)
{
  int iarg = 0;
  if (narg < iarg+9) error->fix_error(FLERR,this,"Illegal fix sph/density/drift/corr command, not enough arguments");

  // parameters
  iarg += 3;

  if (strcmp(arg[iarg],"density") == 0) {
    rhoaveDef = force->numeric(FLERR,arg[iarg+1]);
    iarg += 2;
    if (rhoaveDef <= 0) error->fix_error(FLERR,this,"density <= 0 not allowed");
  } else error->fix_error(FLERR,this,"");

  if (strcmp(arg[iarg],"coeff") == 0) {
    coeff = force->numeric(FLERR,arg[iarg+1]);
    iarg += 2;
    if (coeff <= 0) error->fix_error(FLERR,this,"coeff <= 0 not allowed");
  } else error->fix_error(FLERR,this,"");

  if (strcmp(arg[iarg],"after") == 0) {
    after = force->inumeric(FLERR,arg[iarg+1]);
    iarg += 2;
    if (after < 0) error->fix_error(FLERR,this,"after < 0 not allowed");
  } else error->fix_error(FLERR,this,"");

  while (iarg < narg) {
    // kernel style
    if (strcmp(arg[iarg],"sphkernel") == 0) {
          if (iarg+2 > narg) error->fix_error(FLERR,this,"Illegal fix sph/density/continuity command");

          if(kernel_style) delete []kernel_style;
          kernel_style = new char[strlen(arg[iarg+1])+1];
          strcpy(kernel_style,arg[iarg+1]);

          // check uniqueness of kernel IDs

          int flag = SPH_KERNEL_NS::sph_kernels_unique_id();
          if(flag < 0) error->fix_error(FLERR,this,"Cannot proceed, sph kernels need unique IDs, check all sph_kernel_* files");

          // get kernel id

          kernel_id = SPH_KERNEL_NS::sph_kernel_id(kernel_style);
          if(kernel_id < 0) error->fix_error(FLERR,this,"Illegal fix sph/density/continuity command, unknown sph kernel");

          iarg += 2;

    } else error->fix_error(FLERR,this,"Illegal fix sph/density/continuity command");
  }

  step = 0;
  active = 0;
}

/* ---------------------------------------------------------------------- */

FixSphDensityDriftCorr::~FixSphDensityDriftCorr()
{

}

/* ---------------------------------------------------------------------- */

int FixSphDensityDriftCorr::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSphDensityDriftCorr::init()
{
  FixSph::init();
  //dt = update->dt;
}

/* ---------------------------------------------------------------------- */

void FixSphDensityDriftCorr::pre_force(int)
{
  //template function for using per atom or per atomtype smoothing length
  if (mass_type) pre_force_eval<1>();
  else pre_force_eval<0>();
}

/* ---------------------------------------------------------------------- */

template <int MASSFLAG>
void FixSphDensityDriftCorr::pre_force_eval()
{
  if (active == 0) {
    if (step < after) {
      step += 1;
    } else {
      active = 1;
    }
  }

  if (active == 1) {

    int i,ntotal;
    double rhoavelocal=0.0,rhoavetotal,drho_corr;

    int *mask = atom->mask;
    double *rho = atom->rho;
    double *drho = atom->drho;
    int nlocal = atom->nlocal;

    // sum up mass and density
    for (i = 0; i < nlocal; i++)
    {
      rhoavelocal += rho[i];
    }

    MPI_Allreduce(&nlocal,&ntotal,1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&rhoavelocal,&rhoavetotal,1,MPI_DOUBLE,MPI_SUM,world);

    // average density
    rhoavetotal /= ntotal;

    // density correction
    drho_corr = (rhoaveDef - rhoavetotal) * coeff;

    // correct drho
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        drho[i] += drho_corr;
      }
    }
  }
}
