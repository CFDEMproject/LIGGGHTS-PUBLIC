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

    Copyright 2014-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#include <cmath>
#include <stdlib.h>
#include <string.h>
#include "fix_buoyancy.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "respa.h"
#include "domain.h"
#include "modify.h"
#include "fix_gravity.h"
#include "region.h"
#include "vector_liggghts.h"
#include "math_extra_liggghts.h"
#include "mpi_liggghts.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBuoyancy::FixBuoyancy(LAMMPS *lmp, int narg, char **arg) :
  FixBaseLiggghts(lmp, narg, arg),
  density_(0.),
  dim_(-1),
  direction_(1.),
  fluid_level_(0.),
  fix_gravity_(0),
  force_flag_(false)
{
  vectorZeroize3D(buyoancy_force_total_);

  do_support_respa();

  if (narg < 3) error->fix_error(FLERR,this,"not enough arguments");

  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 1;

  // args

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"density") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments for 'density'");
      density_ = atof(arg[iarg+1]);
      if (density_ <= 0.)
        error->fix_error(FLERR,this,"'density' > 0 required");
      iarg += 2;
    } else if (strcmp(arg[iarg],"level") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments for 'level'");
      fluid_level_ = atof(arg[iarg+1]);
      if (fluid_level_ <= 0.)
        error->fix_error(FLERR,this,"'level' > 0 required");
      iarg += 2;
    } else if (strcmp(arg[iarg],"dim") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments for 'dim'");
      if(0 == strcmp(arg[iarg+1],"x"))
        dim_ = 0;
      else if(0 == strcmp(arg[iarg+1],"y"))
        dim_ = 1;
      else if(0 == strcmp(arg[iarg+1],"z"))
        dim_ = 2;
      else
        error->fix_error(FLERR,this,"expecting 'x' or 'y' or 'z' after 'dim'");
      iarg += 2;
    } else if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments for 'region'");
      process_region(arg[iarg+1]);
      iarg += 2;
    } else error->fix_error(FLERR,this," expecting 'density' or 'region'");
  }

  if(-1 == dim_)
    error->fix_error(FLERR,this," you have to specify 'dim'");

  if(MathExtraLiggghts::compDouble(density_,0.,1e-6))
    error->fix_error(FLERR,this," you have to specify 'density'");

}

/* ---------------------------------------------------------------------- */

FixBuoyancy::~FixBuoyancy()
{
}

/* ---------------------------------------------------------------------- */

int FixBuoyancy::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBuoyancy::init()
{
  FixBaseLiggghts::init();

  if(!atom->density_flag)
    error->fix_error(FLERR,this,"requires per-atom density");

  if(1 != modify->n_fixes_style_strict("gravity"))
      error->fix_error(FLERR,this,"need exactly one fix gravity");

  fix_gravity_ = static_cast<FixGravity*>(modify->find_fix_style_strict("gravity",0));

  // test if dim of fix gravity and this fix coincide
  test_direction();
}

/* ---------------------------------------------------------------------- */

void FixBuoyancy::test_direction()
{
  double gravity[3];
  fix_gravity_->get_gravity(gravity);

  // test if dim of fix gravity and this fix coincide
  double test[3] = {1.,1.,1.};
  test[dim_] = 0.;
  if(MathExtraLiggghts::abs(vectorDot3D(gravity,test)) > 0.)
    error->fix_error(FLERR,this," fix gravity direction is not equal to 'dim'");

  vectorZeroize3D(test);
  test[dim_] = 1.;

  direction_ = vectorDot3D(gravity,test);
  direction_ = direction_ / MathExtraLiggghts::abs(direction_);

}

/* ---------------------------------------------------------------------- */

void FixBuoyancy::setup(int vflag)
{
  FixBaseLiggghts::setup(vflag);

  // error checks on coarsegraining
  if(force->cg_active())
    error->cg(FLERR,this->style);
}

/* ---------------------------------------------------------------------- */

void FixBuoyancy::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixBuoyancy::post_force(int vflag)
{
  // apply drag force to atoms in group
  // direction is opposed to velocity vector
  // magnitude depends on atom type

  double **x = atom->x;
  double **f = atom->f;
  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double buyoancy_force;
  vectorZeroize3D(buyoancy_force_total_);
  force_flag_ = false;

  double gravity[3],gravity_mag;
  fix_gravity_->get_gravity(gravity);
  gravity_mag = vectorMag3D(gravity);

  // test if dim of fix gravity and this fix coincide
  // need to re-do this since orientation of gravity might change
  test_direction();

  double c = 4.*M_PI/3.;
  double pi6 = M_PI/6.;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      if (region_ && !region_->match(x[i][0],x[i][1],x[i][2]))
          continue;

      double h = (x[i][dim_] - fluid_level_)*direction_;
      double hmag = MathExtraLiggghts::abs(h);
      double Vimmersed = 0.;

      // particle not immersed
      if(hmag > radius[i] && h < 0)
      {
          
          continue;
      }
      // particle fully immersed
      else if(hmag > radius[i] && h > 0)
      {
          
          Vimmersed = c*radius[i]*radius[i]*radius[i];
      }
      // larger part of sphere immersed
      else if(h > 0)
      {
          
          h = radius[i] - h;
          double a = sqrt(2.*radius[i]*h - h*h);
          double Vsegment = h*pi6 * (3.*a*a + h*h);
          Vimmersed = (c*radius[i]*radius[i]*radius[i] - Vsegment);
      }
      // smaller part of sphere immersed
      else
      {
          
          h = -h;
          h = radius[i] - h;
          double a = sqrt(2.*radius[i]*h - h*h);
          Vimmersed = h * pi6 * (3.*a*a + h*h);
      }

      buyoancy_force = - direction_*gravity_mag*Vimmersed*(density_);
      f[i][dim_] += buyoancy_force;
      buyoancy_force_total_[dim_] += buyoancy_force;
      
    }
}

/* ---------------------------------------------------------------------- */

void FixBuoyancy::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixBuoyancy::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixBuoyancy::compute_vector(int n)
{
  // only sum across procs one time

  if (!force_flag_) {
    MPI_Sum_Vector(buyoancy_force_total_,3,world);
    force_flag_ = true;
  }
  return buyoancy_force_total_[n];
}
