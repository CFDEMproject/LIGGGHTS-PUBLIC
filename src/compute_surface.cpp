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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "compute_surface.h"
#include "compute_com.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "modify.h"
#include "neighbor.h"
#include "vector_liggghts.h"
#include "math_extra_liggghts.h"
#include "mpi_liggghts.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeSurface::ComputeSurface(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  nmax_(0),
  list_(0),
  is_on_surface_(0),
  use_com_(false),
  angle_(0.),
  cosine_sqr_(0.)
{
    vectorZeroize3D(point_up_);
    vectorZeroize3D(n_vec_up_);
    bool point_up_set = false;
    bool n_vec_up_set = false;

    // parse args

    int iarg = 3;
    if (narg < 3) error->compute_error(FLERR,this,"not enoguh arguments");

    bool hasargs = true;
    while(iarg < narg && hasargs)
    {
      hasargs = false;
      if (strcmp(arg[iarg],"angle") == 0) {
          if (narg < iarg+2)
            error->compute_error(FLERR,this,"not enough arguments for 'angle'");
          iarg++;
          angle_ = force->numeric(FLERR,arg[iarg++]);
          if(angle_ < 0. || angle_ > 45.)
            error->compute_error(FLERR,this,"0° < angle < 45° required");
          // 90 degrees difference in agle definition, so use sin here
          cosine_sqr_ = sin(angle_*M_PI/180.) * sin(angle_*M_PI/180.);
          hasargs = true;
      } else if (strcmp(arg[iarg],"point_up") == 0) {
          if (narg < iarg+4)
            error->compute_error(FLERR,this,"not enough arguments for 'point_up'");
          iarg++;
          point_up_[0] = force->numeric(FLERR,arg[iarg++]);
          point_up_[1] = force->numeric(FLERR,arg[iarg++]);
          point_up_[2] = force->numeric(FLERR,arg[iarg++]);
          point_up_set = true;
          hasargs = true;
      } else if (strcmp(arg[iarg],"n_vec_up") == 0) {
          if (narg < iarg+4)
            error->compute_error(FLERR,this,"not enough arguments for 'n_vec_up'");
          iarg++;
          n_vec_up_[0] = force->numeric(FLERR,arg[iarg++]);
          n_vec_up_[1] = force->numeric(FLERR,arg[iarg++]);
          n_vec_up_[2] = force->numeric(FLERR,arg[iarg++]);
          vectorNormalize3D(n_vec_up_);
          n_vec_up_set = true;
          hasargs = true;
      } else if (strcmp(arg[iarg],"use_com") == 0) {
          if (narg < iarg+2)
            error->compute_error(FLERR,this,"not enough arguments for 'use_com'");
          iarg++;
          if(strcmp(arg[iarg],"yes") == 0)
            use_com_ = true;
          else if(strcmp(arg[iarg],"no") == 0)
            use_com_ = false;
          else
            error->compute_error(FLERR,this,"valid arguments for 'use_com' are 'yes' or 'no'");
          hasargs = true;
          iarg++;
      } else if(strcmp(style,"surface") == 0) {
          char *errmsg = new char[strlen(arg[iarg])+50];
          sprintf(errmsg,"unknown keyword or wrong keyword order: %s", arg[iarg]);
          error->compute_error(FLERR,this,errmsg);
          delete []errmsg;
      }
    }

    if (!point_up_set)
        error->compute_error(FLERR,this,"requires 'point_up' be defined");
    if (!use_com_ && !n_vec_up_set)
        error->compute_error(FLERR,this,"requires either 'use_com = yes' or 'n_vec_up' be defined");
    if (use_com_ && n_vec_up_set)
        error->compute_error(FLERR,this,"must not use both 'use_com = yes' and 'n_vec_up'");
    if (MathExtraLiggghts::compDouble(cosine_sqr_,0.))
        error->compute_error(FLERR,this,"requires 'angle' be defined and 0° < angle < 45° required");

    // other settings

    peratom_flag = 1;
    size_peratom_cols = 0;
    comm_reverse = 1;
    scalar_flag = 1;
}

/* ---------------------------------------------------------------------- */

ComputeSurface::~ComputeSurface()
{
  memory->destroy(is_on_surface_);
}

/* ---------------------------------------------------------------------- */

void ComputeSurface::init()
{
    // error checks

    if (!atom->sphere_flag)
        error->compute_error(FLERR,this,"requires atom style sphere");

    if (force->pair == NULL)
        error->compute_error(FLERR,this,"requires a pair style be defined");

    // need an occasional neighbor list

    int irequest = neighbor->request((void *) this);
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->gran = 1;
    neighbor->requests[irequest]->pair = 0;
    neighbor->requests[irequest]->compute = 1;
    neighbor->requests[irequest]->occasional = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeSurface::init_list(int id, NeighList *ptr)
{
  list_ = ptr;
}

/* ---------------------------------------------------------------------- */

double ComputeSurface::compute_scalar()
{
  invoked_vector = update->ntimestep;

  if(invoked_peratom != update->ntimestep)
    compute_peratom();

  int nlocal = atom->nlocal;

  // sum up
  double n_on_surf = 0.;
  for (int i = 0; i < nlocal; i++)
    n_on_surf += vector_atom[i];

  MPI_Sum_Scalar(n_on_surf,world);

  return n_on_surf;
}

/* ---------------------------------------------------------------------- */

void ComputeSurface::compute_peratom()
{
  int i,j,ii,jj,inum,jnum;
  double xi[3],xi_to_xj[3];
  double com_to_up_unit[3], com[3];
  int *ilist,*jlist,*numneigh,**firstneigh;

  if(invoked_peratom == update->ntimestep)
    return;

  invoked_peratom = update->ntimestep;

  // grow is_on_surface_ array if necessary

  if (atom->nmax > nmax_) {
    memory->destroy(is_on_surface_);
    nmax_ = atom->nmax;
    memory->create(is_on_surface_,nmax_,"surface:is_on_surface_");
    vector_atom = is_on_surface_;
  }

  // invoke neighbor list (will copy or build if necessary)

  neighbor->build_one(list_->index);

  inum = list_->inum;
  ilist = list_->ilist;
  numneigh = list_->numneigh;
  firstneigh = list_->firstneigh;

  // compute number of contacts for each atom in group
  // contact if distance <= sum of radii
  // tally for both I and J

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  // first criterium - must be on right side of "up"

    double com_to_xi[3];

    // calculate groups xcm
    double masstotal = group->mass(igroup);
    group->xcm(igroup,masstotal,com);

    // unit vector from center of mass to "up"
    vectorSubtract3D(point_up_,com,com_to_up_unit);
    vectorNormalize3D(com_to_up_unit);

    for (i = 0; i < nall; i++) {
      vectorSubtract3D(x[i],com,com_to_xi);
      // particle on wrong side of "up" - might be "on surface"
      if(vectorDot3D(com_to_xi,use_com_ ? com_to_up_unit : n_vec_up_) > 0.)
        is_on_surface_[i] = 1.;
      // particle is for sure not "on surface"
      else
        is_on_surface_[i] = 0.;
    }

  // second criterion - neighbor angle
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit) {

      // particle is on right side - check

      vectorCopy3D(x[i],xi);
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        // skip if both particles on the wrong side
        if(MathExtraLiggghts::compDouble(is_on_surface_[i],0.) && MathExtraLiggghts::compDouble(is_on_surface_[j],0.) )
          continue;

        vectorSubtract3D(x[j],xi,xi_to_xj);

        double xi_to_xj_sqr = vectorMag3DSquared(xi_to_xj);
        double dot = vectorDot3D(xi_to_xj,use_com_ ? com_to_up_unit : n_vec_up_);
        double dotsqr = dot*dot;

        // check if particles fulfil angle criterion
        if(dot < 0 && dotsqr > cosine_sqr_*xi_to_xj_sqr)
            is_on_surface_[j] = 0.;
        if(dot > 0 && dotsqr > cosine_sqr_*xi_to_xj_sqr)
            is_on_surface_[i] = 0.;
      }
    }
  }

  // communicate ghost atom counts between neighbor procs if necessary

  if (force->newton_pair) comm->reverse_comm_compute(this);
}

/* ---------------------------------------------------------------------- */

int ComputeSurface::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    buf[m++] = is_on_surface_[i];
  return 1;
}

/* ---------------------------------------------------------------------- */

void ComputeSurface::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    is_on_surface_[j] = MathExtraLiggghts::compDouble(buf[m++],0.) ? 0. : 1.;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeSurface::memory_usage()
{
  double bytes = nmax_ * sizeof(double);
  return bytes;
}
