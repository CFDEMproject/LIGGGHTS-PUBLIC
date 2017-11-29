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
    Alexander Podlozhnyuk (DCS Computing GmbH, Linz)

    Copyright 2015-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#include <cmath>
#include "math_extra.h"
#include "math_extra_liggghts_nonspherical.h"
#include <stdlib.h>
#include <string.h>
#include "atom_vec_superquadric.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
//#include "domain_wedge.h"
#include "modify.h"
#include "force.h"
#include "fix.h"
#include "fix_adapt.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#ifdef SUPERQUADRIC_ACTIVE_FLAG
#include "math_extra_liggghts_superquadric.h"

#define DELTA 10000

/* ---------------------------------------------------------------------- */

AtomVecSuperquadric::AtomVecSuperquadric(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = 0;

  comm_x_only = 0; comm_f_only = 0;
  size_forward = 7;
  size_reverse = 6;
  size_border = 23;
  size_velocity = 6;
  size_data_atom = 8;
  size_data_vel = 7;
  xcol_data = 9;

  atom->superquadric_flag = 1;
  atom->radius_flag = atom->rmass_flag = atom->omega_flag = atom->density_flag =
    atom->torque_flag = atom->angmom_flag = 1;
  nvar_restart = 16 + 16 + 1;
}

void AtomVecSuperquadric::init()
{
  AtomVec::init();
  radvary = 0;
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by DELTA
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecSuperquadric::grow(int n)
{
  if (n == 0) nmax += DELTA;
  else nmax = n;
  atom->nmax = nmax;
  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  tag = memory->grow(atom->tag,nmax,"atom:tag");
  type = memory->grow(atom->type,nmax,"atom:type");
  mask = memory->grow(atom->mask,nmax,"atom:mask");
  image = memory->grow(atom->image,nmax,"atom:image");
  x = memory->grow(atom->x,nmax,3,"atom:x");
  v = memory->grow(atom->v,nmax,3,"atom:v");
  f = memory->grow(atom->f,nmax*comm->nthreads,3,"atom:f");

  radius = memory->grow(atom->radius,nmax,"atom:radius");
  density = memory->grow(atom->density,nmax,"atom:density");
  rmass = memory->grow(atom->rmass,nmax,"atom:rmass");
  omega = memory->grow(atom->omega,nmax,3,"atom:omega");
  torque = memory->grow(atom->torque,nmax*comm->nthreads,3,"atom:torque");
//Superquadric bonus------------------------------------
  shape = memory->grow(atom->shape,nmax,3,"atom:shape");
  blockiness = memory->grow(atom->blockiness,nmax,2,"atom:blockiness");
  inertia = memory->grow(atom->inertia,nmax,3,"atom:inertia");
  volume = memory->grow(atom->volume,nmax,"atom:volume");
  area = memory->grow(atom->area,nmax,"atom:area");
  quaternion = memory->grow(atom->quaternion,nmax,4,"atom:quaternion");
  angmom = memory->grow(atom->angmom,nmax,3,"atom:angmom");
//------------------------------------------------------
  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecSuperquadric::grow_reset()
{
  tag = atom->tag; type = atom->type;
  mask = atom->mask; image = atom->image;
  x = atom->x; v = atom->v; f = atom->f;
  radius = atom->radius; density = atom->density; rmass = atom->rmass;
  omega = atom->omega; torque = atom->torque;
//Superquadric bonus----------------------------------------
  shape = atom->shape; blockiness = atom->blockiness; inertia = atom->inertia;
  volume = atom->volume; area = atom->area; quaternion = atom->quaternion; angmom = atom->angmom;
//----------------------------------------------------------
}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
------------------------------------------------------------------------- */

void AtomVecSuperquadric::copy(int i, int j, int delflag)
{
  tag[j] = tag[i];
  type[j] = type[i];
  mask[j] = mask[i];
  image[j] = image[i];
  x[j][0] = x[i][0];
  x[j][1] = x[i][1];
  x[j][2] = x[i][2];
  v[j][0] = v[i][0];
  v[j][1] = v[i][1];
  v[j][2] = v[i][2];

  radius[j] = radius[i];
  rmass[j] = rmass[i];
  density[j] = density[i];
  omega[j][0] = omega[i][0];
  omega[j][1] = omega[i][1];
  omega[j][2] = omega[i][2];
//Superquadric bonus---------------------------
  shape[j][0] = shape[i][0];
  shape[j][1] = shape[i][1];
  shape[j][2] = shape[i][2];

  blockiness[j][0] = blockiness[i][0];
  blockiness[j][1] = blockiness[i][1];

  inertia[j][0] = inertia[i][0];
  inertia[j][1] = inertia[i][1];
  inertia[j][2] = inertia[i][2];

  volume[j] = volume[i];
  area[j] = area[i];

  quaternion[j][0] = quaternion[i][0];
  quaternion[j][1] = quaternion[i][1];
  quaternion[j][2] = quaternion[i][2];
  quaternion[j][3] = quaternion[i][3];

  angmom[j][0] = angmom[i][0];
  angmom[j][1] = angmom[i][1];
  angmom[j][2] = angmom[i][2];
//----------------------------------------------
  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j,delflag);
}

/* ---------------------------------------------------------------------- */

int AtomVecSuperquadric::pack_comm(int n, int *list, double *buf,
                               int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;
  if(radvary == 1)
    error->one(FLERR,"the case of radvary=1 for superquadrics is not implemented");

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
//Superquadric bonus----------------------------------
      buf[m++] = quaternion[j][0];
      buf[m++] = quaternion[j][1];
      buf[m++] = quaternion[j][2];
      buf[m++] = quaternion[j][3];
//----------------------------------------------------
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
//Superquadric bonus----------------------------------
      buf[m++] = quaternion[j][0];
      buf[m++] = quaternion[j][1];
      buf[m++] = quaternion[j][2];
      buf[m++] = quaternion[j][3];
//----------------------------------------------------
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSuperquadric::pack_comm_vel(int n, int *list, double *buf,
                                   int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;
  m = 0;
  if(radvary == 1)
      error->one(FLERR,"the case of radvary=1 for superquadrics is not implemented");
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
      buf[m++] = omega[j][0];
      buf[m++] = omega[j][1];
      buf[m++] = omega[j][2];
//Superquadric bonus----------------------------------
      buf[m++] = quaternion[j][0];
      buf[m++] = quaternion[j][1];
      buf[m++] = quaternion[j][2];
      buf[m++] = quaternion[j][3];
//----------------------------------------------------
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];

        buf[m++] = omega[j][0];
        buf[m++] = omega[j][1];
        buf[m++] = omega[j][2];
//Superquadric bonus----------------------------------
        buf[m++] = quaternion[j][0];
        buf[m++] = quaternion[j][1];
        buf[m++] = quaternion[j][2];
        buf[m++] = quaternion[j][3];
//-----------------------------------------------------
      }
    } else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
        }
        buf[m++] = omega[j][0];
        buf[m++] = omega[j][1];
        buf[m++] = omega[j][2];
//Superquadric bonus----------------------------------
        buf[m++] = quaternion[j][0];
        buf[m++] = quaternion[j][1];
        buf[m++] = quaternion[j][2];
        buf[m++] = quaternion[j][3];
//----------------------------------------------------
      }
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSuperquadric::pack_comm_hybrid(int n, int *list, double *buf)
{
  error->one(FLERR,"function AtomVecSuperquadric::pack_comm_hybrid is not implemented yet");
  return 0;
}

/* ---------------------------------------------------------------------- */

void AtomVecSuperquadric::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;
  m = 0;
  if(radvary == 1)
      error->one(FLERR,"the case of radvary=1 for superquadrics is not implemented");
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
//Superquadric bonus------------------------
    quaternion[i][0] = buf[m++];
    quaternion[i][1] = buf[m++];
    quaternion[i][2] = buf[m++];
    quaternion[i][3] = buf[m++];
//------------------------------------------
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecSuperquadric::unpack_comm_vel(int n, int first, double *buf)
{
  int i,m,last;
  m = 0;
  if(radvary == 1)
      error->one(FLERR,"the case of radvary=1 for superquadrics is not implemented");
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
    omega[i][0] = buf[m++];
    omega[i][1] = buf[m++];
    omega[i][2] = buf[m++];
//Superquadric bonus----------------------------------
    quaternion[i][0] = buf[m++];
    quaternion[i][1] = buf[m++];
    quaternion[i][2] = buf[m++];
    quaternion[i][3] = buf[m++];
//---------------------------------------------------
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecSuperquadric::unpack_comm_hybrid(int n, int first, double *buf)
{
  error->one(FLERR,"function AtomVecSuperquadric::unpack_comm_hybrid is not implemented yet");
  return 0;
}

/* ---------------------------------------------------------------------- */

int AtomVecSuperquadric::pack_reverse(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = f[i][0];
    buf[m++] = f[i][1];
    buf[m++] = f[i][2];
    buf[m++] = torque[i][0];
    buf[m++] = torque[i][1];
    buf[m++] = torque[i][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSuperquadric::pack_reverse_hybrid(int n, int first, double *buf)
{
  error->one(FLERR,"function AtomVecSuperquadric::pack_reverse_hybrid is not implemented yet");
  return 0;
}

/* ---------------------------------------------------------------------- */

void AtomVecSuperquadric::unpack_reverse(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
    torque[j][0] += buf[m++];
    torque[j][1] += buf[m++];
    torque[j][2] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecSuperquadric::unpack_reverse_hybrid(int n, int *list, double *buf)
{
  error->one(FLERR,"function AtomVecSuperquadric::unpack_reverse_hybrid is not implemented yet");
  return 0;
}

/* ---------------------------------------------------------------------- */

int AtomVecSuperquadric::pack_border(int n, int *list, double *buf,
                                 int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = radius[j];
      buf[m++] = rmass[j];
      buf[m++] = density[j];
//Superquadric bonus----------------------
      buf[m++] = shape[j][0];
      buf[m++] = shape[j][1];
      buf[m++] = shape[j][2];

      buf[m++] = blockiness[j][0];
      buf[m++] = blockiness[j][1];

      buf[m++] = inertia[j][0];
      buf[m++] = inertia[j][1];
      buf[m++] = inertia[j][2];

      buf[m++] = volume[j];
      buf[m++] = area[j];

      buf[m++] = quaternion[j][0];
      buf[m++] = quaternion[j][1];
      buf[m++] = quaternion[j][2];
      buf[m++] = quaternion[j][3];
//----------------------------------------
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = radius[j];
      buf[m++] = rmass[j];
      buf[m++] = density[j];
//Superquadric bonus----------------------
      buf[m++] = shape[j][0];
      buf[m++] = shape[j][1];
      buf[m++] = shape[j][2];

      buf[m++] = blockiness[j][0];
      buf[m++] = blockiness[j][1];

      buf[m++] = inertia[j][0];
      buf[m++] = inertia[j][1];
      buf[m++] = inertia[j][2];

      buf[m++] = volume[j];
      buf[m++] = area[j];

      buf[m++] = quaternion[j][0];
      buf[m++] = quaternion[j][1];
      buf[m++] = quaternion[j][2];
      buf[m++] = quaternion[j][3];
//----------------------------------------
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSuperquadric::pack_border_vel(int n, int *list, double *buf,
                                     int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;

  // error for DomainWedge
  /*if(dynamic_cast<DomainWedge*>(domain))
    error->all(FLERR, "not compatible with domain DomainWedge");*/

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = radius[j];
      buf[m++] = rmass[j];
      buf[m++] = density[j];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
      buf[m++] = omega[j][0];
      buf[m++] = omega[j][1];
      buf[m++] = omega[j][2];
//Superquadric bonus----------------------
      buf[m++] = shape[j][0];
      buf[m++] = shape[j][1];
      buf[m++] = shape[j][2];

      buf[m++] = blockiness[j][0];
      buf[m++] = blockiness[j][1];

      buf[m++] = inertia[j][0];
      buf[m++] = inertia[j][1];
      buf[m++] = inertia[j][2];

      buf[m++] = volume[j];
      buf[m++] = area[j];

      buf[m++] = quaternion[j][0];
      buf[m++] = quaternion[j][1];
      buf[m++] = quaternion[j][2];
      buf[m++] = quaternion[j][3];
//----------------------------------------
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        buf[m++] = radius[j];
        buf[m++] = rmass[j];
        buf[m++] = density[j];
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
        buf[m++] = omega[j][0];
        buf[m++] = omega[j][1];
        buf[m++] = omega[j][2];
//Superquadric bonus----------------------
        buf[m++] = shape[j][0];
        buf[m++] = shape[j][1];
        buf[m++] = shape[j][2];

        buf[m++] = blockiness[j][0];
        buf[m++] = blockiness[j][1];

        buf[m++] = inertia[j][0];
        buf[m++] = inertia[j][1];
        buf[m++] = inertia[j][2];

        buf[m++] = volume[j];
        buf[m++] = area[j];

        buf[m++] = quaternion[j][0];
        buf[m++] = quaternion[j][1];
        buf[m++] = quaternion[j][2];
        buf[m++] = quaternion[j][3];
//----------------------------------------
      }
    } else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        buf[m++] = radius[j];
        buf[m++] = rmass[j];
        buf[m++] = density[j];
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
        }
        buf[m++] = omega[j][0];
        buf[m++] = omega[j][1];
        buf[m++] = omega[j][2];
//Superquadric bonus----------------------
        buf[m++] = shape[j][0];
        buf[m++] = shape[j][1];
        buf[m++] = shape[j][2];

        buf[m++] = blockiness[j][0];
        buf[m++] = blockiness[j][1];

        buf[m++] = inertia[j][0];
        buf[m++] = inertia[j][1];
        buf[m++] = inertia[j][2];

        buf[m++] = volume[j];
        buf[m++] = area[j];

        buf[m++] = quaternion[j][0];
        buf[m++] = quaternion[j][1];
        buf[m++] = quaternion[j][2];
        buf[m++] = quaternion[j][3];
//----------------------------------------
      }
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSuperquadric::pack_border_hybrid(int n, int *list, double *buf)
{
  error->one(FLERR,"function AtomVecSuperquadric::pack_border_hybrid is not implemented yet");
  return 0;
}

/* ---------------------------------------------------------------------- */

void AtomVecSuperquadric::unpack_border(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (int) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    radius[i] = buf[m++];
    rmass[i] = buf[m++];
    density[i] = buf[m++];
//Superquadric bonus----------------------
    shape[i][0] = buf[m++];
    shape[i][1] = buf[m++];
    shape[i][2] = buf[m++];

    blockiness[i][0] = buf[m++];
    blockiness[i][1] = buf[m++];

    inertia[i][0] = buf[m++];
    inertia[i][1] = buf[m++];
    inertia[i][2] = buf[m++];

    volume[i] = buf[m++];
    area[i] = buf[m++];

    quaternion[i][0] = buf[m++];
    quaternion[i][1] = buf[m++];
    quaternion[i][2] = buf[m++];
    quaternion[i][3] = buf[m++];
//----------------------------------------
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

void AtomVecSuperquadric::unpack_border_vel(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (int) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    radius[i] = buf[m++];
    rmass[i] = buf[m++];
    density[i] = buf[m++];
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
    omega[i][0] = buf[m++];
    omega[i][1] = buf[m++];
    omega[i][2] = buf[m++];
//Superquadric bonus----------------------
    shape[i][0] = buf[m++];
    shape[i][1] = buf[m++];
    shape[i][2] = buf[m++];

    blockiness[i][0] = buf[m++];
    blockiness[i][1] = buf[m++];

    inertia[i][0] = buf[m++];
    inertia[i][1] = buf[m++];
    inertia[i][2] = buf[m++];

    volume[i] = buf[m++];
    area[i] = buf[m++];

    quaternion[i][0] = buf[m++];
    quaternion[i][1] = buf[m++];
    quaternion[i][2] = buf[m++];
    quaternion[i][3] = buf[m++];
//----------------------------------------

  }
  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

int AtomVecSuperquadric::unpack_border_hybrid(int n, int first, double *buf)
{
  error->one(FLERR,"function AtomVecSuperquadric::unpack_border_hybrid is not implemented yet");
  return 0;
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVecSuperquadric::pack_exchange(int i, double *buf)
{
  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;

  buf[m++] = radius[i];
  buf[m++] = rmass[i];
  buf[m++] = density[i];
  buf[m++] = omega[i][0];
  buf[m++] = omega[i][1];
  buf[m++] = omega[i][2];

//Superquadric bonus----------------------
  buf[m++] = shape[i][0];
  buf[m++] = shape[i][1];
  buf[m++] = shape[i][2];

  buf[m++] = blockiness[i][0];
  buf[m++] = blockiness[i][1];

  buf[m++] = inertia[i][0];
  buf[m++] = inertia[i][1];
  buf[m++] = inertia[i][2];

  buf[m++] = volume[i];
  buf[m++] = area[i];

  buf[m++] = quaternion[i][0];
  buf[m++] = quaternion[i][1];
  buf[m++] = quaternion[i][2];
  buf[m++] = quaternion[i][3];

  buf[m++] = angmom[i][0];
  buf[m++] = angmom[i][1];
  buf[m++] = angmom[i][2];
//----------------------------------------

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSuperquadric::unpack_exchange(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];
  tag[nlocal] = (int) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (tagint) ubuf(buf[m++]).i;

  radius[nlocal] = buf[m++];
  rmass[nlocal] = buf[m++];
  density[nlocal] = buf[m++];
  omega[nlocal][0] = buf[m++];
  omega[nlocal][1] = buf[m++];
  omega[nlocal][2] = buf[m++];
//Superquadric bonus----------------------
  shape[nlocal][0] = buf[m++];
  shape[nlocal][1] = buf[m++];
  shape[nlocal][2] = buf[m++];

  blockiness[nlocal][0] = buf[m++];
  blockiness[nlocal][1] = buf[m++];

  inertia[nlocal][0] = buf[m++];
  inertia[nlocal][1] = buf[m++];
  inertia[nlocal][2] = buf[m++];

  volume[nlocal] = buf[m++];
  area[nlocal] = buf[m++];

  quaternion[nlocal][0] = buf[m++];
  quaternion[nlocal][1] = buf[m++];
  quaternion[nlocal][2] = buf[m++];
  quaternion[nlocal][3] = buf[m++];

  angmom[nlocal][0] = buf[m++];
  angmom[nlocal][1] = buf[m++];
  angmom[nlocal][2] = buf[m++];
//----------------------------------------

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->
        unpack_exchange(nlocal,&buf[m]);

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   size of restart data for all atoms owned by this proc
   include extra data stored by fixes
------------------------------------------------------------------------- */

int AtomVecSuperquadric::size_restart()
{
  int i;

  int nlocal = atom->nlocal;
  int n = (nvar_restart + 1) * nlocal;

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      for (i = 0; i < nlocal; i++)
        n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);

  return n;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file including extra quantities
   xyz must be 1st 3 values, so that read_restart can test on them
   molecular types may be negative, but write as positive
------------------------------------------------------------------------- */

int AtomVecSuperquadric::pack_restart(int i, double *buf)
{
  int m = 1;
  buf[m++] = x[i][0]; //1
  buf[m++] = x[i][1]; //2
  buf[m++] = x[i][2]; //3
  buf[m++] = ubuf(tag[i]).d;  //4
  buf[m++] = ubuf(type[i]).d; //5
  buf[m++] = ubuf(mask[i]).d; //6
  buf[m++] = ubuf(image[i]).d; //7

  buf[m++] = v[i][0];  //8
  buf[m++] = v[i][1];  //9
  buf[m++] = v[i][2];  //10

  buf[m++] = radius[i];  //11
  buf[m++] = rmass[i];  //12
  buf[m++] = density[i];  //13
  buf[m++] = omega[i][0];  //14
  buf[m++] = omega[i][1];  //15
  buf[m++] = omega[i][2]; //16
//Superquadric bonus-------------------------------
  buf[m++] = quaternion[i][0]; //17
  buf[m++] = quaternion[i][1];  //18
  buf[m++] = quaternion[i][2]; //19
  buf[m++] = quaternion[i][3];  //20

  buf[m++] = shape[i][0]; //21
  buf[m++] = shape[i][1]; //22
  buf[m++] = shape[i][2];  //23

  buf[m++] = blockiness[i][0];  //24
  buf[m++] = blockiness[i][1];  //25

  buf[m++] = inertia[i][0];  //26
  buf[m++] = inertia[i][1];  //27
  buf[m++] = inertia[i][2];  //28

  buf[m++] = volume[i];  //29
  buf[m++] = area[i];  //30

  buf[m++] = angmom[i][0];  //31
  buf[m++] = angmom[i][1];  //32
  buf[m++] = angmom[i][2];  //33
//-------------------------------------------------
  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
------------------------------------------------------------------------- */

int AtomVecSuperquadric::unpack_restart(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) {
    grow(0);
    if (atom->nextra_store)
      memory->grow(atom->extra,nmax,atom->nextra_store,"atom:extra");
  }

  int m = 1;
  x[nlocal][0] = buf[m++]; //1
  x[nlocal][1] = buf[m++]; //2
  x[nlocal][2] = buf[m++]; //3
  tag[nlocal] = (int) ubuf(buf[m++]).i; //4
  type[nlocal] = (int) ubuf(buf[m++]).i; //5
  mask[nlocal] = (int) ubuf(buf[m++]).i; //6
  image[nlocal] = (tagint) ubuf(buf[m++]).i; //7

  v[nlocal][0] = buf[m++]; //8
  v[nlocal][1] = buf[m++]; //9
  v[nlocal][2] = buf[m++]; //10

  radius[nlocal] = buf[m++]; //11
  rmass[nlocal] = buf[m++]; //12
  density[nlocal] = buf[m++]; //13
  omega[nlocal][0] = buf[m++]; //14
  omega[nlocal][1] = buf[m++]; //15
  omega[nlocal][2] = buf[m++]; //16
//Superquadric bonus--------------------------
  quaternion[nlocal][0] = buf[m++]; //17
  quaternion[nlocal][1] = buf[m++]; //18
  quaternion[nlocal][2] = buf[m++]; //19
  quaternion[nlocal][3] = buf[m++]; //20

  shape[nlocal][0] = buf[m++]; //21
  shape[nlocal][1] = buf[m++]; //22
  shape[nlocal][2] = buf[m++]; //23

  blockiness[nlocal][0] = buf[m++]; //24
  blockiness[nlocal][1] = buf[m++]; //25

  inertia[nlocal][0] = buf[m++]; //26
  inertia[nlocal][1] = buf[m++]; //27
  inertia[nlocal][2] = buf[m++]; //28

  volume[nlocal] = buf[m++]; //29
  area[nlocal] = buf[m++]; //30

  angmom[nlocal][0] = buf[m++]; //31
  angmom[nlocal][1] = buf[m++]; //32
  angmom[nlocal][2] = buf[m++]; //33
//---------------------------------------------
  double **extra = atom->extra;
  if (atom->nextra_store) {
    int size = static_cast<int> (buf[0]) - m;

    for (int i = 0; i < size; i++) extra[nlocal][i] = buf[m++];
  }

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   create one atom of itype at coord
   set other values to defaults
------------------------------------------------------------------------- */

void AtomVecSuperquadric::create_atom(int itype, double *coord)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = 0;
  type[nlocal] = itype;
  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];
  mask[nlocal] = 1;
  image[nlocal] = ((tagint) IMGMAX << IMG2BITS) |
    ((tagint) IMGMAX << IMGBITS) | IMGMAX;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  radius[nlocal] = 0.5;
  density[nlocal] = 1.0;

  omega[nlocal][0] = 0.0;
  omega[nlocal][1] = 0.0;
  omega[nlocal][2] = 0.0;
//Superquadric bonus-----------------------
  quatIdentity4D(quaternion[nlocal]);

  shape[nlocal][0] = radius[nlocal];
  shape[nlocal][1] = radius[nlocal];
  shape[nlocal][2] = radius[nlocal];

  blockiness[nlocal][0] = blockiness[nlocal][1] = 2.0;

  angmom[nlocal][0] = angmom[nlocal][1] = angmom[nlocal][2] = 0.0;
//------------------------------------------
  MathExtraLiggghtsNonspherical::volume_superquadric(shape[nlocal], blockiness[nlocal], (volume+nlocal));
  if (domain->dimension == 3)
    rmass[nlocal] = volume[nlocal] * density[nlocal];
  else
    error->one(FLERR,"Superquadrics in 2D are not implemented");
  MathExtraLiggghtsNonspherical::inertia_superquadric(shape[nlocal], blockiness[nlocal], density[nlocal], inertia[nlocal]);
  MathExtraLiggghtsNonspherical::area_superquadric(shape[nlocal], blockiness[nlocal], area+nlocal);
  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecSuperquadric::data_atom(double *coord, tagint imagetmp, char **values)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = atoi(values[0]);
  if (tag[nlocal] <= 0)
    error->one(FLERR,"Invalid atom ID in Atoms section of data file");

  type[nlocal] = atoi(values[1]);
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR,"Invalid atom type in Atoms section of data file");

  density[nlocal] = atof(values[2]);
  if (density[nlocal] <= 0.0)
    error->one(FLERR,"Invalid density in Atoms section of data file");

  shape[nlocal][0] = atof(values[3]);
  if (shape[nlocal][0] <= 0.0)
    error->one(FLERR,"Invalid shape in Atoms section of data file");

  shape[nlocal][1] = atof(values[4]);
  if (shape[nlocal][1] <= 0.0)
      error->one(FLERR,"Invalid shape in Atoms section of data file");

  shape[nlocal][2] = atof(values[5]);
  if (shape[nlocal][2] <= 0.0)
      error->one(FLERR,"Invalid shape in Atoms section of data file");

  blockiness[nlocal][0] = atof(values[6]);
  if (blockiness[nlocal][0] < 2.0)
      error->one(FLERR,"Invalid blockiness in Atoms section of data file");
  blockiness[nlocal][1] = atof(values[7]);
  if (blockiness[nlocal][1] < 2.0)
      error->one(FLERR,"Invalid blockiness in Atoms section of data file");

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  image[nlocal] = imagetmp;

  mask[nlocal] = 1;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;
  omega[nlocal][0] = 0.0;
  omega[nlocal][1] = 0.0;
  omega[nlocal][2] = 0.0;

  angmom[nlocal][0] = angmom[nlocal][1] = angmom[nlocal][2] = 0.0;

  quatIdentity4D(quaternion[nlocal]);
  MathExtraLiggghtsNonspherical::bounding_sphere_radius_superquadric(shape[nlocal], blockiness[nlocal], radius+nlocal);

  MathExtraLiggghtsNonspherical::volume_superquadric(shape[nlocal], blockiness[nlocal], volume+nlocal);
  if (domain->dimension == 3)
    rmass[nlocal] = volume[nlocal] * density[nlocal];
  else
    error->one(FLERR,"Superquadrics in 2D are not implemented");
  MathExtraLiggghtsNonspherical::inertia_superquadric(shape[nlocal], blockiness[nlocal], density[nlocal], inertia[nlocal]);
  MathExtraLiggghtsNonspherical::area_superquadric(shape[nlocal], blockiness[nlocal], area+nlocal);
  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Atoms section of data file
   initialize other atom quantities for this sub-style
------------------------------------------------------------------------- */

int AtomVecSuperquadric::data_atom_hybrid(int nlocal, char **values)
{
  error->one(FLERR,"function AtomVecSuperquadric::data_atom_hybrid is not implemented yet");
  return 0;
}

/* ----------------------------------------------------------------------
   unpack one line from Velocities section of data file
------------------------------------------------------------------------- */

void AtomVecSuperquadric::data_vel(int m, char **values)
{
  v[m][0] = atof(values[0]);
  v[m][1] = atof(values[1]);
  v[m][2] = atof(values[2]);
  omega[m][0] = atof(values[3]);
  omega[m][1] = atof(values[4]);
  omega[m][2] = atof(values[5]);
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Velocities section of data file
------------------------------------------------------------------------- */

int AtomVecSuperquadric::data_vel_hybrid(int m, char **values)
{
  error->one(FLERR,"function AtomVecSuperquadric::data_vel_hybrid is not implemented yet");
  return 0;
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecSuperquadric::pack_data(double **buf)
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = ubuf(type[i]).d;
    buf[i][2] = density[i];
    buf[i][3] = volume[i];
    buf[i][4] = x[i][0];
    buf[i][5] = x[i][1];
    buf[i][6] = x[i][2];
    buf[i][7] = ubuf((image[i] & IMGMASK) - IMGMAX).d;
    buf[i][8] = ubuf((image[i] >> IMGBITS & IMGMASK) - IMGMAX).d;
    buf[i][9] = ubuf((image[i] >> IMG2BITS) - IMGMAX).d;
    buf[i][10] = shape[i][0];
    buf[i][11] = shape[i][1];
    buf[i][12] = shape[i][2];
    buf[i][13] = blockiness[i][0];
    buf[i][14] = blockiness[i][1];
    buf[i][15] = quaternion[i][0];
    buf[i][16] = quaternion[i][1];
    buf[i][17] = quaternion[i][2];
    buf[i][18] = quaternion[i][3];
    buf[i][19] = inertia[i][0];
    buf[i][20] = inertia[i][1];
    buf[i][21] = inertia[i][2];
    buf[i][22] = area[i];
  }
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecSuperquadric::pack_data(double **buf,int tag_offset)
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]+tag_offset).d;
    buf[i][1] = ubuf(type[i]).d;
    buf[i][2] = density[i];
    buf[i][3] = volume[i];
    buf[i][4] = x[i][0];
    buf[i][5] = x[i][1];
    buf[i][6] = x[i][2];
    buf[i][7] = ubuf((image[i] & IMGMASK) - IMGMAX).d;
    buf[i][8] = ubuf((image[i] >> IMGBITS & IMGMASK) - IMGMAX).d;
    buf[i][9] = ubuf((image[i] >> IMG2BITS) - IMGMAX).d;
    buf[i][10] = shape[i][0];
    buf[i][11] = shape[i][1];
    buf[i][12] = shape[i][2];
    buf[i][13] = blockiness[i][0];
    buf[i][14] = blockiness[i][1];
    buf[i][15] = quaternion[i][0];
    buf[i][16] = quaternion[i][1];
    buf[i][17] = quaternion[i][2];
    buf[i][18] = quaternion[i][3];
    buf[i][19] = inertia[i][0];
    buf[i][20] = inertia[i][1];
    buf[i][21] = inertia[i][2];
    buf[i][22] = area[i];
  }
}

/* ----------------------------------------------------------------------
   pack hybrid atom info for data file
------------------------------------------------------------------------- */

int AtomVecSuperquadric::pack_data_hybrid(int i, double *buf)
{
  error->one(FLERR,"function AtomVecSuperquadric::pack_data_hybrid is not implemented yet");
  return 0;
}

/* ----------------------------------------------------------------------
   write atom info to data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecSuperquadric::write_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,"%d %d %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %d %d %d\n",
            (int) ubuf(buf[i][0]).i,(int) ubuf(buf[i][1]).i,
            buf[i][2],buf[i][3],
            buf[i][4],buf[i][5],buf[i][6],
            (int) ubuf(buf[i][7]).i,(int) ubuf(buf[i][8]).i,
            (int) ubuf(buf[i][9]).i);
}

/* ----------------------------------------------------------------------
   write hybrid atom info to data file
------------------------------------------------------------------------- */

int AtomVecSuperquadric::write_data_hybrid(FILE *fp, double *buf)
{
  error->one(FLERR,"function AtomVecSuperquadric::write_data_hybrid is not implemented yet");
  return 0;
}

/* ----------------------------------------------------------------------
   pack velocity info for data file
------------------------------------------------------------------------- */

void AtomVecSuperquadric::pack_vel(double **buf)
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = v[i][0];
    buf[i][2] = v[i][1];
    buf[i][3] = v[i][2];
    buf[i][4] = omega[i][0];
    buf[i][5] = omega[i][1];
    buf[i][6] = omega[i][2];
  }
}

/* ----------------------------------------------------------------------
   pack velocity info for data file
------------------------------------------------------------------------- */

void AtomVecSuperquadric::pack_vel(double **buf,int tag_offset)
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]+tag_offset).d;
    buf[i][1] = v[i][0];
    buf[i][2] = v[i][1];
    buf[i][3] = v[i][2];
    buf[i][4] = omega[i][0];
    buf[i][5] = omega[i][1];
    buf[i][6] = omega[i][2];
  }
}

/* ----------------------------------------------------------------------
   pack hybrid velocity info for data file
------------------------------------------------------------------------- */

int AtomVecSuperquadric::pack_vel_hybrid(int i, double *buf)
{
  error->one(FLERR,"function AtomVecSuperquadric::pack_vel_hybrid is not implemented yet");
  return 0;
}

/* ----------------------------------------------------------------------
   write velocity info to data file
------------------------------------------------------------------------- */

void AtomVecSuperquadric::write_vel(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,"%d %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e\n",
            (int) ubuf(buf[i][0]).i,buf[i][1],buf[i][2],buf[i][3],
            buf[i][4],buf[i][5],buf[i][6]);
}

/* ----------------------------------------------------------------------
   write hybrid velocity info to data file
------------------------------------------------------------------------- */

int AtomVecSuperquadric::write_vel_hybrid(FILE *fp, double *buf)
{
	error->one(FLERR,"function AtomVecSuperquadric::write_vel_hybrid is not implemented yet");
  return 0;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint AtomVecSuperquadric::memory_usage()
{
  bigint bytes = 0;

  if (atom->memcheck("tag")) bytes += memory->usage(tag,nmax);
  if (atom->memcheck("type")) bytes += memory->usage(type,nmax);
  if (atom->memcheck("mask")) bytes += memory->usage(mask,nmax);
  if (atom->memcheck("image")) bytes += memory->usage(image,nmax);
  if (atom->memcheck("x")) bytes += memory->usage(x,nmax,3);
  if (atom->memcheck("v")) bytes += memory->usage(v,nmax,3);
  if (atom->memcheck("f")) bytes += memory->usage(f,nmax*comm->nthreads,3);

  if (atom->memcheck("density")) bytes += memory->usage(density,nmax);
  if (atom->memcheck("radius")) bytes += memory->usage(radius,nmax);
  if (atom->memcheck("rmass")) bytes += memory->usage(rmass,nmax);
  if (atom->memcheck("omega")) bytes += memory->usage(omega,nmax,3);
  if (atom->memcheck("torque"))
    bytes += memory->usage(torque,nmax*comm->nthreads,3);
//Superquadric bonus----------------------------------------------------
  if (atom->memcheck("shape")) bytes += memory->usage(shape,nmax,3);
  if (atom->memcheck("blockiness")) bytes += memory->usage(blockiness,nmax,2);
  if (atom->memcheck("inertia")) bytes += memory->usage(inertia,nmax,3);
  if (atom->memcheck("volume")) bytes += memory->usage(volume,nmax);
  if (atom->memcheck("area")) bytes += memory->usage(area,nmax);
  if (atom->memcheck("quaternion")) bytes += memory->usage(quaternion,nmax,4);
  if (atom->memcheck("angmom")) bytes += memory->usage(angmom,nmax,3);
//----------------------------------------------------

  return bytes;
}

//get pointer to quaternion of i-th particle
double* AtomVecSuperquadric::return_quat_ptr(int i)
{
  return quaternion[i];
}

//set shape parameters - half axes a,b,c for superquadric
void AtomVecSuperquadric::set_shape(int i, double shapex, double shapey, double shapez)
{
  shape[i][0] = shapex;
  shape[i][1] = shapey;
  shape[i][2] = shapez;
}

//set blockiness parameters
void AtomVecSuperquadric::set_blockiness(int i, double blockiness1, double blockiness2)
{
  blockiness[i][0] = blockiness1;
  blockiness[i][1] = blockiness2;
}

#endif
