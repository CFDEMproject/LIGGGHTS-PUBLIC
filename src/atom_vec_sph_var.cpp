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
    Andreas Aigner (JKU Linz))

    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#include <cmath>
#include <stdlib.h>
#include <string.h>
#include "atom_vec_sph_var.h"
#include "atom.h"
#include "domain.h"
#include "modify.h"
#include "force.h"
#include "fix.h"
#include "fix_adapt.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define DELTA 10000

/* ---------------------------------------------------------------------- */

AtomVecSPH2::AtomVecSPH2(LAMMPS *lmp) :
  AtomVec(lmp)
{
  molecular = 0;

  comm_x_only = 0; 
  comm_f_only = 0; 
  size_forward = 6; 
  size_reverse = 5; 
  size_border = 11;  
  size_velocity = 3; 
  size_data_atom = 8; 
  size_data_vel = 4; 
  xcol_data = 6;

  atom->p_flag = atom->rho_flag = atom->e_flag = 1;

  atom->radius_flag = atom->rmass_flag = 1;
  mass_type = 0; // is default!
  radvary = 0;

}

/* ---------------------------------------------------------------------- */

void AtomVecSPH2::init()
{
  AtomVec::init();

  // set radvary if particle diameters are time-varying due to some fix
  for (int i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->rad_mass_vary_flag) {
      radvary = 1;
      size_forward = 8;  
    }

  if(radvary) atom->radvary_flag = 1;
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by DELTA
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecSPH2::grow(int n)
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

  p = memory->grow(atom->p,nmax,"atom:pressure");
  rho = memory->grow(atom->rho,nmax,"atom:rho");
  drho = memory->grow(atom->drho,nmax,"atom:drho");
  e = memory->grow(atom->e,nmax,"atom:e");
  de = memory->grow(atom->e,nmax,"atom:de");

  radius = memory->grow(atom->radius,nmax,"atom:radius");
  rmass = memory->grow(atom->rmass,nmax,"atom:rmass");

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecSPH2::grow_reset()
{
  tag = atom->tag; type = atom->type;
  mask = atom->mask; image = atom->image;
  x = atom->x; v = atom->v; f = atom->f;
  p = atom->p; rho = atom->rho; drho = atom->drho;
  e = atom->e; de = atom->de;
  radius = atom->radius; rmass = atom->rmass;
}

/* ---------------------------------------------------------------------- */

void AtomVecSPH2::copy(int i, int j, int delflag)
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

  p[j] = p[i];
  rho[j] = rho[i];
  drho[j] = drho[i];
  e[j] = e[i];
  de[j] = de[i];

  radius[j] = radius[i];
  rmass[j] = rmass[i];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j,delflag);
}

/* ---------------------------------------------------------------------- */

int AtomVecSPH2::pack_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

  if (radvary == 0) {
    m = 0;
    if (pbc_flag == 0) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0];
        buf[m++] = x[j][1];
        buf[m++] = x[j][2];
        buf[m++] = p[j];
        buf[m++] = rho[j]; 
        buf[m++] = e[j];
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
        buf[m++] = p[j];
        buf[m++] = rho[j]; 
        buf[m++] = e[j];
      }
    }
  } else {
    m = 0;
    if (pbc_flag == 0) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0];
        buf[m++] = x[j][1];
        buf[m++] = x[j][2];
        buf[m++] = p[j];
        buf[m++] = rho[j]; 
        buf[m++] = e[j];
        buf[m++] = radius[j];
        buf[m++] = rmass[j];
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
        buf[m++] = p[j];
        buf[m++] = rho[j]; 
        buf[m++] = e[j];
        buf[m++] = radius[j];
        buf[m++] = rmass[j];
      }
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSPH2::pack_comm_vel(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

  if (radvary == 0) {
    m = 0;
    if (pbc_flag == 0) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0];
        buf[m++] = x[j][1];
        buf[m++] = x[j][2];
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
        buf[m++] = p[j];
        buf[m++] = rho[j]; 
        buf[m++] = e[j];
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
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
        buf[m++] = p[j];
        buf[m++] = rho[j]; 
        buf[m++] = e[j];
      }
    }
  } else {
    m = 0;
    if (pbc_flag == 0) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0];
        buf[m++] = x[j][1];
        buf[m++] = x[j][2];
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
        buf[m++] = p[j];
        buf[m++] = rho[j]; 
        buf[m++] = e[j];
        buf[m++] = radius[j];
        buf[m++] = rmass[j];
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
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
        buf[m++] = p[j];
        buf[m++] = rho[j]; 
        buf[m++] = e[j];
        buf[m++] = radius[j];
        buf[m++] = rmass[j];

      }
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSPH2::pack_comm_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  if (radvary == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = p[j];
      buf[m++] = rho[j]; 
      buf[m++] = e[j];
    }
  } else {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = p[j];
      buf[m++] = rho[j]; 
      buf[m++] = e[j];
      buf[m++] = radius[j];
      buf[m++] = rmass[j];
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecSPH2::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  if (radvary == 0) {
    m = 0;
    last = first + n;
    for (i = first; i < last; i++) {
      x[i][0] = buf[m++];
      x[i][1] = buf[m++];
      x[i][2] = buf[m++];
      p[i] = buf[m++];
      rho[i] = buf[m++];
      e[i] = buf[m++];
    }
  } else {
    m = 0;
    last = first + n;
    for (i = first; i < last; i++) {
      x[i][0] = buf[m++];
      x[i][1] = buf[m++];
      x[i][2] = buf[m++];
      p[i] = buf[m++];
      rho[i] = buf[m++];
      e[i] = buf[m++];
      radius[i] = buf[m++];
      rmass[i] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecSPH2::unpack_comm_vel(int n, int first, double *buf)
{
  int i,m,last;

  if (radvary == 0) {
    m = 0;
    last = first + n;
    for (i = first; i < last; i++) {
      x[i][0] = buf[m++];
      x[i][1] = buf[m++];
      x[i][2] = buf[m++];
      v[i][0] = buf[m++];
      v[i][1] = buf[m++];
      v[i][2] = buf[m++];
      p[i] = buf[m++];
      rho[i] = buf[m++];
      e[i] = buf[m++];
    }
  } else {
    m = 0;
    last = first + n;
    for (i = first; i < last; i++) {
      x[i][0] = buf[m++];
      x[i][1] = buf[m++];
      x[i][2] = buf[m++];
      v[i][0] = buf[m++];
      v[i][1] = buf[m++];
      v[i][2] = buf[m++];
      p[i] = buf[m++];
      rho[i] = buf[m++];
      e[i] = buf[m++];
      radius[i] = buf[m++];
      rmass[i] = buf[m++];

    }
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecSPH2::unpack_comm_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  if (radvary == 0) {
    for (i = first; i < last; i++) {
      p[i] = buf[m++];
      rho[i] = buf[m++];
      e[i] = buf[m++];
    }
  } else {
    for (i = first; i < last; i++) {
      p[i] = buf[m++];
      rho[i] = buf[m++];
      e[i] = buf[m++];
      radius[i] = buf[m++];
      rmass[i] = buf[m++];
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSPH2::pack_reverse(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = f[i][0];
    buf[m++] = f[i][1];
    buf[m++] = f[i][2];
    buf[m++] = drho[i];
    buf[m++] = de[i];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecSPH2::unpack_reverse(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
    drho[j] += buf[m++];
    de[j] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecSPH2::pack_border(int n, int *list, double *buf, int pbc_flag, int *pbc)
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
      buf[m++] = tag[j];
      buf[m++] = type[j];
      buf[m++] = mask[j];
      buf[m++] = p[j];
      buf[m++] = rho[j]; 
      buf[m++] = e[j];
      buf[m++] = radius[j];
      buf[m++] = rmass[j];
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
      buf[m++] = tag[j];
      buf[m++] = type[j];
      buf[m++] = mask[j];
      buf[m++] = p[j];
      buf[m++] = rho[j]; 
      buf[m++] = e[j];
      buf[m++] = radius[j];
      buf[m++] = rmass[j];
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSPH2::pack_border_vel(int n, int *list, double *buf, int pbc_flag, int *pbc)
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
      buf[m++] = tag[j];
      buf[m++] = type[j];
      buf[m++] = mask[j];
      buf[m++] = p[j];
      buf[m++] = rho[j]; 
      buf[m++] = e[j];
      buf[m++] = radius[j];
      buf[m++] = rmass[j];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
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
      buf[m++] = tag[j];
      buf[m++] = type[j];
      buf[m++] = mask[j];
      buf[m++] = p[j];
      buf[m++] = rho[j]; 
      buf[m++] = e[j];
      buf[m++] = radius[j];
      buf[m++] = rmass[j];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSPH2::pack_border_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = p[j];
    buf[m++] = rho[j]; 
    buf[m++] = e[j];
    buf[m++] = radius[j];
    buf[m++] = rmass[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecSPH2::unpack_border(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = static_cast<int> (buf[m++]);
    type[i] = static_cast<int> (buf[m++]);
    mask[i] = static_cast<int> (buf[m++]);
    p[i] = buf[m++];
    rho[i] = buf[m++];
    e[i] = buf[m++];
    radius[i] = buf[m++];
    rmass[i] = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

void AtomVecSPH2::unpack_border_vel(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = static_cast<int> (buf[m++]);
    type[i] = static_cast<int> (buf[m++]);
    mask[i] = static_cast<int> (buf[m++]);
    p[i] = buf[m++];
    rho[i] = buf[m++];
    e[i] = buf[m++];
    radius[i] = buf[m++];
    rmass[i] = buf[m++];
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

int AtomVecSPH2::unpack_border_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    p[i] = buf[m++];
    rho[i] = buf[m++];
    e[i] = buf[m++];
    radius[i] = buf[m++];
    rmass[i] = buf[m++];
  }
  return m;
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVecSPH2::pack_exchange(int i, double *buf)
{
  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];
  buf[m++] = tag[i];
  buf[m++] = type[i];
  buf[m++] = mask[i];
  buf[m] = 0.0;      // for valgrind
  *((tagint *) &buf[m++]) = image[i];

  buf[m++] = p[i];
  buf[m++] = rho[i];
  buf[m++] = e[i];
  buf[m++] = radius[i];
  buf[m++] = rmass[i];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSPH2::unpack_exchange(double *buf)
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
  tag[nlocal] = static_cast<int> (buf[m++]);
  type[nlocal] = static_cast<int> (buf[m++]);
  mask[nlocal] = static_cast<int> (buf[m++]);
  image[nlocal] = *((tagint *) &buf[m++]);

  p[nlocal] = buf[m++];
  rho[nlocal] =  buf[m++];
  e[nlocal] =  buf[m++];
  radius[nlocal] = buf[m++];
  rmass[nlocal] = buf[m++];

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

int AtomVecSPH2::size_restart()
{
  int i;

  int nlocal = atom->nlocal;
  int n = 16 * nlocal; // 11 + p + rho + e + radius + rmass

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

int AtomVecSPH2::pack_restart(int i, double *buf)
{
  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = tag[i];
  buf[m++] = type[i];
  buf[m++] = mask[i];
  buf[m] = 0.0;      // for valgrind
  *((tagint *) &buf[m++]) = image[i];
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];

  buf[m++] = p[i];
  buf[m++] = rho[i];
  buf[m++] = e[i];
  buf[m++] = radius[i];
  buf[m++] = rmass[i];

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
------------------------------------------------------------------------- */

int AtomVecSPH2::unpack_restart(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) {
    grow(0);
    if (atom->nextra_store)
      memory->grow(atom->extra,nmax,atom->nextra_store,"atom:extra");
  }

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  tag[nlocal] = static_cast<int> (buf[m++]);
  type[nlocal] = static_cast<int> (buf[m++]);
  mask[nlocal] = static_cast<int> (buf[m++]);
  image[nlocal] = *((tagint *) &buf[m++]);
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];

  p[nlocal] =  buf[m++];
  rho[nlocal] =  buf[m++];
  e[nlocal] =  buf[m++];
  radius[nlocal] = buf[m++];
  rmass[nlocal] = buf[m++];

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

void AtomVecSPH2::create_atom(int itype, double *coord)
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

  p[nlocal] = 0.0;
  rho[nlocal] = 0.0;
  drho[nlocal] = 0.0;
  e[nlocal] = 0.0;
  de[nlocal] = 0.0;
  radius[nlocal] = 0.5;
  rmass[nlocal] = 1.0;
  // in case of sph an calculation of the mass due density and radius doesn't make any sense

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecSPH2::data_atom(double *coord, int imagetmp, char **values)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = force->inumeric(FLERR,values[0]);
  if (tag[nlocal] <= 0)
    error->one(FLERR,"Invalid atom ID in Atoms section of data file");

  type[nlocal] = force->inumeric(FLERR,values[1]);
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR,"Invalid atom type in Atoms section of data file");

  rho[nlocal] = force->numeric(FLERR,values[2]);
  if (rho[nlocal] <= 0.0)
    error->one(FLERR,"Invalid rho in Atoms section of data file");

  radius[nlocal] = 0.5 * force->numeric(FLERR,values[3]);
  if (radius[nlocal] < 0.0)
    error->one(FLERR,"Invalid radius in Atoms section of data file");

  rmass[nlocal] = force->numeric(FLERR,values[4]);
  if (rmass[nlocal] <= 0.0)
    error->one(FLERR,"Invalid rmass in Atoms section of data file");

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  image[nlocal] = imagetmp;

  mask[nlocal] = 1;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  p[nlocal] = 0.0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Atoms section of data file
   initialize other atom quantities for this sub-style
------------------------------------------------------------------------- */

int AtomVecSPH2::data_atom_hybrid(int nlocal, char **values)
{
  rho[nlocal] = force->numeric(FLERR,values[0]);
  if (rho[nlocal] <= 0.0)
    error->one(FLERR,"Invalid rho in Atoms section of data file");

  radius[nlocal] = 0.5 * force->numeric(FLERR,values[1]);
  if (radius[nlocal] < 0.0)
    error->one(FLERR,"Invalid radius in Atoms section of data file");

  rmass[nlocal] = force->numeric(FLERR,values[2]);
  if (rmass[nlocal] <= 0.0)
    error->one(FLERR,"Invalid rmass in Atoms section of data file");

  return 3;
}

/* ----------------------------------------------------------------------
   unpack one line from Velocities section of data file
------------------------------------------------------------------------- */

void AtomVecSPH2::data_vel(int m, char **values)
{
  v[m][0] = force->numeric(FLERR,values[0]);
  v[m][1] = force->numeric(FLERR,values[1]);
  v[m][2] = force->numeric(FLERR,values[2]);
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Velocities section of data file
------------------------------------------------------------------------- */

int AtomVecSPH2::data_vel_hybrid(int m, char **values)
{
  return 0;
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecSPH2::pack_data(double **buf)
{
  
  error->all(FLERR,"This feature is not supported by SPH");
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = tag[i];
    buf[i][1] = type[i];
    buf[i][2] = 2.0*radius[i];
    if (radius[i] == 0.0) buf[i][3] = rmass[i];
    else
      buf[i][3] = rmass[i] / (4.0*M_PI/3.0 * radius[i]*radius[i]*radius[i]);
    buf[i][4] = x[i][0];
    buf[i][5] = x[i][1];
    buf[i][6] = x[i][2];
    buf[i][7] = (image[i] & IMGMASK) - IMGMAX;
    buf[i][8] = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
    buf[i][9] = (image[i] >> IMG2BITS) - IMGMAX;
  }
}

/* ----------------------------------------------------------------------
   pack hybrid atom info for data file
------------------------------------------------------------------------- */

int AtomVecSPH2::pack_data_hybrid(int i, double *buf)
{
  
  error->all(FLERR,"This feature is not supported by SPH");
  buf[0] = 2.0*radius[i];
  if (radius[i] == 0.0) buf[1] = rmass[i];
  else buf[1] = rmass[i] / (4.0*M_PI/3.0 * radius[i]*radius[i]*radius[i]);
  return 2;
}

/* ----------------------------------------------------------------------
   write atom info to data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecSPH2::write_data(FILE *fp, int n, double **buf)
{
  
  error->all(FLERR,"This feature is not supported by SPH");
  for (int i = 0; i < n; i++)
    fprintf(fp,"%d %d %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %d %d %d\n",
            (int) buf[i][0],(int) buf[i][1],
            buf[i][2],buf[i][3],
            buf[i][4],buf[i][5],buf[i][6],
            (int) buf[i][7],(int) buf[i][8],(int) buf[i][9]);
}

/* ----------------------------------------------------------------------
   write hybrid atom info to data file
------------------------------------------------------------------------- */

int AtomVecSPH2::write_data_hybrid(FILE *fp, double *buf)
{
  
  error->all(FLERR,"This feature is not supported by SPH");
  fprintf(fp," %-1.16e %-1.16e",buf[0],buf[1]);
  return 2;
}

/* ----------------------------------------------------------------------
   pack velocity info for data file
------------------------------------------------------------------------- */

void AtomVecSPH2::pack_vel(double **buf)
{
  
  error->all(FLERR,"This feature is not supported by SPH");
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = tag[i];
    buf[i][1] = v[i][0];
    buf[i][2] = v[i][1];
    buf[i][3] = v[i][2];
  }
}

/* ----------------------------------------------------------------------
   pack hybrid velocity info for data file
------------------------------------------------------------------------- */

int AtomVecSPH2::pack_vel_hybrid(int i, double *buf)
{
  
  error->all(FLERR,"This feature is not supported by SPH");
  return 0;
}

/* ----------------------------------------------------------------------
   write velocity info to data file
------------------------------------------------------------------------- */

void AtomVecSPH2::write_vel(FILE *fp, int n, double **buf)
{
  
  error->all(FLERR,"This feature is not supported by SPH");
  for (int i = 0; i < n; i++)
    fprintf(fp,"%d %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e\n",
            (int) buf[i][0],buf[i][1],buf[i][2],buf[i][3],
            buf[i][4],buf[i][5],buf[i][6]);
}

/* ----------------------------------------------------------------------
   write hybrid velocity info to data file
------------------------------------------------------------------------- */

int AtomVecSPH2::write_vel_hybrid(FILE *fp, double *buf)
{
  
  error->all(FLERR,"This feature is not supported by SPH");
  fprintf(fp," %-1.16e %-1.16e %-1.16e",buf[0],buf[1],buf[2]);
  return 3;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint AtomVecSPH2::memory_usage()
{
  bigint bytes = 0.0;

  if (atom->memcheck("tag")) bytes += memory->usage(tag,nmax);
  if (atom->memcheck("type")) bytes += memory->usage(type,nmax);
  if (atom->memcheck("mask")) bytes += memory->usage(mask,nmax);
  if (atom->memcheck("image")) bytes += memory->usage(image,nmax);
  if (atom->memcheck("x")) bytes += memory->usage(x,nmax,3);
  if (atom->memcheck("v")) bytes += memory->usage(v,nmax,3);
  if (atom->memcheck("f")) bytes += memory->usage(f,nmax*comm->nthreads,3);

  if (atom->memcheck("p")) bytes += memory->usage(p,nmax);
  if (atom->memcheck("rho")) bytes += memory->usage(rho,nmax);
  if (atom->memcheck("drho")) bytes += memory->usage(drho,nmax);
  if (atom->memcheck("e")) bytes += memory->usage(e,nmax);
  if (atom->memcheck("de")) bytes += memory->usage(de,nmax);
  if (atom->memcheck("radius")) bytes += memory->usage(radius,nmax);
  if (atom->memcheck("rmass")) bytes += memory->usage(rmass,nmax);

  return bytes;
}
