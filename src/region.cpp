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
    This file is from LAMMPS, but has been modified. Copyright for
    modification:

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz

    Copyright of original file:
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

#include <cmath>
#include <algorithm>
#include <stdlib.h>
#include <string.h>
#include "region.h"
#include "update.h"
#include "domain.h"
#include "lattice.h"
#include "input.h"
#include "variable.h"
#include "error.h"
#include "force.h"
#include "random_park.h"
#include "vector_liggghts.h" 
#include "mpi_liggghts.h"  
#include "math_extra_liggghts.h" 
#include "comm.h"

#define SMALL 1e-8

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Region::Region(LAMMPS *lmp, int narg, char **arg) : Pointers(lmp)
{
  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  n = strlen(arg[1]) + 1;
  style = new char[n];
  strcpy(style,arg[1]);

  varshape = 0;
  xstr = ystr = zstr = tstr = NULL;
  dx = dy = dz = 0.0;
  lastshape = lastdynamic = -1;

  random = NULL;

  volume_limit_ = 1.e-10;
}

/* ---------------------------------------------------------------------- */

Region::~Region()
{
  delete [] id;
  delete [] style;

  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  delete [] tstr;

  if (random) delete random;
}

/* ---------------------------------------------------------------------- */

void Region::init()
{

  seed = 3012211;

  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0) error->all(FLERR,"Variable name for region does not exist");
    if (!input->variable->equalstyle(xvar))
      error->all(FLERR,"Variable for region is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0) error->all(FLERR,"Variable name for region does not exist");
    if (!input->variable->equalstyle(yvar))
      error->all(FLERR,"Variable for region is not equal style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0) error->all(FLERR,"Variable name for region does not exist");
    if (!input->variable->equalstyle(zvar))
      error->all(FLERR,"Variable for region is not equal style");
  }
  if (tstr) {
    tvar = input->variable->find(tstr);
    if (tvar < 0) error->all(FLERR,"Variable name for region does not exist");
    if (!input->variable->equalstyle(tvar))
      error->all(FLERR,"Variable for region is not equal style");
  }
}

/* ----------------------------------------------------------------------
   return 1 if region is dynamic, 0 if static
   only primitive regions define it here
   union/intersect regions have their own dynamic_check()
------------------------------------------------------------------------- */

int Region::dynamic_check()
{
  return dynamic;
}

/* ----------------------------------------------------------------------
   called before looping over atoms with match() or surface()
   this insures any variables used by region are invoked once per timestep
     also insures variables are invoked by all procs even those w/out atoms
     necessary if equal-style variable invokes global operation
   with MPI_Allreduce, e.g. xcm() or count()
------------------------------------------------------------------------- */

void Region::prematch()
{
  if (varshape) shape_update();
  if (dynamic) pretransform();
}

/* ----------------------------------------------------------------------
   determine if point x,y,z is a match to region volume
   XOR computes 0 if 2 args are the same, 1 if different
   note that inside() returns 1 for points on surface of region
   thus point on surface of exterior region will not match
   if region has variable shape, invoke shape_update() once per timestep
   if region is dynamic, apply inverse transform to x,y,z
     unmove first, then unrotate, so don't have to change rotation point
   caller is responsible for wrapping this call with
     modify->clearstep_compute() and modify->addstep_compute() if needed
------------------------------------------------------------------------- */

int Region::match(double x, double y, double z)
{
  if (varshape && update->ntimestep != lastshape) {
    shape_update();
    lastshape = update->ntimestep;
  }
  if (dynamic) inverse_transform(x,y,z);
  return !(inside(x,y,z) ^ interior);
}

/* ----------------------------------------------------------------------
   generate list of contact points for interior or exterior regions
   if region has variable shape, invoke shape_update() once per timestep
   if region is dynamic:
     before: inverse transform x,y,z (unmove, then unrotate)
     after: forward transform contact point xs,yx,zs (rotate, then move),
            then reset contact delx,dely,delz based on new contact point
            no need to do this if no rotation since delxyz doesn't change
   caller is responsible for wrapping this call with
     modify->clearstep_compute() and modify->addstep_compute() if needed
------------------------------------------------------------------------- */

int Region::surface(double x, double y, double z, double cutoff)
{
  int ncontact;
  double xs,ys,zs;
  double xnear[3],xorig[3]={};

  if (varshape && update->ntimestep != lastshape) {
    shape_update();
    lastshape = update->ntimestep;
  }
  if (dynamic) {
    xorig[0] = x;
    xorig[1] = y;
    xorig[2] = z;
    inverse_transform(x,y,z);
  }

  xnear[0] = x;
  xnear[1] = y;
  xnear[2] = z;

  if (interior) ncontact = surface_interior(xnear,cutoff);
  else ncontact = surface_exterior(xnear,cutoff);

  if (rotateflag && ncontact) {
    for (int i = 0; i < ncontact; i++) {
      xs = xnear[0] - contact[i].delx;
      ys = xnear[1] - contact[i].dely;
      zs = xnear[2] - contact[i].delz;
      forward_transform(xs,ys,zs);
      contact[i].delx = xorig[0] - xs;
      contact[i].dely = xorig[1] - ys;
      contact[i].delz = xorig[2] - zs;
    }
  }

  return ncontact;
}

/* ----------------------------------------------------------------------
   add a single contact at Nth location in contact array
   x = particle position
   xp,yp,zp = region surface point
------------------------------------------------------------------------- */

void Region::add_contact(int n, double *x, double xp, double yp, double zp)
{
  double delx = x[0] - xp;
  double dely = x[1] - yp;
  double delz = x[2] - zp;
  contact[n].r = sqrt(delx*delx + dely*dely + delz*delz);
  contact[n].delx = delx;
  contact[n].dely = dely;
  contact[n].delz = delz;
}

/* ----------------------------------------------------------------------
   pre-compute dx,dy,dz and theta for a moving/rotating region
   called once for the region before per-atom loop, via prematch()
------------------------------------------------------------------------- */

void Region::pretransform()
{
  if (moveflag) {
    if (xstr) dx = input->variable->compute_equal(xvar);
    if (ystr) dy = input->variable->compute_equal(yvar);
    if (zstr) dz = input->variable->compute_equal(zvar);
  }
  if (rotateflag) theta = input->variable->compute_equal(tvar);
}

/* ----------------------------------------------------------------------
   transform a point x,y,z in region space to moved space
   rotate first (around original P), then displace
------------------------------------------------------------------------- */

void Region::forward_transform(double &x, double &y, double &z)
{
  if (rotateflag) {
    if (update->ntimestep != lastdynamic)
      theta = input->variable->compute_equal(tvar);
    rotate(x,y,z,theta);
  }

  if (moveflag) {
    if (update->ntimestep != lastdynamic) {
      if (xstr) dx = input->variable->compute_equal(xvar);
      if (ystr) dy = input->variable->compute_equal(yvar);
      if (zstr) dz = input->variable->compute_equal(zvar);
    }
    x += dx;
    y += dy;
    z += dz;
  }

  lastdynamic = update->ntimestep;
}

/* ----------------------------------------------------------------------
   transform a point x,y,z in moved space back to region space
   undisplace first, then unrotate (around original P)
------------------------------------------------------------------------- */

void Region::inverse_transform(double &x, double &y, double &z)
{
  if (moveflag) {
    if (update->ntimestep != lastdynamic) {
      if (xstr) dx = input->variable->compute_equal(xvar);
      if (ystr) dy = input->variable->compute_equal(yvar);
      if (zstr) dz = input->variable->compute_equal(zvar);
    }
    x -= dx;
    y -= dy;
    z -= dz;
  }

  if (rotateflag) {
    if (update->ntimestep != lastdynamic)
      theta = input->variable->compute_equal(tvar);
    rotate(x,y,z,-theta);
  }

  lastdynamic = update->ntimestep;
}

/* ----------------------------------------------------------------------
   rotate x,y,z by angle via right-hand rule around point and runit normal
   sign of angle determines whether rotating forward/backward in time
   return updated x,y,z
   R = vector axis of rotation
   P = point = point to rotate around
   R0 = runit = unit vector for R
   X0 = x,y,z = initial coord of atom
   D = X0 - P = vector from P to X0
   C = (D dot R0) R0 = projection of D onto R, i.e. Dparallel
   A = D - C = vector from R line to X0, i.e. Dperp
   B = R0 cross A = vector perp to A in plane of rotation, same len as A
   A,B define plane of circular rotation around R line
   new x,y,z = P + C + A cos(angle) + B sin(angle)
------------------------------------------------------------------------- */

void Region::rotate(double &x, double &y, double &z, const double angle)
{
  double a[3],b[3],c[3],d[3],disp[3];

  double sine = sin(angle);
  double cosine = cos(angle);
  d[0] = x - point[0];
  d[1] = y - point[1];
  d[2] = z - point[2];
  double x0dotr = d[0]*runit[0] + d[1]*runit[1] + d[2]*runit[2];
  c[0] = x0dotr * runit[0];
  c[1] = x0dotr * runit[1];
  c[2] = x0dotr * runit[2];
  a[0] = d[0] - c[0];
  a[1] = d[1] - c[1];
  a[2] = d[2] - c[2];
  b[0] = runit[1]*a[2] - runit[2]*a[1];
  b[1] = runit[2]*a[0] - runit[0]*a[2];
  b[2] = runit[0]*a[1] - runit[1]*a[0];
  disp[0] = a[0]*cosine  + b[0]*sine;
  disp[1] = a[1]*cosine  + b[1]*sine;
  disp[2] = a[2]*cosine  + b[2]*sine;
  x = point[0] + c[0] + disp[0];
  y = point[1] + c[1] + disp[1];
  z = point[2] + c[2] + disp[2];

}

/* ----------------------------------------------------------------------
   parse optional parameters at end of region input line
------------------------------------------------------------------------- */

void Region::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal region command");

  // option defaults

  interior = 1;
  scaleflag = 0; 
  moveflag = rotateflag = 0;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal region command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal region command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"side") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal region command");
      if (strcmp(arg[iarg+1],"in") == 0) interior = 1;
      else if (strcmp(arg[iarg+1],"out") == 0) interior = 0;
      else error->all(FLERR,"Illegal region command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"move") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal region command");
      if (strcmp(arg[iarg+1],"NULL") != 0) {
        if (strstr(arg[iarg+1],"v_") != arg[iarg+1])
          error->all(FLERR,"Illegal region command");
        int n = strlen(&arg[iarg+1][2]) + 1;
        xstr = new char[n];
        strcpy(xstr,&arg[iarg+1][2]);
      }
      if (strcmp(arg[iarg+2],"NULL") != 0) {
        if (strstr(arg[iarg+2],"v_") != arg[iarg+2])
          error->all(FLERR,"Illegal region command");
        int n = strlen(&arg[iarg+2][2]) + 1;
        ystr = new char[n];
        strcpy(ystr,&arg[iarg+2][2]);
      }
      if (strcmp(arg[iarg+3],"NULL") != 0) {
        if (strstr(arg[iarg+3],"v_") != arg[iarg+3])
          error->all(FLERR,"Illegal region command");
        int n = strlen(&arg[iarg+3][2]) + 1;
        zstr = new char[n];
        strcpy(zstr,&arg[iarg+3][2]);
      }
      moveflag = 1;
      iarg += 4;

    } else if (strcmp(arg[iarg],"rotate") == 0) {
      if (iarg+8 > narg) error->all(FLERR,"Illegal region command");
      if (strstr(arg[iarg+1],"v_") != arg[iarg+1])
        error->all(FLERR,"Illegal region command");
      int n = strlen(&arg[iarg+1][2]) + 1;
      tstr = new char[n];
      strcpy(tstr,&arg[iarg+1][2]);
      point[0] = force->numeric(FLERR,arg[iarg+2]);
      point[1] = force->numeric(FLERR,arg[iarg+3]);
      point[2] = force->numeric(FLERR,arg[iarg+4]);
      axis[0] = force->numeric(FLERR,arg[iarg+5]);
      axis[1] = force->numeric(FLERR,arg[iarg+6]);
      axis[2] = force->numeric(FLERR,arg[iarg+7]);
      rotateflag = 1;
      iarg += 8;
    
    } else if (strcmp(arg[iarg],"seed") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal region command");
      random = new RanPark(lmp, arg[iarg+1], true);
      iarg += 2;
    } else if (strcmp(arg[iarg],"volume_limit") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal region command");
      volume_limit_ = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    }
     else error->all(FLERR,"Illegal region command");
  }

  if (!random)
      random = new RanPark(lmp, "3012211", true);
  seed = random->getSeed();
  
  // error check

  if ((moveflag || rotateflag) &&
      (strcmp(style,"union") == 0 || strcmp(style,"intersect") == 0))
    error->all(FLERR,"Region union or intersect cannot be dynamic");

  // setup scaling

  if (scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  if (rotateflag) {
    point[0] *= xscale;
    point[1] *= yscale;
    point[2] *= zscale;
  }

  // runit = unit vector along rotation axis

  if (rotateflag) {
    double len = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
    if (len == 0.0)
      error->all(FLERR,"Region cannot have 0 length rotation vector");
    runit[0] = axis[0]/len;
    runit[1] = axis[1]/len;
    runit[2] = axis[2]/len;
  }

  if (moveflag || rotateflag) dynamic = 1;
  else dynamic = 0;
}

/* ---------------------------------------------------------------------- */

void Region::reset_random(int new_seed)
{
    if(comm->me == 0) fprintf(screen,"INFO: Resetting random generator for region %s\n",id);
    random->reset(new_seed + comm->me);
}

/* ---------------------------------------------------------------------- */

inline void Region::rand_bounds(bool subdomain_flag, double *lo, double *hi)
{
    if(!bboxflag) error->one(FLERR,"Impossible to generate random points on region with incomputable bounding box");
    if(subdomain_flag)
    {
        lo[0] = std::max(extent_xlo,domain->sublo[0]);
        lo[1] = std::max(extent_ylo,domain->sublo[1]);
        lo[2] = std::max(extent_zlo,domain->sublo[2]);
        hi[0] = std::min(extent_xhi,domain->subhi[0]);
        hi[1] = std::min(extent_yhi,domain->subhi[1]);
        hi[2] = std::min(extent_zhi,domain->subhi[2]);
        if(lo[0] >= hi[0] || lo[1] >= hi[1] ||lo[2] >= hi[2])
            error->one(FLERR,"Impossible to generate random points on wrong sub-domain");
    }
    else
    {
        vectorConstruct3D(lo,  extent_xlo,extent_ylo,extent_zlo );
        vectorConstruct3D(hi,  extent_xhi,extent_yhi,extent_zhi );
    }
    
}

/* ---------------------------------------------------------------------- */

void Region::generate_random(double *pos,bool subdomain_flag)
{
    double lo[3],hi[3],diff[3];
    rand_bounds(subdomain_flag,lo,hi);
    vectorSubtract3D(hi,lo,diff);
    do
    {
        pos[0] = lo[0] + random->uniform()*diff[0];
        pos[1] = lo[1] + random->uniform()*diff[1];
        pos[2] = lo[2] + random->uniform()*diff[2];
    }
    while(!match(pos[0],pos[1],pos[2]));
    
}

/* ---------------------------------------------------------------------- */

// generates a random point within the region and has a min distance from surface
// i.e. generate random point in region "shrunk" by cut
void Region::generate_random_shrinkby_cut(double *pos,double cut,bool subdomain_flag)
{
    
    double lo[3],hi[3],diff[3];
    rand_bounds(subdomain_flag,lo,hi);
    vectorSubtract3D(hi,lo,diff);

    if((extent_xhi-extent_xlo < 2.*cut) ||
       (extent_yhi-extent_ylo < 2.*cut) ||
       (extent_zhi-extent_zlo < 2.*cut))
        error->one(FLERR,"Impossible to generate random points within region - region too small "
        "(smaller than twice the particle cutoff)");

    do
    {
        pos[0] = lo[0] + random->uniform()*diff[0];
        pos[1] = lo[1] + random->uniform()*diff[1];
        pos[2] = lo[2] + random->uniform()*diff[2];
        
    }
    // pos has to be within region, but not within cut of region surface
    while(!match(pos[0],pos[1],pos[2]) || match_cut(pos,cut));
}

/* ---------------------------------------------------------------------- */

// generate a point inside region OR within a minimum distance from surface
void Region::generate_random_expandby_cut(double *pos,double cut,bool subdomain_flag)
{
    double lo[3],hi[3],diff[3];
    rand_bounds(subdomain_flag,lo,hi);
    vectorSubtract3D(hi,lo,diff);

    do
    {
        pos[0] = lo[0] + random->uniform()*diff[0];
        pos[1] = lo[1] + random->uniform()*diff[1];
        pos[2] = lo[2] + random->uniform()*diff[2];
    }
    // pos has to be within region or within a distance (=cut) to region
    while( !match_expandby_cut(pos,cut) );
}

/* ---------------------------------------------------------------------- */

int Region::match_cut(double *pos,double cut)
{
  double x[3]; 
  vectorCopy3D(pos,x);

  if (dynamic) inverse_transform(x[0],x[1],x[2]);

  if(interior) return surface_interior(x,cut);
  else return surface_exterior(x,cut);
}

/* ---------------------------------------------------------------------- */

int Region::match_expandby_cut(double *pos,double cut)
{
  double x[3]; 
  vectorCopy3D(pos,x);

  if (dynamic) inverse_transform(x[0],x[1],x[2]);

  if(interior) return (match(pos[0],pos[1],pos[2]) || surface_exterior(x,cut));
  else         return (match(pos[0],pos[1],pos[2]) || surface_interior(x,cut));
}

/* ---------------------------------------------------------------------- */

int Region::match_shrinkby_cut(double *pos,double cut)
{
  double x[3]; 
  vectorCopy3D(pos,x);

  if (dynamic) inverse_transform(x[0],x[1],x[2]);

  if(interior) return (match(pos[0],pos[1],pos[2]) && !surface_interior(x,cut) );
  else         return (match(pos[0],pos[1],pos[2]) && !surface_exterior(x,cut) );
}

/* ---------------------------------------------------------------------- */

void Region::volume_mc(int n_test,bool cutflag,double cut,double &vol_global,double &vol_local)
{
    double pos[3],vol_bbox, vol_local_all;
    int n_in_local = 0, n_in_global = 0, n_in_global_all;

    // impossible to calculate volume if bbox non-existent
    if(!bboxflag)
    {
        vol_global = vol_local = 0.;
        error->all(FLERR,"Unable to calculate region volume. Region needs to have existing bounding box");
        return;
    }

    for(int i = 0; i < n_test; i++)
    {
        pos[0] = extent_xlo + random->uniform() * (extent_xhi - extent_xlo);
        pos[1] = extent_ylo + random->uniform() * (extent_yhi - extent_ylo);
        pos[2] = extent_zlo + random->uniform() * (extent_zhi - extent_zlo);

        if(!domain->is_in_domain(pos)) continue;

        // point is in region
        // assume every proc can evaluate this
        
        if(!cutflag)
        {
            if(match(pos[0],pos[1],pos[2]))
            {
                n_in_global++;
                if(domain->is_in_subdomain(pos))
                    n_in_local++;
            }
        }
        else
        {
            if(match(pos[0],pos[1],pos[2]))
            {
                n_in_global++;
                if(domain->is_in_subdomain(pos) && !match_cut(pos,cut) )
                    n_in_local++;
            }
        }
    }

    MPI_Sum_Scalar(n_in_global,n_in_global_all,world);
    if(n_in_global_all == 0)
        error->all(FLERR,"Unable to calculate region volume. Possible sources of error: \n"
                         "   (a) region volume is too small or out of domain (you may want to increase the 'volume_limit' in the input script)\n"
                         "   (b) particles for insertion are too large when using all_in yes\n"
                         "   (c) region is 2d, but should be 3d");

    vol_bbox = (extent_xhi - extent_xlo) * (extent_yhi - extent_ylo) * (extent_zhi - extent_zlo);

    // return calculated values
    vol_global = static_cast<double>(n_in_global_all)/static_cast<double>(n_test*comm->nprocs) * vol_bbox;
    vol_local  = static_cast<double>(n_in_local )/static_cast<double>(n_test) * vol_bbox;

    MPI_Sum_Scalar(vol_local,vol_local_all,world);

    if(vol_local_all < volume_limit_)
        error->all(FLERR,"Unable to calculate region volume. Possible sources of error: \n"
                         "   (a) region volume is too small or out of domain (you may want to decrease the 'volume_limit' in the input script)\n"
                         "   (b) particles for insertion are too large when using all_in yes\n"
                         "   (c) region is 2d, but should be 3d\n");

    vol_local *= (vol_global/vol_local_all);
    
}

/* ---------------------------------------------------------------------- */

int Region::bbox_extends_outside_box()
{
    double min[3],max[3];
    vectorConstruct3D(min,extent_xlo+SMALL,extent_ylo+SMALL,extent_zlo+SMALL);
    vectorConstruct3D(max,extent_xhi-SMALL,extent_yhi-SMALL,extent_zhi-SMALL);
    return (!(domain->is_in_domain(min)) || !(domain->is_in_domain(max)));
}
