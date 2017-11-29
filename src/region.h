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

#ifndef LMP_REGION_H
#define LMP_REGION_H

#include "pointers.h"

namespace LAMMPS_NS {

class Region : protected Pointers {
 public:
  char *id,*style;
  int interior;                     // 1 for interior, 0 for exterior
  int scaleflag;                    // 1 for lattice, 0 for box
  double xscale,yscale,zscale;      // scale factors for box/lattice units
  double extent_xlo,extent_xhi;     // bounding box on region
  double extent_ylo,extent_yhi;
  double extent_zlo,extent_zhi;
  int bboxflag;                     // 1 if bounding box is computable
  int varshape;                     // 1 if region shape changes over time

  // contact = particle near region surface

  struct Contact {
    double r;                 // distance between particle & surf, r > 0.0
    double delx,dely,delz;    // vector from surface pt to particle
  };
  Contact *contact;           // list of contacts
  int cmax;                   // max # of contacts possible with region

  Region(class LAMMPS *, int, char **);
  virtual ~Region();
  virtual void init();
  virtual int dynamic_check();

  // called by other classes to check point versus region

  inline int match(double *point) 
  { return match(point[0],point[1],point[2]);}

  void prematch();
  int match(double, double, double);
  int surface(double, double, double, double);

  // reset random gen - is called out of restart by fix that uses region
  void reset_random(int);

  inline void rand_bounds(bool subdomain_flag, double *lo, double *hi);

  // generates a random point within the region
  virtual void generate_random(double *,bool subdomain_flag);

  // generate a point inside region OR within cut distance from surface
  virtual void generate_random_expandby_cut(double *,double,bool subdomain_flag);

  // generate a point inside region AND further away from surface than cut
  virtual void generate_random_shrinkby_cut(double *,double,bool subdomain_flag);

  // inside region AND within a minimum distance from surface
  int match_cut(double *,double);

  // inside region OR within a minimum distance from surface
  int match_expandby_cut(double *,double);

  // inside region OR within a minimum distance from surface
  int match_shrinkby_cut(double *,double);

  // volume calculation based on MC
  virtual void volume_mc(int n_test,bool cutflag,double cut,double &vol_global,double &vol_local);

  // flag if region bbox extends outside simulation domain
  virtual int bbox_extends_outside_box();

  // implemented by each region, not called by other classes

  virtual int inside(double, double, double) = 0;
  virtual int surface_interior(double *, double) = 0;
  virtual int surface_exterior(double *, double) = 0;
  virtual void shape_update() {}
  virtual void pretransform();

 protected:
  void add_contact(int, double *, double, double, double);
  void options(int, char **);

  int seed;
  class RanPark *random;
  
 private:
  double volume_limit_;             //volume below which error will be thrown
  int dynamic;                      // 1 if region position/orientation changes over time
  int moveflag,rotateflag;          // 1 if position/orientation changes
  double point[3],axis[3],runit[3];
  char *xstr,*ystr,*zstr,*tstr;
  int xvar,yvar,zvar,tvar;
  double dx,dy,dz,theta;
  bigint lastshape,lastdynamic;

  void forward_transform(double &, double &, double &);
  void inverse_transform(double &, double &, double &);
  void rotate(double &, double &, double &, const double);
};

}

#endif

/* ERROR/WARNING messages:

E: Variable name for region does not exist

Self-explanatory.

E: Variable for region is invalid style

Only equal-style variables can be used.

E: Variable for region is not equal style

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region union or intersect cannot be dynamic

The sub-regions can be dynamic, but not the combined region.

E: Region cannot have 0 length rotation vector

Self-explanatory.

U: Use of region with undefined lattice

If units = lattice (the default) for the region command, then a
lattice must first be defined via the lattice command.

*/
