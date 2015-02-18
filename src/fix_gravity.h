/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   This file was modified with respect to the release in LAMMPS
   Modifications are Copyright 2009-2012 JKU Linz
                     Copyright 2012-     DCS Computing GmbH, Linz

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(gravity,FixGravity)

#else

#ifndef LMP_FIX_GRAVITY_H
#define LMP_FIX_GRAVITY_H

#include "fix.h"

namespace LAMMPS_NS {

class FixGravity : public Fix {
  friend class FixPour;

 public:
  FixGravity(class LAMMPS *, int, char **);
  ~FixGravity();
  int setmask();
  void init();
  void setup(int);
  virtual void post_force(int);
  virtual void post_force_respa(int, int, int);
  double compute_scalar();

  void get_gravity(double*); 

 protected:
  int style;
  double magnitude;
  double vert,phi,theta;
  double xdir,ydir,zdir;
  double xgrav,ygrav,zgrav,xacc,yacc,zacc;
  double degree2rad;
  int nlevels_respa;
  int time_origin;
  int eflag;
  double egrav,egrav_all;

  int varflag;
  int mstyle,vstyle,pstyle,tstyle,xstyle,ystyle,zstyle;
  int mvar,vvar,pvar,tvar,xvar,yvar,zvar;
  char *mstr,*vstr,*pstr,*tstr,*xstr,*ystr,*zstr;

  void set_acceleration();
  class FixMultisphere *fm;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Variable name for fix gravity does not exist

Self-explanatory.

E: Variable for fix gravity is invalid style

Only equal-style variables can be used.

*/
