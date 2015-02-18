/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   This file was modified with respect to the release in LAMMPS
   Modifications are Copyright 2009-2012 JKU Linz
                     Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
   the producer of the LIGGGHTS® software and the CFDEM®coupling software
   See http://www.cfdem.com/terms-trademark-policy for details.

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

FixStyle(adapt,FixAdapt)

#else

#ifndef LMP_FIX_ADAPT_H
#define LMP_FIX_ADAPT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAdapt : public Fix {
 public:
  int diamflag;        // 1 if atom diameters will vary, for AtomVecGranular
  int chgflag;

  FixAdapt(class LAMMPS *, int, char **);
  ~FixAdapt();
  int setmask();
  void post_create(); 
  void pre_delete(bool unfixflag); 
  void init();
  void setup_pre_force(int);
  void pre_force(int);
  void post_run();

 private:
  int nadapt,resetflag,scaleflag;
  int anypair;

  struct Adapt {
    int which,ivar;
    char *var;
    char *pstyle,*pparam;
    int ilo,ihi,jlo,jhi;
    int pdim;
    double *scalar,scalar_orig;
    double **array,**array_orig;
    int aparam;
  };

  Adapt *adapt;
  double *kspace_scale;

  void change_settings();
  void restore_settings();

  class FixPropertyAtom *fppat;
  char fixid[100];
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Variable name for fix adapt does not exist

Self-explanatory.

E: Variable for fix adapt is invalid style

Only equal-style variables can be used.

E: Fix adapt pair style does not exist

Self-explanatory

E: Fix adapt pair style param not supported

The pair style does not know about the parameter you specified.

E: Fix adapt type pair range is not valid for pair hybrid sub-style

Self-explanatory.

E: Fix adapt kspace style does not exist

Self-explanatory.

E: Fix adapt requires atom attribute diameter

The atom style being used does not specify an atom diameter.

E: Fix adapt requires atom attribute charge

The atom style being used does not specify an atom charge.

*/
