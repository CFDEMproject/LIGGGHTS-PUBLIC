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

#ifndef LMP_ERROR_H
#define LMP_ERROR_H

#include "pointers.h"

namespace LAMMPS_NS {

class Error : protected Pointers {
 public:
  Error(class LAMMPS *);

  void universe_all(const char *, int, const char *);
  void universe_one(const char *, int, const char *);
  void universe_warn(const char *, int, const char *);

  void all(const char *, int, const char *);
  void fix_error(const char *, int, class Fix*,const char *); 
  void compute_error(const char *, int, class Compute*,const char *); 
  void cg(const char *, int, const char *); 
  void one(const char *, int, const char *);
  void warning(const char *, int, const char *, int = 1);
  void message(const char *, int, const char *, int = 1);
  void done();
};

}

#endif

/* ERROR/WARNING messages:

*/
