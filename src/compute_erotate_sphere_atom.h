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

#ifdef COMPUTE_CLASS

ComputeStyle(erotate/sphere/atom,ComputeErotateSphereAtom)

#else

#ifndef LMP_COMPUTE_EROTATE_SPHERE_ATOM_H
#define LMP_COMPUTE_EROTATE_SPHERE_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeErotateSphereAtom : public Compute {
 public:
  ComputeErotateSphereAtom(class LAMMPS *, int, char **);
  ~ComputeErotateSphereAtom();
  void init();
  void compute_peratom();
  double memory_usage();

 private:
  int nmax;
  double pfactor;
  double *erot;
  class FixMultisphere *fix_ms; 
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute erotate/sphere/atom requires atom style sphere

Self-explanatory.

W: More than one compute erotate/sphere/atom

It is not efficient to use compute erorate/sphere/atom more than once.

*/
