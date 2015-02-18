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

#ifdef COMPUTE_CLASS

ComputeStyle(stress/atom,ComputeStressAtom)

#else

#ifndef LMP_COMPUTE_STRESS_ATOM_H
#define LMP_COMPUTE_STRESS_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeStressAtom : public Compute {
 public:
  ComputeStressAtom(class LAMMPS *, int, char **);
  ~ComputeStressAtom();
  void init() {}
  void compute_peratom();
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();

  int pack_comm(int n, int *list, double *buf, int pbc_flag, int *pbc);
  void unpack_comm(int n, int first, double *buf);

 private:
  int keflag,pairflag,bondflag,angleflag,dihedralflag,improperflag;
  int kspaceflag,fixflag;
  int nmax;
  double **stress;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Per-atom virial was not tallied on needed timestep

You are using a thermo keyword that requires potentials to have
tallied the virial, but they didn't on this timestep.  See the
variable doc page for ideas on how to make this work.

*/
