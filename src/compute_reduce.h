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
    This file is from LAMMPS
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
Contributing authors for SPH:
Andreas Eitzlmayr (Institute for Process and Particle Engineering, TU Graz)
andreas.eitzlmayr@tugraz.at
Christoph Kloss (DCS)
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(reduce,ComputeReduce)
ComputeStyle(reduce/sph,ComputeReduce) 

#else

#ifndef LMP_COMPUTE_REDUCE_H
#define LMP_COMPUTE_REDUCE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeReduce : public Compute {
 public:
  ComputeReduce(class LAMMPS *, int, char **);
  virtual ~ComputeReduce();
  void init();
  double compute_scalar();
  void compute_vector();
  double memory_usage();

 protected:
  int me;
  int mode,nvalues,iregion;
  int *which,*argindex,*flavor,*value2index;
  char **ids;
  double *onevec;
  int *replace,*indices,*owner;
  int index;
  char *idregion;

  int maxatom;
  double *varatom;

  struct Pair {
    double value;
    int proc;
  };
  Pair pairme,pairall;

  virtual double compute_one(int, int);
  virtual bigint count(int);
  void combine(double &, double, int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region ID for compute reduce/region does not exist

Self-explanatory.

E: Compute reduce replace requires min or max mode

Self-explanatory.

E: Invalid replace values in compute reduce

Self-explanatory.

E: Compute ID for compute reduce does not exist

Self-explanatory.

E: Compute reduce compute does not calculate a per-atom vector

Self-explanatory.

E: Compute reduce compute does not calculate a per-atom array

Self-explanatory.

E: Compute reduce compute array is accessed out-of-range

An index for the array is out of bounds.

E: Compute reduce compute does not calculate a local vector

Self-explanatory.

E: Compute reduce compute does not calculate a local array

Self-explanatory.

E: Compute reduce compute calculates global values

A compute that calculates peratom or local values is required.

E: Fix ID for compute reduce does not exist

Self-explanatory.

E: Compute reduce fix does not calculate a per-atom vector

Self-explanatory.

E: Compute reduce fix does not calculate a per-atom array

Self-explanatory.

E: Compute reduce fix array is accessed out-of-range

An index for the array is out of bounds.

E: Compute reduce fix does not calculate a local vector

Self-explanatory.

E: Compute reduce fix does not calculate a local array

Self-explanatory.

E: Compute reduce fix calculates global values

A fix that calculates peratom or local values is required.

E: Variable name for compute reduce does not exist

Self-explanatory.

E: Compute reduce variable is not atom-style variable

Self-explanatory.

E: Fix used in compute reduce not computed at compatible time

Fixes generate their values on specific timesteps.  Compute reduce is
requesting a value on a non-allowed timestep.

*/
