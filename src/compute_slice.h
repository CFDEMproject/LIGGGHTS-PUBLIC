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

#ifdef COMPUTE_CLASS

ComputeStyle(slice,ComputeSlice)

#else

#ifndef LMP_COMPUTE_SLICE_H
#define LMP_COMPUTE_SLICE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeSlice : public Compute {
 public:
  ComputeSlice(class LAMMPS *, int, char **);
  virtual ~ComputeSlice();
  void init();
  void compute_vector();
  void compute_array();

 private:
  int me;
  int nstart,nstop,nskip,nvalues;
  int *which,*argindex,*value2index;
  char **ids;

  void extract_one(int, double *, int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute ID for compute slice does not exist

Self-explanatory.

E: Compute slice compute does not calculate a global array

Self-explanatory.

E: Compute slice compute vector is accessed out-of-range

The index for the vector is out of bounds.

E: Compute slice compute does not calculate a global vector

Self-explanatory.

E: Compute slice compute array is accessed out-of-range

An index for the array is out of bounds.

E: Compute slice compute does not calculate global vector or array

Self-explanatory.

E: Fix ID for compute slice does not exist

Self-explanatory.

E: Compute slice fix does not calculate a global array

Self-explanatory.

E: Compute slice fix vector is accessed out-of-range

The index for the vector is out of bounds.

E: Compute slice fix does not calculate a global vector

Self-explanatory.

E: Compute slice fix array is accessed out-of-range

An index for the array is out of bounds.

E: Compute slice fix does not calculate global vector or array

Self-explanatory.

E: Fix used in compute slice not computed at compatible time

Fixes generate their values on specific timesteps.  Compute slice is
requesting a value on a non-allowed timestep.

*/
