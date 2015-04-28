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

#ifdef FIX_CLASS

FixStyle(ave/spatial,FixAveSpatial)

#else

#ifndef LMP_FIX_AVE_SPATIAL_H
#define LMP_FIX_AVE_SPATIAL_H

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class FixAveSpatial : public Fix {
 public:
  FixAveSpatial(class LAMMPS *, int, char **);
  ~FixAveSpatial();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();
  double compute_array(int,int);
  double memory_usage();
  void reset_timestep(bigint);

 private:
  int me,nvalues;
  int nrepeat,nfreq,irepeat;
  bigint nvalid;
  int ndim,normflag,regionflag,iregion,overwrite;
  char *tstring,*sstring,*idregion;
  int *which,*argindex,*value2index;
  char **ids;
  FILE *fp;
  FILE *fp2; 
  bool write_ts_;
  class Region *region;

  int ave,nwindow,scaleflag;
  int norm,iwindow,window_limit;
  double xscale,yscale,zscale;
  double bin_volume;

  long filepos;
  int dim[3],originflag[3],nlayers[3];
  double origin[3],delta[3];
  double offset[3],invdelta[3];

  int maxvar;
  double *varatom;

  int maxatom;
  int *bin;

  int nbins,maxbin;
  double **coord;
  double *count_one,*count_many,*count_sum;
  double **values_one,**values_many,**values_sum;
  double *count_total,**count_list;
  double **values_total,***values_list;

  int calcStd;
  int lowerLimit, upperLimit;

  void setup_bins();
  void atom2bin1d();
  void atom2bin2d();
  void atom2bin3d();
  bigint nextvalid();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot use fix ave/spatial z for 2 dimensional model

Self-explanatory.

E: Same dimension twice in fix ave/spatial

Self-explanatory.

E: Region ID for fix ave/spatial does not exist

Self-explanatory.

E: Cannot open fix ave/spatial file %s

The specified file cannot be opened.  Check that the path and name are
correct.

E: Compute ID for fix ave/spatial does not exist

Self-explanatory.

E: Fix ave/spatial compute does not calculate per-atom values

A compute used by fix ave/spatial must generate per-atom values.

E: Fix ave/spatial compute does not calculate a per-atom vector

A compute used by fix ave/spatial must generate per-atom values.

E: Fix ave/spatial compute does not calculate a per-atom array

Self-explanatory.

E: Fix ave/spatial compute vector is accessed out-of-range

The index for the vector is out of bounds.

E: Fix ID for fix ave/spatial does not exist

Self-explanatory.

E: Fix ave/spatial fix does not calculate per-atom values

A fix used by fix ave/spatial must generate per-atom values.

E: Fix ave/spatial fix does not calculate a per-atom vector

A fix used by fix ave/spatial must generate per-atom values.

E: Fix ave/spatial fix does not calculate a per-atom array

Self-explanatory.

E: Fix ave/spatial fix vector is accessed out-of-range

The index for the vector is out of bounds.

E: Variable name for fix ave/spatial does not exist

Self-explanatory.

E: Fix ave/spatial variable is not atom-style variable

A variable used by fix ave/spatial must generate per-atom values.

E: Fix ave/spatial for triclinic boxes requires units reduced

Self-explanatory.

E: Fix ave/spatial settings invalid with changing box

If the ave setting is "running" or "window" and the box size/shape
changes during the simulation, then the units setting must be
"reduced", else the number of bins may change.

E: Fix for fix ave/spatial not computed at compatible time

Fixes generate their values on specific timesteps.  Fix ave/spatial is
requesting a value on a non-allowed timestep.

E: Fix ave/spatial missed timestep

You cannot reset the timestep to a value beyond where the fix
expects to next perform averaging.

U: Use of fix ave/spatial with undefined lattice

A lattice must be defined to use fix ave/spatial with units = lattice.

*/
