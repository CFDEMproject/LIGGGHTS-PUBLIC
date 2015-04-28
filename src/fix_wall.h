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

#ifndef LMP_FIX_WALL_H
#define LMP_FIX_WALL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixWall : public Fix {
 public:
  int nwall;
  int wallwhich[6];
  double coord0[6];
  int xflag;           // 1 if any wall position is a variable
  int xstyle[6];
  int xindex[6];
  char *xstr[6];

  FixWall(class LAMMPS *, int, char **);
  virtual ~FixWall();
  int setmask();
  virtual void init();
  void setup(int);
  void min_setup(int);
  void pre_force(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();
  double compute_vector(int);

  virtual void precompute(int) = 0;
  virtual void wall_particle(int, int, double) = 0;

 protected:
  double epsilon[6],sigma[6],cutoff[6];
  double ewall[7],ewall_all[7];
  double xscale,yscale,zscale;
  int estyle[6],sstyle[6],wstyle[6];
  int eindex[6],sindex[6];
  char *estr[6],*sstr[6];
  int varflag;                // 1 if any wall position,epsilon,sigma is a var
  int eflag;                  // per-wall flag for energy summation
  int nlevels_respa;
  double dt;
  int fldflag;
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Wall defined twice in fix wall command

Self-explanatory.

E: Fix wall cutoff <= 0.0

Self-explanatory.

E: Cannot use fix wall zlo/zhi for a 2d simulation

Self-explanatory.

E: Cannot use fix wall in periodic dimension

Self-explanatory.

E: Variable name for fix wall does not exist

Self-explanatory.

E: Variable for fix wall is invalid style

Only equal-style variables can be used.

E: Variable evaluation in fix wall gave bad value

The returned value for epsilon or sigma < 0.0.

U: Use of fix wall with undefined lattice

Must use lattice command with fix wall command if units option is set
to lattice.

*/
