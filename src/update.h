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

#ifndef LMP_UPDATE_H
#define LMP_UPDATE_H

#include "pointers.h"

namespace LAMMPS_NS {

class Update : protected Pointers {
 public:
  double dt;                      // timestep
  double etol,ftol;               // minimizer tolerances on energy/force
  bigint ntimestep;               // current step (dynamics or min iterations)
  int nsteps;                     // # of steps to run (dynamics or min iter)
  int whichflag;                  // 0 for unset, 1 for dynamics, 2 for min
  double atime;                   // simulation time at atime_step
  bigint atimestep;               // last timestep atime was updated
  bigint firststep,laststep;      // 1st & last step of this run
  bigint beginstep,endstep;       // 1st and last step of multiple runs
  int first_update;               // 0 before initial update, 1 after
  int max_eval;                   // max force evaluations for minimizer
  int restrict_output;            // 1 if output should not write dump/restart
  int setupflag;                  // set when setup() is computing forces
  int multireplica;               // 1 if min across replicas, else 0

  bigint eflag_global,eflag_atom;  // timestep global/peratom eng is tallied on
  bigint vflag_global,vflag_atom;  // ditto for virial

  char *unit_style;

  class Integrate *integrate;
  char *integrate_style;

  class Min *minimize;
  char *minimize_style;

  Update(class LAMMPS *);
  ~Update();
  void init();
  void set_units(const char *);
  void create_integrate(int, char **, char *);
  void create_minimize(int, char **);
  void reset_timestep(int, char **);
  void reset_timestep(bigint);
  void update_time();
  bigint memory_usage();

 private:
  void new_integrate(char *, int, char **, char *, int &);

};

}

#endif

/* ERROR/WARNING messages:

E: USER-CUDA mode requires CUDA variant of run style

CUDA mode is enabled, so the run style must include a cuda suffix.

E: USER-CUDA mode requires CUDA variant of min style

CUDA mode is enabled, so the min style must include a cuda suffix.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Illegal integrate style

Self-explanatory.

E: Timestep must be >= 0

Specified timestep is invalid.

E: Too big a timestep

Specified timestep is too large.

E: Cannot reset timestep with a time-dependent fix defined

You cannot reset the timestep when a fix that keeps track of elapsed
time is in place.

E: Cannot reset timestep with a dynamic region defined

Dynamic regions (see the region command) have a time dependence.
Thus you cannot change the timestep when one or more of these
are defined.

*/
