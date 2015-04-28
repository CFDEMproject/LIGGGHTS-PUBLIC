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

#ifdef PAIR_CLASS

PairStyle(hybrid,PairHybrid)

#else

#ifndef LMP_PAIR_HYBRID_H
#define LMP_PAIR_HYBRID_H

#include "stdio.h"
#include "pair.h"

namespace LAMMPS_NS {

class PairHybrid : public Pair {
 public:
  int nstyles;                  // # of sub-styles
  Pair **styles;                // list of Pair style classes
  char **keywords;              // style name of each Pair style
  int *multiple;                // 0 if style used once, else Mth instance

  PairHybrid(class LAMMPS *);
  virtual ~PairHybrid();
  void compute(int, int);
  void settings(int, char **);
  virtual void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double single(int, int, int, int, double, double, double, double &);
  void modify_params(int narg, char **arg);
  double memory_usage();

  void compute_inner();
  void compute_middle();
  void compute_outer(int, int);
  void *extract(const char *, int &);
  void reset_dt();

  int check_ijtype(int, int, char *);

 protected:
  int outerflag;                // toggle compute() when invoked by outer()

  int **nmap;                   // # of sub-styles itype,jtype points to
  int ***map;                   // list of sub-styles itype,jtype points to

  void allocate();
  void flags();
  virtual void modify_requests();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Pair style hybrid cannot have hybrid as an argument

Self-explanatory.

E: Pair style hybrid cannot have none as an argument

Self-explanatory.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair coeff for hybrid has invalid style

Style in pair coeff must have been listed in pair_style command.

E: Pair hybrid sub-style is not used

No pair_coeff command used a sub-style specified in the pair_style
command.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Invoked pair single on pair style none

A command (e.g. a dump) attempted to invoke the single() function on a
pair style none, which is illegal.  You are probably attempting to
compute per-atom quantities with an undefined pair style.

E: Pair hybrid sub-style does not support single call

You are attempting to invoke a single() call on a pair style
that doesn't support it.

E: Coulomb cutoffs of pair hybrid sub-styles do not match

If using a Kspace solver, all Coulomb cutoffs of long pair styles must
be the same.

*/
