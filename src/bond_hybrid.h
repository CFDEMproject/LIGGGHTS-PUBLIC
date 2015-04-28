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

#ifdef BOND_CLASS

BondStyle(hybrid,BondHybrid)

#else

#ifndef LMP_BOND_HYBRID_H
#define LMP_BOND_HYBRID_H

#include "stdio.h"
#include "bond.h"

namespace LAMMPS_NS {

class BondHybrid : public Bond {
  friend class Force;

 public:
  int nstyles;                  // # of different bond styles
  Bond **styles;                // class list for each Bond style
  char **keywords;              // keyword for each Bond style

  BondHybrid(class LAMMPS *);
  ~BondHybrid();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double equilibrium_distance(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double single(int, double, int, int, double &);
  double memory_usage();

 private:
  int *map;                     // which style each bond type points to

  int *nbondlist;               // # of bonds in sub-style bondlists
  int *maxbond;                 // max # of bonds sub-style lists can store
  int ***bondlist;              // bondlist for each sub-style

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Bond style hybrid cannot use same bond style twice

Self-explanatory.

E: Bond style hybrid cannot have hybrid as an argument

Self-explanatory.

E: Bond style hybrid cannot have none as an argument

Self-explanatory.

E: Bond coeff for hybrid has invalid style

Bond style hybrid uses another bond style as one of its coefficients.
The bond style used in the bond_coeff command or read from a restart
file is not recognized.

E: Invoked bond equil distance on bond style none

Self-explanatory.

E: Invoked bond single on bond style none

Self-explanatory.

*/
