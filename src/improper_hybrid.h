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

#ifdef IMPROPER_CLASS

ImproperStyle(hybrid,ImproperHybrid)

#else

#ifndef LMP_IMPROPER_HYBRID_H
#define LMP_IMPROPER_HYBRID_H

#include "stdio.h"
#include "improper.h"

namespace LAMMPS_NS {

class ImproperHybrid : public Improper {
 public:
  int nstyles;                  // # of different improper styles
  Improper **styles;            // class list for each Improper style
  char **keywords;              // keyword for each improper style

  ImproperHybrid(class LAMMPS *);
  ~ImproperHybrid();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double memory_usage();

 private:
  int *map;                     // which style each improper type points to

  int *nimproperlist;           // # of impropers in sub-style improperlists
  int *maximproper;             // max # of impropers sub-style lists can store
  int ***improperlist;          // improperlist for each sub-style

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

E: Improper style hybrid cannot use same improper style twice

Self-explanatory.

E: Improper style hybrid cannot have hybrid as an argument

Self-explanatory.

E: Improper style hybrid cannot have none as an argument

Self-explanatory.

E: Improper coeff for hybrid has invalid style

Improper style hybrid uses another improper style as one of its
coefficients.  The improper style used in the improper_coeff command
or read from a restart file is not recognized.

*/
