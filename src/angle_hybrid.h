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

#ifdef ANGLE_CLASS

AngleStyle(hybrid,AngleHybrid)

#else

#ifndef LMP_ANGLE_HYBRID_H
#define LMP_ANGLE_HYBRID_H

#include "stdio.h"
#include "angle.h"

namespace LAMMPS_NS {

class AngleHybrid : public Angle {
 public:
  int nstyles;                  // # of different angle styles
  Angle **styles;               // class list for each Angle style
  char **keywords;              // keyword for each Angle style

  AngleHybrid(class LAMMPS *);
  ~AngleHybrid();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double equilibrium_angle(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double single(int, int, int, int);
  double memory_usage();

 private:
  int *map;                     // which style each angle type points to

  int *nanglelist;              // # of angles in sub-style anglelists
  int *maxangle;                // max # of angles sub-style lists can store
  int ***anglelist;             // anglelist for each sub-style

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

E: Angle style hybrid cannot use same angle style twice

Self-explanatory.

E: Angle style hybrid cannot have hybrid as an argument

Self-explanatory.

E: Angle style hybrid cannot have none as an argument

Self-explanatory.

E: Angle coeff for hybrid has invalid style

Angle style hybrid uses another angle style as one of its
coefficients.  The angle style used in the angle_coeff command or read
from a restart file is not recognized.

E: Invoked angle equil angle on angle style none

Self-explanatory.

E: Invoked angle single on angle style none

Self-explanatory.

*/
