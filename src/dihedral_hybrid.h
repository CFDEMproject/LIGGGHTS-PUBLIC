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

#ifdef DIHEDRAL_CLASS

DihedralStyle(hybrid,DihedralHybrid)

#else

#ifndef LMP_DIHEDRAL_HYBRID_H
#define LMP_DIHEDRAL_HYBRID_H

#include "stdio.h"
#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralHybrid : public Dihedral {
 public:
  int nstyles;                  // # of different dihedral styles
  Dihedral **styles;            // class list for each Dihedral style
  char **keywords;              // keyword for each dihedral style

  DihedralHybrid(class LAMMPS *);
  ~DihedralHybrid();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  void write_restart(FILE *);
  void read_restart(FILE *);
  double memory_usage();

 private:
  int *map;                     // which style each dihedral type points to

  int *ndihedrallist;           // # of dihedrals in sub-style dihedrallists
  int *maxdihedral;             // max # of dihedrals sub-style lists can store
  int ***dihedrallist;          // dihedrallist for each sub-style

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

E: Dihedral style hybrid cannot use same dihedral style twice

Self-explanatory.

E: Dihedral style hybrid cannot have hybrid as an argument

Self-explanatory.

E: Dihedral style hybrid cannot have none as an argument

Self-explanatory.

E: Dihedral coeff for hybrid has invalid style

Dihedral style hybrid uses another dihedral style as one of its
coefficients.  The dihedral style used in the dihedral_coeff command
or read from a restart file is not recognized.

*/
