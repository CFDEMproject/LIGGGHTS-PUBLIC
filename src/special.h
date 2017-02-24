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

#ifndef LMP_SPECIAL_H
#define LMP_SPECIAL_H

#include "pointers.h"

namespace LAMMPS_NS {

class Special : protected Pointers {
 public:
  Special(class LAMMPS *);
  ~Special();
  void build();

 private:
  int me,nprocs;
  int **onetwo,**onethree,**onefour;
  int dihedral_flag;

  // data used by ring callback methods

  int *count;
  int **dflag;

  void dedup();
  void angle_trim();
  void dihedral_trim();
  void combine();

  // static variable for ring communication callback to access class data
  // callback functions for ring communication

  static Special *sptr;
  static void ring_one(int, char *);
  static void ring_two(int, char *);
  static void ring_three(int, char *);
  static void ring_four(int, char *);
  static void ring_five(int, char *);
  static void ring_six(int, char *);
  static void ring_seven(int, char *);
  static void ring_eight(int, char *);
};

}

#endif

/* ERROR/WARNING messages:

E: 1-3 bond count is inconsistent

An inconsistency was detected when computing the number of 1-3
neighbors for each atom.  This likely means something is wrong with
the bond topologies you have defined.

E: 1-4 bond count is inconsistent

An inconsistency was detected when computing the number of 1-4
neighbors for each atom.  This likely means something is wrong with
the bond topologies you have defined.

*/
