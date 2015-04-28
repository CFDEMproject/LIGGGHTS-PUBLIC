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

#ifndef LMP_DIHEDRAL_H
#define LMP_DIHEDRAL_H

#include "stdio.h"
#include "pointers.h"

namespace LAMMPS_NS {

class Dihedral : protected Pointers {
  friend class ThrOMP;
  friend class FixOMP;
 public:
  int allocated;
  int *setflag;
  int writedata;                     // 1 if writes coeffs to data file
  double energy;                     // accumulated energy
  double virial[6];                  // accumlated virial
  double *eatom,**vatom;             // accumulated per-atom energy/virial
  unsigned int datamask;
  unsigned int datamask_ext;

  Dihedral(class LAMMPS *);
  virtual ~Dihedral();
  virtual void init();
  virtual void init_style() {}
  virtual void compute(int, int) = 0;
  virtual void settings(int, char **) {}
  virtual void coeff(int, char **) = 0;
  virtual void write_restart(FILE *) = 0;
  virtual void read_restart(FILE *) = 0;
  virtual void write_data(FILE *) {}
  virtual double memory_usage();

  virtual unsigned int data_mask() {return datamask;}
  virtual unsigned int data_mask_ext() {return datamask_ext;}

 protected:
  int suffix_flag;             // suffix compatibility flag

  int evflag;
  int eflag_either,eflag_global,eflag_atom;
  int vflag_either,vflag_global,vflag_atom;
  int maxeatom,maxvatom;

  void ev_setup(int, int);
  void ev_tally(int, int, int, int, int, int, double,
                double *, double *, double *, double, double, double,
                double, double, double, double, double, double);
};

}

#endif

/* ERROR/WARNING messages:

E: Dihedral coeffs are not set

No dihedral coefficients have been assigned in the data file or via
the dihedral_coeff command.

E: All dihedral coeffs are not set

All dihedral coefficients must be set in the data file or by the
dihedral_coeff command before running a simulation.

*/
