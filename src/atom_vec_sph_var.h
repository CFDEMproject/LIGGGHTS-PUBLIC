/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
   the producer of the LIGGGHTS® software and the CFDEM®coupling software
   See http://www.cfdem.com/terms-trademark-policy for details.

   LIGGGHTS® is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
Contributing authors:
Andreas Aigner (JKU Linz)
Christoph Kloss (JKU Linz and DCS Computing Gmbh, Linz)
------------------------------------------------------------------------- */

#ifdef ATOM_CLASS

AtomStyle(sph/var,AtomVecSPH2)

#else

#ifndef LMP_ATOM_VEC_SPH2_H
#define LMP_ATOM_VEC_SPH2_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecSPH2 : public AtomVec {
 public:
  AtomVecSPH2(class LAMMPS *);
  ~AtomVecSPH2() {}
  void init();
  void grow(int);
  void grow_reset();
  void copy(int, int, int);

  int pack_comm(int, int *, double *, int, int *);
  int pack_comm_vel(int, int *, double *, int, int *);
  int pack_comm_hybrid(int, int *, double *);
  void unpack_comm(int, int, double *);
  void unpack_comm_vel(int, int, double *);
  int unpack_comm_hybrid(int, int, double *);

  int pack_reverse(int, int, double *);
  void unpack_reverse(int, int *, double *);

  int pack_border(int, int *, double *, int, int *);
  int pack_border_vel(int, int *, double *, int, int *);
  int pack_border_hybrid(int, int *, double *);
  void unpack_border(int, int, double *);
  void unpack_border_vel(int, int, double *);
  int unpack_border_hybrid(int, int, double *);

  int pack_exchange(int, double *);
  int unpack_exchange(double *);

  int size_restart();
  int pack_restart(int, double *);
  int unpack_restart(double *);

  void create_atom(int, double *);
  void data_atom(double *, tagint, char **);
  int data_atom_hybrid(int, char **);
  void data_vel(int, char **);
  int data_vel_hybrid(int, char **);
  void pack_data(double **);
  int pack_data_hybrid(int, double *);
  void write_data(FILE *, int, double **);
  int write_data_hybrid(FILE *, double *);
  void pack_vel(double **);
  int pack_vel_hybrid(int, double *);
  void write_vel(FILE *, int, double **);
  int write_vel_hybrid(FILE *, double *);
  bigint memory_usage();

 private:
  int *tag,*type,*mask,*image;
  double **x,**v,**f;
  double *p,*rho,*drho,*e,*de;
  double *radius,*rmass;
  int radvary;
};

}

#endif
#endif
