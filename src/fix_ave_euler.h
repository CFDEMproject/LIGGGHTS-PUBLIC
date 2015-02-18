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
   Contributing author for triclinic: Andreas Aigner (JKU)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(ave/euler,FixAveEuler)
FixStyle(ave/euler/stress,FixAveEuler)

#else

#ifndef LMP_FIX_AVE_EULER_H
#define LMP_FIX_AVE_EULER_H

#include "stdio.h"
#include "fix.h"
#include "vector_liggghts.h"

namespace LAMMPS_NS {

class FixAveEuler : public Fix {

 public:

  FixAveEuler(class LAMMPS *, int, char **);
  ~FixAveEuler();

  void post_create();
  int setmask();
  void init();
  void setup(int vflag);

  void end_of_step();

  double compute_array(int i, int j);

  // inline access functions for cell based values

  inline int ncells()
  { return ncells_; }

  inline double cell_center(int i, int j)
  { return center_[i][j]; }

  inline double cell_v_av(int i, int j)
  { return v_av_[i][j]; }

  inline double cell_vol_fr(int i)
  { return vol_fr_[i]; }

  inline double cell_radius(int i)
  { return radius_[i]; }

  inline double cell_pressure(int i)
  { return stress_[i][0]; }

  inline double cell_stress(int i,int j)
  { return stress_[i][j+1]; }

 private:

  void setup_bins();
  void bin_atoms();
  void calculate_eu();
  inline int coord2bin(double *x); 

  int exec_every_;
  bool box_change_;
  int triclinic_; 

  // desired cell size over max particle diameter
  double cell_size_ideal_rel_;

  // desired cell size
  double cell_size_ideal_;
  double cell_size_ideal_lamda_[3];

  // number of cells
  int ncells_, ncells_dim_[3];

  // extent of grid in xyz
  double lo_[3],hi_[3];
  double lo_lamda_[3],hi_lamda_[3]; 

  // cell size and inverse size in xyz, cell and volume
  double cell_size_[3];
  double cell_size_inv_[3];
  double cell_volume_;
  double cell_size_lamda_[3]; 
  double cell_size_lamda_inv_[3]; 

  // length of cellhead_, center_, v_av_, vol_fr_ arrays
  int ncells_max_;

  // length of cellptr_ array
  int ncellptr_max_;

  // atom - cell mapping
  int *cellhead_;    // ptr to 1st atom in each cell
  int *cellptr_;       // ptr to next atom in each bin

  /* ---------  DATA  --------- */

  // cell center
  double **center_;

  // cell-based averaged velocity
  double **v_av_;

  // cell-based volume fraction
  double *vol_fr_;

  // cell-based average radius
  double *radius_;

  // cell-based mass
  double *mass_;

  // cell-based stress
  // [0]: pressure
  // [1-3]: 00-11-22
  // [4-6]: 01-02-12
  double **stress_;

  // stress computation
  class ComputeStressAtom *compute_stress_;
};

}

#endif
#endif
