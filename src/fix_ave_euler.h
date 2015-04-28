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
    (if not contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
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

  int ncells_pack();

  // inline access functions for cell based values

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

  inline int ntry_per_cell()
  { return 50; }

  void setup_bins();
  void bin_atoms();
  void calculate_eu();
  void allreduce();
  inline int coord2bin(double *x); 

  bool parallel_;

  int exec_every_;
  bool box_change_size_, box_change_domain_;
  int triclinic_; 

  // desired cell size over max particle diameter
  double cell_size_ideal_rel_;

  // desired cell size
  double cell_size_ideal_;
  double cell_size_ideal_lamda_[3];

  // number of cells, either globally or locally on each proc
  int ncells_, ncells_dim_[3];

  // extent of grid in xyz, either globally or locally on each proc
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

  // region
  char *idregion_;
  class Region *region_;

  /* ---------  DATA  --------- */

  // cell center
  double **center_;

  // cell-based averaged velocity
  double **v_av_;

  // cell-based volume fraction
  double *vol_fr_;

  // cell-based weight for each cell
  
  double *weight_;

  // cell-based average radius
  double *radius_;

  // cell-based number of particles
  int *ncount_;

  // cell-based mass
  double *mass_;

  // cell-based stress
  // [0]: pressure
  // [1-3]: 00-11-22
  // [4-6]: 01-02-12
  double **stress_;

  // stress computation
  class ComputeStressAtom *compute_stress_;

  class RanPark *random_;
};

}

#endif
#endif
