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

#ifdef COMPUTE_CLASS

ComputeStyle(pair/gran/local,ComputePairGranLocal)
ComputeStyle(wall/gran/local,ComputePairGranLocal)

#else

#ifndef LMP_COMPUTE_PAIR_GRAN_LOCAL_H
#define LMP_COMPUTE_PAIR_GRAN_LOCAL_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePairGranLocal : public Compute {

 public:
  ComputePairGranLocal(class LAMMPS *, int, char **);
  ~ComputePairGranLocal();
  void post_create();
  void init();
  void init_cpgl(bool requestflag);
  void init_list(int, class NeighList *);
  void compute_local();
  double memory_usage();
  void reference_deleted();
  void add_pair(int i,int j,double fx,double fy,double fz,double tor1,double tor2,double tor3,double *hist);
  void add_heat(int i,int j,double hf);
  void add_wall_1(int iFMG,int iTri,int iP,double *contact_point,double *v_wall);
  void add_wall_2(int i,double fx,double fy,double fz,double tor1,double tor2,double tor3,double *hist,double rsq);
  void add_heat_wall(int i,double hf);

 private:
  int nvalues;
  int ncount;
  int newton_pair;

  // if 0, pair data is extracted
  // if 1, wall data is extracted
  int wall;

  // flag if compute is working
  // can be set via pair, fix
  int reference_exists;

  // pointers to classes holding the data
  class PairGran *pairgran;
  class FixHeatGranCond *fixheat;
  class FixWallGran *fixwall;

  int ipair;

  int posflag,velflag,idflag,fflag,tflag,hflag,aflag,hfflag;

  int dnum;

  int nmax;
  double *vector;
  double **array;

  class NeighList *list;

  int count_pairs();
  int count_wallcontacts();
  void reallocate(int);
};

}

#endif
#endif
