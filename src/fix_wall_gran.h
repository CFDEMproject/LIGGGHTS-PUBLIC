/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   This file was modified with respect to the release in LAMMPS
   Modifications are Copyright 2009-2012 JKU Linz
                     Copyright 2012-     DCS Computing GmbH, Linz

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

// this class cannot be instantiated

#else

#ifndef LMP_FIX_WALL_GRAN_H
#define LMP_FIX_WALL_GRAN_H

#include "fix_mesh_surface.h"

namespace LAMMPS_NS {

class FixWallGran : public Fix {

 public:
  FixWallGran(class LAMMPS *, int, char **);
  ~FixWallGran();

  /* INHERITED FROM Fix */

  virtual int setmask();
  void post_create();
  void pre_delete(bool unfixflag);
  void init();
  void setup(int);
  void post_force(int);
  void post_force(int,int);
  void post_force_respa(int, int, int);

  int min_type();
  int max_type();

  /* PUBLIC ACCESS FUNCTIONS */

  inline int dnum()
  { return dnum_; }

  inline int n_meshes()
  { return n_FixMesh_; }

  inline class FixMeshSurface ** mesh_list()
  { return FixMesh_list_; }

  inline int atom_type_wall()
  { return atom_type_wall_; }

  inline bool is_mesh_wall()
  { return meshwall_ == 1; }

  int n_contacts_all();
  int n_contacts_all(int);
  int n_contacts_local();
  int n_contacts_local(int);
  int is_moving();

  void register_compute_wall_local(class ComputePairGranLocal *,int&);
  void unregister_compute_wall_local(class ComputePairGranLocal *ptr);

 protected:

  int iarg_, narg_;
  int atom_type_wall_;

  int addflag_;
  class ComputePairGranLocal *cwl_;

  double dt_;
  int shearupdate_;
  int laststep_;

  void set_r0(double _r0)
  { r0_ = _r0; }

  virtual void init_granular() {}
  virtual void init_heattransfer() {}
  bool heattransfer_flag_;

  // virtual functions that allow implementation of the
  // actual physics in the derived classes
  virtual void compute_force(int i,double deltan,double rsq,double meff_wall,
                            double dx,double dy,double dz,double *vwall,
                            double *c_history,double area_ratio) = 0;
  virtual void addHeatFlux(class TriMesh *mesh,int i,double rsq,double area_ratio) {};

  // sets flag that neigh list shall be built
  virtual void pre_neighbor();
  // builds neigh list if necessary
  virtual void pre_force(int vflag);

  // pair style, fix rigid for correct damping
  char *pairstyle_;
  class PairGran *pairgran_;
  class FixRigid *fix_rigid_;
  int *body_;
  double *masstotal_;

 private:

  int nlevels_respa_;

  int shear_, shearDim_, shearAxis_;
  double vshear_;
  double shearAxisVec_[3];

  // number of values for contact history
  int dnum_;

  // flag if mesh wall
  int meshwall_;

  // references to mesh walls
  int n_FixMesh_;
  class FixMeshSurface **FixMesh_list_;

  // flag for stressanalysis
  // true if any of the meshes tracks stresses
  bool stress_flag_;

  class PrimitiveWall *primitiveWall_;

  // class to keep track of wall contacts
  bool rebuildPrimitiveNeighlist_;

  // force storage
  bool store_force_;
  class FixPropertyAtom *fix_wallforce_;

  // class variables for atom properties
  int nlocal_;
  double **x_, **f_, *radius_, *rmass_, **wallforce_, r0_;

  // max neigh cutoff - as in Neighbor
  double cutneighmax_;

  void post_force_mesh(int);
  void post_force_primitive(int);

  inline void post_force_eval_contact(int iPart, double deltan, double *delta, double *v_wall,
                    double *c_history, int iMesh = -1, FixMeshSurface *fix_mesh = 0, TriMesh *mesh = 0, int iTri = 0);
};

}

#endif
#endif
