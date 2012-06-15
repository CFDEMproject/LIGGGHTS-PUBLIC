/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(insert/stream,FixInsertStream)

#else

#ifndef LMP_FIX_INSERT_STREAM_H
#define LMP_FIX_INSERT_STREAM_H

#include "fix_insert.h"

namespace LAMMPS_NS {

class FixInsertStream : public FixInsert {
 public:

  FixInsertStream(class LAMMPS *, int, char **);
  ~FixInsertStream();
  void post_create();
  void pre_delete(bool unfixflag);

  virtual int setmask();
  virtual void init();
  virtual void end_of_step();

  void init_defaults();

 private:

  virtual void calc_insertion_properties();

  void pre_insert();

  int is_nearby(int);
  inline void generate_random(double *pos, double rad);

  void x_v_omega(int ninsert_this,int &ninserted_this, int &ninserted_spheres_this, double &mass_inserted_this);
  double insertion_fraction();
  virtual void finalize_insertion(int);

  // additional insertion settings
  int duration;            //duration for insertion in time-steps

  // stuff for insertion region
  double normalvec[3];
  double extrude_length;
  double p_ref[3];         //reference point on face
  int face_style;
  double v_normal[3];      // insertion velocity projected on face

  // mesh face and bounding box of extruded face
  class TriMesh *ins_face;
  double ins_vol_xmin[3];
  double ins_vol_xmax[3];

  // non-mesh face
  double c_center[3], c_r; // non mesh face - currently only circle

  // per particle release step
  class FixPropertyAtom *fix_release;
  bool i_am_integrator;

};

}

#endif
#endif
