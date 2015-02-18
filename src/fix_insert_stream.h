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

#ifdef FIX_CLASS

FixStyle(insert/stream,FixInsertStream)

#else

#ifndef LMP_FIX_INSERT_STREAM_H
#define LMP_FIX_INSERT_STREAM_H

#include "fix_insert.h"
#include "vector_liggghts.h"

namespace LAMMPS_NS {

class FixInsertStream : public FixInsert {

 public:

  FixInsertStream(class LAMMPS *, int, char **);
  ~FixInsertStream();
  virtual void post_create();
  void pre_delete(bool unfixflag);

  virtual int setmask();
  virtual void init();
  virtual void setup_pre_exchange();
  virtual void end_of_step();

  void init_defaults();

  virtual void reset_timestep(bigint newstep,bigint oldstep);

  void register_tracer_callback(class FixPropertyAtomTracerStream* tr);

  class FixPropertyAtom * fix_prop_release()
  { return fix_release; }

  class TriMesh* face()
  { return ins_face; }

  virtual int release_step_index()
  { return 4; }

  bool vel_is_normal_to_face()
  { return vel_normal_to_face; }

  void vel_normal(double *vn)
  { return vectorCopy3D(v_normal,vn); }

 protected:

  virtual void calc_insertion_properties();

  void pre_insert();

  int is_nearby(int);
  inline void generate_random(double *pos, double rad);
  inline void generate_random_global(double *pos);

  void x_v_omega(int ninsert_this,int &ninserted_this,
                int &ninserted_spheres_this,
                double &mass_inserted_this);

  double insertion_fraction();
  void calc_ins_fraction();
  virtual void finalize_insertion(int);

  virtual void reset_releasedata(bigint newstep,bigint oldstep);

  // additional insertion settings
  int duration;            //duration for insertion in time-steps
  bool parallel;

  // stuff for insertion region
  double normalvec[3];     // points out of extruded volume
  bool vel_normal_to_face;       // points out of extruded volume
  double extrude_length;
  double extrude_length_min, extrude_length_max;
  double p_ref[3];         // reference point on face
  int face_style;
  double v_normal[3];      // insertion velocity projected on face
  double ins_fraction;
  bool do_ins_fraction_calc;

  // mesh face and bounding box of extruded face
  class TriMesh *ins_face;
  double ins_vol_xmin[3];
  double ins_vol_xmax[3];
  int ntry_mc;

  // non-mesh face
  double c_center[3], c_r; // non mesh face - currently only circle

  // per particle release step
  class FixPropertyAtom *fix_release;
  bool i_am_integrator;

 private:

  class FixPropertyAtomTracerStream **tracer;
  int ntracer;
};

}

#endif
#endif
