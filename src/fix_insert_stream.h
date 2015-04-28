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
