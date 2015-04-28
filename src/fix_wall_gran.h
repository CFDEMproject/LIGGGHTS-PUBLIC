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

    Christoph Kloss (DCS Computing GmbH, Linz, JKU Linz)
    Philippe Richard Berger (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
FixStyle(wall/gran,FixWallGran)
#else

#ifndef LMP_FIX_WALL_GRAN_H
#define LMP_FIX_WALL_GRAN_H

#include "fix_mesh_surface.h"
#include "contact_interface.h"
#include "granular_wall.h"
#include <string>
#include <vector>
#include "fix_contact_property_atom_wall.h"
#include "compute_pair_gran_local.h"

namespace LCM = LIGGGHTS::ContactModels;

namespace LAMMPS_NS {

class FixWallGran : public Fix, public LIGGGHTS::IContactHistorySetup {
 public:
  FixWallGran(class LAMMPS *, int, char **);
  ~FixWallGran();

  /* INHERITED FROM Fix */

  virtual int setmask();
  virtual void post_create();
  virtual void pre_delete(bool unfixflag);
  virtual void init();
  virtual void setup(int vflag);
  virtual void post_force(int vflag);
  virtual void post_force_pgl();
  virtual void post_force_respa(int, int, int);

  virtual int min_type();
  virtual int max_type();

  /* PUBLIC ACCESS FUNCTIONS */

  void setSkinDistance(double newSkinDistance)
  { skinDistance_ = newSkinDistance; }

  void setDnum(int newDnum)
  { dnum_ = newDnum; }

  inline int store_force()
  { return store_force_; }

  inline int iarg()
  { return iarg_; }

  int add_history_value(std::string name, std::string newtonflag) {
      return dnum_++;
  }

  inline int dnum()
  { return dnum_; }

  inline int n_meshes()
  { return n_FixMesh_; }

  inline class FixMeshSurface ** mesh_list()
  { return FixMesh_list_; }

  inline int atom_type_wall()
  { return atom_type_wall_; }

  inline bool is_mesh_wall()
  { return 1 == meshwall_; }

  inline bool store_force_contact()
  { return store_force_contact_; }

  class PrimitiveWall* primitiveWall();

  int n_contacts_all();
  int n_contacts_all(int);
  int n_contacts_local();
  int n_contacts_local(int);
  int is_moving();

  void register_compute_wall_local(ComputePairGranLocal *,int&);
  void unregister_compute_wall_local(ComputePairGranLocal *ptr);

  ComputePairGranLocal * compute_pair_gran_local() {
    return cwl_;
  }

  int addflag() const {
    return addflag_;
  }

  int body(int i) {
    return body_[i];
  }

  double masstotal(int i) {
    return masstotal_[i];
  }

  class FixRigid *fix_rigid() {
    return fix_rigid_;
  }

  void add_contactforce_wall(int ip, const LCM::ForceData & i_forces,int idTri)
  {
    // add to fix wallforce contact
    // adds 0 as ID for primitive wall
    double forces_torques_i[6];

    if(!fix_wallforce_contact_->has_partner(ip,idTri))
    {
      vectorCopy3D(i_forces.delta_F,&(forces_torques_i[0]));
      vectorCopy3D(i_forces.delta_torque,&(forces_torques_i[3]));
      fix_wallforce_contact_->add_partner(ip,idTri,forces_torques_i);
    }
  }

  void cwl_add_wall_2(const LCM::SurfacesIntersectData & sidata, const LCM::ForceData & i_forces)
  {
    const double fx = i_forces.delta_F[0];
    const double fy = i_forces.delta_F[1];
    const double fz = i_forces.delta_F[2];
    const double tor1 = i_forces.delta_torque[0]*sidata.area_ratio;
    const double tor2 = i_forces.delta_torque[1]*sidata.area_ratio;
    const double tor3 = i_forces.delta_torque[2]*sidata.area_ratio;
    cwl_->add_wall_2(sidata.i,fx,fy,fz,tor1,tor2,tor3,sidata.contact_history,sidata.rsq);
  }

 protected:

  int iarg_, narg_;
  int atom_type_wall_;

  int computeflag_;

  int addflag_;
  ComputePairGranLocal *cwl_;

  double dt_;
  int shearupdate_;

  // class variables for atom properties
  int nlocal_;
  double **x_, **f_, *radius_, *rmass_, **wallforce_, r0_;

  void set_r0(double _r0)
  { r0_ = _r0; }

  virtual void init_granular() {}

  // heat transfer
  void init_heattransfer();
  bool heattransfer_flag_;
  // model for contact area calculation
  int area_calculation_mode_;

  // mesh and primitive force implementations
  virtual void post_force_mesh(int);
  virtual void post_force_primitive(int);

  // virtual functions that allow implementation of the
  // actual physics in the derived classes
  virtual void compute_force(LCM::SurfacesIntersectData & sidata, double *vwall);
  void addHeatFlux(class TriMesh *mesh,int i,double rsq,double area_ratio);

  // sets flag that neigh list shall be built
  virtual void pre_neighbor();

  // builds neigh list if necessary
  virtual void pre_force(int vflag);

  // references to mesh walls
  int n_FixMesh_;
  class FixMeshSurface **FixMesh_list_;
  class FixRigid *fix_rigid_;
  int *body_;
  double *masstotal_;

  // heat transfer
  class FixPropertyAtom *fppa_T;
  class FixPropertyAtom *fppa_hf;

  double Temp_wall;
  double fixed_contact_area_;
  double Q,Q_add;

  const double *th_cond;
  double const* const* deltan_ratio;

  LIGGGHTS::Walls::IGranularWall * impl;

  // per-contact force storage
  bool store_force_contact_;
  class FixContactPropertyAtomWall *fix_wallforce_contact_;

  int nlevels_respa_;

  int shear_, shearDim_, shearAxis_;
  double vshear_;
  double shearAxisVec_[3];

  // distance in order to calculate interaction with
  // rough wall
  double skinDistance_;

  // number of values for contact history
  int dnum_;

  // flag if mesh wall
  int meshwall_;

  // flag for stressanalysis
  // true if any of the meshes tracks stresses
  bool stress_flag_;

  class PrimitiveWall *primitiveWall_;
  class FixPropertyAtom *fix_history_primitive_;

  // class to keep track of wall contacts
  bool rebuildPrimitiveNeighlist_;

  // force storage
  bool store_force_;
  class FixPropertyAtom *fix_wallforce_;

  // max neigh cutoff - as in Neighbor
  double cutneighmax_;

  virtual void post_force_wall(int vflag);

  inline void post_force_eval_contact(LCM::SurfacesIntersectData & sidata, double * v_wall, int iMesh = -1, FixMeshSurface *fix_mesh = 0, TriMesh *mesh = 0, int iTri = 0);
};

}

#endif
#endif
