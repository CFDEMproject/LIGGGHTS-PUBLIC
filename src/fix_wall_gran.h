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

    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Richard Berger (JKU Linz)
    Arno Mayrhofer (CFDEMresearch GmbH, Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
    Copyright 2016-     CFDEMresearch GmbH, Linz
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

 friend class LIGGGHTS::Walls::IGranularWall;

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

  virtual void createMulticontactData();

  void setDnum(int newDnum)
  { dnum_ = newDnum; }

  inline int store_force() const
  { return store_force_; }

  inline FixPropertyAtom* fix_wallforce() const
  { return fix_wallforce_; }

  inline int iarg() const
  { return iarg_; }

  int add_history_value(std::string name, std::string newtonflag)
  {  return dnum_++; }

  int get_history_offset(const std::string hname)
  {  return impl->get_history_offset(hname);}

  bool contact_match(const std::string mtype, const std::string model)
  {  return impl->contact_match(mtype, model); }

  inline int dnum() const
  { return dnum_; }

  inline int n_meshes() const
  { return n_FixMesh_; }

  inline class FixMeshSurface ** mesh_list() const
  { return FixMesh_list_; }

  inline int atom_type_wall() const
  { return atom_type_wall_; }

  inline bool is_mesh_wall() const
  { return 1 == meshwall_; }

  inline int computeflag() const
  { return computeflag_; }

  inline int shearupdate() const
  { return shearupdate_; }

  inline int heattransfer_flag() const
  { return heattransfer_flag_; }

  inline int stress_flag() const
  { return stress_flag_; }

  inline bool store_force_contact() const
  { return store_force_contact_; }

  inline int store_force_contact_every() const
  { return store_force_contact_every_; }

  inline bool store_force_contact_stress() const
  { return store_force_contact_stress_; }

  inline ComputePairGranLocal * compute_wall_gran_local() const
  { return cwl_; }

  inline int addflag() const
  { return addflag_; }

  inline int body(int i) const
  { return body_[i]; }

  inline double masstotal(int i) const
  { return masstotal_[i]; }

  inline class FixRigid *fix_rigid() const
  { return fix_rigid_; }

  inline void add_contactforce_wall(int ip, const LCM::ForceData & i_forces,int idTri)
  {
    // add to fix wallforce contact
    // adds 0 as ID for primitive wall
    double forces_torques_i[6];

    if(fix_wallforce_contact_->has_partner(ip,idTri) == -1)
    {
      vectorCopy3D(i_forces.delta_F,&(forces_torques_i[0]));
      vectorCopy3D(i_forces.delta_torque,&(forces_torques_i[3]));
      fix_wallforce_contact_->add_partner(ip,idTri,forces_torques_i);
      
    }
  }

  inline void add_contactforce_stress_wall(int ip, const LCM::ForceData & i_forces, const double *const delta, const double *const vwall, int idTri)
  {
    // add to fix wallforce contact
    // adds 0 as ID for primitive wall
    double forces_delta_i[9];

    if(fix_wallforce_contact_stress_->has_partner(ip,idTri) == -1)
    {
      vectorCopy3D(i_forces.delta_F,&(forces_delta_i[0]));
      vectorCopy3D(delta,&(forces_delta_i[3]));
      vectorCopy3D(vwall,&(forces_delta_i[6]));
      fix_wallforce_contact_stress_->add_partner(ip,idTri,forces_delta_i);
    }
  }

  bool store_sum_normal_force() const
  { return fix_sum_normal_force_ != NULL; }

  double * get_sum_normal_force_ptr(const int i)
  { return &(fix_sum_normal_force_->vector_atom[i]); }

  class PrimitiveWall* primitiveWall();

  int n_contacts_all(int &nIntersect);
  int n_contacts_all(int contact_groupbit,int &nIntersect);
  int n_contacts_local(int &nIntersect);
  int n_contacts_local(int contact_groupbit,int &nIntersect);
  int is_moving();

  void register_compute_wall_local(ComputePairGranLocal *,int&);
  void unregister_compute_wall_local(ComputePairGranLocal *ptr);

  void wall_temperature_unique(bool &has_temp,bool &temp_unique, double &temperature_unique);
  void addHeatFlux(class TriMesh *mesh,int i,const double ri,double rsq,double area_ratio);

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

#ifdef SUPERQUADRIC_ACTIVE_FLAG
  double **quat_, **shape_, **blockiness_;
#endif

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
  class FixPropertyAtom *fppa_htcw; 

  double Temp_wall;
  double fixed_contact_area_;
  double Q,Q_add;

  const double *th_cond;
  double const* const* deltan_ratio;

  LIGGGHTS::Walls::IGranularWall * impl;

  // per-contact force storage
  bool store_force_contact_;
  int store_force_contact_every_;
  class FixContactPropertyAtomWall *fix_wallforce_contact_;

  // for stress computation
  bool store_force_contact_stress_;
  class FixContactPropertyAtomWall *fix_wallforce_contact_stress_;

  // storage for per contact data (for multicontact models)
  class FixContactPropertyAtomWall *fix_store_multicontact_data_;

  int nlevels_respa_;

  int shear_, shearDim_, shearAxis_;
  double vshear_;
  double shearAxisVec_[3];

  // distance in order to calculate interaction with
  // rough wall
  //double skinDistance_;

  // number of values for contact history
  int dnum_;

  // flag if mesh wall
  int meshwall_;

  // flag to actiate potential energy calculation
  bool track_energy_;

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

  // storage for simplistic pressure computation via normal forces
  class FixPropertyAtom *fix_sum_normal_force_;

  // max neigh cutoff - as in Neighbor
  double cutneighmax_;

  virtual void post_force_wall(int vflag);
};

}

#endif
#endif
