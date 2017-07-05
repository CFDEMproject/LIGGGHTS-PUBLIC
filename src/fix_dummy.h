#ifndef LMP_FIX_DUMMY_H
#define LMP_FIX_DUMMY_H

#include "fix_rigid.h"
#include "fix_property_atom.h"
#include "custom_value_tracker.h"

namespace LAMMPS_NS {

class Multisphere {

 public:
  int n_body() {return 0;}
  int map(int i) {return 0;}
  int tag(int i) {return 0;}
  int tag_max_body() {return 0;}

  void id_extend() {}
  void generate_map() {}

  inline void x_bound(double *x_bnd,int i){  }
  inline void xcm(double *x_cm,int ibody_local) { }
  inline double volume(int ibody_local) { return 0.;}

  inline double r_bound(int i)
  {
      return 0.;
  }

  inline class CustomValueTracker& prop()
  {
    LAMMPS *dptr = 0;
    return *(new CustomValueTracker(dptr));
  }
  void recalc_n_steps(double dt_ratio) {}
  inline double mass(int i)
  {
      return 0.;
  }

  inline double density(int i)
  {
      return 0.;
  }
  int calc_n_steps(int iatom,double *p_ref,double *normalvec,double *v_normal) {return 0;}
  void* extract(const char*& a, int& b, int& c) {return NULL;}

      void add_body(int nspheres, double *xcm_ins, double *xcm_to_xbound_ins,
                    double r_bound_ins, double *v_ins, double *omega_ins,
                    double mass_ins, double dens_ins, int atomtype_ins, int type_ins,
                    double *inertia_ins, double *ex_space_ins, double *ey_space_ins, double *ez_space_ins,
                    double **displace_ins,bool *fflag, bool *tflag, int start_step_ins = -1, double *v_integrate_ins = NULL) {}

};

class MultisphereParallel : public Multisphere {};

class FixMultisphere : public Fix {

 public:

  FixMultisphere(class LAMMPS * lmp, int narg, char ** arg) : Fix(lmp, narg, arg) {}

  void set_v_integrate(double *v) {}
  int belongs_to(int i) {return -1;}

  void* extract(const char*& a, int& b, int& c) {return NULL;}

  inline double extract_ke()
  { return 0.; }

  inline double extract_rke()
  { return 0.; }
  void release(int,double*,double*) {}

  int calc_n_steps(int iatom,double *p_ref,double *normalvec,double *v_normal)
  { return 0; }

  inline class MultisphereParallel& data()
  { return (*multisphere_);}

  inline int n_body()
  { return data().n_body(); }

  inline int tag_max_body()
  { return data().tag_max_body(); }

  int ntypes() {return 0;}
  double * vclump() {return 0;}

  bool allow_group_and_set()
  { return false; }

  inline void set_v_body_from_atom_index(int iatom,double *vel)
  {}

  inline void set_omega_body_from_atom_index(int iatom,double *omega)
  {}
  class FixPropertyAtom *get_volumeweight() 
      { return 0; }
  class MultisphereParallel *multisphere_;

  void add_body_finalize() {}

  void set_body_displace(int i,double *_displace,int body_id,double volume_weight) {}

};
}

#endif
