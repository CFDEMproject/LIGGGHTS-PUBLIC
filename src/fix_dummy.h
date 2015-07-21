#ifndef LMP_FIX_DUMMY_H
#define LMP_FIX_DUMMY_H

#include "fix_rigid.h"
#include "custom_value_tracker.h"

namespace LAMMPS_NS {

class MultisphereParallel {

 public:
  int n_body() {return 0;}
  int map(int i) {return 0;}
  int tag(int i) {return 0;}
  int tag_max_body() {return 0;}

  void id_extend() {}
  void generate_map() {}

  inline void x_bound(double *x_bnd,int i){  }

  inline double r_bound(int i)
  {
      return 0.;
  }

  inline class CustomValueTracker& prop()
  {
    LAMMPS *dptr = 0;
    return *(new CustomValueTracker(dptr));
  }

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
};

class FixMultisphere : public Fix {

 public:

  FixMultisphere(class LAMMPS * lmp, int narg, char ** arg) : Fix(lmp, narg, arg) {}

  void set_v_integrate(double *v) {}
  int belongs_to(int i) {return -1;}

  void* extract(const char*& a, int& b, int& c) {return NULL;}

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

  class MultisphereParallel *multisphere_;

};
}

#endif
