#ifndef LMP_FIX_DUMMY_H
#define LMP_FIX_DUMMY_H

#include "fix_rigid.h"

namespace LAMMPS_NS {

class FixRigidMultisphere : public FixRigid {

 public:

  int map(int i) {return 0;}
  int tag(int i) {return 0;}
  int tag_max_body() {return 0;}
  int n_bodies() {return 0;}
  void* extract(char*& a, int& b, int& c) {return NULL;}
  void set_v_integrate(double *v) {}
  int belongs_to(int i) {return 0;}
  int calc_n_steps(int iatom,double *p_ref,double *normalvec,double *v_normal) {return 0;}

  inline void x_bound_body(double *x_bnd,int i){  }

  inline double r_bound_body(int i)
  {
      return 0.;
  }

  inline double mass_body(int i)
  {
      return 0.;
  }

  inline double density_body(int i)
  {
      return 0.;
  }

};
}

#endif
