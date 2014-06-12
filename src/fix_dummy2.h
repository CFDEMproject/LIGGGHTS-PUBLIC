#ifndef LMP_FIX_DUMMY2_H
#define LMP_FIX_DUMMY2_H

namespace LAMMPS_NS {

class FixContactPropertyAtom : public Fix {

 public:

  FixContactPropertyAtom(class LAMMPS * lmp, int narg, char ** arg) : Fix(lmp, narg, arg) {}

  inline bool has_partner(int,int)
  {
    return false;
  }

  inline void add_partner(int,int,double*)
  {
  }

  void do_forward_comm() {}
};
}

#endif
