#include "atom_vec.h"
      namespace LAMMPS_NS { class AtomVecConvexHull : public AtomVec { public: int get_ntri_max(){return 0;} int get_ntri(int){return 0;} void get_tri_node(int,int,int, double*) {}  }; } 
