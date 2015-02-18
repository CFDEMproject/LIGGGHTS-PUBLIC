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

/* ----------------------------------------------------------------------
Contributing author for SPH:
Andreas Aigner (CD Lab Particulate Flow Modelling, JKU)
andreas.aigner@jku.at
------------------------------------------------------------------------- */

#ifndef LMP_FIX_SPH
#define LMP_FIX_SPH

#include "fix.h"

namespace LAMMPS_NS {

class FixSph : public Fix {
 public:
  FixSph(class LAMMPS *, int, char **);
  ~FixSph();
  int setmask();
  virtual void updatePtrs();
  void init();
  void init_list(int, class NeighList *);
  virtual void post_integrate() {};
  virtual void post_integrate_respa(int, int);

  int get_kernel_id(){return kernel_id;};
  inline void set_kernel_id(int newid){kernel_id = newid;};

  int kernel_flag;        // 1 if Fix uses sph kernel, 0 if not

 protected:
  inline double interpDist(double disti, double distj) {return 0.5*(disti+distj);};

  class FixPropertyAtom* fppaSl; //smoothing length
  class FixPropertyGlobal* fppaSlType; //per type smoothing length
  double *sl;         // per atom smoothing length
  double **slComType; // common smoothing length in case of mass_type=1

  int kernel_id;
  double kernel_cut;
  char *kernel_style;

  class NeighList *list;
  int nlevels_respa;

  int mass_type; // flag defined in atom_vec*

};

}

#endif
