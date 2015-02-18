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

FixStyle(massflow/mesh,FixMassflowMesh)

#else

#ifndef LMP_FIX_MASSFLOW_MESH_H
#define LMP_FIX_MASSFLOW_MESH_H

#include "fix.h"
#include <vector>

using namespace std;

namespace LAMMPS_NS {

class FixMassflowMesh : public Fix {

 public:

  FixMassflowMesh(class LAMMPS *lmp, int narg, char ** arg);
  ~FixMassflowMesh();

  void post_create();
  void pre_delete(bool unfixflag);

  void init();
  void setup(int vflag);
  int setmask();

  void post_integrate();
  void pre_exchange();

  void write_restart(FILE *fp);
  void restart(char *buf);

  double compute_vector(int index);

 protected:

  // in case particles counted should be deleted or transferred
  bool delete_atoms_;
  vector<int> atom_tags_delete_;
  double mass_deleted_;
  double nparticles_deleted_;

  // true if any given particle is
  // counted only once
  bool once_;

  int  iarg_;
  class FixPropertyAtom* fix_orientation_;

 private:

  class FixMeshSurface *fix_mesh_;
  class FixPropertyAtom *fix_counter_;
  char fixid_[200];
  class FixNeighlistMesh *fix_neighlist_;
  double nvec_[3];
  double pref_[3];
  double sidevec_[3];

  bool   havePointAtOutlet_;
  bool   insideOut_;
  double pointAtOutlet_[3];

  // mass and particles which was counted
  double mass_;
  int nparticles_;

  // additional property to sum
  class FixPropertyAtom *fix_property_;
  double property_sum_;

  // data write
  bool screenflag_;
  FILE *fp_;

  // data for particle and mass flow calculation
  double mass_last_;
  int nparticles_last_;
  double t_count_, delta_t_;
  bool reset_t_count_;

  class FixMultisphere* fix_ms_;
  class MultisphereParallel *ms_;
  class ScalarContainer<int> *ms_counter_;

}; //end class

}
#endif
#endif
