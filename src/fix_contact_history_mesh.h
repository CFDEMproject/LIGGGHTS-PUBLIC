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

FixStyle(contacthistory/mesh,FixContactHistoryMesh) 

#else

#ifndef LMP_FIX_CONTACT_HISTORY_MESH_H
#define LMP_FIX_CONTACT_HISTORY_MESH_H

#include "fix_contact_history.h"
#include "fix_property_atom.h"
#include "my_page.h"
#include "math.h"
#include "vector_liggghts.h"
#include "atom.h"
#include "update.h"
#include "error.h"
#include "tri_mesh.h"

namespace LAMMPS_NS {

class FixContactHistoryMesh : public FixContactHistory {
  friend class Neighbor;
  friend class PairGran;

 public:
  FixContactHistoryMesh(class LAMMPS *, int, char **);
  ~FixContactHistoryMesh();
  virtual int setmask();
  void init();

  void setup_pre_exchange();
  void min_setup_pre_exchange();
  void pre_exchange();

  void setup_pre_force(int dummy);
  void setup_pre_neighbor();
  void min_setup_pre_force(int dummy);
  void pre_force(int dummy);
  void min_pre_force(int dummy);

  void pre_neighbor();

  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int unpack_exchange(int, double *);
  void unpack_restart(int, int);
  void write_restart(FILE *fp);
  double memory_usage();

  // spefific interface for mesh

  bool handleContact(int iPart, int idTri, double *&history);
  void markAllContacts();
  void cleanUpContacts();
  void cleanUpContactJumps();
  
  // OMP interface
  void resetDeletionPage(int tid);
  void markForDeletion(int tid, int ifrom, int ito);
  void cleanUpContacts(int ifrom, int ito);

  void reset_history();

  // return # of contacts
  int n_contacts();
  int n_contacts(int contact_groupbit);

 protected:

  MyPage<int> *ipage1_;        // pages of neighbor tri IDs
  MyPage<double> *dpage1_;     // pages of contact history with neighbors
  MyPage<int> *ipage2_;        // pages of neighbor tri IDs
  MyPage<double> *dpage2_;     // pages of contact history with neighbors
  MyPage<bool> ** keeppage_;   // pages of deletion flags with neighbors

  bool **keepflag_;

  void allocate_pages();

 private:

  // functions specific for mesh - contact management
  bool haveContact(int indexPart, int idTri, double *&history);
  bool coplanarContactAlready(int indexPart, int idTri);
  void checkCoplanarContactHistory(int indexPart, int idTri, double *&history);
  void addNewTriContactToExistingParticle(int indexPart, int idTri, double *&history);

  class TriMesh *mesh_;
  class FixNeighlistMesh *fix_neighlist_mesh_;
  class FixPropertyAtom* fix_nneighs_;
  bool build_neighlist_;
  double *swap_;
  int numpages_;

  void sort_contacts();
  void swap(int ilocal,int ineigh, int jneigh, bool keepflag_swap);
};

// *************************************
#include "fix_contact_history_mesh_I.h"
// *************************************

}

#endif
#endif

/* ERROR/WARNING messages:

E: Pair style granular with history requires atoms have IDs

Atoms in the simulation do not have IDs, so history effects
cannot be tracked by the granular pair potential.

E: Too many touching neighbors - boost MAXTOUCH

A granular simulation has too many neighbors touching one atom.  The
MAXTOUCH parameter in fix_shear_history.cpp must be set larger and
LAMMPS must be re-built.

*/
