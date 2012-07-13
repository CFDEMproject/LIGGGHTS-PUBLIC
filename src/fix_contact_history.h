/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(contacthistory,FixContactHistory) 

FixStyle(contacthistory/mesh,FixContactHistory) 

#else

#ifndef LMP_FIX_CONTACT_HISTORY_H
#define LMP_FIX_CONTACT_HISTORY_H

#include "fix.h"
#include "tri_mesh.h"
#include "atom.h"

namespace LAMMPS_NS {

class FixContactHistory : public Fix {
  friend class Neighbor;
  friend class PairGran;

 public:

  FixContactHistory(class LAMMPS *, int, char **);
  ~FixContactHistory();

  // inherited from Fix

  void post_create();
  int setmask();
  void init();

  void setup_pre_exchange();
  void pre_exchange();

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int);
  void set_arrays(int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();
  void write_restart(FILE *);
  void restart(char *);

  // spefific interface for mesh

  void handleContact(int iPart, int idTri, double *&history);
  void handleNoContact(int iPart, int idTri);
  void cleanUpContacts();

  // return # of contacts
  int n_contacts();
  int n_contacts(int contact_groupbit);

 private:

  // functions specific for pair

  void pre_exchange_pair();

  // functions specific for mesh - contact management

  bool haveContact(int indexPart, int idTri, double *&history);
  bool checkCoplanarContact(int indexPart, int idTri, double *&history);
  void addNewTriContactToExistingParticle(int indexPart, int idTri, double *&history);

  // mem management
  void check_grow();

  // data members

  int *npartner;                // # of touching partners of each atom
  int **partner;                // tags for the partners
  double ***contacthistory;     // history values with the partner
  int maxtouch;                 // max number of partners per atom
  void grow_arrays_maxtouch(int);

  bool **delflag;               

  bool is_pair;
  class PairGran *pair_gran;
  class TriMesh *mesh_;

  int dnum;
  int *newtonflag;
  char **history_id;
  int laststep;
};

// *************************************
#include "fix_contact_history_I.h"
// *************************************

}

#endif
#endif

/* ERROR/WARNING messages:

E: Pair style granular with history requires atoms have IDs

Atoms in the simulation do not have IDs, so history effects
cannot be tracked by the granular pair potential.

*/
