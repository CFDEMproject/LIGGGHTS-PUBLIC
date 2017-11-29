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
    Arno Mayrhofer (CFDEMresearch GmbH, Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
    Copyright 2016-     CFDEMresearch GmbH, Linz
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(contacthistory/mesh,FixContactHistoryMesh) 

#else

#ifndef LMP_FIX_CONTACT_HISTORY_MESH_H
#define LMP_FIX_CONTACT_HISTORY_MESH_H

#include "fix_contact_history.h"
#include "fix_property_atom.h"
#include "my_page.h"
#include <cmath>
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

  bool handleContact(int iPart, int idTri, double *&history, bool intersectflag, bool faceflag);
  void markAllContacts();
  void cleanUpContacts();
  void cleanUpContactJumps();
  
  // OMP interface
  void cleanUpContacts(int ifrom, int ito);
  void reset_history();

  // return # of contacts
  int n_contacts(int & nIntersect);
  int n_contacts(int contact_groupbit, int & nIntersect);

  int get_partner_idTri(const int i, const int j) const
  { return partner_[i][j]; }

  int nneighs(const int iP) const
  { return fix_nneighs_->get_vector_atom_int(iP); }

  int get_contact(const int i, const int j);

 protected:

  MyPage<int> *ipage1_;        // pages of neighbor tri IDs
  MyPage<double> *dpage1_;     // pages of contact history with neighbors
  MyPage<int> *ipage2_;        // pages of neighbor tri IDs
  MyPage<double> *dpage2_;     // pages of contact history with neighbors
  MyPage<bool> ** keeppage_;   // pages of deletion flags with neighbors
  MyPage<bool> ** intersectpage_;   // pages of deletion flags with neighbors

  bool **keepflag_;  
                     
  bool **intersectflag_; 

  void allocate_pages();

 private:

  // functions specific for mesh - contact management
  bool haveContact(int indexPart, int idTri, double *&history, bool intersectflag);
  bool coplanarContactAlready(int indexPart, int idTri);
  void checkCoplanarContactHistory(int indexPart, int idTri, double *&history);
  void addNewTriContactToExistingParticle(int indexPart, int idTri, double *&history, bool intersectflag);

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
