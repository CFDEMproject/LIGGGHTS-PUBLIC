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
    This file is from LAMMPS, but has been modified. Copyright for
    modification:

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz

    Copyright of original file:
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(contacthistory,FixContactHistory) 

#else

#ifndef LMP_FIX_CONTACT_HISTORY_H
#define LMP_FIX_CONTACT_HISTORY_H

#include "fix.h"
#include "my_page.h"
#include "vector_liggghts.h"

namespace LAMMPS_NS {

class FixContactHistory : public Fix {
  friend class Neighbor;
  friend class PairGran;

 public:
  FixContactHistory(class LAMMPS *, int, char **);
  ~FixContactHistory();
  virtual void post_create() {}
  virtual int setmask();
  virtual void init();
  virtual void setup_pre_exchange();
  virtual void setup_pre_neighbor() {}
  virtual void pre_exchange();
  virtual void min_setup_pre_exchange();
  void min_pre_exchange();

  virtual double memory_usage();
  virtual void grow_arrays(int);
  virtual void copy_arrays(int, int, int);
  void set_arrays(int);
  int pack_exchange(int, double *);
  virtual int unpack_exchange(int, double *);
  virtual void write_restart(FILE *fp);
  void restart(char *buf);
  int pack_restart(int, double *);
  virtual void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();

  // inline access
  inline int n_partner(int i)
  { return npartner_[i]; }

  inline int partner(int i,int j)
  { return partner_[i][j]; }

  inline void contacthistory(int i,int j,double *h)
  { vectorCopyN(&(contacthistory_[i][j*dnum_]),h,dnum_); }

  inline double* contacthistory(int i,int j)
  { return &(contacthistory_[i][j*dnum_]); }

 protected:

  int iarg_;

  int dnum_;                      
  char *variablename_;
  int *newtonflag_;
  char **history_id_;
  int index_decide_noncontacting_;

  int *npartner_;                // # of touching partners of each atom
  int **partner_;                // tags for the partners
  double **contacthistory_;     // contact history values with the partner
  int maxtouch_;                 // max # of touching partners for my atoms

  class Pair *pair_gran_;
  int *computeflag_;             // computeflag in PairGranHookeHistory

  int pgsize_,oneatom_;          // copy of settings in Neighbor
  MyPage<int> *ipage_;           // pages of partner atom IDs
  MyPage<double> *dpage_;        // pages of shear history with partners

  virtual void allocate_pages();
};

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
