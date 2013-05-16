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

/* ----------------------------------------------------------------------
   Contributing authors for original version: Leo Silbert (SNL), Gary Grest (SNL)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

#else

#ifndef LMP_PAIR_GRAN_H
#define LMP_PAIR_GRAN_H

#include "pair.h"

namespace LAMMPS_NS {

class PairGran : public Pair {
 public:

  friend class FixWallGran;
  friend class FixCheckTimestepGran;

  PairGran(class LAMMPS *);
  ~PairGran();

  /* INHERITED FROM Pair */

  virtual void compute(int eflag, int vflag);
  virtual void compute_pgl(int eflag, int vflag);
  virtual void settings(int, char **) = 0;
  virtual void coeff(int, char **);
  virtual void init_style();
  virtual void init_granular() {} 
  virtual void init_list(int, class NeighList *);
  virtual double init_one(int, int);
  int pack_comm(int n, int *list,double *buf, int pbc_flag, int *pbc);
  void unpack_comm(int n, int first, double *buf);
  virtual void write_restart(FILE *);
  virtual void read_restart(FILE *);
  virtual void write_restart_settings(FILE *){}
  virtual void read_restart_settings(FILE *){}
  virtual void reset_dt();

  int  cplenable()
  { return cpl_enable; }
  void register_compute_pair_local(class ComputePairGranLocal *,int&);
  void unregister_compute_pair_local(class ComputePairGranLocal *ptr);

  /* PUBLIC ACCESS FUNCTIONS */

  int is_history()
  { return history; }

  int dnum_pair()
  { return dnum_pairgran; }

  class FixRigid* fr_pair()
  { return fix_rigid; }

  class MechParamGran *mpg;

  int fix_extra_dnum_index(class Fix *fix);

  void *extract(const char *str, int &dim);

 protected:

  virtual void history_args(char**) =0;
  virtual void compute_force(int eflag, int vflag,int addflag) = 0;
  virtual bool forceoff() = 0;
  inline int dnum()
  { return dnum_all; }

  virtual void updatePtrs();

  // for parsing settings() args
  int iarg_;

  char * suffix;
  int neighprev;

  // stuff for tracking energy
  int energytrack_enable;
  class FixPropertyAtom* fppaCPEn; //collision potential energy normal
  class FixPropertyAtom* fppaCDEn; //collision dissipation energy normal
  class FixPropertyAtom* fppaCPEt; //collision potential energy tang
  class FixPropertyAtom* fppaCDEVt; //collision dissipation energy viscous tang
  class FixPropertyAtom* fppaCDEFt; //collision dissipation energy friction tang
  class FixPropertyAtom* fppaCTFW; //collision tangential force work
  class FixPropertyAtom* fppaDEH; //dissipation energy of history term (viscous and friction, accumulated over time)
  double *CPEn, *CDEn, *CPEt, *CDEVt, *CDEFt, *CTFW, *DEH;

  // stuff for compute pair gran local
  int cpl_enable;
  class ComputePairGranLocal *cpl;

  class FixRigid* fix_rigid;
  int *body;
  double *masstotal;

  double dt;
  int freeze_group_bit;

  // contact history
  int history;
  int dnum_pairgran;
  class FixContactHistory *fix_history;
  int shearupdate;

  int computeflag;

  double *onerad_dynamic,*onerad_frozen;
  double *maxrad_dynamic,*maxrad_frozen;

  bool needs_neighlist;

  void allocate();

 private:

  // shear history
  int dnum_all;

  int nfix;
  Fix **fix_dnum;
  int *dnum_index;
};

}

#endif
#endif
