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
   Contributing authors for original version: Leo Silbert (SNL), Gary Grest (SNL)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

#else

#ifndef LMP_PAIR_GRAN_H
#define LMP_PAIR_GRAN_H

#include "pair.h"
#include "compute_pair_gran_local.h"
#include "contact_interface.h"
#include <vector>
#include <string>

namespace LCM = LIGGGHTS::ContactModels;

namespace LAMMPS_NS {

class PairGran : public Pair, public LIGGGHTS::IContactHistorySetup {
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
  double memory_usage();

  virtual int64_t hashcode() = 0;

  int  cplenable()
  { return cpl_enable; }

  void register_compute_pair_local(class ComputePairGranLocal *,int&);
  void unregister_compute_pair_local(class ComputePairGranLocal *ptr);

  inline void cpl_add_pair(LCM::CollisionData & cdata, LCM::ForceData & i_forces)
  {
    const double fx = i_forces.delta_F[0];
    const double fy = i_forces.delta_F[1];
    const double fz = i_forces.delta_F[2];
    const double tor1 = i_forces.delta_torque[0];
    const double tor2 = i_forces.delta_torque[1];
    const double tor3 = i_forces.delta_torque[2];
    cpl_->add_pair(cdata.i, cdata.j, fx,fy,fz,tor1,tor2,tor3,cdata.contact_history);
  }

  /* PUBLIC ACCESS FUNCTIONS */

  int is_history()
  { return history; }

  int dnum_pair()
  { return dnum_pairgran; }

  inline int dnum()
  { return dnum_all; }

  inline class ComputePairGranLocal * cpl() {
    return cpl_;
  }

  inline bool storeContactForces() {
    return store_contact_forces_;
  }

  inline int freeze_group_bit() const {
    return freeze_group_bit_;
  }

  inline int computeflag() const {
    return computeflag_;
  }

  inline int shearupdate() const {
    return shearupdate_;
  }

  class FixContactPropertyAtom * fix_contact_forces() {
    return fix_contact_forces_;
  }

  class FixRigid* fr_pair()
  { return fix_rigid; }

  double * mr_pair()
  { return mass_rigid; }

  virtual double stressStrainExponent() = 0;

  class Properties* get_properties()
  {return properties; }

  int fix_extra_dnum_index(class Fix *fix);

  void *extract(const char *str, int &dim);

  int add_history_value(std::string name, std::string newtonflag) {
    int offset = history_arg.size();
    history = true;
    history_arg.push_back(HistoryArg(name, newtonflag));
    dnum_pairgran++;
    return offset;
  }

  void do_store_contact_forces()
  { store_contact_forces_ = true; }

 protected:

  struct HistoryArg {
    std::string name;
    std::string newtonflag;

    HistoryArg(std::string name, std::string newtonflag) : name(name), newtonflag(newtonflag) {}
  };

  std::vector<HistoryArg> history_arg;

  void history_args(char ** args) {
    for(size_t i = 0; i < history_arg.size(); i++) {
      args[2*i] = (char*)history_arg[i].name.c_str();
      args[2*i+1] = (char*)history_arg[i].newtonflag.c_str();
    }
  }

  virtual void compute_force(int eflag, int vflag,int addflag) = 0;
  virtual bool forceoff();

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
  class ComputePairGranLocal *cpl_;

  // storage for per-contact forces
  bool store_contact_forces_;
  class FixContactPropertyAtom *fix_contact_forces_;

  // storage of rigid body masses for use in granular interactions

  class FixRigid *fix_rigid; // ptr to rigid body fix, NULL if none
  double *mass_rigid;        // rigid mass for owned+ghost atoms
  int nmax;                  // allocated size of mass_rigid

  double dt;
  int freeze_group_bit_;

  // contact history
  int history;
  int dnum_pairgran;
  class FixContactHistory *fix_history;
  int shearupdate_;

  int computeflag_;

  double *onerad_dynamic,*onerad_frozen;
  double *maxrad_dynamic,*maxrad_frozen;

  bool needs_neighlist;

  void allocate();

 private:

  class Properties *properties;

  // shear history
  int dnum_all;

  int nfix;
  Fix **fix_dnum;
  int *dnum_index;
};

}

#endif
#endif
