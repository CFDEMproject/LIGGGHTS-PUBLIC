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

FixStyle(transportequation/scalar,FixScalarTransportEquation)

#else

#ifndef LMP_FIX_SCALAR_TRANSPORT_EQUATION_H
#define LMP_FIX_SCALAR_TRANSPORT_EQUATION_H

#include "fix.h"

namespace LAMMPS_NS {

class FixScalarTransportEquation : public Fix {
 public:
  FixScalarTransportEquation(class LAMMPS *, int, char **);
  ~FixScalarTransportEquation();

  virtual int setmask();
  virtual void post_create();
  virtual void pre_delete(bool unfixflag);
  virtual void init();
  virtual int modify_param(int narg, char **arg);
  virtual void updatePtrs();
  virtual void initial_integrate_respa(int,int,int);
  virtual void initial_integrate(int);
  virtual void pre_force(int vflag);
  virtual void final_integrate();
  virtual double compute_scalar();
  bool match_equation_id(const char*);

  double *get_capacity();

 protected:

  int nlevels_respa;

  char *equation_id;

  class FixPropertyAtom* fix_quantity;
  char *quantity_name;
  class FixPropertyAtom* fix_flux;
  char *flux_name;
  class FixPropertyAtom* fix_source;
  char *source_name;

  //storage capacity - would be thermal capacity for heat conduction
  int capacity_flag;
  class FixPropertyGlobal* fix_capacity;
  double *capacity;
  char *capacity_name;

  double quantity_0;  
  double *quantity;   
  double *flux;       
  double *source;     

  // flag if integrate quantity or not
  bool int_flag;

  int  nevery_; //integrate only this many time steps (to avoid round-off issues)
  bool performedIntegrationLastStep_;
};

}

#endif
#endif
