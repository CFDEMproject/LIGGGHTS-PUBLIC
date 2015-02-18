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

FixStyle(couple/cfd,FixCfdCoupling)

#else

#ifndef LMP_FIX_CFD_COUPLING_H
#define LMP_FIX_CFD_COUPLING_H

#include "fix.h"

namespace LAMMPS_NS {

class FixCfdCoupling : public Fix {
 friend class CfdRegionmodel;
 friend class CfdDatacoupling;

 public:

  FixCfdCoupling(class LAMMPS *, int, char **);
  ~FixCfdCoupling();
  void post_create();

  int setmask();
  void init();
  virtual void setup(int);
  virtual void min_setup(int);
  void end_of_step();
  void post_force_respa(int vflag, int ilevel, int iloop);
  void min_post_force(int);

  // pushing and pulling of properties
  //void pull(char *name,char *type,void *&ptr);
  //void push(char *name,char *type,void *&ptr);
  void add_push_property(const char *name, const char *type);
  void add_pull_property(const char *name, const char *type);
  void check_datatransfer();

  int coupleThis() {return couple_this_;}

  class CfdDatacoupling* get_dc(){return dc_;}
  
 protected:

  int iarg_;

  // data transfer is handled by this class
  class CfdDatacoupling *dc_;

 private:

  int couple_this_;

  // couple every couple_nevery_ timesteps
  // not used in case of MPI coupling
  int couple_nevery_,ts_create_;

  // regionmodels
  class CfdRegionmodel *rm_;

  int nlevels_respa;
};

}

#endif
#endif
