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

PairStyle(gran/hooke/history,PairGranHookeHistory)

#else

#ifndef LMP_PAIR_GRAN_HOOKE_HISTORY_H
#define LMP_PAIR_GRAN_HOOKE_HISTORY_H

#include "pair_gran.h"

namespace LAMMPS_NS {

class PairGranHookeHistory : public PairGran {

 friend class FixWallGranHookeHistory;
 friend class FixCheckTimestepGran;

 public:

  PairGranHookeHistory(class LAMMPS *);
  ~PairGranHookeHistory();

  virtual void settings(int, char **);
  virtual void init_granular(); 

  virtual void compute(int, int,int);

  virtual void write_restart_settings(FILE *);
  virtual void read_restart_settings(FILE *);

 protected:

  virtual void history_args(char**);
  void allocate_properties(int);

  class FixPropertyGlobal* Y1; //Youngs Modulus
  class FixPropertyGlobal* v1; //Poisson's ratio
  class FixPropertyGlobal* cohEnergyDens1; //Cohesion energy density

  class FixPropertyGlobal* coeffRest1; //coefficient of restitution
  class FixPropertyGlobal* coeffFrict1; //coefficient of (static) friction
  class FixPropertyGlobal* coeffRollFrict1; //characteristic velocity needed for Linear Spring Model

  int charVelflag;
  class FixPropertyGlobal* charVel1; //characteristic velocity needed for Linear Spring Model

  double **Yeff,**Geff,**betaeff,**veff,**cohEnergyDens,**coeffRestLog,**coeffFrict,charVel,**coeffRollFrict;

  virtual void deriveContactModelParams(int &, int &,double &, double &, double &,double &, double &, double &, double &,double &);
  virtual void addCohesionForce(int &, int &,double &,double &);

  int cohesionflag; 
  int dampflag,rollingflag; 
};

}

#endif
#endif
