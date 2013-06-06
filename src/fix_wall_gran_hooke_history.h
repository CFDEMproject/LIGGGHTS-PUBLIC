/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   This file was modified with respect to the release in LAMMPS
   Modifications are Copyright 2009-2012 JKU Linz
                     Copyright 2012-     DCS Computing GmbH, Linz

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors for original version: Leo Silbert (SNL), Gary Grest (SNL)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(wall/gran/hooke/history,FixWallGranHookeHistory)

#else

#ifndef LMP_FIX_WALL_GRAN_HOOKE_HISTORY_H
#define LMP_FIX_WALL_GRAN_HOOKE_HISTORY_H

#include "fix_wall_gran.h"

namespace LAMMPS_NS {

class FixWallGranHookeHistory : public FixWallGran {
 public:
  FixWallGranHookeHistory(class LAMMPS *, int, char **);
  ~FixWallGranHookeHistory();

 protected:
  virtual void post_create();
  virtual void init_granular();
  virtual void init_heattransfer();

  void addHeatFlux(TriMesh *mesh,int ip, double rsq, double area_ratio);
  virtual void compute_force(int ip, double deltan, double rsq,double meff_wall,
                              double dx, double dy, double dz,double *vwall,
                             double *c_history,double area_ratio);
  virtual void addCohesionForce(int &ip, double &r, double &Fn_coh,double area_ratio);
  template <int ROLLINGFRICTION>
  void addRollingFrictionTorque(int ip, double wr1,double wr2,double wr3,double cr,double ccel,
            double r,double mi,double rmu,double kn,double kt,double dx, double dy, double dz,double rsqinv,double *c_history,double *r_torque);

  virtual void deriveContactModelParams(int ip, double deltan,double meff_wall,
                            double &kn, double &kt, double &gamman, double &gammat,
                            double &xmu,double &rmu,double &vnnr);
  virtual void pre_reset_history(int,double*) {}

  int dampflag,cohesionflag,rollingflag,viscousflag;
  double **Yeff,**Geff,**betaeff,**veff,**cohEnergyDens,**coeffRestLog,**coeffFrict,**coeffRollVisc;
  double charVel,**coeffRollFrict,**coeffMu,**coeffRestMax,**coeffStc;

  // heat transfer

  class FixPropertyAtom *fppa_T;
  class FixPropertyAtom *fppa_hf;

  double Temp_wall;
  double Q,Q_add;

  const double *th_cond;
  double const* const* deltan_ratio;
};

}

#endif
#endif
