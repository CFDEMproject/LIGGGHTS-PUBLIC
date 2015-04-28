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

    Copyright 2013-     TU Graz

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

FixStyle(fiber/spring/simple,FixFiberSpringSimple)

#else

#ifndef LMP_FIX_FIBER_SPRING_SIMPLE_H
#define LMP_FIX_FIBER_SPRING_SIMPLE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixFiberSpringSimple : public Fix {
 public:
  FixFiberSpringSimple(class LAMMPS *, int, char **);
  ~FixFiberSpringSimple();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();
  double compute_vector(int);

 private:
  double xc,yc,zc,r0;
  double k_spring;
  int xflag,yflag,zflag;
  int styleflag;
  char *group2;
  int igroup2,group2bit;
  double masstotal,masstotal2;
  int nlevels_respa;
  double espring,ftotal[4];
  int force_flag;

  class PairGran *pair_gran;

  class FixPropertyGlobal* Y1; //Youngs Modulus
  class FixPropertyGlobal* v1; //Poisson's ratio
  class FixPropertyGlobal* charVel1; //characteristic velocity needed for Linear Spring Model

  double **Yeff,**Geff,**betaeff,**veff,**cohEnergyDens,**coeffRestLog,**coeffFrict,charVel,**coeffRollFrict;

  inline void deriveNormalStiffness(int &ip, int &jp,double &meff, double &kn);

};

}

#endif
#endif

