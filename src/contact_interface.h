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
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */
#ifndef CONTACT_INTERFACE_H_
#define CONTACT_INTERFACE_H_

#include <string>

namespace LIGGGHTS {
namespace ContactModels {

// data available in noCollision() and collision()

struct ContactData {
  double radi;
  double radj;
  double radsum;
  double rsq;
  double delta[3];

  double area_ratio;

  int * touch;
  double * contact_history;

  int i;
  int j;

  bool is_wall;
  bool has_force_update;

  ContactData() : area_ratio(1.0) {}
};

// data available in collision() only

struct CollisionData: ContactData {
  double r;
  double rinv;
  double en[3];
  double * v_i;
  double * v_j;
  double * omega_i;
  double * omega_j;

  double kt;
  double kn;
  double gammat;
  double gamman;

  double Fn;
  double Ft;

  double vn;
  double deltan;
  double cri;
  double crj;
  double wr1;
  double wr2;
  double wr3;

  double vtr1;
  double vtr2;
  double vtr3;

  double mi;
  double mj;
  double meff;

  int computeflag;
  int shearupdate;
  int itype;
  int jtype;

  CollisionData() : Fn(0.0), Ft(0.0) {}
};

struct ForceData {
  double delta_F[3];       // total force acting on particle
  double delta_torque[3];  // torque acting on a particle

  ForceData()
  {
    reset();
  }

  inline void reset() {
    delta_F[0] = 0.0;
    delta_F[1] = 0.0;
    delta_F[2] = 0.0;
    delta_torque[0] = 0.0;
    delta_torque[1] = 0.0;
    delta_torque[2] = 0.0;
  }
};
}

class IContactHistorySetup {
public:
  virtual int add_history_value(std::string name, std::string newtonflag) = 0;
};

}

#endif /* CONTACT_INTERFACE_H_ */
