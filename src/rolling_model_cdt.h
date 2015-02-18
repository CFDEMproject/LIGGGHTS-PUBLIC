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
#ifdef ROLLING_MODEL
ROLLING_MODEL(ROLLING_CDT,cdt,1)
#else
#ifndef ROLLING_MODEL_CDT_H_
#define ROLLING_MODEL_CDT_H_
#include "contact_models.h"
#include <algorithm>
#include "math.h"
#include "math_extra_liggghts.h"

namespace LIGGGHTS {
namespace ContactModels
{
  using namespace LAMMPS_NS;

  template<>
  class RollingModel<ROLLING_CDT> : protected Pointers {
  public:
    static const int MASK = CM_CONNECT_TO_PROPERTIES | CM_COLLISION;

    RollingModel(LAMMPS * lmp, IContactHistorySetup*) : Pointers(lmp), coeffRollFrict(NULL)
    {
      
    }

    void registerSettings(Settings&) {}

    void connectToProperties(PropertyRegistry & registry)
    {
      registry.registerProperty("coeffRollFrict", &MODEL_PARAMS::createCoeffRollFrict);
      registry.connect("coeffRollFrict", coeffRollFrict,"rolling_model cdt");

      // error checks on coarsegraining
      if(force->cg_active())
        error->cg(FLERR,"rolling model cdt");
    }

    void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces) 
    {
      const double rmu= coeffRollFrict[cdata.itype][cdata.jtype];

      double r_torque[3], wr_roll[3];
      vectorZeroize3D(r_torque);

      if(cdata.is_wall){
        const double wr1 = cdata.wr1;
        const double wr2 = cdata.wr2;
        const double wr3 = cdata.wr3;
        const double wrmag = sqrt(wr1*wr1+wr2*wr2+wr3*wr3);
        if (wrmag > 0.)
        {
          const double radius = cdata.radi;
          const double kn = cdata.kn;
          const double enx = cdata.en[0];
          const double eny = cdata.en[1];
          const double enz = cdata.en[2];

          r_torque[0] = rmu*kn*(radius-cdata.r)*wr1/wrmag*cdata.cri;
          r_torque[1] = rmu*kn*(radius-cdata.r)*wr2/wrmag*cdata.cri;
          r_torque[2] = rmu*kn*(radius-cdata.r)*wr3/wrmag*cdata.cri;

          // remove normal (torsion) part of torque
          double rtorque_dot_delta = r_torque[0]*enx+ r_torque[1]*eny + r_torque[2]*enz;
          double r_torque_n[3];
          r_torque_n[0] = enx * rtorque_dot_delta;
          r_torque_n[1] = eny * rtorque_dot_delta;
          r_torque_n[2] = enz * rtorque_dot_delta;
          vectorSubtract3D(r_torque,r_torque_n,r_torque);
        }
      } else {
        vectorSubtract3D(atom->omega[cdata.i],atom->omega[cdata.j],wr_roll);
        const double wr_rollmag = vectorMag3D(wr_roll);

        if(wr_rollmag > 0.)
        {
          const double radi = cdata.radi;
          const double radj = cdata.radj;
          const double enx = cdata.en[0];
          const double eny = cdata.en[1];
          const double enz = cdata.en[2];

          // calculate torque
          const double reff= cdata.is_wall ? radi : (radi*radj/(radi+radj));
          vectorScalarMult3D(wr_roll,rmu*cdata.kn*cdata.deltan*reff/wr_rollmag,r_torque);

          // remove normal (torsion) part of torque
          const double rtorque_dot_delta = r_torque[0]*enx + r_torque[1]*eny + r_torque[2]*enz;
          double r_torque_n[3];
          r_torque_n[0] = enx * rtorque_dot_delta;
          r_torque_n[1] = eny * rtorque_dot_delta;
          r_torque_n[2] = enz * rtorque_dot_delta;
          vectorSubtract3D(r_torque,r_torque_n,r_torque);
        }
      }

      i_forces.delta_torque[0] -= r_torque[0];
      i_forces.delta_torque[1] -= r_torque[1];
      i_forces.delta_torque[2] -= r_torque[2];

      j_forces.delta_torque[0] += r_torque[0];
      j_forces.delta_torque[1] += r_torque[1];
      j_forces.delta_torque[2] += r_torque[2];
    }

    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}
    void noCollision(ContactData&, ForceData&, ForceData&){}

  private:
    double ** coeffRollFrict;
  };
}
}
#endif // ROLLING_MODEL_CDT_H_
#endif
