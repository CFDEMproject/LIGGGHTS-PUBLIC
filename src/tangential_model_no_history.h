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
#ifdef TANGENTIAL_MODEL
TANGENTIAL_MODEL(TANGENTIAL_NO_HISTORY,no_history,0)
#else
#ifndef TANGENTIAL_MODEL_NO_HISTORY_H_
#define TANGENTIAL_MODEL_NO_HISTORY_H_
#include "contact_models.h"
#include <algorithm>
#include "math.h"
#include "global_properties.h"

namespace LIGGGHTS {
namespace ContactModels
{
  using namespace std;

  template<>
  class TangentialModel<TANGENTIAL_NO_HISTORY> : protected Pointers
  {
    double ** coeffFrict;

  public:
    static const int MASK = CM_CONNECT_TO_PROPERTIES | CM_COLLISION;

    TangentialModel(LAMMPS * lmp, IContactHistorySetup*) : Pointers(lmp), coeffFrict(NULL)
    {
      
    }

    inline void registerSettings(Settings&){}

    inline void connectToProperties(PropertyRegistry & registry){
      registry.registerProperty("coeffFrict", &MODEL_PARAMS::createCoeffFrict);
      registry.connect("coeffFrict", coeffFrict,"tangential_model history");
    }

    inline void collision(const CollisionData & cdata, ForceData & i_forces, ForceData & j_forces) {
      const double xmu = coeffFrict[cdata.itype][cdata.jtype];
      const double enx = cdata.en[0];
      const double eny = cdata.en[1];
      const double enz = cdata.en[2];
      const double vrel = sqrt(cdata.vtr1*cdata.vtr1 + cdata.vtr2*cdata.vtr2 + cdata.vtr3*cdata.vtr3);

      // force normalization
      const double Ft_friction = xmu * fabs(cdata.Fn);
      const double Ft_damping = cdata.gammat*vrel;     
      double Ft;

      if (vrel != 0.0) Ft = min(Ft_friction, Ft_damping) / vrel;
      else Ft = 0.0;

      // tangential force due to tangential velocity damping

      const double Ft1 = -Ft*cdata.vtr1;
      const double Ft2 = -Ft*cdata.vtr2;
      const double Ft3 = -Ft*cdata.vtr3;

      // forces & torques

      const double tor1 = (eny*Ft3 - enz*Ft2);
      const double tor2 = (enz*Ft1 - enx*Ft3);
      const double tor3 = (enx*Ft2 - eny*Ft1);

      // return resulting forces
      if(cdata.is_wall) {
        const double area_ratio = cdata.area_ratio;
        i_forces.delta_F[0] += Ft1 * area_ratio;
        i_forces.delta_F[1] += Ft2 * area_ratio;
        i_forces.delta_F[2] += Ft3 * area_ratio;
        i_forces.delta_torque[0] = -cdata.cri * tor1 * area_ratio;
        i_forces.delta_torque[1] = -cdata.cri * tor2 * area_ratio;
        i_forces.delta_torque[2] = -cdata.cri * tor3 * area_ratio;
      } else {
        i_forces.delta_F[0] += Ft1;
        i_forces.delta_F[1] += Ft2;
        i_forces.delta_F[2] += Ft3;
        i_forces.delta_torque[0] = -cdata.cri * tor1;
        i_forces.delta_torque[1] = -cdata.cri * tor2;
        i_forces.delta_torque[2] = -cdata.cri * tor3;

        j_forces.delta_F[0] -= Ft1;
        j_forces.delta_F[1] -= Ft2;
        j_forces.delta_F[2] -= Ft3;
        j_forces.delta_torque[0] = -cdata.crj * tor1;
        j_forces.delta_torque[1] = -cdata.crj * tor2;
        j_forces.delta_torque[2] = -cdata.crj * tor3;
      }
    }

    inline void beginPass(CollisionData&, ForceData&, ForceData&){}
    inline void endPass(CollisionData&, ForceData&, ForceData&){}
    inline void noCollision(ContactData&, ForceData&, ForceData&){}
  };
}
}
#endif // TANGENTIAL_MODEL_NO_HISTORY_H_
#endif
