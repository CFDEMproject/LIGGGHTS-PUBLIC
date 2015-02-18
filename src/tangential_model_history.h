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
TANGENTIAL_MODEL(TANGENTIAL_HISTORY,history,1)
#else
#ifndef TANGENTIAL_MODEL_HISTORY_H_
#define TANGENTIAL_MODEL_HISTORY_H_
#include "contact_models.h"
#include "math.h"
#include "update.h"
#include "global_properties.h"
#include "atom.h"

namespace LIGGGHTS {
namespace ContactModels
{
  template<>
  class TangentialModel<TANGENTIAL_HISTORY> : protected Pointers
  {
    double ** coeffFrict;
    int history_offset;

  public:
    static const int MASK = CM_CONNECT_TO_PROPERTIES | CM_COLLISION | CM_NO_COLLISION;

    TangentialModel(LAMMPS * lmp, IContactHistorySetup * hsetup) : Pointers(lmp),
      coeffFrict(NULL)
    {
      history_offset = hsetup->add_history_value("shearx", "1");
      hsetup->add_history_value("sheary", "1");
      hsetup->add_history_value("shearz", "1");

    }

    inline void registerSettings(Settings&){}

    inline void connectToProperties(PropertyRegistry & registry)
    {
      registry.registerProperty("coeffFrict", &MODEL_PARAMS::createCoeffFrict);
      registry.connect("coeffFrict", coeffFrict,"tangential_model history");
    }

    inline void collision(const CollisionData & cdata, ForceData & i_forces, ForceData & j_forces)
    {
      // normal forces = Hookian contact + normal velocity damping
      const double enx = cdata.en[0];
      const double eny = cdata.en[1];
      const double enz = cdata.en[2];

      // shear history effects
      if(cdata.touch) *cdata.touch |= TOUCH_TANGENTIAL_MODEL;
      double * const shear = &cdata.contact_history[history_offset];

      if (cdata.shearupdate && cdata.computeflag) {
        const double dt = update->dt;
        shear[0] += cdata.vtr1 * dt;
        shear[1] += cdata.vtr2 * dt;
        shear[2] += cdata.vtr3 * dt;

        // rotate shear displacements

        double rsht = shear[0]*enx + shear[1]*eny + shear[2]*enz;
        shear[0] -= rsht * enx;
        shear[1] -= rsht * eny;
        shear[2] -= rsht * enz;
      }

      const double shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] + shear[2]*shear[2]);
      const double kt = cdata.kt;
      const double xmu = coeffFrict[cdata.itype][cdata.jtype];

      // tangential forces = shear + tangential velocity damping
      double Ft1 = -(kt * shear[0]);
      double Ft2 = -(kt * shear[1]);
      double Ft3 = -(kt * shear[2]);

      // rescale frictional displacements and forces if needed
      const double Ft_shear = kt * shrmag; // sqrt(Ft1 * Ft1 + Ft2 * Ft2 + Ft3 * Ft3);
      const double Ft_friction = xmu * fabs(cdata.Fn);

      // energy loss from sliding or damping
      if (Ft_shear > Ft_friction) {
        if (shrmag != 0.0) {
          const double ratio = Ft_friction / Ft_shear;
          Ft1 *= ratio;
          Ft2 *= ratio;
          Ft3 *= ratio;
          shear[0] = -Ft1/kt;
          shear[1] = -Ft2/kt;
          shear[2] = -Ft3/kt;
        }
        else Ft1 = Ft2 = Ft3 = 0.0;
      }
      else
      {
        const double gammat = cdata.gammat;
        Ft1 -= (gammat*cdata.vtr1);
        Ft2 -= (gammat*cdata.vtr2);
        Ft3 -= (gammat*cdata.vtr3);
      }

      // forces & torques
      const double tor1 = eny * Ft3 - enz * Ft2;
      const double tor2 = enz * Ft1 - enx * Ft3;
      const double tor3 = enx * Ft2 - eny * Ft1;

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

        j_forces.delta_F[0] += -Ft1;
        j_forces.delta_F[1] += -Ft2;
        j_forces.delta_F[2] += -Ft3;
        j_forces.delta_torque[0] = -cdata.crj * tor1;
        j_forces.delta_torque[1] = -cdata.crj * tor2;
        j_forces.delta_torque[2] = -cdata.crj * tor3;
      }
    }

    inline void noCollision(ContactData & cdata, ForceData&, ForceData&)
    {
      // unset non-touching neighbors
      // TODO even if shearupdate == false?
      if(cdata.touch) *cdata.touch &= ~TOUCH_TANGENTIAL_MODEL;
      double * const shear = &cdata.contact_history[history_offset];
      shear[0] = 0.0;
      shear[1] = 0.0;
      shear[2] = 0.0;
    }

    inline void beginPass(CollisionData&, ForceData&, ForceData&){}
    inline void endPass(CollisionData&, ForceData&, ForceData&){}
  };
}
}
#endif // TANGENTIAL_MODEL_HISTORY_H_
#endif
