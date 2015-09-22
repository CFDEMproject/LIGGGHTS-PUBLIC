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

    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Richard Berger (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifdef TANGENTIAL_MODEL
TANGENTIAL_MODEL(TANGENTIAL_NO_HISTORY,no_history,1)
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
    static const int MASK = CM_CONNECT_TO_PROPERTIES | CM_SURFACES_INTERSECT;

    TangentialModel(LAMMPS * lmp, IContactHistorySetup*,class ContactModelBase *) : Pointers(lmp), coeffFrict(NULL)
    {
      
    }

    inline void registerSettings(Settings&){}

    inline void connectToProperties(PropertyRegistry & registry){
      registry.registerProperty("coeffFrict", &MODEL_PARAMS::createCoeffFrict);
      registry.connect("coeffFrict", coeffFrict,"tangential_model history");
    }

    inline void surfacesIntersect(const SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces) {
      const double xmu = coeffFrict[sidata.itype][sidata.jtype];
      const double enx = sidata.en[0];
      const double eny = sidata.en[1];
      const double enz = sidata.en[2];
      const double vrel = sqrt(sidata.vtr1*sidata.vtr1 + sidata.vtr2*sidata.vtr2 + sidata.vtr3*sidata.vtr3);

      // force normalization
      const double Ft_friction = xmu * fabs(sidata.Fn);
      const double Ft_damping = sidata.gammat*vrel;     
      double Ft;

      if (vrel != 0.0) Ft = min(Ft_friction, Ft_damping) / vrel;
      else Ft = 0.0;

      // tangential force due to tangential velocity damping

      const double Ft1 = -Ft*sidata.vtr1;
      const double Ft2 = -Ft*sidata.vtr2;
      const double Ft3 = -Ft*sidata.vtr3;

      // forces & torques

      const double tor1 = (eny*Ft3 - enz*Ft2);
      const double tor2 = (enz*Ft1 - enx*Ft3);
      const double tor3 = (enx*Ft2 - eny*Ft1);

      #ifdef SUPERQUADRIC_ACTIVE_FLAG
          double torque_i[3];
          if(sidata.is_non_spherical) {
            double xci[3];
            double Ft_i[3] = { Ft1,  Ft2,  Ft3 };
            vectorSubtract3D(sidata.contact_point, sidata.pos_i, xci);
            vectorCross3D(xci, Ft_i, torque_i);
          } else {
            torque_i[0] = -sidata.cri * tor1;
            torque_i[1] = -sidata.cri * tor2;
            torque_i[2] = -sidata.cri * tor3;
          }
      #endif
      // return resulting forces
      if(sidata.is_wall) {
        const double area_ratio = sidata.area_ratio;
        i_forces.delta_F[0] += Ft1 * area_ratio;
        i_forces.delta_F[1] += Ft2 * area_ratio;
        i_forces.delta_F[2] += Ft3 * area_ratio;
        #ifdef SUPERQUADRIC_ACTIVE_FLAG
                i_forces.delta_torque[0] += torque_i[0] * area_ratio;
                i_forces.delta_torque[1] += torque_i[1] * area_ratio;
                i_forces.delta_torque[2] += torque_i[2] * area_ratio;
        #else
                i_forces.delta_torque[0] = -sidata.cri * tor1 * area_ratio;
                i_forces.delta_torque[1] = -sidata.cri * tor2 * area_ratio;
                i_forces.delta_torque[2] = -sidata.cri * tor3 * area_ratio;
        #endif
      } else {
        i_forces.delta_F[0] += Ft1;
        i_forces.delta_F[1] += Ft2;
        i_forces.delta_F[2] += Ft3;
        j_forces.delta_F[0] -= Ft1;
        j_forces.delta_F[1] -= Ft2;
        j_forces.delta_F[2] -= Ft3;
        #ifdef SUPERQUADRIC_ACTIVE_FLAG
                double torque_j[3];
                if(sidata.is_non_spherical) {
                  double xcj[3];
                  vectorSubtract3D(sidata.contact_point, sidata.pos_j, xcj);
                  double Ft_j[3] = { -Ft1,  -Ft2,  -Ft3 };
                  vectorCross3D(xcj, Ft_j, torque_j);
                } else {
                  torque_j[0] = -sidata.crj * tor1;
                  torque_j[1] = -sidata.crj * tor2;
                  torque_j[2] = -sidata.crj * tor3;
                }
                i_forces.delta_torque[0] += torque_i[0];
                i_forces.delta_torque[1] += torque_i[1];
                i_forces.delta_torque[2] += torque_i[2];

                j_forces.delta_torque[0] += torque_j[0];
                j_forces.delta_torque[1] += torque_j[1];
                j_forces.delta_torque[2] += torque_j[2];
        #else
                i_forces.delta_torque[0] = -sidata.cri * tor1;
                i_forces.delta_torque[1] = -sidata.cri * tor2;
                i_forces.delta_torque[2] = -sidata.cri * tor3;

                j_forces.delta_torque[0] = -sidata.crj * tor1;
                j_forces.delta_torque[1] = -sidata.crj * tor2;
                j_forces.delta_torque[2] = -sidata.crj * tor3;
        #endif
      }
    }

    inline void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    inline void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    inline void surfacesClose(SurfacesCloseData&, ForceData&, ForceData&){}
  };
}
}
#endif // TANGENTIAL_MODEL_NO_HISTORY_H_
#endif
