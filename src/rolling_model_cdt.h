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

    Alexander Podlozhnyuk (DCS Computing GmbH, Linz)
    Andreas Aigner (DCS Computing GmbH, Linz)
    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Richard Berger (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifdef ROLLING_MODEL
ROLLING_MODEL(ROLLING_CDT,cdt,1)
#else
#ifndef ROLLING_MODEL_CDT_H_
#define ROLLING_MODEL_CDT_H_
#include "contact_models.h"
#include "rolling_model_base.h"
#include <algorithm>
#include <cmath>
#include "math_extra_liggghts.h"

namespace LIGGGHTS {
namespace ContactModels
{
  using namespace LAMMPS_NS;

  template<>
  class RollingModel<ROLLING_CDT> : public RollingModelBase {
  public:
    RollingModel(LAMMPS * lmp, IContactHistorySetup * hsetup, class ContactModelBase * c) :
        RollingModelBase(lmp, hsetup, c), coeffRollFrict(NULL)
    {
      
    }

    void registerSettings(Settings& settings)
    {
       settings.registerOnOff("torsionTorque", torsion_torque, false);
    }

    inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb)
    {}

    void connectToProperties(PropertyRegistry & registry)
    {
      registry.registerProperty("coeffRollFrict", &MODEL_PARAMS::createCoeffRollFrict);
      registry.connect("coeffRollFrict", coeffRollFrict,"rolling_model cdt");

      // error checks on coarsegraining
      if(force->cg_active())
        error->cg(FLERR,"rolling model cdt");
    }

    void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces) 
    {
      const double rmu = coeffRollFrict[sidata.itype][sidata.jtype];

      double r_torque[3], wr_roll[3];
      vectorZeroize3D(r_torque);

      const double radi = sidata.radi;
      const double radj = sidata.radj;
      double reff=sidata.is_wall ? radi : (radi*radj/(radi+radj));

#ifdef SUPERQUADRIC_ACTIVE_FLAG
      if(sidata.is_non_spherical && atom->superquadric_flag)
        reff = sidata.reff;
#endif

      if(sidata.is_wall){
        const double wr1 = sidata.wr1;
        const double wr2 = sidata.wr2;
        const double wr3 = sidata.wr3;
        const double wrmag = sqrt(wr1*wr1+wr2*wr2+wr3*wr3);
        if (wrmag > 0.)
        {
          const double kn = sidata.kn;
          const double enx = sidata.en[0];
          const double eny = sidata.en[1];
          const double enz = sidata.en[2];

          const double Fn = sidata.deltan*kn;
          r_torque[0] = rmu*Fn*wr1/wrmag*reff; 
          r_torque[1] = rmu*Fn*wr2/wrmag*reff;
          r_torque[2] = rmu*Fn*wr3/wrmag*reff;

          // remove normal (torsion) part of torque
          if(!torsion_torque)
          {
              double rtorque_dot_delta = r_torque[0]*enx+ r_torque[1]*eny + r_torque[2]*enz;
              double r_torque_n[3];
              r_torque_n[0] = enx * rtorque_dot_delta;
              r_torque_n[1] = eny * rtorque_dot_delta;
              r_torque_n[2] = enz * rtorque_dot_delta;
              vectorSubtract3D(r_torque,r_torque_n,r_torque);
          }
        }
      } else {

        vectorSubtract3D(atom->omega[sidata.i],atom->omega[sidata.j],wr_roll);
        const double wr_rollmag = vectorMag3D(wr_roll);

        if(wr_rollmag > 0.)
        {
          const double enx = sidata.en[0];
          const double eny = sidata.en[1];
          const double enz = sidata.en[2];

          // calculate torque
          vectorScalarMult3D(wr_roll,rmu*sidata.kn*sidata.deltan*reff/wr_rollmag,r_torque);

          // remove normal (torsion) part of torque
          if(!torsion_torque)
          {
              const double rtorque_dot_delta = r_torque[0]*enx + r_torque[1]*eny + r_torque[2]*enz;
              double r_torque_n[3];
              r_torque_n[0] = enx * rtorque_dot_delta;
              r_torque_n[1] = eny * rtorque_dot_delta;
              r_torque_n[2] = enz * rtorque_dot_delta;
              vectorSubtract3D(r_torque,r_torque_n,r_torque);
          }
        }
      }

      i_forces.delta_torque[0] -= r_torque[0];
      i_forces.delta_torque[1] -= r_torque[1];
      i_forces.delta_torque[2] -= r_torque[2];

      j_forces.delta_torque[0] += r_torque[0];
      j_forces.delta_torque[1] += r_torque[1];
      j_forces.delta_torque[2] += r_torque[2];
    }

    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void surfacesClose(SurfacesCloseData&, ForceData&, ForceData&){}

  private:
    double ** coeffRollFrict;
    bool torsion_torque;
  };
}
}
#endif // ROLLING_MODEL_CDT_H_
#endif
