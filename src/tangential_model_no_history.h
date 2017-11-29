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
#include "tangential_model_base.h"
#include <algorithm>
#include <cmath>
#include "global_properties.h"
#include "fix_property_atom.h"

namespace LIGGGHTS {
namespace ContactModels
{
  template<>
  class TangentialModel<TANGENTIAL_NO_HISTORY> : public TangentialModelBase
  {
  public:
    static const int MASK = CM_CONNECT_TO_PROPERTIES | CM_SURFACES_INTERSECT;

    TangentialModel(LAMMPS * lmp, IContactHistorySetup *hsetup, class ContactModelBase *cmb) :
      TangentialModelBase(lmp, hsetup, cmb),
      coeffFrict(NULL),
      disable_when_bonded_(false),
      bond_history_offset_(-1),
      dissipation_history_offset_(-1),
      dissipatedflag_(false),
      fix_dissipated_(NULL)
    {
      
    }

    inline void registerSettings(Settings &settings)
    {
      settings.registerOnOff("disableTangentialWhenBonded", disable_when_bonded_, false);
      settings.registerOnOff("computeDissipatedEnergy", dissipatedflag_, false);
    }

    inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb)
    {
      if (disable_when_bonded_)
      {
        bond_history_offset_ = cmb->get_history_offset("bond_contactflag");
        if (bond_history_offset_ < 0)
          error->one(FLERR, "Could not find bond history offset");
      }
      if (dissipatedflag_)
      {
        if (cmb->is_wall())
        {
            fix_dissipated_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("dissipated_energy_wall", "property/atom", "vector", 0, 0, "dissipated energy"));
            dissipation_history_offset_ = cmb->get_history_offset("dissipation_force");
            if (!dissipation_history_offset_)
                error->one(FLERR, "Internal error: Could not find dissipation history offset");
        }
        else
            fix_dissipated_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("dissipated_energy", "property/atom", "vector", 0, 0, "dissipated energy"));
        if (!fix_dissipated_)
            error->one(FLERR, "Surface model has not registered dissipated_energy fix");
      }
    }

    inline void connectToProperties(PropertyRegistry & registry){
      registry.registerProperty("coeffFrict", &MODEL_PARAMS::createCoeffFrict);
      registry.connect("coeffFrict", coeffFrict,"tangential_model history");
    }

    inline void surfacesIntersect(const SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces) {
      if(sidata.contact_flags) *sidata.contact_flags |= CONTACT_TANGENTIAL_MODEL;
      const double xmu = coeffFrict[sidata.itype][sidata.jtype];
      const double enx = sidata.en[0];
      const double eny = sidata.en[1];
      const double enz = sidata.en[2];
      const double vrel = sqrt(sidata.vtr1*sidata.vtr1 + sidata.vtr2*sidata.vtr2 + sidata.vtr3*sidata.vtr3);

      // force normalization
      const double Ft_friction = xmu * fabs(sidata.Fn);
      double gamma = 0.0;

      if (Ft_friction < sidata.gammat*vrel)
        gamma = Ft_friction/vrel;
      else
        gamma = sidata.gammat;

      // tangential force due to tangential velocity damping

      const double Ft1 = -gamma*sidata.vtr1;
      const double Ft2 = -gamma*sidata.vtr2;
      const double Ft3 = -gamma*sidata.vtr3;

      // forces & torques

      const double tor1 = (eny*Ft3 - enz*Ft2);
      const double tor2 = (enz*Ft1 - enx*Ft3);
      const double tor3 = (enx*Ft2 - eny*Ft1);

      #ifdef NONSPHERICAL_ACTIVE_FLAG
          double torque_i[3];
          if(sidata.is_non_spherical) {
            double xci[3];
            double Ft_i[3] = { Ft1,  Ft2,  Ft3 };
            vectorSubtract3D(sidata.contact_point, atom->x[sidata.i], xci);
            vectorCross3D(xci, Ft_i, torque_i);
          } else {
            torque_i[0] = -sidata.cri * tor1;
            torque_i[1] = -sidata.cri * tor2;
            torque_i[2] = -sidata.cri * tor3;
          }
      #endif
      // return resulting forces
      if (!disable_when_bonded_ || sidata.contact_history[bond_history_offset_] < 0.5)
      {
        // energy balance terms
        if (dissipatedflag_ && sidata.computeflag && sidata.shearupdate)
        {
            const double crj = sidata.is_wall ? 0.0 : sidata.crj;
            // compute increment in dissipated energy
            double * const * const dissipated = fix_dissipated_->array_atom;
            double * const dissipated_i = dissipated[sidata.i];
            double * const dissipated_j = dissipated[sidata.j];
            dissipated_i[1] += -Ft1;
            dissipated_i[2] += -Ft2;
            dissipated_i[3] += -Ft3;
            dissipated_i[4] += sidata.cri*tor1;
            dissipated_i[5] += sidata.cri*tor2;
            dissipated_i[6] += sidata.cri*tor3;
            if (sidata.j < atom->nlocal && !sidata.is_wall)
            {
                dissipated_j[1] -= -Ft1;
                dissipated_j[2] -= -Ft2;
                dissipated_j[3] -= -Ft3;
                dissipated_j[4] += crj*tor1;
                dissipated_j[5] += crj*tor2;
                dissipated_j[6] += crj*tor3;
            }
            else if (sidata.is_wall)
            {
                double * const diss_force = &sidata.contact_history[dissipation_history_offset_];
                diss_force[0] -= -Ft1;
                diss_force[1] -= -Ft2;
                diss_force[2] -= -Ft3;
            }
        }
        if(sidata.is_wall) {
          const double area_ratio = sidata.area_ratio;
          i_forces.delta_F[0] += Ft1 * area_ratio;
          i_forces.delta_F[1] += Ft2 * area_ratio;
          i_forces.delta_F[2] += Ft3 * area_ratio;
          #ifdef NONSPHERICAL_ACTIVE_FLAG
                  i_forces.delta_torque[0] += torque_i[0] * area_ratio;
                  i_forces.delta_torque[1] += torque_i[1] * area_ratio;
                  i_forces.delta_torque[2] += torque_i[2] * area_ratio;
          #else
                  i_forces.delta_torque[0] += -sidata.cri * tor1 * area_ratio;
                  i_forces.delta_torque[1] += -sidata.cri * tor2 * area_ratio;
                  i_forces.delta_torque[2] += -sidata.cri * tor3 * area_ratio;
          #endif
        } else {
          i_forces.delta_F[0] += Ft1;
          i_forces.delta_F[1] += Ft2;
          i_forces.delta_F[2] += Ft3;
          j_forces.delta_F[0] -= Ft1;
          j_forces.delta_F[1] -= Ft2;
          j_forces.delta_F[2] -= Ft3;
          #ifdef NONSPHERICAL_ACTIVE_FLAG
                  double torque_j[3];
                  if(sidata.is_non_spherical) {
                    double xcj[3];
                    vectorSubtract3D(sidata.contact_point, atom->x[sidata.j], xcj);
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
                  i_forces.delta_torque[0] += -sidata.cri * tor1;
                  i_forces.delta_torque[1] += -sidata.cri * tor2;
                  i_forces.delta_torque[2] += -sidata.cri * tor3;

                  j_forces.delta_torque[0] += -sidata.crj * tor1;
                  j_forces.delta_torque[1] += -sidata.crj * tor2;
                  j_forces.delta_torque[2] += -sidata.crj * tor3;
          #endif
        }
      }
    }

    inline void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    inline void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

    inline void surfacesClose(SurfacesCloseData &scdata, ForceData&, ForceData&)
    {
        if(scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_TANGENTIAL_MODEL;
    }

  private:
    double ** coeffFrict;
    bool disable_when_bonded_;
    int bond_history_offset_;
    int dissipation_history_offset_;
    bool dissipatedflag_;
    FixPropertyAtom *fix_dissipated_;
  };
}
}
#endif // TANGENTIAL_MODEL_NO_HISTORY_H_
#endif
