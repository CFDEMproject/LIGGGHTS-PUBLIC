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

    Christoph Kloss (DCS Computing GmbH, Linz, JKU Linz)
    Richard Berger (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifdef NORMAL_MODEL
NORMAL_MODEL(HERTZ_STIFFNESS,hertz/stiffness,4)
#else
#ifndef NORMAL_MODEL_HERTZ_STIFFNESS_H_
#define NORMAL_MODEL_HERTZ_STIFFNESS_H_
#include "contact_models.h"
#include "global_properties.h"
#include "math.h"

namespace LIGGGHTS {
namespace ContactModels
{
  template<>
  class NormalModel<HERTZ_STIFFNESS> : protected Pointers
  {
  public:
    static const int MASK = CM_REGISTER_SETTINGS | CM_CONNECT_TO_PROPERTIES | CM_SURFACES_INTERSECT;

    NormalModel(LAMMPS * lmp, IContactHistorySetup*) : Pointers(lmp),
      k_n(NULL),
      k_t(NULL),
      gamma_n(NULL),
      gamma_t(NULL),
      tangential_damping(false),
      limitForce(false),
      displayedSettings(false)
    {
      
    }

    void registerSettings(Settings & settings)
    {
      settings.registerOnOff("tangential_damping", tangential_damping, true);
      settings.registerOnOff("limitForce", limitForce);
    }

    void connectToProperties(PropertyRegistry & registry)
    {
      registry.registerProperty("k_n", &MODEL_PARAMS::createKn);
      registry.registerProperty("k_t", &MODEL_PARAMS::createKt);
      registry.registerProperty("gamma_n", &MODEL_PARAMS::createGamman);
      registry.registerProperty("gamma_t", &MODEL_PARAMS::createGammat);

      registry.connect("k_n", k_n,"model hertz/stiffness");
      registry.connect("k_t", k_t,"model hertz/stiffness");
      registry.connect("gamma_n", gamma_n,"model hertz/stiffness");
      registry.connect("gamma_t", gamma_t,"model hertz/stiffness");

      // error checks on coarsegraining
      if(force->cg_active())
        error->cg(FLERR,"model hertz/stiffness");
    }

    // effective exponent for stress-strain relationship
    
    inline double stressStrainExponent()
    {
      return 1.5;
    }

    inline void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
      const int itype = sidata.itype;
      const int jtype = sidata.jtype;
      double meff = sidata.meff;
      double reff = sidata.is_wall ? sidata.radi : (sidata.radi*sidata.radj/(sidata.radi+sidata.radj));

      const double polyhertz = sqrt(reff*sidata.deltan);
      double kn = polyhertz*k_n[itype][jtype];
      double kt = polyhertz*k_t[itype][jtype];
      double gamman = polyhertz*meff*gamma_n[itype][jtype];
      double gammat = polyhertz*meff*gamma_t[itype][jtype];

      if(!tangential_damping) gammat = 0.0;

      if(!displayedSettings)
      {
        displayedSettings = true;

        /*
        if(limitForce)
            if(0 == comm->me) fprintf(screen," NormalModel<HERTZ_STIFFNESS>: will limit normal force.\n");
        */
      }
      // convert Kn and Kt from pressure units to force/distance^2
      kn /= force->nktv2p;
      kt /= force->nktv2p;

      const double Fn_damping = -gamman*sidata.vn;
      const double Fn_contact = kn*(sidata.radsum-sidata.r);
      double Fn               = Fn_damping + Fn_contact;

      //limit force to avoid the artefact of negative repulsion force
      if(limitForce && (Fn<0.0) )
      {
          Fn = 0.0;
      }

      sidata.Fn = Fn;
      sidata.kn = kn;
      sidata.kt = kt;
      sidata.gamman = gamman;
      sidata.gammat = gammat;

      // apply normal force
      if(sidata.is_wall) {
        const double Fn_ = Fn * sidata.area_ratio;
        i_forces.delta_F[0] = Fn_ * sidata.en[0];
        i_forces.delta_F[1] = Fn_ * sidata.en[1];
        i_forces.delta_F[2] = Fn_ * sidata.en[2];
      } else {
        i_forces.delta_F[0] = sidata.Fn * sidata.en[0];
        i_forces.delta_F[1] = sidata.Fn * sidata.en[1];
        i_forces.delta_F[2] = sidata.Fn * sidata.en[2];

        j_forces.delta_F[0] = -i_forces.delta_F[0];
        j_forces.delta_F[1] = -i_forces.delta_F[1];
        j_forces.delta_F[2] = -i_forces.delta_F[2];
      }
    }

    void surfacesClose(SurfacesCloseData&, ForceData&, ForceData&){}
    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

  protected:
    double ** k_n;
    double ** k_t;
    double ** gamma_n;
    double ** gamma_t;

    bool tangential_damping;
    bool limitForce;
    bool displayedSettings;
  };
}
}
#endif // NORMAL_MODEL_HERTZ_STIFFNESS_H_
#endif
