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
NORMAL_MODEL(HOOKE_STIFFNESS,hooke/stiffness,1)
#else
#ifndef NORMAL_MODEL_HOOKE_STIFFNESS_H_
#define NORMAL_MODEL_HOOKE_STIFFNESS_H_
#include "contact_models.h"

namespace LIGGGHTS {
namespace ContactModels
{
  template<>
  class NormalModel<HOOKE_STIFFNESS> : protected Pointers
  {
  public:
    static const int MASK = CM_REGISTER_SETTINGS | CM_CONNECT_TO_PROPERTIES | CM_SURFACES_INTERSECT;

    NormalModel(LAMMPS * lmp, IContactHistorySetup*) : Pointers(lmp),
      k_n(NULL),
      k_t(NULL),
      gamma_n(NULL),
      gamma_t(NULL),
      tangential_damping(false),
      absolute_damping(false),
      limitForce(false),
      displayedSettings(false)
    {
      
    }

    void registerSettings(Settings & settings)
    {
      settings.registerOnOff("tangential_damping", tangential_damping, true);
      settings.registerOnOff("absolute_damping", absolute_damping);
      settings.registerOnOff("limitForce", limitForce);
    }

    void connectToProperties(PropertyRegistry & registry) {
      registry.registerProperty("k_n", &MODEL_PARAMS::createKn);
      registry.registerProperty("k_t", &MODEL_PARAMS::createKt);

      registry.connect("k_n", k_n,"model hooke/stiffness");
      registry.connect("k_t", k_t,"model hooke/stiffness");

      if(absolute_damping) {
        registry.registerProperty("gamman_abs", &MODEL_PARAMS::createGammanAbs);
        registry.registerProperty("gammat_abs", &MODEL_PARAMS::createGammatAbs);
        registry.connect("gamman_abs", gamma_n,"model hooke/stiffness");
        registry.connect("gammat_abs", gamma_t,"model hooke/stiffness");
      } else {
        registry.registerProperty("gamman", &MODEL_PARAMS::createGamman);
        registry.registerProperty("gammat", &MODEL_PARAMS::createGammat);
        registry.connect("gamman", gamma_n,"model hooke/stiffness");
        registry.connect("gammat", gamma_t,"model hooke/stiffness");
      }

      // error checks on coarsegraining
      if(force->cg_active())
        error->cg(FLERR,"model hooke/stiffness");
    }

    // effective exponent for stress-strain relationship
    
    inline double stressStrainExponent()
    {
      return 1.;
    }

    inline void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
      const int itype = sidata.itype;
      const int jtype = sidata.jtype;
      double meff=sidata.meff;

      double kn = k_n[itype][jtype];
      double kt = k_t[itype][jtype];
      double gamman, gammat;

      if(!displayedSettings)
      {
        displayedSettings = true;

        /*
        if(limitForce)
            if(0 == comm->me) fprintf(screen," NormalModel<HOOKE_STIFFNESS>: will limit normal force.\n");
        */
      }
      if(absolute_damping)
      {
        gamman = gamma_n[itype][jtype];
        gammat = gamma_t[itype][jtype];
      }
      else
      {
        gamman = meff*gamma_n[itype][jtype];
        gammat = meff*gamma_t[itype][jtype];
      }

      if (!tangential_damping) gammat = 0.0;

      // convert Kn and Kt from pressure units to force/distance^2
      kn /= force->nktv2p;
      kt /= force->nktv2p;

      const double Fn_damping = -gamman*sidata.vn;    
      const double Fn_contact = kn*(sidata.radsum-sidata.r);
      double Fn                       = Fn_damping + Fn_contact;

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
    bool absolute_damping;
    bool limitForce;
    bool displayedSettings;
  };
}
}
#endif // NORMAL_MODEL_HOOKE_STIFFNESS_H_
#endif
