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

NORMAL_MODEL(HOOKE,hooke,0)

#else

#ifndef NORMAL_MODEL_HOOKE_H_
#define NORMAL_MODEL_HOOKE_H_

#include "contact_models.h"
#include "math.h"
#include "atom.h"
#include "force.h"
#include "update.h"

namespace LIGGGHTS {
namespace ContactModels
{
  template<>
  class NormalModel<HOOKE> : protected Pointers
  {
  public:
    static const int MASK = CM_REGISTER_SETTINGS | CM_CONNECT_TO_PROPERTIES | CM_SURFACES_INTERSECT;

    NormalModel(LAMMPS * lmp, IContactHistorySetup*) : Pointers(lmp),
      Yeff(NULL),
      Geff(NULL),
      coeffRestMax(NULL),
      coeffRestLog(NULL),
      coeffMu(NULL),
      coeffStc(NULL),
      charVel(0.0),
      viscous(false),
      tangential_damping(false),
      limitForce(false),
      ktToKn(false),
      displayedSettings(false)
    {
      
    }

    inline void registerSettings(Settings & settings)
    {
      settings.registerOnOff("viscous", viscous);
      settings.registerOnOff("tangential_damping", tangential_damping, true);
      settings.registerOnOff("limitForce", limitForce);
      settings.registerOnOff("ktToKnUser", ktToKn);
    }

    inline void connectToProperties(PropertyRegistry & registry) {
      registry.registerProperty("Yeff", &MODEL_PARAMS::createYeff);
      registry.registerProperty("Geff", &MODEL_PARAMS::createGeff);
      registry.registerProperty("charVel", &MODEL_PARAMS::createCharacteristicVelocity);

      registry.connect("Yeff", Yeff,"model hooke");
      registry.connect("Geff", Geff,"model hooke");
      registry.connect("charVel", charVel,"model hooke");

      if(viscous) {
        registry.registerProperty("coeffMu", &MODEL_PARAMS::createCoeffMu);
        registry.registerProperty("coeffStc", &MODEL_PARAMS::createCoeffStc);
        registry.registerProperty("coeffRestMax", &MODEL_PARAMS::createCoeffRestMax);

        registry.connect("coeffMu", coeffMu,"model hooke viscous");
        registry.connect("coeffStc", coeffStc,"model hooke viscous");
        registry.connect("coeffRestMax", coeffRestMax,"model hooke viscous");
        //registry.connect("log(coeffRestMax)+coeffStc", logRestMaxPlusStc);
      } else {
        registry.registerProperty("coeffRestLog", &MODEL_PARAMS::createCoeffRestLog);

        registry.connect("coeffRestLog", coeffRestLog,"model hooke viscous");
      }

      // error checks on coarsegraining
      if(force->cg_active())
        error->cg(FLERR,"model hooke");
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
      double ri = sidata.radi;
      double rj = sidata.radj;
      double reff=sidata.is_wall ? ri : (ri*rj/(ri+rj));
      double meff=sidata.meff;
      double coeffRestLogChosen;

      const double sqrtval = sqrt(reff);

      if(!displayedSettings)
      {
        displayedSettings = true;
        /*
        if(ktToKn)
            if(0 == comm->me) fprintf(screen," NormalModel<HOOKE>: will use user-modified ktToKn of 2/7.\n");
        if(tangential_damping)
            if(0 == comm->me) fprintf(screen," NormalModel<HOOKE>: will apply tangential damping.\n");
        if(viscous)
            if(0 == comm->me) fprintf(screen," NormalModel<HOOKE>: will apply damping based on Stokes number.\n");
        if(limitForce)
            if(0 == comm->me) fprintf(screen," NormalModel<HOOKE>: will limit normal force.\n");
        */
      }
      if (viscous)  {
         // Stokes Number from MW Schmeeckle (2001)
         const double stokes=sidata.meff*sidata.vn/(6.0*M_PI*coeffMu[itype][jtype]*reff*reff);
         // Empirical from Legendre (2006)
         coeffRestLogChosen=log(coeffRestMax[itype][jtype])+coeffStc[itype][jtype]/stokes;
      } else {
         coeffRestLogChosen=coeffRestLog[itype][jtype];
      }

      double kn = 16./15.*sqrtval*(Yeff[itype][jtype])*pow(15.*meff*charVel*charVel/(16.*sqrtval*Yeff[itype][jtype]),0.2);
      double kt = kn;
      if(ktToKn) kt *= 0.285714286; //2//7
      const double gamman=sqrt(4.*meff*kn/(1.+(M_PI/coeffRestLogChosen)*(M_PI/coeffRestLogChosen)));
      const double gammat = tangential_damping ? gamman : 0.0;

      // convert Kn and Kt from pressure units to force/distance^2
      kn /= force->nktv2p;
      kt /= force->nktv2p;

      const double Fn_damping = -gamman*sidata.vn;
      const double Fn_contact = kn*(sidata.radsum-sidata.r);
      double Fn                         = Fn_damping + Fn_contact;

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
        const double fx = sidata.Fn * sidata.en[0];
        const double fy = sidata.Fn * sidata.en[1];
        const double fz = sidata.Fn * sidata.en[2];

        i_forces.delta_F[0] = fx;
        i_forces.delta_F[1] = fy;
        i_forces.delta_F[2] = fz;

        j_forces.delta_F[0] = -fx;
        j_forces.delta_F[1] = -fy;
        j_forces.delta_F[2] = -fz;
      }
    }

    inline void surfacesClose(SurfacesCloseData&, ForceData&, ForceData&){}
    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

  protected:
    double ** Yeff;
    double ** Geff;
    double ** coeffRestMax;
    double ** coeffRestLog;
    double ** coeffMu;
    double ** coeffStc;
    double charVel;

    bool viscous;
    bool tangential_damping;
    bool limitForce;
    bool ktToKn;
    bool displayedSettings;
  };
}
}
#endif // NORMAL_MODEL_HOOKE_H_
#endif
