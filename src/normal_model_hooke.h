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
    static const int MASK = CM_REGISTER_SETTINGS | CM_CONNECT_TO_PROPERTIES | CM_COLLISION;

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

    inline void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces)
    {
      const int itype = cdata.itype;
      const int jtype = cdata.jtype;
      double ri = cdata.radi;
      double rj = cdata.radj;
      double reff=cdata.is_wall ? ri : (ri*rj/(ri+rj));
      double meff=cdata.meff;
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
         const double stokes=cdata.meff*cdata.vn/(6.0*M_PI*coeffMu[itype][jtype]*reff*reff);
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

      const double Fn_damping = -gamman*cdata.vn;
      const double Fn_contact = kn*(cdata.radsum-cdata.r);
      double Fn                         = Fn_damping + Fn_contact;

      //limit force to avoid the artefact of negative repulsion force
      if(limitForce && (Fn<0.0) )
      {
          Fn = 0.0;
      }
      cdata.Fn = Fn;

      cdata.kn = kn;
      cdata.kt = kt;
      cdata.gamman = gamman;
      cdata.gammat = gammat;

      // apply normal force
      if(cdata.is_wall) {
        const double Fn_ = Fn * cdata.area_ratio;
        i_forces.delta_F[0] = Fn_ * cdata.en[0];
        i_forces.delta_F[1] = Fn_ * cdata.en[1];
        i_forces.delta_F[2] = Fn_ * cdata.en[2];
      } else {
        const double fx = cdata.Fn * cdata.en[0];
        const double fy = cdata.Fn * cdata.en[1];
        const double fz = cdata.Fn * cdata.en[2];

        i_forces.delta_F[0] = fx;
        i_forces.delta_F[1] = fy;
        i_forces.delta_F[2] = fz;

        j_forces.delta_F[0] = -fx;
        j_forces.delta_F[1] = -fy;
        j_forces.delta_F[2] = -fz;
      }
    }

    inline void noCollision(ContactData&, ForceData&, ForceData&){}
    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}

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
