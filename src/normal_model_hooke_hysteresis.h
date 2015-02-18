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
NORMAL_MODEL(HOOKE_HYSTERESIS,hooke/hysteresis,2)
#else
#ifndef NORMAL_MODEL_HOOKE_HYSTERESIS_H_
#define NORMAL_MODEL_HOOKE_HYSTERESIS_H_
#include "contact_models.h"
#include "math.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "global_properties.h"

namespace LIGGGHTS {
namespace ContactModels
{
  template<>
  class NormalModel<HOOKE_HYSTERESIS> : protected NormalModel<HOOKE>
  {
  public:
    static const int MASK = CM_REGISTER_SETTINGS | CM_CONNECT_TO_PROPERTIES | CM_COLLISION | CM_NO_COLLISION;

    NormalModel(LAMMPS * lmp, IContactHistorySetup * hsetup) : NormalModel<HOOKE>(lmp, hsetup),
        kn2k2Max(NULL),
        kn2kc(NULL),
        phiF(NULL)
    {
      history_offset = hsetup->add_history_value("deltaMax", "1");
      
    }

    inline void registerSettings(Settings & settings){
      NormalModel<HOOKE>::registerSettings(settings);
    }

    inline void connectToProperties(PropertyRegistry & registry) {
      NormalModel<HOOKE>::connectToProperties(registry);

      registry.registerProperty("kn2kcMax", &MODEL_PARAMS::createCoeffMaxElasticStiffness);
      registry.registerProperty("kn2kc", &MODEL_PARAMS::createCoeffAdhesionStiffness);
      registry.registerProperty("phiF", &MODEL_PARAMS::createCoeffPlasticityDepth);

      registry.connect("kn2kcMax", kn2k2Max,"model hooke/hysteresis");
      registry.connect("kn2kc", kn2kc,"model hooke/hysteresis");
      registry.connect("phiF", phiF,"model hooke/hysteresis");

      // error checks on coarsegraining
      if(force->cg_active())
        error->cg(FLERR,"model hooke/hysteresis");
    }

    // effective exponent for stress-strain relationship
    
    inline double stressStrainExponent()
    {
      return 1.;
    }

    inline void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces)
    {
      // use these values from HOOKE implementation
      bool & viscous = NormalModel<HOOKE>::viscous;
      double ** & Yeff = NormalModel<HOOKE>::Yeff;
      double & charVel = NormalModel<HOOKE>::charVel;
      bool & tangential_damping = NormalModel<HOOKE>::tangential_damping;
      Force * & force = NormalModel<HOOKE>::force;

      const int itype = cdata.itype;
      const int jtype = cdata.jtype;
      const double deltan = cdata.deltan;
      double ri = cdata.radi;
      double rj = cdata.radj;
      double reff=cdata.is_wall ? cdata.radi : (ri*rj/(ri+rj));
      double meff=cdata.meff;
      double coeffRestLogChosen;

      if (viscous)  {
        double ** & coeffMu = NormalModel<HOOKE>::coeffMu;
        double ** & coeffRestMax = NormalModel<HOOKE>::coeffRestMax;
        double ** & coeffStc = NormalModel<HOOKE>::coeffStc;
        // Stokes Number from MW Schmeeckle (2001)
        const double stokes=cdata.meff*cdata.vn/(6.0*M_PI*coeffMu[itype][jtype]*reff*reff);
        // Empirical from Legendre (2006)
        coeffRestLogChosen=log(coeffRestMax[itype][jtype])+coeffStc[itype][jtype]/stokes;
      } else {
        double ** & coeffRestLog = NormalModel<HOOKE>::coeffRestLog;
        coeffRestLogChosen=coeffRestLog[itype][jtype];
      }

      const double sqrtval = sqrt(reff);
      double kn = 16./15.*sqrtval*(Yeff[itype][jtype])*pow(15.*meff*charVel*charVel/(16.*sqrtval*Yeff[itype][jtype]),0.2);
      double kt = kn;
      const double gamman = sqrt(4.*meff*kn/(1.+(M_PI/coeffRestLogChosen)*(M_PI/coeffRestLogChosen)));
      const double gammat = tangential_damping ? gamman : 0.0;

      // convert Kn and Kt from pressure units to force/distance^2
      kn /= force->nktv2p;
      kt /= force->nktv2p;

      // coefficients
      
      const double k2Max = kn * kn2k2Max[itype][jtype]; 
      const double kc = kn * kn2kc[itype][jtype]; 

      // get the history value -- maximal overlap
      if(cdata.touch) *cdata.touch |= TOUCH_NORMAL_MODEL;
      double * const history = &cdata.contact_history[history_offset];
      double deltaMax; // the 4th value of the history array is deltaMax
      if (deltan > history[0]) {
          history[0] = deltan;
          deltaMax = deltan;
      } else
          deltaMax = history[0];

      // k2 dependent on the maximum overlap
      // this accounts for an increasing stiffness with deformation
      const double deltaMaxLim =(k2Max/(k2Max-kn))*phiF[itype][jtype]*2*reff;
      double k2, fHys;
      if (deltaMax >= deltaMaxLim) // big overlap ... no kn at all
      {
          k2 = k2Max;
          const double fTmp = k2*(deltan-deltaMaxLim)+kn*deltaMaxLim;//k2*(deltan-delta0);
          if (fTmp >= -kc*deltan) { // un-/reloading part (k2)
              fHys = fTmp;
          } else { // cohesion part
              fHys = -kc*deltan;

              const double newDeltaMax = 0.5*(deltan+sqrt(deltan*deltan+4*((kn+kc)*deltan*deltaMaxLim/(k2Max-kn))));
              history[0] = newDeltaMax;
          }
      } else {
          k2 = kn+(k2Max-kn)*deltaMax/deltaMaxLim;
          const double fTmp = k2*(deltan-deltaMax)+kn*deltaMax;//k2*(deltan-delta0);
          if (fTmp >= kn*deltan) { // loading part (kn)
              fHys = kn*deltan;
          } else { // un-/reloading part (k2)
              if (fTmp > -kc*deltan) {
                  fHys = fTmp;
              } else { // cohesion part
                  fHys = -kc*deltan;

                  const double newDeltaMax = 0.5*(deltan+sqrt(deltan*deltan+4*((kn+kc)*deltan*deltaMaxLim/(k2Max-kn))));
                  history[0] = newDeltaMax;
              }
          }
      }

      const double Fn_damping = -gamman*cdata.vn;
      const double Fn = fHys + Fn_damping;

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
        i_forces.delta_F[0] = cdata.Fn * cdata.en[0];
        i_forces.delta_F[1] = cdata.Fn * cdata.en[1];
        i_forces.delta_F[2] = cdata.Fn * cdata.en[2];

        j_forces.delta_F[0] = -i_forces.delta_F[0];
        j_forces.delta_F[1] = -i_forces.delta_F[1];
        j_forces.delta_F[2] = -i_forces.delta_F[2];
      }
    }

    inline void noCollision(ContactData & cdata, ForceData&, ForceData&)
    {
      if(cdata.touch) *cdata.touch &= ~TOUCH_NORMAL_MODEL;
      double * const history = &cdata.contact_history[history_offset];
      history[0] = 0.0;
    }

    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}

  protected:
    double **kn2k2Max;
    double **kn2kc;
    double **phiF;
    int history_offset;
  };
}
}
#endif // NORMAL_MODEL_HOOKE_HYSTERESIS_H_
#endif
