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

#ifdef NORMAL_MODEL
NORMAL_MODEL(HOOKE_HYSTERESIS,hooke/hysteresis,2)
#else
#ifndef NORMAL_MODEL_HOOKE_HYSTERESIS_H_
#define NORMAL_MODEL_HOOKE_HYSTERESIS_H_
#include "contact_models.h"
#include "normal_model_base.h"
#include <cmath>
#include "atom.h"
#include "force.h"
#include "update.h"
#include "global_properties.h"

namespace LIGGGHTS {
namespace ContactModels
{
  template<>
  class NormalModel<HOOKE_HYSTERESIS> : public NormalModel<HOOKE>
  {
  public:
    NormalModel(LAMMPS * lmp, IContactHistorySetup * hsetup,class ContactModelBase *c) :
        NormalModel<HOOKE>(lmp, hsetup,c),
        kn2k2Max(NULL),
        kn2kc(NULL),
        phiF(NULL)
    {
      history_offset = hsetup->add_history_value("deltaMax", "0");
      
    }

    inline void registerSettings(Settings & settings){
      NormalModel<HOOKE>::registerSettings(settings);
    }

    inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb) {}

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

    inline void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
      // use these values from HOOKE implementation
      bool & viscous = NormalModel<HOOKE>::viscous;
      double ** & Yeff = NormalModel<HOOKE>::Yeff;
      double & charVel = NormalModel<HOOKE>::charVel;
      bool & tangential_damping = NormalModel<HOOKE>::tangential_damping;
      Force * & force = NormalModel<HOOKE>::force;

      const int itype = sidata.itype;
      const int jtype = sidata.jtype;
      const double deltan = sidata.deltan;
      const double radi = sidata.radi;
      const double radj = sidata.radj;
      double reff=sidata.is_wall ? radi : (radi*radj/(radi+radj));
#ifdef SUPERQUADRIC_ACTIVE_FLAG
      if(sidata.is_non_spherical && atom->superquadric_flag)
        reff = sidata.reff;
#endif
      double meff=sidata.meff;
      double coeffRestLogChosen;

      if (viscous)  {
        double ** & coeffMu = NormalModel<HOOKE>::coeffMu;
        double ** & coeffRestMax = NormalModel<HOOKE>::coeffRestMax;
        double ** & coeffStc = NormalModel<HOOKE>::coeffStc;
        // Stokes Number from MW Schmeeckle (2001)
        const double stokes=sidata.meff*sidata.vn/(6.0*M_PI*coeffMu[itype][jtype]*reff*reff);
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
      if(sidata.contact_flags) *sidata.contact_flags |= CONTACT_NORMAL_MODEL;
      double * const history = &sidata.contact_history[history_offset];
      if (deltan > history[0]) {
          history[0] = deltan;
      }
      const double deltaMax = history[0]; // the 1st value of the history array is deltaMax

      // k2 dependent on the maximum overlap
      // this accounts for an increasing stiffness with deformation
      const double deltaMaxLim =(k2Max/(k2Max-kn))*phiF[itype][jtype]*2*reff;
      double k2, fHys;
      const bool update_history = sidata.computeflag && sidata.shearupdate;
      if (deltaMax >= deltaMaxLim) // big overlap ... no kn at all
      {
          k2 = k2Max;
          const double fTmp = k2*(deltan-deltaMaxLim)+kn*deltaMaxLim;//k2*(deltan-delta0);
          if (fTmp >= -kc*deltan) { // un-/reloading part (k2)
              fHys = fTmp;
          } else { // cohesion part
              fHys = -kc*deltan;

              const double newDeltaMax = 0.5*(deltan+sqrt(deltan*deltan+4*((kn+kc)*deltan*deltaMaxLim/(k2Max-kn))));
              if (update_history)
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
                  if (update_history)
                      history[0] = newDeltaMax;
              }
          }
      }

      const double Fn_damping = -gamman*sidata.vn;
      const double Fn = fHys + Fn_damping;

      sidata.Fn = Fn;
      sidata.kn = kn;
      sidata.kt = kt;
      sidata.gamman = gamman;
      sidata.gammat = gammat;

      #ifdef NONSPHERICAL_ACTIVE_FLAG
          double Fn_i[3] = { Fn * sidata.en[0], Fn * sidata.en[1], Fn * sidata.en[2]};
          double torque_i[3] = {0.0, 0.0, 0.0}; //initialized here with zeros to avoid compiler warnings
          if(sidata.is_non_spherical) {
            double xci[3];
            vectorSubtract3D(sidata.contact_point, atom->x[sidata.i], xci);
            vectorCross3D(xci, Fn_i, torque_i);
          }
      #endif
      // apply normal force
      if(sidata.is_wall) {
        const double Fn_ = Fn * sidata.area_ratio;
        i_forces.delta_F[0] += Fn_ * sidata.en[0];
        i_forces.delta_F[1] += Fn_ * sidata.en[1];
        i_forces.delta_F[2] += Fn_ * sidata.en[2];
        #ifdef NONSPHERICAL_ACTIVE_FLAG
                if(sidata.is_non_spherical) {
                  //for non-spherical particles normal force can produce torque!
                  i_forces.delta_torque[0] += torque_i[0];
                  i_forces.delta_torque[1] += torque_i[1];
                  i_forces.delta_torque[2] += torque_i[2];
                }
        #endif
      } else {
        i_forces.delta_F[0] += sidata.Fn * sidata.en[0];
        i_forces.delta_F[1] += sidata.Fn * sidata.en[1];
        i_forces.delta_F[2] += sidata.Fn * sidata.en[2];

        j_forces.delta_F[0] += -i_forces.delta_F[0];
        j_forces.delta_F[1] += -i_forces.delta_F[1];
        j_forces.delta_F[2] += -i_forces.delta_F[2];
        #ifdef NONSPHERICAL_ACTIVE_FLAG
                if(sidata.is_non_spherical) {
                  //for non-spherical particles normal force can produce torque!
                  double xcj[3], torque_j[3];
                  double Fn_j[3] = { -Fn_i[0], -Fn_i[1], -Fn_i[2]};
                  vectorSubtract3D(sidata.contact_point, atom->x[sidata.j], xcj);
                  vectorCross3D(xcj, Fn_j, torque_j);

                  i_forces.delta_torque[0] += torque_i[0];
                  i_forces.delta_torque[1] += torque_i[1];
                  i_forces.delta_torque[2] += torque_i[2];

                  j_forces.delta_torque[0] += torque_j[0];
                  j_forces.delta_torque[1] += torque_j[1];
                  j_forces.delta_torque[2] += torque_j[2];
                }
        #endif
      }
    }

    inline void surfacesClose(SurfacesCloseData & scdata, ForceData&, ForceData&)
    {
      if(scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_NORMAL_MODEL;
      double * const history = &scdata.contact_history[history_offset];
      history[0] = 0.0;
    }

    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

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
