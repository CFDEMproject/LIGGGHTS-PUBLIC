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

    Rahul Mohanty (University of Edinburgh, P&G)
------------------------------------------------------------------------- */

#ifdef NORMAL_MODEL
NORMAL_MODEL(LUDING,luding,12)
#else
#ifndef NORMAL_MODEL_LUDING_H_
#define NORMAL_MODEL_LUDING_H_
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
  class NormalModel<LUDING> : public NormalModelBase
  {
  public:
    NormalModel(LAMMPS * lmp, IContactHistorySetup * hsetup,class ContactModelBase *c) :
      NormalModelBase(lmp, hsetup, c),
      K_elastic(NULL),
      CoeffRestLog(NULL),
      kn2k1(NULL),
      kn2kc(NULL),
      phiF(NULL),
      f_adh(NULL),
      limitForce(false)
    {
      history_offset = hsetup->add_history_value("deltaMax", "0");
      kc_offset = hsetup->add_history_value("kc", "1");
      fo_offset = hsetup->add_history_value("fo", "1");
      c->add_history_offset("kc_offset", kc_offset);
      c->add_history_offset("fo_offset", fo_offset);
    }

    inline void registerSettings(Settings & settings){
      settings.registerOnOff("tangential_damping", tangential_damping, true);
      settings.registerOnOff("limitForce", limitForce, true);
    }
    inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb) {}

    inline void connectToProperties(PropertyRegistry & registry) {
      registry.registerProperty("K_elastic", &MODEL_PARAMS::createLoadingStiffness,"model luding");
      registry.registerProperty("CoeffRestLog", &MODEL_PARAMS::createCoeffRestLog,"model luding");
      registry.registerProperty("kn2k1", &MODEL_PARAMS::createUnloadingStiffness,"model luding");
      registry.registerProperty("kn2kc", &MODEL_PARAMS::createCoeffAdhesionStiffness,"model luding");
      registry.registerProperty("phiF", &MODEL_PARAMS::createCoeffPlasticityDepth,"model luding");
      registry.registerProperty("f_adh", &MODEL_PARAMS::createPullOffForce,"model luding");

      registry.connect("K_elastic", K_elastic,"model luding");
      registry.connect("CoeffRestLog", CoeffRestLog,"model luding");
      registry.connect("kn2kc", kn2kc,"model luding");
      registry.connect("kn2k1", kn2k1, "model luding");
      registry.connect("phiF", phiF,"model luding");
      registry.connect("f_adh", f_adh,"model luding");

      // error checks on coarsegraining
      if(force->cg_active())
        error->cg(FLERR,"model luding");
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
      const double deltan = sidata.deltan;
      double ri = sidata.radi;
      double rj = sidata.radj;
      double reff=sidata.is_wall ? sidata.radi : (ri*rj/(ri+rj));
#ifdef SUPERQUADRIC_ACTIVE_FLAG
      if(sidata.is_non_spherical && atom->superquadric_flag)
        reff = sidata.reff;
#endif
      double meff=sidata.meff;
      double kn = K_elastic[itype][jtype];
      double kt = kn;

      // convert Kn and Kt from pressure units to force/distance^2
      kn /= force->nktv2p;
      kt /= force->nktv2p;

      const double k1 = kn;
      const double k2Max = kn * kn2k1[itype][jtype];

      const double kc = kn2kc[itype][jtype] * kn;
      const double f_0 = f_adh[itype][jtype];

      double gamman, gammat;

      gamman = sqrt(4.*meff*kn/(1.+(M_PI/CoeffRestLog[itype][jtype])*(M_PI/CoeffRestLog[itype][jtype])));
      gammat = sqrt(4.*meff*kn/(1.+(M_PI/CoeffRestLog[itype][jtype])*(M_PI/CoeffRestLog[itype][jtype])));

      if (!tangential_damping) gammat = 0.0;

      // get the history value -- maximal overlap
      if(sidata.contact_flags) *sidata.contact_flags |= CONTACT_NORMAL_MODEL;
      double * const history = &sidata.contact_history[history_offset];
      double * const kc_history = &sidata.contact_history[kc_offset];
      double * const fo_history = &sidata.contact_history[fo_offset];
      double deltaMax; // the 4th value of the history array is deltaMax
      if (deltan > history[0]) {
          history[0] = deltan;
          deltaMax = deltan;
      } else{
        deltaMax = history[0];
      }

      // k2 dependent on the maximum overlap
      // this accounts for an increasing stiffness with deformation - to capture nonlinearity
      const double deltaMaxLim =(k2Max/(k2Max-k1))*phiF[itype][jtype]*2*reff;

      double k2, fHys;

      if (deltaMax >= deltaMaxLim) // big overlap ... no kn at all
      {
          k2 = k2Max;
          const double fTmp = k2*(deltan-deltaMaxLim)+k1*deltaMaxLim;//k2*(deltan-delta0);
          if (fTmp >= -kc*deltan) { // un-/reloading part (k2)
              fHys = fTmp;
          } else { // cohesion part
              fHys = -kc*deltan;
              const double newDeltaMax = ((k2 + kc)/(k2-k1))*deltan;
              history[0] = newDeltaMax;
          }
      } else {
          k2 = k1 + (k2Max - k1) * deltaMax/deltaMaxLim;
          const double fTmp = k2*(deltan-deltaMax)+k1*deltaMax;//k2*(deltan-delta0);
          if (fTmp >= k1*deltan) { // loading part (k1)
              fHys = k1*deltan;
          } else { // un-/reloading part (k2)
              if (fTmp > -kc*deltan) {
                  fHys = fTmp;
              } else { // cohesion part
                  fHys = -kc*deltan;
                  const double newDeltaMax = ((k2 + kc)/(k2-k1))*deltan;
                  history[0] = newDeltaMax;
              }
          }
      }

      const double Fn_damping = -gamman*sidata.vn;
      double Fn = fHys + Fn_damping + f_0;

      if(limitForce && (Fn<0.0) && kc == 0 && f_0 == 0.0){
          Fn = 0.0;
      }
      sidata.Fn = Fn;
      sidata.kn = kn;
      sidata.kt = kt;
      kc_history[0] = kc;
      fo_history[0] = f_0;
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
        i_forces.delta_F[0] = Fn_ * sidata.en[0];
        i_forces.delta_F[1] = Fn_ * sidata.en[1];
        i_forces.delta_F[2] = Fn_ * sidata.en[2];
        #ifdef NONSPHERICAL_ACTIVE_FLAG
        if(sidata.is_non_spherical) {
          //for non-spherical particles normal force can produce torque!
          i_forces.delta_torque[0] += torque_i[0];
          i_forces.delta_torque[1] += torque_i[1];
          i_forces.delta_torque[2] += torque_i[2];
        }
        #endif
      } else {
        i_forces.delta_F[0] = sidata.Fn * sidata.en[0];
        i_forces.delta_F[1] = sidata.Fn * sidata.en[1];
        i_forces.delta_F[2] = sidata.Fn * sidata.en[2];

        j_forces.delta_F[0] = -i_forces.delta_F[0];
        j_forces.delta_F[1] = -i_forces.delta_F[1];
        j_forces.delta_F[2] = -i_forces.delta_F[2];
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
    double **K_elastic;
    double **CoeffRestLog;
    double **kn2k1;
    double **kn2kc;
    double **phiF;
    double **f_adh;

    int history_offset;
    int kc_offset;
    int fo_offset;

    bool tangential_damping;
    bool limitForce;
  };
}
}
#endif // NORMAL_MODEL_LUDING_H_
#endif
