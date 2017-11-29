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

   Contributing authors:
   Rahul Mohanty (University of Edinburgh, P&G)
   Tomaz M. Zorec (University of Ljubljana)
------------------------------------------------------------------------- */
#ifdef NORMAL_MODEL
NORMAL_MODEL(EDINBURGH,edinburgh,10)
#else
#ifndef NORMAL_MODEL_EDINBURGH_H_
#define NORMAL_MODEL_EDINBURGH_H_
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
  class NormalModel<EDINBURGH> : public NormalModelBase
  {
  public:
    NormalModel(LAMMPS * lmp, IContactHistorySetup * hsetup,class ContactModelBase *c) :
      NormalModelBase(lmp, hsetup, c),
      Yeff(NULL),
      Geff(NULL),
      CoeffRestLog(NULL),
      betaeff(NULL),
      kn2kc(NULL),
      kn2k1(NULL),
      cex(0.0),
      dex(0.0),
      f_adh(NULL),
      gamma_surf(NULL),
      history_offset(-1),
      kc_offset(-1),
      fo_offset(-1),
      tangential_damping(false),
      limitForce(false),
      fixKc(false)
    {
      history_offset = hsetup->add_history_value("deltaMax", "1");
      hsetup->add_history_value("old_delta", "1");
      kc_offset = hsetup->add_history_value("kc", "1");
      fo_offset = hsetup->add_history_value("fo", "1");
      c->add_history_offset("kc_offset", kc_offset);
      c->add_history_offset("fo_offset", fo_offset);
    }

    void registerSettings(Settings & settings)
    {
      settings.registerOnOff("tangential_damping", tangential_damping, true);
      settings.registerOnOff("limitForce", limitForce,true);
      settings.registerOnOff("fixKc", fixKc);
      //TODO error->one(FLERR,"TODO here also check if right surface model used");
    }

	inline void postSettings(IContactHistorySetup *, ContactModelBase *) {}

    void connectToProperties(PropertyRegistry & registry) {

      registry.registerProperty("Yeff", &MODEL_PARAMS::createYeff,"model edinburgh");                  // only used for non-linear model
      registry.registerProperty("Geff", &MODEL_PARAMS::createGeff,"model edinburgh");                  // only used for non-linear model
      registry.registerProperty("CoeffRestLog", &MODEL_PARAMS::createCoeffRestLog, "model edinburgh");
      registry.registerProperty("betaeff", &MODEL_PARAMS::createBetaEff,"model edinburgh");
      registry.registerProperty("kn2kc", &MODEL_PARAMS::createCoeffAdhesionStiffness,"model edinburgh");
      registry.registerProperty("kn2k1", &MODEL_PARAMS::createUnloadingStiffness, "model edinburgh");
      registry.registerProperty("cex", &MODEL_PARAMS::createAdhesionExponent, "model edinburgh");
      registry.registerProperty("dex", &MODEL_PARAMS::createOverlapExponent, "model edinburgh");
      registry.registerProperty("f_adh", &MODEL_PARAMS::createPullOffForce, "model edinburgh");
      registry.registerProperty("gamma_surf", &MODEL_PARAMS::createSurfaceEnergy, "model edinburgh");

      registry.connect("Yeff", Yeff,"model edinburgh");
      registry.connect("Geff", Geff,"model edinburgh");
      registry.connect("CoeffRestLog", CoeffRestLog, "model edinburgh");
      registry.connect("betaeff", betaeff,"model edinburgh");
      registry.connect("kn2kc", kn2kc,"model edinburgh");
      registry.connect("kn2k1", kn2k1,"model edinburgh");
      registry.connect("cex", cex,"model edinburgh");
      registry.connect("dex", dex,"model edinburgh");
      registry.connect("f_adh", f_adh,"model edinburgh");
      registry.connect("gamma_surf", gamma_surf,"model edinburgh");

    }

    // effective exponent for stress-strain relationship

    inline double stressStrainExponent()
    {
      return 1.5;
    }

    inline double calculate_deltan_p_max (double deltan_p, double * const history, int count_flag, const double f_0, double fTmp, const double k2, double dex, double dex_i, double deltan, double k_adh)
    {
      //calculating the maximum overlap
      double deltan_p_max;

      if (count_flag == 0 ) {
           if (deltan_p > history[0]) {
            deltan_p_max = deltan_p;
            history[0] = deltan_p_max;
           } else {
            deltan_p_max = history[0];
          }
      } else {
        deltan_p_max = pow((pow(history[1], dex) + ((k_adh/k2)*pow(history[1], cex))), dex_i);
        history[0] = deltan_p_max;
      }
      return history[0];
    }

    inline double whichd(double fTmp, double k1, double deltan_e, double d1, double d2)
    {
      //calculation distance between the particles based on branch
      if (fTmp >= k1 * deltan_e){
        return d1;
      } else{
        return d2;
      }
    }

    inline double calculate_k_adh(double d, double risq, double rjsq, const double g_surf, double f_min_lim, double deltan_pe_max, double deltan_p_max, const double k2, const double dex_i, const double cex)
    {
      double new_k_c, delta_min, a, f_min, dsq;
       dsq = d*d;
       a = (1./(2.*d)) * sqrt(4* dsq * risq - ((dsq - rjsq + risq)*(dsq - rjsq + risq)));
       f_min = ( 1.5 * M_PI * g_surf * a);

      if (f_min > f_min_lim)
      {
        f_min = f_min_lim;
        delta_min = 0.5 * deltan_p_max;
      } else {
        delta_min = pow((-f_min + k2 * deltan_pe_max)/k2, dex_i);
      }
      new_k_c = f_min/pow(delta_min, cex);
      return new_k_c;
    }

    inline void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
      const int itype = sidata.itype;
      const int jtype = sidata.jtype;
      const double deltan = sidata.deltan;
      double ri = sidata.radi;
      double rj = sidata.radj;
      double reff=sidata.is_wall ? sidata.radi : ((ri*rj)/(ri+rj));
#ifdef SUPERQUADRIC_ACTIVE_FLAG
        if(sidata.is_non_spherical && atom->superquadric_flag)
            reff = sidata.reff;
#endif
      double meff=sidata.meff;
      const double f_0 = f_adh[itype][jtype];
      double sqrtval = sqrt(reff*sidata.deltan);
      double Sn=2.*Yeff[itype][jtype]*sqrtval;
      double St=8.*Geff[itype][jtype]*sqrtval;
      double kn, kt;

      if(dex==1){
        kn = 2.0*Yeff[itype][jtype]*reff;
        kt = Sn;
      } else {
        kn = (4./3.)*Yeff[itype][jtype]*sqrt(reff);
        kt = St;
        }

      // convert Kn and Kt from pressure units to force/distance^2
      kn /= force->nktv2p;
      kt /= force->nktv2p;

      const double dex_i = 1./dex;
      const double k1 = kn;
      const double k2 = kn*kn2k1[itype][jtype];

      double gamman, gammat;
      gamman = sqrt(4.*meff*kn/(1.+(M_PI/CoeffRestLog[itype][jtype])*(M_PI/CoeffRestLog[itype][jtype])));
      gammat = sqrt(4.*meff*kn/(1.+(M_PI/CoeffRestLog[itype][jtype])*(M_PI/CoeffRestLog[itype][jtype])));

      if(sidata.contact_flags) *sidata.contact_flags |= CONTACT_NORMAL_MODEL;
      double * const history = &sidata.contact_history[history_offset];
      double * const kc_history = &sidata.contact_history[kc_offset];
      double * const fo_history = &sidata.contact_history[fo_offset];

      if (!tangential_damping) gammat = 0.0;

      double Fn_contact, deltan_e, deltan_ce, deltaMax_e;
      double k_adh = 0.0;

      if (fixKc == false) {
        const double lambda = pow((1. - k1/k2), dex_i);
        int count_flag = 0;

        // get the history value -- maximal overlap

        double deltan_p = lambda*deltan;
        double deltan_p_max, deltan_pe_max;

        deltan_e = pow(deltan, dex);
        deltan_ce = pow(deltan, cex);

        double fTmp = k2*(deltan_e - history[0]);
        double  d, d1, d2, r_sum;

        deltan_p_max = calculate_deltan_p_max (deltan_p, &sidata.contact_history[history_offset], count_flag, f_0, fTmp, k2, dex, dex_i, deltan, k_adh);

        // Normal force calculation for Edinburgh model

        const double g_surf = gamma_surf[itype][jtype];

        r_sum = sidata.radsum;
        d1 = sidata.is_wall ? ri : r_sum - deltan;
        d2 = sidata.is_wall ? ri  : r_sum - deltan_p_max;
        double risq = ri*ri;
        double rjsq = rj*rj;

        temp_calc:

        deltan_pe_max = pow(deltan_p_max, dex);
        const double f_min_lim = (k2 * deltan_pe_max) * 0.5;
        fTmp=k2*(deltan_e - deltan_pe_max);

        if (fTmp >= k1 * deltan_e){ // loading
          Fn_contact = k1 * deltan_e;
        }else{
          d = whichd(fTmp, k1, deltan_e, d1, d2);
          k_adh = calculate_k_adh(d, risq, rjsq, g_surf, f_min_lim, deltan_pe_max, deltan_p_max, k2, dex_i, cex);
          kc_history[0] = k_adh;
         if (fTmp > (-k_adh * deltan_ce)){
            Fn_contact = fTmp;
          }else{  // cohesion part

            if (deltan > history[1]){
              count_flag = 1;
              deltan_p_max = calculate_deltan_p_max (deltan_p, &sidata.contact_history[history_offset], count_flag, f_0, fTmp, k2, dex, dex_i, deltan, k_adh);
              goto temp_calc;
            }
            Fn_contact = -k_adh * deltan_ce;
          }
        }
        history[1] = deltan;
      }else{
        const double k_adh = kn2kc[itype][jtype] * kn ;
        kc_history[0] = k_adh;

        double deltaMax; // the 4th value of the history array is deltaMax
        if (deltan > history[0]) {
            history[0] = deltan;
            deltaMax = deltan;
        }else{
          deltaMax = history[0];
        }

        deltan_e = pow(deltan, dex);
        deltan_ce = pow(deltan, cex);
        deltaMax_e = pow(deltaMax, dex);

        const double fTmp = k2*(deltan_e-deltaMax_e)+k1*deltaMax_e;

        if (fTmp >= k1 * deltan_e) // loading
        {
          Fn_contact = k1 * deltan_e;
        } else {
         if (fTmp > (-k_adh * deltan_ce))
          {
            Fn_contact = fTmp;
          } else { // cohesion part
            Fn_contact = -k_adh * deltan_ce;
            const double newDeltaMax = ((k2 + k_adh)/(k2-k1))*deltan;
            history[0] = newDeltaMax;
          }
        }
      }

      const double Fn_damping = -gamman*sidata.vn;
      double Fn = Fn_damping + Fn_contact + f_0;
      fo_history[0]=f_0;

      if(limitForce && (Fn<0.0) && k_adh == 0 && f_0 == 0.0)
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

    void surfacesClose(SurfacesCloseData & scdata, ForceData&, ForceData&)
    {
      if(scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_NORMAL_MODEL;
      double * const history = &scdata.contact_history[history_offset];
      history[0] = 0.0;
      history[1] = 0.0;
    }
    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

  protected:
    double **Yeff;
    double **Geff;
    double **CoeffRestLog;
    double **betaeff;
    double **kn2kc;
    double **kn2k1;
    double cex;
    double dex;
    double **f_adh;
    double **gamma_surf;
    int history_offset;
    int kc_offset;
    int fo_offset;
    bool tangential_damping;
    bool limitForce;
    bool fixKc;
  };

}

}
#endif
#endif
