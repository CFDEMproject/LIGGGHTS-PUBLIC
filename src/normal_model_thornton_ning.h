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
    Tomaz M. Zorec (University of Ljubljana)
------------------------------------------------------------------------- */
#ifdef NORMAL_MODEL
NORMAL_MODEL(THORNTON_NING,thornton_ning,8)
#else
#ifndef NORMAL_THORNTON_NING_H_
#define NORMAL_THORNTON_NING_H_
#include "contact_models.h"
#include "normal_model_base.h"
#include <cmath>
#include "atom.h"
#include "force.h"
#include "update.h"
#include "global_properties.h"
#include "math_extra_liggghts.h"
namespace LIGGGHTS {

namespace ContactModels
{
  template<>
  class NormalModel<THORNTON_NING> : public NormalModelBase
  {
  public:

    NormalModel(LAMMPS * lmp, IContactHistorySetup * hsetup,class ContactModelBase *c) :
      NormalModelBase(lmp, hsetup, c),
      Yeff(NULL),
      Geff(NULL),
      betaeff(NULL),
      limitForce(false),
      displayedSettings(false)
    {
      history_offset = hsetup->add_history_value("tn_virgin_flag", "1");
      hsetup->add_history_value("delta_old", "1");
      hsetup->add_history_value("delta_max", "1");
      hsetup->add_history_value("force_old", "1");
      hsetup->add_history_value("force_max", "1");
      hsetup->add_history_value("adhesion_flag", "1");
      hsetup->add_history_value("detaching_delta", "1");
      hsetup->add_history_value("detaching_flag", "1");
      hsetup->add_history_value("detaching_force", "1");
      hsetup->add_history_value("yielding_flag", "1");
      kc_offset = hsetup->add_history_value("kc", "1");
      fo_offset = hsetup->add_history_value("fo", "1");
      c->add_history_offset("kc_offset", kc_offset);
      c->add_history_offset("fo_offset", fo_offset);

    }

    void registerSettings(Settings & settings)
    {
      settings.registerOnOff("tangential_damping", tangential_damping, true);
      settings.registerOnOff("limitForce", limitForce);
    }

	inline void postSettings(IContactHistorySetup *, ContactModelBase *) {}

    void connectToProperties(PropertyRegistry & registry) {
      registry.registerProperty("Yeff", &MODEL_PARAMS::createYeff,"model thornton_ning");
      registry.registerProperty("Geff", &MODEL_PARAMS::createGeff,"model thornton_ning");
      registry.registerProperty("betaeff", &MODEL_PARAMS::createBetaEff,"model thornton_ning");
      registry.registerProperty("gamma_surf", &MODEL_PARAMS::createSurfaceEnergy, "model thornton_ning");
      registry.registerProperty("yield_ratio", &MODEL_PARAMS::createYieldRatio, "model thornton_ning");
      // registry.registerProperty("coeffRestLog", &MODEL_PARAMS::createCoeffRestLog);

      registry.connect("Yeff", Yeff,"model thornton_ning");
      registry.connect("Geff", Geff,"model thornton_ning");
      registry.connect("betaeff", betaeff,"model thornton_ning");
      registry.connect("gamma_surf", gamma_surf,"model thornton_ning");
      registry.connect("yield_ratio", yield_ratio, "model thornton_ning");
      // registry.connect("coeffRestLog", coeffRestLog,"model thornton_ning_oblique");
    }

    // effective exponent for stress-strain relationship

    inline double stressStrainExponent()
    {
      return 1.5;
    }

    /* ------------------------ CALL FUNCTIONS --------------------------------*/

    inline double calculate_fl(double force_old, double fc, int adhesion_flag)
    {

      double fl;

      if (adhesion_flag == 1)
      {
        fl = force_old + 2*(fc-sqrt(fc * (force_old + fc)));
      } else {
        fl = force_old + 2*(fc+sqrt(fc * (force_old + fc)));
      }

      return fl;
    }

    inline double calculate_elastic_force_differential(double a, double E, double fl, double fc)
    {

      double df;
      if (fc != 0)
      {
        double up = 3*sqrt(fl) - 3*sqrt(fc);
        double down = 3*sqrt(fl) - sqrt(fc);
        df = 2*E*a*up/down;
      } else {
        df = 2*E*a;
      }

      return df;
    }

    inline double calculate_plastic_force_differential(double a_yield, double E, double fl, double fc, double yield_stress, double R)
    {

      double df;
      if (fc != 0)
      {
        double up = 3 * M_PI * R * yield_stress *sqrt(fl) - 2 * a_yield * E *  sqrt(fc);
        double down = 3*sqrt(fl) - sqrt(fc);
        df = up/down;
      } else {
        df = 2*E*a_yield;
      }

      return df;
    }

    inline double calculate_a(double fl, double E, double R, double delta, double gamma_s)
    {
      double a;

      if (gamma_s == 0)
      {
        a = sqrt(R * delta);
      } else {
        a = pow((3. * R * fl)/(4. * E), (1./3.));
      }

      return a;
    }

    inline double calculate_delta_f_ratio(double a, double ac)
    {
      double a_over_ac = a/ac;
      double daf = pow(3, (1./3.))*pow(a_over_ac, 2)*(1-(4./3.)*pow(a_over_ac, -(3./2.)));

      return daf;
    }

    inline void branch_2_force_calc(double force_old, double fc, int adhesion_flag, double E, double rp, double delta, double dn, double gamma_s, double *f_df)
    {
      double fl_r = calculate_fl(force_old, fc, adhesion_flag);
      double a = calculate_a(fl_r, E, rp, delta, gamma_s);

      /* --------------- For checking the model --------------------*/
      // double ac = calculate_a(fc, E, rp, delta, gamma_s);
      // double daf = calculate_delta_f_ratio(a, ac);
      /* --------------- End model checking ------------------------*/

      double df = calculate_elastic_force_differential(a, E, fl_r, fc) * dn;

      f_df[0] = force_old + df;
      f_df[1] = df;

    }

    /* ------------------------- END OF CALL FUNCTIONS ------------------------*/

    inline void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
      const int itype = sidata.itype;
      const int jtype = sidata.jtype;
      double ri = sidata.radi;
      double rj = sidata.radj;
      double reff=sidata.is_wall ? sidata.radi : (ri*rj/(ri+rj));
#ifdef SUPERQUADRIC_ACTIVE_FLAG
      if(sidata.is_non_spherical && atom->superquadric_flag)
        reff = sidata.reff;
#endif
      double meff=sidata.meff;
      double Eeff = Yeff[itype][jtype];
      const double gamma_s = gamma_surf[itype][jtype];
      const double sqrtFiveOverSix = 0.91287092917527685576161630466800355658790782499663875;

      // load history
      if(sidata.contact_flags) *sidata.contact_flags |= CONTACT_NORMAL_MODEL;
      double * const history = &sidata.contact_history[history_offset];
      double * const kc_history = &sidata.contact_history[kc_offset];
      double * const fo_history = &sidata.contact_history[fo_offset];

      int virgin_loading_flag = history[0];
      double delta_old = history[1];
      double delta_max = history[2];
      double force_old = history[3];
      double force_max = history[4];
      int adhesion_flag = history[5];
      double detaching_delta = history[6];
      int detaching_flag = history[7];
      double detaching_force = history[8];
      int yielding_flag = history[9];

      double fc = gamma_s * reff * M_PI * (3./2.);
      fo_history[0] = fc;
      const double delta = sidata.deltan;
      const double a_yield = ri * yield_ratio[itype];
      const double yield_stress = ((2.*Eeff*a_yield)/(M_PI*reff)) - sqrt(2.*gamma_s*Eeff/(M_PI*a_yield));
      const double limit_yield_stress = pow(2.*Eeff*Eeff*gamma_s/(M_PI*M_PI*reff), (1./3.));

      int plastic_from_start = 0;

      if (yield_stress < limit_yield_stress) plastic_from_start = 1;

      if (yield_stress <= 0) error->all(FLERR,"Invalid yield stress, please check surface energy and yield ratio!");

      double f, fl, fl_max;                                                                                                 //forces
      double df;                                                                                                            // differentials
      double a;                                                                                                             // areas, plastic radius

      // Initialization at contact

      double k_t, gammat;
      const double sqrtval = sqrt(reff*sidata.deltan);
      // double coeffRestLogChosen = coeffRestLog[itype][jtype];

      k_t = 8.*Geff[itype][jtype]*sqrtval;
      gammat = 2.*sqrtFiveOverSix*betaeff[itype][jtype]*sqrt(k_t*meff);

      if ((virgin_loading_flag != 1))
      {
        force_old = -(8./9.)*fc;
        delta_old = 0;
        virgin_loading_flag = 1;                                                                                            // update virgin loading flag
        force_max = force_old;
        delta_max = 0;
        adhesion_flag = 0;
        detaching_delta = 0;
        detaching_flag = 0;
        detaching_force = -(5./9.)*fc;
        yielding_flag = 0;
      }

      if (detaching_flag == 1)
        {
          if (delta >= detaching_delta)
          {
            detaching_flag = 0;
            delta_old = detaching_delta;
            force_old = detaching_force;
            goto UNLOADING_RELOADING;

          } else {
            f = 0;
          }
        } else {
          UNLOADING_RELOADING:
          const double dn = delta - delta_old;

          if (delta >= delta_max)                                                                                           // loading or unloading/reloading
          {
            delta_max = delta;

            fl = calculate_fl(force_old, fc, 0);
            a = calculate_a(fl, Eeff, reff, delta, gamma_s);
            if (yielding_flag == 0) if (a >= a_yield) yielding_flag = 1;

            if (plastic_from_start == 1 || yielding_flag == 1)                                                              // if particle already yielded go to plastic, else load elastically
            {

              df = calculate_plastic_force_differential(a_yield, Eeff, fl, fc, yield_stress, reff) * dn;
              f = force_old + df;

            } else {

              double f_df[2];
              branch_2_force_calc(force_old, fc, adhesion_flag, Eeff, reff, delta, dn, gamma_s, f_df);
              f = f_df[0];
              df = f_df[1];

            }
          } else {
            if (yielding_flag == 0)                                                                                           // if particle already yielded, the Reff must be modified, if not we use the original Reff
            {

              double f_df[2];
              branch_2_force_calc(force_old, fc, adhesion_flag, Eeff, reff, delta, dn, gamma_s, f_df);
              f = f_df[0];
              df = f_df[1];

            } else {

              // Calculate modified R and Fc based on that

              fl_max = calculate_fl(force_max, fc, 0);
              reff = reff * fl_max / (force_max + sqrt(4*fc*fl_max));
              fc = (3./2.) * M_PI * gamma_s * reff;
              fo_history[0] = -fc;
              double f_df[2];
              branch_2_force_calc(force_old, fc, adhesion_flag, Eeff, reff, delta, dn, gamma_s, f_df);
              f = f_df[0];
              df = f_df[1];

            }
          }

          if  (gamma_s == 0 && f < 0)                                                                                           // the case of no cohesion
          {
            f = 0;
            detaching_flag = 1;
            detaching_delta = delta;
            detaching_force = 0;
          }

          if (f < -fc) //to prevent stepping into nothingness
          {
            double f_df[2];
            adhesion_flag = 1 - adhesion_flag;
            branch_2_force_calc(force_old, fc, adhesion_flag, Eeff, reff, delta, dn, gamma_s, f_df);
            f = f_df[0];
            df = f_df[1];
          }

          if ((adhesion_flag == 1) && (f + 3*df >= -(5./9.)*fc)) // to make sure particles detach and later reload from same point
          {
            f = 0;

            detaching_flag = 1;
            detaching_delta = delta_old;
            detaching_force = force_old;
          }
        }
      // put force where it belongs
      sidata.Fn = f;
      sidata.kt = k_t;
      kc_history[0] = 0.0;
      sidata.gammat = gammat;
      //sidata.detaching_flag = detaching_flag;

      force_old = f;
      force_max = (f > force_max) ? f : force_max;

      // save history values
      history[0] = virgin_loading_flag;
      history[1] = delta;
      history[2] = delta_max;
      history[3] = f;
      history[4] = force_max;
      history[5] = adhesion_flag;
      history[6] = detaching_delta;
      history[7] = detaching_flag;
      history[8] = detaching_force;
      history[9] = yielding_flag;

      // apply normal force

      if(sidata.is_wall) {
        const double Fn_ = f * sidata.area_ratio;
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

    void surfacesClose(SurfacesCloseData & scdata, ForceData&, ForceData&){
      if(scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_NORMAL_MODEL;
      double * const history = &scdata.contact_history[history_offset];
      history[0] = 0.0;
      history[1] = 0.0;
      history[2] = 0.0;
      history[3] = 0.0;
      history[4] = 0.0;
      history[5] = 0.0;
      history[6] = 0.0;
      history[7] = 0.0;
      history[8] = 0.0;
      history[9] = 0.0;
    }
    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

  protected:
    double ** Yeff;
    double ** Geff;
    double ** betaeff;
    double ** gamma_surf;
    double * yield_ratio;
    double ** coeffRestLog;
    int history_offset;
    int kc_offset;
    int fo_offset;

    bool tangential_damping;
    bool limitForce;
    bool displayedSettings;
    class ContactModelBase *cmb;
  };

}

}
#endif
#endif
