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
NORMAL_MODEL(HERTZ,hertz,3)
#else
#ifndef NORMAL_MODEL_HERTZ_H_
#define NORMAL_MODEL_HERTZ_H_
#include "global_properties.h"
#include <math.h>

namespace LIGGGHTS {

namespace ContactModels
{
  template<>
  class NormalModel<HERTZ> : protected Pointers
  {
  public:
    static const int MASK = CM_REGISTER_SETTINGS | CM_CONNECT_TO_PROPERTIES | CM_SURFACES_INTERSECT;

    NormalModel(LAMMPS * lmp, IContactHistorySetup* hsetup, class ContactModelBase *c) : Pointers(lmp),
      Yeff(NULL),
      Geff(NULL),
      betaeff(NULL),
      limitForce(false),
      displayedSettings(false),
      heating(false),
      heating_track(false),
      elastic_potential_offset_(-1),
      elasticpotflag_(false),
      disable_when_bonded_(false),
      bond_history_offset_(-1),
      cmb(c)
    {
      
    }

    void registerSettings(Settings & settings)
    {
      
      settings.registerOnOff("tangential_damping", tangential_damping, true);
      settings.registerOnOff("limitForce", limitForce);
      settings.registerOnOff("heating_normal_hertz",heating,false);
      settings.registerOnOff("heating_tracking",heating_track,false);
      settings.registerOnOff("computeElasticPotential", elasticpotflag_, false);
      settings.registerOnOff("disableNormalWhenBonded", disable_when_bonded_, false);
      //TODO error->one(FLERR,"TODO here also check if right surface model used");
    }

    inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb)
    {
      if (elasticpotflag_)
      {
        elastic_potential_offset_ = cmb->get_history_offset("elastic_potential");
        if (elastic_potential_offset_ == -1)
        {
          elastic_potential_offset_ = hsetup->add_history_value("elastic_potential", "0");
          cmb->add_history_offset("elastic_potential", elastic_potential_offset_);
        }
      }
      if (disable_when_bonded_)
      {
        bond_history_offset_ = cmb->get_history_offset("bond_contactflag");
        if (bond_history_offset_ < 0)
          error->one(FLERR, "Could not find bond history offset");
      }
    }

    void connectToProperties(PropertyRegistry & registry)
    {
      
      registry.registerProperty("Yeff", &MODEL_PARAMS::createYeff,"model hertz");
      registry.registerProperty("Geff", &MODEL_PARAMS::createGeff,"model hertz");
      registry.registerProperty("betaeff", &MODEL_PARAMS::createBetaEff,"model hertz");

      registry.connect("Yeff", Yeff,"model hertz");
      registry.connect("Geff", Geff,"model hertz");
      registry.connect("betaeff", betaeff,"model hertz");

    }

    // effective exponent for stress-strain relationship
    
    inline double stressStrainExponent()
    {
      return 1.5;
    }

    inline void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {

      //if(!sidata.is_wall) // bond_history_offset_ >= 0)
      //{
          //double * const bond_contact_flag = &sidata.contact_history[bond_history_offset_];
          //if(!MathExtraLiggghts::compDouble(bond_contact_flag[0],0.))
      //       return;
      //}

      const int itype = sidata.itype;
      const int jtype = sidata.jtype;
      const double radi = sidata.radi;
      const double radj = sidata.radj;
      double reff = sidata.is_wall ? radi : (radi*radj/(radi+radj));

      #ifdef SUPERQUADRIC_ACTIVE_FLAG
        if(sidata.is_non_spherical) {
          if(sidata.is_wall)
            reff = MathExtraLiggghtsNonspherical::get_effective_radius_wall(sidata, atom->roundness[sidata.i], error);
          else
            reff = MathExtraLiggghtsNonspherical::get_effective_radius(sidata, atom->roundness[sidata.i], atom->roundness[sidata.j], error);
        }
      #endif
      const double meff=sidata.meff;

      const double sqrtval = sqrt(reff*sidata.deltan);

      const double Sn=2.*Yeff[itype][jtype]*sqrtval;
      const double St=8.*Geff[itype][jtype]*sqrtval;

      double kn=4./3.*Yeff[itype][jtype]*sqrtval;
      double kt=St;
      const double sqrtFiveOverSix = 0.91287092917527685576161630466800355658790782499663875;
      const double gamman=-2.*sqrtFiveOverSix*betaeff[itype][jtype]*sqrt(Sn*meff);
      const double gammat= tangential_damping ? -2.*sqrtFiveOverSix*betaeff[itype][jtype]*sqrt(St*meff) : 0.0;
      
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
      const double Fn_contact = kn*sidata.deltan;
      double Fn = Fn_damping + Fn_contact;

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

      #ifdef NONSPHERICAL_ACTIVE_FLAG
          double torque_i[3] = {0., 0., 0.};
          double Fn_i[3] = { Fn * sidata.en[0], Fn * sidata.en[1], Fn * sidata.en[2]};
          if(sidata.is_non_spherical) {
            double xci[3];
            vectorSubtract3D(sidata.contact_point, atom->x[sidata.i], xci);
            vectorCross3D(xci, Fn_i, torque_i);
          }
          
      #endif

      if(heating)
      {
        sidata.P_diss += fabs(Fn_damping*sidata.vn); //fprintf(screen,"  contrib %f\n",fabs(Fn_damping*sidata.vn));}
        if(heating_track && sidata.is_wall) cmb->tally_pw(fabs(Fn_damping*sidata.vn),sidata.i,jtype,0);
        if(heating_track && !sidata.is_wall) cmb->tally_pp(fabs(Fn_damping*sidata.vn),sidata.i,sidata.j,0);
      }

      // apply normal force
      if (!disable_when_bonded_ || MathExtraLiggghts::compDouble(sidata.contact_history[bond_history_offset_], 0.0, 1e-5))
      {
        // compute increment in elastic potential
        if (elasticpotflag_ && sidata.computeflag && sidata.shearupdate)
          sidata.contact_history[elastic_potential_offset_] += -update->dt*sidata.vn*Fn_contact;

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
    }

    void surfacesClose(SurfacesCloseData&, ForceData&, ForceData&){}
    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

  protected:
    double ** Yeff;
    double ** Geff;
    double ** betaeff;

    bool tangential_damping;
    bool limitForce;
    bool displayedSettings;
    bool heating;
    bool heating_track;
    int elastic_potential_offset_;
    bool elasticpotflag_;
    bool disable_when_bonded_;
    int bond_history_offset_;
    class ContactModelBase *cmb;

  };

}

}
#endif
#endif
