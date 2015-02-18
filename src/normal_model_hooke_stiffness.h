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
    static const int MASK = CM_REGISTER_SETTINGS | CM_CONNECT_TO_PROPERTIES | CM_COLLISION;

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

    inline void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces)
    {
      const int itype = cdata.itype;
      const int jtype = cdata.jtype;
      double meff=cdata.meff;

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

      const double Fn_damping = -gamman*cdata.vn;    
      const double Fn_contact = kn*(cdata.radsum-cdata.r);
      double Fn                       = Fn_damping + Fn_contact;

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
        i_forces.delta_F[0] = cdata.Fn * cdata.en[0];
        i_forces.delta_F[1] = cdata.Fn * cdata.en[1];
        i_forces.delta_F[2] = cdata.Fn * cdata.en[2];

        j_forces.delta_F[0] = -i_forces.delta_F[0];
        j_forces.delta_F[1] = -i_forces.delta_F[1];
        j_forces.delta_F[2] = -i_forces.delta_F[2];
      }
    }

    void noCollision(ContactData&, ForceData&, ForceData&){}
    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}

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
