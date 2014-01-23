/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
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
NORMAL_MODEL(HERTZ_STIFFNESS,hertz/stiffness,4)
#else
#ifndef NORMAL_MODEL_HERTZ_STIFFNESS_H_
#define NORMAL_MODEL_HERTZ_STIFFNESS_H_
#include "contact_models.h"
#include "global_properties.h"
#include "math.h"

namespace ContactModels
{
  template<typename Style>
  class NormalModel<HERTZ_STIFFNESS, Style> : protected Pointers
  {
  public:
    static const int MASK = CM_REGISTER_SETTINGS | CM_CONNECT_TO_PROPERTIES | CM_COLLISION;

    NormalModel(LAMMPS * lmp, IContactHistorySetup*) : Pointers(lmp),
      k_n(NULL),
      k_t(NULL),
      gamma_n(NULL),
      gamma_t(NULL),
      tangential_damping(false)
    {
    }

    void registerSettings(Settings & settings)
    {
      settings.registerOnOff("tangential_damping", tangential_damping, true);
    }

    void connectToProperties(PropertyRegistry & registry)
    {
      registry.registerProperty("k_n", &MODEL_PARAMS::createKn);
      registry.registerProperty("k_t", &MODEL_PARAMS::createKt);
      registry.registerProperty("gamma_n", &MODEL_PARAMS::createGamman);
      registry.registerProperty("gamma_t", &MODEL_PARAMS::createGammat);

      registry.connect("k_n", k_n);
      registry.connect("k_t", k_t);
      registry.connect("gamma_n", gamma_n);
      registry.connect("gamma_t", gamma_t);
    }

    inline void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces)
    {
      const int itype = cdata.itype;
      const int jtype = cdata.jtype;
      double meff=cdata.meff;
      double reff = cdata.is_wall ? cdata.radi : (cdata.radi*cdata.radj/(cdata.radi+cdata.radj));

      const double polyhertz = sqrt(reff*cdata.deltan);
      double kn = polyhertz*k_n[itype][jtype];
      double kt = polyhertz*k_t[itype][jtype];
      double gamman = polyhertz*meff*gamma_n[itype][jtype];
      double gammat = polyhertz*meff*gamma_t[itype][jtype];

      if(!tangential_damping) gammat = 0.0;

      // convert Kn and Kt from pressure units to force/distance^2
      kn /= force->nktv2p;
      kt /= force->nktv2p;

      const double Fn_damping = -gamman*cdata.vn;
      const double Fn_contact = kn*(cdata.radsum-cdata.r);
      const double Fn = Fn_contact + Fn_damping;

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
  };
}
#endif // NORMAL_MODEL_HERTZ_STIFFNESS_H_
#endif
