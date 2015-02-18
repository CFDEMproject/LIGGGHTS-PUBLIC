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
#ifdef COHESION_MODEL
COHESION_MODEL(COHESION_SJKR2,sjkr2,2)
#else
#ifndef COHESION_MODEL_SJKR2_H_
#define COHESION_MODEL_SJKR2_H_

#include "pointers.h"
#include "contact_models.h"
#include "math.h"

namespace LIGGGHTS {
namespace ContactModels {
  using namespace std;
  using namespace LAMMPS_NS;

  template<>
  class CohesionModel<COHESION_SJKR2> : protected Pointers {
  public:
    static const int MASK = CM_CONNECT_TO_PROPERTIES | CM_COLLISION;

    CohesionModel(LAMMPS * lmp, IContactHistorySetup*) : Pointers(lmp), cohEnergyDens(NULL)
    {
      
    }

    void registerSettings(Settings&) {}

    void connectToProperties(PropertyRegistry & registry) {
      registry.registerProperty("cohEnergyDens", &MODEL_PARAMS::createCohesionEnergyDensity);
      registry.connect("cohEnergyDens", cohEnergyDens,"cohesion_model sjkr2");

      // error checks on coarsegraining
      if(force->cg_active())
        error->cg(FLERR,"cohesion model sjkr2");
    }

    void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces) 
    {
      //r is the distance between the sphere's centers
      const double r = cdata.r;
      const double ri = cdata.radi;
      const double rj = cdata.radj;
      double Acont;

      if(cdata.is_wall)
        Acont = M_PI * 2. * ri * (ri - r) * cdata.area_ratio;
      else
        Acont = M_PI * 2. * (2.*ri*rj/(ri+rj)) * (ri + rj - r);

      const double Fn_coh = -cohEnergyDens[cdata.itype][cdata.jtype]*Acont;
      cdata.Fn += Fn_coh;

      if(cdata.touch) *cdata.touch |= TOUCH_COHESION_MODEL;

      // apply normal force
      if(cdata.is_wall) {
        const double Fn_ = Fn_coh * cdata.area_ratio;
        i_forces.delta_F[0] += Fn_ * cdata.en[0];
        i_forces.delta_F[1] += Fn_ * cdata.en[1];
        i_forces.delta_F[2] += Fn_ * cdata.en[2];
      } else {
        const double fx = Fn_coh * cdata.en[0];
        const double fy = Fn_coh * cdata.en[1];
        const double fz = Fn_coh * cdata.en[2];

        i_forces.delta_F[0] += fx;
        i_forces.delta_F[1] += fy;
        i_forces.delta_F[2] += fz;

        j_forces.delta_F[0] -= fx;
        j_forces.delta_F[1] -= fy;
        j_forces.delta_F[2] -= fz;
      }
    }

    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}
    void noCollision(ContactData& cdata, ForceData&, ForceData&)
    {
        if(cdata.touch) *cdata.touch &= ~TOUCH_COHESION_MODEL;
    }

  private:
    double ** cohEnergyDens;
  };
}
}
#endif // COHESION_MODEL_SJKR2_H_
#endif
