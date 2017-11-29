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

#ifdef COHESION_MODEL
COHESION_MODEL(COHESION_SJKR2,sjkr2,2)
#else
#ifndef COHESION_MODEL_SJKR2_H_
#define COHESION_MODEL_SJKR2_H_

#include "pointers.h"
#include "contact_models.h"
#include "cohesion_model_base.h"
#include <cmath>

namespace LIGGGHTS {
namespace ContactModels {
  using namespace std;
  using namespace LAMMPS_NS;

  template<>
  class CohesionModel<COHESION_SJKR2> : public CohesionModelBase {
  public:
    CohesionModel(LAMMPS * lmp, IContactHistorySetup * hsetup, class ContactModelBase * c) :
        CohesionModelBase(lmp, hsetup, c),
        cohEnergyDens(NULL)
    {
        
    }

    void registerSettings(Settings& settings) 
    {
        settings.registerOnOff("tangential_reduce",tangentialReduce_,false);
    }

    inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb) {}

    void connectToProperties(PropertyRegistry & registry)
    {
        registry.registerProperty("cohEnergyDens", &MODEL_PARAMS::createCohesionEnergyDensity);
        registry.connect("cohEnergyDens", cohEnergyDens,"cohesion_model sjkr2");

        // error checks on coarsegraining
        if(force->cg_active())
            error->cg(FLERR,"cohesion model sjkr2");
    }

    void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces) 
    {
      //r is the distance between the sphere's centers
      const double r = sidata.r;
      const double ri = sidata.radi;
      const double rj = sidata.radj;
      double Acont;

      if(sidata.is_wall)
        Acont = M_PI * 2. * ri * (ri - r) * sidata.area_ratio;
      else
        Acont = M_PI * 2. * (2.*ri*rj/(ri+rj)) * (ri + rj - r);

      const double Fn_coh = -cohEnergyDens[sidata.itype][sidata.jtype]*Acont;
      if(tangentialReduce_) sidata.Fn += Fn_coh; 

      if(sidata.contact_flags) *sidata.contact_flags |= CONTACT_COHESION_MODEL;

      // apply normal force
      if(sidata.is_wall) {
        const double Fn_ = Fn_coh * sidata.area_ratio;
        i_forces.delta_F[0] += Fn_ * sidata.en[0];
        i_forces.delta_F[1] += Fn_ * sidata.en[1];
        i_forces.delta_F[2] += Fn_ * sidata.en[2];
      } else {
        const double fx = Fn_coh * sidata.en[0];
        const double fy = Fn_coh * sidata.en[1];
        const double fz = Fn_coh * sidata.en[2];

        i_forces.delta_F[0] += fx;
        i_forces.delta_F[1] += fy;
        i_forces.delta_F[2] += fz;

        j_forces.delta_F[0] -= fx;
        j_forces.delta_F[1] -= fy;
        j_forces.delta_F[2] -= fz;
      }
    }

    inline void endSurfacesIntersect(SurfacesIntersectData &sidata, ForceData&, ForceData&) {}
    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

    void surfacesClose(SurfacesCloseData& scdata, ForceData&, ForceData&)
    {
        if(scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_COHESION_MODEL;
    }

  private:
    double ** cohEnergyDens;
    bool tangentialReduce_;
  };
}
}
#endif // COHESION_MODEL_SJKR2_H_
#endif
