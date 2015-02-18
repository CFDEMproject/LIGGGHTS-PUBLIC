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
#ifdef SURFACE_MODEL
SURFACE_MODEL(SURFACE_DEFAULT,default,0)
#else
#ifndef SURFACE_MODEL_DEFAULT_H_
#define SURFACE_MODEL_DEFAULT_H_
#include "contact_models.h"
#include "math.h"
#include "atom.h"
#include "force.h"
#include "update.h"

namespace LIGGGHTS {
namespace ContactModels
{
  template<>
  class SurfaceModel<SURFACE_DEFAULT> : protected Pointers
  {
  public:
    static const int MASK = CM_COLLISION;

    SurfaceModel(LAMMPS * lmp, IContactHistorySetup*) : Pointers(lmp)
    {
      
    }

    inline void registerSettings(Settings&) {}
    inline void connectToProperties(PropertyRegistry&) {}

    inline void collision(CollisionData & cdata, ForceData&, ForceData&)
    {
      const double enx = cdata.en[0];
      const double eny = cdata.en[1];
      const double enz = cdata.en[2];

      // relative translational velocity
      const double vr1 = cdata.v_i[0] - cdata.v_j[0];
      const double vr2 = cdata.v_i[1] - cdata.v_j[1];
      const double vr3 = cdata.v_i[2] - cdata.v_j[2];

      // normal component
      const double vn = vr1 * enx + vr2 * eny + vr3 * enz;
      const double vn1 = vn * enx;
      const double vn2 = vn * eny;
      const double vn3 = vn * enz;

      // tangential component
      const double vt1 = vr1 - vn1;
      const double vt2 = vr2 - vn2;
      const double vt3 = vr3 - vn3;

      // relative rotational velocity
      const double deltan = cdata.radsum - cdata.r;
      const double dx = cdata.delta[0];
      const double dy = cdata.delta[1];
      const double dz = cdata.delta[2];
      const double rinv = cdata.rinv;
      double wr1, wr2, wr3;

      if(cdata.is_wall) {
        // in case of wall contact, r is the contact radius
        const double cr = cdata.radi - 0.5*cdata.deltan;
        wr1 = cr * cdata.omega_i[0] * rinv;
        wr2 = cr * cdata.omega_i[1] * rinv;
        wr3 = cr * cdata.omega_i[2] * rinv;
        cdata.cri = cr;
      } else {
        const double cri = cdata.radi - 0.5 * deltan;
        const double crj = cdata.radj - 0.5 * deltan;
        wr1 = (cri * cdata.omega_i[0] + crj * cdata.omega_j[0]) * rinv;
        wr2 = (cri * cdata.omega_i[1] + crj * cdata.omega_j[1]) * rinv;
        wr3 = (cri * cdata.omega_i[2] + crj * cdata.omega_j[2]) * rinv;
        cdata.cri = cri;
        cdata.crj = crj;
      }

      // relative velocities
      const double vtr1 = vt1 - (dz * wr2 - dy * wr3);
      const double vtr2 = vt2 - (dx * wr3 - dz * wr1);
      const double vtr3 = vt3 - (dy * wr1 - dx * wr2);

      cdata.vn = vn;
      cdata.deltan = deltan;
      cdata.wr1 = wr1;
      cdata.wr2 = wr2;
      cdata.wr3 = wr3;
      cdata.vtr1 = vtr1;
      cdata.vtr2 = vtr2;
      cdata.vtr3 = vtr3;
    }

    inline void noCollision(ContactData&, ForceData&, ForceData&){}
    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}
  };
}
}
#endif // SURFACE_MODEL_DEFAULT_H_
#endif
