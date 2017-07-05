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
    (if not contributing author is listed, this file has been contributed
    by the core developer)
    Arno Mayrhofer (DCS Computing GmbH, Linz)

    Copyright 2017-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifndef COHESION_MODEL

#ifndef COHESION_MODEL_BASE_H
#define COHESION_MODEL_BASE_H

#include "contact_models.h"
#include "contact_interface.h"
namespace LIGGGHTS
{
namespace ContactModels
{

class CohesionModelBase : protected Pointers
{
  public:
    CohesionModelBase(LAMMPS * lmp, IContactHistorySetup*, class ContactModelBase *) :
        Pointers(lmp)
    {}

    virtual void registerSettings(Settings& settings) = 0;
    virtual void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb) = 0;
    virtual void connectToProperties(PropertyRegistry&) = 0;
    virtual void surfacesIntersect(SurfacesIntersectData & sidata, ForceData&, ForceData&) = 0;
    virtual void endSurfacesIntersect(SurfacesIntersectData &sidata, ForceData &i_forces, ForceData &j_forces) = 0;
    virtual void surfacesClose(SurfacesCloseData &scdata, ForceData&, ForceData&) = 0;
    virtual void beginPass(SurfacesIntersectData&, ForceData&, ForceData&) = 0;
    virtual void endPass(SurfacesIntersectData&, ForceData&, ForceData&) = 0;
};

  template<int Model>
  class CohesionModel : public CohesionModelBase
  {
  public:
    CohesionModel(LAMMPS * lmp, IContactHistorySetup * hsetup, class ContactModelBase *cmb = 0);
    void beginPass(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces);
    void endPass(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces);
    void registerSettings(Settings & settings);
    void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb);
    void connectToProperties(PropertyRegistry & registry);
    void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces);
    void surfacesClose(SurfacesCloseData & scdata, ForceData & i_forces, ForceData & j_forces);
    void endSurfacesIntersect(SurfacesIntersectData &sidata, ForceData&, ForceData&);
  };

  template<>
  class CohesionModel<COHESION_OFF> : public CohesionModelBase
  {
  public:
    CohesionModel(LAMMPS * lmp, IContactHistorySetup * hsetup, class ContactModelBase *c) : CohesionModelBase(lmp, hsetup, c) {}
    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void connectToProperties(PropertyRegistry&){}
    void registerSettings(Settings&){}
    void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb){}
    void surfacesIntersect(SurfacesIntersectData&, ForceData&, ForceData&){}
    void surfacesClose(SurfacesCloseData&, ForceData&, ForceData&){}
    void endSurfacesIntersect(SurfacesIntersectData &sidata, ForceData&, ForceData&) {}
  };

}
}

#endif

#endif
