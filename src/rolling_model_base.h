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

#ifndef ROLLING_MODEL

#ifndef ROLLING_MODEL_BASE_H
#define ROLLING_MODEL_BASE_H

#include "contact_models.h"
#include "contact_interface.h"
namespace LIGGGHTS
{
namespace ContactModels
{

class RollingModelBase : protected Pointers
{
  public:
    RollingModelBase(LAMMPS * lmp, IContactHistorySetup*, class ContactModelBase *) :
        Pointers(lmp)
    {}

    virtual void registerSettings(Settings& settings) = 0;
    virtual void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb) = 0;
    virtual void connectToProperties(PropertyRegistry&) = 0;
    virtual void surfacesIntersect(SurfacesIntersectData & sidata, ForceData&, ForceData&) = 0;
    virtual void surfacesClose(SurfacesCloseData &scdata, ForceData&, ForceData&) = 0;
    virtual void beginPass(SurfacesIntersectData&, ForceData&, ForceData&) = 0;
    virtual void endPass(SurfacesIntersectData&, ForceData&, ForceData&) = 0;
};

  template<int Model>
  class RollingModel : public RollingModelBase
  {
  public:
    RollingModel(LAMMPS * lmp, IContactHistorySetup * hsetup, class ContactModelBase *cmb = 0);
    void beginPass(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces);
    void endPass(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces);
    void registerSettings(Settings & settings);
    void connectToProperties(PropertyRegistry & registry);
    void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces);
    void surfacesClose(SurfacesCloseData & scdata, ForceData & i_forces, ForceData & j_forces);
  };

  template<>
  class RollingModel<ROLLING_OFF> : public RollingModelBase
  {
  public:
    RollingModel(LAMMPS * lmp, IContactHistorySetup * hsetup, class ContactModelBase *c) : RollingModelBase(lmp, hsetup, c) {}
    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void connectToProperties(PropertyRegistry&){}
    void registerSettings(Settings&){}
    void surfacesIntersect(SurfacesIntersectData&, ForceData&, ForceData&){}
    void surfacesClose(SurfacesCloseData&, ForceData&, ForceData&){}
    inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb) {}
  };

}
}

#endif

#endif
