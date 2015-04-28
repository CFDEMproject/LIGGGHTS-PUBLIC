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

    Christoph Kloss (DCS Computing GmbH, Linz, JKU Linz)
    Richard Berger (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef CONTACT_MODELS_H_
#define CONTACT_MODELS_H_

#define STATIC_ASSERT(X)

#include "pointers.h"
#include "lammps.h"
#include "contact_interface.h"
#include "property_registry.h"
#include "pair_gran.h"
#include "settings.h"
#include "contact_model_constants.h"

using namespace LAMMPS_NS;

namespace LIGGGHTS {

namespace ContactModels
{
  static const int CM_REGISTER_SETTINGS     = 1 << 0;
  static const int CM_CONNECT_TO_PROPERTIES = 1 << 1;
  static const int CM_BEGIN_PASS            = 1 << 2;
  static const int CM_END_PASS              = 1 << 3;
  static const int CM_SURFACES_INTERSECT        = 1 << 4;
  static const int CM_SURFACES_CLOSE        = 1 << 5;

  static const int CONTACT_NORMAL_MODEL      = 1 << 0;
  static const int CONTACT_COHESION_MODEL    = 1 << 1;
  static const int CONTACT_TANGENTIAL_MODEL  = 1 << 2;
  static const int CONTACT_ROLLING_MODEL     = 1 << 3;
  static const int CONTACT_SURFACE_MODEL     = 1 << 4;
  static const int CONTACT_FIX               = 1 << 31;

  template
  <
    int M,
    int T = TANGENTIAL_NO_HISTORY,
    int C = COHESION_OFF,
    int R = ROLLING_OFF,
    int S = SURFACE_DEFAULT
  >
  struct GranStyle
  {
  public:
    static const int MODEL = M;
    static const int TANGENTIAL = T;
    static const int COHESION = C;
    static const int ROLLING = R;
    static const int SURFACE = S;
    static const int64_t HASHCODE =
        (((int64_t)M)) |
        (((int64_t)T) << 4) |
        (((int64_t)C) << 8) |
        (((int64_t)R) << 12) |
        (((int64_t)S) << 16);
  };

  int64_t generate_gran_hashcode(int model, int tangential, int cohesion, int rolling, int surface);

  template<int Model>
  class SurfaceModel {
  public:
    SurfaceModel(LAMMPS * lmp, IContactHistorySetup * hsetup);
    inline void beginPass(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces);
    inline void endPass(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces);
    inline void registerSettings(Settings & settings);
    inline void connectToProperties(PropertyRegistry & registry);
    inline bool checkSurfaceIntersect(SurfacesIntersectData & sidata);
    inline void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces);
    inline void surfacesClose(SurfacesCloseData & scdata, ForceData & i_forces, ForceData & j_forces);
  };

  template<int Model>
  class NormalModel {
  public:
    NormalModel(LAMMPS * lmp, IContactHistorySetup * hsetup);
    inline void beginPass(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces);
    inline void endPass(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces);
    inline void registerSettings(Settings & settings);
    inline void connectToProperties(PropertyRegistry & registry);
    inline void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces);
    inline void surfacesClose(SurfacesCloseData & scdata, ForceData & i_forces, ForceData & j_forces);

    inline double stressStrainExponent();
  };

  template<int Model>
  class TangentialModel {
  public:
    TangentialModel(LAMMPS * lmp, IContactHistorySetup * hsetup);
    inline void beginPass(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces);
    inline void endPass(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces);
    inline void registerSettings(Settings & settings);
    inline void connectToProperties(PropertyRegistry & registry);
    inline void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces);
    inline void surfacesClose(SurfacesCloseData & scdata, ForceData & i_forces, ForceData & j_forces);
  };

  template<int Model>
  class CohesionModel {
  public:
    CohesionModel(LAMMPS * lmp, IContactHistorySetup * hsetup);
    void beginPass(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces);
    void endPass(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces);
    void registerSettings(Settings & settings);
    void connectToProperties(PropertyRegistry & registry);
    void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces);
    void surfacesClose(SurfacesCloseData & scdata, ForceData & i_forces, ForceData & j_forces);
  };

  template<int Model>
  class RollingModel {
  public:
    RollingModel(LAMMPS * lmp, IContactHistorySetup * hsetup);
    void beginPass(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces);
    void endPass(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces);
    void registerSettings(Settings & settings);
    void connectToProperties(PropertyRegistry & registry);
    void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces);
    void surfacesClose(SurfacesCloseData & scdata, ForceData & i_forces, ForceData & j_forces);
  };

  template<typename Style>
  class ContactModel {
  private:
    SurfaceModel<Style::SURFACE> surfaceModel;
    NormalModel<Style::MODEL> normalModel;
    CohesionModel<Style::COHESION> cohesionModel;
    TangentialModel<Style::TANGENTIAL> tangentialModel;
    RollingModel<Style::ROLLING> rollingModel;

  public:
    static const int64_t STYLE_HASHCODE = Style::HASHCODE;
    static const int MASK = SurfaceModel<Style::SURFACE>::MASK |
                            NormalModel<Style::MODEL>::MASK |
                            CohesionModel<Style::COHESION>::MASK |
                            TangentialModel<Style::TANGENTIAL>::MASK |
                            RollingModel<Style::ROLLING>::MASK;

    static const int HANDLE_REGISTER_SETTINGS     = MASK & CM_REGISTER_SETTINGS;
    static const int HANDLE_CONNECT_TO_PROPERTIES = MASK & CM_CONNECT_TO_PROPERTIES;
    static const int HANDLE_BEGIN_PASS            = MASK & CM_BEGIN_PASS;
    static const int HANDLE_END_PASS              = MASK & CM_END_PASS;
    static const int HANDLE_SURFACES_INTERSECT    = MASK & CM_SURFACES_INTERSECT;
    static const int HANDLE_SURFACES_CLOSE        = MASK & CM_SURFACES_CLOSE;

    ContactModel(LAMMPS * lmp, IContactHistorySetup * hsetup) :
      surfaceModel(lmp, hsetup),
      normalModel(lmp, hsetup),
      cohesionModel(lmp, hsetup),
      tangentialModel(lmp, hsetup),
      rollingModel(lmp, hsetup)
    {
    }

    int64_t hashcode()
    { return STYLE_HASHCODE; }

    inline void registerSettings(Settings & settings)
    {
      surfaceModel.registerSettings(settings);
      normalModel.registerSettings(settings);
      cohesionModel.registerSettings(settings);
      tangentialModel.registerSettings(settings);
      rollingModel.registerSettings(settings);
    }

    inline void connectToProperties(PropertyRegistry & registry)
    {
      surfaceModel.connectToProperties(registry);
      normalModel.connectToProperties(registry);
      cohesionModel.connectToProperties(registry);
      tangentialModel.connectToProperties(registry);
      rollingModel.connectToProperties(registry);
    }

    inline void beginPass(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
      surfaceModel.beginPass(sidata, i_forces, j_forces);
      normalModel.beginPass(sidata, i_forces, j_forces);
      cohesionModel.beginPass(sidata, i_forces, j_forces);
      tangentialModel.beginPass(sidata, i_forces, j_forces);
      rollingModel.beginPass(sidata, i_forces, j_forces);
    }

    inline void endPass(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
      rollingModel.endPass(sidata, i_forces, j_forces);
      tangentialModel.endPass(sidata, i_forces, j_forces);
      cohesionModel.endPass(sidata, i_forces, j_forces);
      normalModel.endPass(sidata, i_forces, j_forces);
      surfaceModel.endPass(sidata, i_forces, j_forces);
    }

    inline double stressStrainExponent()
    {
      return normalModel.stressStrainExponent();
    }

    inline bool checkSurfaceIntersect(SurfacesIntersectData & sidata)
    {
      return surfaceModel.checkSurfaceIntersect(sidata);
    }

    inline void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
      surfaceModel.surfacesIntersect(sidata, i_forces, j_forces);
      normalModel.surfacesIntersect(sidata, i_forces, j_forces);
      cohesionModel.surfacesIntersect(sidata, i_forces, j_forces);
      tangentialModel.surfacesIntersect(sidata, i_forces, j_forces);
      rollingModel.surfacesIntersect(sidata, i_forces, j_forces);
    }

    inline void surfacesClose(SurfacesCloseData & scdata, ForceData & i_forces, ForceData & j_forces)
    {
      surfaceModel.surfacesClose(scdata, i_forces, j_forces);
      normalModel.surfacesClose(scdata, i_forces, j_forces);
      cohesionModel.surfacesClose(scdata, i_forces, j_forces);
      tangentialModel.surfacesClose(scdata, i_forces, j_forces);
      rollingModel.surfacesClose(scdata, i_forces, j_forces);
    }
  };

  template<>
  class CohesionModel<COHESION_OFF> : protected Pointers
  {
  public:
    static const int MASK = 0;

    CohesionModel(LAMMPS * lmp, IContactHistorySetup*) : Pointers(lmp) {}
    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void connectToProperties(PropertyRegistry&){}
    void registerSettings(Settings&){}
    void surfacesIntersect(SurfacesIntersectData&, ForceData&, ForceData&){}
    void surfacesClose(SurfacesCloseData&, ForceData&, ForceData&){}
  };

  template<>
  class RollingModel<ROLLING_OFF> : protected Pointers
  {
  public:
    static const int MASK = 0;

    RollingModel(LAMMPS * lmp, IContactHistorySetup*) : Pointers(lmp) {}
    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void connectToProperties(PropertyRegistry&){}
    void registerSettings(Settings&){}
    void surfacesIntersect(SurfacesIntersectData&, ForceData&, ForceData&){}
    void surfacesClose(SurfacesCloseData&, ForceData&, ForceData&){}
  };

  class Factory {
    typedef std::map<std::string, int> ModelTable;

    ModelTable surface_models;
    ModelTable normal_models;
    ModelTable tangential_models;
    ModelTable cohesion_models;
    ModelTable rolling_models;

    Factory();
    Factory(const Factory &){}
  public:
    static Factory & instance();
    static int64_t select(int & narg, char ** & args);

    void addNormalModel(const std::string & name, int identifier);
    void addTangentialModel(const std::string & name, int identifier);
    void addCohesionModel(const std::string & name, int identifier);
    void addRollingModel(const std::string & name, int identifier);
    void addSurfaceModel(const std::string & name, int identifier);

  private:
    int64_t select_model(int & narg, char ** & args);
  };

}
}

#endif /* CONTACT_MODELS_H_ */
