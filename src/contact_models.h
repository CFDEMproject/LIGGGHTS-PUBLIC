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
  static const int CM_COLLISION             = 1 << 4;
  static const int CM_NO_COLLISION          = 1 << 5;

  static const int TOUCH_NORMAL_MODEL      = 1 << 0;
  static const int TOUCH_COHESION_MODEL    = 1 << 1;
  static const int TOUCH_TANGENTIAL_MODEL  = 1 << 2;
  static const int TOUCH_ROLLING_MODEL     = 1 << 3;
  static const int TOUCH_SURFACE_MODEL     = 1 << 4;
  static const int TOUCH_FIX               = 1 << 31;

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
    inline void beginPass(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces);
    inline void endPass(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces);
    inline void registerSettings(Settings & settings);
    inline void connectToProperties(PropertyRegistry & registry);
    inline void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces);
    inline void noCollision(ContactData & cdata, ForceData & i_forces, ForceData & j_forces);
  };

  template<int Model>
  class NormalModel {
  public:
    NormalModel(LAMMPS * lmp, IContactHistorySetup * hsetup);
    inline void beginPass(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces);
    inline void endPass(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces);
    inline void registerSettings(Settings & settings);
    inline void connectToProperties(PropertyRegistry & registry);
    inline void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces);
    inline void noCollision(ContactData & cdata, ForceData & i_forces, ForceData & j_forces);

    inline double stressStrainExponent();
  };

  template<int Model>
  class TangentialModel {
  public:
    TangentialModel(LAMMPS * lmp, IContactHistorySetup * hsetup);
    inline void beginPass(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces);
    inline void endPass(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces);
    inline void registerSettings(Settings & settings);
    inline void connectToProperties(PropertyRegistry & registry);
    inline void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces);
    inline void noCollision(ContactData & cdata, ForceData & i_forces, ForceData & j_forces);
  };

  template<int Model>
  class CohesionModel {
  public:
    CohesionModel(LAMMPS * lmp, IContactHistorySetup * hsetup);
    void beginPass(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces);
    void endPass(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces);
    void registerSettings(Settings & settings);
    void connectToProperties(PropertyRegistry & registry);
    void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces);
    void noCollision(ContactData & cdata, ForceData & i_forces, ForceData & j_forces);
  };

  template<int Model>
  class RollingModel {
  public:
    RollingModel(LAMMPS * lmp, IContactHistorySetup * hsetup);
    void beginPass(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces);
    void endPass(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces);
    void registerSettings(Settings & settings);
    void connectToProperties(PropertyRegistry & registry);
    void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces);
    void noCollision(ContactData & cdata, ForceData & i_forces, ForceData & j_forces);
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

    static const int HANDLE_REGISTER_SETTINGS = MASK & CM_REGISTER_SETTINGS;
    static const int HANDLE_CONNECT_TO_PROPERTIES = MASK & CM_CONNECT_TO_PROPERTIES;
    static const int HANDLE_BEGIN_PASS = MASK & CM_BEGIN_PASS;
    static const int HANDLE_END_PASS = MASK & CM_END_PASS;
    static const int HANDLE_COLLISION = MASK & CM_COLLISION;
    static const int HANDLE_NO_COLLISION = MASK & CM_NO_COLLISION;

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

    inline void beginPass(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces)
    {
      surfaceModel.beginPass(cdata, i_forces, j_forces);
      normalModel.beginPass(cdata, i_forces, j_forces);
      cohesionModel.beginPass(cdata, i_forces, j_forces);
      tangentialModel.beginPass(cdata, i_forces, j_forces);
      rollingModel.beginPass(cdata, i_forces, j_forces);
    }

    inline void endPass(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces)
    {
      rollingModel.endPass(cdata, i_forces, j_forces);
      tangentialModel.endPass(cdata, i_forces, j_forces);
      cohesionModel.endPass(cdata, i_forces, j_forces);
      normalModel.endPass(cdata, i_forces, j_forces);
      surfaceModel.endPass(cdata, i_forces, j_forces);
    }

    inline double stressStrainExponent()
    {
      return normalModel.stressStrainExponent();
    }

    inline void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces)
    {
      surfaceModel.collision(cdata, i_forces, j_forces);
      normalModel.collision(cdata, i_forces, j_forces);
      cohesionModel.collision(cdata, i_forces, j_forces);
      tangentialModel.collision(cdata, i_forces, j_forces);
      rollingModel.collision(cdata, i_forces, j_forces);
    }

    inline void noCollision(ContactData & cdata, ForceData & i_forces, ForceData & j_forces)
    {
      surfaceModel.noCollision(cdata, i_forces, j_forces);
      normalModel.noCollision(cdata, i_forces, j_forces);
      cohesionModel.noCollision(cdata, i_forces, j_forces);
      tangentialModel.noCollision(cdata, i_forces, j_forces);
      rollingModel.noCollision(cdata, i_forces, j_forces);
    }
  };

  template<>
  class CohesionModel<COHESION_OFF> : protected Pointers
  {
  public:
    static const int MASK = 0;

    CohesionModel(LAMMPS * lmp, IContactHistorySetup*) : Pointers(lmp) {}
    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}
    void connectToProperties(PropertyRegistry&){}
    void registerSettings(Settings&){}
    void collision(CollisionData&, ForceData&, ForceData&){}
    void noCollision(ContactData&, ForceData&, ForceData&){}
  };

  template<>
  class RollingModel<ROLLING_OFF> : protected Pointers
  {
  public:
    static const int MASK = 0;

    RollingModel(LAMMPS * lmp, IContactHistorySetup*) : Pointers(lmp) {}
    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}
    void connectToProperties(PropertyRegistry&){}
    void registerSettings(Settings&){}
    void collision(CollisionData&, ForceData&, ForceData&){}
    void noCollision(ContactData&, ForceData&, ForceData&){}
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
