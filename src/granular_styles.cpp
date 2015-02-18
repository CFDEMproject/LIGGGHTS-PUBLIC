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
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */
#include "mpi.h"
#include "ctype.h"
#include "float.h"
#include "limits.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "kspace.h"
#include "update.h"
#include "accelerator_cuda.h"
#include "suffix.h"
#include "atom_masks.h"
#include "memory.h"
#include "error.h"

#include "granular_pair_style.h"
#include "granular_wall.h"

#include "contact_models.h"
#include "style_surface_model.h"
#include "style_normal_model.h"
#include "style_tangential_model.h"
#include "style_cohesion_model.h"
#include "style_rolling_model.h"
#include "pair_gran_base.h"
#include "fix_wall_gran_base.h"

namespace LIGGGHTS {

using namespace ContactModels;

// creates global object, which will register its pair styles during static initialization
struct RegisterGranularStyles {
public:
  RegisterGranularStyles() {
    PairStyles::Factory & pair_factory = PairStyles::Factory::instance();
    Walls::Factory & wall_factory = Walls::Factory::instance();

    // select granular variants based on contact models
    pair_factory.addVariantSelector("gran", ContactModels::Factory::select);
    wall_factory.addVariantSelector("gran", ContactModels::Factory::select);

    // register granular pair styles
    #define GRAN_MODEL(MODEL,TANGENTIAL,COHESION,ROLLING,SURFACE) \
    registerPair<MODEL, TANGENTIAL, COHESION, ROLLING, SURFACE>("gran", pair_factory); \
    registerWall<MODEL, TANGENTIAL, COHESION, ROLLING, SURFACE>("gran", wall_factory);
    #include "style_contact_model.h"
    #undef GRAN_MODEL
  }
private:
  template<typename T>
  static PairStyles::IGranularPairStyle * create_pair_style_instance(LAMMPS* lmp, PairGran* parent) {
    return new T(lmp, parent);
  }

  template<typename T>
  static Walls::IGranularWall * create_wall_style_instance(LAMMPS* lmp, FixWallGran* parent) {
    return new T(lmp, parent);
  }

  template<int MODEL, int TANGENTIAL, int COHESION, int ROLLING, int SURFACE>
  void registerPair(const std::string & name, LIGGGHTS::PairStyles::Factory & factory) {
     typedef GranStyle<MODEL,TANGENTIAL,COHESION,ROLLING,SURFACE> Style;
     typedef ContactModel<Style> CModel;
     typedef PairStyles::Granular<CModel> Type;
     int64_t hashcode = Style::HASHCODE;
     factory.addStyle(name, hashcode, &create_pair_style_instance<Type>);
  }

  template<int MODEL, int TANGENTIAL, int COHESION, int ROLLING, int SURFACE>
  void registerWall(const std::string & name, LIGGGHTS::Walls::Factory & factory) {
     typedef GranStyle<MODEL,TANGENTIAL,COHESION,ROLLING,SURFACE> Style;
     typedef ContactModel<Style> CModel;
     typedef Walls::Granular<CModel> Type;
     int64_t hashcode = Style::HASHCODE;
     factory.addStyle(name, hashcode, &create_wall_style_instance<Type>);
  }
} granStyles;
}
