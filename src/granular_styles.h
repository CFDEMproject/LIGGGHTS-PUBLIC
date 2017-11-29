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
    (if no contributing author is listed, this file has been contributed
    by the core developer)

    Christoph Kloss (DCS Computing GmbH, Linz)
    Arno Mayrhofer (DCS Computing GmbH, Linz)
    Richard Berger (JKU Linz)

    Copyright 2015-     DCS Computing GmbH, Linz
    Copyright 2015-     JKU Linz
------------------------------------------------------------------------- */

#ifndef GRANULAR_STYLES_H
#define GRANULAR_STYLES_H

#include <mpi.h>
#include "ctype.h"
#include "float.h"
#include "limits.h"
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

namespace LAMMPS_NS {
    using namespace LIGGGHTS;
    using namespace ContactModels;

    // creates global object, which will register its pair styles during LAMMPS initialization
    class RegisterGranularStyles
    {
    public:
        RegisterGranularStyles()
        {
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
            // default for slow execution
            registerPair<NORMAL_OFF, TANGENTIAL_OFF, COHESION_OFF, ROLLING_OFF, SURFACE_DEFAULT>("gran", pair_factory);
            registerWall<NORMAL_OFF, TANGENTIAL_OFF, COHESION_OFF, ROLLING_OFF, SURFACE_DEFAULT>("gran", wall_factory);
        }

        ~RegisterGranularStyles() {}
    private:
        template<typename T>
        static PairStyles::IGranularPairStyle * create_pair_style_instance(LAMMPS* lmp, PairGran* parent, const int64_t hash)
        {
            return new T(lmp, parent, hash);
        }

        template<typename T>
        static Walls::IGranularWall * create_wall_style_instance(LAMMPS* lmp, FixWallGran* parent, const int64_t hash)
        {
            return new T(lmp, parent, hash);
        }

        template<int MODEL, int TANGENTIAL, int COHESION, int ROLLING, int SURFACE>
        void registerPair(const std::string & name, LIGGGHTS::PairStyles::Factory & factory);

        template<int MODEL, int TANGENTIAL, int COHESION, int ROLLING, int SURFACE>
        void registerWall(const std::string & name, LIGGGHTS::Walls::Factory & factory)
        {
            typedef GranStyle<MODEL, TANGENTIAL, COHESION, ROLLING, SURFACE> Style;
            typedef ContactModel<Style> CModel;
            typedef Walls::Granular<CModel> Type;
            int64_t hashcode = Style::HASHCODE;
            
            factory.addStyle(name, hashcode, &create_wall_style_instance<Type>);
        }
    };

    template<int MODEL, int TANGENTIAL, int COHESION, int ROLLING, int SURFACE>
    inline void RegisterGranularStyles::registerPair(const std::string & name, LIGGGHTS::PairStyles::Factory & factory)
    {
        typedef GranStyle<MODEL, TANGENTIAL, COHESION, ROLLING, SURFACE> Style;
        typedef ContactModel<Style> CModel;
        typedef PairStyles::Granular<CModel> Type;
        int64_t hashcode = Style::HASHCODE;
        
        factory.addStyle(name, hashcode, &create_pair_style_instance<Type>);
    }
}

#endif
