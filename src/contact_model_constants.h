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
    Arno Mayrhofer (DCS Computing GmbH, Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef CONTACT_MODEL_CONSTANTS_H_
#define CONTACT_MODEL_CONSTANTS_H_

namespace LIGGGHTS
{

namespace ContactModels
{
    // This file declares global constant identifiers which are used by contact
    // models for template specialization. Please note that each model constant
    // in a group must be unique.

    static const int CM_REGISTER_SETTINGS      = 1 << 0;
    static const int CM_CONNECT_TO_PROPERTIES  = 1 << 1;
    static const int CM_BEGIN_PASS             = 1 << 2;
    static const int CM_END_PASS               = 1 << 3;
    static const int CM_SURFACES_INTERSECT     = 1 << 4;
    static const int CM_SURFACES_CLOSE         = 1 << 5;

    static const int CONTACT_NORMAL_MODEL      = 1 << 0;
    static const int CONTACT_COHESION_MODEL    = 1 << 1;
    static const int CONTACT_TANGENTIAL_MODEL  = 1 << 2;
    static const int CONTACT_ROLLING_MODEL     = 1 << 3;
    static const int CONTACT_SURFACE_MODEL     = 1 << 4;
    static const int CONTACT_FIX               = 1 << 31;

    // surface models
    #define SURFACE_MODEL(identifier,str,constant) \
    static const int identifier = constant;
    #include "style_surface_model.h"
    #undef SURFACE_MODEL

    // normal models
    static const int NORMAL_OFF = 9;
    #define NORMAL_MODEL(identifier,str,constant) \
    static const int identifier = constant;
    #include "style_normal_model.h"
    #undef NORMAL_MODEL

    // tangential models
    static const int TANGENTIAL_OFF = 0;
    #define TANGENTIAL_MODEL(identifier,str,constant) \
    static const int identifier = constant;
    #include "style_tangential_model.h"
    #undef TANGENTIAL_MODEL

    // cohesion models
    static const int COHESION_OFF = 0;
    #define COHESION_MODEL(identifier,str,constant) \
    static const int identifier = constant;
    #include "style_cohesion_model.h"
    #undef COHESION_MODEL

    // rolling models
    static const int ROLLING_OFF = 0;
    #define ROLLING_MODEL(identifier,str,constant) \
    static const int identifier = constant;
    #include "style_rolling_model.h"
    #undef ROLLING_MODEL
}

}

#endif
