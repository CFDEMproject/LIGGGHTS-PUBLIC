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
namespace LIGGGHTS {

namespace ContactModels
{
  // This file declares global constant identifiers which are used by contact
  // models for template specialization. Please note that each model constant
  // in a group must be unique.

// surface models
#define SURFACE_MODEL(identifier,str,constant) \
static const int identifier = constant;
#include "style_surface_model.h"
#undef SURFACE_MODEL

  // normal models
#define NORMAL_MODEL(identifier,str,constant) \
  static const int identifier = constant;
#include "style_normal_model.h"
#undef NORMAL_MODEL

  // tangential models
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
