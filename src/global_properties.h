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

#ifndef GLOBAL_PROPERTIES_H_
#define GLOBAL_PROPERTIES_H_

#include "property_registry.h"

namespace MODEL_PARAMS
{

  /* -----------------------------------------------------------------------
   * Utility functions
   * ----------------------------------------------------------------------- */

  ScalarProperty* createScalarProperty(PropertyRegistry & registry, const char* name, const char * caller);

  MatrixProperty* createPerTypePairProperty(PropertyRegistry & registry, const char * name, const char * caller);

  /* -----------------------------------------------------------------------
   * Property Creators
   * ----------------------------------------------------------------------- */

  ScalarProperty* createCharacteristicVelocity(PropertyRegistry & registry, const char * caller, bool sanity_checks);
  MatrixProperty* createCohesionEnergyDensity(PropertyRegistry & registry, const char * caller, bool sanity_checks);
  VectorProperty* createYoungsModulus(PropertyRegistry & registry, const char * caller, bool sanity_checks);
  VectorProperty* createPoissonsRatio(PropertyRegistry & registry, const char * caller, bool sanity_checks);
  MatrixProperty* createYeff(PropertyRegistry & registry, const char * caller, bool sanity_checks);
  MatrixProperty* createGeff(PropertyRegistry & registry, const char * caller, bool sanity_checks);
  MatrixProperty* createCoeffRest(PropertyRegistry & registry, const char * caller, bool sanity_checks);
  MatrixProperty* createCoeffRestLog(PropertyRegistry & registry, const char * caller, bool sanity_checks);
  MatrixProperty* createBetaEff(PropertyRegistry & registry, const char * caller, bool sanity_checks);

  MatrixProperty* createCoeffFrict(PropertyRegistry & registry, const char * caller, bool sanity_checks);
  MatrixProperty* createCoeffRollFrict(PropertyRegistry & registry, const char * caller, bool sanity_checks);
  MatrixProperty* createCoeffRollVisc(PropertyRegistry & registry, const char * caller, bool sanity_checks);
  MatrixProperty* createCoeffMu(PropertyRegistry & registry, const char * caller, bool sanity_checks);
  MatrixProperty* createCoeffRestMax(PropertyRegistry & registry, const char * caller, bool sanity_checks);
  MatrixProperty* createCoeffStc(PropertyRegistry & registry, const char * caller, bool sanity_checks);

  ScalarProperty* createLiquidVolume(PropertyRegistry & registry, const char * caller, bool sanity_checks);
  ScalarProperty* createSurfaceTension(PropertyRegistry & registry, const char * caller, bool sanity_checks);
  ScalarProperty* createSwitchModel(PropertyRegistry & registry, const char * caller, bool sanity_checks);
  ScalarProperty* createHistoryIndex(PropertyRegistry & registry, const char * caller, bool sanity_checks);
  VectorProperty* createContactAngle(PropertyRegistry & registry, const char * caller, bool sanity_checks);
  MatrixProperty* createKn(PropertyRegistry & registry, const char * caller, bool sanity_checks);
  MatrixProperty* createKt(PropertyRegistry & registry, const char * caller, bool sanity_checks);
  MatrixProperty* createGamman(PropertyRegistry & registry, const char * caller, bool sanity_checks);
  MatrixProperty* createGammat(PropertyRegistry & registry, const char * caller, bool sanity_checks);
  MatrixProperty* createGammanAbs(PropertyRegistry & registry, const char * caller, bool sanity_checks);
  MatrixProperty* createGammatAbs(PropertyRegistry & registry, const char * caller, bool sanity_checks);
  MatrixProperty* createCoeffMaxElasticStiffness(PropertyRegistry & registry, const char * caller, bool sanity_checks);
  MatrixProperty* createCoeffAdhesionStiffness(PropertyRegistry & registry, const char * caller, bool sanity_checks);
  MatrixProperty* createCoeffPlasticityDepth(PropertyRegistry & registry, const char * caller, bool sanity_checks);
}

#endif /* GLOBAL_PROPERTIES_H_ */
