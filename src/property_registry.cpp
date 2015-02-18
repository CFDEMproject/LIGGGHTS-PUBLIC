/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2013 JKU Linz
   Copyright 2013-     DCS Computing GmbH, Linz

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

#include "pointers.h"
#include "lammps.h"
#include "fix_property_global.h"
#include <map>
#include <set>
#include <string>
#include "properties.h"
#include "error.h"
#include "modify.h"
#include "property_registry.h"

using namespace std;
using namespace LAMMPS_NS;

PropertyRegistry::PropertyRegistry(LAMMPS* lmp) : Pointers(lmp), properties(lmp)
{
}

PropertyRegistry::~PropertyRegistry()
{
  init();
}

int PropertyRegistry::max_type()
{
  return properties.max_type();
}

LAMMPS * PropertyRegistry::getLAMMPS()
{
  return lmp;
}

FixPropertyGlobal* PropertyRegistry::getGlobalProperty(const char *varname, const char *style, const char *svmstyle, int len1, int len2, const char *caller)
{
  return static_cast<FixPropertyGlobal*>(modify->find_fix_property(varname, style, svmstyle, len1, len2, caller));
}

ScalarProperty * PropertyRegistry::getScalarProperty(string varname,const char *caller)
{
  if(scalars.find(varname) == scalars.end()) {
    if(scalar_creators.find(varname) != scalar_creators.end()) {
      scalars[varname] = (*scalar_creators[varname])(*this, caller, use_sanity_checks[varname]);
    } else {
      error->message(FLERR, "unknown scalar property");
    }
  }
  return scalars[varname];
}

VectorProperty * PropertyRegistry::getVectorProperty(string varname,const char *caller)
{
  if(vectors.find(varname) == vectors.end()) {
    if(vector_creators.find(varname) != vector_creators.end()) {
      vectors[varname] = (*vector_creators[varname])(*this, caller, use_sanity_checks[varname]);
    } else {
      error->message(FLERR, "unknown vector property");
    }
  }
  return vectors[varname];
}

MatrixProperty * PropertyRegistry::getMatrixProperty(string varname,const char *caller)
{
  if(matrices.find(varname) == matrices.end()) {
    if(matrix_creators.find(varname) != matrix_creators.end()) {
      matrices[varname] = (*matrix_creators[varname])(*this, caller, use_sanity_checks[varname]);
    } else {
      error->message(FLERR, "unknown matrix property");
    }
  }
  return matrices[varname];
}

void PropertyRegistry::registerProperty(string varname, ScalarPropertyCreator creator, bool sanity_checks)
{
  if(scalar_creators.find(varname) == scalar_creators.end()) {
    scalar_creators[varname] = creator;
    use_sanity_checks[varname] = sanity_checks;
  } else if(scalar_creators[varname] != creator) {
    error->message(FLERR, "property with the same name, but different implementation registered");
  }
}

void PropertyRegistry::registerProperty(string varname, VectorPropertyCreator creator, bool sanity_checks)
{
  if(vector_creators.find(varname) == vector_creators.end()) {
    vector_creators[varname] = creator;
    use_sanity_checks[varname] = sanity_checks;
  } else if(vector_creators[varname] != creator) {
    error->message(FLERR, "property with the same name, but different implementation registered");
  }
}

void PropertyRegistry::registerProperty(string varname, MatrixPropertyCreator creator, bool sanity_checks)
{
  if(matrix_creators.find(varname) == matrix_creators.end()) {
    matrix_creators[varname] = creator;
    use_sanity_checks[varname] = sanity_checks;
  } else if(matrix_creators[varname] != creator) {
    error->message(FLERR, "property with the same name, but different implementation registered");
  }
}

void PropertyRegistry::connect(string varname, double ** & variable, const char *caller)
{
  if(matrices.find(varname) == matrices.end()) {
    if(matrix_creators.find(varname) != matrix_creators.end()) {
      matrices[varname] = (*matrix_creators[varname])(*this, caller, use_sanity_checks[varname]);
    } else {
      // ERROR unknown property
    }
  }
  matrices[varname]->connect(variable);
}

void PropertyRegistry::connect(string varname, double * & variable, const char *caller)
{
  if(vectors.find(varname) == vectors.end()) {
    if(vector_creators.find(varname) != vector_creators.end()) {
      vectors[varname] = (*vector_creators[varname])(*this, caller, use_sanity_checks[varname]);
    } else {
      // ERROR unknown property
    }
  }
  vectors[varname]->connect(variable);
}

void PropertyRegistry::connect(string varname, double & variable, const char *caller)
{
  if(scalars.find(varname) == scalars.end()) {
    if(scalar_creators.find(varname) != scalar_creators.end()) {
      scalars[varname] = (*scalar_creators[varname])(*this, caller, use_sanity_checks[varname]);
    } else {
      // ERROR unknown property
    }
  }
  scalars[varname]->connect(variable);
}

void PropertyRegistry::init()
{
  for(std::map<string,ScalarProperty*>::iterator it = scalars.begin(); it != scalars.end(); ++it) {
      delete it->second;
  }
  for(std::map<string,VectorProperty*>::iterator it = vectors.begin(); it != vectors.end(); ++it) {
      delete it->second;
  }
  for(std::map<string,MatrixProperty*>::iterator it = matrices.begin(); it != matrices.end(); ++it) {
      delete it->second;
  }
  scalars.clear();
  vectors.clear();
  matrices.clear();
}

void PropertyRegistry::print_all(FILE * out)
{
  for(std::map<string,ScalarProperty*>::iterator it = scalars.begin(); it != scalars.end(); ++it) {
      fprintf(out, " %s = ", it->first.c_str());
      it->second->print_value(out);
      fprintf(out, "\n");
  }
  for(std::map<string,VectorProperty*>::iterator it = vectors.begin(); it != vectors.end(); ++it) {
      fprintf(out, " %s = ", it->first.c_str());
      it->second->print_value(out);
      fprintf(out, "\n");
  }
  for(std::map<string,MatrixProperty*>::iterator it = matrices.begin(); it != matrices.end(); ++it) {
      fprintf(out, " %s = ", it->first.c_str());
      it->second->print_value(out);
      fprintf(out, "\n");
  }
}
