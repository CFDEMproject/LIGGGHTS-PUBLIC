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

    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Richard Berger (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
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
      error->message(FLERR, "unknown matrix property");
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
      error->message(FLERR, "unknown vector property");
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
      error->message(FLERR, "unknown scalar property");
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
