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

    Richard Berger (JKU Linz)

    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef UTILS_H
#define UTILS_H

#include "mpi.h"
#include "lmptype.h"
#include "lammps.h"
#include <string>
#include <map>
#include <iostream>
#include "contact_interface.h"

namespace LIGGGHTS {
using namespace LAMMPS_NS;

namespace Utils {
  template<typename Interface>
  class AbstractFactory {
    typedef typename Interface::ParentType ParentType;
    typedef Interface * (*Creator)(class LAMMPS * lmp, ParentType* parent);
    typedef int64_t (*VariantSelector)(int & argc, char ** & argv);
    typedef std::map<std::pair<std::string, int>, Creator> StyleTable;
    typedef std::map<std::string, VariantSelector> VariantSelectorTable;
    StyleTable styleTable;
    VariantSelectorTable variantSelectorTable;
    AbstractFactory(const AbstractFactory &){}
    AbstractFactory& operator=(const AbstractFactory&){}

  protected:
    AbstractFactory() {}

  public:
    Interface * create(const std::string & name, int64_t variant, class LAMMPS * lmp, ParentType* parent) {
      std::pair<std::string, int> key(name, variant);
      if(styleTable.find(key) != styleTable.end()) {
        return styleTable[key](lmp, parent);
      }
      return NULL;
    }

    int64_t selectVariant(const std::string & name, int & argc, char ** & argv) {
      if(variantSelectorTable.find(name) != variantSelectorTable.end()) {
        return variantSelectorTable[name](argc, argv);
      }
      return 0;
    }

    void addStyle(const std::string & name, int variant, Creator create) {
      std::pair<std::string, int> key(name, variant);
      if(styleTable.find(key) != styleTable.end()){
        std::cerr << "WARNING! Style collision detected! Duplicate entry (" << key.first << ", " << key.second <<  ") in style table." << std::endl;
      }
      styleTable[key] = create;
    }

    void addVariantSelector(const std::string & name, VariantSelector selector) {
      if(variantSelectorTable.find(name) != variantSelectorTable.end()){
        std::cerr << "WARNING! VariantSelector collision detected! Duplicate entry '" << name << "' in variant selector table." << std::endl;
      }
      variantSelectorTable[name] = selector;
    }
  };
}

}

#endif // UTILS_H
