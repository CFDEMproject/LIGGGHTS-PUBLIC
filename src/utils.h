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
