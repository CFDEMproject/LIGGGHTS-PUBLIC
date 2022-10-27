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

#include <mpi.h>
#include "lmptype.h"
#include "force.h"
#include <string>
#include <map>
#include <iostream>
#include "contact_interface.h"
#include "contact_model_constants.h"
#include <sstream>

namespace LIGGGHTS {
using namespace LAMMPS_NS;

namespace Utils {

  inline int64_t generate_gran_hashcode(int model, int tangential, int cohesion, int rolling, int surface)
  {
    return (((int64_t)model)           ) |
           (((int64_t)tangential) <<  6) |
           (((int64_t)cohesion)   << 12) |
           (((int64_t)rolling)    << 18) |
           (((int64_t)surface)    << 24) ;
  }

  inline std::string int_to_string(int a)
  {
    // see https://www.cfdem.com/forums/error-non-const-lvalue-reference-type-basicostringstream-cannot-bind-temporary-type
    // return static_cast< std::ostringstream & >(( std::ostringstream() << std::dec << a ) ).str();
    std::ostringstream ss;
    ss << std::dec << a;
    return ss.str();
  }

  inline std::string double_to_string(double dbl)
  {
    std::ostringstream strs;
    strs << dbl;
    std::string str = strs.str();
    return str;
  }

  template <typename T>
  inline T* ptr_reduce(T** &t)
  { return &(t[0][0]); }

  template <typename T>
  inline T* ptr_reduce(T* &t)
  { return t; }

  template <typename T>
  inline T* ptr_reduce(T &t)
  { return &t; }

  template<typename Interface>
  class AbstractFactory {
    typedef typename Interface::ParentType ParentType;
    typedef Interface * (*Creator)(class LAMMPS * lmp, ParentType* parent, int64_t hash);
    typedef int64_t (*VariantSelector)(int & argc, char ** & argv, Custom_contact_models ccm);
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
        return styleTable[key](lmp, parent, variant);
      }
      int64_t default_variant = generate_gran_hashcode(ContactModels::NORMAL_OFF, ContactModels::TANGENTIAL_OFF, ContactModels::COHESION_OFF, ContactModels::ROLLING_OFF, 0);
      std::pair<std::string, int> default_key(name, default_variant);
      if(styleTable.find(default_key) != styleTable.end()) {
        return styleTable[default_key](lmp, parent, variant);
      }
      return NULL;
    }

    int64_t selectVariant(const std::string & name, int & argc, char ** & argv,Custom_contact_models ccm) {
      if(variantSelectorTable.find(name) != variantSelectorTable.end()) {
        return variantSelectorTable[name](argc, argv,ccm);
      }
      return 0;
    }

    void addStyle(const std::string & name, int variant, Creator create) {
      std::pair<std::string, int> key(name, variant);
      if(styleTable.find(key) != styleTable.end()){
        std::cerr << "WARNING! Style collision detected! Duplicate entry (" << key.first << ", " << key.second << ") in style table." << std::endl;
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
