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
