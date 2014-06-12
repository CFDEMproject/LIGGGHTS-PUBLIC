/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
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
#ifndef SETTINGS_H_
#define SETTINGS_H_

#include "pointers.h"
#include <string>
#include <set>

using namespace LAMMPS_NS;
using namespace std;

template<typename T>
class ValuePropagator {
public:
  void registerTarget(T & target)
  {
    targets.insert(&target);
    target = currentValue;
  }

  void setValue(T value)
  {
    currentValue = value;

    for(typename set<T*>::iterator it = targets.begin(); it != targets.end(); ++it) {
      *(*it) = value;
    }
  }

private:
  T currentValue;
  set<T*> targets;
};

class Setting {
public:
  Setting(string name, int num_params) : name(name), num_params(num_params)
  {
  }
  virtual ~Setting(){}

  // subclasses implement this method and return the amount of arguments have
  // been consumed.
  virtual int parseArguments(char ** args) = 0;

  string name;
  int num_params;
  string error_message;
};

template<typename T>
class EnumSetting : public Setting
{
public:
  EnumSetting(string name) : Setting(name, 1)
  {
  }

  virtual ~EnumSetting(){}

  void addOption(string option, T value)
  {
    options[option] = value;
  }

  void setDefault(string option){
    current.setValue(options[option]);
  }

  void registerTarget(T & target) {
    current.registerTarget(target);
  }

  int parseArguments(char ** args) {
    if(name != args[0]) return 0; // argument not consumed
    string selected(args[1]);
    if(options.find(selected) != options.end()){
      current.setValue(options[selected]);
      return 2; // argument consumed
    } else {
      char msg[50];
      sprintf(msg, "while parsing '%s' argument: unknown option or wrong keyword order: '%s'", name.c_str(), args[1]);
      error_message = msg;
    }
    return -1; // error while parsing argument
  }

private:
  ValuePropagator<T> current;
  map<string, T> options;
};

class DoubleSetting : public Setting
{
public:
  DoubleSetting(string name) : Setting(name, 1)
  {
  }

  virtual ~DoubleSetting(){}

  void setDefault(double d){
    current.setValue(d);
  }

  void registerTarget(double & target) {
    current.registerTarget(target);
  }

  int parseArguments(char ** args) {
    if(name != args[0]) return 0;
    current.setValue(atof(args[1]));
    return 2;
  }

private:
  ValuePropagator<double> current;
};

class OnOffSetting : public EnumSetting<bool> {
public:
  OnOffSetting(string name, bool default_value) : EnumSetting<bool>(name)
  {
    addOption("off", default_value);
    addOption("on", !default_value);
    setDefault("off");
  }
};

class Settings : protected Pointers
{
  typedef map<string, OnOffSetting*> OnOffMap;
  typedef map<string, EnumSetting<int>*> EnumMap;
  typedef map<string, DoubleSetting*> DoubleMap;

public:
  std::string error_message;

  Settings(LAMMPS * lmp) : Pointers(lmp) {}
  ~Settings() {
    for(OnOffMap::iterator it = onOffSettings.begin(); it != onOffSettings.end(); ++it) {
      delete it->second;
    }

    for(EnumMap::iterator it = enumSettings.begin(); it != enumSettings.end(); ++it) {
      delete it->second;
    }

    for(DoubleMap::iterator it = doubleSettings.begin(); it != doubleSettings.end(); ++it) {
      delete it->second;
    }
  }

  void registerOnOff(string name, bool & variable, bool default_value = false)
  {
    if(onOffSettings.find(name) == onOffSettings.end()) {
      onOffSettings[name] = new OnOffSetting(name, default_value);
    }
    onOffSettings[name]->registerTarget(variable);
  }

  void registerDoubleSetting(string name, double & variable, double default_value = 0.0)
  {
    if(doubleSettings.find(name) == doubleSettings.end()) {
      doubleSettings[name] = new DoubleSetting(name);
    }
    doubleSettings[name]->setDefault(default_value);
    doubleSettings[name]->registerTarget(variable);
  }

  bool parseArguments(int nargs, char ** args) {
    bool found = false;
    int remaining = nargs;
    char ** remaining_args = args;
    int consumed = 0;

    while(remaining > 0) {
      found = false;

      for(OnOffMap::iterator it = onOffSettings.begin(); it != onOffSettings.end(); ++it) {
        consumed = it->second->parseArguments(remaining_args);
        if(consumed > 0) {
          found = true;
          break;
        } else if(consumed < 0) {
          error_message = it->second->error_message;
          return false;
        }
      }

      if(!found) {
        for(EnumMap::iterator it = enumSettings.begin(); it != enumSettings.end(); ++it) {
          consumed = it->second->parseArguments(remaining_args);
          if(consumed > 0) {
            found = true;
            break;
          } else if(consumed < 0) {
            error_message = it->second->error_message;
            return false;
          }
        }
      }

      if(!found) {
        for(DoubleMap::iterator it = doubleSettings.begin(); it != doubleSettings.end(); ++it) {
          consumed = it->second->parseArguments(remaining_args);
          if(consumed > 0) {
            found = true;
            break;
          } else if(consumed < 0) {
            error_message = it->second->error_message;
            return false;
          }
        }
      }

      if(found) {
        remaining -= consumed;
        remaining_args = &remaining_args[consumed];
      } else {
        char msg[30];
        sprintf(msg, "Unknown argument: '%s'", remaining_args[0]);
        error_message = msg;
        return false;
      }
    }
    return true;
  }

private:
  // map for each type of setting, since RTTI is turned off by default (required by dynamic_cast)
  OnOffMap onOffSettings;
  EnumMap enumSettings;
  DoubleMap doubleSettings;
};

#endif /* SETTINGS_H_ */
