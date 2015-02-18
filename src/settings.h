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
#ifndef SETTINGS_H_
#define SETTINGS_H_

#include "pointers.h"
#include <string>
#include <set>
#include <map>
#include <sstream>

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

  T getValue() const {
    return currentValue;
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

  virtual void print_value(FILE * out) = 0;

  string name;
  int num_params;
  string error_message;
};

template<typename T>
class EnumSetting : public Setting
{
public:
  typedef T value_type;

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
      stringstream ss;
      ss << "while parsing '" << name << "' argument: ";
      ss << "unknown option or wrong keyword order: '" << args[1] << "'";
      error_message = ss.str();
    }
    return -1; // error while parsing argument
  }

  virtual void print_value(FILE* out) {
    T value = current.getValue();
    for(typename map<string, T>::iterator it = options.begin(); it != options.end(); ++it) {
      if(it->second == value) {
        fprintf(out, "%s", it->first.c_str());
        return;
      }
    }
    fprintf(out, "BAD_VALUE");
  }

private:
  ValuePropagator<T> current;
  map<string, T> options;
};

class DoubleSetting : public Setting
{
public:
  typedef double value_type;

  DoubleSetting(string name, double default_value) : Setting(name, 1)
  {
    setDefault(default_value);
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

  void print_value(FILE* out) {
    fprintf(out, "%g", current.getValue());
  }

private:
  ValuePropagator<double> current;
};

class OnOffSetting : public EnumSetting<bool> {
public:
  OnOffSetting(string name, bool default_value) : EnumSetting<bool>(name)
  {
    addOption("off", false);
    addOption("on", true);
    if(default_value) setDefault("on");
    else setDefault("off");
  }

  virtual ~OnOffSetting(){}
};

class YesNoSetting : public EnumSetting<bool> {
public:
  YesNoSetting(string name, bool default_value) : EnumSetting<bool>(name)
  {
    addOption("no", false);
    addOption("yes", true);
    if(default_value) setDefault("yes");
    else setDefault("no");
  }
	virtual ~YesNoSetting(){}
};

class Settings : protected Pointers
{
  typedef map<string, Setting*> SettingMap;
  SettingMap settings;

  template<typename SettingType>
  void registerSetting(string name, typename SettingType::value_type & variable, typename SettingType::value_type default_value) {
    if(settings.find(name) == settings.end()) {
      settings[name] = new SettingType(name, default_value);
    }
    SettingType * setting = dynamic_cast<SettingType*>(settings[name]);
    if(setting) setting->registerTarget(variable);
  }

public:
  std::string error_message;

  Settings(LAMMPS * lmp) : Pointers(lmp) {}
  ~Settings() {
    for(SettingMap::iterator it = settings.begin(); it != settings.end(); ++it) {
      delete it->second;
    }
  }

  void registerOnOff(string name, bool & variable, bool default_value = false)
  {
    registerSetting<OnOffSetting>(name, variable, default_value);
  }

  void registerYesNo(string name, bool & variable, bool default_value = false) {
    registerSetting<YesNoSetting>(name, variable, default_value);
  }

  void registerDoubleSetting(string name, double & variable, double default_value = 0.0)
  {
    registerSetting<DoubleSetting>(name, variable, default_value);
  }

  bool parseArguments(int nargs, char ** args) {
    bool found = false;
    int remaining = nargs;
    char ** remaining_args = args;
    int consumed = 0;

    while(remaining > 0) {
      found = false;

      for(SettingMap::iterator it = settings.begin(); it != settings.end(); ++it) {
        consumed = it->second->parseArguments(remaining_args);
        if(consumed > 0) {
          found = true;
          break;
        } else if(consumed < 0) {
          error_message = it->second->error_message;
          return false;
        }
      }

      if(found) {
        remaining -= consumed;
        remaining_args = &remaining_args[consumed];
      } else {
        stringstream ss;
        ss << "Unknown argument or wrong keyword order: '" << remaining_args[0] << "'";
        error_message = ss.str();
        return false;
      }
    }
    return true;
  }

  void print_all(FILE * out)
  {
    for(SettingMap::iterator it = settings.begin(); it != settings.end(); ++it) {
      fprintf(out, " %s = ", it->first.c_str());
      it->second->print_value(out);
      fprintf(out, "\n");
    }
  }
};

#endif /* SETTINGS_H_ */
