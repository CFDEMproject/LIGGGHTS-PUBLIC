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
    This file is from LAMMPS, but has been modified. Copyright for
    modification:

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz

    Copyright of original file:
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

#ifndef LMP_FORCE_H
#define LMP_FORCE_H

#include "pointers.h"
#include "property_registry.h"
#include <map>
#include <string>
#include <vector>
#include <algorithm>

namespace LAMMPS_NS {

struct Custom_contact_models
{
  std::string custom_surface_model;
  std::string custom_normal_model;
  std::string custom_tangential_model;
  std::string custom_cohesion_model;
  std::string custom_rolling_model;
};

class Force : protected Pointers {
 friend class Coarsegraining;
 friend class StiffnessScaling;

 public:
  double boltz;                      // Boltzmann constant (eng/degree-K)
  double hplanck;                    // Planck's constant (energy-time)
  double mvv2e;                      // conversion of mv^2 to energy
  double ftm2v;                      // conversion of ft/m to velocity
  double mv2d;                       // conversion of mass/volume to density
  double nktv2p;                     // conversion of NkT/V to pressure
  double qqr2e;                      // conversion of q^2/r to energy
  double qe2f;                       // conversion of qE to force
  double vxmu2f;                     // conversion of vx dynamic-visc to force
  double xxt2kmu;                    // conversion of xx/t to kinematic-visc
  double dielectric;                 // dielectric constant
  double qqrd2e;                     // q^2/r to energy w/ dielectric constant
  double e_mass;                     // electron mass
  double hhmrr2e;                    // conversion of (hbar)^2/(mr^2) to energy
  double mvh2r;                      // conversion of mv/hbar to distance
                                     // hbar = h/(2*pi)
  double angstrom;                   // 1 angstrom in native units
  double femtosecond;                // 1 femtosecond in native units
  double qelectron;                  // 1 electron charge abs() in native units

  int newton,newton_pair,newton_bond;   // Newton's 3rd law settings

  class Pair *pair;
  char *pair_style;

  typedef Pair *(*PairCreator)(LAMMPS *);
  std::map<std::string,PairCreator> *pair_map;
  Custom_contact_models custom_contact_models;

  class Bond *bond;
  char *bond_style;

  class Angle *angle;
  char *angle_style;

  class Dihedral *dihedral;
  char *dihedral_style;

  class Improper *improper;
  char *improper_style;

  class KSpace *kspace;
  char *kspace_style;
                             // index [0] is not used in these arrays
  double special_lj[4];      // 1-2, 1-3, 1-4 prefactors for LJ
  double special_coul[4];    // 1-2, 1-3, 1-4 prefactors for Coulombics
  int special_angle;         // 0 if defined angles are ignored
                             // 1 if only weight 1,3 atoms if in an angle
  int special_dihedral;      // 0 if defined dihedrals are ignored
                             // 1 if only weight 1,4 atoms if in a dihedral
  int special_extra;         // extra space for added bonds

  Force(class LAMMPS *);
  ~Force();
  void init();

  void create_pair(const char *, const char *suffix = NULL);
  void create_pair_from_restart(FILE* fp,const char *, const char *suffix = NULL);
  class Pair *new_pair(const char *, const char *, int &);
  class Pair *new_pair_from_restart(FILE* fp, const char *, const char *, int &);

  class Pair *pair_match(const char *, int);

  void create_bond(const char *, const char *suffix = NULL);
  class Bond *new_bond(const char *, const char *, int &);
  class Bond *bond_match(const char *);

  void create_angle(const char *, const char *suffix = NULL);
  class Angle *new_angle(const char *, const char *, int &);

  void create_dihedral(const char *, const char *suffix = NULL);
  class Dihedral *new_dihedral(const char *, const char *, int &);

  void create_improper(const char *, const char *suffix = NULL);
  class Improper *new_improper(const char *, const char *, int &);

  void create_kspace(int, char **, const char *suffix = NULL);
  class KSpace *new_kspace(int, char **, const char *, int &);
  class KSpace *kspace_match(const char *, int);

  void set_special(int, char **);
  void bounds(char *, int, int &, int &, int nmin=1);
  double numeric(const char *, const int, const char *const);
  int inumeric(const char *, const int, const char *const);
  bigint memory_usage();

  bool setCG(double cg)
  {
      bool useTypeSpecific = false;
      if(coarsegraining_>1.0)
      {
        coarsegraining_ = std::max(coarsegraining_,cg); //set maximum we use type-specific CG
        useTypeSpecific = true;
      }
      else
        coarsegraining_ = cg;

      coarsegrainingTypeBased_.push_back(cg);

      return useTypeSpecific;
  }

  void reportCG()
  {
    printf("Force: coarsegrainingfactor: %g.\n", coarsegraining_);
    for(unsigned int it=0;it<coarsegrainingTypeBased_.size();it++)
        printf("Force: type-specific coarsegrainingfactor(%d): %g.\n", it,coarsegrainingTypeBased_[it]);
  }

  inline double cg(int typeID)
  {
    if(typeID<=int(coarsegrainingTypeBased_.size()))
        return coarsegrainingTypeBased_[typeID-1];
    else
        return coarsegraining_;
  }

  inline double cg_max()
  {
    if (coarsegrainingTypeBased_.size() > 0) {
      const double max_cg_type = *(std::max_element(coarsegrainingTypeBased_.begin(),coarsegrainingTypeBased_.end())  );
      if(max_cg_type > coarsegraining_)
        return max_cg_type;
    }

    return coarsegraining_;
  }

  //inline double cg() 
  //{ return coarsegraining_; }

  inline bool cg_active() 
  { return (coarsegraining_ > 1. || coarsegrainingTypeBased_.size() > 0); }

  inline bool error_cg() 
  { return error_coarsegraining_; }

  inline bool warn_cg() 
  { return warn_coarsegraining_; }

  PropertyRegistry registry;

  void set_custom_surface_model(std::string param)
  { custom_contact_models.custom_surface_model = param; }
  std::string get_custom_surface_model()
  { return custom_contact_models.custom_surface_model; }

  void set_custom_normal_model(std::string param)
  { custom_contact_models.custom_normal_model = param; }
  std::string get_custom_normal_model()
  { return custom_contact_models.custom_normal_model; }

  void set_custom_tangential_model(std::string param)
  { custom_contact_models.custom_tangential_model = param; }
  std::string get_custom_tangential_model()
  { return custom_contact_models.custom_tangential_model; }

  void set_custom_cohesion_model(std::string param)
  { custom_contact_models.custom_cohesion_model = param; }
  std::string get_custom_cohesion_model()
  { return custom_contact_models.custom_cohesion_model; }

  void set_custom_rolling_model(std::string param)
  { custom_contact_models.custom_rolling_model = param; }
  std::string get_custom_rolling_model()
  { return custom_contact_models.custom_rolling_model; }

 private:
  template <typename T> static Pair *pair_creator(LAMMPS *);

  double    coarsegraining_;
  std::vector<double> coarsegrainingTypeBased_;
  bool      error_coarsegraining_;
  bool      warn_coarsegraining_;
};

}

#endif

/* ERROR/WARNING messages:

E: Invalid pair style

The choice of pair style is unknown.

E: Invalid bond style

The choice of bond style is unknown.

E: Invalid angle style

The choice of angle style is unknown.

E: Invalid dihedral style

The choice of dihedral style is unknown.

E: Invalid improper style

The choice of improper style is unknown.

E: Invalid kspace style

The choice of kspace style is unknown.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Numeric index is out of bounds

A command with an argument that specifies an integer or range of
integers is using a value that is less than 1 or greater than the
maximum allowed limit.

U: Expected floating point parameter in input script or data file

The quantity being read is an integer on non-numeric value.

U: Expected integer parameter in input script or data file

The quantity being read is a floating point or non-numeric value.

*/
