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

#ifdef FIX_CLASS

FixStyle(property/atom,FixPropertyAtom)

#else

#ifndef LMP_FIX_PROPERTY_ATOM_H
#define LMP_FIX_PROPERTY_ATOM_H

#include "fix.h"
#include "input.h"

namespace LAMMPS_NS {

enum
{
   FIXPROPERTY_ATOM_SCALAR = 0,
   FIXPROPERTY_ATOM_VECTOR = 1
};

class FixPropertyAtom : public Fix {
 friend class Set;
 friend class FixPropertyAtomUpdateFix;
 friend class FixPropertyAtomRandom;
 public:
  FixPropertyAtom(class LAMMPS *, int, char **,bool parse = true);
  ~FixPropertyAtom();
  virtual int setmask();

  void do_forward_comm();
  void do_reverse_comm();

  Fix* check_fix(const char *varname,const char *svmstyle,int len1,int len2,const char *caller,bool errflag);

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int,int);
  void pre_set_arrays();
  virtual void set_arrays(int);

  void set_all(double value);

  void write_restart(FILE *);
  virtual void restart(char *);

  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();
  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double compute_vector(int n);

  virtual void mark_tracers(int ilo, int ihi) { UNUSED(ilo); UNUSED(ihi); }

 protected:
  void parse_args(int narg, char **arg);

 private:
  char *variablename;   // name of the variable (used for identification by other fixes)
  int data_style;            // 0 if a scalar is registered, 1 if vector
  int commGhost;        // 1 if communicated to ghost particles (via pack_comm/unpack_comm), 0 if not
  int commGhostRev;     // 1 if rev communicated from ghost particles (via pack_comm_rev/unpack_comm_rev), 0 if not
  int nvalues;
  double *defaultvalues; // default values at particle creation

  // in case of initialization from property - name of property
  char *propertyname;
  double *property;
}; //end class

}
#endif
#endif
