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

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
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
