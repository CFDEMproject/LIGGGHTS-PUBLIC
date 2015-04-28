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

FixStyle(property/global,FixPropertyGlobal)
FixStyle(property/atomtype,FixPropertyGlobal)
FixStyle(property/atomtypepair,FixPropertyGlobal)

#else

#ifndef LMP_FIX_PROPERTYGLOBAL_H
#define LMP_FIX_PROPERTYGLOBAL_H
#include "fix.h"
#include "input.h"

namespace LAMMPS_NS {

enum
{
        FIXPROPERTY_GLOBAL_SCALAR = 0,
        FIXPROPERTY_GLOBAL_VECTOR = 1,
        FIXPROPERTY_GLOBAL_MATRIX = 2
};

class FixPropertyGlobal : public Fix {
 friend class Modify;
 friend class MechParamGran;
 friend class CfdDatacouplingFile;

 public:
  FixPropertyGlobal(class LAMMPS *, int, char **);
  ~FixPropertyGlobal();
  int setmask();
  void init();
  void pre_delete(bool unfixflag);

  Fix* check_fix(const char *varname,const char *svmstyle,int len1,int len2,const char *caller,bool errflag);

  double memory_usage();
  double compute_scalar();
  double compute_vector(int);
  double compute_vector_modified(int);
  double compute_array(int,int);
  double compute_array_modified(int,int);
  void vector_modify(int,double);
  void array_modify(int,int,double);
  void new_array(int l1,int l2);
  int modify_param(int narg, char **arg);

  //bool checkCorrectness(int,char*,int,int);

  const double* get_values() {return values;}
  const double* get_values_modified() {return values_recomputed;}
  double const* const* get_array() {return array;}
  double const* const* get_array_modified() {return array_recomputed;}

  void grow(int,int);

  void write();

  char *variablename;        // name of the variable (used for identification by other fixes)
  int data_style;            // 0 if a scalar is registered, 1 if vector, 2 if 2d array (matrix)
  int nvalues;
  int nvalues_new_array;

  bool is_symmetric;         // flag if values must be symmetric (only applicable for matrix)
  bool is_atomtype_bound;    // flag if # values is bound to # atom types
  double *values;            // original values to be stored in this fix
  double *values_recomputed; // values to be stored in this fix, recomputed by eg another fix
  double **array;
  double **array_recomputed;

  char *filename;
  char *grpname;
  int me;

}; //end class

}
#endif
#endif
