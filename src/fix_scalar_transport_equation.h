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

FixStyle(transportequation/scalar,FixScalarTransportEquation)

#else

#ifndef LMP_FIX_SCALAR_TRANSPORT_EQUATION_H
#define LMP_FIX_SCALAR_TRANSPORT_EQUATION_H

#include "fix.h"

namespace LAMMPS_NS {

class FixScalarTransportEquation : public Fix {
 public:
  FixScalarTransportEquation(class LAMMPS *, int, char **);
  ~FixScalarTransportEquation();

  virtual int setmask();
  virtual void post_create();
  virtual void pre_delete(bool unfixflag);
  virtual void init();
  virtual int modify_param(int narg, char **arg);
  virtual void updatePtrs();
  virtual void initial_integrate_respa(int,int,int);
  virtual void initial_integrate(int);
  virtual void pre_force(int vflag);
  virtual void final_integrate();
  virtual double compute_scalar();
  bool match_equation_id(const char*);

  double *get_capacity();

 protected:

  int nlevels_respa;

  char *equation_id;

  class FixPropertyAtom* fix_quantity;
  char *quantity_name;
  class FixPropertyAtom* fix_flux;
  char *flux_name;
  class FixPropertyAtom* fix_source;
  char *source_name;

  //storage capacity - would be thermal capacity for heat conduction
  int capacity_flag;
  class FixPropertyGlobal* fix_capacity;
  double *capacity;
  char *capacity_name;

  double quantity_0;  
  double *quantity;   
  double *flux;       
  double *source;     

  // flag if integrate quantity or not
  bool int_flag;

  int  nevery_; //integrate only this many time steps (to avoid round-off issues)
  bool performedIntegrationLastStep_;
};

}

#endif
#endif
