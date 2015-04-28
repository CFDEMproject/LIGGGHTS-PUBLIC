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
    Andreas Aigner (JKU Linz)

    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(sph/density/corr,FixSphDensityCorr)

#else

#ifndef LMP_FIX_SPH_DENSITY_CORR_H
#define LMP_FIX_SPH_DENSITY_CORR_H

#include "fix_sph.h"

namespace LAMMPS_NS {

  enum {CORR_SHEPARD,CORR_MLS};

class FixSphDensityCorr : public FixSph {
 public:
  FixSphDensityCorr(class LAMMPS *, int, char **);
  ~FixSphDensityCorr();
  void pre_delete(bool unfixflag);
  virtual int setmask();
  void updatePtrs();
  void post_create();
  virtual void init();
  virtual void pre_force(int vflag);

 private:
  template <int> void pre_force_eval();

  class FixPropertyAtom* fix_quantity;
  char *quantity_name;

  double quantity_0;
  double *quantity;

  int corrStyle;
  int every;
  int ago;

};

}

#endif
#endif
