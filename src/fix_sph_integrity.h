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
    Andreas Eitzlmayr (TU Graz)

    Copyright 2009-2012 JKU Linz
    Copyright 2013-     TU Graz
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(sph/integrity,FixSphIntegrity)

#else

#ifndef LMP_FIX_SPH_INTEGRITY_H
#define LMP_FIX_SPH_INTEGRITY_H

#include "fix_sph.h"

namespace LAMMPS_NS {

class FixSphIntegrity : public FixSph {
 public:
  FixSphIntegrity(class LAMMPS *, int, char **);
  ~FixSphIntegrity();
  virtual int setmask();
  void post_create();
  void pre_delete(bool unfixflag);
  virtual void init();
  virtual void post_integrate();

 private:
  template <int> void post_integrate_eval();

  class FixPropertyAtom* fix_integrity_;
  double *integrity_;

  int every;
  int ago;
};

}

#endif
#endif
