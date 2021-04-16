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

    Andreas Eitzlmayr (TU Graz)

    Copyright 2013-     TU Graz
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(sph/mixidx,FixSphMixidx)

#else

#ifndef LMP_FIX_SPH_MIXIDX_H
#define LMP_FIX_SPH_MIXIDX_H

#include "fix_sph.h"

namespace LAMMPS_NS {

class FixSphMixidx : public FixSph {
 public:
  FixSphMixidx(class LAMMPS *, int, char **);
  ~FixSphMixidx();
  virtual int setmask();
  void post_create();
  void pre_delete(bool unfixflag);
  virtual void init();
  void post_force(int);

 private:
  template <int> void post_force_eval(int);

  class FixPropertyAtom* fix_dvdx_;
  class FixPropertyAtom* fix_dvdy_;
  class FixPropertyAtom* fix_dvdz_;
  class FixPropertyAtom* fix_gamma_; // shear rate magnitude
  class FixPropertyAtom* fix_omega_; // vorticity magnitude
  class FixPropertyAtom* fix_mixidx_;

  double **dvdx_;
  double **dvdy_;
  double **dvdz_;
  double *gamma_;
  double *omega_;
  double *mixidx_;

  int every;
  int ago;
  int calcgamma;

};

}

#endif
#endif
