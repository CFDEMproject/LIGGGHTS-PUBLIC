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

FixStyle(check/timestep/gran,FixCheckTimestepGran)

#else

#ifndef LMP_FIX_CHECK_TIMESTEP_GRAN_H
#define LMP_FIX_CHECK_TIMESTEP_GRAN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixCheckTimestepGran : public Fix {
 public:
  FixCheckTimestepGran(class LAMMPS *, int, char **);
  int setmask();
  void init();
  void end_of_step();
  double compute_vector(int);

 private:
  class Properties* properties;
  class PairGran* pg;
  class FixWallGran* fwg;
  class FixPropertyGlobal* Y;
  class FixPropertyGlobal* nu;
  void calc_rayleigh_hertz_estims();
  double rayleigh_time,hertz_time;
  double fraction_rayleigh,fraction_hertz,fraction_skin;
  double fraction_rayleigh_lim,fraction_hertz_lim;
  double vmax; //max relative velocity
  double r_min;
  bool warnflag;
  double ** Yeff;
};

}

#endif
#endif
