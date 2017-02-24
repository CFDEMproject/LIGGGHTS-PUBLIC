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

FixStyle(sph/pressure,FixSPHPressure)

#else

#ifndef LMP_FIX_SPH_PRESSURE_H
#define LMP_FIX_SPH_PRESSURE_H

#include "fix_sph.h"

namespace LAMMPS_NS {

  enum {PRESSURESTYLE_ABSOLUT,PRESSURESTYLE_TAIT,PRESSURESTYLE_RELATIV};

class FixSPHPressure : public FixSph {
 public:
  FixSPHPressure(class LAMMPS *, int, char **);
  ~FixSPHPressure();
  int setmask();
  void init();
  void pre_force(int);

  double return_rho0() {
    if (pressureStyle == PRESSURESTYLE_ABSOLUT) return 0;
    else return rho0;
  };

 private:
  int pressureStyle;
  double B,rho0,rho0inv,gamma,P0;
};

}

#endif
#endif
