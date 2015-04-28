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

#ifdef PAIR_CLASS

PairStyle(sph/artVisc/tensCorr,PairSphArtviscTenscorr)

#else

#ifndef LMP_PAIR_SPH_ARTVISC_TENSCORR_H
#define LMP_PAIR_SPH_ARTVISC_TENSCORR_H

#include "pair_sph.h"

namespace LAMMPS_NS {

class PairSphArtviscTenscorr : public PairSph {

 friend class FixSPH;

 public:

  PairSphArtviscTenscorr(class LAMMPS *);
  ~PairSphArtviscTenscorr();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_substyle();
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);

 protected:
  void allocate();
  template <int> void compute_eval(int, int);

  int     artVisc_flag, tensCorr_flag; // flags for additional styles

  class   FixPropertyGlobal* cs;
  class   FixPropertyGlobal* alpha;
  class   FixPropertyGlobal* beta;
  class   FixPropertyGlobal* etaPPG;
  double  **csmean,**alphaMean,**betaMean;
  double  eta;

  class   FixPropertyGlobal* epsilonPPG;
  class   FixPropertyGlobal* deltaP;
  double  **wDeltaPTypeinv;
  double  epsilon; // coeffs for tensile correction

};

}

#endif
#endif
