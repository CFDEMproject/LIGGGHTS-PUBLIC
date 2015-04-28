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

    Christoph Kloss (DCS Computing GmbH, Linz, JKU Linz)
    Richard Berger (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef GRANULAR_PAIR_STYLE_H
#define GRANULAR_PAIR_STYLE_H

#include "utils.h"
#include "pair_gran.h"

namespace LIGGGHTS {

namespace PairStyles {

  /**
   * @brief interface of LIGGGHTS granular pair styles
   */
  class IGranularPairStyle {
  public:
    typedef LAMMPS_NS::PairGran ParentType;

    virtual ~IGranularPairStyle();
    virtual void settings(int nargs, char ** args) = 0;
    virtual void init_granular() = 0;
    virtual void write_restart_settings(FILE * fp) = 0;
    virtual void read_restart_settings(FILE * fp) = 0;
    virtual void compute_force(LAMMPS_NS::PairGran * pg, int eflag, int vflag, int addflag) = 0;

    virtual double stressStrainExponent() = 0;
    virtual int64_t hashcode() = 0;
  };

  /**
   * @brief LIGGGHTS granular pair style factory
   */
  class Factory : public Utils::AbstractFactory<IGranularPairStyle> {
    Factory() {}
  public:
    static Factory & instance();
  };
}

}

#endif // GRANULAR_PAIR_STYLE_H
