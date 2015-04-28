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

    Richard Berger (JKU Linz)
    Christoph Kloss (DCS Computing GmbH, Linz, JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
PairStyle(gran,PairGranProxy)
#else

#ifndef PAIR_GRAN_PROXY_H
#define PAIR_GRAN_PROXY_H

#include "pair_gran.h"
#include "granular_pair_style.h"

namespace LAMMPS_NS
{
class PairGranProxy : public PairGran
{
  LIGGGHTS::PairStyles::IGranularPairStyle * impl;

public:
  PairGranProxy(LAMMPS * lmp);
  virtual ~PairGranProxy();

  virtual void settings(int nargs, char ** args);
  virtual void init_granular();
  virtual void write_restart_settings(FILE * fp);
  virtual void read_restart_settings(FILE * fp);
  virtual void compute_force(int eflag, int vflag, int addflag);

  virtual double stressStrainExponent();
  virtual int64_t hashcode();
};
}

#endif // PAIR_GRAN_PROXY_H

#endif
