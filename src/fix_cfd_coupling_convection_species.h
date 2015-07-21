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

    Copyright 2015-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(couple/cfd/speciesConvection,FixCfdCouplingConvectionSpecies)

#else

#ifndef LMP_FIX_CFD_COUPLING_CONVECTION_SPECIES_H
#define LMP_FIX_CFD_COUPLING_CONVECTION_SPECIES_H

#include "fix_cfd_coupling.h"

namespace LAMMPS_NS {

class FixCfdCouplingConvectionSpecies : public Fix {

 public:
  FixCfdCouplingConvectionSpecies(class LAMMPS *, int, char **);
  ~FixCfdCouplingConvectionSpecies();
  void post_create();
  void pre_delete(bool unfixflag);

  virtual int  setmask();
  virtual void init();
  virtual void post_force(int);

 protected:
  class FixCfdCoupling*  fix_coupling;
  class FixPropertyAtom* fix_speciesConcentration;
  class FixPropertyAtom* fix_convectiveFlux;
  class FixPropertyAtom* fix_totalFlux;

  double species0;
  char   speciesName_[128];
  char   sourceName_[128];
  char   convectiveFluxName_[128];
  char   capacityName_[128];
  char   steName_[128];
  char   totalFluxName_[128];
};

}

#endif
#endif
