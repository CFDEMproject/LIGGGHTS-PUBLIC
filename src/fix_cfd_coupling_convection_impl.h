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

    Christoph Kloss (DCS Computing GmbH, Linz)
    Stefan Radl (TU Graz)

    Copyright 2015-     DCS Computing GmbH, Linz
    Copyright 2015-     TU Graz
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(couple/cfd/convectiveImpl,FixCfdCouplingConvectiveImpl)

#else

#ifndef LMP_FIX_CFD_COUPLING_CONVECTIVE_IMPL_H
#define LMP_FIX_CFD_COUPLING_CONVECTIVE_IMPL_H

#include "fix_cfd_coupling.h"

namespace LAMMPS_NS {

class FixCfdCouplingConvectiveImpl : public Fix {

 public:
  FixCfdCouplingConvectiveImpl(class LAMMPS *, int, char **);
  ~FixCfdCouplingConvectiveImpl();
  void post_create();
  void pre_delete(bool unfixflag);

  virtual int setmask();
  virtual void init();
  virtual void post_force(int);

 protected:
  bool  integrateHeatEqn_;      //set to true to activate integration of heat flux with a scalar transport equation
  bool  forceExplicit_;         //force explicit calculation: i.e., aggregate all fluxes into heatFlux
  class FixCfdCoupling*     fix_coupling;
  class FixPropertyAtom*    fix_heatFluid;
  class FixPropertyAtom*    fix_heatTransCoeff;
  class FixPropertyAtom*    fix_convectiveFlux;

  class FixPropertyAtom*    fix_temperature;  //only needed in case heat integration is performed
  class FixPropertyAtom*    fix_heatFlux;  //only needed in case heat integration is performed
  double T0;
};

}

#endif
#endif
