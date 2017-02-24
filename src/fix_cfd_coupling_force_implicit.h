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

FixStyle(couple/cfd/force/implicit,FixCfdCouplingForceImplicit)

#else

#ifndef LMP_FIX_CFD_COUPLING_FORCE_IMPLICIT_H
#define LMP_FIX_CFD_COUPLING_FORCE_IMPLICIT_H

#include "fix_cfd_coupling_force.h"

namespace LAMMPS_NS {

class FixCfdCouplingForceImplicit : public FixCfdCouplingForce  {
  friend class FixNVEAsphereBase;
 public:
  FixCfdCouplingForceImplicit(class LAMMPS *, int, char **);
  ~FixCfdCouplingForceImplicit();
  void post_create();
  void pre_delete(bool unfixflag);

  int setmask();
  virtual void init();
  void post_force(int);
  void end_of_step();

 protected:
  double deltaT_;

  bool   useCN_;
  double CNalpha_;

  bool   useAM_;
  double CAddRhoFluid_;   //Added mass coefficient times relative fluid density (C_add*rhoFluid/rhoP)
  double onePlusCAddRhoFluid_;

  class FixPropertyAtom* fix_Ksl_;
  class FixPropertyAtom* fix_uf_;
  class FixPropertyAtom* fix_KslRotation_;
  class FixPropertyAtom* fix_ex_;
  class FixPropertyAtom* fix_KslExtra_;
};

}

#endif
#endif
