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
    Philippe Seil (JKU Linz)

    Copyright 2014-     JKU Linz
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(couple/lb/onetoone,FixLbCouplingOnetoone)

#else

#ifndef FIX_LB_COUPLING_ONETOONE
#define FIX_LB_COUPLING_ONETOONE

#include "fix.h"

namespace LAMMPS_NS {
  class FixLbCouplingOnetoone : public Fix {
  public:
    FixLbCouplingOnetoone(class LAMMPS * lmp, int narg, char ** arg);
    ~FixLbCouplingOnetoone();

    virtual int setmask();

    virtual void post_create();
    virtual void pre_delete(bool);
    virtual void init();

    virtual void post_force(int);

    double **get_force_ptr();
    double **get_torque_ptr();
    void comm_force_torque();
  private:
    class FixPropertyAtom* fix_dragforce_;
    class FixPropertyAtom* fix_hdtorque_; // hdtorque = hydrodynamic torque
  }; /* class FixLbCouplingOnetoone */

}; /* LAMMPS_NS */

#endif /* FIX_LB_COUPLING_ONETOONE */
#endif /* FIX_CLASS */

