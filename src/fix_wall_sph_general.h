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

FixStyle(wall/sph/general,FixWallSphGeneral)

#else

#ifndef LMP_FIX_WALL_SPH_GENERAL_H
#define LMP_FIX_WALL_SPH_GENERAL_H

#include "fix_wall_sph_general_base.h"
#include "contact_interface.h"

namespace LCM = LIGGGHTS::ContactModels;

namespace LAMMPS_NS {

class FixWallSphGeneral : public FixWallSphGeneralBase {

   public:
      FixWallSphGeneral(class LAMMPS *, int, char **);
      ~FixWallSphGeneral();
      void compute_density(int ip,double r,double mass);
      void compute_velgrad(int ip,double delx,double dely,double delz,double mass,double *vwall);
      void compute_force(LCM::SurfacesIntersectData & sidata, double *vwall);

    private:
      class FixPropertyAtom* fppaSl; //smoothing length
      class FixPropertyGlobal* fppaSlType; //per type smoothing length
      double *sl;         // per atom smoothing length
      double const*slType; // common smoothing length in case of mass_type=1

      class   FixPropertyGlobal* cs; // speed of sound
      double  const*csValues;

    protected:

      // SPH parameters
      double r0,D,vwallX,vwallY,vwallZ;
      int StaticWall;
      int mass_type; // flag defined in atom_vec*
};

}

#endif
#endif
