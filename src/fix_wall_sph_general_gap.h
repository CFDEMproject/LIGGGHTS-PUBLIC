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

FixStyle(wall/sph/general/gap,FixWallSphGeneralGap)

#else

#ifndef LMP_FIX_WALL_SPH_GENERAL_GAP_H
#define LMP_FIX_WALL_SPH_GENERAL_GAP_H

#include "fix_wall_sph_general_base.h"
#include "contact_interface.h"

namespace LCM = LIGGGHTS::ContactModels;

namespace LAMMPS_NS {

class FixWallSphGeneralGap : public FixWallSphGeneralBase {

   public:
      FixWallSphGeneralGap(class LAMMPS *, int, char **);
      ~FixWallSphGeneralGap();
      virtual void post_create();
      virtual void init();
      void pre_delete(bool unfixflag);
      virtual void post_integrate();
      void compute_density(int ip,double r,double mass);
      void compute_velgrad(int ip,double delx,double dely,double delz,double mass,double *vwall);
      void compute_force(LCM::SurfacesIntersectData & sidata, double *vwall);
      void compute_force_eval_single(int ip,double deltan,double r,double mass,
                                  double dx,double dy,double dz,double *vwall,
                                  double *c_history,double area_ratio);
      void compute_force_eval_gap(int ip,double mass,double r1,
                            double dx1,double dy1,double dz1,double *vwall1,double r2,
                            double dx2,double dy2,double dz2,double *vwall2);

    private:
      class FixPropertyAtom* fppaSl; //smoothing length
      class FixPropertyGlobal* fppaSlType; //per type smoothing length
      double *sl;         // per atom smoothing length
      double const*slType; // common smoothing length in case of mass_type=1

      class   FixPropertyGlobal* cs; // speed of sound
      double  const*csValues;

      // Counter for number of walls (1 -> single wall, 2 -> gap)
      class FixPropertyAtom* fix_wallCount_;
      double *wallCount_;

      class FixPropertyAtom* fix_wallContact2_;
      double **wallContact2_;

      class FixPropertyAtom* fix_wallContact3_;
      double **wallContact3_;

      class FixPropertyAtom* fix_integrity_;
      double *integrity_;

      class FixPropertyAtom* fix_fgradP_;
      double **fgradP_;

      class FixPropertyAtom* fix_wallForce2_;
      class FixPropertyAtom* fix_wallForce3_;

      double **wallForce2_;

      class FixPropertyAtom* fix_usedGapmodel_;
      double *usedGapmodel_;

    protected:

      double r0,D,vwallX,vwallY,vwallZ;
      int StaticWall;
      int mass_type; // flag defined in atom_vec*
      char *fixName;
      double gapWidth;
};

}

#endif
#endif
