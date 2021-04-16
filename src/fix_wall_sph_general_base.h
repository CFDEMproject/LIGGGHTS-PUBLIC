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

// this class cannot be instantiated

#else

#ifndef LMP_FIX_WALL_SPH_GENERAL_BASE_H
#define LMP_FIX_WALL_SPH_GENERAL_BASE_H

#include "fix_wall_gran.h"

namespace LAMMPS_NS {

class FixWallSphGeneralBase : public FixWallGran {

 public:
  FixWallSphGeneralBase(class LAMMPS *, int, char **);
  ~FixWallSphGeneralBase();

  /* INHERITED FROM Fix */

  virtual int setmask();
  virtual void post_create();
  void pre_delete(bool unfixflag);
  virtual void post_integrate();
  void init();
  virtual void pre_force(int vflag);

  virtual void compute_density(int ip,double r,double mass);
  virtual void compute_velgrad(int ip,double delx,double dely,double delz,double mass,double *vwall);

 protected:
  class PairSph *pairsph_;
  int pairStyle;
  double viscosity;
  int kernel_id;
  int densityStyle;
  int modelStyle;
  int firstStep;
  double rho0;
  char *fixName;
  double **wallForce_;
  double **dvdx_;
  double **dvdy_;
  double **dvdz_;
  double *visc_;

  // storage for sph wall contact:
  class FixPropertyAtom* fix_wallContact_;
  double **wallContact_;

  // local viscosity value
  class FixPropertyAtom* fix_visc_;

 private:

  // mesh and primitive force implementations
  virtual void post_force_mesh(int);
  virtual void post_force_primitive(int);

  // force wall->fluidparticle
  class FixPropertyAtom* fix_wallForce_;

  // pressure fix
  class FixSPHPressure* fix_pressure_;

  // velocity gradients
  class FixSphVelgrad* fix_velgrad_;
  class FixPropertyAtom* fix_dvdx_;
  class FixPropertyAtom* fix_dvdy_;
  class FixPropertyAtom* fix_dvdz_;

// class FixSphDensitySumconti* fix_density_sumconti_;
};

}

#endif
#endif
