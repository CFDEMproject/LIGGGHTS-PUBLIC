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
    Alexander Podlozhnyuk (DCS Computing GmbH, Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef CONTACT_INTERFACE_H_
#define CONTACT_INTERFACE_H_

#include <string>
#include "superquadric.h"

namespace LIGGGHTS {
namespace ContactModels {

// data available in noCollision() and collision()

struct SurfacesCloseData {
  double radi;
  double radj;
  double radsum;
  double rsq;
  double delta[3];

  double area_ratio;

  int *contact_flags;
  double *contact_history;

  int i;
  int j;

  bool is_wall;
  bool has_force_update;

#ifdef SUPERQUADRIC_ACTIVE_FLAG
  double *quat_i; //quaternion of i-th particle
  double *quat_j; //quaternion of j-th particle
  double *shape_i; //shape parameters of i-th particle (a,b,c)
  double *shape_j; //shape parameters of j-th particle (a,b,c)
  double *roundness_i; //roundness parameters of i-th particle (eps1, eps2)
  double *roundness_j; //roundness parameters of j-th particle (eps1, eps2)
  double *pos_i;
  double *pos_j;
  SurfacesCloseData() : area_ratio(1.0),
                        quat_i(NULL),
                        quat_j(NULL),
                        shape_i(NULL),
                        shape_j(NULL),
                        roundness_i(NULL),
                        roundness_j(NULL),
                        pos_i(NULL),
                        pos_j(NULL){}
#else
  SurfacesCloseData() : area_ratio(1.0) {}
#endif
};

// data available in collision() only

struct SurfacesIntersectData : SurfacesCloseData {
  double r;
  double rinv;
  double en[3];
  double * v_i;
  double * v_j;
  double * omega_i;
  double * omega_j;

  double kt;
  double kn;
  double gammat;
  double gamman;

  double Fn;
  double Ft;

  double vn;
  double deltan;
  double cri;
  double crj;
  double wr1;
  double wr2;
  double wr3;

  double vtr1;
  double vtr2;
  double vtr3;

  double mi;
  double mj;
  double meff;

  int computeflag;
  int shearupdate;
  int itype;
  int jtype;

  SurfacesIntersectData() : Fn(0.0), Ft(0.0) {}
};

struct ForceData {
  double delta_F[3];       // total force acting on particle
  double delta_torque[3];  // torque acting on a particle

  ForceData()
  {
    reset();
  }

  inline void reset() {
    delta_F[0] = 0.0;
    delta_F[1] = 0.0;
    delta_F[2] = 0.0;
    delta_torque[0] = 0.0;
    delta_torque[1] = 0.0;
    delta_torque[2] = 0.0;
  }
};
}

class IContactHistorySetup {
public:
  virtual int add_history_value(std::string name, std::string newtonflag) = 0;
};

}

#endif /* CONTACT_INTERFACE_H_ */
