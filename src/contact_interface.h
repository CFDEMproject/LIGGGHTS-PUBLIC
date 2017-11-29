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
    Christoph Kloss (JKU Linz)
    Richard Berger (JKU Linz)
    Alexander Podlozhnyuk (DCS Computing GmbH, Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef CONTACT_INTERFACE_H_
#define CONTACT_INTERFACE_H_

#include <string>

// forward declaration
namespace LAMMPS_NS
{
class TriMesh;
class FixMeshSurface;
}

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
  LAMMPS_NS::TriMesh *mesh;
  LAMMPS_NS::FixMeshSurface *fix_mesh;

  int i;
  int j;
  int itype;
  int jtype;

  bool is_wall;
  bool has_force_update;

  double * v_i;
  double * v_j;

  double * omega_i;
  double * omega_j;

  bool is_non_spherical;

#ifdef NONSPHERICAL_ACTIVE_FLAG
  double contact_point[3];
#endif

#ifdef SUPERQUADRIC_ACTIVE_FLAG
  double reff;
#endif

  int computeflag;
  int shearupdate;

  SurfacesCloseData() :
    radi(0.0),
    radj(0.0),
    radsum(0.0),
    rsq(0.0),
    area_ratio(1.0),
    contact_flags(NULL),
    contact_history(NULL),
    mesh(NULL),
    fix_mesh(NULL),
    i(0),
    j(0),
    itype(0),
    jtype(0),
    is_wall(false),
    has_force_update(false),
    v_i(NULL),
    v_j(NULL),
    omega_i(NULL),
    omega_j(NULL),
    is_non_spherical(false),
#ifdef SUPERQUADRIC_ACTIVE_FLAG
    reff(0.0),
#endif
    computeflag(0),
    shearupdate(0)
  {}
};

// data available in collision() only

struct SurfacesIntersectData : SurfacesCloseData {

  double r;         
  double rinv;      
  double en[3];     

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

  mutable double P_diss; 

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
