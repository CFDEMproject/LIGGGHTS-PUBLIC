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
    Alexander Podlozhnyuk, DCS Computing GmbH, Linz
    Christoph Kloss, DCS Computing GmbH, Linz

    Copyright 2015-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

#else

#ifndef LMP_FIX_NVE_ASPHERE_BASE_H
#define LMP_FIX_NVE_ASPHERE_BASE_H

#include "fix_nve.h"
#include "fix_cfd_coupling_force_implicit.h"

namespace LAMMPS_NS {

class FixNVEAsphereBase : public FixNVE {
 public:
  FixNVEAsphereBase(class LAMMPS *, int, char **);
  virtual ~FixNVEAsphereBase() {}
  virtual void init();
  virtual void initial_integrate(int);
  virtual void final_integrate();
  void dynamic_euler(double *wbody, double *tbody, double *inertia, double *result);
  void integrate_dynamic_euler(double dt, double *wbody, double *tbody, double *inertia);
  void integrate_quaternion(double dtq, double *quat, double *wbody);
  void update_hdtorque(int i, double *rotation_matrix, double *omegaOld, double *omegaNew);
  void rotationUpdate(bool updateQuaternion);
  void implicitRotationUpdate(
      double deltaT, double* inertia,
      double *angMom, double *torque, double* KslRot,
      double *omegaNew,
      double *deltaHydrotorquePrime
  );

 private:
  int integration_scheme;
  int couple_fix_id;
  double **ksl_rotation, **orientation, **hdtorque;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix nve/superquadric requires atom style superquadric

Self-explanatory.

*/

