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
    Stefan Radl, TU Graz
    Copyright 2016 - TU Graz
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(massflow/mesh/sieve,FixMassflowMeshSieve)

#else

#ifndef LMP_FIX_MASSFLOW_MESH_SIEVE_H
#define LMP_FIX_MASSFLOW_MESH_SIEVE_H

#include "fix_massflow_mesh.h"
#include <vector>
#define PIOVERFOUR 0.78539816339

namespace LAMMPS_NS {

class FixMassflowMeshSieve : public FixMassflowMesh {

 public:

  FixMassflowMeshSieve(class LAMMPS *lmp, int narg, char ** arg);
  ~FixMassflowMeshSieve();

  void post_create();
//  void pre_delete(bool unfixflag);

  void init();
//  void setup(int vflag);
  int setmask();

  void post_force(int vflag);
  void pre_exchange();

 private:

  bool   sieveMultiSphereCanPass_;
  bool   verbose_;
  double sieveSize_;
  double sieveSpacing_;
  double sieveStiffness_;
  double sieveDamping_;
  double sievePassProbability(double);

 protected:
  char   fixidApplyForce_[200];
  class  FixPropertyAtom *fix_applyForce_;
  class  RanPark *random_;

}; //end class

}
#endif
#endif
