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
    This file is from LAMMPS
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

#ifndef LMP_INTEGRATE_H
#define LMP_INTEGRATE_H

#include "pointers.h"

namespace LAMMPS_NS {

class Integrate : protected Pointers {
 public:
  Integrate(class LAMMPS *, int, char **);
  virtual ~Integrate();
  virtual void init();
  virtual void setup() = 0;
  virtual void setup_minimal(int) = 0;
  virtual void run(int) = 0;
  virtual void cleanup() {}
  virtual void reset_dt() {}
  virtual bigint memory_usage() {return 0;}

 protected:
  int eflag,vflag;                  // flags for energy/virial computation
  int virial_style;                 // compute virial explicitly or implicitly
  int external_force_clear;         // clear forces locally or externally

  int nelist_global,nelist_atom;    // # of PE,virial computes to check
  int nvlist_global,nvlist_atom;
  class Compute **elist_global;     // lists of PE,virial Computes
  class Compute **elist_atom;
  class Compute **vlist_global;
  class Compute **vlist_atom;

  int pair_compute_flag;            // 0 if pair->compute is skipped
  int kspace_compute_flag;          // 0 if kspace->compute is skipped

  void ev_setup();
  void ev_set(bigint);
};

}

#endif
/* ERROR/WARNING messages:

*/
