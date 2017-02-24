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

#ifndef LMP_ACCELERATOR_CUDA_H
#define LMP_ACCELERATOR_CUDA_H

// true interface to USER-CUDA
// used when USER-CUDA is installed

#ifdef LMP_USER_CUDA

#include "cuda.h"
#include "comm_cuda.h"
#include "domain_cuda.h"
#include "neighbor_cuda.h"
#include "modify_cuda.h"
#include "verlet_cuda.h"

#else

// dummy interface to USER-CUDA
// needed for compiling when USER-CUDA is not installed

#include "comm.h"
#include "modify.h"
#include "verlet.h"

namespace LAMMPS_NS {

class Cuda {
 public:
  int cuda_exists;
  int oncpu;

  Cuda(class LAMMPS *) {cuda_exists = 0;}
  ~Cuda() {}
  void accelerator(int, char **) {}
  void evsetup_eatom_vatom(int, int) {}
  void downloadAll() {}
  void uploadAll() {}
};

class CommCuda : public Comm {
 public:
 CommCuda(class LAMMPS *lmp) : Comm(lmp) {}
  ~CommCuda() {}
};

class DomainCuda : public Domain {
 public:
 DomainCuda(class LAMMPS *lmp) : Domain(lmp) {}
  ~DomainCuda() {}
};

class NeighborCuda : public Neighbor {
 public:
 NeighborCuda(class LAMMPS *lmp) : Neighbor(lmp) {}
  ~NeighborCuda() {}
};

class ModifyCuda : public Modify {
 public:
 ModifyCuda(class LAMMPS *lmp) : Modify(lmp) {}
  ~ModifyCuda() {}
};

class VerletCuda : public Verlet {
 public:
 VerletCuda(class LAMMPS *lmp, int narg, char **arg) : Verlet(lmp,narg,arg) {}
  ~VerletCuda() {}
};

}

#endif
#endif
