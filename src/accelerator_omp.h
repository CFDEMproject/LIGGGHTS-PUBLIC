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

// NOTE: this file is *supposed* to be included multiple times

#ifdef LMP_USER_OMP

// true interface to USER-OMP

// this part is used inside the neighbor.h header file to
// add functions to the Neighbor class definition

#ifdef LMP_INSIDE_NEIGHBOR_H

  void half_nsq_no_newton_omp(class NeighList *);
  void half_nsq_no_newton_ghost_omp(class NeighList *);
  void half_nsq_newton_omp(class NeighList *);

  void half_bin_no_newton_omp(class NeighList *);
  void half_bin_no_newton_ghost_omp(class NeighList *);
  void half_bin_newton_omp(class NeighList *);
  void half_bin_newton_tri_omp(class NeighList *);

  void half_multi_no_newton_omp(class NeighList *);
  void half_multi_newton_omp(class NeighList *);
  void half_multi_newton_tri_omp(class NeighList *);

  void full_nsq_omp(class NeighList *);
  void full_nsq_ghost_omp(class NeighList *);
  void full_bin_omp(class NeighList *);
  void full_bin_ghost_omp(class NeighList *);
  void full_multi_omp(class NeighList *);

  void half_from_full_no_newton_omp(class NeighList *);
  void half_from_full_newton_omp(class NeighList *);

  void granular_nsq_no_newton_omp(class NeighList *);
  void granular_nsq_newton_omp(class NeighList *);
  void granular_bin_no_newton_omp(class NeighList *);
  void granular_bin_newton_omp(class NeighList *);
  void granular_bin_newton_tri_omp(class NeighList *);

  void respa_nsq_no_newton_omp(class NeighList *);
  void respa_nsq_newton_omp(class NeighList *);
  void respa_bin_no_newton_omp(class NeighList *);
  void respa_bin_newton_omp(class NeighList *);
  void respa_bin_newton_tri_omp(class NeighList *);

#else /* !LMP_INSIDE_NEIGHBOR_H */

// provide a DomainOMP class with some overrides for Domain
#include "domain.h"

#ifndef LMP_DOMAIN_OMP_H
#define LMP_DOMAIN_OMP_H

namespace LAMMPS_NS {

class DomainOMP : public Domain {
 public:
  DomainOMP(class LAMMPS *lmp) : Domain(lmp) {};
  virtual ~DomainOMP() {};

  // multi-threaded versions
  virtual void pbc();
  virtual void lamda2x(int);
  virtual void lamda2x(double *lamda, double *x) {Domain::lamda2x(lamda,x);}
  virtual void x2lamda(int);
  virtual void x2lamda(double *x, double *lamda) {Domain::x2lamda(x,lamda);}
};
}

#endif /* LMP_DOMAIN_OMP_H */
#endif /* !LMP_INSIDE_NEIGHBOR_H */

#else /* !LMP_USER_OMP */

// needed for compiling Neighbor class when USER-OMP is not installed

#ifdef LMP_INSIDE_NEIGHBOR_H

  void half_nsq_no_newton_omp(class NeighList *) {}
  void half_nsq_no_newton_ghost_omp(class NeighList *) {}
  void half_nsq_newton_omp(class NeighList *) {}

  void half_bin_no_newton_omp(class NeighList *) {}
  void half_bin_no_newton_ghost_omp(class NeighList *) {}
  void half_bin_newton_omp(class NeighList *) {}
  void half_bin_newton_tri_omp(class NeighList *) {}

  void half_multi_no_newton_omp(class NeighList *) {}
  void half_multi_newton_omp(class NeighList *) {}
  void half_multi_newton_tri_omp(class NeighList *) {}

  void full_nsq_omp(class NeighList *) {}
  void full_nsq_ghost_omp(class NeighList *) {}
  void full_bin_omp(class NeighList *) {}
  void full_bin_ghost_omp(class NeighList *) {}
  void full_multi_omp(class NeighList *) {}

  void half_from_full_no_newton_omp(class NeighList *) {}
  void half_from_full_newton_omp(class NeighList *) {}

  void granular_nsq_no_newton_omp(class NeighList *) {}
  void granular_nsq_newton_omp(class NeighList *) {}
  void granular_bin_no_newton_omp(class NeighList *) {}
  void granular_bin_newton_omp(class NeighList *) {}
  void granular_bin_newton_tri_omp(class NeighList *) {}

  void respa_nsq_no_newton_omp(class NeighList *) {}
  void respa_nsq_newton_omp(class NeighList *) {}
  void respa_bin_no_newton_omp(class NeighList *) {}
  void respa_bin_newton_omp(class NeighList *) {}
  void respa_bin_newton_tri_omp(class NeighList *) {}
#endif

#endif /* !LMP_USER_OMP */
