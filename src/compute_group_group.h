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

#ifdef COMPUTE_CLASS

ComputeStyle(group/group,ComputeGroupGroup)

#else

#ifndef LMP_COMPUTE_GROUP_GROUP_H
#define LMP_COMPUTE_GROUP_GROUP_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeGroupGroup : public Compute {
 public:
  ComputeGroupGroup(class LAMMPS *, int, char **);
  ~ComputeGroupGroup();
  void init();
  void init_list(int, class NeighList *);
  double compute_scalar();
  void compute_vector();

 private:
  char *group2;
  int jgroup,jgroupbit,othergroupbit;
  double **cutsq;
  double e_self,e_correction;
  int pairflag,kspaceflag,boundaryflag;
  class Pair *pair;
  class NeighList *list;
  class KSpace *kspace;

  void pair_contribution();
  void kspace_contribution();
  void kspace_correction();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute group/group group ID does not exist

Self-explanatory.

E: No pair style defined for compute group/group

Cannot calculate group interactions without a pair style defined.

E: Pair style does not support compute group/group

The pair_style does not have a single() function, so it cannot be
invokded by the compute group/group command.

E: No Kspace style defined for compute group/group

Self-explanatory.

E: Kspace style does not support compute group/group

Self-explanatory.

W: Both groups in compute group/group have a net charge; the Kspace boundary correction to energy will be non-zero

Self-explantory.

*/
