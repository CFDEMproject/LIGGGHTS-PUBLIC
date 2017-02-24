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

#ifndef LMP_MIN_LSRCH_H
#define LMP_MIN_LSRCH_H

#include "min.h"

namespace LAMMPS_NS {

class MinLineSearch : public Min {
 public:
  MinLineSearch(class LAMMPS *);
  ~MinLineSearch();
  void init();
  void setup_style();
  void reset_vectors();

 protected:
  // vectors needed by linesearch minimizers
  // allocated and stored by fix_minimize
  // x,f are stored by parent or Atom class or Pair class

  double *x0;                 // coords at start of linesearch
  double *g;                  // old gradient vector
  double *h;                  // search direction vector

  double *gextra;             // g,h for extra global dof, x0 is stored by fix
  double *hextra;

  double **x0extra_atom;      // x0,g,h for extra per-atom dof
  double **gextra_atom;
  double **hextra_atom;

  typedef int (MinLineSearch::*FnPtr)(double, double &);
  FnPtr linemin;
  int linemin_backtrack(double, double &);
  int linemin_quadratic(double, double &);
  int linemin_forcezero(double, double &);

  double alpha_step(double, int);
  double compute_dir_deriv(double &);
};

}

#endif
