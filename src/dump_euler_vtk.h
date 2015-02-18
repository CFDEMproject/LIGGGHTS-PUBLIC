/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
   the producer of the LIGGGHTS® software and the CFDEM®coupling software
   See http://www.cfdem.com/terms-trademark-policy for details.

   LIGGGHTS® is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS

DumpStyle(euler/vtk,DumpEulerVTK)

#else

#ifndef LMP_DUMP_EULER_VTK_H
#define LMP_DUMP_EULER_VTK_H

#include "dump.h"

namespace LAMMPS_NS {

class DumpEulerVTK : public Dump {

 public:

  DumpEulerVTK(LAMMPS *, int, char**);
  virtual ~DumpEulerVTK();
  void init_style();

 private:

  class FixAveEuler *fix_euler_;

  int modify_param(int, char **);
  void write_header(bigint ndump);
  int count();
  void pack(int *);
  void write_data(int, double *);

  void write_header_ascii(bigint ndump);
  void write_data_ascii(int n, double *mybuf);

  int n_calls_;

  // buffer for data from all procs
  int n_all_, n_all_max_;
  double *buf_all_;

};

}

#endif
#endif
