/* ----------------------------------------------------------------------
LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
Transfer Simulations

www.liggghts.com | www.cfdem.com
Christoph Kloss, christoph.kloss@cfdem.com

   LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
   the producer of the LIGGGHTS® software and the CFDEM®coupling software
   See http://www.cfdem.com/terms-trademark-policy for details.

   LIGGGHTS® is based on LAMMPS
LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
http://lammps.sandia.gov, Sandia National Laboratories
Steve Plimpton, sjplimp@sandia.gov

Copyright (2003) Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
certain rights in this software. This software is distributed under
the GNU General Public License.

See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS

DumpStyle(decomposition/vtk,DumpDecompositionVTK)

#else

#ifndef LMP_DUMP_DECOMPOSITION_VTK_H
#define LMP_DUMP_DECOMPOSITION_VTK_H

#include "dump.h"

namespace LAMMPS_NS {

class DumpDecompositionVTK : public Dump {
 public:
  DumpDecompositionVTK(LAMMPS *, int, char**);
  ~DumpDecompositionVTK();
  void init_style();

 private:
  int len[3];
  double *xdata, *xdata_all;
  double *ydata, *ydata_all;
  double *zdata, *zdata_all;

  int lasttimestep;

  int modify_param(int, char **);
  void write_header(bigint);
  int count();
  void pack(int *);
  void write_data(int, double *);

  typedef void (DumpDecompositionVTK::*FnPtrHeader)(bigint);
  FnPtrHeader header_choice;           // ptr to write header functions
  void header_item(bigint);
  void footer_item();

  typedef void (DumpDecompositionVTK::*FnPtrPack)();
  FnPtrPack pack_choice;               // ptr to pack functions
  void pack_item();

  typedef void (DumpDecompositionVTK::*FnPtrData)(int, double *);
  FnPtrData write_choice;              // ptr to write data functions
  void write_item(int, double *);

};

}

#endif
#endif
