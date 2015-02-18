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

/* ----------------------------------------------------------------------
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Philippe Seil (JKU Linz)
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS

DumpStyle(mesh/stl,DumpMeshSTL)
DumpStyle(stl,DumpMeshSTL) 

#else

#ifndef LMP_DUMP_MESH_STL_H
#define LMP_DUMP_MESH_STL_H

#include "dump.h"

namespace LAMMPS_NS {

class DumpMeshSTL : public Dump {
 public:
  DumpMeshSTL(LAMMPS *, int, char**);
  virtual ~DumpMeshSTL();
  void init_style();

 private:            // column labels

  int nMesh_;
  class TriMesh **meshList_;
  int dump_what_;

  int n_calls_;

  int writeBinarySTL_;

  // region filter
  int iregion_;

  int modify_param(int, char **);
  void write_header(bigint ndump);
  int count();
  void bounds(int imesh,int &ilo, int &ihi);
  void pack(int *);
  void write_data(int, double *);

  void write_header_ascii(bigint ndump);
  void write_header_binary(bigint ndump);

  void write_data_ascii(int n, double *mybuf);
  void write_data_binary(int n, double *mybuf);

};

}

#endif
#endif
