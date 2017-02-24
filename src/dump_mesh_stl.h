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

    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Philippe Seil (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
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
