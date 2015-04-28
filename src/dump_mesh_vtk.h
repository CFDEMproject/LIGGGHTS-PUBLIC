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
    (if not contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS

DumpStyle(mesh/vtk,DumpMeshVTK)
DumpStyle(mesh/gran/VTK,DumpMeshVTK) 

#else

#ifndef LMP_DUMP_MESH_VTK_H
#define LMP_DUMP_MESH_VTK_H

#include "dump.h"
#include "container.h"

namespace LAMMPS_NS {

class DumpMeshVTK : public Dump {

 public:

  DumpMeshVTK(LAMMPS *, int, char**);
  virtual ~DumpMeshVTK();
  void init_style();

 private:            // column labels

  int dataMode_;

  int nMesh_;
  class TriMesh **meshList_;

  int n_calls_;

  // buffer for data from all procs
  int n_all_, n_all_max_;
  double *buf_all_;

  int dump_what_;

  // properties to be dumped
  // TODO: could make look-up more generic

  // stress
  class ScalarContainer<double> **sigma_n_, **sigma_t_;
  // wear
  class ScalarContainer<double> **wear_;
  // vel
  class MultiVectorContainer<double,3,3> **v_node_;
  // stresscomponents
  class VectorContainer<double,3> **f_node_;
  // temp
  class ScalarContainer<double> **T_;
  // min dist from active edge
  class ScalarContainer<double> **min_active_edge_dist_;

  // general implementation
  class ScalarContainer<double> ***scalar_containers_;
  char **scalar_container_names_;
  int n_scalar_containers_;
  class VectorContainer<double,3> ***vector_containers_;
  char **vector_container_names_;
  int n_vector_containers_;

  char **container_args_;
  int n_container_bases_;

  int modify_param(int, char **);
  void write_header(bigint ndump);
  int count();
  void getRefs();
  void getGeneralRefs();
  void pack(int *);
  void write_data(int, double *);

  void write_header_ascii(bigint ndump);
  void write_data_ascii(int n, double *mybuf);
  void write_data_ascii_point(int n, double *mybuf);
  void write_data_ascii_face(int n, double *mybuf);
};

}

#endif
#endif
