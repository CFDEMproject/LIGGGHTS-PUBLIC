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
    (if no contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
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
