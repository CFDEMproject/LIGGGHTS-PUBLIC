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
    Anton Gladky(TU Bergakademie Freiberg), gladky.anton@gmail.com
------------------------------------------------------------------------- */

#if defined(LAMMPS_VTK) // NP do not use #ifdef here (VS C++ bug)
#ifdef DUMP_CLASS

DumpStyle(atom/vtk,DumpATOMVTK)

#else

#ifndef LMP_DUMP_ATOM_VTK_H
#define LMP_DUMP_ATOM_VTK_H

#include "dump.h"
#include <iostream>
#include <vector>
#include <fstream>
#include "update.h"

namespace LAMMPS_NS {

class DumpATOMVTK : public Dump {
 public:
  DumpATOMVTK(class LAMMPS *, int, char**);

 private:
  void init_style();
  void write_header(bigint);
  int count();
  void pack(int *);
  void write_data(int, double *);

  int n_calls_;
  char * filecurrent;
  void setFileCurrent();

  class V3 {
    public:
      V3() {v[0]=v[1]=v[2]=0.0;}
      V3(double x, double y, double z) {v[0]=x;v[1]=y;v[2]=z;}
      double& operator[](int idx) {return v[idx];}
      const double& operator[](int idx) const {return v[idx];}
      double v[3];
  };

  typedef DumpATOMVTK::V3 V3;

  class DataVTK {
    public:
      V3 _Pos;
      double _Rad;
      double _Mass;
      int _Id;
      int _Type;
      V3 _VelL;
      V3 _VelA;
      V3 _Force;
      int _proc;
      DataVTK(V3, double, double, int, int, V3, V3, V3, int);
      std::string serialize();
  };

  class vtkExportData {
    private:
      std::vector<DumpATOMVTK::DataVTK> vtkData;
      std::ofstream fileVTK;
      const char * _fileName;
      bool _setFileName;
    public:
      vtkExportData();
      void add(DumpATOMVTK::DataVTK &);
      int size();
      void writeSER();
      void setFileName(const char *);
      void show();
      void clear();
  };

  vtkExportData tmpEXP;

};

}

#endif
#endif
#endif
