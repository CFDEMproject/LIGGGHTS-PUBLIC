/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author:
   Anton Gladky(TU Bergakademie Freiberg), gladky.anton@gmail.com
------------------------------------------------------------------------- */

#ifdef LAMMPS_VTK
#ifdef DUMP_CLASS

DumpStyle(atom/vtk,DumpATOMVTK)

#else

#ifndef LMP_DUMP_ATOM_VTK_H
#define LMP_DUMP_ATOM_VTK_H

#include "dump.h"
#include <iostream>
#include <vector>
#include <fstream>

#include<vtkCellArray.h>
#include<vtkFloatArray.h>
#include<vtkDoubleArray.h>
#include<vtkIntArray.h>
#include<vtkPoints.h>
#include<vtkPointData.h>
#include<vtkCellData.h>
#include<vtkSmartPointer.h>
#include<vtkUnstructuredGrid.h>
#include<vtkPolyData.h>
#include<vtkXMLUnstructuredGridWriter.h>
#include<vtkXMLPolyDataWriter.h>
#include<vtkZLibDataCompressor.h>
#include<vtkTriangle.h>
#include<vtkLine.h>
#include<vtkQuad.h>

#include <Eigen/Core>
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

  typedef Eigen::Matrix<double, 3, 1> V3;
  int n_calls_;
  char * filecurrent;
  void setFileCurrent();

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
      ofstream fileVTK;
      const char * _fileName;
      bool _setFileName;
    public:
      vtkExportData();
      void add(DumpATOMVTK::DataVTK &);
      const int size();
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
