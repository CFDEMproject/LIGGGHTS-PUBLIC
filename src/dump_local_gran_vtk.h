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

    Copyright 2014-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#if defined(LAMMPS_VTK) 
#ifdef DUMP_CLASS

DumpStyle(local/gran/vtk,DumpLocalGranVTK)

#else

#ifndef LMP_DUMP_LOCAL_GRAN_VTK_H
#define LMP_DUMP_LOCAL_GRAN_VTK_H

#include "dump.h"
#include <map>
#include <set>
#include <string>

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>

class vtkAbstractArray;
class vtkRectilinearGrid;
class vtkUnstructuredGrid;

namespace LAMMPS_NS {

/**
 * @brief DumpLocalGranVTK class
 *        write gran bond data to vtk files.
 *
 * Similar to the DumpCustomVTK class but uses the vtk library to write data to vtk simple
 * legacy or xml format depending on the filename extension specified. (Since this
 * conflicts with the way binary output is specified, dump_modify allows to set the
 * binary flag for this dump command explicitly).
 * In contrast to DumpCustom class the attributes to be packed are stored in a std::map
 * to avoid duplicate entries and enforce correct ordering of vector components (except
 * for computes and fixes - these have to be given in the right order in the input script).
 * (Note: std::map elements are sorted by their keys.)
 * This dump command does not support compressed files, buffering or custom format strings,
 * multiproc is only supported by the xml formats, multifile option has to be used.
 */
class DumpLocalGranVTK : public Dump {
 public:
  DumpLocalGranVTK(class LAMMPS *, int, char **);
  virtual ~DumpLocalGranVTK();

  virtual void write();

 protected:

  int nevery;                // dump frequency for output
  char *label;               // string for dump file header 
  int iregion;               // -1 if no region, else which region
  char *idregion;            // region ID

  int vtk_file_format;       // which vtk file format to write (vtk, vtp, vtu ...)

  int nchoose;               // # of selected local datums
  int maxlocal;              // size of atom selection and variable arrays

  // private methods

  class ComputePairGranLocal *cpgl_;

  virtual void init_style();
  virtual void write_header(bigint);
  int count();
  void pack(int *);
  virtual void write_data(int, double *);
  bigint memory_usage();

  int parse_fields(int, char **);
  int add_compute(char *);
  int add_fix(char *);
  int add_variable(char *);
  virtual int modify_param(int, char **);

  typedef void (DumpLocalGranVTK::*FnPtrHeader)(bigint);
  FnPtrHeader header_choice;           // ptr to write header functions
  void header_vtk(bigint);

  typedef void (DumpLocalGranVTK::*FnPtrWrite)(int, double *);
  FnPtrWrite write_choice;             // ptr to write data functions
  void write_vtk(int, double *);
  void write_vtp(int, double *);
  void write_vtu(int, double *);

  void define_properties();
  typedef void (DumpLocalGranVTK::*FnPtrPack)(int);
  
  std::map<int, FnPtrPack> pack_choice;  // ptrs to pack functions
  std::map<int, int> vtype;              // data type for each attribute
  std::map<int, std::string> name;       // label for each attribute
  std::set<int> vector_set;              // set of vector attributes; defines which are vectors

  // vtk data containers
  vtkSmartPointer<vtkPoints> points;                            // list of points, 2 points for each line cell
  vtkSmartPointer<vtkCellArray> lineCells;                      // list of line cells
  std::map<int, vtkSmartPointer<vtkAbstractArray> > myarrays;   // list of a list of arrays that is presents data for each atom (x, v,...)
                                                                // is then added to the point cells upon writing
  int n_calls_;

  char *filecurrent;
  char *parallelfilecurrent;
  char *multiname_ex;

  void setFileCurrent();
  void buf2arrays(int, double *); // transfer data from buf array to vtk arrays
  void reset_vtk_data_containers();

  // customize by adding a method prototype
  void pack_x1(int);
  void pack_x2(int);
  void pack_v1(int);
  void pack_v2(int);
  void pack_id1(int);
  void pack_id2(int);
  void pack_id3(int);
  void pack_f(int);
  void pack_fn(int);
  void pack_ft(int);
  void pack_torque(int);
  void pack_torquen(int);
  void pack_torquet(int);
  void pack_area(int);
  void pack_delta(int);
  void pack_heat(int);
};

}

#endif
#endif
#endif

