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
#include "dump_vtk.h"
#include "dump_local_gran.h"
#include <map>
#include <set>
#include <string>

#include <vtkSmartPointer.h>
#include <vtkMultiBlockDataSet.h>

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
class DumpLocalGranVTK : public Dump, public DumpVTK
{
 public:
  DumpLocalGranVTK(class LAMMPS *, int, char **);
  virtual ~DumpLocalGranVTK();

  virtual void write();

 protected:

  int nevery;                // dump frequency for output
  char *label;               // string for dump file header 
  DumpLocalGran *dumpLocalGran; // class that generates the vtk output

  int vtk_file_format;       // which vtk file format to write (vtk, vtp, vtu ...)

  // private methods

  virtual void init_style();
  virtual void write_header(bigint);
  int count();
  void pack(int *) {};
  virtual void write_data(int, double *);
  bigint memory_usage();

  virtual int modify_param(int, char **);

  typedef void (DumpLocalGranVTK::*FnPtrHeader)(bigint);
  FnPtrHeader header_choice;           // ptr to write header functions
  void header_vtk(bigint);

  typedef void (DumpLocalGranVTK::*FnPtrWrite)(int, double *);
  FnPtrWrite write_choice;             // ptr to write data functions
  void write_vtk(int, double *);
  void write_vtp(int, double *);
  void write_vtu(int, double *);

  // vtk data container
  vtkSmartPointer<vtkMultiBlockDataSet> mbSet;

  char *filecurrent;

  void setFileCurrent();
};

}

#endif
#endif
#endif

