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

/* ----------------------------------------------------------------------
   Contributing authors:
   Daniel Queteschiner, daniel.queteschiner@dcs-computing.com
   Alexander Podlozhnyuk, alexander.podlozhnyuk@dcs-computing.com
------------------------------------------------------------------------- */

#if defined(LAMMPS_VTK) 

#ifdef DUMP_CLASS

DumpStyle(custom/vtk,DumpCustomVTK)

#else

#ifndef LMP_DUMP_CUSTOM_VTK_H
#define LMP_DUMP_CUSTOM_VTK_H

#include "dump.h"
#include "dump_vtk.h"
#include "dump_particle.h"
#include <map>
#include <set>
#include <string>

#include <vtkSmartPointer.h>
#include <vtkMultiBlockDataSet.h>

class vtkAbstractArray;
class vtkRectilinearGrid;
class vtkUnstructuredGrid;

namespace LAMMPS_NS {

/**
 * @brief DumpCustomVTK class
 *        write atom data to vtk files.
 *
 * Similar to the DumpCustom class but uses the vtk library to write data to vtk simple
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
class DumpCustomVTK : public Dump, public DumpVTK
{
 public:
  DumpCustomVTK(class LAMMPS *, int, char **);
  virtual ~DumpCustomVTK();

  virtual void write();
 protected:
  char *label;               // string for dump file header 
  DumpParticle *dumpParticle;
  vtkSmartPointer<vtkMultiBlockDataSet> mbSet;

  int nevery;                // dump frequency for output
  int vtk_file_format;       // which vtk file format to write (vtk, vtp, vtu ...)

  int nfield;                // # of keywords listed by user
  int ioptional;             // index of start of optional args

  std::map<int, int> field2index; // which compute,fix,variable calcs this field
  std::map<int, int> argindex;    // index into compute,fix scalar_atom,vector_atom
                                  // 0 for scalar_atom, 1-N for vector_atom values

  // private methods

  virtual void init_style();
  virtual void write_header(bigint);
  int count();
  virtual void write_data(int, double *);
  bigint memory_usage();
  void pack(int *) {}

  virtual int modify_param(int, char **);

  typedef void (DumpCustomVTK::*FnPtrHeader)(bigint);
  FnPtrHeader header_choice;           // ptr to write header functions
  void header_vtk(bigint);

  typedef void (DumpCustomVTK::*FnPtrWrite)(int, double *);
  FnPtrWrite write_choice;             // ptr to write data functions
  void write_vtk(int, double *);
  void write_vtp(int, double *);
  void write_vtu(int, double *);

  void write_domain_vtk();
  void write_domain_vtk_triclinic();
  void write_domain_vtr();
  void write_domain_vtu_triclinic();

                                                                // is then added to the point cells upon writing

  char *filecurrent;
  char *domainfilecurrent;

  void setFileCurrent();
};

}

#endif
#endif
#endif

/* ERROR/WARNING messages:

E: No dump custom arguments specified

The dump custom command requires that atom quantities be specified to
output to dump file.

E: Invalid attribute in dump custom command

Self-explantory.

E: Could not find dump custom compute ID

The compute ID needed by dump custom to compute a per-atom quantity
does not exist.

E: Could not find dump custom fix ID

Self-explanatory.

E: Dump custom and fix not computed at compatible times

The fix must produce per-atom quantities on timesteps that dump custom
needs them.

E: Could not find dump custom variable name

Self-explanatory.

E: Region ID for dump custom does not exist

Self-explanatory.

E: Threshhold for an atom property that isn't allocated

A dump threshhold has been requested on a quantity that is
not defined by the atom style used in this simulation.

E: Dumping an atom property that isn't allocated

The chosen atom style does not define the per-atom quantity being
dumped.

E: Dumping an atom quantity that isn't allocated

Only per-atom quantities that are defined for the atom style being
used are allowed.

E: Dump custom compute does not compute per-atom info

Self-explanatory.

E: Dump custom compute does not calculate per-atom vector

Self-explanatory.

E: Dump custom compute does not calculate per-atom array

Self-explanatory.

E: Dump custom compute vector is accessed out-of-range

Self-explanatory.

E: Dump custom fix does not compute per-atom info

Self-explanatory.

E: Dump custom fix does not compute per-atom vector

Self-explanatory.

E: Dump custom fix does not compute per-atom array

Self-explanatory.

E: Dump custom fix vector is accessed out-of-range

Self-explanatory.

E: Dump custom variable is not atom-style variable

Only atom-style variables generate per-atom quantities, needed for
dump output.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Dump_modify region ID does not exist

Self-explanatory.

E: Dump modify element names do not match atom types

Number of element names must equal number of atom types.

E: Invalid attribute in dump modify command

Self-explantory.

E: Could not find dump modify compute ID

Self-explanatory.

E: Dump modify compute ID does not compute per-atom info

Self-explanatory.

E: Dump modify compute ID does not compute per-atom vector

Self-explanatory.

E: Dump modify compute ID does not compute per-atom array

Self-explanatory.

E: Dump modify compute ID vector is not large enough

Self-explanatory.

E: Could not find dump modify fix ID

Self-explanatory.

E: Dump modify fix ID does not compute per-atom info

Self-explanatory.

E: Dump modify fix ID does not compute per-atom vector

Self-explanatory.

E: Dump modify fix ID does not compute per-atom array

Self-explanatory.

E: Dump modify fix ID vector is not large enough

Self-explanatory.

E: Could not find dump modify variable name

Self-explanatory.

E: Dump modify variable is not atom-style variable

Self-explanatory.

E: Invalid dump_modify threshhold operator

Operator keyword used for threshold specification in not recognized.

*/
