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

#ifdef LAMMPS_VTK

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "dump_local_gran_vtk.h"
#include "atom.h"
#include "force.h"
#include "domain.h"
#include "region.h"
#include "group.h"
#include "input.h"
#include "variable.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "compute_pair_gran_local.h"
#include "fix.h"
#include "memory.h"
#include "error.h"
#include "sort_buffer.h"
#include <vector>
#include <sstream>
#include <vtkVersion.h>
#ifndef VTK_MAJOR_VERSION
#include <vtkConfigure.h>
#endif
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLPPolyDataWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>

using namespace LAMMPS_NS;

enum{VTK,VTP,VTU,PVTP,PVTU};   // file formats

/* ---------------------------------------------------------------------- */

DumpLocalGranVTK::DumpLocalGranVTK(LAMMPS *lmp, int narg, char **arg) :
    Dump(lmp, narg, arg),
    DumpVTK(lmp)
{
  //if (narg == 5) error->all(FLERR,"No dump custom/vtk arguments specified");

  clearstep = 1;

  nevery = force->inumeric(FLERR,arg[3]);

  if(narg < 6)
    error->all(FLERR,"dump local/gran/vtk requires 6 arguments");

  label = NULL; 

  {
    // parallel vtp/vtu requires proc number to be preceded by underscore '_'
    multiname_ex = NULL;
    char *ptr = strchr(filename,'%');
    if (ptr) {
      multiname_ex = new char[strlen(filename) + 16];
      *ptr = '\0';
      sprintf(multiname_ex,"%s_%d%s",filename,me,ptr+1);
      *ptr = '%';
    }
  }

  vtk_file_format = VTK;

  char *suffix = filename + strlen(filename) - strlen(".vtp");
  if (suffix > filename && strcmp(suffix,".vtp") == 0) {
    if (multiproc) vtk_file_format = PVTP;
    else           vtk_file_format = VTP;
  } else if (suffix > filename && strcmp(suffix,".vtu") == 0) {
    if (multiproc) vtk_file_format = PVTU;
    else           vtk_file_format = VTU;
  }

  if (vtk_file_format == VTK) { // no multiproc support for legacy vtk format
    if (me != 0) filewriter = 0;
    fileproc = 0;
    multiproc = 0;
    nclusterprocs = nprocs;
  }

  filecurrent = NULL;
  parallelfilecurrent = NULL;

  dumpLocalGran = new DumpLocalGran(lmp, igroup, nclusterprocs, multiproc, nevery, filewriter, fileproc);
  dumpLocalGran->parse_parameters(narg-5, &arg[5], true);
}

/* ---------------------------------------------------------------------- */

DumpLocalGranVTK::~DumpLocalGranVTK()
{
  delete [] filecurrent;
  delete [] parallelfilecurrent;
  delete [] multiname_ex;
  delete [] label;
  delete dumpLocalGran;
}

/* ---------------------------------------------------------------------- */

void DumpLocalGranVTK::init_style()
{
  dumpLocalGran->init_style();

  // setup function ptrs

  header_choice = &DumpLocalGranVTK::header_vtk;

  if (vtk_file_format == VTP || vtk_file_format == PVTP)
    write_choice = &DumpLocalGranVTK::write_vtp;
  else if (vtk_file_format == VTU || vtk_file_format == PVTU)
    write_choice = &DumpLocalGranVTK::write_vtu;
  else
    write_choice = &DumpLocalGranVTK::write_vtk;
}

/* ---------------------------------------------------------------------- */

void DumpLocalGranVTK::write_header(bigint /*ndump*/)
{
}

/* ---------------------------------------------------------------------- */

void DumpLocalGranVTK::header_vtk(bigint)
{
}

/* ---------------------------------------------------------------------- */

int DumpLocalGranVTK::count()
{
  return 0;
}

/* ---------------------------------------------------------------------- */

void DumpLocalGranVTK::write()
{
  
  // nme = # of dump lines this proc contributes to dump
  mbSet = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  bool usePolyData = false;
  if (vtk_file_format == VTP || vtk_file_format == PVTP)
    usePolyData = true;
#ifndef UNSTRUCTURED_GRID_VTK
  if (vtk_file_format == VTK)
    usePolyData = true;
#endif
  dumpLocalGran->prepare_mbSet(mbSet, usePolyData);
  if (filewriter)
      write_data(0, NULL);

}

/* ---------------------------------------------------------------------- */

void DumpLocalGranVTK::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpLocalGranVTK::setFileCurrent() {
  delete [] filecurrent;
  filecurrent = NULL;

  char *filestar = filename;
  if (multiproc) {
    if (multiproc > 1) { // if dump_modify fileper or nfile was used
      delete [] multiname_ex;
      multiname_ex = NULL;
      char *ptr = strchr(filename,'%');
      if (ptr) {
        int id;
        if (me + nclusterprocs == nprocs) // last filewriter
          id = multiproc -1;
        else
          id = me/nclusterprocs;
        multiname_ex = new char[strlen(filename) + 16];
        *ptr = '\0';
        sprintf(multiname_ex,"%s_%d%s",filename,id,ptr+1);
        *ptr = '%';
      }
    } // else multiname_ex built in constructor is OK
    filestar = multiname_ex;
  }

  if (multifile == 0) {
    filecurrent = new char[strlen(filestar) + 1];
    strcpy(filecurrent, filestar);
  } else {
    filecurrent = new char[strlen(filestar) + 16];
    char *ptr = strchr(filestar,'*');
    *ptr = '\0';
    if (padflag == 0) {
      sprintf(filecurrent,"%s" BIGINT_FORMAT "%s",
              filestar,update->ntimestep,ptr+1);
    } else {
      char bif[8],pad[16];
      strcpy(bif,BIGINT_FORMAT);
      sprintf(pad,"%%s%%0%d%s%%s",padflag,&bif[1]);
      sprintf(filecurrent,pad,filestar,update->ntimestep,ptr+1);
    }
    *ptr = '*';
  }

  // filename of parallel file
  if (multiproc && me == 0) {
    delete [] parallelfilecurrent;
    parallelfilecurrent = NULL;

    // remove '%' character and add 'p' to file extension
    // -> string length stays the same
    char *ptr = strchr(filename,'%');
    filestar = new char[strlen(filename) + 1];
    *ptr = '\0';
    sprintf(filestar,"%s%s",filename,ptr+1);
    *ptr = '%';
    ptr = strrchr(filestar,'.');
    ptr++;
    *ptr++='p';
    *ptr++='v';
    *ptr++='t';
    *ptr++= (vtk_file_format == PVTP)?'p':'u';
    *ptr++= 0;

    if (multifile == 0) {
      parallelfilecurrent = new char[strlen(filestar) + 1];
      strcpy(parallelfilecurrent, filestar);
    } else {
      parallelfilecurrent = new char[strlen(filestar) + 16];
      char *ptr = strchr(filestar,'*');
      *ptr = '\0';
      if (padflag == 0) {
        sprintf(parallelfilecurrent,"%s" BIGINT_FORMAT "%s",
                filestar,update->ntimestep,ptr+1);
      } else {
        char bif[8],pad[16];
        strcpy(bif,BIGINT_FORMAT);
        sprintf(pad,"%%s%%0%d%s%%s",padflag,&bif[1]);
        sprintf(parallelfilecurrent,pad,filestar,update->ntimestep,ptr+1);
      }
      *ptr = '*';
    }
    delete [] filestar;
    filestar = NULL;
  }
}

/* ---------------------------------------------------------------------- */

void DumpLocalGranVTK::write_vtk(int n, double *mybuf)
{
  setFileCurrent();

  {
#ifdef UNSTRUCTURED_GRID_VTK

    vtkSmartPointer<vtkDataObject> unstructuredGrid = mbSet->GetBlock(0);
    vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();

#else
    vtkSmartPointer<vtkDataObject> polyData = mbSet->GetBlock(0);
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
#endif

    if(label) writer->SetHeader(label);
    else      writer->SetHeader("Generated by LIGGGHTS");

    setVtkWriterOptions(vtkDataWriter::SafeDownCast(writer));

#ifdef UNSTRUCTURED_GRID_VTK
  #if VTK_MAJOR_VERSION < 6
    writer->SetInput(unstructuredGrid);
  #else
    writer->SetInputData(unstructuredGrid);
  #endif
#else
  #if VTK_MAJOR_VERSION < 6
    writer->SetInput(polyData);
  #else
    writer->SetInputData(polyData);
  #endif
#endif
    writer->SetFileName(filecurrent);
    writer->Write();
  }
}

/* ---------------------------------------------------------------------- */

void DumpLocalGranVTK::write_vtp(int n, double *mybuf)
{
  setFileCurrent();

  {
    vtkSmartPointer<vtkDataObject> polyData = mbSet->GetBlock(0);

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    setVtkWriterOptions(vtkXMLWriter::SafeDownCast(writer));

#if VTK_MAJOR_VERSION < 6
    writer->SetInput(polyData);
#else
    writer->SetInputData(polyData);
#endif
    writer->SetFileName(filecurrent);
    writer->Write();

    if (me == 0) {
      if (multiproc) {
        vtkSmartPointer<vtkXMLPPolyDataWriter> pwriter = vtkSmartPointer<vtkXMLPPolyDataWriter>::New();
        pwriter->SetFileName(parallelfilecurrent);
        pwriter->SetNumberOfPieces((multiproc > 1)?multiproc:nprocs);
        setVtkWriterOptions(vtkXMLWriter::SafeDownCast(pwriter));

#if VTK_MAJOR_VERSION < 6
        pwriter->SetInput(polyData);
#else
        pwriter->SetInputData(polyData);
#endif
        pwriter->Write();
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpLocalGranVTK::write_vtu(int n, double *mybuf)
{
  setFileCurrent();

  {
    vtkSmartPointer<vtkDataObject> unstructuredGrid = mbSet->GetBlock(0);
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    setVtkWriterOptions(vtkXMLWriter::SafeDownCast(writer));

#if VTK_MAJOR_VERSION < 6
    writer->SetInput(unstructuredGrid);
#else
    writer->SetInputData(unstructuredGrid);
#endif
    writer->SetFileName(filecurrent);
    writer->Write();

    if (me == 0) {
      if (multiproc) {
        vtkSmartPointer<vtkXMLPUnstructuredGridWriter> pwriter = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
        setVtkWriterOptions(vtkXMLWriter::SafeDownCast(pwriter));
        pwriter->SetFileName(parallelfilecurrent);
        pwriter->SetNumberOfPieces((multiproc > 1)?multiproc:nprocs);

#if VTK_MAJOR_VERSION < 6
        pwriter->SetInput(unstructuredGrid);
#else
        pwriter->SetInputData(unstructuredGrid);
#endif
        pwriter->Write();
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

int DumpLocalGranVTK::modify_param(int narg, char **arg)
{
    if (strcmp(arg[0],"label") == 0)
    {
        if (narg < 2)
            error->all(FLERR,"Illegal dump_modify command [label]");
        delete [] label;
        int n = strlen(arg[1]) + 1;
        label = new char[n];
        strcpy(label,arg[1]);
        return 2;
    }

    const int mvtk = DumpVTK::modify_param(narg, arg);
    if (mvtk > 0)
        return mvtk;

    return dumpLocalGran->modify_param(narg, arg);
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory in buf, choose, variable arrays
------------------------------------------------------------------------- */

bigint DumpLocalGranVTK::memory_usage()
{
  bigint bytes = Dump::memory_usage();
  bytes += dumpLocalGran->memory_usage();
  return bytes;
}

#endif
