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

#include <cmath>
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
#include <vtkMPIController.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

DumpLocalGranVTK::DumpLocalGranVTK(LAMMPS *lmp, int narg, char **arg) :
    Dump(lmp, narg, arg),
    DumpVTK(lmp),
    filecurrent(NULL)
{
    clearstep = 1;

    nevery = force->inumeric(FLERR,arg[3]);

    if(narg < 6)
        error->all(FLERR,"dump local/gran/vtk requires 6 arguments");

    label = NULL; 

    std::list<int> allowed_extensions;
    allowed_extensions.push_back(VTK_FILE_FORMATS::VTK);
    allowed_extensions.push_back(VTK_FILE_FORMATS::VTU);
    allowed_extensions.push_back(VTK_FILE_FORMATS::PVTU);
    allowed_extensions.push_back(VTK_FILE_FORMATS::VTP);
    allowed_extensions.push_back(VTK_FILE_FORMATS::PVTP);
    vtk_file_format = DumpVTK::identify_file_type(filename, allowed_extensions, style, multiproc, nclusterprocs, filewriter, fileproc, world, clustercomm);

    // ensure no old % format is used
    // this is set in dump.cpp
    if (multiproc)
        error->all(FLERR, "dump local/gran/vtk no longer allow parallel writing by setting the \% character. Instead use a filename with suffix .pvtX (X = {u, p}).");
    // find last dot
    char *suffix = strrchr(filename, '.');
    if (strlen(suffix) == 5)
    {
        multiproc = 1;
        nclusterprocs = 1;
        filewriter = 1;
        fileproc = comm->me;
        MPI_Comm_split(world,me,0,&clustercomm);
        if (strcmp(suffix, ".pvtp") == 0)
            vtk_file_format = VTK_FILE_FORMATS::PVTP;
        else if (strcmp(suffix, ".pvtu") == 0)
            vtk_file_format = VTK_FILE_FORMATS::PVTU;
        else
            error->all(FLERR, "dump local/gran/vtk only allows .pvtu or .pvtp for parallel writing");
    }
    else if (strlen(suffix) == 4)
    {
        if (strcmp(suffix, ".vtp") == 0)
            vtk_file_format = VTK_FILE_FORMATS::VTP;
        else if (strcmp(suffix, ".vtu") == 0)
            vtk_file_format = VTK_FILE_FORMATS::VTU;
        else if (strcmp(suffix, ".vtk") == 0)
            vtk_file_format = VTK_FILE_FORMATS::VTK;
        else
            error->warning(FLERR, "Unknown suffix in dump local/gran/vtk, expected .vt{u,p,k}. Writing legacy VTK file.");
    }
    else
        error->warning(FLERR, "Unknown suffix in dump local/gran/vtk, expected .vt{u,p,k}. Writing legacy VTK file.");

    filecurrent = NULL;

    dumpLocalGran = new DumpLocalGran(lmp, igroup, nclusterprocs, multiproc, nevery, filewriter, fileproc);
    dumpLocalGran->parse_parameters(narg-5, &arg[5], true);

    if (!vtkMultiProcessController::GetGlobalController())
    {
        vtkMPIController *vtkController = vtkMPIController::New();
        vtkController->Initialize();
        vtkMultiProcessController::SetGlobalController(vtkController);
    }
}

/* ---------------------------------------------------------------------- */

DumpLocalGranVTK::~DumpLocalGranVTK()
{
  delete [] filecurrent;
  delete [] label;
  delete dumpLocalGran;
}

/* ---------------------------------------------------------------------- */

void DumpLocalGranVTK::init_style()
{
  dumpLocalGran->init_style();

  // setup function ptrs

  header_choice = &DumpLocalGranVTK::header_vtk;

  if (vtk_file_format == VTK_FILE_FORMATS::VTP || vtk_file_format == VTK_FILE_FORMATS::PVTP)
    write_choice = &DumpLocalGranVTK::write_vtp;
  else if (vtk_file_format == VTK_FILE_FORMATS::VTU || vtk_file_format == VTK_FILE_FORMATS::PVTU)
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
  if (vtk_file_format == VTK_FILE_FORMATS::VTP || vtk_file_format == VTK_FILE_FORMATS::PVTP)
    usePolyData = true;
#ifndef UNSTRUCTURED_GRID_VTK
  if (vtk_file_format == VTK_FILE_FORMATS::VTK)
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

void DumpLocalGranVTK::setFileCurrent()
{
    DumpVTK::setFileCurrent(filecurrent, filename, multifile, padflag);
}

/* ---------------------------------------------------------------------- */

void DumpLocalGranVTK::write_vtk(int n, double *mybuf)
{
    setFileCurrent();
#ifdef UNSTRUCTURED_GRID_VTK
    vtkSmartPointer<vtkDataObject> unstructuredGrid = mbSet->GetBlock(0);
    DumpVTK::write_vtk_unstructured_grid(unstructuredGrid, vtk_file_format, filecurrent);
#else
    vtkSmartPointer<vtkDataObject> polyData = mbSet->GetBlock(0);
    DumpVTK::write_vtk_poly(polyData, vtk_file_format, filecurrent);
#endif
}

/* ---------------------------------------------------------------------- */

void DumpLocalGranVTK::write_vtp(int n, double *mybuf)
{
    setFileCurrent();

    vtkSmartPointer<vtkDataObject> polyData = mbSet->GetBlock(0);
    DumpVTK::write_vtp(polyData, vtk_file_format, filecurrent);
}

/* ---------------------------------------------------------------------- */

void DumpLocalGranVTK::write_vtu(int n, double *mybuf)
{
    setFileCurrent();

    vtkSmartPointer<vtkDataObject> unstructuredGrid = mbSet->GetBlock(0);
    DumpVTK::write_vtu(unstructuredGrid, vtk_file_format, filecurrent);
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
