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
   Daniel Queteschiner, DCS Computing daniel.queteschiner@dcs-computing.com
   Alexander Podlozhnyuk, DCS Computing alexander.podlozhnyuk@dcs-computing.com
   Christoph Kloss, DCS Computing
------------------------------------------------------------------------- */

#ifdef LAMMPS_VTK

#include <cmath>
#include "math_extra_liggghts.h"
#include <stdlib.h>
#include <string.h>
#include "dump_custom_vtk.h"

#ifdef CONVEX_ACTIVE_FLAG
#include "atom_vec_convexhull.h"
#endif

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
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkStringArray.h>
#include <vtkPolyData.h>
#include <vtkTriangle.h>
#include <vtkRectilinearGrid.h>
#include <vtkHexahedron.h>
#include <vtkUnstructuredGrid.h>
#include <vtkMPIController.h>
#include <vtkMPI.h>
#include <vtkMPICommunicator.h>

// For compatibility with new VTK generic data arrays (VTK >= 7.0)
#ifdef vtkGenericDataArray_h
#define InsertNextTupleValue InsertNextTypedTuple
#endif

using namespace LAMMPS_NS;

// customize by
// * adding an enum constant (add vector components in consecutive order)
// * adding a pack_*(int) function for the value
// * adjusting parse_fields function to add the pack_* function to pack_choice
//   (in case of vectors, adjust identify_vectors as well)
// * adjusting thresh part in modify_param and count functions

enum{X,Y,Z, // required for vtk, must come first
     POINTS_CONVEXHULL, 
     ID,MOL,TYPE,ELEMENT,MASS,
     XS,YS,ZS,XSTRI,YSTRI,ZSTRI,XU,YU,ZU,XUTRI,YUTRI,ZUTRI,
     XSU,YSU,ZSU,XSUTRI,YSUTRI,ZSUTRI,
     IX,IY,IZ,
     VX,VY,VZ,FX,FY,FZ,
     Q, MUX,MUY,MUZ,MU,RADIUS,DIAMETER,
     OMEGAX,OMEGAY,OMEGAZ,ANGMOMX,ANGMOMY,ANGMOMZ,
     TQX,TQY,TQZ,SPIN,ERADIUS,ERVEL,ERFORCE,
     DENSITY, RHO, P, 
     VARIABLE,COMPUTE,FIX,
     SHAPEX, SHAPEY, SHAPEZ,
     QUAT1, QUAT2, QUAT3, QUAT4,
     EXTRA1, EXTRA2, TENSOR,
     BLOCKINESS1, BLOCKINESS2, 
     ATTRIBUTES}; // must come last
enum{LT,LE,GT,GE,EQ,NEQ};
enum{INT,DOUBLE,STRING,TENSOR_DOUBLE};    // same as in DumpCFG
enum{VTK,VTP,VTU,PVTP,PVTU}; // file formats

/* ---------------------------------------------------------------------- */

DumpCustomVTK::DumpCustomVTK(LAMMPS *lmp, int narg, char **arg) :
    Dump(lmp, narg, arg),
    DumpVTK(lmp),
    filecurrent(NULL),
    domainfilecurrent(NULL)
{
    if (narg == 5)
        error->all(FLERR,"No dump custom/vtk arguments specified");

    clearstep = 1;

    nevery = force->inumeric(FLERR,arg[3]);

    // size_one may be shrunk below if additional optional args exist

    size_one = nfield = narg - 5;

    // atom selection arrays

    label = NULL; 

    std::list<int> allowed_extensions;
    allowed_extensions.push_back(VTK_FILE_FORMATS::VTK);
    allowed_extensions.push_back(VTK_FILE_FORMATS::VTU);
    allowed_extensions.push_back(VTK_FILE_FORMATS::PVTU);
    allowed_extensions.push_back(VTK_FILE_FORMATS::VTP);
    allowed_extensions.push_back(VTK_FILE_FORMATS::PVTP);
    vtk_file_format = DumpVTK::identify_file_type(filename, allowed_extensions, style, multiproc, nclusterprocs, filewriter, fileproc, world, clustercomm);

    filecurrent = NULL;
    domainfilecurrent = NULL;

    dumpParticle = new DumpParticle(lmp, igroup, nclusterprocs, multiproc, nevery, filewriter, fileproc);
    dumpParticle->parse_parameters(narg-5, &arg[5], true);

    if (!vtkMultiProcessController::GetGlobalController())
    {
        vtkMPICommunicatorOpaqueComm vtkWorldOpaqueComm(&world);
        vtkMPICommunicator * vtkWorldComm = vtkMPICommunicator::New();
        vtkWorldComm->InitializeExternal(&vtkWorldOpaqueComm);
        vtkMPIController *vtkController = vtkMPIController::New();
        vtkController->SetCommunicator(vtkWorldComm);
        vtkMultiProcessController::SetGlobalController(vtkController);
    }
}

/* ---------------------------------------------------------------------- */

DumpCustomVTK::~DumpCustomVTK()
{
  delete [] filecurrent;
  delete [] domainfilecurrent;

  delete dumpParticle;
  delete [] label;
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::init_style()
{
  dumpParticle->init_style();

  // setup function ptrs

  header_choice = &DumpCustomVTK::header_vtk;

  if (vtk_file_format == VTK_FILE_FORMATS::VTP || vtk_file_format == VTK_FILE_FORMATS::PVTP)
    write_choice = &DumpCustomVTK::write_vtp;
  else if (vtk_file_format == VTK_FILE_FORMATS::VTU || vtk_file_format == VTK_FILE_FORMATS::PVTU)
    write_choice = &DumpCustomVTK::write_vtu;
  else
    write_choice = &DumpCustomVTK::write_vtk;
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::write_header(bigint /*ndump*/)
{
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::header_vtk(bigint)
{
}

/* ---------------------------------------------------------------------- */

int DumpCustomVTK::count()
{
    return dumpParticle->count();
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::write()
{
  nme = count();
  mbSet = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  bool usePolyData = false;
  if (vtk_file_format == VTK_FILE_FORMATS::VTP || vtk_file_format == VTK_FILE_FORMATS::PVTP)
    usePolyData = true;
#ifndef UNSTRUCTURED_GRID_VTK
  if (vtk_file_format == VTK_FILE_FORMATS::VTK)
    usePolyData = true;
#endif
  dumpParticle->prepare_mbSet(mbSet, usePolyData);
  if (filewriter)
    write_data(0, NULL);
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::setFileCurrent()
{
    DumpVTK::setFileCurrent(filecurrent, filename, multifile, padflag);

    // filename of domain box data file
    delete [] domainfilecurrent;
    domainfilecurrent = NULL;
    domainfilecurrent = new char[strlen(filecurrent) + 14];
    char *tmp = new char[strlen(filecurrent)+1];
    strcpy(tmp, filecurrent);
    char * bbInsert;
    if (multifile == 0)
        bbInsert = strrchr(filecurrent, '.');
    else
        bbInsert = filecurrent + (strchr(filename, '*') - filename);
    char *ptr = tmp + (bbInsert - filecurrent);
    *ptr = '\0';
    if (multifile == 0)
        sprintf(domainfilecurrent,"%s_boundingBox%s",tmp,bbInsert);
    else
        sprintf(domainfilecurrent,"%sboundingBox_%s",tmp,bbInsert);
    *ptr = '.';
    delete [] tmp;
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::write_domain_vtk()
{
    vtkSmartPointer<vtkDataObject> rgrid = mbSet->GetBlock(1);
    DumpVTK::write_vtk_rectilinear_grid(rgrid, vtk_file_format, domainfilecurrent, label);
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::write_domain_vtk_triclinic()
{
    vtkSmartPointer<vtkDataObject> hexahedronGrid = mbSet->GetBlock(1);
    DumpVTK::write_vtk_unstructured_grid(hexahedronGrid, vtk_file_format, domainfilecurrent, label);
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::write_domain_vtr()
{
    vtkSmartPointer<vtkDataObject> rgrid = mbSet->GetBlock(1);
    DumpVTK::write_vtr(rgrid, VTK_FILE_FORMATS::VTR, domainfilecurrent);
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::write_domain_vtu_triclinic()
{
    vtkSmartPointer<vtkDataObject> hexahedronGrid = mbSet->GetBlock(1);
    DumpVTK::write_vtu(hexahedronGrid, VTK_FILE_FORMATS::VTU, domainfilecurrent);
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::write_vtk(int n, double *mybuf)
{
    setFileCurrent();

#ifdef UNSTRUCTURED_GRID_VTK
    vtkSmartPointer<vtkDataObject> unstructuredGrid = mbSet->GetBlock(0);
    DumpVTK::write_vtk_unstructured_grid(unstructuredGrid, vtk_file_format, filecurrent, label);
#else
    vtkSmartPointer<vtkDataObject> polyData = mbSet->GetBlock(0);
    DumpVTK::write_vtk_poly(polyData, vtk_file_format, filecurrent, label);
#endif

    if (domain->triclinic == 0)
        write_domain_vtk();
    else
        write_domain_vtk_triclinic();
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::write_vtp(int n, double *mybuf)
{
    setFileCurrent();

    vtkSmartPointer<vtkDataObject> polyData = mbSet->GetBlock(0);

    DumpVTK::write_vtp(polyData, vtk_file_format, filecurrent);

    if (me == 0)
    {
        if (domain->triclinic == 0)
        {
            domainfilecurrent[strlen(domainfilecurrent)-1] = 'r'; // adjust filename extension
            write_domain_vtr();
        }
        else
        {
            domainfilecurrent[strlen(domainfilecurrent)-1] = 'u'; // adjust filename extension
            write_domain_vtu_triclinic();
        }
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::write_vtu(int n, double *mybuf)
{
    setFileCurrent();

    vtkSmartPointer<vtkDataObject> unstructuredGrid = mbSet->GetBlock(0);

    DumpVTK::write_vtu(unstructuredGrid, vtk_file_format, filecurrent);

    if (me == 0)
    {
        if (domain->triclinic == 0)
        {
            domainfilecurrent[strlen(domainfilecurrent)-1] = 'r'; // adjust filename extension
            write_domain_vtr();
        }
        else
            write_domain_vtu_triclinic();
    }
}

/* ---------------------------------------------------------------------- */

int DumpCustomVTK::modify_param(int narg, char **arg)
{
    if (strcmp(arg[0],"label") == 0) { 
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

    return dumpParticle->modify_param(narg, arg);
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint DumpCustomVTK::memory_usage()
{
  bigint bytes = Dump::memory_usage();
  bytes += dumpParticle->memory_usage();
  return bytes;
}

#endif
