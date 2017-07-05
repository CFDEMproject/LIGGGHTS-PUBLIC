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

#include <math.h>
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
#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLPPolyDataWriter.h>
#include <vtkRectilinearGrid.h>
#include <vtkRectilinearGridWriter.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkHexahedron.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkMPIController.h>

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
     ROUNDNESS1, ROUNDNESS2, 
     ATTRIBUTES}; // must come last
enum{LT,LE,GT,GE,EQ,NEQ};
enum{INT,DOUBLE,STRING,TENSOR_DOUBLE};    // same as in DumpCFG
enum{VTK,VTP,VTU,PVTP,PVTU}; // file formats

/* ---------------------------------------------------------------------- */

DumpCustomVTK::DumpCustomVTK(LAMMPS *lmp, int narg, char **arg) :
  Dump(lmp, narg, arg),
  DumpVTK(lmp)
{
  if (narg == 5) error->all(FLERR,"No dump custom/vtk arguments specified");

  clearstep = 1;

  nevery = force->inumeric(FLERR,arg[3]);

  // size_one may be shrunk below if additional optional args exist

  size_one = nfield = narg - 5;

  // atom selection arrays

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
  domainfilecurrent = NULL;
  parallelfilecurrent = NULL;

  dumpParticle = new DumpParticle(lmp, igroup, nclusterprocs, multiproc, nevery, filewriter, fileproc);
  dumpParticle->parse_parameters(narg-5, &arg[5], true);

    if (!vtkMultiProcessController::GetGlobalController())
    {
        vtkMPIController *vtkController = vtkMPIController::New();
        vtkController->Initialize();
        vtkMultiProcessController::SetGlobalController(vtkController);
    }
}

/* ---------------------------------------------------------------------- */

DumpCustomVTK::~DumpCustomVTK()
{
  delete [] filecurrent;
  delete [] domainfilecurrent;
  delete [] parallelfilecurrent;
  delete [] multiname_ex;

  delete dumpParticle;
  delete [] label; 
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::init_style()
{
  dumpParticle->init_style();

  // setup function ptrs

  header_choice = &DumpCustomVTK::header_vtk;

  if (vtk_file_format == VTP || vtk_file_format == PVTP)
    write_choice = &DumpCustomVTK::write_vtp;
  else if (vtk_file_format == VTU || vtk_file_format == PVTU)
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
  if (vtk_file_format == VTP || vtk_file_format == PVTP)
    usePolyData = true;
#ifndef UNSTRUCTURED_GRID_VTK
  if (vtk_file_format == VTK)
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

void DumpCustomVTK::setFileCurrent() {
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

  char * bbInsert = NULL;
  if (multifile == 0) {
    filecurrent = new char[strlen(filestar) + 1];
    strcpy(filecurrent, filestar);
    bbInsert = strrchr(filecurrent, '.');
  } else {
    filecurrent = new char[strlen(filestar) + 16];
    char *ptr = strchr(filestar,'*');
    *ptr = '\0';
    if (padflag == 0)
    {
      sprintf(filecurrent,"%s" BIGINT_FORMAT "%s",
              filestar,update->ntimestep,ptr+1);
    } else {
      char bif[8],pad[16];
      strcpy(bif,BIGINT_FORMAT);
      sprintf(pad,"%%s%%0%d%s%%s",padflag,&bif[1]);
      sprintf(filecurrent,pad,filestar,update->ntimestep,ptr+1);
    }
    bbInsert = filecurrent + (ptr - filestar);
    *ptr = '*';
  }

  // filename of domain box data file
  delete [] domainfilecurrent;
  domainfilecurrent = NULL;
  if (multiproc) {
    // remove '%' character
    char *ptr = strchr(filename,'%');
    domainfilecurrent = new char[strlen(filename)];
    *ptr = '\0';
    sprintf(domainfilecurrent,"%s%s",filename,ptr+1);
    *ptr = '%';

    if (multifile == 0)
    {
      // insert "boundingBox_" string
      ptr = strrchr(domainfilecurrent,'.');
      filestar = new char[strlen(domainfilecurrent)+15];
      *ptr = '\0';
      sprintf(filestar,"%s_boundingBox.%s",domainfilecurrent,ptr+1);
      delete [] domainfilecurrent;

      domainfilecurrent = new char[strlen(filestar) + 1];
      strcpy(domainfilecurrent, filestar);
    }
    else
    {
      // insert "boundingBox_" string
      ptr = strrchr(domainfilecurrent,'*');
      filestar = new char[strlen(domainfilecurrent)+15];
      *ptr = '\0';
      sprintf(filestar,"%sboundingBox_*%s",domainfilecurrent,ptr+1);
      delete [] domainfilecurrent;

      domainfilecurrent = new char[strlen(filestar) + 16];
      char *ptr = strchr(filestar,'*');
      *ptr = '\0';
      if (padflag == 0) {
        sprintf(domainfilecurrent,"%s" BIGINT_FORMAT "%s",
                filestar,update->ntimestep,ptr+1);
      } else {
        char bif[8],pad[16];
        strcpy(bif,BIGINT_FORMAT);
        sprintf(pad,"%%s%%0%d%s%%s",padflag,&bif[1]);
        sprintf(domainfilecurrent,pad,filestar,update->ntimestep,ptr+1);
      }
      *ptr = '*';
    }
    delete [] filestar;
    filestar = NULL;
  } else {
    domainfilecurrent = new char[strlen(filecurrent) + 14];
    char *tmp = new char[strlen(filecurrent)+1];
    strcpy(tmp, filecurrent);
    char *ptr = tmp + (bbInsert - filecurrent);
    *ptr = '\0';
    sprintf(domainfilecurrent,"%sboundingBox_%s",tmp,bbInsert);
    *ptr = '.';
  }

  // filename of parallel file
  if (multiproc) {
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

void DumpCustomVTK::write_domain_vtk()
{
  vtkSmartPointer<vtkDataObject> rgrid = mbSet->GetBlock(1);
  vtkSmartPointer<vtkRectilinearGridWriter> gwriter = vtkSmartPointer<vtkRectilinearGridWriter>::New();

  if(label) gwriter->SetHeader(label);
  else      gwriter->SetHeader("Generated by LIGGGHTS");

  setVtkWriterOptions(vtkDataWriter::SafeDownCast(gwriter));

#if VTK_MAJOR_VERSION < 6
  gwriter->SetInput(rgrid);
#else
  gwriter->SetInputData(rgrid);
#endif
  gwriter->SetFileName(domainfilecurrent);
  gwriter->Write();
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::write_domain_vtk_triclinic()
{
  vtkSmartPointer<vtkUnstructuredGridWriter> gwriter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();

  if(label) gwriter->SetHeader(label);
  else      gwriter->SetHeader("Generated by LIGGGHTS");

  setVtkWriterOptions(vtkDataWriter::SafeDownCast(gwriter));

  vtkSmartPointer<vtkDataObject> hexahedronGrid = mbSet->GetBlock(1);
#if VTK_MAJOR_VERSION < 6
  gwriter->SetInput(hexahedronGrid);
#else
  gwriter->SetInputData(hexahedronGrid);
#endif
  gwriter->SetFileName(domainfilecurrent);
  gwriter->Write();
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::write_domain_vtr()
{
    vtkSmartPointer<vtkXMLRectilinearGridWriter> gwriter = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();

    setVtkWriterOptions(vtkXMLWriter::SafeDownCast(gwriter));

    vtkSmartPointer<vtkDataObject> rgrid = mbSet->GetBlock(1);
#if VTK_MAJOR_VERSION < 6
    gwriter->SetInput(rgrid);
#else
    gwriter->SetInputData(rgrid);
#endif
    gwriter->SetFileName(domainfilecurrent);
    gwriter->Write();
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::write_domain_vtu_triclinic()
{
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> gwriter = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

  setVtkWriterOptions(vtkXMLWriter::SafeDownCast(gwriter));

  vtkSmartPointer<vtkDataObject> hexahedronGrid = mbSet->GetBlock(1);
#if VTK_MAJOR_VERSION < 6
  gwriter->SetInput(hexahedronGrid);
#else
  gwriter->SetInputData(hexahedronGrid);
#endif
  gwriter->SetFileName(domainfilecurrent);
  gwriter->Write();
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::write_vtk(int n, double *mybuf)
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

    if (domain->triclinic == 0)
      write_domain_vtk();
    else
      write_domain_vtk_triclinic();
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::write_vtp(int n, double *mybuf)
{
    setFileCurrent();

    vtkSmartPointer<vtkDataObject> polyData = mbSet->GetBlock(0);

    if (multiproc)
    {
        vtkSmartPointer<vtkXMLPPolyDataWriter> pwriter = vtkSmartPointer<vtkXMLPPolyDataWriter>::New();
        pwriter->SetFileName(parallelfilecurrent);

        setVtkWriterOptions(vtkXMLWriter::SafeDownCast(pwriter));

#if VTK_MAJOR_VERSION < 6
        pwriter->SetInput(polyData);
#else
        pwriter->SetInputData(polyData);
#endif

        pwriter->SetNumberOfPieces((multiproc > 1)?multiproc:nprocs);
        pwriter->SetStartPiece(me);
        pwriter->SetEndPiece(me);
        pwriter->Write();
    }
    else
    {
        vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

        setVtkWriterOptions(vtkXMLWriter::SafeDownCast(writer));

#if VTK_MAJOR_VERSION < 6
        writer->SetInput(polyData);
#else
        writer->SetInputData(polyData);
#endif

        writer->SetFileName(filecurrent);
        writer->Write();
    }

    if (me == 0) {
        if (domain->triclinic == 0) {
            domainfilecurrent[strlen(domainfilecurrent)-1] = 'r'; // adjust filename extension
            write_domain_vtr();
        } else {
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

    if (multiproc) {
        vtkSmartPointer<vtkXMLPUnstructuredGridWriter> pwriter = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
        pwriter->SetFileName(parallelfilecurrent);

        setVtkWriterOptions(vtkXMLWriter::SafeDownCast(pwriter));

#if VTK_MAJOR_VERSION < 6
        pwriter->SetInput(unstructuredGrid);
#else
        pwriter->SetInputData(unstructuredGrid);
#endif

        pwriter->SetNumberOfPieces((multiproc > 1)?multiproc:nprocs);
        pwriter->SetStartPiece(me);
        pwriter->SetEndPiece(me);
        pwriter->Write();
    }
    else
    {
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        setVtkWriterOptions(vtkXMLWriter::SafeDownCast(writer));

#if VTK_MAJOR_VERSION < 6
    writer->SetInput(unstructuredGrid);
#else
    writer->SetInputData(unstructuredGrid);
#endif

        writer->SetFileName(filecurrent);
        writer->Write();
    }

    if (me == 0) {
        if (domain->triclinic == 0) {
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
