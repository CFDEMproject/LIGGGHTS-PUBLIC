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
#include <vector>
#include <sstream>
#include <vtkVersion.h>
#ifndef VTK_MAJOR_VERSION
#include <vtkConfigure.h>
#endif
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkLine.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkStringArray.h>
#include <vtkPolyData.h>
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

using namespace LAMMPS_NS;

enum{INT,DOUBLE,STRING};       // same as in DumpCFG
enum{VTK,VTP,VTU,PVTP,PVTU};   // file formats
enum{X1,X2,V1,V2,ID1,ID2,F,FN,FT,TORQUE,TORQUEN,TORQUET,AREA,DELTA,HEAT}; // dumps positions, force, normal and tangential forces, torque, normal and tangential torque

#define INVOKED_VECTOR 2
#define INVOKED_ARRAY 4
#define INVOKED_PERATOM 8
#define INVOKED_LOCAL 16

/* ---------------------------------------------------------------------- */

DumpLocalGranVTK::DumpLocalGranVTK(LAMMPS *lmp, int narg, char **arg) :
  Dump(lmp, narg, arg)
{
  //if (narg == 5) error->all(FLERR,"No dump custom/vtk arguments specified");

  clearstep = 1;

  nevery = force->inumeric(FLERR,arg[3]);

  if(narg < 6)
    error->all(FLERR,"dump local/gran/vtk requires 6 arguments");

  Compute *comp = modify->find_compute_id(arg[5]);
  if(!comp || !dynamic_cast<ComputePairGranLocal*>(comp))
    error->all(FLERR,"dump local/gran/vtk requires a valid ID of a compute pair/gran/local to be provided");

  cpgl_ = static_cast<ComputePairGranLocal*>(comp);

  // error if compute does not write pos
  if(cpgl_->offset_x1() < 0 || cpgl_->offset_x2() < 0)
    error->all(FLERR,"dump local/gran/vtk requires a valid ID of a compute pair/gran/local that writes the positions");

  // do stuff which needs the cpgl_ ptr
  size_one = cpgl_->get_nvalues();

  pack_choice.clear();
  vtype.clear();
  name.clear();

  // fill data into containers
  define_properties();

  iregion = -1;
  idregion = NULL;

  // computes, fixes, variables which the dump accesses

  myarrays.clear();
  n_calls_ = 0;

  if (filewriter) reset_vtk_data_containers();

  // atom selection arrays

  maxlocal = 0;

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
}

/* ---------------------------------------------------------------------- */

DumpLocalGranVTK::~DumpLocalGranVTK()
{
  delete [] filecurrent;
  delete [] parallelfilecurrent;
  delete [] multiname_ex;

  delete [] idregion;
  delete [] label; 
}

/* ---------------------------------------------------------------------- */

void DumpLocalGranVTK::init_style()
{

  // setup function ptrs

  header_choice = &DumpLocalGranVTK::header_vtk;

  if (vtk_file_format == VTP || vtk_file_format == PVTP)
    write_choice = &DumpLocalGranVTK::write_vtp;
  else if (vtk_file_format == VTU || vtk_file_format == PVTU)
    write_choice = &DumpLocalGranVTK::write_vtu;
  else
    write_choice = &DumpLocalGranVTK::write_vtk;

  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for dump custom/vtk does not exist");
  }
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
  
  n_calls_ = 0;

  //TODO generalize
  //also check if have same length
  cpgl_->compute_local();
  cpgl_->invoked_flag |= INVOKED_LOCAL;

  nchoose = cpgl_->get_ncount();

  return nchoose;
}

/* ---------------------------------------------------------------------- */

void DumpLocalGranVTK::write()
{
  
  // nme = # of dump lines this proc contributes to dump

  nme = count();

  // ntotal = total # of dump lines in snapshot
  // nmax = max # of dump lines on any proc

  bigint bnme = nme;
  MPI_Allreduce(&bnme,&ntotal,1,MPI_LMP_BIGINT,MPI_SUM,world);

  int nmax;
  if (multiproc != nprocs) MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,world);
  else nmax = nme;

  // write timestep header
  // for multiproc,
  //   nheader = # of lines in this file via Allreduce on clustercomm

  bigint nheader = ntotal;
  if (multiproc)
    MPI_Allreduce(&bnme,&nheader,1,MPI_LMP_BIGINT,MPI_SUM,clustercomm);

  if (filewriter) write_header(nheader);

  // insure buf is sized for packing and communicating
  // use nmax to insure filewriter proc can receive info from others
  // limit nmax*size_one to int since used as arg in MPI calls

  if (nmax > maxbuf) {
    if ((bigint) nmax * size_one > MAXSMALLINT)
      error->all(FLERR,"Too much per-proc info for dump");
    maxbuf = nmax;
    memory->destroy(buf);
    memory->create(buf,maxbuf*size_one,"dump:buf");
    
  }

  // insure ids buffer is sized for sorting

  if (sort_flag && sortcol == 0 && nmax > maxids) {
    maxids = nmax;
    memory->destroy(ids);
    memory->create(ids,maxids,"dump:ids");
  }

  // pack my data into buf
  // if sorting on IDs also request ID list from pack()
  // sort buf as needed

  if (sort_flag && sortcol == 0) pack(ids);
  else pack(NULL);
  if (sort_flag) sort();

  // filewriter = 1 = this proc writes to file
  //   ping each proc in my cluster, receive its data, write data to file
  // else wait for ping from fileproc, send my data to fileproc

  int tmp,nlines;
  MPI_Status status;
  MPI_Request request;

  // comm and output buf of doubles

  if (filewriter) {
    for (int iproc = 0; iproc < nclusterprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(buf,maxbuf*size_one,MPI_DOUBLE,me+iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_DOUBLE,&nlines);
        nlines /= size_one;
      } else nlines = nme;

      write_data(nlines,buf);
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,&status);
    MPI_Rsend(buf,nme*size_one,MPI_DOUBLE,fileproc,0,world);
  }

}

/* ---------------------------------------------------------------------- */

void DumpLocalGranVTK::pack(int *ids)
{
  int n = 0;
  for (std::map<int,FnPtrPack>::iterator it = pack_choice.begin(); it != pack_choice.end(); ++it)
  {
      (this->*(it->second))(n);

      // increase n by length of data
      if(vector_set.find(it->first) != vector_set.end())
        n += 3;
      else
        n++;

  }

  // similar to dump local, IDs are not used here (since are atom IDS
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

void DumpLocalGranVTK::buf2arrays(int n, double *mybuf)
{
  
  for (int idata=0; idata < n; ++idata) {

    // stores the ID of newly added points
    vtkIdType pid[2];

    pid[0] = points->InsertNextPoint(mybuf[idata*size_one],mybuf[idata*size_one+1],mybuf[idata*size_one+2]);
    pid[1] = points->InsertNextPoint(mybuf[idata*size_one+3],mybuf[idata*size_one+4],mybuf[idata*size_one+5]);

    // define the line going from point pid[0] to pid[1]
    vtkSmartPointer<vtkLine> line0 = vtkSmartPointer<vtkLine>::New();
    line0->GetPointIds()->SetId(0,pid[0]);
    line0->GetPointIds()->SetId(1,pid[1]);

    lineCells->InsertNextCell(line0);

    int j = 6; // 0,1,2,3,4,5 = 2 * (x,y,z) handled just above
    for (std::map<int, vtkSmartPointer<vtkAbstractArray> >::iterator it=myarrays.begin(); it!=myarrays.end(); ++it) {

      vtkAbstractArray *paa = it->second;
      if (it->second->GetNumberOfComponents() == 3) {
        switch (vtype[it->first]) {
          case INT:
            {
              int iv3[3] = { static_cast<int>(mybuf[idata*size_one+j  ]),
                             static_cast<int>(mybuf[idata*size_one+j+1]),
                             static_cast<int>(mybuf[idata*size_one+j+2]) };
              vtkIntArray *pia = static_cast<vtkIntArray*>(paa);
              pia->InsertNextTupleValue(iv3);
              break;
            }
          case DOUBLE:
            {
              vtkDoubleArray *pda = static_cast<vtkDoubleArray*>(paa);
              pda->InsertNextTupleValue(&mybuf[idata*size_one+j]);
              break;
            }
        }
        j+=3;
      } else {
        switch (vtype[it->first]) {
          case INT:
            {
              vtkIntArray *pia = static_cast<vtkIntArray*>(paa);
              pia->InsertNextValue(mybuf[idata*size_one+j]);
              break;
            }
          case DOUBLE:
            {
              vtkDoubleArray *pda = static_cast<vtkDoubleArray*>(paa);
              pda->InsertNextValue(mybuf[idata*size_one+j]);
              break;
            }
        }
        ++j;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpLocalGranVTK::write_vtk(int n, double *mybuf)
{
  ++n_calls_;

  buf2arrays(n, mybuf);

  if (n_calls_ < nclusterprocs)
    return; // multiple processors but only proc 0 is a filewriter (-> nclusterprocs procs contribute to the filewriter's output data)

  setFileCurrent();

  {
#ifdef UNSTRUCTURED_GRID_VTK

    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructuredGrid->SetPoints(points);
    unstructuredGrid->SetCells(VTK_LINE, lineCells);

    for (std::map<int, vtkSmartPointer<vtkAbstractArray> >::iterator it=myarrays.begin(); it!=myarrays.end(); ++it) {
      unstructuredGrid->GetCellData()->AddArray(it->second);
    }

    vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();

#else
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetLines(lineCells);

    for (std::map<int, vtkSmartPointer<vtkAbstractArray> >::iterator it=myarrays.begin(); it!=myarrays.end(); ++it) {
      polyData->GetCellData()->AddArray(it->second);
    }

    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
#endif

    if(label) writer->SetHeader(label);
    else      writer->SetHeader("Generated by LIGGGHTS");

    if (binary) writer->SetFileTypeToBinary();
    else        writer->SetFileTypeToASCII();

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

  reset_vtk_data_containers();
}

/* ---------------------------------------------------------------------- */

void DumpLocalGranVTK::write_vtp(int n, double *mybuf)
{
  ++n_calls_;

  buf2arrays(n, mybuf);

  if (n_calls_ < nclusterprocs)
    return; // multiple processors but not all are filewriters (-> nclusterprocs procs contribute to the filewriter's output data)

  setFileCurrent();

  {
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();

    polyData->SetPoints(points);
    polyData->SetLines(lineCells);

    for (std::map<int, vtkSmartPointer<vtkAbstractArray> >::iterator it=myarrays.begin(); it!=myarrays.end(); ++it) {
      polyData->GetCellData()->AddArray(it->second);
    }

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    if (binary) writer->SetDataModeToBinary();
    else        writer->SetDataModeToAscii();

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
        if (binary) pwriter->SetDataModeToBinary();
        else        pwriter->SetDataModeToAscii();

#if VTK_MAJOR_VERSION < 6
        pwriter->SetInput(polyData);
#else
        pwriter->SetInputData(polyData);
#endif
        pwriter->Write();
      }
    }
  }

  reset_vtk_data_containers();
}

/* ---------------------------------------------------------------------- */

void DumpLocalGranVTK::write_vtu(int n, double *mybuf)
{
  ++n_calls_;

  buf2arrays(n, mybuf);

  if (n_calls_ < nclusterprocs)
    return; // multiple processors but not all are filewriters (-> nclusterprocs procs contribute to the filewriter's output data)

  setFileCurrent();

  {
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    unstructuredGrid->SetPoints(points);
    unstructuredGrid->SetCells(VTK_LINE, lineCells);

    for (std::map<int, vtkSmartPointer<vtkAbstractArray> >::iterator it=myarrays.begin(); it!=myarrays.end(); ++it) {
      unstructuredGrid->GetCellData()->AddArray(it->second);
    }

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    if (binary) writer->SetDataModeToBinary();
    else        writer->SetDataModeToAscii();

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
        pwriter->SetFileName(parallelfilecurrent);
        pwriter->SetNumberOfPieces((multiproc > 1)?multiproc:nprocs);
        if (binary) pwriter->SetDataModeToBinary();
        else        pwriter->SetDataModeToAscii();

#if VTK_MAJOR_VERSION < 6
        pwriter->SetInput(unstructuredGrid);
#else
        pwriter->SetInputData(unstructuredGrid);
#endif
        pwriter->Write();
      }
    }
  }

  reset_vtk_data_containers();
}

/* ---------------------------------------------------------------------- */

void DumpLocalGranVTK::reset_vtk_data_containers()
{
  points = vtkSmartPointer<vtkPoints>::New();
  lineCells = vtkSmartPointer<vtkCellArray>::New();

  std::map<int,int>::iterator it=vtype.begin();

  ++it;
  ++it;

  for (; it!=vtype.end(); ++it) {

    // part 1: add VTK array to myarrays
    switch(vtype[it->first]) {
      case INT:
        myarrays[it->first] = vtkSmartPointer<vtkIntArray>::New();
        break;
      case DOUBLE:
        myarrays[it->first] = vtkSmartPointer<vtkDoubleArray>::New();
        break;
      case STRING:
        myarrays[it->first] = vtkSmartPointer<vtkStringArray>::New();
        break;
    }

    // part 2: if vector, set length to 3; set name
    if (vector_set.find(it->first) != vector_set.end()) {
      myarrays[it->first]->SetNumberOfComponents(3);
      myarrays[it->first]->SetName(name[it->first].c_str());
    } else {
      myarrays[it->first]->SetName(name[it->first].c_str());
    }
  }
}

/* ---------------------------------------------------------------------- */

// customize here

void DumpLocalGranVTK::define_properties()
{
  pack_choice[X1] = &DumpLocalGranVTK::pack_x1;
  vtype[X1] = DOUBLE;
  name[X1] = "pos1";
  vector_set.insert(X1);

  pack_choice[X2] = &DumpLocalGranVTK::pack_x2;
  vtype[X2] = DOUBLE;
  name[X2] = "pos2";
  vector_set.insert(X2);

  if(cpgl_->offset_v1() > 0)
  {
      pack_choice[V1] = &DumpLocalGranVTK::pack_v1;
      vtype[V1] = DOUBLE;
      name[V1] = "vel1";
      vector_set.insert(V1);
  }

  if(cpgl_->offset_v2() > 0)
  {
      pack_choice[V2] = &DumpLocalGranVTK::pack_v2;
      vtype[V2] = DOUBLE;
      name[V2] = "vel2";
      vector_set.insert(V2);
  }

  if(cpgl_->offset_id1() > 0)
  {
      pack_choice[ID1] = &DumpLocalGranVTK::pack_id1;
      vtype[ID1] = DOUBLE;
      name[ID1] = "id1";
      //scalar
  }

  if(cpgl_->offset_id2() > 0)
  {
      pack_choice[ID2] = &DumpLocalGranVTK::pack_id2;
      vtype[ID2] = DOUBLE;
      name[ID2] = "id2";
      //scalar
  }

  if(cpgl_->offset_f() > 0)
  {
      pack_choice[F] = &DumpLocalGranVTK::pack_f;
      vtype[F] = DOUBLE;
      name[F] = "force";
      vector_set.insert(F);
  }

  if(cpgl_->offset_fn() > 0)
  {
      pack_choice[FN] = &DumpLocalGranVTK::pack_fn;
      vtype[FN] = DOUBLE;
      name[FN] = "force_normal";
      vector_set.insert(FN);
  }

  if(cpgl_->offset_ft() > 0)
  {
      pack_choice[FT] = &DumpLocalGranVTK::pack_ft;
      vtype[FT] = DOUBLE;
      name[FT] = "force_tangential";
      vector_set.insert(FT);
  }

  if(cpgl_->offset_torque() > 0)
  {
      pack_choice[TORQUE] = &DumpLocalGranVTK::pack_torque;
      vtype[TORQUE] = DOUBLE;
      name[TORQUE] = "torque";
      vector_set.insert(TORQUE);
  }

  if(cpgl_->offset_torquen() > 0)
  {
      pack_choice[TORQUEN] = &DumpLocalGranVTK::pack_torquen;
      vtype[TORQUEN] = DOUBLE;
      name[TORQUEN] = "torque_normal";
      vector_set.insert(TORQUEN);
  }

  if(cpgl_->offset_torquet() > 0)
  {
      pack_choice[TORQUET] = &DumpLocalGranVTK::pack_torquet;
      vtype[TORQUET] = DOUBLE;
      name[TORQUET] = "torque_tangential";
      vector_set.insert(TORQUET);
  }

  if(cpgl_->offset_area() > 0)
  {
      
      pack_choice[AREA] = &DumpLocalGranVTK::pack_area;
      vtype[AREA] = DOUBLE;
      name[AREA] = "contact_area";
      //scalar
  }

  if(cpgl_->offset_delta() > 0)
  {
      pack_choice[DELTA] = &DumpLocalGranVTK::pack_delta;
      vtype[DELTA] = DOUBLE;
      name[DELTA] = "delta";
      //scalar
  }

  if(cpgl_->offset_heat() > 0)
  {
      pack_choice[HEAT] = &DumpLocalGranVTK::pack_heat;
      vtype[HEAT] = DOUBLE;
      name[HEAT] = "heat_flux";
      //scalar
  }
}

/* ---------------------------------------------------------------------- */

int DumpLocalGranVTK::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"region") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    if (strcmp(arg[1],"none") == 0) iregion = -1;
    else {
      iregion = domain->find_region(arg[1]);
      if (iregion == -1)
        error->all(FLERR,"Dump_modify region ID does not exist");
      delete [] idregion;
      int n = strlen(arg[1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[1]);
    }
    return 2;
  }

  if (strcmp(arg[0],"label") == 0) { 
     if (narg < 2) error->all(FLERR,"Illegal dump_modify command [label]");
     delete [] label;
     int n = strlen(arg[1]) + 1;
     label = new char[n];
     strcpy(label,arg[1]);
     return 2;
   }

  if (strcmp(arg[0],"binary") == 0) {
     if (narg < 2) error->all(FLERR,"Illegal dump_modify command [binary]");
     if (strcmp(arg[1],"yes") == 0) binary = 1;
     else if (strcmp(arg[1],"no") == 0) binary = 0;
     else error->all(FLERR,"Illegal dump_modify command [binary]");
     return 2;
  }

  if (strcmp(arg[0],"element") == 0) {
    error->all(FLERR,"Dump local/gran/vtk does not support dump_modify 'element' ");
    return 0;
  }

  if (strcmp(arg[0],"thresh") == 0) {
    error->all(FLERR,"Dump local/gran/vtk does not support dump_modify 'thresh' ");
  }

  return 0;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory in buf, choose, variable arrays
------------------------------------------------------------------------- */

bigint DumpLocalGranVTK::memory_usage()
{
  bigint bytes = Dump::memory_usage();
  return bytes;
}

/* ----------------------------------------------------------------------
   pack properties
   TODO generalize function, and length of packing
------------------------------------------------------------------------- */

// customize here

void DumpLocalGranVTK::pack_x1(int n)
{
  int offset = cpgl_->offset_x1();

  for (int i = 0; i < nchoose; i++) {
    vectorCopy3D(&cpgl_->get_data()[i][offset],&buf[n]);
    n += size_one;
  }
}

void DumpLocalGranVTK::pack_x2(int n)
{
  int offset = cpgl_->offset_x2();

  for (int i = 0; i < nchoose; i++) {
    vectorCopy3D(&cpgl_->get_data()[i][offset],&buf[n]);
    n += size_one;
  }
}

void DumpLocalGranVTK::pack_v1(int n)
{
  int offset = cpgl_->offset_v1();

  for (int i = 0; i < nchoose; i++) {
    vectorCopy3D(&cpgl_->get_data()[i][offset],&buf[n]);
    n += size_one;
  }
}

void DumpLocalGranVTK::pack_v2(int n)
{
  int offset = cpgl_->offset_v2();

  for (int i = 0; i < nchoose; i++) {
    vectorCopy3D(&cpgl_->get_data()[i][offset],&buf[n]);
    n += size_one;
  }
}

void DumpLocalGranVTK::pack_id1(int n)
{
  int offset = cpgl_->offset_id1();

  for (int i = 0; i < nchoose; i++) {
    buf[n] = cpgl_->get_data()[i][offset];
    n += size_one;
  }
}

void DumpLocalGranVTK::pack_id2(int n)
{
  int offset = cpgl_->offset_id2();

  for (int i = 0; i < nchoose; i++) {
    buf[n] = cpgl_->get_data()[i][offset];
    n += size_one;
  }
}

void DumpLocalGranVTK::pack_f(int n)
{
  int offset = cpgl_->offset_f();

  for (int i = 0; i < nchoose; i++) {
    vectorCopy3D(&cpgl_->get_data()[i][offset],&buf[n]);
    n += size_one;
  }
}

void DumpLocalGranVTK::pack_fn(int n)
{
  int offset = cpgl_->offset_fn();

  for (int i = 0; i < nchoose; i++) {
    vectorCopy3D(&cpgl_->get_data()[i][offset],&buf[n]);
    n += size_one;
  }
}

void DumpLocalGranVTK::pack_ft(int n)
{
  int offset = cpgl_->offset_ft();

  for (int i = 0; i < nchoose; i++) {
    vectorCopy3D(&cpgl_->get_data()[i][offset],&buf[n]);
    n += size_one;
  }
}

void DumpLocalGranVTK::pack_torque(int n)
{
  int offset = cpgl_->offset_torque();

  for (int i = 0; i < nchoose; i++) {
    vectorCopy3D(&cpgl_->get_data()[i][offset],&buf[n]);
    n += size_one;
  }
}

void DumpLocalGranVTK::pack_torquen(int n)
{
  int offset = cpgl_->offset_torquen();

  for (int i = 0; i < nchoose; i++) {
    vectorCopy3D(&cpgl_->get_data()[i][offset],&buf[n]);
    n += size_one;
  }
}

void DumpLocalGranVTK::pack_torquet(int n)
{
  int offset = cpgl_->offset_torquet();

  for (int i = 0; i < nchoose; i++) {
    vectorCopy3D(&cpgl_->get_data()[i][offset],&buf[n]);
    n += size_one;
  }
}

void DumpLocalGranVTK::pack_area(int n)
{
  int offset = cpgl_->offset_area();

  for (int i = 0; i < nchoose; i++) {
    buf[n] = cpgl_->get_data()[i][offset];
    n += size_one;
  }
}

void DumpLocalGranVTK::pack_delta(int n)
{
  int offset = cpgl_->offset_delta();

  for (int i = 0; i < nchoose; i++) {
    buf[n] = cpgl_->get_data()[i][offset];
    n += size_one;
  }
}

void DumpLocalGranVTK::pack_heat(int n)
{
  int offset = cpgl_->offset_heat();

  for (int i = 0; i < nchoose; i++) {
    buf[n] = cpgl_->get_data()[i][offset];
    n += size_one;
  }
}

#endif
