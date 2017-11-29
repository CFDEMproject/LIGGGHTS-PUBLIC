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

    Arno Mayrhofer (DCS Computing GmbH, Linz)

    Copyright 2016-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifdef LAMMPS_VTK

#include <cmath>
#include <stdlib.h>
#include <string.h>
#include "dump_local_gran.h"
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
#include <vtkRectilinearGrid.h>
#include <vtkHexahedron.h>
#include <vtkUnstructuredGrid.h>
#include <vtkInformation.h>

#ifdef vtkGenericDataArray_h
#define InsertNextTupleValue InsertNextTypedTuple
#endif

using namespace LAMMPS_NS;

enum{INT,DOUBLE,STRING}; // same as in DumpCFG
enum{X1,X2,CP,V1,V2,ID1,ID2,F,FN,FT,TORQUE,TORQUEN,TORQUET,AREA,DELTA,HEAT,MSID1,MSID2}; // dumps positions, force, normal and tangential forces, torque, normal and tangential torque

/* ---------------------------------------------------------------------- */

DumpLocalGran::DumpLocalGran(LAMMPS *lmp, int _igroup, int _nclusterprocs, int _multiproc, int _nevery, int _filewriter, int _fileproc) :
    Pointers(lmp),
    nevery(_nevery),
    nclusterprocs(_nclusterprocs),
    multiproc(_multiproc),
    filewriter(_filewriter),
    fileproc(_fileproc),
    iregion(-1),
    idregion(NULL),
    igroup(_igroup),
    groupbit(group->bitmask[igroup]),
    nchoose(0),
    sortBuffer(NULL),
    maxbuf(0),
    buf(NULL),
    size_one(0),
    cpgl_(0),
    n_calls_(0)
{
    pack_choice.clear();
    vtype.clear();
    name.clear();
    myarrays.clear();
}

/* ---------------------------------------------------------------------- */

DumpLocalGran::~DumpLocalGran()
{
    delete [] idregion;
    if (sortBuffer)
        delete sortBuffer;
}

/* ---------------------------------------------------------------------- */

int DumpLocalGran::parse_parameters(const int narg, const char *const *const arg, bool optional_keyword, std::list<std::string> keyword_list)
{
    int iarg = 0;

    if (narg < 1)
        error->all(FLERR, "dump local/gran is missing arguments");

    if (strcmp(arg[iarg], "local_gran") != 0)
    {
        if (!optional_keyword)
            error->all(FLERR, "Missing keyword \"local_gran\" in dump local/gran");
    }
    else
        iarg++;

    if (narg < 1+iarg)
        error->all(FLERR, "dump local/gran is missing arguments");

    Compute *comp = modify->find_compute_id(arg[iarg++]);
    if(!comp || !dynamic_cast<ComputePairGranLocal*>(comp))
        error->all(FLERR,"dump local/gran requires a valid ID of a compute pair/gran/local to be provided");

    cpgl_ = static_cast<ComputePairGranLocal*>(comp);

    // error if compute does not write pos
    if(cpgl_->offset_x1() < 0 || cpgl_->offset_x2() < 0)
        error->all(FLERR,"dump local/gran requires a valid ID of a compute pair/gran/local that writes the positions");

    // do stuff which needs the cpgl_ ptr
    size_one = cpgl_->get_nvalues();

    // fill data into containers
    define_properties();

    if (filewriter) reset_vtk_data_containers();

    if (narg > iarg)
    {
        std::list<std::string>::iterator it;
        bool found_keyword = false;
        for (it = keyword_list.begin(); it != keyword_list.end(); it++)
        {
            if (it->compare(arg[iarg]) == 0)
            {
                found_keyword = true;
                break;
            }
        }
        if (!found_keyword)
            error->all(FLERR, "Could not find follow-up keyword in local/gran");
    }

    return iarg;
}

/* ---------------------------------------------------------------------- */

void DumpLocalGran::init_style()
{
    // set index and check validity of region

    if (iregion >= 0) {
        iregion = domain->find_region(idregion);
        if (iregion == -1)
            error->all(FLERR,"Region ID for dump custom/vtk does not exist");
    }
}

/* ---------------------------------------------------------------------- */

int DumpLocalGran::count()
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

void DumpLocalGran::prepare_mbSet(vtkSmartPointer<vtkMultiBlockDataSet> mbSet, bool use_poly_data)
{
    
    // nme = # of dump lines this proc contributes to dump

    int nme = count();

    // ntotal = total # of dump lines in snapshot
    // nmax = max # of dump lines on any proc

    bigint bnme = nme;
    bigint ntotal = nme;
    MPI_Allreduce(&bnme,&ntotal,1,MPI_LMP_BIGINT,MPI_SUM,world);

    int nmax;
    if (multiproc != comm->nprocs) MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,world);
    else nmax = nme;

    // ensure buf is sized for packing and communicating
    // use nmax to ensure filewriter proc can receive info from others
    // limit nmax*size_one to int since used as arg in MPI calls

    if (nmax > maxbuf) {
        if ((bigint) nmax * size_one > MAXSMALLINT)
            error->all(FLERR,"Too much per-proc info for dump");
        maxbuf = nmax;
        memory->destroy(buf);
        memory->create(buf,maxbuf*size_one,"dump:buf");
        
    }

    // ensure ids buffer is sized for sorting
    if (sortBuffer)
            sortBuffer->realloc_ids(nmax);

    // pack my data into buf
    // if sorting on IDs also request ID list from pack()
    // sort buf as needed

    if (sortBuffer)
            pack(sortBuffer->get_ids());
    else
            pack(NULL);
    if (sortBuffer)
            sortBuffer->sort(buf, nme, maxbuf, size_one, ntotal);

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
                MPI_Irecv(buf,maxbuf*size_one,MPI_DOUBLE,comm->me+iproc,0,world,&request);
                MPI_Send(&tmp,0,MPI_INT,comm->me+iproc,0,world);
                MPI_Wait(&request,&status);
                MPI_Get_count(&status,MPI_DOUBLE,&nlines);
                nlines /= size_one;
            } else nlines = nme;

            write_data(nlines, buf, mbSet, use_poly_data);
        }
    } else {
        MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,&status);
        MPI_Rsend(buf,nme*size_one,MPI_DOUBLE,fileproc,0,world);
    }

}

/* ---------------------------------------------------------------------- */

void DumpLocalGran::pack(int *ids)
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

void DumpLocalGran::buf2arrays(int n, double *mybuf)
{
    
    const bool have_cp = cpgl_->offset_contact_point() >= 0;

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
        if(have_cp)
        {
            vtkIdType pidCP[1];
            pidCP[0] = points->InsertNextPoint(mybuf[idata*size_one+6],mybuf[idata*size_one+7],mybuf[idata*size_one+8]);
            line0->GetPointIds()->SetId(0,pid[0]);
            line0->GetPointIds()->SetId(1,pidCP[0]);
            lineCells->InsertNextCell(line0);
            line0->GetPointIds()->SetId(0,pid[1]);
            line0->GetPointIds()->SetId(1,pidCP[0]);
            lineCells->InsertNextCell(line0);
        }

        int j = 6; // 0,1,2,3,4,5 = 2 * (x,y,z) handled just above
        if(have_cp)
            j += 3;
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
                            if (have_cp)
                            {
                                pia->InsertNextTupleValue(iv3);
                                pia->InsertNextTupleValue(iv3);
                            }
                            break;
                        }
                    case DOUBLE:
                        {
                            vtkDoubleArray *pda = static_cast<vtkDoubleArray*>(paa);
                            pda->InsertNextTupleValue(&mybuf[idata*size_one+j]);
                            if (have_cp)
                            {
                                pda->InsertNextTupleValue(&mybuf[idata*size_one+j]);
                                pda->InsertNextTupleValue(&mybuf[idata*size_one+j]);
                            }
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
                            if (have_cp)
                            {
                                pia->InsertNextValue(mybuf[idata*size_one+j]);
                                pia->InsertNextValue(mybuf[idata*size_one+j]);
                            }
                            break;
                        }
                    case DOUBLE:
                        {
                            vtkDoubleArray *pda = static_cast<vtkDoubleArray*>(paa);
                            pda->InsertNextValue(mybuf[idata*size_one+j]);
                            if (have_cp)
                            {
                                pda->InsertNextValue(mybuf[idata*size_one+j]);
                                pda->InsertNextValue(mybuf[idata*size_one+j]);
                            }
                            break;
                        }
                }
                ++j;
            }
        }
    }
}

/* ---------------------------------------------------------------------- */

void DumpLocalGran::write_data(int n, double *mybuf, vtkSmartPointer<vtkMultiBlockDataSet> mbSet, bool use_poly_data)
{
    ++n_calls_;
    int cur_block = mbSet->GetNumberOfBlocks();

    buf2arrays(n, mybuf);

    if (n_calls_ < nclusterprocs)
        return; // multiple processors but only proc 0 is a filewriter (-> nclusterprocs procs contribute to the filewriter's output data)

    if (!use_poly_data)
    {
        vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        unstructuredGrid->SetPoints(points);
        unstructuredGrid->SetCells(VTK_LINE, lineCells);

        for (std::map<int, vtkSmartPointer<vtkAbstractArray> >::iterator it=myarrays.begin(); it!=myarrays.end(); ++it) {
            unstructuredGrid->GetCellData()->AddArray(it->second);
        }
        mbSet->SetBlock(cur_block++, unstructuredGrid);

    }
    else
    {
        vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
        polyData->SetPoints(points);
        polyData->SetLines(lineCells);

        for (std::map<int, vtkSmartPointer<vtkAbstractArray> >::iterator it=myarrays.begin(); it!=myarrays.end(); ++it) {
            polyData->GetCellData()->AddArray(it->second);
        }
        mbSet->SetBlock(cur_block++, polyData);
    }
    std::string name = "local_gran_";
    name.append(cpgl_->id);
    mbSet->GetMetaData(cur_block-1)->Set(mbSet->NAME(), name.c_str());

    reset_vtk_data_containers();
}

/* ---------------------------------------------------------------------- */

void DumpLocalGran::reset_vtk_data_containers()
{
    points = vtkSmartPointer<vtkPoints>::New();
    lineCells = vtkSmartPointer<vtkCellArray>::New();

    std::map<int,int>::iterator it=vtype.begin();

    ++it;
    ++it;
    if(cpgl_->offset_contact_point() >= 0)
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

void DumpLocalGran::define_properties()
{
    pack_choice[X1] = &DumpLocalGran::pack_x1;
    vtype[X1] = DOUBLE;
    name[X1] = "pos1";
    vector_set.insert(X1);

    pack_choice[X2] = &DumpLocalGran::pack_x2;
    vtype[X2] = DOUBLE;
    name[X2] = "pos2";
    vector_set.insert(X2);

    if(cpgl_->offset_contact_point() >= 0)
    {
            pack_choice[CP] = &DumpLocalGran::pack_contact_point;
            vtype[CP] = DOUBLE;
            name[CP] = "contact_point";
            vector_set.insert(CP);
    }

    if(cpgl_->offset_v1() >= 0)
    {
            pack_choice[V1] = &DumpLocalGran::pack_v1;
            vtype[V1] = DOUBLE;
            name[V1] = "vel1";
            vector_set.insert(V1);
    }

    if(cpgl_->offset_v2() >= 0)
    {
            pack_choice[V2] = &DumpLocalGran::pack_v2;
            vtype[V2] = DOUBLE;
            name[V2] = "vel2";
            vector_set.insert(V2);
    }

    if(cpgl_->offset_id1() >= 0)
    {
            pack_choice[ID1] = &DumpLocalGran::pack_id1;
            vtype[ID1] = DOUBLE;
            name[ID1] = "id1";
            //scalar
    }

    if(cpgl_->offset_id2() >= 0)
    {
            pack_choice[ID2] = &DumpLocalGran::pack_id2;
            vtype[ID2] = DOUBLE;
            name[ID2] = "id2";
            //scalar
    }

    if(cpgl_->offset_f() >= 0)
    {
            pack_choice[F] = &DumpLocalGran::pack_f;
            vtype[F] = DOUBLE;
            name[F] = "force";
            vector_set.insert(F);
    }

    if(cpgl_->offset_fn() >= 0)
    {
            pack_choice[FN] = &DumpLocalGran::pack_fn;
            vtype[FN] = DOUBLE;
            name[FN] = "force_normal";
            vector_set.insert(FN);
    }

    if(cpgl_->offset_ft() >= 0)
    {
            pack_choice[FT] = &DumpLocalGran::pack_ft;
            vtype[FT] = DOUBLE;
            name[FT] = "force_tangential";
            vector_set.insert(FT);
    }

    if(cpgl_->offset_torque() >= 0)
    {
            pack_choice[TORQUE] = &DumpLocalGran::pack_torque;
            vtype[TORQUE] = DOUBLE;
            name[TORQUE] = "torque";
            vector_set.insert(TORQUE);
    }

    if(cpgl_->offset_torquen() >= 0)
    {
            pack_choice[TORQUEN] = &DumpLocalGran::pack_torquen;
            vtype[TORQUEN] = DOUBLE;
            name[TORQUEN] = "torque_normal";
            vector_set.insert(TORQUEN);
    }

    if(cpgl_->offset_torquet() >= 0)
    {
            pack_choice[TORQUET] = &DumpLocalGran::pack_torquet;
            vtype[TORQUET] = DOUBLE;
            name[TORQUET] = "torque_tangential";
            vector_set.insert(TORQUET);
    }

    if(cpgl_->offset_area() >= 0)
    {
            
            pack_choice[AREA] = &DumpLocalGran::pack_area;
            vtype[AREA] = DOUBLE;
            name[AREA] = "contact_area";
            //scalar
    }

    if(cpgl_->offset_delta() >= 0)
    {
            pack_choice[DELTA] = &DumpLocalGran::pack_delta;
            vtype[DELTA] = DOUBLE;
            name[DELTA] = "delta";
            //scalar
    }

    if(cpgl_->offset_heat() >= 0)
    {
            pack_choice[HEAT] = &DumpLocalGran::pack_heat;
            vtype[HEAT] = DOUBLE;
            name[HEAT] = "heat_flux";
            //scalar
    }

    if(cpgl_->offset_ms_id1() >= 0)
    {
            pack_choice[MSID1] = &DumpLocalGran::pack_ms_id1;
            vtype[MSID1] = DOUBLE;
            name[MSID1] = "ms_id1";
            //scalar
    }

    if(cpgl_->offset_ms_id2() >= 0)
    {
            pack_choice[MSID2] = &DumpLocalGran::pack_ms_id2;
            vtype[MSID2] = DOUBLE;
            name[MSID2] = "ms_id2";
            //scalar
    }
}

/* ---------------------------------------------------------------------- */

int DumpLocalGran::modify_param(int narg, char **arg)
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

    if (strcmp(arg[0],"element") == 0) {
        error->all(FLERR,"Dump local/gran does not support dump_modify 'element' ");
        return 0;
    }

    if (strcmp(arg[0],"thresh") == 0) {
        error->all(FLERR,"Dump local/gran does not support dump_modify 'thresh' ");
    }

    if (!sortBuffer)
        sortBuffer = new SortBuffer(lmp, false);

    return sortBuffer->modify_param(narg, arg);
}

/* ----------------------------------------------------------------------
     return # of bytes of allocated memory in buf, choose, variable arrays
------------------------------------------------------------------------- */

bigint DumpLocalGran::memory_usage()
{
    bigint bytes = memory->usage(buf,maxbuf*size_one);
    if (sortBuffer)
        bytes += sortBuffer->memory_usage(size_one);
    return bytes;
}

/* ----------------------------------------------------------------------
     pack properties
     TODO generalize function, and length of packing
------------------------------------------------------------------------- */

// customize here

void DumpLocalGran::pack_x1(int n)
{
    int offset = cpgl_->offset_x1();

    for (int i = 0; i < nchoose; i++) {
        vectorCopy3D(&cpgl_->get_data()[i][offset],&buf[n]);
        n += size_one;
    }
}

void DumpLocalGran::pack_x2(int n)
{
    int offset = cpgl_->offset_x2();

    for (int i = 0; i < nchoose; i++) {
        vectorCopy3D(&cpgl_->get_data()[i][offset],&buf[n]);
        n += size_one;
    }
}

void DumpLocalGran::pack_v1(int n)
{
    int offset = cpgl_->offset_v1();

    for (int i = 0; i < nchoose; i++) {
        vectorCopy3D(&cpgl_->get_data()[i][offset],&buf[n]);
        n += size_one;
    }
}

void DumpLocalGran::pack_v2(int n)
{
    int offset = cpgl_->offset_v2();

    for (int i = 0; i < nchoose; i++) {
        vectorCopy3D(&cpgl_->get_data()[i][offset],&buf[n]);
        n += size_one;
    }
}

void DumpLocalGran::pack_id1(int n)
{
    int offset = cpgl_->offset_id1();

    for (int i = 0; i < nchoose; i++) {
        buf[n] = cpgl_->get_data()[i][offset];
        n += size_one;
    }
}

void DumpLocalGran::pack_id2(int n)
{
    int offset = cpgl_->offset_id2();

    for (int i = 0; i < nchoose; i++) {
        buf[n] = cpgl_->get_data()[i][offset];
        n += size_one;
    }
}

void DumpLocalGran::pack_f(int n)
{
    int offset = cpgl_->offset_f();

    for (int i = 0; i < nchoose; i++) {
        vectorCopy3D(&cpgl_->get_data()[i][offset],&buf[n]);
        n += size_one;
    }
}

void DumpLocalGran::pack_fn(int n)
{
    int offset = cpgl_->offset_fn();

    for (int i = 0; i < nchoose; i++) {
        vectorCopy3D(&cpgl_->get_data()[i][offset],&buf[n]);
        n += size_one;
    }
}

void DumpLocalGran::pack_ft(int n)
{
    int offset = cpgl_->offset_ft();

    for (int i = 0; i < nchoose; i++) {
        vectorCopy3D(&cpgl_->get_data()[i][offset],&buf[n]);
        n += size_one;
    }
}

void DumpLocalGran::pack_torque(int n)
{
    int offset = cpgl_->offset_torque();

    for (int i = 0; i < nchoose; i++) {
        vectorCopy3D(&cpgl_->get_data()[i][offset],&buf[n]);
        n += size_one;
    }
}

void DumpLocalGran::pack_torquen(int n)
{
    int offset = cpgl_->offset_torquen();

    for (int i = 0; i < nchoose; i++) {
        vectorCopy3D(&cpgl_->get_data()[i][offset],&buf[n]);
        n += size_one;
    }
}

void DumpLocalGran::pack_torquet(int n)
{
    int offset = cpgl_->offset_torquet();

    for (int i = 0; i < nchoose; i++) {
        vectorCopy3D(&cpgl_->get_data()[i][offset],&buf[n]);
        n += size_one;
    }
}

void DumpLocalGran::pack_area(int n)
{
    int offset = cpgl_->offset_area();

    for (int i = 0; i < nchoose; i++) {
        buf[n] = cpgl_->get_data()[i][offset];
        n += size_one;
    }
}

void DumpLocalGran::pack_delta(int n)
{
    int offset = cpgl_->offset_delta();

    for (int i = 0; i < nchoose; i++) {
        buf[n] = cpgl_->get_data()[i][offset];
        n += size_one;
    }
}

void DumpLocalGran::pack_heat(int n)
{
    int offset = cpgl_->offset_heat();

    for (int i = 0; i < nchoose; i++) {
        buf[n] = cpgl_->get_data()[i][offset];
        n += size_one;
    }
}
void DumpLocalGran::pack_contact_point(int n)
{
    int offset = cpgl_->offset_contact_point();

    for (int i = 0; i < nchoose; i++) {
      vectorCopy3D(&cpgl_->get_data()[i][offset],&buf[n]);
        n += size_one;
    }
}

void DumpLocalGran::pack_ms_id1(int n)
{
    int offset = cpgl_->offset_ms_id1();

    for (int i = 0; i < nchoose; i++) {
        buf[n] = cpgl_->get_data()[i][offset];
        n += size_one;
    }
}

void DumpLocalGran::pack_ms_id2(int n)
{
    int offset = cpgl_->offset_ms_id2();

    for (int i = 0; i < nchoose; i++) {
        buf[n] = cpgl_->get_data()[i][offset];
        n += size_one;
    }
}

#endif
