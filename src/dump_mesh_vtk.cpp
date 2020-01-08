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
    (if not contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz

    Contributing author for interpolated output:
    Felix Kleinfeldt (OVGU Magdeburg)
------------------------------------------------------------------------- */

#ifdef LAMMPS_VTK

#include <string.h>
#include "dump_mesh_vtk.h"
#include "tri_mesh.h"
#include "domain.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "error.h"
#include "fix.h"
#include "fix_mesh_surface.h"
#include "modify.h"
#include "comm.h"
#include <stdint.h>
#include <vtkDataSet.h>
#include <vtkDataSetAttributes.h>
#include <vtkAbstractArray.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkType.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCellDataToPointData.h>
#include <vtkPointDataToCellData.h>
#include <vtkMPIController.h>
#include <vtkMPI.h>
#include <vtkMPICommunicator.h>

// For compatibility with new VTK generic data arrays (VTK >= 7.0)
#ifdef vtkGenericDataArray_h
#define InsertNextTupleValue InsertNextTypedTuple
#endif

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

DumpMeshVTK::DumpMeshVTK(LAMMPS *lmp, int narg, char **arg) :
    Dump(lmp, narg, arg),
    DumpVTK(lmp),
    filecurrent(NULL),
    dumpMesh_(NULL),
    vtk_file_format_(VTK_FILE_FORMATS::VTK),
    dataMode_(0)
{
    if (narg < 5)
        error->all(FLERR,"Illegal dump mesh/vtk command");

    //INFO: CURRENTLY ONLY PROC 0 writes

    format_default = NULL;

    int iarg = 5;

    char **dump_mesh_args(NULL);
    char **dump_mesh_properties(NULL);

    // +3 because we have the 'meshes' 'all' and 'mesh_properties' keywords that do not appear
    // by default in the mesh/vtk argument list but we need them for DumpMesh
    dump_mesh_args = new char*[narg+3];
    dump_mesh_properties = new char*[narg];
    int narg_dump_mesh = 0;
    int nproperties_dump_mesh = 0;

    dump_mesh_args[narg_dump_mesh++] = (char*)"meshes";

    bool hasargs = true;
    while (iarg < narg && hasargs)
    {
        hasargs = false;
        if(strcmp(arg[iarg],"output")==0)
        {
            if (iarg+2 > narg)
                error->all(FLERR,"Dump mesh/vtk: not enough arguments for 'interpolate'");
            if(strcmp(arg[iarg+1],"face")==0)
                dataMode_ = 0;
            else if(strcmp(arg[iarg+1],"interpolate")==0)
                dataMode_ = 1;
            else if (strcmp(arg[iarg+1], "original") == 0)
                dataMode_ = 2;
            else
                error->all(FLERR,"Dump mesh/vtk: wrong argument for 'output'");
            iarg += 2;
            hasargs = true;
        }
        else
        {
            // check if input is a mesh
            int ifix = modify->find_fix(arg[iarg]);
            FixMeshSurface *fms = (ifix < 0) ? 0 : dynamic_cast<FixMeshSurface*>(modify->fix[ifix]);
            if (fms)
                dump_mesh_args[narg_dump_mesh++] = arg[iarg];
            else
                dump_mesh_properties[nproperties_dump_mesh++] = arg[iarg];
            iarg++;
            hasargs = true;
        }
    }
    // if we did not find any mesh keywords dump all
    if (narg_dump_mesh == 1)
        dump_mesh_args[narg_dump_mesh++] = (char*)"all";
    // check if we dump at least one property
    if (nproperties_dump_mesh == 0)
        error->all(FLERR,"Dump mesh/vtk: No dump quantity selected");
    dump_mesh_args[narg_dump_mesh++] = (char*)"mesh_properties";
    // copy list of properties to args
    for (int i = 0; i < nproperties_dump_mesh; i++)
        dump_mesh_args[narg_dump_mesh++] = dump_mesh_properties[i];

    if (!vtkMultiProcessController::GetGlobalController())
    {
        vtkMPICommunicatorOpaqueComm vtkWorldOpaqueComm(&world);
        vtkMPICommunicator * vtkWorldComm = vtkMPICommunicator::New();
        vtkWorldComm->InitializeExternal(&vtkWorldOpaqueComm);
        vtkMPIController *vtkController = vtkMPIController::New();
        vtkController->SetCommunicator(vtkWorldComm);
        vtkMultiProcessController::SetGlobalController(vtkController);
    }
    vtkMPIController * controller = getLocalController();

    std::list<int> allowed_extensions;
    allowed_extensions.push_back(VTK_FILE_FORMATS::VTK);
    allowed_extensions.push_back(VTK_FILE_FORMATS::VTP);
    allowed_extensions.push_back(VTK_FILE_FORMATS::PVTP);
    vtk_file_format_= DumpVTK::identify_file_type(filename, allowed_extensions, style, multiproc, nclusterprocs, filewriter, fileproc, world, clustercomm);

    if (multiproc && dataMode_ != 2)
        error->all(FLERR, "Parallel writing does not allow interpolation on meshes. It is advised to do this in post-processing");

    dumpMesh_ = new DumpMesh(lmp, nclusterprocs, multiproc, filewriter, fileproc, controller);
    int ioptional = dumpMesh_->parse_parameters(narg_dump_mesh, dump_mesh_args);

    if (ioptional < narg_dump_mesh)
        error->all(FLERR,"Invalid attribute in dump mesh/vtm command");

    delete [] dump_mesh_args;
    delete [] dump_mesh_properties;
}

/* ---------------------------------------------------------------------- */

DumpMeshVTK::~DumpMeshVTK()
{
    if (filecurrent)
        delete [] filecurrent;
    delete dumpMesh_;
}

/* ---------------------------------------------------------------------- */

void DumpMeshVTK::init_style()
{
    size_one = dumpMesh_->init_style();
}

/* ---------------------------------------------------------------------- */

int DumpMeshVTK::modify_param(int narg, char **arg)
{
    const int mvtk = DumpVTK::modify_param(narg, arg);
    if (mvtk > 0)
        return mvtk;

    return 0;
}

/* ---------------------------------------------------------------------- */

int DumpMeshVTK::count()
{
    return 0;
}

/* ---------------------------------------------------------------------- */

void DumpMeshVTK::setFileCurrent()
{
    DumpVTK::setFileCurrent(filecurrent, filename, multifile, padflag);
}

/* ---------------------------------------------------------------------- */

void DumpMeshVTK::write()
{
    write_data(0, NULL);
}

/* ---------------------------------------------------------------------- */

void DumpMeshVTK::write_data(int n, double *mybuf)
{
    setFileCurrent();

    vtkSmartPointer<vtkMultiBlockDataSet> mbSet_ = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    dumpMesh_->prepare_mbSet(mbSet_);

    if (!filewriter)
        return;

    const unsigned int nblocks = mbSet_->GetNumberOfBlocks();
    vtkSmartPointer<vtkDataObject> polyData;
    if (nblocks == 1)
        polyData = mbSet_->GetBlock(0);
    else
    {
        std::list<std::string> point_attributes, cell_attributes;
        // VTK_INT or VTK_DOUBLE expected
        std::list<int> point_types, cell_types;
        // 1 = scalar, 3 = vector
        std::list<int> point_ncomp, cell_ncomp;
        int npoints = 0;
        for (unsigned int i = 0; i < nblocks; i++)
        {
            vtkSmartPointer<vtkDataObject> mesh = mbSet_->GetBlock(i);
            if (!mesh->IsA("vtkDataSet"))
                error->one(FLERR, "Internal error");
            npoints += static_cast<vtkDataSet*>(mesh.GetPointer())->GetNumberOfPoints();
            vtkSmartPointer<vtkDataSetAttributes> pointData = mesh->GetAttributes(vtkDataSet::POINT);
            vtkSmartPointer<vtkDataSetAttributes> cellData = mesh->GetAttributes(vtkDataSet::CELL);
            const int nPointArrays = pointData->GetNumberOfArrays();
            const int nCellArrays = cellData->GetNumberOfArrays();
            for (int j = 0; j < nPointArrays; j++)
            {
                vtkSmartPointer<vtkAbstractArray> data = pointData->GetAbstractArray(j);
                point_attributes.push_back(std::string(data->GetName()));
                point_types.push_back(data->GetDataType());
                point_ncomp.push_back(data->GetNumberOfComponents());
            }
            for (int j = 0; j < nCellArrays; j++)
            {
                vtkSmartPointer<vtkAbstractArray> data = cellData->GetAbstractArray(j);
                cell_attributes.push_back(std::string(data->GetName()));
                cell_types.push_back(data->GetDataType());
                cell_ncomp.push_back(data->GetNumberOfComponents());
            }
        }
        point_attributes.unique();
        cell_attributes.unique();

        // points & cells
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        points->SetNumberOfPoints(npoints);
        vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
        int offset = 0;
        for (unsigned int i = 0; i < nblocks; i++)
        {
            // this is allowed because we checked above
            vtkSmartPointer<vtkDataSet> mesh = static_cast<vtkDataSet*>(mbSet_->GetBlock(i));
            int niPoints = mesh->GetNumberOfPoints();
            for (int j = 0; j < niPoints; j++)
            {
                points->SetPoint(j+offset, mesh->GetPoint(j));
                if (j%3 == 2)
                {
                    cells->InsertNextCell(3);
                    cells->InsertCellPoint(j+offset-2);
                    cells->InsertCellPoint(j+offset-1);
                    cells->InsertCellPoint(j+offset-0);
                }
            }
            offset += niPoints;
        }

        // data
        std::list<vtkSmartPointer<vtkAbstractArray> > cell_data;
        std::list<vtkSmartPointer<vtkAbstractArray> > point_data;

        std::list<std::string>::iterator data = point_attributes.begin();
        std::list<int>::iterator type = point_types.begin();
        std::list<int>::iterator ncomp = point_ncomp.begin();
        for (; data != point_attributes.end() && type != point_types.end() && ncomp != point_ncomp.end(); data++, type++, ncomp++)
        {
            if (*type == VTK_INT)
            {
                point_data.push_back(static_cast<vtkSmartPointer<vtkAbstractArray> >(vtkSmartPointer<vtkIntArray>::New()));
                vtkIntArray *array = static_cast<vtkIntArray*>(point_data.back().GetPointer());
                array->SetName(data->c_str());
                if (*ncomp != 1)
                    error->one(FLERR, "Int arrays can only be scalars");
                for (unsigned int i = 0; i < nblocks; i++)
                {
                    // this is allowed because we checked above
                    vtkSmartPointer<vtkDataSet> mesh = static_cast<vtkDataSet*>(mbSet_->GetBlock(i));
                    vtkSmartPointer<vtkDataSetAttributes> pointData = mesh->GetAttributes(vtkDataSet::POINT);
                    int arrayId = -1;
                    for (int j = 0; j < pointData->GetNumberOfArrays(); j++)
                    {
                        if (data->compare(std::string(pointData->GetArrayName(j))) == 0)
                        {
                            arrayId = j;
                            break;
                        }
                    }
                    if (arrayId != -1)
                    {
                        vtkSmartPointer<vtkAbstractArray> old_array = pointData->GetArray(arrayId);
                        const int ntuples = old_array->GetNumberOfTuples();
                        for (int j = 0; j < ntuples; j++)
                            array->InsertNextTuple(j, old_array);
                    }
                    else
                    {
                        const int ntuples = mesh->GetNumberOfPoints();
                        for (int j = 0; j < ntuples; j++)
                            array->InsertNextValue(0);
                    }
                }
            }
            else if (*type == VTK_DOUBLE)
            {
                point_data.push_back(static_cast<vtkSmartPointer<vtkAbstractArray> >(vtkSmartPointer<vtkDoubleArray>::New()));
                vtkDoubleArray *array = static_cast<vtkDoubleArray*>(point_data.back().GetPointer());
                array->SetNumberOfComponents(*ncomp);
                array->SetName(data->c_str());
                for (unsigned int i = 0; i < nblocks; i++)
                {
                    // this is allowed because we checked above
                    vtkSmartPointer<vtkDataSet> mesh = static_cast<vtkDataSet*>(mbSet_->GetBlock(i));
                    vtkSmartPointer<vtkDataSetAttributes> pointData = mesh->GetAttributes(vtkDataSet::POINT);
                    int arrayId = -1;
                    for (int j = 0; j < pointData->GetNumberOfArrays(); j++)
                    {
                        if (data->compare(std::string(pointData->GetArrayName(j))) == 0)
                        {
                            arrayId = j;
                            break;
                        }
                    }
                    if (arrayId != -1)
                    {
                        vtkSmartPointer<vtkAbstractArray> old_array = pointData->GetArray(arrayId);
                        const int ntuples = old_array->GetNumberOfTuples();
                        for (int j = 0; j < ntuples; j++)
                            array->InsertNextTuple(j, old_array);
                    }
                    else
                    {
                        const int ntuples = mesh->GetNumberOfPoints();
                        double * dum = new double[*ncomp];
                        for (int j = 0; j < *ncomp; j++)
                            dum[j] = 0.0;
                        for (int j = 0; j < ntuples; j++)
                            array->InsertNextTuple(dum);
                        delete [] dum;
                    }
                }
            }
            else
                error->one(FLERR, "Only vtk types int and double allowd");
        }

        data = cell_attributes.begin();
        type = cell_types.begin();
        ncomp = cell_ncomp.begin();
        for (; data != cell_attributes.end() && type != cell_types.end() && ncomp != cell_ncomp.end(); data++, type++, ncomp++)
        {
            if (*type == VTK_INT)
            {
                cell_data.push_back(static_cast<vtkSmartPointer<vtkAbstractArray> >(vtkSmartPointer<vtkIntArray>::New()));
                vtkIntArray *array = static_cast<vtkIntArray*>(cell_data.back().GetPointer());
                array->SetName(data->c_str());
                if (*ncomp != 1)
                    error->one(FLERR, "Int arrays can only be scalars");
                for (unsigned int i = 0; i < nblocks; i++)
                {
                    // this is allowed because we checked above
                    vtkSmartPointer<vtkDataSet> mesh = static_cast<vtkDataSet*>(mbSet_->GetBlock(i));
                    vtkSmartPointer<vtkDataSetAttributes> cellData = mesh->GetAttributes(vtkDataSet::CELL);
                    int arrayId = -1;
                    for (int j = 0; j < cellData->GetNumberOfArrays(); j++)
                    {
                        if (data->compare(std::string(cellData->GetArrayName(j))) == 0)
                        {
                            arrayId = j;
                            break;
                        }
                    }
                    if (arrayId != -1)
                    {
                        vtkSmartPointer<vtkAbstractArray> old_array = cellData->GetArray(arrayId);
                        const int ntuples = old_array->GetNumberOfTuples();
                        for (int j = 0; j < ntuples; j++)
                            array->InsertNextTuple(j, old_array);
                    }
                    else
                    {
                        const int ntuples = mesh->GetNumberOfCells();
                        for (int j = 0; j < ntuples; j++)
                            array->InsertNextValue(0);
                    }
                }
            }
            else if (*type == VTK_DOUBLE)
            {
                cell_data.push_back(static_cast<vtkSmartPointer<vtkAbstractArray> >(vtkSmartPointer<vtkDoubleArray>::New()));
                vtkDoubleArray *array = static_cast<vtkDoubleArray*>(cell_data.back().GetPointer());
                array->SetNumberOfComponents(*ncomp);
                array->SetName(data->c_str());
                for (unsigned int i = 0; i < nblocks; i++)
                {
                    // this is allowed because we checked above
                    vtkSmartPointer<vtkDataSet> mesh = static_cast<vtkDataSet*>(mbSet_->GetBlock(i));
                    vtkSmartPointer<vtkDataSetAttributes> cellData = mesh->GetAttributes(vtkDataSet::CELL);
                    int arrayId = -1;
                    for (int j = 0; j < cellData->GetNumberOfArrays(); j++)
                    {
                        if (data->compare(std::string(cellData->GetArrayName(j))) == 0)
                        {
                            arrayId = j;
                            break;
                        }
                    }
                    if (arrayId != -1)
                    {
                        vtkSmartPointer<vtkAbstractArray> old_array = cellData->GetArray(arrayId);
                        const int ntuples = old_array->GetNumberOfTuples();
                        for (int j = 0; j < ntuples; j++)
                            array->InsertNextTuple(j, old_array);
                    }
                    else
                    {
                        const int ntuples = mesh->GetNumberOfCells();
                        double * dum = new double[*ncomp];
                        for (int j = 0; j < *ncomp; j++)
                            dum[j] = 0.0;
                        for (int j = 0; j < ntuples; j++)
                            array->InsertNextTuple(dum);
                        delete [] dum;
                    }
                }
            }
            else
                error->one(FLERR, "Only vtk types int and double allowd");
        }

        vtkSmartPointer<vtkPolyData> new_pData = vtkSmartPointer<vtkPolyData>::New();
        new_pData->SetPoints(points);
        new_pData->SetPolys(cells);

        std::list<vtkSmartPointer<vtkAbstractArray> >::iterator attributes;
        // cell attributes
        for (attributes = cell_data.begin(); attributes != cell_data.end(); attributes++)
            new_pData->GetCellData()->AddArray(*attributes);
        // point attributes
        for (attributes = point_data.begin(); attributes != point_data.end(); attributes++)
            new_pData->GetPointData()->AddArray(*attributes);
        polyData = new_pData;
    }

    // interpolate everything from points to cells
    if (dataMode_ == 0)
    {
        vtkSmartPointer<vtkPointDataToCellData> converter = vtkSmartPointer<vtkPointDataToCellData>::New();
        converter->SetPassPointData(0);

        #if VTK_MAJOR_VERSION <= 5
        converter->SetInputConnection(polyData->GetProducerPort());
        #else
        converter->SetInputData(polyData);
        #endif

        converter->Update();
        polyData = vtkPolyData::SafeDownCast(converter->GetOutput());
    }
    // interpolate everything from cells to points
    else if (dataMode_ == 1)
    {
        vtkSmartPointer<vtkCellDataToPointData> converter = vtkSmartPointer<vtkCellDataToPointData>::New();
        converter->SetPassCellData(0);

        #if VTK_MAJOR_VERSION <= 5
        converter->SetInputConnection(polyData->GetProducerPort());
        #else
        converter->SetInputData(polyData);
        #endif

        converter->Update();
        polyData = vtkPolyData::SafeDownCast(converter->GetOutput());
    }

    if (vtk_file_format_ == VTK_FILE_FORMATS::PVTP || vtk_file_format_ == VTK_FILE_FORMATS::VTP)
        DumpVTK::write_vtp(polyData, vtk_file_format_, filecurrent);
    else
        DumpVTK::write_vtk_poly(polyData, vtk_file_format_, filecurrent);
}

#endif
