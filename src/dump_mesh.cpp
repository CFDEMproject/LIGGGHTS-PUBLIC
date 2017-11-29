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
    Arno Mayrhofer (DCS Computing GmbH, Linz)

    Copyright 2016-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifdef LAMMPS_VTK

#include <string.h>
#include "dump_mesh.h"
#include "domain.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "error.h"
#include "fix.h"
#include "fix_mesh_surface.h"
#include "modify.h"
#include "comm.h"
#include "dump_vtk.h"
#include <stdint.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkStringArray.h>
#include <vtkPolyData.h>
#include <vtkInformation.h>
#include <vtkCellArray.h>
#include <vtkMPIController.h>
#include <vtkAppendPolyData.h>

// For compatibility with new VTK generic data arrays (VTK >= 7.0)
#ifdef vtkGenericDataArray_h
#define InsertNextTupleValue InsertNextTypedTuple
#endif

using namespace LAMMPS_NS;

enum
{
    DUMP_STRESS               = 1<<0,
    DUMP_STRESSCOMPONENTS     = 1<<1,
    DUMP_ID                   = 1<<2,
    DUMP_VEL                  = 1<<3,
    DUMP_WEAR                 = 1<<4,
    DUMP_TEMP                 = 1<<5,
    DUMP_OWNER                = 1<<6,
    DUMP_AREA                 = 1<<7,
    DUMP_AEDGES               = 1<<8,
    DUMP_ACORNERS             = 1<<9,
    DUMP_INDEX                = 1<<10,
    DUMP_NNEIGHS              = 1<<11,
    DUMP_MIN_ACTIVE_EDGE_DIST = 1<<12,
    DUMP_LIQUID_CONTENT       = 1<<13
};

/* ---------------------------------------------------------------------- */

DumpMesh::DumpMesh(LAMMPS *lmp, int _nclusterprocs, int _multiproc, int _filewriter, int _fileproc, vtkMPIController *controller) :
    Pointers(lmp),
    dump_what_(0),
    mbSet(NULL),
    localController(controller),
    n_calls_(0),
    nclusterprocs(_nclusterprocs),
    multiproc(_multiproc),
    filewriter(_filewriter),
    fileproc(_fileproc),
    n_all_(0),
    n_all_max_(0),
    buf_all_(NULL),
    sigma_n_(NULL),
    sigma_t_(NULL),
    wear_(NULL),
    v_node_(NULL),
    f_node_(NULL),
    T_(NULL),
    min_active_edge_dist_(NULL),
    liquid_content_(NULL),
    scalar_containers_(NULL),
    scalar_container_names_(NULL),
    n_scalar_containers_(0),
    vector_containers_(NULL),
    vector_container_names_(NULL),
    n_vector_containers_(0),
    container_args_(NULL),
    n_container_bases_(0)
{}

int DumpMesh::parse_parameters(const int narg, const char *const *const arg, std::list<std::string> keyword_list)
{
    int iarg = 0;

    if (narg < 2)
        error->all(FLERR, "Illegal dump mesh command");

    if (strcmp(arg[iarg++], "meshes") != 0)
        error->all(FLERR, "dump mesh expected keyword meshes");

    bool hasargs = true;
    bool first_arg = true;
    bool all_set = false;
    bool mesh_properties_follow = false;
    while (iarg < narg && hasargs)
    {
        hasargs = false;
        if (strcmp(arg[iarg], "all") == 0)
        {
            if (!first_arg)
                error->all(FLERR, "Keyword 'all' can only appear directly after meshes");
            first_arg = false;
            all_set = true;
            // dump all meshes so search for them
            int nMesh = modify->n_fixes_style("mesh/surface");

            for (int iMesh = 0; iMesh < nMesh; iMesh++)
            {
                FixMeshSurface* fms = static_cast<FixMeshSurface*>(modify->find_fix_style("mesh/surface",iMesh));
                meshList_.push_back(fms->triMesh());
                fms->dumpAdd();
            }

            if (meshList_.empty())
                error->warning(FLERR,"Dump mesh cannot find any fix of type 'mesh/surface' to dump");

            iarg++;
            hasargs = true;
        }
        else if (strcmp(arg[iarg], "mesh_properties") == 0)
        {
            first_arg = false;
            if (meshList_.empty())
                error->all(FLERR, "No mesh was added to the dump, use either keyword 'all' or one or several mesh ids.");
            mesh_properties_follow = true;
            iarg++;
            break;
        }
        else
        {
            // check if we found a keyword
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
            if (found_keyword)
                break;
            if (all_set)
                error->all(FLERR, "Keyword \"all\" was set in meshes and thus no additional mesh can be added");
            // If keyword is not all a mesh id needs to follow
            int ifix = modify->find_fix(arg[iarg++]);
            FixMeshSurface *fms = dynamic_cast<FixMeshSurface*>(modify->fix[ifix]);
            if (!fms)
                error->all(FLERR, "dump mesh could not find mesh id after 'meshes' keyword");
            meshList_.push_back(fms->triMesh());
            fms->dumpAdd();
            hasargs = true;
            first_arg = false;
        }
    }

    if (meshList_.empty())
        error->warning(FLERR,"Dump mesh included no mesh list");

    if (mesh_properties_follow)
    {
        bool hasargs = true;
        while (iarg < narg && hasargs)
        {
            hasargs = false;
            if(strcmp(arg[iarg],"stress")==0)
            {
                dump_what_ |= DUMP_STRESS;
                iarg++;
                hasargs = true;
            }
            else if(strcmp(arg[iarg],"stresscomponents")==0)
            {
                dump_what_ |= DUMP_STRESSCOMPONENTS;
                iarg++;
                hasargs = true;
            }
            else if(strcmp(arg[iarg],"id")==0)
            {
                dump_what_ |= DUMP_ID;
                iarg++;
                hasargs = true;
            }
            else if(strcmp(arg[iarg],"vel")==0)
            {
                dump_what_ |= DUMP_VEL;
                iarg++;
                hasargs = true;
            }
            else if(strcmp(arg[iarg],"wear")==0)
            {
                dump_what_ |= DUMP_WEAR;
                iarg++;
                hasargs = true;
            }
            else if(strcmp(arg[iarg],"temp")==0)
            {
                dump_what_ |= DUMP_TEMP;
                iarg++;
                hasargs = true;
            }
            else if(strcmp(arg[iarg],"owner")==0)
            {
                dump_what_ |= DUMP_OWNER;
                iarg++;
                hasargs = true;
            }
            else if(strcmp(arg[iarg],"area")==0)
            {
                dump_what_ |= DUMP_AREA;
                iarg++;
                hasargs = true;
            }
            else if(strcmp(arg[iarg],"aedges")==0)
            {
                dump_what_ |= DUMP_AEDGES;
                iarg++;
                hasargs = true;
            }
            else if(strcmp(arg[iarg],"acorners")==0)
            {
                dump_what_ |= DUMP_ACORNERS;
                iarg++;
                hasargs = true;
            }
            else if(strcmp(arg[iarg],"index")==0)
            {
                dump_what_ |= DUMP_INDEX;
                iarg++;
                hasargs = true;
            }
            else if(strcmp(arg[iarg],"nneighs")==0)
            {
                dump_what_ |= DUMP_NNEIGHS;
                iarg++;
                hasargs = true;
            }
            else if(strcmp(arg[iarg],"liquid")==0)
            {
                dump_what_ |= DUMP_LIQUID_CONTENT;
                iarg++;
                hasargs = true;
            }
            else if(strcmp(arg[iarg],"distaa")==0)
            {
                dump_what_ |= DUMP_MIN_ACTIVE_EDGE_DIST;
                iarg++;
                hasargs = true;
            }
            else
            {
                // check if we found a keyword
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
                if (found_keyword)
                    break;
                // if not we assume that it's some container
                n_container_bases_++;
                memory->grow(container_args_,n_container_bases_,100,"container_args_");
                strcpy(container_args_[n_container_bases_-1],arg[iarg++]);
                hasargs = true;
            }
        }
    }

    //INFO: CURRENTLY ONLY PROC 0 writes

    int nMesh = meshList_.size();
    // allocate arrays
    sigma_n_ = new ScalarContainer<double>*[nMesh];
    sigma_t_ = new ScalarContainer<double>*[nMesh];
    wear_ = new ScalarContainer<double>*[nMesh];
    v_node_ = new MultiVectorContainer<double,3,3>*[nMesh];
    f_node_ = new VectorContainer<double,3>*[nMesh];
    T_ = new ScalarContainer<double>*[nMesh];
    min_active_edge_dist_ = new ScalarContainer<double>*[nMesh];
    liquid_content_ = new ScalarContainer<double>*[nMesh];

    scalar_containers_ = new ScalarContainer<double>**[n_container_bases_];
    scalar_container_names_ = new char*[n_container_bases_];
    vector_containers_ = new VectorContainer<double,3>**[n_container_bases_];
    vector_container_names_ = new char*[n_container_bases_];
    for(int i = 0; i < n_container_bases_; i++)
    {
        scalar_containers_[i] = new ScalarContainer<double>*[nMesh];
        vector_containers_[i] = new VectorContainer<double,3>*[nMesh];
        scalar_container_names_[i] = new char[200];
        vector_container_names_[i] = new char[200];

        for(int j = 0; j < nMesh; j++)
        {
            scalar_containers_[i][j] = 0;
            vector_containers_[i][j] = 0;
        }
    }

    if(dump_what_ == 0 && n_container_bases_ == 0 && mesh_properties_follow)
        error->all(FLERR,"Dump mesh/vtk: No dump quantity selected");

    return iarg;
}

/* ---------------------------------------------------------------------- */

DumpMesh::~DumpMesh()
{
    std::list<TriMesh*>::iterator mesh;
    for (mesh = meshList_.begin(); mesh != meshList_.end(); mesh++)
    {
        static_cast<FixMeshSurface*>(modify->find_fix_id((*mesh)->mesh_id()))->dumpRemove();
    }

    meshList_.clear();
    memory->destroy(buf_all_);

    delete [] sigma_n_;
    delete [] sigma_t_;
    delete [] wear_;
    delete [] v_node_;
    delete [] f_node_;
    delete [] T_;
    delete [] min_active_edge_dist_;
    delete [] liquid_content_;

    for(int i = 0; i < n_container_bases_; i++)
    {
        delete [] scalar_containers_[i];
        delete [] vector_containers_[i];
        delete [] scalar_container_names_[i];
        delete [] vector_container_names_[i];
    }
    delete [] scalar_containers_;
    delete [] vector_containers_;
    delete [] scalar_container_names_;
    delete [] vector_container_names_;
}

/* ---------------------------------------------------------------------- */

int DumpMesh::init_style()
{
    // nodes
    int size_one = 9;

    // add sizes and get references to properties - some may stay NULL
    if(dump_what_ & DUMP_STRESS)
        size_one += 2;
    if(dump_what_ & DUMP_STRESSCOMPONENTS)
        size_one += 3;
    if(dump_what_ & DUMP_ID)
        size_one += 1;
    if(dump_what_ & DUMP_VEL)
        size_one += 3;
    if(dump_what_ & DUMP_WEAR)
        size_one += 1;
    if(dump_what_ & DUMP_TEMP)
      size_one += 1;
    if(dump_what_ & DUMP_OWNER)
      size_one += 1;
    if(dump_what_ & DUMP_AREA)
      size_one += 1;
    if(dump_what_ & DUMP_AEDGES)
      size_one += 1;
    if(dump_what_ & DUMP_ACORNERS)
      size_one += 1;
    if(dump_what_ & DUMP_INDEX)
      size_one += 1;
    if(dump_what_ & DUMP_NNEIGHS)
      size_one += 1;
    if(dump_what_ & DUMP_MIN_ACTIVE_EDGE_DIST)
      size_one += 1;
    if(dump_what_ & DUMP_LIQUID_CONTENT)
      size_one += 1;

    getGeneralRefs();

    size_one += n_scalar_containers_*1;
    size_one += n_vector_containers_*3;

    return size_one;
}

/* ---------------------------------------------------------------------- */

void DumpMesh::getGeneralRefs()
{
    
    n_scalar_containers_ = 0;
    n_vector_containers_ = 0;
    char cid[200];

    for(int ib = 0; ib < n_container_bases_; ib++)
    {
        bool found_scalar = false, found_vector = false;

        std::list<TriMesh*>::iterator mesh;
        int i = 0;
        for (mesh = meshList_.begin(); mesh != meshList_.end(); mesh++, i++)
        {
            if((*mesh)->prop().getElementProperty<ScalarContainer<double> >(container_args_[ib]))
            {
                found_scalar = true;
                scalar_containers_[n_scalar_containers_][i] = (*mesh)->prop().getElementProperty<ScalarContainer<double> >(container_args_[ib]);
                scalar_containers_[n_scalar_containers_][i]->id(cid);
                strcpy(scalar_container_names_[n_scalar_containers_],cid);
            }

            if((*mesh)->prop().getElementProperty<VectorContainer<double,3> >(container_args_[ib]))
            {
                found_vector = true;
                vector_containers_[n_vector_containers_][i] = (*mesh)->prop().getElementProperty<VectorContainer<double,3> >(container_args_[ib]);
                vector_containers_[n_vector_containers_][i]->id(cid);
                strcpy(vector_container_names_[n_vector_containers_],cid);
            }
        }

        if(found_scalar)
          n_scalar_containers_++;
        if(found_vector)
          n_vector_containers_++;

        if(!found_scalar && !found_vector)
          error->all(FLERR,"Illegal dump mesh/vtk command, unknown keyword or mesh");
    }
}

/* ---------------------------------------------------------------------- */

int DumpMesh::modify_param(int narg, char **arg)
{
    error->warning(FLERR,"dump_modify keyword is not supported by 'dump mesh' and is thus ignored");
    return 0;
}

/* ---------------------------------------------------------------------- */

int DumpMesh::count()
{
    int numTri = 0;

    n_calls_ = 0;
    n_all_ = 0;

    getRefs();

    std::list<TriMesh*>::iterator mesh;
    for (mesh = meshList_.begin(); mesh != meshList_.end(); mesh++)
    {
      if(!(*mesh)->isParallel() && 0 != comm->me)
          continue;
      numTri += (*mesh)->sizeLocal();
      
    }

    return numTri;
}

/* ---------------------------------------------------------------------- */

void DumpMesh::getRefs()
{
    // add sizes and get references to properties - some may stay NULL
    if(dump_what_ & DUMP_STRESS)
    {
        std::list<TriMesh*>::iterator mesh;
        int i = 0;
        for (mesh = meshList_.begin(); mesh != meshList_.end(); mesh++, i++)
        {
            sigma_n_[i] = (*mesh)->prop().getElementProperty<ScalarContainer<double> >("sigma_n");
            sigma_t_[i] = (*mesh)->prop().getElementProperty<ScalarContainer<double> >("sigma_t");
            //if(0 == comm->me && (!sigma_n_[i] || !sigma_t_[i]))
            //  error->warning(FLERR,"Trying to dump stress for mesh which does not calculate stress, will dump '0' instead");
        }
    }
    if(dump_what_ & DUMP_STRESSCOMPONENTS)
    {
        std::list<TriMesh*>::iterator mesh;
        int i = 0;
        for (mesh = meshList_.begin(); mesh != meshList_.end(); mesh++, i++)
        {
            f_node_[i] = (*mesh)->prop().getElementProperty<VectorContainer<double,3> >("f");
        }
    }
    if(dump_what_ & DUMP_VEL)
    {
        std::list<TriMesh*>::iterator mesh;
        int i = 0;
        for (mesh = meshList_.begin(); mesh != meshList_.end(); mesh++, i++)
        {
            v_node_[i] = (*mesh)->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v");
            
        }
    }
    if(dump_what_ & DUMP_WEAR)
    {
        std::list<TriMesh*>::iterator mesh;
        int i = 0;
        for (mesh = meshList_.begin(); mesh != meshList_.end(); mesh++, i++)
        {
            wear_[i] = (*mesh)->prop().getElementProperty<ScalarContainer<double> >("wear");
            //if(0 == comm->me && !wear_[i])
            //  error->warning(FLERR,"Trying to dump wear for mesh which does not calculate wear, will dump '0' instead");
        }
    }
    if(dump_what_ & DUMP_TEMP)
    {
        temp_per_element_.clear();
        std::list<TriMesh*>::iterator mesh;
        int i = 0;
        for (mesh = meshList_.begin(); mesh != meshList_.end(); mesh++, i++)
        {
            if( (*mesh)->prop().getElementProperty<ScalarContainer<double> >("Temp"))
            {
                T_[i] = (*mesh)->prop().getElementProperty<ScalarContainer<double> >("Temp");
                temp_per_element_.push_back(true);
            }
            else
            {
                temp_per_element_.push_back(true);
                T_[i] = (*mesh)->prop().getGlobalProperty<ScalarContainer<double> >("Temp");
            }
            //if(0 == comm->me && !T_[i])
            //  error->warning(FLERR,"Trying to dump temperature for mesh which does not calculate temperature, will dump '0' instead");
        }
    }
    if(dump_what_ & DUMP_MIN_ACTIVE_EDGE_DIST)
    {
        std::list<TriMesh*>::iterator mesh;
        int i = 0;
        for (mesh = meshList_.begin(); mesh != meshList_.end(); mesh++, i++)
        {
            min_active_edge_dist_[i] = (*mesh)->prop().getElementProperty<ScalarContainer<double> >("minActiveEdgeDist");
        }
    }
    if(dump_what_ & DUMP_LIQUID_CONTENT)
    {
        std::list<TriMesh*>::iterator mesh;
        int i = 0;
        for (mesh = meshList_.begin(); mesh != meshList_.end(); mesh++, i++)
        {
            liquid_content_[i] = (*mesh)->prop().getElementProperty<ScalarContainer<double> >("LiquidContent");
        }
    }
}
/* ---------------------------------------------------------------------- */

void DumpMesh::pack(int *ids)
{ }

/* ---------------------------------------------------------------------- */

void DumpMesh::write_data(int n, double *mybuf)
{ }

/* ---------------------------------------------------------------------- */

void DumpMesh::prepare_mbSet(vtkSmartPointer<vtkMultiBlockDataSet> mbSet)
{
    if (!localController)
        error->all(FLERR, "VTK MPI Controller not found");

    count();

    int cur_block = mbSet->GetNumberOfBlocks();

    std::list<TriMesh*>::iterator mesh;
    int iMesh = 0;
    for (mesh = meshList_.begin(); mesh != meshList_.end(); mesh++, iMesh++)
    {
        if(!(*mesh)->isParallel() && 0 != comm->me)
        {
            // empty mesh
            vtkSmartPointer<vtkPolyData> poly_data = vtkSmartPointer<vtkPolyData>::New();
            int status = localController->Send(poly_data.GetPointer(), fileproc, comm->me-fileproc);
            if (status != 1)
            {
                
                error->one(FLERR, "MPI send error");
            }

            mbSet->SetBlock(cur_block++, poly_data);
            mbSet->GetMetaData(cur_block-1)->Set(mbSet->NAME(), (*mesh)->mesh_id());
            continue;
        }

        int nlocal = (*mesh)->sizeLocal();

        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        
        points->SetNumberOfPoints(nlocal*3);
        vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
        for (int iTri = 0; iTri < nlocal; iTri++)
        {
            // Position information
            for(int j=0;j<3;j++)
            {
                double node[3];
                (*mesh)->node(iTri,j,node);
                points->SetPoint(3*iTri+j, node[0], node[1], node[2]);
            }
            // Construct triangles
            cells->InsertNextCell(3);
            cells->InsertCellPoint(iTri*3 + 0);
            cells->InsertCellPoint(iTri*3 + 1);
            cells->InsertCellPoint(iTri*3 + 2);
        }

        std::list<vtkSmartPointer<vtkAbstractArray> > cell_attributes;
        std::list<vtkSmartPointer<vtkAbstractArray> > point_attributes;
        // Cell attributes
        if (dump_what_ & DUMP_STRESS)
        {
            if (sigma_n_[iMesh])
            {
                cell_attributes.push_back(static_cast<vtkSmartPointer<vtkAbstractArray> >(vtkSmartPointer<vtkDoubleArray>::New()));
                vtkDoubleArray *sigma_n = static_cast<vtkDoubleArray*>(cell_attributes.back().GetPointer());
                sigma_n->SetName("normal_stress_average");
                for (int iTri = 0; iTri < nlocal; iTri++)
                    sigma_n->InsertNextValue(sigma_n_[iMesh]->get(iTri));
            }
            if (sigma_t_[iMesh])
            {
                cell_attributes.push_back(static_cast<vtkSmartPointer<vtkAbstractArray> >(vtkSmartPointer<vtkDoubleArray>::New()));
                vtkDoubleArray *sigma_t = static_cast<vtkDoubleArray*>(cell_attributes.back().GetPointer());
                sigma_t->SetName("shear_stress_average");
                for (int iTri = 0; iTri < nlocal; iTri++)
                    sigma_t->InsertNextValue(sigma_t_[iMesh]->get(iTri));
            }
        }
        if(dump_what_ & DUMP_STRESSCOMPONENTS)
        {
            if (f_node_[iMesh])
            {
                cell_attributes.push_back(static_cast<vtkSmartPointer<vtkAbstractArray> >(vtkSmartPointer<vtkDoubleArray>::New()));
                vtkDoubleArray *f_node = static_cast<vtkDoubleArray*>(cell_attributes.back().GetPointer());
                f_node->SetNumberOfComponents(3);
                f_node->SetName("stress");
                for (int iTri = 0; iTri < nlocal; iTri++)
                {
                    const double invArea = 1./(*mesh)->areaElem(iTri);
                    double f[3];
                    f_node_[iMesh]->get(iTri, f);
                    f[0] *= invArea;
                    f[1] *= invArea;
                    f[2] *= invArea;
                    f_node->InsertNextTupleValue(f);
                }
            }
        }
        if(dump_what_ & DUMP_ID)
        {
            cell_attributes.push_back(static_cast<vtkSmartPointer<vtkAbstractArray> >(vtkSmartPointer<vtkIntArray>::New()));
            vtkIntArray *id = static_cast<vtkIntArray*>(cell_attributes.back().GetPointer());
            id->SetName("meshid");
            for (int iTri = 0; iTri < nlocal; iTri++)
                id->InsertNextValue((*mesh)->id(iTri));
        }
        if(dump_what_ & DUMP_VEL)
        {
            if (v_node_[iMesh])
            {
                point_attributes.push_back(static_cast<vtkSmartPointer<vtkAbstractArray> >(vtkSmartPointer<vtkDoubleArray>::New()));
                vtkDoubleArray *v_node = static_cast<vtkDoubleArray*>(point_attributes.back().GetPointer());
                v_node->SetNumberOfComponents(3);
                v_node->SetName("v");
                double **v;
                memory->create<double>(v,3,3,"DumpMeshVTK:v");
                for (int iTri = 0; iTri < nlocal; iTri++)
                {
                    v_node_[iMesh]->get(iTri, v);
                    v_node->InsertNextTupleValue(v[0]);
                    v_node->InsertNextTupleValue(v[1]);
                    v_node->InsertNextTupleValue(v[2]);
                }
                memory->destroy<double>(v);
            }
        }
        if(dump_what_ & DUMP_WEAR)
        {
            if (wear_[iMesh])
            {
                cell_attributes.push_back(static_cast<vtkSmartPointer<vtkAbstractArray> >(vtkSmartPointer<vtkDoubleArray>::New()));
                vtkDoubleArray *wear = static_cast<vtkDoubleArray*>(cell_attributes.back().GetPointer());
                wear->SetName("wear");
                for (int iTri = 0; iTri < nlocal; iTri++)
                    wear->InsertNextValue(wear_[iMesh]->get(iTri));
            }
        }
        if(dump_what_ & DUMP_TEMP)
        {
            if (T_[iMesh])
            {
                cell_attributes.push_back(static_cast<vtkSmartPointer<vtkAbstractArray> >(vtkSmartPointer<vtkDoubleArray>::New()));
                vtkDoubleArray *T = static_cast<vtkDoubleArray*>(cell_attributes.back().GetPointer());
                T->SetName("Temp");
                if (temp_per_element_[iMesh])
                {
                    for (int iTri = 0; iTri < nlocal; iTri++)
                        T->InsertNextValue(T_[iMesh]->get(iTri));
                }
                else
                {
                    for (int iTri = 0; iTri < nlocal; iTri++)
                        T->InsertNextValue(T_[iMesh]->get(0));
                }
            }
        }
        if(dump_what_ & DUMP_OWNER)
        {
            cell_attributes.push_back(static_cast<vtkSmartPointer<vtkAbstractArray> >(vtkSmartPointer<vtkIntArray>::New()));
            vtkIntArray *id = static_cast<vtkIntArray*>(cell_attributes.back().GetPointer());
            id->SetName("owner");
            int me = comm->me;
            for (int iTri = 0; iTri < nlocal; iTri++)
                id->InsertNextValue(me);
        }
        if(dump_what_ & DUMP_AREA)
        {
            cell_attributes.push_back(static_cast<vtkSmartPointer<vtkAbstractArray> >(vtkSmartPointer<vtkDoubleArray>::New()));
            vtkDoubleArray *area = static_cast<vtkDoubleArray*>(cell_attributes.back().GetPointer());
            area->SetName("area");
            for (int iTri = 0; iTri < nlocal; iTri++)
                area->InsertNextValue((*mesh)->areaElem(iTri));
        }
        if(dump_what_ & DUMP_AEDGES)
        {
            cell_attributes.push_back(static_cast<vtkSmartPointer<vtkAbstractArray> >(vtkSmartPointer<vtkIntArray>::New()));
            vtkIntArray *aedges = static_cast<vtkIntArray*>(cell_attributes.back().GetPointer());
            aedges->SetName("active_edges");
            for (int iTri = 0; iTri < nlocal; iTri++)
                aedges->InsertNextValue((*mesh)->n_active_edges(iTri));
        }
        if(dump_what_ & DUMP_ACORNERS)
        {
            cell_attributes.push_back(static_cast<vtkSmartPointer<vtkAbstractArray> >(vtkSmartPointer<vtkIntArray>::New()));
            vtkIntArray *acorners = static_cast<vtkIntArray*>(cell_attributes.back().GetPointer());
            acorners->SetName("active_corners");
            for (int iTri = 0; iTri < nlocal; iTri++)
                acorners->InsertNextValue((*mesh)->n_active_corners(iTri));
        }
        if(dump_what_ & DUMP_INDEX)
        {
            cell_attributes.push_back(static_cast<vtkSmartPointer<vtkAbstractArray> >(vtkSmartPointer<vtkIntArray>::New()));
            vtkIntArray *index = static_cast<vtkIntArray*>(cell_attributes.back().GetPointer());
            index->SetName("index");
            for (int iTri = 0; iTri < nlocal; iTri++)
                index->InsertNextValue(iTri);
        }
        if(dump_what_ & DUMP_NNEIGHS)
        {
            cell_attributes.push_back(static_cast<vtkSmartPointer<vtkAbstractArray> >(vtkSmartPointer<vtkIntArray>::New()));
            vtkIntArray *nneighs = static_cast<vtkIntArray*>(cell_attributes.back().GetPointer());
            nneighs->SetName("nneighs");
            for (int iTri = 0; iTri < nlocal; iTri++)
                nneighs->InsertNextValue((*mesh)->nNeighs(iTri));
        }
        if(dump_what_ & DUMP_MIN_ACTIVE_EDGE_DIST)
        {
            if (min_active_edge_dist_[iMesh])
            {
                cell_attributes.push_back(static_cast<vtkSmartPointer<vtkAbstractArray> >(vtkSmartPointer<vtkDoubleArray>::New()));
                vtkDoubleArray *min_active_edge_dist = static_cast<vtkDoubleArray*>(cell_attributes.back().GetPointer());
                min_active_edge_dist->SetName("min_active_edge_dist");
                for (int iTri = 0; iTri < nlocal; iTri++)
                    min_active_edge_dist->InsertNextValue(min_active_edge_dist_[iMesh]->get(iTri));
            }
        }
        if(dump_what_ & DUMP_LIQUID_CONTENT)
        {
            if (liquid_content_[iMesh])
            {
                cell_attributes.push_back(static_cast<vtkSmartPointer<vtkAbstractArray> >(vtkSmartPointer<vtkDoubleArray>::New()));
                vtkDoubleArray *liquid_content = static_cast<vtkDoubleArray*>(cell_attributes.back().GetPointer());
                liquid_content->SetName("liquid_content");
                for (int iTri = 0; iTri < nlocal; iTri++)
                    liquid_content->InsertNextValue(liquid_content_[iMesh]->get(iTri));
            }
        }
        // scalar containers
        ScalarContainer<double> *scalar_cont_;
        for(int ib = 0; ib < n_scalar_containers_; ib++)
        {
           scalar_cont_ = scalar_containers_[ib][iMesh];
           if (scalar_cont_)
           {
                cell_attributes.push_back(static_cast<vtkSmartPointer<vtkAbstractArray> >(vtkSmartPointer<vtkDoubleArray>::New()));
                vtkDoubleArray *scalar_cont = static_cast<vtkDoubleArray*>(cell_attributes.back().GetPointer());
                scalar_cont->SetName(scalar_container_names_[ib]);
                for (int iTri = 0; iTri < nlocal; iTri++)
                    scalar_cont->InsertNextValue((*scalar_cont_)(iTri));
           }
        }
        // vector containers
        VectorContainer<double,3> *vector_cont_;
        for(int ib = 0; ib < n_vector_containers_; ib++)
        {
           vector_cont_ = vector_containers_[ib][iMesh];
           if (vector_cont_)
           {
                cell_attributes.push_back(static_cast<vtkSmartPointer<vtkAbstractArray> >(vtkSmartPointer<vtkDoubleArray>::New()));
                vtkDoubleArray *vector_cont = static_cast<vtkDoubleArray*>(cell_attributes.back().GetPointer());
                vector_cont->SetNumberOfComponents(3);
                vector_cont->SetName(vector_container_names_[ib]);
                for (int iTri = 0; iTri < nlocal; iTri++)
                    vector_cont->InsertNextTupleValue((*vector_cont_)(iTri));
           }
        }
        // Construct polygonal mesh
        vtkSmartPointer<vtkPolyData> poly_data = vtkSmartPointer<vtkPolyData>::New();
        poly_data->SetPoints(points);
        poly_data->SetPolys(cells);

        std::list<vtkSmartPointer<vtkAbstractArray> >::iterator attributes;
        // cell attributes
        for (attributes = cell_attributes.begin(); attributes != cell_attributes.end(); attributes++)
            poly_data->GetCellData()->AddArray(*attributes);
        // point attributes
        for (attributes = point_attributes.begin(); attributes != point_attributes.end(); attributes++)
            poly_data->GetPointData()->AddArray(*attributes);

        if (nclusterprocs > 1)
        {
            if (filewriter)
            {
                std::list<vtkSmartPointer<vtkPolyData> > all_poly_data;
                vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
                // self
                all_poly_data.push_back(vtkSmartPointer<vtkPolyData>::New());
                all_poly_data.back()->DeepCopy(poly_data);
                #if VTK_MAJOR_VERSION <= 5
                // TODO check if this works
                appendFilter->AddInputConnection(all_poly_data.back()->GetProducerPort());
                #else
                appendFilter->AddInputData(all_poly_data.back());
                #endif
                // others
                for (int iproc = 1; iproc < nclusterprocs; iproc++)
                {
                    all_poly_data.push_back(vtkSmartPointer<vtkPolyData>::New());
                    vtkSmartPointer<vtkPolyData> remote_poly_data = vtkSmartPointer<vtkPolyData>::New();
                    int status = localController->Receive(all_poly_data.back().GetPointer(), comm->me+iproc, iproc);
                    if (status != 1)
                    {
                        
                        error->one(FLERR, "MPI receive error");
                    }
                    #if VTK_MAJOR_VERSION <= 5
                    // TODO check if this works
                    appendFilter->AddInputConnection(all_poly_data.back()->GetProducerPort());
                    #else
                    appendFilter->AddInputData(all_poly_data.back());
                    #endif
                }
                appendFilter->Update();
                poly_data = appendFilter->GetOutput();
            }
            else
            {
                int status = localController->Send(poly_data.GetPointer(), fileproc, comm->me-fileproc);
                if (status != 1)
                {
                    
                    error->one(FLERR, "MPI send error");
                }
            }
        }

        mbSet->SetBlock(cur_block++, poly_data);
        mbSet->GetMetaData(cur_block-1)->Set(mbSet->NAME(), (*mesh)->mesh_id());
    }

    return;
}

#endif //LAMMPS_VTK
