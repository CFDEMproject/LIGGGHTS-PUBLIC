/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
   the producer of the LIGGGHTS® software and the CFDEM®coupling software
   See http://www.cfdem.com/terms-trademark-policy for details.

   LIGGGHTS® is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author for interpolated output:
   Felix Kleinfeldt (OVGU Magdeburg)
------------------------------------------------------------------------- */

#include "string.h"
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

using namespace LAMMPS_NS;

enum
{
    DUMP_STRESS = 1,
    DUMP_STRESSCOMPONENTS = 2,
    DUMP_ID = 4,
    DUMP_VEL = 8,
    DUMP_WEAR = 16,
    DUMP_TEMP = 32,
    DUMP_OWNER = 64,
    DUMP_AREA = 128,
    DUMP_AEDGES = 256,
    DUMP_ACORNERS = 512,
    DUMP_INDEX = 1024,
    DUMP_NNEIGHS = 2048,
    DUMP_MIN_ACTIVE_EDGE_DIST = 2048
};

enum
{
    DUMP_POINTS = 1,
    DUMP_FACE = 2
};

/* ---------------------------------------------------------------------- */

DumpMeshVTK::DumpMeshVTK(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg),
  dataMode_(0),
  nMesh_(0),
  meshList_(0),
  n_calls_(0),
  n_all_(0),
  n_all_max_(0),
  buf_all_(0),
  dump_what_(0),
  sigma_n_(0),
  sigma_t_(0),
  wear_(0),
  v_node_(0),
  f_node_(0),
  T_(0),
  min_active_edge_dist_(0),
  scalar_containers_(0),
  scalar_container_names_(0),
  n_scalar_containers_(0),
  vector_containers_(0),
  vector_container_names_(0),
  n_vector_containers_(0),
  container_args_(0),
  n_container_bases_(0)
{
  if (narg < 5)
    error->all(FLERR,"Illegal dump mesh/vtk command");

  //INFO: CURRENTLY ONLY PROC 0 writes

  format_default = NULL;

  nMesh_ = 0;

  int iarg = 5;

  bool hasargs = true;
  while (iarg < narg && hasargs)
  {
      hasargs = false;
      if(strcmp(arg[iarg],"output")==0)
      {
          if (iarg+2 > narg) error->all(FLERR,"Dump mesh/vtk: not enough arguments for 'interpolate'");
          if(strcmp(arg[iarg+1],"face")==0) dataMode_ = 0;
          else if(strcmp(arg[iarg+1],"interpolate")==0) dataMode_ = 1;
          else error->all(FLERR,"Dump mesh/vtk: wrong arrgument for 'output'");
          iarg += 2;
          hasargs = true;
      }
      else if(strcmp(arg[iarg],"stress")==0)
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
      else if(strcmp(arg[iarg],"distaa")==0)
      {
          dump_what_ |= DUMP_MIN_ACTIVE_EDGE_DIST;
          iarg++;
          hasargs = true;
      }
      else
      {
          // assume it's a mesh
          TriMesh **meshListNew = new TriMesh*[nMesh_+1];
          for(int i = 0; i < nMesh_; i++)
            meshListNew[i] = meshList_[i];
          if(meshList_) delete[] meshList_;
          meshList_ = meshListNew;

          int ifix = modify->find_fix(arg[iarg++]);

          FixMeshSurface *fms = (ifix < 0) ? 0 : dynamic_cast<FixMeshSurface*>(modify->fix[ifix]);

          if(!fms)
          {
              n_container_bases_++;
              memory->grow(container_args_,n_container_bases_,100,"container_args_");
              strcpy(container_args_[n_container_bases_-1],arg[iarg-1]);
              
          }
          else
          {
              meshList_[nMesh_] = fms->triMesh();
              fms->dumpAdd();
              nMesh_++;
          }
          hasargs = true;
      }
  }

  // in case meshes not specified explicitly, take all meshes
  if (nMesh_ == 0)
  {
      nMesh_ = modify->n_fixes_style("mesh/surface");

      meshList_ = new TriMesh*[nMesh_];
      for (int iMesh = 0; iMesh < nMesh_; iMesh++)
      {
          
          meshList_[iMesh] = static_cast<FixMeshSurface*>(modify->find_fix_style("mesh/surface",iMesh))->triMesh();
      }

      if (nMesh_ == 0)
          error->warning(FLERR,"Dump mesh/vtk cannot find any fix of type 'mesh/surface' to dump");
  }

  // allocate arrays
  sigma_n_ = new ScalarContainer<double>*[nMesh_];
  sigma_t_ = new ScalarContainer<double>*[nMesh_];
  wear_ = new ScalarContainer<double>*[nMesh_];
  v_node_ = new MultiVectorContainer<double,3,3>*[nMesh_];
  f_node_ = new VectorContainer<double,3>*[nMesh_];
  T_ = new ScalarContainer<double>*[nMesh_];
  min_active_edge_dist_ = new ScalarContainer<double>*[nMesh_];

  scalar_containers_ = new ScalarContainer<double>**[n_container_bases_];
  scalar_container_names_ = new char*[n_container_bases_];
  vector_containers_ = new VectorContainer<double,3>**[n_container_bases_];
  vector_container_names_ = new char*[n_container_bases_];
  for(int i = 0; i < n_container_bases_; i++)
  {
      scalar_containers_[i] = new ScalarContainer<double>*[nMesh_];
      vector_containers_[i] = new VectorContainer<double,3>*[nMesh_];
      scalar_container_names_[i] = new char[200];
      vector_container_names_[i] = new char[200];

      for(int j = 0; j < nMesh_; j++)
      {
          scalar_containers_[i][j] = 0;
          vector_containers_[i][j] = 0;
      }
  }

  if(dump_what_ == 0 && n_container_bases_ == 0)
    error->all(FLERR,"Dump mesh/vtk: No dump quantity selected");
}

/* ---------------------------------------------------------------------- */

DumpMeshVTK::~DumpMeshVTK()
{
  for (int iMesh = 0; iMesh < nMesh_; iMesh++)
  {
      static_cast<FixMeshSurface*>(modify->find_fix_id(meshList_[iMesh]->mesh_id()))->dumpRemove();
  }

  delete[] meshList_;
  memory->destroy(buf_all_);

  delete [] sigma_n_;
  delete [] sigma_t_;
  delete [] wear_;
  delete [] v_node_;
  delete [] f_node_;
  delete [] T_;
  delete [] min_active_edge_dist_;

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

void DumpMeshVTK::init_style()
{
  
  //multifile=1;             // 0 = one big file, 1 = one file per timestep
  //multiproc=0;             // 0 = proc 0 writes for all, 1 = one file/proc
  if (multifile != 1)
    error->all(FLERR,"You should use a filename like 'dump*.vtk' for the 'dump mesh/vtk' command to produce one file per time-step");
  if (multiproc != 0)
    error->all(FLERR,"Your 'dump mesh/vtk' command is writing one file per processor, where all the files contain the same data");

  if (domain->triclinic == 1)
    error->all(FLERR,"Can not dump VTK files for triclinic box");
  if (binary)
    error->all(FLERR,"Can not dump VTK files in binary mode");

  if(n_container_bases_ && dataMode_ == 1)
    error->all(FLERR,"Dump mesh/vtk: Can not dump interpolated general properties");

  // nodes
  size_one = 9;

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

  getGeneralRefs();

  size_one += n_scalar_containers_*1;
  size_one += n_vector_containers_*3;

  delete [] format;
}

/* ---------------------------------------------------------------------- */

int DumpMeshVTK::modify_param(int narg, char **arg)
{
  error->warning(FLERR,"dump_modify keyword is not supported by 'dump mesh/vtk' and is thus ignored");
  return 0;
}

/* ---------------------------------------------------------------------- */

void DumpMeshVTK::write_header(bigint ndump)
{
  write_header_ascii(ndump);
}

void DumpMeshVTK::write_header_ascii(bigint ndump)
{
  if (comm->me!=0) return;
  fprintf(fp,"# vtk DataFile Version 2.0\nLIGGGHTS mesh/VTK export\nASCII\n");
}

/* ---------------------------------------------------------------------- */

int DumpMeshVTK::count()
{
  int numTri = 0;

  n_calls_ = 0;
  n_all_ = 0;

  getRefs();

  for(int i = 0; i < nMesh_; i++)
  {
    if(!meshList_[i]->isParallel() && 0 != comm->me)
        continue;
    numTri += meshList_[i]->sizeLocal();
    
  }

  return numTri;
}

/* ---------------------------------------------------------------------- */

void DumpMeshVTK::getRefs()
{
  // add sizes and get references to properties - some may stay NULL
  if(dump_what_ & DUMP_STRESS)
  {
      for(int i = 0; i < nMesh_; i++)
      {
          sigma_n_[i] = meshList_[i]->prop().getElementProperty<ScalarContainer<double> >("sigma_n");
          sigma_t_[i] = meshList_[i]->prop().getElementProperty<ScalarContainer<double> >("sigma_t");
          if(0 == comm->me && (!sigma_n_[i] || !sigma_t_[i]))
            error->warning(FLERR,"Trying to dump stress for mesh which does not calculate stress, will dump '0' instead");
      }
  }
  if(dump_what_ & DUMP_STRESSCOMPONENTS)
  {
      for(int i = 0; i < nMesh_; i++)
      {
          f_node_[i] = meshList_[i]->prop().getElementProperty<VectorContainer<double,3> >("f");
      }
  }
  if(dump_what_ & DUMP_VEL)
  {
      for(int i = 0; i < nMesh_; i++)
      {
          v_node_[i] = meshList_[i]->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v");
          
      }
  }
  if(dump_what_ & DUMP_WEAR)
  {
      for(int i = 0; i < nMesh_; i++)
      {
          wear_[i] = meshList_[i]->prop().getElementProperty<ScalarContainer<double> >("wear");
          if(0 == comm->me && !wear_[i])
            error->warning(FLERR,"Trying to dump wear for mesh which does not calculate wear, will dump '0' instead");
      }
  }
  if(dump_what_ & DUMP_TEMP)
  {
      for(int i = 0; i < nMesh_; i++)
      {
          T_[i] = meshList_[i]->prop().getGlobalProperty<ScalarContainer<double> >("Temp");
          if(0 == comm->me && !T_[i])
            error->warning(FLERR,"Trying to dump temperature for mesh which does not calculate temperature, will dump '0' instead");
      }
  }
  if(dump_what_ & DUMP_MIN_ACTIVE_EDGE_DIST)
  {
      for(int i = 0; i < nMesh_; i++)
      {
          min_active_edge_dist_[i] = meshList_[i]->prop().getElementProperty<ScalarContainer<double> >("minActiveEdgeDist");
      }
  }
}

/* ---------------------------------------------------------------------- */

void DumpMeshVTK::getGeneralRefs()
{
  
  n_scalar_containers_ = 0;
  n_vector_containers_ = 0;
  char cid[200];

  for(int ib = 0; ib < n_container_bases_; ib++)
  {
      bool found_scalar = false, found_vector = false;

      for(int i = 0; i < nMesh_; i++)
      {
          if(meshList_[i]->prop().getElementProperty<ScalarContainer<double> >(container_args_[ib]))
          {
              
              found_scalar = true;
              scalar_containers_[n_scalar_containers_][i] = meshList_[i]->prop().getElementProperty<ScalarContainer<double> >(container_args_[ib]);
              scalar_containers_[n_scalar_containers_][i]->id(cid);
              strcpy(scalar_container_names_[n_scalar_containers_],cid);
          }
          
          if(meshList_[i]->prop().getElementProperty<VectorContainer<double,3> >(container_args_[ib]))
          {
              found_vector = true;
              vector_containers_[n_vector_containers_][i] = meshList_[i]->prop().getElementProperty<VectorContainer<double,3> >(container_args_[ib]);
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

void DumpMeshVTK::pack(int *ids)
{
  int m = 0;
  double node[3];
  TriMesh *mesh;
  ScalarContainer<double> *cb1;
  VectorContainer<double,3> *cb2;

  double **tmp;
  memory->create<double>(tmp,3,3,"DumpMeshVTK:tmp");

  // have to stick with this order (all per-element props)
  // as multiple procs pack
  for(int iMesh = 0;iMesh < nMesh_; iMesh++)
  {
    mesh = meshList_[iMesh];

    if(!mesh->isParallel() && 0 != comm->me)
        continue;

    int nlocal = meshList_[iMesh]->sizeLocal();

    for(int iTri = 0; iTri < nlocal; iTri++)
    {
        for(int j=0;j<3;j++)
        {
            mesh->node(iTri,j,node);
            for(int k=0;k<3;k++)
                buf[m++] = node[k];
        }

        if(dump_what_ & DUMP_STRESS)
        {
            buf[m++] = sigma_n_[iMesh] ? sigma_n_[iMesh]->get(iTri) : 0.;
            buf[m++] = sigma_t_[iMesh] ? sigma_t_[iMesh]->get(iTri) : 0.;
        }

        if(dump_what_ & DUMP_STRESSCOMPONENTS)
        {
            double invArea = 1./mesh->areaElem(iTri);
            double f[3];
            if(f_node_[iMesh])
                f_node_[iMesh]->get(iTri,f);
            else
                vectorZeroize3D(f);
            buf[m++] = f[0]*invArea;
            buf[m++] = f[1]*invArea;
            buf[m++] = f[2]*invArea;
        }

        if(dump_what_ & DUMP_ID)
        {
            buf[m++] = static_cast<double>(mesh->id(iTri));
        }

        if(dump_what_ & DUMP_VEL)
        {
            double avg[3];

            if(v_node_[iMesh])
            {
                // get vel for element, copy it to tmp
                v_node_[iMesh]->get(iTri,tmp);

                // calculate average
                vectorZeroize3D(avg);
                vectorCopy3D(tmp[0],avg);
                vectorAdd3D(tmp[1],avg,avg);
                vectorAdd3D(tmp[2],avg,avg);
                vectorScalarDiv3D(avg,3.);
            }
            else
                vectorZeroize3D(avg);

            // push to buffer
            buf[m++] = avg[0];
            buf[m++] = avg[1];
            buf[m++] = avg[2];
        }

        if(dump_what_ & DUMP_WEAR)
        {
            buf[m++] = wear_[iMesh] ? wear_[iMesh]->get(iTri) : 0.;
        }

        if(dump_what_ & DUMP_TEMP)
        {
            buf[m++] = T_[iMesh] ? T_[iMesh]->get(0) : 0.;
        }

        if(dump_what_ & DUMP_OWNER)
        {
            int me = comm->me;
            buf[m++] = static_cast<double>(me);
        }

        if(dump_what_ & DUMP_AREA)
        {
            buf[m++] = mesh->areaElem(iTri);
        }

        if(dump_what_ & DUMP_AEDGES)
        {
            buf[m++] = static_cast<double>(mesh->n_active_edges(iTri));
        }

        if(dump_what_ & DUMP_ACORNERS)
        {
            buf[m++] = static_cast<double>(mesh->n_active_corners(iTri));
        }
        if(dump_what_ & DUMP_INDEX)
        {
            buf[m++] = static_cast<double>(iTri);
        }
        if(dump_what_ & DUMP_NNEIGHS)
        {
            buf[m++] = static_cast<double>(mesh->nNeighs(iTri));
        }
        if(dump_what_ & DUMP_MIN_ACTIVE_EDGE_DIST)
        {
            buf[m++] = min_active_edge_dist_[iMesh] ? min_active_edge_dist_[iMesh]->get(iTri) : 0.;
        }

        for(int ib = 0; ib < n_scalar_containers_; ib++)
        {
           cb1 = scalar_containers_[ib][iMesh];
           buf[m++] = cb1 ? (*cb1)(iTri) : 0.;
        }

        for(int ib = 0; ib < n_vector_containers_; ib++)
        {
           cb2 = vector_containers_[ib][iMesh];
           buf[m++] = cb2 ? (*cb2)(iTri)[0] : 0.;
           buf[m++] = cb2 ? (*cb2)(iTri)[1] : 0.;
           buf[m++] = cb2 ? (*cb2)(iTri)[2] : 0.;
        }
    }
  }

  memory->destroy<double>(tmp);

  return;
}

/* ---------------------------------------------------------------------- */

void DumpMeshVTK::write_data(int n, double *mybuf)
{
    if (comm->me != 0) return;

    n_calls_++;

    // grow buffer if necessary
    if(n_all_+n*size_one > n_all_max_)
    {
        n_all_max_ = n_all_ + n*size_one;
        memory->grow(buf_all_,n_all_max_,"DumpMeshVTK:buf_all_");
    }

    // copy to buffer
    vectorCopyN(mybuf,&(buf_all_[n_all_]),n*size_one);
    n_all_ += n*size_one;

    // write on last call
    if(n_calls_ == comm->nprocs)
        write_data_ascii(n_all_/size_one,buf_all_);
}

/* ---------------------------------------------------------------------- */

void DumpMeshVTK::write_data_ascii(int n, double *mybuf)
{
    if(dataMode_ == 0)
        write_data_ascii_face(n,mybuf);
    else if(dataMode_ == 1)
        write_data_ascii_point(n,mybuf);
    return;
}

/* ---------------------------------------------------------------------- */

void DumpMeshVTK::write_data_ascii_point(int n, double *mybuf)
{
    int m, buf_pos;

    // n is the number of tri

    m = 0;
    buf_pos = 0;

    ScalarContainer<double> points("DumpMeshVTK::points");
    ScalarContainer<int> tri_points("DumpMeshVTK::tri_points");
    bool add;
    for (int i = 0; i < n; i++)
      {
        for (int j=0;j<3;j++)
        {
        add = true;
        for (int k=0; k < points.size(); k+=3)
        {
            if (points.get(k)==mybuf[m+0] && points.get(k+1)==mybuf[m+1] && points.get(k+2)==mybuf[m+2])
            {
            add=false;
            tri_points.add(k/3);
            break;
            }
        }
        if (add)
        {
            for (int k=0;k<3;k++)
            points.add(mybuf[m+k]);
            int max=tri_points.max();
            if (max<0) max=-1;
            tri_points.add(max+1);
        }
        m+=3;
        }
    m += size_one-9;
    }

    buf_pos += 9;

    // points_neightri
    class ScalarContainer<int> **points_neightri;
    points_neightri = new ScalarContainer<int>*[(int)points.size()/3];
    for (int i=0; i < points.size()/3; i++)
        points_neightri[i] = new ScalarContainer<int>("DumpMeshVTK::points_neightri");
    for (int i=0; i < 3*n; i+=3)
    {
        for (int j=0; j<3;j++)
        points_neightri[tri_points.get(i+j)]->add(i/3);
    }

    // write point data
      fprintf(fp,"DATASET UNSTRUCTURED_GRID\nPOINTS %d float\n", points.size()/3);
    for (int i=0; i < points.size(); i+=3)
        fprintf(fp,"%f %f %f\n",points.get(i+0),points.get(i+1),points.get(i+2));

    // write polygon data
      fprintf(fp,"CELLS %d %d\n",n,4*n);
      for (int i = 0; i < 3*n; i+=3) fprintf(fp,"%d %d %d %d\n",3,tri_points.get(i+0),tri_points.get(i+1),tri_points.get(i+2));
    // write cell types
      fprintf(fp,"CELL_TYPES %d\n",n);
      for (int i = 0; i < n; i++)
        fprintf(fp,"5\n");

    // write point data header
      fprintf(fp,"POINT_DATA %d\n",points.size()/3);

    // write cell data

      if(dump_what_ & DUMP_STRESS)
      {
        // write pressure and shear stress
        fprintf(fp,"SCALARS pressure float 1\nLOOKUP_TABLE default\n");
            m = buf_pos;
        for (int i = 0; i < points.size()/3; i++)
        {
        double helper=0;
        for (int j=0; j < points_neightri[i]->size();j++) helper += mybuf[m + points_neightri[i]->get(j)*size_one];
        helper /= points_neightri[i]->size();
        fprintf(fp,"%f\n",helper);
        }
        buf_pos++;

        // write shear stress
        fprintf(fp,"SCALARS shearstress float 1\nLOOKUP_TABLE default\n");
        m = buf_pos;
        for (int i = 0; i < points.size()/3; i++)
        {
        double helper=0;
        for (int j=0; j < points_neightri[i]->size();j++) helper += mybuf[m + points_neightri[i]->get(j)*size_one];
        helper /= points_neightri[i]->size();
        fprintf(fp,"%f\n",helper);
        }
        buf_pos++;
      }
      if(dump_what_ & DUMP_STRESSCOMPONENTS)
      {
        //write x y z stress component
             fprintf(fp,"VECTORS stress float\n");
        m = buf_pos;
        for (int i = 0; i < points.size()/3; i++)
        {
        double helper1=0, helper2=0, helper3=0;
        for (int j=0; j < points_neightri[i]->size();j++)
        {
            helper1 += mybuf[m + points_neightri[i]->get(j)*size_one];
            helper2 += mybuf[m+1 + points_neightri[i]->get(j)*size_one];
             helper3 += mybuf[m+2 + points_neightri[i]->get(j)*size_one];
        }
        helper1 /= points_neightri[i]->size();
        helper2 /= points_neightri[i]->size();
        helper3 /= points_neightri[i]->size();
        fprintf(fp,"%f %f %f\n",helper1,helper2,helper3);
        }
        buf_pos += 3;
     }

      if(dump_what_ & DUMP_ID)
      {
        // write id
        fprintf(fp,"SCALARS meshid float 1\nLOOKUP_TABLE default\n");
        m = buf_pos;
        for (int i = 0; i < points.size()/3; i++)
        {
        double helper=0;
        for (int j=0; j < points_neightri[i]->size();j++) helper += mybuf[m + points_neightri[i]->get(j)*size_one];
        helper /= points_neightri[i]->size();
        fprintf(fp,"%f\n",helper);
        }
        buf_pos++;
      }

      if(dump_what_ & DUMP_VEL)
      {
        //write vel data
        fprintf(fp,"VECTORS v float\n");
        m = buf_pos;
        for (int i = 0; i < points.size()/3; i++)
        {
        double helper1=0, helper2=0, helper3=0;
        for (int j=0; j < points_neightri[i]->size();j++)
        {
            helper1 += mybuf[m + points_neightri[i]->get(j)*size_one];
            helper2 += mybuf[m+1 + points_neightri[i]->get(j)*size_one];
            helper3 += mybuf[m+2 + points_neightri[i]->get(j)*size_one];
        }
        helper1 /= points_neightri[i]->size();
        helper2 /= points_neightri[i]->size();
        helper3 /= points_neightri[i]->size();
        fprintf(fp,"%f %f %f\n",helper1,helper2,helper3);
        }
        buf_pos += 3;
      }

      if(dump_what_ & DUMP_WEAR)
      {
        //write wear data
        fprintf(fp,"SCALARS wear float 1\nLOOKUP_TABLE default\n");
        m = buf_pos;
        for (int i = 0; i < points.size()/3; i++)
        {
        double helper=0;
        for (int j=0; j < points_neightri[i]->size();j++) helper += mybuf[m + points_neightri[i]->get(j)*size_one];
        helper /= points_neightri[i]->size();
        fprintf(fp,"%f\n",helper);
        }
        buf_pos++;
      }

      if(dump_what_ & DUMP_TEMP)
      {
        //write wear data
        fprintf(fp,"SCALARS Temp float 1\nLOOKUP_TABLE default\n");
        m = buf_pos;
        for (int i = 0; i < points.size()/3; i++)
        {
        double helper=0;
        for (int j=0; j < points_neightri[i]->size();j++) helper += mybuf[m + points_neightri[i]->get(j)*size_one];
        helper /= points_neightri[i]->size();
        fprintf(fp,"%f\n",helper);
        }
        buf_pos++;
      }

      if(dump_what_ & DUMP_OWNER)
      {
        //write owner data
        fprintf(fp,"SCALARS owner float 1\nLOOKUP_TABLE default\n");
        m = buf_pos;
        for (int i = 0; i < points.size()/3; i++)
        {
        double helper=0;
        for (int j=0; j < points_neightri[i]->size();j++) helper += mybuf[m + points_neightri[i]->get(j)*size_one];
        helper /= points_neightri[i]->size();
        fprintf(fp,"%f\n",helper);
        }
        buf_pos++;
      }

      if(dump_what_ & DUMP_AREA)
      {
        //write area data
        fprintf(fp,"SCALARS area float 1\nLOOKUP_TABLE default\n");
        m = buf_pos;
        for (int i = 0; i < points.size()/3; i++)
        {
        double helper=0;
        for (int j=0; j < points_neightri[i]->size();j++) helper += mybuf[m + points_neightri[i]->get(j)*size_one];
        helper /= points_neightri[i]->size();
        fprintf(fp,"%f\n",helper);
        }
        buf_pos++;
      }

      // not able to interpolate this
      if(dump_what_ & DUMP_AEDGES)
      {
        buf_pos++;
      }

      // not able to interpolate this
      if(dump_what_ & DUMP_ACORNERS)
      {
        buf_pos++;
      }
      // not able to interpolate this
      if(dump_what_ & DUMP_INDEX)
      {
        buf_pos++;
      }
      // not able to interpolate this
      if(dump_what_ & DUMP_NNEIGHS)
      {
        buf_pos++;
      }
      // not able to interpolate this
      if(dump_what_ & DUMP_MIN_ACTIVE_EDGE_DIST)
      {
        buf_pos++;
      }
    return;
}

/* ---------------------------------------------------------------------- */

void DumpMeshVTK::write_data_ascii_face(int n, double *mybuf)
{
  int k, m, buf_pos;

  // n is the number of elements

  // write point data
  fprintf(fp,"DATASET UNSTRUCTURED_GRID\nPOINTS %d float\n",3*n);
  m = 0;
  buf_pos = 0;
  for (int i = 0; i < n; i++)
  {
      fprintf(fp,"%f %f %f\n",mybuf[m+0],mybuf[m+1],mybuf[m+2]);
      fprintf(fp,"%f %f %f\n",mybuf[m+3],mybuf[m+4],mybuf[m+5]);
      fprintf(fp,"%f %f %f\n",mybuf[m+6],mybuf[m+7],mybuf[m+8]);
      m += size_one;
  }
  buf_pos += 9;

  // write polygon data
  fprintf(fp,"CELLS %d %d\n",n,4*n);
  k = 0;
  for (int i = 0; i < n; i++)
  {
      fprintf(fp,"%d %d %d %d\n",3,k,k+1,k+2);
      k += 3;
  }

  // write cell data header
  fprintf(fp,"CELL_TYPES %d\n",n);
  for (int i = 0; i < n; i++)
      fprintf(fp,"5\n");

  // write cell data header
  fprintf(fp,"CELL_DATA %d\n",n);

  // write cell data

  if(dump_what_ & DUMP_STRESS)
  {
      // write pressure and shear stress
      fprintf(fp,"SCALARS pressure float 1\nLOOKUP_TABLE default\n");
      m = buf_pos;
      for (int i = 0; i < n; i++)
      {
          fprintf(fp,"%f\n",mybuf[m]);
          m += size_one;
      }
      buf_pos++;

      // write shear stress
      fprintf(fp,"SCALARS shearstress float 1\nLOOKUP_TABLE default\n");
      m = buf_pos;
      for (int i = 0; i < n; i++)
      {
          fprintf(fp,"%f\n",mybuf[m]);
          m += size_one;
      }
      buf_pos++;
  }

  if(dump_what_ & DUMP_STRESSCOMPONENTS)
  {
      //write x y z stress component
      fprintf(fp,"VECTORS stress float\n");
      m = buf_pos;
      for (int i = 0; i < n; i++)
      {
          fprintf(fp,"%f %f %f\n",mybuf[m],mybuf[m+1],mybuf[m+2]);
          m += size_one;
      }
      buf_pos += 3;
  }

  if(dump_what_ & DUMP_ID)
  {
      // write id
      fprintf(fp,"SCALARS meshid float 1\nLOOKUP_TABLE default\n");
      m = buf_pos;
      for (int i = 0; i < n; i++)
      {
          fprintf(fp,"%f\n",mybuf[m]);
          m += size_one;
      }
      buf_pos++;
  }

  if(dump_what_ & DUMP_VEL)
  {
      //write vel data
      fprintf(fp,"VECTORS v float\n");
      m = buf_pos;
      for (int i = 0; i < n; i++)
      {
          fprintf(fp,"%f %f %f\n",mybuf[m],mybuf[m+1],mybuf[m+2]);
          m += size_one;
      }
      buf_pos += 3;
  }

  if(dump_what_ & DUMP_WEAR)
  {
      //write wear data
      fprintf(fp,"SCALARS wear float 1\nLOOKUP_TABLE default\n");
      m = buf_pos;
      for (int i = 0; i < n; i++)
      {
          fprintf(fp,"%f\n",mybuf[m]);
          m += size_one;
      }
      buf_pos++;
  }

  if(dump_what_ & DUMP_TEMP)
  {
      //write wear data
      fprintf(fp,"SCALARS Temp float 1\nLOOKUP_TABLE default\n");
      m = buf_pos;
      for (int i = 0; i < n; i++)
      {
          fprintf(fp,"%f\n",mybuf[m]);
          m += size_one;
      }
      buf_pos++;
  }

  if(dump_what_ & DUMP_OWNER)
  {
      //write owner data
      fprintf(fp,"SCALARS owner float 1\nLOOKUP_TABLE default\n");
      m = buf_pos;
      for (int i = 0; i < n; i++)
      {
          fprintf(fp,"%f\n",mybuf[m]);
          m += size_one;
      }
      buf_pos++;
  }

  if(dump_what_ & DUMP_AREA)
  {
      //write owner data
      fprintf(fp,"SCALARS area float 1\nLOOKUP_TABLE default\n");
      m = buf_pos;
      for (int i = 0; i < n; i++)
      {
          fprintf(fp,"%f\n",mybuf[m]);
          m += size_one;
      }
      buf_pos++;
  }

  if(dump_what_ & DUMP_AEDGES)
  {
      //write owner data
      fprintf(fp,"SCALARS active_edges float 1\nLOOKUP_TABLE default\n");
      m = buf_pos;
      for (int i = 0; i < n; i++)
      {
          fprintf(fp,"%f\n",mybuf[m]);
          m += size_one;
      }
      buf_pos++;
  }
  if(dump_what_ & DUMP_ACORNERS)
  {
      //write owner data
      fprintf(fp,"SCALARS active_corners float 1\nLOOKUP_TABLE default\n");
      m = buf_pos;
      for (int i = 0; i < n; i++)
      {
          fprintf(fp,"%f\n",mybuf[m]);
          m += size_one;
      }
      buf_pos++;
  }
  if(dump_what_ & DUMP_INDEX)
  {
      //write owner data
      fprintf(fp,"SCALARS index float 1\nLOOKUP_TABLE default\n");
      m = buf_pos;
      for (int i = 0; i < n; i++)
      {
          fprintf(fp,"%f\n",mybuf[m]);
          m += size_one;
      }
      buf_pos++;
  }
  if(dump_what_ & DUMP_NNEIGHS)
  {
      //write owner data
      fprintf(fp,"SCALARS nneighs float 1\nLOOKUP_TABLE default\n");
      m = buf_pos;
      for (int i = 0; i < n; i++)
      {
          fprintf(fp,"%f\n",mybuf[m]);
          m += size_one;
      }
      buf_pos++;
  }
  if(dump_what_ & DUMP_NNEIGHS)
  {
      //write owner data
      fprintf(fp,"SCALARS min_active_edge_dist float 1\nLOOKUP_TABLE default\n");
      m = buf_pos;
      for (int i = 0; i < n; i++)
      {
          fprintf(fp,"%f\n",mybuf[m]);
          m += size_one;
      }
      buf_pos++;
  }

  for(int ib = 0; ib < n_scalar_containers_; ib++)
  {
      //write owner data
      fprintf(fp,"SCALARS %s float 1\nLOOKUP_TABLE default\n",scalar_container_names_[ib]);
      m = buf_pos;
      
      for (int i = 0; i < n; i++)
      {
          fprintf(fp,"%f\n",mybuf[m]);
          m += size_one;
      }
      buf_pos++;
  }

  for(int ib = 0; ib < n_vector_containers_; ib++)
  {
      //write vector data
      fprintf(fp,"VECTORS %s float\n",vector_container_names_[ib]);
      m = buf_pos;
      for (int i = 0; i < n; i++)
      {
          fprintf(fp,"%f %f %f\n",mybuf[m],mybuf[m+1],mybuf[m+2]);
          m += size_one;
      }
      buf_pos += 3;
  }

  // footer not needed
  // if would be needed, would do like in dump stl
}
