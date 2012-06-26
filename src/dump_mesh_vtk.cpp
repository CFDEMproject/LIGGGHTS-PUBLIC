/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
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

enum{
    DUMP_STRESS = 1,
    DUMP_ID = 2,
    DUMP_WEAR = 4,
    DUMP_VEL = 8,
    DUMP_STRESSCOMPONENTS = 16,
    DUMP_TEMP = 32,
    DUMP_OWNER = 64,
    };

/* ---------------------------------------------------------------------- */

DumpMeshVTK::DumpMeshVTK(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg),
  nMesh_(0), meshList_(0), dump_what_(0),
  n_calls_(0), n_all_(0), n_all_max_(0),buf_all_(0)
{
  if (narg < 5)
    error->all(FLERR,"Illegal dump mesh/vtk command");

  //INFO: CURRENTLY ONLY PROC 0 writes

  format_default = NULL;

  nMesh_ = modify->n_fixes_style("mesh/surface");
  
  if (nMesh_ == 0)
    error->warning(FLERR,"Dump mesh/vtk cannot find any fix of type 'mesh/gran' to dump");

  meshList_ = new TriMesh*[nMesh_];
  for (int iMesh = 0; iMesh < nMesh_; iMesh++)
  {
      
      meshList_[iMesh] = static_cast<FixMeshSurface*>(modify->find_fix_style("mesh/surface",iMesh))->triMesh();
  }

  int iarg = 5;
  dump_what_ = 0;

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
      else if(strcmp(arg[iarg],"wear")==0)
      {
          dump_what_ |= DUMP_WEAR;
          iarg++;
          hasargs = true;
      }
      else if(strcmp(arg[iarg],"vel")==0)
      {
          dump_what_ |= DUMP_VEL;
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
  }

  if(dump_what_ == 0)
    error->all(FLERR,"Dump mesh/vtk: No dump quantity selected");
}

/* ---------------------------------------------------------------------- */

DumpMeshVTK::~DumpMeshVTK()
{
  delete[] meshList_;
  memory->destroy(buf_all_);
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

  // nodes
  size_one = 9;

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

  delete [] format;
}

/* ---------------------------------------------------------------------- */

int DumpMeshVTK::modify_param(int narg, char **arg)
{
  error->warning(FLERR,"dump_modify keyword is not supported by 'dump stl' and is thus ignored");
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

  for(int i = 0; i < nMesh_; i++)
    numTri += meshList_[i]->sizeLocal();

  return numTri;
}

/* ---------------------------------------------------------------------- */

void DumpMeshVTK::pack(int *ids)
{
  int m = 0;
  double node[3];
  TriMesh *mesh;

  double **tmp;
  memory->create<double>(tmp,3,3,"DumpMeshVTK:tmp");

  // have to stick with this order (all per-element props)
  // as multiple procs pack
  for(int iMesh = 0;iMesh < nMesh_; iMesh++)
  {
    mesh = meshList_[iMesh];
    int nlocal = meshList_[iMesh]->sizeLocal();

    for(int iTri=0;iTri<nlocal;iTri++)
    {
        for(int j=0;j<3;j++)
        {
            mesh->node(iTri,j,node);
            for(int k=0;k<3;k++)
                buf[m++] = node[k];
        }
        if(dump_what_ & DUMP_STRESS)
        {
            buf[m++] = /*TODO P*/ 0.;
            buf[m++] = /*TODO SHEARSTRESS*/ 0.;
        }
        if(dump_what_ & DUMP_STRESSCOMPONENTS)
        {
            buf[m++] = /*TODO sx*/ 0.;
            buf[m++] = /*TODO sy*/ 0.;
            buf[m++] = /*TODO sz*/ 0.;
        }
        if(dump_what_ & DUMP_ID)
        {
            buf[m++] = static_cast<double>(mesh->id(iTri));
        }
        if(dump_what_ & DUMP_VEL)
        {
            MultiVectorContainer<double,3,3> *v_node;
            double avg[3];

            v_node = mesh->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v");

            if(v_node)
            {
                // get vel for element, copy it to tmp
                v_node->get(iTri,tmp);

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
            buf[m++] = /*TODO wear*/ 0.;
        }
        if(dump_what_ & DUMP_TEMP)
        {
            buf[m++] = /*TODO temp*/ 0.;
        }
        if(dump_what_ & DUMP_OWNER)
        {
            int me = comm->me;
            buf[m++] = static_cast<double>(me);
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

void DumpMeshVTK::write_data_ascii(int n, double *mybuf)
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

  // footer not needed
  // if would be needed, would do like in dump stl
}
