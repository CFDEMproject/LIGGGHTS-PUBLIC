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

    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Philippe Seil (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#include "string.h"
#include "dump_mesh_stl.h"
#include "tri_mesh.h"
#include "domain.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "error.h"
#include "fix.h"
#include "fix_mesh_surface.h"
#include "region_mesh_tet.h"
#include "region.h"
#include "modify.h"
#include "comm.h"
#include <stdint.h>

using namespace LAMMPS_NS;

enum
{
    NLOCAL,
    NGHOST,
    NALL
};

/* ---------------------------------------------------------------------- */

DumpMeshSTL::DumpMeshSTL(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg),
  nMesh_(0),
  meshList_(0),
  iregion_(-1)
{
  if (narg < 5)
    error->all(FLERR,"Illegal dump mesh/stl command");

  //INFO: CURRENTLY ONLY PROC 0 writes

  format_default = NULL;
  binary = 0;
  dump_what_ =  NLOCAL;

  int iarg = 5;

  while(iarg < narg){
    if(strncmp(arg[iarg],"binary",6) == 0){
      binary = 1;
      iarg++;
    } else if(strcmp(arg[iarg],"region") == 0){
      if (narg < iarg+2)
        error->all(FLERR,"Illegal dump mesh/stl command, not enough arguments");
      iregion_ = domain->find_region(arg[iarg+1]);
      if (iregion_ == -1)
        error->all(FLERR,"Illegal dump mesh/stl command, region ID does not exist");
      iarg += 2;
    } else if(strcmp(arg[iarg],"all") == 0){
      dump_what_ = NALL;
      iarg++;
    } else if(strcmp(arg[iarg],"local") == 0){
      dump_what_ = NLOCAL;
      iarg++;
    } else if(strcmp(arg[iarg],"ghost") == 0){
      dump_what_ = NGHOST;
      iarg++;
    } else {

      TriMesh **meshListNew = new TriMesh*[nMesh_+1];
      for(int i = 0; i < nMesh_; i++)
        meshListNew[i] = meshList_[i];
      delete[] meshList_;
      meshList_ = meshListNew;

      // assume it's a mesh
      int ifix = modify->find_fix(arg[iarg]);
      if(ifix >= 0)
      {
          FixMeshSurface *fms = static_cast<FixMeshSurface*>(modify->fix[ifix]);
          meshList_[nMesh_] = fms->triMesh();
          fms->dumpAdd();
      }
      // assume it's a surface of a tet mesh
      else
      {
         int ireg = domain->find_region(arg[iarg]);

         if(ireg < 0)
            error->all(FLERR,"Illegal dump mesh/stl command, unknown keyword or mesh");

         // only mesh/tet style valid
         if(strcmp(domain->regions[ireg]->style,"mesh/tet"))
            error->all(FLERR,"Illegal dump mesh/stl command, region found, but not of style 'mesh/tet'");

         // region is of right style
         meshList_[nMesh_] = static_cast<RegTetMesh*>(domain->regions[ireg])->get_tri_mesh();
      }
      iarg++;
      nMesh_++;
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
          static_cast<FixMeshSurface*>(modify->find_fix_style("mesh/surface",iMesh))->dumpAdd();
      }

      if (nMesh_ == 0)
          error->warning(FLERR,"Dump mesh/stl cannot find any fix of type 'mesh/surface' to dump");
  }
}

/* ---------------------------------------------------------------------- */

DumpMeshSTL::~DumpMeshSTL()
{
  for (int iMesh = 0; iMesh < nMesh_; iMesh++)
  {
      if(meshList_[iMesh]->mesh_id())
      {
          Fix *f = modify->find_fix_id(meshList_[iMesh]->mesh_id());
          if(f)
            static_cast<FixMeshSurface*>(f)->dumpRemove();
      }
  }

  delete[] meshList_;
}

/* ---------------------------------------------------------------------- */

void DumpMeshSTL::init_style()
{
  //multifile=1;             // 0 = one big file, 1 = one file per timestep
  //multiproc=0;             // 0 = proc 0 writes for all, 1 = one file/proc
  if (multifile != 1)
    error->all(FLERR,"You should use a filename like 'dump*.stl' for the 'mesh/stl' command to produce one file per time-step");

/*
  if (multiproc != 0)
    error->all(FLERR,"Your 'dump mesh/stl' command is writing one file per processor, where all the files contain the same data");
*/

  size_one = 12;

  delete [] format;
  format = new char[150];
  strcpy(format,"  facet normal %g %g %g\n");
  strcat(format,"    outer loop\n");
  strcat(format,"      vertex %g %g %g\n");
  strcat(format,"      vertex %g %g %g\n");
  strcat(format,"      vertex %g %g %g\n");
  strcat(format,"    endloop\n");
  strcat(format,"  endfacet\n");
}

/* ---------------------------------------------------------------------- */

int DumpMeshSTL::modify_param(int narg, char **arg)
{
  error->warning(FLERR,"dump_modify keyword is not supported by 'dump stl' and is thus ignored");
  return 0;
}

/* ---------------------------------------------------------------------- */

void DumpMeshSTL::write_header(bigint ndump)
{
  if(binary) write_header_binary(ndump);
  else write_header_ascii(ndump);
}

void DumpMeshSTL::write_header_binary(bigint ndump)
{
  // ndump = # of dump lines this proc will contribute to dump
  if(!multiproc && comm->me != 0) return;
  char *header;
  header = new char[81];

  strcpy(header,"STL File written by LIGGGHTS");
  for(int i=strlen(header);i<80;i++)
    strcat(header," ");
  fwrite((void*)header, 1, 80, fp);

  uint32_t nTri = nme;
  fwrite((void*)&nTri,1,4,fp);

  delete[] header;
}

void DumpMeshSTL::write_header_ascii(bigint ndump)
{
  if (!multiproc && comm->me != 0) return;
  fprintf(fp,"solid LIGGGHTS_STL_EXPORT\n");
}

/* ---------------------------------------------------------------------- */

int DumpMeshSTL::count()
{
  int ilo, ihi;
  int numTri = 0;
  double center[3];

  n_calls_ = 0;

  for(int imesh = 0; imesh < nMesh_; imesh++)
  {
      bounds(imesh,ilo,ihi);
      if(iregion_ == -1)
      {
          numTri += ihi-ilo;
      }
      else
      {
          for(int i = ilo; i < ihi; i++)
          {
              meshList_[imesh]->center(i,center);
              if(domain->regions[iregion_]->match(center[0],center[1],center[2]))
                numTri++;
          }
      }
  }

  return numTri;
}

/* ---------------------------------------------------------------------- */

void DumpMeshSTL::bounds(int imesh,int &ilo, int &ihi)
{
    if(!meshList_[imesh]->isParallel() && 0 != comm->me)
    {
        ilo = ihi = 0;
        return;
    }

    if(dump_what_ == NLOCAL)
    {
        ilo = 0;
        ihi = meshList_[imesh]->sizeLocal();
    }
    else if(dump_what_ == NGHOST)
    {
        ilo = meshList_[imesh]->sizeLocal();
        ihi = meshList_[imesh]->sizeLocal() + meshList_[imesh]->sizeGhost();
    }
    else if(dump_what_ == NALL)
    {
        ilo = 0;
        ihi = meshList_[imesh]->sizeLocal() + meshList_[imesh]->sizeGhost();
    }
    return;
}

/* ---------------------------------------------------------------------- */

void DumpMeshSTL::pack(int *ids)
{
  int ilo,ihi;
  int m = 0;
  double node[3],surfaceNorm[3], center[3];
  TriMesh *mesh;

  for(int iMesh = 0;iMesh < nMesh_; iMesh++)
  {
    mesh = meshList_[iMesh];
    bounds(iMesh,ilo,ihi);

    for(int iTri = ilo; iTri < ihi; iTri++)
    {
        if(iregion_ >= 0)
        {
            mesh->center(iTri,center);
            if(!domain->regions[iregion_]->match(center[0],center[1],center[2]))
                continue;
        }

        mesh->surfaceNorm(iTri,surfaceNorm);
        for(int j=0;j<3;j++)
            buf[m++] = surfaceNorm[j];

        for(int j=0;j<3;j++)
        {
            mesh->node(iTri,j,node);
            for(int k=0;k<3;k++)
                buf[m++] = node[k];
        }
    }
  }
  return;
}

/* ---------------------------------------------------------------------- */

void DumpMeshSTL::write_data(int n, double *mybuf)
{
  
  if(!multiproc && comm->me != 0) return;

  if(binary) write_data_binary(n,mybuf);
  else write_data_ascii(n,mybuf);
}

void DumpMeshSTL::write_data_binary(int n, double *mybuf)
{
  double *cursor = mybuf;
  for(int i=0;i<n;i++){
    for(int j=0;j<size_one;j++,cursor++){
      float dummy = *cursor;
      fwrite((void*)&dummy,4,1,fp);
    }
    uint16_t x = 0;
    fwrite((void*)&x,2,1,fp);

  }
}

void DumpMeshSTL::write_data_ascii(int n, double *mybuf)
{
  
  n_calls_++;

  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp,format,
            mybuf[m],mybuf[m+1],mybuf[m+2],mybuf[m+3],mybuf[m+4],mybuf[m+5],
            mybuf[m+6],mybuf[m+7],mybuf[m+8],mybuf[m+9],mybuf[m+10],mybuf[m+11]
    );
    m += size_one;
  }
  //write footer
  
  if(multiproc || n_calls_ == comm->nprocs)
    fprintf(fp,"endsolid LIGGGHTS_STL_EXPORT\n");
}
