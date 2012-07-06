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

/* ----------------------------------------------------------------------
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Philippe Seil (JKU Linz)
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
  nMesh_(0), meshList_(0), writeBinarySTL_(0)
{
  if (narg < 5)
    error->all(FLERR,"Illegal dump mesh/stl command");

  //INFO: CURRENTLY ONLY PROC 0 writes

  format_default = NULL;
  writeBinarySTL_ = 0;
  dump_what_ =  NLOCAL;

  int iarg = 5;

  while(iarg < narg){
    if(strncmp(arg[iarg],"binary",6) == 0){
      writeBinarySTL_ = 1;
      iarg++;
    } else if(strcmp(arg[iarg],"all") == 0){
      dump_what_ = NALL;
      iarg++;
    } else if(strcmp(arg[iarg],"local") == 0){
      dump_what_ = NLOCAL;
      iarg++;
    } else if(strcmp(arg[iarg],"ghost") == 0){
      dump_what_ = NGHOST;
      iarg++;
    } else error->all(FLERR,"Illegal dump mesh/stl command, unknown keyword");
  }

  nMesh_ = modify->n_fixes_style("mesh/surface");
  
  if (nMesh_ == 0)
    error->warning(FLERR,"Dump mesh/stl cannot find any fix of type 'mesh/gran' to dump");

  meshList_ = new TriMesh*[nMesh_];
  for (int iMesh = 0; iMesh < nMesh_; iMesh++)
  {
      
      meshList_[iMesh] =static_cast<FixMeshSurface*>(modify->find_fix_style("mesh/surface",iMesh))->triMesh();
  }
}

/* ---------------------------------------------------------------------- */

DumpMeshSTL::~DumpMeshSTL()
{
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
  if(writeBinarySTL_) write_header_binary(ndump);
  else write_header_ascii(ndump);
}

void DumpMeshSTL::write_header_binary(bigint ndump)
{
  // ndump = # of dump lines this proc will contribute to dump
  if(comm->me != 0) return;
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

  n_calls_ = 0;

  for(int i = 0; i < nMesh_; i++)
  {
      bounds(i,ilo,ihi);
      numTri += ihi-ilo;
  }

  return numTri;
}

/* ---------------------------------------------------------------------- */

void DumpMeshSTL::bounds(int imesh,int &ilo, int &ihi)
{
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
  double node[3],surfaceNorm[3];
  TriMesh *mesh;

  for(int iMesh = 0;iMesh < nMesh_; iMesh++)
  {
    mesh = meshList_[iMesh];
    bounds(iMesh,ilo,ihi);

    for(int iTri = ilo; iTri < ihi; iTri++)
    {
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

  if(writeBinarySTL_) write_data_binary(n,mybuf);
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
  
  if(n && n_calls_ == comm->nprocs)
    fprintf(fp,"endsolid LIGGGHTS_STL_EXPORT\n");
}
