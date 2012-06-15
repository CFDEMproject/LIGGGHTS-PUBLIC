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

#include "tri_mesh.h"
#include <stdio.h>
#include <cstring>
#include <cstdlib>
#include "memory.h"
#include "myvector.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   constructor, destructor
------------------------------------------------------------------------- */

TriMesh::TriMesh(LAMMPS *lmp) : SurfaceMesh<3>(lmp)
{
}

TriMesh::~TriMesh()
{
}

/* ----------------------------------------------------------------------
   add a new triangle to mesh
------------------------------------------------------------------------- */

void TriMesh::addTriangle(double *a, double *b, double *c)
{
    double **nodeTmp = create<double>(nodeTmp,3,3);
    for(int i=0;i<3;i++){
      nodeTmp[0][i] = a[i];
      nodeTmp[1][i] = b[i];
      nodeTmp[2][i] = c[i];
    }
    addElement(nodeTmp);
    destroy<double>(nodeTmp);
}

/* ----------------------------------------------------------------------
   delete a triangle
------------------------------------------------------------------------- */

void TriMesh::deleteTriangle(int n)
{
    SurfaceMesh<3>::deleteElement(n);
    // triangle-specific code comes here
}

/* ----------------------------------------------------------------------
   stub
------------------------------------------------------------------------- */

void TriMesh::calcTriProperties()
{

}

