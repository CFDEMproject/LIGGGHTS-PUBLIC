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
   Philippe Seil (JKU Linz)
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
------------------------------------------------------------------------- */

#include "fix_neighlist_mesh.h"
#include "fix_mesh_surface.h"
#include "modify.h"
#include "container.h"
#include "bounding_box.h"
#include "neighbor.h"
#include "atom.h"
#include "vector_liggghts.h"
#include "update.h"
#include <stdio.h>

using namespace LAMMPS_NS;
using namespace FixConst;

#define SMALL_DELTA skin/(70.*M_PI)

FixNeighlistMesh::FixNeighlistMesh(LAMMPS *lmp, int narg, char **arg)
: Fix(lmp,narg,arg),
  contactList("contactList"),
  numContacts("numContacts"),
  buildNeighList(false),
  movingMesh(false)
{
    if(!modify->find_fix_id(arg[3]) || !dynamic_cast<FixMeshSurface*>(modify->find_fix_id(arg[3])))
        error->fix_error(FLERR,this,"illegal caller");

    caller_ = static_cast<FixMeshSurface*>(modify->find_fix_id(arg[3]));
    mesh_ = caller_->triMesh();
}

/* ---------------------------------------------------------------------- */

FixNeighlistMesh::~FixNeighlistMesh()
{
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::post_create()
{
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::initializeNeighlist()
{
    // remove old lists, init new ones
    
    for(int i = 0; i < numContacts.size(); i++)
        numContacts.del(i);
    for(int i = 0; i < contactList.size(); i++)
        contactList.del(i);

    int nall = mesh_->sizeLocal()+mesh_->sizeGhost();

    for(int iTri = 0; iTri < nall; iTri++)
        numContacts.add(0);

    contactList.add(0);
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::setup_pre_force(int foo)
{
    
    pre_neighbor();
    pre_force(0);
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::pre_delete(bool unfixflag)
{

}

/* ---------------------------------------------------------------------- */

int FixNeighlistMesh::setmask()
{
    int mask = 0;
    mask |= PRE_NEIGHBOR;
    mask |= PRE_FORCE;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::pre_neighbor()
{
    buildNeighList = true;
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::pre_force(int vflag)
{
    if(!buildNeighList) return;

    movingMesh = mesh_->isMoving();

    buildNeighList = false;

    contactList.empty();
    numContacts.empty();
    numAllContacts_ = 0;

    x = atom->x;
    r = atom->radius;

    if(neighbor->style != 1)
        error->all(FLERR,"Please use style 'bin' in the 'neighbor' command together with triangular walls");

    double rmax = 0.5*(neighbor->cutneighmax - neighbor->skin);
    if(movingMesh)
    {
      skin = neighbor->skin;
      distmax = neighbor->cutneighmax + SMALL_DELTA;
    }
    else
    {
      skin = 0.5*neighbor->skin;
      distmax = neighbor->cutneighmax - rmax + SMALL_DELTA;
    }

    mbinx = neighbor->mbinx;
    mbiny = neighbor->mbiny;
    mbinz = neighbor->mbinz;
    bins = neighbor->bins;
    binhead = neighbor->binhead;
    maxhead = neighbor->maxhead;

    int nall = mesh_->sizeLocal() + mesh_->sizeGhost();

    for(int iTri = 0; iTri < nall; iTri++)
      handleTriangle(iTri);

    MPI_Sum_Scalar(numAllContacts_,world);
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::handleTriangle(int iTri)
{
    int *mask = atom->mask;

    // get bounding box of element on this subdomain
    
    BoundingBox b = mesh_->getElementBoundingBoxOnSubdomain(iTri);

    int ixMin(0),ixMax(0),iyMin(0),iyMax(0),izMin(0),izMax(0);
    int nlocal = atom->nlocal;
    double lo[3],hi[3];

    // extend bbox by cutneighmax and get bin boundaries
    getBinBoundariesFromBoundingBox(b,ixMin,ixMax,iyMin,iyMax,izMin,izMax);

    int numContTmp = 0;

    // only do this if I own particles
    if(nlocal)
    {
        for(int ix=ixMin;ix<=ixMax;ix++)
          for(int iy=iyMin;iy<=iyMax;iy++)
            for(int iz=izMin;iz<=izMax;iz++)
            {
              int iBin = iz*mbiny*mbinx + iy*mbinx + ix;
              if(iBin < 0 || iBin >= maxhead) continue;

              int iAtom = binhead[iBin];
              while(iAtom != -1 && iAtom < nlocal)
              {
                if(! (mask[iAtom] & groupbit))
                {
                    if(bins) iAtom = bins[iAtom];
                    else iAtom = -1;
                    continue;
                }

                if(mesh_->resolveTriSphereNeighbuild(iTri,r ? r[iAtom] : 0. ,x[iAtom],r ? skin : (distmax+skin) ))
                {
                  
                  numContTmp++;
                  contactList.add(iAtom);
                }
                if(bins) iAtom = bins[iAtom];
                else iAtom = -1;
              }
            }
    }

    numContacts.add(numContTmp);
    numAllContacts_ += numContTmp;
    
    return;
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::getBinBoundariesFromBoundingBox(BoundingBox &b,
      int &ixMin,int &ixMax,int &iyMin,int &iyMax,int &izMin,int &izMax)
{
    double delta = distmax;
    double tri_xmin[3] = {b.xLo-delta,b.yLo-delta,b.zLo-delta};
    double tri_xmax[3] = {b.xHi+delta,b.yHi+delta,b.zHi+delta};

    int binmin = neighbor->coord2bin(tri_xmin,ixMin,iyMin,izMin);
    
    int binmax= neighbor->coord2bin(tri_xmax,ixMax,iyMax,izMax);
    
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::getPointers(int *&cList, int *&nContact)
{
    cList = contactList.begin();
    nContact = numContacts.begin();
}
