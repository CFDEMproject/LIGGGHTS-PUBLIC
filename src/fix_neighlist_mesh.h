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

#ifdef FIX_CLASS

FixStyle(neighlist/mesh,FixNeighlistMesh)

#else

#ifndef LMP_FIX_NEIGHLIST_MESH_H
#define LMP_FIX_NEIGHLIST_MESH_H

#include "fix.h"
#include "container.h"

namespace LAMMPS_NS
{

class FixNeighlistMesh : public Fix
{
  public:

    FixNeighlistMesh(LAMMPS *lmp, int narg, char **arg);
    virtual
    ~FixNeighlistMesh();
    virtual int setmask();
    virtual void post_create();
    virtual void pre_delete(bool unfixflag);

    void initializeNeighlist();
    virtual void setup_pre_force(int);

    virtual void pre_neighbor(); 
    virtual void pre_force(int vflag); 

    void getPointers(int *&cList, int *&nContact);

    int getSizeNumContacts() {return numContacts.size();}
    int getTotalNumContacts() {return numAllContacts_;}

  private:

    void handleTriangle(int iTri);
    void getBinBoundariesFromBoundingBox(class BoundingBox &b,
        int &ixMin,int &ixMax,int &iyMin,int &iyMax,int &izMin,int &izMax);

    class FixMeshSurface *caller_;
    class TriMesh *mesh_;

    bool buildNeighList;

    ScalarContainer<int> contactList, numContacts;

    int numAllContacts_;

    int mbinx,mbiny,mbinz,maxhead, *bins, *binhead;
    double skin;

    // max distance from triangle to particle center
    double distmax;

    double **x, *r;

    bool movingMesh;
};

} /* namespace LAMMPS_NS */

#endif /* FIX_MESH_NEIGHLIST_H_ */
#endif /* FIX_CLASS */
