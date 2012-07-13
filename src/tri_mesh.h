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

#ifndef LMP_TRI_MESH_H
#define LMP_TRI_MESH_H

#include "surface_mesh.h"
#include "math_extra_liggghts.h"
#include <fstream>

namespace LAMMPS_NS
{

  class TriMesh : public SurfaceMesh<3>
  {
      public:

        TriMesh(LAMMPS *lmp);
        virtual ~TriMesh();

        void addTriangle(double *a, double *b, double *c);
        void addTriangleComplete();
        void deleteTriangle(int n);

        double resolveTriSphereContact(int nTri, double rSphere, double *cSphere, double *delta);
        double resolveTriSphereContactBary(int nTri, double rSphere, double *cSphere, double *contactPoint,double *bary);
        bool resolveTriSphereNeighbuild(int nTri, double rSphere, double *cSphere, double treshold);

        int generateRandomOwnedGhost(double *pos);
        int generateRandomSubbox(double *pos);
        int generateRandomSubboxWithin(double *pos,double delta);

      protected:

        double calcArea(int n);
        bool isInElement(double *pos,int i);

      private:

        double calcDist(double *cs, double *closestPoint, double *en0);
        double calcDistToPlane(double *p, double *pPlane, double *nPlane);
        double resolveEdgeContact(int iTri, int iEdge, double *p, double *pPlane, double *en0, double *bary);
  };

  // *************************************
  #include "tri_mesh_I.h"
  // *************************************

} /* LAMMPS_NS */
#endif /* TRIMESH_H_ */
