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

#ifndef LMP_SURFACE_MESH_H
#define LMP_SURFACE_MESH_H

#include "tracking_mesh.h"
#include "container.h"
#include "container.h"
#include "random_park.h"
#include "vector_liggghts.h"
#include "random_park.h"
#include "memory_ns.h"
#include "mpi_liggghts.h"
#include "comm.h"
#include <cmath>
#include "math_extra_liggghts.h"

#define EPSILON_CURVATURE 0.0001

using namespace LAMMPS_MEMORY_NS;

namespace LAMMPS_NS{

template<int NUM_NODES>
class SurfaceMesh : public TrackingMesh<NUM_NODES>
{
      public:

        //virtual double distToElem(int nElem, double *p) = 0;

        void setCurvature(double _curvature);
        void useAsInsertionMesh();
        void addElement(double **nodeToAdd);

        bool isPlanar();
        bool areCoplanar(int tag_i, int tag_j);
        bool isOnSurface(double *pos);

        void buildNeighbours();
        void growSurface(int iSrf, double by = 1e-13);

        void move(double *vecTotal, double *vecIncremental);
        void move(double *vecIncremental);
        void scale(double factor);

        virtual int generateRandomOwnedGhost(double *pos) = 0;
        virtual int generateRandomSubbox(double *pos) = 0;
        virtual int generateRandomSubboxWithin(double *pos,double delta) = 0;

        // public inline access

        // area of total mesh - all elements (all processes)
        
        inline double areaMeshGlobal()
        { return areaMesh_(0);}

        // area of owned elements
        inline double areaMeshOwned()
        { return areaMesh_(1);}

        // area of ghost elements
        inline double areaMeshGhost()
        { return areaMesh_(2);}

        // area of owned and ghost elements in my subdomain
        inline double areaMeshSubdomain()
        { return areaMesh_(3);}

        inline void surfaceNorm(int i,double *sn)
        { vectorCopy3D(surfaceNorm(i),sn); }

        inline double areaElem(int i)
        { return (area_)(i); }

      protected:

        SurfaceMesh(LAMMPS *lmp);
        virtual ~SurfaceMesh();

        void deleteElement(int n);

        void refreshOwned(int setupFlag);
        void refreshGhosts(int setupFlag);

        // (re)calculate properties
        inline void recalcVectors(); 
        inline void recalcLocalSurfProperties();
        inline void recalcGhostSurfProperties();

        // returns true if surfaces share an edge
        // called with local index
        // iEdge, jEdge return indices of first shared edge
        bool shareEdge(int i, int j, int &iEdge, int &jEdge);

        // mesh topology functions, called with local index
        void handleSharedEdge(int iSrf, int iEdge, int jSrf, int jEdge, bool coplanar);
        void handleSharedNode(int iSrf, int iNode, int jSrf, int jNode, bool coplanar);

        // calc properties
        void calcSurfPropertiesOfNewElement();
        void calcSurfPropertiesOfElement(int n);
        virtual double calcArea(int nElem) =0;
        void calcEdgeVecLen(int nElem, double *len, double **vec);
        void calcEdgeLen(int nElem, double *edgeLen);
        void calcSurfaceNorm(int nElem, double *surfNorm);
        void calcEdgeNormals(int nElem, double **edgeNorm);

        virtual bool isInElement(double *pos, int i) = 0;

        int randomOwnedGhostElement();

        void rotate(double *totalQ, double *dQ,double *totalDispl,double *dDispl);
        void rotate(double *dQ,double *dDispl);

        // inline access
        inline double&  area(int i)         {return (area_)(i);}
        inline double&  areaAcc(int i)      {return (areaAcc_)(i);}
        inline double*  edgeLen(int i)      {return (edgeLen_)(i);}
        inline double** edgeVec(int i)      {return (edgeVec_)(i);}
        inline double** edgeNorm(int i)     {return (edgeNorm_)(i);}
        inline double*  surfaceNorm(int i)  {return (surfaceNorm_)(i);}
        inline bool*    edgeActive(int i)   {return (edgeActive_)(i);}
        inline bool*    cornerActive(int i) {return (cornerActive_)(i);}
        inline bool*    hasNonCoplanarSharedNode(int i) {return (hasNonCoplanarSharedNode_)(i);}

      private:

        int searchElementByAreaAcc(double area,int lo, int hi);

        // flag indicating usage as insertion mesh
        bool isInsertionMesh_;

        // mesh properties
        double curvature_;
        ScalarContainer<double>& areaMesh_; 

        // per-element properties
        // add more properties as they are needed, such as
        // edge length, edge vectors, surface normal ....

        // surface properties
        ScalarContainer<double>& area_;
        ScalarContainer<double>& areaAcc_;                    // accumulated area of owned and ghost particles
        VectorContainer<double,NUM_NODES>& edgeLen_;          // len of edgeVec
        MultiVectorContainer<double,NUM_NODES,3>& edgeVec_;   // unit vec from node0 to node1 etc
        MultiVectorContainer<double,NUM_NODES,3>& edgeNorm_;
        VectorContainer<double,3>& surfaceNorm_;

        // neighbor topology
        
        ScalarContainer<int>& nNeighs_;
        VectorContainer<int,NUM_NODES>& neighFaces_;
        VectorContainer<bool,NUM_NODES>& hasNonCoplanarSharedNode_;
        VectorContainer<bool,NUM_NODES>& edgeActive_;
        VectorContainer<bool,NUM_NODES>& cornerActive_;
};

// *************************************
#include "surface_mesh_I.h"
// *************************************

} /* LAMMPS_NS */
#endif /* SURFACEMESH_H_ */
