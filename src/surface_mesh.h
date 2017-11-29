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

#ifndef LMP_SURFACE_MESH_H
#define LMP_SURFACE_MESH_H

#include "tracking_mesh.h"
#include "container.h"
#include "container.h"
#include "random_park.h"
#include "vector_liggghts.h"
#include "random_park.h"
#include "memory_ns.h"
#include "region_neighbor_list.h"
#include "mpi_liggghts.h"
#include "comm.h"
#include "math_extra_liggghts.h"
#include <cmath>
#include <algorithm>

#define EPSILON_CURVATURE 0.00001

using namespace LAMMPS_MEMORY_NS;

namespace LAMMPS_NS{

template<int NUM_NODES, int NUM_NEIGH_MAX>
class SurfaceMesh : public TrackingMesh<NUM_NODES>
{
      public:

        //virtual double distToElem(int nElem, double *p) = 0;

        void setCurvature(double _curvature);
        void setCurvatureTolerant(bool _tol);
        
        bool addElement(double **nodeToAdd,int lineNumb);

        bool isPlanar();
        bool areCoplanar(int tag_i, int tag_j);
        bool areCoplanarNodeNeighs(int tag_i, int tag_j);
        bool areCoplanarNeighs(int tag_i, int tag_j);
        bool isOnSurface(double *pos);

        // returns true if surfaces share an edge
        // called with local index
        // iEdge, jEdge return indices of first shared edge
        bool shareEdge(int i, int j, int &iEdge, int &jEdge);

        void move(const double * const vecTotal, const double * const vecIncremental);
        void move(const double * const vecIncremental);
        void scale(double factor);

        using MultiNodeMesh<NUM_NODES>::rotate;
        void rotate(const double * const totalQ, const double * const dQ, const double * const origin);
        void rotate(const double * const dQ, const double * const origin);

        virtual int generateRandomOwnedGhost(double *pos) = 0;
        virtual int generateRandomOwnedGhostWithin(double *pos,double delta) = 0;
        virtual int generateRandomSubbox(double *pos) = 0;

        inline double edgeEdgeDist(int iSrf, int iEdge, int jSrf, int jEdge);
        inline double edgeNodeDist(int iSrf, int iEdge, int jSrf, int jNode);
        inline double edgePointDist(int iSrf, int iEdge, double *point);

        int n_active_edges(int i);
        int n_active_corners(int i);

        void extrudePlanarMesh(const double length, double * &extrusion_tri_nodes, int &extrusion_tri_count);

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

        inline void edgeLen(int i,double *el)
        { vectorCopy3D(edgeLen_(i),el); }

        inline double edgeLen(int i,int iEdge)
        { return edgeLen_(i)[iEdge]; }

        inline void edgeVec(int i,int j,double *ev)
        { vectorCopy3D(edgeVec_(i)[j],ev); }

        inline void surfaceNorm(int i,double *sn)
        { vectorCopy3D(surfaceNorm(i),sn); }

        inline double areaElem(int i)
        { return (area_)(i); }

        inline int nNeighs(int i)
        { return nNeighs_(i); }

        inline int neighFaces(int i,int j)
        { return neighFaces_(i)[j]; }

        inline bool edgeActive(int i,int j)
        { return (edgeActive_)(i)[j]; }

        inline bool cornerActive(int i,int j)
        { return (cornerActive_)(i)[j]; }

        static const int NO_OBTUSE_ANGLE = -1;

      protected:

        SurfaceMesh(LAMMPS *lmp);
        virtual ~SurfaceMesh();

        void deleteElement(int n);
        void qualityCheck();

        void buildNeighbours();
        virtual void parallelCorrectionActiveInactive();
        virtual void parallelCorrectionNeighs();

        bool edgeVecsColinear(double *v,double *w);
        bool coplanarNeighsOverlap(int iSrf,int iEdge,int jSrf,int jEdge);

        // mesh topology functions, called with local index
        void handleSharedEdge(int iSrf, int iEdge, int jSrf, int jEdge,
                          bool coplanar,bool neighflag = true);
        virtual int handleCorner(int iSrf, int iNode,int *idListVisited,int *idListHasNode,
                          double **edgeList,double **edgeEndPoint);
        void checkNodeRecursive(int iSrf,double *nodeToCheck,int &nIdListVisited,
                          int *idListVisited,int &nIdListHasNode,int *idListHasNode,
                          double **edgeList,double **edgeEndPoint,bool &anyActiveEdge);

        #include "surface_mesh_feature_remove.h"

        // (re) calc properties
        void calcSurfPropertiesOfNewElement();
        void calcSurfPropertiesOfElement(int n);
        void refreshOwned(int setupFlag);
        void refreshGhosts(int setupFlag);
        inline void recalcLocalSurfProperties();
        inline void recalcGhostSurfProperties();

        virtual double calcArea(int nElem) = 0;
        void calcEdgeVecLen(int nElem, double *len, double **vec);
        void calcSurfaceNorm(int nElem, double *surfNorm);
        void calcEdgeNormals(int nElem, double **edgeNorm);
        void calcObtuseAngleIndex(int nElem, int iNode, double &dot);

        virtual bool isInElement(double *pos, int i) = 0;

        int randomOwnedGhostElement();

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
        inline int& obtuseAngleIndex(int i) {return obtuseAngleIndex_(i); }

        inline int nBelowAngleSoftLimit()   {return nBelowAngle_softLimit_;}
        inline double angleSoftLimit()      {return acos(minAngle_softLimit_)*180./M_PI;}
        inline double angleHardLimit()      {return acos(minAngle_hardLimit_)*180./M_PI;}
        inline int nTooManyNeighs()         {return nTooManyNeighs_;}
        inline int nOverlapping()           {return nOverlapping_;}

      private:

        int searchElementByAreaAcc(double area,int lo, int hi);

        void growSurface(int iSrf, double by = 1e-13);

        void extrudeEdge(const int nElem, const int edge, const double * const extrudeVec, int &count, double * extrusion_tri_nodes);

        // mesh properties
        double curvature_;
        bool curvature_tolerant_;
        double minAngle_softLimit_;         
        double minAngle_hardLimit_;         
        ScalarContainer<double>& areaMesh_; 

        int nBelowAngle_softLimit_;
        int nTooManyNeighs_;
        int nOverlapping_;

        // per-element properties
        // add more properties as they are needed, such as
        // edge length, edge vectors, surface normal ....

        // surface properties
        ScalarContainer<double>& area_;
        ScalarContainer<double>& areaAcc_;                    // accumulated area of owned and ghost particles
        VectorContainer<double,NUM_NODES>& edgeLen_;          // len of edgeVec
        MultiVectorContainer<double,NUM_NODES,3>& edgeVec_;   // unit vec from node0 to node1, node1 to node2 etc
        MultiVectorContainer<double,NUM_NODES,3>& edgeNorm_;
        VectorContainer<double,3>& surfaceNorm_;
        ScalarContainer<int>& obtuseAngleIndex_;

        // neighbor topology
        
        ScalarContainer<int>& nNeighs_;
        VectorContainer<int,NUM_NEIGH_MAX>& neighFaces_;
        VectorContainer<bool,NUM_NODES>& hasNonCoplanarSharedNode_;
        VectorContainer<bool,NUM_NODES>& edgeActive_;
        VectorContainer<bool,NUM_NODES>& cornerActive_;

        // for overlap check on element insertion
        
        RegionNeighborList<interpolate_no> &neighList_;
};

// *************************************
#include "surface_mesh_I.h"

#include "surface_mesh_feature_remove_I.h"

// *************************************

} /* LAMMPS_NS */
#endif /* SURFACEMESH_H_ */
