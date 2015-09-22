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

#ifndef LMP_SURFACE_MESH_I_H
#define LMP_SURFACE_MESH_I_H

#define NTRY_MC_SURFACE_MESH_I_H 30000
#define NITER_MC_SURFACE_MESH_I_H 5
#define TOLERANCE_MC_SURFACE_MESH_I_H 0.05

/* ----------------------------------------------------------------------
   constructors, destructor
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::SurfaceMesh(LAMMPS *lmp)
:   TrackingMesh<NUM_NODES>(lmp),
    curvature_(1.-EPSILON_CURVATURE),
    curvature_tolerant_(false),
    minAngle_(cos(MIN_ANGLE_MESH*M_PI/180.)),

    // TODO should keep areaMeshSubdomain up-to-date more often for insertion faces
    areaMesh_     (*this->prop().template addGlobalProperty   < ScalarContainer<double> >                 ("areaMesh",     "comm_none","frame_trans_rot_invariant","restart_no",2)),

    nBelowAngle_(0),
    nTooManyNeighs_(0),
    nOverlapping_(0),

    area_         (*this->prop().template addElementProperty< ScalarContainer<double> >                   ("area",         "comm_none","frame_trans_rot_invariant", "restart_no",2)),
    areaAcc_      (*this->prop().template addElementProperty< ScalarContainer<double> >                   ("areaAcc",      "comm_none","frame_trans_rot_invariant", "restart_no",2)),
    edgeLen_      (*this->prop().template addElementProperty< VectorContainer<double,NUM_NODES> >         ("edgeLen",      "comm_none","frame_trans_rot_invariant", "restart_no")),
    edgeVec_      (*this->prop().template addElementProperty< MultiVectorContainer<double,NUM_NODES,3> >  ("edgeVec",      "comm_none","frame_scale_trans_invariant","restart_no")),
    edgeNorm_     (*this->prop().template addElementProperty< MultiVectorContainer<double,NUM_NODES,3> >  ("edgeNorm",     "comm_none","frame_scale_trans_invariant","restart_no")),
    surfaceNorm_  (*this->prop().template addElementProperty< VectorContainer<double,3> >                 ("surfaceNorm",  "comm_none","frame_scale_trans_invariant","restart_no")),
    obtuseAngleIndex_   (*this->prop().template addElementProperty< ScalarContainer<int> >                ("obtuseAngleIndex","comm_exchange_borders","frame_invariant","restart_no")),
    nNeighs_      (*this->prop().template addElementProperty< ScalarContainer<int> >                      ("nNeighs",      "comm_exchange_borders","frame_invariant","restart_no")),
    neighFaces_   (*this->prop().template addElementProperty< VectorContainer<int,NUM_NEIGH_MAX> >        ("neighFaces",   "comm_exchange_borders","frame_invariant","restart_no")),
    hasNonCoplanarSharedNode_(*this->prop().template addElementProperty< VectorContainer<bool,NUM_NODES> >("hasNonCoplanarSharedNode","comm_exchange_borders","frame_invariant", "restart_no")),
    edgeActive_   (*this->prop().template addElementProperty< VectorContainer<bool,NUM_NODES> >           ("edgeActive",   "comm_exchange_borders","frame_invariant","restart_no")),
    cornerActive_ (*this->prop().template addElementProperty< VectorContainer<bool,NUM_NODES> >           ("cornerActive", "comm_exchange_borders","frame_invariant","restart_no"))

{
    
    areaMesh_.add(0.);
    areaMesh_.add(0.);
    areaMesh_.add(0.);
    
}

template<int NUM_NODES, int NUM_NEIGH_MAX>
SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::~SurfaceMesh()
{}

/* ----------------------------------------------------------------------
   set mesh curvature, used for mesh topology
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::setCurvature(double _curvature)
{
    curvature_ = _curvature;
}

/* ----------------------------------------------------------------------
   set mesh curvature tolerance
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::setCurvatureTolerant(bool _tol)
{
    curvature_tolerant_ = _tol;
}

/* ----------------------------------------------------------------------
   add and delete an element
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
bool SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::addElement(double **nodeToAdd,int lineNumb)
{
    if(TrackingMesh<NUM_NODES>::addElement(nodeToAdd,lineNumb))
    {

        calcSurfPropertiesOfNewElement();
        return true;
    }
    return false;
}

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::deleteElement(int n)
{
    TrackingMesh<NUM_NODES>::deleteElement(n);
}

/* ----------------------------------------------------------------------
   recalculate properties on setup (on start and during simulation)
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::refreshOwned(int setupFlag)
{
    TrackingMesh<NUM_NODES>::refreshOwned(setupFlag);
    // (re)calculate all properties for owned elements
    
    recalcLocalSurfProperties();
}

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::refreshGhosts(int setupFlag)
{
    TrackingMesh<NUM_NODES>::refreshGhosts(setupFlag);

    recalcGhostSurfProperties();
}

/* ----------------------------------------------------------------------
   recalculate properties of local elements
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::recalcLocalSurfProperties()
{
    
    // areaMeshGlobal [areaMesh_(0)] and areaMeshOwned [areaMesh_(1)]
    // calculated here

    areaMesh_(0) = 0.;
    areaMesh_(1) = 0.;

    int nlocal = this->sizeLocal();

    for(int i = 0; i < nlocal; i++)
    {
      calcEdgeVecLen(i, edgeLen(i), edgeVec(i));
      calcSurfaceNorm(i, surfaceNorm(i));
      calcEdgeNormals(i, edgeNorm(i));
      for(int j=0;j<NUM_NODES;j++)
      {
          double dot;
          calcObtuseAngleIndex(i,j,dot);
      }

      area(i) = calcArea(i);
      areaAcc(i) = area(i);
      if(i > 0) areaAcc(i) += areaAcc(i-1);

      // add to local area
      areaMesh_(1) += area(i);
      
    }

    // mesh area must be summed up
    MPI_Sum_Scalar(areaMesh_(1),areaMesh_(0),this->world);

}

/* ----------------------------------------------------------------------
   recalculate properties of ghost elements
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::recalcGhostSurfProperties()
{
    int nlocal = this->sizeLocal();
    int nall = this->sizeLocal()+this->sizeGhost();

    // areaMeshGhost [areaMesh_(2)] calculated here

    // accumulated area includes owned and ghosts
    areaMesh_(2) = 0.;
    for(int i = nlocal; i < nall; i++)
    {
      calcEdgeVecLen(i, edgeLen(i), edgeVec(i));
      calcSurfaceNorm(i, surfaceNorm(i));
      calcEdgeNormals(i, edgeNorm(i));

      for(int j=0;j<NUM_NODES;j++)
      {
          double dot;
          calcObtuseAngleIndex(i,j,dot);
      }

      area(i) = calcArea(i);
      areaAcc(i) = area(i);
      if(i > 0) areaAcc(i) += areaAcc(i-1);

      // add to ghost area
      areaMesh_(2) += area(i);
    }

}

/* ----------------------------------------------------------------------
   generate a random Element by areaAcc
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
inline int SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::randomOwnedGhostElement()
{
    
    if(!this->isInsertionMesh())
        this->error->one(FLERR,"Illegal call for non-insertion mesh");

    double area = areaMeshOwned()+areaMeshGhost();

    double r = this->random_->uniform() * area;
    
    int first = 0;
    int last = this->sizeLocal()+this->sizeGhost()-1;

    return searchElementByAreaAcc(r,first,last);
}

template<int NUM_NODES, int NUM_NEIGH_MAX>
inline int SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::searchElementByAreaAcc(double area,int lo, int hi)
{
    
    if( (lo < 1 || area > areaAcc(lo-1)) && (area <= areaAcc(lo)) )
        return lo;
    if( (hi < 1 || area > areaAcc(hi-1)) && (area <= areaAcc(hi)) )
        return hi;

    int mid = static_cast<int>((lo+hi)/2);
    if(area > areaAcc(mid))
        return searchElementByAreaAcc(area,mid,hi);
    else
        return searchElementByAreaAcc(area,lo,mid);
}

/* ----------------------------------------------------------------------
   calculate surface properties of new element
   only called once on import
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::calcSurfPropertiesOfNewElement()
{
    
    int n = this->sizeLocal()-1;

    double *vecTmp3,*vecTmpNumNodes,**nodeTmp;
    create<double>(vecTmp3,3);
    create<double>(vecTmpNumNodes,NUM_NODES);
    create<double>(nodeTmp,NUM_NODES,3);

    // calculate edge vectors and lengths
    calcEdgeVecLen(n,vecTmpNumNodes,nodeTmp);
    edgeLen_.set(n,vecTmpNumNodes);
    edgeVec_.set(n,nodeTmp);

    // calc surface normal
    calcSurfaceNorm(n,vecTmp3);
    surfaceNorm_.set(n,vecTmp3);

    // calc edge normal in plane pointing outwards of area_
    // should be (edgeVec_ cross surfaceNormal)
    calcEdgeNormals(n,nodeTmp);
    edgeNorm_.set(n,nodeTmp);

    obtuseAngleIndex_.set(n,NO_OBTUSE_ANGLE);

    bool hasSmallAngle = false;

    for(int i=0;i<NUM_NODES;i++){
      double dot;
      calcObtuseAngleIndex(n,i,dot);
      if(-dot > minAngle_)
        hasSmallAngle = true;
    }

    if(hasSmallAngle)
    {
        if(TrackingMesh<NUM_NODES>::verbose() && 0 == this->comm->me)
            fprintf(this->screen,"Mesh %s: elements %d (line %d) has high aspect ratio (angle < %f °) \n",
                    this->mesh_id_,n,TrackingMesh<NUM_NODES>::lineNo(n),this->angleLimit());
        nBelowAngle_++;
    }

    // calc area_ from previously obtained values and add to container
    // calcArea is pure virtual and implemented in derived class(es)
    
    double area_elem = calcArea(n);
    areaMesh_(0) += area_elem;
    area_(n) = area_elem;
    areaAcc_(n) = area_elem;
    if(n > 0) areaAcc_(n) += areaAcc_(n-1);

    // cannot calc areaMesh_(1), areaMesh_(2), areaMesh_(3) here since
    // not parallelized at this point
    
    destroy<double>(nodeTmp);
    destroy<double>(vecTmpNumNodes);
    destroy<double>(vecTmp3);

}

/* ----------------------------------------------------------------------
   sub-functions needed to calculate mesh properties
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::calcEdgeVecLen(int nElem, double *len, double **vec)
{
    for(int i=0;i<NUM_NODES;i++)
    {
      vectorSubtract3D(
        MultiNodeMesh<NUM_NODES>::node_(nElem)[(i+1)%NUM_NODES],
        MultiNodeMesh<NUM_NODES>::node_(nElem)[i],vec[i]);
      len[i] = vectorMag3D(vec[i]);
      vectorScalarDiv3D(vec[i],len[i]);
    }
}

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::calcSurfaceNorm(int nElem, double *surfNorm)
{
    vectorCross3D(edgeVec(nElem)[0],edgeVec(nElem)[1],surfNorm);
    vectorScalarDiv3D(surfNorm, vectorMag3D(surfNorm));
}

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::calcEdgeNormals(int nElem, double **edgeNorm)
{
    for(int i=0;i<NUM_NODES;i++){
      vectorCross3D(edgeVec(nElem)[i],surfaceNorm(nElem),edgeNorm[i]);
      vectorScalarDiv3D(edgeNorm[i],vectorMag3D(edgeNorm[i]));
    }
}

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::calcObtuseAngleIndex(int nElem, int iNode, double &dot)
{
    
    dot = vectorDot3D(edgeVec_(nElem)[iNode],edgeVec_(nElem)[(iNode-1+NUM_NODES)%NUM_NODES]);
    
    if(dot > 0.)
    {
        
        obtuseAngleIndex_.set(nElem,iNode);
    }
    else
        obtuseAngleIndex_.set(nElem,NO_OBTUSE_ANGLE);
}

/* ----------------------------------------------------------------------
   build neighlist, generate mesh topology, check (in)active edges and nodes
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::buildNeighbours()
{
    
    int nall = this->sizeLocal()+this->sizeGhost();

    bool t[NUM_NODES], f[NUM_NODES];
    int neighs[NUM_NEIGH_MAX];

    for(int i = 0; i < NUM_NODES; i++)
    {
        t[i] = true;
        f[i] = false;
    }
    for(int i = 0; i < NUM_NEIGH_MAX; i++)
        neighs[i] = -1;

    for(int i = 0; i < nall; i++)
    {
        nNeighs_.set(i,0);
        neighFaces_.set(i,neighs);
        edgeActive_.set(i,t);
        cornerActive_.set(i,t);
        hasNonCoplanarSharedNode_.set(i,f);
    }

    // build neigh topology and edge activity, ~n*n/2
    for(int i = 0; i < nall; i++)
    {
      for(int j = i+1; j < nall; j++)
      {
        
        int iEdge(0), jEdge(0);

        if(shareEdge(i,j,iEdge,jEdge))
          handleSharedEdge(i,iEdge,j,jEdge, areCoplanar(TrackingMesh<NUM_NODES>::id(i),TrackingMesh<NUM_NODES>::id(j)));
      }
    }

    int *idListVisited = new int[nall];
    int *idListHasNode = new int[nall];
    double **edgeList,**edgeEndPoint;
    this->memory->create(edgeList,2*nall,3,"SurfaceMesh:edgeList");
    this->memory->create(edgeEndPoint,2*nall,3,"SurfaceMesh:edgeEndPoint");

    // recursively handle corner activity, ~n
    for(int i = 0; i < nall; i++)
    {
        for(int iNode = 0; iNode < NUM_NODES; iNode++)
            handleCorner(i,iNode,idListVisited,idListHasNode,edgeList,edgeEndPoint);
    }

    delete []idListVisited;
    delete []idListHasNode;
    this->memory->destroy(edgeList);
    this->memory->destroy(edgeEndPoint);
    
    // correct edge and corner activation/deactivation in parallel
    
    parallelCorrection();
}

/* ----------------------------------------------------------------------
   quality check for surface mesh
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::qualityCheck()
{
    // iterate over surfaces
    
    int nlocal = this->sizeLocal();
    int nall = this->sizeLocal()+this->sizeGhost();
    int me = this->comm->me;

    // check duplicate elements, n^2/2 operation
    
    for(int i = 0; i < nlocal; i++)
    {
        for(int j = i+1; j < nall; j++)
        {
            if(this->nSharedNodes(i,j) == NUM_NODES)
            {
                fprintf(this->screen,"ERROR: Mesh %s: elements %d and %d (lines %d and %d) are duplicate\n",
                        this->mesh_id_,TrackingMesh<NUM_NODES>::id(i),TrackingMesh<NUM_NODES>::id(j),
                        TrackingMesh<NUM_NODES>::lineNo(i),TrackingMesh<NUM_NODES>::lineNo(j));
                if(!this->removeDuplicates())
                    this->error->one(FLERR,"Fix mesh: Bad mesh, cannot continue. You can try re-running with 'heal auto_remove_duplicates'");
                else
                    this->error->one(FLERR,"Fix mesh: Bad mesh, cannot continue. The mesh probably reached the precision you defined. "
                                           "You can try re-running with a lower value for 'precision'");
            }
        }
    }

    for(int i = 0; i < nlocal; i++)
    {
      for(int iNode = 0; iNode < NUM_NODES; iNode++)
      {
        double dot;
        calcObtuseAngleIndex(i,iNode,dot);
        if(-dot > curvature_)
        {
            fprintf(this->screen,"%s: Mesh %s: The minumum angle of mesh element %d (line %d) is lower than the specified curvature. "
                                  "Increase mesh quality or decrease curvature (currently %f°)\n",
                                  curvature_tolerant_?"WARNING:":"ERROR",this->mesh_id_,TrackingMesh<NUM_NODES>::id(i),
                                  TrackingMesh<NUM_NODES>::lineNo(i),acos(curvature_)*180./M_PI);
            if(!curvature_tolerant_)
                this->error->one(FLERR,"Fix mesh: Bad mesh, cannot continue. You can try setting 'curvature' to 1e-5 or lower or use 'curvature_tolerant yes'");
        }
      }
    }

    if(this->nBelowAngle() > 0 && 0 == me)
    {
        fprintf(this->screen,"Mesh %s: %d elements have high aspect ratio (angle < %f °)\n",
                this->mesh_id_,this->nBelowAngle(),this->angleLimit());
        this->error->warning(FLERR,"Fix mesh: Mesh contains highly skewed element, moving mesh (if used) will not parallelize well");
    }

    if(this->nTooManyNeighs() > 0 && 0 == me)
    {
        
        fprintf(this->screen,"Mesh %s: %d mesh elements have more than %d neighbors \n",
                this->mesh_id_,this->nTooManyNeighs(),NUM_NEIGH_MAX);
        this->error->one(FLERR,"Fix mesh: Bad mesh, cannot continue. Possibly corrupt elements with too many neighbors.\n"
                                "If you know what you're doing, you can try to change the definition of SurfaceMeshBase in tri_mesh.h and recompile");
    }

    if(nOverlapping() > 0)
    {
        fprintf(this->screen,"WARNING: Mesh %s: proc %d has %d element pairs that are coplanar, "
                "share an edge and overlap (but are not duplicate)\n",
                this->mesh_id_,me,nOverlapping());
    }
}

/* ----------------------------------------------------------------------
   correct edge and corner activation/deactivation in parallel
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::parallelCorrection()
{
    
    int iGlobal,iLocal;
    int mine = this->sizeLocal()+this->sizeGhost();
    int sizeGlob = this->sizeGlobal();
    int len = NUM_NODES*sizeGlob;

    int *edgea = new int[len];
    int *cornera = new int[len];
    vectorInitializeN(edgea,len,2);
    vectorInitializeN(cornera,len,2);

    for(int i = 0; i < mine; i++)
    {
        iGlobal = TrackingMesh<NUM_NODES>::id(i);

        for(int j = 0; j < NUM_NODES; j++)
        {
            edgea[iGlobal*NUM_NODES+j] = edgeActive(i)[j]?1:0;
            cornera[iGlobal*NUM_NODES+j] = cornerActive(i)[j]?1:0;
        }
    }

    MPI_Min_Vector(edgea,len,this->world);
    MPI_Min_Vector(cornera,len,this->world);

    for(int i = 0; i < sizeGlob; i++)
    {
        iLocal = this->map(i);
        if(iLocal >= 0)
        {
            for(int j = 0; j < NUM_NODES; j++)
            {
                if(edgea[i*NUM_NODES+j] == 0)
                    edgeActive(iLocal)[j] = false;
                else if(edgea[i*NUM_NODES+j] == 1)
                    edgeActive(iLocal)[j] = true;
                else
                    this->error->one(FLERR,"Illegal situation in SurfaceMesh::parallelCorrection()");
                if(cornera[i*NUM_NODES+j] == 0)
                    cornerActive(iLocal)[j] = false;
                else if(cornera[i*NUM_NODES+j] == 1)
                    cornerActive(iLocal)[j] = true;
                else
                    this->error->one(FLERR,"Illegal situation in SurfaceMesh::parallelCorrection()");
            }
        }
    }

    delete []edgea;
    delete []cornera;
}

/* ----------------------------------------------------------------------
   functions to generate mesh topology
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
bool SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::areCoplanar(int tag_a, int tag_b)
{
    int a = this->map(tag_a);
    int b = this->map(tag_b);

    if(a < 0 || b < 0)
        this->error->one(FLERR,"Internal error: Illegal call to SurfaceMesh::areCoplanar()");

    // check if two faces are coplanar
    
    double dot = vectorDot3D(surfaceNorm(a),surfaceNorm(b));
    
    // need fabs in case surface normal is other direction
    if(fabs(dot) >= curvature_) return true;
    else return false;
}

/* ---------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
bool SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::areCoplanarNeighs(int tag_a, int tag_b)
{
    bool areNeighs = false;
    int a = this->map(tag_a);
    int b = this->map(tag_b);

    if(a < 0 || b < 0)
        this->error->one(FLERR,"Internal error: Illegal call to SurfaceMesh::areCoplanarNeighs()");

    // check if two faces are coplanar
    
    // must be neighs, otherwise not considered coplanar
    for(int i = 0; i < nNeighs_(a); i++)
        if(neighFaces_(a)[i] == tag_b)
            areNeighs = true;

    if(!areNeighs) return false;

    double dot = vectorDot3D(surfaceNorm(a),surfaceNorm(b));
    
    // need fabs in case surface normal is other direction
    if(fabs(dot) > curvature_) return true;
    else return false;
}

/* ---------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
bool SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::areCoplanarNodeNeighs(int tag_a, int tag_b)
{
    bool areNeighs = false;
    int a = this->map(tag_a);
    int b = this->map(tag_b);

    if(a < 0 || b < 0)
        this->error->one(FLERR,"Internal error: Illegal call to SurfaceMesh::areCoplanarNeighs()");

    // check if two faces are coplanar
    
    // must be neighs, otherwise not considered coplanar
    for(int i = 0; i < nNeighs_(a); i++)
        if(neighFaces_(a)[i] == tag_b)
            areNeighs = true;

    if(!areNeighs && MultiNodeMesh<NUM_NODES>::nSharedNodes(a,b) == 0) return false;

    double dot = vectorDot3D(surfaceNorm(a),surfaceNorm(b));
    
    // need fabs in case surface normal is other direction
    if(fabs(dot) > curvature_) return true;
    else return false;
}

/* ---------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
bool SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::coplanarNeighsOverlap(int iSrf,int iEdge,int jSrf,int jEdge)
{
    
    double vecI[3],vecJ[3], pRef[3], edgeN[3], dot1, dot2;

    vectorCopy3D(MultiNodeMesh<NUM_NODES>::node_(iSrf)[iEdge],pRef);
    vectorCopy3D(edgeNorm(iSrf)[iEdge],edgeN);

    vectorSubtract3D(MultiNodeMesh<NUM_NODES>::node_(iSrf)[(iEdge+2)%NUM_NODES],pRef,vecI);
    vectorSubtract3D(MultiNodeMesh<NUM_NODES>::node_(jSrf)[(jEdge+2)%NUM_NODES],pRef,vecJ);

    dot1 = vectorDot3D(vecI,edgeN);
    dot2 = vectorDot3D(vecJ,edgeN);

    if(dot1*dot2 > 0.)
    {
        if(TrackingMesh<NUM_NODES>::verbose())
        {
            
            int nlocal = this->sizeLocal();
            fprintf(this->screen,"WARNING: Mesh %s: elements %d and %d are coplanar, "
                    "share an edge and overlap (but are not duplicate)\n",
                    this->mesh_id_,TrackingMesh<NUM_NODES>::id(iSrf),TrackingMesh<NUM_NODES>::id(jSrf));
            if(iSrf < nlocal)
                fprintf(this->screen,"INFO: Mesh %s: element %d corresponds to line # %d\n",
                    this->mesh_id_,TrackingMesh<NUM_NODES>::id(iSrf),TrackingMesh<NUM_NODES>::lineNo(iSrf));
            if(jSrf < nlocal)
                fprintf(this->screen,"INFO: Mesh %s: element %d corresponds to line # %d\n",
                    this->mesh_id_,TrackingMesh<NUM_NODES>::id(jSrf),TrackingMesh<NUM_NODES>::lineNo(jSrf));
        }

        nOverlapping_++;
        
        //this->error->warning(FLERR,"Fix mesh: Check overlapping mesh elements");
        return true;
    }
    else return false;
}

/* ---------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
bool SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::edgeVecsColinear(double *v,double *w)
{
    // need normalized vectors
    double dot = vectorDot3D(v,w);
    // need fabs in case vectors are in different direction
    if(fabs(dot) > curvature_) return true;
    else return false;
}

/* ---------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::growSurface(int iSrf, double by)
{
    double *tmp = new double[3];
    for(int i=0;i<NUM_NODES;i++)
    {
      vectorSubtract3D(MultiNodeMesh<NUM_NODES>::node(iSrf)[i],this->center_(iSrf),tmp);
      vectorScalarMult3D(tmp,by);
      vectorAdd3D(MultiNodeMesh<NUM_NODES>::node(iSrf)[i],
                      tmp,MultiNodeMesh<NUM_NODES>::node(iSrf)[i]);
    }
    delete[] tmp;
    return;
}

/* ---------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
bool SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::shareEdge(int iSrf, int jSrf, int &iEdge, int &jEdge)
{
    int iNode1=0,jNode1=0,iNode2,jNode2;
    if(this->share2Nodes(iSrf,jSrf,iNode1,jNode1,iNode2,jNode2)){
      // following implementation of shareNode(), the only remaining option to
      // share an edge is that the next node of iSrf is equal to the next or previous
      // node if jSrf
      
      if(2 == iNode1+iNode2)
        iEdge = 2;
      
      else
        iEdge = MathExtraLiggghts::min(iNode1,iNode2);

      if(2 == jNode1+jNode2)
        jEdge = 2;
      else
        jEdge = MathExtraLiggghts::min(jNode1,jNode2);

      return true;
    }
    
    iEdge = -1; jEdge = -1;
    return false;
}

/* ---------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::handleSharedEdge(int iSrf, int iEdge, int jSrf, int jEdge,
                                            bool coplanar, bool neighflag)
{
    
    if(neighflag)
    {
        
        if(nNeighs_(iSrf) == NUM_NEIGH_MAX || nNeighs_(jSrf) == NUM_NEIGH_MAX)
        {
            int ii;
            if(nNeighs_(iSrf) == NUM_NEIGH_MAX)
                ii = iSrf;
            else
                ii = jSrf;

            nTooManyNeighs_++;
            fprintf(this->screen,"Mesh %s: element id %d (line %d) has %d neighs, but only %d expected\n",
                    this->mesh_id_,TrackingMesh<NUM_NODES>::id(ii),TrackingMesh<NUM_NODES>::lineNo(ii),nNeighs_(ii)+1,NUM_NEIGH_MAX);
            if(MultiNodeMesh<NUM_NODES>::elementExclusionList())
                fprintf(MultiNodeMesh<NUM_NODES>::elementExclusionList(),"%d\n",TrackingMesh<NUM_NODES>::lineNo(ii));
            
        }

        // set neighbor topology
        if(nNeighs_(iSrf) < NUM_NEIGH_MAX)
            neighFaces_(iSrf)[nNeighs_(iSrf)] = TrackingMesh<NUM_NODES>::id(jSrf);
        if(nNeighs_(jSrf) < NUM_NEIGH_MAX)
            neighFaces_(jSrf)[nNeighs_(jSrf)] = TrackingMesh<NUM_NODES>::id(iSrf);
        nNeighs_(iSrf)++;
        nNeighs_(jSrf)++;
    }

    // deactivate one egde
    // other as well if coplanar
    
    if(!coplanar || coplanarNeighsOverlap(iSrf,iEdge,jSrf,jEdge))
    {
        if(TrackingMesh<NUM_NODES>::id(iSrf) < TrackingMesh<NUM_NODES>::id(jSrf))
        {
            
            edgeActive(iSrf)[iEdge] = false;
            edgeActive(jSrf)[jEdge] = true;
        }
        else
        {
            
            edgeActive(iSrf)[iEdge] = true;
            edgeActive(jSrf)[jEdge] = false;
        }
    }
    else // coplanar
    {
        if(!coplanar) this->error->one(FLERR,"internal error");
        
        edgeActive(iSrf)[iEdge] = false;
        edgeActive(jSrf)[jEdge] = false;
    }
}

/* ---------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
int SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::handleCorner(int iSrf, int iNode,
        int *idListVisited,int *idListHasNode,double **edgeList,double **edgeEndPoint)
{
    double nodeToCheck[3];
    bool hasTwoColinearEdges, anyActiveEdge;
    int nIdListVisited = 0, nIdListHasNode = 0, maxId = -1, nEdgeList;

    this->node(iSrf,iNode,nodeToCheck);
    anyActiveEdge = false;
    checkNodeRecursive(iSrf,nodeToCheck,nIdListVisited,idListVisited,
        nIdListHasNode,idListHasNode,edgeList,edgeEndPoint,anyActiveEdge);

    // each element that shares the node contributes two edges
    nEdgeList = 2*nIdListHasNode;

    // get max ID
    for(int i = 0; i < nIdListHasNode; i++)
        maxId = MathExtraLiggghts::max(maxId,idListHasNode[i]);

    // check if any 2 edges coplanar
    
    hasTwoColinearEdges = false;
    for(int i = 0; i < nEdgeList; i++)
    {
        for(int j = i+1; j < nEdgeList; j++)
        {
            
            if(edgeVecsColinear(edgeList[i],edgeList[j]) && !this->nodesAreEqual(edgeEndPoint[i],edgeEndPoint[j]))
                hasTwoColinearEdges = true;
        }
    }

    if(hasTwoColinearEdges || !anyActiveEdge)
        cornerActive(iSrf)[iNode] = false;
    
    else if(TrackingMesh<NUM_NODES>::id(iSrf) == maxId)
        cornerActive(iSrf)[iNode] = true;
    else
        cornerActive(iSrf)[iNode] = false;

    return nIdListHasNode;
}

/* ---------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::checkNodeRecursive(int iSrf,double *nodeToCheck,
        int &nIdListVisited,int *idListVisited,int &nIdListHasNode,int *idListHasNode,
        double **edgeList,double **edgeEndPoint,bool &anyActiveEdge)
{
    int idNeigh, iNeigh, nEdgeList = 2*nIdListHasNode, nEdgeEndPoint = 2*nIdListHasNode;

    // check if I have been here already
    for(int i = 0; i < nIdListVisited; i++)
        if(idListVisited[i] == TrackingMesh<NUM_NODES>::id(iSrf)) return;

    // add to visited list
    idListVisited[nIdListVisited++] = TrackingMesh<NUM_NODES>::id(iSrf);

    // if contains node, add to list and call neighbors
    int iNode = this->containsNode(iSrf, nodeToCheck);
    if(iNode >= 0)
    {
        
        idListHasNode[nIdListHasNode++] = TrackingMesh<NUM_NODES>::id(iSrf);
        // node iNode is associated with edge iNode and iNode-1
        vectorCopy3D(edgeVec(iSrf)[iNode],edgeList[nEdgeList++]);
        vectorCopy3D(edgeVec(iSrf)[(iNode-1+NUM_NODES)%NUM_NODES],edgeList[nEdgeList++]);
        vectorCopy3D(this->node_(iSrf)[(iNode+1)%NUM_NODES],edgeEndPoint[nEdgeEndPoint++]);
        vectorCopy3D(this->node_(iSrf)[(iNode-1+NUM_NODES)%NUM_NODES],edgeEndPoint[nEdgeEndPoint++]);
        if(edgeActive(iSrf)[iNode]) anyActiveEdge = true;
        else if(edgeActive(iSrf)[(iNode-1+NUM_NODES)%NUM_NODES]) anyActiveEdge = true;

        // only call recursive if have neighbor and if I have data of neigh element (own or ghost)

        for(int iN = 0; iN < nNeighs_(iSrf); iN++)
        {
            idNeigh = neighFaces_(iSrf)[iN];
            if(idNeigh < 0) return;
            iNeigh = this->map(idNeigh);
            if(iNeigh >= 0)
                checkNodeRecursive(iNeigh,nodeToCheck,nIdListVisited,idListVisited,nIdListHasNode,
                                    idListHasNode,edgeList,edgeEndPoint,anyActiveEdge);
        }
    }
    
}

/* ----------------------------------------------------------------------
   move mesh
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::move(double *vecTotal, double *vecIncremental)
{
    TrackingMesh<NUM_NODES>::move(vecTotal,vecIncremental);
}

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::move(double *vecIncremental)
{
    TrackingMesh<NUM_NODES>::move(vecIncremental);
}

/* ----------------------------------------------------------------------
   scale mesh
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::scale(double factor)
{
    TrackingMesh<NUM_NODES>::scale(factor);

}

/* ----------------------------------------------------------------------
   rotate mesh
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::rotate(double *totalQ, double *dQ,double *origin)
{
    TrackingMesh<NUM_NODES>::rotate(totalQ,dQ,origin);

    // find out if rotating every property is cheaper than
    // re-calculating them from the new nodes
    
}

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::rotate(double *dQ,double *origin)
{
    TrackingMesh<NUM_NODES>::rotate(dQ,origin);

    // find out if rotating every property is cheaper than
    // re-calculating them from the new nodes
    
}

/* ----------------------------------------------------------------------
   check if faces is planar
   used to check if a face can be used for particle insertion
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
bool SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::isPlanar()
{
    int id_j;
    int flag = 0;

    int nlocal = this->sizeLocal();

    for(int i = 0; i < nlocal; i++)
    {
        if(flag) break;

        for(int ineigh = 0; ineigh < nNeighs_(i); ineigh++)
        {
            id_j = neighFaces_(i)[ineigh];
            if(!areCoplanarNeighs(TrackingMesh<NUM_NODES>::id(i),id_j))
                flag = 1;
        }
    }

    MPI_Max_Scalar(flag,this->world);

    if(flag) return false;
    return true;
}

/* ----------------------------------------------------------------------
   check if point on surface - only valid if pos is in my subbox
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
bool SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::isOnSurface(double *pos)
{
    bool on_surf = false;

    int nall = this->sizeLocal()+this->sizeGhost();

    // brute force
    // loop over ghosts as well as they might overlap my subbox
    for(int i = 0; i < nall; i++)
    {
        on_surf = on_surf || isInElement(pos,i);
        if(on_surf) break;
    }

    return on_surf;
}

/* ----------------------------------------------------------------------
   return number of active edges and corners for debugging
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
int SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::n_active_edges(int i)
{
    int n = 0;
    if(i > this->size()) return n;

    if(edgeActive(i)[0]) n++;
    if(edgeActive(i)[1]) n++;
    if(edgeActive(i)[2]) n++;
    return n;
}

template<int NUM_NODES, int NUM_NEIGH_MAX>
int SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::n_active_corners(int i)
{
    int n = 0;
    if(i > this->size()) return n;

    if(cornerActive(i)[0]) n++;
    if(cornerActive(i)[1]) n++;
    if(cornerActive(i)[2]) n++;
    return n;
}

/* ----------------------------------------------------------------------
   edge-edge, edge-node, edge-point distance
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
double SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::edgeEdgeDist(int iSrf, int iEdge, int jSrf, int jEdge)
{
    double d1,d2,d3,d4;
    d1 = edgeNodeDist(iSrf,iEdge,jSrf,jEdge);
    d2 = edgeNodeDist(iSrf,iEdge,jSrf,(jEdge+1)%NUM_NODES);
    d3 = edgeNodeDist(jSrf,jEdge,iSrf,(iEdge+1)%NUM_NODES);
    d4 = edgeNodeDist(jSrf,jEdge,iSrf,(iEdge+1)%NUM_NODES);
    return MathExtraLiggghts::min(d1,d2,d3,d4);
}

template<int NUM_NODES, int NUM_NEIGH_MAX>
double SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::edgeNodeDist(int iSrf, int iEdge, int jSrf, int jNode)
{
    return edgePointDist(iSrf, iEdge, MultiNodeMesh<NUM_NODES>::node_(jSrf)[jNode]);
}

template<int NUM_NODES, int NUM_NEIGH_MAX>
double SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::edgePointDist(int iSrf, int iEdge, double *point)
{
        double nodeToP[3], dot;

        vectorSubtract3D(point,MultiNodeMesh<NUM_NODES>::node_(iSrf)[iEdge],nodeToP);
        dot = vectorDot3D(edgeVec(iSrf)[iEdge],nodeToP);

        if(dot < 0)
            return vectorMag3D(nodeToP);
        
        else if(dot > edgeLen(iSrf)[iEdge])
        {
            vectorSubtract3D(point,MultiNodeMesh<NUM_NODES>::node_(iSrf)[(iEdge+1)%NUM_NODES],nodeToP);
            return vectorMag3D(nodeToP);
        }
        
        else
            return MathExtraLiggghts::abs(vectorDot3D(edgeNorm(iSrf)[iEdge],nodeToP));
}

#endif
