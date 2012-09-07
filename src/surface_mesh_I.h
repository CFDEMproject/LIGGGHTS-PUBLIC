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

#ifndef LMP_SURFACE_MESH_I_H
#define LMP_SURFACE_MESH_I_H

#define NTRY_MC_SURFACE_MESH_I_H 30000
#define NITER_MC_SURFACE_MESH_I_H 5
#define TOLERANCE_MC_SURFACE_MESH_I_H 0.05

/* ----------------------------------------------------------------------
   constructors, destructor
------------------------------------------------------------------------- */

template<int NUM_NODES>
SurfaceMesh<NUM_NODES>::SurfaceMesh(LAMMPS *lmp)
:   TrackingMesh<NUM_NODES>(lmp),
    isInsertionMesh_(false),
    isShallowGlobalMesh_(false),
    curvature_(1.-EPSILON_CURVATURE),

    // TODO should keep areaMeshSubdomain up-to-date more often for insertion faces
    areaMesh_     (*this->prop().template addGlobalProperty   < ScalarContainer<double> >                 ("areaMesh",     "comm_none","frame_trans_rot_invariant","restart_no",2)),

    area_         (*this->prop().template addElementProperty< ScalarContainer<double> >                   ("area",         "comm_none","frame_trans_rot_invariant", "restart_no",2)),
    areaAcc_      (*this->prop().template addElementProperty< ScalarContainer<double> >                   ("areaAcc",      "comm_none","frame_trans_rot_invariant", "restart_no",2)),
    edgeLen_      (*this->prop().template addElementProperty< VectorContainer<double,NUM_NODES> >         ("edgeLen",      "comm_none","frame_trans_rot_invariant", "restart_no")),
    edgeVec_      (*this->prop().template addElementProperty< MultiVectorContainer<double,NUM_NODES,3> >  ("edgeVec",      "comm_none","frame_scale_trans_invariant","restart_no")),
    edgeNorm_     (*this->prop().template addElementProperty< MultiVectorContainer<double,NUM_NODES,3> >  ("edgeNorm",     "comm_none","frame_scale_trans_invariant","restart_no")),
    surfaceNorm_  (*this->prop().template addElementProperty< VectorContainer<double,3> >                 ("surfaceNorm",  "comm_none","frame_scale_trans_invariant","restart_no")),
    edgeActive_   (*this->prop().template addElementProperty< VectorContainer<bool,NUM_NODES> >           ("edgeActive",   "comm_exchange_borders","frame_invariant",            "restart_no")),
    cornerActive_ (*this->prop().template addElementProperty< VectorContainer<bool,NUM_NODES> >           ("cornerActive", "comm_exchange_borders","frame_invariant",            "restart_no")),
    hasNonCoplanarSharedNode_(*this->prop().template addElementProperty< VectorContainer<bool,NUM_NODES> >("hasNonCoplanarSharedNode","comm_exchange_borders","frame_invariant", "restart_no")),
    nNeighs_      (*this->prop().template addElementProperty< ScalarContainer<int> >                      ("nNeighs",      "comm_exchange_borders","frame_invariant",            "restart_no")),
    
    neighFaces_   (*this->prop().template addElementProperty< VectorContainer<int,NUM_NODES> >            ("neighFaces",   "comm_exchange_borders","frame_invariant",            "restart_no")),
    obtuseAngleIndex_   (*this->prop().template addElementProperty< ScalarContainer<int> >            ("obtuseAngleIndex",   "comm_exchange_borders","frame_invariant",            "restart_no"))
{
    
    areaMesh_.add(0.);
    areaMesh_.add(0.);
    areaMesh_.add(0.);
    areaMesh_.add(0.);
    
}

template<int NUM_NODES>
SurfaceMesh<NUM_NODES>::~SurfaceMesh()
{}

/* ----------------------------------------------------------------------
   set mesh curvature, used for mesh topology
------------------------------------------------------------------------- */

template<int NUM_NODES>
void SurfaceMesh<NUM_NODES>::setCurvature(double _curvature)
 {
    curvature_ = _curvature;
}

/* ----------------------------------------------------------------------
   set flag if used as insertion mesh
------------------------------------------------------------------------- */

template<int NUM_NODES>
void SurfaceMesh<NUM_NODES>::useAsInsertionMesh()
{
    
    isInsertionMesh_ = true;
}

template<int NUM_NODES>
void SurfaceMesh<NUM_NODES>::useAsShallowGlobalMesh()
{
    
    isShallowGlobalMesh_ = true;
}

/* ----------------------------------------------------------------------
   add and delete an element
------------------------------------------------------------------------- */

template<int NUM_NODES>
void SurfaceMesh<NUM_NODES>::addElement(double **nodeToAdd)
{
    TrackingMesh<NUM_NODES>::addElement(nodeToAdd);

    calcSurfPropertiesOfNewElement();
}

template<int NUM_NODES>
void SurfaceMesh<NUM_NODES>::deleteElement(int n)
{
    TrackingMesh<NUM_NODES>::deleteElement(n);
}

/* ----------------------------------------------------------------------
   recalculate properties on setup (on start and during simulation)
------------------------------------------------------------------------- */

template<int NUM_NODES>
void SurfaceMesh<NUM_NODES>::refreshOwned(int setupFlag)
{
    TrackingMesh<NUM_NODES>::refreshOwned(setupFlag);
    // (re)calculate all properties for owned elements
    
    recalcLocalSurfProperties();
}

template<int NUM_NODES>
void SurfaceMesh<NUM_NODES>::refreshGhosts(int setupFlag)
{
    TrackingMesh<NUM_NODES>::refreshGhosts(setupFlag);

    recalcGhostSurfProperties();
}

/* ----------------------------------------------------------------------
   recalculate properties of local elements
------------------------------------------------------------------------- */

template<int NUM_NODES>
void SurfaceMesh<NUM_NODES>::recalcLocalSurfProperties()
{
    
    // areaMeshGlobal [areaMesh_(0)] and areaMeshOwned [areaMesh_(1)]
    // calculated here

    double areaAccOff;

    areaMesh_(0) = 0.;
    areaMesh_(1) = 0.;

    int nlocal = this->sizeLocal();

    for(int i = 0; i < nlocal; i++)
    {
      calcEdgeVecLen(i, edgeLen(i), edgeVec(i));
      calcSurfaceNorm(i, surfaceNorm(i));
      calcEdgeNormals(i, edgeNorm(i));
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

template<int NUM_NODES>
void SurfaceMesh<NUM_NODES>::recalcGhostSurfProperties()
{
    double pos[3], areaCheck;
    int n_succ, n_iter;
    int nlocal = this->sizeLocal();
    int nall = this->sizeLocal()+this->sizeGhost();

    // areaMeshGhost [areaMesh_(2)] and areaMeshSubdomain [areaMesh_(3)]
    // calculated here

    // accumulated area includes owned and ghosts
    areaMesh_(2) = 0.;
    for(int i = nlocal; i < nall; i++)
    {

      calcEdgeVecLen(i, edgeLen(i), edgeVec(i));
      calcSurfaceNorm(i, surfaceNorm(i));
      calcEdgeNormals(i, edgeNorm(i));
      area(i) = calcArea(i);
      areaAcc(i) = area(i);
      if(i > 0) areaAcc(i) += areaAcc(i-1);

      // add to ghost area
      areaMesh_(2) += area(i);
    }

    // calc area of owned and ghost elements in my subdomain
    
    areaMesh_(3) = 0.;
    areaCheck = 0.;

    if(isInsertionMesh_)
    {
        n_succ = 0;
        n_iter = 0;

        // iterate long enough so MC has the desired tolerance
        while( (n_iter < NITER_MC_SURFACE_MESH_I_H) &&
               (fabs((areaCheck-areaMeshGlobal()))/areaMeshGlobal() > TOLERANCE_MC_SURFACE_MESH_I_H) )
        {
            // only generate random positions if I have any mesh elements
            if(nall)
            {
                for(int i = 0; i < NTRY_MC_SURFACE_MESH_I_H; i++)
                {
                    // pick a random position on owned or ghost element
                    if((generateRandomOwnedGhost(pos)) >= 0 && (this->domain->is_in_extended_subdomain(pos)))
                        n_succ++;
                }
            }
            n_iter++;
            areaMesh_(3) = static_cast<double>(n_succ)/static_cast<double>(NTRY_MC_SURFACE_MESH_I_H*n_iter) * (areaMeshOwned()+areaMeshGhost());

            MPI_Sum_Scalar(areaMesh_(3),areaCheck,this->world);
            
        }

        if(fabs((areaCheck-areaMeshGlobal()))/areaMeshGlobal() > TOLERANCE_MC_SURFACE_MESH_I_H)
        {
            
            this->error->all(FLERR,"Local mesh area calculation failed, try boosting NITER_MC_SURFACE_MESH_I_H");
        }

        // correct so sum of all owned areas is equal to global area
        areaMesh_(3) *= areaMeshGlobal()/areaCheck;

    }

}

/* ----------------------------------------------------------------------
   recalculate some of the properties
------------------------------------------------------------------------- */

template<int NUM_NODES>
void SurfaceMesh<NUM_NODES>::recalcVectors()
{
    for(int i=0;i<this->size();i++)
    {
      calcEdgeVecLen(i, edgeLen_(i), edgeVec(i));
      calcSurfaceNorm(i, surfaceNorm(i));
      calcEdgeNormals(i, edgeNorm(i));
    }
}

/* ----------------------------------------------------------------------
   generate a random Element by areaAcc
------------------------------------------------------------------------- */

template<int NUM_NODES>
inline int SurfaceMesh<NUM_NODES>::randomOwnedGhostElement()
{
    
    if(!isInsertionMesh_ && !isShallowGlobalMesh_) this->error->one(FLERR,"Illegal call for non-insertion mesh");

    double area = isInsertionMesh_?(areaMeshOwned()+areaMeshGhost()):areaMeshGlobal();

    double r = this->random_->uniform() * area;
    
    int first = 0;
    int last = this->sizeLocal()+this->sizeGhost()-1;

    return searchElementByAreaAcc(r,first,last);
}

template<int NUM_NODES>
inline int SurfaceMesh<NUM_NODES>::searchElementByAreaAcc(double area,int lo, int hi)
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

template<int NUM_NODES>
void SurfaceMesh<NUM_NODES>::calcSurfPropertiesOfNewElement()
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

    for(int i=0;i<NUM_NODES;i++){
      
      double angle = vectorDot3D(edgeVec_(n)[i],edgeVec_(n)[(i+1)%NUM_NODES]);
      if(angle > 0.){
        obtuseAngleIndex_.set(n,i);
        break;
      }
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

template<int NUM_NODES>
void SurfaceMesh<NUM_NODES>::calcEdgeVecLen(int nElem, double *len, double **vec)
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

template<int NUM_NODES>
void SurfaceMesh<NUM_NODES>::calcEdgeLen(int nElem, double *edgeLen)
{
    for(int i=0;i<NUM_NODES;i++)
      edgeLen[i] = vectorMag3D(edgeVec(nElem)[i]);
}

template<int NUM_NODES>
void SurfaceMesh<NUM_NODES>::calcSurfaceNorm(int nElem, double *surfNorm)
{
    vectorCross3D(edgeVec(nElem)[0],edgeVec(nElem)[1],surfNorm);
    vectorScalarDiv3D(surfNorm, vectorMag3D(surfNorm));
}

template<int NUM_NODES>
void SurfaceMesh<NUM_NODES>::calcEdgeNormals(int nElem, double **edgeNorm)
{
    for(int i=0;i<NUM_NODES;i++){
      vectorCross3D(edgeVec(nElem)[i],surfaceNorm(nElem),edgeNorm[i]);
      vectorScalarDiv3D(edgeNorm[i],vectorMag3D(edgeNorm[i]));
    }
}

/* ----------------------------------------------------------------------
   build neighlist, generate mesh topology, check (in)active edges and nodes
------------------------------------------------------------------------- */

template<int NUM_NODES>
void SurfaceMesh<NUM_NODES>::buildNeighbours()
{
    // iterate over all surfaces, over ghosts as well
    int nall = this->sizeLocal()+this->sizeGhost();

    // inititalize neigh topology - reset to default, ~n
    bool t[NUM_NODES], f[NUM_NODES];
    int neighs[NUM_NODES];
    for(int i=0;i<NUM_NODES;i++)
    {
        neighs[i] = -1;
        t[i] = true;
        f[i] = false;
    }
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
        
        int iNode(0), jNode(0), iEdge(0), jEdge(0);
        if(!this->shareNode(i,j,iNode,jNode)) continue;

        if(shareEdge(i,j,iEdge,jEdge))
          handleSharedEdge(i,iEdge,j,jEdge, areCoplanar(this->id(i),this->id(j)));
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
   correct edge and corner activation/deactivation in parallel
------------------------------------------------------------------------- */

template<int NUM_NODES>
void SurfaceMesh<NUM_NODES>::parallelCorrection()
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
        iGlobal = this->id(i);

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

template<int NUM_NODES>
bool SurfaceMesh<NUM_NODES>::areCoplanar(int tag_a, int tag_b)
{
    int a = this->map(tag_a);
    int b = this->map(tag_b);

    if(a < 0 || b < 0)
        this->error->one(FLERR,"Internal error: Illegal call to SurfaceMesh::areCoplanar()");

    // check if two faces are coplanar
    // eg used to transfer shear history btw planar faces

    double dot = vectorDot3D(surfaceNorm(a),surfaceNorm(b));
    
    // need fabs in case surface normal is other direction
    if(fabs(dot) > curvature_) return true;
    else return false;
}

/* ---------------------------------------------------------------------- */

template<int NUM_NODES>
bool SurfaceMesh<NUM_NODES>::edgeVecsColinear(double *v,double *w)
{
    // need normalized vectors
    double dot = vectorDot3D(v,w);
    // need fabs in case vectors are in different direction
    if(fabs(dot) > curvature_) return true;
    else return false;
}

/* ---------------------------------------------------------------------- */

template<int NUM_NODES>
void SurfaceMesh<NUM_NODES>::growSurface(int iSrf, double by)
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

template<int NUM_NODES>
bool SurfaceMesh<NUM_NODES>::shareEdge(int iSrf, int jSrf, int &iEdge, int &jEdge)
{
    int i,j;
    if(this->shareNode(iSrf,jSrf,i,j)){
      // following implementation of shareNode(), the only remaining option to
      // share an edge is that the next node of iSrf is equal to the previous
      // node if jSrf
      if(i==0 && MultiNodeMesh<NUM_NODES>::nodesAreEqual(iSrf,NUM_NODES-1,jSrf,(j+1)%NUM_NODES)){
        iEdge = NUM_NODES-1;
        jEdge = j;
        return true;
      }
      if(MultiNodeMesh<NUM_NODES>::nodesAreEqual
            (iSrf,(i+1)%NUM_NODES,jSrf,(j-1+NUM_NODES)%NUM_NODES)){
        iEdge = i;//(ii-1+NUM_NODES)%NUM_NODES;
        jEdge = (j-1+NUM_NODES)%NUM_NODES;
        return true;
      }
    }
    iEdge = -1; jEdge = -1;
    return false;
}

/* ---------------------------------------------------------------------- */

template<int NUM_NODES>
void SurfaceMesh<NUM_NODES>::handleSharedEdge(int iSrf, int iEdge, int jSrf, int jEdge, bool coplanar)
{
    // set neighbor topology
    neighFaces_(iSrf)[nNeighs_(iSrf)] = this->id(jSrf);
    neighFaces_(jSrf)[nNeighs_(jSrf)] = this->id(iSrf);
    nNeighs_(iSrf)++;
    nNeighs_(jSrf)++;

    // deactivate one egde
    // other as well if coplanar
    
    if(coplanar)
    {
        edgeActive(iSrf)[iEdge] = false;
        edgeActive(jSrf)[jEdge] = false;
    }
    else
    {
        if(this->id(iSrf) < this->id(jSrf))
            edgeActive(iSrf)[iEdge] = false;
        else
            edgeActive(jSrf)[jEdge] = false;
    }
}

/* ---------------------------------------------------------------------- */

template<int NUM_NODES>
void SurfaceMesh<NUM_NODES>::handleCorner(int iSrf, int iNode,
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
        for(int j = i+1; j < nEdgeList; j++)
        {
            
            if(edgeVecsColinear(edgeList[i],edgeList[j]) && !this->nodesAreEqual(edgeEndPoint[i],edgeEndPoint[j]))
                hasTwoColinearEdges = true;
        }

    // deactivate all
    if(hasTwoColinearEdges || !anyActiveEdge)
        cornerActive(iSrf)[iNode] = false;
    // let the highest ID live
    else if(this->id(iSrf) == maxId)
        cornerActive(iSrf)[iNode] = true;
    else
        cornerActive(iSrf)[iNode] = false;
}

/* ---------------------------------------------------------------------- */

template<int NUM_NODES>
void SurfaceMesh<NUM_NODES>::checkNodeRecursive(int iSrf,double *nodeToCheck,
        int &nIdListVisited,int *idListVisited,int &nIdListHasNode,int *idListHasNode,
        double **edgeList,double **edgeEndPoint,bool &anyActiveEdge)
{
    int idNeigh, iNeigh, nEdgeList = 2*nIdListHasNode, nEdgeEndPoint = 2*nIdListHasNode;

    // check if I have been here already
    for(int i = 0; i < nIdListVisited; i++)
        if(idListVisited[i] == this->id(iSrf)) return;

    // add to visited list
    idListVisited[nIdListVisited++] = this->id(iSrf);

    // if contains node, add to list and call neighbors
    int iNode = this->containsNode(iSrf, nodeToCheck);
    if(iNode >= 0)
    {
        
        idListHasNode[nIdListHasNode++] = this->id(iSrf);
        // node iNode is associated with edge iNode and iNode-1
        vectorCopy3D(edgeVec(iSrf)[iNode],edgeList[nEdgeList++]);
        vectorCopy3D(edgeVec(iSrf)[(iNode-1+NUM_NODES)%NUM_NODES],edgeList[nEdgeList++]);
        vectorCopy3D(this->node_(iSrf)[(iNode+1)%NUM_NODES],edgeEndPoint[nEdgeEndPoint++]);
        vectorCopy3D(this->node_(iSrf)[(iNode-1+NUM_NODES)%NUM_NODES],edgeEndPoint[nEdgeEndPoint++]);
        if(edgeActive(iSrf)[iNode]) anyActiveEdge = true;
        else if(edgeActive(iSrf)[(iNode-1+NUM_NODES)%NUM_NODES]) anyActiveEdge = true;

        if(nNeighs_(iSrf) > NUM_NODES)
        {
            fprintf(this->screen,"mesh %s: element id %d has %d neighs, but only %d expected\n",
                    this->mesh_id_,this->id(iSrf),nNeighs_(iSrf),NUM_NODES);
            //for(int iN = 0; iN < nNeighs_(iSrf); iN++) (fprintf(this->screen,"   neigh %d\n",neighFaces_(iSrf)[iN]));
            this->error->warning(FLERR,"more mesh element neighbors than expected - corrupt mesh?");
        }

        // only call recursive if have neighbor and if I have neigh element (own or ghost)

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

template<int NUM_NODES>
void SurfaceMesh<NUM_NODES>::move(double *vecTotal, double *vecIncremental)
{
    TrackingMesh<NUM_NODES>::move(vecTotal,vecIncremental);
}

template<int NUM_NODES>
void SurfaceMesh<NUM_NODES>::move(double *vecIncremental)
{
    TrackingMesh<NUM_NODES>::move(vecIncremental);
}

/* ----------------------------------------------------------------------
   scale mesh
------------------------------------------------------------------------- */

template<int NUM_NODES>
void SurfaceMesh<NUM_NODES>::scale(double factor)
{
    TrackingMesh<NUM_NODES>::scale(factor);

}

/* ----------------------------------------------------------------------
   rotate mesh
------------------------------------------------------------------------- */

template<int NUM_NODES>
void SurfaceMesh<NUM_NODES>::rotate(double *totalQ, double *dQ,double *origin)
{
    TrackingMesh<NUM_NODES>::rotate(totalQ,dQ,origin);

    // find out if rotating every property is cheaper than
    // re-calculating them from the new nodes
    
}

template<int NUM_NODES>
void SurfaceMesh<NUM_NODES>::rotate(double *dQ,double *origin)
{
    TrackingMesh<NUM_NODES>::rotate(dQ,origin);

    // find out if rotating every property is cheaper than
    // re-calculating them from the new nodes
    
}

/* ----------------------------------------------------------------------
   check if faces is planar
   used to check if a face can be used for particle insertion
------------------------------------------------------------------------- */

template<int NUM_NODES>
bool SurfaceMesh<NUM_NODES>::isPlanar()
{
    int id_j;
    int flag = 0;

    int nlocal = this->sizeLocal();

    for(int i = 0; i < this->sizeLocal(); i++)
    {
        if(flag) break;

        for(int ineigh = 0; ineigh < nNeighs_(i); ineigh++)
        {
            id_j = neighFaces_(i)[ineigh];
            if(!areCoplanar(this->id(i),id_j))
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

template<int NUM_NODES>
bool SurfaceMesh<NUM_NODES>::isOnSurface(double *pos)
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

template<int NUM_NODES>
int SurfaceMesh<NUM_NODES>::n_active_edges(int i)
{
    int n = 0;
    if(i > this->size()) return n;

    if(edgeActive(i)[0]) n++;
    if(edgeActive(i)[1]) n++;
    if(edgeActive(i)[2]) n++;
    return n;
}

template<int NUM_NODES>
int SurfaceMesh<NUM_NODES>::n_active_corners(int i)
{
    int n = 0;
    if(i > this->size()) return n;

    if(cornerActive(i)[0]) n++;
    if(cornerActive(i)[1]) n++;
    if(cornerActive(i)[2]) n++;
    return n;
}

#endif
