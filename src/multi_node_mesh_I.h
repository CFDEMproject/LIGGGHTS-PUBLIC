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

#ifndef LMP_MULTI_NODE_MESH_I_H
#define LMP_MULTI_NODE_MESH_I_H

  /* ----------------------------------------------------------------------
   consturctors
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  MultiNodeMesh<NUM_NODES>::MultiNodeMesh(LAMMPS *lmp)
  : AbstractMesh(lmp),
    node_orig_(0),
    precision_(EPSILON_PRECISION),
    autoRemoveDuplicates_(false),
    nMove_(0),
    stepLastReset_(-1),
    nScale_(0),
    nTranslate_(0),
    nRotate_(0),
    random_(new RanPark(lmp,179424799)), // big prime #
    mesh_id_(0)
  {
  }

  /* ----------------------------------------------------------------------
   destructor
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  MultiNodeMesh<NUM_NODES>::~MultiNodeMesh()
  {
      if(node_orig_) delete node_orig_;
      delete random_;
      if(mesh_id_) delete []mesh_id_;
  }

  /* ----------------------------------------------------------------------
   set ID and mesh accuracy (latter used for mesh topology)
   set healing options
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMesh<NUM_NODES>::setMeshID(const char *_id)
  {
      if(mesh_id_) delete []mesh_id_;
      mesh_id_ = new char[strlen(_id)+1];
      strcpy(mesh_id_,_id);
  }

  template<int NUM_NODES>
  void MultiNodeMesh<NUM_NODES>::setPrecision(double _precision)
  {
      precision_ = _precision;
  }

  template<int NUM_NODES>
  void MultiNodeMesh<NUM_NODES>::autoRemoveDuplicates()
  {
      autoRemoveDuplicates_ = true;
  }

  /* ----------------------------------------------------------------------
   add an element - only called at mesh construction
   i.e. only used to construct local elements
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  bool MultiNodeMesh<NUM_NODES>::addElement(double **nodeToAdd)
  {
    
    double avg[3];

    // add node
    node_.add(nodeToAdd);

    int n = sizeLocal();

    // calculate center
    vectorZeroize3D(avg);
    for(int i = 0; i < NUM_NODES; i++)
        vectorAdd3D(nodeToAdd[i],avg,avg);
    vectorScalarDiv3D(avg,static_cast<double>(NUM_NODES));
    center_.add(avg);

    // extend bbox
    this->extendToElem(n);

    // calculate rounding radius
    double rb = 0.;
    double vec[3];
    for(int i = 0; i < NUM_NODES; i++)
    {
        vectorSubtract3D(center_(n),node_(n)[i],vec);
        rb = MathExtraLiggghts::max(rb,vectorMag3D(vec));
    }
    rBound_.add(rb);

    if(autoRemoveDuplicates_)
    {
        for(int i = 0; i < n; i++)
        {
            
            if(this->nSharedNodes(i,n) == NUM_NODES)
            {
                
                node_.del(n);
                center_.del(n);
                return false;
            }
        }
    }

    return true;
  }

  /* ----------------------------------------------------------------------
   delete an element - may delete an owned or ghost
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMesh<NUM_NODES>::deleteElement(int n)
  {
    node_.del(n);
    if(node_orig_) node_orig_->del(n);
    center_.del(n);
    rBound_.del(n);

    // do not re-calc bbox here
  }

  /* ----------------------------------------------------------------------
   recalculate properties on setup (on start and during simulation)
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMesh<NUM_NODES>::refreshOwned(int setupFlag)
  {
      int ilo = 0, ihi = sizeLocal();

      if(isDeforming())
          updateCenterRbound(ilo,ihi);

      storeNodePosRebuild();

      if(node_orig_ && setupFlag)
        storeNodePosOrig(ilo,ihi);

      // nothing more to do here, necessary initialitation done in addElement()
  }

  template<int NUM_NODES>
  void MultiNodeMesh<NUM_NODES>::refreshGhosts(int setupFlag)
  {
      int ilo = sizeLocal(), ihi = sizeLocal()+sizeGhost();

      if(isDeforming())
          updateCenterRbound(ilo,ihi);

      if(node_orig_ && setupFlag)
        storeNodePosOrig(ilo,ihi);

      // nothing more to do here, necessary initialitation done in addElement()
  }

  /* ----------------------------------------------------------------------
   comparison
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  bool MultiNodeMesh<NUM_NODES>::nodesAreEqual(int iElem, int iNode, int jElem, int jNode)
  {
    for(int i=0;i<3;i++)
      if(!MathExtraLiggghts::compDouble(node_(iElem)[iNode][i],node_(jElem)[jNode][i],precision_))
        return false;
    return true;
  }

  template<int NUM_NODES>
  bool MultiNodeMesh<NUM_NODES>::nodesAreEqual(double *nodeToCheck1,double *nodeToCheck2)
  {
    for(int i=0;i<3;i++)
      if(!MathExtraLiggghts::compDouble(nodeToCheck1[i],nodeToCheck2[i],precision_))
        return false;
    return true;
  }

  template<int NUM_NODES>
  int MultiNodeMesh<NUM_NODES>::containsNode(int iElem, double *nodeToCheck)
  {
      for(int iNode = 0; iNode < NUM_NODES; iNode++)
      {
          if(MathExtraLiggghts::compDouble(node_(iElem)[iNode][0],nodeToCheck[0],precision_) &&
             MathExtraLiggghts::compDouble(node_(iElem)[iNode][1],nodeToCheck[1],precision_) &&
             MathExtraLiggghts::compDouble(node_(iElem)[iNode][2],nodeToCheck[2],precision_))
                return iNode;
      }
      return -1;
  }

  /* ----------------------------------------------------------------------
   return the lowest iNode/jNode combination that is shared
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  bool MultiNodeMesh<NUM_NODES>::shareNode(int iElem, int jElem, int &iNode, int &jNode)
  {
    // broad phase
    double dist[3], radsum;
    vectorSubtract3D(center_(iElem),center_(jElem),dist);
    radsum = rBound_(iElem) + rBound_(jElem);
    if(vectorMag3DSquared(dist) > radsum*radsum)
    {
        iNode = -1; jNode = -1;
        
        return false;
    }

    // narrow phase
    for(int i=0;i<NUM_NODES;i++){
      for(int j=0;j<NUM_NODES;j++){
        if(MultiNodeMesh<NUM_NODES>::nodesAreEqual(iElem,i,jElem,j)){
          iNode = i; jNode = j;
          
          return true;
        }
      }
    }
    iNode = -1; jNode = -1;
    
    return false;
  }

  /* ----------------------------------------------------------------------
   return the number of shared nodes
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMesh<NUM_NODES>::nSharedNodes(int iElem, int jElem)
  {
    double dist[3], radsum;
    int nShared = 0;

    // broad phase
    vectorSubtract3D(center_(iElem),center_(jElem),dist);
    radsum = rBound_(iElem) + rBound_(jElem);
    if(vectorMag3DSquared(dist) > radsum*radsum)
        return 0;

    // narrow phase
    for(int i=0;i<NUM_NODES;i++)
      for(int j=0;j<NUM_NODES;j++)
        if(MultiNodeMesh<NUM_NODES>::nodesAreEqual(iElem,i,jElem,j))
          nShared++;

    return nShared;
  }

  /* ----------------------------------------------------------------------
   register and unregister mesh movement
   on registration, return bool staing if this is first mover on this mesh
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  bool MultiNodeMesh<NUM_NODES>::registerMove(bool _scale, bool _translate, bool _rotate)
  {
      bool isFirst = true;
      if(nMove_ > 0)
        isFirst = false;

      nMove_ ++;
      if(_scale) nScale_++;
      if(_translate) nTranslate_++;
      if(_rotate) nRotate_++;

      if(isFirst)
      {
          int nall = sizeLocal()+sizeGhost();
          
          double **tmp;
          this->memory->template create<double>(tmp,NUM_NODES,3,"MultiNodeMesh:tmp");

          if(node_orig_)
            error->one(FLERR,"Illegal situation in MultiNodeMesh<NUM_NODES>::registerMove");

          node_orig_ = new MultiVectorContainer<double,NUM_NODES,3>;
          for(int i = 0; i < nall; i++)
          {
            for(int j = 0; j < NUM_NODES; j++)
              vectorCopy3D(node_(i)[j],tmp[j]);

            node_orig_->add(tmp);
          }

          this->memory->destroy<double>(tmp);
      }

      return isFirst;
  }

  template<int NUM_NODES>
  void MultiNodeMesh<NUM_NODES>::unregisterMove(bool _scale, bool _translate, bool _rotate)
  {
      nMove_ --;
      if(_scale) nScale_--;
      if(_translate) nTranslate_--;
      if(_rotate) nRotate_--;

      bool del = true;
      if(nMove_ > 0)
        del = false;

      if(del)
      {
          delete node_orig_;
          node_orig_ = NULL;
      }
  }

  /* ----------------------------------------------------------------------
   store current node position as original node position for use by moving mesh
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMesh<NUM_NODES>::storeNodePosOrig(int ilo, int ihi)
  {
    if(!node_orig_)
        error->one(FLERR,"Internal error in MultiNodeMesh<NUM_NODES>::storeNodePosOrig");

    for(int i = ilo; i < ihi; i++)
        for(int j = 0; j < NUM_NODES; j++)
        {
            vectorCopy3D(node_(i)[j],node_orig(i)[j]);
            
        }
  }

  /* ----------------------------------------------------------------------
   reset mesh nodes to original position, done before movements are added
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  bool MultiNodeMesh<NUM_NODES>::resetToOrig()
  {
    if(!node_orig_)
        error->all(FLERR,"Internal error in MultiNodeMesh<NUM_NODES>::resetNodesToOrig");

    int ntimestep = update->ntimestep;

    if(stepLastReset_ < ntimestep)
    {
        int nall = sizeLocal() + sizeGhost();
        stepLastReset_ = ntimestep;
        for(int i = 0; i < nall; i++)
            for(int j = 0; j < NUM_NODES; j++)
                vectorCopy3D(node_orig(i)[j],node_(i)[j]);

        return true;
    }
    return false;
  }

  /* ----------------------------------------------------------------------
   move mesh by amount vecTotal, starting from original position
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMesh<NUM_NODES>::move(double *vecTotal, double *vecIncremental)
  {
    if(!isTranslating())
        this->error->all(FLERR,"Illegal call, need to register movement first");

    resetToOrig();

    int n = sizeLocal() + sizeGhost();

    for(int i = 0; i < n; i++)
    {
        vectorZeroize3D(center_(i));

        for(int j = 0; j < NUM_NODES; j++)
        {
            vectorAdd3D(node_(i)[j],vecTotal,node_(i)[j]);
            vectorAdd3D(node_(i)[j],center_(i),center_(i));
        }
        vectorScalarDiv3D(center_(i),static_cast<double>(NUM_NODES));
    }

    updateGlobalBoundingBox();
  }

  /* ----------------------------------------------------------------------
   move mesh incrementally by amount vecIncremental
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMesh<NUM_NODES>::move(double *vecIncremental)
  {
    
    int n = sizeLocal() + sizeGhost();

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < NUM_NODES; j++)
            vectorAdd3D(node_(i)[j],vecIncremental,node_(i)[j]);

        vectorAdd3D(center_(i),vecIncremental,center_(i));
    }

    updateGlobalBoundingBox();
  }
  /* ----------------------------------------------------------------------
   move mesh incrementally by amount vecIncremental
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMesh<NUM_NODES>::moveElement(int i,double *vecIncremental)
  {
    for(int j = 0; j < NUM_NODES; j++)
            vectorAdd3D(node_(i)[j],vecIncremental,node_(i)[j]);

    vectorAdd3D(center_(i),vecIncremental,center_(i));

    extendToElem(bbox_,i);
  }

  /* ----------------------------------------------------------------------
   rotate mesh interface, takes both total angle and dAngle in rad
   assumes axis stays the same over time
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMesh<NUM_NODES>::rotate(double totalAngle, double dAngle, double *axis, double *p)
  {
    double totalQ[4],dQ[4], axisNorm[3], origin[3];

    // rotates around axis through p

    // normalize axis
    vectorCopy3D(axis,axisNorm);
    vectorScalarDiv3D(axisNorm,vectorMag3D(axisNorm));

    // quat for total rotation from original position
    totalQ[0] = cos(totalAngle*0.5);
    for(int i=0;i<3;i++)
      totalQ[i+1] = axis[i]*sin(totalAngle*0.5);

    // quat for rotation since last time-step
    dQ[0] = cos(dAngle*0.5);
    for(int i = 0; i < 3; i++)
      dQ[i+1] = axis[i]*sin(dAngle*0.5);

    vectorCopy3D(p,origin);

    // apply rotation around center axis + displacement
    // = rotation around axis through p
    rotate(totalQ,dQ,origin);
  }

  template<int NUM_NODES>
  void MultiNodeMesh<NUM_NODES>::rotate(double *totalQ, double *dQ,double *origin)
  {
    if(!isRotating())
        this->error->all(FLERR,"Illegal call, need to register movement first");

    resetToOrig();

    int n = sizeLocal() + sizeGhost();

    bool trans = vectorMag3DSquared(origin) > 0.;

    // perform total rotation for data in this class
    for(int i = 0; i < n; i++)
    {
      vectorZeroize3D(center_(i));

      for(int j = 0; j < NUM_NODES; j++)
      {
        if(trans) vectorSubtract3D(node_(i)[j],origin,node_(i)[j]);
        MathExtraLiggghts::vec_quat_rotate(node_(i)[j], totalQ, node_(i)[j]);
        if(trans) vectorAdd3D(node_(i)[j],origin,node_(i)[j]);
        vectorAdd3D(node_(i)[j],center_(i),center_(i));
      }
      vectorScalarDiv3D(center_(i),static_cast<double>(NUM_NODES));
    }

    updateGlobalBoundingBox();
  }

  /* ----------------------------------------------------------------------
   rotate mesh interface, takes only dAngle
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMesh<NUM_NODES>::rotate(double dAngle, double *axis, double *p)
  {
    double dQ[4], axisNorm[3], origin[3];

    // rotates around axis through p

    // normalize axis
    vectorCopy3D(axis,axisNorm);
    vectorScalarDiv3D(axisNorm,vectorMag3D(axisNorm));

    // quat for rotation since last time-step
    dQ[0] = cos(dAngle*0.5);
    for(int i = 0; i < 3; i++)
      dQ[i+1] = axisNorm[i]*sin(dAngle*0.5);

    vectorCopy3D(p,origin);

    // apply rotation around center axis + displacement
    // = rotation around axis through p
    rotate(dQ,origin);
  }

  template<int NUM_NODES>
  void MultiNodeMesh<NUM_NODES>::rotate(double *dQ, double *origin)
  {
    
    int n = sizeLocal() + sizeGhost();

    bool trans = vectorMag3DSquared(origin) > 0.;

    // perform total rotation for data in this class
    
    for(int i = 0; i < n; i++)
    {
      vectorZeroize3D(center_(i));
      for(int j = 0; j < NUM_NODES; j++)
      {
        if(trans) vectorSubtract3D(node_(i)[j],origin,node_(i)[j]);
        MathExtraLiggghts::vec_quat_rotate(node_(i)[j], dQ,node_(i)[j]);
        if(trans) vectorAdd3D(node_(i)[j],origin,node_(i)[j]);
        vectorAdd3D(node_(i)[j],center_(i),center_(i));
      }
      vectorScalarDiv3D(center_(i),static_cast<double>(NUM_NODES));
    }

    updateGlobalBoundingBox();
  }

  /* ----------------------------------------------------------------------
   scale mesh
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMesh<NUM_NODES>::scale(double factor)
  {
    
    int n = sizeLocal() + sizeGhost();

    for(int i = 0; i < n; i++)
    {
      vectorZeroize3D(center_(i));
      for(int j = 0; j < NUM_NODES; j++)
      {
        node_(i)[j][0] *= factor;
        node_(i)[j][1] *= factor;
        node_(i)[j][2] *= factor;
        vectorAdd3D(node_(i)[j],center_(i),center_(i));
      }
      vectorScalarDiv3D(center_(i),static_cast<double>(NUM_NODES));

      // calculate rounding radius
      double rb = 0.;
      double vec[3];
      for(int j = 0; j < NUM_NODES; j++)
      {
         vectorSubtract3D(center_(i),node_(i)[j],vec);
         rb = MathExtraLiggghts::max(rb,vectorMag3D(vec));
      }
      rBound_(i) = rb;
    }

    updateGlobalBoundingBox();
  }

  /* ----------------------------------------------------------------------
   update center and rbound from node data
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMesh<NUM_NODES>::updateCenterRbound(int ilo, int ihi)
  {
    for(int i = ilo; i < ihi; i++)
    {
      vectorZeroize3D(center_(i));
      for(int j = 0; j < NUM_NODES; j++)
        vectorAdd3D(node_(i)[j],center_(i),center_(i));
      vectorScalarDiv3D(center_(i),static_cast<double>(NUM_NODES));

      // calculate rounding radius
      double rb = 0.;
      double vec[3];
      for(int j = 0; j < NUM_NODES; j++)
      {
         vectorSubtract3D(center_(i),node_(i)[j],vec);
         rb = MathExtraLiggghts::max(rb,vectorMag3D(vec));
      }
      rBound_(i) = rb;
    }

    updateGlobalBoundingBox();
  }

  /* ----------------------------------------------------------------------
   bounding box funtions
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  BoundingBox MultiNodeMesh<NUM_NODES>::getElementBoundingBoxOnSubdomain(int const n)
  {
    BoundingBox ret;
    extendToElem(ret,n);
    ret.shrinkToSubbox(this->domain->sublo,this->domain->subhi);
    return ret;
  }

  template<int NUM_NODES>
  BoundingBox MultiNodeMesh<NUM_NODES>::getGlobalBoundingBox() const
  {
    return bbox_;
  }

  template<int NUM_NODES>
  void MultiNodeMesh<NUM_NODES>::updateGlobalBoundingBox()
  {
    bbox_.reset();
    
    int n = sizeLocal();

    for(int i = 0; i < n; i++)
      extendToElem(bbox_,i);
    bbox_.extendToParallel(this->world);
  }

  template<int NUM_NODES>
  void MultiNodeMesh<NUM_NODES>::extendToElem(int const nElem)
  {
    for(int i = 0; i < NUM_NODES; ++i)
      bbox_.extendToContain(node_(nElem)[i]);
  }

  template<int NUM_NODES>
  void MultiNodeMesh<NUM_NODES>::extendToElem(BoundingBox &box, int const nElem)
  {
    for(int i = 0; i < NUM_NODES; ++i)
      box.extendToContain(node_(nElem)[i]);
  }

  /* ----------------------------------------------------------------------
   decide if any node has moved far enough to trigger re-build
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  bool MultiNodeMesh<NUM_NODES>::decideRebuild()
  {
    // just return for non-moving mesh
    if(!isMoving() && !isDeforming()) return false;

    double ***node = node_.begin();
    double ***old = nodesLastRe_.begin();
    int flag = 0;
    int nlocal = sizeLocal();
    double triggersq = 0.25*this->neighbor->skin*this->neighbor->skin;

    if(nlocal != nodesLastRe_.size())
        this->error->one(FLERR,"Internal error in MultiNodeMesh::decide_rebuild()");

    for(int iTri = 0; iTri < nlocal; iTri++)
    {
      for(int iNode = 0; iNode < NUM_NODES; iNode++)
      {
        double deltaX[3];
        vectorSubtract3D(node[iTri][iNode],old[iTri][iNode],deltaX);
        double distSq = deltaX[0]*deltaX[0] + deltaX[1]*deltaX[1] + deltaX[2]*deltaX[2];
        if(distSq > triggersq){
          
          flag = 1;
        }
      }
      if (flag) break;
    }

    // allreduce result
    MPI_Max_Scalar(flag,this->world);

    if(flag) return true;
    else     return false;
}

  /* ----------------------------------------------------------------------
   store node pos at last re-build
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMesh<NUM_NODES>::storeNodePosRebuild()
  {
    // just return for non-moving mesh
    if(!isMoving() && !isDeforming()) return;

    int nlocal = sizeLocal();
    double ***node = node_.begin();

    nodesLastRe_.empty();
    for(int i = 0; i < nlocal; i++)
        nodesLastRe_.add(node[i]);
  }

#endif
