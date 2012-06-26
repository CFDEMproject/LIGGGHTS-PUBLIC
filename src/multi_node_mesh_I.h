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
  : AbstractMesh(lmp), node_(), node_orig_(0), nMove_(0), stepLastReset_(-1),
    nScale_(0), nTranslate_(0), nRotate_(0),
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
   set ID
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMesh<NUM_NODES>::setMeshID(const char *_id)
  {
      if(mesh_id_) delete []mesh_id_;
      mesh_id_ = new char[strlen(_id)+1];
      strcpy(mesh_id_,_id);
  }

  /* ----------------------------------------------------------------------
   add an element - only called at mesh construction
   i.e. only used to construct local elements
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMesh<NUM_NODES>::addElement(double **nodeToAdd)
  {
    double avg[3];

    // add node
    node_.add(nodeToAdd);

    // calculate center
    vectorZeroize3D(avg);
    for(int i = 0; i < NUM_NODES; i++)
        vectorAdd3D(nodeToAdd[i],avg,avg);
    vectorScalarDiv3D(avg,static_cast<double>(NUM_NODES));
    center_.add(avg);

    int n = sizeLocal();

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
      if(node_orig_ && setupFlag)
        storeNodePos(0,sizeLocal());

      // nothing more to do here, necessary initialitation done in addElement()
  }

  template<int NUM_NODES>
  void MultiNodeMesh<NUM_NODES>::refreshGhosts(int setupFlag)
  {
      if(node_orig_ && setupFlag)
        storeNodePos(sizeLocal(),sizeLocal()+sizeGhost());

      // nothing more to do here, necessary initialitation done in addElement()
  }

  /* ----------------------------------------------------------------------
   comparison
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  bool MultiNodeMesh<NUM_NODES>::nodesAreEqual(int iGrp, int iNode, int jGrp, int jNode)
  {
    for(int i=0;i<3;i++)
      if(!MathExtraLiggghts::compDouble(node_(iGrp)[iNode][i],node_(jGrp)[jNode][i],1e-8))
        return false;
    return true;
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
          this->memory->create<double>(tmp,NUM_NODES,3,"MultiNodeMesh:tmp");

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
   store current node position for use by moving mesh
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMesh<NUM_NODES>::storeNodePos(int ilo, int ihi)
  {
    if(!node_orig_)
        error->one(FLERR,"Internal error in MultiNodeMesh<NUM_NODES>::storeNodePos");

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
  void MultiNodeMesh<NUM_NODES>::resetNodesToOrig()
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
    }
  }

  /* ----------------------------------------------------------------------
   move mesh by amount vecTotal, starting from original position
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMesh<NUM_NODES>::move(double *vecTotal, double *vecIncremental)
  {
    
    resetNodesToOrig();

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
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMesh<NUM_NODES>::rotate(double totalAngle, double dAngle, double *axis, double *p)
  {
    double totalQ[4],dQ[4], axisNorm[3];
    double p_rot[3], totalDispl[3], dDispl[3];

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

    // calc displacement caused by total rotation around axis through center
    MathExtraLiggghts::vec_quat_rotate(p, totalQ,p_rot);
    vectorSubtract3D(p,p_rot,totalDispl);

    // calc displacement caused by incremental rotation around axis through center
    MathExtraLiggghts::vec_quat_rotate(p, dQ,p_rot);
    vectorSubtract3D(p,p_rot,dDispl);

    // apply rotation around center axis + displacement
    // = rotation around axis through p
    rotate(totalQ,dQ,totalDispl,dDispl);
  }

  template<int NUM_NODES>
  void MultiNodeMesh<NUM_NODES>::rotate(double *totalQ, double *dQ,double *totalDispl, double *dDispl)
  {
    
    resetNodesToOrig();

    int n = sizeLocal() + sizeGhost();

    // perform total rotation for data in this class
    for(int i = 0; i < n; i++)
    {
      vectorZeroize3D(center_(i));

      for(int j = 0; j < NUM_NODES; j++)
      {
        MathExtraLiggghts::vec_quat_rotate(node_(i)[j], totalQ, node_(i)[j]);
        vectorAdd3D(node_(i)[j],totalDispl,node_(i)[j]);
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
    double dQ[4], axisNorm[3];
    double p_rot[3], dDispl[3];

    // rotates around axis through p

    // normalize axis
    vectorCopy3D(axis,axisNorm);
    vectorScalarDiv3D(axisNorm,vectorMag3D(axisNorm));

    // quat for rotation since last time-step
    dQ[0] = cos(dAngle*0.5);
    for(int i = 0; i < 3; i++)
      dQ[i+1] = axis[i]*sin(dAngle*0.5);

    // calc displacement caused by incremental rotation around axis through center
    MathExtraLiggghts::vec_quat_rotate(p, dQ,p_rot);
    vectorSubtract3D(p,p_rot,dDispl);

    // apply rotation around center axis + displacement
    // = rotation around axis through p
    rotate(dQ,dDispl);
  }

  template<int NUM_NODES>
  void MultiNodeMesh<NUM_NODES>::rotate(double *dQ, double *dDispl)
  {
    
    int n = sizeLocal() + sizeGhost();

    // perform total rotation for data in this class
    for(int i = 0; i < n; i++)
    {
      vectorZeroize3D(center_(i));
      for(int j = 0; j < NUM_NODES; j++)
      {
        MathExtraLiggghts::vec_quat_rotate(node_(i)[j], dQ,node_(i)[j]);
        vectorAdd3D(node_(i)[j],dDispl,node_(i)[j]);
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

#endif
