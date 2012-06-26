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

#ifndef LMP_TRACKING_MESH_I_H
#define LMP_TRACKING_MESH_I_H

  /* ----------------------------------------------------------------------
   constructor
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  TrackingMesh<NUM_NODES>::TrackingMesh(LAMMPS *lmp)
  : MultiNodeMeshParallel<NUM_NODES>(lmp),
    customValues_(*(new CustomValueTracker(lmp,*this))),
    mapArray_(0),
    mapTagMax_(0),
    id_ (*this->prop().template addElementProperty< ScalarContainer<int> >("id","comm_none"/*ID does never change*/,"frame_invariant","restart_yes"))
  {
  }

  /* ----------------------------------------------------------------------
   destructor
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  TrackingMesh<NUM_NODES>::~TrackingMesh()
  {
     delete &customValues_;

     // deallocate map memory if exists
      if(mapArray_) clearMap();
  }

  /* ----------------------------------------------------------------------
   add / delete element
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void TrackingMesh<NUM_NODES>::addElement(double **nodeToAdd)
  {
    // this function is always called in serial mode
    
    MultiNodeMeshParallel<NUM_NODES>::addElement(nodeToAdd);

    // tracking mesh add memory
    
    customValues_.grow(this->sizeLocal());

    // set ID for element
    // ID starts from 0
    id_(this->sizeLocal()-1) = this->sizeLocal()-1;
  }

  template<int NUM_NODES>
  void TrackingMesh<NUM_NODES>::deleteElement(int n)
  {
    MultiNodeMeshParallel<NUM_NODES>::deleteElement(n);

    // tracking mesh delete code
    customValues_.deleteElement(n);
  }

  /* ----------------------------------------------------------------------
   recalculate properties on setup (on start and during simulation)
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void TrackingMesh<NUM_NODES>::refreshOwned(int setupFlag)
  {
    MultiNodeMeshParallel<NUM_NODES>::refreshOwned(setupFlag);
  }

  template<int NUM_NODES>
  void TrackingMesh<NUM_NODES>::refreshGhosts(int setupFlag)
  {
    MultiNodeMeshParallel<NUM_NODES>::refreshGhosts(setupFlag);
  }

  /* ----------------------------------------------------------------------
   clear and generate a global map for global-local lookup
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void TrackingMesh<NUM_NODES>::clearMap()
  {
      // deallocate old memory
      this->memory->destroy(mapArray_);
      mapArray_ = NULL;
  }

  template<int NUM_NODES>
  void TrackingMesh<NUM_NODES>::generateMap()
  {
      // deallocate old memory if exists
      if(mapArray_) clearMap();

      // get max ID of all proc
      int idmax = id_.max();
     MPI_Max_Scalar(idmax,mapTagMax_,this->world);

      // alocate and initialize new array
      // IDs start at 0, so have to use mapTagMax_+1
      this->memory->create(mapArray_,mapTagMax_+1,"TrackingMesh:mapArray_");
      for(int i = 0; i < mapTagMax_; i++)
        mapArray_[i] = -1;

      // build map for owned and ghost particles
      for (int i = this->sizeLocal()+this->sizeGhost()-1; i >= 0 ; i--)
      {
          
          mapArray_[id_(i)] = i;
      }
  }

  /* ----------------------------------------------------------------------
   clear ghost data that is communicated via forward comm - called in forw comm
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void TrackingMesh<NUM_NODES>::clearGhostForward(bool scale,bool translate,bool rotate)
  {
      MultiNodeMeshParallel<NUM_NODES>::clearGhostForward(scale,translate,rotate);

      // delete ghost data from container classes
      // delete only data that is communicated afterwards
      for(int i = this->sizeLocal()+this->sizeGhost()-1; i >= this->sizeLocal(); i--)
          customValues_.deleteForwardElement(i,scale,translate,rotate);

  }

  /* ----------------------------------------------------------------------
   push / pop functions for a list of
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int TrackingMesh<NUM_NODES>::elemListBufSize(int n,int operation,bool scale,bool translate,bool rotate)
  {
    int buf_size = 0;
    buf_size += MultiNodeMeshParallel<NUM_NODES>::elemListBufSize(n,operation,scale,translate,rotate);
    buf_size += customValues_.elemListBufSize(n,operation,scale,translate,rotate);
    return buf_size;
  }

  template<int NUM_NODES>
  int TrackingMesh<NUM_NODES>::pushElemListToBuffer(int n, int *list, double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    int nsend = 0;
    nsend += MultiNodeMeshParallel<NUM_NODES>::pushElemListToBuffer(n,list,&buf[nsend],operation,scale,translate,rotate);
    nsend += customValues_.pushElemListToBuffer(n,list,&buf[nsend],operation,scale,translate,rotate);
    return nsend;
  }

  template<int NUM_NODES>
  int TrackingMesh<NUM_NODES>::popElemListFromBuffer(int first, int n,double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    int nrecv = 0;
    nrecv += MultiNodeMeshParallel<NUM_NODES>::popElemListFromBuffer(first,n,&buf[nrecv],operation,scale,translate,rotate);
    nrecv += customValues_.popElemListFromBuffer(first,n,&buf[nrecv],operation,scale,translate,rotate);
    return nrecv;
  }

  /* ----------------------------------------------------------------------
   push / pop functions for a single element
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int TrackingMesh<NUM_NODES>::elemBufSize(int operation,bool scale,bool translate,bool rotate)
  {
    int buf_size = 0;
    buf_size += MultiNodeMeshParallel<NUM_NODES>::elemBufSize(operation,scale,translate,rotate);
    
    buf_size += customValues_.elemBufSize(operation,scale,translate,rotate);
    
    return buf_size;
  }

  template<int NUM_NODES>
  int TrackingMesh<NUM_NODES>::pushElemToBuffer(int n, double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    int nsend = 0;
    nsend += MultiNodeMeshParallel<NUM_NODES>::pushElemToBuffer(n,&buf[nsend],operation,scale,translate,rotate);
    nsend += customValues_.pushElemToBuffer(n,&buf[nsend],operation,scale,translate,rotate);
    return nsend;
  }

  template<int NUM_NODES>
  int TrackingMesh<NUM_NODES>::popElemFromBuffer(double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    int nrecv = 0;
    nrecv += MultiNodeMeshParallel<NUM_NODES>::popElemFromBuffer(&buf[nrecv],operation,scale,translate,rotate);
    nrecv += customValues_.popElemFromBuffer(&buf[nrecv],operation,scale,translate,rotate);
    return nrecv;
  }

  /* ----------------------------------------------------------------------
   push / pop functions for mesh properties
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int TrackingMesh<NUM_NODES>::meshPropsBufSize(int operation,bool scale,bool translate,bool rotate)
  {
    int buf_size = 0;
    buf_size += customValues_.meshPropsBufSize(operation,scale,translate,rotate);
    
    return buf_size;
  }

  template<int NUM_NODES>
  int TrackingMesh<NUM_NODES>::pushMeshPropsToBuffer(double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    int nsend = 0;
    nsend += customValues_.pushMeshPropsToBuffer(&buf[nsend],operation,scale,translate,rotate);
    return nsend;
  }

  template<int NUM_NODES>
  int TrackingMesh<NUM_NODES>::popMeshPropsFromBuffer(double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    int nrecv = 0;
    nrecv += customValues_.popMeshPropsFromBuffer(&buf[nrecv],operation,scale,translate,rotate);
    return nrecv;
  }

  /* ----------------------------------------------------------------------
   move / rotate  / scale
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void TrackingMesh<NUM_NODES>::move(double *vecTotal, double *vecIncremental)
  {
    
    MultiNodeMesh<NUM_NODES>::move(vecTotal, vecIncremental);
    customValues_.move(vecIncremental);
  }

  template<int NUM_NODES>
  void TrackingMesh<NUM_NODES>::move(double *vecIncremental)
  {
    
    MultiNodeMesh<NUM_NODES>::move(vecIncremental);
    customValues_.move(vecIncremental);
  }

  template<int NUM_NODES>
  void TrackingMesh<NUM_NODES>::moveElement(int i,double *vecIncremental)
  {
    MultiNodeMesh<NUM_NODES>::moveElement(i,vecIncremental);
    customValues_.moveElement(i,vecIncremental);
  }

  template<int NUM_NODES>
  void TrackingMesh<NUM_NODES>::rotate(double *totalQ, double *dQ,double *totalDispl, double *dDispl)
  {
    MultiNodeMesh<NUM_NODES>::rotate(totalQ,dQ,totalDispl,dDispl);

    customValues_.rotate(dQ);
    customValues_.move(dDispl);
  }

  template<int NUM_NODES>
  void TrackingMesh<NUM_NODES>::rotate(double *dQ,double *dDispl)
  {
    MultiNodeMesh<NUM_NODES>::rotate(dQ,dDispl);

    customValues_.rotate(dQ);
    customValues_.move(dDispl);
  }

  template<int NUM_NODES>
  void TrackingMesh<NUM_NODES>::scale(double factor)
  {
    
    MultiNodeMesh<NUM_NODES>::scale(factor);
    customValues_.scale(factor);
  }

#endif
