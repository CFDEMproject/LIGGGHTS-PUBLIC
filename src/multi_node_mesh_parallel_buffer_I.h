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

#ifndef LMP_MULTI_NODE_MESH_PARALLEL_BUFFER_I_H
#define LMP_MULTI_NODE_MESH_PARALLEL_BUFFER_I_H

  /* ----------------------------------------------------------------------
   push / pop for exchange
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::pushExchange(int dim,double *buf)
  {
      // scale translate rotate not needed here
      bool dummy = false, is_in_subdom;
      double checklo,checkhi;

      checklo = this->domain->sublo[dim];
      if(this->domain->subhi[dim] == this->domain->boxhi[dim])
        checkhi = this->domain->boxhi[dim] + SMALL_DMBRDR;
      else
        checkhi = this->domain->subhi[dim];

      int nsend = 0, nsend_this = 0;
      int i = 0;
      while(i < nLocal_)
      {
          if(!(this->center_(i)[dim] >= checklo && this->center_(i)[dim] < checkhi))
          {
              nsend_this = pushElemToBuffer(i,&(buf[nsend+1]),OPERATION_COMM_EXCHANGE,dummy,dummy,dummy);
              buf[nsend] = static_cast<double>(nsend_this+1);
              nsend += (nsend_this+1);
              
              if (nsend > maxsend_)
                  grow_send(nsend,1);
              this->deleteElement(i); 
          }
          else i++;
      }
      return nsend;
  }

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::popExchange(int nrecv,double *buf)
  {
      double center_elem[3];
      int m = 0, nrecv_this;

      // scale translate rotate not needed here
      bool dummy = false;

      while (m < nrecv)
      {
          // number of values is first in buffer
          nrecv_this = static_cast<int>(buf[m]);

          // center is next in buffer, test it
          vectorCopy3D(&(buf[m+1]),center_elem);

          if(this->domain->is_in_subdomain(center_elem))
          {
            popElemFromBuffer(&(buf[m+1]),OPERATION_COMM_EXCHANGE,dummy,dummy,dummy);
            nLocal_++;
            
          }
          
          m += nrecv_this;
      }
  }

  /* ----------------------------------------------------------------------
   return required buffer size for a list of elements for borders(),forwardComm()
   must match push / pop implementation
   depending on operation and if mesh scales, translates or rotates,
   different properties are communicated
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::listBufSize(int n,int operation,bool scale,bool translate,bool rotate)
  {
      return n*elemBufSize(operation,scale,translate,rotate);
  }

  /* ----------------------------------------------------------------------
   push a list of elements for borders(), forwardComm()
   depending on operation and if mesh scales, translates or rotates,
   different properties are communicated
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::pushListToBuffer(int n, int *list, double *buf,int operation,bool scale,bool translate,bool rotate)
  {
      
      int nsend = 0;

      if(operation == OPERATION_COMM_REVERSE)
      {
        this->error->one(FLERR,"TODO here");
        return nsend;
      }

      if(operation == OPERATION_RESTART)
      {
          this->error->one(FLERR,"TODO here");
          nsend += MultiNodeMesh<NUM_NODES>::node_.pushListToBuffer(n,list,&(buf[nsend]),operation);
          if(this->node_orig_)
            nsend += this->node_orig_->pushListToBuffer(n,list,&(buf[nsend]),operation);

          return nsend;
      }

      if(operation == OPERATION_COMM_EXCHANGE || operation == OPERATION_COMM_BORDERS)
      {
          
          nsend += MultiNodeMesh<NUM_NODES>::center_.pushListToBuffer(n,list,&(buf[nsend]),operation);
          nsend += MultiNodeMesh<NUM_NODES>::node_.pushListToBuffer(n,list,&(buf[nsend]),operation);
          nsend += MultiNodeMesh<NUM_NODES>::rBound_.pushListToBuffer(n,list,&(buf[nsend]),operation);
          if(this->node_orig_)
              nsend += this->node_orig_->pushListToBuffer(n,list,&(buf[nsend]),operation);
          return nsend;
      }

      if(operation == OPERATION_COMM_FORWARD)
      {
          
          // node_orig cannot change during a run
          //if(translate || rotate || scale)
          //  nsend += MultiNodeMesh<NUM_NODES>::node_.pushListToBuffer(n,list,&(buf[nsend]),operation);
          return nsend;
      }

      this->error->one(FLERR,"Illegal operation in MultiNodeMeshParallel<NUM_NODES>::pushElemToBuffer");
      return 0;
  }

  /* ----------------------------------------------------------------------
   pop a list of elements for borders(), forwardComm()
   depending on operation and if mesh scales, translates or rotates,
   different properties are communicated
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::popListFromBuffer(int first, int n, double *buf,int operation,bool scale,bool translate,bool rotate)
  {
      int nrecv = 0;

      if(operation == OPERATION_COMM_REVERSE)
      {
        this->error->one(FLERR,"TODO here");
        return nrecv;
      }

      if(operation == OPERATION_RESTART)
      {
          this->error->one(FLERR,"TODO here - add, recalc properties etc");
          nrecv += MultiNodeMesh<NUM_NODES>::node_.popListFromBuffer(first,n,&(buf[nrecv]),operation);
          if(MultiNodeMesh<NUM_NODES>::node_orig_)
            nrecv += MultiNodeMesh<NUM_NODES>::node_orig_->popListFromBuffer(first,n,&(buf[nrecv]),operation);

          return nrecv;
      }

      if(operation == OPERATION_COMM_EXCHANGE || operation == OPERATION_COMM_BORDERS)
      {
          nrecv += MultiNodeMesh<NUM_NODES>::center_.popListFromBuffer(first,n,&(buf[nrecv]),operation);
          nrecv += MultiNodeMesh<NUM_NODES>::node_.popListFromBuffer(first,n,&(buf[nrecv]),operation);
          nrecv += MultiNodeMesh<NUM_NODES>::rBound_.popListFromBuffer(first,n,&(buf[nrecv]),operation);
          if(MultiNodeMesh<NUM_NODES>::node_orig_)
            nrecv += MultiNodeMesh<NUM_NODES>::node_orig_->popListFromBuffer(first,n,&(buf[nrecv]),operation);

          return nrecv;
      }

      if(operation == OPERATION_COMM_FORWARD)
      {
          // node_orig cannot change during a run
          //if(translate || rotate || scale)
          //    nrecv += MultiNodeMesh<NUM_NODES>::node_.popListFromBuffer(first,n,&(buf[nrecv]),operation);

          return nrecv;
      }

      this->error->one(FLERR,"Illegal operation in MultiNodeMeshParallel<NUM_NODES>::popElemFromBuffer");
      return 0;
  }

  /* ----------------------------------------------------------------------
   return required buffer size for one element for exchange()
   must match push / pop implementation
   depending on operation and if mesh scales, translates or rotates,
   different properties are communicated
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::elemBufSize(int operation,bool scale,bool translate,bool rotate)
  {
      int size_buf = 0;

      if(operation == OPERATION_COMM_REVERSE)
      {
        this->error->one(FLERR,"TODO here");
        return size_buf;
      }

      if(operation == OPERATION_RESTART)
      {
          this->error->one(FLERR,"TODO here - add, recalc properties etc");
          this->error->one(FLERR,"TODO PUSH and POP random seed, maybe do this via fix mesh?");
          size_buf += MultiNodeMesh<NUM_NODES>::node_.elemBufSize();
          if(MultiNodeMesh<NUM_NODES>::node_orig_)
            size_buf += MultiNodeMesh<NUM_NODES>::node_orig_->elemBufSize();

          return size_buf;
      }

      if(operation == OPERATION_COMM_EXCHANGE || operation == OPERATION_COMM_BORDERS)
      {
          size_buf += MultiNodeMesh<NUM_NODES>::center_.elemBufSize();
          size_buf += MultiNodeMesh<NUM_NODES>::node_.elemBufSize();
          size_buf += MultiNodeMesh<NUM_NODES>::rBound_.elemBufSize();
          if(MultiNodeMesh<NUM_NODES>::node_orig_)
            size_buf += MultiNodeMesh<NUM_NODES>::node_orig_->elemBufSize();

          return size_buf;
      }

      if(operation == OPERATION_COMM_FORWARD)
      {
          // node_orig cannot change during a run
          //if(translate || rotate || scale)
          //  size_buf += MultiNodeMesh<NUM_NODES>::node_.elemBufSize();

          return size_buf;
      }

      this->error->one(FLERR,"Illegal operation in MultiNodeMeshParallel<NUM_NODES>::elemBufSize");
      return 0;
  }

  /* ----------------------------------------------------------------------
   push one element for exchange()
   depending on operation and if mesh scales, translates or rotates,
   different properties are communicated
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::pushElemToBuffer(int i, double *buf,int operation,bool scale,bool translate,bool rotate)
  {
      int nsend = 0;

      if(operation == OPERATION_COMM_REVERSE)
      {
        this->error->one(FLERR,"TODO here");
        return nsend;
      }

      if(operation == OPERATION_RESTART)
      {
          this->error->one(FLERR,"TODO here");
          nsend += MultiNodeMesh<NUM_NODES>::node_.pushElemToBuffer(i,&(buf[nsend]),operation);
          if(this->node_orig_)
            nsend += this->node_orig_->pushElemToBuffer(i,&(buf[nsend]),operation);

          return nsend;
      }

      if(operation == OPERATION_COMM_EXCHANGE || operation == OPERATION_COMM_BORDERS)
      {
          
          nsend += MultiNodeMesh<NUM_NODES>::center_.pushElemToBuffer(i,&(buf[nsend]),operation);
          nsend += MultiNodeMesh<NUM_NODES>::node_.pushElemToBuffer(i,&(buf[nsend]),operation);
          nsend += MultiNodeMesh<NUM_NODES>::rBound_.pushElemToBuffer(i,&(buf[nsend]),operation);
          if(this->node_orig_)
              nsend += this->node_orig_->pushElemToBuffer(i,&(buf[nsend]),operation);
          return nsend;
      }

      if(operation == OPERATION_COMM_FORWARD)
      {
          
          // node_orig cannot change during a run
          //if(translate || rotate || scale)
          //    nsend += MultiNodeMesh<NUM_NODES>::node_.pushElemToBuffer(i,&(buf[nsend]),operation);

          return nsend;
      }

      this->error->one(FLERR,"Illegal operation in MultiNodeMeshParallel<NUM_NODES>::pushElemToBuffer");
      return 0;
  }

  /* ----------------------------------------------------------------------
   pop one element for exchange()
   depending on operation and if mesh scales, translates or rotates,
   different properties are communicated
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::popElemFromBuffer(double *buf,int operation,bool scale,bool translate,bool rotate)
  {
      int nrecv = 0;

      if(operation == OPERATION_COMM_REVERSE)
      {
        this->error->one(FLERR,"TODO here");
        return nrecv;
      }

      if(operation == OPERATION_RESTART)
      {
          this->error->one(FLERR,"TODO here - add, recalc properties etc");
          nrecv += MultiNodeMesh<NUM_NODES>::node_.popElemFromBuffer(&(buf[nrecv]),operation);
          if(MultiNodeMesh<NUM_NODES>::node_orig_)
            nrecv += MultiNodeMesh<NUM_NODES>::node_orig_->popElemFromBuffer(&(buf[nrecv]),operation);

          return nrecv;
      }

      if(operation == OPERATION_COMM_EXCHANGE || operation == OPERATION_COMM_BORDERS)
      {
          nrecv += MultiNodeMesh<NUM_NODES>::center_.popElemFromBuffer(&(buf[nrecv]),operation);
          nrecv += MultiNodeMesh<NUM_NODES>::node_.popElemFromBuffer(&(buf[nrecv]),operation);
          nrecv += MultiNodeMesh<NUM_NODES>::rBound_.popElemFromBuffer(&(buf[nrecv]),operation);
          if(MultiNodeMesh<NUM_NODES>::node_orig_)
            nrecv += MultiNodeMesh<NUM_NODES>::node_orig_->popElemFromBuffer(&(buf[nrecv]),operation);

          return nrecv;
      }

      if(operation == OPERATION_COMM_FORWARD)
      {
          // node_orig cannot change during a run
          //if(translate || rotate || scale)
          //    nrecv += MultiNodeMesh<NUM_NODES>::node_.popElemFromBuffer(&(buf[nrecv]),operation);

          return nrecv;
      }

      this->error->one(FLERR,"Illegal operation in MultiNodeMeshParallel<NUM_NODES>::popElemFromBuffer");
      return 0;
  }

#endif
