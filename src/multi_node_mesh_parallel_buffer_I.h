/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
   the producer of the LIGGGHTS® software and the CFDEM®coupling software
   See http://www.cfdem.com/terms-trademark-policy for details.

   LIGGGHTS® is based on LAMMPS
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
  int MultiNodeMeshParallel<NUM_NODES>::pushExchange(int dim)
  {
      
      // scale translate rotate not needed here
      bool dummy = false;
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
              nsend_this = pushElemToBuffer(i,&(buf_send_[nsend+1]),OPERATION_COMM_EXCHANGE,dummy,dummy,dummy);
              buf_send_[nsend] = static_cast<double>(nsend_this+1);
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
  void MultiNodeMeshParallel<NUM_NODES>::popExchange(int nrecv,int dim,double *buf)
  {
      double center_elem[3];
      double checklo,checkhi;
      int m = 0, nrecv_this;

      // scale translate rotate not needed here
      bool dummy = false;

      checklo = this->domain->sublo[dim];
      if(this->domain->subhi[dim] == this->domain->boxhi[dim])
        checkhi = this->domain->boxhi[dim] + SMALL_DMBRDR;
      else
        checkhi = this->domain->subhi[dim];

      while (m < nrecv)
      {
          // number of values is first in buffer
          nrecv_this = static_cast<int>(buf[m]);

          // center is next in buffer, test it
          vectorCopy3D(&(buf[m+1]),center_elem);

          if(center_elem[dim] >= checklo && center_elem[dim] < checkhi)
          {
            popElemFromBuffer(&(buf[m+1]),OPERATION_COMM_EXCHANGE,dummy,dummy,dummy);
            nLocal_++;
            
          }
          
          m += nrecv_this;
      }
  }

  /* ----------------------------------------------------------------------
   restart functionality - write all required data into restart buffer
   executed on all processes, but only proc 0 writes into writebuf
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::writeRestart(FILE *fp)
  {
      int size_this;

      // # elements
      int nlocal = this->sizeLocal();
      int nglobal = sizeGlobal();

      // buffer sizes
      int sizeMesh, sizeElements, sizeElements_all;

      sizeMesh = sizeRestartMesh();
      sizeElements = nlocal * (sizeRestartElement() + 1); 

      double *bufMesh = NULL, *sendbufElems = NULL, *recvbufElems = NULL;
      bool dummy = false;

      // pack global data into buffer
      // do this only on proc 0
      if(this->comm->me == 0)
      {
          this->memory->create(bufMesh,sizeMesh,"MultiNodeMeshParallel::writeRestart:bufMesh");
          pushMeshPropsToBuffer(bufMesh, OPERATION_RESTART,dummy,dummy,dummy);
      }

      // allocate send buffer and pack element data
      // all local elements are in list
      this->memory->create(sendbufElems,sizeElements,"MultiNodeMeshParallel::writeRestart:sendbuf");
      sizeElements = 0;
      for(int i = 0; i < nlocal; i++)
      {
          size_this = pushElemToBuffer(i,&(sendbufElems[sizeElements+1]),OPERATION_RESTART,dummy,dummy,dummy);
          sendbufElems[sizeElements] = static_cast<double>(size_this+1);
          sizeElements += (size_this+1);
      }

      // gather the per-element data
      
      sizeElements_all = MPI_Gather0_Vector(sendbufElems,sizeElements,recvbufElems,this->world);

      // actually write data to restart file
      // do this only on proc 0
      if(this->comm->me == 0)
      {
        double nG = static_cast<double>(nglobal);

        // for error check
        double sE = static_cast<double>(sizeRestartElement());
        double sM = static_cast<double>(sizeRestartMesh());

        // size with 3 extra values
        int size = (sizeMesh+sizeElements_all+3) * sizeof(double);

        // write size
        fwrite(&size,sizeof(int),1,fp);

        // write 3 extra values
        fwrite(&nG,sizeof(double),1,fp);
        fwrite(&sE,sizeof(double),1,fp);
        fwrite(&sM,sizeof(double),1,fp);

        // write per-element and mesh data
        fwrite(recvbufElems,sizeof(double),sizeElements_all,fp);
        fwrite(bufMesh,sizeof(double),sizeMesh,fp);
      }

      // free mem

      if(bufMesh)
        this->memory->destroy(bufMesh);

      this->memory->destroy(sendbufElems);

      if(recvbufElems)
        delete []recvbufElems;
  }

  /* ----------------------------------------------------------------------
   restart functionality - read all required data from restart buffer
   executed on all processes
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::restart(double *list)
  {
      int m, nglobal, nrecv_this, sE, sM;
      bool dummy = false;

      m = 0;

      nglobal = static_cast<int> (list[m++]);
      sE = static_cast<int> (list[m++]);
      sM = static_cast<int> (list[m++]);

      if(sE != sizeRestartElement() || sM != sizeRestartMesh())
          this->error->all(FLERR,"Incompatible mesh restart file - mesh has different properties in restarted simulation");

      for(int i = 0; i < nglobal; i++)
      {
          nrecv_this = static_cast<int>(list[m]);
          
          popElemFromBuffer(&(list[m+1]),OPERATION_RESTART,dummy,dummy,dummy);
          m += nrecv_this;
      }

      this->prop().deleteRestartGlobal(dummy,dummy,dummy);
      popMeshPropsFromBuffer(&list[m],OPERATION_RESTART,dummy,dummy,dummy);
  }

  /* ----------------------------------------------------------------------
   size of restart data
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::sizeRestartMesh()
  {
      
      bool dummy = false;
      return meshPropsBufSize(OPERATION_RESTART,dummy,dummy,dummy);
  }

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::sizeRestartElement()
  {
      
      bool dummy = false;
      return elemBufSize(OPERATION_RESTART,dummy,dummy,dummy);
  }

  /* ----------------------------------------------------------------------
   return required buffer size for a list of elements for borders(),forwardComm()
   must match push / pop implementation
   depending on operation and if mesh scales, translates or rotates,
   different properties are communicated
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::elemListBufSize(int n,int operation,bool scale,bool translate,bool rotate)
  {
      return n*elemBufSize(operation,scale,translate,rotate);
  }

  /* ----------------------------------------------------------------------
   push a list of elements for borders(), forwardComm()
   depending on operation and if mesh scales, translates or rotates,
   different properties are communicated
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::pushElemListToBuffer(int n, int *list, double *buf, int operation, bool , bool, bool)
  {
      
      int nsend = 0;

      if(OPERATION_COMM_EXCHANGE == operation || OPERATION_COMM_BORDERS == operation)
      {
          
          nsend += MultiNodeMesh<NUM_NODES>::center_.pushElemListToBuffer(n,list,&(buf[nsend]),operation);
          nsend += MultiNodeMesh<NUM_NODES>::node_.pushElemListToBuffer(n,list,&(buf[nsend]),operation);
          nsend += MultiNodeMesh<NUM_NODES>::rBound_.pushElemListToBuffer(n,list,&(buf[nsend]),operation);
          if(this->node_orig_)
              nsend += this->node_orig_->pushElemListToBuffer(n,list,&(buf[nsend]),operation);
          return nsend;
      }

      if(OPERATION_COMM_FORWARD == operation)
      {
          
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
  int MultiNodeMeshParallel<NUM_NODES>::popElemListFromBuffer(int first, int n, double *buf, int operation, bool, bool, bool)
  {
      int nrecv = 0;

      if(OPERATION_COMM_EXCHANGE == operation || OPERATION_COMM_BORDERS == operation)
      {
          nrecv += MultiNodeMesh<NUM_NODES>::center_.popElemListFromBuffer(first,n,&(buf[nrecv]),operation);
          nrecv += MultiNodeMesh<NUM_NODES>::node_.popElemListFromBuffer(first,n,&(buf[nrecv]),operation);
          nrecv += MultiNodeMesh<NUM_NODES>::rBound_.popElemListFromBuffer(first,n,&(buf[nrecv]),operation);
          if(MultiNodeMesh<NUM_NODES>::node_orig_)
            nrecv += MultiNodeMesh<NUM_NODES>::node_orig_->popElemListFromBuffer(first,n,&(buf[nrecv]),operation);
          return nrecv;
      }

      if(OPERATION_COMM_FORWARD == operation)
      {
          
          //    nrecv += MultiNodeMesh<NUM_NODES>::node_.popListFromBuffer(first,n,&(buf[nrecv]),operation);
          return nrecv;
      }

      this->error->one(FLERR,"Illegal operation in MultiNodeMeshParallel<NUM_NODES>::popElemFromBuffer");
      return 0;
  }

  /* ----------------------------------------------------------------------
   push a list of elements for reverseComm()
   depending on operation and if mesh scales, translates or rotates,
   different properties are communicated
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::pushElemListToBufferReverse(int, int, double*, int operation, bool, bool, bool)
  {
      int nsend = 0;

      if(OPERATION_COMM_REVERSE == operation)
      {
        
        return nsend;
      }

      this->error->one(FLERR,"Illegal operation in MultiNodeMeshParallel<NUM_NODES>::popElemFromBuffer");
      return 0;
  }

  /* ----------------------------------------------------------------------
   pop a list of elements for reverseComm()
   depending on operation and if mesh scales, translates or rotates,
   different properties are communicated
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::popElemListFromBufferReverse(int, int*, double*, int operation, bool, bool, bool)
  {
      int nrecv = 0;

      if(OPERATION_COMM_REVERSE == operation)
      {
        
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
  int MultiNodeMeshParallel<NUM_NODES>::elemBufSize(int operation, bool, bool, bool)
  {
      int size_buf = 0;

      if(OPERATION_RESTART == operation)
      {
          size_buf += MultiNodeMesh<NUM_NODES>::node_.elemBufSize();
          return size_buf;
      }

      if(OPERATION_COMM_EXCHANGE == operation || OPERATION_COMM_BORDERS == operation)
      {
          size_buf += MultiNodeMesh<NUM_NODES>::center_.elemBufSize();
          size_buf += MultiNodeMesh<NUM_NODES>::node_.elemBufSize();
          size_buf += MultiNodeMesh<NUM_NODES>::rBound_.elemBufSize();
          if(MultiNodeMesh<NUM_NODES>::node_orig_)
            size_buf += MultiNodeMesh<NUM_NODES>::node_orig_->elemBufSize();
          return size_buf;
      }

      if(OPERATION_COMM_FORWARD == operation)
      {
          
          return size_buf;
      }

      if(OPERATION_COMM_REVERSE == operation)
      {
        
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
  int MultiNodeMeshParallel<NUM_NODES>::pushElemToBuffer(int i, double *buf, int operation, bool, bool, bool)
  {
      int nsend = 0;

      if(OPERATION_RESTART == operation)
      {
          nsend += MultiNodeMesh<NUM_NODES>::node_.pushElemToBuffer(i,&(buf[nsend]),operation);

          return nsend;
      }

      if(OPERATION_COMM_EXCHANGE == operation || OPERATION_COMM_BORDERS == operation)
      {
          
          nsend += MultiNodeMesh<NUM_NODES>::center_.pushElemToBuffer(i,&(buf[nsend]),operation);
          nsend += MultiNodeMesh<NUM_NODES>::node_.pushElemToBuffer(i,&(buf[nsend]),operation);
          nsend += MultiNodeMesh<NUM_NODES>::rBound_.pushElemToBuffer(i,&(buf[nsend]),operation);
          if(this->node_orig_)
              nsend += this->node_orig_->pushElemToBuffer(i,&(buf[nsend]),operation);
          return nsend;
      }

      this->error->one(FLERR,"Illegal operation in MultiNodeMeshParallel<NUM_NODES>::pushElemToBuffer");
      return 0;
  }

  /* ----------------------------------------------------------------------
   pop one element for exchange, restart or restart
   depending on operation and if mesh scales, translates or rotates,
   different properties are communicated
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::popElemFromBuffer(double *buf, int operation, bool, bool, bool)
  {
      int nrecv = 0;

      if(OPERATION_RESTART == operation)
      {
          MultiVectorContainer<double,NUM_NODES,3> nodeTmp("nodeTmp");
          
          nrecv += nodeTmp.popElemFromBuffer(&(buf[nrecv]),operation);
          this->addElement(nodeTmp.begin()[0],-1);

          this->prop().deleteRestartElement(nLocal_-1,false,false,false);

          return nrecv;
      }

      if(OPERATION_COMM_EXCHANGE == operation || OPERATION_COMM_BORDERS == operation)
      {
          nrecv += MultiNodeMesh<NUM_NODES>::center_.popElemFromBuffer(&(buf[nrecv]),operation);
          nrecv += MultiNodeMesh<NUM_NODES>::node_.popElemFromBuffer(&(buf[nrecv]),operation);
          nrecv += MultiNodeMesh<NUM_NODES>::rBound_.popElemFromBuffer(&(buf[nrecv]),operation);
          if(MultiNodeMesh<NUM_NODES>::node_orig_)
            nrecv += MultiNodeMesh<NUM_NODES>::node_orig_->popElemFromBuffer(&(buf[nrecv]),operation);
          return nrecv;
      }

      this->error->one(FLERR,"Illegal operation in MultiNodeMeshParallel<NUM_NODES>::popElemFromBuffer");
      return 0;
  }

#endif
