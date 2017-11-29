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
    (if not contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef LMP_MULTI_NODE_MESH_PARALLEL_I_H
#define LMP_MULTI_NODE_MESH_PARALLEL_I_H

#define BIG_MNMP 1.0e20
#define BUFFACTOR_MNMP 1.5
#define BUFMIN_MNMP 2000
#define BUFEXTRA_MNMP 2000

  /* ----------------------------------------------------------------------
   consturctors
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  MultiNodeMeshParallel<NUM_NODES>::MultiNodeMeshParallel(LAMMPS *lmp)
  : MultiNodeMesh<NUM_NODES>(lmp),
    doParallellization_(true),
    nLocal_(0), nGhost_(0), nGlobal_(0), nGlobalOrig_(0),
    isParallel_(false),
    isInsertionMesh_(false),
    maxsend_(0), maxrecv_(0),
    buf_send_(0), buf_recv_(0),
    half_atom_cut_(0.),
    size_exchange_(0),
    size_forward_(0),
    size_border_(0),
    maxforward_(0),maxreverse_(0),
    nswap_(0),
    maxswap_(0),
    sendnum_(0),recvnum_(0),
    firstrecv_(0),
    sendproc_(0),recvproc_(0),
    size_forward_recv_(0),
    size_reverse_recv_(0),
    slablo_(0),slabhi_(0),
    sendlist_(0),
    sendwraplist_(0),
    maxsendlist_(0),
    pbc_flag_(0),
    pbc_(0)
  {
      // initialize comm buffers & exchange memory
      
      maxsend_ = BUFMIN_MNMP;
      this->memory->create(buf_send_,maxsend_+BUFEXTRA_MNMP,"MultiNodeMeshParallel:buf_send");
      maxrecv_ = BUFMIN_MNMP;
      this->memory->create(buf_recv_,maxrecv_,"MultiNodeMeshParallel:buf_recv");

      maxswap_ = 6;
      allocate_swap(maxswap_);

      sendlist_ = (int **) this->memory->smalloc(maxswap_*sizeof(int *),"MultiNodeMeshParallel:sendlist");
      sendwraplist_ = (int **) this->memory->smalloc(maxswap_*sizeof(int *),"MultiNodeMeshParallel:sendlist");
      this->memory->create(maxsendlist_,maxswap_,"MultiNodeMeshParallel:maxsendlist");
      for (int i = 0; i < maxswap_; i++) {
        maxsendlist_[i] = BUFMIN_MNMP;
        this->memory->create(sendlist_[i],BUFMIN_MNMP,"MultiNodeMeshParallel:sendlist[i]");
        this->memory->create(sendwraplist_[i],BUFMIN_MNMP,"MultiNodeMeshParallel:sendlist[i]");
      }
  }

  /* ----------------------------------------------------------------------
   destructor
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  MultiNodeMeshParallel<NUM_NODES>::~MultiNodeMeshParallel()
  {
      free_swap();

      if (sendlist_)
        for (int i = 0; i < maxswap_; i++)
            this->memory->destroy(sendlist_[i]);
      if (sendwraplist_)
        for (int i = 0; i < maxswap_; i++)
            this->memory->destroy(sendwraplist_[i]);

      this->memory->sfree(sendlist_);
      this->memory->sfree(sendwraplist_);
      this->memory->destroy(maxsendlist_);

      this->memory->destroy(buf_send_);
      this->memory->destroy(buf_recv_);
  }

  /* ----------------------------------------------------------------------
   add and delete elements
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  bool MultiNodeMeshParallel<NUM_NODES>::addElement(double **nodeToAdd)
  {
    
    if(MultiNodeMesh<NUM_NODES>::addElement(nodeToAdd))
    {
        nLocal_++;
        return true;
    }
    return false;
  }

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::deleteElement(int n)
  {
    
    if(n < nLocal_ && nGhost_ != 0)
        this->error->one(FLERR,"Illegal call to MultiNodeMeshParallel<NUM_NODES>::deleteElement");

    MultiNodeMesh<NUM_NODES>::deleteElement(n);

    if(n >= nLocal_)
        nGhost_--;
    else
        nLocal_--;
  }

  /* ----------------------------------------------------------------------
   recalculate properties on setup (on start and during simulation)
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::refreshOwned(int setupFlag)
  {
    MultiNodeMesh<NUM_NODES>::refreshOwned(setupFlag);
  }

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::refreshGhosts(int setupFlag)
  {
    MultiNodeMesh<NUM_NODES>::refreshGhosts(setupFlag);
  }

  /* ----------------------------------------------------------------------
   completely clear ghosts - called in borders()
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::clearGhosts()
  {
      // delete ghost data from container classes

      while(nGhost_ > 0)
      {
          
          deleteElement(nLocal_);
      }
  }

  /* ----------------------------------------------------------------------
   clear ghost data that is communicated via forward comm - called in forw comm
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::clearGhostForward(bool scale,bool translate,bool rotate)
  {
      // delete ghost data from container classes
      // delete only data that is communicated afterwards

      for(int i = this->sizeLocal()+this->sizeGhost()-1; i >= this->sizeLocal(); i--)
      {
          // clear ghost data that belongs to this class
          // must match push/pop implementation for forward comm in this class
          if(translate || rotate || scale)
          {
            this->node_.del(i);
            this->center_.del(i);
          }
          if(scale)
            this->rBound_.del(i);
      }
  }

  /* ----------------------------------------------------------------------
   check if all elements are in domain
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  bool MultiNodeMeshParallel<NUM_NODES>::allNodesInsideSimulationBox()
  {
    int flag = 0;
    for(int i=0;i<sizeLocal();i++)
      for(int j=0;j<NUM_NODES;j++)
      {
        
        if(!this->domain->is_in_domain(this->node_(i)[j]))
        {
            flag = 1;
            break;
        }
      }

    MPI_Max_Scalar(flag,this->world);
    if(flag) return false;
    else return true;
  }

  /* ----------------------------------------------------------------------
   set flag if used as insertion mesh
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::useAsInsertionMesh(bool parallelflag)
  {
    
    isInsertionMesh_ = true;
    
    if(!parallelflag)
    {
        if(isParallel())
            this->error->all(FLERR,"If a run command is between the fix mesh/surface and the "
                             "fix insert command, you have to use fix mesh/surface/planar for "
                             "the insertion mesh");
        doParallellization_ = false;
    }
  }

  /* ----------------------------------------------------------------------
   setup of communication
  ------------------------------------------------------------------------- */

   template<int NUM_NODES>
   void MultiNodeMeshParallel<NUM_NODES>::setup()
   {
       if(!doParallellization_) return;

       double sublo[3],subhi[3], extent_acc;
       double rBound_max, cut_ghost;
       double **sublo_all, **subhi_all;

       int nprocs = this->comm->nprocs;
       int myloc[3], loc_dim, nextproc, need_this;

       // get required size of communication per element
       bool scale = this->isScaling();
       bool translate = this->isTranslating();
       bool rotate = this->isRotating();

       size_exchange_ = elemBufSize(OPERATION_COMM_EXCHANGE, NULL, scale,translate,rotate) + 1;
       size_border_ = elemBufSize(OPERATION_COMM_BORDERS, NULL, scale,translate,rotate);
       size_forward_ = elemBufSize(OPERATION_COMM_FORWARD, NULL, scale,translate,rotate);
       size_reverse_ = elemBufSize(OPERATION_COMM_REVERSE, NULL, scale,translate,rotate);

       // maxforward = # of datums in largest forward communication
       // maxreverse = # of datums in largest reverse communication

       maxforward_ = MathExtraLiggghts::max(size_exchange_,size_border_,size_forward_);
       maxreverse_ = size_reverse_;

       // copy comm and domain data
       vectorCopy3D(this->comm->myloc,myloc);
       vectorCopy3D(this->domain->sublo,sublo);
       vectorCopy3D(this->domain->subhi,subhi);

       this->memory->create(sublo_all,nprocs,3,"MultiNodeMeshParallel::setup() sublo_all");
       this->memory->create(subhi_all,nprocs,3,"MultiNodeMeshParallel::setup() subhi_all");

       // ghost elements are for computing interaction with owned particles
       // so need to aquire ghost elements that overlap my subbox extened by
       // half neigh cutoff
       
       half_atom_cut_ = this->neighbor->cutneighmax / 2.;

       if(this->isMoving())
         half_atom_cut_+= this->neighbor->skin / 2.;

       // calculate maximum bounding radius of elements across all procs
       rBound_max = 0.;
       for(int i = 0; i < sizeLocal(); i++)
           rBound_max = std::max(this->rBound_(i),rBound_max);
       MPI_Max_Scalar(rBound_max,this->world);

       // mesh element ghost cutoff is element bounding radius plus half atom neigh cut
       cut_ghost = rBound_max + half_atom_cut_;

       // set up maxneed_, sendneed_
       // account for non-uniform boundaries due to load-balancing
       // so aquire sub-box bounds from all processors

       MPI_Allgather(sublo,3,MPI_DOUBLE,&(sublo_all[0][0]),3,MPI_DOUBLE,this->world);
       MPI_Allgather(subhi,3,MPI_DOUBLE,&(subhi_all[0][0]),3,MPI_DOUBLE,this->world);

       // set up maxneed_ and sendneed_
       // assume element with max bound radius is in my subbox
       
       for(int dim = 0; dim < 3; dim++)
       {
           bool is_x = dim == 0 ? true : false;
           bool is_y = dim == 1 ? true : false;
           bool is_z = dim == 2 ? true : false;

           // go each direction (N-S-E-W-UP-DN)
           
           maxneed_[dim] = 0;
           for(int way = -1; way <= 1; way += 2)
           {
               // start from location of myself
               // reset accumulated extent
               loc_dim = myloc[dim];
               extent_acc = 0.;
               need_this = 0;
               sendneed_[dim][way == -1 ? 0 : 1] = 0;

               while(extent_acc < cut_ghost)
               {
                   // increase or decrease location
                   loc_dim += way;

                   // break if at dead end and non-pbc
                   if( (loc_dim < 0 && !this->domain->periodicity[dim]) ||
                       (loc_dim > this->comm->procgrid[dim]-1 && !this->domain->periodicity[dim]) )
                           break;

                   // wrap around PBCs
                   if(loc_dim < 0 && this->domain->periodicity[dim])
                      loc_dim = this->comm->procgrid[dim]-1;

                   if(loc_dim > this->comm->procgrid[dim]-1)
                      loc_dim = 0;

                   // increase counters
                   need_this++;
                   sendneed_[dim][way == -1 ? 0 : 1]++;

                   // go to next proc in proc grid and add its extent
                   nextproc = this->comm->grid2proc[is_x ? loc_dim : myloc[0]]
                                                   [is_y ? loc_dim : myloc[1]]
                                                   [is_z ? loc_dim : myloc[2]];
                   extent_acc += subhi_all[nextproc][dim] - sublo_all[nextproc][dim];
               }

               maxneed_[dim] = std::max(maxneed_[dim],need_this);
           }

           // limit maxneed for non-pbc
           
           if(maxneed_[dim] > this->comm->procgrid[dim]-1 && !this->domain->periodicity[dim])
               maxneed_[dim] = this->comm->procgrid[dim]-1;
       }

       // maxneed_ summed accross all processors
       MPI_Max_Vector(maxneed_,3,this->world);

       destroy(sublo_all);
       destroy(subhi_all);

       // allocate comm memory
       
       nswap_ = 2 * (maxneed_[0]+maxneed_[1]+maxneed_[2]);
       if (nswap_ > maxswap_) grow_swap(nswap_);

       // setup parameters for each exchange:
       //   slablo_/slabhi_ = boundaries for slab of elements to send at each swap
       //   use -BIG/midpt/BIG to insure all elements included even if round-off occurs
       //   if round-off, atoms elements across PBC can be < or > than subbox boundary
       //   note that borders() only loops over subset of elements during each swap

       // treat all as PBC here, non-PBC is handled in borders() via r/s need[][]
       // pbc_flag_: 0 = nothing across a boundary, 1 = something across a boundary
       // pbc_ = -1/0/1 for PBC factor in each of 3/6 orthogonal/triclinic dirs
       // 1st part of if statement is sending to the west/south/down
       // 2nd part of if statement is sending to the east/north/up

       int dim,ineed;

       int iswap = 0;
       for (dim = 0; dim < 3; dim++)
       {
         for (ineed = 0; ineed < 2*maxneed_[dim]; ineed++)
         {
           pbc_flag_[iswap] = 0;
           vectorZeroizeN(pbc_[iswap],6);

           // send left, receive right
           if (ineed % 2 == 0)
           {
               sendproc_[iswap] = this->comm->procneigh[dim][0];
               recvproc_[iswap] = this->comm->procneigh[dim][1];

               if (ineed < 2) slablo_[iswap] = -BIG_MNMP;
               else slablo_[iswap] = 0.5 * (this->domain->sublo[dim] + this->domain->subhi[dim]);

               // use half cut here, since rBound is used (added) in checkBorderElement()
               slabhi_[iswap] = this->domain->sublo[dim] + half_atom_cut_;

               if (myloc[dim] == 0)
               {
                   pbc_flag_[iswap] = 1;
                   pbc_[iswap][dim] = 1;
               }
           }
           // send right, receive left
           else
           {
               sendproc_[iswap] = this->comm->procneigh[dim][1];
               recvproc_[iswap] = this->comm->procneigh[dim][0];

               // use half cut here, since rBound is used (added) in checkBorderElement()
               slablo_[iswap] = this->domain->subhi[dim] -  half_atom_cut_;
               if (ineed < 2) slabhi_[iswap] = BIG_MNMP;
               else slabhi_[iswap] = 0.5 * (this->domain->sublo[dim] + this->domain->subhi[dim]);

               if (myloc[dim] == this->comm->procgrid[dim]-1)
               {
                   pbc_flag_[iswap] = 1;
                   pbc_[iswap][dim] = -1;
               }
           }
           iswap++;
         }
       }
   }

/* ----------------------------------------------------------------------
   realloc the buffers needed for communication and swaps
------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::grow_swap(int n)
  {
      free_swap();
      allocate_swap(n);

      sendlist_ = (int **)
        this->memory->srealloc(sendlist_,n*sizeof(int *),"MultiNodeMeshParallel:sendlist_");
      sendwraplist_ = (int **)
        this->memory->srealloc(sendwraplist_,n*sizeof(int *),"MultiNodeMeshParallel:sendwraplist_");
      this->memory->grow(maxsendlist_,n,"MultiNodeMeshParallel:maxsendlist_");
      for (int i = maxswap_; i < n; i++)
      {
        maxsendlist_[i] = BUFMIN_MNMP;
        this->memory->create(sendlist_[i],BUFMIN_MNMP,"MultiNodeMeshParallel:sendlist_[i]");
        this->memory->create(sendwraplist_[i],BUFMIN_MNMP,"MultiNodeMeshParallel:sendwraplist_[i]");
      }
      maxswap_ = n;
  }

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::allocate_swap(int n)
  {
      this->memory->create(sendnum_,n,"MultiNodeMeshParallel:sendnum_");
      this->memory->create(recvnum_,n,"MultiNodeMeshParallel:recvnum_");
      this->memory->create(sendproc_,n,"MultiNodeMeshParallel:sendproc_");
      this->memory->create(recvproc_,n,"MultiNodeMeshParallel:recvproc_");
      this->memory->create(size_forward_recv_,n,"MultiNodeMeshParallel:size");
      this->memory->create(size_reverse_recv_,n,"MultiNodeMeshParallel:size");
      this->memory->create(slablo_,n,"MultiNodeMeshParallel:slablo_");
      this->memory->create(slabhi_,n,"MultiNodeMeshParallel:slabhi_");
      this->memory->create(firstrecv_,n,"MultiNodeMeshParallel:firstrecv");
      this->memory->create(pbc_flag_,n,"MultiNodeMeshParallel:pbc_flag_");
      this->memory->create(pbc_,n,6,"MultiNodeMeshParallel:pbc_");
  }

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::free_swap()
  {
      this->memory->destroy(sendnum_);
      this->memory->destroy(recvnum_);
      this->memory->destroy(sendproc_);
      this->memory->destroy(recvproc_);
      this->memory->destroy(size_forward_recv_);
      this->memory->destroy(size_reverse_recv_);
      this->memory->destroy(slablo_);
      this->memory->destroy(slabhi_);
      this->memory->destroy(firstrecv_);
      this->memory->destroy(pbc_flag_);
      this->memory->destroy(pbc_);
  }

  /* ----------------------------------------------------------------------
   realloc the size of the send buffer as needed with BUFFACTOR & BUFEXTRA
   if flag = 1, realloc
   if flag = 0, don't need to realloc with copy, just free/malloc
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::grow_send(int n, int flag)
  {
      maxsend_ = static_cast<int> (BUFFACTOR_MNMP * n);
      
      if (flag)
        this->memory->grow(buf_send_,(maxsend_+BUFEXTRA_MNMP),"MultiNodeMeshParallel:buf_send");
      else {
        this->memory->destroy(buf_send_);
        this->memory->create(buf_send_,maxsend_+BUFEXTRA_MNMP,"MultiNodeMeshParallel:buf_send");
      }
  }

  /* ----------------------------------------------------------------------
   free/malloc the size of the recv buffer as needed with BUFFACTOR
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::grow_recv(int n)
  {
      maxrecv_ = static_cast<int> (BUFFACTOR_MNMP * n);
      this->memory->destroy(buf_recv_);
      this->memory->create(buf_recv_,maxrecv_,"MultiNodeMeshParallel:buf_recv");
  }

  /* ----------------------------------------------------------------------
   realloc the size of the iswap sendlist as needed with BUFFACTOR
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::grow_list(int iswap, int n)
  {
    maxsendlist_[iswap] = static_cast<int> (BUFFACTOR_MNMP * n)+1;
    this->memory->grow(sendlist_[iswap],maxsendlist_[iswap],"MultiNodeMeshParallel:sendlist[iswap]");
    this->memory->grow(sendwraplist_[iswap],maxsendlist_[iswap],"MultiNodeMeshParallel:sendlist[iswap]");
  }

  /* ----------------------------------------------------------------------
   parallelization -
   initially, all processes have read the whole data
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::initialSetup()
  {
      nGlobalOrig_ = sizeLocal();

      // check for possible round-off isues
      
      double span = this->node_.max_scalar()-this->node_.min_scalar();
      if(span < 1e-4)
        this->error->all(FLERR,"Mesh error - root causes: (a) mesh empty or (b) dimensions too small - use different unit system");

      double comBefore[3];
      this->center_of_mass(comBefore);
      
      // delete all elements that do not belong to this processor
      
      deleteUnowned();

      if(sizeGlobal() != sizeGlobalOrig())
      {
        
        char errstr[1024];

        if(0 == sizeGlobal())
        {
            sprintf(errstr,"Mesh (id %s): All %d mesh elements have been lost / left the domain. \n"
                           "Please use 'boundary m m m' or scale/translate/rotate the mesh or change its dynamics\n"
                           "FYI: center of mass of mesh including scale/tranlate/rotate is %f / %f / %f\n"
                           "     simulation box x from %f to %f y  from %f to %f z from %f to %f\n"
                           "     (gives indication about changes in scale/tranlate/rotate necessary to make simulation run)\n",
                       this->mesh_id_,sizeGlobalOrig()-sizeGlobal(),comBefore[0],comBefore[1],comBefore[2],
                       this->domain->boxlo[0],this->domain->boxhi[0],this->domain->boxlo[1],this->domain->boxhi[1],this->domain->boxlo[2],this->domain->boxhi[2]);
        }
        else
        {
            double comAfter[3];
            this->center_of_mass(comAfter);

            sprintf(errstr,"Mesh (id %s): %d mesh elements have been lost / left the domain. \n"
                           "Please use 'boundary m m m' or scale/translate/rotate the mesh or change its dynamics\n"
                           "FYI: center of mass of mesh including scale/tranlate/rotate before cutting out elements is %f / %f / %f\n"
                           "     simulation box x from %f to %f y  from %f to %f z from %f to %f\n"
                           "     center of mass of mesh after cutting out elements outside simulation box is is        %f / %f / %f\n"
                           "     (gives indication about changes in scale/tranlate/rotate necessary to make simulation run)\n",
                       this->mesh_id_,sizeGlobalOrig()-sizeGlobal(),comBefore[0],comBefore[1],comBefore[2],
                       this->domain->boxlo[0],this->domain->boxhi[0],this->domain->boxlo[1],this->domain->boxhi[1],this->domain->boxlo[2],this->domain->boxhi[2],
                       comAfter[0],comAfter[1],comAfter[2]);
        }
        this->error->all(FLERR,errstr);
      }

      // perform operations that should be done before initial setup
      
      preInitialSetup();

      // set-up mesh parallelism
      
      setup();

      // re-calculate properties for owned particles
      
      refreshOwned(1);

      // identify elements that are near borders
      // forward communicate them
      
      borders();

      // re-calculate properties for ghost particles
      
      refreshGhosts(1);

      // build mesh topology and neigh list
      
      buildNeighbours();

      // perform quality check on the mesh
      
      qualityCheck();

      if(doParallellization_) isParallel_ = true;

      postInitialSetup();

      // stuff that should be done before resuming simulation
      
      postBorders();
      
  }

  /* ----------------------------------------------------------------------
   parallelization - aggregates pbc, exchange and borders
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::pbcExchangeBorders(int setupFlag)
  {
      // need not do this during simulation for non-moving mesh and non-changing simulation box
      
      if(setupFlag) this->reset_stepLastReset();

      // perform operations that should be done before setting up parallellism and exchanging elements
      preSetup();

      if(!setupFlag && !this->isMoving() && !this->isDeforming() && !this->domain->box_change) return;

      // set-up mesh parallelism
      setup();

      // enforce pbc
      pbc();

      // communicate particles
      exchange();

      if(sizeGlobal() != sizeGlobalOrig())
      {
        
        //this->error->all(FLERR,"Mesh elements have been lost");
        char errstr[500];
        sprintf(errstr,"Mesh (id %s): Mesh elements have been lost / left the domain. Please use "
                       "'boundary m m m' or scale/translate/rotate the mesh or change its dynamics",
                       this->mesh_id_);
        this->error->all(FLERR,errstr);
      }

      // re-calculate properties for owned particles
      
      refreshOwned(setupFlag);

      // identify elements that are near borders
      // forward communicate them
      
      borders();

      // re-calculate properties for ghosts
      refreshGhosts(setupFlag);

      // stuff that should be done before resuming simulation
      postBorders();

  }

  /* ----------------------------------------------------------------------
   parallelization - clear data of reverse comm properties
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::clearReverse()
  {
      // nothing to do here
  }

  /* ----------------------------------------------------------------------
   delete all particles which are not owned on this proc
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::deleteUnowned()
  {
      
      int i = 0;

      if(doParallellization_)
      {

          while(i < nLocal_)
          {
              if(!this->domain->is_in_subdomain(this->center_(i)))
                  this->deleteElement(i);
              else i++;
          }

          // calculate nGlobal for the first time
          MPI_Sum_Scalar(nLocal_,nGlobal_,this->world);
      }
      else
        nGlobal_ = nLocal_;

  }

  /* ----------------------------------------------------------------------
   enforce periodic boundary conditions
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::pbc()
  {
      if(!doParallellization_) return;

      double centerNew[3], delta[3];

      for(int i = 0; i < this->sizeLocal(); i++)
      {
          vectorCopy3D(this->center_(i),centerNew);
          this->domain->remap(centerNew);
          vectorSubtract3D(centerNew,this->center_(i),delta);

          // move element i incremental
          if(vectorMag3DSquared(delta) > 1e-9)
            this->moveElement(i,delta);
      }
  }

  /* ----------------------------------------------------------------------
   exchange elements with nearby processors
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::exchange()
  {
      if(!doParallellization_) return;

      int nrecv, nsend = 0;
      int nrecv1,nrecv2;
      double *buf;
      MPI_Request request;
      MPI_Status status;
      MPI_Comm world = this->world;

      //int nprocs = this->comm->nprocs;
      int *procgrid = this->comm->procgrid;
      int procneigh[3][2];

      // clear global->local map for owned and ghost atoms
      
      clearMap();

      // clear old ghosts
      
      clearGhosts();

      // copy procneigh
      for (int i = 0; i < 3; i++)
        for( int j = 0; j < 2; j++)
            procneigh[i][j] = this->comm->procneigh[i][j];

      for (int dim = 0; dim < 3; dim++)
      {
          // push data to buffer
          
          nsend = pushExchange(dim);

          // send/recv in both directions
          // if 1 proc in dimension, no send/recv, set recv buf to send buf
          // if 2 procs in dimension, single send/recv
          // if more than 2 procs in dimension, send/recv to both neighbors

          if (procgrid[dim] == 1)
          {
            nrecv = nsend;
            buf = buf_send_;
          }
          else
          {
            MPI_Sendrecv(&nsend,1,MPI_INT,procneigh[dim][0],0,&nrecv1,1,MPI_INT,procneigh[dim][1],0,world,&status);
            nrecv = nrecv1;

            if (this->comm->procgrid[dim] > 2)
            {
                MPI_Sendrecv(&nsend,1,MPI_INT,procneigh[dim][1],0,&nrecv2,1,MPI_INT,procneigh[dim][0],0,world,&status);
                nrecv += nrecv2;
            }

            if (nrecv > maxrecv_) grow_recv(nrecv);

            MPI_Irecv(buf_recv_,nrecv1,MPI_DOUBLE,procneigh[dim][1],0,world,&request);
            MPI_Send(buf_send_,nsend,MPI_DOUBLE,procneigh[dim][0],0,world);
            MPI_Wait(&request,&status);

            if (procgrid[dim] > 2)
            {
                MPI_Irecv(&buf_recv_[nrecv1],nrecv2,MPI_DOUBLE,procneigh[dim][0],0,world,&request);
                MPI_Send(buf_send_,nsend,MPI_DOUBLE,procneigh[dim][1],0,world);
                MPI_Wait(&request,&status);
            }

            buf = buf_recv_;
          }

          // check incoming elements to see if they are in my box
          // if so, add on this proc

          popExchange(nrecv,dim, buf);
          
      }

      // re-calculate nGlobal as some element might have been lost
     MPI_Sum_Scalar(nLocal_,nGlobal_,world);
  }

  /* ----------------------------------------------------------------------
   generate ghost elements, refresh global map
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::borders()
  {
      if(doParallellization_)
      {
          int iswap, twoneed, nfirst, nlast, n, nsend, nrecv, smax, rmax;
          bool sendflag, dummy = false;
          double lo,hi;
          MPI_Request request;
          MPI_Status status;

          nfirst = 0;
          iswap = 0;
          smax = rmax = 0;

          for (int dim = 0; dim < 3; dim++)
          {
              nlast = 0;

              // need to go left and right in each dim
              twoneed = 2*maxneed_[dim];
              for (int ineed = 0; ineed < twoneed; ineed++)
              {
                  lo = slablo_[iswap];
                  hi = slabhi_[iswap];

                  // find elements within slab boundaries lo/hi using <= and >=
                  
                  if (ineed % 2 == 0)
                  {
                      nfirst = nlast;
                      nlast = sizeLocal() + sizeGhost();
                  }

                  nsend = 0;

                  // sendflag = 0 if I do not send on this swap
                  
                  sendflag = true;
                  int wrap = 0;
                  
                  if(ineed % 2 == 0 && this->comm->myloc[dim] == 0)
                  {
                      if(this->domain->periodicity[dim] && !this->domain->triclinic && !dynamic_cast<DomainWedge*>(this->domain))
                          wrap = 1;
                      else
                          sendflag = false;
                  }

                  if(ineed % 2 == 1 && this->comm->myloc[dim] == this->comm->procgrid[dim]-1)
                  {
                      if(this->domain->periodicity[dim] && !this->domain->triclinic && !dynamic_cast<DomainWedge*>(this->domain))
                          wrap = -1;
                      else
                          sendflag = false;
                  }

                  // find send elements
                  if(sendflag)
                  {
                      
                      for (int i = nfirst; i < nlast; i++)
                      {
                          int type = checkBorderElement(ineed, i, dim, lo, hi);
                          if(type != NOT_GHOST)
                          {
                              if (nsend >= maxsendlist_[iswap])
                                  grow_list(iswap,nsend);
                              sendlist_[iswap][nsend] = i;
                              if (wrap == 1)
                              {
                                  switch (dim)
                                  {
                                  case 0:
                                      type = IS_GHOST_WRAP_DIM_0_POS;
                                      break;
                                  case 1:
                                      type = IS_GHOST_WRAP_DIM_1_POS;
                                      break;
                                  case 2:
                                      type = IS_GHOST_WRAP_DIM_2_POS;
                                      break;
                                  }
                              }
                              else if (wrap == -1)
                              {
                                  switch (dim)
                                  {
                                  case 0:
                                      type = IS_GHOST_WRAP_DIM_0_NEG;
                                      break;
                                  case 1:
                                      type = IS_GHOST_WRAP_DIM_1_NEG;
                                      break;
                                  case 2:
                                      type = IS_GHOST_WRAP_DIM_2_NEG;
                                      break;
                                  }
                              }
                              sendwraplist_[iswap][nsend] = type;
                              nsend++;

                          }
                      }
                  }

                  // pack up list of border elements

                  if(nsend*size_border_ > maxsend_)
                    grow_send(nsend*size_border_,0);

                  n = pushElemListToBuffer(nsend, sendlist_[iswap], sendwraplist_[iswap], buf_send_, OPERATION_COMM_BORDERS, NULL, this->domain->boxlo, this->domain->boxhi,dummy,dummy,dummy);

                  // swap atoms with other proc
                  // no MPI calls except SendRecv if nsend/nrecv = 0
                  // put incoming ghosts at end of my atom arrays
                  // if swapping with self, simply copy, no messages

                  double *buf = NULL;

                  if (sendproc_[iswap] != this->comm->me)
                  {
                      MPI_Sendrecv(&nsend,1,MPI_INT,sendproc_[iswap],0,&nrecv,1,MPI_INT,recvproc_[iswap],0,this->world,&status);
                      if (nrecv*size_border_ > maxrecv_)
                          grow_recv(nrecv*size_border_);
                      if (nrecv)
                          MPI_Irecv(buf_recv_,nrecv*size_border_,MPI_DOUBLE,recvproc_[iswap],0,this->world,&request);

                      if (n)
                          MPI_Send(buf_send_,n,MPI_DOUBLE,sendproc_[iswap],0,this->world);

                      if (nrecv)
                          MPI_Wait(&request,&status);

                      buf = buf_recv_;
                  }
                  else
                  {
                      nrecv = nsend;
                      buf = buf_send_;
                  }

                  // unpack buffer

                  n = popElemListFromBuffer(nLocal_+nGhost_, nrecv, buf, OPERATION_COMM_BORDERS, NULL, dummy,dummy,dummy);

                  // set pointers & counters

                  smax = MAX(smax,nsend);
                  rmax = MAX(rmax,nrecv);
                  sendnum_[iswap] = nsend;
                  recvnum_[iswap] = nrecv;
                  size_forward_recv_[iswap] = nrecv*size_forward_;
                  size_reverse_recv_[iswap] = nsend*size_reverse_;
                  firstrecv_[iswap] = nLocal_+nGhost_;
                  nGhost_ += nrecv;
                  iswap++;
              }
          }

          // insure send/recv buffers are long enough for all forward & reverse comm
          int max = MAX(maxforward_*smax,maxreverse_*rmax);
          if (max > maxsend_) grow_send(max,0);
          max = MAX(maxforward_*rmax,maxreverse_*smax);
          if (max > maxrecv_) grow_recv(max);
      }

      // build global-local map
      this->generateMap();
  }

  /* ----------------------------------------------------------------------
   check if element qualifies as ghost
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  inline int MultiNodeMeshParallel<NUM_NODES>::checkBorderElement(const int ineed, const int i, const int dim, const double lo, const double hi) const
  {
      if (ineed % 2 == 0)
          return checkBorderElementLeft(i,dim,lo,hi);
      else
          return checkBorderElementRight(i,dim,lo,hi);
  }
  
  template<int NUM_NODES>
  inline int MultiNodeMeshParallel<NUM_NODES>::checkBorderElementLeft(const int i, const int dim, const double lo, const double hi) const
  {
      
      // center of triangle
      const double pos = this->center_(i)[dim];
      // hi is extended by the bounding radius
      const double hi_extended = hi + this->rBound_(i);
      // check whether center is inside interval
      if (pos >= lo && pos <= hi_extended)
        return IS_GHOST;

      return NOT_GHOST;
  }

  template<int NUM_NODES>
  inline int MultiNodeMeshParallel<NUM_NODES>::checkBorderElementRight(const int i, const int dim, const double lo, const double hi) const
  {
      // center of triangle
      const double pos = this->center_(i)[dim];
      // lo is extended by the bounding radius
      const double lo_extended = lo - this->rBound_(i);
      // check whether center is inside interval
      if (pos >= lo_extended && pos <= hi)
        return IS_GHOST;

      return NOT_GHOST;
  }

  /* ----------------------------------------------------------------------
   communicate properties to ghost elements
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::forwardComm(std::string property)
  {
      std::list<std::string> properties (1, property);
      forwardComm(&properties);
  }

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::forwardComm(std::list<std::string> * properties)
  {
      int n;
      MPI_Request request;
      MPI_Status status;
      int me = this->comm->me;

      bool scale = this->isScaling();
      bool translate = this->isTranslating();
      bool rotate = this->isRotating();

      // exit here if no forward communication at all
      
      if(size_forward_ == 0)
        return;

      const int size_this = properties ? elemBufSize(OPERATION_COMM_REVERSE, properties, scale, translate, rotate) : 1;

      // exchange data with another proc
      // if other proc is self, just copy

      for (int iswap = 0; iswap < nswap_; iswap++)
      {
          if (sendproc_[iswap] != me)
          {
                if (size_forward_recv_[iswap] && size_this)
                {
                    int nrecv = size_forward_recv_[iswap];
                    if (properties)
                    {
                        // size forward is the size of all forward buffers
                        nrecv /= size_forward_;
                        // size_this is the size of the buffers listed in properties
                        nrecv *= size_this;
                    }
                    MPI_Irecv(buf_recv_, nrecv, MPI_DOUBLE,recvproc_[iswap],0,this->world,&request);
                }

                n = pushElemListToBuffer(sendnum_[iswap],sendlist_[iswap], sendwraplist_[iswap],buf_send_,OPERATION_COMM_FORWARD, properties, this->domain->boxlo, this->domain->boxhi,scale,translate,rotate);
                
                if (n)
                    MPI_Send(buf_send_,n,MPI_DOUBLE,sendproc_[iswap],0,this->world);

                if (size_forward_recv_[iswap] && size_this)
                    MPI_Wait(&request,&status);

                n = popElemListFromBuffer(firstrecv_[iswap],recvnum_[iswap],buf_recv_,OPERATION_COMM_FORWARD, properties, scale,translate,rotate);
                
          }
          else
          {
              n = pushElemListToBuffer(sendnum_[iswap], sendlist_[iswap], sendwraplist_[iswap], buf_send_, OPERATION_COMM_FORWARD, properties, this->domain->boxlo, this->domain->boxhi, scale, translate, rotate);

              // note buf_recv_ not used in this case (just use buf_send_ as receive buffer
              n = popElemListFromBuffer(firstrecv_[iswap], recvnum_[iswap], buf_send_, OPERATION_COMM_FORWARD, properties, scale, translate, rotate);
          }
      }
  }

  /* ----------------------------------------------------------------------
   reverse communication of properties on atoms every timestep
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::reverseComm(std::string property)
  {
      std::list<std::string> properties (1, property);
      reverseComm(&properties);
  }

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::reverseComm(std::list<std::string> * properties)
  {
      int n;
      MPI_Request request;
      MPI_Status status;
      int me = this->comm->me;

      bool scale = this->isScaling();
      bool translate = this->isTranslating();
      bool rotate = this->isRotating();

      const int size_this = properties ? elemBufSize(OPERATION_COMM_REVERSE, properties, scale, translate, rotate) : 1;

      // exchange data with another proc
      // if other proc is self, just copy

      for (int iswap = nswap_-1; iswap >= 0; iswap--)
      {
          if (sendproc_[iswap] != me)
          {
              if (size_reverse_recv_[iswap] && size_this)
              {
                  int nrecv = size_reverse_recv_[iswap];
                  if (properties)
                  {
                      // size reverse is the size of all reverse buffers
                      nrecv /= size_reverse_;
                      // size_this is the size of the buffers listed in properties
                      nrecv *= size_this;
                  }
                  MPI_Irecv(buf_recv_, nrecv, MPI_DOUBLE, sendproc_[iswap], 0, this->world, &request);
              }

              n = pushElemListToBufferReverse(firstrecv_[iswap], recvnum_[iswap], buf_send_, OPERATION_COMM_REVERSE, properties, scale, translate, rotate);

              if (n)
                  MPI_Send(buf_send_,n,MPI_DOUBLE,recvproc_[iswap],0,this->world);

              if (size_reverse_recv_[iswap] && size_this)
                  MPI_Wait(&request,&status);

              n = popElemListFromBufferReverse(sendnum_[iswap], sendlist_[iswap], buf_recv_, OPERATION_COMM_REVERSE, properties, scale, translate, rotate);
          }
          else
          {
              n = pushElemListToBufferReverse(firstrecv_[iswap], recvnum_[iswap], buf_send_, OPERATION_COMM_REVERSE, properties, scale, translate, rotate);
              // note buf_recv_ not used in this case (just use buf_send_ as receive buffer
              n = popElemListFromBufferReverse(sendnum_[iswap], sendlist_[iswap], buf_send_, OPERATION_COMM_REVERSE, properties, scale, translate, rotate);
          }
      }
  }

#endif
