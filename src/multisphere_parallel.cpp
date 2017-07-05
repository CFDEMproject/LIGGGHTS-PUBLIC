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

#define DELTA 10000

#include "multisphere_parallel.h"
#include "atom.h"
#include "atom_vec.h"
#include "vector_liggghts.h"
#include "domain.h"
#include "memory.h"

#define BUFFACTOR 1.5
#define BUFMIN 1000
#define BUFEXTRA 1000

/* ----------------------------------------------------------------------
   constructor / destructor
------------------------------------------------------------------------- */

MultisphereParallel::MultisphereParallel(LAMMPS *lmp) :
  Multisphere(lmp),

  // initialize comm buffers & exchange memory
  maxsend_(BUFMIN),
  maxrecv_(BUFMIN),
  buf_send_((double *) memory->smalloc((maxsend_+BUFEXTRA)*sizeof(double),"frm:buf_send_")),
  buf_recv_((double *) memory->smalloc((maxsend_+BUFEXTRA)*sizeof(double),"frm:buf_send_"))
{

}

MultisphereParallel::~MultisphereParallel()
{
    memory->sfree(buf_send_);
    memory->sfree(buf_recv_);
}

/* ----------------------------------------------------------------------
   realloc the size of the send buffer as needed with BUFFACTOR & BUFEXTRA
   if flag = 1, realloc
   if flag = 0, don't need to realloc with copy, just free/malloc
------------------------------------------------------------------------- */

void MultisphereParallel::grow_send(int n, int flag)
{
  maxsend_ = static_cast<int> (BUFFACTOR * n);
  if (flag)
    buf_send_ = (double *) memory->srealloc(buf_send_,(maxsend_+BUFEXTRA)*sizeof(double),"comm:buf_send_");
  else {
    memory->sfree(buf_send_);
    buf_send_ = (double *) memory->smalloc((maxsend_+BUFEXTRA)*sizeof(double),"comm:buf_send_");
  }
}

/* ----------------------------------------------------------------------
   free/malloc the size of the recv buffer as needed with BUFFACTOR
------------------------------------------------------------------------- */

void MultisphereParallel::grow_recv(int n)
{
  maxrecv_ = static_cast<int> (BUFFACTOR * n);
  memory->sfree(buf_recv_);
  buf_recv_ = (double *) memory->smalloc(maxrecv_*sizeof(double),        "comm:buf_recv_");
}

/* ----------------------------------------------------------------------
   exchange bodies with neighbor procs
------------------------------------------------------------------------- */

void MultisphereParallel::exchange()
{
  int i,m,nsend,nrecv,nrecv1,nrecv2;
  double lo,hi,value;
  double x[3];
  double *sublo,*subhi,*buf;
  MPI_Request request;
  MPI_Status status;

  // subbox bounds for orthogonal
  // triclinic not implemented

  sublo = domain->sublo;
  subhi = domain->subhi;

  // loop over dimensions

  for (int dim = 0; dim < 3; dim++) {

    // fill buffer with atoms leaving my box, using < and >=
    // when atom is deleted, fill it in with last atom

    lo = sublo[dim];
    hi = subhi[dim];
    i = nsend = 0;

    while (i < nbody_) {

          MathExtraLiggghts::local_coosys_to_cartesian(x,xcm_to_xbound_(i),ex_space_(i),ey_space_(i),ez_space_(i));
          vectorAdd3D(xcm_(i),x,x);

          if (x[dim] < lo || x[dim] >= hi)
          {
            if (nsend > maxsend_)
                grow_send(nsend,1);
            nsend += pack_exchange_rigid(i,&buf_send_[nsend]);
            
            remove_body(i);
            
          }
          else i++;
    }

    // send/recv atoms in both directions
    // if 1 proc in dimension, no send/recv, set recv buf to send buf
    // if 2 procs in dimension, single send/recv
    // if more than 2 procs in dimension, send/recv to both neighbors

    int procneigh[3][2];
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 2; j++)
            procneigh[i][j] = comm->procneigh[i][j];

    int *procgrid = comm->procgrid;

    if (procgrid[dim] == 1) {
      nrecv = nsend;
      buf = buf_send_;

    }
    else
    {
          MPI_Sendrecv(&nsend,1,MPI_INT,procneigh[dim][0],0,&nrecv1,1,MPI_INT,procneigh[dim][1],0,world,&status);
          nrecv = nrecv1;
          if (procgrid[dim] > 2) {
             MPI_Sendrecv(&nsend,1,MPI_INT,procneigh[dim][1],0,&nrecv2,1,MPI_INT,procneigh[dim][0],0,world,&status);
             nrecv += nrecv2;
          }

          if (nrecv > maxrecv_) grow_recv(nrecv);

          MPI_Irecv(buf_recv_,nrecv1,MPI_DOUBLE,procneigh[dim][1],0,world,&request);
          MPI_Send(buf_send_,nsend,MPI_DOUBLE,procneigh[dim][0],0,world);
          MPI_Wait(&request,&status);

          if (procgrid[dim] > 2) {
            MPI_Irecv(&buf_recv_[nrecv1],nrecv2,MPI_DOUBLE,procneigh[dim][0],0,world,&request);
            MPI_Send(buf_send_,nsend,MPI_DOUBLE,procneigh[dim][1],0,world);
            MPI_Wait(&request,&status);
          }

          buf = buf_recv_;
    }

    // check incoming atoms to see if they are in my box
    // if so, add to my list

    m = 0;

    while (m < nrecv) {
      value = buf[m+dim+1];
      
      if (value >= lo && value < hi) m += unpack_exchange_rigid(&buf[m]);
      else m += static_cast<int> (buf[m]);
    }
  }

}

/* ----------------------------------------------------------------------
   restart functionality - write all required data into restart buffer
   executed on all processes, but only proc 0 writes into writebuf
------------------------------------------------------------------------- */

void MultisphereParallel::writeRestart(FILE *fp)
{
    double *sendbuf = 0, *recvbuf = 0;
    double xbnd[3];
    bool dummy = false;
    double nba = static_cast<double>(n_body_all());

    int sizeLocal = n_body() * (customValues_.elemBufSize(OPERATION_RESTART, NULL, dummy,dummy,dummy) + 4);
    int sizeGlobal = 0, sizeOne = 0;

    // allocate send buffer and pack element data
    // all local elements are in list
    
    memory->create(sendbuf,sizeLocal,"MultiNodeMeshParallel::writeRestart:sendbuf");
    sizeLocal = 0;
    for(int i = 0; i < n_body(); i++)
    {
        x_bound(xbnd,i);
        sizeOne = customValues_.pushElemToBuffer(i,&(sendbuf[sizeLocal+4]),OPERATION_RESTART,dummy,dummy,dummy);
        sendbuf[sizeLocal] = static_cast<double>(sizeOne+4);
        sendbuf[sizeLocal+1] = xbnd[0];
        sendbuf[sizeLocal+2] = xbnd[1];
        sendbuf[sizeLocal+3] = xbnd[2];
        
        sizeLocal += (sizeOne+4);
    }

    // gather the per-element data
    
    sizeGlobal = MPI_Gather0_Vector(sendbuf,sizeLocal,recvbuf,world);

    // write data to file
    if(comm->me == 0)
    {
        
        // size with 1 extra value (nba)
        int size = (1+sizeGlobal) * sizeof(double);

        // write size
        fwrite(&size,sizeof(int),1,fp);

        // write extra value
        fwrite(&nba,sizeof(double),1,fp);

        // write per-element data
        fwrite(recvbuf,sizeof(double),sizeGlobal,fp);
    }

    // clean up

    memory->destroy(sendbuf);
    
    if(recvbuf)
      delete []recvbuf;
}

/* ----------------------------------------------------------------------
   restart functionality - read all required data from restart buffer
   executed on all processes
------------------------------------------------------------------------- */

void MultisphereParallel::restart(double *list)
{
    bool dummy = false;
    int m = 0, nrecv_this;

    int nbody_all_old = static_cast<int> (list[m++]);

    nbody_ = nbody_all_ = 0;

    for(int i = 0; i < nbody_all_old; i++)
    {
        nrecv_this = static_cast<int>(list[m]);
        
        double *x_bnd = &(list[m+1]);
        
        if(domain->is_in_subdomain(x_bnd))
        {
            
            customValues_.addZeroElement();
            customValues_.deleteRestartElement(nbody_,dummy,dummy,dummy);
            customValues_.popElemFromBuffer(&(list[m+4]),OPERATION_RESTART,dummy,dummy,dummy);

            nbody_++;
        }
        m += nrecv_this;
    }

    // do initialization tasks

    MPI_Sum_Scalar(nbody_,nbody_all_,world);
    generate_map();
    reset_forces(true);

}
