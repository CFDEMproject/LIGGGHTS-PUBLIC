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

#ifndef LMP_COMM_I_H
#define LMP_COMM_I_H

#include "atom.h"
#include "domain_wedge.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   decide if use comm optimizations for granular systems
   don't use for triclinic (current implementation not valid for triclinic,
   would have to translate radius into triclinic coordinates)
------------------------------------------------------------------------- */

inline bool Comm::use_gran_opt()
{
    return (0 == domain->triclinic && atom->radius);
}

/* ----------------------------------------------------------------------
   decide if border element, optimization for granular
------------------------------------------------------------------------- */

inline bool Comm::decide(int i,int dim,double lo,double hi,int ineed)
{
    double **x = atom->x;
    double *radius = atom->radius;

    if( ((ineed % 2 == 0) && x[i][dim] >= lo && x[i][dim] <= (hi + (use_gran_opt()? (radius[i]) : 0.)) ) ||
        ((ineed % 2 == 1) && x[i][dim] >= (lo - (use_gran_opt()? radius[i] : 0.)) && x[i][dim] <= hi )   )
        return true;

    return false;
}

/* ----------------------------------------------------------------------
   decide if border element for wedge case, optimization for granular
------------------------------------------------------------------------- */

inline bool Comm::decide_wedge(int i,int dim,double lo,double hi,int ineed)
{
    double **x = atom->x;
    double *radius = atom->radius;
    double coo[2],d[2];
    coo[0] = x[i][iphi];
    coo[1] = x[i][(iphi+1)%3];

    if (ineed % 2 == 0)
    {
        vectorSubtract2D(coo,pleft,d);
        if(vectorDot2D(d,nleft) >= -(use_gran_opt()? radius[i] : 0.))
        {
            
            return true;
        }
    }
    
    else if (ineed % 2 == 1)
    {
        vectorSubtract2D(coo,pright,d);
        if(vectorDot2D(d,nright) >= -(use_gran_opt()? radius[i] : 0.))
        {
            
            return true;
        }
    }
    
    return false;
}

/* ----------------------------------------------------------------------
   routine to determine Exchange Events
   similar to Comm:Exchange()
------------------------------------------------------------------------- */

void Comm::exchangeEventsRecorder()
{
    if(!exchangeEvents || nprocs==1 )
        return;

    exchangeEventsLocalId.clear();
    exchangeEventsReceivingProcess.clear();
    exchangeEventsGlobalProblemIds.clear();

    bool verbose = false; //true; //Developer to set here if needed for debugging
    int i;
    double **x;
    int    *tag;
    double *sublo,*subhi;
    double subloNeigh[3],subhiNeigh[3]; //boarders of 1st neighbor process (0=left if two), in each dim

    MPI_Status  status;
    MPI_Request request;

    //determine process boundaries
    //since called in Comm::Exchange, triclinic boundaries are already correct
    if (triclinic == 0) {
      sublo = domain->sublo;
      subhi = domain->subhi;
    } else {
      sublo = domain->sublo_lamda;
      subhi = domain->subhi_lamda;
    }

    //Exchange info on sublo/subhi
    for(int dim=0;dim<3;dim++)
    {
       MPI_Send (&(sublo[dim]),      1, MPI_DOUBLE, procneigh[dim][1],0,world);          //to right
       MPI_Irecv(&(subloNeigh[dim]), 1, MPI_DOUBLE, procneigh[dim][0],0,world,&request); //from left
       MPI_Wait (&request,&status); //wait to receive

       MPI_Send (&(subhi[dim]),      1, MPI_DOUBLE, procneigh[dim][1],0,world);          //to right
       MPI_Irecv(&(subhiNeigh[dim]), 1, MPI_DOUBLE, procneigh[dim][0],0,world,&request); //from left
       MPI_Wait (&request,&status); //wait to receive
    }
    if(verbose)
       printf("[%d/%d]:Comm::exchangeEventsRecorder(); lo/hi: %g %g %g / %g %g %g, exchanging sublo/hi information with process %d %d %d: %g %g %g / %g %g %g \n",
               me, nprocs,
               sublo[0],sublo[1],sublo[2],subhi[0],subhi[1],subhi[2],
               procneigh[0][0], procneigh[1][0], procneigh[2][0],
               subloNeigh[0],subloNeigh[1],subloNeigh[2],
               subhiNeigh[0],subhiNeigh[1],subhiNeigh[2]);

    // loop over dimensions,
    x   = atom->x;
    tag = atom->tag;
    i = 0;
    while (i < atom->nlocal)
    {
        bool  willExchange    = false;
        bool  problemDetected = false;
//        if(verbose)
//           printf("[%d/%d]:Comm::exchangeEventsRecorder() looping  atom %d with pos: %g %g %g\n",
//                       me, nprocs, i, x[i][0], x[i][1], x[i][2]);

        for (int dim = 0; dim < 3; dim++)
        {
            if( procgrid[dim] == 1 ) //nothing to do if only one processor in this direction
              continue;

            double lo = sublo[dim];
            double hi = subhi[dim];
            double position = x[i][dim];
            bool   willExchangeThisDim = false;

            if(verbose)
              printf("[%d/%d]:Comm::exchangeEventsRecorder() looping dim: %d, atom %d with pos: %g\n",
                       me, nprocs, dim, i, position);

            if (position < lo || position >= hi)
            {
                if(willExchange) //exchange over multiple dims, has been already registered!
                {
                    exchangeEventsReceivingProcess.back() = -1; //reset receiving processor
                    problemDetected = true;
                    continue;
                }

                //Append results from this direction to overall event list
                //can have only one per local ID!
                willExchange        = true;
                willExchangeThisDim = true;
                exchangeEventsLocalId.push_back(i);
                exchangeEventsReceivingProcess.push_back(procneigh[dim][0]); //by default, use first neighbor, correct later

                if(verbose)
                printf("[%d/%d]:Comm::exchangeEventsRecorder() recording exchange for atom %d with pos: %g\n",
                       me, nprocs, i, position );
            }

            if( (procgrid[dim]>2) && willExchangeThisDim )
            {
                if(verbose)
                    printf("[%d/%d]: checking dim: %d because procgrid[dim] = %d \n",
                       me, nprocs, dim, procgrid[dim]);

                if( position < subloNeigh[dim] || position >= subhiNeigh[dim] ) //not on left neighbor
                {
                    exchangeEventsReceivingProcess.back() = procneigh[dim][1];
                    if(verbose)
                        printf("[%d/%d]: checking dim: %d, detected transfer to right process with id %d \n",
                               me, nprocs, dim, exchangeEventsReceivingProcess.back());
                }
            } //handle more than 2 processors in this dimension
        } //loop over dimensions

        //if problem detected, regist globalId for later determination of receiving process
        if(problemDetected)
           exchangeEventsGlobalProblemIds.push_back(tag[i]);

        i++;

    } //loop over atoms
}

/* ----------------------------------------------------------------------
   routine to determine the id of the process a particle has been
   transferred to in case of multi process dimensions crossing
------------------------------------------------------------------------- */

void Comm::exchangeEventsCorrector()
{
    if( !exchangeEvents || nprocs==1 )
        return;

    //Determine global count of problemIds, exit if none
    bool verbose = false; //true; //Developer to set here if needed for debugging
    int  global_sum = 0;
    int  local_sum = exchangeEventsGlobalProblemIds.size();
    MPI_Barrier(world);
    MPI_Allreduce(&local_sum, &global_sum,
                  1, MPI_INT, MPI_SUM, world
                 );
    if( global_sum==0 )
        return;

    //Collect information on problem IDs
    int * exchangeReceiveCounts;
    int * exchangeReceiveDisplacement;
    int * exchangeGlobalProblems;
    int * exchangeOwningProcess;
    memory->create(exchangeReceiveCounts,nprocs,"comm:exchangeReceiveCounts"); //TODO: do just once
    memory->create(exchangeReceiveDisplacement,nprocs,"comm:exchangeReceiveDisplacement");  //TODO
    memory->create(exchangeGlobalProblems,global_sum,"comm:exchangeGlobalProblems");
    memory->create(exchangeOwningProcess, global_sum,"comm:exchangeOwningProcess");

  	MPI_Allgather(&local_sum,            1, MPI_INT,
   	               exchangeReceiveCounts, 1, MPI_INT,
                   world
                 );

    exchangeReceiveDisplacement[0] = 0;
    for(int iPro=1; iPro<nprocs; iPro++)
	    exchangeReceiveDisplacement[iPro] = exchangeReceiveCounts[iPro-1]
                                          + exchangeReceiveDisplacement[iPro-1];

    //Inform all CPUs about problem IDs
    MPI_Barrier(world);
    MPI_Allgatherv(&(exchangeEventsGlobalProblemIds.front()), local_sum, MPI_INT,
                   exchangeGlobalProblems, exchangeReceiveCounts, exchangeReceiveDisplacement,
                   MPI_INT, world
                  );

    if(verbose && me==0)
    {
        printf("**exchangeEventsCorrector: globalProblems %d with IDs: ",
                global_sum
              );
        for(int k=0;k<global_sum;k++)
          printf(" %d, ",exchangeGlobalProblems[k] );

        printf("; exchangeReceiveCounts:");

        for(int k=0;k<nprocs;k++)
          printf(" %d, ",
                exchangeReceiveCounts[k]
              );
        printf("\n");
    }

    //Loop through list of global problem IDs and detect if CPU owns this atom
    for(int iter=0; iter<global_sum; iter++)
    {
        exchangeOwningProcess[iter] = -1;
        int currLocal     = atom->map(exchangeGlobalProblems[iter]);
        if( (currLocal<0) || (currLocal >= atom->nlocal) ) //not owned
            continue;
        else
            exchangeOwningProcess[iter] = me;
    }

    MPI_Barrier(world);
    MPI_Allreduce(exchangeOwningProcess,
                  exchangeGlobalProblems, //re-use this container, is now (received) owning process id
                  global_sum, MPI_INT, MPI_MAX, world
                 );

    //Check if all problem IDs have been detected
    for(int j=0;j<global_sum;j++)
    {
       if(verbose && me==0)
           printf("**exchangeEventsCorrector: receivedOwningProcess[%d]: %d \n",
                   j, exchangeGlobalProblems[j]
                 );

       if(exchangeGlobalProblems[j]<0)
            error->one(FLERR,"Comm::exchangeEventsCorrector: Could not find particle. Must have left the domain. You must upgrade this check in order to handle this situation.");
    }

    //Fill in the process ids
    unsigned int checkCounter=0;
    for(int iter =  exchangeReceiveDisplacement[me];  //loop global list to get correct iter
            iter < (exchangeReceiveDisplacement[me]+exchangeReceiveCounts[me]);
            iter++
       )
       for(unsigned int j=0; j<exchangeEventsReceivingProcess.size(); j++)
         if(exchangeEventsReceivingProcess[j]==-1) //have invalid receiving process id
         {
            exchangeEventsReceivingProcess[j]=exchangeGlobalProblems[iter];
            checkCounter++;
         }

    if( checkCounter!=exchangeEventsGlobalProblemIds.size() )
        error->all(FLERR,"Comm::exchangeEventsCorrector: Problem when fill in corrected process ids! This is fatal.");

    //Clear memory
    memory->destroy(exchangeReceiveCounts);
    memory->destroy(exchangeReceiveDisplacement);
    memory->destroy(exchangeGlobalProblems);
    memory->destroy(exchangeOwningProcess);
}

#endif
