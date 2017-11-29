/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com


   Copyright (C): 2014 DCS Computing GmbH (www.dcs-computing.com), Linz, Austria
                  2014 Graz University of Technology (ippt.tugraz.at), Graz, Austria

   LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
   the producer of the LIGGGHTS® software and the CFDEM®coupling software
   See http://www.cfdem.com/terms-trademark-policy for details.

   LIGGGHTS® is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   Parts of the code were developped in the frame of the NanoSim project funded
   by the European Commission through FP7 Grant agreement no. 604656.
------------------------------------------------------------------------- */

#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "pascal.h"          // these are PASCAL include files
#include "fix_pascal_couple.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "modify.h"
#include "force.h"
#include "output.h"
#include "group.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "fix_property_atom.h"
#include "cfd_datacoupling_simple.h"

using namespace LAMMPS_NS;
using namespace PASCAL_NS;
using namespace FixConst;

/* ----------------------------------------------------------------------
   Constructor, calls ParScale's constructor
------------------------------------------------------------------------- */
FixParScaleCouple::FixParScaleCouple(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  exchangeEventsLocalId(NULL),
  exchangeEventsReceivingProcess(NULL),
  dc_(0),
  verbose_(false),
  reneighbor_at_least_every_(0),
  couple_this_step_(false),
  pascal_setup_(false),
  prePostRun_(false),
  time_(0.),
  pasc_(0)
{

  iarg_ = 3;
  dc_ = NULL;
  map_copy = NULL;

  int me, nprocs;
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

  if (narg < 7)
    error->fix_error(FLERR,this,"not enough arguments");

  if (strcmp(arg[iarg_],"reneighbor_at_least_every"))
    error->fix_error(FLERR,this,"expecting keyword 'reneighbor_at_least_every'");

  iarg_++;
  reneighbor_at_least_every_ = force->inumeric(FLERR,arg[iarg_]);
  if (reneighbor_at_least_every_ < 0)
    error->fix_error(FLERR,this,"'reneighbor_at_least_every' >= 0 required");

  iarg_++;
  if (strcmp(arg[iarg_],"couple_every"))
    error->fix_error(FLERR,this,"expecting keyword 'couple_every'");

  iarg_++;
  couple_every_ = force->inumeric(FLERR,arg[iarg_]);
  if (couple_every_ != 1)
    error->fix_error(FLERR,this,"'couple_every' = 1 required");


  iarg_++;
  if (narg > 7)
  {
    if (strcmp(arg[iarg_],"verbose")==0)
        verbose_ = true;
  }
  if (narg > 8)
  {
    if (strcmp(arg[iarg_],"prePostRun")==0)
        prePostRun_ = true;
  }

  nevery = 1;

  //Create datacoupling
  dc_ = new CfdDatacouplingSimple(lmp,iarg_+1,narg,arg,this);

  // set next reneighbor
  force_reneighbor = 1;

  // create ParScale instance
  int pascal;
  if (me < nprocs) pascal = 1;
  else pascal = MPI_UNDEFINED;

  MPI_Comm comm_pascal;
  MPI_Comm_split(MPI_COMM_WORLD,pascal,0,&comm_pascal);

  // Open ParScale input script and create ParScale Object
  if(0 == me)
    fprintf(screen, "\n...creating ParScale object... \n");
  if(pascal == 1) pasc_ = new PASCAL_NS::ParScale(0, NULL, comm_pascal,lmp);
  pascal_setup_ = true;

  char *runDirectory = new char[128];
  int   nStrLength = 0;
  if (me == 0)
  {
      sprintf(runDirectory,"%s", "pascal"); //TODO: add option to specify via parsing
      nStrLength = strlen(runDirectory) + 1;
  }
  MPI_Bcast(&nStrLength,1,MPI_INT,0,MPI_COMM_WORLD);
  if (nStrLength > 0)
  {
        MPI_Bcast(runDirectory,nStrLength,MPI_CHAR,0,MPI_COMM_WORLD);
        if (pascal == 1) pasc_->set_dir(runDirectory);
  }

  char *pascalFile = new char[256];
  nStrLength = 0;
  if (me == 0)
  {
      sprintf(pascalFile, "./%s/%s", runDirectory, "in.pascal");  //TODO: add option to specify via parsing
      nStrLength = strlen(pascalFile) + 1;
  }
  MPI_Bcast(&nStrLength,1,MPI_INT,0,MPI_COMM_WORLD);
  if (nStrLength > 0)
  {
        MPI_Bcast(pascalFile,nStrLength,MPI_CHAR,0,MPI_COMM_WORLD);
        if (pascal == 1) pasc_->set_input(pascalFile);
  }

  pasc_->input();

  //Finalized initialization
  delete [] pascalFile;
  delete [] runDirectory;

  ts_create_ = update->ntimestep;
  fprintf(screen, "...ParScale object initialized! \n\n");
}

/* ----------------------------------------------------------------------
   free all memory for ParScale
------------------------------------------------------------------------- */

FixParScaleCouple::~FixParScaleCouple()
{
  delete pasc_;
  delete dc_;
}

/* ---------------------------------------------------------------------- */

void FixParScaleCouple::post_create()
{
  if(dc_)
        dc_->post_create();
  else
    error->all(FLERR,"internal error");


  if(verbose_) fprintf(screen, "ParScale::post_create()!\n");
  // register fixes for quantities to be saved to disk
  // see fix_property_atom.cpp for meaning of fixargs
}

/* ---------------------------------------------------------------------- */

void FixParScaleCouple::updatePtrs()
{
    //TODO
}

/* ---------------------------------------------------------------------- */

int FixParScaleCouple::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixParScaleCouple::init()
{
    // TODO: warn if something is wrong
    dc_->init();

    if(0 == atom->map_style)
      error->fix_error(FLERR,this,"requires an 'atom_modify map' command to allocate an atom map");

    //TODO
//    error->all(FLERR,"TODO: add more properties here to be pushed/pulled");
//    error->all(FLERR,"TODO: separate framework and model; put this in derived class");

    //  values to be transfered to OF
    dc_->add_push_property("id","scalar-atom");
    dc_->add_push_property("radius","scalar-atom");

    //Set points and activate exchange event recording
    comm->exchangeEvents            = true;
    exchangeEventsLocalId           = &(comm->exchangeEventsLocalId);
    exchangeEventsReceivingProcess  = &(comm->exchangeEventsReceivingProcess);

}

/* ----------------------------------------------------------------------
   make some setup calls to ParScale if necessary
   set-up is
------------------------------------------------------------------------- */

void FixParScaleCouple::setup(int vflag)
{

}

/* ----------------------------------------------------------------------
   Get the Global-to-Local index mapping
------------------------------------------------------------------------- */

int* FixParScaleCouple::get_liggghts_map(int &length)
{
    int size_map = atom->get_map_size();
    length       = size_map-1; //since LIGGGHTS' GLOBAL indexing starts with 1!

    memory->destroy(map_copy);
    memory->create<int>(map_copy, length,"FixParScaleCouple:map");
    memcpy(map_copy,&(atom->get_map_array()[1]),length*sizeof(int)); //offset by one in LIGGGHTS!

    if(verbose_)
    {
        printf("[%d/%d]:FixParScaleCouple::get_liggghts_map. nlocal: %d, nghost: %d \n",
               comm->me, comm->nprocs, atom->nlocal, atom->nghost);
        printf("atom->map_style: %d, size_map: %d, length of map_copy: %d \n",
                atom->map_style, size_map, length);

        printf("[%d/%d]:original_map[-1]: %d. \n",
               comm->me, comm->nprocs,atom->get_map_array()[0]);

        for(int iGlobal=0; iGlobal < length; iGlobal++)
        {
            printf("[%d/%d]:original_map[%d]: %d, map_copy[%d]: %d, particle r: %g\n",
                    comm->me, comm->nprocs, 
                    iGlobal, atom->get_map_array()[iGlobal+1], //offset by one in LIGGGHTS
                    iGlobal, map_copy[iGlobal],
                    atom->radius[map_copy[iGlobal]]);
        }
    }

    //Check map
    int validMapEntry = 0;
    for(int iGlobal=0; iGlobal < length; iGlobal++)
        if(map_copy[iGlobal]>-1)
            validMapEntry += 1;

    if(verbose_)
        printf("[%d/%d]: CouplingModelLiggghts: validMapEntry %d, atom->nlocal: %d atom->nghost: %d \n",
              comm->me, comm->nprocs, 
              validMapEntry,atom->nlocal,atom->nghost);

    return map_copy;
}

/* ----------------------------------------------------------------------
   detect neigh list build / exchange and trigger coupling
   schedule next coupling
------------------------------------------------------------------------- */

void FixParScaleCouple::pre_exchange()
{



}

/* ----------------------------------------------------------------------
   call ParScale to catch up with LIGGGHTS
------------------------------------------------------------------------- */
void FixParScaleCouple::end_of_step()
{
    int ts = update->ntimestep;
    time_ += update->dt;

    couple_this_step_ = false;
    
    // only execute if pushing or pulling is desired
    // do not execute on step of creation
    if( ts%couple_every_ ) return;
    couple_this_step_ = true;
    
//    if(screen && comm->me == 0)
//        fprintf(screen,"ParScale Coupling established at step %d\n",ts);
    if(logfile && comm->me == 0)
        fprintf(logfile,"ParScale Coupling established at step %d\n",ts);

    // assemble command and run in ParScale
    // init upon first use of ParScale, but not afterwards

    char commandstr[200];
    sprintf(commandstr,"control run %g init %s",
            time_,
            pascal_setup_?"yes":"no"
//            prePostRun_?"yes":"no"    //TODO: implement silent run in ParScale
            );
    pasc_->runCommand(commandstr);

    // reset flags and time counter

    pascal_setup_ = false;
    time_ = 0.;


#if 0
   int currAtom=1;
   fprintf(screen, "currAtom: %d xcm %g %g %g,vcm %g %g %g ,omega %g %g %g, torque  %g %g %g, fcm  %g %g %g\n",
            currAtom,
             xcm[currAtom][0], xcm[currAtom][1], xcm[currAtom][2],
             vcm[currAtom][0], vcm[currAtom][1], vcm[currAtom][2],
             omega[currAtom][0], omega[currAtom][1], omega[currAtom][2],
             torque[currAtom][0],torque[currAtom][1],torque[currAtom][2],
             fcm[currAtom][0],fcm[currAtom][1],fcm[currAtom][2]);
#endif
}

//////////////////////////////////////////////////////////////
void* LAMMPS_NS::FixParScaleCouple::find_pull_property(const char *name, const char *type, int &len1, int &len2)
{
    return dc_->find_pull_property(name,type,len1,len2);
}

void* LAMMPS_NS::FixParScaleCouple::find_push_property(const char *name, const char *type, int &len1, int &len2)
{
    return dc_->find_push_property(name,type,len1,len2);
}
