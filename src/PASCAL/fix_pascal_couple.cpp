/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com


   Copyright (C): 2014 DCS Computing GmbH (www.dcs-computing.com), Linz, Austria
                  2014 Graz University of Technology (ippt.tugraz.at), Graz, Austria

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   Parts of the code were developped in the frame of the NanoSim project funded
   by the European Commission through FP7 Grant agreement no. 604656.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "stdio.h"
#include "string.h"
#include "stdlib.h"
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
   Constructor, calls PaScal's constructor
------------------------------------------------------------------------- */
FixPaScalCouple::FixPaScalCouple(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  fix_shellTemperature_(0),
  fix_shellHeatFlux_(0),
  dc_(0),
  verbose_(false),
  couple_at_least_every_(0),
  couple_this_step_(false),
  pascal_setup_(0),
  pasc_(0),
  time_(0.)
{
  iarg_ = 3;

  int me, nprocs;
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

  if (narg < 5)
    error->fix_error(FLERR,this,"not enough arguments");

  if (strcmp(arg[iarg_++],"couple_at_least_every"))
    error->fix_error(FLERR,this,"expecting keyword 'couple_at_least_every'");

  couple_at_least_every_ = force->inumeric(FLERR,arg[iarg_++]);
  if (couple_at_least_every_ <= 0)
    error->fix_error(FLERR,this,"'couple_at_least_every' > 0 required");

  nevery = 1;

  //Create datacoupling
  dc_ = new CfdDatacouplingSimple(lmp,iarg_+1,narg,arg,this);

  // set next reneighbor
  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + couple_at_least_every_;

  // TODO: meed parsing here!
  // TODO: need to parse filename with pascal input

  // create PaScal instance
  int pascal;
  if (me < nprocs) pascal = 1;
  else pascal = MPI_UNDEFINED;

  MPI_Comm comm_pascal;
  MPI_Comm_split(MPI_COMM_WORLD,pascal,0,&comm_pascal);

  // Open PaScal input script and create PaScal Object
  if(0 == me)
    fprintf(screen, "\n...creating PaScal object... \n");
  if(pascal == 1) pasc_ = new PASCAL_NS::PaScal(0, NULL, comm_pascal,lmp);
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

  fprintf(screen, "...PaScal object initialized! \n\n");
}

/* ----------------------------------------------------------------------
   free all memory for PaScal
------------------------------------------------------------------------- */

FixPaScalCouple::~FixPaScalCouple()
{
  delete pasc_;
}

/* ---------------------------------------------------------------------- */

void FixPaScalCouple::post_create()
{
  if(dc_)
        dc_->post_create();
  else
    error->all(FLERR,"internal error");


  if(verbose_) fprintf(screen, "PaScal::post_create()!\n");
  // register fixes for quantities to be saved to disk
  // see fix_property_atom.cpp for meaning of fixargs
  if(!fix_shellTemperature_)
  {
        const char* fixarg[9];
        fixarg[0]="shellT";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="shellT";
        fixarg[4]="scalar";
        fixarg[5]="no";
        fixarg[6]="yes";
        fixarg[7]="no";
        fixarg[8]="0.";
	    fix_shellTemperature_ =
	                modify->add_fix_property_atom(9,
	                                              const_cast<char**>(fixarg),
	                                              style);
  }
}

/* ---------------------------------------------------------------------- */

void FixPaScalCouple::updatePtrs()
{
    //TODO
}

/* ---------------------------------------------------------------------- */

int FixPaScalCouple::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPaScalCouple::init()
{
    // TODO: warn if something is wrong

    dc_->init();

    if(0 == atom->map_style)
      error->fix_error(FLERR,this,"requires an 'atom_modify map' command to allocate an atom map");

    //TODO
//    error->all(FLERR,"TODO: add more properties here to be pushed/pulled");
//    error->all(FLERR,"TODO: separate framework and model; put this in derived class");

    //  values to be transfered to OF

    dc_->add_push_property("x","vector-atom");
    dc_->add_push_property("v","vector-atom");
    dc_->add_push_property("radius","scalar-atom");
    dc_->add_push_property("id","scalar-atom");
}

/* ----------------------------------------------------------------------
   make some setup calls to PaScal if necessary
   set-up is
------------------------------------------------------------------------- */

void FixPaScalCouple::setup(int vflag)
{

}

/* ----------------------------------------------------------------------
   detect neigh list build / exchange and trigger coupling
   schedule next coupling
------------------------------------------------------------------------- */

int* FixPaScalCouple::get_liggghts_map(int &length)
{

    int size_map = atom->get_map_size();
    length       = size_map;

    int *map_copy = NULL;
    map_copy = memory->create<int>(map_copy, size_map,"map");
    memcpy(map_copy,atom->get_map_array(),size_map);
    return map_copy;
}

/* ----------------------------------------------------------------------
   detect neigh list build / exchange and trigger coupling
   schedule next coupling
------------------------------------------------------------------------- */

void FixPaScalCouple::pre_exchange()
{
    couple_this_step_ = true;

    if (next_reneighbor != update->ntimestep)
        fprintf(screen,"'premature' LIGGGHTS-PaScal coupling because of high flow dynamics\n");

    next_reneighbor = update->ntimestep + couple_at_least_every_;
}

/* ----------------------------------------------------------------------
   call PaScal to catch up with LIGGGHTS
------------------------------------------------------------------------- */
void FixPaScalCouple::end_of_step()
{
    time_ += update->dt;

    if(!couple_this_step_)
        return;

    // assemble command and run in PaScal
    // init upon first use of PaScal, but not afterwards

    char commandstr[200];
    sprintf(commandstr,"control run %f init %s",time_,pascal_setup_?"yes":"no");
    fprintf(screen,"PaScal::runCommand()!\n");
    pasc_->runCommand(commandstr);

    // reset flags and time counter

    couple_this_step_ = false;
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
void* LAMMPS_NS::FixPaScalCouple::find_pull_property(const char *name, const char *type, int &len1, int &len2)
{ 
    return dc_->find_pull_property(name,type,len1,len2); 
}

void* LAMMPS_NS::FixPaScalCouple::find_push_property(const char *name, const char *type, int &len1, int &len2)
{ 
    return dc_->find_push_property(name,type,len1,len2); 
}
