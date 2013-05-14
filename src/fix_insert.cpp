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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "domain.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"
#include "fix_multisphere.h"
#include "fix_particledistribution_discrete.h"
#include "fix_template_sphere.h"
#include "fix_insert.h"
#include "math_extra_liggghts.h"
#include "mpi_liggghts.h"
#include "vector_liggghts.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define EPSILON 0.001

#define LMP_DEBUGMODE_FIXINSERT false //(667 == update->ntimestep)//  true
#define LMP_DEBUG_OUT_FIXINSERT screen

/* ---------------------------------------------------------------------- */

FixInsert::FixInsert(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) error->fix_error(FLERR,this,"not enough arguments");

  time_depend = 1; 
  restart_global = 1;

  setup_flag = false;

  fix_distribution = NULL;
  fix_multisphere = NULL;
  multisphere = NULL;

  // required args
  iarg = 3;

  if(strcmp(arg[iarg++],"seed")) error->fix_error(FLERR,this,"expecting keyword 'seed'");
  seed = atoi(arg[iarg++]) + comm->me;
  if (seed <= 0) error->fix_error(FLERR,this,"illegal seed");

  // random number generator, seed depends on proc
  random = new RanPark(lmp,seed);

  // set defaults
  init_defaults();

  xnear = NULL;

  // parse args
  
  bool hasargs = true;
  while(iarg < narg && hasargs)
  {
    hasargs = false;
    if(strcmp(arg[iarg],"distributiontemplate") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      int ifix = modify->find_fix(arg[iarg+1]);
      if(ifix < 0 || strcmp(modify->fix[ifix]->style,"particledistribution/discrete"))
        error->fix_error(FLERR,this,"Fix insert requires you to define a valid ID for a fix of type particledistribution/discrete");
      fix_distribution = static_cast<FixParticledistributionDiscrete*>(modify->fix[ifix]);
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"maxattempt") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      maxattempt = atoi(arg[iarg+1]);
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"nparticles") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      if(strcmp(arg[iarg+1],"INF") == 0)
        ninsert_exists = 0;
      else ninsert = atof(arg[iarg+1]);
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"mass") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      if(strcmp(arg[iarg+1],"INF") == 0)
        ninsert_exists = 0;
      else massinsert = atof(arg[iarg+1]);
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"massrate") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      massflowrate = atof(arg[iarg+1]);
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"particlerate") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      nflowrate = atof(arg[iarg+1]);
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"insert_every") == 0 || strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      if(strcmp(arg[iarg+1],"once") == 0) insert_every = 0;
      else insert_every = atoi(arg[iarg+1]);
      if(insert_every < 0) error->fix_error(FLERR,this,"insert_every must be >= 0");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"start") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      first_ins_step = atoi(arg[iarg+1]);
      if(first_ins_step < update->ntimestep + 1 && !modify->fix_restart_in_progress())
        error->fix_error(FLERR,this,"'start' step can not be before current step");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"overlapcheck") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      if(strcmp(arg[iarg+1],"yes")==0) check_ol_flag = 1;
      else if(strcmp(arg[iarg+1],"no")==0) check_ol_flag = 0;
      else error->fix_error(FLERR,this,"");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"all_in") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      if(strcmp(arg[iarg+1],"yes")==0) all_in_flag = 1;
      else if(strcmp(arg[iarg+1],"no")==0) all_in_flag = 0;
      else error->fix_error(FLERR,this,"");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"random_distribute") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      if(strcmp(arg[iarg+1],"uncorrelated")==0) exact_number = 0;
      else if(strcmp(arg[iarg+1],"exact")==0) exact_number = 1;
      else error->fix_error(FLERR,this,"");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"verbose") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      if(strcmp(arg[iarg+1],"no")==0) print_stats_during_flag = 0;
      else if(strcmp(arg[iarg+1],"yes")==0) print_stats_during_flag = 1;
      else error->fix_error(FLERR,this,"");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"vel") == 0) {
      if (iarg+5 > narg) error->fix_error(FLERR,this,"not enough keyword for 'vel'");
      if (strcmp(arg[iarg+1],"constant") == 0)  {
          v_insert[0] = atof(arg[iarg+2]);
          v_insert[1] = atof(arg[iarg+3]);
          v_insert[2] = atof(arg[iarg+4]);
          iarg += 5;
      } else if (strcmp(arg[iarg+1],"uniform") == 0) {
          if (iarg+8 > narg) error->fix_error(FLERR,this,"not enough keyword for 'uniform'");
          v_randomSetting = 1; //switch 1...distribute with equal prop.
          v_insert[0] = atof(arg[iarg+2]);
          v_insert[1] = atof(arg[iarg+3]);
          v_insert[2] = atof(arg[iarg+4]);
          v_insertFluct[0] = atof(arg[iarg+5]);
          v_insertFluct[1] = atof(arg[iarg+6]);
          v_insertFluct[2] = atof(arg[iarg+7]);
          iarg += 8;
      } else if (strcmp(arg[iarg+1],"gaussian") == 0) {
          if (iarg+8 > narg) error->fix_error(FLERR,this,"not enough keyword for 'gaussian'");
          v_randomSetting = 2; //switch 2...distribute with gaussian distrib.
          v_insert[0] = atof(arg[iarg+2]);
          v_insert[1] = atof(arg[iarg+3]);
          v_insert[2] = atof(arg[iarg+4]);
          v_insertFluct[0] = atof(arg[iarg+5]);
          v_insertFluct[1] = atof(arg[iarg+6]);
          v_insertFluct[2] = atof(arg[iarg+7]);
          iarg += 8;
      } else
          error->fix_error(FLERR,this,"expecting keyword 'constant' or 'uniform' or 'gaussian' after keyword 'vel'");
      hasargs = true;
    } else if (strcmp(arg[iarg],"omega") == 0) {
      if (iarg+5 > narg) error->fix_error(FLERR,this,"");
      if (strcmp(arg[iarg+1],"constant") == 0)
      {
          omega_insert[0] = atof(arg[iarg+2]);
          omega_insert[1] = atof(arg[iarg+3]);
          omega_insert[2] = atof(arg[iarg+4]);
      } else error->fix_error(FLERR,this,"expecting keyword 'constant' after keyword 'omega'");
      iarg += 5;
      hasargs = true;
    } else if (strcmp(arg[iarg],"quat") == 0) {
      if (iarg+5 > narg) error->fix_error(FLERR,this,"");
      if (strcmp(arg[iarg+1],"constant") == 0)
      {
          quat_insert[0] = atof(arg[iarg+2]);
          quat_insert[1] = atof(arg[iarg+3]);
          quat_insert[2] = atof(arg[iarg+4]);
          quat_insert[3] = atof(arg[iarg+5]);
      } else error->fix_error(FLERR,this,"expecting keyword 'constant' after keyword 'quat'");
      iarg += 6;
      hasargs = true;
    }
    
    else if(strcmp(style,"insert") == 0) error->fix_error(FLERR,this,"unknown keyword");
  }

  // memory not allocated initially
  ninsert_this_max_local = 0;

  // check for missing or contradictory settings
  sanity_check();

  //min/max type to be inserted, need that to check if material properties defined for all materials
  type_max = fix_distribution->max_type();
  type_min = fix_distribution->min_type();

  // allgather arrays
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  recvcounts = new int[nprocs];
  displs = new int[nprocs];

  // set next reneighbor
  force_reneighbor = 1;
  next_reneighbor = first_ins_step;
  most_recent_ins_step = -1;

  vector_flag = 1;
  size_vector = 2;

  print_stats_start_flag = 1;

  // calc max insertion radius
  int ntypes = atom->ntypes;
  maxrad = 0.;
  minrad = 1000.;
  for(int i = 1; i <= ntypes; i++)
  {
     maxrad = MathExtraLiggghts::max(maxrad,max_rad(i));
     minrad = MathExtraLiggghts::min(minrad,min_rad(i));
  }
}

/* ---------------------------------------------------------------------- */

FixInsert::~FixInsert()
{
  delete random;
  delete [] recvcounts;
  delete [] displs;
}

/* ---------------------------------------------------------------------- */

void FixInsert::setup(int vflag)
{
  
  // do this only once
  if(setup_flag) return;
  else setup_flag = true;

  // calculate ninsert, insert_every, ninsert_per
  calc_insertion_properties();

  // calc last step of insertion
  if(ninsert_exists)
  {
      if(ninsert < ninsert_per)
        final_ins_step = first_ins_step;
      else
        final_ins_step = first_ins_step +
                static_cast<int>(static_cast<double>(ninsert)/ninsert_per *  static_cast<double>(insert_every));

      if(final_ins_step < 0)
        error->fix_error(FLERR,this,"Particle insertion: Overflow - need too long for particle insertion. "
                                    "Please decrease # particles to insert or increase insertion rate");
      if(ninsert < 0)
        error->fix_error(FLERR,this,"Particle insertion: Overflow - too many particles for particle insertion. "
                                    "Please decrease # particles to insert.");
  }
  else
    final_ins_step = -1;

  // print statistics
  print_stats_start();

}

/* ---------------------------------------------------------------------- */

void FixInsert::init_defaults()
{
  // default is that total # of particles to insert by this command is known
  ninsert_exists = 1;

  ninsert = ninserted = 0;
  massinsert = massinserted = 0.;
  nflowrate = massflowrate = 0.;

  insert_every = -1;
  ninsert_per = 0.;

  // 1st insertion on next timestep is default
  first_ins_step = update->ntimestep + 1;

  maxattempt = 50;

  check_ol_flag = 1;
  all_in_flag = 0;

  exact_number = 1;

  v_randomSetting = 0;
  vectorZeroize3D(v_insert);
  vectorZeroize3D(v_insertFluct);
  vectorZeroize3D(omega_insert);

  quatUnitize4D(quat_insert);

  print_stats_during_flag = 1;
}

/* ---------------------------------------------------------------------- */

void FixInsert::sanity_check()
{
    if(fix_distribution == NULL) error->fix_error(FLERR,this,"have to define a 'distributiontemplate'");
    if(vectorMag4DSquared(quat_insert) != 1.) error->fix_error(FLERR,this,"quaternion not valid");

    if(ninsert > 0 && massinsert > 0.) error->fix_error(FLERR,this,"must not define both 'nparticles' and 'mass'");
    if(nflowrate > 0. && massflowrate > 0.) error->fix_error(FLERR,this,"must not define both 'particlerate' and 'massrate'");

    if(insert_every == 0 && (massflowrate > 0. || nflowrate > 0.)) error->fix_error(FLERR,this,"must not define 'particlerate' or 'massrate' for 'insert_every' = 0");
}

/* ---------------------------------------------------------------------- */

void FixInsert::print_stats_start()
{
  if (me == 0 && print_stats_start_flag) {

    if(ninsert_exists)
    {
        if (screen)
            fprintf(screen ,"INFO: Particle insertion %s: %f particles every %d steps - particle rate %f  (mass rate %f)\n"
                            "      %d particles (mass %f) within %d steps\n",
                id,ninsert_per,insert_every,nflowrate,massflowrate,ninsert,massinsert,final_ins_step-first_ins_step);

        if (logfile)
            fprintf(logfile,"INFO: Particle insertion %s: %f particles every %d steps - particle rate %f, (mass rate %f)\n"
                            "      %d particles (mass %f) within %d steps\n",
                id,ninsert_per,insert_every,nflowrate,massflowrate,ninsert,massinsert,final_ins_step-first_ins_step);
    }
    else
    {
        if (screen)
            fprintf(screen ,"INFO: Particle insertion %s: %f particles every %d steps - particle rate %f  (mass rate %f)\n",
                id,ninsert_per,insert_every,nflowrate,massflowrate);

        if (logfile)
            fprintf(logfile,"INFO: Particle insertion %s: %f particles every %d steps - particle rate %f, (mass rate %f)\n",
                id,ninsert_per,insert_every,nflowrate,massflowrate);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixInsert::print_stats_during(int ninsert_this, double mass_inserted_this)
{
  int step = update->ntimestep;

  if (me == 0 && print_stats_during_flag)
  {
    if (screen)
      fprintf(screen ,"INFO: Particle insertion %s: inserted %d particle templates (mass %f) at step %d\n - a total of %d particle templates (mass %f) inserted so far.\n",
              id,ninsert_this,mass_inserted_this,step,ninserted,massinserted);

    if (logfile)
      fprintf(logfile,"INFO: Particle insertion %s: inserted %d particle templates (mass %f) at step %d\n - a total of %d particle templates (mass %f) inserted so far.\n",
              id,ninsert_this,mass_inserted_this,step,ninserted,massinserted);
  }
}

/* ---------------------------------------------------------------------- */

int FixInsert::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixInsert::init()
{
    int ntimestep = update->ntimestep;

    if (!atom->radius_flag || !atom->rmass_flag)
        error->fix_error(FLERR,this,"Fix insert requires atom attributes radius, rmass");
    if (domain->triclinic)
        error->fix_error(FLERR,this,"Cannot use with triclinic box");
    if (domain->dimension != 3)
        error->fix_error(FLERR,this,"Can use fix insert for 3d simulations only");
    
    fix_multisphere = static_cast<FixMultisphere*>(modify->find_fix_style("multisphere", 0));
    if(!fix_multisphere) multisphere = NULL;
    else multisphere = &fix_multisphere->data();

    if(fix_multisphere && fix_multisphere->igroup != igroup)
        error->fix_error(FLERR,this,"Fix insert command and fix multisphere command are not compatible, must be same group");

    // in case of new fix insert in a restarted simulation, have to add current time-step
    if(next_reneighbor > 0 && next_reneighbor < ntimestep)
        error->fix_error(FLERR,this,"'start' step can not be before current step");
}

/* ---------------------------------------------------------------------- */

int FixInsert::min_type()
{
    return type_min;
}

/* ---------------------------------------------------------------------- */

int FixInsert::max_type()
{
    return type_max;
}

/* ---------------------------------------------------------------------- */

double FixInsert::max_rad(int type)
{
    return fix_distribution->max_rad(type);
}

/* ---------------------------------------------------------------------- */

double FixInsert::min_rad(int type)
{
    return fix_distribution->min_rad(type);
}

/* ---------------------------------------------------------------------- */

double FixInsert::max_r_bound()
{
    return fix_distribution->max_r_bound();
}

/* ---------------------------------------------------------------------- */

double FixInsert::extend_cut_ghost()
{
    
    return 2.*fix_distribution->max_r_bound();
}

/* ---------------------------------------------------------------------- */

int FixInsert::calc_ninsert_this()
{
  if(ninsert_per == 0.) error->fix_error(FLERR,this,"ninsert_per == 0.");

  // number of bodies to insert this timestep
  int ninsert_this = static_cast<int>(ninsert_per + random->uniform());
  if (ninsert_exists && ninserted + ninsert_this > ninsert) ninsert_this = ninsert - ninserted;

  return ninsert_this;
}

/* ----------------------------------------------------------------------
   perform particle insertion
------------------------------------------------------------------------- */

void FixInsert::pre_exchange()
{
  
  int ninsert_this, ninsert_this_local; // global and local # bodies to insert this time-step

  // just return if should not be called on this timestep
  
  if (next_reneighbor != update->ntimestep || most_recent_ins_step == update->ntimestep) return;
  most_recent_ins_step = update->ntimestep;

  // things to be done before inserting new particles
  pre_insert();

  // number of particles to insert this timestep
  ninsert_this = calc_ninsert_this();
  
  // limit to max number of particles that shall be inserted
  // to avoid that max # may be slightly exceeded by random processes
  // in fix_distribution->randomize_list, set exact_number to 1
  if (ninsert_exists && ninserted + ninsert_this > ninsert)
  {
      ninsert_this = ninsert - ninserted;
      exact_number = 1;
  }

  // distribute ninsert_this across processors
  ninsert_this_local = distribute_ninsert_this(ninsert_this);
  
  // re-allocate list if necessary
  
  if(ninsert_this_local > ninsert_this_max_local)
  {
      fix_distribution->random_init_list(ninsert_this_local);
      ninsert_this_max_local = ninsert_this_local;
  }

  // generate list of insertions
  // number of inserted particles can change if exact_number = 0
  
  ninsert_this_local = fix_distribution->randomize_list(ninsert_this_local,groupbit,exact_number);
  
  MPI_Sum_Scalar(ninsert_this_local,ninsert_this,world);

  if(ninsert_this == 0)
  {
      // warn if flowrate should be fulfilled
      if((nflowrate > 0. || massflowrate > 0.) && comm->me == 0)
        error->warning(FLERR,"Particle insertion: Inserting no particle - check particle insertion settings");

      // schedule next insertion
      if (insert_every && (!ninsert_exists || ninserted < ninsert))
        next_reneighbor += insert_every;

      return;
  }
  else if(ninsert_this < 0)
  {
      error->one(FLERR,"Particle insertion: Internal error");
  }

  // warn if max # insertions exceeded by random processes
  if (ninsert_exists && ninserted + ninsert_this > ninsert)
  {
      error->warning(FLERR,"INFO: Particle insertion: Number of particles to insert was slightly exceeded by random process");
  }

  // fill xnear array with particles to check overlap against
  
  // add particles in insertion volume to xnear list
  nspheres_near = 0;
  xnear = NULL;
  if(check_ol_flag)
      nspheres_near = load_xnear(ninsert_this_local);

  // insertion counters in this step
  int ninserted_this = 0, ninserted_spheres_this = 0;
  int ninserted_this_local = 0, ninserted_spheres_this_local = 0;
  double mass_inserted_this = 0.;
  double mass_inserted_this_local = 0.;

  // randomize insertion positions and set v, omega
  // also performs overlap check via xnear if requested
  // returns # bodies and # spheres that could actually be inserted
  x_v_omega(ninsert_this_local,ninserted_this_local,ninserted_spheres_this_local,mass_inserted_this_local);

  // actual particle insertion
  
  ninserted_spheres_this_local = fix_distribution->insert(ninserted_this_local);

  // warn if max # insertions exceeded by random processes
  if (ninsert_exists && ninserted + ninsert_this > ninsert)
  {
      error->warning(FLERR,"INFO: Particle insertion: Number of particles to insert was slightly exceeded by random process");
  }

  // set tag # of new particles beyond all previous atoms, reset global natoms
  // if global map exists, reset it now instead of waiting for comm
  // since deleting atoms messes up ghosts

  if (atom->tag_enable)
  {
    atom->tag_extend();
    atom->natoms += static_cast<double>(ninserted_spheres_this);
    if (atom->map_style)
    {
      atom->nghost = 0;
      atom->map_init();
      atom->map_set();
    }
  }

  // give particle distributions the chance to do some wrap-up
  
  fix_distribution->finalize_insertion();

  // give derived classes the chance to do some wrap-up
  
  finalize_insertion(ninserted_spheres_this_local);

  // tally stats
  MPI_Sum_Scalar(ninserted_this_local,ninserted_this,world);
  ninserted += ninserted_this;
  MPI_Sum_Scalar(mass_inserted_this_local,mass_inserted_this,world);
  massinserted += mass_inserted_this;
  print_stats_during(ninserted_this,mass_inserted_this);

  if(ninserted_this < ninsert_this && comm->me == 0)
      error->warning(FLERR,"Particle insertion: Less insertions than requested");

  // free local memory
  if(xnear) memory->destroy(xnear);

  // next timestep to insert
  if (insert_every && (!ninsert_exists || ninserted < ninsert)) next_reneighbor += insert_every;
  else next_reneighbor = 0;

}

/* ----------------------------------------------------------------------
   distribute insertions across processors
------------------------------------------------------------------------- */

//not standard c, VS complains about unknown function tztzt
double round(double d)
{
  return floor(d + 0.5);
}

int FixInsert::distribute_ninsert_this(int ninsert_this)
{
    int me, nprocs, ngap, ninsert_this_local, *ninsert_this_local_all;
    double fraction_local, fraction_local_all_sum, *fraction_local_all, *remainder, r, rsum;

    me = comm->me;
    nprocs = comm->nprocs;

    fraction_local = insertion_fraction();
    
    if(!exact_number)
        return static_cast<int>(fraction_local*static_cast<double>(ninsert_this) + random->uniform());

    // for exact_number==1, have to allgather to exactly match ninsert_this

    fraction_local_all = new double[nprocs];
    remainder = new double[nprocs];
    ninsert_this_local_all = new int[nprocs];

    // allgather local fractions
    MPI_Allgather(&fraction_local,1,MPI_DOUBLE,fraction_local_all,1,MPI_DOUBLE,world);

    // proc0 calculates ninsert_this_local for all processes
    
    if(me == 0)
    {
        // remove fractions < 2% / nprocs
        // have to normalize so not all portions get cancelled away for higher proc counts
        // normalize fraction_local_all so sum across processors is 1

        double lower_thresh = 0.02 / static_cast<double>(nprocs);

        fraction_local_all_sum = 0.;
        for(int iproc = 0; iproc < nprocs; iproc++)
        {
            if(fraction_local_all[iproc] < lower_thresh)
                fraction_local_all[iproc] = 0.;
            fraction_local_all_sum += fraction_local_all[iproc];
        }

        if(fraction_local_all_sum == 0.)
            error->one(FLERR,"Internal error distributing particles to processes");

        for(int iproc = 0; iproc < nprocs; iproc++)
            fraction_local_all[iproc] /= fraction_local_all_sum;

        rsum = 0.;
        for(int iproc = 0; iproc < nprocs; iproc++)
        {
            ninsert_this_local_all[iproc] = static_cast<int>(fraction_local_all[iproc]*static_cast<double>(ninsert_this));
            remainder[iproc] = fraction_local_all[iproc]*static_cast<double>(ninsert_this) - ninsert_this_local_all[iproc];
            rsum += remainder[iproc];
            
        }

        ngap = round(rsum);
        
        for(int i = 0; i < ngap; i++)
        {
            r = random->uniform() * static_cast<double>(ngap);
            int iproc = 0;
            rsum = remainder[iproc];

            while(iproc < (nprocs-1) && rsum < r)
            {
                iproc++;
                rsum += remainder[iproc];
            }
            ninsert_this_local_all[iproc]++;
        }
    }

    // Bcast the result
    MPI_Bcast(ninsert_this_local_all,nprocs, MPI_INT,0,world);
    ninsert_this_local = ninsert_this_local_all[me];

    delete []fraction_local_all;
    delete []remainder;
    delete []ninsert_this_local_all;

    return ninsert_this_local;
}

/* ----------------------------------------------------------------------
   count # of particles that could overlap
   must loop local + ghost particles
------------------------------------------------------------------------- */

int FixInsert::count_nnear()
{
    int nall = atom->nlocal + atom->nghost;
    int ncount = 0;

    for(int i = 0; i < nall; i++)
        ncount += is_nearby(i);

    return ncount;
}

/* ----------------------------------------------------------------------
   fill xnear with nearby particles
------------------------------------------------------------------------- */

int FixInsert::load_xnear(int ninsert_this_local)
{
  // count nearby spheres
  // setup for allgatherv
  int nspheres_near_local = count_nnear();

  // data size per particle: x and radius
  int n = 4*nspheres_near;

  // xnear is for my atoms + atoms to be inserted
  
  memory->create(xnear,nspheres_near_local + ninsert_this_local*fix_distribution->max_nspheres(), 4, "FixInsert::xnear");

  // load up xnear array with local and ghosts

  double **x = atom->x;
  double *radius = atom->radius;
  int nall = atom->nlocal + atom->nghost;

  int ncount = 0;
  for (int i = 0; i < nall; i++)
  {
    if (is_nearby(i))
    {
      xnear[ncount][0] = x[i][0];
      xnear[ncount][1] = x[i][1];
      xnear[ncount][2] = x[i][2];
      xnear[ncount][3] = radius[i];
      ncount++;
    }
  }

  return nspheres_near_local;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixInsert::write_restart(FILE *fp)
{
  int n = 0;
  double list[5];
  list[n++] = static_cast<double>(random->state());
  list[n++] = static_cast<double>(ninserted);
  list[n++] = static_cast<double>(first_ins_step);
  list[n++] = static_cast<double>(next_reneighbor);
  list[n++] = massinserted;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixInsert::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;
  double next_reneighbor_re;

  seed = static_cast<int> (list[n++]) + comm->me;
  ninserted = static_cast<int> (list[n++]);
  first_ins_step = static_cast<int> (list[n++]);
  next_reneighbor_re = static_cast<int> (list[n++]);
  massinserted = list[n++];

  random->reset(seed);

  // in order to be able to continue pouring with increased number of particles
  // if insert was already finished in run to be restarted
  if(next_reneighbor_re != 0 && ninserted < ninsert) next_reneighbor = next_reneighbor_re;
}

/* ----------------------------------------------------------------------
   output
------------------------------------------------------------------------- */

double FixInsert::compute_vector(int index)
{
    if(index == 0) return static_cast<double>(ninserted);
    if(index == 1) return massinserted;
}
