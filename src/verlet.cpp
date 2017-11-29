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
    This file is from LAMMPS, but has been modified. Copyright for
    modification:

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz

    Copyright of original file:
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

#include <string.h>
#include <stdio.h>
#include <time.h>
#include "verlet.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "atom.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "output.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
#include "signal_handling.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Verlet::Verlet(LAMMPS *lmp, int narg, char **arg) :
  Integrate(lmp, narg, arg) {}

/* ----------------------------------------------------------------------
   initialization before run
------------------------------------------------------------------------- */

void Verlet::init()
{
  Integrate::init();

  // warn if no fixes

  if (modify->nfix == 0 && comm->me == 0)
    error->warning(FLERR,"No fixes defined, atoms won't move");

  // virial_style:
  // 1 if computed explicitly by pair->compute via sum over pair interactions
  // 2 if computed implicitly by pair->virial_fdotr_compute via sum over ghosts

  if (force->newton_pair) virial_style = 2;
  else virial_style = 1;

  // setup lists of computes for global and per-atom PE and pressure

  ev_setup();

  // detect if fix omp is present for clearing force arrays

  int ifix = modify->find_fix("package_omp");
  if (ifix >= 0) external_force_clear = 1;

  // set flags for what arrays to clear in force_clear()
  // need to clear additionals arrays if they exist

  torqueflag = 0;
  if (atom->torque_flag) torqueflag = 1;
  erforceflag = 0;
  if (atom->erforce_flag) erforceflag = 1;
  e_flag = 0;
  if (atom->e_flag) e_flag = 1;
  rho_flag = 0;
  if (atom->rho_flag) rho_flag = 1;

  // orthogonal vs triclinic simulation box

  triclinic = domain->triclinic;
}

/* ----------------------------------------------------------------------
   setup before run
------------------------------------------------------------------------- */

void Verlet::setup()
{
  time_t curtime;
  time(&curtime);

  if (comm->me == 0 && screen) fprintf(screen,"Setting up run at %s\n",ctime(&curtime));
  if (comm->me == 0 && logfile) fprintf(logfile,"Setting up run at %s\n",ctime(&curtime));

  update->setupflag = 1;

  // setup domain, communication and neighboring
  // acquire ghosts
  // build neighbor lists
  
  atom->setup();
  modify->setup_pre_exchange();
  
  if (triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  
  comm->setup();
  if (neighbor->style) neighbor->setup_bins();
  
  comm->exchange();
  if (atom->sortfreq > 0) atom->sort();
  
  comm->borders();
  if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
  
  modify->setup_pre_neighbor(); 
  
  neighbor->build();
  neighbor->ncalls = 0;

  // compute all forces

  ev_set(update->ntimestep);
  force_clear();
  
  modify->setup_pre_force(vflag);

  if (pair_compute_flag) force->pair->compute(eflag,vflag);
  else if (force->pair) force->pair->compute_dummy(eflag,vflag);

  if (atom->molecular) {
    if (force->bond) force->bond->compute(eflag,vflag);
    if (force->angle) force->angle->compute(eflag,vflag);
    if (force->dihedral) force->dihedral->compute(eflag,vflag);
    if (force->improper) force->improper->compute(eflag,vflag);
  }

  if (force->kspace) {
    force->kspace->setup();
    if (kspace_compute_flag) force->kspace->compute(eflag,vflag);
    else force->kspace->compute_dummy(eflag,vflag);
  }

  if (force->newton) comm->reverse_comm();

  modify->setup(vflag);
  
  output->setup();
  
  update->setupflag = 0;
}

/* ----------------------------------------------------------------------
   setup without output
   flag = 0 = just force calculation
   flag = 1 = reneighbor and force calculation
------------------------------------------------------------------------- */

void Verlet::setup_minimal(int flag)
{
  update->setupflag = 1;

  // setup domain, communication and neighboring
  // acquire ghosts
  // build neighbor lists

  if (flag) {
    modify->setup_pre_exchange();
    if (triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    domain->reset_box();
    comm->setup();
    if (neighbor->style) neighbor->setup_bins();
    comm->exchange();
    comm->borders();
    if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
    domain->image_check();
    domain->box_too_small_check();
    modify->setup_pre_neighbor();
    neighbor->build();
    neighbor->ncalls = 0;
  }

  // compute all forces

  ev_set(update->ntimestep);
  force_clear();
  modify->setup_pre_force(vflag);

  if (pair_compute_flag) force->pair->compute(eflag,vflag);
  else if (force->pair) force->pair->compute_dummy(eflag,vflag);

  if (atom->molecular) {
    if (force->bond) force->bond->compute(eflag,vflag);
    if (force->angle) force->angle->compute(eflag,vflag);
    if (force->dihedral) force->dihedral->compute(eflag,vflag);
    if (force->improper) force->improper->compute(eflag,vflag);
  }

  if (force->kspace) {
    force->kspace->setup();
    if (kspace_compute_flag) force->kspace->compute(eflag,vflag);
    else force->kspace->compute_dummy(eflag,vflag);
  }

  if (force->newton) comm->reverse_comm();

  modify->setup(vflag);
  update->setupflag = 0;
}

/* ----------------------------------------------------------------------
   run for N steps
------------------------------------------------------------------------- */

void Verlet::run(int n)
{
  bigint ntimestep;
  int nflag,sortflag;

  const int n_pre_initial_integrate = modify->n_pre_initial_integrate;
  const int n_post_integrate = modify->n_post_integrate;
  const int n_pre_exchange = modify->n_pre_exchange;
  const int n_pre_neighbor = modify->n_pre_neighbor;
  const int n_pre_force = modify->n_pre_force;
  const int n_post_force = modify->n_post_force;
  const int n_pre_final_integrate = modify->n_pre_final_integrate;
  const int n_end_of_step = modify->n_end_of_step;

  if (atom->sortfreq > 0) sortflag = 1;
  else sortflag = 0;

  for (int i = 0; i < n; i++) {

    ntimestep = ++update->ntimestep;
    
    ev_set(ntimestep);

    // pre-integration step

    if (n_pre_initial_integrate) modify->pre_initial_integrate();

    // initial time integration

    modify->initial_integrate(vflag);
    
    if (n_post_integrate) modify->post_integrate();

    // regular communication vs neighbor list rebuild

    nflag = neighbor->decide();

    if (nflag == 0) {
      timer->stamp();
      comm->forward_comm();
      timer->stamp(TIME_COMM);
    } else {
      
      if (n_pre_exchange) modify->pre_exchange();
      
      if (triclinic) domain->x2lamda(atom->nlocal);
      domain->pbc();
      if (domain->box_change) {
        domain->reset_box();
        
        comm->setup();
        if (neighbor->style) neighbor->setup_bins();
      }
      timer->stamp();
      
      comm->exchange();
      
      if (sortflag && ntimestep >= atom->nextsort) atom->sort();
      comm->borders();
      if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
      timer->stamp(TIME_COMM);
      
      if (n_pre_neighbor) modify->pre_neighbor();
      
      neighbor->build();
      timer->stamp(TIME_NEIGHBOR);
    }

    // force computations
    // important for pair to come before bonded contributions
    // since some bonded potentials tally pairwise energy/virial
    // and Pair:ev_tally() needs to be called before any tallying

    force_clear();
    if (n_pre_force) modify->pre_force(vflag);

    timer->stamp();

    if (pair_compute_flag) {
      force->pair->compute(eflag,vflag);
      timer->stamp(TIME_PAIR);
    }

    if (atom->molecular) {
      if (force->bond) force->bond->compute(eflag,vflag);
      if (force->angle) force->angle->compute(eflag,vflag);
      if (force->dihedral) force->dihedral->compute(eflag,vflag);
      if (force->improper) force->improper->compute(eflag,vflag);
      timer->stamp(TIME_BOND);
    }

    if (kspace_compute_flag) {
      force->kspace->compute(eflag,vflag);
      timer->stamp(TIME_KSPACE);
    }

    // reverse communication of forces
    
    if (force->newton) {
      comm->reverse_comm();
      timer->stamp(TIME_COMM);
    }

    // force modifications, final time integration, diagnostics
    
    if (n_post_force) modify->post_force(vflag);
    
    if (n_pre_final_integrate) modify->pre_final_integrate();
    
    modify->final_integrate();
    
    if (n_end_of_step) modify->end_of_step();

    // all output
    
    if (ntimestep == output->next) {
      timer->stamp();
      output->write(ntimestep);
      timer->stamp(TIME_OUTPUT);
    }
    
    if (SignalHandler::request_quit && !SignalHandler::request_write_restart)
        break;
  }
}

/* ---------------------------------------------------------------------- */

void Verlet::cleanup()
{
  modify->post_run();
  domain->box_too_small_check();
  update->update_time();
}

/* ----------------------------------------------------------------------
   clear force on own & ghost atoms
   clear other arrays as needed
------------------------------------------------------------------------- */

void Verlet::force_clear()
{
  int i;

  if (external_force_clear) return;

  // clear force on all particles
  // if either newton flag is set, also include ghosts
  // when using threads always clear all forces.

  if (neighbor->includegroup == 0) {
    int nall;
    
    nall = atom->nlocal + atom->nghost;
    //if (force->newton) nall = atom->nlocal + atom->nghost;
    //else nall = atom->nlocal;

    size_t nbytes = sizeof(double) * nall;

    if (nbytes) {
      memset(&(atom->f[0][0]),0,3*nbytes);
      if (torqueflag)  memset(&(atom->torque[0][0]),0,3*nbytes);
      if (erforceflag) memset(&(atom->erforce[0]),  0,  nbytes);
      if (e_flag)      memset(&(atom->de[0]),       0,  nbytes);
      if (rho_flag)    memset(&(atom->drho[0]),     0,  nbytes);
    }

  // neighbor includegroup flag is set
  // clear force only on initial nfirst particles
  // if either newton flag is set, also include ghosts

  } else {
    int nall = atom->nfirst;

    double **f = atom->f;
    for (i = 0; i < nall; i++) {
      f[i][0] = 0.0;
      f[i][1] = 0.0;
      f[i][2] = 0.0;
    }

    if (torqueflag) {
      double **torque = atom->torque;
      for (i = 0; i < nall; i++) {
        torque[i][0] = 0.0;
        torque[i][1] = 0.0;
        torque[i][2] = 0.0;
      }
    }

    if (erforceflag) {
      double *erforce = atom->erforce;
      for (i = 0; i < nall; i++) erforce[i] = 0.0;
    }

    if (e_flag) {
      double *de = atom->de;
      for (i = 0; i < nall; i++) de[i] = 0.0;
    }

    if (rho_flag) {
      double *drho = atom->drho;
      for (i = 0; i < nall; i++) drho[i] = 0.0;
    }

    if (force->newton) {
      nall = atom->nlocal + atom->nghost;

      for (i = atom->nlocal; i < nall; i++) {
        f[i][0] = 0.0;
        f[i][1] = 0.0;
        f[i][2] = 0.0;
      }

      if (torqueflag) {
        double **torque = atom->torque;
        for (i = atom->nlocal; i < nall; i++) {
          torque[i][0] = 0.0;
          torque[i][1] = 0.0;
          torque[i][2] = 0.0;
        }
      }

      if (erforceflag) {
        double *erforce = atom->erforce;
        for (i = atom->nlocal; i < nall; i++) erforce[i] = 0.0;
      }

      if (e_flag) {
        double *de = atom->de;
        for (i = 0; i < nall; i++) de[i] = 0.0;
      }

      if (rho_flag) {
        double *drho = atom->drho;
        for (i = 0; i < nall; i++) drho[i] = 0.0;
      }
    }
  }
}
