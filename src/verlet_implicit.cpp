/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   This file was modified with respect to the release in LAMMPS
   Modifications are Copyright 2009-2012 JKU Linz
                     Copyright 2012-     DCS Computing GmbH, Linz

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include "string.h"
#include "verlet_implicit.h"
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

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

VerletImplicit::VerletImplicit(LAMMPS *lmp, int narg, char **arg) :
  Verlet(lmp, narg, arg) {}

/* ----------------------------------------------------------------------
   run for N steps iteratively
------------------------------------------------------------------------- */

void VerletImplicit::run(int n)
{
  bigint ntimestep;
  int nflag,sortflag;

  int n_post_integrate = modify->n_post_integrate;
  int n_pre_exchange = modify->n_pre_exchange;
  int n_pre_neighbor = modify->n_pre_neighbor;
  int n_pre_force = modify->n_pre_force;
  int n_post_force = modify->n_post_force;
  int n_end_of_step = modify->n_end_of_step;

  if (atom->sortfreq > 0) sortflag = 1;
  else sortflag = 0;

  for (int i = 0; i < n; i++) {

    ntimestep = ++update->ntimestep;
    
    // neigh list build if required

    nflag = neighbor->decide();

    if (nflag == 1) {
          
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

    do
    {
        ev_set(ntimestep);

        // initial time integration

        modify->initial_integrate(vflag);
        
        if (n_post_integrate) modify->post_integrate();

        // regular communication here

        timer->stamp();
        comm->forward_comm();
        timer->stamp(TIME_COMM);

        // force computations
        
        force_clear();
        if (n_pre_force) modify->pre_force(vflag);

        timer->stamp();

        if (force->pair) {
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
        
        modify->final_integrate();
    }
    while(modify->iterate_implicitly());

    if (n_end_of_step) modify->end_of_step();

    // all output
    
    if (ntimestep == output->next) {
      timer->stamp();
      output->write(ntimestep);
      timer->stamp(TIME_OUTPUT);
    }
    
  }
}
