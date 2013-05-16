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

/* ----------------------------------------------------------------------
   Contributing authors for original version: Leo Silbert (SNL), Gary Grest (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_gran_hertz_history_simple.h"
#include "atom.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "error.h"
#include "fix_property_global.h"
#include "mech_param_gran.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairGranHertzHistorySimple::PairGranHertzHistorySimple(LAMMPS *lmp) : PairGranHookeHistorySimple(lmp)
{}

/* ---------------------------------------------------------------------- */

PairGranHertzHistorySimple::~PairGranHertzHistorySimple()
{}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairGranHertzHistorySimple::settings(int narg, char **arg) 
{
    PairGranHookeHistory::settings(narg,arg);

    // set defaults
    damp_massflag = 1;

    // parse args

    bool hasargs = true;
    while(iarg_ < narg && hasargs)
    {
        // no args to parse here
    }
}

/* ----------------------------------------------------------------------
   init specific to this granular substyle
------------------------------------------------------------------------- */

void PairGranHertzHistorySimple::init_granular()
{
  PairGranHookeHistorySimple::init_granular();

  // error checks on coarsegraining
  if(force->cg_active())
    error->cg(FLERR,force->pair_style);
}

/* ----------------------------------------------------------------------
 return appropriate params
------------------------------------------------------------------------- */

inline void PairGranHertzHistorySimple::deriveContactModelParams(int &ip, int &jp,double &meff,double &deltan, double &kn, double &kt, double &gamman, double &gammat, double &xmu, double &rmu,double &vnnr)
{
    int itype = atom->type[ip];
    int jtype = atom->type[jp];
    double rj = atom->radius[jp];
    double ri = atom->radius[ip];
    double polyhertz = sqrt(ri*rj/(ri+rj)*deltan);
    kn = polyhertz*k_n[itype][jtype];
    kt = polyhertz*k_t[itype][jtype];
    gamman = polyhertz*meff*gamma_n[itype][jtype];
    gammat = polyhertz*meff*gamma_t[itype][jtype];

    xmu=coeffFrict[itype][jtype];
    if(rollingflag)rmu=coeffRollFrict[itype][jtype];
    if (dampflag == 0) gammat = 0.0;

    // convert Kn and Kt from pressure units to force/distance^2
    kn /= force->nktv2p;
    kt /= force->nktv2p;
    return;
}
