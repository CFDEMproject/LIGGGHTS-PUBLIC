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
#include "stdlib.h"
#include "string.h"
#include "fix_wall_gran_hertz_history_simple.h"
#include "pair_gran_hertz_history_simple.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "respa.h"
#include "memory.h"
#include "error.h"
#include "fix_rigid.h"
#include "fix_property_global.h"
#include "mech_param_gran.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define BIG 1.0e20

#define MIN(A,B) (((A) < (B)) ? (A) : (B))
#define MAX(A,B) (((A) > (B)) ? (A) : (B))

/* ---------------------------------------------------------------------- */

FixWallGranHertzHistorySimple::FixWallGranHertzHistorySimple(LAMMPS *lmp, int narg, char **arg) :
  FixWallGranHookeHistorySimple(lmp, narg, arg)
{}

/* ---------------------------------------------------------------------- */

inline void FixWallGranHertzHistorySimple::deriveContactModelParams(int ip, double deltan, double meff_wall, double &kn, double &kt, double &gamman, double &gammat, double &xmu,double &rmu,double &vnnr)  
{
    double polyhertz = sqrt(atom->radius[ip]*deltan);
    int itype = atom->type[ip];

    kn = polyhertz*k_n[itype][atom_type_wall_];
    kt = polyhertz*k_t[itype][atom_type_wall_];
    gamman = polyhertz*meff_wall*gamma_n[itype][atom_type_wall_];
    gammat = polyhertz*meff_wall*gamma_t[itype][atom_type_wall_];

    xmu=coeffFrict[itype][atom_type_wall_];
    if(rollingflag)rmu=coeffRollFrict[itype][atom_type_wall_];

    if (dampflag == 0) gammat = 0.0;

    // convert Kn and Kt from pressure units to force/distance^2
    kn /= force->nktv2p;
    kt /= force->nktv2p;
    return;
}
