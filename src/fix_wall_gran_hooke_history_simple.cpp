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
#include "fix_wall_gran_hooke_history_simple.h"
#include "pair_gran_hooke_history_simple.h"
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

FixWallGranHookeHistorySimple::FixWallGranHookeHistorySimple(LAMMPS *lmp, int narg, char **arg) :
  FixWallGranHookeHistory(lmp, narg, arg)
{
    k_n1 = k_t1 = gamma_n1 = gamma_t1 = NULL;
}

/* ---------------------------------------------------------------------- */

void FixWallGranHookeHistorySimple::init_granular()
{
  k_n = ((PairGranHookeHistorySimple*)pairgran_)->k_n;
  k_t = ((PairGranHookeHistorySimple*)pairgran_)->k_t;
  gamma_n = ((PairGranHookeHistorySimple*)pairgran_)->gamma_n;
  gamma_t = ((PairGranHookeHistorySimple*)pairgran_)->gamma_t;
  coeffFrict = ((PairGranHookeHistorySimple*)pairgran_)->coeffFrict;
  coeffRollFrict = ((PairGranHookeHistorySimple*)pairgran_)->coeffRollFrict;

  damp_massflag = ((PairGranHookeHistorySimple*)pairgran_)->damp_massflag;

  if(viscousflag)
    error->fix_error(FLERR,this,"cannot use 'stokes' model for fix wall/gran model without coefficient of restitution");

  //need to check properties for rolling friction and cohesion energy density here
  //since these models may not be active in the pair style
  
  FixPropertyGlobal *coeffRollFrict1, *cohEnergyDens1;
  int max_type = pairgran_->mpg->max_type();
  if(rollingflag)
    coeffRollFrict1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("coefficientRollingFriction","property/global","peratomtypepair",max_type,max_type,style));
  if(cohesionflag)
    cohEnergyDens1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("cohesionEnergyDensity","property/global","peratomtypepair",max_type,max_type,style));

  //pre-calculate parameters for possible contact material combinations
  for(int i=1;i< max_type+1; i++)
  {
      for(int j=1;j<max_type+1;j++)
      {
          if(rollingflag) coeffRollFrict[i][j] = coeffRollFrict1->compute_array(i-1,j-1);
          if(cohesionflag) cohEnergyDens[i][j] = cohEnergyDens1->compute_array(i-1,j-1);
      }
  }

  // error checks on coarsegraining
  if((rollingflag || cohesionflag) && force->cg_active())
    error->cg(FLERR,"Granular model with rolling friction and / or cohesion");
}

/* ----------------------------------------------------------------------
   contact model parameters derived for hertz model 
------------------------------------------------------------------------- */

inline void FixWallGranHookeHistorySimple::deriveContactModelParams(int ip, double deltan,double meff_wall, double &kn, double &kt, double &gamman, double &gammat, double &xmu,double &rmu,double &vnnr)  
{
    int itype = atom->type[ip];

    kn = k_n[itype][atom_type_wall_];
    kt = k_t[itype][atom_type_wall_];

    if(damp_massflag)
    {
        gamman = meff_wall*gamma_n[itype][atom_type_wall_];
        gammat = meff_wall*gamma_t[itype][atom_type_wall_];
    }
    else
    {
        gamman = gamma_n[itype][atom_type_wall_];
        gammat = gamma_t[itype][atom_type_wall_];
    }

    xmu=coeffFrict[itype][atom_type_wall_];
    if(rollingflag)rmu=coeffRollFrict[itype][atom_type_wall_];

    if (dampflag == 0) gammat = 0.0;

    // convert Kn and Kt from pressure units to force/distance^2
    kn /= force->nktv2p;
    kt /= force->nktv2p;
    return;
}
