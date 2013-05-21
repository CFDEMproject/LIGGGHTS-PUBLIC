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
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_gran_hooke_simple.h"
#include "atom.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "error.h"
#include "fix_property_global.h"
#include "mech_param_gran.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairGranHookeSimple::PairGranHookeSimple(LAMMPS *lmp) : PairGranHooke(lmp)
{
    k_n = k_t = gamma_n = gamma_t = NULL;
}

/* ---------------------------------------------------------------------- */

PairGranHookeSimple::~PairGranHookeSimple()
{
    memory->destroy(k_n);
    memory->destroy(k_t);
    memory->destroy(gamma_n);
    memory->destroy(gamma_t);

    // do not destroy coeffFrict, coeffRollFrict, cohEnergyDens
    // since are destroyed in ~PairGranHookeHistory()
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairGranHookeSimple::settings(int narg, char **arg) 
{
    PairGranHooke::settings(narg,arg);

    // set defaults
    damp_massflag = 1;

    // parse args

    bool hasargs = true;
    while(iarg_ < narg && hasargs)
    {
        hasargs = false;
        if (strcmp(arg[iarg_],"absolute_damping") == 0) {
            if (narg < iarg_+2) error->all(FLERR,"Pair gran: not enough arguments for 'absolute_damping'");
            iarg_++;
            if(strcmp(arg[iarg_],"on") == 0)
                damp_massflag = 0;
            else if(strcmp(arg[iarg_],"off") == 0)
                damp_massflag = 1;
            else
                error->all(FLERR,"Illegal pair_style gran command, expecting 'on' or 'off' after keyword 'absolute_damping'");
            iarg_++;
            hasargs = true;
        }
    }
}

/* ----------------------------------------------------------------------
   init specific to this granular substyle
------------------------------------------------------------------------- */

void PairGranHookeSimple::init_granular()
{
  int max_type = mpg->max_type();
  allocate_properties(max_type);

  //Get pointer to the fixes that have the material properties

  k_n1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("kn","property/global","peratomtypepair",max_type,max_type,force->pair_style));
  k_t1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("kt","property/global","peratomtypepair",max_type,max_type,force->pair_style));

  // can be either absolute damping value or relative (which is multiplied by effectivemass afterwards)
  if(damp_massflag == 1)
  {
    gamma_n1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("gamman","property/global","peratomtypepair",max_type,max_type,force->pair_style));
    gamma_t1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("gammat","property/global","peratomtypepair",max_type,max_type,force->pair_style));
  }
  else if(damp_massflag == 0)
  {
    gamma_n1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("gamman_abs","property/global","peratomtypepair",max_type,max_type,force->pair_style));
    gamma_t1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("gammat_abs","property/global","peratomtypepair",max_type,max_type,force->pair_style));
  }

  coeffFrict1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("coefficientFriction","property/global","peratomtypepair",max_type,max_type,force->pair_style));
  if(rollingflag)
    coeffRollFrict1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("coefficientRollingFriction","property/global","peratomtypepair",max_type,max_type,force->pair_style));

  if(cohesionflag)
    cohEnergyDens1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("cohesionEnergyDensity","property/global","peratomtypepair",max_type,max_type,force->pair_style));

  //pre-calculate parameters for possible contact material combinations
  for(int i=1;i< max_type+1; i++)
  {
      for(int j=1;j<max_type+1;j++)
      {
          k_n[i][j] = force->cg()*k_n1->compute_array(i-1,j-1);
          k_t[i][j] = k_t1->compute_array(i-1,j-1);

          // decide on coarse graining for damping depending on formulation

          if(damp_massflag == 0)
            gamma_n[i][j] = force->cg()*force->cg()*gamma_n1->compute_array(i-1,j-1);
          else if(damp_massflag == 1)
            gamma_n[i][j] = (1./force->cg())*gamma_n1->compute_array(i-1,j-1);
          gamma_t[i][j] = gamma_t1->compute_array(i-1,j-1);

          coeffFrict[i][j] = coeffFrict1->compute_array(i-1,j-1);
          if(rollingflag) coeffRollFrict[i][j] = coeffRollFrict1->compute_array(i-1,j-1);

          if(cohesionflag) cohEnergyDens[i][j] = cohEnergyDens1->compute_array(i-1,j-1);
      }
  }

  // error checks on coarsegraining
  if((rollingflag || cohesionflag) && force->cg_active())
    error->cg(FLERR,"Granular model with rolling friction and / or cohesion");
}

/* ----------------------------------------------------------------------
  allocate per-type and per-type pair properties
------------------------------------------------------------------------- */

void PairGranHookeSimple::allocate_properties(int size)
{
    memory->destroy(k_n);
    memory->destroy(k_t);
    memory->destroy(gamma_n);
    memory->destroy(gamma_t);

    memory->destroy(coeffFrict);
    memory->destroy(coeffRollFrict);

    memory->destroy(cohEnergyDens);

    memory->create(k_n,size+1,size+1,"kn");
    memory->create(k_t,size+1,size+1,"kt");
    memory->create(gamma_n,size+1,size+1,"gamman");
    memory->create(gamma_t,size+1,size+1,"gammat");

    memory->create(coeffFrict,size+1,size+1,"coeffFrict");
    memory->create(coeffRollFrict,size+1,size+1,"coeffRollFrict");

    memory->create(cohEnergyDens,size+1,size+1,"cohEnergyDens");
}

/* ----------------------------------------------------------------------
 return appropriate params
------------------------------------------------------------------------- */

inline void PairGranHookeSimple::deriveContactModelParams(int &ip, int &jp,double &meff,double &deltan, double &kn, double &kt, double &gamman, double &gammat, double &xmu, double &rmu,double &vnnr)
{
    int itype = atom->type[ip];
    int jtype = atom->type[jp];

    kn = k_n[itype][jtype];
    kt = k_t[itype][jtype];

    if(damp_massflag)
    {
        gamman = meff*gamma_n[itype][jtype];
        gammat = meff*gamma_t[itype][jtype];
    }
    else
    {
        gamman = gamma_n[itype][jtype];
        gammat = gamma_t[itype][jtype];
    }

    xmu=coeffFrict[itype][jtype];
    if(rollingflag)rmu=coeffRollFrict[itype][jtype];
    if (dampflag == 0) gammat = 0.0;

    // convert Kn and Kt from pressure units to force/distance^2
    kn /= force->nktv2p;
    kt /= force->nktv2p;
    return;
}
