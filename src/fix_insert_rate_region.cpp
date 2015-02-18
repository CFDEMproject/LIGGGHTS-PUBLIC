/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
   the producer of the LIGGGHTS® software and the CFDEM®coupling software
   See http://www.cfdem.com/terms-trademark-policy for details.

   LIGGGHTS® is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_insert_rate_region.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "region.h"
#include "domain.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"
#include "fix_particledistribution_discrete.h"
#include "fix_template_sphere.h"
#include "vector_liggghts.h"
#include "particleToInsert.h"

#define SEED_OFFSET 12

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixInsertRateRegion::FixInsertRateRegion(LAMMPS *lmp, int narg, char **arg) :
  FixInsertPack(lmp, narg, arg)
{
  bool hasargs = true;
  while(iarg < narg && hasargs)
  {
    hasargs = false;
    if (strcmp(arg[iarg],"some_arg") == 0) {
    } else if(strcmp(style,"insert/rate/region") == 0)
    error->fix_error(FLERR,this,"unknown keyword or wrong keyword order");
  }
}

/* ---------------------------------------------------------------------- */

FixInsertRateRegion::~FixInsertRateRegion()
{
}

/* ---------------------------------------------------------------------- */

void FixInsertRateRegion::calc_insertion_properties()
{
    double dt = update->dt;

    // error check on region
    if(!ins_region)
        error->fix_error(FLERR,this,"must define an insertion region");
    ins_region->reset_random(seed + SEED_OFFSET);
    ins_region->volume_mc(ntry_mc,all_in_flag==0?false:true,fix_distribution->max_r_bound(),
                          region_volume,region_volume_local);
    if(region_volume <= 0. || region_volume_local < 0. || region_volume_local > region_volume)
        error->one(FLERR,"Fix insert: Region volume calculation with MC failed");

    if(ins_region->dynamic_check())
        error->fix_error(FLERR,this,"dynamic regions are not allowed");

    // error check on insert_every
    if(insert_every < 0)
        error->fix_error(FLERR,this,"must define 'insert_every'");
    if(insert_every == 0)
        error->fix_error(FLERR,this,"'insert_every' = once not allowed");

    // some error checks
    if(nflowrate > 0. && massflowrate > 0.)
        error->fix_error(FLERR,this,"both 'nflowrate' and 'massflowrate' not allowed");
    if(ninsert > 0 && massinsert > 0.)
        error->fix_error(FLERR,this,"must not define both 'nparticles' and 'mass'");

    // ninsert - either defined defined directly or calculated
    if(ninsert == 0&& ninsert_exists)
    {
        if(massinsert > 0.) ninsert = static_cast<int>(massinsert / fix_distribution->mass_expect());
        else error->fix_error(FLERR,this,"must define either 'nparticles' or 'mass'");
    }

    // flow rate, ninsert_per
    if(nflowrate == 0.)
    {
        if(massflowrate == 0.) error->fix_error(FLERR,this,"must define either 'massrate' or 'particlerate'");
        nflowrate = massflowrate / fix_distribution->mass_expect();
    }
    else massflowrate = nflowrate * fix_distribution->mass_expect();

    ninsert_per = nflowrate*(static_cast<double>(insert_every)*dt);
    if(ninsert_exists) massinsert = static_cast<double>(ninsert) * fix_distribution->mass_expect();

}

/* ---------------------------------------------------------------------- */

int FixInsertRateRegion::calc_ninsert_this()
{
  // check if region extends outside simulation box
  // if so, throw error if boundary setting is "f f f"

  if(ins_region->bbox_extends_outside_box())
  {
      for(int idim = 0; idim < 3; idim++)
        for(int iface = 0; iface < 2; iface++)
            if(domain->boundary[idim][iface] == 1)
                error->fix_error(FLERR,this,"Insertion region extends outside simulation box and a fixed boundary is used."
                            "Please use non-fixed boundaries in this case only");
  }

  return FixInsert::calc_ninsert_this();
}

/* ----------------------------------------------------------------------
   calc # of maximum tries - directly linked to number of particles to insert
------------------------------------------------------------------------- */

int FixInsertRateRegion::calc_maxtry(int ninsert_this_local)
{
    return ninsert_this_local * maxattempt;
}
