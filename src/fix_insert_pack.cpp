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
    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Richard Berger (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2015 JKU Linz
------------------------------------------------------------------------- */

#include <cmath>
#include <stdlib.h>
#include <string.h>
#include "fix_insert_pack.h"
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
#include "mpi_liggghts.h"
#include "particleToInsert.h"
#include "fix_multisphere.h"
#include "math_extra_liggghts.h"

#define SEED_OFFSET 12

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixInsertPack::FixInsertPack(LAMMPS *lmp, int narg, char **arg) :
  FixInsert(lmp, narg, arg)
{
  // set defaults first, then parse args
  init_defaults();

  bool hasargs = true;
  while(iarg < narg && hasargs)
  {
    hasargs = false;
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      int iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->fix_error(FLERR,this,"region ID does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      ins_region = domain->regions[iregion];
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"volumefraction_region") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      volumefraction_region = atof(arg[iarg+1]);
      if(volumefraction_region < 0. || volumefraction_region > 1.)
        error->fix_error(FLERR,this,"Invalid volumefraction");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"particles_in_region") == 0) {
      if (iarg+2 > narg)
        error->fix_error(FLERR,this,"");
      ntotal_region = atoi(arg[iarg+1]);
      if(ntotal_region <= 0)
        error->fix_error(FLERR,this,"'ntotal_region' > 0 required");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"mass_in_region") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      masstotal_region = atof(arg[iarg+1]);
      if(masstotal_region <= 0.)
        error->fix_error(FLERR,this,"'masstotal_region' > 0 required");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"ntry_mc") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      ntry_mc = atoi(arg[iarg+1]);
      if(ntry_mc < 1000)
        error->fix_error(FLERR,this,"ntry_mc must be > 1000");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"warn_region") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      if(strcmp(arg[iarg+1],"yes") == 0)
        warn_region = true;
      else if(strcmp(arg[iarg+1],"no") == 0)
        warn_region = false;
      else
        error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'warn_region'");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"check_dist_from_subdomain_border") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      if(strcmp(arg[iarg+1],"yes") == 0)
        check_dist_from_subdomain_border_ = true;
      else if(strcmp(arg[iarg+1],"no") == 0)
        check_dist_from_subdomain_border_ = false;
      else
        error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'check_dist_from_subdomain_border'");
      iarg += 2;
      hasargs = true;
    } else if(strcmp(style,"insert/pack") == 0)
        error->fix_error(FLERR,this,"unknown keyword");
  }

  // no fixed total number of particles inserted by this fix exists
  if(strcmp(style,"insert/pack") == 0)
    ninsert_exists = 0;
}

/* ---------------------------------------------------------------------- */

FixInsertPack::~FixInsertPack()
{
    if(idregion) delete []idregion;
}

/* ---------------------------------------------------------------------- */

void FixInsertPack::init_defaults()
{
      ins_region = NULL;
      idregion = 0;
      ntry_mc = 100000;

      volumefraction_region = 0.0;
      ntotal_region = 0;
      masstotal_region = 0.0;

      region_volume = region_volume_local = 0.;

      insertion_ratio = 0.;

      check_dist_from_subdomain_border_ = true;

      warn_region = true;
}

/* ---------------------------------------------------------------------- */

void FixInsertPack::init()
{
    FixInsert::init();

    if (ins_region)
    {
        int iregion = domain->find_region(idregion);
        if (iregion == -1)
            error->fix_error(FLERR,this,"regions used by this command must not be deleted");
        ins_region = domain->regions[iregion];
    }
}

/* ----------------------------------------------------------------------
   perform error checks
------------------------------------------------------------------------- */

void FixInsertPack::calc_insertion_properties()
{
    // error check on region
    if(!ins_region)
        error->fix_error(FLERR,this,"must define an insertion region");
    ins_region->reset_random(seed + SEED_OFFSET);

    calc_region_volume_local();
    //if(comm->me == 0)
    //  printf("Total volume of the insertion region: %e\n",region_volume);
    if(region_volume <= 0. || region_volume_local < 0. || (region_volume_local - region_volume)/region_volume > 1e-3 )
        error->one(FLERR,"Fix insert: Region volume calculation with MC failed");

    if(ins_region->dynamic_check())
        error->fix_error(FLERR,this,"dynamic regions are not allowed");

    // error check on insert_every
    if(insert_every < 0)
        error->fix_error(FLERR,this,"must define 'insert_every'");

    // error checks to disallow args from FixInsert
    if(ninsert > 0 || massinsert > 0.)
        error->fix_error(FLERR,this,"specifying 'nparticles' or 'mass' not allowed");
    if(nflowrate > 0. || massflowrate > 0.)
        error->fix_error(FLERR,this,"specifying 'nflowrate' or 'massflowrate' not allowed");

    // error check if exactly one target is specified
    int n_defined = 0;
    if(volumefraction_region > 0.) n_defined++;
    if(ntotal_region > 0) n_defined++;
    if(masstotal_region > 0.) n_defined++;

    if(n_defined != 1)
        error->fix_error(FLERR,this,"must define exactly one keyword out of 'volumefraction_region', 'particles_in_region', and 'mass_in_region'");

}

/* ----------------------------------------------------------------------
   calculate volume of region on my subbox
   has to be called at initialization and before every insertion in case
   box is changing
------------------------------------------------------------------------- */

void FixInsertPack::calc_region_volume_local()
{
    ins_region->volume_mc(ntry_mc,all_in_flag==0?false:true,fix_distribution->max_r_bound(),
                          region_volume,region_volume_local);
    
}

/* ----------------------------------------------------------------------
   number of particles to insert this timestep
   depends on number of particles in region already
------------------------------------------------------------------------- */

int FixInsertPack::calc_ninsert_this()
{
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double *rmass = atom->rmass;
  double *radius = atom->radius;

  int ninsert_this = 0;

  // check if region extends outside simulation box
  // if so, throw error if boundary setting is "f f f"

  if(warn_region && ins_region->bbox_extends_outside_box())
  {
      for(int idim = 0; idim < 3; idim++)
        for(int iface = 0; iface < 2; iface++)
            if(domain->boundary[idim][iface] == 1)
                error->fix_error(FLERR,this,"Insertion region extends outside simulation box and a fixed boundary is used."
                            "Please use non-fixed boundaries in this case only");
  }

  // get number of particles, masss and occupied volume in insertion region
  // use all particles, not only those in the fix group

  int np_region = 0;
  double vol_region = 0., mass_region = 0.;
  double _4Pi3 = 4.*M_PI/3.;
  for(int i = 0; i < nlocal; i++)
  {
      
      if(fix_multisphere && fix_multisphere->belongs_to(i) >= 0) continue;
      if(ins_region->match(x[i][0],x[i][1],x[i][2]))
        
      {
          np_region++;
          vol_region += _4Pi3*radius[i]*radius[i]*radius[i];
          mass_region += rmass[i];
      }
  }

  int nbody;
  double x_bound_body[3], mass_body, density_body;
  if(multisphere)
  {
      nbody = multisphere->n_body();

      for(int ibody = 0; ibody < nbody; ibody++)
      {
          multisphere->x_bound(x_bound_body,ibody);
          //r_bound_body = multisphere->r_bound(ibody); // DEAD CODE? Side Effects?
          multisphere->r_bound(ibody); // DEAD CODE? Side Effects?
          if(ins_region->match(x_bound_body[0],x_bound_body[1],x_bound_body[2]))
            
          {
              np_region++;
              mass_body = multisphere->mass(ibody);
              density_body = multisphere->density(ibody);
              vol_region += mass_body/density_body;
              mass_region += mass_body;
          }
      }
  }

  // calculate and return number of particles that is missing

  if(volumefraction_region > 0.)
  {
      MPI_Sum_Scalar(vol_region,world);
      ninsert_this = static_cast<int>((volumefraction_region*region_volume - vol_region) / fix_distribution->vol_expect() + random->uniform());
      insertion_ratio = vol_region / (volumefraction_region*region_volume);
      
  }
  else if(ntotal_region > 0)
  {
      MPI_Sum_Scalar(np_region,world);
      ninsert_this = ntotal_region - np_region;
      insertion_ratio = static_cast<double>(np_region) / static_cast<double>(ntotal_region);
      
  }
  else if(masstotal_region > 0.)
  {
      MPI_Sum_Scalar(mass_region,world);
      ninsert_this = static_cast<int>((masstotal_region - mass_region) / fix_distribution->mass_expect() + random->uniform());
      insertion_ratio = mass_region / masstotal_region;
  }
  else error->one(FLERR,"Internal error in FixInsertPack::calc_ninsert_this()");
  //if(comm->me == 0)
  //  printf("Target number of particles to be inserted: %d\n", ninsert_this);

  // can be < 0 due to overflow, round-off etc
  if(ninsert_this < -200000)
    error->fix_error(FLERR,this,"overflow in particle number calculation: inserting too many particles in one step");
  if(ninsert_this < 0) ninsert_this = 0;

  if(insertion_ratio < 0.) insertion_ratio = 0.;
  if(insertion_ratio > 1.) insertion_ratio = 1.;

  return ninsert_this;
}

/* ---------------------------------------------------------------------- */

double FixInsertPack::insertion_fraction()
{
    // have to re-calculate region_volume_local in case simulation box is changing
    if(domain->box_change)
        calc_region_volume_local();

    return region_volume_local/region_volume;
}

/* ---------------------------------------------------------------------- */

inline int FixInsertPack::is_nearby(int i)
{
    double pos[3], rad, cut;

    vectorCopy3D(atom->x[i],pos);
    rad = atom->radius[i];

    // choose right distance depending on all_in_flag

    if(all_in_flag) cut = maxrad;
    else cut = rad + maxrad;

    if(ins_region->match_expandby_cut(pos,cut)) return 1;
    return 0;
}

/* ---------------------------------------------------------------------- */

BoundingBox FixInsertPack::getBoundingBox() {
  BoundingBox bb(ins_region->extent_xlo, ins_region->extent_xhi,
                 ins_region->extent_ylo, ins_region->extent_yhi,
                 ins_region->extent_zlo, ins_region->extent_zhi);

  double extend = 2*maxrad /*cut*/ + this->extend_cut_ghost(); 
  bb.shrinkToSubbox(domain->sublo, domain->subhi);
  bb.extendByDelta(extend);
  return bb;
}

/* ----------------------------------------------------------------------
   calc # of maximum tries
   propertional to total desired # of particles to insert on this
   subdomain to ensure insertion "does not give up too early" if a low
   remaining # of insertions is to be performed
------------------------------------------------------------------------- */

int FixInsertPack::calc_maxtry(int ninsert_this_local)
{
    if(insertion_ratio >= 1.) return ninsert_this_local * maxattempt;
    else return static_cast<int>( static_cast<double>(ninsert_this_local*maxattempt) / (1.-insertion_ratio));
}

/* ----------------------------------------------------------------------
   generate random positions within insertion volume
   perform overlap check via xnear if requested
   returns # bodies and # spheres that could actually be inserted
------------------------------------------------------------------------- */

void FixInsertPack::x_v_omega(int ninsert_this_local,int &ninserted_this_local, int &ninserted_spheres_this_local, double &mass_inserted_this_local)
{
    ninserted_this_local = ninserted_spheres_this_local = 0;
    mass_inserted_this_local = 0.;

    double pos[3];
    ParticleToInsert *pti;

    int ntry = 0;
    int maxtry = calc_maxtry(ninsert_this_local);

    double v_toInsert[3];
    vectorZeroize3D(v_toInsert);

    // no overlap check
    if(!check_ol_flag)
    {
        for(int itotal = 0; itotal < ninsert_this_local; itotal++)
        {
            pti = fix_distribution->pti_list[ninserted_this_local];
            double rbound = pti->r_bound_ins;

            if(screen && print_stats_during_flag && (ninsert_this_local >= 10) && (0 == itotal % (ninsert_this_local/10)))
                fprintf(screen,"insertion: proc %d at %d %%\n",comm->me,10*itotal/(ninsert_this_local/10));

            if(all_in_flag) ins_region->generate_random_shrinkby_cut(pos,rbound,true);
            else ins_region->generate_random(pos,true);

            // randomize vel, omega, quat here
            vectorCopy3D(v_insert,v_toInsert);
            // could ramdonize vel, omega, quat here
            generate_random_velocity(v_toInsert);

            if(quat_random_)
                MathExtraLiggghts::random_unit_quat(random,quat_insert);

            if(pos[0] == 0. && pos[1] == 0. && pos[2] == 0.)
                error->one(FLERR,"FixInsertPack::x_v_omega() illegal position");
            ninserted_spheres_this_local += pti->set_x_v_omega(pos,v_toInsert,omega_insert,quat_insert);
            mass_inserted_this_local += pti->mass_ins;
            ninserted_this_local++;

        }
    }
    // overlap check
    // account for maxattempt
    // pti checks against xnear and adds self contributions
    else
    {
        
        while(ntry < maxtry && ninserted_this_local < ninsert_this_local)
        {
            
            pti = fix_distribution->pti_list[ninserted_this_local];
            double rbound = pti->r_bound_ins;

            if(screen && print_stats_during_flag && (ninsert_this_local >= 10) && (0 == ninserted_this_local % (ninsert_this_local/10)) )
                fprintf(screen,"insertion: proc %d at %d %%\n",comm->me,10*ninserted_this_local/(ninsert_this_local/10));

            int nins = 0;
            while(nins == 0 && ntry < maxtry)
            {
                do
                {
                    
                    if(all_in_flag) ins_region->generate_random_shrinkby_cut(pos,rbound,true);
                    else ins_region->generate_random(pos,true);
                    ntry++;

                }
                
                while((check_dist_from_subdomain_border_) && (ntry < maxtry && domain->dist_subbox_borders(pos) < rbound));

                if(ntry == maxtry) break;

                // randomize vel, omega, quat here
                vectorCopy3D(v_insert,v_toInsert);

                // could ramdonize vel, omega, quat here
                generate_random_velocity(v_toInsert);

                if(quat_random_)
                    MathExtraLiggghts::random_unit_quat(random,quat_insert);

                nins = pti->check_near_set_x_v_omega(pos,v_toInsert,omega_insert,quat_insert,neighList);

            }

            if(nins > 0)
            {
                ninserted_spheres_this_local += nins;
                mass_inserted_this_local += pti->mass_ins;
                ninserted_this_local++;
            }
        }
    }
    
}

/* ---------------------------------------------------------------------- */

void FixInsertPack::restart(char *buf)
{
    FixInsert::restart(buf);

    ins_region->reset_random(seed + SEED_OFFSET);
}
