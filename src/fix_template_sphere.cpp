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
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "fix_template_sphere.h"
#include "atom.h"
#include "atom_vec.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "random_park.h"
#include "particleToInsert.h"
#include "region.h"
#include "domain.h"
#include "force.h"
#include "comm.h"
#include "vector_liggghts.h"
#include "fix_region_variable.h"

using namespace LAMMPS_NS;
using namespace LMP_PROBABILITY_NS;
using namespace FixConst;

#define LMP_DEBUGMODE_SPHERE false

/* ---------------------------------------------------------------------- */

FixTemplateSphere::FixTemplateSphere(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (domain->dimension != 3)
    error->fix_error(FLERR,this,"this fix is for 3D simulations only");

  restart_global = 1;

  // random number generator, same for all procs
  if (narg < 4) error->fix_error(FLERR,this,"not enough arguments");
  seed = atoi(arg[3]) + comm->me;
  random = new RanPark(lmp,seed);

  iarg = 4;

  // set default values
  atom_type = 1;
  vol_limit = 1e-12;

  pdf_radius = NULL;
  pdf_density = NULL;

  pti = new ParticleToInsert(lmp);

  n_pti_max = 0;
  pti_list = NULL;

  reg = NULL;
  reg_var = NULL;

  //parse further args
  bool hasargs = true;
  while (iarg < narg && hasargs)
  {
    hasargs = false;
    if (strcmp(arg[iarg],"atom_type") == 0)
    {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments");
      atom_type=atoi(arg[iarg+1]);
      if (atom_type < 1) error->fix_error(FLERR,this,"invalid atom type (must be >=1)");
      hasargs = true;
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"region") == 0)
    {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments");
      int ireg = domain->find_region(arg[iarg+1]);
      if (ireg < 0) error->fix_error(FLERR,this,"illegal region");
      reg = domain->regions[ireg];
      hasargs = true;
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"region_variable") == 0)
    {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments");
      int ifix = modify->find_fix(arg[iarg+1]);
      if (ifix < 0) error->fix_error(FLERR,this,"illegal region/variable fix");
      reg_var = static_cast<FixRegionVariable*>(modify->fix[ifix]);
      hasargs = true;
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"volume_limit") == 0)
    {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments for 'volume_limit'");
      vol_limit = atof(arg[iarg+1]);
      if(vol_limit <= 0)
        error->fix_error(FLERR,this,"volume_limit > 0 required");
      hasargs = true;
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"radius") == 0) {
      hasargs = true;
      if(strcmp(this->style,"particletemplate/sphere"))
        error->fix_error(FLERR,this,"keyword radius only valid for particletemplate/sphere");
      if (iarg+3 > narg)
        error->all(FLERR,"Illegal fix particletemplate/sphere command, not enough arguments");
      pdf_radius = new PDF(error);
      if (strcmp(arg[iarg+1],"constant") == 0)
      {
          double value = atof(arg[iarg+2])*force->cg();
          if( value <= 0.)
            error->all(FLERR,"Illegal fix particletemplate/sphere command, radius must be >= 0");
          pdf_radius->set_params<RANDOM_CONSTANT>(value);
          iarg += 3;
      }
      else if (strcmp(arg[iarg+1],"uniform") == 0)
      {
          if(iarg+5 > narg)
              error->fix_error(FLERR,this,"not enough arguments for uniform");
          if(strcmp(arg[iarg+2],"mass") == 0)
              pdf_radius->activate_mass_shift();
          else if(strcmp(arg[iarg+2],"number"))
              error->fix_error(FLERR,this,"expecting 'mass' or 'number'");
          if(force->cg() > 1.)
              error->fix_error(FLERR,this,"cannot use distribution with coarse-graining.");
          double min = atof(arg[iarg+3]);
          double max = atof(arg[iarg+4]);
          if( min <= 0. || max <= 0. || min >= max)
            error->fix_error(FLERR,this,"illegal min or max value for radius");
          pdf_radius->set_params<RANDOM_UNIFORM>(min,max);
          iarg += 5;
      }
      else if (strcmp(arg[iarg+1],"gaussian") == 0)
      {
          if(iarg+5 > narg)
              error->fix_error(FLERR,this,"not enough arguments for gaussian");
          if(strcmp(arg[iarg+2],"mass") == 0)
              pdf_radius->activate_mass_shift();
          else if(strcmp(arg[iarg+2],"number"))
              error->fix_error(FLERR,this,"expecting 'mass' or 'number'");
          if(force->cg() > 1.)
              error->fix_error(FLERR,this,"cannot use distribution with coarse-graining.");
          double mu = atof(arg[iarg+3]);
          double sigma = atof(arg[iarg+4]);
          if( mu <= 0. ) error->fix_error(FLERR,this,"illegal mu value for radius");
          if( sigma <= 0. ) error->fix_error(FLERR,this,"illegal sigma value for radius");
          pdf_radius->set_params<RANDOM_GAUSSIAN>(mu,sigma);
          iarg += 5;
      }
      else if (strcmp(arg[iarg+1],"lognormal") == 0)
      {
          if(iarg+5 > narg)
              error->fix_error(FLERR,this,"not enough arguments for lognormal");
          if(strcmp(arg[iarg+2],"mass") == 0)
              pdf_radius->activate_mass_shift();
          else if(strcmp(arg[iarg+2],"number"))
              error->fix_error(FLERR,this,"expecting 'mass' or 'number'");
          if(force->cg() > 1.)
              error->fix_error(FLERR,this,"cannot use distribution with coarse-graining.");
          double mu = atof(arg[iarg+3]);
          double sigma = atof(arg[iarg+4]);
          if( sigma <= 0. ) error->fix_error(FLERR,this,"illegal sigma value for radius");
          pdf_radius->set_params<RANDOM_LOGNORMAL>(mu,sigma);
          iarg += 5;
      }
      else error->fix_error(FLERR,this,"invalid radius random style");
    }
    else if (strcmp(arg[iarg],"density") == 0) {
      hasargs = true;
      if (iarg+3 > narg) error->fix_error(FLERR,this,"not enough arguments");
      pdf_density = new PDF(error);
      if (strcmp(arg[iarg+1],"constant") == 0)
      {
          double value = atof(arg[iarg+2]);
          if( value <= 0.) error->fix_error(FLERR,this,"density must be >= 0");
          pdf_density->set_params<RANDOM_CONSTANT>(value);
          iarg += 3;
      }
      else if (strcmp(arg[iarg+1],"uniform") == 0)
      {
          if (iarg+4 > narg) error->fix_error(FLERR,this,"not enough arguments");
          double min = atof(arg[iarg+2]);
          double max = atof(arg[iarg+3]);
          if( min <= 0. || max <= 0. || min >= max)
            error->fix_error(FLERR,this,"illegal min or max value for density");
          pdf_density->set_params<RANDOM_UNIFORM>(min,max);
          iarg += 4;
      }
      else if (strcmp(arg[iarg+1],"gaussian") == 0)
      {
          if (iarg+4 > narg)
            error->fix_error(FLERR,this,"not enough arguments");
          double mu = atof(arg[iarg+2]);
          double sigma = atof(arg[iarg+3]);
          if( mu <= 0. )
            error->fix_error(FLERR,this,"illegal mu value for density");
          if( sigma <= 0. ) error->all(FLERR,"illegal sigma value for density");
          pdf_density->set_params<RANDOM_GAUSSIAN>(mu,sigma);
          iarg += 4;
      }
      else if (strcmp(arg[iarg+1],"lognormal") == 0)
      {
          if (iarg+4 > narg)
            error->fix_error(FLERR,this,"not enough arguments");
          double mu = atof(arg[iarg+2]);
          double sigma = atof(arg[iarg+3]);
          if( sigma <= 0. ) error->all(FLERR,"illegal sigma value for density");
          pdf_density->set_params<RANDOM_LOGNORMAL>(mu,sigma);
          iarg += 4;
      }
      else error->fix_error(FLERR,this,"invalid density random style");
      
    }
    else if(strcmp(style,"particletemplate/sphere") == 0)
        error->fix_error(FLERR,this,"unrecognized keyword");
  }

  if(pdf_density == NULL) error->fix_error(FLERR,this,"have to define 'density'");

  // end here for derived classes
  if(strcmp(this->style,"particletemplate/sphere"))return;

  if(pdf_radius == NULL) error->fix_error(FLERR,this,"have to define 'radius'");

  // set mass and volume expectancy
  volume_expect = cubic_expectancy(pdf_radius)*4.*M_PI/3.;
  mass_expect = expectancy(pdf_density) * volume_expect;
}

/* ---------------------------------------------------------------------- */

FixTemplateSphere::~FixTemplateSphere()
{
    delete random;

    delete pdf_density;
    delete pdf_radius;

    if(strcmp(style,"particletemplate/sphere") == 0)
    {
        delete pti;
        if(pti_list) delete_ptilist();
    }
}

/* ----------------------------------------------------------------------*/

int FixTemplateSphere::setmask()
{
  int mask = 0;
  return mask;
}

/* ----------------------------------------------------------------------*/

Region* FixTemplateSphere::region()
{
    if(reg_var) return reg_var->region();
    else return reg;
}

/* ----------------------------------------------------------------------*/

void FixTemplateSphere::randomize_single()
{
    
    pti->atom_type = atom_type;

    // randomize radius
    double radius = rand(pdf_radius,random);
    pti->radius_ins[0] = pti->r_bound_ins = radius;

    // randomize density
    pti->density_ins = rand(pdf_density,random);

    // calculate volume and mass
    pti->volume_ins = radius * radius * radius * 4.*M_PI/3.;
    pti->mass_ins = pti->density_ins*pti->volume_ins;

    // init insertion position
    vectorZeroize3D(pti->x_ins[0]);

    pti->groupbit = groupbit;

}

/* ----------------------------------------------------------------------*/

void FixTemplateSphere::init_ptilist(int n_random_max)
{
    if(pti_list) error->one(FLERR,"invalid FixTemplateSphere::init_list()");
    n_pti_max = n_random_max;
    pti_list = (ParticleToInsert**) memory->smalloc(n_pti_max*sizeof(ParticleToInsert*),"pti_list");
    for(int i = 0; i < n_pti_max; i++)
       pti_list[i] = new ParticleToInsert(lmp);
}

/* ----------------------------------------------------------------------*/

void FixTemplateSphere::delete_ptilist()
{
    if(n_pti_max == 0) return;

    for(int i = 0; i < n_pti_max; i++)
       delete pti_list[i];

    memory->sfree(pti_list);
    pti_list = NULL;
    n_pti_max = 0;
}

/* ----------------------------------------------------------------------*/

void FixTemplateSphere::randomize_ptilist(int n_random,int distribution_groupbit)
{
    for(int i = 0; i < n_random; i++)
    {
        
        pti_list[i]->atom_type = atom_type;

        // randomize radius
        double radius = rand(pdf_radius,random);

        pti_list[i]->radius_ins[0] = pti_list[i]->r_bound_ins = radius;

        // randomize density
        pti_list[i]->density_ins = rand(pdf_density,random);

        // calculate volume and mass
        pti_list[i]->volume_ins = radius * radius * radius * 4.*M_PI/3.;
        pti_list[i]->mass_ins = pti_list[i]->density_ins*pti_list[i]->volume_ins;

        // init insertion position
        vectorZeroize3D(pti_list[i]->x_ins[0]);
        vectorZeroize3D(pti_list[i]->v_ins);
        vectorZeroize3D(pti_list[i]->omega_ins);

        pti_list[i]->groupbit = groupbit | distribution_groupbit; 
    }
    
}

/* ----------------------------------------------------------------------*/

double FixTemplateSphere::min_rad()
{
    
    return pdf_min(pdf_radius);
}

/* ----------------------------------------------------------------------*/

double FixTemplateSphere::max_rad()
{
    
    return pdf_max(pdf_radius);
}

/* ----------------------------------------------------------------------*/

double FixTemplateSphere::max_r_bound()
{
    return pdf_max(pdf_radius);
}

/* ----------------------------------------------------------------------*/

double FixTemplateSphere::volexpect()
{
    if(volume_expect < vol_limit)
    {
        
        error->fix_error(FLERR,this,"Volume expectancy too small. Change 'volume_limit' "
        "if you are sure you know what you're doing");
    }
    return volume_expect;
}

/* ----------------------------------------------------------------------*/

double FixTemplateSphere::massexpect()
{
    return mass_expect;
}

/* ----------------------------------------------------------------------*/

int FixTemplateSphere::number_spheres()
{
    return 1;
}

/* ----------------------------------------------------------------------*/

int FixTemplateSphere::maxtype()
{
    return atom_type;
}

/* ----------------------------------------------------------------------*/

int FixTemplateSphere::mintype()
{
    return atom_type;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixTemplateSphere::write_restart(FILE *fp)
{
  int n = 0;
  double list[1];
  list[n++] = static_cast<int>(random->state());

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixTemplateSphere::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  seed = static_cast<int> (list[n++]) + comm->me;

  random->reset(seed);
}
