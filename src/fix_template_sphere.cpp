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
    (if not contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <climits>
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
#include "mpi_liggghts.h"
#include "fix_region_variable.h"
#include "group.h"

using namespace LAMMPS_NS;
using namespace LMP_PROBABILITY_NS;
using namespace FixConst;

#define LMP_DEBUGMODE_SPHERE false

/* ---------------------------------------------------------------------- */

FixTemplateSphere::FixTemplateSphere(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg),
    iarg(0),
    reg(NULL),
    reg_var(NULL),
    random_insertion(NULL),
    random_mc(NULL),
    seed_insertion(0),
    seed_mc(0),
    seed_orig(0),
    atom_type(0),
    pdf_radius(NULL),
    pdf_density(NULL),
    volume_expect(0.0),
    mass_expect(0.0),
    vol_limit(0.0),
    relative(false)
{
  if (domain->dimension != 3)
    error->fix_error(FLERR,this,"this fix is for 3D simulations only");

  restart_global = 1;

  // random number generator, same for all procs
  if (narg < 4) error->fix_error(FLERR,this,"not enough arguments");
  random_insertion = new RanPark(lmp, arg[3], true);
  random_mc = new RanPark(lmp, arg[3], false);
  seed_insertion = random_insertion->getSeed();
  seed_mc = random_mc->getSeed();
  seed_orig = seed_insertion;

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

  relative = false;

  //parse further args
  bool hasargs = true;
  while (iarg < narg && hasargs)
  {
    hasargs = false;
    if (strcmp(arg[iarg],"atom_type") == 0)
    {
        if (iarg+2 > narg)
            error->fix_error(FLERR,this,"not enough arguments");
        atom_type=atoi(arg[iarg+1]);
        if (atom_type < 1)
            error->fix_error(FLERR,this,"invalid atom type (must be >=1)");
        hasargs = true;
        iarg += 2;
    }
    else if (strcmp(arg[iarg],"region") == 0)
    {
        if (iarg+2 > narg)
              error->fix_error(FLERR,this,"not enough arguments");
        int ireg = domain->find_region(arg[iarg+1]);
        if (ireg < 0)
            error->fix_error(FLERR,this,"illegal region");
        reg = domain->regions[ireg];
        hasargs = true;
        iarg += 2;
    }
    else if (strcmp(arg[iarg],"relative") == 0)
    {
        if (iarg+2 > narg)
            error->fix_error(FLERR,this,"not enough arguments for 'relative'");
        if(0 == strcmp(arg[iarg+1],"yes"))
            relative = true;
        else if(0 == strcmp(arg[iarg+1],"no"))
            relative = false;
        else
            error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'relative'");
        hasargs = true;
        iarg += 2;
    }
    else if (strcmp(arg[iarg],"region_variable") == 0)
    {
        if (iarg+2 > narg)
            error->fix_error(FLERR,this,"not enough arguments");
        int ifix = modify->find_fix(arg[iarg+1]);
        if (ifix < 0)
            error->fix_error(FLERR,this,"illegal region/variable fix");
        reg_var = static_cast<FixRegionVariable*>(modify->fix[ifix]);
        hasargs = true;
        iarg += 2;
    }
    else if (strcmp(arg[iarg],"volume_limit") == 0)
    {
        if (iarg+2 > narg)
            error->fix_error(FLERR,this,"not enough arguments for 'volume_limit'");
        vol_limit = atof(arg[iarg+1]);
        if(vol_limit <= 0)
            error->fix_error(FLERR,this,"volume_limit > 0 required");
        hasargs = true;
        iarg += 2;
    }
    else if (strcmp(arg[iarg],"radius") == 0)
    {
        hasargs = true;
        if(strcmp(this->style,"particletemplate/sphere"))
            error->fix_error(FLERR,this,"keyword radius only valid for particletemplate/sphere");
        if (iarg+3 > narg)
            error->all(FLERR,"Illegal fix particletemplate/sphere command, not enough arguments");
        pdf_radius = new PDF(error);
        if (strcmp(arg[iarg+1],"constant") == 0)
        {
            double value = atof(arg[iarg+2])*force->cg(atom_type);
            if( value <= 0.)
              error->all(FLERR,"Illegal fix particletemplate/sphere command, radius must be >= 0");
            pdf_radius->set_params<RANDOM_CONSTANT>(value);
            iarg += 3;
        }
        else
            error->fix_error(FLERR,this,"invalid radius random style");
    }
    else if (strcmp(arg[iarg],"density") == 0)
    {
        hasargs = true;
        if (iarg+3 > narg)
            error->fix_error(FLERR,this,"not enough arguments");
        pdf_density = new PDF(error);
        if (strcmp(arg[iarg+1],"constant") == 0)
        {
            double value = atof(arg[iarg+2]);
            if (value <= 0.)
                error->fix_error(FLERR,this,"density must be >= 0");
            pdf_density->set_params<RANDOM_CONSTANT>(value);
            iarg += 3;
        }
        else
            error->fix_error(FLERR,this,"invalid density random style");
    }
    else if(strcmp(arg[iarg],"additional_group") == 0)
    {
        if (iarg+2 > narg)
            error->fix_error(FLERR,this,"not enough arguments for 'additional_group'");
        iarg++;
        const int iAddGroup = group->find(arg[iarg++]);
        if (iAddGroup == -1) error->all(FLERR,"Could not find additional fix group ID");
        groupbit = groupbit | group->bitmask[iAddGroup];
        
        hasargs = true;
    }
    else if(strcmp(style,"particletemplate/sphere") == 0)
        error->fix_error(FLERR,this,"unrecognized keyword");
  }

  if(pdf_density == NULL)
      error->fix_error(FLERR,this,"have to define 'density'");

  // end here for derived classes
  if(strcmp(this->style,"particletemplate/sphere"))
      return;

  if(pdf_radius == NULL)
      error->fix_error(FLERR,this,"have to define 'radius'");

  // set mass and volume expectancy
  volume_expect = cubic_expectancy(pdf_radius)*4.*M_PI/3.;
  mass_expect = expectancy(pdf_density) * volume_expect;
}

/* ---------------------------------------------------------------------- */

FixTemplateSphere::~FixTemplateSphere()
{
    delete random_insertion;
    delete random_mc;

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
    double radius = rand(pdf_radius,random_insertion);
    pti->radius_ins[0] = pti->r_bound_ins = radius;

    // randomize density
    pti->density_ins = rand(pdf_density,random_insertion);

    // calculate volume and mass
    pti->volume_ins = radius * radius * radius * 4.*M_PI/3.;
    pti->mass_ins = pti->density_ins*pti->volume_ins;

    // init insertion position
    vectorZeroize3D(pti->x_ins[0]);

    pti->groupbit = groupbit;

}

/* ----------------------------------------------------------------------*/

void FixTemplateSphere::init_ptilist(int n_random_max, const bool enforce_single, FixPropertyAtom * const fix_release)
{
    if(pti_list) error->one(FLERR,"invalid FixTemplateSphere::init_list()");
    n_pti_max = n_random_max;
    pti_list = (ParticleToInsert**) memory->smalloc(n_pti_max*sizeof(ParticleToInsert*),"pti_list");
    for(int i = 0; i < n_pti_max; i++)
       pti_list[i] = new ParticleToInsert(lmp, 1, fix_release);
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

void FixTemplateSphere::randomize_ptilist(int n_random,int distribution_groupbit,int distorder)
{
    for(int i = 0; i < n_random; i++)
    {
        
        pti_list[i]->atom_type = atom_type;

        // randomize radius
        double radius = rand(pdf_radius,random_insertion);

        pti_list[i]->radius_ins[0] = pti_list[i]->r_bound_ins = radius;

        // randomize density
        pti_list[i]->density_ins = rand(pdf_density,random_insertion);

        // calculate volume and mass
        pti_list[i]->volume_ins = radius * radius * radius * 4.*M_PI/3.;
        pti_list[i]->mass_ins = pti_list[i]->density_ins*pti_list[i]->volume_ins;

        // init insertion position
        vectorZeroize3D(pti_list[i]->x_ins[0]);
        vectorZeroize3D(pti_list[i]->v_ins);
        vectorZeroize3D(pti_list[i]->omega_ins);

        pti_list[i]->groupbit = groupbit | distribution_groupbit; 

        pti_list[i]->distorder = distorder;
    }
    
}

/* ----------------------------------------------------------------------*/

void FixTemplateSphere::direct_set_ptlist(const int i, const void * const data, const int distribution_groupbit, const int distorder)
{
    const PARTICLE_PACKING::Sphere * const sphere = static_cast<const PARTICLE_PACKING::Sphere * const>(data);
    pti_list[i]->atom_type = atom_type;
    const double radius = sphere->get_radius();
    pti_list[i]->radius_ins[0] = radius;
    pti_list[i]->density_ins = sphere->get_density();
    pti_list[i]->volume_ins = radius*radius*radius*4.1887902047863909;
    pti_list[i]->mass_ins = pti_list[i]->density_ins*pti_list[i]->volume_ins;
    pti_list[i]->id_ins = sphere->get_id();

    // set fix_property_atom
    if (pti_list[i]->fix_property || pti_list[i]->fix_property_value)
        error->one(FLERR, "Ensure that set_property is not used in fix insert");
    if (sphere->n_fix_properties() > 0)
    {
        const int n = sphere->n_fix_properties();
        pti_list[i]->n_fix_property = n;
        pti_list[i]->fix_property = new FixPropertyAtom*[n];
        pti_list[i]->fix_property_value = new double*[n];
        for (int j = 0; j < n; j++)
        {
            pti_list[i]->fix_property[j] = sphere->get_fix_property(j);
            const int m = sphere->fix_property_nentries(j);
            pti_list[i]->fix_property_value[j] = new double[m];
            for (int k = 0; k < m; k++)
                pti_list[i]->fix_property_value[j][k] = sphere->fix_property_value(j, k);
        }
    }

    // init insertion position
    vectorZeroize3D(pti_list[i]->x_ins[0]);
    vectorZeroize3D(pti_list[i]->v_ins);
    vectorZeroize3D(pti_list[i]->omega_ins);

    pti_list[i]->groupbit = groupbit | distribution_groupbit; 

    pti_list[i]->distorder = distorder;
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
  int nprocs = comm->nprocs;
  double *list = new double[2+nprocs];

  int *state_insertion_all_0 = 0;
  int state_insertion_me = random_insertion->state();
  int size_state_insertion_all_0 = MPI_Gather0_Vector<int>(&state_insertion_me,1,state_insertion_all_0,world);

  if (comm->me == 0) {
    if(size_state_insertion_all_0 != nprocs)
        error->one(FLERR,"internal error");

    list[n++] = static_cast<double>(nprocs);
    for(int iproc = 0; iproc < nprocs; iproc++)
      list[n++] = static_cast<double>(state_insertion_all_0[iproc]);
    list[n++] = static_cast<double>(random_mc->state());

    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }

  delete []list;
  if(state_insertion_all_0) delete []state_insertion_all_0;
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixTemplateSphere::restart(char *buf)
{
  double *list = (double *) buf;

  int nprocs = comm->nprocs;
  int me = comm->me;
  int nprocs_old = static_cast<double>(list[0]);

  seed_mc = static_cast<int> (list[nprocs_old+1]);

  if(nprocs <= nprocs_old)
  {
    seed_insertion = static_cast<int> (list[me+1]);
  }
  else
  {
    if(me < nprocs_old)
      seed_insertion = static_cast<int> (list[me+1]);
    else
      seed_insertion = static_cast<int> (list[(me%nprocs_old)+1]);
  }

  random_insertion->reset(seed_insertion);
}

/* ----------------------------------------------------------------------
   generate hash to identify this template
------------------------------------------------------------------------- */

unsigned int FixTemplateSphere::generate_hash()
{
    unsigned int hash = 0;
    unsigned int start = seed_orig*420001; // it's magic
    add_hash_value(atom_type, start, hash);
    add_hash_value(pdf_radius->rand_style(), start, hash);
    add_hash_value(expectancy(pdf_radius), start, hash);
    add_hash_value(cubic_expectancy(pdf_radius), start, hash);
    add_hash_value(pdf_density->rand_style(), start, hash);
    add_hash_value(expectancy(pdf_density), start, hash);
    add_hash_value(cubic_expectancy(pdf_density), start, hash);
    return hash;
}

void FixTemplateSphere::add_hash_value(const int value, unsigned int &start, unsigned int &hash)
{
    if (value >= 0)
        hash = hash*start + (unsigned int)value;
    else
        hash = hash*start + (unsigned int)(-value) + INT_MAX;
    start = start*seed_orig;
}

void FixTemplateSphere::add_hash_value(double value, unsigned int &start, unsigned int &hash)
{
    int ivalue = 0;
    if (value < 0)
        value *= -1.0;
    if (value > 1e-50)
    {
        while (value > 1e6)
            value *= 1e-6;

        while (value < 1e0)
            value *= 1e6;
        int high = (int) value;
        double remainder = (value - (double)high)*1e6;
        int low = (int) remainder;
        ivalue = high + low;
    }
    add_hash_value(ivalue, start, hash);
}
