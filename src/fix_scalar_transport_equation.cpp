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
#include "stdlib.h"
#include "string.h"
#include "fix_scalar_transport_equation.h"
#include "atom.h"
#include "domain.h"
#include "group.h"
#include "force.h"
#include "update.h"
#include "error.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair.h"
#include "math_extra.h"
#include "fix_property_global.h"
#include "fix_property_atom.h"
#include "respa.h"
#include "mech_param_gran.h"
#include "pair_gran.h"
#include "mpi_liggghts.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define SMALL 1e-8

/* ---------------------------------------------------------------------- */

FixScalarTransportEquation::FixScalarTransportEquation(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  int iarg = 3;

  if (narg < 15) error->all(FLERR,"Illegal fix transportequation/scalar command, not enough arguments");

  if(strcmp(arg[iarg++],"equation_id")) error->all(FLERR,"Fix transportequation/scalar: expecting keyword 'equation_id'");
  equation_id = new char[strlen(arg[iarg])+1];
  strcpy(equation_id,arg[iarg++]);

  if(strcmp(arg[iarg++],"quantity")) error->all(FLERR,"Fix transportequation/scalar: expecting keyword 'quantity'");
  quantity_name = new char[strlen(arg[iarg])+1];
  strcpy(quantity_name,arg[iarg++]);

  if(strcmp(arg[iarg++],"default_value")) error->all(FLERR,"Fix transportequation/scalar: expecting keyword 'default_value'");
  quantity_0 = atof(arg[iarg++]);

  if(strcmp(arg[iarg++],"flux_quantity")) error->all(FLERR,"Fix transportequation/scalar: expecting keyword 'flux_quantity'");
  flux_name = new char[strlen(arg[iarg])+1];
  strcpy(flux_name,arg[iarg++]);

  if(strcmp(arg[iarg++],"source_quantity")) error->all(FLERR,"Fix transportequation/scalar: expecting keyword 'source_quantity'");
  source_name = new char[strlen(arg[iarg])+1];
  strcpy(source_name,arg[iarg++]);

  capacity_name = NULL;
  capacity_flag = 0;

  if(strcmp(arg[iarg++],"capacity_quantity")) error->all(FLERR,"Fix transportequation/scalar: expecting keyword 'capacity_quantity'");
  if(strcmp(arg[iarg],"none"))
  {
      capacity_flag = 1;
      capacity_name = new char[strlen(arg[iarg])+1];
      strcpy(capacity_name,arg[iarg++]);
  }

  fix_quantity = fix_flux = fix_source = NULL; fix_capacity = NULL;
  capacity = NULL;

  peratom_flag = 1;              
  size_peratom_cols = 0;         
  peratom_freq = 1;

  scalar_flag = 1; 
  global_freq = 1; 
}

/* ---------------------------------------------------------------------- */

FixScalarTransportEquation::~FixScalarTransportEquation()
{
    delete []quantity_name;
    delete []flux_name;
    delete []source_name;
    delete []capacity_name;
    delete []equation_id;

    if(capacity) delete []capacity;
    
}

void FixScalarTransportEquation::pre_delete(bool unfixflag)
{
    //unregister property/atom fixes
    if (fix_quantity) modify->delete_fix(quantity_name);
    if (fix_flux) modify->delete_fix(flux_name);
    if (fix_source) modify->delete_fix(source_name);
}

/* ---------------------------------------------------------------------- */

int FixScalarTransportEquation::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= INITIAL_INTEGRATE;
  mask |= PRE_FORCE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixScalarTransportEquation::updatePtrs()
{
  quantity = fix_quantity->vector_atom;
  flux = fix_flux->vector_atom;
  source = fix_source->vector_atom;

  vector_atom = quantity; 
}

/* ---------------------------------------------------------------------- */

void FixScalarTransportEquation::post_create()
{
  char **fixarg;
  fixarg=new char*[9];
  for (int kk=0;kk<9;kk++) fixarg[kk] = new char[30];

  if (fix_quantity==NULL) {
  //register Temp as property/atom
    strcpy(fixarg[0],quantity_name);
    fixarg[1]="all";
    fixarg[2]="property/atom";
    strcpy(fixarg[3],quantity_name);
    fixarg[4]="scalar"; 
    fixarg[5]="yes";    
    fixarg[6]="yes";    
    fixarg[7]="no";    
    sprintf(fixarg[8],"%e",quantity_0);
    modify->add_fix(9,fixarg);
    fix_quantity=static_cast<FixPropertyAtom*>(modify->find_fix_property(quantity_name,"property/atom","scalar",0,0,style));
  }

  if (fix_flux==NULL){
    //register heatFlux as property/atom
    strcpy(fixarg[0],flux_name);
    fixarg[1]="all";
    fixarg[2]="property/atom";
    strcpy(fixarg[3],flux_name);
    fixarg[4]="scalar"; 
    fixarg[5]="yes";    
    fixarg[6]="no";    
    fixarg[7]="yes";    
    fixarg[8]="0.";     
    modify->add_fix(9,fixarg);
    fix_flux=static_cast<FixPropertyAtom*>(modify->find_fix_property(flux_name,"property/atom","scalar",0,0,style));
  }

  if (fix_source==NULL){
    //register heatSource as property/atom
    strcpy(fixarg[0],source_name);
    fixarg[1]="all";
    fixarg[2]="property/atom";
    strcpy(fixarg[3],source_name);
    fixarg[4]="scalar"; 
    fixarg[5]="yes";    
    fixarg[6]="yes";    
    fixarg[7]="no";    
    fixarg[8]="0.";     
    modify->add_fix(9,fixarg);
    fix_source=static_cast<FixPropertyAtom*>(modify->find_fix_property(source_name,"property/atom","scalar",0,0,style));
  }

  delete []fixarg;

  updatePtrs();
}

/* ---------------------------------------------------------------------- */

void FixScalarTransportEquation::init()
{
  if (!atom->rmass_flag) error->all(FLERR,"Please use an atom style that defines per-particle mass for fix transportequation/scalar");

  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  if(capacity_flag)
  {
      if(!force->pair_match("gran", 0)) error->all(FLERR,"Please use a granular pair style for fix transportequation/scalar with capacityflag");
      PairGran* pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));

      int max_type = pair_gran->mpg->max_type();

      if(capacity) delete []capacity;
      capacity = new double[max_type+1];

      fix_capacity = static_cast<FixPropertyGlobal*>(modify->find_fix_property(capacity_name,"property/global","peratomtype",max_type,0,style));

      //pre-calculate parameters for possible contact material combinations
      for(int i=1;i< max_type+1; i++)
          for(int j=1;j<max_type+1;j++)
              capacity[i] = fix_capacity->compute_vector(i-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixScalarTransportEquation::initial_integrate(int vflag)
{
  
  updatePtrs();

  //reset heat flux
  //sources are not reset
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
  {
       if (mask[i] & groupbit)
           flux[i]=0.;
  }

  fix_quantity->do_forward_comm();
}

/* ---------------------------------------------------------------------- */

void FixScalarTransportEquation::pre_force(int vflag)
{
    
    if(neighbor->ago == 0)
        fix_quantity->do_forward_comm();
}

/* ---------------------------------------------------------------------- */

void FixScalarTransportEquation::final_integrate()
{
    double dt = update->dt;
    int nlocal = atom->nlocal;
    double *rmass = atom->rmass;
    int *type = atom->type;
    int *mask = atom->mask;
    double capacity;

    updatePtrs();

    fix_source->do_forward_comm();

    if(capacity_flag)
    {
        for (int i = 0; i < nlocal; i++)
        {
           if (mask[i] & groupbit){
              capacity = fix_capacity->compute_vector(type[i]-1);
              if(fabs(capacity) > SMALL) quantity[i] += (flux[i] + source[i]) * dt / (rmass[i]*capacity);
           }
        }
    }
    else
    {
        for (int i = 0; i < nlocal; i++)
        {
           if (mask[i] & groupbit){
              quantity[i] += (flux[i] + source[i]) * dt;
           }
        }
    }
}

/* ---------------------------------------------------------------------- */

void FixScalarTransportEquation::initial_integrate_respa(int vflag, int ilevel, int flag)
{
  // outermost level - update
  // all other levels - nothing

  if (ilevel == nlevels_respa-1) initial_integrate(vflag);
}

/* ---------------------------------------------------------------------- */

double FixScalarTransportEquation::compute_scalar()
{
    double *rmass = atom->rmass;
    int *type = atom->type;
    int nlocal = atom->nlocal;
    double capacity;

    updatePtrs();
    
    double quantity_sum = 0.;

    if(capacity_flag)
    {
        for (int i = 0; i < nlocal; i++)
        {
           capacity = fix_capacity->compute_vector(type[i]-1);
           quantity_sum += capacity * rmass[i] * quantity[i];
        }
    }
    else
    {
        for (int i = 0; i < nlocal; i++)
        {
           quantity_sum += quantity[i];
        }
    }

   MPI_Sum_Scalar(quantity_sum,world);
    return quantity_sum;
}

/* ---------------------------------------------------------------------- */

bool FixScalarTransportEquation::match_equation_id(const char* id)
{
    if(strcmp(id,equation_id)) return false;
    return true;
}
