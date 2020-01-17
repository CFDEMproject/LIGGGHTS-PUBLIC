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

    Implementation of implicit update algorithm
    Copyright 2016      TU Graz, Stefan Radl
    Copyright 2016      DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#include <cmath>
#include <stdlib.h>
#include <string.h>
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
#include "properties.h"
#include "pair_gran.h"
#include "mpi_liggghts.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define SMALL 1e-8

/* ---------------------------------------------------------------------- */

FixScalarTransportEquation::FixScalarTransportEquation(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  fix_quantity = fix_flux = fix_source = NULL; fix_capacity = NULL;
  fix_fluidQty_ = fix_transCoeffQty_ = NULL;
  fluidQty_0_       = 0;
  transCoeffQty_0_  = 0;

  capacity = NULL;

  int_flag = true;

  nevery_  = 1;
  performedIntegrationLastStep_ = true; //ensure flux is reset at the very first time step

  peratom_flag = 1;              
  size_peratom_cols = 0;         
  peratom_freq = 1;

  scalar_flag = 1; 
  global_freq = 1; 

  capacity_name = NULL;
  capacity_flag = 0;

  crankNicholsonFactor_ = 0.5; //0 ... fully explicit; 1 ... fully implicit
  implicitMode_ = false;
  advanceQty=&FixScalarTransportEquation::advanceQtyExplicit;

  if(strcmp(arg[2],"transportequation/scalar"))
    return;

  int iarg = 3;

  if (narg < 15)
    error->fix_error(FLERR,this,"not enough arguments");

  if(strcmp(arg[iarg++],"equation_id"))
    error->fix_error(FLERR,this,"expecting keyword 'equation_id'");
  equation_id = new char[strlen(arg[iarg])+1];
  strcpy(equation_id,arg[iarg++]);

  if(strcmp(arg[iarg++],"quantity"))
    error->fix_error(FLERR,this,"expecting keyword 'quantity'");
  quantity_name = new char[strlen(arg[iarg])+1];
  strcpy(quantity_name,arg[iarg++]);

  if(strcmp(arg[iarg++],"default_value"))
    error->fix_error(FLERR,this,"expecting keyword 'default_value'");
  quantity_0 = atof(arg[iarg++]);

  if(strcmp(arg[iarg++],"flux_quantity"))
    error->fix_error(FLERR,this,"expecting keyword 'flux_quantity'");
  flux_name = new char[strlen(arg[iarg])+1];
  strcpy(flux_name,arg[iarg++]);

  if(strcmp(arg[iarg++],"source_quantity"))
    error->fix_error(FLERR,this,"expecting keyword 'source_quantity'");
  source_name = new char[strlen(arg[iarg])+1];
  strcpy(source_name,arg[iarg++]);

  if(strcmp(arg[iarg++],"capacity_quantity"))
    error->fix_error(FLERR,this,"expecting keyword 'capacity_quantity'");
  if(strcmp(arg[iarg],"none"))
  {
      capacity_flag = 1;
      capacity_name = new char[strlen(arg[iarg])+1];
      strcpy(capacity_name,arg[iarg++]);
  }
}

/* ---------------------------------------------------------------------- */

FixScalarTransportEquation::~FixScalarTransportEquation()
{
    delete []quantity_name;
    delete []flux_name;
    delete []source_name;
    delete []capacity_name;
    delete []equation_id;

    if(implicitMode_)
    {
        delete []fluid_name_;
        delete []transCoeff_name_;
    }

    if(capacity) delete []capacity;
    
}

/* ---------------------------------------------------------------------- */
void FixScalarTransportEquation::register_implicit_fixes(char* fluidName, double fluid0, char* transCoeffName, double transCoeff0)
{
  const char *fixarg[9];

  implicitMode_ = true;
  advanceQty=&FixScalarTransportEquation::advanceQtyImplicit;

  fluid_name_ = new char[strlen(fluidName)+1];
  strcpy(fluid_name_,fluidName);
  fluidQty_0_ = fluid0;

  transCoeff_name_ = new char[strlen(transCoeffName)+1];
  strcpy(transCoeff_name_,transCoeffName);
  transCoeffQty_0_ = transCoeff0;

  fix_fluidQty_=static_cast<FixPropertyAtom*>(modify->find_fix_property(fluid_name_,"property/atom","scalar",0,0,style));
  if (fix_fluidQty_==NULL) {
        //register fluid quantity as property/atom
        fixarg[0]=fluid_name_;
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]=fluid_name_;
        fixarg[4]="scalar";
        fixarg[5]="yes";    //restart
        fixarg[6]="no";     //commGhost
        fixarg[7]="yes";    //commGhostRev
        char arg8[30];
        sprintf(arg8,"%e", fluidQty_0_);
        fixarg[8]=arg8;
        modify->add_fix(9,const_cast<char**>(fixarg));
        fix_fluidQty_=static_cast<FixPropertyAtom*>(modify->find_fix_property(fluid_name_,"property/atom","scalar",0,0,style));
  }

  fix_transCoeffQty_=static_cast<FixPropertyAtom*>(modify->find_fix_property(transCoeff_name_,"property/atom","scalar",0,0,style));
  if (fix_transCoeffQty_==NULL){
        //register transfer coefficient as property/atom
        fixarg[0]=transCoeff_name_;
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]=transCoeff_name_;
        fixarg[4]="scalar";
        fixarg[5]="yes";    //restart
        fixarg[6]="no";     //commGhost
        fixarg[7]="yes";    //commGhostRev
        char arg8[30];
        sprintf(arg8,"%e", transCoeffQty_0_);
        fixarg[8]=arg8;
        modify->add_fix(9,const_cast<char**>(fixarg));
        fix_transCoeffQty_=static_cast<FixPropertyAtom*>(modify->find_fix_property(transCoeff_name_,"property/atom","scalar",0,0,style));
  }

  updatePtrsImpl();

  return;
}

/* ---------------------------------------------------------------------- */

void FixScalarTransportEquation::pre_delete(bool unfixflag)
{
    //unregister property/atom fixes
    if(unfixflag)
    {
        if (fix_quantity) modify->delete_fix(quantity_name);
        if (fix_flux) modify->delete_fix(flux_name);
        if (fix_source) modify->delete_fix(source_name);
        if(implicitMode_)
        {
            if (fix_fluidQty_) modify->delete_fix(fluid_name_);
            if (fix_transCoeffQty_) modify->delete_fix(transCoeff_name_);
        }
    }
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

  if(implicitMode_)
    updatePtrsImpl();

  vector_atom = quantity; 
}

/* ---------------------------------------------------------------------- */
void FixScalarTransportEquation::updatePtrsImpl()
{
  fluidQty_         = fix_fluidQty_->vector_atom;
  transCoeffQty_    = fix_transCoeffQty_->vector_atom;
}

/* ---------------------------------------------------------------------- */

void FixScalarTransportEquation::post_create()
{
  const char *fixarg[9];

  if (fix_quantity==NULL) {
    //register Temp as property/atom
    fixarg[0]=quantity_name;
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]=quantity_name;
    fixarg[4]="scalar"; 
    fixarg[5]="yes";    
    fixarg[6]="yes";    
    fixarg[7]="no";    
    char arg8[30];
    sprintf(arg8,"%e",quantity_0);
    fixarg[8]=arg8;
    modify->add_fix(9,const_cast<char**>(fixarg));
    fix_quantity=static_cast<FixPropertyAtom*>(modify->find_fix_property(quantity_name,"property/atom","scalar",0,0,style));
  }

  if (fix_flux==NULL){
    //register heatFlux as property/atom
    fixarg[0]=flux_name;
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]=flux_name;
    fixarg[4]="scalar"; 
    fixarg[5]="yes";    
    fixarg[6]="no";    
    fixarg[7]="yes";    
    fixarg[8]="0.";     
    modify->add_fix(9,const_cast<char**>(fixarg));
    fix_flux=static_cast<FixPropertyAtom*>(modify->find_fix_property(flux_name,"property/atom","scalar",0,0,style));
  }

  if (fix_source==NULL){
    //register heatSource as property/atom
    fixarg[0]=source_name;
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]=source_name;
    fixarg[4]="scalar"; 
    fixarg[5]="yes";    
    fixarg[6]="yes";    
    fixarg[7]="no";    
    fixarg[8]="0.";     
    modify->add_fix(9,const_cast<char**>(fixarg));
    fix_source=static_cast<FixPropertyAtom*>(modify->find_fix_property(source_name,"property/atom","scalar",0,0,style));
  }

  updatePtrs();
}

/* ---------------------------------------------------------------------- */

double* FixScalarTransportEquation::get_capacity()
{
    return capacity;
}

/* ---------------------------------------------------------------------- */

void FixScalarTransportEquation::init()
{
  if (!atom->rmass_flag) error->all(FLERR,"Please use an atom style that defines per-particle mass for fix transportequation/scalar");

  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  if(capacity_flag)
  {
      int max_type = atom->get_properties()->max_type();

      if(capacity) delete []capacity;
      capacity = new double[max_type+1];

      fix_capacity = static_cast<FixPropertyGlobal*>(modify->find_fix_property(capacity_name,"property/global","peratomtype",max_type,0,style));

      //pre-calculate parameters for possible contact material combinations
      for(int i=1;i< max_type+1; i++)
          capacity[i] = fix_capacity->compute_vector(i-1);
  }
}

/* ---------------------------------------------------------------------- */

int FixScalarTransportEquation::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"integrate") == 0) {
    if (narg < 2) error->fix_error(FLERR,this,"not enough arguments for fix_modify 'integrate'");

    if (strcmp(arg[1],"start") == 0) {
      int_flag = true;
    } else if (strcmp(arg[1],"stop") == 0) {
      int_flag = false;
    } else
      error->fix_error(FLERR,this,"wrong argument for fix_modify 'integrate'");
    return 2;
  }

  if (strcmp(arg[0],"every") == 0) {
    if (narg < 2) error->fix_error(FLERR,this,"not enough arguments for fix_modify 'every'");

    nevery_ = force->inumeric(FLERR,arg[1]);
    
    return 1;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

void FixScalarTransportEquation::initial_integrate(int vflag)
{
  
  updatePtrs();

  //Skip in case there was NO Integration the last time step (to keep flux in mem)
  if(!performedIntegrationLastStep_)
        return;

  //reset heat flux
  //sources are not reset
  //int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
  {
       
       //if (mask[i] & groupbit)
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

    // skip if integration turned off
    if(!int_flag)
        return;

    // skip if integration not wanted at this timestep
    if (update->ntimestep % nevery_)
    {
        performedIntegrationLastStep_ = false;
        return;
    }

    updatePtrs();

    fix_source->do_forward_comm();

    (this->*advanceQty)();

    performedIntegrationLastStep_ = true;
}

/* ---------------------------------------------------------------------- */
void FixScalarTransportEquation::advanceQtyExplicit()
{
    double  dt = update->dt;
    int     nlocal = atom->nlocal;
    double  *rmass = atom->rmass;
    int     *type = atom->type;
    int     *mask = atom->mask;
    double  capacity;

    if(capacity_flag)
    {
        for (int i = 0; i < nlocal; i++)
        {
           if (mask[i] & groupbit){
              capacity = fix_capacity->compute_vector(type[i]-1);
              if(fabs(capacity) > SMALL) quantity[i] += (
                                                            flux[i]
                                                          + source[i]*double(nevery_) //multiply source to account for missing steps
                                                        ) * dt
                                                      / (rmass[i]*capacity);
           }
        }
    }
    else
    {
        for (int i = 0; i < nlocal; i++)
        {
           if (mask[i] & groupbit){
              quantity[i] += (
                                  flux[i]
                                + source[i]*double(nevery_) //multiply source to account for missing steps
                             ) * dt ;
           }
        }
    }
}

/* ---------------------------------------------------------------------- */
void FixScalarTransportEquation::advanceQtyImplicit()
{

    double  dt = update->dt;
    int     nlocal = atom->nlocal;
    double  *rmass = atom->rmass;
    int     *type = atom->type;
    int     *mask = atom->mask;
    double  *radius  = atom->radius;
    double  capacity;
    double  OneMinusCN = 1.0 - crankNicholsonFactor_;

    if(capacity_flag)
    {
        for (int i = 0; i < nlocal; i++)
        {
           if (mask[i] & groupbit){
              capacity = fix_capacity->compute_vector(type[i]-1);
              if(fabs(capacity) > SMALL)
              {
                     double currRadius = radius[i];
                     double termM  = dt / (rmass[i]*capacity);
                     double termMP = termM
                                   * transCoeffQty_[i]
                                   * currRadius * currRadius * 12.5663706144 ; //surface area
                     quantity[i]  = (   quantity[i]
                                      * (1.0 - termMP * OneMinusCN)
                                      + termMP * fluidQty_[i]
                                      + termM
                                      * (  flux[i]
                                         + source[i]*double(nevery_)
                                        )
                                    ) / (1.0 + termMP*crankNicholsonFactor_);
              }
           }
        }
    }
    else
    {
        for (int i = 0; i < nlocal; i++)
        {
           if (mask[i] & groupbit)
           {
                     double currRadius = radius[i];
                     double termMP = dt
                                   * transCoeffQty_[i]
                                   * currRadius * currRadius * 12.5663706144 ; //surface area
                     quantity[i]  = (   quantity[i]
                                      * (1.0 - termMP * OneMinusCN)
                                      + termMP * fluidQty_[i]
                                      + dt
                                      * (  flux[i]
                                         + source[i]*double(nevery_)
                                        )
                                    ) / (1.0 + termMP*crankNicholsonFactor_);
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
