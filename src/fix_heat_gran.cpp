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
#include "fix_heat_gran.h"
#include "atom.h"
#include "domain.h"
#include "group.h"
#include "force.h"
#include "comm.h"
#include "update.h"
#include "error.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair_gran.h"
#include "math_extra.h"
#include "fix_property_global.h"
#include "fix_property_atom.h"
#include "fix_scalar_transport_equation.h"
#include "mech_param_gran.h"
#include "respa.h"
#include "compute_pair_gran_local.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define SMALL 1e-8

/* ---------------------------------------------------------------------- */

FixHeatGran::FixHeatGran(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if ((!atom->radius_flag)||(!atom->rmass_flag)) error->all(FLERR,"Fix heat/gran needs per particle radius and mass");

  if (narg < 5)
    error->fix_error(FLERR,this,"not enough arguments");

  int iarg = 3;

  if(strcmp(arg[iarg++],"initial_temperature"))
    error->fix_error(FLERR,this,"expecting keyword 'initial_temperature'");
  T0 = atof(arg[iarg++]);

  area_correction_flag = 0;

  bool hasargs = true;
  while(iarg < narg && hasargs)
  {
    hasargs = false;
    if(strcmp(arg[iarg],"area_correction") == 0) {
        if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments for keyword 'area_correction'");
        if(strcmp(arg[iarg+1],"yes") == 0)
          area_correction_flag = 1;
        else if(strcmp(arg[iarg+1],"no") == 0)
          area_correction_flag = 0;
        else error->fix_error(FLERR,this,"");
        iarg += 2;
        hasargs = true;
    } else if(strcmp(style,"heat/gran") == 0)
        error->fix_error(FLERR,this,"unknown keyword");
  }

  fix_temp = fix_heatFlux = fix_heatSource = NULL;
  fix_conductivity = NULL;
  fix_ste = NULL;

  conductivity = NULL;

  peratom_flag = 1;              
  size_peratom_cols = 0;         
  peratom_freq = 1;
  time_depend = 1;

  scalar_flag = 1; 
  global_freq = 1; 

  cpl = NULL;

}

/* ---------------------------------------------------------------------- */

FixHeatGran::~FixHeatGran()
{
    
    if(conductivity) delete []conductivity;
}

/* ---------------------------------------------------------------------- */

void FixHeatGran::post_create()
{
    fix_ste = modify->find_fix_scalar_transport_equation("heattransfer");

    if(!fix_ste)
    {
        char **newarg = new char*[15];
        newarg[0] = (char *) "ste_heattransfer";
        newarg[1] = group->names[igroup];
        newarg[2] = (char *) "transportequation/scalar";
        newarg[3] = (char *) "equation_id";
        newarg[4] = (char *) "heattransfer";
        newarg[5] = (char *) "quantity";
        newarg[6] = (char *) "Temp";
        newarg[7] = (char *) "default_value";
        newarg[8] = new char[30];
        sprintf(newarg[8],"%f",T0);
        newarg[9] = (char *) "flux_quantity";
        newarg[10] = (char *) "heatFlux";
        newarg[11] = (char *) "source_quantity";
        newarg[12] = (char *) "heatSource";
        newarg[13] = (char *) "capacity_quantity";
        newarg[14] = (char *) "thermalCapacity";
        modify->add_fix(15,newarg);
        delete [] newarg[8];
        delete [] newarg;
    }
}

/* ---------------------------------------------------------------------- */

void FixHeatGran::pre_delete(bool unfixflag)
{

    // tell cpl that this fix is deleted
    if(cpl && unfixflag) cpl->reference_deleted();
}

/* ---------------------------------------------------------------------- */

int FixHeatGran::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */
void FixHeatGran::updatePtrs()
{
  Temp = fix_temp->vector_atom;
  vector_atom = Temp; 

  heatFlux = fix_heatFlux->vector_atom;
  heatSource = fix_heatSource->vector_atom;
}

/* ---------------------------------------------------------------------- */

void FixHeatGran::init()
{
  const double *Y,*nu,*Y_orig;
  double expo,Yeff_ij,Yeff_orig_ij,ratio;
  Fix* ymo_fix;

  if (!atom->radius_flag || !atom->rmass_flag) error->all(FLERR,"Please use a granular atom style for fix heat/gran");

  // check if a fix of this style already exists
  if(modify->n_fixes_style(style) > 1)
    error->fix_error(FLERR,this,"cannot have more than one fix of this style");

  if(!force->pair_match("gran", 0)) error->all(FLERR,"Please use a granular pair style for fix heat/gran");

  pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));
  history_flag = pair_gran->is_history();

  fix_ste = modify->find_fix_scalar_transport_equation("heattransfer");
  if(!fix_ste) error->all(FLERR,"Fix heat/gran needs a fix transportequation/scalar to work with");

  fix_temp = static_cast<FixPropertyAtom*>(modify->find_fix_property("Temp","property/atom","scalar",0,0,style));
  fix_heatFlux = static_cast<FixPropertyAtom*>(modify->find_fix_property("heatFlux","property/atom","scalar",0,0,style));
  fix_heatSource = static_cast<FixPropertyAtom*>(modify->find_fix_property("heatSource","property/atom","scalar",0,0,style));

  int max_type = pair_gran->mpg->max_type();

  if(conductivity) delete []conductivity;
  conductivity = new double[max_type];
  fix_conductivity = static_cast<FixPropertyGlobal*>(modify->find_fix_property("thermalConductivity","property/global","peratomtype",max_type,0,style));

  // pre-calculate conductivity for possible contact material combinations
  for(int i=1;i< max_type+1; i++)
      for(int j=1;j<max_type+1;j++)
      {
          conductivity[i-1] = fix_conductivity->compute_vector(i-1);
          if(conductivity[i-1] < 0.) error->all(FLERR,"Fix heat/gran: Thermal conductivity must not be < 0");
      }

  // calculate heat transfer correction

  ymo_fix = NULL;
  if(area_correction_flag)
  {
      ymo_fix = modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",0,0,style);

      if(force->pair_match("gran/hooke",0)) expo = 1.;
      else if(force->pair_match("gran/hertz",0)) expo = 2./3.;
      else error->all(FLERR,"Fix heat/gran with area correction could not identify the granular pair style you are using, supported are hooke and hertz types");

      Y = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulus","property/global","peratomtype",max_type,0,style))->get_values();
      nu = static_cast<FixPropertyGlobal*>(modify->find_fix_property("poissonsRatio","property/global","peratomtype",max_type,0,style))->get_values();
      Y_orig = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0,style))->get_values();

      // allocate a new array within youngsModulusOriginal
      static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0,style))->new_array(max_type,max_type);

      // feed deltan_ratio into this array
      for(int i = 1; i < max_type+1; i++)
      {
          for(int j = 1; j < max_type+1; j++)
          {
              Yeff_ij      = 1./((1.-pow(nu[i-1],2.))/Y[i-1]     +(1.-pow(nu[j-1],2.))/Y[j-1]);
              Yeff_orig_ij = 1./((1.-pow(nu[i-1],2.))/Y_orig[i-1]+(1.-pow(nu[j-1],2.))/Y_orig[j-1]);
              ratio = pow(Yeff_ij/Yeff_orig_ij,expo);
              
              static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0,style))->array_modify(i-1,j-1,ratio);
          }
      }

      // get reference to deltan_ratio
      deltan_ratio = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0,style))->get_array_modified();
  }

  updatePtrs();
}

/* ---------------------------------------------------------------------- */

void FixHeatGran::post_force(int vflag)
{
    //template function for using touchflag or not
    if(history_flag == 0) post_force_eval<0>(vflag,0);
    if(history_flag == 1) post_force_eval<1>(vflag,0);
}

/* ---------------------------------------------------------------------- */

void FixHeatGran::cpl_evaluate(ComputePairGranLocal *caller)
{
    if(caller != cpl) error->all(FLERR,"Illegal situation in FixHeatGran::cpl_evaluate");
    if(history_flag == 0) post_force_eval<0>(0,1);
    if(history_flag == 1) post_force_eval<1>(0,1);
}

/* ---------------------------------------------------------------------- */

template <int HISTFLAG>
void FixHeatGran::post_force_eval(int vflag,int cpl_flag)
{
  double hc,contactArea,delta_n,flux;
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double radi,radj,radsum,rsq,r,rinv,rsqinv,tcoi,tcoj;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;

  int newton_pair = force->newton_pair;

  if (strcmp(force->pair_style,"hybrid")==0)
    error->warning(FLERR,"Fix heat/gran implementation may not be valid for pair style hybrid");
  if (strcmp(force->pair_style,"hybrid/overlay")==0)
    error->warning(FLERR,"Fix heat/gran implementation may not be valid for pair style hybrid/overlay");

  inum = pair_gran->list->inum;
  ilist = pair_gran->list->ilist;
  numneigh = pair_gran->list->numneigh;
  firstneigh = pair_gran->list->firstneigh;
  if(HISTFLAG) firsttouch = pair_gran->listgranhistory->firstneigh;

  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  updatePtrs();

  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    if(HISTFLAG) touch = firsttouch[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      if (!(mask[i] & groupbit) && !(mask[j] & groupbit)) continue;

      if(!HISTFLAG)
      {
          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx*delx + dely*dely + delz*delz;
          radj = radius[j];
          radsum = radi + radj;
      }

      if (HISTFLAG && touch[jj] || !HISTFLAG && (rsq < radsum*radsum)) {  //contact
         
         if(HISTFLAG)
         {
              delx = xtmp - x[j][0];
              dely = ytmp - x[j][1];
              delz = ztmp - x[j][2];
              rsq = delx*delx + dely*dely + delz*delz;
              radj = radius[j];
              radsum = radi + radj;
         }

         r = sqrt(rsq);

         if(area_correction_flag)
         {
             delta_n = radsum - r;
             delta_n *= deltan_ratio[type[i]-1][type[j]-1];
             r = radsum - delta_n;
         }

         contactArea = - M_PI/4 * ( (r-radi-radj)*(r+radi-radj)*(r-radi+radj)*(r+radi+radj) )/(r*r); //contact area of the two spheres

         tcoi = conductivity[type[i]-1];
         tcoj = conductivity[type[j]-1];
         if (tcoi < SMALL || tcoj < SMALL) hc = 0.;
         else hc = 4.*tcoi*tcoj/(tcoi+tcoj)*sqrt(contactArea);

         flux = (Temp[j]-Temp[i])*hc;

         if(!cpl_flag)
         {
             heatFlux[i] += flux;
             if (newton_pair || j < nlocal) heatFlux[j] -= flux;
         }

         if(cpl_flag && cpl) cpl->add_heat(i,j,flux);
      }
    }
  }

  if(newton_pair) fix_heatFlux->do_reverse_comm();
}

/* ---------------------------------------------------------------------- */
double FixHeatGran::compute_scalar()
{
    return fix_ste->compute_scalar();
}

/* ----------------------------------------------------------------------
   register and unregister callback to compute
------------------------------------------------------------------------- */

void FixHeatGran::register_compute_pair_local(ComputePairGranLocal *ptr)
{
   
   if(cpl != NULL)
      error->all(FLERR,"Fix heat/gran allows only one compute of type pair/local");
   cpl = ptr;
}

void FixHeatGran::unregister_compute_pair_local(ComputePairGranLocal *ptr)
{
   
   if(cpl != ptr)
       error->all(FLERR,"Illegal situation in FixHeatGran::unregister_compute_pair_local");
   cpl = NULL;
}
