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
    This file is from LAMMPS, but has been modified. Copyright for
    modification:

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
    Copyright 2016-     CFDEMresearch GmbH, Linz

    Copyright of original file:
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "modify.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "fix_contact_history.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "properties.h"
#include "fix_rigid.h"
#include "fix_particledistribution_discrete.h"
#include "fix_property_global.h"
#include "fix_property_atom.h"
#include "fix_contact_property_atom.h"
#include "compute_pair_gran_local.h"
#include "pair_gran.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairGran::PairGran(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  no_virial_fdotr_compute = 1;

  suffix = NULL;
  neighprev = 0;

  history = 0;
  dnum_pairgran = 0;
  dnum_all = 0;
  fix_history = NULL;

  nfix = 0;
  fix_dnum = NULL;
  dnum_index = NULL;

  mass_rigid = NULL;

  nmax = 0;

  cpl_enable = 1;
  cpl_ = NULL;

  energytrack_enable = 0;
  fppaCPEn = fppaCDEn = fppaCPEt = fppaCDEVt = fppaCDEFt = fppaCTFW = fppaDEH = NULL;
  CPEn = CDEn = CPEt = CDEVt = CDEFt = CTFW = DEH = NULL;

  computeflag_ = 0;

  needs_neighlist = true;

  fix_contact_forces_ = NULL;
  store_contact_forces_ = false;
  store_contact_forces_every_ = 1;
  fix_contact_forces_stress_ = NULL;
  store_contact_forces_stress_ = false;
  fix_store_multicontact_data_ = NULL;
  store_multicontact_data_ = false;

  fix_relax_ = NULL;

  if(modify->n_fixes_style("multisphere/advanced"))
    do_store_contact_forces();
}

/* ---------------------------------------------------------------------- */

PairGran::~PairGran()
{
  if (fix_history) modify->delete_fix("contacthistory");
  if (fix_contact_forces_) modify->delete_fix(fix_contact_forces_->id);
  if (fix_contact_forces_stress_) modify->delete_fix(fix_contact_forces_stress_->id);
  if (fix_store_multicontact_data_) modify->delete_fix(fix_store_multicontact_data_->id);

  if (suffix) delete[] suffix;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete [] onerad_dynamic;
    delete [] onerad_frozen;
    delete [] maxrad_dynamic;
    delete [] maxrad_frozen;
  }

  // tell cpl that pair gran is deleted
  if(cpl_) cpl_->reference_deleted();

  //unregister energy terms as property/atom
  if (fppaCPEn) modify->delete_fix("CPEn");
  if (fppaCDEn) modify->delete_fix("CDEn");
  if (fppaCPEt) modify->delete_fix("CPEt");
  if (fppaCDEVt) modify->delete_fix("CDEVt");
  if (fppaCDEFt) modify->delete_fix("CDEFt");
  if (fppaCTFW) modify->delete_fix("CTFW");
  if (fppaDEH) modify->delete_fix("DEH");

  if(fix_dnum) delete []fix_dnum;
  if(dnum_index) delete []dnum_index;
}

/* ---------------------------------------------------------------------- */

void PairGran::updatePtrs()
{
   if(fppaCPEn) CPEn = fppaCPEn->vector_atom;
   if(fppaCDEn) CDEn = fppaCDEn->vector_atom;
   if(fppaCPEt) CPEt = fppaCPEt->vector_atom;
   if(fppaCDEVt) CDEVt = fppaCDEVt->vector_atom;
   if(fppaCDEFt) CDEFt = fppaCDEFt->vector_atom;
   if(fppaCTFW) CTFW = fppaCTFW->vector_atom;
   if(fppaDEH) DEH = fppaDEH->vector_atom;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairGran::allocate()
{
  allocated = 1;
  int n = atom->ntypes; //ensured elsewhere that this is high enough

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  onerad_dynamic = new double[n+1];
  onerad_frozen = new double[n+1];
  maxrad_dynamic = new double[n+1];
  maxrad_frozen = new double[n+1];
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairGran::coeff(int narg, char **arg)
{
  if (narg > 2) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairGran::init_style()
{
  int i;

  // error and warning checks

  if(0 == comm->me && 0 != neighbor->delay)
    error->warning(FLERR,"It is heavily recommended to use 'neigh_modify delay 0' with granular pair styles");

  if(strcmp(update->unit_style,"metal") ==0 || strcmp(update->unit_style,"real") == 0)
    error->all(FLERR,"Cannot use a non-consistent unit system with pair gran. Please use si,cgs or lj.");

  if (!atom->sphere_flag && !atom->superquadric_flag && !atom->shapetype_flag)
    error->all(FLERR,"Pair granular requires atom style sphere, superquadric or convexhull");
  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair granular requires ghost atoms store velocity");

  // check if a fix rigid is registered - important for damp
  
  fix_rigid = static_cast<FixRigid*>(modify->find_fix_style_strict("rigid",0));

  if(modify->n_fixes_style("rigid") > 1)
    error->warning(FLERR,"Pair gran does currently not support more than one fix rigid. This may result in under-damping.");

  dt = update->dt;

  // if shear history is stored:
  // check if newton flag is valid
  // if first init, create Fix needed for storing shear history

  if (history && force->newton_pair == 1)
    error->all(FLERR,"Pair granular with shear history requires newton pair off");

  dnum_all = dnum_pairgran;

  int dnum_extra = 0;
  nfix = modify->nfix;

  if(fix_dnum) delete []fix_dnum;
  if(dnum_index) delete []dnum_index;
  fix_dnum = new Fix*[nfix];
  dnum_index = new int[nfix];

  for(int ifix = 0; ifix < nfix; ifix++)
  {
      fix_dnum[ifix] = modify->fix[ifix];
      dnum_index[ifix] = dnum_pairgran + dnum_extra;
      
      dnum_extra += fix_dnum[ifix]->n_history_extra();
  }

  if(dnum_extra && !history)
    error->all(FLERR,"Cannot use extra dnum with non-history pair style");
  else
    dnum_all += dnum_extra;

  // init contact history
  if(history)
  {
    if (!fix_history)
    {
        
        fix_history = static_cast<FixContactHistory*>(modify->find_fix_style_strict("contacthistory",0));
        if(fix_history)
            error->all(FLERR,"Can not use more than one pair style that uses contact history");

        char **fixarg = new char*[4+2*dnum_all];
        fixarg[3] = new char[3];

        fixarg[0] = (char *) "contacthistory";
        fixarg[1] = (char *) "all";
        fixarg[2] = (char *) "contacthistory";
        sprintf(fixarg[3],"%d",dnum_all);
        history_args(&fixarg[4]);

        // add history args from fixes
        for(int ifix = 0; ifix < nfix; ifix++)
        {
            if(fix_dnum[ifix]->n_history_extra())
                if(!fix_dnum[ifix]->history_args(&fixarg[4+2*dnum_index[ifix]]))
                    error->all(FLERR,"Missing history_args() implementation for fix that requests extra history");
        }

        modify->add_fix(4+2*dnum_all,fixarg);

        if(modify->n_fixes_style_strict("contacthistory") != 1) error->all(FLERR,"Pair granular with shear history requires exactly one fix of style contacthistory");
        fix_history = static_cast<FixContactHistory*>(modify->find_fix_style_strict("contacthistory",0));

        delete [] fixarg[3];
        delete [] fixarg;
    }
  }

  // register per-particle properties for energy tracking
  if(energytrack_enable)
  {
      if(comm->nprocs > 1) error->all(FLERR,"check communication settings");
      char **fixarg = new char*[9];
      if (fppaCPEn == NULL) {
        fixarg[0]=(char *) "CPEn";
        fixarg[1]=(char *) "all";
        fixarg[2]=(char *) "property/atom";
        fixarg[3]=(char *) "CPEn";
        fixarg[4]=(char *) "scalar";
        fixarg[5]=(char *) "yes";
        fixarg[6]=(char *) "no";
        fixarg[7]=(char *) "no";
        fixarg[8]=(char *) "0.0";
        modify->add_fix(9,fixarg);
        fppaCPEn=static_cast<FixPropertyAtom*>(modify->find_fix_property("CPEn","property/atom","scalar",0,0,force->pair_style));
      }
      if (fppaCDEn == NULL) {
        fixarg[0]=(char *) "CDEn";
        fixarg[1]=(char *) "all";
        fixarg[2]=(char *) "property/atom";
        fixarg[3]=(char *) "CDEn";
        fixarg[4]=(char *) "scalar";
        fixarg[5]=(char *) "yes";
        fixarg[6]=(char *) "no";
        fixarg[7]=(char *) "no";
        fixarg[8]=(char *) "0.0";
        modify->add_fix(9,fixarg);
        fppaCDEn=static_cast<FixPropertyAtom*>(modify->find_fix_property("CDEn","property/atom","scalar",0,0,force->pair_style));
      }
      if (fppaCPEt == NULL) {
          fixarg[0]=(char *) "CPEt";
          fixarg[1]=(char *) "all";
          fixarg[2]=(char *) "property/atom";
          fixarg[3]=(char *) "CPEt";
          fixarg[4]=(char *) "scalar";
          fixarg[5]=(char *) "yes";
          fixarg[6]=(char *) "no";
          fixarg[7]=(char *) "no";
          fixarg[8]=(char *) "0.0";
          modify->add_fix(9,fixarg);
          fppaCPEt=static_cast<FixPropertyAtom*>(modify->find_fix_property("CPEt","property/atom","scalar",0,0,force->pair_style));
       }
       if (fppaCDEVt == NULL) {
         fixarg[0]=(char *) "CDEVt";
         fixarg[1]=(char *) "all";
         fixarg[2]=(char *) "property/atom";
         fixarg[3]=(char *) "CDEVt";
         fixarg[4]=(char *) "scalar";
         fixarg[5]=(char *) "yes";
         fixarg[6]=(char *) "no";
         fixarg[7]=(char *) "no";
         fixarg[8]=(char *) "0.0";
         modify->add_fix(9,fixarg);
         fppaCDEVt=static_cast<FixPropertyAtom*>(modify->find_fix_property("CDEVt","property/atom","scalar",0,0,force->pair_style));
       }
       if (fppaCDEFt == NULL) {
         fixarg[0]=(char *) "CDEFt";
         fixarg[1]=(char *) "all";
         fixarg[2]=(char *) "property/atom";
         fixarg[3]=(char *) "CDEFt";
         fixarg[4]=(char *) "scalar";
         fixarg[5]=(char *) "yes";
         fixarg[6]=(char *) "no";
         fixarg[7]=(char *) "no";
         fixarg[8]=(char *) "0.0";
         modify->add_fix(9,fixarg);
         fppaCDEFt=static_cast<FixPropertyAtom*>(modify->find_fix_property("CDEFt","property/atom","scalar",0,0,force->pair_style));
       }
       if (fppaCTFW == NULL) {
         fixarg[0]=(char *) "CTFW";
         fixarg[1]=(char *) "all";
         fixarg[2]=(char *) "property/atom";
         fixarg[3]=(char *) "CTFW";
         fixarg[4]=(char *) "scalar";
         fixarg[5]=(char *) "yes";
         fixarg[6]=(char *) "no";
         fixarg[7]=(char *) "no";
         fixarg[8]=(char *) "0.0";
         modify->add_fix(9,fixarg);
         fppaCTFW=static_cast<FixPropertyAtom*>(modify->find_fix_property("CTFW","property/atom","scalar",0,0,force->pair_style));
      }
      if (fppaDEH == NULL) {
       fixarg[0]=(char *) "DEH";
       fixarg[1]=(char *) "all";
       fixarg[2]=(char *) "property/atom";
       fixarg[3]=(char *) "DEH";
       fixarg[4]=(char *) "scalar";
       fixarg[5]=(char *) "yes";
       fixarg[6]=(char *) "no";
       fixarg[7]=(char *) "no";
       fixarg[8]=(char *) "0.0";
       modify->add_fix(9,fixarg);
       fppaDEH=static_cast<FixPropertyAtom*>(modify->find_fix_property("DEH","property/atom","scalar",0,0,force->pair_style));
     }
    delete []fixarg;

    if(force->newton_pair == 1) error->all(FLERR,"Have to implement rev comm of energy terms");
  }

  // register per-particle properties for energy tracking
  if(store_contact_forces_ && fix_contact_forces_ == NULL)
  {
    char **fixarg = new char*[19];
    fixarg[0]=(char *) "contactforces";
    fixarg[1]=(char *) "all";
    fixarg[2]=(char *) "contactproperty/atom";
    fixarg[3]=(char *) "contactforces";
    fixarg[4]=(char *) "6";
    fixarg[5]=(char *) "fx";
    fixarg[6]=(char *) "0";
    fixarg[7]=(char *) "fy";
    fixarg[8]=(char *) "0";
    fixarg[9]=(char *) "fz";
    fixarg[10]=(char *) "0";
    fixarg[11]=(char *) "tx";
    fixarg[12]=(char *) "0";
    fixarg[13]=(char *) "ty";
    fixarg[14]=(char *) "0";
    fixarg[15]=(char *) "tz";
    fixarg[16]=(char *) "0";
    fixarg[17]=(char *) "reset";
    fixarg[18]=(char *) "yes";
    modify->add_fix(19,fixarg);
    fix_contact_forces_ = static_cast<FixContactPropertyAtom*>(modify->find_fix_id("contactforces"));
    delete []fixarg;
  }

  // register per-particle properties for energy tracking
  if(store_contact_forces_stress_ && fix_contact_forces_stress_ == NULL)
  {
    char **fixarg = new char*[15];
    fixarg[0]=(char *) "contactforces_stress_";
    fixarg[1]=(char *) "all";
    fixarg[2]=(char *) "contactproperty/atom";
    fixarg[3]=(char *) "contactforces_stress_";
    fixarg[4]=(char *) "4";
    fixarg[5]=(char *) "fx";
    fixarg[6]=(char *) "0";
    fixarg[7]=(char *) "fy";
    fixarg[8]=(char *) "0";
    fixarg[9]=(char *) "fz";
    fixarg[10]=(char *) "0";
    fixarg[11]=(char *) "j";
    fixarg[12]=(char *) "0";
    fixarg[13]=(char *) "reset";
    fixarg[14]=(char *) "yes";
    modify->add_fix(15,fixarg);
    fix_contact_forces_stress_ = static_cast<FixContactPropertyAtom*>(modify->find_fix_id("contactforces_stress_"));
    delete []fixarg;
  }

    if(store_multicontact_data_ && fix_store_multicontact_data_ == NULL)
    {
        // create a new per contact property which will contain the data for the computation according to Brodu et. al. 2016
        // surfPosIJ will contain the position of the contact surface ij, realtive to position i
        // normalForce will contain the normal component of the contact force
        char **fixarg = new char*[15];
        fixarg[0]=(char *) "multicontactData_";
        fixarg[1]=(char *) "all";
        fixarg[2]=(char *) "contactproperty/atom";
        fixarg[3]=(char *) "multicontactData_";
        fixarg[4]=(char *) "4";
        fixarg[5]=(char *) "surfPosIJ_x";
        fixarg[6]=(char *) "0";
        fixarg[7]=(char *) "surfPosIJ_y";
        fixarg[8]=(char *) "0";
        fixarg[9]=(char *) "surfPosIJ_z";
        fixarg[10]=(char *) "0";
        fixarg[11]=(char *) "normalForce";
        fixarg[12]=(char *) "0";
        fixarg[13]=(char *) "reset";
        fixarg[14]=(char *) "no";
        modify->add_fix(15,fixarg);
        fix_store_multicontact_data_ = static_cast<FixContactPropertyAtom*>(modify->find_fix_id("multicontactData_"));
    delete []fixarg;
  }

  fix_sum_normal_force_ =
      static_cast<FixPropertyAtom*>
      (
          modify->find_fix_property("sum_normal_force_","property/atom","scalar",0,0, "pair/gran", false)
      );

  // need a gran neigh list and optionally a granular history neigh list

  if(needs_neighlist)
  {
      int irequest = neighbor->request(this);
      neighbor->requests[irequest]->pairgran_hashcode = hashcode();
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->gran = 1;

      if (history) {
        irequest = neighbor->request(this);
        neighbor->requests[irequest]->pairgran_hashcode = hashcode();
        neighbor->requests[irequest]->id = 1; 
        neighbor->requests[irequest]->half = 0;
        neighbor->requests[irequest]->granhistory = 1;
        neighbor->requests[irequest]->dnum = dnum_all;
      }
  }

  // check for freeze Fix and set freeze_group_bit_

  for (i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"freeze") == 0) break;
  if (i < modify->nfix) freeze_group_bit_ = modify->fix[i]->groupbit;
  else freeze_group_bit_ = 0;

  // set maxrad_dynamic and maxrad_frozen for each type
  for (i = 1; i <= atom->ntypes; i++)
  onerad_dynamic[i] = onerad_frozen[i] = 0.0;

  // include future particles as dynamic

  for (i = 0; i < modify->nfix; i++)
  {
    
    if(!modify->fix[i]->use_rad_for_cut_neigh_and_ghost())
        continue;

    for(int j=1;j<=atom->ntypes;j++)
    {
        int pour_type = 0;
        double pour_maxrad = 0.0;
        pour_type = j;
        pour_maxrad = modify->fix[i]->max_rad(pour_type);
        onerad_dynamic[pour_type] = MAX(onerad_dynamic[pour_type],pour_maxrad);
    }
  }

  // further dynamic and frozen

  double *radius = atom->radius;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & freeze_group_bit_)
      onerad_frozen[type[i]] = MAX(onerad_frozen[type[i]],radius[i]);
    else
      onerad_dynamic[type[i]] = MAX(onerad_dynamic[type[i]],radius[i]);

  MPI_Allreduce(&onerad_dynamic[1],&maxrad_dynamic[1],atom->ntypes,MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(&onerad_frozen[1],&maxrad_frozen[1],atom->ntypes,MPI_DOUBLE,MPI_MAX,world);

  init_granular();
}

/* ----------------------------------------------------------------------
   register and unregister callback to compute
------------------------------------------------------------------------- */

void PairGran::register_compute_pair_local(ComputePairGranLocal *ptr,int &dnum_compute)
{
   if(cpl_ != NULL) error->all(FLERR,"Pair gran allows only one compute of type pair/local");
   cpl_ = ptr;
   dnum_compute = dnum_pairgran; //history values
}

void PairGran::unregister_compute_pair_local(ComputePairGranLocal *ptr)
{
   if(cpl_ != ptr) error->all(FLERR,"Illegal situation in PairGran::unregister_compute_pair_local");
   cpl_ = NULL;
}

/* ----------------------------------------------------------------------
   return index for extra dnum
------------------------------------------------------------------------- */

int PairGran::fix_extra_dnum_index(class Fix *fix)
{
    for(int ifix = 0; ifix < nfix; ifix++)
    {
        if(fix == fix_dnum[ifix])
            return dnum_index[ifix];
    }
    error->all(FLERR,"Internal error - illegal fix_extra_dnum_index() call");
    return 0;
}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   optional granular history list
------------------------------------------------------------------------- */

void PairGran::init_list(int id, NeighList *ptr)
{
  
  if (id == 0) list = ptr;
  else if (id == 1) listgranhistory = ptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairGran::init_one(int i, int j)
{
  if (!allocated) allocate();

  // cutoff = sum of max I,J radii for
  // dynamic/dynamic & dynamic/frozen interactions, but not frozen/frozen

  double cutoff = maxrad_dynamic[i]+maxrad_dynamic[j];
  cutoff = MAX(cutoff,maxrad_frozen[i]+maxrad_dynamic[j]);
  cutoff = MAX(cutoff,maxrad_dynamic[i]+maxrad_frozen[j]);

  return cutoff * neighbor->contactDistanceFactor;
}

/* ----------------------------------------------------------------------
   compute as called via force
------------------------------------------------------------------------- */

void PairGran::compute(int eflag, int vflag)
{
   if(forceoff()) return;

  // update rigid body info for owned & ghost atoms if using FixRigid masses
  // body[i] = which body atom I is in, -1 if none
  // mass_body = mass of each rigid body

  if (fix_rigid && neighbor->ago == 0) {
    int tmp;
    int *body = (int *) fix_rigid->extract("body",tmp);
    double *mass_body = (double *) fix_rigid->extract("masstotal",tmp);
    if (atom->nmax > nmax) {
      memory->destroy(mass_rigid);
      nmax = atom->nmax;
      memory->create(mass_rigid,nmax,"pair:mass_rigid");
    }
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++)
      if (body[i] >= 0) mass_rigid[i] = mass_body[body[i]];
      else mass_rigid[i] = 0.0;
    comm->forward_comm_pair(this);
  }

   computeflag_ = 1;
   shearupdate_ = 1;
   if (update->setupflag) shearupdate_ = 0;

   compute_force(eflag,vflag,0);
}

/* ----------------------------------------------------------------------
   compute as called via compute pair gran local
------------------------------------------------------------------------- */

void PairGran::compute_pgl(int eflag, int vflag)
{
  
  // update rigid body info for owned & ghost atoms if using FixRigid masses
  // body[i] = which body atom I is in, -1 if none
  // mass_body = mass of each rigid body

  if (fix_rigid && neighbor->ago == 0) {
    int tmp;
    int *body = (int *) fix_rigid->extract("body",tmp);
    double *mass_body = (double *) fix_rigid->extract("masstotal",tmp);
    if (atom->nmax > nmax) {
      memory->destroy(mass_rigid);
      nmax = atom->nmax;
      memory->create(mass_rigid,nmax,"pair:mass_rigid");
    }
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++)
      if (body[i] >= 0) mass_rigid[i] = mass_body[body[i]];
      else mass_rigid[i] = 0.0;
    comm->forward_comm_pair(this);
  }

  bool reset_computeflag = (computeflag_ == 1) ? true : false;

  computeflag_ = 0;
  shearupdate_ = 0;

  compute_force(eflag,vflag,1);

  if(reset_computeflag)
    computeflag_ = 1;
}

/* ---------------------------------------------------------------------- */

int PairGran::pack_comm(int n, int *list,double *buf, int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = mass_rigid[j];
  }
  return 1;
}

/* ---------------------------------------------------------------------- */

void PairGran::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    mass_rigid[i] = static_cast<int> (buf[m++]);
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairGran::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++)
      fwrite(&setflag[i][j],sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairGran::read_restart(FILE *fp, const int major, const int minor)
{
  read_restart_settings(fp, major, minor);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
    }
}

/* ---------------------------------------------------------------------- */

void PairGran::reset_dt()
{
  dt = update->dt;
}

/* ---------------------------------------------------------------------- */

void *PairGran::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"computeflag") == 0) return (void *) &computeflag_;
  return NULL;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairGran::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}

bool PairGran::forceoff() {
  return false;
}

void PairGran::cpl_add_pair(LCM::SurfacesIntersectData & sidata, LCM::ForceData & i_forces)
{
    const double fx = i_forces.delta_F[0];
    const double fy = i_forces.delta_F[1];
    const double fz = i_forces.delta_F[2];
    const double tor1 = i_forces.delta_torque[0];
    const double tor2 = i_forces.delta_torque[1];
    const double tor3 = i_forces.delta_torque[2];
    const double * const contact_point =
#ifdef NONSPHERICAL_ACTIVE_FLAG
        atom->shapetype_flag ? sidata.contact_point : NULL;
#else
        NULL;
#endif
    cpl_->add_pair(sidata.i, sidata.j, fx,fy,fz,tor1,tor2,tor3,sidata.contact_history, contact_point);
}

void PairGran::cpl_pair_finalize()
{
    cpl_->pair_finalize();
}
