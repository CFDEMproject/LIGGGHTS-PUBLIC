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
Contributing author for SPH:
Andreas Aigner (CD Lab Particulate Flow Modelling, JKU)
andreas.aigner@jku.at
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_sph.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "force.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "fix.h"
#include "integrate.h"
#include "respa.h"
#include "memory.h"
#include "error.h"
#include "sph_kernels.h"
#include "fix_property_atom.h"
#include "fix_property_global.h"
#include "timer.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSph::PairSph(LAMMPS *lmp) : Pair(lmp)
{
    single_enable = 0;
 //   no_virial_compute = 1;

    calcMode_ = 0;

    kernel_style = NULL;

    fppaSl = NULL;
    fppaSlType = NULL;
    sl = NULL;
    slComType = NULL;

    mass_type = atom->avec->mass_type; // get flag for mass per type
}

/* ---------------------------------------------------------------------- */

PairSph::~PairSph()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    if (mass_type) memory->destroy(slComType);
  }

  delete [] maxrad;
  delete [] onerad;

  if(kernel_style) delete []kernel_style;
  if(fppaSl) modify->delete_fix("sl");
//  if(fppaSlType) modify->delete_fix("sl");
}

/* ----------------------------------------------------------------------
    set kernel id and smoothing length
------------------------------------------------------------------------- */

void PairSph::setKernelAndLength(int narg, char **arg)
{
  int iarg = 0;

  if (iarg+2 > narg) error->all(FLERR, "Illegal pair_style sph command");

  // kernel style

  if(kernel_style) delete []kernel_style;
  kernel_style = new char[strlen(arg[iarg])+1];
  strcpy(kernel_style,arg[iarg]);

  // check uniqueness of kernel IDs

  int flag = SPH_KERNEL_NS::sph_kernels_unique_id();
  if(flag < 0) error->all(FLERR, "Cannot proceed, sph kernels need unique IDs, check all sph_kernel_* files");

  // get kernel id

  kernel_id = SPH_KERNEL_NS::sph_kernel_id(kernel_style);
  if(kernel_id < 0) error->all(FLERR, "Illegal pair_style sph command, unknown sph kernel");

  // smoothing length

  //if(strcmp(arg[iarg+1],"default_value")) error->all(FLERR, "Fix transportequation/scalar: expecting keyword 'default_value'");
  sl_0 = force->numeric(arg[iarg+1]); // in case of mass_type == 1, this is only a dummy!
  iarg += 2;

}

/* ---------------------------------------------------------------------- */

void PairSph::updatePtrs()
{
  if(fppaSl) sl = fppaSl->vector_atom;
  if(fppaSlType) sl = fppaSlType->values;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairSph::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  if (mass_type) memory->create(slComType,n+1,n+1,"pair:slComType");

  onerad = new double[n+1];
  maxrad = new double[n+1];
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairSph::coeff(int narg, char **arg)
{
  if (narg > 2) error->all(FLERR,"Incorrect args for pair sph coefficients");
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

  if (count == 0) error->all(FLERR,"Incorrect args for pair sph coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSph::init_style()
{
  int i,j= 0;
  int ntypes = atom->ntypes;

  // check for mass, rho
  if(!atom->rho_flag || !atom->p_flag) error->all(FLERR,"Pair sph requires atom_style sph");

  // request neighbor list
  // if mass_type -> regular (half) list
  // if !mass_type -> granular list (neighborlist build by radius)
  int irequest = neighbor->request(this);
  if (mass_type) {
    neighbor->requests[irequest]->half = 1; //default
  } else {
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->gran = 1;
  }

  // check for fixes
  int ifix_p = -1, ifix_d = -1;
  for (int ifix = 0; ifix < modify->nfix; ifix++)
  {
    if (strncmp("sph/density",modify->fix[ifix]->style,11) == 0) ifix_d = ifix;
    if (strcmp("sph/pressure",modify->fix[ifix]->style) == 0) ifix_p = ifix;
  }
  if (ifix_d == -1) error->all(FLERR,"Pair sph requires a fix sph/density");
  if (ifix_p == -1) error->all(FLERR,"Pair sph requires a fix sph/pressure");

  // init individual for mass_type 0/1
  // mass_type = 1 ... get fppaSlType, the perAtomType smoothing length
  //                   define cutsq und cutforce
  // mass_type = 0 ... get fppaSl, the perAtom smoothing length
  //                   define maxrad for all types
  if (mass_type) {
    // get pointer to the smoothing length pointer

    fppaSlType=static_cast<FixPropertyGlobal*>(modify->find_fix_property("sl","property/global","peratomtype",ntypes,0,force->pair_style));
    if(!fppaSlType) error->all(FLERR, "Pairstyle sph only works with a fix property/global that defines sl");

    for (i = 1; i <= ntypes; i++)
      for (j = i; j <= ntypes; j++) {
        double sli = fppaSlType->compute_vector(i-1);
        double slj = fppaSlType->compute_vector(j-1);

        slComType[i][j] = slComType[j][i] = interpDist(sli,slj);;
}

  } else {
    // register per-particle property smoothing length

    if (fppaSl == NULL) {
      char **fixarg = new char*[9];
      for (int kk=0;kk<9;kk++) fixarg[kk]=new char[30];

      fixarg[0]=(char *) "sl";
      fixarg[1]=(char *) "all";
      fixarg[2]=(char *) "property/atom";
      fixarg[3]=(char *) "sl";
      fixarg[4]=(char *) "scalar";
      fixarg[5]=(char *) "yes"; // restart_peratom = 1
      fixarg[6]=(char *) "yes"; // commGhost = 1
      fixarg[7]=(char *) "no";  // commGhostRev = 0
      sprintf(fixarg[8],"%f",sl_0); // default value
      modify->add_fix(9,fixarg);
      delete []fixarg;

      fppaSl=static_cast<FixPropertyAtom*>(modify->find_fix_property("sl","property/atom","scalar",0,0,force->pair_style));
    }
    if(!fppaSl) error->all(FLERR, "Pairstyle sph only works with a internal fix property/atom that defines sl.");

    // communicate to other procs
    // TODO: Is this necessary? All procs have the same sl_0.
    timer->stamp();
    fppaSl->do_forward_comm();
    timer->stamp(TIME_COMM);

    // update pointers
    updatePtrs();

    // update radius
    updateRadius();

    // set maxrad_dynamic and maxrad_frozen for each type
    for (i = 1; i <= atom->ntypes; i++)  onerad[i] = 0.0;

    // include future Fix pour particles as dynamic

    for (i = 0; i < modify->nfix; i++){
      for(int j=1;j<=atom->ntypes;j++)
      {
          int pour_type = 0;
          double pour_maxrad = 0.0;
          pour_type = j;
          pour_maxrad = modify->fix[i]->max_rad(pour_type);
          onerad[pour_type] = MAX(onerad[pour_type],pour_maxrad);
      }
    }

    //further dynamic and frozen

    double *radius = atom->radius;
    int *mask = atom->mask;
    int *type = atom->type;
    int nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++)
      if (mask[i])
        onerad[type[i]] = MAX(onerad[type[i]],radius[i]);

    MPI_Allreduce(&onerad[1],&maxrad[1],atom->ntypes,MPI_DOUBLE,MPI_MAX,world);
  }

  // proceed with initialisation of the substyle
  init_substyle();

}

/* ---------------------------------------------------------------------- */

void PairSph::init_list(int id, NeighList *ptr)
{
  if (id == 0) list = ptr;
  else error->all(FLERR,"Internal error in PairSph::init_list");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairSph::init_one(int i, int j)
{
  if (!allocated) allocate();

  if (setflag[i][j] == 0) {
    error->all(FLERR,"PairSph: Illegal pair_style sph command,");
  }

  if (mass_type) {
    return slComType[i][j]*SPH_KERNEL_NS::sph_kernel_cut(kernel_id); //return cutoff
  } else {
    return interpDist(maxrad[i],maxrad[j]); //return cutoff .. until now dummy!
  }
}

/* ----------------------------------------------------------------------
   calculate radius
   radius[i] = kernel_cut * sl[i];
   only valid for atom styles with mass_type = 0
------------------------------------------------------------------------- */

void PairSph::updateRadius()
{
  int i;

  // TODO: include future Fix pour particles

  //find local maximum smoothing length
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // set radius for each particle
  for (i = 0; i < nlocal; i++)
    if (mask[i])
      atom->radius[i] = sl[i]*SPH_KERNEL_NS::sph_kernel_cut(kernel_id);

  // update ghost particles
  timer->stamp();
  comm->forward_comm();
  timer->stamp(TIME_COMM);
}

/* ----------------------------------------------------------------------
   return common radius for types with different smoothing lengths
   atm: simple arithmetic mean (compare Morris)
------------------------------------------------------------------------- */
// XXX: .. inserted in pair_sph.h
/*
double PairSph::interpDist(double disti, double distj)
{
  return 0.5*(disti+distj);
}
*/
