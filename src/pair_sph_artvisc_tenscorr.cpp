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
#include "pair_sph_artvisc_tenscorr.h"
#include "fix_property_global.h"
#include "fix_property_atom.h"
#include "atom.h"
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
#include "timer.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSphArtviscTenscorr::PairSphArtviscTenscorr(LAMMPS *lmp) : PairSph(lmp)
{
  respa_enable = 0;
  single_enable = 0;
  pairStyle_ = 1;

  csmean = NULL;
  wDeltaPTypeinv = NULL;
}

/* ---------------------------------------------------------------------- */

PairSphArtviscTenscorr::~PairSphArtviscTenscorr()
{
  if (allocated) {
    memory->destroy(csmean);
    if (mass_type && tensCorr_flag) memory->destroy(wDeltaPTypeinv);
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairSphArtviscTenscorr::allocate()
{
  PairSph::allocate();

  int n = atom->ntypes;

  memory->create(csmean,n+1,n+1,"pair:csmean");
  if (mass_type && tensCorr_flag) memory->create(wDeltaPTypeinv,n+1,n+1,"pair:wDeltaPTypeinv");

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSphArtviscTenscorr::settings(int narg, char **arg)
{

  int iarg = 0;

  // The first two input arguments are reserved for the kernel style and default smoothing length
  PairSph::setKernelAndLength(narg, arg);

/*
  if (iarg+1 > narg) error->all(FLERR, "Illegal pair_style sph command");

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
*/

  iarg += 2;

  // optional parameters

  artVisc_flag = tensCorr_flag = 0;
  alpha = beta = eta = epsilon = 0.0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"artVisc") == 0) {
      // parameters for artifical viscosity
      if (iarg+4 > narg) error->all(FLERR, "Illegal pair_style sph command");
      artVisc_flag = 1;
      alpha = force->numeric(arg[iarg+1]);
      viscosity_ = alpha;
      beta = force->numeric(arg[iarg+2]);
      eta = force->numeric(arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"tensCorr") == 0) {
      // parameters for tensile correction
      if (iarg+3 > narg) error->all(FLERR, "Illegal pair_style sph command");
      tensCorr_flag = 1;
      epsilon = force->numeric(arg[iarg+1]);
      deltaP = force->numeric(arg[iarg+2]);
      iarg += 3;
    } else error->all(FLERR, "Illegal pair_style sph command");
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairSphArtviscTenscorr::coeff(int narg, char **arg)
{
  if (narg > 3) error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  //if (screen) fprintf(screen,"ilo= %d, ihi= %d, jlo= %d, jhi= %d, ntypes= %d\n",ilo,ihi,jlo,jhi,atom->ntypes);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      setflag[i][j] = 1;
      count++;
      //cut[i][j] = cut[j][i] = cut_global; // sl_0?
      //cutsq[i][j] = cutsq[j][i] = cut_global*cut_global; // sl_0?
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");

}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSphArtviscTenscorr::init_substyle()
{
  int max_type = atom->ntypes;
  int i,j;

  //create wDeltaPTypeInv
  if (mass_type && tensCorr_flag) {

    double slCom,slComInv;

    //pre-calculate common smoothing length
    for(i = 1; i < max_type+1; i++) {
      for(j = 1; j < max_type+1; j++) {

        slCom = slComType[i][j];
        slComInv = 1./slCom;
        wDeltaPTypeinv[i][j] = 1./SPH_KERNEL_NS::sph_kernel(kernel_id,deltaP * slComInv,slCom,slComInv);
      }
    }
  }

  //Get pointer to the fixes that have the material properties

  if (artVisc_flag) {
    cs=static_cast<FixPropertyGlobal*>(modify->find_fix_property("speedOfSound","property/global","peratomtype",max_type,0,force->pair_style));
    if(!cs) error->all(FLERR, "Pairstyle sph/artVisc/tensCorr only works with a fix property/global that defines speedOfSound");

    double csi,csj;

    //pre-calculate parameters for possible contact material combinations
    for(i=1;i< max_type+1; i++)
    {
      for(j=1;j<max_type+1;j++)
      {
        csi=cs->compute_vector(i-1);
        csj=cs->compute_vector(j-1);

        csmean[i][j] = 0.5*(csi+csj);

      }
    }
  }
}

/* ----------------------------------------------------------------------
  allocate per-type and per-type pair properties
------------------------------------------------------------------------- */
/*
void PairSphArtviscTenscorr::allocate_properties(int size)
{
    memory->destroy_2d_double_array(csmean);
    csmean = memory->create_2d_double_array(size+1,size+1,"pair:csmean");
}
*/
/* ---------------------------------------------------------------------- */

void PairSphArtviscTenscorr::compute(int eflag, int vflag)
{
  if (mass_type) compute_eval<1>(eflag,vflag);
  else compute_eval<0>(eflag,vflag);
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSphArtviscTenscorr::write_restart(FILE *fp)
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

void PairSphArtviscTenscorr::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  PairSph::allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSphArtviscTenscorr::write_restart_settings(FILE *fp)
{
  fwrite(&artVisc_flag,sizeof(int),1,fp);
  fwrite(&tensCorr_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSphArtviscTenscorr::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&artVisc_flag,sizeof(int),1,fp);
    fread(&tensCorr_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&artVisc_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tensCorr_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   artificial viscosity according to Monaghan (1992)
   XXX: function call is slower! may be the double **v call?
------------------------------------------------------------------------- */
/*
double PairSphArtviscTenscorr::artificialViscosity(int ip, int jp, int itype, int jtype,
    double rsq, double slCom, double delx, double dely, double delz, double rhoi, double rhoj, double **v)
{
  // artifical viscosity [Monaghan, 1992]
  // alpha ... shear viscosity
  // beta  ... bulk viscosity
  // eta   ... avoid singularities ~ 0.01*h*h

  double muAB, rhoMeanInv,dotDelVDelR;

  //double **v = atom->v;

  dotDelVDelR = ( (v[ip][0]-v[jp][0])*delx + (v[ip][1]-v[jp][1])*dely + (v[ip][2]-v[jp][2])*delz );

  if ( dotDelVDelR < 0.0 ) {
    muAB = slCom * dotDelVDelR / (rsq + eta);
    rhoMeanInv = 2/(rhoi+rhoj);
    return ((- alpha * csmean[itype][jtype] * muAB + beta * muAB * muAB) * rhoMeanInv);
  } else return 0.0;
}
*/
/* ----------------------------------------------------------------------
   tensile correction according to Monaghan (2000)
   calculates the values for rAB and fAB4
------------------------------------------------------------------------- */
/*
template <int MASSFLAG>
void PairSphArtviscTenscorr::tensileCorrection(int itype, int jtype, double rhoi, double rhoj,
    double pi, double pj, double s, double slCom, double slComInv, double &rAB, double &fAB4)
{
  double rA,rB,fAB,fAB2;
  double wDeltaPinv;

  // repulsive term for tensile instability [Monaghan, 2000]
  if (pi > 0.0 && pj > 0.0) {
    rAB = 0.01 * (pi / (rhoi * rhoi) + pj / (rhoj * rhoj));
  } else {
    if (pi < 0.0) rA = epsilon * -1.0 * pi / (rhoi * rhoi);
    else rA = 0;
    if (pj < 0.0) rB = epsilon * -1.0 * pj / (rhoj * rhoj);
    else rB = 0;
    rAB = rA+rB;
  }

  if (MASSFLAG) {
    wDeltaPinv = wDeltaPTypeinv[itype][jtype];
  } else {
    wDeltaPinv = 1./SPH_KERNEL_NS::sph_kernel(kernel_id,deltaP * slComInv,slCom,slComInv);
  }

  //TODO: Is fAB4 in this form ok?!
  fAB =  SPH_KERNEL_NS::sph_kernel(kernel_id,s,slCom,slComInv) * wDeltaPinv;
  fAB2 = fAB * fAB;
  fAB4 = fAB2 * fAB2;
}
*/
/* ----------------------------------------------------------------------
   template compute
------------------------------------------------------------------------- */

template <int MASSFLAG>
void PairSphArtviscTenscorr::compute_eval(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double xtmp,ytmp,ztmp,delx,dely,delz,r,rsq,rinv,s;
  double gradWmag,fpair;

  double rhoi,rhoj,pi,pj;
  double sli,slj,slCom,slComInv,imass,jmass;

  double artVisc,fAB4,rAB;
  double muAB, rhoMeanInv,dotDelVDelR;
  double rA,rB,fAB,fAB2;
  double wDeltaPinv;

  double radi,radj,rcom;

  double **x = atom->x;
  double **v = atom->v;
  double *p = atom->p;
  double *rho = atom->rho;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  // individual properties
  double *mass = atom->mass;
  double *radius = atom->radius;
  double *rmass = atom->rmass;

  // TODO: Use this eflag/vflag .. what is it for?
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  // depend on mass_type

  // no communication of sl, in case of perAtomType mass
  if (!MASSFLAG) {
    // update smoothing length
//    timer->stamp(); // no timer step inside of a pair style!!
    fppaSl->do_forward_comm();
//    timer->stamp(TIME_COMM);

    updatePtrs(); // get sl
  }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    rhoi = rho[i];
    pi = p[i];

    if (MASSFLAG) {
      imass = mass[itype];
    } else {
      sli = sl[i];
      radi = radius[i];
      imass = rmass[i];
    }

    // derivative of kernel must be 0 at s = 0
    // so particle itself is not contributing

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      jtype = type[j];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (!MASSFLAG) {
        radj = radius[j];
        rcom = interpDist(radi,radj);
      }

      if ((MASSFLAG && rsq < cutsq[itype][jtype]) || (!MASSFLAG && rsq < rcom*rcom)) {

        if (MASSFLAG) {
          jmass = mass[jtype];
          slCom = slComType[itype][jtype];
        } else {
          jmass = rmass[j];
          slj = sl[j];
          slCom = interpDist(sli,slj);
        }

        pj = p[j];
        rhoj = rho[j];
        slComInv = 1./slCom;
        //cut = slCom*SPH_KERNEL_NS::sph_kernel_cut(kernel_id);

        // get distance and normalized distance
        r = sqrt(rsq);
        rinv = 1./r;
        s = r * slComInv;

        // calculate value for magnitude of grad W
        gradWmag = SPH_KERNEL_NS::sph_kernel_der(kernel_id,s,slCom,slComInv);

        // artificial viscosity
        artVisc = 0.0;
        if (artVisc_flag) {
          // artifical viscosity [Monaghan, 1992]
          // alpha ... shear viscosity
          // beta  ... bulk viscosity
          // eta   ... avoid singularities ~ 0.01*h*h

          dotDelVDelR = ( (v[i][0]-v[j][0])*delx + (v[i][1]-v[j][1])*dely + (v[i][2]-v[j][2])*delz );

          if ( dotDelVDelR < 0.0 ) {
            muAB = slCom * dotDelVDelR / (rsq + eta);
            rhoMeanInv = 2/(rhoi+rhoj);
            artVisc = ((- alpha * csmean[itype][jtype] * muAB + beta * muAB * muAB) * rhoMeanInv);
          }

          //version with function
          //artVisc = artificialViscosity(i,j,itype,jtype,r,slCom,delx,dely,delz,rhoi,rhoj,v); //XXX: function call is slower :-/
        }

        // tensile correction
        rAB = fAB4 = 0.0;
        if (tensCorr_flag) {
          // repulsive term for tensile instability [Monaghan, 2000]
          if (pi > 0.0 && pj > 0.0) {
            rAB = 0.01 * (pi / (rhoi * rhoi) + pj / (rhoj * rhoj));
          } else {
            if (pi < 0.0) rA = epsilon * -1.0 * pi / (rhoi * rhoi);
            else rA = 0;
            if (pj < 0.0) rB = epsilon * -1.0 * pj / (rhoj * rhoj);
            else rB = 0;
            rAB = rA+rB;
          }

          if (MASSFLAG) {
            wDeltaPinv = wDeltaPTypeinv[itype][jtype];
          } else {
            wDeltaPinv = 1./SPH_KERNEL_NS::sph_kernel(kernel_id,deltaP * slComInv,slCom,slComInv);
          }

          //TODO: Is fAB4 in this form ok?!
          fAB =  SPH_KERNEL_NS::sph_kernel(kernel_id,s,slCom,slComInv) * wDeltaPinv;
          fAB2 = fAB * fAB;
          fAB4 = fAB2 * fAB2;

          // version with function
          //tensileCorrection<MASSFLAG>(itype,jtype,rhoi,rhoj,pi,pj,s,slCom,slComInv,rAB,fAB4); //XXX: function call is slower :-/
        }

        // calculate the force
        fpair = - rinv * imass * jmass * (pi/(rhoi*rhoi) + pj/(rhoj*rhoj) + rAB*fAB4 + artVisc) * gradWmag; // mass[i] for integration.. check fix_nve.cpp

        // apply the force

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;

        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,0.0,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();

}
