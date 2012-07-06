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

#include "domain.h"
#include "lattice.h"

using namespace LAMMPS_NS;

#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

PairSPH::PairSPH(LAMMPS *lmp) : Pair(lmp)
{
    kernel_style = NULL;
    single_enable = 0;
}

/* ---------------------------------------------------------------------- */

PairSPH::~PairSPH()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
  }
  if(kernel_style) delete []kernel_style;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairSPH::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut,n+1,n+1,"pair:cut");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSPH::settings(int narg, char **arg)
{

  int iarg = 0;

  if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style sph command");

  // kernel style

  if(kernel_style) delete []kernel_style;
  kernel_style = new char[strlen(arg[iarg])+1];
  strcpy(kernel_style,arg[iarg]);

  // check uniqueness of kernel IDs

  int flag = SPH_KERNEL_NS::sph_kernels_unique_id();
  if(flag < 0) error->all(FLERR,"Cannot proceed, sph kernels need unique IDs, check all sph_kernel_* files");

  // get kernel id

  kernel_id = SPH_KERNEL_NS::sph_kernel_id(kernel_style);
  if(kernel_id < 0) error->all(FLERR,"Illegal pair_style sph command, unknown sph kernel");

  // get h

  h = force->numeric(arg[iarg+1]);
  hinv = 1./h;

  // get cutoff

  cut_global = SPH_KERNEL_NS::sph_kernel_cut(kernel_id) * h;
  //if (screen) fprintf(screen,"cut_global %f\n",cut_global);

  iarg += 2;

  // optional parameters

  artVisc_flag = tensCorr_flag = 0;
  alpha = beta = cAB = eta = epsilon = 0.0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"artVisc") == 0) {
      // parameters for artifical viscosity
      if (iarg+5 > narg) error->all(FLERR,"Illegal pair_style sph command");
      artVisc_flag = 1;
      alpha = force->numeric(arg[iarg+1]);
      beta = force->numeric(arg[iarg+2]);
      cAB = force->numeric(arg[iarg+3]);
      eta = force->numeric(arg[iarg+4]);
      iarg += 5;
    } else if (strcmp(arg[iarg],"tensCorr") == 0) {
      // parameters for tensile correction
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style sph command");
      tensCorr_flag = 1;
      epsilon = force->numeric(arg[iarg+1]);
      iarg += 2;
    } else iarg++;
  }

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairSPH::coeff(int narg, char **arg)
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
      cut[i][j] = cut_global;
      cutsq[i][j] = cut_global*cut_global;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair sph coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSPH::init_style()
{
  // check for mass, density

  if(!atom->density_flag || !atom->q_flag) error->all(FLERR,"Pair sph requires atom_style sph");

  // request regular neighbor list
  // is a half list
  int irequest = neighbor->request(this);

  // check for fixes
  int ifix_p = -1, ifix_d = -1;
  for (int ifix = 0; ifix < modify->nfix; ifix++)
  {
    if (strncmp("sph/density",modify->fix[ifix]->style,11) == 0) ifix_d = ifix;
    if (strcmp("sph/pressure",modify->fix[ifix]->style) == 0) ifix_p = ifix;
  }
  if (ifix_d == -1) error->all(FLERR,"Pair sph requires a fix sph/density");
  if (ifix_p == -1) error->all(FLERR,"Pair sph requires a fix sph/pressure");

  // set once WdeltaPinv
  WdeltaPinv = 1./SPH_KERNEL_NS::sph_kernel(kernel_id,domain->lattice->xlattice * hinv, h, hinv);
  if (WdeltaPinv < 0. || std::isnan(WdeltaPinv)) error->all(FLERR,"Wrong kernel calculation for tensile correction. Check inital lattice.");

}

/* ---------------------------------------------------------------------- */

void PairSPH::init_list(int id, NeighList *ptr)
{
  if (id == 0) list = ptr;
  else error->all(FLERR,"Internal error in PairSPH::init_list");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairSPH::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    error->all(FLERR,"PairSPH: Illegal pair_style sph command,");
  }
  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

void PairSPH::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,r,rsq,rinv,s;
  double gradWmag,fpair;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double delvx,delvy,delvz,densAB,mhu,artVisc;

  double fAB,fAB2,fAB4,rA,rB;

  double **x = atom->x;
  double **v = atom->v;
  double *q = atom->q;
  double *density = atom->density;
  double *mass = atom->mass;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    // derivative of kernel must be 0 at s = 0
    // so particle itself is not contributing

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {

        // get distance and normalized distance

        r = sqrt(rsq);
        rinv = 1./r;
        s = r * hinv;

        artVisc = 0.0;
        if (artVisc_flag) {
          // artifical viscosity [Monaghan, 1992]
          // alpha ... shear viscosity
          // beta  ... bulk viscosity
          // eta   ... avoid singularities ~ 0.01*h*h
          delvx = v[i][0] - v[j][0];
          delvy = v[i][1] - v[j][1];
          delvz = v[i][2] - v[j][2];

          if ( (delvx*delx+delvy*dely+delvz*delz) < 0.0 ) {
            mhu = h * (delvx*delx+delvy*dely+delvz*delz) / (rsq + eta);
            densAB = 0.5*(density[i]+density[j]);
            artVisc = (- alpha * cAB * mhu + beta * mhu * mhu)/densAB;
          }
        }

        rA = rB = 0.0;
        fAB = 1.0;
        if (tensCorr_flag) {
          // repulsive term for tensile instability [Monaghan, 2000]
          if (q[i] < 0.0) rA = epsilon * -1.0 * q[i] / (density[i] * density[i]);
          else rA = 0.01 * q[i] / (density[i] * density[i]);
          if (q[j] < 0.0) rB = epsilon * -1.0 * q[j] / (density[j] * density[j]);
          else rB = 0.01 * q[j] / (density[j] * density[j]);
          fAB =  SPH_KERNEL_NS::sph_kernel(kernel_id,s,h,hinv) * WdeltaPinv;
          fAB2 = fAB * fAB;
          fAB4 = fAB2 * fAB2;
        }

        // calculate value for magnitude of grad W

        gradWmag = SPH_KERNEL_NS::sph_kernel_der(kernel_id,s,h,hinv);

        // calculate the force

        fpair = - rinv * mass[itype] * mass[jtype] * (q[i]/(density[i]*density[i]) + q[j]/(density[j]*density[j]) + (rA+rB)*fAB4 + artVisc) * gradWmag; // mass[type[i]] for integration.. check fix_nve.cpp

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
/*
  for (i = 0; i < nlocal; i++) {
    fprintf(screen,"particle %d: force %f %f %f\n",i,f[i][0],f[i][1],f[i][2]);
  }*/
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSPH::write_restart(FILE *fp)
{
  write_restart_settings(fp);
/*
  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }*/
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSPH::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();
/*
  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&epsilon[i][j],sizeof(double),1,fp);
          fread(&sigma[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }*/
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSPH::write_restart_settings(FILE *fp)
{/*
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);*/
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSPH::read_restart_settings(FILE *fp)
{/*
  int me = comm->me;
  if (me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);*/
}
