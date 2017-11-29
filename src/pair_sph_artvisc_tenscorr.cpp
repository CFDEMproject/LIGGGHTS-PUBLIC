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

    Andreas Aigner (JKU Linz)

    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

PairSphArtviscTenscorr::PairSphArtviscTenscorr(LAMMPS *lmp) : PairSph(lmp),
    artVisc_flag(false),
    tensCorr_flag(false),
    cs(NULL),
    alpha(NULL),
    beta(NULL),
    etaPPG(NULL),
    csmean(NULL),
    alphaMean(NULL),
    betaMean(NULL),
    eta(0.),
    epsilonPPG(NULL),
    deltaP(NULL),
    wDeltaPTypeinv(NULL),
    epsilon(0.)
{
  respa_enable = 0;
  single_enable = 0;
  pairStyle_ = 1;
}

/* ---------------------------------------------------------------------- */

PairSphArtviscTenscorr::~PairSphArtviscTenscorr()
{
  if (allocated) {
    memory->destroy(csmean);
    memory->destroy(alphaMean);
    memory->destroy(betaMean);
    memory->destroy(wDeltaPTypeinv);
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairSphArtviscTenscorr::allocate()
{
  PairSph::allocate();

  int n = atom->ntypes;

  if (artVisc_flag) {
    memory->create(csmean,n+1,n+1,"pair:csmean");
    memory->create(alphaMean,n+1,n+1,"pair:alphaMean");
    memory->create(betaMean,n+1,n+1,"pair:betaMean");
  }

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

  iarg += 2;

  // optional parameters

  artVisc_flag = tensCorr_flag = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"artVisc") == 0) {
      // parameters for artifical viscosity
      if (iarg+1 > narg) error->all(FLERR, "Illegal pair_style sph command");
      artVisc_flag = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg],"tensCorr") == 0) {
      // parameters for tensile correction
      if (iarg+1 > narg) error->all(FLERR, "Illegal pair_style sph command");
      tensCorr_flag = 1;
      iarg += 1;
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
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");

}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSphArtviscTenscorr::init_substyle()
{
  const int max_type = atom->ntypes;

  //create wDeltaPTypeInv
  if (mass_type && tensCorr_flag) {

    deltaP=static_cast<FixPropertyGlobal*>(modify->find_fix_property("tensCorrDeltaP","property/global","peratomtype",max_type,0,force->pair_style));
    if(!deltaP) error->all(FLERR, "Pairstyle sph/artVisc/tensCorr only works with a fix property/global that defines tensCorrDeltaP");

    epsilonPPG=static_cast<FixPropertyGlobal*>(modify->find_fix_property("tensCorrEpsilon","property/global","scalar",0,0,force->pair_style));
    if(!epsilonPPG) error->all(FLERR, "Pairstyle sph/artVisc/tensCorr only works with a fix property/global that defines tensCorrEpsilon");
    epsilon = epsilonPPG->compute_scalar(); // const for all types

    //pre-calculate common smoothing length
    for(int i = 1; i < max_type+1; i++) {
      for(int j = 1; j < max_type+1; j++) {
        const double deltaPi = deltaP->compute_vector(i-1);
        const double deltaPj = deltaP->compute_vector(j-1);
        const double meanDeltaP = 0.5*(deltaPi+deltaPj);

        const double slCom = slComType[i][j];
        const double slComInv = 1./slCom;
        wDeltaPTypeinv[i][j] = 1./SPH_KERNEL_NS::sph_kernel(kernel_id,meanDeltaP * slComInv,slCom,slComInv);
      }
    }
  }

  //Get pointer to the fixes that have the material properties

  if (artVisc_flag) {
    cs=static_cast<FixPropertyGlobal*>(modify->find_fix_property("speedOfSound","property/global","peratomtype",max_type,0,force->pair_style));
    if(!cs) error->all(FLERR, "Pairstyle sph/artVisc/tensCorr only works with a fix property/global that defines speedOfSound");

    alpha=static_cast<FixPropertyGlobal*>(modify->find_fix_property("artViscAlpha","property/global","peratomtype",max_type,0,force->pair_style));
    if(!alpha) error->all(FLERR, "Pairstyle sph/artVisc/tensCorr only works with a fix property/global that defines artViscAlpha");

    beta=static_cast<FixPropertyGlobal*>(modify->find_fix_property("artViscBeta","property/global","peratomtype",max_type,0,force->pair_style));
    if(!beta) error->all(FLERR, "Pairstyle sph/artVisc/tensCorr only works with a fix property/global that defines artViscBeta");

    etaPPG=static_cast<FixPropertyGlobal*>(modify->find_fix_property("artViscEta","property/global","scalar",0,0,force->pair_style));
    if(!etaPPG) error->all(FLERR, "Pairstyle sph/artVisc/tensCorr only works with a fix property/global that defines artViscEta");
    eta = etaPPG->compute_scalar(); 

    viscosity_ = 1; 
    //viscosity_ = alpha;

    //pre-calculate parameters for possible contact material combinations
    for(int i=1;i< max_type+1; i++)
    {
      for(int j=1;j<max_type+1;j++)
      {
        const double csi=cs->compute_vector(i-1);
        const double csj=cs->compute_vector(j-1);

        const double alphai=alpha->compute_vector(i-1);
        const double alphaj=alpha->compute_vector(j-1);

        const double betai=beta->compute_vector(i-1);
        const double betaj=beta->compute_vector(j-1);

        csmean[i][j] = 0.5*(csi+csj);
        alphaMean[i][j] = 0.5*(alphai+alphaj);
        betaMean[i][j] = 0.5*(betai+betaj);

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

void PairSphArtviscTenscorr::read_restart(FILE *fp, const int major, const int minor)
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

void PairSphArtviscTenscorr::read_restart_settings(FILE *fp, const int major, const int minor)
{
  if (comm->me == 0) {
    fread(&artVisc_flag,sizeof(int),1,fp);
    fread(&tensCorr_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&artVisc_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tensCorr_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   template compute
------------------------------------------------------------------------- */

template <int MASSFLAG>
void PairSphArtviscTenscorr::compute_eval(int eflag, int vflag)
{
  double sli,slCom,imass,jmass;
  double artVisc,fAB4,rAB;
  double rA,rB;
  double wDeltaPinv;

  double radi,rcom;

  double **x = atom->x;
  double **v = atom->vest;
  double *p = atom->p;
  double *rho = atom->rho;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  const int newton_pair = force->newton_pair;

  // individual properties
  double *mass = atom->mass;
  double *radius = atom->radius;
  double *rmass = atom->rmass;

  // TODO: Use this eflag/vflag .. what is it for?
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  const int inum = list->inum;
  int * const ilist = list->ilist;
  int * const numneigh = list->numneigh;
  int ** const firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  // depend on mass_type

  // no communication of sl, in case of perAtomType mass
  if (!MASSFLAG) {
    // update smoothing length
    fppaSl->do_forward_comm();
    updatePtrs(); // get sl
  }

  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    const int itype = type[i];
    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];
    int * const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    const double rhoi = rho[i];
    const double pi = p[i];

    if (MASSFLAG) {
      imass = mass[itype];
    } else {
      sli = sl[i];
      radi = radius[i];
      imass = rmass[i];
    }

    // derivative of kernel must be 0 at s = 0
    // so particle itself is not contributing

    for (int jj = 0; jj < jnum; jj++) {
      const int j = jlist[jj];
      const int jtype = type[j];

      const double delx = xtmp - x[j][0];
      const double dely = ytmp - x[j][1];
      const double delz = ztmp - x[j][2];
      const double rsq = delx*delx + dely*dely + delz*delz;

      if (!MASSFLAG) {
        const double radj = radius[j];
        rcom = interpDist(radi,radj);
      }

      if ((MASSFLAG && rsq < cutsq[itype][jtype]) || (!MASSFLAG && rsq < rcom*rcom)) {

        if (MASSFLAG) {
          jmass = mass[jtype];
          slCom = slComType[itype][jtype];
        } else {
          jmass = rmass[j];
          const double slj = sl[j];
          slCom = interpDist(sli,slj);
        }

        const double pj = p[j];
        const double rhoj = rho[j];
        const double slComInv = 1./slCom;
        //cut = slCom*SPH_KERNEL_NS::sph_kernel_cut(kernel_id);

        // get distance and normalized distance
        const double r = sqrt(rsq);
        if (r == 0.) {
          printf("Particle %i and %i are at same position (%f, %f, %f)",i,j,xtmp,ytmp,ztmp);
          error->one(FLERR,"Zero distance between SPH particles!");
        }
        const double rinv = 1./r;
        const double s = r * slComInv;

        // calculate value for magnitude of grad W
        const double gradWmag = SPH_KERNEL_NS::sph_kernel_der(kernel_id,s,slCom,slComInv);

        // artificial viscosity
        artVisc = 0.0;
        if (artVisc_flag) {
          // artifical viscosity [Monaghan, 1992]
          // alpha ... shear viscosity
          // beta  ... bulk viscosity
          // eta   ... avoid singularities ~ 0.01*h*h

          const double dotDelVDelR = ( (v[i][0]-v[j][0])*delx + (v[i][1]-v[j][1])*dely + (v[i][2]-v[j][2])*delz );

          if ( dotDelVDelR < 0.0 ) { 
            const double muAB = slCom * dotDelVDelR / (rsq + eta);
            const double rhoMeanInv = 2/(rhoi+rhoj);
            artVisc = ((- alphaMean[itype][jtype] * csmean[itype][jtype] * muAB + betaMean[itype][jtype] * muAB * muAB) * rhoMeanInv);
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
            // assumption that deltaP = sl / 1.2
            const double deltaPOne = slCom/1.2;
            wDeltaPinv = 1./SPH_KERNEL_NS::sph_kernel(kernel_id,deltaPOne * slComInv,slCom,slComInv);
          }

          //TODO: Is fAB4 in this form ok?!
          const double fAB =  SPH_KERNEL_NS::sph_kernel(kernel_id,s,slCom,slComInv) * wDeltaPinv;
          const double fAB2 = fAB * fAB;
          fAB4 = fAB2 * fAB2;
        }

        // calculate the force
        const double fpair = - rinv * imass * jmass * (pi/(rhoi*rhoi) + pj/(rhoj*rhoj) + rAB*fAB4 + artVisc) * gradWmag; // mass[i] for integration.. check fix_nve.cpp

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
