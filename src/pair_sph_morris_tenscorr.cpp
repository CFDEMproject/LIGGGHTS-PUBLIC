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
    Andreas Eitzlmayr (TU Graz)

    Copyright 2009-2012 JKU Linz
    Copyright 2013-     TU Graz
------------------------------------------------------------------------- */

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_sph_morris_tenscorr.h"
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

PairSphMorrisTenscorr::PairSphMorrisTenscorr(LAMMPS *lmp) : PairSph(lmp)
{
  respa_enable = 0;
  single_enable = 0;
  pairStyle_ = 2;

  wDeltaPTypeinv = NULL;
}

/* ---------------------------------------------------------------------- */

PairSphMorrisTenscorr::~PairSphMorrisTenscorr()
{
	if (allocated) {
    if (mass_type && tensCorr_flag) memory->destroy(wDeltaPTypeinv);
	}

  if (modelStyle > 1) {
    if(fix_gamma_)
        modify->delete_fix(fix_gamma_->id);
    if(fix_visc_)
        modify->delete_fix(fix_visc_->id);
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairSphMorrisTenscorr::allocate()
{
	PairSph::allocate();

	int n = atom->ntypes;

	if (mass_type && tensCorr_flag) memory->create(wDeltaPTypeinv,n+1,n+1,"pair:wDeltaPTypeinv");

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSphMorrisTenscorr::settings(int narg, char **arg)
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

  tensCorr_flag = 0;
  modelStyle = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"newton") == 0) {
      // dynamic viscosity [Pas]
      if (iarg+2 > narg) error->all(FLERR, "Illegal pair_style sph command");
      dynVisc = force->numeric(FLERR,arg[iarg+1]);
      viscosity_ = dynVisc;
      if (modelStyle == 0) modelStyle = 1;
      else error->all(FLERR, "Illegal pair_style sph command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"power") == 0) {
      // power law parameters (K, n, etaMax, etaMin)
      if (iarg+5 > narg) error->all(FLERR, "Illegal pair_style sph command");
      conIdx = force->numeric(FLERR,arg[iarg+1]);
      powIdx = force->numeric(FLERR,arg[iarg+2]);
      etaMax = force->numeric(FLERR,arg[iarg+3]);
      etaMin = force->numeric(FLERR,arg[iarg+4]);
      OneMpowIdx = 1 - powIdx;
      if (modelStyle == 0) modelStyle = 2;
      else error->all(FLERR, "Illegal pair_style sph command");
      iarg += 5;
    } else if (strcmp(arg[iarg],"carreau") == 0) {
      // carreau parameters (eta0, etaInf, lambda, a, n)
      if (iarg+6 > narg) error->all(FLERR, "Illegal pair_style sph command");
      eta0 = force->numeric(FLERR,arg[iarg+1]);
      etaInf = force->numeric(FLERR,arg[iarg+2]);
      eta0_inf = eta0 - etaInf;
      lambda = force->numeric(FLERR,arg[iarg+3]);
      aExp = force->numeric(FLERR,arg[iarg+4]);
      powIdx = force->numeric(FLERR,arg[iarg+5]);
      OneMpowIdx_a = (1 - powIdx)/aExp;
      if (modelStyle == 0) modelStyle = 3;
      else error->all(FLERR, "Illegal pair_style sph command");
      iarg += 6;
    } else if (strcmp(arg[iarg],"tensCorr") == 0) {
      // parameters for tensile correction
      if (iarg+3 > narg) error->all(FLERR, "Illegal pair_style sph command");
      tensCorr_flag = 1;
      epsilon = force->numeric(FLERR,arg[iarg+1]);
      deltaP = force->numeric(FLERR,arg[iarg+2]);
      iarg += 3;
    } else error->all(FLERR, "Illegal pair_style sph command");
  }

  if (modelStyle > 1) {
    const char *fixarg[9];
    fixarg[0]="gamma";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="gamma";
    fixarg[4]="scalar";
    fixarg[5]="yes";
    fixarg[6]="yes";
    fixarg[7]="yes";
    fixarg[8]="0";
    fix_gamma_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),"PairSphMorrisTenscorr");

    fixarg[0]="viscosity";
    fixarg[3]="viscosity";
    fix_visc_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),"PairSphMorrisTenscorr");
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairSphMorrisTenscorr::coeff(int narg, char **arg)
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

void PairSphMorrisTenscorr::init_substyle()
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

  if (modelStyle > 1) {

    // check if there is a fix sph/velgrad present
    int ifix_velgrad = -1;
    for(int i = 0; i < modify->nfix; i++)
    {
      if(strncmp("sph/velgrad",modify->fix[i]->style,11) == 0) {
        ifix_velgrad = i;
      }
    }
    if(ifix_velgrad == -1) error->all(FLERR,"Requires to define a fix sph/velgrad also\n");

    fix_dvdx_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("dvdx","property/atom","vector",0,0,"PairSphMorrisTenscorr",false));
    fix_dvdy_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("dvdy","property/atom","vector",0,0,"PairSphMorrisTenscorr",false));
    fix_dvdz_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("dvdz","property/atom","vector",0,0,"PairSphMorrisTenscorr",false));
  }

}

/* ----------------------------------------------------------------------
  allocate per-type and per-type pair properties
------------------------------------------------------------------------- */
/*
void PairSphMorrisTenscorr::allocate_properties(int size)
{

}
*/
/* ---------------------------------------------------------------------- */

void PairSphMorrisTenscorr::compute(int eflag, int vflag)
{
  if (mass_type) compute_eval<1>(eflag,vflag);
  else compute_eval<0>(eflag,vflag);
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSphMorrisTenscorr::write_restart(FILE *fp)
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

void PairSphMorrisTenscorr::read_restart(FILE *fp, const int major, const int minor)
{
  read_restart_settings(fp, major, minor);
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

void PairSphMorrisTenscorr::write_restart_settings(FILE *fp)
{
  fwrite(&tensCorr_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSphMorrisTenscorr::read_restart_settings(FILE *fp, const int major, const int minor)
{
  if (comm->me == 0) {
    fread(&tensCorr_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&tensCorr_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   template compute
------------------------------------------------------------------------- */

template <int MASSFLAG>
void PairSphMorrisTenscorr::compute_eval(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double xtmp,ytmp,ztmp,delx,dely,delz,r,rsq,rinv,s;
  double gradWmag,prefac,fgradP,fgradPX,fgradPY,fgradPZ,fpairX,fpairY,fpairZ;
  double D11,D22,D33,D12,D13,D23,gamma2,gamma;
  double A, B, C;

  double rhoi,rhoj,pi,pj,etai,etaj;
  double sli,slj,slCom,slComInv,imass,jmass;

  double morrisVisc,fAB4,rAB;
  double rA,rB,fAB,fAB2;
  double wDeltaPinv;

  double radi,radj,rcom;

  double **x = atom->x;
  double **v = atom->vest;
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

  fgradP_ = fix_fgradP_->array_atom;

  if (modelStyle > 1) {
    dvdx_ = fix_dvdx_->array_atom;
    dvdy_ = fix_dvdy_->array_atom;
    dvdz_ = fix_dvdz_->array_atom;
    gamma_ = fix_gamma_->vector_atom;
    visc_ = fix_visc_->vector_atom;
  }

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

  // Reset fgradP
  for (i = 0; i < nlocal; i++) {
    fgradP_[i][0] = 0;
    fgradP_[i][1] = 0;
    fgradP_[i][2] = 0;

    if (modelStyle > 1) {
      D11 = dvdx_[i][0];
      D22 = dvdy_[i][1];
      D33 = dvdz_[i][2];
      D12 = 0.5 * (dvdy_[i][0] + dvdx_[i][1]);
      D13 = 0.5 * (dvdz_[i][0] + dvdx_[i][2]);
      D23 = 0.5 * (dvdz_[i][1] + dvdy_[i][2]);

      // gamma = sqrt(2*D:D) = sqrt(2*Frobeniusnorm(D)), similarly in e.g.,
      // A. Ficarella, M. Milanese, D. Laforgia, J. Food Eng. 72 (2006) 179-188.
      // S. Shao, E.Y.M. Lo, Adv. Water Res. 26 (2003) 787-800.
      // H.A. Ardakani, E. Mitsoulis, S.G. Hatzikiriakos, J. Non-Newt. Fluid Mech. 166 (2011) 1261-1271.

      gamma2 = 2*(D11*D11 + D22*D22 + D33*D33 + 2*D12*D12 + 2*D13*D13 + 2*D23*D23);
      gamma_[i] = sqrt(gamma2);
    }

    if (modelStyle == 2) { // power law
      gamma = gamma_[i];
      if (gamma <= 0) gamma = 1e-10;
      visc_[i] = conIdx / pow(gamma,OneMpowIdx);
      if (visc_[i] > etaMax) visc_[i] = etaMax;
      if (visc_[i] < etaMin) visc_[i] = etaMin;
    }

    if (modelStyle == 3) { // carreau
      A = gamma_[i]*lambda;
      B = 1 + pow(A,aExp);
      C = pow(B,OneMpowIdx_a);
      visc_[i] = etaInf + eta0_inf / C;
    }
  }

  if (modelStyle > 1) {
    // send to ghosts
    fix_visc_->do_forward_comm();
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
        if (r == 0.) {
          printf("Particle %i and %i are at same position (%f, %f, %f)\n",i,j,xtmp,ytmp,ztmp);
          error->one(FLERR,"Zero distance between SPH particles!");
        }
        rinv = 1./r;
        s = r * slComInv;

        // calculate value for magnitude of grad W
        gradWmag = SPH_KERNEL_NS::sph_kernel_der(kernel_id,s,slCom,slComInv);

        // viscosity
        if (modelStyle == 1) {// Newtonian
          etai = dynVisc;
          etaj = dynVisc;
        } else { // non-Newtonian
          etai = visc_[i];
          etaj = visc_[j];
        }

				// Morris viscosity term
				morrisVisc = imass * jmass * (etai + etaj) * rinv * gradWmag / (rhoi * rhoj);

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
        }

        // calculate the force
        prefac = - rinv * imass * jmass * gradWmag;
        fgradP = prefac * (pi/(rhoi*rhoi) + pj/(rhoj*rhoj) + rAB*fAB4);
						// mass[i] for integration.. check fix_nve.cpp

        fgradPX = fgradP*delx;
        fgradPY = fgradP*dely;
        fgradPZ = fgradP*delz;

				fpairX = fgradPX + morrisVisc * (v[i][0]-v[j][0]);
				fpairY = fgradPY + morrisVisc * (v[i][1]-v[j][1]);
				fpairZ = fgradPZ + morrisVisc * (v[i][2]-v[j][2]);

        // apply the force

        fgradP_[i][0] += fgradPX;
        fgradP_[i][1] += fgradPY;
        fgradP_[i][2] += fgradPZ;

        f[i][0] += fpairX;
        f[i][1] += fpairY;
        f[i][2] += fpairZ;

        if (newton_pair || j < nlocal) {
          fgradP_[j][0] -= fgradPX;
          fgradP_[j][1] -= fgradPY;
          fgradP_[j][2] -= fgradPZ;

          f[j][0] -= fpairX;
          f[j][1] -= fpairY;
          f[j][2] -= fpairZ;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,0.0,0.0,fgradP,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();

}
