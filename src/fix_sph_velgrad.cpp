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

    Andreas Eitzlmayr (TU Graz)

    Copyright 2013-     TU Graz
------------------------------------------------------------------------- */

#include <cmath>
#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include "fix_sph_velgrad.h"
#include "update.h"
#include "respa.h"
#include "atom.h"
#include "force.h"
#include "modify.h"
#include "pair.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "sph_kernels.h"
#include "fix_property_atom.h"
#include "timer.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSphVelgrad::FixSphVelgrad(LAMMPS *lmp, int narg, char **arg) :
  FixSph(lmp, narg, arg)
{
  int iarg = 0;

  if (iarg+3 > narg) error->fix_error(FLERR,this,"Not enough input arguments");

  iarg += 3;

  every = 1; // optional argument, default value = 1
  ago = 0;
  velgrad_flag = 1;

  while (iarg < narg) {
    // kernel style
    if (strcmp(arg[iarg],"sphkernel") == 0) {
          if (iarg+2 > narg) error->fix_error(FLERR,this,"Illegal use of keyword 'sphkernel'. Not enough input arguments");

          if(kernel_style) delete []kernel_style;
          kernel_style = new char[strlen(arg[iarg+1])+1];
          strcpy(kernel_style,arg[iarg+1]);

          // check uniqueness of kernel IDs

          int flag = SPH_KERNEL_NS::sph_kernels_unique_id();
          if(flag < 0) error->fix_error(FLERR,this,"Cannot proceed, sph kernels need unique IDs, check all sph_kernel_* files");

          // get kernel id

          kernel_id = SPH_KERNEL_NS::sph_kernel_id(kernel_style);
          if(kernel_id < 0) error->fix_error(FLERR,this,"Unknown sph kernel");

          iarg += 2;

    } else if (strcmp(arg[iarg],"every") == 0) {
      every = force->inumeric(FLERR,arg[iarg+1]);
      if (every <= 0) error->fix_error(FLERR,this,"every <= 0 not allowed");
      iarg += 2;

    } else error->fix_error(FLERR,this,"Wrong keyword.");
  }

  fix_dvdx_ = NULL;
  fix_dvdy_ = NULL;
  fix_dvdz_ = NULL;
}

/* ---------------------------------------------------------------------- */

FixSphVelgrad::~FixSphVelgrad()
{

}

/* ---------------------------------------------------------------------- */

void FixSphVelgrad::post_create()
{
  const char *fixarg[11];
  fixarg[0]="dvdx";
  fixarg[1]="all";
  fixarg[2]="property/atom";
  fixarg[3]="dvdx";
  fixarg[4]="vector";
  fixarg[5]="yes";
  fixarg[6]="yes";
  fixarg[7]="yes";
  fixarg[8]="0";
  fixarg[9]="0";
  fixarg[10]="0";
  fix_dvdx_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);

  fixarg[0]="dvdy";
  fixarg[3]="dvdy";
  fix_dvdy_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);

  fixarg[0]="dvdz";
  fixarg[3]="dvdz";
  fix_dvdz_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
}

/* ---------------------------------------------------------------------- */

void FixSphVelgrad::pre_delete(bool unfixflag)
{
    if(unfixflag && fix_dvdx_)
        modify->delete_fix(fix_dvdx_->id);

    if(unfixflag && fix_dvdy_)
        modify->delete_fix(fix_dvdy_->id);

    if(unfixflag && fix_dvdz_)
        modify->delete_fix(fix_dvdz_->id);
}

/* ---------------------------------------------------------------------- */

int FixSphVelgrad::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSphVelgrad::init()
{
  FixSph::init();

  // check if there is an sph/density fix present

  int dens = -1;

  for(int i = 0; i < modify->nfix; i++)
  {
    if(strncmp("sph/density",modify->fix[i]->style,11) == 0) {
      dens = i;
    }
  }

  if(dens == -1) error->fix_error(FLERR,this,"Requires to define a fix sph/density also \n");
}

/* ---------------------------------------------------------------------- */

void FixSphVelgrad::pre_force(int vflag)
{
  //template function for using per atom or per atomtype smoothing length
  if (mass_type) pre_force_eval<1>(vflag);
  else pre_force_eval<0>(vflag);
}

/* ---------------------------------------------------------------------- */

template <int MASSFLAG>
void FixSphVelgrad::pre_force_eval(int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,r,s,gradWmag,m_rhoGradWmag_r,delvx,delvy,delvz;
  double sli,slj,slCom,slComInv,cut,imass,jmass,irho,jrho;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  double **v = atom->vest;
  int *mask = atom->mask;
  double *rho = atom->rho;
  int newton_pair = force->newton_pair;

  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  ago++;

  if (every > 1) {
    if (ago == 1) {
      dvdx_ = fix_dvdx_->array_atom;
      dvdy_ = fix_dvdy_->array_atom;
      dvdz_ = fix_dvdz_->array_atom;

      // reset values
      velgrad_flag = 0;

      //int nlocal = atom->nlocal;
      for (i = 0; i < nlocal; i++) {
        dvdx_[i][0] = 0;
        dvdx_[i][1] = 0;
        dvdx_[i][2] = 0;
        dvdy_[i][0] = 0;
        dvdy_[i][1] = 0;
        dvdy_[i][2] = 0;
        dvdz_[i][0] = 0;
        dvdz_[i][1] = 0;
        dvdz_[i][2] = 0;
      }
    }
  }

  if (ago % every == 0) {
    ago = 0;
    velgrad_flag = 1;

    updatePtrs(); // get sl

    dvdx_ = fix_dvdx_->array_atom;
    dvdy_ = fix_dvdy_->array_atom;
    dvdz_ = fix_dvdz_->array_atom;

    // reset values

    for (i = 0; i < nlocal; i++) {
      dvdx_[i][0] = 0;
      dvdx_[i][1] = 0;
      dvdx_[i][2] = 0;
      dvdy_[i][0] = 0;
      dvdy_[i][1] = 0;
      dvdy_[i][2] = 0;
      dvdz_[i][0] = 0;
      dvdz_[i][1] = 0;
      dvdz_[i][2] = 0;
    }

    // need updated ghost positions and self contributions
    timer->stamp();
    comm->forward_comm();
    timer->stamp(TIME_COMM);

    // loop over neighbors of my atoms

    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];

      if (!(mask[i] & groupbit)) continue;
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      if (MASSFLAG) {
        itype = type[i];
        imass = mass[itype];
      } else {
        imass = rmass[i];
        sli = sl[i];
      }
      irho = rho[i];

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];

        if (!(mask[j] & groupbit)) continue;

        if (MASSFLAG) {
          jtype = type[j];
          jmass = mass[jtype];
          slCom = slComType[itype][jtype];
        } else {
          jmass = rmass[j];
          slj = sl[j];
          slCom = interpDist(sli,slj);
        }

        jrho = rho[j];
        slComInv = 1./slCom;
        cut = slCom*kernel_cut;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;

        if (rsq >= cut*cut) continue;
        // calculate distance and normalized distance

        r = sqrt(rsq);
        slComInv = 1./slCom;
        s = r*slComInv;

        // calculate value for magnitude of grad W
        gradWmag = SPH_KERNEL_NS::sph_kernel_der(kernel_id,s,slCom,slComInv);

        m_rhoGradWmag_r = jmass / jrho * gradWmag / r;

        delvx = v[j][0] - v[i][0];
        delvy = v[j][1] - v[i][1];
        delvz = v[j][2] - v[i][2];

        // dvdx_[i][0] is the x-derivative of the x-velocity component (dvx/dx)
        // dvdx_[i][1] is the x-derivative of the y-velocity component (dvy/dx)
        // analogous for the others ...

        // calculation of the velocity derivatives as shown e.g., in:
        // Y. Amini, H. Emdad, M. Farid, Europ. J. Mech. B/Fluids 30 (2011) 184-194.
        // J. Fang, A. Parriaux, M. Rentschler, C. Ancey, Appl. Num. Math. 59 (2009) 251-271.
        // M. Ellero, M. Kröger, S. Hess, J. Non-Newt. Fluid Mech. 105 (2002) 35-51.

        dvdx_[i][0] += m_rhoGradWmag_r * delvx * delx;
        dvdx_[i][1] += m_rhoGradWmag_r * delvy * delx;
        dvdx_[i][2] += m_rhoGradWmag_r * delvz * delx;
        dvdy_[i][0] += m_rhoGradWmag_r * delvx * dely;
        dvdy_[i][1] += m_rhoGradWmag_r * delvy * dely;
        dvdy_[i][2] += m_rhoGradWmag_r * delvz * dely;
        dvdz_[i][0] += m_rhoGradWmag_r * delvx * delz;
        dvdz_[i][1] += m_rhoGradWmag_r * delvy * delz;
        dvdz_[i][2] += m_rhoGradWmag_r * delvz * delz;

        if (newton_pair || j < nlocal)
          m_rhoGradWmag_r = imass / irho * gradWmag / r;
        dvdx_[j][0] += m_rhoGradWmag_r * delvx * delx;
        dvdx_[j][1] += m_rhoGradWmag_r * delvy * delx;
        dvdx_[j][2] += m_rhoGradWmag_r * delvz * delx;
        dvdy_[j][0] += m_rhoGradWmag_r * delvx * dely;
        dvdy_[j][1] += m_rhoGradWmag_r * delvy * dely;
        dvdy_[j][2] += m_rhoGradWmag_r * delvz * dely;
        dvdz_[j][0] += m_rhoGradWmag_r * delvx * delz;
        dvdz_[j][1] += m_rhoGradWmag_r * delvy * delz;
        dvdz_[j][2] += m_rhoGradWmag_r * delvz * delz;
      }
    }

    // dvdx, dvdy, dvdz are now correct, send to ghosts
    timer->stamp();
    comm->forward_comm();
    timer->stamp(TIME_COMM);
  }
}
