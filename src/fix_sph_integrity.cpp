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
    Copyrigth 2013-     TU Graz
------------------------------------------------------------------------- */

#include <cmath>
#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include "fix_sph_integrity.h"
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

FixSphIntegrity::FixSphIntegrity(LAMMPS *lmp, int narg, char **arg) :
  FixSph(lmp, narg, arg)
{
  int iarg = 0;

  if (iarg+3 > narg) error->fix_error(FLERR,this,"Not enough input arguments");

  iarg += 3;

  every = 1;

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

  fix_integrity_ = NULL;
}

/* ---------------------------------------------------------------------- */

FixSphIntegrity::~FixSphIntegrity()
{

}

/* ---------------------------------------------------------------------- */

void FixSphIntegrity::post_create()
{
  const char *fixarg[9];
  fixarg[0]="int";
  fixarg[1]="all";
  fixarg[2]="property/atom";
  fixarg[3]="int";
  fixarg[4]="scalar";
  fixarg[5]="yes";
  fixarg[6]="yes";
  fixarg[7]="yes";
  fixarg[8]="1";
  fix_integrity_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
}

/* ---------------------------------------------------------------------- */

void FixSphIntegrity::pre_delete(bool unfixflag)
{
    if(unfixflag && fix_integrity_)
        modify->delete_fix(fix_integrity_->id);
}

/* ---------------------------------------------------------------------- */

int FixSphIntegrity::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSphIntegrity::init()
{
  FixSph::init();

  // check if there is an sph/pressure fix present
  // must come before me, because
  int dens = -1;
  int me = -1;
  for(int i = 0; i < modify->nfix; i++)
  {
    if(strcmp("sph/integrity",modify->fix[i]->style)) {
      me = i;
    }
    if(strncmp("sph/density",modify->fix[i]->style,11) == 0) {
      dens = i;
    }
  }

  if(dens > me) error->fix_error(FLERR,this,"Fix sph/density has to be defined before sph/integrity \n");
  if(dens == -1) error->fix_error(FLERR,this,"Requires to define a fix sph/density also \n");
}

/* ---------------------------------------------------------------------- */

void FixSphIntegrity::post_integrate()
{
  //template function for using per atom or per atomtype smoothing length
  if (mass_type) post_integrate_eval<1>();
  else post_integrate_eval<0>();
}

/* ---------------------------------------------------------------------- */

template <int MASSFLAG>
void FixSphIntegrity::post_integrate_eval()
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,r,s=0.0,W;
  double sli,sliInv,slj,slCom,slComInv,cut,imass,jmass,irho,jrho;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
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
      integrity_ = fix_integrity_->vector_atom;

      //int nlocal = atom->nlocal;
      for (i = 0; i < nlocal; i++) {
        integrity_[i] += 0;
      }
    }
  }

  if (ago % every == 0) {
    ago = 0;

  updatePtrs(); // get sl

  integrity_ = fix_integrity_->vector_atom;

  // reset and add rho contribution of self
  for (i = 0; i < nlocal; i++) {
    if (MASSFLAG) {
      itype = type[i];
      sli = sl[itype-1];
      imass = mass[itype];
    } else {
      sli = sl[i];
      imass = rmass[i];
    }

    sliInv = 1./sli;
    irho = rho[i];

    // this gets a value for W at self, perform error check

    W = SPH_KERNEL_NS::sph_kernel(kernel_id,0.,sli,sliInv);
    if (W < 0.)
    {
      fprintf(screen,"s = %f, W = %f\n",s,W);
      error->one(FLERR,"Illegal kernel used, W < 0");
    }

    // add contribution of self
    integrity_[i] = W * imass / irho;
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

      // this gets a value for W at self, perform error check

      W = SPH_KERNEL_NS::sph_kernel(kernel_id,s,slCom,slComInv);
      if (W < 0.)
      {
        fprintf(screen,"s = %f, W = %f\n",s,W);
        error->one(FLERR,"Illegal kernel used, W < 0");
      }

      // add contribution of neighbor
      // have a half neigh list, so do it for both if necessary

      // What we called "integrity" here, is the sum of (mass/rho * W),
      // sometimes called "the kernel estimate of unity" (e.g., P.W. Randles, L.D. Libersky,

      // add contribution of neighbor
      // have a half neigh list, so do it for both if necessary

      integrity_[i] += W * jmass / jrho;

      if (newton_pair || j < nlocal)
        integrity_[j] += W * imass / irho;
    }
  }

    // rho is now correct, send to ghosts
    timer->stamp();
    comm->forward_comm();
    timer->stamp(TIME_COMM);
  }
}
