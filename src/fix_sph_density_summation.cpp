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
#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include "fix_sph_density_summation.h"
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

FixSPHDensitySum::FixSPHDensitySum(LAMMPS *lmp, int narg, char **arg) :
  FixSph(lmp, narg, arg)
{
  int iarg = 0;

  if (iarg+3 > narg) error->fix_error(FLERR,this,"Not enough input arguments");

  iarg += 3;

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

    } else error->fix_error(FLERR,this,"Wrong keyword.");
  }
}

/* ---------------------------------------------------------------------- */

FixSPHDensitySum::~FixSPHDensitySum()
{

}

/* ---------------------------------------------------------------------- */

int FixSPHDensitySum::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSPHDensitySum::init()
{
  FixSph::init();

  // check if there is an sph/pressure fix present
  // must come before me, because
  // a - it needs the rho for the pressure calculation
  // b - i have to do the forward comm to have updated ghost properties

  int pres = -1;
  int me = -1;
  for(int i = 0; i < modify->nfix; i++)
  {
    if(strcmp("sph/density/summation",modify->fix[i]->style)) {
      me = i;
    }
    if(strncmp("sph/pressure",modify->fix[i]->style,12) == 0) {
      pres = i;
      break;
    }
  }

  if(me == -1 && pres >= 0) error->fix_error(FLERR,this,"Fix sph/pressure has to be defined after sph/density/summation \n");
  if(pres == -1) error->fix_error(FLERR,this,"Requires to define a fix sph/pressure also \n");
}

/* ---------------------------------------------------------------------- */

void FixSPHDensitySum::post_integrate()
{
  //template function for using per atom or per atomtype smoothing length
  if (mass_type) post_integrate_eval<1>();
  else post_integrate_eval<0>();
}

/* ---------------------------------------------------------------------- */

template <int MASSFLAG>
void FixSPHDensitySum::post_integrate_eval()
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,r,s=0.0,W;
  double sli,sliInv,slj,slCom,slComInv,cut,imass,jmass;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  int *mask = atom->mask;
  double *rho = atom->rho;
  int newton_pair = force->newton_pair;

  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;

  updatePtrs(); // get sl

  // reset and add rho contribution of self

  int nlocal = atom->nlocal;
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

    // this gets a value for W at self, perform error check

    W = SPH_KERNEL_NS::sph_kernel(kernel_id,0.,sli,sliInv);
    if (W < 0.)
    {
      fprintf(screen,"s = %f, W = %f\n",s,W);
      error->one(FLERR,"Illegal kernel used, W < 0");
    }

    // add contribution of self
    rho[i] = W * imass;
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

      rho[i] += W * jmass;

      if (newton_pair || j < nlocal)
        rho[j] += W * imass;
    }
  }

  // rho is now correct, send to ghosts
  timer->stamp();
  comm->forward_comm();
  timer->stamp(TIME_COMM);

}
