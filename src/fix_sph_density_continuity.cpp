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
    Andreas Aigner (JKU Linz))

    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#include <cmath>
#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include "fix_sph_density_continuity.h"
#include "update.h"
#include "atom.h"
#include "force.h"
#include "modify.h"
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

FixSphDensityContinuity::FixSphDensityContinuity(LAMMPS *lmp, int narg, char **arg) :
  FixSph(lmp, narg, arg)
{
  int iarg = 0;

  if (iarg+3 > narg) error->all(FLERR,"Illegal fix sph/density/continuity command");

  iarg += 3;

  while (iarg < narg) {
    // kernel style
    if (strcmp(arg[iarg],"sphkernel") == 0) {
          if (iarg+2 > narg) error->all(FLERR,"Illegal fix sph/density/continuity command");

          if(kernel_style) delete []kernel_style;
          kernel_style = new char[strlen(arg[iarg+1])+1];
          strcpy(kernel_style,arg[iarg+1]);

          // check uniqueness of kernel IDs

          int flag = SPH_KERNEL_NS::sph_kernels_unique_id();
          if(flag < 0) error->all(FLERR,"Cannot proceed, sph kernels need unique IDs, check all sph_kernel_* files");

          // get kernel id

          kernel_id = SPH_KERNEL_NS::sph_kernel_id(kernel_style);
          if(kernel_id < 0) error->all(FLERR,"Illegal fix sph/density/continuity command, unknown sph kernel");

          iarg += 2;

    } else error->all(FLERR,"Illegal fix sph/density/continuity command");
  }

  time_depend = 0; // only time step is used, but not a relative time

}

/* ---------------------------------------------------------------------- */

FixSphDensityContinuity::~FixSphDensityContinuity()
{

}

/* ---------------------------------------------------------------------- */

int FixSphDensityContinuity::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSphDensityContinuity::init()
{
  FixSph::init();

  // check if there is an nve/sph fix present

  int idx_integ = -1;
  for(int i = 0; i < modify->nfix; i++)
  {
    if(strncmp("nve/sph",modify->fix[i]->style,7) == 0) {
      idx_integ = i;
      break;
    }
    if(strncmp("nve/xsph",modify->fix[i]->style,8) == 0) {
      idx_integ = i;
      break;
    }
  }

  if(idx_integ == -1) error->fix_error(FLERR,this,"Requires to define a fix nve/sph also \n");
}

/* ---------------------------------------------------------------------- */

void FixSphDensityContinuity::pre_force(int vflag)
{
  //template function for using per atom or per atomtype smoothing length
  if (mass_type) pre_force_eval<1>(vflag);
  else pre_force_eval<0>(vflag);
}

/* ---------------------------------------------------------------------- */

template <int MASSFLAG>
void FixSphDensityContinuity::pre_force_eval(int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,r,rinv,s,gradWmag;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double sli,slj,slCom,slComInv,cut,delVDotDelR,imass,jmass;

  double **x = atom->x;
  double **v = atom->vest;
  double *drho = atom->drho;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int *type = atom->type;       // if MASSFLAG
  double *mass = atom->mass;    // if MASSFLAG
  double *rmass = atom->rmass;  // if !MASSFLAG

  int newton_pair = force->newton_pair;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // need updated ghost positions and self contributions
  timer->stamp();
  
  if (!MASSFLAG) fppaSl->do_forward_comm(); 
  timer->stamp(TIME_COMM);

  if (!MASSFLAG) updatePtrs(); // get sl

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

      cut = slCom*kernel_cut;//SPH_KERNEL_NS::sph_kernel_cut(kernel_id);

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      // TODO: Cutsq from pair?
      if (rsq >= cut*cut) continue;

      // calculate normalized distance

      r = sqrt(rsq);
      if (r == 0.) {
        fprintf(screen,"Particle %i and %i are at same position (%f, %f, %f)\n",i,j,xtmp,ytmp,ztmp);
        error->one(FLERR,"Zero distance between SPH particles!");
      }
      rinv = 1./r;
      slComInv = 1./slCom;
      s = r * slComInv;

      //    scalar product of delV and delR/R
      delVDotDelR = rinv * ( delx*(v[i][0]-v[j][0]) + dely*(v[i][1]-v[j][1]) + delz*(v[i][2]-v[j][2]) );

      // calculate value for magnitude of grad W
      gradWmag = SPH_KERNEL_NS::sph_kernel_der(kernel_id,s,slCom,slComInv);

      // add contribution of neighbor
      // have a half neigh list, so do it for both if necessary

      drho[i] += jmass*gradWmag*delVDotDelR;

      if (newton_pair || j < nlocal) {
        drho[j] += imass*gradWmag*delVDotDelR;
      }
    }
  }

}
