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
#include "fix_sph_mixidx.h"
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

FixSphMixidx::FixSphMixidx(LAMMPS *lmp, int narg, char **arg) :
  FixSph(lmp, narg, arg)
{
  int iarg = 0;

  if (iarg+3 > narg) error->fix_error(FLERR,this,"Not enough input arguments");

  iarg += 3;

  every = 1; // optional argument, default value = 1
  ago = 0;

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
  fix_gamma_ = NULL;
  fix_omega_ = NULL;
  fix_mixidx_ = NULL;
}

/* ---------------------------------------------------------------------- */

FixSphMixidx::~FixSphMixidx()
{

}

/* ---------------------------------------------------------------------- */

void FixSphMixidx::post_create()
{
  fix_dvdx_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("dvdx","property/atom","vector",0,0,"FixSphMixidx",false));
  fix_dvdy_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("dvdy","property/atom","vector",0,0,"FixSphMixidx",false));
  fix_dvdz_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("dvdz","property/atom","vector",0,0,"FixSphMixidx",false));

  const char *fixarg[9];
  fixarg[0]="omega";
  fixarg[1]="all";
  fixarg[2]="property/atom";
  fixarg[3]="omega";
  fixarg[4]="scalar";
  fixarg[5]="yes";
  fixarg[6]="yes";
  fixarg[7]="yes";
  fixarg[8]="0";
  fix_omega_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);

  fixarg[0]="mixidx";
  fixarg[3]="mixidx";
  fix_mixidx_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  char* fixID;
  int ifix_gamma = -1;

  for (int ifix = 0; ifix < modify->nfix; ifix++)
  {
    fixID = modify->fix[ifix]->id;

    if (strcmp("gamma",fixID) == 0) {
      ifix_gamma = ifix;
    }
  }

  if (ifix_gamma == -1) { // create fix for gamma
    fixarg[0]="gamma";
    fixarg[3]="gamma";
    fix_gamma_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    calcgamma = 1;
  } else { // gamma already existing (from pair_sph_morris_tenscorr)
    fix_gamma_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("gamma","property/atom","scalar",0,0,"FixSphMixidx",false));
    calcgamma = 0;
  }
}

/* ---------------------------------------------------------------------- */

void FixSphMixidx::pre_delete(bool unfixflag)
{
    if(unfixflag && fix_dvdx_)
        modify->delete_fix(fix_dvdx_->id);

    if(unfixflag && fix_dvdy_)
        modify->delete_fix(fix_dvdy_->id);

    if(unfixflag && fix_dvdz_)
        modify->delete_fix(fix_dvdz_->id);

    if (calcgamma == 1) {
      if(unfixflag && fix_gamma_)
        modify->delete_fix(fix_gamma_->id);
    }

    if(unfixflag && fix_omega_)
        modify->delete_fix(fix_omega_->id);

    if(unfixflag && fix_mixidx_)
        modify->delete_fix(fix_mixidx_->id);
}

/* ---------------------------------------------------------------------- */

int FixSphMixidx::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSphMixidx::init()
{
  FixSph::init();

  // check if there is an sph/velgrad fix present

  int ifix_velgrad = -1;

  for(int i = 0; i < modify->nfix; i++)
  {
    if(strncmp("sph/velgrad",modify->fix[i]->style,11) == 0) {
      ifix_velgrad = i;
    }
  }

  if(ifix_velgrad == -1) error->fix_error(FLERR,this,"Requires to define a fix sph/velgrad also \n");
}

/* ---------------------------------------------------------------------- */

void FixSphMixidx::post_force(int vflag)
{
  //template function for using per atom or per atomtype smoothing length
  if (mass_type) post_force_eval<1>(vflag);
  else post_force_eval<0>(vflag);
}

/* ---------------------------------------------------------------------- */

template <int MASSFLAG>
void FixSphMixidx::post_force_eval(int vflag)
{
  double D11,D22,D33,D12,D13,D23,W12,W13,W23,gamma2,omega2;
  int nlocal = atom->nlocal;
  int i;

  ago++;
  if (ago % every == 0) {
    ago = 0;

    dvdx_ = fix_dvdx_->array_atom;
    dvdy_ = fix_dvdy_->array_atom;
    dvdz_ = fix_dvdz_->array_atom;
    gamma_ = fix_gamma_->vector_atom;
    omega_ = fix_omega_->vector_atom;
    mixidx_ = fix_mixidx_->vector_atom;

    // calculate shear rate and vorticity magnitude, mixing index
    for (i = 0; i < nlocal; i++) {

      // Components of rate-of-strain tensor D and vorticity tensor W,
      // see e.g., G. Boehme, Stroemungsmechanik nicht-newtonscher Fluide,
      // B.G. Teubner, Stuttgart, 1981.

      if (calcgamma == 1) {
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

      W12 = 0.5 * (dvdy_[i][0] - dvdx_[i][1]);
      W13 = 0.5 * (dvdz_[i][0] - dvdx_[i][2]);
      W23 = 0.5 * (dvdz_[i][1] - dvdy_[i][2]);

      // omega = = norm(nabla x v) = sqrt(2*W:W) = sqrt(2*Frobeniusnorm(W))
      omega2 = 4*(W12*W12 + W13*W13 + W23*W23);
      omega_[i] = sqrt(omega2);

      // The mixing Index (called similarly in Fluent) describes the local flow type
      // (0 ... rotational flow / 0.5 ... shear flow / 1 ... elongational flow)
      // similarly shown by:
      // I. Manas-Zloczower, Mixing and Compounding of Polymers, Carl Hanser Verlag, Munich, 2009.
      // Ansys Polyflow User's Guide, October 2012, Eq. 28.10.

      mixidx_[i] = gamma_[i] / (gamma_[i] + omega_[i] + 1e-10);
    }

    // dvdx, dvdy, dvdz are now correct, send to ghosts
    timer->stamp();
    comm->forward_comm();
    timer->stamp(TIME_COMM);
  }
}
