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
#include "fix_sph.h"
#include "update.h"
#include "respa.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "pair_sph.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "sph_kernels.h"
#include "modify.h"
#include "fix_property_atom.h"
#include "fix_property_global.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSph::FixSph(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  kernel_flag = 1;  // default: kernel is used
  kernel_id = -1;   // default value
  kernel_cut = -1;
  kernel_style = NULL;

  fppaSl = NULL;
  fppaSlType = NULL;
  sl = NULL;
  slComType = NULL;
}

/* ---------------------------------------------------------------------- */

FixSph::~FixSph()
{
  if(kernel_style) delete []kernel_style;
  if(mass_type) memory->destroy(slComType);
}

/* ---------------------------------------------------------------------- */

int FixSph::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSph::init()
{
  mass_type = atom->avec->mass_type;
  int ntypes = atom->ntypes;
  // need a half neighbor list, built when ever re-neighboring occurs

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;

  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  // check if kernel id is set
  if (kernel_flag && kernel_id < 0) error->all(FLERR,"No sph kernel for fixes is set.");
  // set kernel_cut
  kernel_cut = SPH_KERNEL_NS::sph_kernel_cut(kernel_id);

  // get the fix_property containing the smoothing length
  if (mass_type) {
    if (fppaSlType == NULL) {
    fppaSlType=static_cast<FixPropertyGlobal*>(modify->find_fix_property("sl","property/global","peratomtype",ntypes,0,force->pair_style));
    }
    if (!fppaSlType) error->all(FLERR,"Fix sph only works with a fix property/global that defines sl");

    // allocate memory for per atom-type property
    //TODO: copy slComType from pair?
    if (!slComType) memory->create(slComType,ntypes+1,ntypes+1,"fix:slComType");

    for (int i = 1; i <= ntypes; i++)
      for (int j = i; j <= ntypes; j++) {
        double sli = fppaSlType->compute_vector(i-1);
        double slj = fppaSlType->compute_vector(j-1);

        slComType[i][j] = slComType[j][i] = interpDist(sli,slj);;
      }

  } else {
    if (fppaSl == NULL) {
      fppaSl=static_cast<FixPropertyAtom*>(modify->find_fix_property("sl","property/atom","scalar",0,0,"FixSph",false));
    }
    if(!fppaSl) error->all(FLERR,"Fix sph only works with a fix property/atom that defines sl. Internal error!");
  }
}

/* ---------------------------------------------------------------------- */

void FixSph::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixSph::post_integrate_respa(int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_integrate();
}

/* ---------------------------------------------------------------------- */

void FixSph::updatePtrs()
{
  if (fppaSl) sl = fppaSl->vector_atom;
  if (fppaSlType) sl = fppaSlType->values;
}

/* ----------------------------------------------------------------------
   return common radius for types with different smoothing lengths
   atm: simple arithmetic mean (compare Morris)
------------------------------------------------------------------------- */
/*
inline double FixSph::interpDist(double disti, double distj)
{
  return 0.5*(disti+distj);
}
*/
