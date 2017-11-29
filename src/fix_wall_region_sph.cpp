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
#include <stdlib.h>
#include <string.h>
#include "fix_wall_region.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "region.h"
#include "lattice.h"
#include "update.h"
#include "output.h"
#include "respa.h"
#include "error.h"

#include "fix_wall_region_sph.h"
#include "sph_kernels.h"
#include "force.h"
#include "pair.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixWallRegionSph::FixWallRegionSph(LAMMPS *lmp, int narg, char **arg) :
  FixSph(lmp, narg, arg)
{
  if (narg != 6) error->all(FLERR,"Illegal fix wall/region/sph command");

  int iarg = 3;

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;

  // parse args

  iregion = domain->find_region(arg[iarg]);
  if (iregion == -1) error->all(FLERR,"Fix wall/region/sph region ID does not exist");

  r0 = force->numeric(FLERR,arg[iarg+1]);
  D  = force->numeric(FLERR,arg[iarg+2]);

  iarg += 3;

  //only for force!
  cutoff = r0;
  if (cutoff <= 0.0) error->all(FLERR,"Fix wall/region/sph cutoff <= 0.0");

  eflag = 0;
  ewall[0] = ewall[1] = ewall[2] = ewall[3] = 0.0;
}

/* ---------------------------------------------------------------------- */

int FixWallRegionSph::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallRegionSph::init()
{

  FixSph::init();

  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixWallRegionSph::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixWallRegionSph::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixWallRegionSph::post_force(int vflag)
{
  //template function for using per atom or per atomtype smoothing length
  if (mass_type) post_force_eval<1>(vflag);
  else post_force_eval<0>(vflag);
}

/* ---------------------------------------------------------------------- */

void FixWallRegionSph::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixWallRegionSph::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   energy of wall interaction
------------------------------------------------------------------------- */

double FixWallRegionSph::compute_scalar()
{
  // only sum across procs one time

  if (eflag == 0) {
    MPI_Allreduce(ewall,ewall_all,4,MPI_DOUBLE,MPI_SUM,world);
    eflag = 1;
  }
  return ewall_all[0];
}

/* ----------------------------------------------------------------------
   components of force on wall
------------------------------------------------------------------------- */

double FixWallRegionSph::compute_vector(int n)
{
  // only sum across procs one time

  if (eflag == 0) {
    MPI_Allreduce(ewall,ewall_all,4,MPI_DOUBLE,MPI_SUM,world);
    eflag = 1;
  }
  return ewall_all[n+1];
}

/* ---------------------------------------------------------------------- */

template <int MASSFLAG>
void FixWallRegionSph::post_force_eval(int vflag)
{
  int i,m,n;
  double fx,fy,fz,sli,imass;

  eflag = 0;
  ewall[0] = ewall[1] = ewall[2] = ewall[3] = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;

  Region *region = domain->regions[iregion];
  int onflag = 0;

  updatePtrs(); // get sl

  // TODO: mass_type dependent declaration?
  int itype;
  int *type = atom->type;
  double *mass = atom->mass;

  double *rmass = atom->rmass;

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      if (MASSFLAG) {
        itype = type[i];
        sli = sl[itype-1];
        imass = mass[itype];
      } else {
        sli = sl[i];
        imass = rmass[i];
      }

      // region->match() insures particle is in region or on surface, else error
      // if returned contact dist r = 0, is on surface, also an error
      if (!region->match(x[i][0],x[i][1],x[i][2])) {
        onflag = 1;
        fprintf(screen,"Particle %d with the Coordinates x= %f, y= %f, z= %f is on or inside fix wall/region/sph surface. \n",tag[i],x[i][0],x[i][1],x[i][2]); //TEST OUTPUT
        continue;
      }

      // number of wall contacts
      n = region->surface(x[i][0],x[i][1],x[i][2],cutoff);

      for (m = 0; m < n; m++) {

        // check wall distance
        if (region->contact[m].r <= 0.0) {
          onflag = 1;
          fprintf(screen,"Particle %d with the Coordinates x= %f, y= %f, z= %f has zero distance. \n",tag[i],x[i][0],x[i][1],x[i][2]); // TEST OUTPUT
          continue;
        }

        // add fwall, due to self influence
        fwall = selfInfluenceForce(i,region->contact[m].r,sli,imass);

        // add fwall, due to repulsive force
        fwall += repulsivSph(region->contact[m].r);

        ewall[0] += eng;
        fx = fwall * region->contact[m].delx;
        fy = fwall * region->contact[m].dely;
        fz = fwall * region->contact[m].delz;
        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
        ewall[1] -= fx;
        ewall[2] -= fy;
        ewall[3] -= fz;

      }
    }
  }

  if (onflag) {
    error->one(FLERR,"Particle on or inside fix wall/region/sph surface \n");
  }

}

/* ----------------------------------------------------------------------
   repulsiv force for sph simulation
   compute eng and fwall = magnitude of wall force
------------------------------------------------------------------------- */

double FixWallRegionSph::repulsivSph(double r)
{
  double rinv,frac,frac2;

  if (r <= r0) {
    rinv = 1.0/r;
    frac = r0*rinv;
    frac2 = frac*frac;
    eng = 0;
    return (D * frac2 * (frac2 - 1) * rinv);
  } else {
    eng = 0;
    return 0;
  }
}

/* ----------------------------------------------------------------------
   force due to self influence
------------------------------------------------------------------------- */

double FixWallRegionSph::selfInfluenceForce(int ip, double r, double sl, double mass)
{
  double isl,ir,s,gradWmag;
  double *rho = atom->rho;
  double *p = atom->p;

  // fwall, due to self influence
  // a small perturbation to avoid stacking particles
  ir = 1./r;
  isl = 1./sl;
  s = r*isl;

  // calculate value for magnitude of grad W
  gradWmag = SPH_KERNEL_NS::sph_kernel_der(kernel_id,s,sl,isl);

  // zero order approximation - properties at wall assumed to be equal to properties at particle (thuis factor 2)
  return (- ir * mass * mass * (2.*p[ip]/(rho[ip]*rho[ip])) * gradWmag);
}
