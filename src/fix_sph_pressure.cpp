/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
   the producer of the LIGGGHTS® software and the CFDEM®coupling software
   See http://www.cfdem.com/terms-trademark-policy for details.

   LIGGGHTS® is based on LAMMPS
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
-------------------------------------------------------------------------
Contributing author for SPH:
Andreas Eitzlmayr (Institute for Process and Particle Engineering, TU Graz)
andreas.eitzlmayr@tugraz.at
------------------------------------------------------------------------- */

#include "math.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "fix_sph_pressure.h"
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
#include "timer.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSPHPressure::FixSPHPressure(LAMMPS *lmp, int narg, char **arg) :
  FixSph(lmp, narg, arg)
{
    //Check args
    int iarg = 3;
    if (narg < iarg+1) error->fix_error(FLERR,this,"Not enough arguments \n");

    if (strcmp(arg[iarg],"absolut") == 0)
    {
      if (narg < iarg+2) error->fix_error(FLERR,this,"Not enough arguments for 'absolut' pressure style \n");
      B = force->numeric(FLERR,arg[iarg+1]);
      pressureStyle = PRESSURESTYLE_ABSOLUT;
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"Tait") == 0)
    {
      if (narg < iarg+4) error->fix_error(FLERR,this,"Not enough arguments for 'Tait' pressure style \n");
      B = force->numeric(FLERR,arg[iarg+1]);
      rho0 = force->numeric(FLERR,arg[iarg+2]);
      if (rho0 > 0) rho0inv = 1./rho0;
      else error->fix_error(FLERR,this," rho0 is zero or negativ \n");
      gamma = force->numeric(FLERR,arg[iarg+3]);
      if (narg < iarg+5)
      {
        P0 = 0;
        iarg += 4;
      }
      else
      {
        P0 = force->numeric(FLERR,arg[iarg+4]);
        iarg += 5;
      }
      pressureStyle = PRESSURESTYLE_TAIT;
    }
    else if (strcmp(arg[iarg],"relativ") == 0)
    {
      if (narg < iarg+3) error->fix_error(FLERR,this,"Not enough arguments for 'relativ' pressure style \n");
      B = force->numeric(FLERR,arg[iarg+1]);
      rho0 = force->numeric(FLERR,arg[iarg+2]);
      if (narg < iarg+4)
      {
        P0 = 0;
        iarg += 3;
      }
      else
      {
        P0 = force->numeric(FLERR,arg[iarg+3]);
        iarg += 4;
      }
      pressureStyle = PRESSURESTYLE_RELATIV;
    }
    else error->fix_error(FLERR,this,"Unknown style. Valid styles are 'absolut' or 'Tait' \n");

    kernel_flag = 0; // does not need any kernel
}

/* ---------------------------------------------------------------------- */

FixSPHPressure::~FixSPHPressure()
{

}

/* ---------------------------------------------------------------------- */

int FixSPHPressure::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void FixSPHPressure::init()
{
  FixSph::init();

  // check if there is an sph/density fix present
  // must come before me, because -- (sph/density has post_integrate routine... it will be called first)
  // a - need the pressure for the rho
  // b - does the forward comm for me to have updated ghost properties
  // check is done by fix sph/density itself

  int dens = -1;
  for(int i = 0; i < modify->nfix; i++)
  {
    if(strncmp("sph/density",modify->fix[i]->style,11) == 0) {
      dens = i;
      break;
    }
  }

  if(dens == -1) error->fix_error(FLERR,this,"Requires to define a fix sph/density also \n");
}

/* ---------------------------------------------------------------------- */

void FixSPHPressure::pre_force(int vflag)
{
  int *mask = atom->mask;
  double *rho = atom->rho;
  double *p = atom->p;
  int nlocal = atom->nlocal;

  // already have updated ghost positions due to regular communication

  // set pressure

  if (pressureStyle == PRESSURESTYLE_TAIT)
  {
    for (int i = 0; i < nlocal; i++)
    {

      if (mask[i] & groupbit)
      {

      p[i] = B*(pow(rho[i]*rho0inv,gamma) - 1) + P0; // Tait's equation
      // Added a background pressure P0 (see e.g., S. Adami, X.Y. Hu, N.A. Adams,
      // J. Comput. Phys. 241 (2013) 292-307)
    }
    }
  }
  else if (pressureStyle == PRESSURESTYLE_RELATIV)
  {
    for (int i = 0; i < nlocal; i++)
    {
      if (mask[i] & groupbit)
      {
        p[i] = B * (rho[i] - rho0) + P0;
      }
    }
  }
  else if (pressureStyle == PRESSURESTYLE_ABSOLUT)
  {
    for (int i = 0; i < nlocal; i++)
    {
      if (mask[i] & groupbit)
      {
        p[i] = B * B * rho[i];
      }
    }
  }
}
