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

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "comm.h"
#include "update.h"
#include "region.h"
#include "domain.h"
#include "fix_region_variable.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixRegionVariable::FixRegionVariable(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  
  restart_global = 1;
  time_depend = 1;

  if (narg < 4) error->all(FLERR,"Illegal fix region/variable command, not enough arguments");
  iarg = 3;

  if(strcmp("n_regions",arg[iarg++])) error->all(FLERR,"Illegal fix region/variable command, expecting keyword 'n_regions'");
  n_regions = atoi(arg[iarg++]);
  if(n_regions < 1) error->all(FLERR,"Illegal fix region/variable command, expecting n_regions >= 1");

  steps = new double[n_regions];
  times = new double[n_regions];
  regions = (Region **) memory->smalloc(n_regions*sizeof(Region*),"FixRegionVariable:regions");

  step_start = update->ntimestep;

  bool t_set = false, r_set = false;

  // parse further args
  bool hasargs = true;
  while (iarg < narg && hasargs)
  {
    hasargs = false;
    if (strcmp(arg[iarg],"times") == 0) {
      for(int i = 0; i < n_regions; i++)
      {
          times[i] = atof(arg[iarg+1+i]);
          if (times[i] <= 0) error->all(FLERR,"Illegal fix region/variable command, time must be > 0");
      }
      hasargs = true;
      t_set = true;
      iarg += 1+n_regions;
    }
    else if (strcmp(arg[iarg],"regions") == 0) {
      if (iarg+1+n_regions > narg) error->all(FLERR,"Illegal fix region/variable command, not enough arguments for regions");

      for(int i = 0; i < n_regions; i++)
      {
          int ireg = domain->find_region(arg[iarg+1+i]);
          if (ireg < 0) error->all(FLERR,"Illegal fix region/variable command, illegal region");
          regions[i] = domain->regions[ireg];
      }
      hasargs = true;
      r_set = true;
      iarg += 1+n_regions;
    }
    else if(strcmp(style,"region/variable") == 0) error->all(FLERR,"Illegal fix region/variable command, unrecognized keyword");
  }

  if(!r_set || !t_set) error->all(FLERR,"Illegal fix region/variable command, must set values for times and regions");
}

/* ---------------------------------------------------------------------- */

FixRegionVariable::~FixRegionVariable()
{
  delete []steps;
  delete []times;
}

/* ----------------------------------------------------------------------*/

int FixRegionVariable::setmask()
{
  int mask = 0;
  return mask;
}

/* ----------------------------------------------------------------------*/

void FixRegionVariable::init()
{
    dt = update->dt;

    for(int i = 0; i < n_regions; i++)
    {
        steps[i] = times[i] / dt;
    }
}

/* ----------------------------------------------------------------------*/

Region* FixRegionVariable::region()
{
    int i_chosen;
    double steps_evolved = static_cast<double>(update->ntimestep - step_start);

    // choose region, assuming periodicity
    i_chosen = -1;
    while (steps_evolved > 0)
    {
        i_chosen++;
        if(i_chosen == n_regions) i_chosen = 0;
        steps_evolved -= steps[i_chosen];
    }
    return regions[i_chosen];
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixRegionVariable::write_restart(FILE *fp)
{
  int n = 0;
  double list[2];
  list[n++] = dt;
  list[n++] = static_cast<int>(step_start);

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixRegionVariable::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  dt = list[n++];
  step_start = static_cast<int> (list[n++]);
  if(dt != update->dt)
  {
      if(comm->me == 0) fprintf(screen,"Fix region/variable used a time-step of %f, you are now using %f\n",dt,update->dt);
      error->all(FLERR,"This is fatal");
  }
}
