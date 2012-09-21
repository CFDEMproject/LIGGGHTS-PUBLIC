/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include "string.h"
#include "atom.h"
#include "mpi.h"
#include "math.h"
#include "modify.h"
#include "mech_param_gran.h"
#include "error.h"
#include "fix_property_global.h"
#include "memory.h"

using namespace LAMMPS_NS;

MechParamGran::MechParamGran(LAMMPS *lmp): Pointers(lmp)
{
}

MechParamGran::~MechParamGran()
{
}

int MechParamGran::max_type()
{
  //loop over all particles to check how many atom types are present
  mintype=100000;
  maxtype=1;

  for (int i=0;i<atom->nlocal;i++)
  {
      if (atom->type[i]<mintype) mintype=atom->type[i];
      if (atom->type[i]>maxtype) maxtype=atom->type[i];
  }

  // check all fixes
  // such as fix insert, fix change/type, fix wall, fix pour
  for(int i=0;i<modify->nfix;i++)
  {
      // checks
      Fix *fix = modify->fix[i];
      if(fix->min_type() > 0 &&  fix->min_type() < mintype)
        mintype = fix->min_type();
      if(fix->max_type() > 0 &&  fix->max_type() > maxtype)
        maxtype = fix->max_type();
  }

  //Get min/max from other procs
  int mintype_all,maxtype_all;
  MPI_Allreduce(&mintype,&mintype_all, 1, MPI_INT, MPI_MIN, world);
  MPI_Allreduce(&maxtype,&maxtype_all, 1, MPI_INT, MPI_MAX, world);
  mintype=mintype_all;
  maxtype=maxtype_all;

  //error check
  if(mintype != 1) error->all(FLERR,"Atom types must start from 1 for granular simulations");
  if(maxtype > atom->ntypes) error->all(FLERR,"Please increase the number of atom types in the 'create_box' command to match the number of atom types you use in the simulation");

  return maxtype;
}
