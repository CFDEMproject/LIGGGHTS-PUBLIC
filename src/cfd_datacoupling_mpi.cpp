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

#include "string.h"
#include "stdlib.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "memory.h"
#include "comm.h"
#include "modify.h"
#include "math.h"
#include "vector_liggghts.h"
#include "fix_cfd_coupling.h"
#include "fix_multisphere.h"
#include "cfd_datacoupling_mpi.h"

using namespace LAMMPS_NS;

#define DELTA 10000

/* ---------------------------------------------------------------------- */

CfdDatacouplingMPI::CfdDatacouplingMPI(LAMMPS *lmp,int iarg, int narg, char **arg,FixCfdCoupling* fc) :
  CfdDatacoupling(lmp, iarg, narg, arg,fc)
{
  liggghts_is_active = false;

  if(!atom->tag_enable) error->one(FLERR,"CFD-DEM coupling via MPI requires particles to have tags");

  this->fc_ = fc;

  len_allred_double = 0;
  allred_double = NULL;
  len_allred_int = 0;
  allred_int = NULL;

  if(comm->me == 0) error->message(FLERR,"nevery as specified in LIGGGHTS is overriden by calling external program",1);

}

CfdDatacouplingMPI::~CfdDatacouplingMPI()
{
    memory->sfree(allred_double);
    memory->sfree(allred_int);
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingMPI::exchange()
{
    // does nothing since done by OF
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingMPI::pull(const char *name,const char *type,void *&from,const char *datatype)
{
    CfdDatacoupling::pull(name,type,from,datatype);

    if(strcmp(datatype,"double") == 0)
        pull_mpi<double>(name,type,from);
    else if(strcmp(datatype,"int") == 0)
        pull_mpi<int>(name,type,from);
    else error->one(FLERR,"Illegal call to CfdDatacouplingMPI::pull, valid datatypes are 'int' and double'");
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingMPI::push(const char *name,const char *type,void *&to,const char *datatype)
{
    CfdDatacoupling::push(name,type,to,datatype);

    if(strcmp(datatype,"double") == 0)
        push_mpi<double>(name,type,to);
    else if(strcmp(datatype,"int") == 0)
        push_mpi<int>(name,type,to);
    else error->one(FLERR,"Illegal call to CfdDatacouplingMPI::pull, valid datatypes are 'int' and double'");
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingMPI::allocate_external(int **&data, int len2,int len1,int initvalue)
{
  if(len1 < 1 || len2 < 1)
    error->one(FLERR,"Illegal length used in CfdDatacouplingMPI::allocate_external");

  memory->grow(data, len1,len2, "CfdDatacouplingMPI:data");
  for (int i = 0; i < len1; i++)
    for (int j = 0; j < len2; j++)
      data[i][j] = initvalue;
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingMPI::allocate_external(int    **&data, int len2,char *keyword,int initvalue)
{
  int len1 = 0;
  MultisphereParallel *ms_data = properties_->ms_data();

  if(strcmp(keyword,"nparticles") == 0) len1 = atom->tag_max();
  else if(strcmp(keyword,"nbodies") == 0)
  {
      if(ms_data)
        len1 = ms_data->tag_max_body();
      else error->one(FLERR,"CFD datacoupling keyword 'nbodies' may only be used with multisphere model in LIGGGHTS");
  }
  else error->one(FLERR,"Illegal keyword used in CfdDatacouplingMPI::allocate_external");
  if(len1 < 1 || len2 < 1)
   len1 = len2 = 1;

  memory->grow(data, len1,len2, "CfdDatacouplingMPI:data");
  for (int i = 0; i < len1; i++)
    for (int j = 0; j < len2; j++)
      data[i][j] = initvalue;
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingMPI::allocate_external(double **&data, int len2,int len1,double initvalue)
{
  if(len1 < 1 || len2 < 1)
    error->one(FLERR,"Illegal length used in CfdDatacouplingMPI::allocate_external");

  memory->grow(data, len1,len2, "CfdDatacouplingMPI:data");
  for (int i = 0; i < len1; i++)
    for (int j = 0; j < len2; j++)
        data[i][j] = initvalue;
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingMPI::allocate_external(double **&data, int len2,char *keyword,double initvalue)
{
  int len1 = 0;
  MultisphereParallel *ms_data = properties_->ms_data();

  if(strcmp(keyword,"nparticles") == 0) len1 = atom->tag_max();
  else if(strcmp(keyword,"nbodies") == 0)
  {
      if(ms_data)
        len1 = ms_data->tag_max_body();
      else error->one(FLERR,"CFD datacoupling keyword 'nbodies' may only be used with multisphere model in LIGGGHTS");
  }
  else error->one(FLERR,"Illegal keyword used in CfdDatacouplingMPI::allocate_external");
  if(len1 < 1 || len2 < 1)
    len1 = len2 = 1;

  memory->grow(data, len1,len2, "CfdDatacouplingMPI:data");
  for (int i = 0; i < len1; i++)
    for (int j = 0; j < len2; j++)
        data[i][j] = initvalue;
}

/* ---------------------------------------------------------------------- */
