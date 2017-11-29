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
    (if not contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#include <string.h>
#include <stdlib.h>
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "memory.h"
#include "comm.h"
#include "modify.h"
#include <cmath>
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

void CfdDatacouplingMPI::allocate_external(int    **&data, int len2,const char *keyword,int initvalue)
{
  int len1 = 0;
  Multisphere *ms_data = properties_->ms_data();

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

void CfdDatacouplingMPI::allocate_external(double **&data, int len2,const char *keyword,double initvalue)
{
  int len1 = 0;
  Multisphere *ms_data = properties_->ms_data();

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
