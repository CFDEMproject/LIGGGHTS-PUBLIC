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
   CFD-DEM Coupling Stuff
------------------------------------------------------------------------- */

#include "mpi.h"
#include "string.h"
#include "library_cfd_coupling.h"
#include "lammps.h"
#include "input.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "fix_cfd_coupling.h"
#include "fix_multisphere.h"
#include "cfd_regionmodel.h"
#include "memory.h"
#include "error.h"
#include "comm.h"
#include "cfd_datacoupling.h"

using namespace LAMMPS_NS;

#define LMP_GROW_DELTA 11000

/* ---------------------------------------------------------------------- */

int liggghts_get_maxtag(void *ptr)
{
  LAMMPS *lmp = (LAMMPS *) ptr;
  return lmp->atom->tag_max();
}

/* ---------------------------------------------------------------------- */

int liggghts_get_maxtag_ms(void *ptr)
{
  // currently no possibility to delete multisphere bodies
  // so just return # of bodies

  LAMMPS *lmp = (LAMMPS *) ptr;
  FixMultisphere *fix_ms = static_cast<FixMultisphere*>(lmp->modify->find_fix_style_strict("multisphere",0));
  if(!fix_ms) return 0;
  return fix_ms->tag_max_body();
}

/* ---------------------------------------------------------------------- */

int liggghts_get_ntypes_ms(void *ptr)
{
  // currently no possibility to delete multisphere bodies
  // so just return # of bodies

  LAMMPS *lmp = (LAMMPS *) ptr;
  FixMultisphere *fix_ms = static_cast<FixMultisphere*>(lmp->modify->find_fix_style_strict("multisphere",0));
  if(!fix_ms) return 0;
  return fix_ms->ntypes();
}

/* ---------------------------------------------------------------------- */

double* liggghts_get_vclump_ms(void *ptr)
{
  // currently no possibility to delete multisphere bodies
  // so just return # of bodies

  LAMMPS *lmp = (LAMMPS *) ptr;
  FixMultisphere *fix_ms = static_cast<FixMultisphere*>(lmp->modify->find_fix_style_strict("multisphere",0));
  if(!fix_ms) return 0;
  return fix_ms->vclump();
}

/* ---------------------------------------------------------------------- */

void* locate_coupling_fix(void *ptr)
{
    LAMMPS *lmp = (LAMMPS *) ptr;
    int ifix = -1;
    for(int i=0;i<lmp->modify->nfix;i++)
      if(strcmp(lmp->modify->fix[i]->style,"couple/cfd") == 0)
        ifix = i;

    if(ifix ==-1) lmp->error->all(FLERR,"No fix of style 'couple/cfd' found, aborting.");

    return ((void*)lmp->modify->fix[ifix]);
}

/* ---------------------------------------------------------------------- */

void data_liggghts_to_of(char *name,char *type,void *ptr,void *&data,char* datatype)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->push(name,type,data,datatype);
}

/* ---------------------------------------------------------------------- */

void data_of_to_liggghts(char *name,char *type,void *ptr,void *data,char* datatype)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->pull(name,type,data,datatype);
}

/* ---------------------------------------------------------------------- */

void update_rm(void *ptr)
{
    LAMMPS *lmp = (LAMMPS *) ptr;
    //FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    locate_coupling_fix(ptr);
    //CfdRegionmodel *rm = fcfd->rm;

    //if(rm) rm->rm_update();
    lmp->error->all(FLERR,"Region model update not implemented aborting.");
}

/* ---------------------------------------------------------------------- */

void allocate_external_int(int    **&data, int len2,int len1,int    initvalue,void *ptr)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->allocate_external(data,len2,len1,initvalue);
}
/* ---------------------------------------------------------------------- */

void allocate_external_int(int    **&data, int len2,char *keyword,int    initvalue,void *ptr)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->allocate_external(data,len2,keyword,initvalue);
}

/* ---------------------------------------------------------------------- */

void allocate_external_double(double **&data, int len2,int len1,double initvalue,void *ptr)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->allocate_external(data,len2,len1,initvalue);
}

/* ---------------------------------------------------------------------- */

void allocate_external_double(double **&data, int len2,char* keyword,double initvalue,void *ptr)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->allocate_external(data,len2,keyword,initvalue);
}

/* ---------------------------------------------------------------------- */

void check_datatransfer(void *ptr)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->check_datatransfer();
}
