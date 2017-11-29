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
    (if no contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2015-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#include <string.h>
#include "atom.h"
#include <mpi.h>
#include <cmath>
#include "modify.h"
#include "properties.h"
#include "error.h"
#include "memory.h"
#include "fix_multisphere.h"
#include "fix_property_atom.h"
#include "fix_property_global.h"

#define BIG 1e20
#define SMALL 1e-6

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Properties::Properties(LAMMPS *lmp): Pointers(lmp),
  ms_(0),
  ms_data_(0),
  mintype_(-1),
  maxtype_(-1),
  allow_soft_particles_(false),
  allow_hard_particles_(false)
{
}

/* ---------------------------------------------------------------------- */

Properties::~Properties()
{
}

/* ----------------------------------------------------------------------
   get max type used in the simulation
   is not necessarly equal to ntypes -1 as defined by create_box
   since not all atom types have to be used in the simulation
   error check so that atom types start with 1
------------------------------------------------------------------------- */

int Properties::max_type()
{
  // loop over all particles to check how many atom types are present
  mintype_ = 100000;
  maxtype_ = 1;

  for (int i=0;i<atom->nlocal;i++)
  {
      if (atom->type[i]<mintype_)
        mintype_=atom->type[i];
      if (atom->type[i]>maxtype_)
        maxtype_=atom->type[i];
  }

  // check all fixes
  // such as fix insert, fix change/type, fix wall, fix pour
  for(int i=0;i<modify->nfix;i++)
  {
      // checks
      Fix *fix = modify->fix[i];
      if(fix->min_type() > 0 &&  fix->min_type() < mintype_)
        mintype_ = fix->min_type();
      if(fix->max_type() > 0 &&  fix->max_type() > maxtype_)
        maxtype_ = fix->max_type();
  }

  //Get min/max from other procs
  int mintype_all,maxtype_all;
  MPI_Allreduce(&mintype_,&mintype_all, 1, MPI_INT, MPI_MIN, world);
  MPI_Allreduce(&maxtype_,&maxtype_all, 1, MPI_INT, MPI_MAX, world);
  mintype_ = mintype_all;
  maxtype_ = maxtype_all;

  //error check
  if(!lmp->wb)
  {
      if(mintype_ != 1)
        error->all(FLERR,"Atom types must start from 1 for granular simulations");
      if(maxtype_ > atom->ntypes)
        error->all(FLERR,"Please increase the number of atom types in the 'create_box' command to match the number of atom types you use in the simulation");
  }
  else
  {
      if(mintype_ != 1)
        error->all(FLERR,"Materials defined but not used in the simulation as particle or wall material must be the last materials defined");
      if(maxtype_ > atom->ntypes)
        error->all(FLERR,"Please increase the number of atom types in the 'create_box' command to match the number of atom types you use in the simulation");
  }
  return maxtype_;
}

/* ----------------------------------------------------------------------
   get minimum of radii used in the simulation
------------------------------------------------------------------------- */

double Properties::min_radius()
{
  const double maxtype = max_type();
  double minRadius = BIG;

  // check local particles
  for (int i=0;i<atom->nlocal;i++)
  {
    const double irad = atom->radius[i];
    // minimum
    if (irad < minRadius)
      minRadius = irad;
  }

  // check all fixes
  // such as fix insert, fix change/type, fix wall, fix pour
  for(int i=0;i<modify->nfix;i++)
  {
      // checks
      Fix *fix = modify->fix[i];

      if(!fix->use_rad_for_cut_neigh_and_ghost())
          continue;

      // loop over all types since min_rad(int) and max_rad(int) need a type
      for (int j=1;j<maxtype+1;j++)
      {
        const double f_minrad = fix->min_rad(j);
        if(f_minrad > SMALL &&  f_minrad < minRadius)
          minRadius = f_minrad;
      }
  }

  //Get min/max from other procs
  double minRadius_all;
  MPI_Allreduce(&minRadius,&minRadius_all, 1, MPI_DOUBLE, MPI_MIN, world);
  minRadius = minRadius_all;

  //error check
  if(minRadius <= SMALL)
    error->all(FLERR,"Atom radius must be bigger than zero for granular simulations");

  return  minRadius;
}

/* ----------------------------------------------------------------------
   get maximum of radii used in the simulation
------------------------------------------------------------------------- */

double Properties::max_radius()
{
  const double maxtype = max_type();
  double maxRadius = -1.0;

  // check local particles
  for (int i=0;i<atom->nlocal;i++)
  {
    const double irad = atom->radius[i];
    // maximum
    if (irad > maxRadius)
      maxRadius = irad;
  }

  // check all fixes
  // such as fix insert, fix change/type, fix wall, fix pour
  for(int i=0;i<modify->nfix;i++)
  {
      // checks
      Fix *fix = modify->fix[i];

      if(!fix->use_rad_for_cut_neigh_and_ghost())
          continue;

      // loop over all types since min_rad(int) and max_rad(int) need a type
      for (int j=1;j<maxtype+1;j++)
      {
        const double f_maxrad = fix->max_rad(j);
        if(f_maxrad > SMALL &&  f_maxrad > maxRadius)
          maxRadius = f_maxrad;
      }
  }

  //Get min/max from other procs
  double maxRadius_all;
  MPI_Allreduce(&maxRadius,&maxRadius_all, 1, MPI_DOUBLE, MPI_MAX, world);
  maxRadius = maxRadius_all;

  //error check
  if(maxRadius <= SMALL)
    error->all(FLERR,"Atom radius must be bigger than zero for granular simulations");

  return  maxRadius;
}

/* ----------------------------------------------------------------------
   find a property that was requested
   called e.g. from CFD data exchange model
   property may come from atom class, from a fix property, or fix ms
   last 2 args are the data length and are used for all data
------------------------------------------------------------------------- */

void* Properties::find_property(const char *name, const char *type, int &len1, int &len2)
{
    
    void *ptr = NULL;

    // possiblility 1
    // may be atom property - look up in atom class

    ptr = atom->extract(name,len2);
    // if nlocal == 0 and property found atom->extract returns NULL and len2 >= 0
    if(ptr || len2 >= 0)
    {
        len1 = atom->tag_max();
        // check if length correct
        if(((strcmp(type,"scalar-atom") == 0) && (len2 != 1)) || ((strcmp(type,"vector-atom") == 0) && (len2 != 3)))
            return NULL;
        else if(ptr && strstr(type,"multisphere"))
        {
            error->one(FLERR,"mismatch of data found and type specified");
            return NULL;
        }
        return ptr;
    }

    // possiblility 2
    // may come from a fix multisphere
    // also handles scalar-multisphere and vector-multisphere

    ms_ = static_cast<FixMultisphere*>(modify->find_fix_style("multisphere",0));
    if(ms_) ms_data_ = &ms_->data();

    if(ms_)
    {
        ptr = ms_->extract(name,len1,len2);
        if(((strcmp(type,"scalar-multisphere") == 0) && (len2 != 1)) || ((strcmp(type,"vector-multisphere") == 0) && (len2 != 3)))
            return NULL;
        
        else if(ptr && (strcmp(name,"body") && strstr(type,"atom")))
        {
            error->one(FLERR,"mismatch of data found and type specified");
            return NULL;
        }

        if(ptr || ((len1 >= 0) && (len2 >= 0)))
            return ptr;
    }

    // possiblility 3
    // may be fix property per atom - look up in modify class

    Fix *fix = NULL;

    if(strcmp(type,"scalar-atom") == 0)
    {
       fix = modify->find_fix_property(name,"property/atom","scalar",0,0,"cfd coupling",false);
       if(fix)
       {
           len1 = atom->tag_max();
           len2 = 1;
           return (void*) static_cast<FixPropertyAtom*>(fix)->vector_atom;
       }
    }
    else if(strcmp(type,"vector-atom") == 0)
    {
       fix = modify->find_fix_property(name,"property/atom","vector",0,0,"cfd coupling",false);
       if(fix)
       {
           len1 = atom->tag_max();
           len2 = 3;
           return (void*) static_cast<FixPropertyAtom*>(fix)->array_atom;
       }
    }
    else if(strcmp(type,"scalar-global") == 0)
    {
       fix = modify->find_fix_property(name,"property/global","scalar",0,0,"cfd coupling",false);
       len1 = len2 = 1;
       if(fix) return (void*) static_cast<FixPropertyGlobal*>(fix)->values;
    }
    else if(strcmp(type,"vector-global") == 0)
    {
       fix = modify->find_fix_property(name,"property/global","vector",0,0,"cfd coupling",false);
       if(fix)
       {
           len1 = static_cast<FixPropertyGlobal*>(fix)->nvalues;
           len2 = 1;
           return (void*) static_cast<FixPropertyGlobal*>(fix)->values;
       }
    }
    else if(strcmp(type,"matrix-global") == 0)
    {
       fix = modify->find_fix_property(name,"property/global","matrix",0,0,"cfd coupling",false);
       if(fix)
       {
           len1  = static_cast<FixPropertyGlobal*>(fix)->size_array_rows;
           len2  = static_cast<FixPropertyGlobal*>(fix)->size_array_cols;
           return (void*) static_cast<FixPropertyGlobal*>(fix)->array;
       }
    }
    else if(strcmp(name,"ex") == 0) 
    {
        // possiblility 4A - Dipole is specified as atom property (requires DIPOLE package)
        ptr = atom->extract("mu",len2);
        printf("len2 of mu: %d \n", len2);
        if(ptr)
        {
            len1 = atom->tag_max();
            // check if length correct
            if( (strcmp(type,"vector-atom") == 0) && (len2 != 3) )
                return NULL;
            return ptr;
        }

        // possiblility 4B - Quaternion is specified as atom property (requires ASPHERE package)
        // requires a fix that computes orientation data from quaternion,
        // or another fix that provides orientationEx (e.g., from POEMS)
        //TODO: Write fix that computes orientation from quaternion
        //TODO: Check if correct data is drawn in case a FixPOEMS is used
        fix = modify->find_fix_property("orientationEx","property/atom","vector",0,0,"cfd coupling",false);
        if(fix)
        {
               len1 = atom->tag_max();
               len2 = 3;
               return (void*) static_cast<FixPropertyAtom*>(fix)->array_atom;
        }
        else printf("WARNING: Fix with name 'orientationEx' not found that stores orientation information. \n");
    }
    return NULL;
}
