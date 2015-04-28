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

#include "string.h"
#include "atom.h"
#include "mpi.h"
#include "math.h"
#include "modify.h"
#include "properties.h"
#include "error.h"
#include "memory.h"
#include "fix_multisphere.h"
#include "multisphere.h"
#include "fix_property_atom.h"
#include "fix_property_global.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Properties::Properties(LAMMPS *lmp): Pointers(lmp),
  ms_(0),
  ms_data_(0)
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
  mintype = 100000;
  maxtype = 1;

  for (int i=0;i<atom->nlocal;i++)
  {
      if (atom->type[i]<mintype)
        mintype=atom->type[i];
      if (atom->type[i]>maxtype)
        maxtype=atom->type[i];
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
  mintype = mintype_all;
  maxtype = maxtype_all;

  //error check
  if(mintype != 1)
    error->all(FLERR,"Atom types must start from 1 for granular simulations");
  if(maxtype > atom->ntypes)
    error->all(FLERR,"Please increase the number of atom types in the 'create_box' command to match the number of atom types you use in the simulation");

  return maxtype;
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
        return ptr;
    }

    // possiblility 2
    // may come from a fix multisphere
    // also handles scalar-multisphere and vector-multisphere

    ms_ = static_cast<FixMultisphere*>(modify->find_fix_style_strict("multisphere",0));
    if(ms_) ms_data_ = &ms_->data();

    if(ms_)
    {
        ptr = ms_->extract(name,len1,len2);
        if(((strcmp(type,"scalar-multisphere") == 0) && (len2 != 1)) || ((strcmp(type,"vector-multisphere") == 0) && (len2 != 3)))
            return NULL;

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
