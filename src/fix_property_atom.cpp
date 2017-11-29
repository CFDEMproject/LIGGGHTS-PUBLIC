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

#include <cmath>
#include <stdlib.h>
#include <string.h>
#include "fix_property_atom.h"
#include "atom.h"
#include "memory.h"
#include "error.h"

#include "pair_gran.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "group.h"
#include "timer.h"
#include "neighbor.h"

#include "mpi_liggghts.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define EPSILON 0.001

/* ---------------------------------------------------------------------- */

FixPropertyAtom::FixPropertyAtom(LAMMPS *lmp, int narg, char **arg, bool parse) :
  Fix(lmp, narg, arg),
  propertyname(0),
  property(0),
  internal(false)
{
    
    if(parse) parse_args(narg,arg);
}

void FixPropertyAtom::parse_args(int narg, char **arg)
{
    // Check args
    if (narg < 9) error->all(FLERR,"Illegal fix property/atom command, not enough arguments");
    if (narg > 29) error->warning(FLERR,"Vector length in fix property/atom larger than 20. Are you sure you want that?");

    // Read args
    
    int n = strlen(arg[3]) + 1;
    variablename = new char[n];
    strcpy(variablename,arg[3]);

    bool vector_with_one_entry = false;
    if (strcmp(arg[4],"scalar") == 0) data_style = FIXPROPERTY_ATOM_SCALAR;
    else if (strcmp(arg[4],"vector") == 0) data_style = FIXPROPERTY_ATOM_VECTOR;
    // This vector style allows for a vector to have only one entry. Under normal circumstances the scalar style should be chosen instead.
    else if (strcmp(arg[4],"vector_one_entry") == 0)
    {
        vector_with_one_entry = true;
        data_style = FIXPROPERTY_ATOM_VECTOR;
    }
    else error->all(FLERR,"Unknown style for fix property/atom. Valid styles are 'scalar', 'vector' or 'vector_one_entry'");

    if (strcmp(arg[5],"yes") == 0)
    {
            restart_peratom = 1;
            restart_global = 1;
    }
    else if (strcmp(arg[5],"no") == 0)
    {
         restart_peratom = 0;
         restart_global = 0;
    }
    else error->all(FLERR,"Unknown restart style for fix property/atom. Valid styles are 'yes' or 'no'");

    if (strcmp(arg[6],"yes") == 0) commGhost = 1;
    else if (strcmp(arg[6],"no") == 0) commGhost = 0;
    else error->all(FLERR,"Unknown communicate_ghost style for fix property/atom. Valid styles are 'yes' or 'no'");

    if (strcmp(arg[7],"yes") == 0) commGhostRev = 1;
    else if (strcmp(arg[7],"no") == 0) commGhostRev = 0;
    else error->all(FLERR,"Unknown communicate_reverse_ghost style for fix property/atom. Valid styles are 'yes' or 'no'");

    nvalues = narg - 8;
    
    if ((nvalues == 1) && !(data_style == FIXPROPERTY_ATOM_SCALAR || (vector_with_one_entry && data_style == FIXPROPERTY_ATOM_VECTOR)))
      error->all(FLERR,"Error in fix property/atom: Number of default values provided not consistent with vector style. Provide more than 1 value or use style 'scalar'");

    if ((nvalues >1) && (data_style != FIXPROPERTY_ATOM_VECTOR))
      error->all(FLERR,"Error in fix property/atom: Number of default values provided not consistent with scalar style. Provide 1 value or use style 'vector'");

    defaultvalues = new double[nvalues];

    // fix handles properties that need to be initialized at particle creation
    create_attribute = 1;

    // special case: only implemented for scalar default value from scalar
    // initialize property from another existing scalar, which is found via Property class
    // since this fix must exist already, it is initialized before this fix property
    propertyname = 0;
    if(FIXPROPERTY_ATOM_SCALAR == data_style)
    {
        char *prop = arg[8];
        int n = strlen(prop);
        bool is_digit = false;
        for (int i = 0; i < n; i++)
            if (isdigit(prop[i])) is_digit = true;

        if(!is_digit)
        {
            // scalar property found
            int len1,len2;
            if(atom->get_properties()->find_property(prop,"scalar-atom",len1,len2))
            {
                propertyname = new char[n+1];
                strcpy(propertyname,prop);
            }
        }
    }

    // set default values
    if(!propertyname)
    {
        for (int j = 0; j < nvalues; j++)
        {
            // if any of the values is none, this fix will not init properties
            if(strcmp(arg[8+j],"none") == 0)
            {
                create_attribute = 0;
                continue;
            }
            defaultvalues[j] = force->numeric(FLERR,arg[8+j]);
        }
    }

    if (data_style) size_peratom_cols = nvalues;
    else size_peratom_cols = 0;

    peratom_flag=1; 
    peratom_freq=1;
    extvector=0; 

    if (commGhost) comm_forward = nvalues;
    if (commGhostRev) comm_reverse = nvalues;

    // perform initial allocation of atom-based array
    // register with Atom class
    vector_atom = NULL; array_atom = NULL;
    grow_arrays(atom->nmax); 
    atom->add_callback(0); 
    if (restart_peratom) atom->add_callback(1); 

    // init all arrays since dump may access it on timestep 0
    // or a variable may access it before first run
    
    int nlocal = atom->nlocal;
    if(create_attribute)
    {
        if(propertyname)
        {
            // not implemented
            if(data_style)
                error->all(FLERR,"internal error");
            
            pre_set_arrays();
            for (int i = 0; i < nlocal; i++)
                vector_atom[i] = property[i];
        }
        else
        {
            for (int i = 0; i < nlocal; i++)
            {
              if (data_style)
              {
                for (int m = 0; m < nvalues; m++)
                    array_atom[i][m] = defaultvalues[m];
              }
              else vector_atom[i] = defaultvalues[0];
            }
        }
    }

    // check if there is already a fix that tries to register a property with the same name
    
    for (int ifix = 0; ifix < modify->nfix; ifix++)
        if ((modify->fix[ifix]) && (strcmp(modify->fix[ifix]->style,style) == 0) && (strcmp(((FixPropertyAtom*)(modify->fix[ifix]))->variablename,variablename)==0) )
            error->fix_error(FLERR,this,"there is already a fix that registers a variable of the same name");

    // flags for vector output
    //vector_flag = 1;
    size_vector = nvalues;
    global_freq = 1;
    extvector = 1;
}

/* ---------------------------------------------------------------------- */

FixPropertyAtom::~FixPropertyAtom()
{
  // unregister callbacks to this fix from Atom class
  atom->delete_callback(id,0);
  if (restart_peratom) atom->delete_callback(id,1);

  // delete locally stored arrays
  delete []variablename;
  delete []defaultvalues;
  if(propertyname) delete []propertyname;

  if (data_style) memory->destroy(array_atom);
  else memory->destroy(vector_atom);
}

/* ---------------------------------------------------------------------- */

Fix* FixPropertyAtom::check_fix(const char *varname,const char *svmstyle,int len1,int len2,const char *caller,bool errflag)
{
    char errmsg[400];
    
    if(strcmp(varname,variablename) == 0)
    {
        if(strcmp(svmstyle,"scalar") == 0) len1 = 1;

        // check variable style
        if(
            (strcmp(svmstyle,"scalar") == 0 && data_style != FIXPROPERTY_ATOM_SCALAR) ||
            (strcmp(svmstyle,"vector") == 0 && data_style != FIXPROPERTY_ATOM_VECTOR) ||
            (strcmp(svmstyle,"vector2D") == 0 && data_style != FIXPROPERTY_ATOM_VECTOR2D) ||
            (strcmp(svmstyle,"quaternion") == 0 && data_style != FIXPROPERTY_ATOM_QUATERNION)
        )
        {
            if(errflag)
            {
                sprintf(errmsg,"%s style required for fix property/atom variable %s for usage with caller %s",
                        svmstyle,varname,caller);
                error->all(FLERR,errmsg);
            }
            else return NULL;
        }

        // check length
        if(len1 > nvalues)
        {
            if(errflag)
            {
                sprintf(errmsg,"Fix property/atom variable %s has wrong length (length is %d but length %d expected) for usage with caller %s",
                        varname,nvalues,len1,caller);
                error->all(FLERR,errmsg);
            }
            else return NULL;
        }

        // success
        return static_cast<Fix*>(this);
    }
    return NULL;
}

/* ---------------------------------------------------------------------- */

int FixPropertyAtom::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ----------------------------------------------------------------------
   forward and backward comm to be used by other fixes as needed
------------------------------------------------------------------------- */

void FixPropertyAtom::do_forward_comm()
{
    timer->stamp();
    if (commGhost) comm->forward_comm_fix(this);
    else error->all(FLERR,"FixPropertyAtom: Faulty implementation - forward_comm invoked, but not registered");
    timer->stamp(TIME_COMM);
}

void FixPropertyAtom::do_reverse_comm()
{
   timer->stamp();
   if (commGhostRev)  comm->reverse_comm_fix(this);
   else error->all(FLERR,"FixPropertyAtom: Faulty implementation - reverse_comm invoked, but not registered");
   timer->stamp(TIME_COMM);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixPropertyAtom::memory_usage()
{
    int nmax = atom->nmax;
    double bytes = nmax * nvalues * sizeof(double);
    return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixPropertyAtom::grow_arrays(int nmax)
{
    if (data_style) memory->grow(array_atom,nmax,nvalues,"FixPropertyAtom:array_atom");
    else memory->grow(vector_atom, nmax, "FixPropertyAtom:vector_atom");
    
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixPropertyAtom::copy_arrays(int i, int j, int delflag)
{
    if (data_style) for(int k=0;k<nvalues;k++) array_atom[j][k] = array_atom[i][k];
    else vector_atom[j]=vector_atom[i];
}

/* ----------------------------------------------------------------------
   called before set_arrays is called for each atom
------------------------------------------------------------------------- */

void FixPropertyAtom::pre_set_arrays()
{
    
    property = 0;
    if(propertyname)
    {
        
        int len1,len2;
        
        property = (double*) atom->get_properties()->find_property(propertyname,"scalar-atom",len1,len2);
        if(!property)
        {
            char errstr[200];
            sprintf(errstr,"Property %s not found",propertyname);
            
            error->fix_error(FLERR,this,errstr);
        }

    }
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixPropertyAtom::set_arrays(int i)
{
    
    if (data_style)
        for(int k=0;k<nvalues;k++)
            array_atom[i][k] = defaultvalues[k];
    else vector_atom[i] = property ? property[i] : defaultvalues[0];
}

/* ----------------------------------------------------------------------
   set all atoms values
------------------------------------------------------------------------- */

void FixPropertyAtom::set_all(double value, bool ghost)
{
    
    int nall;

    if(!ghost)
        nall = atom->nlocal;
    else
        nall = atom->nlocal + atom->nghost;
    if (data_style)
    {
        for(int i = 0; i < nall; i++)
        {
            for(int k=0;k<nvalues;k++)
                array_atom[i][k] = value;
        }
    }
    else
    {
        for(int i = 0; i < nall; i++)
            vector_atom[i] = value;
    }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixPropertyAtom::pack_exchange(int i, double *buf)
{
    if (data_style) for(int k=0;k<nvalues;k++) buf[k] = array_atom[i][k];
    else buf[0] = vector_atom[i];
    return nvalues;
}

/* ----------------------------------------------------------------------
   unpack values into local atom-based arrays after exchange
------------------------------------------------------------------------- */

int FixPropertyAtom::unpack_exchange(int nlocal, double *buf)
{
    if (data_style) for(int k=0;k<nvalues;k++) array_atom[nlocal][k] = buf[k];
    else vector_atom[nlocal]=buf[0];
    return nvalues;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixPropertyAtom::pack_restart(int i, double *buf)
{
  buf[0] = static_cast<double>(nvalues+1);
  if (data_style) for(int k=0;k<nvalues;k++) buf[k+1] = array_atom[i][k];
  else buf[1] = vector_atom[i];

  return (nvalues+1);
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixPropertyAtom::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  if (data_style) for(int k=0;k<nvalues;k++) array_atom[nlocal][k] = extra[nlocal][m++];
  else vector_atom[nlocal] = extra[nlocal][m++];
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixPropertyAtom::maxsize_restart()
{
  return nvalues+1;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixPropertyAtom::size_restart(int nlocal)
{
  return nvalues+1;
}

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

int FixPropertyAtom::pack_comm(int n, int *list, double *buf,
                             int pbc_flag, int *pbc)
{
    int i,j;
    //we dont need to account for pbc here
    int m = 0;
    for (i = 0; i < n; i++) {
      j = list[i];
      if (data_style) for(int k=0;k<nvalues;k++) buf[m++] = array_atom[j][k];
      else buf[m++] = vector_atom[j];
    }
    return nvalues;
}

/* ---------------------------------------------------------------------- */

void FixPropertyAtom::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
      if (data_style) for(int k=0;k<nvalues;k++) array_atom[i][k]=buf[m++];
      else vector_atom[i]=buf[m++];
  }

}

/* ---------------------------------------------------------------------- */

int FixPropertyAtom::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (data_style) for(int k=0;k<nvalues;k++) buf[m++] = array_atom[i][k];
    else buf[m++] = vector_atom[i];
  }
  return nvalues;
}

/* ---------------------------------------------------------------------- */

void FixPropertyAtom::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    if (data_style) for(int k=0;k<nvalues;k++) array_atom[j][k]+=buf[m++];
    else vector_atom[j]+=buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void FixPropertyAtom::write_restart(FILE *fp)
{
  int n = 0;
  double list[1];
  list[n++] = static_cast<double>(nvalues);

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixPropertyAtom::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;
  int nvalues_re;

  nvalues_re = static_cast<int> (list[n++]);

  if(nvalues_re != nvalues)
    error->fix_error(FLERR,this,"restarted fix has incompatible data size");
}

/* ----------------------------------------------------------------------
   return components of property sum on fix group, n = 0..nvalues-1
------------------------------------------------------------------------- */

double FixPropertyAtom::compute_vector(int n)
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  double value = 0.;

  for (int i = 0; i < nlocal; i++)
  {
      if (mask[i] & groupbit)
      {
          if (data_style) value += array_atom[i][n];
          else value += vector_atom[i];
      }
  }

  MPI_Sum_Scalar(value,world);
  return value;
}
