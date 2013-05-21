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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "group.h"
#include "timer.h"
#include "neighbor.h"
#include "fix_property_atom.h"
#include "mpi_liggghts.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define EPSILON 0.001
#define myAtof force->numeric 

/* ---------------------------------------------------------------------- */

FixPropertyAtom::FixPropertyAtom(LAMMPS *lmp, int narg, char **arg, bool parse) :
  Fix(lmp, narg, arg)
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

    if (strcmp(arg[4],"scalar") == 0) data_style = FIXPROPERTY_ATOM_SCALAR;
    else if (strcmp(arg[4],"vector") == 0) data_style = FIXPROPERTY_ATOM_VECTOR;
    else error->all(FLERR,"Unknown style for fix property/atom. Valid styles are 'scalar' or 'vector'");

    if (strcmp(arg[5],"yes") == 0) restart_peratom = 1;
    else if (strcmp(arg[5],"no") == 0) restart_peratom = 0;
    else error->all(FLERR,"Unknown restart style for fix property/atom. Valid styles are 'yes' or 'no'");

    if (strcmp(arg[6],"yes") == 0) commGhost = 1;
    else if (strcmp(arg[6],"no") == 0) commGhost = 0;
    else error->all(FLERR,"Unknown communicate_ghost style for fix property/atom. Valid styles are 'yes' or 'no'");

    if (strcmp(arg[7],"yes") == 0) commGhostRev = 1;
    else if (strcmp(arg[7],"no") == 0) commGhostRev = 0;
    else error->all(FLERR,"Unknown communicate_reverse_ghost style for fix property/atom. Valid styles are 'yes' or 'no'");

    nvalues = narg - 8;
    if ((nvalues == 1) && (data_style != FIXPROPERTY_ATOM_SCALAR))
      error->all(FLERR,"Error in fix property/atom: Number of default values provided not consistent with vector style. Provide more than 1 value or use style 'scalar'");

    defaultvalues = new double[nvalues];

    // fix handles properties that need to be initialized at particle creation
    create_attribute = 1;
    for (int j = 0; j < nvalues; j++)
    {
        // if any of the values is none, this fix will not init properties
        if(strcmp(arg[8+j],"none") == 0)
        {
            create_attribute = 0;
            continue;
        }
        defaultvalues[j] = myAtof(arg[8+j]);
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

    // check if there is already a fix that tries to register a property with the same name
    for (int ifix = 0; ifix < modify->nfix; ifix++)
        if ((strcmp(modify->fix[ifix]->style,style) == 0) && (strcmp(((FixPropertyAtom*)(modify->fix[ifix]))->variablename,variablename)==0) )
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
  delete[] variablename;
  delete[] defaultvalues;

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
            (strcmp(svmstyle,"vector") == 0 && data_style != FIXPROPERTY_ATOM_VECTOR)
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

void FixPropertyAtom::copy_arrays(int i, int j)
{
    if (data_style) for(int k=0;k<nvalues;k++) array_atom[j][k] = array_atom[i][k];
    else vector_atom[j]=vector_atom[i];
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixPropertyAtom::set_arrays(int i)
{
    
    if (data_style)
        for(int k=0;k<nvalues;k++)
            array_atom[i][k] = defaultvalues[k];
    else vector_atom[i] = defaultvalues[0];
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
