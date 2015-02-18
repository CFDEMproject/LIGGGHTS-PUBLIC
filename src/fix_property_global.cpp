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
#include "neighbor.h"
#include "fix_property_global.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define EPSILON 0.001

/* ---------------------------------------------------------------------- */

FixPropertyGlobal::FixPropertyGlobal(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
    //Check args
    if (narg < 6) error->all(FLERR,"Illegal fix property/global command, not enough arguments");

    //Read args
    int n = strlen(arg[3]) + 1;
    variablename = new char[n];
    strcpy(variablename,arg[3]);
    is_symmetric = false;
    is_atomtype_bound = false;

    if (strcmp(arg[4],"scalar") == 0)
        data_style = FIXPROPERTY_GLOBAL_SCALAR;
    else if (strcmp(arg[4],"vector") == 0)
        data_style = FIXPROPERTY_GLOBAL_VECTOR;
    else if (strcmp(arg[4],"peratomtype") == 0 || strcmp(arg[4],"atomtype") == 0)
    {
        data_style = FIXPROPERTY_GLOBAL_VECTOR;
        is_atomtype_bound = true;
    }
    else if (strcmp(arg[4],"matrix") == 0)
        data_style = FIXPROPERTY_GLOBAL_MATRIX;
    else if (strcmp(arg[4],"peratomtypepair") == 0 || strcmp(arg[4],"atomtypepair") == 0)
    {
        data_style = FIXPROPERTY_GLOBAL_MATRIX;
        is_symmetric = true;
        is_atomtype_bound = true;
    }
    else error->fix_error(FLERR,this,"Unknown style. Valid styles are scalar, vector, atomtype/peratomtype, matrix, or atomtypepair/peratomtypepair");

    int darg = 0;
    if (data_style == FIXPROPERTY_GLOBAL_MATRIX) darg = 1;

    //assign values
    nvalues = narg - 5 - darg;
    nvalues_new_array = 0;
    
    values = (double*) memory->smalloc(nvalues*sizeof(double),"values");
    values_recomputed = (double*) memory->smalloc(nvalues*sizeof(double),"values");

    if(narg < 5+darg+nvalues) error->fix_error(FLERR,this,"not enough arguments");

    for (int j = 0; j < nvalues; j++)
        values[j] = force->numeric(FLERR,arg[5+darg+j]);

    if (data_style == FIXPROPERTY_GLOBAL_SCALAR)
        scalar_flag = 1;
    else if (data_style==FIXPROPERTY_GLOBAL_VECTOR) {
        vector_flag = 1;
        size_vector = nvalues;
    }
    else if (data_style == FIXPROPERTY_GLOBAL_MATRIX) {
        array_flag = 1;
        size_array_cols = force->inumeric(FLERR,arg[5]);
        if (fmod(static_cast<double>(nvalues),size_array_cols) != 0.)
          error->fix_error(FLERR,this,"the number of default values must be a multiple of nCols.");
        size_array_rows = static_cast<int>(static_cast<double>(nvalues)/size_array_cols);
    }

    extvector=0; 

    filename = 0;
    grpname = 0;

    //check if there is already a fix that tries to register a property with the same name
    
    for (int ifix = 0; ifix < modify->nfix; ifix++)
        if ((modify->fix[ifix]) && (strcmp(modify->fix[ifix]->style,style) == 0) && (strcmp(((FixPropertyGlobal*)(modify->fix[ifix]))->variablename,variablename)==0) )
            error->fix_error(FLERR,this,"There is already a fix that registers a variable of the same name");

    array = NULL;
    array_recomputed = NULL;
    if(data_style == FIXPROPERTY_GLOBAL_MATRIX)
    {
        array = (double**)memory->smalloc(size_array_rows*sizeof(double**),"FixPropGlob:array");
        array_recomputed = (double**)memory->smalloc(size_array_rows*sizeof(double**),"FixPropGlob:array_recomputed");
        for(int i = 0; i < size_array_rows; i++) array[i] = &values[i*size_array_cols];
        for(int i = 0; i < size_array_rows; i++) array_recomputed[i] = &values_recomputed[i*size_array_cols];
    }

    // error check if matrix is symmetric (if required)
    if(is_symmetric)
    {
        if(size_array_rows != size_array_cols)
            error->fix_error(FLERR,this,"per-atomtype property matrix must be symmetric, i.e. N atom types "
                                        "require you to define N columns and N rows with N*N total values");

        int sflag = true;
        for(int i = 0; i < size_array_rows; i++)
            for(int j = 0; j < size_array_cols; j++)
                if(array[i][j] != array[j][i])
                    sflag = false;

        if(!sflag)
            error->fix_error(FLERR,this,"per-atomtype property matrix must be symmetric");
    }
}

/* ---------------------------------------------------------------------- */

FixPropertyGlobal::~FixPropertyGlobal()
{
  // delete locally stored arrays
  delete[] variablename;

  if(filename) delete[] filename;
  if(grpname) delete[] grpname;

  memory->sfree(values);
  memory->sfree(values_recomputed);

  if(array)            memory->sfree(array);
  if(array_recomputed) memory->sfree(array_recomputed);
}

/* ---------------------------------------------------------------------- */

void FixPropertyGlobal::pre_delete(bool unfixflag)
{
    if(filename) write();
}

/* ---------------------------------------------------------------------- */

Fix* FixPropertyGlobal::check_fix(const char *varname,const char *svmstyle,int len1,int len2,const char *caller,bool errflag)
{
    char errmsg[400];

    if(strcmp(varname,variablename) == 0)
    {
        if(strcmp(svmstyle,"scalar") == 0) len1 = 1;

        // check variable style
        if(
            (strcmp(svmstyle,"scalar") == 0 && data_style != FIXPROPERTY_GLOBAL_SCALAR) ||
            ((strcmp(svmstyle,"vector") == 0 || strcmp(svmstyle,"peratomtype") == 0) && data_style != FIXPROPERTY_GLOBAL_VECTOR) ||
            ((strcmp(svmstyle,"matrix") == 0 || strcmp(svmstyle,"peratomtypepair") == 0) && data_style != FIXPROPERTY_GLOBAL_MATRIX)
        )
        {
            if(errflag)
            {
                sprintf(errmsg,"%s style required for fix property/global variable %s for usage with %s",svmstyle,varname,caller);
                error->fix_error(FLERR,this,errmsg);
            }
            else return NULL;
        }

        // check length
        if((nvalues < len1) && ((data_style != FIXPROPERTY_GLOBAL_MATRIX) || ((data_style == FIXPROPERTY_GLOBAL_MATRIX) && (size_array_cols < len2))))
        {
            if(errflag)
            {
                sprintf(errmsg,"Length not sufficient for variable %s for usage with %s",varname,caller);
                error->fix_error(FLERR,this,errmsg);
            }
            else return NULL;
        }

        // success
        return static_cast<Fix*>(this);
    }
    return NULL;
}

/* ---------------------------------------------------------------------- */

void FixPropertyGlobal::init()
{
    me = comm->me;

    char errmsg[300];
    int ntypes = atom->ntypes;

    if(FIXPROPERTY_GLOBAL_VECTOR == data_style && is_atomtype_bound && nvalues != ntypes)
    {
        
        sprintf(errmsg,"Fix property/global: Length not correct for variable %s, length should be equal to %d (= number of atom types)",
                variablename,ntypes);
        error->fix_error(FLERR,this,errmsg);
    }
    if(FIXPROPERTY_GLOBAL_MATRIX == data_style && is_atomtype_bound && nvalues != ntypes*ntypes)
    {
        sprintf(errmsg,"Fix property/global: Length not correct for variable %s, length should be equal to %d ( = number of atom types * number of atom types)",
                variablename,ntypes*ntypes);
        error->fix_error(FLERR,this,errmsg);
    }
}

/* ---------------------------------------------------------------------- */

void FixPropertyGlobal::grow(int len1, int len2)
{
    if(data_style == FIXPROPERTY_GLOBAL_SCALAR) error->fix_error(FLERR,this,"Can not grow global property of type scalar");
    else if(data_style == FIXPROPERTY_GLOBAL_VECTOR && len1 > nvalues)
    {
        memory->grow(values,len1,"FixPropertyGlobal:values");
    }
    else if(data_style == FIXPROPERTY_GLOBAL_MATRIX && len1*len2 > nvalues)
    {
        values = (double*) memory->srealloc(values,len1*len2*sizeof(double),"FixPropertyGlobal:values");
        size_array_rows = len1;
        size_array_cols = len2;
        nvalues = len1*len2;
        array = (double**)memory->srealloc(array,size_array_rows*sizeof(double**),"FixPropGlob:array");
        for(int i = 0; i < size_array_rows; i++) array[i] = &values[i*size_array_cols];
    }
}

/* ---------------------------------------------------------------------- */

double FixPropertyGlobal::compute_scalar()
{
  return values[0];
}

/* ---------------------------------------------------------------------- */

double FixPropertyGlobal::compute_vector(int i)
{
    if (i>(nvalues-1))error->fix_error(FLERR,this,"Trying to access vector, but index out of bounds");
    return values[i];
}

void FixPropertyGlobal::vector_modify(int i,double val)
{
    if (i>(nvalues-1))error->fix_error(FLERR,this,"Trying to access vector, but index out of bounds");
    values_recomputed[i] = val;
}

double FixPropertyGlobal::compute_vector_modified(int i)
{
    if (i>(nvalues-1))error->fix_error(FLERR,this,"Trying to access vector, but index out of bounds");
    return values_recomputed[i];
}

/* ---------------------------------------------------------------------- */

double FixPropertyGlobal::compute_array(int i, int j) //i is row, j is column
{
    if (i>(size_array_rows-1))error->fix_error(FLERR,this,"Trying to access matrix, but row index out of bounds");
    if (j>(size_array_cols-1))error->fix_error(FLERR,this,"Trying to access matrix, but column index out of bounds");

    return array[i][j];
}

void FixPropertyGlobal::array_modify(int i, int j,double val) //i is row, j is column
{
    if (i>(size_array_rows-1))error->fix_error(FLERR,this,"Trying to access matrix, but row index out of bounds");
    if (j>(size_array_cols-1))error->fix_error(FLERR,this,"Trying to access matrix, but column index out of bounds");

    array_recomputed[i][j] = val;
}

double FixPropertyGlobal::compute_array_modified(int i, int j) //i is row, j is column
{
    if (i>(size_array_rows-1))error->fix_error(FLERR,this,"Trying to access matrix, but row index out of bounds");
    if (j>(size_array_cols-1))error->fix_error(FLERR,this,"Trying to access matrix, but column index out of bounds");

    return array_recomputed[i][j];
}

/* ---------------------------------------------------------------------- */

int FixPropertyGlobal::setmask()
{
  int mask = 0;
  return mask;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixPropertyGlobal::memory_usage()
{
  double bytes = nvalues * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

void FixPropertyGlobal::new_array(int l1,int l2)
{
    
    if (data_style == FIXPROPERTY_GLOBAL_MATRIX) error->fix_error(FLERR,this,"Can not allocate extra array for matrix style");
    array_flag = 1;
    size_array_rows = l1;
    size_array_cols = l2;
    nvalues_new_array = l1*l2;

    memory->create(array,size_array_rows,size_array_cols,"FixPropGlob:array");
    memory->create(array_recomputed,size_array_rows,size_array_cols,"FixPropGlob:array_recomputed");
}

/* ----------------------------------------------------------------------
   write out command
------------------------------------------------------------------------- */

void FixPropertyGlobal::write()
{
    
    if(0 != me)
        return;

    FILE *file = fopen(filename,"w");

    if(!file)
        error->one(FLERR,"Fix property/global cannot open file");

    // fix id group style variablename
    fprintf(file,"fix %s %s %s %s ",id,grpname,style,variablename);

    // datatype
    char *datatyp;
    if(0 == data_style) datatyp = (char*) "scalar";
    if(1 == data_style) datatyp = (char*) "vector";
    if(2 == data_style && is_symmetric) datatyp = (char*) "atomtypepair";
    else if(2 == data_style) datatyp            = (char*) "matrix";
    fprintf(file,"%s ",datatyp);

    // size_array_cols if required
    if(2 == data_style) fprintf(file,"%d ",size_array_cols);

    // values
    for(int i = 0; i < nvalues; i++)
        fprintf(file,"%f ",values[i]);

    fprintf(file,"\n");
    fclose(file);
}

/* ---------------------------------------------------------------------- */

int FixPropertyGlobal::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"file") == 0) {
    if (narg < 2) error->fix_error(FLERR,this,"not enough arguments for fix_modify 'file'");

    filename = new char[strlen(arg[1])+1];
    strcpy(filename,arg[1]);
    grpname = new char[strlen(group->names[igroup])+1];
    strcpy(grpname,group->names[igroup]);
    return 2;
  }

  return 0;
}
