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

#include "sys/stat.h"
#include "string.h"
#include "stdlib.h"
#include "atom.h"
#include "comm.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "math.h"
#include "vector_liggghts.h"
#include "fix_property_atom.h"
#include "fix_property_global.h"
#include "fix_cfd_coupling.h"
#include "cfd_datacoupling_file.h"
#include <iostream>
#include <fstream>
#include <unistd.h>

#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#define sleep Sleep
#endif

using namespace LAMMPS_NS;
using namespace std;

/* ---------------------------------------------------------------------- */

CfdDatacouplingFile::CfdDatacouplingFile(LAMMPS *lmp, int iarg,int narg, char **arg,FixCfdCoupling *fc)  :
  CfdDatacoupling(lmp, iarg, narg, arg,fc)
{
    iarg_ = iarg;
    int n_arg = narg - iarg_;
    
    if(n_arg < 1) error->all(FLERR,"Cfd file coupling: wrong # arguments");

    liggghts_is_active = true;
    firstexec = true;

    this->fc_ = fc;
    is_parallel = false;

    filepath = new char[strlen(arg[iarg_])+2];
    strcpy(filepath,arg[iarg_]);

    t0 = -1;

    iarg_++;

    append = 1;
}

/* ---------------------------------------------------------------------- */

CfdDatacouplingFile::~CfdDatacouplingFile()
{
    delete []filepath;
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingFile::post_create()
{
    if(!is_parallel && comm->nprocs > 1)  error->all(FLERR,"Fix couple/cfd with file coupling is for serial computation only");
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingFile::exchange()
{
    void *dummy = NULL;

    // write to file
    for(int i = 0; i < npush_; i++)
    {
       push(pushnames_[i],pushtypes_[i],dummy,"");
    }

    // read from files
    for(int i = 0; i < npull_; i++)
    {
       pull(pullnames_[i],pulltypes_[i],dummy,"");
    }
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingFile::pull(const char *name, const char *type, void *&from, const char *datatype)
{
    CfdDatacoupling::pull(name,type,from,datatype);

    int len1 = -1, len2 = -1;

    void * to = find_pull_property(name,type,len1,len2);

    if(to && strcmp(type,"scalar-atom") == 0)
    {
        readScalarData(name,(double*)to);
    }

    else if(to && strcmp(type,"vector-atom") == 0)
    {
        readVectorData(name,(double**)to);
    }

    else if(to &&  strcmp(type,"vector-global") == 0)
    {
        readGlobalVectorData(name,(double*)from,len1);
    }

    else if(to && strcmp(type,"array-global") == 0)
    {
        readGlobalArrayData(name,(double**)from,len1,len2);
    }

    else
    {
        if(screen) fprintf(screen,"LIGGGHTS could not find property %s to write data from calling program to.\n",name);
        lmp->error->all(FLERR,"This error is fatal");
    }
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingFile::push(const char *name, const char *type, void *&to, const char *datatype)
{
    CfdDatacoupling::push(name,type,to,datatype);

    int len1 = -1, len2 = -1;

    if(t0 == -1) t0 = update->ntimestep;
    if(update->ntimestep > t0) firstexec = false;

    void * from = find_push_property(name,type,len1,len2);

    if(from && strcmp(type,"scalar-atom") == 0)
    {
        writeScalarData(name,(double*)from);
    }

    else if(from && strcmp(type,"vector-atom") == 0)
    {
        writeVectorData(name,(double**)from);
    }

    else if(from && strcmp(type,"vector-global") == 0)
    {
        writeGlobalVectorData(name,(double*)from,len1);
    }

    else if(from && strcmp(type,"array-global") == 0)
    {
        writeGlobalArrayData(name,(double**)from,len1,len2);
    }

    else
    {
        if(screen) fprintf(screen,"LIGGGHTS could not find property %s to write to calling program.\n",name);
        lmp->error->all(FLERR,"This error is fatal");
    }
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

char * CfdDatacouplingFile::getFilePath(const char *name, bool flag)
{
    
    if(!append) {
      
      char *file = new char[strlen(name)+1];
      strcpy(file,name);
      return file;
    }

    char *file = new char[strlen(filepath)+strlen(name)+3];
    strcpy(file,filepath);
    strcat(file,name);
    if(flag) strcat(file,"0");
    else strcat(file,"1");
    return file;
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingFile::op_complete(const char *name)
{
    if(!append) return;
    char *oldfile = getFilePath(name,true);
    char *newfile = getFilePath(name,false);
    
    rename(oldfile,newfile);
    delete []oldfile;
    delete []newfile;
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingFile::readVectorData(const char *name, double ** field)
{
    // get output path
    char *file = getFilePath(name,true);

    fprintf(screen,"Fix couple/cfd/file: waiting for file: %s\n",file);
    struct stat st;
    while (stat(file,&st)) sleep(10);

    // set file pointer
    ifstream inputPtr(file);

    // skip lines starting with #
    while(inputPtr.peek() == '#')  inputPtr.ignore(1000,'\n');

    // write data to variable
    int numberOfParticles;
    inputPtr >> numberOfParticles;
    
    if(atom->nlocal!=numberOfParticles) error->all(FLERR,"Fix couple/cfd/file: Data corruption: # particles in file does not match # particles in LIGGGHTS.\n"
                                                   "Note that file-based coupling currently does not support inserting or deleting particles during a coupled run.");

    for(int index = 0;index < numberOfParticles; ++index)
    {
        for(int i=0;i<3;i++) inputPtr >> field[index][i];
    }

    // clean up inputStream
    delete []file;
    op_complete(name);
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingFile::readScalarData(const char* name, double *field)
{
    // get output path
    char *file = getFilePath(name,true);

    fprintf(screen,"Fix couple/cfd/file: waiting for file: %s\n",file);
    struct stat st;
    while (stat(file,&st)) sleep(10);

    // set file pointer
    ifstream inputPtr(file);

    // skip lines starting with #
    while(inputPtr.peek() == '#')  inputPtr.ignore(1000,'\n');

    // write data to variable
    int numberOfParticles;
    inputPtr >> numberOfParticles;

    if(atom->nlocal!=numberOfParticles) error->all(FLERR,"Fix couple/cfd/file: Data corruption: # particles in file does not match # particles in LIGGGHTS.\n"
                                                   "Note that file-based coupling currently does not support inserting or deleting particles during a coupled run.");

    // write data to variable
    for(int index = 0;index < numberOfParticles; ++index)
    {
        inputPtr >> field[index];
    }

    // clean up inputStream
    delete []file;
    op_complete(name);
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingFile::readGlobalArrayData(const char *name, double ** field, int &len1, int &len2)
{
    // get output path
    char *file = getFilePath(name,true);

    fprintf(screen,"Fix couple/cfd/file: waiting for file: %s\n",file);
    struct stat st;
    while (stat(file,&st)) sleep(10);

    // set file pointerfrom
    ifstream inputPtr(file);

    // skip lines starting with #
    while(inputPtr.peek() == '#')  inputPtr.ignore(1000,'\n');

    // write data to variable
    int l1,l2;
    inputPtr >> l1;
    inputPtr >> l2;

    if(l1 != len1 || l2 != len2)
        error->one(FLERR,"Global array received has different length than the corresponding global array in LIGGGHTS");

    for(int index = 0; index < len1; ++index)
    {
        for(int i = 0; i < len2; i++)
        {
            if(inputPtr.eof())
                error->one(FLERR,"Global array received has different length than the corresponding global array in LIGGGHTS");
            inputPtr >> field[index][i];
            
        }
    }

    // clean up inputStream
    delete []file;
    op_complete(name);
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingFile::readGlobalVectorData(const char* name, double *field, int &len)
{
    // get output path
    char *file = getFilePath(name,true);

    fprintf(screen,"Fix couple/cfd/file: waiting for file: %s\n",file);
    struct stat st;
    while (stat(file,&st)) sleep(10);

    // set file pointer
    int l1;
    ifstream inputPtr(file);

    // skip lines starting with #
    while(inputPtr.peek() == '#')  inputPtr.ignore(1000,'\n');

    inputPtr >> l1;

    if(l1 != len) error->all(FLERR,"Global vector received has different length than the corresponding global array in LIGGGHTS");

    // write data to variable
    for(int index = 0;index < len; ++index)
    {
        inputPtr >> field[index];
    }

    // clean up inputStream
    delete []file;
    op_complete(name);
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingFile::writeVectorData(const char *name,  double ** field)
{
    // get output path
    char *file = getFilePath(name,true);

    if(!firstexec)
    {
      fprintf(screen,"Fix couple/cfd/file: waiting for file: %s\n",file);
       struct stat st;
       while (stat(file,&st)) sleep(10);
    }

    // set file pointer
    ofstream outputPtr(file);

    // write data to file
    int numberOfParticles = atom->nlocal;
    outputPtr << numberOfParticles << endl;
    for(int index = 0;index < numberOfParticles; ++index)
    {
        for(int i=0;i<3;i++) outputPtr << field[index][i] << " ";
        outputPtr << endl;
    }

    // clean up outputStream and rename file
    delete []file;
    op_complete(name);
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingFile::writeScalarData(const char* name, double * field)
{
    // get output path
    char *file = getFilePath(name,true);

    if(!firstexec)
    {
      fprintf(screen,"Fix couple/cfd/file: waiting for file: %s\n",file);
       struct stat st;
       while (stat(file,&st)) sleep(10);
    }

    // set file pointer
    ofstream outputPtr(file);

    // write data to file
    int numberOfParticles = atom->nlocal;
    outputPtr << numberOfParticles << endl;
    for(int index = 0;index < numberOfParticles; ++index)
    {
        outputPtr << field[index] << endl;
    }

    // clean up outputStream and rename file
    delete []file;
    op_complete(name);
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingFile::writeGlobalVectorData(const char *name,  double *field, int len)
{
    if(len < 0) error->all(FLERR,"Internal error in CfdDatacouplingFile");

    // get output path
    char *file = getFilePath(name,true);

    if(!firstexec)
    {
      fprintf(screen,"Fix couple/cfd/file: waiting for file: %s\n",file);
       struct stat st;
       while (stat(file,&st)) sleep(10);
    }

    // set file pointer
    ofstream outputPtr(file);

    // write data to file
    outputPtr << len << endl;
    for(int index = 0;index < len; ++index)
    {
        outputPtr << field[index];
        outputPtr << endl;
    }

    // clean up outputStream and rename file
    delete []file;
    op_complete(name);
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingFile::writeGlobalArrayData(const char* name, double **field, int len1, int len2)
{
    if(len1 < 0 || len2 < 0) error->all(FLERR,"Internal error in CfdDatacouplingFile");

    // get output path
    char *file = getFilePath(name,true);

    if(!firstexec)
    {
      fprintf(screen,"Fix couple/cfd/file: waiting for file: %s\n",file);
       struct stat st;
       while (stat(file,&st)) sleep(10);
    }

    // set file pointer
    ofstream outputPtr(file);

    // write data to file
    outputPtr << len1 << endl;
    outputPtr << len2 << endl;
    for(int index = 0;index < len1; ++index)
    {
        for(int i=0;i<len2;i++) outputPtr << field[index][i] << " ";
        outputPtr << endl;
    }

    // clean up outputStream and rename file
    delete []file;
    op_complete(name);
}
