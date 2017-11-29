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
    This file is from LAMMPS
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

#include "lmptype.h"
#include <mpi.h>
#include <stdlib.h>
#include <string>
#include <stdio.h>
#include "dump.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "group.h"
#include "output.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#if !defined(_WINDOWS) && !defined(__MINGW32__)
#include <sys/stat.h>
#endif

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Dump::Dump(LAMMPS *lmp, int narg, char **arg) :
    Pointers(lmp),
    sortBuffer(NULL)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  igroup = group->find(arg[1]);
  groupbit = group->bitmask[igroup];

  n = strlen(arg[2]) + 1;
  style = new char[n];
  strcpy(style,arg[2]);

  n = strlen(arg[4]) + 1;
  filename = new char[n];
  strcpy(filename,arg[4]);

  // check whether the folder is accessible, not available on windows
#if !defined(_WINDOWS) && !defined(__MINGW32__)
    std::string fname(filename);
    std::size_t last_slash = fname.rfind("/");
    // check if we use directories at all
    if (last_slash != std::string::npos)
    {
        std::size_t next_slash = fname.find("/", 1);
        while (next_slash != std::string::npos)
        {
            std::string curdir = fname.substr(0, next_slash);
            struct stat statbuf;
            const bool exists = (stat(curdir.c_str(), &statbuf) != -1) && S_ISDIR(statbuf.st_mode);
            if (!exists)
                mkdir(curdir.c_str(), S_IRWXU | S_IRGRP | S_IXGRP);
            next_slash = fname.find("/", next_slash+1);
        }
    }
#endif

  comm_forward = comm_reverse = 0;

  first_flag = 0;
  flush_flag = 1;
  format = NULL;
  format_user = NULL;
  format_default = NULL;
  clearstep = 0;
  append_flag = 0;
  buffer_allow = 0;
  buffer_flag = 0;
  padflag = 0;

  maxbuf = 0;
  buf = NULL;

  size_one = 0;

  maxsbuf = 0;
  sbuf = NULL;

  // parse filename for special syntax
  // if contains '%', write one file per proc and replace % with proc-ID
  // if contains '*', write one file per timestep and replace * with timestep
  // check file suffixes
  //   if ends in .bin = binary file
  //   else if ends in .gz = gzipped text file
  //   else ASCII text file

  fp = NULL;
  singlefile_opened = 0;
  compressed = 0;
  binary = 0;
  multifile = 0;

  multiproc = 0;
  nclusterprocs = nprocs;
  filewriter = 0;
  if (me == 0) filewriter = 1;
  fileproc = 0;
  multiname = NULL;

  char *ptr;
  if ((ptr = strchr(filename,'%'))) { 
    multiproc = 1;
    nclusterprocs = 1;
    filewriter = 1;
    fileproc = me;
    MPI_Comm_split(world,me,0,&clustercomm);
    multiname = new char[strlen(filename) + 16];
    *ptr = '\0';
    sprintf(multiname,"%s%d%s",filename,me,ptr+1);
    *ptr = '%';
  }

  if (strchr(filename,'*')) multifile = 1;

  char *suffix = filename + strlen(filename) - strlen(".bin");
  if (suffix > filename && strcmp(suffix,".bin") == 0) binary = 1;
  suffix = filename + strlen(filename) - strlen(".gz");
  if (suffix > filename && strcmp(suffix,".gz") == 0) compressed = 1;
}

/* ---------------------------------------------------------------------- */

Dump::~Dump()
{
  delete [] id;
  delete [] style;
  delete [] filename;
  delete [] multiname;

  delete [] format;
  delete [] format_default;
  delete [] format_user;

  memory->destroy(buf);

  memory->destroy(sbuf);

  if (multiproc) MPI_Comm_free(&clustercomm);

  // XTC style sets fp to NULL since it closes file in its destructor

  if (multifile == 0 && fp != NULL) {
    if (compressed) {
      if (filewriter) pclose(fp);
    } else {
      if (filewriter) fclose(fp);
    }
  }
}

/* ---------------------------------------------------------------------- */

void Dump::init()
{
    init_style();

    if (sortBuffer)
    {
        sortBuffer->init(igroup);

        if (sortBuffer->sort_set())
        {
            if (multiproc > 1)
                error->all(FLERR,
                           "Cannot dump sort when multiple procs write the dump file");
            if (sortBuffer->get_sortcol() == 0 && atom->tag_enable == 0)
                error->all(FLERR,"Cannot dump sort on atom IDs with no atom IDs defined");
            if (sortBuffer->get_sortcol() && sortBuffer->get_sortcol() > size_one)
                error->all(FLERR,"Dump sort column is invalid");
        }
    }
}

/* ---------------------------------------------------------------------- */

int Dump::count()
{
  if (igroup == 0) return atom->nlocal;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int m = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) m++;
  return m;
}

/* ---------------------------------------------------------------------- */

void Dump::write()
{
  // if file per timestep, open new file

  if (multifile) openfile();

  // simulation box bounds

  if (domain->triclinic == 0) {
    boxxlo = domain->boxlo[0];
    boxxhi = domain->boxhi[0];
    boxylo = domain->boxlo[1];
    boxyhi = domain->boxhi[1];
    boxzlo = domain->boxlo[2];
    boxzhi = domain->boxhi[2];
  } else {
    boxxlo = domain->boxlo_bound[0];
    boxxhi = domain->boxhi_bound[0];
    boxylo = domain->boxlo_bound[1];
    boxyhi = domain->boxhi_bound[1];
    boxzlo = domain->boxlo_bound[2];
    boxzhi = domain->boxhi_bound[2];
    boxxy = domain->xy;
    boxxz = domain->xz;
    boxyz = domain->yz;
  }

  // nme = # of dump lines this proc contributes to dump

  nme = count();

  // ntotal = total # of dump lines in snapshot
  // nmax = max # of dump lines on any proc

  bigint bnme = nme;
  MPI_Allreduce(&bnme,&ntotal,1,MPI_LMP_BIGINT,MPI_SUM,world);

  int nmax;
  if (multiproc != nprocs) MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,world);
  else nmax = nme;

  // write timestep header
  // for multiproc,
  //   nheader = # of lines in this file via Allreduce on clustercomm

  bigint nheader = ntotal;
  if (multiproc)
    MPI_Allreduce(&bnme,&nheader,1,MPI_LMP_BIGINT,MPI_SUM,clustercomm);

  if (filewriter) write_header(nheader);

  // ensure buf is sized for packing and communicating
  // use nmax to ensure filewriter proc can receive info from others
  // limit nmax*size_one to int since used as arg in MPI calls

  if (nmax > maxbuf) {
    if ((bigint) nmax * size_one > MAXSMALLINT)
      error->all(FLERR,"Too much per-proc info for dump");
    maxbuf = nmax;
    memory->destroy(buf);
    memory->create(buf,maxbuf*size_one,"dump:buf");
  }

  // ensure ids buffer is sized for sorting

  if (sortBuffer)
    sortBuffer->realloc_ids(nmax);

  // pack my data into buf
  // if sorting on IDs also request ID list from pack()
  // sort buf as needed

  if (sortBuffer)
    pack(sortBuffer->get_ids());
  else
    pack(NULL);
  if (sortBuffer)
    sortBuffer->sort(buf, nme, maxbuf, size_one, ntotal);

  // if buffering, convert doubles into strings
  // ensure sbuf is sized for communicating
  // cannot buffer if output is to binary file

  if (buffer_flag && !binary) {
    nsme = convert_string(nme,buf);
    int nsmin,nsmax;
    MPI_Allreduce(&nsme,&nsmin,1,MPI_INT,MPI_MIN,world);
    if (nsmin < 0) error->all(FLERR,"Too much buffered per-proc info for dump");
    if (multiproc != nprocs)
      MPI_Allreduce(&nsme,&nsmax,1,MPI_INT,MPI_MAX,world);
    else nsmax = nsme;
    if (nsmax > maxsbuf) {
      maxsbuf = nsmax;
      memory->grow(sbuf,maxsbuf,"dump:sbuf");
    }
  }

  // filewriter = 1 = this proc writes to file
  //   ping each proc in my cluster, receive its data, write data to file
  // else wait for ping from fileproc, send my data to fileproc

  int tmp,nlines,nchars;
  MPI_Status status;
  MPI_Request request;

  // comm and output buf of doubles

  if (buffer_flag == 0 || binary)
  {
    if (filewriter)
    {
        for (int iproc = 0; iproc < nclusterprocs; iproc++)
        {
            if (iproc)
            {
                MPI_Irecv(buf,maxbuf*size_one,MPI_DOUBLE,me+iproc,0,world,&request);
                MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
                MPI_Wait(&request,&status);
                MPI_Get_count(&status,MPI_DOUBLE,&nlines);
                nlines /= size_one;
            }
            else
                nlines = nme;

            write_data(nlines,buf);
        }
        if (flush_flag)
            fflush(fp);

    }
    else
    {
        MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,&status);
        MPI_Rsend(buf,nme*size_one,MPI_DOUBLE,fileproc,0,world);
    }

  // comm and output sbuf = one big string of formatted values per proc

  } else {
    if (filewriter) {
      for (int iproc = 0; iproc < nclusterprocs; iproc++) {
        if (iproc) {
          MPI_Irecv(sbuf,maxsbuf,MPI_CHAR,me+iproc,0,world,&request);
          MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
          MPI_Wait(&request,&status);
          MPI_Get_count(&status,MPI_CHAR,&nchars);
        } else nchars = nsme;

        write_data(nchars,(double *) sbuf);
      }
      if (flush_flag) fflush(fp);

    } else {
      MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,&status);
      MPI_Rsend(sbuf,nsme,MPI_CHAR,fileproc,0,world);
    }
  }

  // if file per timestep, close file if I am filewriter

  if (multifile) {
    if (compressed) {
      if (filewriter) pclose(fp);
    } else {
      if (filewriter) fclose(fp);
    }
  }
}

/* ----------------------------------------------------------------------
   generic opening of a dump file
   ASCII or binary or gzipped
   some derived classes override this function
------------------------------------------------------------------------- */

void Dump::openfile()
{
  // single file, already opened, so just return

  if (singlefile_opened) return;
  if (multifile == 0) singlefile_opened = 1;

  // if one file per timestep, replace '*' with current timestep

  char *filecurrent = filename;
  if (multiproc) filecurrent = multiname;

  if (multifile) {
    char *filestar = filecurrent;
    filecurrent = new char[strlen(filestar) + 16];
    char *ptr = strchr(filestar,'*');
    *ptr = '\0';
    if (padflag == 0)
      sprintf(filecurrent,"%s" BIGINT_FORMAT "%s",
              filestar,update->ntimestep,ptr+1);
    else {
      char bif[8],pad[16];
      strcpy(bif,BIGINT_FORMAT);
      sprintf(pad,"%%s%%0%d%s%%s",padflag,&bif[1]);
      sprintf(filecurrent,pad,filestar,update->ntimestep,ptr+1);
    }
    *ptr = '*';
  }

  // each proc with filewriter = 1 opens a file

  if (filewriter) {
    if (compressed) {
#ifdef LAMMPS_GZIP
      char gzip[128];
      sprintf(gzip,"gzip -6 > %s",filecurrent);
#ifdef _WIN32
      fp = _popen(gzip,"wb");
#else
      fp = popen(gzip,"w");
#endif
#else
      error->one(FLERR,"Cannot open gzipped file");
#endif
    } else if (binary) {
      fp = fopen(filecurrent,"wb");
      
    } else if (append_flag) {
      fp = fopen(filecurrent,"a");
    } else {
      fp = fopen(filecurrent,"w");
      
    }

    if (fp == NULL) error->one(FLERR,"Cannot open dump file");
  } else fp = NULL;

  // delete string with timestep replaced

  if (multifile) delete [] filecurrent;
}

/* ----------------------------------------------------------------------
   process params common to all dumps here
   if unknown param, call modify_param specific to the dump
------------------------------------------------------------------------- */

void Dump::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all(FLERR,"Illegal dump_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"append") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) append_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) append_flag = 0;
      else error->all(FLERR,"Illegal dump_modify command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"buffer") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) buffer_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) buffer_flag = 0;
      else error->all(FLERR,"Illegal dump_modify command");
      if (buffer_flag && buffer_allow == 0)
        error->all(FLERR,"Dump_modify buffer yes not allowed for this style");
      iarg += 2;

    } else if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      int idump;
      for (idump = 0; idump < output->ndump; idump++)
        if (strcmp(id,output->dump[idump]->id) == 0) break;
      int n;
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        delete [] output->var_dump[idump];
        n = strlen(&arg[iarg+1][2]) + 1;
        output->var_dump[idump] = new char[n];
        strcpy(output->var_dump[idump],&arg[iarg+1][2]);
        n = 0;
      } else {
        n = force->inumeric(FLERR,arg[iarg+1]);
        if (n <= 0) error->all(FLERR,"Illegal dump_modify command");
      }
      output->every_dump[idump] = n;
      iarg += 2;

    } else if (strcmp(arg[iarg],"first") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) first_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) first_flag = 0;
      else error->all(FLERR,"Illegal dump_modify command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"fileper") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      if (!multiproc)
        error->all(FLERR,"Cannot use dump_modify fileper "
                   "without % in dump file name");
      int nper = force->inumeric(FLERR,arg[iarg+1]);
      if (nper <= 0) error->all(FLERR,"Illegal dump_modify command");

      multiproc = nprocs/nper;
      if (nprocs % nper) multiproc++;
      fileproc = me/nper * nper;
      int fileprocnext = MIN(fileproc+nper,nprocs);
      nclusterprocs = fileprocnext - fileproc;
      if (me == fileproc) filewriter = 1;
      else filewriter = 0;
      int icluster = fileproc/nper;

      MPI_Comm_free(&clustercomm);
      MPI_Comm_split(world,icluster,0,&clustercomm);

      delete [] multiname;
      multiname = new char[strlen(filename) + 16];
      char *ptr = strchr(filename,'%');
      if (ptr)
      {
          *ptr = '\0';
          sprintf(multiname,"%s%d%s",filename,icluster,ptr+1);
          *ptr = '%';
      }
      iarg += 2;

    } else if (strcmp(arg[iarg],"flush") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) flush_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) flush_flag = 0;
      else error->all(FLERR,"Illegal dump_modify command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"format") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      delete [] format_user;
      format_user = NULL;
      if (strcmp(arg[iarg+1],"none")) {
        int n = strlen(arg[iarg+1]) + 1;
        format_user = new char[n];
        strcpy(format_user,arg[iarg+1]);
      }
      iarg += 2;

    } else if (strcmp(arg[iarg],"nfile") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      if (!multiproc)
        error->all(FLERR,"Cannot use dump_modify nfile "
                   "without % in dump file name");
      int nfile = force->inumeric(FLERR,arg[iarg+1]);
      if (nfile <= 0) error->all(FLERR,"Illegal dump_modify command");
      nfile = MIN(nfile,nprocs);

      multiproc = nfile;
      int icluster = static_cast<int> ((bigint) me * nfile/nprocs);
      fileproc = static_cast<int> ((bigint) icluster * nprocs/nfile);
      int fcluster = static_cast<int> ((bigint) fileproc * nfile/nprocs);
      if (fcluster < icluster) fileproc++;
      int fileprocnext =
        static_cast<int> ((bigint) (icluster+1) * nprocs/nfile);
      fcluster = static_cast<int> ((bigint) fileprocnext * nfile/nprocs);
      if (fcluster < icluster+1) fileprocnext++;
      nclusterprocs = fileprocnext - fileproc;
      if (me == fileproc) filewriter = 1;
      else filewriter = 0;

      MPI_Comm_free(&clustercomm);
      MPI_Comm_split(world,icluster,0,&clustercomm);

      delete [] multiname;
      multiname = new char[strlen(filename) + 16];
      char *ptr = strchr(filename,'%');
      if (ptr)
      {
          *ptr = '\0';
          sprintf(multiname,"%s%d%s",filename,icluster,ptr+1);
          *ptr = '%';
      }
      iarg += 2;

    } else if (strcmp(arg[iarg],"pad") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      padflag = force->inumeric(FLERR,arg[iarg+1]);
      if (padflag < 0) error->all(FLERR,"Illegal dump_modify command");
      iarg += 2;

    } else {
      int n = modify_param(narg-iarg,&arg[iarg]);
      if (n == 0)
      {
        if (!sortBuffer)
            sortBuffer = new SortBuffer(lmp, false);
        n = sortBuffer->modify_param(narg-iarg, &arg[iarg]);
        if (n == 0)
            error->all(FLERR,"Illegal dump_modify command");
      }
      iarg += n;
    }
  }
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint Dump::memory_usage()
{
  bigint bytes = 0;
  bytes += memory->usage(buf,size_one*maxbuf);
  bytes += memory->usage(sbuf,maxsbuf);
  if (sortBuffer) {
    bytes += sortBuffer->memory_usage(size_one);
  }
  return bytes;
}
