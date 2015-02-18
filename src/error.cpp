/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   This file was modified with respect to the release in LAMMPS
   Modifications are Copyright 2009-2012 JKU Linz
                     Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
   the producer of the LIGGGHTS® software and the CFDEM®coupling software
   See http://www.cfdem.com/terms-trademark-policy for details.

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdlib.h"
#include "error.h"
#include "universe.h"
#include "output.h"
#include "fix.h" 
#include "force.h" 
#include "compute.h" 
#include "string.h" 

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Error::Error(LAMMPS *lmp) : Pointers(lmp) {}

/* ----------------------------------------------------------------------
   called by all procs in universe
   close all output, screen, and log files in world and universe
   no abort, so insure all procs in universe call, else will hang
------------------------------------------------------------------------- */

void Error::universe_all(const char *file, int line, const char *str)
{
  MPI_Barrier(universe->uworld);

  if (universe->me == 0) {
    if (universe->uscreen) fprintf(universe->uscreen,
                                   "ERROR: %s (%s:%d)\n",str,file,line);
    if (universe->ulogfile) fprintf(universe->ulogfile,
                                    "ERROR: %s (%s:%d)\n",str,file,line);
  }

  if (output) delete output;
  if (universe->nworlds > 1) {
    if (screen && screen != stdout) fclose(screen);
    if (logfile) fclose(logfile);
  }
  if (universe->ulogfile) fclose(universe->ulogfile);

  MPI_Finalize();
  exit(1);
}

/* ----------------------------------------------------------------------
   called by one proc in universe
   forces abort of entire universe if any proc in universe calls
------------------------------------------------------------------------- */

void Error::universe_one(const char *file, int line, const char *str)
{
  if (universe->uscreen)
    fprintf(universe->uscreen,"ERROR on proc %d: %s (%s:%d)\n",
            universe->me,str,file,line);
  MPI_Abort(universe->uworld,1);
}

/* ----------------------------------------------------------------------
   called by one proc in universe
   prints a warning message to the screen
------------------------------------------------------------------------- */

void Error::universe_warn(const char *file, int line, const char *str)
{
  if (universe->uscreen)
    fprintf(universe->uscreen,"WARNING on proc %d: %s (%s:%d)\n",
            universe->me,str,file,line);
}

/* ----------------------------------------------------------------------
   called by all procs in one world
   close all output, screen, and log files in world
   insure all procs in world call, else will hang
   force MPI_Abort if running in multi-partition mode
------------------------------------------------------------------------- */

void Error::all(const char *file, int line, const char *str)
{
  MPI_Barrier(world);

  int me;
  MPI_Comm_rank(world,&me);

  if (me == 0) {
    if (screen) fprintf(screen,"ERROR: %s (%s:%d)\n",str,file,line);
    if (logfile) fprintf(logfile,"ERROR: %s (%s:%d)\n",str,file,line);
  }

  if (output) delete output;
  if (screen && screen != stdout) fclose(screen);
  if (logfile) fclose(logfile);

  if (universe->nworlds > 1) MPI_Abort(universe->uworld,1);
  MPI_Finalize();
  exit(1);
}

/* ----------------------------------------------------------------------
   similar to Error::all
   called by by fix constructors so fixes identify themselves correctly
   even if derived class
------------------------------------------------------------------------- */

void Error::fix_error(const char *file, int line, Fix *fix,const char *str)
{
  MPI_Barrier(world);

  int me;
  MPI_Comm_rank(world,&me);

  if (me == 0) {
    if(strlen(str) > 2)
    {
        if (screen) fprintf(screen,"ERROR: Fix %s (id %s): %s (%s:%d)\n",fix->style,fix->id,str,file,line);
        if (logfile) fprintf(logfile,"ERROR: Fix %s (id %s): %s (%s:%d)\n",fix->style,fix->id,str,file,line);
    }
    else
    {
        if (screen) fprintf(screen,"ERROR: Illegal fix %s (id %s) command (%s:%d)\n",fix->style,fix->id,file,line);
        if (logfile) fprintf(logfile,"ERROR: Illegal fix %s (id %s) command (%s:%d)\n",fix->style,fix->id,file,line);
    }
  }

  if (output) delete output;
  if (screen && screen != stdout) fclose(screen);
  if (logfile) fclose(logfile);

  if (universe->nworlds > 1) MPI_Abort(universe->uworld,1);
  MPI_Finalize();
  exit(1);
}

void Error::compute_error(const char *file, int line, Compute *compute,const char *str)
{
  MPI_Barrier(world);

  int me;
  MPI_Comm_rank(world,&me);

  if (me == 0) {
    if(strlen(str) > 2)
    {
        if (screen) fprintf(screen,"ERROR: Compute %s (id %s): %s (%s:%d)\n",compute->style,compute->id,str,file,line);
        if (logfile) fprintf(logfile,"ERROR: Compute %s (id %s): %s (%s:%d)\n",compute->style,compute->id,str,file,line);
    }
    else
    {
        if (screen) fprintf(screen,"ERROR: Illegal compute %s (id %s) command (%s:%d)\n",compute->style,compute->id,file,line);
        if (logfile) fprintf(logfile,"ERROR: Illegal compute %s (id %s) command (%s:%d)\n",compute->style,compute->id,file,line);
    }
  }

  if (output) delete output;
  if (screen && screen != stdout) fclose(screen);
  if (logfile) fclose(logfile);

  if (universe->nworlds > 1) MPI_Abort(universe->uworld,1);
  MPI_Finalize();
  exit(1);
}

/* ----------------------------------------------------------------------
   cg error and warnign message
------------------------------------------------------------------------- */

void Error::cg(const char *file, int line, const char *str)
{
    char *catstr = new char[strlen(str)+1+100];
    strcpy(catstr,"The following model does not yield consistent results with coarse-graining: ");
    strcat(catstr,str);
    if(force->error_cg())
      all(file,line,catstr);
    else
      warning(file,line,catstr,1);
    delete []catstr;
}

/* ----------------------------------------------------------------------
   called by one proc in world
   write to world screen only if non-NULL on this proc
   always write to universe screen
   forces abort of entire world (and universe) if any proc in world calls
------------------------------------------------------------------------- */

void Error::one(const char *file, int line, const char *str)
{
  int me;
  MPI_Comm_rank(world,&me);
  if (screen) fprintf(screen,"ERROR on proc %d: %s (%s:%d)\n",
                      me,str,file,line);
  if (universe->nworlds > 1)
    if (universe->uscreen)
      fprintf(universe->uscreen,"ERROR on proc %d: %s (%s:%d)\n",
              universe->me,str,file,line);
  MPI_Abort(world,1);
}

/* ----------------------------------------------------------------------
   called by one proc in world
   only write to screen if non-NULL on this proc since could be file
------------------------------------------------------------------------- */

void Error::warning(const char *file, int line, const char *str, int logflag)
{
  if (screen) fprintf(screen,"WARNING: %s (%s:%d)\n",str,file,line);
  if (logflag && logfile) fprintf(logfile,"WARNING: %s (%s:%d)\n",
                                  str,file,line);
}

/* ----------------------------------------------------------------------
   called by one proc in world, typically proc 0
   write message to screen and logfile (if logflag is set)
------------------------------------------------------------------------- */

void Error::message(const char *file, int line, const char *str, int logflag)
{
  if (screen) fprintf(screen,"%s (%s:%d)\n",str,file,line);
  if (logflag && logfile) fprintf(logfile,"%s (%s:%d)\n",str,file,line);
}

/* ----------------------------------------------------------------------
   shutdown LAMMPS
   called by all procs in one world
   close all output, screen, and log files in world
   no abort, so insure all procs in world call, else will hang
------------------------------------------------------------------------- */

void Error::done()
{
  MPI_Barrier(world);

  if (output) delete output;
  if (screen && screen != stdout) fclose(screen);
  if (logfile) fclose(logfile);

  MPI_Finalize();
  exit(1);
}
