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
#include "stdlib.h"
#include "string.h"
#include "rerun.h"
#include "read_dump.h"
#include "domain.h"
#include "update.h"
#include "integrate.h"
#include "modify.h"
#include "output.h"
#include "finish.h"
#include "timer.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Rerun::Rerun(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void Rerun::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Rerun command before simulation box is defined");

  if (narg < 2) error->all(FLERR,"Illegal rerun command");

  // list of dump files = args until a keyword

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"first") == 0) break;
    if (strcmp(arg[iarg],"last") == 0) break;
    if (strcmp(arg[iarg],"every") == 0) break;
    if (strcmp(arg[iarg],"skip") == 0) break;
    if (strcmp(arg[iarg],"start") == 0) break;
    if (strcmp(arg[iarg],"stop") == 0) break;
    if (strcmp(arg[iarg],"dump") == 0) break;
    iarg++;
  }
  int nfile = iarg;
  if (nfile == 0 || nfile == narg) error->all(FLERR,"Illegal rerun command");

  // parse optional args up until "dump"
  // user MAXBIGINT -1 so Output can add 1 to it and still be a big int

  bigint first = 0;
  bigint last = MAXBIGINT - 1;
  int nevery = 0;
  int nskip = 1;
  int startflag = 0;
  int stopflag = 0;
  bigint start=0,stop=0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"first") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal rerun command");
      first = ATOBIGINT(arg[iarg+1]);
      if (first < 0) error->all(FLERR,"Illegal rerun command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"last") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal rerun command");
      last = ATOBIGINT(arg[iarg+1]);
      if (last < 0) error->all(FLERR,"Illegal rerun command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal rerun command");
      nevery = force->inumeric(FLERR,arg[iarg+1]);
      if (nevery < 0) error->all(FLERR,"Illegal rerun command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"skip") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal rerun command");
      nskip = force->inumeric(FLERR,arg[iarg+1]);
      if (nskip <= 0) error->all(FLERR,"Illegal rerun command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"start") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal rerun command");
      startflag = 1;
      start = ATOBIGINT(arg[iarg+1]);
      if (start < 0) error->all(FLERR,"Illegal rerun command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"stop") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal rerun command");
      stopflag = 1;
      stop = ATOBIGINT(arg[iarg+1]);
      if (stop < 0) error->all(FLERR,"Illegal rerun command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"dump") == 0) {
      break;
    } else error->all(FLERR,"Illegal rerun command");
  }

  int nremain = narg - iarg - 1;
  if (nremain <= 0) error->all(FLERR,"Illegal rerun command");
  if (first > last) error->all(FLERR,"Illegal rerun command");
  if (startflag && stopflag && start > stop)
    error->all(FLERR,"Illegal rerun command");

  // pass list of filenames to ReadDump
  // along with post-"dump" args and post-"format" args

  ReadDump *rd = new ReadDump(lmp);

  rd->store_files(nfile,arg);
  if (nremain)
    nremain = rd->fields_and_keywords(nremain,&arg[narg-nremain]);
  else nremain = rd->fields_and_keywords(0,NULL);
  if (nremain) rd->setup_reader(nremain,&arg[narg-nremain]);
  else rd->setup_reader(0,NULL);

  // perform the psuedo run
  // invoke lmp->init() only once
  // read all relevant snapshots
  // uset setup_minimal() since atoms are already owned by correct procs
  // addstep_compute_all() insures energy/virial computed on every snapshot

  update->whichflag = 1;

  if (startflag) update->beginstep = update->firststep = start;
  else update->beginstep = update->firststep = first;
  if (stopflag) update->endstep = update->laststep = stop;
  else update->endstep = update->laststep = last;

  int firstflag = 1;
  int ndump = 0;

  lmp->init();

  timer->init();
  timer->barrier_start(TIME_LOOP);

  bigint ntimestep = rd->seek(first,0);
  if (ntimestep < 0)
    error->all(FLERR,"Rerun dump file does not contain requested snapshot");

  while (1) {
    ndump++;
    rd->header(firstflag);
    update->reset_timestep(ntimestep);
    rd->atoms();

    modify->init();
    update->integrate->setup_minimal(1);
    modify->end_of_step();
    if (firstflag) output->setup();
    else if (output->next) output->write(ntimestep);

    firstflag = 0;
    ntimestep = rd->next(ntimestep,last,nevery,nskip);
    if (ntimestep < 0) break;
  }

  // insure thermo output on last dump timestep

  output->next_thermo = update->ntimestep;
  output->write(update->ntimestep);

  timer->barrier_stop(TIME_LOOP);

  update->integrate->cleanup();

  // set update->nsteps to ndump for Finish stats to print

  update->nsteps = ndump;

  Finish finish(lmp);
  finish.end(1);

  update->whichflag = 0;
  update->firststep = update->laststep = 0;
  update->beginstep = update->endstep = 0;

  // clean-up

  delete rd;
}
