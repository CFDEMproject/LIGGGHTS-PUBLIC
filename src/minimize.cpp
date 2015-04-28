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

#include "stdlib.h"
#include "minimize.h"
#include "domain.h"
#include "update.h"
#include "min.h"
#include "finish.h"
#include "timer.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Minimize::Minimize(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void Minimize::command(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR,"Illegal minimize command");

  if (domain->box_exist == 0)
    error->all(FLERR,"Minimize command before simulation box is defined");

  update->etol = force->numeric(FLERR,arg[0]);
  update->ftol = force->numeric(FLERR,arg[1]);
  update->nsteps = force->inumeric(FLERR,arg[2]);
  update->max_eval = force->inumeric(FLERR,arg[3]);

  if (update->etol < 0.0 || update->ftol < 0.0)
    error->all(FLERR,"Illegal minimize command");

  update->whichflag = 2;
  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->firststep + update->nsteps;
  if (update->laststep < 0 || update->laststep > MAXBIGINT)
    error->all(FLERR,"Too many iterations");

  lmp->init();
  update->minimize->setup();

  timer->init();
  timer->barrier_start(TIME_LOOP);
  update->minimize->run(update->nsteps);
  timer->barrier_stop(TIME_LOOP);

  update->minimize->cleanup();

  Finish finish(lmp);
  finish.end(1);

  update->whichflag = 0;
  update->firststep = update->laststep = 0;
  update->beginstep = update->endstep = 0;
}
