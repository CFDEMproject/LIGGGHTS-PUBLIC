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

    Arno Mayrhofer (DCS Computing GmbH, Linz)

    This file is from LAMMPS, but has been modified. Copyright for
    modification:

    Copyright 2017-     DCS Computing GmbH, Linz

    Copyright of original file:
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include <mpi.h>
#include "lammps.h"
#include "input.h"
#include <string.h>
#include <signal.h>
#include "signal_handling.h"

#ifdef LIGGGHTS_DEBUG
#include "fenv.h"
#endif

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   main program to drive LAMMPS
------------------------------------------------------------------------- */

int main(int argc, char **argv)
{
  #if !defined(_WINDOWS) && !defined(__MINGW32__)
  struct sigaction int_action, usr1_action, term_action;
  memset(&int_action, 0, sizeof(struct sigaction));
  memset(&usr1_action, 0, sizeof(struct sigaction));
  memset(&term_action, 0, sizeof(struct sigaction));

  int_action.sa_handler = SignalHandler::int_handler;
  sigaction(SIGINT, &int_action, NULL);
  // SIGTERM is handled the same way as sigint. Note that OpenMPI (and possibly other flavours)
  // convert a SIGINT to mpirun to a SIGTERM to its children. That's why we need to catch it too.
  sigaction(SIGTERM, &int_action, NULL);
  usr1_action.sa_handler = SignalHandler::usr1_handler;
  sigaction(SIGUSR1, &usr1_action, NULL);
  #else  // _WINDOWS or __MINGW32__
  signal(SIGINT, SignalHandler::int_handler);
  signal(SIGTERM, SignalHandler::int_handler);
  // no SIGUSR1 treatment because Windows
  #endif // _WINDOWS

  MPI_Init(&argc,&argv);
  #ifdef LIGGGHTS_DEBUG
  feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO);
  #endif

  LAMMPS *lammps = new LAMMPS(argc,argv,MPI_COMM_WORLD);
  lammps->input->file();
  delete lammps;

  MPI_Finalize();
}
