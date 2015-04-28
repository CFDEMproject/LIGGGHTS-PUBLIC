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

#ifndef LMP_UNIVERSE_H
#define LMP_UNIVERSE_H

#include "stdio.h"
#include "pointers.h"

namespace LAMMPS_NS {

class Universe : protected Pointers {
 public:
  char *version;          // LAMMPS version string = date

  MPI_Comm uworld;        // communicator for entire universe
  int me,nprocs;          // my place in universe

  FILE *uscreen;          // universe screen output
  FILE *ulogfile;         // universe logfile
                          
  FILE *uthermofile;      // universe thermofile

  int existflag;          // 1 if universe exists due to -partition flag
  int nworlds;            // # of worlds in universe
  int iworld;             // which world I am in
  int *procs_per_world;   // # of procs in each world
  int *root_proc;         // root proc in each world

  MPI_Comm uorig;         // original communicator passed to LAMMPS instance
  int *uni2orig;          // proc I in universe uworld is
                          // proc uni2orig[I] in original communicator

  char *universe_id;       

  Universe(class LAMMPS *, MPI_Comm);
  ~Universe();
  void reorder(char *, char *);
  void add_world(char *);
  int consistent();
  void id(char *id);       
};

}

#endif

/* ERROR/WARNING messages:

E: Invalid -reorder N value

Self-explanatory.

E: Nprocs not a multiple of N for -reorder

Self-explanatory.

E: Cannot open -reorder file

Self-explanatory.

E: Unexpected end of -reorder file

Self-explanatory.

E: Invalid entry in -reorder file

Self-explanatory.

E: Invalid command-line argument

One or more command-line arguments is invalid.  Check the syntax of
the command you are using to launch LAMMPS.

*/
