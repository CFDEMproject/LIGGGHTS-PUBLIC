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

/* ----------------------------------------------------------------------
   Contributing author:  Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "string.h"
#include "write_dump.h"
#include "style_dump.h"
#include "dump.h"
#include "dump_image.h"
#include "atom.h"
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

void WriteDump::command(int narg, char **arg)
{

  if (narg < 3) error->all(FLERR,"Illegal write_dump command");

  // modindex = index in args of "modify" keyword
  // will be narg if "modify" is not present

  int modindex;
  for (modindex = 0; modindex < narg; modindex++)
    if (strcmp(arg[modindex],"modify") == 0) break;

  // create the Dump instance
  // create dump command line with extra required args

  Dump *dump = NULL;

  char **dumpargs = new char*[modindex+2];
  dumpargs[0] = (char *) "WRITE_DUMP"; // dump id
  dumpargs[1] = arg[0];                // group
  dumpargs[2] = arg[1];                // dump style
  dumpargs[3] = (char *) "0";          // dump frequency

  for (int i = 2; i < modindex; ++i)
    dumpargs[i+2] = arg[i];

  if (0) return;         // dummy line to enable else-if macro expansion

#define DUMP_CLASS
#define DumpStyle(key,Class) \
  else if (strcmp(arg[1],#key) == 0) dump = new Class(lmp,modindex+2,dumpargs);
#include "style_dump.h"
#undef DUMP_CLASS

  else error->all(FLERR,"Invalid dump style");

  if (modindex < narg) dump->modify_params(narg-modindex-1,&arg[modindex+1]);

  // write out one frame and then delete the dump again
  // set multifile_override for DumpImage so that filename needs no "*"

  if (strcmp(arg[1],"image") == 0)
    ((DumpImage *) dump)->multifile_override = 1;

  dump->init();
  dump->write();

  // delete the Dump instance and local storage

  delete dump;
  delete[] dumpargs;
}

