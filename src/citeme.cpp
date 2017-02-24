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

#include "citeme.h"
#include "version.h"
#include "universe.h"
#include "error.h"

using namespace LAMMPS_NS;

static const char cite_header[] =
  "This LAMMPS simulation made specific use of work described in the\n"
  "following references.  See http://lammps.sandia.gov/cite.html\n"
  "for details.\n\n";

static const char cite_nagline[] = "\nPlease see the log.cite file "
  "for references relevant to this simulation\n\n";

/* ---------------------------------------------------------------------- */

CiteMe::CiteMe(LAMMPS *lmp) : Pointers(lmp)
{
  fp = NULL;
  cs = new citeset();
}

/* ----------------------------------------------------------------------
   write out nag-line at the end of the regular output and clean up
------------------------------------------------------------------------- */

CiteMe::~CiteMe()
{
  if (universe->me || cs->size() == 0) {
    delete cs;
    return;
  }

  delete cs;

  if (screen) fprintf(screen,cite_nagline);
  if (logfile) fprintf(logfile,cite_nagline);

  if (fp) fclose(fp);
}

/* ----------------------------------------------------------------------
   write out and register a citation so it will be written only once
------------------------------------------------------------------------- */

void CiteMe::add(const char *ref)
{
  if (universe->me) return;
  if (cs->find(ref) != cs->end()) return;
  cs->insert(ref);

  if (!fp) {
    fp = fopen("log.cite","w");
    if (!fp) error->universe_one(FLERR,"Could not open log.cite file");
    fputs(cite_header,fp);
    fflush(fp);
  }

  fputs(ref,fp);
  fflush(fp);
}
