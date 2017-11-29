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

// Park/Miller RNG

#include <cmath>
#include "random_park.h"
#include "error.h"

using namespace LAMMPS_NS;

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836

/* ---------------------------------------------------------------------- */

RanPark::RanPark(LAMMPS *lmp, const char * seed_char, bool proc_shift, int multiplier) :
    Random(lmp, seed_char, proc_shift, multiplier)
{
  if (seed <= 0)
    error->one(FLERR,"Invalid seed for Park random # generator");
  save = 0;
}

/* ----------------------------------------------------------------------
   uniform RN
------------------------------------------------------------------------- */

double RanPark::uniform()
{
  const int k = seed/IQ;
  seed = IA*(seed-k*IQ) - IR*k;
  if (seed < 0) seed += IM;
  return AM*seed;
}

/* ----------------------------------------------------------------------
   gaussian RN
------------------------------------------------------------------------- */

double RanPark::gaussian()
{
  double first,v1,v2,rsq,fac;

  if (!save) {
    int again = 1;
    while (again) {
      v1 = 2.0*uniform()-1.0;
      v2 = 2.0*uniform()-1.0;
      rsq = v1*v1 + v2*v2;
      if (rsq < 1.0 && rsq != 0.0) again = 0;
    }
    fac = sqrt(-2.0*log(rsq)/rsq);
    second = v1*fac;
    first = v2*fac;
    save = 1;
  } else {
    first = second;
    save = 0;
  }
  return first;
}

/* ---------------------------------------------------------------------- */

void RanPark::reset(int seed_init)
{
  if (seed_init <= 0)
    error->all(FLERR,"Invalid seed for Park random # generator");
  seed = seed_init;
  save = 0;
}

/* ----------------------------------------------------------------------
   reset the seed based on atom coords and ibase = caller seed
   use hash function, treating user seed and coords as sequence of input ints
   this is Jenkins One-at-a-time hash, see Wikipedia entry on hash tables
------------------------------------------------------------------------- */

void RanPark::reset(int ibase, double *coord)
{
  char *str = (char *) &ibase;
  int n = sizeof(int);

  unsigned int hash = 0;
  for (int i = 0; i < n; i++) {
    hash += str[i];
    hash += (hash << 10);
    hash ^= (hash >> 6);
  }

  str = (char *) coord;
  n = 3 * sizeof(double);
  for (int i = 0; i < n; i++) {
    hash += str[i];
    hash += (hash << 10);
    hash ^= (hash >> 6);
  }

  hash += (hash << 3);
  hash ^= (hash >> 11);
  hash += (hash << 15);

  // keep 31 bits of unsigned int as new seed
  // do not allow seed = 0, since will cause hang in gaussian()

  seed = hash & 0x7ffffff;
  if (!seed) seed = 1;

  // warm up the RNG

  for (int i = 0; i < 5; i++) uniform();
  save = 0;
}

/* ---------------------------------------------------------------------- */

int RanPark::state()
{
  return seed;
}
