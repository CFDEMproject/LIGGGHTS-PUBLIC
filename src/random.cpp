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
    (if no contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2016-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#include <cmath>
#include "random.h"
#include "error.h"
#include "comm.h"
#include "input.h"
#include "math_extra_liggghts.h"
#include <stdlib.h>
#include <string>
#include <climits>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Random::Random(LAMMPS *lmp, const char * seed_char, bool proc_shift, int multiplier) : Pointers(lmp)
{
    if (!seed_char)
        error->all(FLERR, "Internal error: NULL seed_char");
    long seedl = atol(seed_char);
    seed = atoi(seed_char);
    if ((long)seed != seedl)
    {
        char errstr[1024];

        sprintf(errstr,"Seed %ld is larger than INT_MAX (%d)\n", seedl, INT_MAX);
        error->all(FLERR, errstr);
    }

    const int offset = proc_shift ? multiplier*comm->me : 0;
    if (seedl + (long)offset > (long)INT_MAX)
    {
        char errstr[1024];

        sprintf(errstr,"Seed %ld + %d (offset) is larger than INT_MAX (%d)\n", seedl, offset, INT_MAX);
        error->all(FLERR, errstr);
    }
    seed += offset;

    if(0 == comm->me)
    {
        if(seed < 10000 || !MathExtraLiggghts::isPrime(seed))
        {
            char errstr[1024];

            sprintf(errstr,"Random number generation: It is required that the random seed value is > 10000 and a prime number.\n"
                           "The random seed used was %d\n"
                           "  Hint 1: start with 'liggghts -echo both < in.script' to find out which command caused this\n"
                           "  Hint 2: possible valid seeds would be the following numbers:\n"
                           "          15485863, 15485867, 32452843, 32452867, 49979687, 49979693, 67867967, 67867979, 86028121, 86028157",
                           seed);

            if(input->seed_check_throw_error())
                error->one(FLERR,errstr);
            else
                error->warning(FLERR,errstr);
        }
    }

}
