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

#include <mpi.h>
#include "compute_ke.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "group.h"
#include "error.h"
#include "modify.h" 
#include "fix_multisphere.h" 

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeKE::ComputeKE(LAMMPS *lmp, int &iarg, int narg, char **arg) :
    Compute(lmp, iarg, narg, arg),
    pfactor(1.0),
    halfstep(false),
    fix_ms(NULL)
{
    if (narg < iarg || narg > iarg+1) error->all(FLERR,"Illegal compute ke command");

    if (narg == iarg+1)
    {
        if (strcmp(arg[iarg++], "halfstep") == 0)
            halfstep = true;
        else
            error->all(FLERR,"Illegal compute ke option");
    }

    scalar_flag = 1;
    extscalar = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeKE::init()
{
  pfactor = 0.5 * force->mvv2e;
  fix_ms =  static_cast<FixMultisphere*>(modify->find_fix_style("multisphere",0)); 
}

/* ---------------------------------------------------------------------- */

double ComputeKE::compute_scalar()
{
    invoked_scalar = update->ntimestep;

    double **v = atom->v;
    double **f = atom->f;
    double *rmass = atom->rmass;
    double *mass = atom->mass;
    int *mask = atom->mask;
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double ke = 0.0;

    if (rmass) {
        for (int i = 0; i < nlocal; i++)
        {
            if (mask[i] & groupbit && (!fix_ms || fix_ms->belongs_to(i) < 0))
            {
                const double mult = halfstep ? update->dt*0.5/rmass[i] : 0.0;
                const double v0 = v[i][0] + mult*f[i][0];
                const double v1 = v[i][1] + mult*f[i][1];
                const double v2 = v[i][2] + mult*f[i][2];
                ke += rmass[i] * (v0*v0 + v1*v1 + v2*v2);
            }
        }
    }
    else
    {
        for (int i = 0; i < nlocal; i++)
            if (mask[i] & groupbit)
                ke += mass[type[i]] *
                      (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
    }

    MPI_Allreduce(&ke,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
    scalar *= pfactor;

    // for multispheres we need to get the kinetic energy from the multisphere fix
    if (fix_ms)
        scalar += fix_ms->extract_ke();

    return scalar;
}
