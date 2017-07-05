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
    This file is from LAMMPS, but has been modified. Copyright for
    modification:

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz

    Copyright of original file:
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

#include <mpi.h>
#include "compute_erotate_sphere.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "modify.h" 
#include "fix_multisphere.h" 
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

#define INERTIA 0.4          // moment of inertia prefactor for sphere

/* ---------------------------------------------------------------------- */

ComputeERotateSphere::ComputeERotateSphere(LAMMPS *lmp, int &iarg, int narg, char **arg) :
    Compute(lmp, iarg, narg, arg),
    pfactor(1.0),
    halfstep(false),
    fix_ms(NULL)
{
    if (narg < iarg || narg > iarg+1) error->all(FLERR,"Illegal compute erotate/sphere command");

    if (narg == iarg+1)
    {
        if (strcmp(arg[iarg++], "halfstep") == 0)
            halfstep = true;
        else
            error->all(FLERR,"Illegal compute erotate/sphere option");
    }

    scalar_flag = 1;
    extscalar = 1;

    // error check

    if (!atom->sphere_flag)
      error->all(FLERR,"Compute erotate/sphere requires atom style sphere");
}

/* ---------------------------------------------------------------------- */

void ComputeERotateSphere::init()
{
  pfactor = 0.5 * force->mvv2e * INERTIA;

  fix_ms =  static_cast<FixMultisphere*>(modify->find_fix_style("multisphere",0)); 
}

/* ---------------------------------------------------------------------- */

double ComputeERotateSphere::compute_scalar()
{
    if (invoked_scalar == update->ntimestep)
        return scalar;
    invoked_scalar = update->ntimestep;

    double **omega = atom->omega;
    double **torque = atom->torque;
    double *radius = atom->radius;
    double *rmass = atom->rmass;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    // sum rotational energy for each particle
    // point particles will not contribute, due to radius = 0.0

    double erotate = 0.0;
    for (int i = 0; i < nlocal; i++)
    {
        if (mask[i] & groupbit && (!fix_ms || fix_ms->belongs_to(i) < 0)) 
        {
            const double mult = halfstep ? update->dt*0.5/(INERTIA*radius[i]*radius[i]*rmass[i]) : 0.0;
            const double omega0 = omega[i][0] + mult*torque[i][0];
            const double omega1 = omega[i][1] + mult*torque[i][1];
            const double omega2 = omega[i][2] + mult*torque[i][2];
            erotate += (omega0*omega0 + omega1*omega1 + omega2*omega2) * radius[i]*radius[i]*rmass[i];
        }
    }

    MPI_Allreduce(&erotate,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
    scalar *= pfactor;
    
    return scalar;
}
