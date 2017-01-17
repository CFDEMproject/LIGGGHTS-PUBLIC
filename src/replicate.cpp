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

#include <stdlib.h>
#include <string.h>
#include "replicate.h"
#include "atom.h"
#include "atom_vec.h"
#include "atom_vec_hybrid.h"
#include "force.h"
#include "domain.h"
#include "comm.h"
#include "special.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define LB_FACTOR 1.1
#define EPSILON   1.0e-6

/* ---------------------------------------------------------------------- */

Replicate::Replicate(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void Replicate::command(int narg, char **arg)
{
    if (domain->box_exist == 0)
        error->all(FLERR,"Replicate command before simulation box is defined");
    if (narg <= 7)
        error->all(FLERR,"Illegal replicate command");

    int me = comm->me;
    if (me == 0 && screen)
        fprintf(screen,"Replicating atoms ...\n");

    const int nx = force->inumeric(FLERR,arg[0]);
    const int ny = force->inumeric(FLERR,arg[1]);
    const int nz = force->inumeric(FLERR,arg[2]);
    if (strcmp(arg[3], "offset") != 0)
        error->all(FLERR, "Could not find offset keyword in replicate command");
    const double size_x = force->numeric(FLERR,arg[4]);
    const double size_y = force->numeric(FLERR,arg[5]);
    const double size_z = force->numeric(FLERR,arg[6]);
    double shift_x = 0;
    double shift_y = 0;
    double shift_z = 0;
    if (narg >= 8 && strcmp(arg[7], "shift") == 0)
    {
        if (narg != 11)
            error->all(FLERR, "Invalid replicate command");
        shift_x = force->numeric(FLERR,arg[8]);
        shift_y = force->numeric(FLERR,arg[9]);
        shift_z = force->numeric(FLERR,arg[10]);
    }

    // error and warning checks

    if (nx <= 0 || ny <= 0 || nz <= 0)
        error->all(FLERR,"Illegal replicate command");
    if (domain->dimension == 2 && nz != 1)
        error->all(FLERR,"Cannot replicate 2d simulation in z dimension");
    if ((nx > 1 && domain->xperiodic == 0) ||
        (ny > 1 && domain->yperiodic == 0) ||
        (nz > 1 && domain->zperiodic == 0))
    {
        if (comm->me == 0)
            error->warning(FLERR,"Replicating in a non-periodic dimension");
    }

    int triclinic = domain->triclinic;
    double old_xprd = domain->xprd;
    double old_yprd = domain->yprd;
    double old_zprd = domain->zprd;

    // setup new simulation box

    domain->print_box("  old: ");
    domain->boxlo[0] = domain->boxlo[0] + shift_x;
    domain->boxlo[1] = domain->boxlo[1] + shift_y;
    domain->boxlo[2] = domain->boxlo[2] + shift_z;
    domain->boxhi[0] = domain->boxlo[0] + nx*old_xprd;
    domain->boxhi[1] = domain->boxlo[1] + ny*old_yprd;
    domain->boxhi[2] = domain->boxlo[2] + nz*old_zprd;
    if (triclinic) {
        domain->xy *= ny;
        domain->xz *= nz;
        domain->yz *= nz;
    }

    // new problem setup using new box boundaries

    domain->print_box("  new: ");
    domain->set_initial_box();
    domain->set_global_box();
    comm->set_proc_grid();
    domain->set_local_box();
    domain->print_box("  fin: ");

    const int initial_natoms = atom->nlocal;
    for (int i=0; i<nx; i++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int k=0; k<nz; k++)
            {
                int offset = i+j+k==0 ? 0 : atom->nlocal;
                for (int l=0; l < initial_natoms; l++)
                {
                    if (i + j + k != 0)
                    {
                        double zero[3] = {0.0, 0.0, 0.0};
                        atom->avec->create_atom(atom->type[l], zero);
                        atom->avec->copy(l, l+offset, 0);
                        atom->x[l+offset][0] += i*size_x;
                        atom->x[l+offset][1] += j*size_y;
                        atom->x[l+offset][2] += k*size_z;
                        atom->tag[l+offset] += offset;
                    }
                    else
                    {
                        atom->x[l][0] += shift_x;
                        atom->x[l][1] += shift_y;
                        atom->x[l][2] += shift_z;
                    }
                }
            }
        }
    }
    if (atom->map_style)
    {
        atom->nghost = 0;
        atom->map_init();
        atom->map_set();
    }
    if (atom->molecular)
    {
        Special special(lmp);
        special.build();
    }
}
