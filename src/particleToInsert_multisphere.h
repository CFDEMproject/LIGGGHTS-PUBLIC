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
    (if not contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef LMP_PARTICLE_TO_INSERT_MULTISPHERE_H
#define LMP_PARTICLE_TO_INSERT_MULTISPHERE_H

#include "particleToInsert.h"

using namespace LAMMPS_NS;

namespace LAMMPS_NS {
    class ParticleToInsertMultisphere : public ParticleToInsert
    {
        public:

           ParticleToInsertMultisphere(LAMMPS* lmp,int ns);
           virtual ~ParticleToInsertMultisphere();

           // per-particle displace in body coordinates
           double **displace;

           // vector to center of bounding sphere in body coos
           double xcm_to_xbound[3];

           // center of mass, should be 0/0/0
           double xcm_ins[3];

           double quat_ins[4];

           // body principal axes in space coords
           // = coordinate axis for body co sys
           double ex_space[3],ey_space[3],ez_space[3];

           double inertia[3];

           bool fflag[3],tflag[3];

           int type_ms;

           int insert();
           int check_near_set_x_v_omega(double *x,double *v, double *omega, double *quat, RegionNeighborList & neighList);//, double **xnear, int &nnear)
           int set_x_v_omega(double *,double *,double *, double *);

           void random_rotate(double,double,double);

           virtual void scale_pti(double r_scale);
    };
}

#endif
