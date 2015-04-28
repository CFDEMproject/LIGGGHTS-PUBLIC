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

#ifndef LMP_PARTICLE_TO_INSERT_H
#define LMP_PARTICLE_TO_INSERT_H

#include "memory.h"
#include "pointers.h"

using namespace LAMMPS_NS;

namespace LAMMPS_NS {
    class ParticleToInsert : protected Pointers
    {
     public:

        ParticleToInsert(LAMMPS* lmp,int ns = 1);

        virtual ~ParticleToInsert();

        // insertion properties
        int nspheres;
        int groupbit;
        int atom_type;
        double density_ins;
        double volume_ins;
        double mass_ins;
        double r_bound_ins;

        // per-sphere radius, position
        // if atom_type_vector exists, each sphere has different type
        double *radius_ins;
        double **x_ins;
        bool atom_type_vector_flag;
        int *atom_type_vector;

        // center of bounding sphere
        
        double x_bound_ins[3];

        // velocity and omega at insertion
        
        double v_ins[3];
        double omega_ins[3];

        // value of a fix property/atom at insertion
        class FixPropertyAtom *fix_property;
        double fix_property_value;

        virtual int insert();
        virtual int check_near_set_x_v_omega(double *x,double *v, double *omega, double *quat, double **xnear, int &nnear);
        virtual int check_near_set_x_v_omega_ms(double *x,double *v, double *omega, double *quat, double **xnear, int &nnear);
        virtual int set_x_v_omega(double *,double *,double *,double *);

        virtual void scale_pti(double r_scale);
    };

}

#endif
