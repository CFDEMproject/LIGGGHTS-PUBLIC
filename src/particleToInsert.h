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
    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Richard Berger (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2015 JKU Linz
------------------------------------------------------------------------- */

#ifndef LMP_PARTICLE_TO_INSERT_H
#define LMP_PARTICLE_TO_INSERT_H

#include "memory.h"
#include "pointers.h"
#include "region_neighbor_list.h"

using namespace LAMMPS_NS;

namespace LAMMPS_NS {
    class ParticleToInsert : protected Pointers
    {
     public:

        ParticleToInsert(LAMMPS* lmp, int ns = 1, class FixPropertyAtom *_fix_release = NULL);

        virtual ~ParticleToInsert();

        // insertion properties
        int nparticles; 
        int groupbit;
        int atom_type;
        double density_ins;
        double volume_ins;
        double area_ins;
        double mass_ins;
        double r_bound_ins;

        int distorder;

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

        // a custom id for the history writing
        // in fix insert/stream/predefined
        int id_ins;
        class FixPropertyAtom * const fix_release;

        // value of a fix property/atoms at insertion
        class FixPropertyAtom **fix_property;
        int n_fix_property;
        int *fix_property_nentry;
        double **fix_property_value;

        virtual int insert();
        virtual int check_near_set_x_v_omega(double *x,double *v, double *omega, double *quat, RegionNeighborList<interpolate_no> & neighList);
        virtual int check_near_set_x_v_omega_ms(double *x,double *v, double *omega, double *quat, RegionNeighborList<interpolate_no> & neighList);
        //virtual int check_near_set_x_v_omega(double *x,double *v, double *omega, double *quat, double **xnear, int &nnear);
        //virtual int check_near_set_x_v_omega_ms(double *x,double *v, double *omega, double *quat, double **xnear, int &nnear);

        virtual int set_x_v_omega(double *,double *,double *,double *);

        virtual void scale_pti(double r_scale);

        void setFixTemplate(FixPropertyAtom* fix_template)
        { fix_template_ = fix_template; }

    private:

        // Fix property atom that identifies the template
        FixPropertyAtom *fix_template_;
    };

}

#endif
