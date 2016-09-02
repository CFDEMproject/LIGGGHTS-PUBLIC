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
    Arno Mayrhofer (CFDEMresearch GmbH, Linz)

    Copyright 2016-     CFDEMresearch GmbH, Linz
------------------------------------------------------------------------- */

/*
 * Fix for computing the stress tensor according to
 * Goldhirsch 2010 and/or the strain tensor according
 * to Zhang et al. 2010
 * Boundary influences to the stress tensor were included
 * according to Weinhart et al. 2012.
 */
#ifdef FIX_CLASS

FixStyle(continuum/weighted,FixContinuumWeighted)

#else

#ifndef LMP_FIX_CONTINUUM_WEIGHTED_H
#define LMP_FIX_CONTINUUM_WEIGHTED_H

#include "fix.h"
#include <vector>

namespace LAMMPS_NS {

enum kernel_type_t {
    TOP_HAT,
    GAUSSIAN,
    WENDLAND
};

class FixContinuumWeighted : public Fix {
  public:
    FixContinuumWeighted(class LAMMPS *, int, char **);
    ~FixContinuumWeighted();
    void post_create();

    int setmask();
    void init();

    void post_integrate();

    double get_phi(const double r);
    double get_grad_phi(const double r);

  private:
    double kernel_radius_;
    double kernel_sqRadius_;
    class PairGran *pairgran_;
    class FixPropertyAtom *fix_stress_;
    class FixPropertyAtom *fix_strain_;
    class FixPropertyAtom *fix_cont_vars_;
    class FixContactPropertyAtom *fix_contact_forces_;
    std::vector<FixContactPropertyAtom *> fix_wall_contact_forces_vector_;
    double integrate_phi(const double *const xij, const double *const nkj, const double a, const double b);
    template<kernel_type_t kernel_type> double weightingFunction(const double r);
    template<kernel_type_t kernel_type> double gradWeightingFunction(const double r);
    inline double compute_line_sphere_intersection(const double *const xij, const double *const xjk);
    bool compute_stress;
    bool compute_strain;
    kernel_type_t kernel_type;
};

}

#endif
#endif
