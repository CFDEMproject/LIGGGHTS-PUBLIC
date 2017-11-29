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
    Arno Mayrhofer (DCS Computing GmbH, Linz)
    Stefan Radl (TU Graz)

    Copyright 2017 - DCS Computing GmbH, Linz
    Copyright 2016 - TU Graz
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(multisphere/break,FixMultisphereBreak)

#else

#ifndef LMP_FIX_MULTISPHERE_BREAK_H
#define LMP_FIX_MULTISPHERE_BREAK_H

#include "fix_multisphere.h"

namespace LAMMPS_NS
{

class FixMultisphereBreak : public FixMultisphere
{
public:

    FixMultisphereBreak(class LAMMPS *, int, char **);
    virtual ~FixMultisphereBreak();

    void init();
    void final_integrate();
    void pre_neighbor();
    void calc_force(bool setupflag);

protected:

    char*               triggerFixName_;
    FixPropertyAtom*    triggerFix_;
    int                 triggerIdx_;
    double *            triggerArray_;
    int                 maxatom_;
    int                 triggerType_;
    char*               triggerName_;
    int                 triggerIndex_;
    double              triggerThreshold_;
    int                 triggerTimeStep_;
};

}

#endif
#endif
