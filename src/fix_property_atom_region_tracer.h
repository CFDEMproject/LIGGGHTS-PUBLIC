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

    Copyright 2013-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(property/atom/timetracer,FixPropertyAtomRegionTracer)
FixStyle(property/atom/regiontracer/time,FixPropertyAtomRegionTracer)

#else

#ifndef LMP_FIX_PROPERTY_ATOM_REGION_TRACER_H
#define LMP_FIX_PROPERTY_ATOM_REGION_TRACER_H

#include "fix_property_atom.h"
#include <vector>

namespace LAMMPS_NS {

class FixPropertyAtomRegionTracer : public FixPropertyAtom {

public:

    FixPropertyAtomRegionTracer(class LAMMPS *, int, char **, bool parse = true);
    ~FixPropertyAtomRegionTracer();

    virtual void init();
    virtual int setmask();
    void end_of_step();

    virtual double compute_vector(int n);

    inline int n_values()
    { return 0 == iregion_.size() ? 1: iregion_.size(); }

    // makes this a time tracer by default
    virtual inline double trace_multiplier(const int i) const
    { return 1.0; }

protected:

    int iarg_;

    // params for region marker
    std::vector<int>   iregion_;
    std::vector<char*> idregion_;

    int check_every_;
}; //end class

}
#endif
#endif
