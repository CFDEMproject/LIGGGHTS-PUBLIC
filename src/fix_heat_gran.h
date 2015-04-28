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

#ifndef LMP_FIX_HEATGRAN_ABSTRACT_H
#define LMP_FIX_HEATGRAN_ABSTRACT_H

#include "fix.h"

#define SMALL 1e-8

namespace LAMMPS_NS {

  class FixHeatGran : public Fix {
  public:
    FixHeatGran(class LAMMPS *, int, char **);
    ~FixHeatGran(){};
    virtual void post_create();
    virtual void pre_delete(bool unfixflag){ UNUSED(unfixflag); };

    void initial_integrate(int vflag);

    virtual double compute_scalar();
    virtual int setmask();
    virtual void init();

    // per default these three methods throw errors.
    virtual void cpl_evaluate(class ComputePairGranLocal *);
    virtual void register_compute_pair_local(class ComputePairGranLocal *);
    virtual void unregister_compute_pair_local(class ComputePairGranLocal *);
    void updatePtrs();

  protected:
    class ComputePairGranLocal *cpl;
    class FixPropertyAtom* fix_heatFlux;
    class FixPropertyAtom* fix_heatSource;
    class FixPropertyAtom* fix_temp;
    class FixScalarTransportEquation *fix_ste;
    class FixPropertyAtom* fix_directionalHeatFlux;

    double *heatFlux;   
    double *heatSource; 
    double *Temp;       
    double T0;          
    double **directionalHeatFlux;

    class PairGran *pair_gran;
    int history_flag;
  };

}

#endif
