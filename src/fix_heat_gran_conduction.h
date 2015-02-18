/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
   the producer of the LIGGGHTS® software and the CFDEM®coupling software
   See http://www.cfdem.com/terms-trademark-policy for details.

   LIGGGHTS® is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(heat/gran/conduction,FixHeatGranCond)
FixStyle(heat/gran,FixHeatGranCond)

#else

#ifndef LMP_FIX_HEATGRAN_CONDUCTION_H
#define LMP_FIX_HEATGRAN_CONDUCTION_H

#include "fix_heat_gran.h"

namespace LAMMPS_NS {

  class FixHeatGranCond : public FixHeatGran {
  public:
    FixHeatGranCond(class LAMMPS *, int, char **);
    ~FixHeatGranCond();
    virtual void post_create();
    void pre_delete(bool);

    int setmask();
    void init();
    virtual void post_force(int);

    void cpl_evaluate(class ComputePairGranLocal *);
    void register_compute_pair_local(ComputePairGranLocal *);
    void unregister_compute_pair_local(ComputePairGranLocal *);

  protected:
    int iarg_;

  private:
    template <int,int> void post_force_eval(int,int);

    class FixPropertyGlobal* fix_conductivity_;
    double *conductivity_;

    // model for contact area calculation
    int area_calculation_mode_;

    double fixed_contact_area_;

    // for heat transfer area correction
    int area_correction_flag_;
    double const* const* deltan_ratio_;
  };

}

#endif
#endif

