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
