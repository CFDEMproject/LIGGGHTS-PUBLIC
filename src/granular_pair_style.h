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

/* ----------------------------------------------------------------------
   Contributing authors:
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */

#ifndef GRANULAR_PAIR_STYLE_H
#define GRANULAR_PAIR_STYLE_H

#include "utils.h"
#include "pair_gran.h"

namespace LIGGGHTS {

namespace PairStyles {

  /**
   * @brief interface of LIGGGHTS granular pair styles
   */
  class IGranularPairStyle {
  public:
    typedef LAMMPS_NS::PairGran ParentType;

    virtual ~IGranularPairStyle();
    virtual void settings(int nargs, char ** args) = 0;
    virtual void init_granular() = 0;
    virtual void write_restart_settings(FILE * fp) = 0;
    virtual void read_restart_settings(FILE * fp) = 0;
    virtual void compute_force(LAMMPS_NS::PairGran * pg, int eflag, int vflag, int addflag) = 0;

    virtual double stressStrainExponent() = 0;
    virtual int64_t hashcode() = 0;
  };

  /**
   * @brief LIGGGHTS granular pair style factory
   */
  class Factory : public Utils::AbstractFactory<IGranularPairStyle> {
    Factory() {}
  public:
    static Factory & instance();
  };
}

}

#endif // GRANULAR_PAIR_STYLE_H
