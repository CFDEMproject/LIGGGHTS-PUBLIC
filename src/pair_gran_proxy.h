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

#ifdef PAIR_CLASS
PairStyle(gran,PairGranProxy)
#else

#ifndef PAIR_GRAN_PROXY_H
#define PAIR_GRAN_PROXY_H

#include "pair_gran.h"
#include "granular_pair_style.h"

namespace LAMMPS_NS
{
class PairGranProxy : public PairGran
{
  LIGGGHTS::PairStyles::IGranularPairStyle * impl;

public:
  PairGranProxy(LAMMPS * lmp);
  virtual ~PairGranProxy();

  virtual void settings(int nargs, char ** args);
  virtual void init_granular();
  virtual void write_restart_settings(FILE * fp);
  virtual void read_restart_settings(FILE * fp);
  virtual void compute_force(int eflag, int vflag, int addflag);

  virtual double stressStrainExponent();
  virtual int64_t hashcode();
};
}

#endif // PAIR_GRAN_PROXY_H

#endif
