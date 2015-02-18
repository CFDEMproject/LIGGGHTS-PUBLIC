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

#ifdef COMPUTE_CLASS

ComputeStyle(nparticles/tracer/region,ComputeNparticlesTracerRegion)

#else

#ifndef LMP_COMPUTE_NPARTICLES_TRACER_REGION_H
#define LMP_COMPUTE_NPARTICLES_TRACER_REGION_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeNparticlesTracerRegion : public Compute {

 public:

  ComputeNparticlesTracerRegion(class LAMMPS *, int, char **);
  ~ComputeNparticlesTracerRegion();

  void init();
  void compute_vector();

 private:

  template<bool IMAGE>
  void compute_vector_eval(bool, double&, double&);

  // image stuff
  int image_dim_, image_no_;
  bool reset_marker_;

  // params for regions where to mark and where to count
  int iregion_count_;
  char *idregion_count_;

  class FixPropertyAtomTracer *fix_tracer_;
  char *fix_tracer_name_;
};

}

#endif
#endif
