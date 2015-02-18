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
Contributing author for SPH:
Andreas Aigner (CD Lab Particulate Flow Modelling, JKU)
andreas.aigner@jku.at
------------------------------------------------------------------------- */

/*-----------------------------------------------------------------------
2D Spiky Kernel by Markus Schörgenhumer, mkschoe@gmail.com
-------------------------------------------------------------------------*/

#ifdef SPH_KERNEL_CLASS

    // kernel identifier (a unique integer >= 0)
    // a name for the kernel
    // name of the functions for the kernel, its derivative, and the cutoff are defined
    SPHKernel
    (
        3,
        spiky2d,
        sph_kernel_spiky2d,
        sph_kernel_spiky2d_der,
        sph_kernel_spiky2d_cut
    )

#else

#ifndef LMP_SPH_KERNEL_SPIKY2D
#define LMP_SPH_KERNEL_SPIKY2D

namespace SPH_KERNEL_NS {
  inline double sph_kernel_spiky2d(double s, double h, double hinv);
  inline double sph_kernel_spiky2d_der(double s, double h, double hinv);
  inline double sph_kernel_spiky2d_cut();
}

/* ----------------------------------------------------------------------
   Spiky SPH kernel 2D
   h is kernel parameter
   s is distance normalized by h
   0.09947183943 is 5 over 16*pi
------------------------------------------------------------------------- */

inline double SPH_KERNEL_NS::sph_kernel_spiky2d(double s, double, double hinv)
{
  if (s < 2.)
  {
      return (0.09947183943*hinv*hinv * (2.-s)*(2.-s)*(2.-s));
  }
  else
  {
      return 0;
  }
}

/* ----------------------------------------------------------------------
   Derivative of Spiky SPH kernel 2D
   is equal to grad W if multiplied with radial unit vector
   h is kernel parameter
   s is distance normalized by h
   0.298415518297304 is 15 over 16*pi
------------------------------------------------------------------------- */

inline double SPH_KERNEL_NS::sph_kernel_spiky2d_der(double s, double, double hinv)
{
    if (s < 2.)
    {
        return (-0.298415518297304*hinv*hinv*hinv * (2.-s)*(2.-s));
    }
    else
    {
        return 0;
    }
}

/* ----------------------------------------------------------------------
   Definition of Spiky SPH 2D kernel cutoff in terms of s
   s is normalized distance
------------------------------------------------------------------------- */

inline double SPH_KERNEL_NS::sph_kernel_spiky2d_cut()
{
  return 2.;
}

#endif
#endif

