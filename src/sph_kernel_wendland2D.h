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

#ifdef SPH_KERNEL_CLASS

    // kernel identifier (a unique integer >= 0)
    // a name for the kernel
    // name of the functions for the kernel, its derivative, and the cutoff are defined
    SPHKernel
    (
        5,
        wendland2d,
        sph_kernel_wendland2d,
        sph_kernel_wendland2d_der,
        sph_kernel_wendland2d_cut
    )

#else

#ifndef LMP_SPH_KERNEL_WENDLAND2D
#define LMP_SPH_KERNEL_WENDLAND2D

namespace SPH_KERNEL_NS {
  inline double sph_kernel_wendland2d(double s, double h, double hinv);
  inline double sph_kernel_wendland2d_der(double s, double h, double hinv);
  inline double sph_kernel_wendland2d_cut();
}

/* ----------------------------------------------------------------------
   Wendland SPH kernel
   h is kernel parameter
   s is distance normalized by h
   0.557042301 is 7 over 4pi
------------------------------------------------------------------------- */

inline double SPH_KERNEL_NS::sph_kernel_wendland2d(double s, double, double hinv)
{
    return (0.557042301*hinv*hinv * (1.-0.5*s)*(1.-0.5*s)*(1.-0.5*s)*(1.-0.5*s) * (2.*s+1));
}

/* ----------------------------------------------------------------------
   Derivative of Wendland SPH kernel
   is equal to grad W if multiplied with radial unit vector
   h is kernel parameter
   s is distance normalized by h
   1.114084602 is 14 over 4pi
------------------------------------------------------------------------- */

inline double SPH_KERNEL_NS::sph_kernel_wendland2d_der(double s, double, double hinv)
{
    return (1.114084602*hinv*hinv*hinv * ((1-0.5*s)*(1-0.5*s)*(1-0.5*s)*((1-0.5*s)-(2.*s+1))));
}

/* ----------------------------------------------------------------------
   Definition of Wendland SPH kernel cutoff in terms of s
   s is normalized distance
------------------------------------------------------------------------- */

inline double SPH_KERNEL_NS::sph_kernel_wendland2d_cut()
{
    return 2.;
}

#endif
#endif

