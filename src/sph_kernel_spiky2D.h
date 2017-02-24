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
    Andreas Aigner (JKU Linz)
    Markus Schoergenhumer (JKU Linz, 2D kernel)

    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

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

