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

    Andreas Aigner (JKU Linz))

    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifdef SPH_KERNEL_CLASS

    // kernel identifier (a unique integer >= 0)
    // a name for the kernel
    // name of the functions for the kernel, its derivative, and the cutoff are defined
    SPHKernel
    (
        4,
        spiky,
        sph_kernel_spiky,
        sph_kernel_spiky_der,
        sph_kernel_spiky_cut
    )

#else

#ifndef LMP_SPH_KERNEL_SPIKY
#define LMP_SPH_KERNEL_SPIKY

namespace SPH_KERNEL_NS {
  inline double sph_kernel_spiky(double s, double h, double hinv);
  inline double sph_kernel_spiky_der(double s, double h, double hinv);
  inline double sph_kernel_spiky_cut();
}

/* ----------------------------------------------------------------------
   Cubic spline SPH kernel
   h is kernel parameter
   s is distance normalized by h
   4.7746 is 15 over pi
   0.07460388 is 15 over 64*pi
------------------------------------------------------------------------- */

inline double SPH_KERNEL_NS::sph_kernel_spiky(double s, double, double hinv)
{
  if (s < 2.)
  {
      return (0.07460388*hinv*hinv*hinv * (2.-s)*(2.-s)*(2.-s));
  }
  else
  {
      return 0;
  }
  /*
    if (s < 1.)
    {
        return (4.7746*hinv*hinv*hinv * (1.-s)*(1.-s)*(1.-s));
    }
    else
    {
        return 0;
    }
  */
}

/* ----------------------------------------------------------------------
   Derivative of cubic spline SPH kernel
   is equal to grad W if multiplied with radial unit vector
   h is kernel parameter
   s is distance normalized by h
   14.323944876 is 45 over pi
   0.223811639 is 45 over 64*pi
------------------------------------------------------------------------- */

inline double SPH_KERNEL_NS::sph_kernel_spiky_der(double s, double, double hinv)
{
    if (s < 2.)
    {
        return (-0.223811639*hinv*hinv*hinv*hinv * (2.-s)*(2.-s));
    }
    else
    {
        return 0;
    }
  /*
    if (s < 1.)
    {
        return (-14.324*hinv*hinv*hinv*hinv * (1.-s)*(1.-s));
    }
    else
    {
        return 0;
    }
  */
}

/* ----------------------------------------------------------------------
   Definition of cubic spline SPH kernel cutoff in terms of s
   s is normalized distance
------------------------------------------------------------------------- */

inline double SPH_KERNEL_NS::sph_kernel_spiky_cut()
{
  return 2.;
  /*
    return 1.;
  */
}

#endif
#endif

