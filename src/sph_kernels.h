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

    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef LMP_SPH_KERNELS
#define LMP_SPH_KERNELS

#include "style_sph_kernel.h"

namespace SPH_KERNEL_NS {
  inline int sph_kernels_unique_id();
  inline int sph_kernel_id(char *style);
  inline double sph_kernel(int id,double s,double h,double hinv);
  inline double sph_kernel_der(int id,double s,double h,double hinv);
  inline double sph_kernel_cut(int id);
}

/* ---------------------------------------------------------------------- */

inline int SPH_KERNEL_NS::sph_kernels_unique_id()
{
  int ids[50];
  int nkernels = 0;

  if (0) return 0.;
  #define SPH_KERNEL_CLASS
  #define SPHKernel(kernel_id,kernelstyle,SPHKernelCalculation,SPHKernelCalculationDer,SPHKernelCalculationCut) \
  ids[nkernels++] = kernel_id;
  #include "style_sph_kernel.h"
  #undef SPH_KERNEL_CLASS
  #undef SPHKernel

  // check if double ids
  for (int i = 0; i < nkernels; i++)
      for (int j = i+1; j < nkernels; j++)
         if(ids[i] == ids[j]) return -1;

  return 0;
}

/* ---------------------------------------------------------------------- */

inline int SPH_KERNEL_NS::sph_kernel_id(char *style)
{
  if (0) return 0.;
  #define SPH_KERNEL_CLASS
  #define SPHKernel(kernel_id,kernelstyle,SPHKernelCalculation,SPHKernelCalculationDer,SPHKernelCalculationCut) \
  else if (strcmp(style,#kernelstyle) == 0) return kernel_id;
  #include "style_sph_kernel.h"
  #undef SPH_KERNEL_CLASS
  #undef SPHKernel
  return -1;
}

/* ---------------------------------------------------------------------- */

inline double SPH_KERNEL_NS::sph_kernel(int id,double s,double h,double hinv)
{
  if (0) return 0.;
  #define SPH_KERNEL_CLASS
  #define SPHKernel(kernel_id,kernelstyle,SPHKernelCalculation,SPHKernelCalculationDer,SPHKernelCalculationCut) \
  else if (kernel_id == id) return SPH_KERNEL_NS::SPHKernelCalculation(s,h,hinv);
  #include "style_sph_kernel.h"
  #undef SPH_KERNEL_CLASS
  #undef SPHKernel
  return 0.;
}

/* ---------------------------------------------------------------------- */

inline double SPH_KERNEL_NS::sph_kernel_der(int id,double s,double h,double hinv)
{
  if (0) return 0.;
  #define SPH_KERNEL_CLASS
  #define SPHKernel(kernel_id,kernelstyle,SPHKernelCalculation,SPHKernelCalculationDer,SPHKernelCalculationCut) \
  else if (kernel_id == id) return SPH_KERNEL_NS::SPHKernelCalculationDer(s,h,hinv);
  #include "style_sph_kernel.h"
  #undef SPH_KERNEL_CLASS
  #undef SPHKernel
  return 0.;
}

/* ---------------------------------------------------------------------- */

inline double SPH_KERNEL_NS::sph_kernel_cut(int id)
{
  if (0) return 0.;
  #define SPH_KERNEL_CLASS
  #define SPHKernel(kernel_id,kernelstyle,SPHKernelCalculation,SPHKernelCalculationDer,SPHKernelCalculationCut) \
  else if (kernel_id == id) return SPH_KERNEL_NS::SPHKernelCalculationCut();
  #include "style_sph_kernel.h"
  #undef SPH_KERNEL_CLASS
  #undef SPHKernel
  return 0.;
}

#endif
