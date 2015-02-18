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
