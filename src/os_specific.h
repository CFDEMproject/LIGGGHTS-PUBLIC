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
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */
#include <stdlib.h>

#if defined(_WIN32) || defined(_WIN64)
#include <malloc.h>
static inline void * aligned_malloc(size_t alignment, size_t size)
{
  return _aligned_malloc(size, alignment);
}

static inline void aligned_free(void * ptr)
{
  _aligned_free(ptr);
}

template<typename T>
static inline T * aligned_malloc(size_t alignment)
{
  return (T*)_aligned_malloc(sizeof(T), alignment);
}
#else
static inline void * aligned_malloc(size_t alignment, size_t size)
{
  void * ptr = NULL;
  posix_memalign(&ptr, alignment, size);
  return ptr;
}

static inline void aligned_free(void * ptr)
{
  free(ptr);
}

template<typename T>
static inline T * aligned_malloc(size_t alignment)
{
  void * ptr = NULL;
  posix_memalign(&ptr, alignment, sizeof(T));
  return (T*)ptr;
}
#endif
