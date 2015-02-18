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
   Philippe Seil (JKU Linz)
------------------------------------------------------------------------- */

#ifndef LMP_VECTOR_CONTAINER
#define LMP_VECTOR_CONTAINER

#include "general_container.h"
#include "memory.h"

namespace LAMMPS_NS
{
  template<typename T, int LEN_VEC>
  class VectorContainer : public GeneralContainer <T, 1, LEN_VEC>
  {
    public:
          VectorContainer(const char *_id);
          VectorContainer(const char *_id, const char *_comm, const char *_ref, const char *_restart, int _scalePower = 1);
          VectorContainer(VectorContainer<T,LEN_VEC> const &orig);
          virtual ~VectorContainer();

          void add(T* elem);
          void get(int n, T* elem);
          void set(int n, T* elem);

          T max_elem(int n);
          T min_elem(int n);

          //void setAll(T def);
          T*& operator() (int n);
          T* const& operator() (int n) const;

          T** begin();
          void* begin_slow_dirty();
  };

  /* ----------------------------------------------------------------------
   constructors
  ------------------------------------------------------------------------- */

  template<typename T, int LEN_VEC>
  VectorContainer<T,LEN_VEC>::VectorContainer(const char *_id)
  : GeneralContainer<T,1,LEN_VEC>(_id)
  {

  }

  template<typename T, int LEN_VEC>
  VectorContainer<T,LEN_VEC>::VectorContainer(const char *_id, const char *_comm, const char *_ref, const char *_restart, int _scalePower)
  : GeneralContainer<T,1,LEN_VEC>(_id, _comm, _ref, _restart, _scalePower)
  {

  }

  template<typename T, int LEN_VEC>
  VectorContainer<T,LEN_VEC>::VectorContainer(VectorContainer<T,LEN_VEC> const &orig)
  : GeneralContainer<T,1,LEN_VEC>(orig)
  {

  }

  /* ----------------------------------------------------------------------
   destructor
  ------------------------------------------------------------------------- */

  template<typename T, int LEN_VEC>
  VectorContainer<T,LEN_VEC>::~VectorContainer()
  {

  }

  /* ----------------------------------------------------------------------
   add element
  ------------------------------------------------------------------------- */

  template<typename T, int LEN_VEC>
  void VectorContainer<T,LEN_VEC>::add(T* elem)
  {
          if(GeneralContainer<T,1,LEN_VEC>::numElem_ == GeneralContainer<T,1,LEN_VEC>::maxElem_)
          {
                  grow(GeneralContainer<T,1,LEN_VEC>::arr_,GeneralContainer<T,1,LEN_VEC>::maxElem_+GROW,1,LEN_VEC);
                  GeneralContainer<T,1,LEN_VEC>::maxElem_ += GROW;
          }
          for(int i=0;i<LEN_VEC;i++)
                  GeneralContainer<T,1,LEN_VEC>::arr_[GeneralContainer<T,1,LEN_VEC>::numElem_][0][i] = elem[i];

          GeneralContainer<T,1,LEN_VEC>::numElem_++;
  }

  /* ----------------------------------------------------------------------
   access
  ------------------------------------------------------------------------- */

  template<typename T, int LEN_VEC>
  T*& VectorContainer<T,LEN_VEC>::operator() (int n)
  {
          return GeneralContainer<T,1,LEN_VEC>::arr_[n][0];
  }

  template<typename T, int LEN_VEC>
  T* const& VectorContainer<T,LEN_VEC>::operator() (int n) const
  {
          return GeneralContainer<T,1,LEN_VEC>::arr_[n][0];
  }

  template<typename T, int LEN_VEC>
  void VectorContainer<T,LEN_VEC>::get(int n, T* elem)
  {
          for(int i = 0; i < LEN_VEC; i++)
                  elem[i] = GeneralContainer<T,1,LEN_VEC>::arr_[n][0][i];
  }

  template<typename T, int LEN_VEC>
  void VectorContainer<T,LEN_VEC>::set(int n, T* elem)
  {
          for(int i = 0; i < LEN_VEC; i++)
                  GeneralContainer<T,1,LEN_VEC>::arr_[n][0][i] = elem[i];
  }

  template<typename T, int LEN_VEC>
  T VectorContainer<T,LEN_VEC>::max_elem(int n)
  {
          T max = GeneralContainer<T,1,LEN_VEC>::arr_[n][0][0];

          for(int i = 1; i < LEN_VEC; i++)
                  if(GeneralContainer<T,1,LEN_VEC>::arr_[n][0][i] > max)
                    max = GeneralContainer<T,1,LEN_VEC>::arr_[n][0][i];

          return max;
  }

  template<typename T, int LEN_VEC>
  T VectorContainer<T,LEN_VEC>::min_elem(int n)
  {
          T min = GeneralContainer<T,1,LEN_VEC>::arr_[n][0][0];

          for(int i = 1; i < LEN_VEC; i++)
                  if(GeneralContainer<T,1,LEN_VEC>::arr_[n][0][i] < min)
                    min = GeneralContainer<T,1,LEN_VEC>::arr_[n][0][i];

          return min;
  }
/*
  template<typename T, int LEN_VEC>
  void VectorContainer<T,LEN_VEC>::setAll(T def)
  {
      int len = this->size();
      for(int n = 0; n < len; n++)
          for(int i = 0; i < LEN_VEC; i++)
                  GeneralContainer<T,1,LEN_VEC>::arr_[n][0][i] = def;
  }
*/
  template<typename T, int LEN_VEC>
  T** VectorContainer<T,LEN_VEC>::begin()
  {
          return &(GeneralContainer<T,1,LEN_VEC>::arr_[0][0]);
  }

  template<typename T, int LEN_VEC>
  void* VectorContainer<T,LEN_VEC>::begin_slow_dirty()
  {
          return (void*) &(GeneralContainer<T,1,LEN_VEC>::arr_[0][0]);
  }

} /* LAMMPS_NS */
#endif /* VECTORCONTAINER_H_ */
