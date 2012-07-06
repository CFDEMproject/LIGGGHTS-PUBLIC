/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
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

#ifndef LMP_SCALAR_CONTAINER
#define LMP_SCALAR_CONTAINER

#include "general_container.h"
#include "memory.h"

namespace LAMMPS_NS
{
  template<typename T>
  class ScalarContainer : public GeneralContainer <T, 1, 1>
  {
    public:
          ScalarContainer();
          ScalarContainer(char *_id, char *_comm, char *_ref, char *_restart, int _scalePower = 1);
          ScalarContainer(ScalarContainer<T> const &orig);
          virtual ~ScalarContainer();

          void add(T elem);
          T get(int n);
          void set(int n, T elem);
          void setAll(T def);
          T* begin();
          T& operator() (int n);
          T const& operator() (int n) const;
          T max();
  };

  /* ----------------------------------------------------------------------
   constructors
  ------------------------------------------------------------------------- */

  template<typename T>
  ScalarContainer<T>::ScalarContainer()
  : GeneralContainer<T,1,1>()
  {

  }

  template<typename T>
  ScalarContainer<T>::ScalarContainer(char *_id, char *_comm, char *_ref, char *_restart, int _scalePower)
  : GeneralContainer<T,1,1>(_id, _comm, _ref, _restart, _scalePower)
  {

  }

  template<typename T>
  ScalarContainer<T>::ScalarContainer(ScalarContainer<T> const &orig)
  : GeneralContainer<T,1,1>(orig)
  {

  }

  /* ----------------------------------------------------------------------
   destructor
  ------------------------------------------------------------------------- */

  template<typename T>
  ScalarContainer<T>::~ScalarContainer()
  {

  }

  /* ----------------------------------------------------------------------
   add element
  ------------------------------------------------------------------------- */

  template<typename T>
  void ScalarContainer<T>::add(T elem)
  {
          if(this->numElem_ == this->maxElem_)
          {
                  grow<T>(this->arr_,this->maxElem_+GROW,1,1);
                  this->maxElem_ += GROW;
          }
          this->arr_[this->numElem_][0][0] = elem;
          this->numElem_++;
  }

  /* ----------------------------------------------------------------------
   access
  ------------------------------------------------------------------------- */

  template<typename T>
  T& ScalarContainer<T>::operator() (int n)
  {
          return this->arr_[n][0][0];
  }

  template<typename T>
  T const& ScalarContainer<T>::operator() (int n) const
  {
          return this->arr_[n][0][0];
  }

  template<typename T>
  T ScalarContainer<T>::get(int n)
  {
          return this->arr_[n][0][0];
  }

  template<typename T>
  void ScalarContainer<T>::set(int n, T elem)
  {
          this->arr_[n][0][0] = elem;
  }

  template<typename T>
  void ScalarContainer<T>::setAll(T def)
  {
      for(int n=0;n<this->size();n++)
          this->arr_[n][0][0] = def;
  }

  template<typename T>
  T* ScalarContainer<T>::begin()
  {
          return &(this->arr_[0][0][0]);
  }

  template<typename T>
  T ScalarContainer<T>::max()
  {
      T maxim = this->arr_[0][0][0];

      for(int n=0;n<this->size();n++)
          if(this->arr_[n][0][0] > maxim)
            maxim = this->arr_[n][0][0];

      return maxim;
  }
} /* LAMMPS_NS */
#endif /* SCALARCONTAINER_H_ */
