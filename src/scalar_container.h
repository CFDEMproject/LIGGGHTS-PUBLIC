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

#ifndef LMP_SCALAR_CONTAINER
#define LMP_SCALAR_CONTAINER

#include "general_container.h"
#include "memory.h"
#include <limits>

namespace LAMMPS_NS
{
  template<typename T>
  class ScalarContainer : public GeneralContainer <T, 1, 1>
  {
    public:
          ScalarContainer();
          ScalarContainer(const char *_id);
          ScalarContainer(const char *_id, const char *_comm, const char *_ref, const char *_restart, int _scalePower = 1);
          ScalarContainer(ScalarContainer<T> const &orig);
          virtual ~ScalarContainer();

          void add(T elem);
          T get(int n);
          void set(int n, T elem);
          void setAll(T def);
          T* begin();
          void* begin_slow_dirty();
          T& operator() (int n);
          T const& operator() (int n) const;
          T max();
          T max(int);
          int pushToBuffer_plain(double *buf);
          int pullFromBuffer_plain(double *buf);
  };

  /* ----------------------------------------------------------------------
   constructors
  ------------------------------------------------------------------------- */

  template<typename T>
  ScalarContainer<T>::ScalarContainer()
  : GeneralContainer<T,1,1>(0)
  {

  }

  template<typename T>
  ScalarContainer<T>::ScalarContainer(const char *_id)
  : GeneralContainer<T,1,1>(_id)
  {

  }

  template<typename T>
  ScalarContainer<T>::ScalarContainer(const char *_id, const char *_comm, const char *_ref, const char *_restart, int _scalePower)
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
      int len = this->size();
      for(int n = 0; n < len; n++)
          this->arr_[n][0][0] = def;
  }

  template<typename T>
  T* ScalarContainer<T>::begin()
  {
          return &(this->arr_[0][0][0]);
  }

  template<typename T>
  void* ScalarContainer<T>::begin_slow_dirty()
  {
          return (void*) &(this->arr_[0][0][0]);
  }

  template<typename T>
  T ScalarContainer<T>::max()
  {
      int len = this->size();

      if(len == 0)
        return  (std::numeric_limits<T>::min)();

      T maxim = this->arr_[0][0][0];

      for(int n = 0; n < len; n++)
          if(this->arr_[n][0][0] > maxim)
            maxim = this->arr_[n][0][0];

      return maxim;
  }

  template<typename T>
  T ScalarContainer<T>::max(int to)
  {
      T maxim;

      int nn = MathExtraLiggghts::min(to,this->size());

      if(nn == 0)
        return  (std::numeric_limits<T>::min)();

      maxim = this->arr_[0][0][0];

      for(int n = 1; n < nn; n++)
          if(this->arr_[n][0][0] > maxim)
            maxim = this->arr_[n][0][0];

      return maxim;
  }

  template<typename T>
  int ScalarContainer<T>::pushToBuffer_plain(double *buf)
  {
      int len = this->size();

      int m = 0;

      for(int i = 0; i < len; i++)
        buf[m++] = static_cast<double>(this->arr_[i][0][0]);

      return len;
  }

  template<typename T>
  int ScalarContainer<T>::pullFromBuffer_plain(double *buf)
  {
      int len = this->size();

      int m = 0;

      for(int i = 0; i < len; i++)
        this->arr_[i][0][0] = static_cast<T>(buf[m++]);

      return len;
  }
} /* LAMMPS_NS */
#endif /* SCALARCONTAINER_H_ */
