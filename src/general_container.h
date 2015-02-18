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

#ifndef LMP_GENERAL_CONTAINER
#define LMP_GENERAL_CONTAINER

#include "container_base.h"
#include "memory_ns.h"
#include "math_extra_liggghts.h"
#include <string.h>

#define GROW 100

using namespace LAMMPS_MEMORY_NS;

namespace LAMMPS_NS
{

  template<typename T, int NUM_VEC, int LEN_VEC>
  class GeneralContainer : public ContainerBase
  {
      public:

          bool isDoubleData();
          bool isIntData();

          void add(T** elem);
          void addZero();

          void copy(int from,int to);
          void del(int n);
          void delForward(int n,bool scale,bool translate,bool rotate);
          void delRestart(int n,bool scale,bool translate,bool rotate);
          void delRestart(bool scale,bool translate,bool rotate);
          void clearReverse(bool scale,bool translate,bool rotate);

          void get(int n, T** elem);

          void setToDefault(int n);
          void setAll(T def);
          void setAll(int to, T def);
          void set(int i, T** elem);
          void set(int i, int j, T* elem);

          bool setFromContainer(ContainerBase *cont);

          T max_scalar();
          T min_scalar();

          T**& operator()(int n);
          T** const& operator()(int n) const;
          T*** begin();
          virtual void* begin_slow_dirty();

          inline void scale(double factor);
          inline void move(double *dx);
          inline void moveElement(int i,double *dx);
          inline void rotate(double *dQ);

          // all push and pop functions return number of bytes taken from / added to buf
          // all push and pop functions expect buf to point to first element with usable data

          // push / pop all elements
          
          inline int bufSize(int operation = OPERATION_UNDEFINED,
                            bool scale=false,bool translate=false, bool rotate=false);
          inline int pushToBuffer(double *buf, int operation,
                           bool scale=false,bool translate=false, bool rotate=false);
          inline int popFromBuffer(double *buf, int operation,
                           bool scale=false,bool translate=false, bool rotate=false);

          // push / pop a list elements
          
          inline int elemListBufSize(int n, int operation = OPERATION_UNDEFINED,
                            bool scale=false,bool translate=false, bool rotate=false);
          inline int pushElemListToBuffer(int n, int *list, double *buf, int operation,
                           bool scale=false,bool translate=false, bool rotate=false);
          inline int popElemListFromBuffer(int first, int n, double *buf, int operation,
                           bool scale=false,bool translate=false, bool rotate=false);
          inline int pushElemListToBufferReverse(int first, int n, double *buf, int operation,
                           bool scale=false,bool translate=false, bool rotate=false);
          inline int popElemListFromBufferReverse(int n, int *list, double *buf, int operation,
                           bool scale=false,bool translate=false, bool rotate=false);

          // push / pop one single element
          
          inline int elemBufSize(int operation = OPERATION_UNDEFINED,
                            bool scale=false,bool translate=false, bool rotate=false);
          inline int pushElemToBuffer(int i, double *buf, int operation,
                           bool scale=false,bool translate=false, bool rotate=false);
          inline int popElemFromBuffer(double *buf, int operation,
                           bool scale=false,bool translate=false, bool rotate=false);

          void addUninitialized(int n);

          inline int size()
          { return numElem_; }

          inline int nVec()
          { return NUM_VEC; }

          inline int lenVec()
          { return LEN_VEC; }

          inline int capacity()
          { return maxElem_; }

          inline void clearContainer()
          { numElem_ = 0; }

          inline void setDefaultValue(T val)
          { defaultValue_ = val; useDefault_ = true; }

      protected:

          GeneralContainer(const char *_id);
          GeneralContainer(const char *_id, const char *_comm, const char *_ref, const char *_restart, int _scalePower = 1);
          GeneralContainer(GeneralContainer<T,NUM_VEC,LEN_VEC> const &orig);
          virtual ~GeneralContainer();

          // shall return the size of an entry in bytes
          int getElemSize();

          int numElem_, maxElem_;

          T*** arr_;

          T defaultValue_;
  };

  // *************************************
  #include "general_container_I.h"
  // *************************************

} /* LAMPPS_NS */
#endif /* CONTAINERBASE_H_ */
