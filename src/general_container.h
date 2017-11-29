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

    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Philippe Seil (JKU Linz)
    Andreas Aigner (DCS Computing GmbH, Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef LMP_GENERAL_CONTAINER
#define LMP_GENERAL_CONTAINER

#include "container_base.h"
#include "memory_ns.h"
#include "math_extra_liggghts.h"
#include "domain.h"
#include <limits>
#include <cmath>
#include <algorithm>

inline int GROW_CONTAINER()
{ return 10000; }

using namespace LAMMPS_MEMORY_NS;

namespace LAMMPS_NS
{

  template<typename T, int NUM_VEC, int LEN_VEC>
  class GeneralContainer : public ContainerBase
  {
      public:

          bool isDoubleData();
          bool isIntData();

          bool subtract (GeneralContainer<T,NUM_VEC,LEN_VEC> const &A,
                         GeneralContainer<T,NUM_VEC,LEN_VEC> const &minusB);

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
          void setAllToZero()
          { setAll(static_cast<T>(0)); }

          void set(int i, T** elem);
          void set(int i, int j, T* elem);

          bool setFromContainer(ContainerBase *cont);

          bool calcAvgFromContainer();
          bool calcMeanSquareFromContainer();
          bool calcSumFromContainer();

          T max_scalar();
          T min_scalar();

          T**& operator()(int n);
          T** const& operator()(int n) const;
          T*** begin();
          virtual void* begin_slow_dirty();

          inline void scale(double factor);
          inline void move(const double * const dx);
          inline void moveElement(const int i, const double * const dx);
          inline void rotate(const double * const dQ);

          // all push and pop functions return number of bytes taken from / added to buf
          // all push and pop functions expect buf to point to first element with usable data

          // push / pop all elements
          
          inline int bufSize(int operation = OPERATION_UNDEFINED,
                            bool scale=false,bool translate=false, bool rotate=false) const;
          inline int pushToBuffer(double *buf, int operation = OPERATION_UNDEFINED,
                           bool scale=false,bool translate=false, bool rotate=false);
          inline int popFromBuffer(double *buf, int operation = OPERATION_UNDEFINED,
                           bool scale=false,bool translate=false, bool rotate=false);

          // push / pop a list elements
          
          inline int elemListBufSize(int n, int operation = OPERATION_UNDEFINED,
                            bool scale=false,bool translate=false, bool rotate=false);
          inline int pushElemListToBuffer(int n, int *list, int *wraplist, double *buf, int operation, double *dlo, double *dhi,
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

          inline int size() const
          { return numElem_; }

          inline int nVec() const
          { return NUM_VEC; }

          inline int lenVec() const
          { return LEN_VEC; }

          inline int capacity() const
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
