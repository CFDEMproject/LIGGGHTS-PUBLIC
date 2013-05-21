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

#ifndef LMP_ASSOCIATIVE_POINTER_ARRAY_H
#define LMP_ASSOCIATIVE_POINTER_ARRAY_H

#include <string.h>
#include "memory.h"

namespace LAMMPS_NS
{
  #define ID_LEN 100

template<typename T>
class AssociativePointerArray
{
      public:
        AssociativePointerArray();
        ~AssociativePointerArray();

        template <typename U>
        U* add(char *_id, char* _comm, char* _ref, char *_restart,int _scalePower = 1);

        void remove(char *_id);

        template <typename U>
        U* getPointerById(char *_id);

        T* getBasePointerById(char *_id);

        template <typename U>
        U* getPointerByIndex(int i);

        T* getBasePointerByIndex(int i);

        void grow(int to);

        int size();

        inline void copyElement(int from, int to);
        inline void addZeroElement();
        inline void deleteElement(int n);
        inline void deleteForwardElement(int n,bool scale,bool translate,bool rotate);
        inline void deleteRestartElement(int n,bool scale,bool translate,bool rotate);

        inline void clearReverse(bool scale,bool translate,bool rotate);

        inline void storeOrig(class AssociativePointerArray &orig);
        inline void storeOrig(char *_id,class AssociativePointerArray &orig);
        inline bool reset(class AssociativePointerArray &orig);
        inline bool reset(char *_id,class AssociativePointerArray &orig);

        void rotate(double *dQ);
        void move(double *delta);
        void moveElement(int i,double *delta);
        void scale(double factor);

        inline int bufSize(int operation,bool scale,bool translate,bool rotate);
        inline int pushToBuffer(double *buf, int operation,bool scale,bool translate, bool rotate);
        inline int popFromBuffer(double *buf, int operation,bool scale,bool translate, bool rotate);

        inline int elemListBufSize(int n,int operation,bool scale,bool translate,bool rotate);
        inline int pushElemListToBuffer(int n, int *list, double *buf, int operation,bool scale,bool translate, bool rotate);
        inline int popElemListFromBuffer(int first, int n, double *buf, int operation,bool scale,bool translate, bool rotate);
        inline int pushElemListToBufferReverse(int first, int n, double *buf, int operation,bool scale,bool translate, bool rotate);
        inline int popElemListFromBufferReverse(int n, int *list, double *buf, int operation,bool scale,bool translate, bool rotate);

        inline int elemBufSize(int operation,bool scale,bool translate,bool rotate);
        inline int pushElemToBuffer(int n, double *buf, int operation,bool scale,bool translate, bool rotate);
        inline int popElemFromBuffer(double *buf, int operation,bool scale,bool translate, bool rotate);

        int idToIndex(char *_id);
        void indexToId(int index, char *_id);

      private:

        T **content_;
        int numElem_, maxElem_;

        void growArrays();
};

  // *************************************
  #include "associative_pointer_array_I.h"
  // *************************************

} /* LAMMPS_NS */
#endif /* ASSOCIATIVEPOINTERARRAY_H_ */
