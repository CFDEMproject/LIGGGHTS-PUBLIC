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

#ifndef LMP_CUSTOM_VALUE_TRACKER_H
#define LMP_CUSTOM_VALUE_TRACKER_H

#include "associative_pointer_array.h"
#include "container.h"
#include "abstract_mesh.h"

namespace LAMMPS_NS
{
  class CustomValueTracker : protected Pointers
  {
      public:
        CustomValueTracker(LAMMPS *lmp);
        CustomValueTracker(LAMMPS *lmp,AbstractMesh *_owner);
        ~CustomValueTracker();

        // per-element properties

        template<typename T>
        T* addElementProperty(char *_id, char* _comm, char* _ref,char *_restart,int _scalePower = 1);

        template<typename T>
        T* getElementProperty(char *_id);

        inline ContainerBase* getElementPropertyBase(char *_id);

        inline int getElementPropertyIndex(char *_id);

        template<typename T, typename U>
        void setElementProperty(char *_id, U def);

        void removeElementProperty(char *_id);

        // global (e.g. mesh) properties

        template<typename T>
        T* addGlobalProperty(char *_id, char* _comm, char* _ref,char *_restart, int _scalePower = 1);

        template<typename T>
        T* getGlobalProperty(char *_id);

        template<typename T, typename U>
        void setGlobalProperty(char *_id, U def);

        void removeGlobalProperty(char *_id);

        // operation with
        // per-element properties

        inline void copyElement(int from, int to);
        inline void addZeroElement();
        inline void deleteElement(int i);
        inline void deleteForwardElement(int i,bool scale,bool translate,bool rotate);
        inline void deleteRestartElement(int i,bool scale,bool translate,bool rotate);
        void clearReverse(bool scale,bool translate,bool rotate);

        void storeOrig();
        void resetToOrig();

        inline void storeGlobalPropOrig(char *_id);
        inline void resetGlobalPropToOrig(char *_id);

        inline void moveElement(int i, double *delta);
        void move(double *vecTotal, double *vecIncremental);
        void move(double *vecIncremental);
        void rotate(double *totalQ, double *dQ);
        void rotate(double *dQ);
        void scale(double factor);

        // buffer operations

        inline int elemListBufSize(int n,int operation,bool scale,bool translate,bool rotate);
        inline int pushElemListToBuffer(int n, int *list, double *buf, int operation,bool scale,bool translate, bool rotate);
        inline int popElemListFromBuffer(int first, int n, double *buf, int operation,bool scale,bool translate, bool rotate);
        inline int pushElemListToBufferReverse(int first, int n, double *buf, int operation,bool scale,bool translate, bool rotate);
        inline int popElemListFromBufferReverse(int n, int *list, double *buf, int operation,bool scale,bool translate, bool rotate);

        inline int elemBufSize(int operation,bool scale,bool translate,bool rotate);
        inline int pushElemToBuffer(int n, double *buf, int operation,bool scale,bool translate, bool rotate);
        inline int popElemFromBuffer(double *buf, int operation,bool scale,bool translate, bool rotate);

        inline int globalPropsBufSize(int operation,bool scale,bool translate,bool rotate);
        inline int pushGlobalPropsToBuffer(double *buf, int operation,bool scale,bool translate, bool rotate);
        inline int popGlobalPropsFromBuffer(double *buf, int operation,bool scale,bool translate, bool rotate);

        // mem managenement

        int getCapacity();
        inline void grow(int to);

      private:

        class AbstractMesh *ownerMesh_;

        int capacityElement_; 
        class AssociativePointerArray<ContainerBase> elementProperties_;
        class AssociativePointerArray<ContainerBase> globalProperties_;
        class AssociativePointerArray<ContainerBase> globalProperties_orig_;
  };

  // *************************************
  #include "custom_value_tracker_I.h"
  // *************************************

} /* LAMMPS_NS */
#endif
