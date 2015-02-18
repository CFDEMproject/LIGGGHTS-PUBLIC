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
        T* addElementProperty(const char *_id, const char* _comm, const char* _ref, const char *_restart, int _scalePower = 1, int _init_len = 0);

        template<typename T>
        T* getElementProperty(const char *_id);

        inline ContainerBase* getElementPropertyBase(const char *_id);
        inline ContainerBase* getElementPropertyBase(int i);

        inline int getElementPropertyIndex(const char *_id);

        template<typename T, typename U>
        void setElementProperty(const char *_id, U def);

        void removeElementProperty(const char *_id);

        void check_element_property_consistency(int _len);

        // global (e.g. mesh) properties

        template<typename T>
        T* addGlobalProperty(const char *_id, const char* _comm, const char* _ref, const char *_restart, int _scalePower = 1);

        template<typename T>
        T* getGlobalProperty(const char *_id);

        template<typename T, typename U>
        void setGlobalProperty(const char *_id, U def);

        void removeGlobalProperty(const char *_id);

        // operation with
        // per-element properties

        inline void copyElement(int from, int to);
        inline void addUninitializedElement();
        inline void addZeroElement();
        inline void deleteElement(int i);
        inline void deleteForwardElement(int i,bool scale,bool translate,bool rotate);
        inline void deleteRestartElement(int i,bool scale,bool translate,bool rotate);
        inline void deleteRestartGlobal(bool scale,bool translate,bool rotate);
        void clearReverse(bool scale,bool translate,bool rotate);

        void storeOrig();
        void resetToOrig();

        inline void storeGlobalPropOrig(const char *_id);
        inline void resetGlobalPropToOrig(const char *_id);

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
