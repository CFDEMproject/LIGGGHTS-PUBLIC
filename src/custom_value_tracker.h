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
        CustomValueTracker(LAMMPS *lmp,AbstractMesh &_owner);
        ~CustomValueTracker();

        // per mesh-element properties

        template<typename T>
        T* addElementProperty(char *_id, char* _comm, char* _ref,char *_restart,int _scalePower = 1);

        template<typename T>
        T* getElementProperty(char *_id);

        template<typename T, typename U>
        void setElementProperty(char *_id, U def);

        void removeElementProperty(char *_id);

        // operation with
        // per mesh-element properties

        inline void deleteElement(int i);
        inline void deleteForwardElement(int i,bool scale,bool translate,bool rotate);
        inline void deleteRestartElement(int i,bool scale,bool translate,bool rotate);

        void move(double *delta);
        inline void moveElement(int i, double *delta);
        void rotate(double *dQ);
        void scale(double factor);

        // buffer operations

        inline int elemListBufSize(int n,int operation,bool scale,bool translate,bool rotate);
        inline int pushElemListToBuffer(int n, int *list, double *buf, int operation,bool scale,bool translate, bool rotate);
        inline int popElemListFromBuffer(int first, int n, double *buf, int operation,bool scale,bool translate, bool rotate);

        inline int elemBufSize(int operation,bool scale,bool translate,bool rotate);
        inline int pushElemToBuffer(int n, double *buf, int operation,bool scale,bool translate, bool rotate);
        inline int popElemFromBuffer(double *buf, int operation,bool scale,bool translate, bool rotate);

        inline int meshPropsBufSize(int operation,bool scale,bool translate,bool rotate);
        inline int pushMeshPropsToBuffer(double *buf, int operation,bool scale,bool translate, bool rotate);
        inline int popMeshPropsFromBuffer(double *buf, int operation,bool scale,bool translate, bool rotate);

        // scalar mesh properties

        template<typename T>
        T* addMeshProperty(char *_id, char* _comm, char* _ref,char *_restart, int _scalePower = 1);

        template<typename T>
        T* getMeshProperty(char *_id);

        template<typename T, typename U>
        void setMeshProperty(char *_id, U def);

        void removeMeshProperty(char *_id);

        int getCapacity();
        inline void grow(int to);

      private:

        class AbstractMesh &owner_;
        int capacityElement_; 
        class AssociativePointerArray<ContainerBase> elementProperties_;
        class AssociativePointerArray<ContainerBase> meshProperties_;
  };

  // *************************************
  #include "custom_value_tracker_I.h"
  // *************************************

} /* LAMMPS_NS */
#endif /* CUSTOMVALUETRACKER_H_FOO */
