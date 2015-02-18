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

#ifndef LMP_CUSTOM_VALUE_TRACKER_I_H
#define LMP_CUSTOM_VALUE_TRACKER_I_H

  /* ----------------------------------------------------------------------
   add property
  ------------------------------------------------------------------------- */

  template<typename T>
  T* CustomValueTracker::addElementProperty(const char *_id, const char* _comm, const char* _ref, const char *_restart, int _scalePower, int _init_len)
  {
     // error if property exists already
     if(elementProperties_.getPointerById<T>(_id))
     {
         char *errmsg = new char[strlen(_id)+200];
         sprintf(errmsg,"Illegal command, features are incompatible - element property '%s' exists already",_id);
         error->all(FLERR,errmsg);
         delete []errmsg;
     }

     // add property
     elementProperties_.add<T>(_id,_comm,_ref,_restart,_scalePower);

     // check if properties were set correctly
     // error here since ContainerBase not derived from Pointers
     if(!elementProperties_.getPointerById<T>(_id)->propertiesSetCorrectly())
     {
         char *errmsg = new char[strlen(_id)+200];
         sprintf(errmsg,"Illegal element property, comm or frame property not set correctly for property '%s'",_id);
         error->all(FLERR,errmsg);
         delete []errmsg;
     }

     // allocate memory and initialize
     
     if(ownerMesh_)
     {
        elementProperties_.getPointerById<T>(_id)->addUninitialized(ownerMesh_->sizeLocal()+ownerMesh_->sizeGhost());
     }
     if(_init_len > 0)
        elementProperties_.getPointerById<T>(_id)->addUninitialized(_init_len);

     elementProperties_.getPointerById<T>(_id)->setAll(0);

     // return pointer
     return elementProperties_.getPointerById<T>(_id);
  }

  template<typename T>
  T* CustomValueTracker::addGlobalProperty(const char *_id, const char* _comm, const char* _ref, const char *_restart, int _scalePower)
  {
     // error if property exists already
     if(globalProperties_.getPointerById<T>(_id))
     {
         char *errmsg = new char[strlen(_id)+200];
         sprintf(errmsg,"Illegal command, features are incompatible - global property '%s' already exists",_id);
         error->all(FLERR,errmsg);
         delete []errmsg;
     }

     // add property
     globalProperties_.add<T>(_id,_comm,_ref,_restart,_scalePower);
     globalProperties_orig_.add<T>(_id,_comm,_ref,_restart,_scalePower);

     // check if properties were set correctly
     // error here since ContainerBase not derived from Pointers
     if(!globalProperties_.getPointerById<T>(_id)->propertiesSetCorrectly())
     {
         char *errmsg = new char[strlen(_id)+200];
         sprintf(errmsg,"Illegal global property, comm or frame property not set correctly for property '%s'",_id);
         error->all(FLERR,errmsg);
         delete []errmsg;
     }

     // allocate memory
     //globalProperties_.getPointerById<T>(_id)->addUninitialized(capacityElement_);

     // return pointer
     return globalProperties_.getPointerById<T>(_id);
  }

  /* ----------------------------------------------------------------------
   mem management
  ------------------------------------------------------------------------- */

  void CustomValueTracker::grow(int to)
  {
      elementProperties_.grow(to);
      capacityElement_ = to;
  }

  /* ----------------------------------------------------------------------
   get reference
  ------------------------------------------------------------------------- */

  template<typename T>
  T* CustomValueTracker::getElementProperty(const char *_id)
  {
     return elementProperties_.getPointerById<T>(_id);
  }

  inline ContainerBase* CustomValueTracker::getElementPropertyBase(const char *_id)
  {
     return elementProperties_.getBasePointerById(_id);
  }

  inline ContainerBase* CustomValueTracker::getElementPropertyBase(int i)
  {
     return elementProperties_.getBasePointerByIndex(i);
  }

  inline int CustomValueTracker::getElementPropertyIndex(const char *_id)
  {
     return elementProperties_.idToIndex(_id);
  }

  template<typename T>
  T* CustomValueTracker::getGlobalProperty(const char *_id)
  {
     return globalProperties_.getPointerById<T>(_id);
  }

  /* ----------------------------------------------------------------------
   set property
  ------------------------------------------------------------------------- */

  template<typename T, typename U>
  void CustomValueTracker::setElementProperty(const char *_id, U def)
  {
     elementProperties_.getPointerById<T>(_id)->set(def);
  }

  template<typename T, typename U>
  void CustomValueTracker::setGlobalProperty(const char *_id, U def)
  {
     
     if(globalProperties_.getPointerById<T>(_id)->size() == 0)
        globalProperties_.getPointerById<T>(_id)->addUninitialized(1);
     globalProperties_.getPointerById<T>(_id)->set(0,def);

     if(globalProperties_orig_.getPointerById<T>(_id)->size() == 0)
        globalProperties_orig_.getPointerById<T>(_id)->addUninitialized(1);
     globalProperties_orig_.getPointerById<T>(_id)->set(0,def);
  }

  /* ----------------------------------------------------------------------
   store global property orig - only needs to be done manually for
   special cases, eg moving mesh ref points
  ------------------------------------------------------------------------- */

  inline void CustomValueTracker::storeGlobalPropOrig(const char *_id)
  {
      globalProperties_.storeOrig(_id,globalProperties_orig_);
  }

  /* ----------------------------------------------------------------------
   reset global property to orig - only needs to be done manually for
   special cases, eg moving mesh ref points
  ------------------------------------------------------------------------- */

  inline void CustomValueTracker::resetGlobalPropToOrig(const char *_id)
  {
      globalProperties_.reset(_id,globalProperties_orig_);
  }

  /* ----------------------------------------------------------------------
   copy data from element from to element to
  ------------------------------------------------------------------------- */

  void CustomValueTracker::copyElement(int from, int to)
  {
      elementProperties_.copyElement(from,to);
  }

  /* ----------------------------------------------------------------------
   add an element and initialize its properties with 0
  ------------------------------------------------------------------------- */

  void CustomValueTracker::addZeroElement()
  {
      elementProperties_.addZeroElement();
  }

  /* ----------------------------------------------------------------------
   add an uninitialized element
  ------------------------------------------------------------------------- */

  void CustomValueTracker::addUninitializedElement()
  {
      elementProperties_.addUninitializedElement();
  }

  /* ----------------------------------------------------------------------
   delete element i
  ------------------------------------------------------------------------- */

  void CustomValueTracker::deleteElement(int i)
  {
      elementProperties_.deleteElement(i);
  }

  /* ----------------------------------------------------------------------
   delete forward comm properties of element i
  ------------------------------------------------------------------------- */

  void CustomValueTracker::deleteForwardElement(int i,bool scale,bool translate,bool rotate)
  {
      elementProperties_.deleteForwardElement(i,scale,translate,rotate);
  }

  /* ----------------------------------------------------------------------
   delete restart properties of element i
  ------------------------------------------------------------------------- */

  void CustomValueTracker::deleteRestartElement(int i,bool scale,bool translate,bool rotate)
  {
      elementProperties_.deleteRestartElement(i,scale,translate,rotate);
  }

  /* ----------------------------------------------------------------------
   delete global restart properties
  ------------------------------------------------------------------------- */

  void CustomValueTracker::deleteRestartGlobal(bool scale,bool translate,bool rotate)
  {
      globalProperties_.deleteRestartGlobal(scale,translate,rotate);
      globalProperties_orig_.deleteRestartGlobal(scale,translate,rotate);
  }

  /* ----------------------------------------------------------------------
   move element i
  ------------------------------------------------------------------------- */

  void CustomValueTracker::moveElement(int i, double *delta)
  {
      elementProperties_.moveElement(i,delta);
  }

  /* ----------------------------------------------------------------------
   push / pop for list of elements
  ------------------------------------------------------------------------- */

  int CustomValueTracker::elemListBufSize(int n,int operation,bool scale,bool translate,bool rotate)
  {
    return elementProperties_.elemListBufSize(n,operation,scale,translate,rotate);
  }

  int CustomValueTracker::pushElemListToBuffer(int n, int *list, double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    return elementProperties_.pushElemListToBuffer(n,list,buf,operation,scale,translate,rotate);
  }

  int CustomValueTracker::popElemListFromBuffer(int first, int n, double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    return elementProperties_.popElemListFromBuffer(first,n,buf,operation,scale,translate,rotate);
  }

  int CustomValueTracker::pushElemListToBufferReverse(int first, int n, double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    return elementProperties_.pushElemListToBufferReverse(first,n,buf,operation,scale,translate,rotate);
  }

  int CustomValueTracker::popElemListFromBufferReverse(int n, int *list, double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    return elementProperties_.popElemListFromBufferReverse(n,list,buf,operation,scale,translate,rotate);
  }

  /* ----------------------------------------------------------------------
   push / pop for element i
  ------------------------------------------------------------------------- */

  int CustomValueTracker::elemBufSize(int operation,bool scale,bool translate,bool rotate)
  {
    
    return elementProperties_.elemBufSize(operation,scale,translate,rotate);
  }

  int CustomValueTracker::pushElemToBuffer(int i, double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    return elementProperties_.pushElemToBuffer(i,buf,operation,scale,translate,rotate);
  }

  int CustomValueTracker::popElemFromBuffer(double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    return elementProperties_.popElemFromBuffer(buf,operation,scale,translate,rotate);
  }

  /* ----------------------------------------------------------------------
   push / pop for global properties
  ------------------------------------------------------------------------- */

  int CustomValueTracker::globalPropsBufSize(int operation,bool scale,bool translate,bool rotate)
  {
    int n = 0;
    n += globalProperties_.bufSize(operation,scale,translate,rotate);
    n += globalProperties_orig_.bufSize(operation,scale,translate,rotate);
    return n;
  }

  int CustomValueTracker::pushGlobalPropsToBuffer(double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    int n = 0;
    n += globalProperties_.pushToBuffer(&(buf[n]),operation,scale,translate,rotate);
    n += globalProperties_orig_.pushToBuffer(&(buf[n]),operation,scale,translate,rotate);
    return n;
  }

  int CustomValueTracker::popGlobalPropsFromBuffer(double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    int n = 0;
    n += globalProperties_.popFromBuffer(&(buf[n]),operation,scale,translate,rotate);
    n += globalProperties_orig_.popFromBuffer(&(buf[n]),operation,scale,translate,rotate);
    return n;
  }

#endif
