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

#ifndef LMP_ASSOCIATIVE_POINTER_ARRAY_I_H
#define LMP_ASSOCIATIVE_POINTER_ARRAY_I_H

  /* ----------------------------------------------------------------------
   constructors, destructor
  ------------------------------------------------------------------------- */

  template<typename T>
  AssociativePointerArray<T>::AssociativePointerArray()
   : content_(0), numElem_(0), maxElem_(1)
  {
    content_ = new T*[1];
    content_[0] = 0;
  }

  template<typename T>
  AssociativePointerArray<T>::~AssociativePointerArray()
  {
    for(int i = 0; i < numElem_; i++)
      delete content_[i];

    delete[] content_;
  }

  /* ----------------------------------------------------------------------
   add for per-element and pre-mesh properties
  ------------------------------------------------------------------------- */

  template<typename T> template<typename U>
  U* AssociativePointerArray<T>::add(char *_id, char* _comm, char* _ref,char *_restart,int _scalePower)
  {
    if(numElem_ == maxElem_)
      growArrays();

    content_[numElem_] = static_cast<T*>(new U(_id,_comm,_ref,_restart,_scalePower));

    numElem_++;

    //for(int i=0;i<numElem_;i++)
    //  printf("%d %s\n",i,id_[i]);
    return static_cast<U*>(content_[numElem_-1]);
  }

  /* ----------------------------------------------------------------------
   delete properties
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::remove(char *_id)
  {
    int index = idToIndex(_id);
    if(index == -1) return;

    numElem_--;

    delete content_[index];

    if(numElem_ > 0)
        content_[index] = content_[numElem_];
  }

  /* ----------------------------------------------------------------------
   get pointer to property
  ------------------------------------------------------------------------- */

  template<typename T> template<typename U>
  U* AssociativePointerArray<T>::getPointerById(char *_id)
  {
    int ind = idToIndex(_id);
    if(ind != -1)
      return dynamic_cast<U*>(content_[ind]);
    else
      return 0;
  }

  template<typename T> template<typename U>
  U* AssociativePointerArray<T>::getPointerByIndex(int i)
  {
    if(i >= size() || i < 0) return 0;
    else return dynamic_cast<U>(content_[i]);
  }

  template<typename T>
  T* AssociativePointerArray<T>::getBasePointerByIndex(int i)
  {
    if(i >= size() || i < 0) return 0;
    else return content_[i];
  }

  /* ----------------------------------------------------------------------
   memory management
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::growArrays()
  {

    // for(int i=0;i<numElem_+1;i++)
    //  printf("%d %s %d\n",i,id_[i], strcmp(id_[i],"v"));

    T *tmp[maxElem_];

    for(int i = 0; i < maxElem_; i++)
        tmp[i] = content_[i];

    delete[] content_;

    maxElem_++;
    content_ = new T*[maxElem_];

    for(int i = 0; i < numElem_; i++)
        content_[i] = tmp[i];

    //for(int i=0;i<numElem_+1;i++)
    //  printf("%d %s %d\n",i,id_[i], strcmp(id_[i],"v"));
  }

  template<typename T>
  void AssociativePointerArray<T>::grow(int to)
   {
      int by;
      for(int i = 0; i < maxElem_; i++)
      {
          by = to - getBasePointerByIndex(i)->size();
          if(by > 0)
            getBasePointerByIndex(i)->addUninitialized(by);
      }
  }

  template<typename T>
  int AssociativePointerArray<T>::size()
  {
    return numElem_;
  }

  /* ----------------------------------------------------------------------
   copy data from element from to element to
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::copyElement(int from, int to)
  {
      for(int i=0;i<numElem_;i++)
        content_[i]->copy(from,to);
  }

  /* ----------------------------------------------------------------------
   delete element n
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::deleteElement(int n)
  {
      for(int i=0;i<numElem_;i++)
        content_[i]->del(n);
  }

  /* ----------------------------------------------------------------------
   delete forward properties of element i
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::deleteForwardElement(int n,bool scale,bool translate,bool rotate)
  {
      for(int i=0;i<numElem_;i++)
        content_[i]->delForward(n,scale,translate,rotate);
  }

  /* ----------------------------------------------------------------------
   delete restart properties of element i
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::deleteRestartElement(int n,bool scale,bool translate,bool rotate)
  {
      for(int i=0;i<numElem_;i++)
        content_[i]->delRestart(n,scale,translate,rotate);
  }

  /* ----------------------------------------------------------------------
   clear reverse properties, i.e. reset all of them to 0
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::clearReverse(bool scale,bool translate,bool rotate)
  {
      for(int i=0;i<numElem_;i++)
        content_[i]->clearReverse(scale,translate,rotate);
  }

  /* ----------------------------------------------------------------------
   id 2 index
  ------------------------------------------------------------------------- */

  template<typename T>
  int AssociativePointerArray<T>::idToIndex(char *_id)
  {
    for(int i=0;i<numElem_;i++)
      if(content_[i]->matches_id(_id))
        return i;
    return -1;
  }

  template<typename T>
  void AssociativePointerArray<T>::indexToId(int index, char *_id)
  {
      content_[index]->id(_id);
  }

  /* ----------------------------------------------------------------------
   move, rotate scale all properties
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::rotate(double *dQ)
  {
      for(int i=0;i<numElem_;i++)
        content_[i]->rotate(dQ);
  }

  template<typename T>
  void AssociativePointerArray<T>::scale(double factor)
  {
      for(int i=0;i<numElem_;i++)
        content_[i]->scale(factor);
  }

  template<typename T>
  void AssociativePointerArray<T>::move(double *delta)
  {
      for(int i=0;i<numElem_;i++)
        content_[i]->move(delta);
  }

  template<typename T>
  void AssociativePointerArray<T>::moveElement(int n,double *delta)
  {
      for(int i=0;i<numElem_;i++)
        content_[i]->moveElement(n,delta);
  }

  /* ----------------------------------------------------------------------
   buf size, push, pop for all elements
  ------------------------------------------------------------------------- */

  template<typename T>
  int AssociativePointerArray<T>::bufSize(int operation,bool scale,bool translate,bool rotate)
  {
    int buf_size = 0;
    for(int i=0;i<numElem_;i++)
      buf_size += getBasePointerByIndex(i)->bufSize(operation,scale,translate,rotate);
    return buf_size;
  }

  template<typename T>
  int AssociativePointerArray<T>::pushToBuffer(double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    int nsend = 0;
    for(int i=0;i<numElem_;i++)
      nsend += getBasePointerByIndex(i)->pushToBuffer(&buf[nsend],operation,scale,translate,rotate);
    return nsend;
  }

  template<typename T>
  int AssociativePointerArray<T>::popFromBuffer(double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    int nrecv = 0;
    for(int i=0;i<numElem_;i++)
      nrecv += getBasePointerByIndex(i)->popFromBuffer(&buf[nrecv],operation,scale,translate,rotate);
    return nrecv;
  }

  /* ----------------------------------------------------------------------
   buf size, push, pop for list of elements
  ------------------------------------------------------------------------- */

  template<typename T>
  int AssociativePointerArray<T>::elemListBufSize(int n,int operation,bool scale,bool translate,bool rotate)
  {
    int buf_size = 0;
    for(int i=0;i<numElem_;i++)
      buf_size += getBasePointerByIndex(i)->elemListBufSize(n,operation,scale,translate,rotate);
    return buf_size;
  }

  template<typename T>
  int AssociativePointerArray<T>::pushElemListToBuffer(int n, int *list, double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    int nsend = 0;
    for(int i=0;i<numElem_;i++)
      nsend += getBasePointerByIndex(i)->pushElemListToBuffer(n,list,&buf[nsend],operation,scale,translate,rotate);
    return nsend;
  }

  template<typename T>
  int AssociativePointerArray<T>::popElemListFromBuffer(int first, int n, double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    int nrecv = 0;
    for(int i=0;i<numElem_;i++)
      nrecv += getBasePointerByIndex(i)->popElemListFromBuffer(first,n,&buf[nrecv],operation,scale,translate,rotate);
    return nrecv;
  }

  template<typename T>
  int AssociativePointerArray<T>::pushElemListToBufferReverse(int first, int n, double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    int nrecv = 0;
    for(int i=0;i<numElem_;i++)
      nrecv += getBasePointerByIndex(i)->pushElemListToBufferReverse(first,n,&buf[nrecv],operation,scale,translate,rotate);
    return nrecv;
  }

  template<typename T>
  int AssociativePointerArray<T>::popElemListFromBufferReverse(int n, int *list, double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    int nsend = 0;
    for(int i=0;i<numElem_;i++)
      nsend += getBasePointerByIndex(i)->popElemListFromBufferReverse(n,list,&buf[nsend],operation,scale,translate,rotate);
    return nsend;
  }

  /* ----------------------------------------------------------------------
   buf size, push, pop for single element
  ------------------------------------------------------------------------- */

  template<typename T>
  int AssociativePointerArray<T>::elemBufSize(int operation,bool scale,bool translate,bool rotate)
  {
    int buf_size = 0;
    for(int i=0;i<numElem_;i++)
      buf_size += getBasePointerByIndex(i)->elemBufSize(operation,scale,translate,rotate);
    return buf_size;
  }

  template<typename T>
  int AssociativePointerArray<T>::pushElemToBuffer(int n, double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    int nsend = 0;
    for(int i=0;i<numElem_;i++)
      nsend += getBasePointerByIndex(i)->pushElemToBuffer(n,&buf[nsend],operation,scale,translate,rotate);
    return nsend;
  }

  template<typename T>
  int AssociativePointerArray<T>::popElemFromBuffer(double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    int nrecv = 0;
    for(int i=0;i<numElem_;i++)
      nrecv += getBasePointerByIndex(i)->popElemFromBuffer(&buf[nrecv],operation,scale,translate,rotate);
    return nrecv;
  }

#endif
