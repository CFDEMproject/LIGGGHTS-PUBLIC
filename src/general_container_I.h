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

#ifndef LMP_GENERAL_CONTAINER_I_H
#define LMP_GENERAL_CONTAINER_I_H

  /* ----------------------------------------------------------------------
   constructors
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  GeneralContainer<T,NUM_VEC,LEN_VEC>::GeneralContainer()
  : ContainerBase(),
    numElem_(0),
    maxElem_(GROW)
  {
          create<T>(arr_,GROW,NUM_VEC,LEN_VEC);
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  GeneralContainer<T,NUM_VEC,LEN_VEC>::GeneralContainer(char *_id, char *_comm, char *_ref, char *_restart, int _scalePower)
  : ContainerBase(_id, _comm, _ref, _restart, _scalePower),
    numElem_(0),
    maxElem_(GROW)
  {
          create<T>(arr_,GROW,NUM_VEC,LEN_VEC);
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  GeneralContainer<T,NUM_VEC,LEN_VEC>::GeneralContainer(GeneralContainer<T,NUM_VEC,LEN_VEC> const &orig)
  : ContainerBase(orig),
    numElem_(orig.numElem_),
    maxElem_(orig.numElem_)
  {
          create<T>(arr_,maxElem_,NUM_VEC,LEN_VEC);
          for(int i=0;i<maxElem_;i++)
                  for(int ii=0;ii<NUM_VEC;ii++)
                          for(int jj=0;jj<LEN_VEC;jj++)
                                  arr_[i][ii][jj] = orig.arr_[i][ii][jj];
  }

  /* ----------------------------------------------------------------------
   destructor
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  GeneralContainer<T,NUM_VEC,LEN_VEC>::~GeneralContainer()
  {
          destroy<T>(arr_);
  }

  /* ----------------------------------------------------------------------
   add element(s)
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  void GeneralContainer<T,NUM_VEC,LEN_VEC>::add(T** elem)
  {
          if(numElem_ == maxElem_)
          {
                  grow<T>(arr_,maxElem_+GROW,NUM_VEC,LEN_VEC);
                  maxElem_ += GROW;
          }
          for(int i=0;i<NUM_VEC;i++)
                  for(int j=0;j<LEN_VEC;j++)
                          arr_[numElem_][i][j] = elem[i][j];
          numElem_++;
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  void GeneralContainer<T,NUM_VEC,LEN_VEC>::addUninitialized(int n)
  {
        numElem_ += n;
        if(numElem_ >= maxElem_)
        {
            grow(arr_,numElem_+GROW,NUM_VEC,LEN_VEC);
            maxElem_ = numElem_ + GROW;
        }
  }

  /* ----------------------------------------------------------------------
   delete an element
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  void GeneralContainer<T,NUM_VEC,LEN_VEC>::del(int n)
  {
          numElem_--;
          if(numElem_ == n) return;
          for(int i=0;i<NUM_VEC;i++)
                  for(int j=0;j<LEN_VEC;j++)
                          arr_[n][i][j] = arr_[numElem_][i][j];
  }

  /* ----------------------------------------------------------------------
   copy element data
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  void GeneralContainer<T,NUM_VEC,LEN_VEC>::copy(int from,int to)
  {
          for(int i=0;i<NUM_VEC;i++)
                  for(int j=0;j<LEN_VEC;j++)
                          arr_[to][i][j] = arr_[from][i][j];
  }

  /* ----------------------------------------------------------------------
   delete an element
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  void GeneralContainer<T,NUM_VEC,LEN_VEC>::delForward(int n,bool scale,bool translate,bool rotate)
  {
          // do only delete property if it is a forward comm property
          if(!decideBufferOperation(OPERATION_COMM_FORWARD, scale, translate, rotate))
            return;

          numElem_--;
          if(numElem_ == n) return;
          for(int i=0;i<NUM_VEC;i++)
                  for(int j=0;j<LEN_VEC;j++)
                          arr_[n][i][j] = arr_[numElem_][i][j];
  }

  /* ----------------------------------------------------------------------
   delete an element
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  void GeneralContainer<T,NUM_VEC,LEN_VEC>::delRestart(int n,bool scale,bool translate,bool rotate)
  {
          // do only delete property if it is a forward comm property
          if(!decideBufferOperation(OPERATION_RESTART, scale, translate, rotate))
            return;

          numElem_--;
          if(numElem_ == n) return;
          for(int i=0;i<NUM_VEC;i++)
                  for(int j=0;j<LEN_VEC;j++)
                          arr_[n][i][j] = arr_[numElem_][i][j];
  }

  /* ----------------------------------------------------------------------
   get an element
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  void GeneralContainer<T,NUM_VEC,LEN_VEC>::get(int n, T** elem)
  {
          for(int i=0;i<NUM_VEC;i++)
                  for(int j=0;j<LEN_VEC;j++)
                          elem[i][j] = arr_[n][i][j];
  }

  /* ----------------------------------------------------------------------
   operator()
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  T**& GeneralContainer<T,NUM_VEC,LEN_VEC>::operator() (int n)
  {
          return arr_[n];
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  T** const& GeneralContainer<T,NUM_VEC,LEN_VEC>::operator() (int n) const
  {
          return arr_[n];
  }

  /* ---------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  void GeneralContainer<T,NUM_VEC,LEN_VEC>::set(int n, T** elem)
  {
          for(int i=0;i<NUM_VEC;i++)
                          for(int j=0;j<LEN_VEC;j++)
                                  arr_[n][i][j] = elem[i][j];
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  void GeneralContainer<T,NUM_VEC,LEN_VEC>::setAll(T def)
  {
      for(int n=0;n<size();n++)
          for(int i=0;i<NUM_VEC;i++)
                          for(int j=0;j<LEN_VEC;j++)
                                  arr_[n][i][j] = def;
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  T*** GeneralContainer<T,NUM_VEC,LEN_VEC>::begin()
  {
          return arr_;
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  int GeneralContainer<T,NUM_VEC,LEN_VEC>::getElemSize()
  {
          return NUM_VEC*LEN_VEC*sizeof(T);
  }

  /* ----------------------------------------------------------------------
   translate, rotate, scale
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  void GeneralContainer<T,NUM_VEC,LEN_VEC>::scale(double factor)
  {
      if(isScaleInvariant()) return;

      double factorApplied = 1.;
      for(int i = 0; i < scalePower_; i++)
        factorApplied *= factor;

      for(int i=0;i<size();i++)
            for(int j=0;j<NUM_VEC;j++)
                for(int k=0;k<LEN_VEC;k++)
                    arr_[i][j][k] *= factorApplied;
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  void GeneralContainer<T,NUM_VEC,LEN_VEC>::move(double *delta)
  {
      if(isTranslationInvariant()) return;

      for(int i=0;i<size();i++)
            for(int j=0;j<NUM_VEC;j++)
                for(int k=0;k<LEN_VEC;k++)
                    arr_[i][j][k] += delta[k];
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  void GeneralContainer<T,NUM_VEC,LEN_VEC>::moveElement(int i,double *delta)
  {
      if(isTranslationInvariant()) return;

            for(int j=0;j<NUM_VEC;j++)
                for(int k=0;k<LEN_VEC;k++)
                    arr_[i][j][k] += delta[k];
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  void GeneralContainer<T,NUM_VEC,LEN_VEC>::rotate(double *dQ)
  {
      if(isRotationInvariant()) return;

      // ATTENTION: only correct for 3D vectors
      for(int i=0;i<size();i++)
            for(int j=0;j<NUM_VEC;j++)
              MathExtraLiggghts::vec_quat_rotate(arr_[i][j],dQ);
  }

  /* ----------------------------------------------------------------------
   buffer size for all elements, push / pop for all elements
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  int GeneralContainer<T,NUM_VEC,LEN_VEC>::bufSize(int operation,bool scale,bool translate,bool rotate)
  {
      if(!this->decideBufferOperation(operation,scale,translate,rotate))
            return 0;

      return (1 + size()*NUM_VEC*LEN_VEC);
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  int GeneralContainer<T,NUM_VEC,LEN_VEC>::pushToBuffer(double *buf,int operation,bool scale,bool translate, bool rotate)
  {
          //TODO throw error if sizeof(T) > sizeof(double)

          int m = 0;

          if(!this->decideBufferOperation(operation,scale,translate,rotate))
            return 0;

          buf[m++] = static_cast<double>(size());

          for(int i=0;i<size();i++)
            for(int j=0;j<NUM_VEC;j++)
                for(int k=0;k<LEN_VEC;k++)
                    buf[m++] = static_cast<double>(arr_[i][j][k]);

          return (1 + size()*NUM_VEC*LEN_VEC);
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  int GeneralContainer<T,NUM_VEC,LEN_VEC>::popFromBuffer(double *buf,int operation,bool scale,bool translate, bool rotate)
  {
          int nNew, m = 0;

          if(!this->decideBufferOperation(operation,scale,translate,rotate))
            return 0;

          T** tmp;
          create<T>(tmp,NUM_VEC,LEN_VEC);

          nNew = static_cast<int>(buf[m++]);

          for(int i=0;i<nNew;i++)
          {
            for(int j=0;j<NUM_VEC;j++)
                for(int k=0;k<LEN_VEC;k++)
                    tmp[j][k] = static_cast<T>(buf[m++]);
            add(tmp);
          }

          destroy<T>(tmp);

          return (1 + nNew*NUM_VEC*LEN_VEC);
  }

  /* ----------------------------------------------------------------------
   buffer size for a list of elements, push / pop a list of elements
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  int GeneralContainer<T,NUM_VEC,LEN_VEC>::elemListBufSize(int n,int operation,bool scale,bool translate,bool rotate)
  {
      if(!this->decideBufferOperation(operation,scale,translate,rotate))
            return 0;

      return (n*NUM_VEC*LEN_VEC);
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  int GeneralContainer<T,NUM_VEC,LEN_VEC>::pushElemListToBuffer(int n, int *list,double *buf,int operation,bool scale,bool translate, bool rotate)
  {
        int i,m = 0;

        if(!this->decideBufferOperation(operation,scale,translate,rotate))
            return 0;

        for(int ii=0;ii<n;ii++)
        {
            i = list[ii];
            for(int j=0;j<NUM_VEC;j++)
                for(int k=0;k<LEN_VEC;k++)
                    buf[m++] = static_cast<double>(arr_[i][j][k]);
        }

        return (n*NUM_VEC*LEN_VEC);
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  int GeneralContainer<T,NUM_VEC,LEN_VEC>::popElemListFromBuffer(int first, int n, double *buf,int operation,bool scale,bool translate, bool rotate)
  {
        int m = 0;

        if(!this->decideBufferOperation(operation,scale,translate,rotate))
            return 0;

        T** tmp;
        create<T>(tmp,NUM_VEC,LEN_VEC);

        for(int i=first;i<first+n;i++)
        {
            for(int j=0;j<NUM_VEC;j++)
                for(int k=0;k<LEN_VEC;k++)
                    tmp[j][k] = static_cast<T>(buf[m++]);

            add(tmp);
        }

        destroy<T>(tmp);

        return (n*NUM_VEC*LEN_VEC);
  }

  /* ----------------------------------------------------------------------
   buffer size for a single element, push / pop a single element
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  int GeneralContainer<T,NUM_VEC,LEN_VEC>::elemBufSize(int operation,bool scale,bool translate,bool rotate)
  {
      if(!this->decideBufferOperation(operation,scale,translate,rotate))
            return 0;

      return (NUM_VEC*LEN_VEC);
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  int GeneralContainer<T,NUM_VEC,LEN_VEC>::pushElemToBuffer(int i, double *buf,int operation,bool scale,bool translate, bool rotate)
  {
        int m = 0;

        if(!this->decideBufferOperation(operation,scale,translate,rotate))
            return 0;

        for(int j=0;j<NUM_VEC;j++)
            for(int k=0;k<LEN_VEC;k++)
                buf[m++] = static_cast<double>(arr_[i][j][k]);

        return (NUM_VEC*LEN_VEC);
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  int GeneralContainer<T,NUM_VEC,LEN_VEC>::popElemFromBuffer(double *buf,int operation,bool scale,bool translate, bool rotate)
  {
        int m = 0;

        if(!this->decideBufferOperation(operation,scale,translate,rotate))
            return 0;

        T** tmp;
        create<T>(tmp,NUM_VEC,LEN_VEC);

        for(int j=0;j<NUM_VEC;j++)
            for(int k=0;k<LEN_VEC;k++)
                tmp[j][k] = static_cast<T>(buf[m++]);
        add(tmp);

        destroy<T>(tmp);

        return (NUM_VEC*LEN_VEC);
  }

#endif
