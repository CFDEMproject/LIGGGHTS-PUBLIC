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

#ifndef LMP_GENERAL_CONTAINER_I_H
#define LMP_GENERAL_CONTAINER_I_H

  /* ----------------------------------------------------------------------
   constructors
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  GeneralContainer<T,NUM_VEC,LEN_VEC>::GeneralContainer(const char *_id)
  : ContainerBase(_id),
    numElem_(0),
    maxElem_(GROW),
    defaultValue_(0)
  {
          create<T>(arr_,GROW,NUM_VEC,LEN_VEC);
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  GeneralContainer<T,NUM_VEC,LEN_VEC>::GeneralContainer(const char *_id, const char *_comm, const char *_ref, const char *_restart, int _scalePower)
  : ContainerBase(_id, _comm, _ref, _restart, _scalePower),
    numElem_(0),
    maxElem_(GROW),
    defaultValue_(0)
  {
          create<T>(arr_,GROW,NUM_VEC,LEN_VEC);
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  GeneralContainer<T,NUM_VEC,LEN_VEC>::GeneralContainer(GeneralContainer<T,NUM_VEC,LEN_VEC> const &orig)
  : ContainerBase(orig),
    numElem_(orig.numElem_),
    maxElem_(orig.numElem_),
    defaultValue_(orig.defaultValue_)
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
   check if data is of type double
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  bool GeneralContainer<T,NUM_VEC,LEN_VEC>::isDoubleData()
  {
      // partial templatization does not work
      // std::is_same<T,double>::value is from C++11
      // this is work-around

      if(sizeof(T) == sizeof(double))
        return true;
      else
        return false;
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  bool GeneralContainer<T,NUM_VEC,LEN_VEC>::isIntData()
  {
      // partial templatization does not work
      // std::is_same<T,double>::value is from C++11
      // this is work-around

      if(sizeof(T) == sizeof(int))
        return true;
      else
        return false;
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
  void GeneralContainer<T,NUM_VEC,LEN_VEC>::addZero()
  {
          if(numElem_ == maxElem_)
          {
                  grow<T>(arr_,maxElem_+GROW,NUM_VEC,LEN_VEC);
                  maxElem_ += GROW;
          }
          for(int i=0;i<NUM_VEC;i++)
                  for(int j=0;j<LEN_VEC;j++)
                          arr_[numElem_][i][j] = static_cast<T>(0);
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
          if(!decidePackUnpackOperation(OPERATION_COMM_FORWARD, scale, translate, rotate))
            return;

          numElem_--;
          if(numElem_ == n) return;
          for(int i=0;i<NUM_VEC;i++)
                  for(int j=0;j<LEN_VEC;j++)
                          arr_[n][i][j] = arr_[numElem_][i][j];
  }

  /* ----------------------------------------------------------------------
   clear reverse properties, i.e. reset all of them to 0
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  void GeneralContainer<T,NUM_VEC,LEN_VEC>::clearReverse(bool scale,bool translate,bool rotate)
  {
      // do only reset property if it is a reverse comm property
      if(!decidePackUnpackOperation(OPERATION_COMM_REVERSE, scale, translate, rotate))
        return;

      int len = size();
      for(int i = 0; i < len; i++)
            for(int j = 0; j < NUM_VEC; j++)
                for(int k = 0; k < LEN_VEC; k++)
                    arr_[i][j][k] = 0.;
  }

  /* ----------------------------------------------------------------------
   delete an element if restart
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  void GeneralContainer<T,NUM_VEC,LEN_VEC>::delRestart(int n,bool scale,bool translate,bool rotate)
  {
          // do only delete property if it is a restart property
          if(!decidePackUnpackOperation(OPERATION_RESTART, scale, translate, rotate))
            return;

          numElem_--;
          if(numElem_ == n) return;
          for(int i=0;i<NUM_VEC;i++)
                  for(int j=0;j<LEN_VEC;j++)
                          arr_[n][i][j] = arr_[numElem_][i][j];
  }

  /* ----------------------------------------------------------------------
   delete all elements if restart
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  void GeneralContainer<T,NUM_VEC,LEN_VEC>::delRestart(bool scale,bool translate,bool rotate)
  {
          // do only delete property if it is a restart property
          if(!decidePackUnpackOperation(OPERATION_RESTART, scale, translate, rotate))
            return;

          numElem_ = 0;
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

  /* ----------------------------------------------------------------------
   set all data by copy from other container
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  bool GeneralContainer<T,NUM_VEC,LEN_VEC>::setFromContainer(ContainerBase *cont)
  {
      GeneralContainer<T,NUM_VEC,LEN_VEC> *gcont = static_cast<GeneralContainer<T,NUM_VEC,LEN_VEC>* >(cont);

      if(size() != gcont->size() || nVec() != gcont->nVec() || lenVec() != gcont->lenVec())
        return false;

      int len = size();
      for(int n = 0; n < len; n++)
          for(int i=0;i<NUM_VEC;i++)
                  for(int j=0;j<LEN_VEC;j++)
                  {
                          arr_[n][i][j] = gcont->arr_[n][i][j];
                          
                  }

      return true;
  }

  /* ---------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  void GeneralContainer<T,NUM_VEC,LEN_VEC>::setToDefault(int n)
  {
    
          for(int i = 0; i < NUM_VEC; i++)
                          for(int j = 0; j < LEN_VEC; j++)
                                  arr_[n][i][j] = defaultValue_;
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  void GeneralContainer<T,NUM_VEC,LEN_VEC>::set(int n, T** elem)
  {
          for(int i = 0; i < NUM_VEC; i++)
                          for(int j = 0; j < LEN_VEC; j++)
                                  arr_[n][i][j] = elem[i][j];
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  void GeneralContainer<T,NUM_VEC,LEN_VEC>::set(int n, int m, T* elem)
  {
      for(int j = 0; j < LEN_VEC; j++)
          arr_[n][m][j] = elem[j];
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  void GeneralContainer<T,NUM_VEC,LEN_VEC>::setAll(T def)
  {
      int len = size();
      for(int n = 0; n < len; n++)
          for(int i = 0; i < NUM_VEC; i++)
                          for(int j = 0; j < LEN_VEC; j++)
                                  arr_[n][i][j] = def;
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  void GeneralContainer<T,NUM_VEC,LEN_VEC>::setAll(int to,T def)
  {
      int len = MathExtraLiggghts::min(to,size());
      for(int n = 0; n < len; n++)
          for(int i = 0; i < NUM_VEC; i++)
                          for(int j = 0; j < LEN_VEC; j++)
                                  arr_[n][i][j] = def;
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  T*** GeneralContainer<T,NUM_VEC,LEN_VEC>::begin()
  {
          return arr_;
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  void* GeneralContainer<T,NUM_VEC,LEN_VEC>::begin_slow_dirty()
  {
          return (void*) arr_;
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  int GeneralContainer<T,NUM_VEC,LEN_VEC>::getElemSize()
  {
          return NUM_VEC*LEN_VEC*sizeof(T);
  }

  /* ----------------------------------------------------------------------
   min,max
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  T GeneralContainer<T,NUM_VEC,LEN_VEC>::max_scalar()
  {
      T max = arr_[0][0][0];

      int len = size();
      for(int i = 0; i < len; i++)
            for(int j = 0; j < NUM_VEC; j++)
                for(int k = 0; k < LEN_VEC; k++)
                    if(arr_[i][j][k] > max)
                        max = arr_[i][j][k];

      return max;
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  T GeneralContainer<T,NUM_VEC,LEN_VEC>::min_scalar()
  {
      T min = arr_[0][0][0];

      int len = size();
      for(int i = 0; i < len; i++)
            for(int j = 0; j < NUM_VEC; j++)
                for(int k = 0; k < LEN_VEC; k++)
                    if(arr_[i][j][k] < min)
                        min = arr_[i][j][k];

      return min;
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

      int len = size();
      for(int i = 0; i < len; i++)
            for(int j = 0; j < NUM_VEC;j++)
                for(int k = 0; k < LEN_VEC; k++)
                    arr_[i][j][k] *= factorApplied;
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  void GeneralContainer<T,NUM_VEC,LEN_VEC>::move(double *delta)
  {
      if(isTranslationInvariant()) return;

      int len = size();

      for(int i = 0; i < len; i++)
            for(int j = 0; j < NUM_VEC; j++)
                for(int k = 0; k < LEN_VEC; k++)
                    arr_[i][j][k] += delta[k];
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  void GeneralContainer<T,NUM_VEC,LEN_VEC>::moveElement(int i,double *delta)
  {
      if(isTranslationInvariant()) return;

            for(int j = 0; j < NUM_VEC; j++)
                for(int k = 0; k < LEN_VEC; k++)
                    arr_[i][j][k] += delta[k];
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  void GeneralContainer<T,NUM_VEC,LEN_VEC>::rotate(double *dQ)
  {
      if(isRotationInvariant()) return;

      // ATTENTION: only correct for 3D vectors
      int len = size();
      for(int i = 0; i < len; i++)
            for(int j = 0; j < NUM_VEC; j++)
              MathExtraLiggghts::vec_quat_rotate(arr_[i][j],dQ);
  }

  /* ----------------------------------------------------------------------
   buffer size for all elements, push / pop for all elements
   used for global properties
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  int GeneralContainer<T,NUM_VEC,LEN_VEC>::bufSize(int operation,bool scale,bool translate,bool rotate)
  {
      if(!this->decidePackUnpackOperation(operation,scale,translate,rotate))
            return 0;

      if(!this->decideCommOperation(operation))
            return 0;

      return (1 + size()*NUM_VEC*LEN_VEC);
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  int GeneralContainer<T,NUM_VEC,LEN_VEC>::pushToBuffer(double *buf,int operation,bool scale,bool translate, bool rotate)
  {
          //TODO throw error if sizeof(T) > sizeof(double)

          int m = 0;

          if(!this->decidePackUnpackOperation(operation,scale,translate,rotate))
            return 0;

          int len = size();

          buf[m++] = static_cast<double>(len);

          for(int i = 0; i < len; i++)
            for(int j = 0; j < NUM_VEC; j++)
                for(int k = 0; k < LEN_VEC; k++)
                    buf[m++] = static_cast<double>(arr_[i][j][k]);

          return (1 + len*NUM_VEC*LEN_VEC);
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  int GeneralContainer<T,NUM_VEC,LEN_VEC>::popFromBuffer(double *buf,int operation,bool scale,bool translate, bool rotate)
  {
          int nNew, m = 0;

          if(!this->decidePackUnpackOperation(operation,scale,translate,rotate))
            return 0;

          if(decideCreateNewElements(operation))
          {
              T** tmp;
              create<T>(tmp,NUM_VEC,LEN_VEC);

              nNew = static_cast<int>(buf[m++]);

              for(int i = 0; i < nNew; i++)
              {
                for(int j = 0; j < NUM_VEC; j++)
                    for(int k = 0; k < LEN_VEC; k++)
                        tmp[j][k] = static_cast<T>(buf[m++]);
                add(tmp);
              }

              destroy<T>(tmp);

              return (1 + nNew*NUM_VEC*LEN_VEC);
          }
          else return 0;
  }

  /* ----------------------------------------------------------------------
   buffer size for a list of elements, push / pop a list of elements
   used for borders, fw and rev comm for element properties
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  int GeneralContainer<T,NUM_VEC,LEN_VEC>::elemListBufSize(int n,int operation,bool scale,bool translate,bool rotate)
  {
      if(!this->decidePackUnpackOperation(operation,scale,translate,rotate))
            return 0;

      if(!this->decideCommOperation(operation))
            return 0;

      return (n*NUM_VEC*LEN_VEC);
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  int GeneralContainer<T,NUM_VEC,LEN_VEC>::pushElemListToBuffer(int n, int *list,double *buf,int operation,bool scale,bool translate, bool rotate)
  {
        int i,m = 0;

        if(!this->decidePackUnpackOperation(operation,scale,translate,rotate))
            return 0;

        if(!this->decideCommOperation(operation))
            return 0;

        for(int ii = 0; ii < n; ii++)
        {
            i = list[ii];
            for(int j = 0; j < NUM_VEC; j++)
                for(int k = 0; k < LEN_VEC; k++)
                    buf[m++] = static_cast<double>(arr_[i][j][k]);
        }

        return (n*NUM_VEC*LEN_VEC);
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  int GeneralContainer<T,NUM_VEC,LEN_VEC>::popElemListFromBuffer(int first, int n, double *buf,int operation,bool scale,bool translate, bool rotate)
  {
        int m = 0;

        if(!this->decidePackUnpackOperation(operation,scale,translate,rotate))
            return 0;

        bool pullBuf = decideCommOperation(operation);

        bool createElem = decideCreateNewElements(operation);

        T** tmp;
        create<T>(tmp,NUM_VEC,LEN_VEC);

        for(int i = first; i < first+n; i++)
        {
            for(int j = 0; j < NUM_VEC; j++)
                for(int k = 0; k < LEN_VEC; k++)
                    (createElem ? tmp[j][k] : arr_[i][j][k]) = (pullBuf ? static_cast<T>(buf[m++]) : static_cast<T>(0));

            if(createElem) add(tmp);
        }

        destroy<T>(tmp);

        return m;
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  int GeneralContainer<T,NUM_VEC,LEN_VEC>::pushElemListToBufferReverse(int first, int n, double *buf,int operation,bool scale,bool translate, bool rotate)
  {
        int m = 0;

        if(!this->decidePackUnpackOperation(operation,scale,translate,rotate))
            return 0;

        for(int i = first; i < first+n; i++)
        {
            for(int j = 0; j < NUM_VEC; j++)
                for(int k = 0; k < LEN_VEC; k++)
                    buf[m++] = static_cast<double>(arr_[i][j][k]);
        }

        return (n*NUM_VEC*LEN_VEC);
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  int GeneralContainer<T,NUM_VEC,LEN_VEC>::popElemListFromBufferReverse(int n, int *list,double *buf,int operation,bool scale,bool translate, bool rotate)
  {
        int i,m = 0;

        if(!this->decidePackUnpackOperation(operation,scale,translate,rotate))
            return 0;

        for(int ii = 0; ii < n; ii++)
        {
            i = list[ii];
            for(int j = 0; j < NUM_VEC; j++)
                for(int k = 0; k < LEN_VEC; k++)
                    arr_[i][j][k] += static_cast<T>(buf[m++]);
        }

        return (n*NUM_VEC*LEN_VEC);
  }

  /* ----------------------------------------------------------------------
   buffer size for a single element, push / pop a single element
   used for exchange of single elements
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  int GeneralContainer<T,NUM_VEC,LEN_VEC>::elemBufSize(int operation,bool scale,bool translate,bool rotate)
  {
      
      if(!this->decidePackUnpackOperation(operation,scale,translate,rotate))
            return 0;

      if(!this->decideCommOperation(operation))
            return 0;
      
      return (NUM_VEC*LEN_VEC);
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  int GeneralContainer<T,NUM_VEC,LEN_VEC>::pushElemToBuffer(int i, double *buf,int operation,bool scale,bool translate, bool rotate)
  {
        int m = 0;

        if(!this->decidePackUnpackOperation(operation,scale,translate,rotate))
            return 0;

        if(!this->decideCommOperation(operation))
            return 0;

        for(int j = 0; j < NUM_VEC; j++)
            for(int k = 0; k < LEN_VEC; k++)
                buf[m++] = static_cast<double>(arr_[i][j][k]);

        return m;
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  int GeneralContainer<T,NUM_VEC,LEN_VEC>::popElemFromBuffer(double *buf,int operation,bool scale,bool translate, bool rotate)
  {
        int m = 0;

        if(!this->decidePackUnpackOperation(operation,scale,translate,rotate))
            return 0;

        bool pullBuf = decideCommOperation(operation);

        T** tmp;
        create<T>(tmp,NUM_VEC,LEN_VEC);

        for(int j = 0; j < NUM_VEC; j++)
            for(int k = 0; k < LEN_VEC; k++)
                tmp[j][k] = pullBuf ? static_cast<T>(buf[m++]) : static_cast<T>(0);

        add(tmp);
        destroy<T>(tmp);

        return m;
  }

#endif
