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

#ifndef LMP_GENERAL_CONTAINER_I_H
#define LMP_GENERAL_CONTAINER_I_H

/* ----------------------------------------------------------------------
constructors
------------------------------------------------------------------------- */

template<typename T, int NUM_VEC, int LEN_VEC>
GeneralContainer<T,NUM_VEC,LEN_VEC>::GeneralContainer(const char *_id)
: ContainerBase(_id),
numElem_(0),
maxElem_(GROW_CONTAINER()),
defaultValue_(0)
{
    create<T>(arr_,GROW_CONTAINER(),NUM_VEC,LEN_VEC);
}

template<typename T, int NUM_VEC, int LEN_VEC>
GeneralContainer<T,NUM_VEC,LEN_VEC>::GeneralContainer(const char *_id, const char *_comm, const char *_ref, const char *_restart, int _scalePower)
: ContainerBase(_id, _comm, _ref, _restart, _scalePower),
numElem_(0),
maxElem_(GROW_CONTAINER()),
defaultValue_(0)
{
    create<T>(arr_,GROW_CONTAINER(),NUM_VEC,LEN_VEC);
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
bool GeneralContainer<T,NUM_VEC,LEN_VEC>::subtract(GeneralContainer<T,NUM_VEC,LEN_VEC> const &A,
                                                   GeneralContainer<T,NUM_VEC,LEN_VEC> const &minusB)
{
    int len = size();
    int lenA = A.size();
    int lenB = minusB.size();

    if(lenA != lenB)
        return false;

    if(len < lenA)
        addUninitialized(lenA-len);

    for(int i = 0; i < len; i++)
        for(int j = 0; j < NUM_VEC; j++)
            for(int k = 0; k < LEN_VEC; k++)
                arr_[i][j][k] = A(i)[j][k] - minusB(i)[j][k];

    return true;
}

/* ----------------------------------------------------------------------
add element(s)
------------------------------------------------------------------------- */

template<typename T, int NUM_VEC, int LEN_VEC>
void GeneralContainer<T,NUM_VEC,LEN_VEC>::add(T** elem)
{
    if(numElem_ == maxElem_)
    {
        grow<T>(arr_,maxElem_+GROW_CONTAINER(),NUM_VEC,LEN_VEC);
        maxElem_ += GROW_CONTAINER();
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
              grow<T>(arr_,maxElem_+GROW_CONTAINER(),NUM_VEC,LEN_VEC);
              maxElem_ += GROW_CONTAINER();
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
        T init_val = static_cast<T>(0);
        grow(arr_,numElem_+GROW_CONTAINER(),NUM_VEC,LEN_VEC);
        for(int i = numElem_; i < numElem_+GROW_CONTAINER(); i++)
            for(int j=0;j<NUM_VEC;j++)
                for(int k=0;k<LEN_VEC;k++)
                      arr_[i][j][k] = init_val;
        maxElem_ = numElem_ + GROW_CONTAINER();
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

/* ----------------------------------------------------------------------
average from other container
------------------------------------------------------------------------- */

template<typename T, int NUM_VEC, int LEN_VEC>
bool GeneralContainer<T,NUM_VEC,LEN_VEC>::calcAvgFromContainer()
{
    
    GeneralContainer<T,NUM_VEC,LEN_VEC> *gcont = static_cast<GeneralContainer<T,NUM_VEC,LEN_VEC>* >(container_statistics_raw_data_);
    GeneralContainer<T,1,1> *gscale = dynamic_cast<GeneralContainer<T,1,1>* >(container_statistics_scale_data_);
    GeneralContainer<T,1,1> *gscaleAvg = dynamic_cast<GeneralContainer<T,1,1>* >(container_statistics_scale_average_data_);

    // source has to be defined
    if (!gcont)
        return false;

    // only use if identical dimensions
    if(size() != gcont->size() || nVec() != gcont->nVec() || lenVec() != gcont->lenVec())
        return false;

    const int len = size();

    T epsilon = std::numeric_limits<T>::epsilon();

    if (enable_favre_)
    {
        for(int n = 0; n < len; n++)
        {
            const double scale = (*gscaleAvg)(n)[0][0] < epsilon ? 0.0 : (*gscale)(n)[0][0]/(*gscaleAvg)(n)[0][0];
            for(int i=0;i<NUM_VEC;i++)
                for(int j=0;j<LEN_VEC;j++)
                {
                    const T contribution = gcont->arr_[n][i][j];

                    if(MathExtraLiggghts::abs(arr_[n][i][j]) < epsilon)
                        arr_[n][i][j] = contribution;
                    else
                        arr_[n][i][j] = (1.-weighting_factor_*scale)*arr_[n][i][j] +
                                        weighting_factor_*scale*contribution;
                }
        }
    }
    else
    {
        for(int n = 0; n < len; n++)
            for(int i=0;i<NUM_VEC;i++)
                for(int j=0;j<LEN_VEC;j++)
                {
                    const T contribution = gcont->arr_[n][i][j];

                    if(MathExtraLiggghts::abs(arr_[n][i][j]) < epsilon)
                        arr_[n][i][j] = contribution;
                    else
                        arr_[n][i][j] = (1.-weighting_factor_)*arr_[n][i][j] +
                                        weighting_factor_*contribution;
                }
    }
    return true;
}

/* ----------------------------------------------------------------------
mean square from other container
------------------------------------------------------------------------- */

template<typename T, int NUM_VEC, int LEN_VEC>
bool GeneralContainer<T,NUM_VEC,LEN_VEC>::calcMeanSquareFromContainer()
{
    
    GeneralContainer<T,NUM_VEC,LEN_VEC> *gcont = static_cast<GeneralContainer<T,NUM_VEC,LEN_VEC>* >(container_statistics_raw_data_);
    GeneralContainer<T,1,1> *gscale = dynamic_cast<GeneralContainer<T,1,1>* >(container_statistics_scale_data_);
    GeneralContainer<T,1,1> *gscaleAvg = dynamic_cast<GeneralContainer<T,1,1>* >(container_statistics_scale_average_data_);

    // at least source has to be defined
    if (!gcont)
        return false;

    // only copy if identical
    if(size() != gcont->size() || nVec() != gcont->nVec() || lenVec() != gcont->lenVec())
        return false;

    const int len = size();

    T epsilon = std::numeric_limits<T>::epsilon();

    if (enable_favre_)
    {
        for(int n = 0; n < len; n++)
        {
            const double scale = (*gscaleAvg)(n)[0][0] < epsilon ? 0.0 : (*gscale)(n)[0][0]/(*gscaleAvg)(n)[0][0];
            for(int i=0;i<NUM_VEC;i++)
                for(int j=0;j<LEN_VEC;j++)
                {
                    const T contribution = gcont->arr_[n][i][j];

                    if(MathExtraLiggghts::abs(arr_[n][i][j]) < epsilon)
                        arr_[n][i][j] = contribution*contribution;
                    else
                        arr_[n][i][j] = (1.-weighting_factor_*scale)*arr_[n][i][j] +
                                        weighting_factor_*scale*contribution*contribution;
                }
                
        }
    }
    else
    {
        for(int n = 0; n < len; n++)
            for(int i=0;i<NUM_VEC;i++)
                for(int j=0;j<LEN_VEC;j++)
                {
                    const T contribution = gcont->arr_[n][i][j];

                    if(MathExtraLiggghts::abs(arr_[n][i][j]) < epsilon)
                        arr_[n][i][j] = contribution*contribution;
                    else
                        arr_[n][i][j] = (1.-weighting_factor_)*arr_[n][i][j] +
                                        weighting_factor_*contribution*contribution;
                }
                
    }

    return true;
}

// This is the averaging for the scaling arrays
template<typename T, int NUM_VEC, int LEN_VEC>
bool GeneralContainer<T,NUM_VEC,LEN_VEC>::calcSumFromContainer()
{
    
    GeneralContainer<T,NUM_VEC,LEN_VEC> *gcont = static_cast<GeneralContainer<T,NUM_VEC,LEN_VEC>* >(container_statistics_raw_data_);

    // at least source has to be defined
    if (!gcont)
        return false;

    // only copy if identical
    if(size() != gcont->size() || nVec() != gcont->nVec() || lenVec() != gcont->lenVec())
        return false;

    const int len = size();
    for(int n = 0; n < len; n++)
        for(int i=0;i<NUM_VEC;i++)
            for(int j=0;j<LEN_VEC;j++)
            {
                arr_[n][i][j] = (1.-weighting_factor_)*arr_[n][i][j] +
                                weighting_factor_*(*gcont)(n)[i][j];
                
                if (arr_[n][i][j] < std::numeric_limits<T>::epsilon())
                    arr_[n][i][j] = 0;
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
    int len = std::min(to,size());
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
void GeneralContainer<T,NUM_VEC,LEN_VEC>::move(const double * const delta)
{
  if(isTranslationInvariant()) return;

  int len = size();

  for(int i = 0; i < len; i++)
        for(int j = 0; j < NUM_VEC; j++)
            for(int k = 0; k < LEN_VEC; k++)
                arr_[i][j][k] += delta[k];
}

template<typename T, int NUM_VEC, int LEN_VEC>
void GeneralContainer<T,NUM_VEC,LEN_VEC>::moveElement(const int i, const double * const delta)
{
  if(isTranslationInvariant()) return;

        for(int j = 0; j < NUM_VEC; j++)
            for(int k = 0; k < LEN_VEC; k++)
                arr_[i][j][k] += delta[k];
}

template<typename T, int NUM_VEC, int LEN_VEC>
void GeneralContainer<T,NUM_VEC,LEN_VEC>::rotate(const double * const dQ)
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
int GeneralContainer<T,NUM_VEC,LEN_VEC>::bufSize(int operation,bool scale,bool translate,bool rotate) const
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
int GeneralContainer<T,NUM_VEC,LEN_VEC>::pushElemListToBuffer(int n, int *list, int *wraplist, double *buf,int operation, double *dlo, double *dhi, bool scale,bool translate, bool rotate)
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
            {
                buf[m] = static_cast<double>(arr_[i][j][k]);
                if (wrapPeriodic())
                {
                    const int wrap = wraplist[ii];
                    if (wrap != IS_GHOST)
                    {
                        if      ((k == 0 && wrap == IS_GHOST_WRAP_DIM_0_NEG) ||
                                 (k == 1 && wrap == IS_GHOST_WRAP_DIM_1_NEG) ||
                                 (k == 2 && wrap == IS_GHOST_WRAP_DIM_2_NEG)   )
                            buf[m] -= dhi[k] - dlo[k];
                        else if ((k == 0 && wrap == IS_GHOST_WRAP_DIM_0_POS) ||
                                 (k == 1 && wrap == IS_GHOST_WRAP_DIM_1_POS) ||
                                 (k == 2 && wrap == IS_GHOST_WRAP_DIM_2_POS)   )
                            buf[m] += dhi[k] - dlo[k];
                    }
                }
                m++;
            }
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

    if(COMM_TYPE_REVERSE == this->communicationType())
    {
        
        for(int ii = 0; ii < n; ii++)
        {
            i = list[ii];
            for(int j = 0; j < NUM_VEC; j++)
                for(int k = 0; k < LEN_VEC; k++)
                    arr_[i][j][k] += static_cast<T>(buf[m++]);
        }
    }
    else if(sizeof(int) == sizeof(T) && COMM_TYPE_REVERSE_BITFIELD == this->communicationType())
    {
        
        for(int ii = 0; ii < n; ii++)
        {
            i = list[ii];
            for(int j = 0; j < NUM_VEC; j++)
                for(int k = 0; k < LEN_VEC; k++)
                    arr_[i][j][k] = (T) (static_cast<int>(arr_[i][j][k]) | static_cast<int>(buf[m++]));
        }
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
