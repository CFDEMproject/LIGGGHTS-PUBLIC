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
    (if not contributing author is listed, this file has been contributed
    by the core developer)

    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Philippe Seil (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef LMP_CUSTOM_VALUE_TRACKER_I_H
#define LMP_CUSTOM_VALUE_TRACKER_I_H

  /* ----------------------------------------------------------------------
   add property
  ------------------------------------------------------------------------- */

template<typename T>
T* CustomValueTracker::addElementProperty(const char *_id,
                                          const char* _comm,
                                          const char* _ref,
                                          const char *_restart,
                                          int _scalePower,
                                          int _init_len,
                                          const char *_statistics,
                                          const double _weighting_factor,
                                          ScalarContainer<double> * const _scale,
                                          ScalarContainer<double> * const _scaleAvg,
                                          const bool _enable_favre)
{
    // error if property exists already
    if(elementProperties_.getPointerById<T>(_id))
    {
        char *errmsg = new char[strlen(_id)+200];
        sprintf(errmsg,"Illegal command, features are incompatible - element property '%s' exists already",_id);
        error->all(FLERR,errmsg);
        delete []errmsg;
    }

    std::vector<std::string> id_list;
    std::string id_string(_id);

    // add property
    ContainerBase * cb_new = elementProperties_.add<T>(_id,_comm,_ref,_restart,_scalePower);
    id_list.push_back(id_string);

    // check if properties were set correctly
    // error here since ContainerBase not derived from Pointers
    if(!elementProperties_.getPointerById<T>(_id)->propertiesSetCorrectly())
    {
        char *errmsg = new char[strlen(_id)+200];
        sprintf(errmsg,"Illegal element property, comm or frame property not set correctly for property '%s'",_id);
        error->all(FLERR,errmsg);
        delete []errmsg;
    }

    // add to statistics if applicable
    if(_statistics)
    {
        if(strstr(_statistics,ContainerBase::AVERAGESUFFIX))
        {
            std::string id_string_ave = id_string;
            id_string_ave.append(ContainerBase::AVERAGESUFFIX);
            T* cb_average = elementProperties_.add<T>(id_string_ave.c_str(),_comm,_ref,_restart,_scalePower);
            cb_average->setContainerStatistics(_weighting_factor, cb_new, _scale, _scaleAvg, _enable_favre);
            id_list.push_back(id_string_ave);

            if(strstr(_statistics,"avgVar"))
            {
                // this one uses the average as reference field
                // TODO: Hard-coded higher weighting factor for second stage statistics
                std::string id_string_avg_avg = id_string_ave;
                id_string_avg_avg.append(ContainerBase::AVERAGESUFFIX);
                elementProperties_.add<T>(id_string_avg_avg.c_str(),_comm,_ref,_restart,_scalePower)->setContainerStatistics(5*_weighting_factor, cb_average,0,0,_enable_favre );
                id_list.push_back(id_string_avg_avg);

                std::string id_string_avg_mean_square = id_string_ave;
                id_string_avg_mean_square.append(ContainerBase::MEANSQUARESUFFIX);
                elementProperties_.add<T>(id_string_avg_mean_square.c_str(),_comm,_ref,_restart,_scalePower)->setContainerStatistics(5*_weighting_factor, cb_average,0,0,_enable_favre );
                id_list.push_back(id_string_avg_mean_square);
            }
        }
        if(strstr(_statistics,ContainerBase::MEANSQUARESUFFIX))
        {
            std::string id_string_mean_square = id_string;
            id_string_mean_square.append(ContainerBase::MEANSQUARESUFFIX);
            elementProperties_.add<T>(id_string_mean_square.c_str(),_comm,_ref,_restart,_scalePower)->setContainerStatistics(_weighting_factor, cb_new, _scale, _scaleAvg, _enable_favre );
            id_list.push_back(id_string_mean_square);
        }
    }

    // allocate memory and initialize
    
    for(size_t inew = 0; inew < id_list.size(); inew++)
    {
        T * const idPtr = elementProperties_.getPointerById<T>(id_list[inew].c_str());
        if(ownerMesh_)
            idPtr->addUninitialized(ownerMesh_->sizeLocal()+ownerMesh_->sizeGhost());

        if(_init_len > 0)
            idPtr->addUninitialized(_init_len);

        idPtr->setAll(0);
    }

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
   get reference
  ------------------------------------------------------------------------- */

  template<typename T>
  T* CustomValueTracker::getElementProperty(const char *_id)
  {
     return elementProperties_.getPointerById<T>(_id);
  }

  template<typename T>
  T* CustomValueTracker::getElementProperty(int _i)
  {
     return elementProperties_.getPointerByIndex<T>(_i);
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

  template<typename T>
  T* CustomValueTracker::getAvgElementProperty(const char *_id)
  {
      std::string id_string(_id);
      id_string.append(ContainerBase::AVERAGESUFFIX);
      return getElementProperty<T>(id_string.c_str());
  }

  template<typename T>
  T* CustomValueTracker::getMeanSquareElementProperty(const char *_id)
  {
      std::string id_string(_id);
      id_string.append(ContainerBase::MEANSQUARESUFFIX);
      return getElementProperty<T>(id_string.c_str());
  }

  template<typename T>
  T* CustomValueTracker::getAvgAvgElementProperty(const char *_id)
  {
      std::string id_string(_id);
      id_string.append(ContainerBase::AVERAGESUFFIX).append(ContainerBase::AVERAGESUFFIX);
      return getElementProperty<T>(id_string.c_str());
  }

  template<typename T>
  T* CustomValueTracker::getAvgMeanSquareElementProperty(const char *_id)
  {
      std::string id_string(_id);
      id_string.append(ContainerBase::AVERAGESUFFIX).append(ContainerBase::MEANSQUARESUFFIX);
      return getElementProperty<T>(id_string.c_str());
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
   delete all elements
  ------------------------------------------------------------------------- */

  void CustomValueTracker::deleteAllElements()
  {
      elementProperties_.deleteAllElements();
  }

  /* ----------------------------------------------------------------------
   delete all elements
  ------------------------------------------------------------------------- */

  void CustomValueTracker::deleteRestart(bool scale,bool translate,bool rotate)
  {
      elementProperties_.deleteRestart(scale,translate,rotate);
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

  void CustomValueTracker::moveElement(const int i, const double * const delta)
  {
      elementProperties_.moveElement(i,delta);
  }

  /* ----------------------------------------------------------------------
   push / pop for all lements
  ------------------------------------------------------------------------- */

  int CustomValueTracker::allElemBufSize(int operation,bool scale,bool translate,bool rotate) const
  {
    return elementProperties_.bufSize(operation,scale,translate,rotate);
  }

  int CustomValueTracker::pushAllElemToBuffer(double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    return elementProperties_.pushToBuffer(buf,operation,scale,translate,rotate);
  }

  int CustomValueTracker::popAllElemFromBuffer(double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    return elementProperties_.popFromBuffer(buf, operation,scale,translate,rotate);
  }

  /* ----------------------------------------------------------------------
   push / pop for list of elements
  ------------------------------------------------------------------------- */

  int CustomValueTracker::elemListBufSize(int n,int operation,bool scale,bool translate,bool rotate)
  {
    return elementProperties_.elemListBufSize(n,operation,scale,translate,rotate);
  }

  int CustomValueTracker::pushElemListToBuffer(int n, int *list, int *wraplist, double *buf, int operation, std::list<std::string> * properties, double *dlo, double *dhi, bool scale,bool translate, bool rotate)
  {
    return elementProperties_.pushElemListToBuffer(n,list, wraplist, buf,operation, properties, dlo, dhi, scale,translate,rotate);
  }

  int CustomValueTracker::popElemListFromBuffer(int first, int n, double *buf, int operation, std::list<std::string> * properties, bool scale,bool translate, bool rotate)
  {
    return elementProperties_.popElemListFromBuffer(first,n,buf,operation, properties, scale,translate,rotate);
  }

  int CustomValueTracker::pushElemListToBufferReverse(int first, int n, double *buf, int operation, std::list<std::string> * properties, bool scale,bool translate, bool rotate)
  {
    return elementProperties_.pushElemListToBufferReverse(first,n,buf,operation, properties, scale,translate,rotate);
  }

  int CustomValueTracker::popElemListFromBufferReverse(int n, int *list, double *buf, int operation, std::list<std::string> * properties, bool scale,bool translate, bool rotate)
  {
    return elementProperties_.popElemListFromBufferReverse(n,list,buf,operation, properties, scale,translate,rotate);
  }

  /* ----------------------------------------------------------------------
   push / pop for element i
  ------------------------------------------------------------------------- */

  int CustomValueTracker::elemBufSize(int operation, std::list<std::string> * properties, bool scale,bool translate,bool rotate)
  {
    
    return elementProperties_.elemBufSize(operation, properties, scale,translate,rotate);
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
