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

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef LMP_CONTAINER_BASE_H
#define LMP_CONTAINER_BASE_H

#include <string>
#include <list>
#include <string.h>

namespace LAMMPS_NS
{
  // buffer operation types (for push and pop)

  enum{ OPERATION_COMM_EXCHANGE,
        OPERATION_COMM_BORDERS,
        OPERATION_COMM_FORWARD,
        OPERATION_COMM_REVERSE,
        OPERATION_RESTART,
        OPERATION_UNDEFINED};

  /* ----------------------------------------------------------------------
   definition of reference frames and comm types
  ------------------------------------------------------------------------- */

  // reference frame types
  // invariant: invariant to scaling, translation, rotation
  // trans invariant: invariant to translation, not invariant to scaling, rotation
  // trans+rot invariant: invariant to translation, rotation, not invariant to scaling
  // general: not invariant to scaling, translation, rotation

  enum{ REF_FRAME_UNDEFINED,
        REF_FRAME_INVARIANT,
        REF_FRAME_SCALE_TRANS_INVARIANT,
        REF_FRAME_TRANS_ROT_INVARIANT,
        REF_FRAME_TRANS_INVARIANT,
        REF_FRAME_GENERAL};

  // communication types

  enum{ // communication invoked manually
        COMM_TYPE_MANUAL,
        // only exchange and borders comm
        COMM_EXCHANGE_BORDERS,
        // forward comm every step
        COMM_TYPE_FORWARD,
        // forward comm based on reference frame setting
        // ie if mesh rotates, egdeVecs are communicated
        
        COMM_TYPE_FORWARD_FROM_FRAME,
        // reverse comm every step
        
        COMM_TYPE_REVERSE,
        // reverse comm every step
        
        COMM_TYPE_REVERSE_BITFIELD,
        // no comm at all
        
        COMM_TYPE_NONE,
        // undefined state for error check
        COMM_TYPE_UNDEFINED};  // communication types

  // restart types

  enum{ RESTART_TYPE_UNDEFINED,
        RESTART_TYPE_YES,
        RESTART_TYPE_NO};

    // return status of checkBorderElement{Left,Right} in multi_node_mesh_parallel
    // also enables container properties to be wrapped around periodic boundaries
    // e.g. node positions of triangles
    enum
    {
        NOT_GHOST               = 0,
        IS_GHOST                = 1<<0,
        IS_GHOST_WRAP_DIM_0_NEG = 1<<1,
        IS_GHOST_WRAP_DIM_0_POS = 1<<2,
        IS_GHOST_WRAP_DIM_1_NEG = 1<<3,
        IS_GHOST_WRAP_DIM_1_POS = 1<<4,
        IS_GHOST_WRAP_DIM_2_NEG = 1<<5,
        IS_GHOST_WRAP_DIM_2_POS = 1<<6
    };

  /* ----------------------------------------------------------------------
   class definitions
  ------------------------------------------------------------------------- */

  class ContainerBase
  {
      public:

          ContainerBase(const char *_id);

          virtual ~ContainerBase();

          void setProperties(const char *_id, const char* _comm, const char* _ref, const char *_restart,int _scalePower = 1);
          bool propertiesSetCorrectly();

          void setContainerStatistics(const double _weighting_factor, class ContainerBase *_cb_stat, class ContainerBase * const _cb_scale,
                                      class ContainerBase * const _cb_scale_avg = 0, const bool _enable_favre = false);

          inline const char* id()
          {return id_; }

          inline void setDoNotReset(bool _doNotReset)
          { doNotReset_ = _doNotReset; }

          inline bool doNotReset()
          { return doNotReset_; }

          inline void setWrapPeriodic(bool wrap)
          { wrapPeriodic_ = wrap; }

          inline bool wrapPeriodic()
          { return wrapPeriodic_; }

          inline void id(char *_id);
          inline bool matches_id(const char *_id);
          inline bool matches_any_id(std::list<std::string> * ids);

          virtual bool isDoubleData() = 0;
          virtual bool isIntData() = 0;

          virtual void addZero() = 0;
          virtual void addUninitialized(int n) = 0;
          virtual int size() const = 0;
          virtual int capacity() const = 0;
          virtual int nVec() const = 0;
          virtual int lenVec() const = 0;
          virtual void* begin_slow_dirty() = 0;

          virtual void clearContainer() = 0;

          virtual void copy(int from,int to) = 0;
          virtual void del(int n) = 0;
          virtual void delForward(int n,bool scale,bool translate,bool rotate) = 0;
          virtual void delRestart(int n,bool scale,bool translate,bool rotate) = 0;
          virtual void delRestart(bool scale,bool translate,bool rotate) = 0;
          virtual void clearReverse(bool scale,bool translate,bool rotate) = 0;

          virtual bool setFromContainer(ContainerBase *cont) = 0;

          bool isStatisticsContainer()
          { return (container_statistics_raw_data_!=0); }
          bool calcStatistics();
          bool updateScalingContainer();
          virtual bool calcAvgFromContainer() = 0;
          virtual bool calcMeanSquareFromContainer() = 0;
          virtual bool calcSumFromContainer() = 0;

          virtual void scale(double factor) = 0;
          virtual void move(const double * const dx) = 0;
          virtual void moveElement(const int i, const double * const dx) = 0;
          virtual void rotate(const double * const dQ) = 0;

          virtual void setToDefault(int n) = 0;
          virtual void setAllToZero() = 0;

          inline bool useDefault()
          { return useDefault_ ; }

          inline int getStatLevel() const
          { return statLevel_; }

          bool isScalingContainer() const
          { return scalingContainer_; }

          void setScalingContainer(bool _value)
          { scalingContainer_ = _value; }

          inline void setWeightingFactor(double _value)
          { weighting_factor_ = _value; }

          inline int communicationType() const
          { return communicationType_; }

          // buffer functions for parallelization

          virtual int bufSize(int operation = OPERATION_UNDEFINED,
                            bool scale=false,bool translate=false, bool rotate=false) const = 0;
          virtual int popFromBuffer(double *buf,int operation,
                            bool scale=false,bool translate=false, bool rotate=false) = 0;
          virtual int pushToBuffer(double *buf,int operation,
                            bool scale=false,bool translate=false, bool rotate=false) = 0;

          virtual int elemListBufSize(int n, int operation = OPERATION_UNDEFINED,
                            bool scale=false,bool translate=false, bool rotate=false) = 0;
          virtual int pushElemListToBuffer(int n, int *list, int *wraplist, double *buf, int operation, double *dlo, double *dhi,
                           bool scale=false,bool translate=false, bool rotate=false) = 0;
          virtual int popElemListFromBuffer(int first, int n, double *buf, int operation,
                           bool scale=false,bool translate=false, bool rotate=false) = 0;
          virtual int pushElemListToBufferReverse(int first, int n, double *buf, int operation,
                           bool scale=false,bool translate=false, bool rotate=false) = 0;
          virtual int popElemListFromBufferReverse(int n, int *list, double *buf, int operation,
                           bool scale=false,bool translate=false, bool rotate=false) = 0;

          virtual int elemBufSize(int operation = OPERATION_UNDEFINED,
                            bool scale=false,bool translate=false, bool rotate=false) = 0;
          virtual int pushElemToBuffer(int n, double *buf,int operation,
                            bool scale=false,bool translate=false, bool rotate=false) = 0;
          virtual int popElemFromBuffer(double *buf,int operation,
                            bool scale=false,bool translate=false, bool rotate=false) = 0;

          // static elements
          static const char * AVERAGESUFFIX;
          static const char * MEANSQUARESUFFIX;

     protected:

          ContainerBase(const char *_id, const char* _comm, const char* _ref, const char *_restart,int _scalePower);
          ContainerBase(ContainerBase const &orig);

          inline bool isScaleInvariant() const;
          inline bool isTranslationInvariant() const;
          inline bool isRotationInvariant() const;

          inline bool decidePackUnpackOperation(int operation,bool scale,bool translate, bool rotate) const;

          inline bool decideCommOperation(int operation) const;

          inline bool decideCreateNewElements(int operation);

          char *id_;
          int communicationType_;
          int refFrame_;
          int restartType_;
          int scalePower_;

          bool useDefault_;

          bool doNotReset_;

          class ContainerBase *container_statistics_raw_data_;

          class ContainerBase *container_statistics_scale_data_;
          class ContainerBase *container_statistics_scale_average_data_;

          int statLevel_;
          double weighting_factor_;

          bool scalingContainer_;

          // Enable favre averaged time averaging
          bool enable_favre_;

          // ignore zero values for averaging
          // default is false

          // Decides whether a property is shifted when a periodic boundary is encountered
          // in pushElemListToBuffer
          bool wrapPeriodic_;

     private:

         ContainerBase();
  };

  // *************************************
  #include "container_base_I.h"
  // *************************************

} /* LAMPPS_NS */
#endif /* CONTAINERBASE_H_ */
