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

#include "string.h"

namespace LAMMPS_NS
{
  // buffer operation types (for push and pop)

  enum{ OPERATION_COMM_EXCHANGE,
        OPERATION_COMM_BORDERS,
        OPERATION_COMM_FORWARD,
        OPERATION_COMM_REVERSE,
        OPERATION_RESTART,
        OPERATION_UNDEFINED};

  class ContainerBase
  {
      public:

          ContainerBase(const char *_id);

          virtual ~ContainerBase();

          void setProperties(const char *_id, const char* _comm, const char* _ref, const char *_restart,int _scalePower = 1);
          bool propertiesSetCorrectly();

          inline void id(char *_id);
          inline bool matches_id(const char *_id);

          virtual bool isDoubleData() = 0;
          virtual bool isIntData() = 0;

          virtual void addZero() = 0;
          virtual void addUninitialized(int n) = 0;
          virtual int size() = 0;
          virtual int nVec() = 0;
          virtual int lenVec() = 0;
          virtual void* begin_slow_dirty() = 0;

          virtual void copy(int from,int to) = 0;
          virtual void del(int n) = 0;
          virtual void delForward(int n,bool scale,bool translate,bool rotate) = 0;
          virtual void delRestart(int n,bool scale,bool translate,bool rotate) = 0;
          virtual void delRestart(bool scale,bool translate,bool rotate) = 0;
          virtual void clearReverse(bool scale,bool translate,bool rotate) = 0;

          virtual bool setFromContainer(ContainerBase *cont) = 0;

          virtual void scale(double factor) = 0;
          virtual void move(double *dx) = 0;
          virtual void moveElement(int i,double *dx) = 0;
          virtual void rotate(double *dQ) = 0;

          virtual void setToDefault(int n) = 0;

          inline bool useDefault()
          { return useDefault_ ; }

          // buffer functions for parallelization

          virtual int bufSize(int operation = OPERATION_UNDEFINED,
                            bool scale=false,bool translate=false, bool rotate=false) = 0;
          virtual int popFromBuffer(double *buf,int operation,
                            bool scale=false,bool translate=false, bool rotate=false) = 0;
          virtual int pushToBuffer(double *buf,int operation,
                            bool scale=false,bool translate=false, bool rotate=false) = 0;

          virtual int elemListBufSize(int n, int operation = OPERATION_UNDEFINED,
                            bool scale=false,bool translate=false, bool rotate=false) = 0;
          virtual int pushElemListToBuffer(int n, int *list, double *buf, int operation,
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

     protected:

          ContainerBase(const char *_id, const char* _comm, const char* _ref, const char *_restart,int _scalePower);
          ContainerBase(ContainerBase const &orig);

          inline bool isScaleInvariant();
          inline bool isTranslationInvariant();
          inline bool isRotationInvariant();

          inline bool decidePackUnpackOperation(int operation,bool scale,bool translate, bool rotate);

          inline bool decideCommOperation(int operation);

          inline bool decideCreateNewElements(int operation);

          char *id_;
          int communicationType_;
          int refFrame_;
          int restartType_;
          int scalePower_;

          bool useDefault_;

     private:

         ContainerBase();
  };

  // *************************************
  #include "container_base_I.h"
  // *************************************

} /* LAMPPS_NS */
#endif /* CONTAINERBASE_H_ */
