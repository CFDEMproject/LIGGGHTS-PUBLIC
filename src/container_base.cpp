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

#include "container_base.h"
#include <string.h>

#define GROW 100

using namespace LAMMPS_NS;

  /* ----------------------------------------------------------------------
   constructor
  ------------------------------------------------------------------------- */

  ContainerBase::ContainerBase()
  : id_(0),
    communicationType_(COMM_TYPE_MANUAL),
    refFrame_(REF_FRAME_UNDEFINED),
    restartType_(RESTART_TYPE_UNDEFINED),
    scalePower_(-1),
    useDefault_(false)
  {
  }

  ContainerBase::ContainerBase(const char *_id)
  : id_(0),
    communicationType_(COMM_TYPE_MANUAL),
    refFrame_(REF_FRAME_UNDEFINED),
    restartType_(RESTART_TYPE_UNDEFINED),
    scalePower_(-1),
    useDefault_(false)
  {
      if(_id)
      {
        id_ = new char[strlen(_id)+1];
        strcpy(id_,_id);
      }
  }

  ContainerBase::ContainerBase(const char *_id, const char* _comm, const char* _ref, const char *_restart,int _scalePower)
  : id_(0),
    communicationType_(COMM_TYPE_MANUAL),
    refFrame_(REF_FRAME_UNDEFINED),
    restartType_(RESTART_TYPE_UNDEFINED),
    scalePower_(-1),
    useDefault_(false)
  {
          setProperties(_id, _comm, _ref,_restart,_scalePower);
  }

  ContainerBase::ContainerBase(ContainerBase const &orig)
  :  id_(0),
     communicationType_(orig.communicationType_),
     refFrame_(orig.refFrame_),
     restartType_(orig.restartType_),
     scalePower_(orig.scalePower_),
     useDefault_(orig.useDefault_)
  {

  }

  ContainerBase::~ContainerBase()
  {
      if(id_) delete []id_;
  }

  /* ----------------------------------------------------------------------
   set comm and reference properties
  ------------------------------------------------------------------------- */

  void ContainerBase::setProperties(const char *_id, const char* _comm, const char* _ref, const char *_restart, int _scalePower)
  {
      id_ = new char[strlen(_id)+1];
      strcpy(id_,_id);

      if      (strcmp(_comm,"comm_forward") == 0) communicationType_ = COMM_TYPE_FORWARD;
      else if (strcmp(_comm,"comm_forward_from_frame") == 0) communicationType_ = COMM_TYPE_FORWARD_FROM_FRAME;
      else if (strcmp(_comm,"comm_reverse") == 0) communicationType_ = COMM_TYPE_REVERSE;
      else if (strcmp(_comm,"comm_exchange_borders") == 0) communicationType_ = COMM_EXCHANGE_BORDERS;
      else if (strcmp(_comm,"comm_none") == 0) communicationType_ = COMM_TYPE_NONE;
      else if (strcmp(_comm,"comm_manual") == 0) communicationType_ = COMM_TYPE_MANUAL;
      else communicationType_ = COMM_TYPE_UNDEFINED;

      if      (strcmp(_ref,"frame_invariant") == 0) refFrame_ = REF_FRAME_INVARIANT;
      else if (strcmp(_ref,"frame_trans_rot_invariant") == 0) refFrame_ = REF_FRAME_TRANS_ROT_INVARIANT;
      else if (strcmp(_ref,"frame_scale_trans_invariant") == 0) refFrame_ = REF_FRAME_SCALE_TRANS_INVARIANT;
      else if (strcmp(_ref,"frame_trans_invariant") == 0) refFrame_ = REF_FRAME_TRANS_INVARIANT;
      else if (strcmp(_ref,"frame_general") == 0) refFrame_ = REF_FRAME_GENERAL;
      else refFrame_ = REF_FRAME_UNDEFINED;

      if      (strcmp(_restart,"restart_yes") == 0) restartType_ = RESTART_TYPE_YES;
      else if (strcmp(_restart,"restart_no") == 0) restartType_ = RESTART_TYPE_NO;
      else restartType_ = RESTART_TYPE_UNDEFINED;

      scalePower_ = _scalePower;
  }

  bool ContainerBase::propertiesSetCorrectly()
  {
      if(refFrame_ == REF_FRAME_UNDEFINED ||
        communicationType_ == COMM_TYPE_UNDEFINED ||
        restartType_ == RESTART_TYPE_UNDEFINED ||
        scalePower_ < 0)
            return false;

      return true;
  }
