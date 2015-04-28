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

    Christoph Kloss (DCS Computing GmbH, Linz, JKU Linz)
    Philippe Seil (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef LMP_CONTAINER_BASE_I_H
#define LMP_CONTAINER_BASE_I_H

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
        // no comm at all
        
        COMM_TYPE_NONE,
        // undefined state for error check
        COMM_TYPE_UNDEFINED};  // communication types

  // restart types

  enum{ RESTART_TYPE_UNDEFINED,
        RESTART_TYPE_YES,
        RESTART_TYPE_NO};

  /* ----------------------------------------------------------------------
   decide if property is pushed or pulled at all
  ------------------------------------------------------------------------- */

  inline bool ContainerBase::decidePackUnpackOperation(int operation,bool scale,bool translate, bool rotate)
  {
      // return true for manual communication, such as for node_, node_orig_
      // etc in MultiNodeMeshParallel
      if(COMM_TYPE_MANUAL == communicationType_)
        return true;

      if(OPERATION_RESTART == operation)
      {
          if(restartType_ == RESTART_TYPE_YES)
            return true;
          return false;
      }

      if(OPERATION_COMM_BORDERS == operation ||
         OPERATION_COMM_EXCHANGE == operation )
        return true;

      if(COMM_TYPE_NONE == communicationType_)
        return false;

      if(OPERATION_COMM_REVERSE == operation &&
         COMM_TYPE_REVERSE == communicationType_)
        return true;

      if(OPERATION_COMM_FORWARD == operation &&
         COMM_TYPE_FORWARD == communicationType_)
        return true;

      if(OPERATION_COMM_FORWARD == operation &&
         COMM_TYPE_FORWARD_FROM_FRAME == communicationType_)
      {
         if(scale && !isScaleInvariant())
           return true;
         if(translate && !isTranslationInvariant())
           return true;
         if(rotate && !isRotationInvariant())
           return true;

         return false;
      }

      // default
      return false;
  }

  /* ----------------------------------------------------------------------
   decide if operation performs data communication
  ------------------------------------------------------------------------- */

  inline bool ContainerBase::decideCommOperation(int operation)
  {
      
      if(operation == OPERATION_RESTART)
          return true;

      if(operation == OPERATION_COMM_FORWARD ||
         operation == OPERATION_COMM_REVERSE )
        return true;

      if(operation == OPERATION_COMM_BORDERS ||
         operation == OPERATION_COMM_EXCHANGE )
      {
          
          if(communicationType_ == COMM_TYPE_NONE ||
             communicationType_ == COMM_TYPE_REVERSE)
             return false;

          return true;
      }

      // default
      return true;
  }

  /* ----------------------------------------------------------------------
   decide if unpack creates new element or overwrites existing data
  ------------------------------------------------------------------------- */

  inline bool ContainerBase::decideCreateNewElements(int operation)
  {
      
      if(operation == OPERATION_RESTART)
          return true;

      if(operation == OPERATION_COMM_BORDERS ||
         operation == OPERATION_COMM_EXCHANGE )
        return true;

      if(operation == OPERATION_COMM_FORWARD ||
         operation == OPERATION_COMM_REVERSE )
        return false;

      // default
      return false;
  }

  /* ----------------------------------------------------------------------
   fast test for reference frame
   note that rotation is only carried out for LEN_VEC==3
  ------------------------------------------------------------------------- */

    bool ContainerBase::isScaleInvariant()
    {
       return ( refFrame_ == REF_FRAME_INVARIANT ||
                refFrame_ == REF_FRAME_SCALE_TRANS_INVARIANT);
    }

    bool ContainerBase::isTranslationInvariant()
    {
        return ( refFrame_ == REF_FRAME_INVARIANT ||
                 refFrame_ == REF_FRAME_TRANS_ROT_INVARIANT ||
                 refFrame_ == REF_FRAME_SCALE_TRANS_INVARIANT ||
                 refFrame_ == REF_FRAME_TRANS_INVARIANT);
    }

    bool ContainerBase::isRotationInvariant()
    {
        return ( refFrame_ == REF_FRAME_INVARIANT ||
                 refFrame_ == REF_FRAME_TRANS_ROT_INVARIANT ||
                 lenVec() != 3);
    }

  /* ----------------------------------------------------------------------
   ID operations
  ------------------------------------------------------------------------- */

  inline void ContainerBase::id(char *_id)
  {
      strcpy(_id,id_);
  }

  inline bool ContainerBase::matches_id(const char *_id)
  {
      if(strcmp(_id,id_) == 0) return true;
      return false;
  }

#endif
