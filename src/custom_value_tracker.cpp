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

#include "custom_value_tracker.h"

using namespace LAMMPS_NS;

  /* ----------------------------------------------------------------------
   constructor, destructor
  ------------------------------------------------------------------------- */

  CustomValueTracker::CustomValueTracker(LAMMPS *lmp,AbstractMesh *_ownerMesh)
   : Pointers(lmp),
     ownerMesh_(_ownerMesh),
     capacityElement_(0)
  {
  }

  CustomValueTracker::CustomValueTracker(LAMMPS *lmp)
   : Pointers(lmp),
     ownerMesh_(NULL),
     capacityElement_(0)
  {
  }

  CustomValueTracker::~CustomValueTracker()
  {
  }

  /* ----------------------------------------------------------------------
   memory management
  ------------------------------------------------------------------------- */

  int CustomValueTracker::getCapacity()
  {
    return capacityElement_;
  }

  /* ----------------------------------------------------------------------
   check if all containers have same length
  ------------------------------------------------------------------------- */

  void CustomValueTracker::check_element_property_consistency(int _len)
  {
    if (!elementProperties_.sameLength(_len))
    {
         error->one(FLERR,"all element properties must have the same length.\n"
                          "For meshes, all elem properties with restart must be added in post_create_pre_restart().\n"
                          "For meshes, all elem properties without restart must be added after the constructor()\n");
    }
  }

  /* ----------------------------------------------------------------------
   remove property
  ------------------------------------------------------------------------- */

  void CustomValueTracker::removeElementProperty(const char *_id)
  {
     elementProperties_.remove(_id);
  }

  void CustomValueTracker::removeGlobalProperty(const char *_id)
  {
     globalProperties_.remove(_id);
     globalProperties_orig_.remove(_id);
  }

  /* ----------------------------------------------------------------------
   store current values of global properties as orig
  ------------------------------------------------------------------------- */

  void CustomValueTracker::storeOrig()
  {
      
      globalProperties_.storeOrig(globalProperties_orig_);
      
  }

  /* ----------------------------------------------------------------------
   reset global properties to orig
  ------------------------------------------------------------------------- */

  void CustomValueTracker::resetToOrig()
  {
      
      globalProperties_.reset(globalProperties_orig_);
      
  }

  /* ----------------------------------------------------------------------
   rotate all properties, applies to vector and multivector only
  ------------------------------------------------------------------------- */

  void CustomValueTracker::rotate(double *totalQ,double *dQ)
  {
      
      elementProperties_.rotate(dQ);
      globalProperties_.rotate(totalQ);
  }

  void CustomValueTracker::rotate(double *dQ)
  {
      
      elementProperties_.rotate(dQ);
      globalProperties_.rotate(dQ);
  }

  /* ----------------------------------------------------------------------
   scale all properties, applies to vectors and multivectors only
  ------------------------------------------------------------------------- */

  void CustomValueTracker::scale(double factor)
  {
      
      elementProperties_.scale(factor);
      globalProperties_.scale(factor);
  }

  /* ----------------------------------------------------------------------
   move all properties
  ------------------------------------------------------------------------- */

  void CustomValueTracker::move(double *vecTotal, double *vecIncremental)
  {
      
      elementProperties_.move(vecIncremental);
      globalProperties_.move(vecTotal);
  }

  void CustomValueTracker::move(double *vecIncremental)
  {
      
      elementProperties_.move(vecIncremental);
      globalProperties_.move(vecIncremental);
  }

  /* ----------------------------------------------------------------------
   clear reverse properties, i.e. reset all of them to 0
  ------------------------------------------------------------------------- */

  void CustomValueTracker::clearReverse(bool scale,bool translate,bool rotate)
  {
      
      elementProperties_.clearReverse(scale,translate,rotate);
  }
