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

    Christoph Kloss (DCS Computing GmbH, Linz, JKU Linz)
    Philippe Seil (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef LMP_MULTI_VECTOR_CONTAINER
#define LMP_MULTI_VECTOR_CONTAINER

#include "general_container.h"

namespace LAMMPS_NS{
  template<typename T, int NUM_VEC, int LEN_VEC>
  class MultiVectorContainer : public GeneralContainer <T, NUM_VEC, LEN_VEC>
  {
      public:
          MultiVectorContainer(const char *_id);
          MultiVectorContainer(const char *_id, const char *_comm, const char *_ref, const char *_restart, int _scalePower = 1);
          MultiVectorContainer(MultiVectorContainer<T,NUM_VEC,LEN_VEC> const &orig);
          virtual ~MultiVectorContainer();
  };

  /* ----------------------------------------------------------------------
   constructors
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  MultiVectorContainer<T,NUM_VEC,LEN_VEC>::MultiVectorContainer(const char *_id)
  : GeneralContainer<T,NUM_VEC,LEN_VEC>(_id)
  {

  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  MultiVectorContainer<T,NUM_VEC,LEN_VEC>::MultiVectorContainer(const char *_id, const char* _comm, const char* _ref, const char *_restart, int _scalePower)
  : GeneralContainer<T,NUM_VEC,LEN_VEC>(_id, _comm, _ref,_restart,_scalePower)
  {

  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  MultiVectorContainer<T,NUM_VEC,LEN_VEC>::MultiVectorContainer(MultiVectorContainer<T,NUM_VEC,LEN_VEC> const &orig)
  : GeneralContainer<T,NUM_VEC,LEN_VEC>(orig)
  {

  }

  /* ----------------------------------------------------------------------
   destructor
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  MultiVectorContainer<T,NUM_VEC,LEN_VEC>::~MultiVectorContainer()
  {

  }

} /* LAMMPS_NS */
#endif /* MULTIVECTORCONTAINER_H_ */
