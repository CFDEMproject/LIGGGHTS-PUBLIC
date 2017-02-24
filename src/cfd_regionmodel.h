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

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef LMP_CFD_REGIONMODEL_H
#define LMP_CFD_REGIONMODEL_H

#include "pointers.h"

namespace LAMMPS_NS {

class CfdRegionmodel : protected Pointers {
 public:
  CfdRegionmodel(class LAMMPS *lmp, int jarg, int narg, char **arg,class FixCfdCoupling* fc) : Pointers(lmp)
  {
    UNUSED(narg);
    UNUSED(jarg);
    UNUSED(arg);
      this->fc = fc;
  }
  ~CfdRegionmodel() {}

  int get_iarg() {return iarg;}
  bool liggghts_is_active;

  virtual void init() {};
  virtual void rm_update() {};

 protected:
  int iarg;
  class FixCfdCoupling *fc;

  void add_push_property(const char *name, const char *type)
  {
     fc->add_push_property(name,type);
  }

  void add_pull_property(const char *name, const char *type)
  {
     fc->add_pull_property(name,type);
  }
};

}

#endif
