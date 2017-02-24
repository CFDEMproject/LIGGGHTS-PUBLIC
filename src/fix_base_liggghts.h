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
    (if no contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2015-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

//non-constructable

#else

#ifndef LMP_FIX_BASE_LIGGGHTS_H
#define LMP_FIX_BASE_LIGGGHTS_H

#include "fix.h"
#include "container.h"
#include "region.h"

namespace LAMMPS_NS {

class FixBaseLiggghts : public Fix {

 public:

  FixBaseLiggghts(class LAMMPS *, int, char **);
  virtual ~FixBaseLiggghts();

  void process_region(char *regid);

  virtual void init();
  virtual void setup(int vflag);

  void do_support_multisphere()
  { support_ms_ = true; }

  void do_support_respa()
  { support_respa_ = true; }

  void do_not_need_radius()
  { do_need_radius_ = false; }

  void do_not_need_mass()
  { do_need_mass_ = false; }

  inline Region *region()
  {return region_;}

 protected:

  void count_eligible(double &mass_counted,double &volume_counted, int &nparticles_counted);

  bool support_respa_;
  int nlevels_respa_;

  bool do_need_radius_;
  bool do_need_mass_;

  bool support_ms_;
  class FixMultisphere *fix_ms_;
  class Multisphere *ms_;

  // region to be used for replacement
  char *idregion_;
  class Region *region_;
  int iregion_;
};

}

#endif
#endif
