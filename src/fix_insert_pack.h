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

#ifdef FIX_CLASS

FixStyle(insert/pack,FixInsertPack)

#else

#ifndef LMP_FIX_INSERT_PACK_H
#define LMP_FIX_INSERT_PACK_H

#include "fix_insert.h"

namespace LAMMPS_NS {

class FixInsertPack : public FixInsert {
 public:

  FixInsertPack(class LAMMPS *, int, char **);
  ~FixInsertPack();

  void init();
  virtual void restart(char *);

 protected:

  virtual void calc_insertion_properties();
  void init_defaults();

  void calc_region_volume_local();

  virtual int calc_ninsert_this();
  virtual int calc_maxtry(int);
  void x_v_omega(int,int&,int&,double&);
  double insertion_fraction();

  int is_nearby(int);
  int is_nearby_body(int);

  // region to be used for insertion
  class Region *ins_region;
  char *idregion;
  double region_volume,region_volume_local;
  int ntry_mc;

  // target that region should fulfil after each insertion
  double volumefraction_region;
  int ntotal_region;
  double masstotal_region;

  // ratio how many particles have been inserted
  double insertion_ratio;

  // warn if region extends outside box
  bool warn_region;

};

}

#endif
#endif
