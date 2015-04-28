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

FixStyle(couple/cfd,FixCfdCoupling)

#else

#ifndef LMP_FIX_CFD_COUPLING_H
#define LMP_FIX_CFD_COUPLING_H

#include "fix.h"

namespace LAMMPS_NS {

class FixCfdCoupling : public Fix {
 friend class CfdRegionmodel;
 friend class CfdDatacoupling;

 public:

  FixCfdCoupling(class LAMMPS *, int, char **);
  ~FixCfdCoupling();
  void post_create();

  int setmask();
  void init();
  virtual void setup(int);
  virtual void min_setup(int);
  void end_of_step();
  void post_force_respa(int vflag, int ilevel, int iloop);
  void min_post_force(int);

  // pushing and pulling of properties
  //void pull(char *name,char *type,void *&ptr);
  //void push(char *name,char *type,void *&ptr);
  void add_push_property(const char *name, const char *type);
  void add_pull_property(const char *name, const char *type);
  void check_datatransfer();

  int coupleThis() {return couple_this_;}

  class CfdDatacoupling* get_dc(){return dc_;}
  
 protected:

  int iarg_;

  // data transfer is handled by this class
  class CfdDatacoupling *dc_;

 private:

  int couple_this_;

  // couple every couple_nevery_ timesteps
  // not used in case of MPI coupling
  int couple_nevery_,ts_create_;

  // regionmodels
  class CfdRegionmodel *rm_;

  int nlevels_respa;
};

}

#endif
#endif
