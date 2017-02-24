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
    Andreas Aigner (JKU Linz))

    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef LMP_FIX_SPH
#define LMP_FIX_SPH

#include "fix.h"

namespace LAMMPS_NS {

class FixSph : public Fix {
 public:
  FixSph(class LAMMPS *, int, char **);
  ~FixSph();
  int setmask();
  virtual void updatePtrs();
  void init();
  void init_list(int, class NeighList *);
  virtual void post_integrate() {};
  virtual void post_integrate_respa(int, int);

  int get_kernel_id(){return kernel_id;};
  inline void set_kernel_id(int newid){kernel_id = newid;};

  int kernel_flag;        // 1 if Fix uses sph kernel, 0 if not

 protected:
  inline double interpDist(double disti, double distj) {return 0.5*(disti+distj);};

  class FixPropertyAtom* fppaSl; //smoothing length
  class FixPropertyGlobal* fppaSlType; //per type smoothing length
  double *sl;         // per atom smoothing length
  double **slComType; // common smoothing length in case of mass_type=1

  int kernel_id;
  double kernel_cut;
  char *kernel_style;

  class NeighList *list;
  int nlevels_respa;

  int mass_type; // flag defined in atom_vec*

};

}

#endif
