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
    Andreas Aigner (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

#else

#ifndef LMP_PAIR_SPH_H
#define LMP_PAIR_SPH_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSph : public Pair {
 public:

  friend class FixSPH;

  PairSph(class LAMMPS *);
  ~PairSph();

  /* INHERITED FROM Pair */

  virtual void compute(int, int) = 0;
  virtual void settings(int, char **) = 0;
  virtual void setKernelAndLength(int, char **);
  virtual void coeff(int, char **) = 0;
  virtual void init_style();
  virtual void init_substyle() = 0;
  virtual void init_list(int, class NeighList *);
  virtual double init_one(int, int);
  virtual void write_restart(FILE *){}
  virtual void read_restart(FILE *){}
  virtual void write_restart_settings(FILE *){}
  virtual void read_restart_settings(FILE *){}
  //virtual void reset_dt();

  /* PUBLIC ACCESS FUNCTIONS */

  int sph_kernel_id(){return kernel_id;}
  int returnPairStyle(){return pairStyle_; };
  double returnViscosity() {return viscosity_; };

 protected:

  void allocate();
  virtual void updatePtrs();
  virtual void updateRadius();
  //virtual double interpDist(double, double);
  inline double interpDist(double disti, double distj) {return 0.5*(disti+distj);}

  class FixPropertyAtom* fppaSl; //fix for smoothing length
  class FixPropertyGlobal* fppaSlType; //fix for per type smoothing length
  double *sl;         // per atom smoothing length
  double **slComType; // common smoothing length in case of mass_type=1
  double sl_0;

  int kernel_id;
  char *kernel_style;

  double *onerad;
  double *maxrad;

  int mass_type; // flag defined in atom_vec*

  int pairStyle_;
  double viscosity_;

  // storage for force part caused by pressure gradient (grad P / rho):
  class FixPropertyAtom* fix_fgradP_;
  double **fgradP_;
};

}

#endif
#endif
