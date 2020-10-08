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

    Andreas Aigner (JKU Linz)
    Andreas Eitzlmayr (TU Graz)

    Copyright 2009-2012 JKU Linz
    Copyright 2013-     TU Graz
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(sph/morris/tensCorr,PairSphMorrisTenscorr)

#else

#ifndef LMP_PAIR_SPH_MORRIS_TENSCORR_H
#define LMP_PAIR_SPH_MORRIS_TENSCORR_H

#include "pair_sph.h"

namespace LAMMPS_NS {

class PairSphMorrisTenscorr : public PairSph {

 friend class FixSPH;

 public:

  PairSphMorrisTenscorr(class LAMMPS *);
  ~PairSphMorrisTenscorr();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_substyle();
  void write_restart(FILE *);
  void read_restart(FILE *, const int major, const int minor);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *, const int major, const int minor);

 private:
  class FixPropertyAtom* fix_dvdx_;
  class FixPropertyAtom* fix_dvdy_;
  class FixPropertyAtom* fix_dvdz_;
  class FixPropertyAtom* fix_gamma_; // shear rate magnitude
  class FixPropertyAtom* fix_visc_; // local viscosity value

  double **dvdx_;
  double **dvdy_;
  double **dvdz_;
  double *gamma_;
  double *visc_;

 protected:
  void allocate();
//  void allocate_properties(int);
//  double artificialViscosity(int, int, int, int, double, double, double, double, double, double, double, double **);
//  template <int> void tensileCorrection(int, int, double, double, double, double, double, double, double, double &, double &);
  template <int> void compute_eval(int, int);

  double  **wDeltaPTypeinv;

  int     tensCorr_flag; // flags for additional styles

  int     modelStyle; // 1 ... Newtonian, 2 ... power law, 3 ... Carreau model

  double  dynVisc; // dynamic viscosity (Newtonian)
  double  conIdx; // consistency index (power law)
  double  powIdx; // power law index
  double  OneMpowIdx; // 1 - powIdx
  double  etaMax; // maximum viscosity
  double  etaMin; // minimum viscosity

  double  eta0; // viscosity at zero shear rate (Carreau)
  double  etaInf; // viscosity at infinite shear rate (Carreau)
  double  eta0_inf; // eta0 - etainf
  double  lambda; // critical shear rate (Carreau)
  double  aExp; // exponent a (Carreau)
  double  OneMpowIdx_a; // (1 - powIdx)/a

  double  epsilon,deltaP; // coeffs for tensile correction

};

}

#endif
#endif
