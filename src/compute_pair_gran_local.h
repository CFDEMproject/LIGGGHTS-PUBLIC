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

#ifdef COMPUTE_CLASS

ComputeStyle(pair/gran/local,ComputePairGranLocal)
ComputeStyle(wall/gran/local,ComputePairGranLocal)

#else

#ifndef LMP_COMPUTE_PAIR_GRAN_LOCAL_H
#define LMP_COMPUTE_PAIR_GRAN_LOCAL_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePairGranLocal : public Compute {

 public:
  ComputePairGranLocal(class LAMMPS *, int, char **);
  ~ComputePairGranLocal();
  void post_create();
  void pre_delete(bool uncomputeflag);
  void init();
  virtual void init_cpgl(bool requestflag);
  void init_list(int, class NeighList *);
  void compute_local();
  double memory_usage();
  void reference_deleted();
  virtual void add_pair(int i,int j,double fx,double fy,double fz,double tor1,double tor2,double tor3,double *hist);
  virtual void add_heat(int i,int j,double hf);
  virtual void add_wall_1(int iFMG,int iTri,int iP,double *contact_point,double *v_wall);
  virtual void add_wall_2(int i,double fx,double fy,double fz,double tor1,double tor2,double tor3,double *hist,double rsq);
  virtual void add_heat_wall(int i,double hf);

 protected:

  int nvalues;
  int ncount;
  int newton_pair;

  // if 0, pair data is extracted
  // if 1, wall data is extracted
  int wall;

  // flag if compute is working
  // can be set via pair, fix
  int reference_exists;

  // pointers to classes holding the data
  class PairGran *pairgran;
  class FixHeatGranCond *fixheat;
  class FixWallGran *fixwall;

  int ipair;

  int posflag,velflag,idflag,fflag,fnflag,ftflag,tflag,hflag,aflag,dflag,hfflag;

  bool   verbose;
  double extraSurfDistance;

  int dnum;

  int nmax;
  double *vector;
  double **array;

  class NeighList *list;

  virtual int count_pairs(int &nCountWithOverlap);
  int count_wallcontacts(int &nCountWithOverlap);
  void reallocate(int);
};

}

#endif
#endif
