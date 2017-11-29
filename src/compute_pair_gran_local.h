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
//#include "pair_gran.h"
//#include "pair_gran_proxy.h"

namespace LAMMPS_NS {

class PairGran;

class ComputePairGranLocal : public Compute {

 public:
  ComputePairGranLocal(class LAMMPS *, int &iarg, int, char **);
  ~ComputePairGranLocal();
  void post_create();
  void pre_delete(bool uncomputeflag);
  void init();
  virtual void init_cpgl(bool requestflag);
  void init_list(int, class NeighList *);
  void compute_local();
  double memory_usage();
  void reference_deleted();
  virtual void add_pair(int i,int j,double fx,double fy,double fz,double tor1,double tor2,double tor3,double *hist, const double * const contact_point);
  virtual void add_heat(int i,int j,double hf);
  virtual void add_wall_1(int iFMG,int iTri,int iP,double *contact_point,double *v_wall);
  virtual void add_wall_2(int i,double fx,double fy,double fz,double tor1,double tor2,double tor3,double *hist,double rsq, double *normal_);
  virtual void add_heat_wall(int i,double hf);

  virtual void pair_finalize();
  int get_history_offset(const char * const name);

  /* inline access */

  virtual bool decide_add(double *hist, double * &contact_pos)
  { return true; }

  inline int get_nvalues()
  { return nvalues; }

  inline double** get_data()
  { return array; }

  inline int get_ncount()
  { return ncount_added_via_pair; }

  virtual int offset_x1()
  { return (posflag > 0 ? 0 : -1); }

  virtual int offset_x2()
  { return (posflag > 0 ? 3 : -1); }

  virtual int offset_v1()
  { return (velflag > 0 ? posflag*6 : -1); }

  virtual int offset_v2()
  { return (velflag > 0 ? posflag*6+3 : -1); }

  virtual int offset_id1()
  { return (idflag > 0 ? posflag*6+velflag*6 : -1);}

  virtual int offset_id2()
  { return (idflag > 0 ? posflag*6+velflag*6+1 : -1);}

  virtual int offset_f()
  { return (fflag > 0 ? posflag*6+velflag*6+idflag*3 : -1);}

  virtual int offset_fn()
  { return (fnflag > 0 ? posflag*6+velflag*6+idflag*3+fflag*3 : -1);}

  virtual int offset_ft()
  { return (ftflag > 0 ? posflag*6+velflag*6+idflag*3+fflag*3+fnflag*3 : -1);}

  virtual int offset_torque()
  { return (torqueflag > 0 ? posflag*6+velflag*6+idflag*3+fflag*3+fnflag*3+ftflag*3 : -1);}

  virtual int offset_torquen()
  { return (torquenflag > 0 ? posflag*6+velflag*6+idflag*3+fflag*3+fnflag*3+ftflag*3+torqueflag*3 : -1);}

  virtual int offset_torquet()
  { return (torquetflag > 0 ? posflag*6+velflag*6+idflag*3+fflag*3+fnflag*3+ftflag*3+torqueflag*3+torquenflag*3 : -1);}

  virtual int offset_history()
  { return (histflag > 0 ? posflag*6+velflag*6+idflag*3+fflag*3+fnflag*3+ftflag*3+torqueflag*3+torquenflag*3+torquetflag*3 : -1);}

  virtual int offset_area()
  { return (areaflag > 0 ? posflag*6+velflag*6+idflag*3+fflag*3+fnflag*3+ftflag*3+torqueflag*3+torquenflag*3+torquetflag*3+histflag*dnum : -1);}

  virtual int offset_delta()
  { return (deltaflag > 0 ? posflag*6+velflag*6+idflag*3+fflag*3+fnflag*3+ftflag*3+torqueflag*3+torquenflag*3+torquetflag*3+histflag*dnum+areaflag*1 : -1);}

  virtual int offset_heat()
  { return (heatflag > 0 ? posflag*6+velflag*6+idflag*3+fflag*3+fnflag*3+ftflag*3+torqueflag*3+torquenflag*3+torquetflag*3+histflag*dnum+areaflag*1 + deltaflag*1 : -1);}

  virtual int offset_contact_point()
  {
    return (cpflag > 0 ?   posflag*6+velflag*6+idflag*3+fflag*3+fnflag*3+ftflag*3+torqueflag*3+torquenflag*3+torquetflag*3+histflag*dnum+areaflag*1 + deltaflag*1 + heatflag*1 : -1);
    //return offset_history() + get_history_offset("contact_point_offset");
  }

  virtual int offset_ms_id1()
  {
    return (msidflag > 0 ?   posflag*6+velflag*6+idflag*3+fflag*3+fnflag*3+ftflag*3+torqueflag*3+torquenflag*3+torquetflag*3+histflag*dnum+areaflag*1 + deltaflag*1 + heatflag*1 + cpflag*3 : -1);
  }

  virtual int offset_ms_id2()
  {
    return (msidflag > 0 ?   posflag*6+velflag*6+idflag*3+fflag*3+fnflag*3+ftflag*3+torqueflag*3+torquenflag*3+torquetflag*3+histflag*dnum+areaflag*1 + deltaflag*1 + heatflag*1 + cpflag*3 +1 : -1);
  }

 protected:

  int nvalues;      // number of double values per entry
  int ncount;       // count of eligible pair - all who are eligible for surfacesIntersect or surfacesClose

  int ncount_added_via_pair; // count actually added via call from pair_gran
                             // might be lower than ncount because is based on hasForceUpdate occurrences
                             // not all the surfacesClose calls have hasForceUpdate=true

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
  class FixMultisphere *fix_ms;

  int ipair;

  int posflag,velflag,idflag,fflag,fnflag,ftflag,torqueflag,torquenflag,torquetflag,histflag,areaflag,deltaflag,heatflag,cpflag,msidflag;

  bool   verbose;

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
