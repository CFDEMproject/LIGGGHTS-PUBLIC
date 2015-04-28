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

#else

#ifndef LMP_FIX_INSERT_H
#define LMP_FIX_INSERT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixInsert : public Fix {
 public:
  FixInsert(class LAMMPS *, int, char **);
  ~FixInsert();

  virtual int setmask();
  virtual void init();
  virtual void setup_pre_exchange() {}
  void setup(int vflag);
  virtual double extend_cut_ghost();
  void pre_exchange();
  virtual void end_of_step() {}

  void write_restart(FILE *);
  virtual void restart(char *);

  virtual void reset_timestep(bigint) {}

  double compute_vector(int index);

  virtual double min_rad(int);  
  virtual double max_rad(int);  
  virtual double max_r_bound();  
  int min_type();
  int max_type();

  int ins_every()
  { return insert_every; }

 protected:

  int iarg;

  /*---TARGET QUANTITIES---how much will be inserted*/

  // flag if total # of particle to insert by this command exists
  int ninsert_exists;

  // insertion target - can bei either number or mass
  // in case of multi-sphere, ninsert is number of bodies to insert
  int ninsert;
  double massinsert;

  // insertion rate - can bei either number rate or mass rate
  // in particle per time and mass per time units - both can be 0
  double nflowrate;
  double massflowrate;

  // how much we have inserted to far
  int ninserted;
  double massinserted;

  // max allocation for insertion
  int ninsert_this_max_local;

  // first, most recent, final insertion step
  int first_ins_step,most_recent_ins_step, final_ins_step;

  /*---INSERTION QUANTITIES---what, where and how exactly will we insert*/

  //particle distribution
  class FixParticledistributionDiscrete *fix_distribution;

  // insert 'ninsert_per' particles every 'insert_every' steps
  // 'ninsert_per' is the default, actual # of inserted particles
  // to insert is defined in calc_ninsert_this()
  int insert_every;
  double ninsert_per;

  // determines how particle distributions are interpreted
  // total number of particles inserted for each vary with calc_ninsert_this() result
  // if exact_number = 0, total # inserted particles will deviate from calc_ninsert_this() result
  //                      distribution will be fulfilled statistically
  // if exact_number = 1, total # inserted particles will match result from calc_ninsert_this()
  int exact_number;

  // max and min type to be inserted
  int type_max,type_min;

  // minmum/maximum radius to be inserted
  double minrad, maxrad;

  // flag if overlap is checked upon insertion (via all-to-all comm)
  int check_ol_flag;

  // if flag is 1, particles are generated to be in the region as a whole
  // if flag is 0, particles centers are in region
  int all_in_flag;

  int maxattempt;

  // positions generated, and for overlap check
  int nspheres_near;
  double **xnear;

  // velocity and ang vel distribution
  // currently constant for omega - could also be a distribution
  int    v_randomSetting;
  double v_insert[3];
  double v_insertFluct[3];
  double omega_insert[3];
  double quat_insert[4];
  bool quat_random_;

  /*---FURTHER THINGS THAT WE NEED---*/

  // random generator
  class RanPark *random;
  int seed;

  // all2all comm
  int me,nprocs;
  int *recvcounts,*displs;

  // determine if print stats
  int print_stats_start_flag;
  int print_stats_during_flag;

  // warn if box extent too small for insertion
  bool warn_boxentent;

  class FixMultisphere *fix_multisphere;
  class MultisphereParallel *multisphere;

  /*---FUNCTION MEMBERS---*/

  virtual void print_stats_start();
  virtual void print_stats_during(int,double);

  virtual void init_defaults();
  virtual void sanity_check();
  virtual void calc_insertion_properties() = 0;

  virtual void pre_insert() {};
  virtual int calc_ninsert_this();
  virtual int load_xnear(int);
  virtual int count_nnear();
  virtual int is_nearby(int) = 0;

  virtual void x_v_omega(int,int&,int&,double&) = 0;
  virtual double insertion_fraction() = 0;

  virtual void finalize_insertion(int){};

 private:

  char *property_name;
  class FixPropertyAtom *fix_property;
  double fix_property_value;

  bool setup_flag;

  virtual int distribute_ninsert_this(int);
};

}

#endif
#endif
