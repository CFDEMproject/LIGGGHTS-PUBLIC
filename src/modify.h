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
    This file is from LAMMPS, but has been modified. Copyright for
    modification:

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz

    Copyright of original file:
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

#ifndef LMP_MODIFY_H
#define LMP_MODIFY_H

#include <stdio.h>
#include "pointers.h"
#include "fix.h"
#include <map>
#include <string>

namespace LAMMPS_NS {

class Modify : protected Pointers {
  friend class Info;
 public:
  int nfix,maxfix;
  int n_pre_initial_integrate, n_initial_integrate,n_post_integrate,n_pre_exchange,n_pre_neighbor;
  int n_pre_force,n_post_force;
  int n_iterate_implicitly, n_pre_final_integrate; 
  int n_final_integrate,n_end_of_step,n_thermo_energy;
  int n_initial_integrate_respa,n_post_integrate_respa;
  int n_pre_force_respa,n_post_force_respa,n_final_integrate_respa;
  int n_min_pre_exchange,n_min_pre_neighbor;
  int n_min_pre_force,n_min_post_force,n_min_energy;

  int restart_pbc_any;       // 1 if any fix sets restart_pbc
  int nfix_restart_global;   // stored fix global info from restart file
  int nfix_restart_peratom;  // stored fix peratom info from restart file

  class Fix **fix;           // list of fixes
  int *fmask;                // bit mask for when each fix is applied

  int timing;                // 1 if fix calls are timed

  int ncompute,maxcompute;   // list of computes
  class Compute **compute;

  Modify(class LAMMPS *);
  virtual ~Modify();
  virtual void init();
  virtual void setup(int);
  virtual void setup_pre_exchange();
  virtual void setup_pre_neighbor();
  virtual void setup_pre_force(int);
  virtual void pre_initial_integrate();
  virtual void initial_integrate(int);
  virtual void post_integrate();
  void pre_decide();
  virtual void pre_exchange();
  virtual void pre_neighbor();
  virtual void pre_force(int);
  virtual void post_force(int);
  virtual void pre_final_integrate();
  virtual void final_integrate();
  virtual bool iterate_implicitly();  
  virtual void end_of_step();
  virtual double thermo_energy();
  virtual void post_run();

  void setup_pre_force_respa(int, int);
  void initial_integrate_respa(int, int, int);
  void post_integrate_respa(int, int);
  void pre_force_respa(int, int, int);
  void post_force_respa(int, int, int);
  void final_integrate_respa(int, int);

  void min_pre_exchange();
  void min_pre_neighbor();
  void min_pre_force(int);
  void min_post_force(int);

  double min_energy(double *);
  void min_store();
  void min_step(double, double *);
  void min_clearstore();
  void min_pushstore();
  void min_popstore();
  int min_reset_ref();
  double max_alpha(double *);
  int min_dof();

  void add_fix(int, char **, char *suffix = NULL);
  void modify_fix(int, char **);
  void delete_fix(const char *,bool unfixflag = false); 
  int find_fix(const char *);

  class FixPropertyGlobal* add_fix_property_global(int narg,char **arg,const char *);
  class FixPropertyAtom* add_fix_property_atom(int narg,char **arg,const char *);
  class Fix* find_fix_property(const char *,const char *,const char *,int ,int,const char * );
  class Fix* find_fix_property(const char *,const char *,const char *,int ,int,const char *,bool );
  class Fix* find_fix_id(const char *id);
  class Compute* find_compute_id(const char *id);
  class Fix* find_fix_id_style(const char *id,const char *style);
  class Fix* find_fix_style(const char *style, int rank);
  class Fix* find_fix_style_strict(const char *style, int rank);
  class Compute* find_compute_style_strict(const char *style, int rank);
  int n_fixes_style(const char *style); 
  int n_computes_style(const char *style); 
  int n_fixes_style_strict(const char *style); 
  int n_fixes_property_atom();
  class FixPropertyAtom* find_fix_property_atom(int rank);
  int n_fixes_property_atom_not_internal();
  int dump_size_fixes_property_atom_not_internal();
  class FixPropertyAtom* find_fix_property_atom_not_internal(int rank);
  bool i_am_first_of_style(class Fix *fix_to_check); 
  int index_first_fix_of_style(const char *style); 
  int index_last_fix_of_style(const char *style); 
  int my_index(class Fix *fixptr);
  int index_first_fix_with_function(const int FUNCTION, bool integrate=false); 
  class FixScalarTransportEquation* find_fix_scalar_transport_equation(const char *equation_id); 
  class FixScalarTransportEquation* find_fix_scalar_transport_equation_strict(const char *equation_id); 
  void box_extent(double &xlo,double &xhi,double &ylo,double &yhi,double &zlo,double &zhi); 

  void add_compute(int, char **, char *suffix = NULL);
  void modify_compute(int, char **);
  void delete_compute(const char *,bool uncomputeflag = false); 
  int find_compute(const char *);
  void clearstep_compute();
  void addstep_compute(bigint);
  void addstep_compute_all(bigint);

  void write_restart(FILE *);
  int read_restart(FILE *);
  void restart_deallocate();

  bigint memory_usage();

  int fix_restart_in_progress();
  bool have_restart_data(Fix *f);
  bool have_restart_data_style(const char* _style);
  int n_restart_data_global_style(const char* _style);
  char* id_restart_data_global_style(const char* _style,int _rank);
  void max_min_rad(double &maxrad,double &minrad); 

  void forceMeshExchange();

  // updates all computes at the end of a run
  void update_computes_on_run_end();

 protected:

  // lists of fixes to apply at different stages of timestep

  int *list_pre_initial_integrate, *list_initial_integrate,*list_post_integrate;
  int *list_pre_exchange,*list_pre_neighbor;
  int *list_pre_force,*list_post_force;
  int *list_iterate_implicitly, *list_pre_final_integrate; 
  int *list_final_integrate,*list_end_of_step,*list_thermo_energy;
  int *list_initial_integrate_respa,*list_post_integrate_respa;
  int *list_pre_force_respa,*list_post_force_respa;
  int *list_final_integrate_respa;
  int *list_min_pre_exchange,*list_min_pre_neighbor;
  int *list_min_pre_force,*list_min_post_force;
  int *list_min_energy;

  int *end_of_step_every;

  int n_timeflag;            // list of computes that store time invocation
  int *list_timeflag;

  char **id_restart_global;           // stored fix global info
  char **style_restart_global;        // from read-in restart file
  char **state_restart_global;

  char **id_restart_peratom;          // stored fix peratom info
  char **style_restart_peratom;       // from read-in restart file
  int *index_restart_peratom;

  int index_permanent;        // fix/compute index returned to library call

  void list_init(int, int &, int *&);
  void list_init_pre_exchange(int, int &, int *&);
  void list_init_end_of_step(int, int &, int *&);
  void list_init_thermo_energy(int, int &, int *&);
  void list_init_compute();

private:
  inline void call_method_on_fixes(FixMethod method);
  inline void call_method_on_fixes(FixMethod method, int *& ilist, int & inum);
  inline void call_method_on_fixes(FixMethodWithVFlag method, int vflag);
  inline void call_method_on_fixes(FixMethodWithVFlag method, int vflag, int *& ilist, int & inum);

  inline void call_respa_method_on_fixes(FixMethodRESPA2 method, int arg1, int arg2, int *& ilist, int & inum);
  inline void call_respa_method_on_fixes(FixMethodRESPA3 method, int arg1, int arg2, int arg3, int *& ilist, int & inum);

  typedef Compute *(*ComputeCreator)(LAMMPS *, int, int, char **);
  std::map<std::string,ComputeCreator> *compute_map;

  typedef Fix *(*FixCreator)(LAMMPS *, int, char **);
  std::map<std::string,FixCreator> *fix_map;

  template <typename T> static Compute *compute_creator(LAMMPS *, int, int, char **);
  template <typename T> static Fix *fix_creator(LAMMPS *, int, char **);
};

}

#endif

/* ERROR/WARNING messages:

W: One or more atoms are time integrated more than once

This is probably an error since you typically do not want to
advance the positions or velocities of an atom more than once
per timestep.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix command before simulation box is defined

The fix command cannot be used before a read_data, read_restart, or
create_box command.

E: Could not find fix group ID

A group ID used in the fix command does not exist.

E: Replacing a fix, but new style != old style

A fix ID can be used a 2nd time, but only if the style matches the
previous fix.  In this case it is assumed you with to reset a fix's
parameters.  This error may mean you are mistakenly re-using a fix ID
when you do not intend to.

W: Replacing a fix, but new group != old group

The ID and style of a fix match for a fix you are changing with a fix
command, but the new group you are specifying does not match the old
group.

E: Invalid fix style

The choice of fix style is unknown.

E: Could not find fix_modify ID

A fix ID used in the fix_modify command does not exist.

E: Could not find fix ID to delete

Self-explanatory.

E: Reuse of compute ID

A compute ID cannot be used twice.

E: Invalid compute style

Self-explanatory.

E: Could not find compute_modify ID

Self-explanatory.

E: Could not find compute ID to delete

Self-explanatory.

*/
