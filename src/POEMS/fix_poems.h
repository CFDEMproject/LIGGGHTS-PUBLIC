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
    This file is from LAMMPS
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.

    POEMS and the POEMS fix has been re-worked by Stefan Radl and 
    Mingqiu WU (TU Graz) to be integrated with LIGGGHTS-TUG's fibre modules
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(poems,FixPOEMS)

#else

#ifndef LMP_FIX_POEMS_H
#define LMP_FIX_POEMS_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPOEMS : public Fix  {
 public:
  FixPOEMS(class LAMMPS *, int narg, char **arg);
  ~FixPOEMS();
  virtual void post_create();
  int setmask();
  int model_switch_value();
  void init();
  void setup(int);
  void initial_integrate(int);
  void post_force(int);
  void final_integrate();
  void initial_integrate_respa(int, int, int);
  void post_force_respa(int, int, int);
  void final_integrate_respa(int, int);

  void grow_arrays(int);
  void copy_arrays(int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  double memory_usage();

  void pre_neighbor();
  int dof(int);
  void deform(int);
  void reset_dt();

  void updatePtrs();

 protected:
  class FixPropertyAtom* fix_xcm;
  class FixPropertyAtom* fix_orientation;

 private:
  bool mydebug;
  int me;
  double dtv,dtf,dthalf;
  double *step_respa;
  int nlevels_respa;
  double total_ke;

  // atom assignment to rigid segments
  // double count joint atoms as being in multiple segments

  int *natom2segment;         // # of segments each atom is part of
  int **atom2segment;         // list of segments each atom is part of
  double **displace;       // atom displace in segment coords for 1st segment it's in

  // rigid segment properties
  // only nrigid double counts joint atoms as being in multiple segments
  // other quantities only count a joint atom as being in 1st segment
  int nsystem;               // # of systems 
  int nsegment;                // # of rigid segments
  int *nrigid;              // # of atoms in each rigid segment
  double *masstotal;        // total mass of each rigid segment
  double **xcm;             // coords of center-of-mass of each rigid segment
  double **vcm;             // velocity of center-of-mass of each rigid segment
  double **fcm;             // force on center-of-mass of each rigid segment
  double **inertia;         // 3 inertia components of each rigid segment (xx,yy,zz,xy,yz,xz)
  double **ex_space,**ey_space,**ez_space;
                            // orientation of each rigid segments' principal axes in space coords
  double **angmom;          // angular momentum of each rigid segment in space coords
  double **omega;           // angular velocity of each rigid segment in space coords
  double **torque;          // torque on each rigid segment in space coords
  double **sum,**all;       // work vectors

  // joint attributes between pairs of rigid segments

  int ncluster;             // # of independent clusters of coupled segments
  int njoint;               // # of interbody joints
  int **jointsegment;          // indices of 2 rigid segments in each joint (1-N)
  double **xjoint;          // coords of each joint point, NOT updated!
  int nfree;                // # of isolated unconnected segments
  int *freelist;            // indices of isolated segments (1-N)

  // POEMS object

  class Workspace *poems;

  // internal class functions

  void readfile(char *);
  int model_switch_flag;
  int readline(FILE *, char **, int *);
  void jointbuild();
  void sortlist(int, int **);
  int loopcheck(int, int, int **);
  int jacobi(double **, double *, double **);
  void rotate(double **, int, int, int, int, double, double);
  void omega_from_mq(double *, double *, double *, double *,
                     double *, double *);
  void set_v();
  void set_xv();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Could not find fix poems group ID

A group ID used in the fix poems command does not exist.

E: Must use a molecular atom style with fix poems molecule

Self-explanatory.

E: No rigid segments defined

The fix specification did not end up defining any rigid segments.

E: Atom in too many rigid segments - boost MAXSEGMENT

Fix poems has a parameter MAXSEGMENT (in fix_poems.cpp) which determines
the maximum number of rigid segments a single atom can belong to (i.e. a
multibody joint).  The segments you have defined exceed this limit.

E: One or zero atoms in rigid segment

Any rigid segment defined by the fix rigid command must contain 2 or more
atoms.

W: More than one fix poems

It is not efficient to use fix poems more than once.

E: POEMS fix must come before NPT/NPH fix

NPT/NPH fix must be defined in input script after all poems fixes,
else the fix contribution to the pressure virial is incorrect.

E: Insufficient Jacobi rotations for POEMS segment

Eigensolve for rigid segment was not sufficiently accurate.

E: Rigid segment has degenerate moment of inertia

Fix poems will only work with segments (collections of atoms) that have
non-zero principal moments of inertia.  This means they must be 3 or
more non-collinear atoms, even with joint atoms removed.

E: Bad principal moments

Fix rigid did not compute the principal moments of inertia of a rigid
group of atoms correctly.

E: Cannot open fix poems file %s

The specified file cannot be opened.  Check that the path and name are
correct.

W: No joints between rigid segments, use fix rigid instead

The segments defined by fix poems are not connected by joints.  POEMS
will integrate the segment motion, but it would be more efficient to use
fix rigid.

E: Cyclic loop in joint connections

Fix poems cannot (yet) work with coupled segments whose joints connect
the segments in a ring (or cycle).

E: Tree structure in joint connections

Fix poems cannot (yet) work with coupled segments whose joints connect
the segments in a tree structure.

*/
