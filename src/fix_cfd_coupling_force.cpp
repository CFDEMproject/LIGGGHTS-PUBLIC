/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include "string.h"
#include "stdlib.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "comm.h"
#include "math.h"
#include "vector_liggghts.h"
#include "mpi_liggghts.h"
#include "fix_cfd_coupling_force.h"
#include "fix_property_atom.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCfdCouplingForce::FixCfdCouplingForce(LAMMPS *lmp, int narg, char **arg) : Fix(lmp,narg,arg)
{
    fix_dragforce_ = NULL;
    fix_coupling_ = NULL;
    fix_volumeweight_ = NULL;

    // flags for vector output
    vector_flag = 1;
    size_vector = 3;
    global_freq = 1;
    extvector = 1;
}

/* ---------------------------------------------------------------------- */

FixCfdCouplingForce::~FixCfdCouplingForce()
{

}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForce::post_create()
{
    // register dragforce
    if(!fix_dragforce_)
    {
        char* fixarg[11];
        fixarg[0]="dragforce";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="dragforce";
        fixarg[4]="vector"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fixarg[9]="0.";
        fixarg[10]="0.";
        fix_dragforce_ = modify->add_fix_property_atom(11,fixarg,style);
    }

    // register volume weight for volume fraction calculation if not present
    // is 1 per default
    fix_volumeweight_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("volumeweight","property/atom","scalar",0,0,style,false));
    if(!fix_volumeweight_)
    {
        char* fixarg[9];
        fixarg[0]="volumeweight";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="volumeweight";
        fixarg[4]="scalar"; // 1 vector per particle to be registered
        fixarg[5]="no";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="1.";
        fix_volumeweight_ = modify->add_fix_property_atom(9,fixarg,style);
    }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForce::pre_delete(bool unfixflag)
{
    if(fix_dragforce_) modify->delete_fix("dragforce");
}

/* ---------------------------------------------------------------------- */

int FixCfdCouplingForce::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForce::init()
{
    // make sure there is only one fix of this style
    if(modify->n_fixes_style(style) != 1)
      error->fix_error(FLERR,this,"More than one fix of style couple/cfd/force is not allowed");

    // find coupling fix
    fix_coupling_ = static_cast<FixCfdCoupling*>(modify->find_fix_style_strict("couple/cfd",0));
    if(!fix_coupling_)
      error->fix_error(FLERR,this,"Fix couple/cfd/force needs a fix of type couple/cfd");

    //  values to be transfered to OF

    fix_coupling_->add_push_property("x","vector-atom");
    fix_coupling_->add_push_property("v","vector-atom");
    fix_coupling_->add_push_property("radius","scalar-atom");
    fix_coupling_->add_push_property("volumeweight","scalar-atom");

    // values to come from OF
    fix_coupling_->add_pull_property("dragforce","vector-atom");

    vectorZeroize3D(dragforce_total);
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForce::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double **dragforce = fix_dragforce_->array_atom;

  vectorZeroize3D(dragforce_total);

  // add dragforce to force vector
  for (int i = 0; i < nlocal; i++)
  {
    if (mask[i] & groupbit)
    {
        vectorAdd3D(f[i],dragforce[i],f[i]);
        vectorAdd3D(dragforce_total,dragforce[i],dragforce_total);
    }
  }
}

/* ----------------------------------------------------------------------
   return components of total force on fix group
------------------------------------------------------------------------- */

double FixCfdCouplingForce::compute_vector(int n)
{
  MPI_Sum_Vector(dragforce_total,3,world);
  return dragforce_total[n];
}
