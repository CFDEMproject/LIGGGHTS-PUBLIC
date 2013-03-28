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
#include "fix_cfd_coupling_force_implicit.h"
#include "fix_property_atom.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCfdCouplingForceImplicit::FixCfdCouplingForceImplicit(LAMMPS *lmp, int narg, char **arg) :
    FixCfdCouplingForce(lmp,narg,arg),
    fix_Ksl_(0),
    fix_uf_(0)
{

}

/* ---------------------------------------------------------------------- */

FixCfdCouplingForceImplicit::~FixCfdCouplingForceImplicit()
{

}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForceImplicit::post_create()
{
    // do mother class init w/o dragforce
    FixCfdCouplingForce::post_create();

    // register Ksl
    if(!fix_Ksl_)
    {
        char* fixarg[9];
        fixarg[0]="Ksl";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="Ksl";
        fixarg[4]="scalar"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fix_Ksl_ = modify->add_fix_property_atom(9,fixarg,style);
    }

    // register uf
    if(!fix_uf_)
    {
        char* fixarg[11];
        fixarg[0]="uf";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="uf";
        fixarg[4]="vector"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fixarg[9]="0.";
        fixarg[10]="0.";
        fix_uf_ = modify->add_fix_property_atom(11,fixarg,style);
    }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForceImplicit::pre_delete(bool unfixflag)
{
    if(unfixflag && fix_Ksl_) modify->delete_fix("Ksl");
    if(unfixflag && fix_uf_) modify->delete_fix("uf");
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForceImplicit::init()
{
    FixCfdCouplingForce::init();

    // values to come from OF
    fix_coupling_->add_pull_property("Ksl","scalar-atom");
    fix_coupling_->add_pull_property("uf","vector-atom");
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForceImplicit::post_force(int vflag)
{
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *Ksl = fix_Ksl_->vector_atom;
  double **uf = fix_uf_->array_atom;
  double **dragforce = fix_dragforce_->array_atom;
  double frc[3];

  vectorZeroize3D(dragforce_total);

  // add dragforce to force vector
  for (int i = 0; i < nlocal; i++)
  {
    if (mask[i] & groupbit)
    {
        // calc force
        vectorSubtract3D(uf[i],v[i],frc);
        vectorScalarMult3D(frc,Ksl[i]);

        // add force
        vectorAdd3D(f[i],frc,f[i]);
        vectorAdd3D(f[i],dragforce[i],f[i]);

        // add up forces for post-proc
        vectorAdd3D(dragforce_total,frc,dragforce_total);
        vectorAdd3D(dragforce_total,dragforce[i],dragforce_total);
    }
  }
}
