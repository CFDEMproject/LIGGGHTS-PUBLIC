/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
   the producer of the LIGGGHTS® software and the CFDEM®coupling software
   See http://www.cfdem.com/terms-trademark-policy for details.

   LIGGGHTS® is based on LAMMPS
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

FixCfdCouplingForce::FixCfdCouplingForce(LAMMPS *lmp, int narg, char **arg) : Fix(lmp,narg,arg),
    fix_coupling_(0),
    fix_dragforce_(0),
    fix_hdtorque_(0),
    fix_volumeweight_(0),
    use_force_(true),
    use_torque_(true),
    use_dens_(false),
    use_type_(false),
    use_property_(false)
{
    int iarg = 3;

    bool hasargs = true;
    while(iarg < narg && hasargs)
    {
        hasargs = false;

        if(strcmp(arg[iarg],"transfer_density") == 0) {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'transfer_density'");
            iarg++;
            if(strcmp(arg[iarg],"yes") == 0)
                use_dens_ = true;
            else if(strcmp(arg[iarg],"no") == 0)
                use_dens_ = false;
            else
                error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'transfer_density'");
            iarg++;
            hasargs = true;
        }
        else if(strcmp(arg[iarg],"transfer_torque") == 0)
        {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'transfer_torque'");
            iarg++;
            if(strcmp(arg[iarg],"yes") == 0)
                use_torque_ = true;
            else if(strcmp(arg[iarg],"no") == 0)
                use_torque_ = false;
            else
                error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'transfer_torque'");
            iarg++;
            hasargs = true;
        }
        else if(strcmp(arg[iarg],"transfer_type") == 0)
        {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'transfer_type'");
            iarg++;
            if(strcmp(arg[iarg],"yes") == 0)
                use_type_ = true;
            else if(strcmp(arg[iarg],"no") == 0)
                use_type_ = false;
            else
                error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'transfer_type'");
            iarg++;
            hasargs = true;
        } else if(strcmp(arg[iarg],"transfer_property") == 0) {
            if(narg < iarg+5)
                error->fix_error(FLERR,this,"not enough arguments for 'transfer_type'");
            iarg++;
            use_property_ = true;
            if(strcmp(arg[iarg++],"name"))
                error->fix_error(FLERR,this,"expecting 'name' after 'transfer_property'");
            sprintf(property_name,"%s",arg[iarg++]);
            if(strcmp(arg[iarg++],"type"))
                error->fix_error(FLERR,this,"expecting 'type' after property name");
            sprintf(property_type,"%s",arg[iarg++]);
            iarg++;
            hasargs = true;
        } else if (strcmp(this->style,"couple/cfd/force") == 0) {
            error->fix_error(FLERR,this,"unknown keyword");
        }
    }

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
    if(!fix_dragforce_ && use_force_)
    {
        const char* fixarg[11];
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
        fix_dragforce_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
    }

    // register hydrodynamic torque
    if(!fix_hdtorque_ && use_torque_)
    {
        const char* fixarg[11];
        fixarg[0]="hdtorque";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="hdtorque";
        fixarg[4]="vector"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fixarg[9]="0.";
        fixarg[10]="0.";
        fix_hdtorque_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
    }

    // register volume weight for volume fraction calculation if not present
    // is 1 per default
    fix_volumeweight_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("volumeweight","property/atom","scalar",0,0,style,false));
    if(!fix_volumeweight_)
    {
        const char* fixarg[9];
        fixarg[0]="volumeweight";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="volumeweight";
        fixarg[4]="scalar"; // 1 vector per particle to be registered
        fixarg[5]="no";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="1.";
        fix_volumeweight_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForce::pre_delete(bool unfixflag)
{
    if(unfixflag && fix_dragforce_) modify->delete_fix("dragforce");
    if(unfixflag && fix_hdtorque_) modify->delete_fix("hdtorque");
    if(unfixflag && fix_volumeweight_) modify->delete_fix("volumeweight");
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
      error->fix_error(FLERR,this,"More than one fix of this style is not allowed");

    // find coupling fix
    fix_coupling_ = static_cast<FixCfdCoupling*>(modify->find_fix_style_strict("couple/cfd",0));
    if(!fix_coupling_)
      error->fix_error(FLERR,this,"Fix couple/cfd/force needs a fix of type couple/cfd");

    //  values to be transfered to OF

    fix_coupling_->add_push_property("x","vector-atom");
    fix_coupling_->add_push_property("v","vector-atom");
    fix_coupling_->add_push_property("radius","scalar-atom");
    if(use_type_) fix_coupling_->add_push_property("type","scalar-atom");
    if(use_dens_) fix_coupling_->add_push_property("density","scalar-atom");
    if(use_torque_) fix_coupling_->add_push_property("omega","vector-atom");
    fix_coupling_->add_push_property("volumeweight","scalar-atom");

    if(use_property_) fix_coupling_->add_push_property(property_name,property_type);

    // values to come from OF
    if(use_force_) fix_coupling_->add_pull_property("dragforce","vector-atom");
    if(use_torque_) fix_coupling_->add_pull_property("hdtorque","vector-atom");

    vectorZeroize3D(dragforce_total);
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForce::post_force(int)
{
  double **f = atom->f;
  double **torque = atom->torque;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double **dragforce = fix_dragforce_->array_atom;
  double **hdtorque = fix_hdtorque_->array_atom;

  vectorZeroize3D(dragforce_total);

  // add dragforce to force vector
  
  for (int i = 0; i < nlocal; i++)
  {
    if (mask[i] & groupbit)
    {
        if(use_force_) vectorAdd3D(f[i],dragforce[i],f[i]);
        if(use_torque_) vectorAdd3D(torque[i],hdtorque[i],torque[i]);
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
