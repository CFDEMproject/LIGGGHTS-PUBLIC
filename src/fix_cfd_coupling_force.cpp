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

#include <string.h>
#include <stdlib.h>
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "comm.h"
#include <cmath>
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
    fix_dispersionTime_(0),
    fix_dispersionVel_(0),
    fix_UrelOld_(0),
    use_force_(true),
    use_torque_(true),
    use_dens_(false),
    use_type_(false),
    use_stochastic_(false),
    use_virtualMass_(false),
    use_superquadric_(false),
    use_id_(false),
    use_property_(false),
    use_fiber_topo_(false),
    fix_fiber_axis_(0),
    fix_fiber_ends_(0)
{
    iarg = 3;

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
        }
        else if(strcmp(arg[iarg],"transfer_id") == 0)
        {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'transfer_type'");
            iarg++;
            if(strcmp(arg[iarg],"yes") == 0)
                use_id_ = true;
            else if(strcmp(arg[iarg],"no") == 0)
                use_id_ = false;
            else
                error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'transfer_id'");
            iarg++;
            hasargs = true;
        }
        else if(strcmp(arg[iarg],"transfer_stochastic") == 0)
        {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'transfer_stochastic'");
            iarg++;
            if(strcmp(arg[iarg],"yes") == 0)
                use_stochastic_ = true;
            else if(strcmp(arg[iarg],"no") == 0)
                use_stochastic_ = false;
            else
                error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'transfer_stochastic'");
            iarg++;
            hasargs = true;
        }
        else if(strcmp(arg[iarg],"transfer_virtualMass") == 0)
        {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'transfer_virtualMass'");
            iarg++;
            if(strcmp(arg[iarg],"yes") == 0)
                use_virtualMass_ = true;
            else if(strcmp(arg[iarg],"no") == 0)
                use_virtualMass_ = false;
            else
                error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'transfer_virtualMass'");
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
        } else if(strcmp(arg[iarg],"transfer_fiber_topology") == 0) {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'transfer_fiber_topology'");
            if(strcmp(arg[iarg],"yes") == 0)
                use_fiber_topo_ = true;
            else if(strcmp(arg[iarg],"no") == 0)
                use_fiber_topo_ = false;
            else
                error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'transfer_fiber_topology'");
            iarg++;
            hasargs = true;
        } else if(strcmp(arg[iarg],"transfer_superquadric") == 0) {
            if(narg < iarg+2)
              error->fix_error(FLERR,this,"not enough arguments for 'transfer_superquadric'");
            iarg++;
            if(strcmp(arg[iarg],"yes") == 0) {
              use_superquadric_ = true;
              use_torque_ = true;
              use_force_ = true;
            }
            else if(strcmp(arg[iarg],"no") == 0) {
              use_superquadric_ = false;
            }
            else
              error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'transfer_superquadric'");
            iarg++;
            hasargs = true;
        } else if (strcmp(this->style,"couple/cfd/force") == 0) {
            error->fix_error(FLERR,this,"unknown keyword");
        }
    }

    // flags for vector output
    vector_flag = 1;
    size_vector = 6;
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

    if(!fix_dispersionTime_ && use_stochastic_)
    {
        const char* fixarg[9];
        fixarg[0]="dispersionTime";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="dispersionTime";
        fixarg[4]="scalar"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="1e12";
        fix_dispersionTime_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    if(!fix_dispersionVel_ && use_stochastic_)
    {
        const char* fixarg[11];
        fixarg[0]="dispersionVel";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="dispersionVel";
        fixarg[4]="vector"; // vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0";
        fixarg[9]="0";
        fixarg[10]="0";
        fix_dispersionVel_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
    }

    if(!fix_UrelOld_ && use_virtualMass_)
    {
        const char* fixarg[11];
        fixarg[0]="UrelOld";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="UrelOld";
        fixarg[4]="vector"; // vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0";
        fixarg[9]="0";
        fixarg[10]="0";
        fix_dispersionVel_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
    }

    if(use_fiber_topo_)
    {
        const char *fixarg[] = {
              "topo",       // fix id
              "all",        // fix group
              "bond/fiber/topology" // fix style
        };
        modify->add_fix(3,const_cast<char**>(fixarg));
    }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForce::pre_delete(bool unfixflag)
{
    if(unfixflag && fix_dragforce_) modify->delete_fix("dragforce");
    if(unfixflag && fix_hdtorque_) modify->delete_fix("hdtorque");
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
    if(use_superquadric_) {
      fix_coupling_->add_push_property("volume","scalar-atom");
      fix_coupling_->add_push_property("area","scalar-atom");
      fix_coupling_->add_push_property("shape","vector-atom");
      fix_coupling_->add_push_property("blockiness","vector2D-atom");
      fix_coupling_->add_push_property("quaternion","quaternion-atom");
    }
    if(use_type_) fix_coupling_->add_push_property("type","scalar-atom");
    if(use_dens_) fix_coupling_->add_push_property("density","scalar-atom");
    if(use_torque_) fix_coupling_->add_push_property("omega","vector-atom");
    if(use_id_) fix_coupling_->add_push_property("id","scalar-atom");

    if(use_property_) fix_coupling_->add_push_property(property_name,property_type);

    // values to come from OF
    if(use_force_) fix_coupling_->add_pull_property("dragforce","vector-atom");
    if(use_torque_) fix_coupling_->add_pull_property("hdtorque","vector-atom");

    if(use_stochastic_)
    {
        fix_coupling_->add_pull_property("dispersionTime","scalar-atom");
        fix_coupling_->add_pull_property("dispersionVel","vector-atom");
    }

    if(use_fiber_topo_)
    {
        fix_coupling_->add_pull_property("fiber_axis","vector-atom");
        fix_coupling_->add_pull_property("fiber_ends","vector-atom");
    }

    vectorZeroize3D(dragforce_total);
    vectorZeroize3D(hdtorque_total);

    if (strcmp(update->integrate_style,"respa") == 0)
       error->fix_error(FLERR,this,"'run_style respa' not supported.");

}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForce::setup(int vflag)
{
    if (strstr(update->integrate_style,"verlet"))
        post_force(vflag);
    else
        error->fix_error(FLERR,this,"only 'run_style verlet' supported.");
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
  vectorZeroize3D(hdtorque_total);

  // add dragforce to force vector
  
  for (int i = 0; i < nlocal; i++)
  {
    if (mask[i] & groupbit)
    {
        if(use_force_)
        {
            vectorAdd3D(f[i],dragforce[i],f[i]);
            vectorAdd3D(dragforce_total,dragforce[i],dragforce_total);
        }
        if(use_torque_)
        {
            vectorAdd3D(torque[i],hdtorque[i],torque[i]);
            vectorAdd3D(hdtorque_total,hdtorque[i],hdtorque_total);
        }
    }
  }
}

/* ----------------------------------------------------------------------
   return components of total force on fix group
------------------------------------------------------------------------- */

double FixCfdCouplingForce::compute_vector(int n)
{
  if(n < 3)
  {
    double dragtotal = dragforce_total[n];
    MPI_Sum_Scalar(dragtotal,world);
    return dragtotal;
  }

  double hdtorque = hdtorque_total[n-3];
  MPI_Sum_Scalar(hdtorque,world);
  return hdtorque;
}
