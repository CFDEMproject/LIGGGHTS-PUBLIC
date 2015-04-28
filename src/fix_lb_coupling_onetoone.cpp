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
    Philippe Seil (JKU Linz)

    Copyright 2004-     JKU Linz
------------------------------------------------------------------------- */

#include "fix_lb_coupling_onetoone.h"

#include "lammps.h"
#include "modify.h"
#include "error.h"
#include "atom.h"
#include "comm.h"

#include "fix_property_atom.h"

#include <iomanip>
#include <iostream>

namespace LAMMPS_NS {
  FixLbCouplingOnetoone::FixLbCouplingOnetoone(LAMMPS *lmp, int narg, char **arg)
    : Fix(lmp,narg,arg), fix_dragforce_(0), fix_hdtorque_(0)
  {
  }

  FixLbCouplingOnetoone::~FixLbCouplingOnetoone()
  {

  }

  int FixLbCouplingOnetoone::setmask()
  {
    int mask = 0;
    mask |= FixConst::POST_FORCE;
    return mask;
  }

  void FixLbCouplingOnetoone::init()
  {
    // make sure there is only one fix of this style
    if(modify->n_fixes_style(style) != 1)
      error->fix_error(FLERR,this,"More than one fix of this style is not allowed");

  }

  void FixLbCouplingOnetoone::post_create()
  {
    // register dragforce
    if(!fix_dragforce_)
      {
        const char *fixarg[] = {
          "dragforce", // fix id
          "all",       // fix group
          "property/atom", // fix style: property/atom
          "dragforce",     // property name
          "vector", // 1 vector per particle
          "yes",    // restart
          "no",     // communicate ghost
          "yes",    // communicate rev
          "0.","0.","0." // default values
        };
        fix_dragforce_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
      }

    // register hydrodynamic torque
    if(!fix_hdtorque_)
      {
        const char *fixarg[] = {
          "hdtorque", // fix id
          "all",      // fix group
          "property/atom", // fix style: property/atom
          "hdtorque",      // property name
          "vector", // 1 vector per particle
          "yes",    // restart
          "no",     // communicate ghost
          "yes",    // communicate rev
          "0.","0.","0."
        };
        fix_hdtorque_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
      }
  }

  void FixLbCouplingOnetoone::pre_delete(bool unfixflag)
  {
    if(!unfixflag) return;

    if(fix_dragforce_)
      modify->delete_fix(fix_dragforce_->id);
    if(fix_hdtorque_)
      modify->delete_fix(fix_hdtorque_->id);
  }

  void FixLbCouplingOnetoone::post_force(int)
  {
    double **f_ext = fix_dragforce_->array_atom;
    //double **t_ext = fix_hdtorque_->array_atom;
    double **t_ext = fix_hdtorque_->array_atom;
    double **f = atom->f;
    //double **t = atom->torque;

    // for(int i=0;i<atom->nlocal;i++)
    //   std::cout << comm->me << " force_liggghts "
    //             << std::setprecision(12) << f_ext[i][2] << std::endl;

    // std::cout << "before" << std::endl;
    // fix_dragforce_->do_reverse_comm();
    // fix_hdtorque_->do_reverse_comm();
    // std::cout << "after" << std::endl;

    // for(int i=0;i<atom->nlocal;i++)
    //   std::cout << comm->me << " force_liggghts_after "
    //             << std::setprecision(12) << f_ext[i][2] << " x " << atom->x[0][2] << std::endl;
    double **t = atom->torque;

    for(int i=0;i<atom->nlocal;i++){
      for(int j=0;j<3;j++){
        f[i][j] += f_ext[i][j];
        // t[i][j] += t_ext[i][j];
        t[i][j] += t_ext[i][j];
      }
    }

  }

  double **FixLbCouplingOnetoone::get_force_ptr()
  {
    return fix_dragforce_->array_atom;
  }
  double **FixLbCouplingOnetoone::get_torque_ptr()
  {
    return fix_hdtorque_->array_atom;
  }
  void FixLbCouplingOnetoone::comm_force_torque()
  {
    fix_dragforce_->do_reverse_comm();
    fix_hdtorque_->do_reverse_comm();
  }

}; /* LAMMPS_NS */
