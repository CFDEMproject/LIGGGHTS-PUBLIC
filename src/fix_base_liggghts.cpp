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
    (if no contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2015-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#include <cmath>
#include <string.h>
#include <stdlib.h>
#include "fix_base_liggghts.h"
#include "atom.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "update.h"
#include "respa.h"
#include "region.h"
#include "domain.h"
#include "modify.h"
#include "atom_vec.h"
#include "comm.h"
#include "fix_multisphere.h"
#include "mpi_liggghts.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBaseLiggghts::FixBaseLiggghts(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg),
  support_respa_(false),
  nlevels_respa_(0),
  do_need_radius_(true),
  do_need_mass_(true),
  support_ms_(false),
  fix_ms_(0),
  ms_(0),
  idregion_(0),
  region_(0),
  iregion_(-1)
{
}

/* ---------------------------------------------------------------------- */

FixBaseLiggghts::~FixBaseLiggghts()
{
    if(idregion_) delete []idregion_;
}

/* ---------------------------------------------------------------------- */

void FixBaseLiggghts::process_region(char *regid)
{
    iregion_ = domain->find_region(regid);
    if (iregion_ == -1)
        error->fix_error(FLERR,this,"Region ID does not exist");
    int n = strlen(regid) + 1;
    idregion_ = new char[n];
    strcpy(idregion_,regid);
    region_ = domain->regions[iregion_];
}

/* ---------------------------------------------------------------------- */

void FixBaseLiggghts::init()
{
    // error checks
    if (do_need_radius_ && !atom->radius_flag)
        error->fix_error(FLERR,this,"requires atom attribute radius (per-particle)");

    // error checks
    if (do_need_mass_ && !atom->rmass_flag)
        error->fix_error(FLERR,this,"requires atom attribute mass (per-particle)");

    // check validity of region
    iregion_ = -1;
    if (idregion_)
    {
      iregion_ = domain->find_region(idregion_);
      if (iregion_ == -1)
        error->fix_error(FLERR,this,"Region ID does not exist");
      region_ = domain->regions[iregion_];
    }

    // multisphere support

    fix_ms_ = static_cast<FixMultisphere*>(modify->find_fix_style("multisphere",0));
    if(modify->n_fixes_style("multisphere") > 1)
      error->fix_error(FLERR,this,"does not support more than one fix multisphere.");
    if(fix_ms_)
      ms_ = &(fix_ms_->data());
    else
      ms_ = 0;

    if(!support_ms_ && ms_)
        error->fix_error(FLERR,this,"does not support multi-sphere");

    // currently no groups for MS
    int igrp_all = group->find("all");
    if(fix_ms_ && (igrp_all!= igroup))
       error->fix_error(FLERR,this,"does only support fix group 'all' when multi-sphere particles present");

    if (strstr(update->integrate_style,"respa"))
        nlevels_respa_ = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixBaseLiggghts::setup(int vflag)
{
    if (strstr(update->integrate_style,"verlet"))
        post_force(vflag);
    else
    {
        if(!support_respa_)
            error->fix_error(FLERR,this,"does nor support run_style respa");
        ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa_-1);
        post_force_respa(vflag,nlevels_respa_-1,0);
        ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa_-1);
    }
}

/* ---------------------------------------------------------------------- */

void FixBaseLiggghts::count_eligible(double &mass_counted,double &volume_counted, int &nparticles_counted)
{
    int nlocal = atom->nlocal;
    int nbody_local = ms_? ms_->n_body() : 0;
    int *mask = atom->mask;
    double **x = atom->x;
    double *rmass = atom->rmass;
    double *radius = atom->radius;
    mass_counted = 0.;
    volume_counted = 0.;
    nparticles_counted = 0.;

    double _4_pi_over_3 = 4.*M_PI/3.;

    // spheres contribution

    for(int ilocal = 0; ilocal < nlocal; ilocal++)
    {
        if (mask[ilocal] & groupbit)
        {
            // do not handle multi-sphere case
            // skip if not in region
            if (
                 (!fix_ms_ || fix_ms_->belongs_to(ilocal) < 0) &&
                 (!region_ || region_->match(x[ilocal][0],x[ilocal][1],x[ilocal][2]))
               )
            {
                mass_counted += rmass[ilocal];
                volume_counted += _4_pi_over_3*radius[ilocal]*radius[ilocal]*radius[ilocal];
                nparticles_counted += 1;
            }
        }
    }

    // multi-spheres contribution

    for(int ibody_local = 0; ibody_local < nbody_local ; ibody_local++)
    {
        double xcm[3];
        ms_->xcm(xcm,ibody_local);

        // skip if not in region; group not accounted for here
        if (!region_ || region_->match(xcm[0],xcm[1],xcm[2]))
        {
            mass_counted += ms_->mass(ibody_local);
            volume_counted += ms_->volume(ibody_local);
            nparticles_counted += 1;
        }
    }

    double vec_mpi[3];
    vec_mpi[0] = mass_counted;
    vec_mpi[1] = volume_counted;
    vec_mpi[2] = static_cast<double>(nparticles_counted);
    MPI_Sum_Vector(vec_mpi,3,world);
    mass_counted = vec_mpi[0];
    volume_counted = vec_mpi[1];
    nparticles_counted = static_cast<int>(vec_mpi[2]);
}
