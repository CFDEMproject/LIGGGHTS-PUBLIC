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
    Arno Mayrhofer (DCS Computing GmbH, Linz)

    Copyright 2017-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#include "compute_erotate.h"
#include "error.h"
#include "fix.h"
#include "modify.h"
#include "atom.h"
#include "update.h"
#include <cstring>

using namespace LAMMPS_NS;

ComputeERotate::ComputeERotate(LAMMPS * lmp, int &iarg, int narg, char **arg) :
    Compute(lmp, iarg, narg, arg),
    csphere_(NULL),
    cms_(NULL),
#ifdef SUPERQUADRIC_ACTIVE_FLAG
    csq_(NULL),
#endif
    id_group_(NULL)
{
    if (narg < iarg)
        error->all(FLERR, "Illegal compute erotate command");

    scalar_flag = 1;
    extscalar = 1;

    bool found = true;
    for (int i = 0; found; i++)
    {
        csphere_ = static_cast<ComputeERotateSphere *>(modify->find_compute_style_strict("erotate/sphere", i));
        if (csphere_)
        {
            found = true;
            if (csphere_->igroup != igroup)
                csphere_ = NULL;
        }
        else
            found = false;
    }

    found = true;
    for (int i = 0; found; i++)
    {
        cms_ = static_cast<ComputeERotateMultisphere *>(modify->find_compute_style_strict("erotate/multisphere", i));
        if (cms_)
        {
            found = true;
            if (cms_->igroup != igroup)
                cms_ = NULL;
        }
        else
            found = false;
    }

#ifdef SUPERQUADRIC_ACTIVE_FLAG
    found = true;
    for (int i = 0; found; i++)
    {
        csq_ = static_cast<ComputeERotateSuperquadric *>(modify->find_compute_style_strict("erotate/superquadric", i));
        if (csq_)
        {
            found = true;
            if (csq_->igroup != igroup)
                csq_ = NULL;
        }
        else
            found = false;
    }
#endif

    const int len = strlen(arg[1])+1;
    id_group_ = new char[len];
    strncpy(id_group_, arg[1], len);
    printf("igrpu: %s\n", id_group_);
}

void ComputeERotate::post_create()
{
    if (!csphere_ && atom->sphere_flag)
    {
        const int cnargs = update_on_run_end_ ? 5 : 3;
        char **carg = new char * [cnargs];
        carg[0] = (char *) "erotate_sphere_";
        carg[1] = id_group_;
        carg[2] = (char *) "erotate/sphere";
        if (update_on_run_end_)
        {
            carg[3] = (char *) "update_on_run_end";
            carg[4] = (char *) "yes";
        }
        modify->add_compute(cnargs, carg);
        csphere_ = static_cast<ComputeERotateSphere *>(modify->compute[modify->find_compute("erotate_sphere_")]);
        delete [] carg;
    }

    Fix * fix_ms =  modify->find_fix_style("multisphere",0);
    if (!cms_ && fix_ms)
    {
        const int cnargs = update_on_run_end_ ? 5 : 3;
        char **carg = new char * [cnargs];
        carg[0] = (char *) "erotate_multisphere_";
        carg[1] = id_group_;
        carg[2] = (char *) "erotate/multisphere";
        if (update_on_run_end_)
        {
            carg[3] = (char *) "update_on_run_end";
            carg[4] = (char *) "yes";
        }
        modify->add_compute(cnargs, carg);
        cms_ = static_cast<ComputeERotateMultisphere *>(modify->compute[modify->find_compute("erotate_multisphere_")]);
        delete [] carg;
    }

#ifdef SUPERQUADRIC_ACTIVE_FLAG
    if (!csq_ && atom->superquadric_flag)
    {
        const int cnargs = update_on_run_end_ ? 5 : 3;
        char **carg = new char * [cnargs];
        carg[0] = (char *) "erotate_superquadric_";
        carg[1] = id_group_;
        carg[2] = (char *) "erotate/superquadric";
        if (update_on_run_end_)
        {
            carg[3] = (char *) "update_on_run_end";
            carg[4] = (char *) "yes";
        }
        modify->add_compute(cnargs, carg);
        csq_ = static_cast<ComputeERotateSuperquadric *>(modify->compute[modify->find_compute("erotate_superquadric_")]);
        delete [] carg;
    }
#endif

}

void ComputeERotate::init()
{}

double ComputeERotate::compute_scalar()
{
    if (invoked_scalar == update->ntimestep)
        return scalar;
    invoked_scalar = update->ntimestep;

    scalar = 0.0;
    if (csphere_)
        scalar += csphere_->compute_scalar();

    if (cms_)
        scalar += cms_->compute_scalar();

#ifdef SUPERQUADRIC_ACTIVE_FLAG
    if (csq_)
        scalar += csq_->compute_scalar();
#endif

    return scalar;
}
