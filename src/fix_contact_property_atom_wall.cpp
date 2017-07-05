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

    Arno Mayrhofer (CFDEMresearch GmbH, Linz)
    Christoph Kloss (DCS Computing GmbH, Linz)

    Copyright 2016-     CFDEMresearch GmbH, Linz
    Copyright 2014-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#include <string.h>
#include <stdio.h>
#include "fix_contact_property_atom_wall.h"
#include "atom.h"
#include "neighbor.h"
#include "fix_wall_gran.h"
#include "fix_mesh_surface.h"
#include "pair_gran.h"
#include "force.h"
#include "primitive_wall.h"
#include "pair.h"
#include "modify.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixContactPropertyAtomWall::FixContactPropertyAtomWall(LAMMPS *lmp, int narg, char **arg) :
  FixContactPropertyAtom(lmp, narg, arg),
  fix_mesh_surface_(0),
  fix_nneighs_(0),
  primitive_wall_(0)
{
    // read style: primitive or mesh
    if(narg-iarg_ < 2)
     error->fix_error(FLERR,this,"not enough parameters");

    // parse further args
    if (strcmp(arg[iarg_],"primitive") == 0)
    {
        iarg_++;
        FixWallGran *fwg = static_cast<FixWallGran*>(modify->find_fix_id(arg[iarg_++]));
        if(!fwg)
            error->fix_error(FLERR,this,"illegal FixWallGran id");
        primitive_wall_ = fwg->primitiveWall();
    }
    else if (strcmp(arg[iarg_],"mesh") == 0)
    {
        iarg_++;
        fix_mesh_surface_ = static_cast<FixMeshSurface*>(modify->find_fix_id(arg[iarg_++]));
        if(!fix_mesh_surface_)
            error->fix_error(FLERR,this,"illegal FixMeshSurface id");

        fix_nneighs_ = fix_mesh_surface_->meshNeighlist()->fix_nneighs();
    }
    else
        error->fix_error(FLERR,this,"expecting 'primitive' or 'mesh'");
}

/* ---------------------------------------------------------------------- */

FixContactPropertyAtomWall::~FixContactPropertyAtomWall()
{

}

/* ----------------------------------------------------------------------
   allocate storage
------------------------------------------------------------------------- */

void FixContactPropertyAtomWall::clear()
{
    int nall = atom->nlocal+atom->nghost;

    // reset number of partners every time-step
    vectorZeroizeN(npartner_,nall);

    // other stuff to do only upon neigh list rebuild
    if(!build_neighlist_)
        return;
    build_neighlist_ = false;

    int nneighs_next;
    double half_skin = neighbor->skin * 0.5;

    ipage_->reset();
    dpage_->reset();

    // allocate for owned and ghost
    for (int i = 0; i < nall; i++)
    {
        // just one possible contact with primitive wall
        nneighs_next = fix_nneighs_ ? fix_nneighs_->get_vector_atom_int(i) : primitive_wall_->isNear(i,half_skin);

        npartner_[i] = 0;

        partner_[i] = ipage_->get(nneighs_next);
        vectorInitializeN(partner_[i],nneighs_next,-1);
        contacthistory_[i] = dpage_->get(nneighs_next*dnum_);
        vectorZeroizeN(contacthistory_[i],nneighs_next*dnum_);

   }
}

/* ---------------------------------------------------------------------- */

bool FixContactPropertyAtomWall::haveContact(const int iP, const int idTri, double *&history) const
{
    int *tri = partner_[iP];
    const double half_skin = neighbor->skin * 0.5;
    const int nneighs = fix_nneighs_ ? fix_nneighs_->get_vector_atom_int(iP) : primitive_wall_->isNear(iP, half_skin);

    for(int i = 0; i < nneighs; i++)
    {
        if(tri[i] == idTri)
        {
            if(dnum_ > 0) history = &(contacthistory_[iP][i*dnum_]);
            return true;
        }
    }
    return false;
}

FixMeshSurface *FixContactPropertyAtomWall::getMesh() const
{
    return fix_mesh_surface_;
}

