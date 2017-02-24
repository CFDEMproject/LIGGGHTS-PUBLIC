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

    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Philippe Seil (JKU Linz)
    Richard Berger (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#include "fix_move_mesh.h"
#include "fix_mesh.h"
#include "tri_mesh.h"
#include "modify.h"
#include "error.h"
#include "force.h"
#include "update.h"
#include "mesh_mover.h"
#include "container.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMoveMesh::FixMoveMesh(LAMMPS *lmp, int narg, char **arg)
  : Fix(lmp,narg,arg),
    fix_mesh_(0),
    move_(0),
    mesh_(0),
    time_(0),
    time_since_setup_(0)
{
    vectorZeroize3D(reference_point_);

    if(narg < 6)
      error->all(FLERR,"Illegal fix move/mesh command, you need to specify a mesh");

    int iarg = 3;

    if(strcmp(arg[iarg++],"mesh"))
      error->all(FLERR,"Illegal fix move/mesh command, expecting keyword 'mesh'");

    fix_mesh_ = dynamic_cast<FixMesh*>(modify->find_fix_id(arg[iarg++]));
    if(fix_mesh_ == 0)
        error->all(FLERR,"Illegal fix move/mesh command, illegal mesh ID provided");

    mesh_ = fix_mesh_->mesh();
    move_ = createMeshMover(lmp,mesh_,this,&arg[iarg],narg-iarg);

    if(move_ == 0)
      error->all(FLERR,"Illegal fix move/mesh command, illegal arguments");

    // not compatible because surface velocity is just set once
    // and for moving mesh it is set every step
    if(fix_mesh_->surfaceVel())
      error->all(FLERR,"Illegal fix move/mesh command, cannot apply move to a mesh using keywords 'velocity' or 'angular_velocity'");

    restart_global = 1;
}

/* ---------------------------------------------------------------------- */

void FixMoveMesh:: post_create()
{
    
    move_->post_create();
}

/* ---------------------------------------------------------------------- */

void FixMoveMesh::pre_delete(bool unfixflag)
{
    // check if another fix move operates on the same mesh
    // which came after me in the imput script
    // if so, it is illegal to delete this command
    // without first deleting the other

    if(unfixflag)
    {
        int nmove = modify->n_fixes_style("move/mesh");

        for(int imove = 0; imove < nmove; imove++)
        {
            FixMoveMesh* fix_move_mesh = static_cast<FixMoveMesh*>(modify->find_fix_style("move/mesh",imove));
            if(fix_move_mesh != this && fix_move_mesh->fixMesh() == fixMesh() && move_->isFirst())
                error->all(FLERR,"Illegal deletion of a fix move/mesh. There is another fix move/mesh command active on the same mesh. "
                           "Superposed fix move/mesh commands must be unfixed in reverse order of creation");
        }

        move_->pre_delete();

        // do not delete property v, as a dump command may still refer to it
        // set velocity to zero
        
        MultiVectorContainer<double,3,3> *v;
        v = mesh_->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v");
        if(v) v->setAll(0.);

        // remove reference point if have one
        char refpt_id[200];
        sprintf(refpt_id, "REFPT_%s",id);

        if(mesh_->prop().getGlobalProperty<   VectorContainer<double,3> >(refpt_id))
           mesh_->prop().removeGlobalProperty(refpt_id);

    }

    delete move_;
}

/* ---------------------------------------------------------------------- */

FixMoveMesh::~FixMoveMesh()
{
    
}

/* ---------------------------------------------------------------------- */

int FixMoveMesh::setmask()
{
    int mask = 0;
    mask |= INITIAL_INTEGRATE;
    mask |= FINAL_INTEGRATE;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixMoveMesh::setup(int vflag)
{
    
    time_since_setup_ = 0.;
    
    reset_reference_point();

    if(!mesh_->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v"))
    {
        
        mesh_->prop().addElementProperty<MultiVectorContainer<double,3,3> >("v","comm_exchange_borders","frame_invariant","restart_no");
        
    }

    move_->setup();

    // do set-up in case velocity is dumped (called via set-up)
    if(move_->isFirst())
    {
        MultiVectorContainer<double,3,3> *v;
        v = mesh_->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v");
        v->setAll(0.);
    }
}

/* ---------------------------------------------------------------------- */

void FixMoveMesh::initial_integrate(int dummy)
{
    MultiVectorContainer<double,3,3> *v;

    double dt = update->dt;
    time_ += update->dt;
    time_since_setup_ += update->dt;

    if(move_->isFirst())
    {
        v = mesh_->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v");
        v->setAll(0.);
    }

    // integration
    
    move_->initial_integrate(time_,time_since_setup_,dt);
}

/* ---------------------------------------------------------------------- */

void FixMoveMesh::final_integrate()
{
    double dt = update->dt;

    // useful only if accelerations are known
    
    move_->final_integrate(time_,time_since_setup_,dt);
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixMoveMesh::write_restart(FILE *fp)
{
  int n = 0;
  double * list = new double[1 + move_->n_restart()];
  list[n++] = time_;
  
  move_->write_restart(&(list[n]));
  n += move_->n_restart();

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }

  delete[] list;
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixMoveMesh::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  time_ = static_cast<double> (list[n++]);
  move_->read_restart(&(list[n]));

}

/* ----------------------------------------------------------------------
   called by mesh mover
------------------------------------------------------------------------- */

void FixMoveMesh::add_reference_point(double *point)
{
    char refpt_id[200];
    sprintf(refpt_id, "REFPT_%s",id);

    if(mesh_->prop().getGlobalProperty<VectorContainer<double,3> >(refpt_id))
        error->fix_error(FLERR,this,"only one reference point allowed");

    vectorCopy3D(point,reference_point_);

    mesh_->prop().addGlobalProperty<VectorContainer<double,3> >(refpt_id,"comm_none","frame_general","restart_no");
    mesh_->prop().setGlobalProperty<VectorContainer<double,3> >(refpt_id,point);
}

/* ---------------------------------------------------------------------- */

void FixMoveMesh::get_reference_point(double *point)
{
    VectorContainer<double,3> *refpt;
    char refpt_id[200];

    sprintf(refpt_id, "REFPT_%s",id);
    refpt = mesh_->prop().getGlobalProperty<VectorContainer<double,3> >(refpt_id);

    if(!refpt)
        error->fix_error(FLERR,this,"internal error");

    if(move_->isFirst())
        mesh_->prop().resetGlobalPropToOrig(refpt_id);

    refpt->get(0,point);
    vectorCopy3D(point,reference_point_);
    
}

/* ---------------------------------------------------------------------- */

void FixMoveMesh::reset_reference_point()
{
    
    VectorContainer<double,3> *refpt;
    char refpt_id[200];

    sprintf(refpt_id, "REFPT_%s",id);
    refpt = mesh_->prop().getGlobalProperty<VectorContainer<double,3> >(refpt_id);

    // no error since not all moves have reference points
    if(!refpt)
        return;

    // set value for property
    refpt->set(0,reference_point_);

    // set orig value for property
    mesh_->prop().storeGlobalPropOrig(refpt_id);

}
