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

/* ----------------------------------------------------------------------
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Philippe Seil (JKU Linz)
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
    mesh_(0),
    move_(0),
    neighListFresh_(false),
    time_(0),
    time_since_setup_(0)
{
    if(narg < 6)
      error->all(FLERR,"Illegal fix move/mesh command, you need to specify a mesh");

    int iarg = 3;

    if(strcmp(arg[iarg++],"mesh"))
      error->all(FLERR,"Illegal fix move/mesh command, expecting keyword 'mesh'");

    fix_mesh_ = static_cast<FixMesh*>(modify->find_fix_id(arg[iarg++]));
    if(fix_mesh_ == 0)
        error->all(FLERR,"Illegal fix move/mesh command, illegal mesh ID provided");

    mesh_ = fix_mesh_->mesh();
    move_ = createMeshMover(lmp,mesh_,&arg[iarg],narg-iarg);

    if(move_ == 0)
      error->all(FLERR,"Illegal fix move/mesh command, illegal arguments");

    // not compatible because surface velocity is just set once
    // and for moving mesh it is set every step
    if(fix_mesh_->surfaceVel())
      error->all(FLERR,"Illegal fix move/mesh command, cannot apply move to a mesh using keywords 'velocity' or 'angular_velocity'");

    restart_global = 1;

    force_reneighbor = 1;

    next_reneighbor = -1;
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

        if(nmove > 1 && move_->isFirst())
        error->all(FLERR,"Illegal deletion of a fix move/mesh. There is another fix move/mesh command active on the same mesh. "
                           "Superposed fix move/mesh commands must be unfixed in reverse order of creation");

        move_->pre_delete();

        // do not delete property v, as a dump command may still refer to it
        // set velocity to zero
        
        MultiVectorContainer<double,3,3> *v;
        v = mesh_->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v");
        if(v) v->setAll(0.);
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
    mask |= PRE_FORCE;
    mask |= FINAL_INTEGRATE;
    mask |= PRE_NEIGHBOR;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixMoveMesh::setup_pre_force(int vflag)
{
    time_since_setup_ = 0.;

    pre_neighbor();
    pre_force(0);
}

/* ---------------------------------------------------------------------- */

void FixMoveMesh::setup(int vflag)
{
    if(!mesh_->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v"))
    {
        mesh_->prop().addElementProperty<MultiVectorContainer<double,3,3> >("v","comm_none","frame_general","restart_no");
        
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
    move_->initial_integrate(time_since_setup_,dt);
}

/* ---------------------------------------------------------------------- */

void FixMoveMesh::pre_neighbor()
{
    neighListFresh_ = true;
}

/* ---------------------------------------------------------------------- */

void FixMoveMesh::pre_force(int vflag)
{
    // only first move on this mesh handles neigh list build
    if(!move_->isFirst()) return;

    // neigh list built on this step
    
    if(neighListFresh_)
    {
        store_node_pos();

    }
    // check for re-build of neigh list
    else
    {
        if(decide_rebuild())
            next_reneighbor = update->ntimestep + 1;
    }
}

/* ---------------------------------------------------------------------- */

void FixMoveMesh::final_integrate()
{
    double dt = update->dt;

    // useful only if accelerations are known
    move_->final_integrate(time_,dt);
}
/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixMoveMesh::write_restart(FILE *fp)
{
  int n = 0;
  double list[1];
  list[n++] = time_;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixMoveMesh::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  time_ = static_cast<int> (list[n++]);
}

/* ----------------------------------------------------------------------
   decide if any node has moved far enough to trigger re-build
------------------------------------------------------------------------- */

bool FixMoveMesh::decide_rebuild()
{
    
    double ***node = mesh_->nodePtr();
    double ***old = oldNodes.begin();
    int flag = 0;
    int nlocal = mesh_->sizeLocal();
    double triggersq = 0.25*neighbor->skin*neighbor->skin;

    for(int iTri = 0; iTri < nlocal; iTri++)
    {
      for(int iNode = 0; iNode < 3; iNode++)
      {
        double deltaX[3];
        vectorSubtract3D(node[iTri][iNode],old[iTri][iNode],deltaX);
        double distSq = deltaX[0]*deltaX[0] + deltaX[1]*deltaX[1] + deltaX[2]*deltaX[2];
        if(distSq > triggersq){
          
          flag = 1;
        }
      }
      if (flag) break;
    }

    // allreduce result
   MPI_Max_Scalar(flag,this->world);

    if(flag) return true;
    else     return false;
}

/* ----------------------------------------------------------------------
   store node pos at lat re-build
------------------------------------------------------------------------- */

void FixMoveMesh::store_node_pos()
{
    int nlocal = mesh_->sizeLocal();
    double ***node = mesh_->nodePtr();

    oldNodes.empty();
    for(int i = 0; i < nlocal; i++)
        oldNodes.add(node[i]);
}
