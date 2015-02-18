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

/* ----------------------------------------------------------------------
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Philippe Seil (JKU Linz)
   Evan Smuts (U Cape Town, surface velocity rotation)
------------------------------------------------------------------------- */

#include "fix_mesh.h"
#include <stdio.h>
#include <string.h>
#include "error.h"
#include "force.h"
#include "bounding_box.h"
#include "input_mesh_tri.h"
#include "fix_contact_history.h"
#include "fix_neighlist_mesh.h"
#include "tri_mesh_deform.h"
#include "tri_mesh_planar.h"
#include "modify.h"
#include "comm.h"
#include "math_extra.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define EPSILON_V 0.00001

FixMesh::FixMesh(LAMMPS *lmp, int narg, char **arg)
: Fix(lmp, narg, arg),
  atom_type_mesh_(-1),
  mesh_(NULL),
  setupFlag_(false),
  pOpFlag_(false),
  manipulated_(false),
  verbose_(false),
  autoRemoveDuplicates_(false),
  precision_(0.)
{
    if(narg < 5)
      error->fix_error(FLERR,this,"not enough arguments - at least keyword 'file' and a filename are required.");

    restart_global = 1;

    force_reneighbor = 1;
    next_reneighbor = -1;

    // parse args

    iarg_ = 3;

    char mesh_fname[256];
    if(strcmp(arg[iarg_++],"file"))
        error->fix_error(FLERR,this,"expecting keyword 'file'");
    strcpy(mesh_fname,arg[iarg_++]);

    // parse args

    bool hasargs = true;
    while(iarg_ < narg && hasargs)
    {
        hasargs = false;
        if(strcmp(arg[iarg_],"type") == 0) {
            iarg_++;
            atom_type_mesh_ = force->inumeric(FLERR,arg[iarg_++]);
            if(atom_type_mesh_ < 1)
                error->fix_error(FLERR,this,"'type' > 0 required");
            hasargs = true;
        } else if(strcmp(arg[iarg_],"verbose") == 0) {
            if(narg < iarg_+2)
                error->fix_error(FLERR,this,"not enough arguments for 'verbose'");
            if(strcmp(arg[iarg_+1],"yes") == 0)
              verbose_ = true;
            else if(strcmp(arg[iarg_+1],"no"))
                error->fix_error(FLERR,this,"expecing 'yes' or 'no' for 'verbose'");
            iarg_ += 2;
            hasargs = true;
        } else if(strcmp(arg[iarg_],"heal") == 0) {
            if(narg < iarg_+2)
                error->fix_error(FLERR,this,"not enough arguments for 'heal'");
            if(strcmp(arg[iarg_+1],"auto_remove_duplicates") == 0)
                autoRemoveDuplicates_ = true;
            else if(strcmp(arg[iarg_+1],"no"))
                error->fix_error(FLERR,this,"expecing 'auto_remove_duplicates' or 'no' for 'heal'");
            iarg_ += 2;
            hasargs = true;
        } else if (strcmp(arg[iarg_],"precision") == 0) {
            if (narg < iarg_+2) error->fix_error(FLERR,this,"not enough arguments");
            iarg_++;
            precision_ = force->numeric(FLERR,arg[iarg_++]);
            if(precision_ < 0. || precision_ > 0.001)
              error->fix_error(FLERR,this,"0 < precision < 0.001 required");
            hasargs = true;
        }
    }

    // construct a mesh - can be surface or volume mesh
    // just create object and return if reading data from restart file
    
    if(modify->have_restart_data(this)) create_mesh_restart();
    else create_mesh(mesh_fname);

    // parse further args

    hasargs = true;
    while(iarg_ < narg && hasargs)
    {
      hasargs = false;
      if(strcmp(arg[iarg_],"move") == 0) {
          if (narg < iarg_+4) error->fix_error(FLERR,this,"not enough arguments");
          moveMesh(force->numeric(FLERR,arg[iarg_+1]),force->numeric(FLERR,arg[iarg_+2]),force->numeric(FLERR,arg[iarg_+3]));
          manipulated_ = true;
          iarg_ += 4;
          hasargs = true;
      } else if(strcmp(arg[iarg_],"rotate") == 0) {
          if (narg < iarg_+7)
              error->fix_error(FLERR,this,"not enough arguments");
          if(strcmp(arg[iarg_+1],"axis"))
              error->fix_error(FLERR,this,"expecting keyword 'axis' after keyword 'rotate'");
          if(strcmp(arg[iarg_+5],"angle"))
              error->fix_error(FLERR,this,"expecting keyword 'angle' after axis definition");
          rotateMesh(force->numeric(FLERR,arg[iarg_+2]),force->numeric(FLERR,arg[iarg_+3]),force->numeric(FLERR,arg[iarg_+4]),
                   force->numeric(FLERR,arg[iarg_+6]));
          manipulated_ = true;
          iarg_ += 7;
          hasargs = true;
      } else if(strcmp(arg[iarg_],"scale") == 0) {
          if (narg < iarg_+2) error->fix_error(FLERR,this,"not enough arguments");
          scaleMesh(force->numeric(FLERR,arg[iarg_+1]));
          manipulated_ = true;
          iarg_ += 2;
          hasargs = true;
      } else if (strcmp(arg[iarg_],"temperature") == 0) {
          iarg_++;
          double Temp_mesh = atof(arg[iarg_++]);
          mesh_->prop().addGlobalProperty< ScalarContainer<double> >("Temp","comm_none","frame_invariant","restart_yes");
          mesh_->prop().setGlobalProperty< ScalarContainer<double> >("Temp",Temp_mesh);
          mesh_->prop().addGlobalProperty< ScalarContainer<double> >("heatFlux","comm_none","frame_invariant","restart_no");
          mesh_->prop().setGlobalProperty< ScalarContainer<double> >("heatFlux",0.);
          mesh_->prop().addGlobalProperty< ScalarContainer<double> >("heatFluxTotal","comm_none","frame_invariant","restart_yes");
          mesh_->prop().setGlobalProperty< ScalarContainer<double> >("heatFluxTotal",0.);
          
          hasargs = true;
      }
    }
}

/* ---------------------------------------------------------------------- */

FixMesh::~FixMesh()
{
    delete mesh_;
}

/* ---------------------------------------------------------------------- */

void FixMesh::post_create()
{
    // check if all element property container have same length
    // could potentially be whacked by adding element properties
    // at the wrong place in code
    mesh_->check_element_property_consistency();
}

/* ---------------------------------------------------------------------- */

void FixMesh::create_mesh(char *mesh_fname)
{
    
    if(strncmp(style,"mesh/surface",12) == 0)
    {
        if(strcmp(style,"mesh/surface/stress/deform") == 0)
            mesh_ = new TriMeshDeformable(lmp);
        else if(strcmp(style,"mesh/surface/planar") == 0)
            mesh_ = new TriMeshPlanar(lmp);
        else
            mesh_ = new TriMesh(lmp);

        // set properties that are important for reading
        mesh_->setMeshID(id);
        if(verbose_) mesh_->setVerbose();
        if(autoRemoveDuplicates_) mesh_->autoRemoveDuplicates();
        if(precision_ > 0.) mesh_->setPrecision(precision_);

        // read file
        // can be from STL file or VTK file
        InputMeshTri *mesh_input = new InputMeshTri(lmp,0,NULL);
        
        mesh_input->meshtrifile(mesh_fname,static_cast<TriMesh*>(mesh_),verbose_);
        
        delete mesh_input;
    }
    else error->one(FLERR,"Illegal implementation of create_mesh();");
}

/* ---------------------------------------------------------------------- */

void FixMesh::create_mesh_restart()
{
    
    if(strcmp(style,"mesh/surface/stress/deform") == 0)
        mesh_ = new TriMeshDeformable(lmp);
    else if(strcmp(style,"mesh/surface/planar") == 0)
        mesh_ = new TriMeshPlanar(lmp);
    else if(strncmp(style,"mesh/surface",12) == 0)
        mesh_ = new TriMesh(lmp);
    else error->one(FLERR,"Illegal implementation of create_mesh();");

    // set properties that are important for reading
    mesh_->setMeshID(id);
    if(verbose_) mesh_->setVerbose();
    if(autoRemoveDuplicates_) mesh_->autoRemoveDuplicates();
    if(precision_ > 0.) mesh_->setPrecision(precision_);
}

/* ---------------------------------------------------------------------- */

void FixMesh::pre_delete(bool unfixflag)
{
    // error if moving mesh is operating on a mesh to be deleted
    
    if(unfixflag)
    {
        if(mesh_->isMoving() && modify->n_fixes_style("move/mesh") > 0)
            error->fix_error(FLERR,this,
                    "illegal unfix command, may not unfix a moving mesh while a fix move is applied."
                    "Unfix the fix move/mesh first");
    }
}

/* ---------------------------------------------------------------------- */

int FixMesh::setmask()
{
    int mask = 0;
    mask |= PRE_EXCHANGE;
    mask |= PRE_FORCE;
    mask |= FINAL_INTEGRATE;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixMesh::setup_pre_force(int vflag)
{
    // first-time set-up
    
    if(!setupFlag_)
    {
        initialSetup();
        setupFlag_ = true;
    }
    // if mesh already set-up and parallelized
    
    else
    {
        mesh_->pbcExchangeBorders(1);
    }

    // clear reverse comm properties
    mesh_->clearReverse();

    pOpFlag_ = false;

}

/* ---------------------------------------------------------------------- */

void FixMesh::initialSetup()
{
    
    mesh_->initalSetup();

    // warn if there are elements that extend outside box
    if(!mesh_->allNodesInsideSimulationBox() && 0 == comm->me)
       error->warning(FLERR,"Not all nodes of fix mesh inside simulation box, "
                            "elements will be deleted or wrapped around periodic boundary conditions");

    if(comm->me == 0)
       fprintf(screen,"Import and parallelization of mesh %s containing %d triangle(s) successful\n",
               id,mesh_->sizeGlobal());
}

/* ----------------------------------------------------------------------
   invoke parallelism
------------------------------------------------------------------------- */

void FixMesh::pre_exchange()
{
    
    pOpFlag_ = true;
}

/* ----------------------------------------------------------------------
   forward comm for mesh
------------------------------------------------------------------------- */

void FixMesh::pre_force(int vflag)
{
    // case re-neigh step
    if(pOpFlag_)
    {
        
        mesh_->pbcExchangeBorders(0);

        pOpFlag_ = false;
    }
    // case regular step
    else
    {
        mesh_->forwardComm();

        if(mesh_->decideRebuild())
        {
            
            next_reneighbor = update->ntimestep + 1;
        }
    }

    // clear reverse comm properties
    mesh_->clearReverse();
}

/* ----------------------------------------------------------------------
   reverse comm for mesh
------------------------------------------------------------------------- */

void FixMesh::final_integrate()
{
    
    mesh_->reverseComm();
}

/* ---------------------------------------------------------------------- */

int FixMesh::min_type()
{
    if(atom_type_mesh_ == -1) return 1;
    return atom_type_mesh_;
}

/* ---------------------------------------------------------------------- */

int FixMesh::max_type()
{
    if(atom_type_mesh_ == -1) return 1;
    return atom_type_mesh_;
}

/* ----------------------------------------------------------------------
   moves the mesh by a vector - (dx dy dz) is the displacement vector
------------------------------------------------------------------------- */

void FixMesh::moveMesh(double const dx, double const dy, double const dz)
{
    if (comm->me == 0)
    {
      //fprintf(screen,"moving mesh by ");
      //fprintf(screen,"%f %f %f\n", dx,dy,dz);
    }

    double arg[3] = {dx,dy,dz};
    mesh_->move(arg);
}

/* ----------------------------------------------------------------------
   rotates the mesh around an axis through the origin
   phi the rotation angle
   (axisX axisY axisZ) vector in direction of the axis
------------------------------------------------------------------------- */

void FixMesh::rotateMesh(double const axisX, double const axisY, double const axisZ, double const phi)
{
    double axis[3] = {axisX,axisY,axisZ}, p[3] = {0.,0.,0.};

    if(vectorMag3D(axis) < 1e-5)
        error->fix_error(FLERR,this,"illegal magnitude of rotation axis");

    if (comm->me == 0)
    {
      //fprintf(screen,"rotate ");
      //fprintf(screen,"%f %f %f %f\n", phi, axisX, axisY, axisZ);
    }

    mesh_->rotate(phi*3.14159265/180.0,axis,p);
}

/* ----------------------------------------------------------------------
   scales the mesh by a factor in x, y and z direction
   can also be used to distort a mesh
   (factorX factorY factorZ) vector contains the factors to scale the
   mesh in the xyz direction
------------------------------------------------------------------------- */

void FixMesh::scaleMesh(double const factor)
{
    if (comm->me == 0){
      //fprintf(screen,"scale ");
      //fprintf(screen,"%f \n", factor);
    }
    mesh_->scale(factor);
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixMesh::write_restart(FILE *fp)
{
    mesh_->writeRestart(fp);
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixMesh::restart(char *buf)
{
    double *list = (double *) buf;
    mesh_->restart(list);
}

/* ----------------------------------------------------------------------
   change box extent due to mesh node position
------------------------------------------------------------------------- */

void FixMesh::box_extent(double &xlo,double &xhi,double &ylo,double &yhi,double &zlo,double &zhi)
{
    double node[3];
    int size = mesh()->size();
    int numNodes = mesh()->numNodes();

    for(int i = 0; i < size; i++)
    {
        
        for(int j = 0; j < numNodes; j++)
        {
            mesh()->node_slow(i,j,node);
            xlo = MathExtraLiggghts::min(xlo,node[0]);
            xhi = MathExtraLiggghts::max(xhi,node[0]);
            ylo = MathExtraLiggghts::min(ylo,node[1]);
            yhi = MathExtraLiggghts::max(yhi,node[1]);
            zlo = MathExtraLiggghts::min(zlo,node[2]);
            zhi = MathExtraLiggghts::max(zhi,node[2]);
        }
    }
}
