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
#include "multi_node_mesh.h"
#include "modify.h"
#include "comm.h"
#include "math_extra.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define EPSILON_V 0.00001

FixMesh::FixMesh(LAMMPS *lmp, int narg, char **arg)
: Fix(lmp, narg, arg),
  fix_contact_history_(0),
  fix_mesh_neighlist_(0),
  surfaceVel_(false),
  setupFlag_(false),
  pOpFlag_(false)
{
    if(narg < 7)
      error->fix_error(FLERR,this,"not enough arguments - filename and walltype are required.");

    iarg_ = 3;

    char mesh_fname[256];
    if(strcmp(arg[iarg_++],"file"))
        error->fix_error(FLERR,this,"expecting keyword 'file'");

    strcpy(mesh_fname,arg[iarg_++]);
    if(strcmp(arg[iarg_++],"type"))
        error->fix_error(FLERR,this,"expecting keyword 'type'");

    atom_type_wall_ = force->inumeric(arg[iarg_++]);

    mesh_ = new TriMesh(lmp);
    mesh_->setMeshID(id);

    // read file
    // can be from STL file or VTK file
    InputMeshTri *mesh_input = new InputMeshTri(lmp,0,NULL);
    mesh_input->meshtrifile(mesh_fname,mesh_);
    delete mesh_input;

    // parse further args

    bool hasargs = true;
    while(iarg_ < narg && hasargs)
    {
      hasargs = false;
      if(strcmp(arg[iarg_],"move") == 0){
          if (narg < iarg_+4) error->fix_error(FLERR,this,"not enough arguments");
          moveMesh(force->numeric(arg[iarg_+1]),force->numeric(arg[iarg_+2]),force->numeric(arg[iarg_+3]));
          iarg_ += 4;
          hasargs = true;
      } else if(strcmp(arg[iarg_],"rotate") == 0){
          if (narg < iarg_+7)
              error->fix_error(FLERR,this,"not enough arguments");
          if(strcmp(arg[iarg_+1],"axis"))
              error->fix_error(FLERR,this,"expecting keyword 'axis' after keyword 'rotate'");
          if(strcmp(arg[iarg_+5],"angle"))
              error->fix_error(FLERR,this,"expecting keyword 'angle' after axis definition");
          rotateMesh(force->numeric(arg[iarg_+2]),force->numeric(arg[iarg_+3]),force->numeric(arg[iarg_+4]),
                   force->numeric(arg[iarg_+6]));
          iarg_ += 7;
          hasargs = true;
      } else if(strcmp(arg[iarg_],"scale") == 0){
          if (narg < iarg_+2) error->fix_error(FLERR,this,"not enough arguments");
          scaleMesh(force->numeric(arg[iarg_+1]));
          iarg_ += 2;
          hasargs = true;
      } else if (strcmp(arg[iarg_],"temperature") == 0) {
          iarg_++;
          double Temp_wall = atof(arg[iarg_++]);
          mesh_->prop().addMeshProperty< ScalarContainer<double> >("Temp","comm_none","frame_invariant");
          mesh_->prop().setMeshProperty< ScalarContainer<double> >("Temp",Temp_wall);
          mesh_->prop().addMeshProperty< ScalarContainer<double> >("heatFlux","comm_none","frame_invariant");
          mesh_->prop().setMeshProperty< ScalarContainer<double> >("heatFlux",0);
          mesh_->prop().addMeshProperty< ScalarContainer<double> >("heatFluxTotal","comm_none","frame_invariant");
          mesh_->prop().setMeshProperty< ScalarContainer<double> >("heatFluxTotal",0);
          
          hasargs = true;
      } else if (strcmp(arg[iarg_],"surface_vel") == 0) {
          if (narg < iarg_+4) error->fix_error(FLERR,this,"not enough arguments");
          iarg_++;
          double vel[3];
          vel[0] = force->numeric(arg[iarg_++]);
          vel[1] = force->numeric(arg[iarg_++]);
          vel[2] = force->numeric(arg[iarg_++]);
          initVel(vel);
          hasargs = true;
      } else if (strcmp(arg[iarg_],"surface_ang_vel") == 0) {
          if (narg < iarg_+11) error->fix_error(FLERR,this,"not enough arguments");
          iarg_++;
          if(strcmp(arg[iarg_++],"origin"))
              error->fix_error(FLERR,this,"expecting keyword 'origin' after 'rotation'");
          double origin[3];
          origin[0] = force->numeric(arg[iarg_++]);
          origin[1] = force->numeric(arg[iarg_++]);
          origin[2] = force->numeric(arg[iarg_++]);
          if(strcmp(arg[iarg_++],"axis"))
              error->fix_error(FLERR,this,"expecting keyword 'axis' after definition of 'origin'");
          double axis[3];
          axis[0] = force->numeric(arg[iarg_++]);
          axis[1] = force->numeric(arg[iarg_++]);
          axis[2] = force->numeric(arg[iarg_++]);
          if(strcmp(arg[iarg_++],"omega"))
              error->fix_error(FLERR,this,"expecting keyword 'omega' after definition of 'axis'");
          // positive omega give anti-clockwise (CCW) rotation
          double omega = force->numeric(arg[iarg_++]);
          initAngVel(origin,axis,omega);
          hasargs = true;
      } else if (strcmp(arg[iarg_],"curvature") == 0) {
          if (narg < iarg_+2) error->fix_error(FLERR,this,"not enough arguments");
          iarg_++;
          double curvature = force->numeric(arg[iarg_++]);
          if(curvature < 0. || curvature > 60)
            error->fix_error(FLERR,this,"0° < curvature < 60° required");
          curvature = cos(curvature*M_PI/180.);
          mesh_->setCurvature(curvature);
      } else if(strcmp(style,"mesh") == 0 || strcmp(style,"mesh/gran") == 0) {
          char *errmsg = new char[strlen(arg[iarg_])+20];
          sprintf(errmsg,"unknown keyword: %s", arg[iarg_]);
          error->fix_error(FLERR,this,errmsg);
          delete []errmsg;
      }
    }
}

/* ---------------------------------------------------------------------- */

FixMesh::~FixMesh()
{
    if(mesh_->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v"))
        mesh_->prop().removeElementProperty("v");

    delete mesh_;
}

/* ---------------------------------------------------------------------- */

void FixMesh::post_create()
{

}

/* ---------------------------------------------------------------------- */

void FixMesh::pre_delete(bool unfixflag)
{
    // error if moving mesh is operating on a mesh to be deleted
    
    if(unfixflag)
    {
        if(mesh()->isMoving())
            error->fix_error(FLERR,this,
                    "illegal unfix command, may not unfix a mesh while a fix move is applied."
                    "Unfix the fix move/mesh first");
    }

    // contact tracker and neighlist are created via fix wall/gran
    if(fix_contact_history_)
      modify->delete_fix(fix_contact_history_->id);
    if(fix_mesh_neighlist_)
      modify->delete_fix(fix_mesh_neighlist_->id);
}

/* ---------------------------------------------------------------------- */

int FixMesh::setmask()
{
    int mask = 0;
    mask |= PRE_EXCHANGE;
    mask |= PRE_FORCE;
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
        mesh_->pbcExchangeBorders();
    }

    pOpFlag_ = false;

    // create neigh list for owned and local elements
    
    if(meshNeighlist())
        meshNeighlist()->initializeNeighlist();

}

/* ---------------------------------------------------------------------- */
void FixMesh::initialSetup()
{
    
    mesh_->initalSetup();

    // warn if there are elements that extend outside box
    if(!mesh_->allNodesInsideSimulationBox())
       error->warning(FLERR,"Not all nodes of fix mesh inside simulation box, "
                            "elements will be deleted or wrapped around periodic boundary conditions");

    if(comm->me == 0)
       fprintf(screen,"Import and parallelization of mesh %s containing %d triangle(s) successful\n",
               id,mesh_->sizeGlobal());
}

/* ---------------------------------------------------------------------- */

int FixMesh::min_type()
{
    return atom_type_wall_;
}

/* ---------------------------------------------------------------------- */

int FixMesh::max_type()
{
    return atom_type_wall_;
}

/* ----------------------------------------------------------------------
   called from fix wall/gran out of post_create()
------------------------------------------------------------------------- */

void FixMesh::createNeighList()
{
    if(fix_mesh_neighlist_) return;
    char *neighlist_name = new char[strlen(id)+1+10];
    sprintf(neighlist_name,"neighlist_%s",id);

    char **fixarg = new char*[4];
    fixarg[0]= neighlist_name;
    fixarg[1]= (char *) "all";
    fixarg[2]= (char *) "neighlist/mesh";
    fixarg[3]= id;
    modify->add_fix(4,fixarg);

    fix_mesh_neighlist_ =
        static_cast<FixNeighlistMesh*>(modify->find_fix_id(neighlist_name));

    delete []fixarg;
    delete []neighlist_name;

    /*
    
    fix_mesh_neighlist_->pre_neighbor();
    fix_mesh_neighlist_->pre_force(0);
    */
}

/* ----------------------------------------------------------------------
   called from fix wall/gran out of post_create()
------------------------------------------------------------------------- */

void FixMesh::createContactHistory(int dnum)
{
    if(fix_contact_history_) return;

    // create a contact tracker for the mesh
    char *contacthist_name = new char[strlen(id)+1+8];
    char *my_id = new char[strlen(id)+1];
    sprintf(contacthist_name,"tracker_%s",id);
    sprintf(my_id,"%s",id);

    char dnum_char[10];
    sprintf(dnum_char,"%d",dnum);

    char **fixarg = new char*[6];
    fixarg[0] = contacthist_name;
    fixarg[1] = (char *) "all";
    fixarg[2] = (char *) "contacthistory/mesh";
    fixarg[3] = (char *) "mesh";
    fixarg[4] = my_id;
    fixarg[5]= dnum_char;

    modify->add_fix(6,fixarg);

    fix_contact_history_ = static_cast<FixContactHistory*>(modify->find_fix_id(contacthist_name));

    delete []fixarg;
    delete []contacthist_name;
    delete []my_id;
}

/* ----------------------------------------------------------------------
   invoke parallelism
------------------------------------------------------------------------- */

void FixMesh::pre_exchange()
{
    // flag parallel operations on this step
    
    pOpFlag_ = true;
}

/* ----------------------------------------------------------------------
   forward comm for mesh, currently no reverse comm invoked
------------------------------------------------------------------------- */

void FixMesh::pre_force(int vflag)
{
    
    // case re-neigh step
    if(pOpFlag_)
    {
        mesh_->pbcExchangeBorders();
        pOpFlag_ = false;
    }
    // case regular step
    else
        mesh_->forwardComm();
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

    if (comm->me == 0)
    {
      //fprintf(screen,"rotate ");
      //fprintf(screen,"%f %f %f %f\n", phi, axisX, axisY, axisZ);
    }

    static_cast<MultiNodeMesh<3>*>(mesh_)->rotate(phi*3.14159265/180.0,axis,p);
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
   sets mesh velocity for conveyor model
------------------------------------------------------------------------- */

void FixMesh::initVel(double *conv_vel)
{
    int size, nVec;
    double scp, tmp[3], facenormal[3], ***v_node;
    double conv_vel_mag = vectorMag3D(conv_vel);

    surfaceVel_ = true;

    // register new mesh property v
    mesh_->prop().addMeshProperty< VectorContainer<double,3> >("v","comm_none","frame_invariant");
    mesh_->prop().setMeshProperty< VectorContainer<double,3> >("v",conv_vel);

    // register new element property v - error if exists already via moving mesh
    mesh_->prop().addElementProperty<MultiVectorContainer<double,3,3> >("v","comm_none","frame_invariant");
    size = mesh_->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v")->size();
    nVec = mesh_->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v")->nVec();

    v_node = mesh_->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v")->begin();

    // set mesh velocity
    for (int i = 0; i < size; i++)
    {
        mesh_->surfaceNorm(i,facenormal);
        scp = vectorDot3D(conv_vel,facenormal);
        vectorScalarMult3D(facenormal,scp,tmp);
        for(int j = 0; j < nVec; j++)
        {
            vectorSubtract3D(conv_vel,tmp,v_node[i][j]);
            if(vectorMag3D(v_node[i][j]) > 0.)
            {
                vectorScalarDiv3D(v_node[i][j],vectorMag3D(v_node[i][j]));
                vectorScalarMult3D(v_node[i][j],conv_vel_mag);
            }
        }
    }
}

/* ----------------------------------------------------------------------
   sets mesh velocity for rotation model
------------------------------------------------------------------------- */

void FixMesh::initAngVel(double *rot_origin,double *rot_axis,double rot_omega)
{
    int size, nVec;
    double tmp[3], scp, unitAxis[3], tangComp[3], Utang[3], surfaceV[3];
    double node[3],facenormal[3], magAxis, magUtang, ***v_node;

    surfaceVel_ = true;

    // register new element property v - error if exists already via moving mesh
    mesh_->prop().addElementProperty<MultiVectorContainer<double,3,3> >("v","comm_none","frame_invariant");
    size = mesh_->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v")->size();
    nVec = mesh_->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v")->nVec();

    size = mesh_->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v")->size();
    nVec = mesh_->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v")->nVec();
    v_node = mesh_->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v")->begin();

    // calculate unit vector of rotation axis
    magAxis = vectorMag3D(rot_axis);
    vectorScalarDiv3D(rot_axis,magAxis,unitAxis);

    for (int i = 0; i < size; i++)
    {
        mesh_->surfaceNorm(i,facenormal);
        // number of nodes per face (3)
        for(int j = 0; j < nVec; j++)
        {
            // position of node - origin of rotation (to get lever arm)
            mesh_->node(i,j,node);
            vectorSubtract3D(node,rot_origin,surfaceV);

    	    // lever arm X rotational axis = tangential component
    	    vectorCross3D(surfaceV,unitAxis,tangComp);

    	    // multiplying by omega scales the tangential component to give tangential velocity
    	    vectorScalarMult3D(tangComp,-rot_omega,Utang);

            // EPSILON is wall velocity, not rotational omega
            if(vectorMag3D(Utang) < EPSILON_V)
                error->fix_error(FLERR,this,"Rotation velocity too low");

            scp = vectorDot3D(Utang, facenormal);
            vectorScalarMult3D(facenormal,scp,tmp);

            // removes components normal to wall
            vectorSubtract3D(Utang,tmp,v_node[i][j]);
    	    magUtang = vectorMag3D(Utang);

            if(vectorMag3D(v_node[i][j]) > 0.)
            {
               vectorScalarDiv3D(v_node[i][j],vectorMag3D(v_node[i][j]));
               vectorScalarMult3D(v_node[i][j],magUtang);
            }
        }
    }
}
