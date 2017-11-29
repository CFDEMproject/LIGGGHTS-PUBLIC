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
    Niels Dallinger (TU Chemnitz, viblin and vibrot)
    Christian Richter (OVGU Magdeburg, linear variable and rotate/variable)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
    Copyright 2013-2015 TU Chemnitz
    Copyright 2013      OVGU Magdeburg
------------------------------------------------------------------------- */

#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include "mesh_mover_rotation.h"
#include "vector_liggghts.h"
#include "math_extra_liggghts.h"
#include "input.h"
#include "modify.h"
#include "variable.h"
#include "fix_mesh.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   MeshMoverRotate
------------------------------------------------------------------------- */

MeshMoverRotate::MeshMoverRotate(LAMMPS *lmp,AbstractMesh *_mesh, FixMoveMesh *_fix_move_mesh,
                                 const char * const * const arg, const int narg) :
    MeshMover(lmp,_mesh,_fix_move_mesh)
{
    if (narg < 11)
        error->all(FLERR, "Not enough arguments for movement type rotate");
    if (narg > 11)
        error->warning(FLERR, "Movement type rotate requires only 11 arguments, excess arguments will be ignored");

    if (strcmp(arg[1], "origin"))
        error->all(FLERR, "Expected keyword 'origin'");
    point_[0] = force->numeric(FLERR, arg[2]);
    point_[1] = force->numeric(FLERR, arg[3]);
    point_[2] = force->numeric(FLERR, arg[4]);

    if (strcmp(arg[5], "axis"))
        error->all(FLERR, "Expected keyword 'axis'");
    axis_[0] = force->numeric(FLERR, arg[6]);
    axis_[1] = force->numeric(FLERR, arg[7]);
    axis_[2] = force->numeric(FLERR, arg[8]);
    vectorNormalize3D(axis_);

    if (strcmp(arg[9], "period"))
        error->all(FLERR, "Expected keyword 'period'");
    omega_ = 2.*M_PI/force->numeric(FLERR, arg[10]);

    add_reference_point(point_);
}

void MeshMoverRotate::post_create()
{
    isFirst_ = mesh_->registerMove(false,true,true);
}

void MeshMoverRotate::pre_delete()
{
    mesh_->unregisterMove(false,true,true);
}

MeshMoverRotate::~MeshMoverRotate()
{ }

/* ---------------------------------------------------------------------- */

void MeshMoverRotate::initial_integrate(double dTAbs,double dTSetup,double dt)
{
    double node[3],vRot[3],omegaVec[3],rPA[3];
    double reference_point[3];
    double incrementalPhi = omega_*dt;

    get_reference_point(reference_point);

    int size = mesh_->size();
    int numNodes = mesh_->numNodes();
    double ***v_node = get_v();
    double ***nodes = get_nodes();

    // rotate the mesh
    fix_move_mesh_->fixMesh()->rotate(incrementalPhi,axis_,reference_point, fix_move_mesh_);

    // set mesh velocity, w x rPA
    vectorScalarMult3D(axis_,omega_,omegaVec);
    for(int i = 0; i < size; i++)
    {
        for(int iNode = 0; iNode < numNodes; iNode++)
        {
            vectorCopy3D(nodes[i][iNode],node);
            vectorSubtract3D(node,reference_point,rPA);
            vectorCross3D(omegaVec,rPA,vRot);
            vectorAdd3D(v_node[i][iNode],vRot,v_node[i][iNode]);
        }
    }
}

/* ----------------------------------------------------------------------
   MeshMoverRotateVariable
------------------------------------------------------------------------- */

MeshMoverRotateVariable::MeshMoverRotateVariable(LAMMPS *lmp,AbstractMesh *_mesh, FixMoveMesh *_fix_move_mesh,
                                 const char * const * const arg, const int narg) :
    MeshMover(lmp,_mesh,_fix_move_mesh)
{
    if (narg < 11)
        error->all(FLERR, "Not enough arguments for movement type rotate/variable");
    if (narg > 11)
        error->warning(FLERR, "Movement type rotate/variable requires only 11 arguments, excess arguments will be ignored");

    if (strcmp(arg[1], "origin"))
        error->all(FLERR, "Expected keyword 'origin'");
    point_[0] = force->numeric(FLERR, arg[2]);
    point_[1] = force->numeric(FLERR, arg[3]);
    point_[2] = force->numeric(FLERR, arg[4]);

    if (strcmp(arg[5], "axis"))
        error->all(FLERR, "Expected keyword 'axis'");
    axis_[0] = force->numeric(FLERR, arg[6]);
    axis_[1] = force->numeric(FLERR, arg[7]);
    axis_[2] = force->numeric(FLERR, arg[8]);
    vectorNormalize3D(axis_);

    if (strcmp(arg[9], "omega"))
        error->all(FLERR, "Expected keyword 'omega'");
    int n;
    n = strlen(&arg[10][2]) + 1;
    var1str_ = new char[n];
    strcpy(var1str_,&arg[10][2]);
    myvar1_ = input->variable->find(var1str_);
    if (myvar1_ < 0)
        error->all(FLERR,"Variable name 1 for fix move/mesh rotate/variable does not exist");

    omega_ = input->variable->compute_equal(myvar1_);
    totalPhi_ = 0.;

    add_reference_point(point_);
}

void MeshMoverRotateVariable::post_create()
{
    isFirst_ = mesh_->registerMove(false,true,true);
}

void MeshMoverRotateVariable::pre_delete()
{
    mesh_->unregisterMove(false,true,true);
}

MeshMoverRotateVariable::~MeshMoverRotateVariable()
{
    delete []var1str_;
}

/* ---------------------------------------------------------------------- */

void MeshMoverRotateVariable::setup()
{
    // check if variable still exists
    myvar1_ = input->variable->find(var1str_);
    if (myvar1_ < 0)
        error->all(FLERR,"Variable name 1 for fix move/mesh rotate dynamic does not exist");
}

/* ---------------------------------------------------------------------- */

void MeshMoverRotateVariable::initial_integrate(double,double,double dt)
{
    double node[3],vRot[3],omegaVec[3],rPA[3];
    double reference_point[3];
    double incrementalPhi;

    int size = mesh_->size();
    int numNodes = mesh_->numNodes();
    double ***v_node = get_v();
    double ***nodes = get_nodes();

    modify->clearstep_compute();

    // re-evaluation of omega (global,private)
    omega_ = input->variable->compute_equal(myvar1_);

    modify->addstep_compute(update->ntimestep + 1);

    get_reference_point(reference_point);

    incrementalPhi = omega_*dt;

    // rotate the mesh
    fix_move_mesh_->fixMesh()->rotate(incrementalPhi,axis_,reference_point, fix_move_mesh_);

    // set mesh velocity, w x rPA
    vectorScalarMult3D(axis_,omega_,omegaVec);
    for(int i = 0; i < size; i++)
    {
        for(int iNode = 0; iNode < numNodes; iNode++)
        {
            vectorCopy3D(nodes[i][iNode],node);
            vectorSubtract3D(node,reference_point,rPA);
            vectorCross3D(omegaVec,rPA,vRot);
            vectorAdd3D(v_node[i][iNode],vRot,v_node[i][iNode]);
        }
    }
}

/* ----------------------------------------------------------------------
   MeshMoverRiggle
------------------------------------------------------------------------- */

MeshMoverRiggle::MeshMoverRiggle(LAMMPS *lmp,AbstractMesh *_mesh, FixMoveMesh *_fix_move_mesh,
                                 const char * const * const arg, const int narg) :
    MeshMover(lmp,_mesh,_fix_move_mesh)
{
    if (narg < 13)
        error->all(FLERR, "Not enough arguments for movement type riggle");
    if (narg > 13)
        error->warning(FLERR, "Movement type riggle requires only 13 arguments, excess arguments will be ignored");

    if (strcmp(arg[1], "origin"))
        error->all(FLERR, "Expected keyword 'origin'");
    point_[0] = force->numeric(FLERR, arg[2]);
    point_[1] = force->numeric(FLERR, arg[3]);
    point_[2] = force->numeric(FLERR, arg[4]);

    if (strcmp(arg[5], "axis"))
        error->all(FLERR, "Expected keyword 'axis'");
    axis_[0] = force->numeric(FLERR, arg[6]);
    axis_[1] = force->numeric(FLERR, arg[7]);
    axis_[2] = force->numeric(FLERR, arg[8]);
    vectorNormalize3D(axis_);

    if (strcmp(arg[9], "period"))
        error->all(FLERR, "Expected keyword 'period'");
    omega_ = 2.*M_PI/force->numeric(FLERR, arg[10]);

    if (strcmp(arg[11], "amplitude"))
        error->all(FLERR, "Expected keyword 'amplitude'");
    amplitude_ = M_PI*force->numeric(FLERR, arg[12])/180.;
}

void MeshMoverRiggle::post_create()
{
    isFirst_ = mesh_->registerMove(false,true,true);
}

void MeshMoverRiggle::pre_delete()
{
    mesh_->unregisterMove(false,true,true);
}

MeshMoverRiggle::~MeshMoverRiggle()
{ }

/* ---------------------------------------------------------------------- */

void MeshMoverRiggle::initial_integrate(double dTAbs,double dTSetup,double dt)
{
    double node[3],vRot[3],omegaVec[3],rPA[3];

    double vel_prefactor = omega_*amplitude_*cos(omega_ * dTAbs);

    int size = mesh_->size();
    int numNodes = mesh_->numNodes();
    double ***v_node = get_v();
    double ***nodes = get_nodes();

    // calculate total and incremental angle
    double incrementalPhi = vel_prefactor*dt;

    // rotate the mesh
    fix_move_mesh_->fixMesh()->rotate(incrementalPhi,axis_,point_, fix_move_mesh_);

    // set mesh velocity, vel_prefactor * w/|w| x rPA
    vectorScalarMult3D(axis_,vel_prefactor,omegaVec);
    for(int i = 0; i < size; i++)
    {
        for(int iNode = 0; iNode < numNodes; iNode++)
        {
            vectorCopy3D(nodes[i][iNode],node);
            vectorSubtract3D(node,point_,rPA);
            vectorCross3D(omegaVec,rPA,vRot);
            vectorAdd3D(v_node[i][iNode],vRot,v_node[i][iNode]);
        }
    }
}

/* ----------------------------------------------------------------------
   MeshMoverVibRot
------------------------------------------------------------------------- */

MeshMoverVibRot::MeshMoverVibRot(LAMMPS *lmp,AbstractMesh *_mesh, FixMoveMesh *_fix_move_mesh,
                                 const char * const * const arg, const int narg) :
    MeshMover(lmp,_mesh,_fix_move_mesh)
{
    if (narg < 11)
        error->all(FLERR, "Not enough arguments for movement type vibrot");
    if (strcmp(arg[9], "order"))
        error->all(FLERR, "Expected keyword 'order'");
    ord = force->inumeric(FLERR, arg[10]);
    if (ord > 30 || ord < 1)
        error->all(FLERR, "order can be at most 30 and must be greater 0");
    if (narg < 14+2*ord)
        error->all(FLERR, "Not enough arguments for movement type vibrot");
    if (narg > 14+2*ord)
        error->warning(FLERR, "Movement type vibrot requires only (14 + 2*$order) arguments, excess arguments will be ignored");

    if (strcmp(arg[1], "origin"))
        error->all(FLERR, "Expected keyword 'origin'");
    p_[0] = force->numeric(FLERR, arg[2]);
    p_[1] = force->numeric(FLERR, arg[3]);
    p_[2] = force->numeric(FLERR, arg[4]);

    if (strcmp(arg[5], "axis"))
        error->all(FLERR, "Expected keyword 'axis'");
    axis_[0] = force->numeric(FLERR, arg[6]);
    axis_[1] = force->numeric(FLERR, arg[7]);
    axis_[2] = force->numeric(FLERR, arg[8]);
    vectorNormalize3D(axis_);

    if (strcmp(arg[11], "amplitude"))
        error->all(FLERR, "Expected keyword 'amplitude'");
    if (strcmp(arg[12+ord], "phase"))
        error->all(FLERR, "Expected keyword 'phase'");
    if (strcmp(arg[13+2*ord], "period"))
        error->all(FLERR, "Expected keyword 'period'");
    //array transfer
    for (int j=0;j<ord; j++)
    {
       ampl[j] = force->numeric(FLERR, arg[12+j]);
       phi[j] = force->numeric(FLERR, arg[13+ord+j]);
       omega[j] = 2.*M_PI/force->numeric(FLERR, arg[14+2*ord+j]);
    }
}

void MeshMoverVibRot::post_create()
{
    isFirst_ = mesh_->registerMove(false,true,true);
}

void MeshMoverVibRot::pre_delete()
{
    mesh_->unregisterMove(false,true,true);
}

MeshMoverVibRot::~MeshMoverVibRot()
{ }

/* ---------------------------------------------------------------------- */

void MeshMoverVibRot::initial_integrate(double dTAbs,double dTSetup,double dt)
{
    double node[3],omegaVec[3],rPA[3],vRot[3];

    double arg = 0;
    double vR = 0;

    for (int j=0;j<ord; j++)
    {
        arg = arg+ampl[j]* ( cos(omega[j]*dTAbs+phi[j]) - cos(omega[j]*(dTAbs-dTSetup)+phi[j]) );
        vR = vR-ampl[j]*omega[j]*sin(omega[j]*dTAbs+phi[j]);
    }

    int size = mesh_->size();
    int numNodes = mesh_->numNodes();
    double ***v_node = get_v();
    double ***nodes = get_nodes();

    double incrementalPhi = vR*dt;

    // rotate the mesh
    fix_move_mesh_->fixMesh()->rotate(incrementalPhi,axis_,p_, fix_move_mesh_);

    // set mesh velocity, vel_prefactor * w/|w| x rPA
    vectorScalarMult3D(axis_,vR,omegaVec);

    for(int i = 0; i < size; i++)
    {
        for(int iNode = 0; iNode < numNodes; iNode++)
        {
            vectorCopy3D(nodes[i][iNode],node);
            vectorSubtract3D(node,p_,rPA);
            vectorCross3D(omegaVec,rPA,vRot);
            vectorAdd3D(v_node[i][iNode],vRot,v_node[i][iNode]);
        }
    }
}
