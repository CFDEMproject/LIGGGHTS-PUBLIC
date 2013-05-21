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
   Christian Richter (OVGU Magdeburg, linear/variable and rotate/variable)
   Niels Dallinger (TU Chemnitz, viblin and vibrot)
------------------------------------------------------------------------- */

#include "mesh_mover.h"
#include "math.h"
#include "vector_liggghts.h"
#include "math_extra_liggghts.h"
#include "input.h"
#include "modify.h"
#include "variable.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   MeshMoverLinear
------------------------------------------------------------------------- */

MeshMoverLinear::MeshMoverLinear(LAMMPS *lmp,AbstractMesh *_mesh,
    FixMoveMesh *_fix_move_mesh,double vx, double vy, double vz)
  : MeshMover(lmp,_mesh,_fix_move_mesh)
{
    vel_[0] = vx;
    vel_[1] = vy;
    vel_[2] = vz;

    isFirst_ = mesh_->registerMove(false,true,false);
}

void MeshMoverLinear::pre_delete()
{
    mesh_->unregisterMove(false,true,false);
}

MeshMoverLinear::~MeshMoverLinear()
{

}

/* ---------------------------------------------------------------------- */

void MeshMoverLinear::initial_integrate(double dTAbs,double dTSetup,double dt)
{
    double dX[3],dx[3];

    int size = mesh_->size();
    int numNodes = mesh_->numNodes();
    double ***v_node = get_v();

    // calculate total and incremental displacement
    vectorScalarMult3D(vel_,dTSetup,dX);
    vectorScalarMult3D(vel_,dt,dx);

    // apply move
    mesh_->move(dX,dx);

    // set mesh velocity
    for (int i = 0; i < size; i++)
        for(int j = 0; j < numNodes; j++)
            vectorAdd3D(v_node[i][j],vel_,v_node[i][j]);
}

/* ----------------------------------------------------------------------
   MeshMoverLinearVariable
------------------------------------------------------------------------- */

MeshMoverLinearVariable::MeshMoverLinearVariable(LAMMPS *lmp,AbstractMesh *_mesh,
    FixMoveMesh *_fix_move_mesh, char* var1, char* var2, char* var3)
  : MeshMover(lmp,_mesh,_fix_move_mesh)
{
      int n;
      n = strlen(&var1[2]) + 1;
      var1str_ = new char[n];
      strcpy(var1str_,&var1[2]);
      myvar1_ = input->variable->find(var1str_);

      n = strlen(&var1[2]) + 1;
      var2str_ = new char[n];
      strcpy(var2str_,&var2[2]);
      myvar2_ = input->variable->find(var2str_);

      n = strlen(&var1[2]) + 1;
      var3str_ = new char[n];
      strcpy(var3str_,&var3[2]);
      myvar3_ = input->variable->find(var3str_);

      if (myvar1_ < 0) error->all(FLERR,"Variable name 1 for fix move/mesh linear/variable does not exist");
      if (myvar2_ < 0) error->all(FLERR,"Variable name 2 for fix move/mesh linear/variable does not exist");
      if (myvar3_ < 0) error->all(FLERR,"Variable name 3 for fix move/mesh linear/variable does not exist");

      vectorZeroize3D(dX_);

      vel_[0] = input->variable->compute_equal(myvar1_);
      vel_[1] = input->variable->compute_equal(myvar2_);
      vel_[2] = input->variable->compute_equal(myvar3_);

      isFirst_ = mesh_->registerMove(false,true,false);
}

void MeshMoverLinearVariable::pre_delete()
{
    mesh_->unregisterMove(false,true,false);
}

MeshMoverLinearVariable::~MeshMoverLinearVariable()
{
    delete []var1str_;
    delete []var2str_;
    delete []var3str_;
}

/* ---------------------------------------------------------------------- */

int MeshMoverLinearVariable::n_restart()
{
    return 3;
}

void MeshMoverLinearVariable::write_restart(double *buf)
{
    int n = 0;
    buf[n++] = dX_[0];
    buf[n++] = dX_[1];
    buf[n++] = dX_[2];
}

void MeshMoverLinearVariable::read_restart(double *buf)
{
    int n = 0;
    dX_[0] = buf[n++];
    dX_[1] = buf[n++];
    dX_[2] = buf[n++];
}

/* ---------------------------------------------------------------------- */

void MeshMoverLinearVariable::setup()
{
    
    vectorZeroize3D(dX_);
}

/* ---------------------------------------------------------------------- */

void MeshMoverLinearVariable::initial_integrate(double dTAbs,double dTSetup,double dt)
{
    double dx[3];

    int size = mesh_->size();
    int numNodes = mesh_->numNodes();
    double ***v_node = get_v();

    modify->clearstep_compute();

    // evaluate variable
    
    vel_[0] = input->variable->compute_equal(myvar1_);
    vel_[1] = input->variable->compute_equal(myvar2_);
    vel_[2] = input->variable->compute_equal(myvar3_);
    
    modify->addstep_compute(update->ntimestep + 1);

    // calculate total and incremental displacement
    vectorScalarMult3D(vel_,dt,dx);
    vectorAdd3D(dX_,dx,dX_);

    // apply move
    mesh_->move(dX_,dx);

    // set mesh velocity
    for (int i = 0; i < size; i++)
        for(int j = 0; j < numNodes; j++)
            vectorAdd3D(v_node[i][j],vel_,v_node[i][j]);
}

/* ----------------------------------------------------------------------
   MeshMoverWiggle
------------------------------------------------------------------------- */

MeshMoverWiggle::MeshMoverWiggle(LAMMPS *lmp,AbstractMesh *_mesh,
    FixMoveMesh *_fix_move_mesh, double ax, double ay, double az,double T)
  : MeshMover(lmp,_mesh,_fix_move_mesh),
    omega_(2.*M_PI/T)
{
    amplitude_[0] = ax;
    amplitude_[1] = ay;
    amplitude_[2] = az;

    isFirst_ = mesh_->registerMove(false,true,false);
}

void MeshMoverWiggle::pre_delete()
{
    mesh_->unregisterMove(false,true,false);
}

MeshMoverWiggle::~MeshMoverWiggle()
{

}

/* ---------------------------------------------------------------------- */

void MeshMoverWiggle::initial_integrate(double dTAbs,double dTSetup,double dt)
{

    double dX[3],dx[3],vNode[3];
    double sine =   sin(omega_ * dTAbs) - sin(omega_ * (dTAbs-dTSetup));
    double cosine = cos(omega_ * dTAbs) - cos(omega_ * (dTAbs-dTSetup));

    int size = mesh_->size();
    int numNodes = mesh_->numNodes();
    double ***v_node = get_v();

    // calculate velocity, same for all nodes
    vectorScalarMult3D(amplitude_,omega_*cos(omega_ * dTAbs),vNode);

    // calculate total and incremental displacement
    vectorScalarMult3D(amplitude_,sine,dX);
    vectorScalarMult3D(vNode,dt,dx);

    // apply move
    mesh_->move(dX,dx);

    // set mesh velocity
    for (int i = 0; i < size; i++)
        for(int j = 0; j < numNodes; j++)
            vectorAdd3D(v_node[i][j],vNode,v_node[i][j]);
}

/* ----------------------------------------------------------------------
   MeshMoverRotate
------------------------------------------------------------------------- */

MeshMoverRotate::MeshMoverRotate(LAMMPS *lmp,AbstractMesh *_mesh,
                                 FixMoveMesh *_fix_move_mesh,
                                 double px, double py, double pz,
                                 double axisX, double axisY, double axisZ,
                                 double T)
  : MeshMover(lmp,_mesh,_fix_move_mesh),
    omega_(2.*M_PI/T)
{
    axis_[0] = axisX;
    axis_[1] = axisY;
    axis_[2] = axisZ;

    vectorScalarDiv3D(axis_,vectorMag3D(axis_));

    point_[0] = px;
    point_[1] = py;
    point_[2] = pz;

    add_reference_point(point_);

    isFirst_ = mesh_->registerMove(false,true,true);
}

void MeshMoverRotate::pre_delete()
{
    mesh_->unregisterMove(false,true,true);
}

MeshMoverRotate::~MeshMoverRotate()
{

}

/* ---------------------------------------------------------------------- */

void MeshMoverRotate::initial_integrate(double dTAbs,double dTSetup,double dt)
{
    double xOld[3],node[3],vRot[3],omegaVec[3],rPA[3];
    double reference_point[3];
    double totalPhi = omega_*dTSetup;
    double incrementalPhi = omega_*dt;

    get_reference_point(reference_point);

    int size = mesh_->size();
    int numNodes = mesh_->numNodes();
    double ***v_node = get_v();
    double ***nodes = get_nodes();

    // rotate the mesh
    mesh_->rotate(totalPhi,incrementalPhi,axis_,reference_point);

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

MeshMoverRotateVariable::MeshMoverRotateVariable(LAMMPS *lmp,AbstractMesh *_mesh,
                                 FixMoveMesh *_fix_move_mesh,
                                 double px, double py, double pz,
                                 double axisX, double axisY, double axisZ,
                                 char* var1) // omega
  : MeshMover(lmp,_mesh,_fix_move_mesh)
{
    int n;
    n = strlen(&var1[2]) + 1;
    var1str_ = new char[n];
    strcpy(var1str_,&var1[2]);
    myvar1_ = input->variable->find(var1str_);
    if (myvar1_ < 0)
        error->all(FLERR,"Variable name 1 for fix move/mesh rotate dynamic does not exist");

    omega_ = input->variable->compute_equal(myvar1_);
    totalPhi_ = 0.;

    axis_[0] = axisX;
    axis_[1] = axisY;
    axis_[2] = axisZ;

    vectorScalarDiv3D(axis_,vectorMag3D(axis_));

    point_[0] = px;
    point_[1] = py;
    point_[2] = pz;

    add_reference_point(point_);

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

int MeshMoverRotateVariable::n_restart()
{
    return 1;
}

void MeshMoverRotateVariable::write_restart(double *buf)
{
    int n = 0;
    buf[n++] = totalPhi_;
}

void MeshMoverRotateVariable::read_restart(double *buf)
{
    int n = 0;
    totalPhi_ = buf[n++];
}

/* ---------------------------------------------------------------------- */

void MeshMoverRotateVariable::setup()
{
    
    totalPhi_ = 0.;
}

/* ---------------------------------------------------------------------- */

void MeshMoverRotateVariable::initial_integrate(double dTAbs,double dTSetup,double dt)
{
    double xOld[3],node[3],vRot[3],omegaVec[3],rPA[3];
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
    totalPhi_ += incrementalPhi;

    // rotate the mesh
    mesh_->rotate(totalPhi_,incrementalPhi,axis_,reference_point);

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

MeshMoverRiggle::MeshMoverRiggle(LAMMPS *lmp,AbstractMesh *_mesh,
                                 FixMoveMesh *_fix_move_mesh,
                                 double px, double py, double pz,
                                 double axisX, double axisY, double axisZ,
                                 double T, double ampl)
  : MeshMover(lmp,_mesh,_fix_move_mesh),
    omega_(2.*M_PI/T),
    amplitude_(ampl*M_PI/180.)
{
    axis_[0] = axisX;
    axis_[1] = axisY;
    axis_[2] = axisZ;

    vectorScalarDiv3D(axis_,vectorMag3D(axis_));

    point_[0] = px;
    point_[1] = py;
    point_[2] = pz;

    isFirst_ = mesh_->registerMove(false,true,true);
}

void MeshMoverRiggle::pre_delete()
{
    mesh_->unregisterMove(false,true,true);
}

MeshMoverRiggle::~MeshMoverRiggle()
{

}

/* ---------------------------------------------------------------------- */

void MeshMoverRiggle::initial_integrate(double dTAbs,double dTSetup,double dt)
{
    double xOld[3],node[3],vRot[3],omegaVec[3],rPA[3];

    double sine =   amplitude_*(sin(omega_ * dTAbs)-sin(omega_ * (dTAbs-dTSetup)));
    double vel_prefactor = omega_*amplitude_*cos(omega_ * dTAbs);

    int size = mesh_->size();
    int numNodes = mesh_->numNodes();
    double ***v_node = get_v();
    double ***nodes = get_nodes();

    // calculate total and incremental angle
    double totalPhi = sine;
    double incrementalPhi = vel_prefactor*dt;

    // rotate the mesh
    mesh_->rotate(totalPhi,incrementalPhi,axis_,point_);

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
   MeshMoverVibLin
------------------------------------------------------------------------- */

MeshMoverVibLin::MeshMoverVibLin(LAMMPS *lmp,AbstractMesh *_mesh,
                                 FixMoveMesh *_fix_move_mesh,
                                 double axisX, double axisY, double axisZ,
                                 int order, double amplitude[10], double phase[10],
                                 double T)
  : MeshMover(lmp,_mesh,_fix_move_mesh), omega_(2.*M_PI/T)
{
    axis_[0] = axisX;
    axis_[1] = axisY;
    axis_[2] = axisZ;
    ord = order;

    vectorScalarDiv3D(axis_,vectorMag3D(axis_));
    //array transfer
    for (int j=0;j<order; j++) {
       phi[j] = phase[j];
       ampl[j] = amplitude[j];
     }

    isFirst_ = mesh_->registerMove(false,true,false);
}

void MeshMoverVibLin::pre_delete()
{
    mesh_->unregisterMove(false,true,false);
}

MeshMoverVibLin::~MeshMoverVibLin()
{

}

/* ---------------------------------------------------------------------- */
void MeshMoverVibLin::initial_integrate(double dTAbs,double dTSetup,double dt)
{
    double dX[3],dx[3],vNode[3];
    int size = mesh_->size();
    int numNodes = mesh_->numNodes();
    double ***v_node = get_v();

    double arg = 0;
    double vA = 0;

    for (int j=0;j<ord; j++)
    {
        arg = arg+ampl[j]*(cos(omega_*(j+1)*dTAbs+phi[j]) - cos(omega_*(j+1)*(dTAbs-dTSetup)+phi[j]));
        vA= vA-ampl[j]*(j+1)*omega_*sin(omega_*(j+1)*dTAbs+phi[j]);
    }

    // calculate velocity, same for all nodes
    vectorScalarMult3D(axis_,vA,vNode);
    // calculate total and incremental displacement
    vectorScalarMult3D(axis_,arg,dX); //total
    vectorScalarMult3D(vNode,dt,dx); //incremental

    // apply linear move
    mesh_->move(dX,dx);

    // set mesh velocity
    for (int i = 0; i < size; i++)
        for(int j = 0; j < numNodes; j++)
            vectorAdd3D(v_node[i][j],vNode,v_node[i][j]);

}

/* ----------------------------------------------------------------------
   MeshMoverVibRot
------------------------------------------------------------------------- */

MeshMoverVibRot::MeshMoverVibRot(LAMMPS *lmp,AbstractMesh *_mesh,
                                 FixMoveMesh *_fix_move_mesh,
                                 double px, double py, double pz,
                                 double axisX, double axisY, double axisZ,
                                 int order, double amplitude[10], double phase[10],
                                 double T)
  : MeshMover(lmp,_mesh,_fix_move_mesh), omega_(2.*M_PI/T)
{
    axis_[0] = axisX;
    axis_[1] = axisY;
    axis_[2] = axisZ;

    vectorScalarDiv3D(axis_,vectorMag3D(axis_));

    p_[0] = px;
    p_[1] = py;
    p_[2] = pz;
    ord = order;

    //array transfer
    for (int j=0;j<order; j++) {
       phi[j] = phase[j];
       ampl[j] = amplitude[j];
     }

    isFirst_ = mesh_->registerMove(false,true,true);
}
void MeshMoverVibRot::pre_delete()
{
    mesh_->unregisterMove(false,true,true);
}

MeshMoverVibRot::~MeshMoverVibRot()
{

}

/* ---------------------------------------------------------------------- */

void MeshMoverVibRot::initial_integrate(double dTAbs,double dTSetup,double dt)
{
    double xOld[3],node[3],omegaVec[3],rPA[3],vRot[3];

    double arg = 0;
    double vR = 0;

    for (int j=0;j<ord; j++)
    {
        arg = arg+ampl[j]* ( cos(omega_*(j+1)*dTAbs+phi[j]) - cos(omega_*(j+1)*(dTAbs-dTSetup)+phi[j]) );
        vR = vR-ampl[j]*(j+1)*omega_*sin(omega_*(j+1)*dTAbs+phi[j]);
    }

    int size = mesh_->size();
    int numNodes = mesh_->numNodes();
    double ***v_node = get_v();
    double ***nodes = get_nodes();

    double totalPhi = arg;
    double incrementalPhi = vR*dt;

    // rotate the mesh
    mesh_->rotate(totalPhi,incrementalPhi,axis_,p_);

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

