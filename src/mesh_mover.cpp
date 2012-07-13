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

#include "mesh_mover.h"
#include "math.h"
#include "vector_liggghts.h"
#include "math_extra_liggghts.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   MeshMoverLinear
------------------------------------------------------------------------- */

MeshMoverLinear::MeshMoverLinear(LAMMPS *lmp,AbstractMesh *_mesh, double vx, double vy, double vz)
  : MeshMover(lmp,_mesh)
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

void MeshMoverLinear::initial_integrate(double dT,double dt)
{
    double dX[3],dx[3];

    int size = mesh_->size();
    int numNodes = mesh_->numNodes();
    double ***v_node = get_v();

    // calculate total and incremental displacement
    vectorScalarMult3D(vel_,dT,dX);
    vectorScalarMult3D(vel_,dt,dx);

    // apply move
    mesh_->move(dX,dx);

    // set mesh velocity
    for (int i = 0; i < size; i++)
        for(int j = 0; j < numNodes; j++)
            vectorAdd3D(v_node[i][j],vel_,v_node[i][j]);
}

/* ----------------------------------------------------------------------
   MeshMoverWiggle
------------------------------------------------------------------------- */

MeshMoverWiggle::MeshMoverWiggle(LAMMPS *lmp,AbstractMesh *_mesh,
    double ax, double ay, double az,double T)
  : MeshMover(lmp,_mesh), omega(2.*M_PI/T)
{
    amplitude[0] = ax;
    amplitude[1] = ay;
    amplitude[2] = az;

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

void MeshMoverWiggle::initial_integrate(double dT,double dt)
{
    double dX[3],dx[3],vNode[3];
    double sine = sin(omega * dT);
    double cosine = cos(omega * dT);

    int size = mesh_->size();
    int numNodes = mesh_->numNodes();
    double ***v_node = get_v();

    // calculate velocity, same for all nodes
    vectorScalarMult3D(amplitude,omega*cosine,vNode);

    // calculate total and incremental displacement
    vectorScalarMult3D(amplitude,sine,dX);
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
                                 double px, double py, double pz,
                                 double axisX, double axisY, double axisZ,
                                 double T)
  : MeshMover(lmp,_mesh), omega(2.*M_PI/T)
{
    axis[0] = axisX;
    axis[1] = axisY;
    axis[2] = axisZ;

    vectorScalarDiv3D(axis,vectorMag3D(axis));

    p[0] = px;
    p[1] = py;
    p[2] = pz;

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

void MeshMoverRotate::initial_integrate(double dT,double dt)
{
    double xOld[3],node[3],vRot[3],omegaVec[3],rPA[3];
    double totalPhi = omega*dT;
    double incrementalPhi = omega*dt;
    double oneOverDT = 1./dT;

    int size = mesh_->size();
    int numNodes = mesh_->numNodes();
    double ***v_node = get_v();
    double ***nodes = get_nodes();

    // rotate the mesh
    mesh_->rotate(totalPhi,incrementalPhi,axis,p);

    // set mesh velocity, w x rPA
    vectorScalarMult3D(axis,omega,omegaVec);
    for(int i = 0; i < size; i++)
    {
      for(int iNode = 0; iNode < numNodes; iNode++)
      {
          vectorCopy3D(nodes[i][iNode],node);
          vectorSubtract3D(node,p,rPA);
          vectorCross3D(omegaVec,rPA,vRot);
          vectorAdd3D(v_node[i][iNode],vRot,v_node[i][iNode]);
      }
    }
}

/* ----------------------------------------------------------------------
   MeshMoverRiggle
------------------------------------------------------------------------- */

MeshMoverRiggle::MeshMoverRiggle(LAMMPS *lmp,AbstractMesh *_mesh,
                                 double px, double py, double pz,
                                 double axisX, double axisY, double axisZ,
                                 double T, double ampl)
  : MeshMover(lmp,_mesh), omega(2.*M_PI/T), amplitude(ampl*M_PI/180.)
{
    axis[0] = axisX;
    axis[1] = axisY;
    axis[2] = axisZ;

    vectorScalarDiv3D(axis,vectorMag3D(axis));

    p[0] = px;
    p[1] = py;
    p[2] = pz;

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

void MeshMoverRiggle::initial_integrate(double dT,double dt)
{
    double xOld[3],node[3],vRot[3],omegaVec[3],rPA[3];

    double oneOverDT = 1./dT;

    double sine = amplitude*sin(omega * dT);
    double cosine = amplitude*cos(omega * dT);

    double vel_prefactor = cosine*omega;

    int size = mesh_->size();
    int numNodes = mesh_->numNodes();
    double ***v_node = get_v();
    double ***nodes = get_nodes();

    // calculate total and incremental angle
    double totalPhi = sine;
    double incrementalPhi = cosine*omega*dt;

    // rotate the mesh
    mesh_->rotate(totalPhi,incrementalPhi,axis,p);

    // set mesh velocity, vel_prefactor * w/|w| x rPA
    vectorScalarMult3D(axis,vel_prefactor,omegaVec);
    for(int i = 0; i < size; i++)
    {
      for(int iNode = 0; iNode < numNodes; iNode++)
      {
          vectorCopy3D(nodes[i][iNode],node);
          vectorSubtract3D(node,p,rPA);
          vectorCross3D(omegaVec,rPA,vRot);
          vectorAdd3D(v_node[i][iNode],vRot,v_node[i][iNode]);
      }
    }
}
