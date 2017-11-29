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
#include "mesh_mover_linear.h"
#include "vector_liggghts.h"
#include "math_extra_liggghts.h"
#include "input.h"
#include "modify.h"
#include "variable.h"
#include "fix_mesh.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   MeshMoverLinear
------------------------------------------------------------------------- */

MeshMoverLinear::MeshMoverLinear(LAMMPS *lmp, AbstractMesh *_mesh, FixMoveMesh *_fix_move_mesh,
                                 const char * const * const arg, const int narg) :
    MeshMover(lmp,_mesh,_fix_move_mesh)
{
    if (narg < 4)
        error->all(FLERR, "Not enough arguments for movement type linear");
    if (narg > 4)
        error->warning(FLERR, "Movement type linear requires only 4 arguments, excess arguments will be ignored");
    vel_[0] = force->numeric(FLERR, arg[1]);
    vel_[1] = force->numeric(FLERR, arg[2]);
    vel_[2] = force->numeric(FLERR, arg[3]);
}

void MeshMoverLinear::post_create()
{
    isFirst_ = mesh_->registerMove(false,true,false);
}

void MeshMoverLinear::pre_delete()
{
    mesh_->unregisterMove(false,true,false);
}

MeshMoverLinear::~MeshMoverLinear()
{ }

/* ---------------------------------------------------------------------- */

void MeshMoverLinear::initial_integrate(double dTAbs,double dTSetup,double dt)
{
    double dx[3];

    int size = mesh_->size();
    int numNodes = mesh_->numNodes();
    double ***v_node = get_v();

    // calculate total and incremental displacement
    vectorScalarMult3D(vel_,dt,dx);

    // apply move
    fix_move_mesh_->fixMesh()->move(dx, fix_move_mesh_);

    // set mesh velocity
    for (int i = 0; i < size; i++)
        for(int j = 0; j < numNodes; j++)
            vectorAdd3D(v_node[i][j],vel_,v_node[i][j]);
}

/* ----------------------------------------------------------------------
   MeshMoverLinearVariable
------------------------------------------------------------------------- */

MeshMoverLinearVariable::MeshMoverLinearVariable(LAMMPS *lmp,AbstractMesh *_mesh, FixMoveMesh *_fix_move_mesh,
                                 const char * const * const arg, const int narg) :
    MeshMover(lmp,_mesh,_fix_move_mesh)
{
    if (narg < 4)
        error->all(FLERR, "Not enough arguments for movement type linear/variable");
    if (narg > 4)
        error->warning(FLERR, "Movement type linear/variable requires only 4 arguments, excess arguments will be ignored");

    int n;
    n = strlen(&arg[1][2]) + 1;
    var1str_ = new char[n];
    strcpy(var1str_,&arg[1][2]);
    myvar1_ = input->variable->find(var1str_);

    n = strlen(&arg[2][2]) + 1;
    var2str_ = new char[n];
    strcpy(var2str_,&arg[2][2]);
    myvar2_ = input->variable->find(var2str_);

    n = strlen(&arg[3][2]) + 1;
    var3str_ = new char[n];
    strcpy(var3str_,&arg[3][2]);
    myvar3_ = input->variable->find(var3str_);

    if (myvar1_ < 0)
        error->all(FLERR,"Variable name 1 for fix move/mesh linear/variable does not exist");
    if (myvar2_ < 0)
        error->all(FLERR,"Variable name 2 for fix move/mesh linear/variable does not exist");
    if (myvar3_ < 0)
        error->all(FLERR,"Variable name 3 for fix move/mesh linear/variable does not exist");

    vel_[0] = 0.; //input->variable->compute_equal(myvar1_);
    vel_[1] = 0.; //input->variable->compute_equal(myvar2_);
    vel_[2] = 0.; //input->variable->compute_equal(myvar3_);
}

void MeshMoverLinearVariable::post_create()
{
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

void MeshMoverLinearVariable::setup()
{
    // check if variable still exists
    myvar1_ = input->variable->find(var1str_);
    myvar2_ = input->variable->find(var2str_);
    myvar3_ = input->variable->find(var3str_);

    if (myvar1_ < 0)
        error->all(FLERR,"Variable name 1 for fix move/mesh linear/variable does not exist");
    if (myvar2_ < 0)
        error->all(FLERR,"Variable name 2 for fix move/mesh linear/variable does not exist");
    if (myvar3_ < 0)
        error->all(FLERR,"Variable name 3 for fix move/mesh linear/variable does not exist");

    vel_[0] = input->variable->compute_equal(myvar1_);
    vel_[1] = input->variable->compute_equal(myvar2_);
    vel_[2] = input->variable->compute_equal(myvar3_);
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

    // apply move
    fix_move_mesh_->fixMesh()->move(dx, fix_move_mesh_);

    // set mesh velocity
    for (int i = 0; i < size; i++)
        for(int j = 0; j < numNodes; j++)
            vectorAdd3D(v_node[i][j],vel_,v_node[i][j]);
}

/* ----------------------------------------------------------------------
   MeshMoverWiggle
------------------------------------------------------------------------- */

MeshMoverWiggle::MeshMoverWiggle(LAMMPS *lmp,AbstractMesh *_mesh, FixMoveMesh *_fix_move_mesh,
                                 const char * const * const arg, const int narg) :
    MeshMover(lmp,_mesh,_fix_move_mesh)
{
    if (narg < 7)
        error->all(FLERR, "Not enough arguments for movement type wiggle");
    if (narg > 7)
        error->warning(FLERR, "Movement type wiggle requires only 7 arguments, excess arguments will be ignored");
    if (strcmp(arg[1], "amplitude"))
        error->all(FLERR, "Expected keyword 'amplitude'");
    amplitude_[0] = force->numeric(FLERR, arg[2]);
    amplitude_[1] = force->numeric(FLERR, arg[3]);
    amplitude_[2] = force->numeric(FLERR, arg[4]);
    if (strcmp(arg[5], "period"))
        error->all(FLERR, "Expected keyword 'period'");
    omega_ = 2.*M_PI/force->numeric(FLERR, arg[6]);
}

void MeshMoverWiggle::post_create()
{
    isFirst_ = mesh_->registerMove(false,true,false);
}

void MeshMoverWiggle::pre_delete()
{
    mesh_->unregisterMove(false,true,false);
}

MeshMoverWiggle::~MeshMoverWiggle()
{ }

/* ---------------------------------------------------------------------- */

void MeshMoverWiggle::initial_integrate(double dTAbs,double dTSetup,double dt)
{

    double dx[3],vNode[3];
    //double cosine = cos(omega_ * dTAbs) - cos(omega_ * (dTAbs-dTSetup));

    int size = mesh_->size();
    int numNodes = mesh_->numNodes();
    double ***v_node = get_v();

    // calculate velocity, same for all nodes
    vectorScalarMult3D(amplitude_,omega_*cos(omega_ * dTAbs),vNode);

    // calculate total and incremental displacement
    vectorScalarMult3D(vNode,dt,dx);

    // apply move
    fix_move_mesh_->fixMesh()->move(dx, fix_move_mesh_);

    // set mesh velocity
    for (int i = 0; i < size; i++)
        for(int j = 0; j < numNodes; j++)
            vectorAdd3D(v_node[i][j],vNode,v_node[i][j]);
}

/* ----------------------------------------------------------------------
   MeshMoverVibLin
------------------------------------------------------------------------- */

MeshMoverVibLin::MeshMoverVibLin(LAMMPS *lmp,AbstractMesh *_mesh, FixMoveMesh *_fix_move_mesh,
                                 const char * const * const arg, const int narg) :
    MeshMover(lmp,_mesh,_fix_move_mesh)
{
    if (narg < 7)
        error->all(FLERR, "Not enough arguments for movement type viblin");
    if (strcmp(arg[5], "order"))
        error->all(FLERR, "Expected keyword 'order'");
    ord = force->inumeric(FLERR, arg[6]);
    if (narg < 10+2*ord)
        error->all(FLERR, "Not enough arguments for movement type viblin");
    if (narg > 10+2*ord)
        error->warning(FLERR, "Movement type wiggle requires only (10 + 2*$order) arguments, excess arguments will be ignored");
    if (ord > 30 || ord < 1)
        error->all(FLERR, "order can be at most 30 and must be greater 0");

    if (strcmp(arg[1], "axis"))
        error->all(FLERR, "Expected keyword 'axis'");
    axis_[0] = force->numeric(FLERR, arg[2]);
    axis_[1] = force->numeric(FLERR, arg[3]);
    axis_[2] = force->numeric(FLERR, arg[4]);
    vectorNormalize3D(axis_);

    if (strcmp(arg[7], "amplitude"))
        error->all(FLERR, "Expected keyword 'amplitude'");
    if (strcmp(arg[8+ord], "phase"))
        error->all(FLERR, "Expected keyword 'phase'");
    if (strcmp(arg[9+2*ord], "period"))
        error->all(FLERR, "Expected keyword 'period'");
    //array transfer
    for (int j=0;j<ord; j++)
    {
       ampl[j]   = force->numeric(FLERR, arg[8+j]);
       phi[j]  = force->numeric(FLERR, arg[9+ord+j]);
       omega[j] = 2.*M_PI/force->numeric(FLERR, arg[10+2*ord+j]);
    }
}

void MeshMoverVibLin::post_create()
{
    isFirst_ = mesh_->registerMove(false,true,false);
}

void MeshMoverVibLin::pre_delete()
{
    mesh_->unregisterMove(false,true,false);
}

MeshMoverVibLin::~MeshMoverVibLin()
{ }

/* ---------------------------------------------------------------------- */
void MeshMoverVibLin::initial_integrate(double dTAbs,double dTSetup,double dt)
{
    double dx[3],vNode[3];
    int size = mesh_->size();
    int numNodes = mesh_->numNodes();
    double ***v_node = get_v();

    double arg = 0;
    double vA = 0;

    for (int j=0;j<ord; j++)
    {
        arg = arg+ampl[j]*(cos(omega[j]*dTAbs+phi[j]) - cos(omega[j]*(dTAbs-dTSetup)+phi[j]));
        vA= vA-ampl[j]*omega[j]*sin(omega[j]*dTAbs+phi[j]);
    }

    // calculate velocity, same for all nodes
    vectorScalarMult3D(axis_,vA,vNode);
    // calculate total and incremental displacement
    vectorScalarMult3D(vNode,dt,dx); //incremental

    // apply linear move
    fix_move_mesh_->fixMesh()->move(dx, fix_move_mesh_);

    // set mesh velocity
    for (int i = 0; i < size; i++)
        for(int j = 0; j < numNodes; j++)
            vectorAdd3D(v_node[i][j],vNode,v_node[i][j]);

}
