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

#ifdef MESHMOVER_CLASS

MeshMoverStyle(rotate,MeshMoverRotate)
MeshMoverStyle(rotate/variable,MeshMoverRotateVariable)
MeshMoverStyle(riggle,MeshMoverRiggle)
MeshMoverStyle(vibrot,MeshMoverVibRot)

#else

#ifndef LMP_MESH_MOVER_ROTATION_H
#define LMP_MESH_MOVER_ROTATION_H

#include <vector>
#include "tri_mesh.h"
#include "fix_move_mesh.h"
#include "mesh_mover.h"
#include "force.h"

namespace LAMMPS_NS
{

/* ---------------------------------------------------------------------- */

class MeshMoverRotate : public MeshMover {

public:

    MeshMoverRotate(LAMMPS *lmp,AbstractMesh *_mesh,FixMoveMesh *_fix_move_mesh,
                    const char * const * const arg, const int narg);
    virtual ~MeshMoverRotate();

    void initial_integrate(double dTAbs,double dTSetup,double dt);
    void final_integrate(double dTAbs,double dTSetup,double dt) {}
    void pre_delete();
    void post_create();

private:

    double axis_[3], point_[3], omega_;
};

/* ---------------------------------------------------------------------- */

class MeshMoverRotateVariable : public MeshMover {

public:

    MeshMoverRotateVariable(LAMMPS *lmp,AbstractMesh *_mesh,FixMoveMesh *_fix_move_mesh,
                            const char * const * const arg, const int narg);
    virtual ~MeshMoverRotateVariable();

    void pre_delete();
    void post_create();
    void setup();

    void initial_integrate(double dTAbs,double dTSetup,double dt);
    void final_integrate(double dTAbs,double dTSetup,double dt) {}

private:

    char *var1str_;
    int myvar1_;
    double axis_[3], point_[3], omega_, totalPhi_;
};

/* ---------------------------------------------------------------------- */

class MeshMoverRiggle : public MeshMover {

public:

    MeshMoverRiggle(LAMMPS *lmp,AbstractMesh *_mesh,FixMoveMesh *_fix_move_mesh,
                    const char * const * const arg, const int narg);
    virtual ~MeshMoverRiggle();

    void initial_integrate(double dTAbs,double dTSetup,double dt);
    void final_integrate(double dTAbs,double dTSetup,double dt) {}
    void pre_delete();
    void post_create();

private:

    double axis_[3], point_[3], omega_, amplitude_;
};

/* ---------------------------------------------------------------------- */

class MeshMoverVibRot : public MeshMover {

public:
    MeshMoverVibRot(LAMMPS *lmp,AbstractMesh *_mesh,FixMoveMesh *_fix_move_mesh,
                    const char * const * const arg, const int narg);
    virtual ~MeshMoverVibRot();

    void initial_integrate(double dTAbs,double dTSetup,double dt);
    void final_integrate(double dTAbs,double dTSetup,double dt) {}
    void pre_delete();
    void post_create();

private:
    double axis_[3], ampl[30], phi[30], p_[3], omega[30];
    int ord;

};

} /* LAMMPS_NS */

#endif /* MESHMOVER_ROTATION_H_ */

#endif
