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

MeshMoverStyle(linear,MeshMoverLinear)
MeshMoverStyle(linear/variable,MeshMoverLinearVariable)
MeshMoverStyle(wiggle,MeshMoverWiggle)
MeshMoverStyle(viblin,MeshMoverVibLin)

#else

#ifndef LMP_MESH_MOVER_LINEAR_H
#define LMP_MESH_MOVER_LINEAR_H

#include <vector>
#include "tri_mesh.h"
#include "fix_move_mesh.h"
#include "mesh_mover.h"
#include "force.h"

namespace LAMMPS_NS
{

/* ---------------------------------------------------------------------- */

class MeshMoverLinear : public MeshMover{

public:

    MeshMoverLinear(LAMMPS *lmp,AbstractMesh *_mesh,FixMoveMesh *_fix_move_mesh,
                    const char * const * const arg, const int narg);
    virtual ~MeshMoverLinear();

    void initial_integrate(double dTAbs,double dTSetup,double dt);
    void final_integrate(double dTAbs,double dTSetup,double dt) {}
    void pre_delete();
    void post_create();

private:

    double vel_[3];
};

/* ----------------------------------------------------------------------- */

class MeshMoverLinearVariable : public MeshMover{

public:

    MeshMoverLinearVariable(LAMMPS *lmp,AbstractMesh *_mesh,FixMoveMesh *_fix_move_mesh,
                            const char * const * const arg, const int narg);
    virtual ~MeshMoverLinearVariable();

    void pre_delete();
    void post_create();
    void setup();

    void initial_integrate(double dTAbs,double dTSetup,double dt);
    void final_integrate(double dTAbs,double dTSetup,double dt) {}

private:

    char *var1str_,*var2str_,*var3str_;
    int myvar1_,myvar2_,myvar3_;
    double vel_[3];
};

/* ---------------------------------------------------------------------- */

class MeshMoverWiggle : public MeshMover{

public:

    MeshMoverWiggle(LAMMPS *lmp,AbstractMesh *_mesh,FixMoveMesh *_fix_move_mesh,
                    const char * const * const arg, const int narg);
    virtual ~MeshMoverWiggle();

    void initial_integrate(double dTAbs,double dTSetup,double dt);
    void final_integrate(double dTAbs,double dTSetup,double dt) {}
    void pre_delete();
    void post_create();

private:

    double amplitude_[3],omega_;
};

/* ---------------------------------------------------------------------- */
class MeshMoverVibLin : public MeshMover {

public:
    MeshMoverVibLin(LAMMPS *lmp,AbstractMesh *_mesh,FixMoveMesh *_fix_move_mesh,
                    const char * const * const arg, const int narg);
    virtual ~MeshMoverVibLin();

    void initial_integrate(double dTAbs,double dTSetup,double dt);
    void final_integrate(double dTAbs,double dTSetup,double dt) {}
    void pre_delete();
    void post_create();

private:
    double axis_[3], omega[30], ampl[30], phi[30];
    int ord;

};

} /* LAMMPS_NS */

#endif /* MESHMOVER_LINEAR_H_ */

#endif
