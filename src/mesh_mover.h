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

#ifndef LMP_MESH_MOVER_H
#define LMP_MESH_MOVER_H

#include <vector>
#include "tri_mesh.h"
#include "fix_move_mesh.h"
#include "force.h"

namespace LAMMPS_NS
{

  class MeshMover : protected Pointers
  {
      public:

        MeshMover(LAMMPS * lmp,AbstractMesh *_mesh,FixMoveMesh *_fix_move_mesh) :
            Pointers(lmp),
            mesh_(_mesh),
            fix_move_mesh_(_fix_move_mesh),
            isFirst_(false),
            has_reference_point_(false),
            last_reset_(0)
        {
            vectorZeroize3D(reference_point_);
        }

        virtual ~MeshMover()
        {}

        virtual void post_create() = 0;
        virtual void pre_delete() = 0;
        virtual void setup() {};
        virtual void init()
        {}

        virtual void initial_integrate(double dTAbs,double dTSetup,double dt) = 0;
        virtual void final_integrate(double dTAbs,double dTSetup,double dt) {};

        inline bool isFirst()
        { return isFirst_; }

        virtual int n_restart()
        { return 0; }

        virtual void write_restart(double *buf) {}
        virtual void read_restart(double *buf) {}

        void add_reference_point(double *point)
        {
            if (has_reference_point_)
                error->all(FLERR, "Internal error: Mesh mover can only have one reference point");
            vectorCopy3D(point, reference_point_);
            has_reference_point_ = true;
        }

        void get_reference_point(double *point)
        { vectorCopy3D(reference_point_, point); }

        double ***get_nodes()
        { return mesh_->nodePtr(); }

        double ***get_v()
        {
            double ***ptr = NULL;
            if(mesh_->numNodes() == 3)
                ptr = mesh_->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v")->begin();
            else if(mesh_->numNodes() == 4)
                return mesh_->prop().getElementProperty<MultiVectorContainer<double,4,3> >("v")->begin();
            if(!ptr)
                error->one(FLERR,"Illegal call to MeshMover::get_v");
            return ptr;
        }

        virtual void move(const double * const dx)
        {
            if (has_reference_point_)
                vectorAdd3D(reference_point_, dx, reference_point_);
        }

        virtual void rotate(const double dphi, const double * const axis, const double * const center)
        {
            if (has_reference_point_)
            {
                double tmp[3], quat[4];
                vectorSubtract3D(reference_point_, center, tmp);
                // build rotation quaternion
                phiToQuat(dphi, axis, quat);
                // only rotate if center != centerOfRotation
                if (vectorMag3DSquared(tmp) > 1e-20*vectorMag3DSquared(reference_point_))
                {
                    MathExtraLiggghts::vec_quat_rotate(tmp, quat);
                    vectorAdd3D(tmp, center, reference_point_);
                }
            }
        }

        AbstractMesh *mesh_;
        FixMoveMesh *fix_move_mesh_;
        bool isFirst_;
        bool has_reference_point_;
        double reference_point_[3];
        bigint last_reset_;
  };

} /* LAMMPS_NS */

#endif /* MESHMOVER_H_ */
