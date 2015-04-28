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
    (if not contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef LMP_TET_MESH_H
#define LMP_TET_MESH_H

#include "volume_mesh.h"

namespace LAMMPS_NS
{

  class TetMesh : public VolumeMesh<4,4,3>
  {
      public:

        TetMesh(LAMMPS *lmp);
        virtual ~TetMesh();

        int generateRandomOwnedGhost(double *pos);
        int generateRandomSubbox(double *pos);
        int generateRandomSubboxWithin(double *pos,double delta);

      protected:

        double calcVol(int n);
        double calcCenter(int n);
        bool isInside(int iTet,double *pos);

        bool shareFace(int i, int j, int &iFace, int &jFace);

      private:

        double calcTetVol(double* v0, double* v1, double* v2, double* v3);
        void baryToCart(int iTet,double *bary_coo,double *pos);
  };

  // *************************************
  #include "tet_mesh_I.h"
  // *************************************

} /* LAMMPS_NS */
#endif /* TRIMESH_H_ */
