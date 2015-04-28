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

#ifndef LMP_TRI_MESH_PLANAR_H
#define LMP_TRI_MESH_PLANAR_H

#include "tri_mesh.h"
#include "math_extra_liggghts.h"
#include "vector_liggghts.h"

#define MAXTRY_WITHIN 500

namespace LAMMPS_NS
{
  class TriMeshPlanar : public TriMesh
  {
      public:

        TriMeshPlanar(LAMMPS *lmp);
        virtual ~TriMeshPlanar();

        static const int NUM_NODES = 3;

        int generateRandomOwnedGhostWithin(double *pos,double delta);
        bool locatePosition(double *pos,int &triID,double *bary,double &distance);
        bool constructPositionFromBary(int triID,double *bary,double *pos);

      protected:

        void postInitialSetup();

      private:

        void buildEdgeLists();

        // store nearest active edges
        
        VectorContainer<int,2*NUM_NODES> &nearestActiveEdgeID_;
        VectorContainer<int,2*NUM_NODES> &nearestActiveEdgeIndex_;

        ScalarContainer<double> &minActiveEdgeDist_;
  };

  // *************************************
  #include "tri_mesh_planar_I.h"
  // *************************************

} /* LAMMPS_NS */
#endif
