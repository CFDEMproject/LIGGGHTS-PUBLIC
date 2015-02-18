/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
   the producer of the LIGGGHTS® software and the CFDEM®coupling software
   See http://www.cfdem.com/terms-trademark-policy for details.

   LIGGGHTS® is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
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
