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

#ifndef LMP_ABSTRACT_MESH_H
#define LMP_ABSTRACT_MESH_H

#include "pointers.h"

namespace LAMMPS_NS
{
  class AbstractMesh : protected Pointers
  {
      friend class FixMoveMesh;
      friend class MeshMover;
      friend class FixMeshSurfaceStressServo;

      public:

        virtual void setMeshID(const char *_mesh_id) = 0;

        virtual void setPrecision(double _precision) = 0;

        virtual void autoRemoveDuplicates() = 0;

        // scale mesh
        virtual void scale(double factor) = 0;

        // linear move w/ total and incremental displacement
        virtual void move(double *vecTotal, double *vecIncremental) = 0;

        // linear move w/ incremental displacement
        virtual void move(double *vecIncremental) = 0;

        // rotation w/ total and incremental displacement
        //   calls rotate(double *totalQuat,double *dQuat,double *displacement)
        virtual void rotate(double totalAngle, double dAngle, double *axis, double *p) = 0;

        // rotation w/ incremental displacement
        //   calls rotate(double *dQuat,double *displacement)
        virtual void rotate(double dAngle, double *axis, double *p) = 0;

        // rotation using quaternions
        virtual void rotate(double *totalQ, double *dQ,double *origin) = 0;

        // initialize movement
        virtual bool registerMove(bool _scale, bool _translate, bool _rotate) = 0;
        virtual void unregisterMove(bool _scale, bool _translate, bool _rotate) = 0;

        virtual bool isMoving() = 0;
        virtual int nMove() = 0;

        // get node j of element i
        
        virtual void node_slow(int i,int j,double *node) = 0;

        // neigh list stuff for moving mesh
        virtual bool decideRebuild() = 0;

        virtual void initalSetup() = 0;
        virtual void pbcExchangeBorders(int setupFlag) = 0;
        virtual void clearReverse() = 0;
        virtual void forwardComm() = 0;
        virtual void reverseComm() = 0;

        virtual void writeRestart(FILE *fp) = 0;
        virtual void restart(double *list) = 0;

        virtual bool allNodesInsideSimulationBox() = 0;

        virtual int numNodes() = 0;

        virtual class CustomValueTracker& prop() = 0;

        /*
        virtual ContainerBase* container(double type,int lenVec) = 0;
        virtual ContainerBase* container(int type,int lenVec) = 0;
        virtual ContainerBase* container(bool type,int lenVec) = 0;
        */
        virtual int id_slow(int i) = 0;

        virtual void setVerbose() = 0;

        virtual void check_element_property_consistency() = 0;

        virtual int nBelowAngle() = 0;
        virtual double angleLimit() = 0;
        virtual int nTooManyNeighs() = 0;

        // size includes owned and ghost elements
        inline int size()
        { return sizeLocal()+sizeGhost(); }

        // virtual functions for size
        // parallelism implemented in derived class
        virtual int sizeLocal() = 0;
        virtual int sizeGhost() = 0;
        virtual int sizeGlobal() = 0;

        virtual ~AbstractMesh()
        {}

      protected:

          AbstractMesh(LAMMPS *lmp)
          : Pointers(lmp)
          {}

      private:

         virtual double*** nodePtr() = 0;
  };

} /* LAMMPS_NS */
#endif
