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

        virtual void setElementExclusionList(FILE *_file) = 0;

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
