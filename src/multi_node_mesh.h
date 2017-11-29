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

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef LMP_MULTI_NODE_MESH_H
#define LMP_MULTI_NODE_MESH_H

#include "abstract_mesh.h"
#include "domain.h"
#include "vector_liggghts.h"
#include "math_extra_liggghts.h"
#include "update.h"
#include "container.h"
#include "bounding_box.h"
#include "random_park.h"
#include <stdio.h>
#include <cmath>
#include <algorithm>

#define EPSILON_PRECISION 1e-8

namespace LAMMPS_NS
{
  template<int NUM_NODES>
  class MultiNodeMesh : public AbstractMesh
  {
      friend class FixMoveMesh;
      friend class MeshMover;

      public:

        void setMeshID(const char *_mesh_id);

        void setPrecision(double _precision);

        void setMinFeatureLength(double _min_feature_length);

        void setElementExclusionList(FILE *_file);

        void autoRemoveDuplicates();

        // scale mesh
        virtual void scale(double factor);

        // linear move w/ total and incremental displacement
        virtual void move(const double * const vecTotal, const double * const vecIncremental);

        // linear move w/ incremental displacement
        virtual void move(const double * const vecIncremental);

        // rotation w/ total and incremental displacement
        //   calls rotate(double *totalQuat,double *dQuat,double *displacement)
        void rotate(const double totalAngle, const double dAngle, const double * const axis, const double * const p);

        // rotation w/ incremental displacement
        //   calls rotate(double *dQuat,double *displacement)
        void rotate(const double dAngle, const double * const axis, const double * const p);

        void updateCenterRbound(int ilo,int ihi);

        // initialize movement
        bool registerMove(bool _scale, bool _translate, bool _rotate);
        void unregisterMove(bool _scale, bool _translate, bool _rotate);

        // flags for saving global linear and angular velocity
        void set_store_vel()
        { store_vel++; }
        void set_store_omega()
        { store_omega++; }
        void unset_store_vel()
        { store_vel = store_vel > 1 ? store_vel-1 : 0; }
        void unset_store_omega()
        { store_omega = store_omega > 1 ? store_omega-1 : 0; }
        bool is_set_store_vel()
        { return store_vel ? true : false; }
        bool is_set_store_omega()
        { return store_omega ? true : false; }

        // functions to obtain global linear and angular velocity
        void get_global_vel(double * vel);
        void get_global_omega(double * omega);

        // bbox stuff
        BoundingBox getGlobalBoundingBox() const;
        BoundingBox getElementBoundingBoxOnSubdomain(int const n);
        void updateGlobalBoundingBox();

        // neigh list stuff for moving mesh
        bool decideRebuild();
        void storeNodePosRebuild();

        // function for mesh import debugging
        void center_of_mass(double *_com);

        // inline access

        inline bool isMoving()
        { return nMove_ > 0; }

        inline int nMove()
        { return nMove_; }

        inline bool isScaling()
        { return nScale_ > 0; }

        inline bool isTranslating()
        { return nTranslate_ > 0; }

        inline bool isRotating()
        { return nRotate_ > 0; }

        inline void node(int i,int j,double *node)
        { vectorCopy3D(node_(i)[j],node);}

        void node_slow(int i,int j,double *node)
        { vectorCopy3D(node_(i)[j],node);}

        inline void center(int i,double *center)
        { vectorCopy3D(center_(i),center);}

        inline int numNodes()
        { return NUM_NODES; }

        inline char* mesh_id()
        { return mesh_id_; }

        inline bool removeDuplicates()
        { return autoRemoveDuplicates_; }

        // virtual functions for size
        // parallelism implemented in derived class
        virtual int sizeLocal() = 0;
        virtual int sizeGhost() = 0;
        virtual int sizeGlobal() = 0;

        virtual bool isDeforming()
        { return false; }

      protected:
        MultiNodeMesh(LAMMPS *lmp);
        virtual ~MultiNodeMesh();

        virtual bool addElement(double **nodeToAdd);
        virtual void deleteElement(int n);

        virtual void refreshOwned(int setupFlag);
        virtual void refreshGhosts(int setupFlag);

        bool nodesAreEqual(int iSurf, int iNode, int jSurf, int jNode);
        bool nodesAreEqual(double *nodeToCheck1,double *nodeToCheck2);

        // returns true if surfaces share 2 nodes
        // called with local element indices
        // returns indices of nodes with iNode1 < iNode2
        bool share2Nodes(int iElem, int jElem, int &iNode1, int &jNode1, int &iNode2, int &jNode2);

        // returns number of shared nodes
        // called with local index
        int nSharedNodes(int iElem, int jElem);

        // returns node index if iElem contains nodeToCheck
        int containsNode(int iElem, double *nodeToCheck);

        void extendToElem(int const nElem) const;

        // linear move of single element w/ incremental displacement
        virtual void moveElement(const int i, const double * const vecIncremental);

        // rotation using quaternions
        
        virtual void rotate(const double * const totalQ, const double * const dQ, const double * const origin);
        virtual void rotate(const double * const dQ, const double * const origin);

        // mesh nodes
        MultiVectorContainer<double,NUM_NODES,3> node_;

        // original mesh node_ position, used for moving mesh
        MultiVectorContainer<double,NUM_NODES,3> *node_orig_;
        double** node_orig(int i) {return (*node_orig_)(i);}

        // node pos stored by external trigger at neigh build
        MultiVectorContainer<double,NUM_NODES,3> nodesLastRe_;

        // mesh center
        VectorContainer<double,3> center_; 
        ScalarContainer<double> rBound_; 

        // global bounding box for mesh across all processors
        BoundingBox bbox_;

        // random generator
        RanPark *random_;

        // mesh ID - same as fix mesh ID
        char *mesh_id_;

        // for overlap check on element insertion
        
        inline void reset_stepLastReset()
        { stepLastReset_ = -1; }

        // reset mesh nodes to original position
        // called via mesh move functions
        
        virtual bool resetToOrig();

        inline double precision()
        { return precision_; }

        inline double minFeatureLength()
        { return min_feature_length_; }

        inline FILE* elementExclusionList()
        { return element_exclusion_list_; }

      private:

        // mesh precision
        double precision_;

        // ignore features smaller than this size
        double min_feature_length_;

        FILE *element_exclusion_list_;

        // state if elements should be automatically removed if duplicate
        bool autoRemoveDuplicates_;

        // flags stating how many move operations are performed on the mesh
        int nMove_;
        int nScale_,nTranslate_,nRotate_;

        // counters to decide whether mesh linear and angular velocity are stored
        int store_vel, store_omega;

        // step when velocities where last stored
        int step_store_vel, step_store_omega;

        // storage for global mesh linear and angular velocity
        double global_vel[3], global_quaternion[4], prev_quaternion[4];

        // store current node position for use by moving mesh
        void storeNodePosOrig(int ilo, int ihi);

        // step when nodes have been reset the last time
        // only relevant for moving mesh
        int stepLastReset_;

        // extends a given bbox to include element number nElem
        void extendToElem(BoundingBox &box, int const nElem);
        void extendToElem(int const nElem);

        inline double*** nodePtr()
        { return node_.begin(); }
  };

  // *************************************
  #include "multi_node_mesh_I.h"
  // *************************************

} /* LAMMPS_NS */
#endif /* MULTINODEMESH_H_ */
