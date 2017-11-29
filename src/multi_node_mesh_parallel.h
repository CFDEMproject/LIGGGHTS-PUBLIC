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

#ifndef LMP_MULTI_NODE_MESH_PARALLEL_H
#define LMP_MULTI_NODE_MESH_PARALLEL_H

#include "mpi_liggghts.h"
#include "multi_node_mesh.h"
#include "comm.h"
#include "error.h"
#include "vector_liggghts.h"
#include "neighbor.h"
#include "math_extra_liggghts.h"
#include "container_base.h"
#include "domain_wedge.h"
#include <string>
#include <list>
#include <cmath>
#include <algorithm>

namespace LAMMPS_NS
{

  template<int NUM_NODES>
  class MultiNodeMeshParallel : public MultiNodeMesh<NUM_NODES>
  {
      public:

        void initialSetup();
        void pbcExchangeBorders(int setupFlag);
        void clearReverse();
        void forwardComm(std::string property);
        void forwardComm(std::list<std::string> * properties = NULL);
        void reverseComm(std::string property);
        void reverseComm(std::list<std::string> * properties = NULL);

        void writeRestart(FILE *fp);
        void restart(double *list);

        bool allNodesInsideSimulationBox();
        void useAsInsertionMesh(bool parallel);

        inline bool isInsertionMesh()
        { return isInsertionMesh_; }

        inline int sizeLocal()
        { return nLocal_; }

        inline int sizeGhost()
        { return nGhost_; }

        inline int sizeGlobal()
        { return nGlobal_; }

        inline int sizeGlobalOrig()
        { return nGlobalOrig_; }

        inline bool isParallel()
        { return isParallel_; }

        virtual int id(int i) = 0;

      protected:

        MultiNodeMeshParallel(LAMMPS *lmp);
        virtual ~MultiNodeMeshParallel();

        virtual bool addElement(double **nodeToAdd);
        
        virtual bool addElement(double **nodeToAdd,int lineNumb) = 0;
        virtual void deleteElement(int n);

        virtual void buildNeighbours() = 0;

        virtual void refreshOwned(int setupFlag);
        virtual void refreshGhosts(int setupFlag);

        virtual void qualityCheck() = 0;

        virtual void preSetup() {}
        virtual void preInitialSetup() {}
        virtual void postInitialSetup() {}
        virtual void postBorders() {}

        virtual void clearMap() = 0;
        virtual void generateMap() = 0;
        virtual void clearGhostForward(bool scale,bool translate,bool rotate);

        // lo-level parallelization also used by derived classes

        virtual int elemListBufSize(int n,int operation,bool scale,bool translate,bool rotate);
        virtual int pushElemListToBuffer(int n, int *list, int *wraplist, double *buf, int operation, std::list<std::string> * properties, double *dlo, double *dhi, bool scale,bool translate, bool rotate);
        virtual int popElemListFromBuffer(int first, int n, double *buf, int operation, std::list<std::string> * properties, bool scale,bool translate, bool rotate);
        virtual int pushElemListToBufferReverse(int first, int n, double *buf, int operation, std::list<std::string> * properties, bool scale,bool translate, bool rotate);
        virtual int popElemListFromBufferReverse(int n, int *list, double *buf, int operation, std::list<std::string> * properties, bool scale,bool translate, bool rotate);

        virtual int elemBufSize(int operation, std::list<std::string> * properties, bool scale,bool translate,bool rotate);
        virtual int pushElemToBuffer(int i, double *buf,int operation,bool scale,bool translate,bool rotate);
        virtual int popElemFromBuffer(double *buf,int operation,bool scale,bool translate,bool rotate);

        virtual int meshPropsBufSize(int operation,bool scale,bool translate,bool rotate) = 0;
        virtual int pushMeshPropsToBuffer(double *buf, int operation,bool scale,bool translate, bool rotate) = 0;
        virtual int popMeshPropsFromBuffer(double *buf, int operation,bool scale,bool translate, bool rotate) = 0;

        // flags if mesh should be parallelized
        bool doParallellization_;

      private:

        // parallelization functions

        void setup();
        void deleteUnowned();
        void pbc();
        void exchange();
        void borders();
        void clearGhosts();

        int checkBorderElement (const int, const int, const int, const double, const double) const;
        int checkBorderElementLeft (const int, const int, const double, const double) const;
        int checkBorderElementRight(const int, const int, const double, const double) const;

        // lo-level parallelization
        int pushExchange(int dim);
        void popExchange(int nrecv,int dim,double *buf);

        int sizeRestartMesh();
        int sizeRestartElement();

        // number of local and ghost elements
        int nLocal_;
        int nGhost_;
        int nGlobal_;

        // initial number of elements
        int nGlobalOrig_;

        // flags if mesh is parallelized
        bool isParallel_;

        // flag indicating usage as insertion mesh
        bool isInsertionMesh_;

        // *************************************
        // comm stuff - similar to Comm class
        // *************************************

        void grow_swap(int n);
        void allocate_swap(int n);
        void free_swap();
        void grow_send(int,int);
        void grow_recv(int);
        void grow_list(int, int);

        // communication buffers
        int maxsend_,maxrecv_;       // current size of send/recv buffer
        double *buf_send_, *buf_recv_;

        // comm
        int sendneed_[3][2];         // # of procs away I send elements to
        int maxneed_[3];             // max procs away any proc needs, per dim
        double half_atom_cut_;       // half atom neigh cut

        int size_exchange_;          // # of per-element datums in exchange
        int size_forward_;           // # of per-element datums in forward comm
        int size_reverse_;           // # of datums in reverse comm
        int size_border_;            // # of datums in forward border comm

        int maxforward_,maxreverse_; // max # of datums in forward/reverse comm

        // comm swaps
        
        int nswap_;                  // # of swaps to perform = sum of maxneed
        int maxswap_;                // # of swaps where mem was allocated
        int *sendnum_,*recvnum_;     // # of atoms to send/recv in each swap
        int *firstrecv_;             // where to put 1st recv atom in each swap
        int *sendproc_,*recvproc_;   // proc to send/recv to/from at each swap
        int *size_forward_recv_;     // # of values to recv in each forward comm
        int *size_reverse_recv_;     // # to recv in each reverse comm

        double *slablo_,*slabhi_;
        int **sendlist_;             // list of elements to send in each swap
        int **sendwraplist_;         // whether an element needs to be wrapped or not
        int *maxsendlist_;           // max size of send list for each swap

        int *pbc_flag_;              // general flag for sending atoms thru PBC
        int **pbc_;                  // dimension flags for PBC adjustments
  };

  // *************************************
  #include "multi_node_mesh_parallel_I.h"
  #include "multi_node_mesh_parallel_buffer_I.h"
  // *************************************

} /* LAMMPS_NS */
#endif /* MULTINODEMESH_H_ */
