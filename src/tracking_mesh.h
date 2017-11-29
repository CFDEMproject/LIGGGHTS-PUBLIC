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

/* ----------------------------------------------------------------------
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Philippe Seil (JKU Linz)
------------------------------------------------------------------------- */

#ifndef LMP_TRACKING_MESH_H
#define LMP_TRACKING_MESH_H

#include "multi_node_mesh_parallel.h"
#include "custom_value_tracker.h"
#include "container.h"
#include "memory.h"
#include <map>
#include <vector>
#include <list>

namespace LAMMPS_NS{

  template<int NUM_NODES>
  class TrackingMesh : public MultiNodeMeshParallel<NUM_NODES>
  {
      public:

        void clearReverse();

        void setVerbose();

        virtual void buildNeighbours() = 0;

        virtual void move(const double * const vecTotal, const double * const vecIncremental);
        virtual void move(const double * const vecIncremental);

        virtual void scale(double factor);

        virtual int generateRandomOwnedGhost(double *pos) = 0;
        virtual int generateRandomOwnedGhostWithin(double *pos,double delta) = 0;
        virtual int generateRandomSubbox(double *pos) = 0;

        // inline access

        inline CustomValueTracker& prop()
        {return customValues_;}

        // global to local ID mapping
        inline int map(const int global, const int j)
        {
            if (mapArray_.empty() || mapArray_.find(global) == mapArray_.end())
                return -1;
            else
            {
                if (j < (int)mapArray_[global].size())
                    return mapArray_[global][j];
                else
                    return -1;
            }
        }

        inline int map_size(const int global)
        {
            if (mapArray_.empty() || mapArray_.find(global) == mapArray_.end())
                return 0;
            else
                return mapArray_[global].size();
        }

        //AM-TODO check if still in use
        inline int tag_max()
        { return mapTagMax_; }

        int id_slow(int i)
        { return id_(i); }

        inline int id(int i)
        { return id_(i); }

        inline int lineNo(int i)
        { return (lineNo_?(*lineNo_)(i):-1); }

        inline bool verbose()
        { return verbose_; }

        inline void check_element_property_consistency()
        { customValues_.check_element_property_consistency(this->sizeLocal()+this->sizeGhost()); }

      protected:

        TrackingMesh(LAMMPS *lmp);
        virtual ~TrackingMesh();

        virtual bool addElement(double **nodeToAdd,int lineNumb);
        virtual void deleteElement(int n);

        void clearGhostForward(bool scale,bool translate,bool rotate);

        // called via mesh move functions
        bool resetToOrig();

        void postInitialSetup();

        virtual void refreshOwned(int setupFlag);
        virtual void refreshGhosts(int setupFlag);

        void clearMap();
        void generateMap();

        virtual void moveElement(const int i, const double * const vecIncremental);

        virtual void rotate(const double * const totalQ, const double * const dQ, const double * const origin);
        virtual void rotate(const double * const dQ, const double * const origin);

        // buffer operations

        inline int elemListBufSize(int n,int operation,bool scale,bool translate,bool rotate);
        inline int pushElemListToBuffer(int n, int *list, int *wraplist, double *buf, int operation, std::list<std::string> * properties, double *dlo, double *dhi, bool scale,bool translate, bool rotate);
        inline int popElemListFromBuffer(int first, int n, double *buf, int operation, std::list<std::string> * properties, bool scale,bool translate, bool rotate);
        inline int pushElemListToBufferReverse(int first, int n, double *buf, int operation, std::list<std::string> * properties, bool scale,bool translate, bool rotate);
        inline int popElemListFromBufferReverse(int n, int *list, double *buf, int operation, std::list<std::string> * properties, bool scale,bool translate, bool rotate);

        inline int elemBufSize(int operation, std::list<std::string> * properties, bool scale,bool translate,bool rotate);
        inline int pushElemToBuffer(int n, double *buf, int operation,bool scale,bool translate,bool rotate);
        inline int popElemFromBuffer(double *buf, int operation,bool scale,bool translate,bool rotate);

        int meshPropsBufSize(int operation,bool scale,bool translate,bool rotate);
        int pushMeshPropsToBuffer(double *buf, int operation,bool scale,bool translate, bool rotate);
        int popMeshPropsFromBuffer(double *buf, int operation,bool scale,bool translate, bool rotate);

      private:

        // class holding fields
        CustomValueTracker &customValues_;

        // ID of element
        ScalarContainer<int> &id_;

        // line where element was read
        
        ScalarContainer<int> *lineNo_;

        // global-local lookup from ID to local index
        int mapTagMax_;
        std::map<int, std::vector<int> > mapArray_;

        bool verbose_;

  };

  // *************************************
  #include "tracking_mesh_I.h"
  // *************************************

} /* LAMMPS_NS */
#endif /* TRACKINGMESH_H_ */
