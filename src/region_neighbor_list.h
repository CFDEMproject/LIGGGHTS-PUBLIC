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

    Richard Berger (JKU Linz)
    Christoph Kloss (DCS Computing GmbH)
    Alexander Podlozhnyuk (DCS Computing GmbH)

    Copyright 2014-2015 JKU Linz
    Copyright 2015-     DCS Computing GmbH
------------------------------------------------------------------------- */

#ifndef REGION_NEIGHBOR_LIST_H
#define REGION_NEIGHBOR_LIST_H

#include "region_neighbor_list_base.h"
#include "region_neighbor_list_definitions.h"
#include "pointers.h"
#include "vector_liggghts.h"
#include "bounding_box.h"
#include "lmptype.h"
#include "bounding_box.h"
#include "error.h"
#include "region.h"

#include <vector>
#include <mpi.h>
#include <limits>
#include <algorithm>

#ifdef SUPERQUADRIC_ACTIVE_FLAG
#include "math_extra_liggghts_superquadric.h"
#endif

namespace LAMMPS_NS {

/**
 * @brief A neighbor list of of a certain region
 *
 * This implementation uses the same binning approach used in LAMMPS neighbor lists.
 * Instead of accessing internal data structures directly, manipulations and queries
 * can only occur through the given interface.
 * This class is based on LOCAL binning like in the Neighbor class, i.e. each proc
 * only allocates bins for his own sub-box
 */

template<bool INTERPOLATE>
class RegionNeighborList : public IRegionNeighborList, protected Pointers
{
  friend class FixAddforceSteadystate;
  friend class FixAddforceSteadystateExperimental;

  public:
    RegionNeighborList(LAMMPS_NS::LAMMPS *lmp);
    virtual ~RegionNeighborList() {}

    bool hasOverlap(double * x, double radius) const;
    bool hasOverlapWith(double * x, double radius, std::vector<int> &overlap_list) const;
    void insert(double * x, double radius,int index = -1);
#ifdef SUPERQUADRIC_ACTIVE_FLAG
    bool hasOverlap_superquadric(double * x, double radius, double *quaternion, double *shape, double *blockiness) const;
    void insert_superquadric(double * x, double radius, double *quaternion, double *shape, double *blockiness, int index = -1);
    void set_obb_flag(int check_obb_flag_) {check_obb_flag = check_obb_flag_;}
#endif

    size_t count() const;
    virtual void clear();
    virtual void reset();

    int getSizeOne() const;
    int pushBinToBuffer(int i, double *buf) const;

    inline void setBoundingBox_calc_interpolation_stencil(Bin<INTERPOLATE> &it,int ibin,int ix,int iy, int iz) const;

    virtual bool setBoundingBox(LAMMPS_NS::BoundingBox & bb, double maxrad, bool extend = true, bool failsafe = false);
    virtual BoundingBox setBoundingBoxRegion(const Region &region, double maxrad, bool extend = true, bool failsafe = false);

    bool isInBoundingBox(double *pos) const;

    inline void coord2bin_calc_interpolation_weights(double *x,int ibin,int ix,int iy, int iz,int &quadrant,double &wx,double &wy,double &wz) const;

    int coord2bin(double *x,int &quadrant,double &wx,double &wy,double &wz) const;

    inline int coord2bin(double *x) const
    { int quadrant; double wx,wy,wz; return coord2bin(x,quadrant,wx,wy,wz); }

    bool boundingBoxSet() const
    {return bbox_set; }

    std::vector<LAMMPS_NS::Bin<INTERPOLATE> > getBins() const
    {return bins; }

    double invBinVolume() const
    { return bininvx*bininvy*bininvz; }

    inline int nbins() const
    { return nbinx*nbiny*nbinz; }

    inline int mbins() const
    { return mbinx*mbiny*mbinz; }

    inline void getDimensions(int *dims) const
    { dims[0] = mbinx; dims[1] = mbiny; dims[2] = mbinz; }

    inline void getLocalDimensions(int *dims) const
    {
        // NP TODO make this on local subproc
        dims[0] = 0;
        dims[1] = mbinx;
        dims[2] = 0;
        dims[3] = mbiny;
        dims[4] = 0;
        dims[5] = mbinz;
    }

    inline void getBinSize(double *binsize) const
    { binsize[0] = binsizex; binsize[1] = binsizey; binsize[2] = binsizez; }

    inline void getOrigin(double *origin) const
    { vectorCopy3D(bboxlo, origin); }

  protected:

    std::vector<Bin<INTERPOLATE> > bins;// list of particle bins
    std::vector<int> stencil;           // stencil used to check bins for collisions
    size_t ncount;                      // total number of particles in neighbor list

    bool bbox_set;

    double bboxlo[3];               // lowest point of bounding box
    double bboxhi[3];               // highest point of bounding box

    int nbinx,nbiny,nbinz;          // # of global bins
    int mbinx,mbiny,mbinz;          // # of bins in each dimension
    int mbinxlo,mbinylo,mbinzlo;    // offsets of (0,0,0) bin

    double binsizex,binsizey,binsizez;  // bin sizes
    double bininvx,bininvy,bininvz;     // inverse of bin sizes

    double bin_distance(int i, int j, int k);

#ifdef SUPERQUADRIC_ACTIVE_FLAG
  int check_obb_flag;
#endif
};

/*
INCLUDE INLINE HEADER FILE
*/
#include "region_neighbor_list_I.h"

}

#endif // REGION_NEIGHBOR_LIST_H
