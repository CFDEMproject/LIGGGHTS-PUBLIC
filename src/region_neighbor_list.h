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

#include <vector>
#include "vector_liggghts.h"
#include "pointers.h"
#include "bounding_box.h"
#include "superquadric_flag.h"

namespace LAMMPS_NS {

/**
 * @brief A small particle structure
 */
struct Particle {
  int index;
  double x[3];
  double radius;
#ifdef SUPERQUADRIC_ACTIVE_FLAG
  double shape[3];
  double quaternion[4];
#endif

  Particle(int i,double * pos, double rad) {
    index = i;
    LAMMPS_NS::vectorCopy3D(pos, x);
    radius = rad;
#ifdef SUPERQUADRIC_ACTIVE_FLAG
    quaternion[0] = 1.0;
    quaternion[1] = quaternion[2] = quaternion[3] = 0.0;
    shape[0] = shape[1] = shape[2] = radius;
#endif
  }
  Particle(double * pos, double rad) {
    index = -1;
    LAMMPS_NS::vectorCopy3D(pos, x);
    radius = rad;
#ifdef SUPERQUADRIC_ACTIVE_FLAG
    quaternion[0] = 1.0;
    quaternion[1] = quaternion[2] = quaternion[3] = 0.0;
    shape[0] = shape[1] = shape[2] = radius;
#endif
  }
#ifdef SUPERQUADRIC_ACTIVE_FLAG
  Particle(double * pos, double rad, double *quaternion_, double *shape_) {
      index = -1;
      LAMMPS_NS::vectorCopy3D(pos, x);
      radius = rad;
      LAMMPS_NS::vectorCopy4D(quaternion_, quaternion);
      LAMMPS_NS::vectorCopy3D(shape_, shape);
    }
#endif
};

typedef std::vector<Particle> ParticleBin;

struct Bin {
    double center[3];
    ParticleBin p_array;
};

/**
 * @brief A neighbor list of of a certain region
 *
 * This implementation uses the same binning approach used in LAMMPS neighbor lists.
 * Instead of accessing internal data structures directly, manipulations and queries
 * can only occur through the given interface.
 * This class is based on LOCAL binning like in the Neighbor class, i.e. each proc
 * only allocates bins for his own sub-box
 */

class RegionNeighborList : protected LAMMPS_NS::Pointers
{
public:
    RegionNeighborList(LAMMPS_NS::LAMMPS *lmp);

    bool hasOverlap(double * x, double radius) const;
    bool hasOverlapWith(double * x, double radius, std::vector<int> &overlap_list) const ;
    void insert(double * x, double radius,int index = -1);
#ifdef SUPERQUADRIC_ACTIVE_FLAG
    bool hasOverlap_superquadric(double * x, double radius, double *quaternion, double *shape) const;
    void insert_superquadric(double * x, double radius, double *quaternion, double *shape);
    void set_obb_flag(int check_obb_flag_) {check_obb_flag = check_obb_flag_;}
#endif

    size_t count() const;
    void clear();
    virtual bool setBoundingBox(LAMMPS_NS::BoundingBox & bb, double maxrad, bool extend = true, bool failsafe = false);
    bool isInBoundingBox(double *pos);

    bool boundingBoxSet()
    {return bbox_set; }

    int coord2binGlobal(double *x) const;
    std::vector<Bin> getBins() const
    {return bins; }

    double invBinVolume() const
    { return bininvx*bininvy*bininvz; }

protected:

    std::vector<Bin> bins;          // list of particle bins
    std::vector<int> stencil;       // stencil used to check bins for collisions
    size_t ncount;                  // total number of particles in neighbor list

    bool bbox_set;

    double bboxlo[3];               // lowest point of bounding box
    double bboxhi[3];               // highest point of bounding box

    int nbinx,nbiny,nbinz;          // # of global bins
    int mbinx,mbiny,mbinz;          // # of bins in each dimension
    int mbinxlo,mbinylo,mbinzlo;    // offsets of (0,0,0) bin

    double binsizex,binsizey,binsizez;  // bin sizes
    double bininvx,bininvy,bininvz;     // inverse of bin sizes

    double bin_distance(int i, int j, int k);
    int coord2binLocal(double *x) const;

#ifdef SUPERQUADRIC_ACTIVE_FLAG
  int check_obb_flag;
#endif
};

}

#endif // REGION_NEIGHBOR_LIST_H
