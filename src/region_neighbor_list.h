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

    Philippe Richard Berger (JKU Linz)

    Copyright 2014-2015 JKU Linz
------------------------------------------------------------------------- */

#ifndef REGION_NEIGHBOR_LIST_H
#define REGION_NEIGHBOR_LIST_H

#include <vector>
#include "vector_liggghts.h"
#include "pointers.h"
#include "bounding_box.h"

namespace LAMMPS_NS {

/**
 * @brief A small particle structure
 */
struct Particle {
  double x[3];
  double radius;

  Particle(double * pos, double rad) {
    LAMMPS_NS::vectorCopy3D(pos, x);
    radius = rad;
  }
};

/**
 * @brief A neighbor list of of a certain region
 *
 * This implementation uses the same binning approach used in LAMMPS neighbor lists.
 * Instead of accessing internal data structures directly, manipulations and queries
 * can only occur through the given interface.
 */
class RegionNeighborList : protected LAMMPS_NS::Pointers
{
  typedef std::vector<Particle> ParticleBin;

  std::vector<ParticleBin> bins;  // list of particle bins
  std::vector<int> stencil;       // stencil used to check bins for collisions
  size_t ncount;                  // total number of particles in neighbor list

  double bboxlo[3];               // lowest point of bounding box
  double bboxhi[3];               // highest point of bounding box

  int nbinx,nbiny,nbinz;          // # of global bins
  int mbinx,mbiny,mbinz;          // # of bins in each dimension
  int mbinxlo,mbinylo,mbinzlo;    // offsets of (0,0,0) bin

  double binsizex,binsizey,binsizez;  // bin sizes
  double bininvx,bininvy,bininvz;     // inverse of bin sizes

  double bin_distance(int i, int j, int k);
  int coord2bin(double *x) const;

public:
    RegionNeighborList(LAMMPS_NS::LAMMPS *lmp);

    bool hasOverlap(double * x, double radius) const;
    bool hasOverlapWith(double * x, double radius, std::vector<int> &overlap_list) const ;
    void insert(double * x, double radius);
    size_t count() const;
    void reset();
    bool setBoundingBox(LAMMPS_NS::BoundingBox & bb, double maxrad);
};

}

#endif // REGION_NEIGHBOR_LIST_H
