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

#include "lmptype.h"
#include "mpi.h"
#include "bounding_box.h"
#include "error.h"
#include "region_neighbor_list.h"
#include <limits>
#include <algorithm>
#ifdef SUPERQUADRIC_ACTIVE_FLAG
#include "math_extra_liggghts_superquadric.h"
#endif

static const double SMALL = 1.0e-6;
static const double BIG = 1.0e20;

using namespace LAMMPS_NS;
using namespace std;

/**
 * @brief Default constructor which will create an empty neighbor list
 */
RegionNeighborList::RegionNeighborList(LAMMPS *lmp) :
    Pointers(lmp),
    bbox_set(false)
{
}

/**
 * @brief Determine if the given particle overlaps with any particle in this neighbor list
 * @param x        position of particle to check
 * @param radius   radius of particle to check
 * @return true if particle has an overlap with a particle in this neighbor list, false otherwise
 */
bool RegionNeighborList::hasOverlap(double * x, double radius) const {
  int ibin = coord2binLocal(x);

  for(std::vector<int>::const_iterator it = stencil.begin(); it != stencil.end(); ++it) {
    const int offset = *it;
    if((ibin+offset < 0) || ((size_t)(ibin+offset) >= bins.size()))
    {
        
        error->one(FLERR,"assertion failed");
    }
    const ParticleBin & bin = bins[ibin+offset].p_array;

    for(ParticleBin::const_iterator pit = bin.begin(); pit != bin.end(); ++pit) {
      const Particle & p = *pit;
      double del[3];
      vectorSubtract3D(x, p.x, del);
      const double rsq = vectorMag3DSquared(del);
      const double radsum = radius + p.radius;
      if (rsq <= radsum*radsum) return true;
    }
  }

  return false;
}

#ifdef SUPERQUADRIC_ACTIVE_FLAG
//the same for superquadrics
bool RegionNeighborList::hasOverlap_superquadric(double * x, double radius, double *quaternion, double *shape) const {
  int ibin = coord2bin(x);

  for(std::vector<int>::const_iterator it = stencil.begin(); it != stencil.end(); ++it) {
    const int offset = *it;
    if((ibin+offset < 0) || ((size_t)(ibin+offset) >= bins.size()))
    {

        error->one(FLERR,"assertion failed");
    }
    const ParticleBin & bin = bins[ibin+offset];

    double roundness[2] = {2.0, 2.0};
    Superquadric particle1(x, quaternion, shape, roundness);
    for(ParticleBin::const_iterator pit = bin.begin(); pit != bin.end(); ++pit) {
      const Particle & p = *pit;
      double del[3];
      vectorSubtract3D(x, p.x, del);
      const double rsq = vectorMag3DSquared(del);
      const double radsum = radius + p.radius;
      if(check_obb_flag) {
        double x_copy[3], quaternion_copy[4], shape_copy[3];
        vectorCopy3D(p.x, x_copy);
        vectorCopy4D(p.quaternion, quaternion_copy);
        vectorCopy3D(p.shape, shape_copy);
        Superquadric particle2(x_copy, quaternion_copy, shape_copy, roundness);

        if (rsq <= radsum*radsum and MathExtraLiggghtsSuperquadric::obb_intersect(&particle1, &particle2)) return true;
      } else
        if (rsq <= radsum*radsum) return true;
    }
  }

  return false;
}

#endif
/**
 * @brief Determine if the given particle overlaps with any particle in this neighbor list
 * @param x        position of particle to check
 * @param radius   radius of particle to check
 * @param overlap_list  list of overlaps, to be populated by this function
 * @return true if particle has an overlap with a particle in this neighbor list, false otherwise
 */

bool RegionNeighborList::hasOverlapWith(double * x, double radius, std::vector<int> &overlap_list) const {
  int ibin = coord2binLocal(x);

  bool overlap = false;

  for(std::vector<int>::const_iterator it = stencil.begin(); it != stencil.end(); ++it) {
    const int offset = *it;
    if((ibin+offset < 0) || ((size_t)(ibin+offset) >= bins.size()))
    {
        
        error->one(FLERR,"assertion failed");
    }
    const ParticleBin & bin = bins[ibin+offset].p_array;

    for(ParticleBin::const_iterator pit = bin.begin(); pit != bin.end(); ++pit) {
      const Particle & p = *pit;
      double del[3];
      vectorSubtract3D(x, p.x, del);
      const double rsq = vectorMag3DSquared(del);
      const double radsum = radius + p.radius;
      if (rsq <= radsum*radsum)
      {
        overlap_list.push_back(p.index);
        overlap = true;
      }
    }
  }

  return overlap;
}

/**
 * @brief Insert a new particle into neighbor list
 * @param x        position in 3D
 * @param radius   particle radius
 */
void RegionNeighborList::insert(double * x, double radius,int index) {
  int ibin = coord2binLocal(x);
  if((ibin < 0) || ((size_t)(ibin) >= bins.size()))
  {
      
      error->one(FLERR,"assertion failed");
  }

  bins[ibin].p_array.push_back(Particle(index,x, radius));
  ++ncount;
}

#ifdef SUPERQUADRIC_ACTIVE_FLAG
//the same for superquadrics
void RegionNeighborList::insert_superquadric(double * x, double radius, double *quaternion, double *shape) {
  int ibin = coord2bin(x);
  if((ibin < 0) || ((size_t)(ibin) >= bins.size()))
  {

      error->one(FLERR,"assertion failed");
  }

  bins[ibin].push_back(Particle(x, radius, quaternion, shape));
  ++ncount;
}
#endif

/**
 * @brief Clears neighbor list and brings it into initial state
 */
void RegionNeighborList::clear() {
  bins.clear();
  stencil.clear();
  ncount = 0;
}

/**
 * @brief Returns the number of particles inserted into the neighbor list
 * @return number of particles in neighbor list
 */
size_t RegionNeighborList::count() const {
  return ncount;
}

/**
 * @brief Update the region bounding box
 *
 * This will update internal data structures to ensure they can handle the new
 * region, which is defined by its bounding box.
 *
 * @param bb        bounding box of region
 * @param maxrad    largest particle radius
 * @return true if bounding box was set successfully, false bounding box could not
 * be set and neighbor list is not usable
 */
bool RegionNeighborList::setBoundingBox(BoundingBox & bb, double maxrad, bool extend, bool failsafe) {
  double extent[3];
  bb.getExtent(extent);

  if(extent[0] <= 0.0 || extent[1] <= 0.0 || extent[2] <= 0.0) {
    // empty or invalid region
    bins.clear();
    stencil.clear();
    return false;
  }

  bb.getBoxBounds(bboxlo, bboxhi);

  // testing code
  double binsize_optimal = 4*maxrad;
  double binsizeinv = 1.0/binsize_optimal;

  // test for too many global bins in any dimension due to huge global domain or small maxrad
  const int max_small_int = std::numeric_limits<int>::max();

  if (extent[0]*binsizeinv > max_small_int || extent[1]*binsizeinv > max_small_int ||
      extent[2]*binsizeinv > max_small_int)
  {

    if(failsafe)
    {
        binsizeinv = 1./ (vectorMax3D(extent) / 100.);
    }
    else
    {
        printf("ERROR: too many bins for this domain\n");
        return false;
    }
  }

  // create actual bins
  nbinx = static_cast<int>(extent[0]*binsizeinv);
  nbiny = static_cast<int>(extent[1]*binsizeinv);
  nbinz = static_cast<int>(extent[2]*binsizeinv);

  if (nbinx == 0) nbinx = 1;
  if (nbiny == 0) nbiny = 1;
  if (nbinz == 0) nbinz = 1;

  binsizex = extent[0]/nbinx;
  binsizey = extent[1]/nbiny;
  binsizez = extent[2]/nbinz;

  bininvx = 1.0 / binsizex;
  bininvy = 1.0 / binsizey;
  bininvz = 1.0 / binsizez;

  // mbinlo/hi = lowest and highest global bins my ghost atoms could be in
  // coord = lowest and highest values of coords for my ghost atoms
  // static_cast(-1.5) = -1, so subract additional -1
  // add in SMALL for round-off safety
  double bsubboxlo[3], bsubboxhi[3];
  bb.getBoxBounds(bsubboxlo, bsubboxhi);

  // list is extended for ghost atoms or
  // the list covers just exactly the region
  if (extend)
  {

      double coord = bsubboxlo[0] - SMALL*extent[0];
      mbinxlo = static_cast<int> ((coord-bboxlo[0])*bininvx);
      if (coord < bboxlo[0]) mbinxlo = mbinxlo - 1;
      coord = bsubboxhi[0] + SMALL*extent[0];
      int mbinxhi = static_cast<int> ((coord-bboxlo[0])*bininvx);

      coord = bsubboxlo[1] - SMALL*extent[1];
      mbinylo = static_cast<int> ((coord-bboxlo[1])*bininvy);
      if (coord < bboxlo[1]) mbinylo = mbinylo - 1;
      coord = bsubboxhi[1] + SMALL*extent[1];
      int mbinyhi = static_cast<int> ((coord-bboxlo[1])*bininvy);

      coord = bsubboxlo[2] - SMALL*extent[2];
      mbinzlo = static_cast<int> ((coord-bboxlo[2])*bininvz);
      if (coord < bboxlo[2]) mbinzlo = mbinzlo - 1;
      coord = bsubboxhi[2] + SMALL*extent[2];
      int mbinzhi = static_cast<int> ((coord-bboxlo[2])*bininvz);

      // extend bins by 1 to insure stencil extent is included

      mbinxlo = mbinxlo - 1;
      mbinxhi = mbinxhi + 1;
      mbinylo = mbinylo - 1;
      mbinyhi = mbinyhi + 1;
      mbinzlo = mbinzlo - 1;
      mbinzhi = mbinzhi + 1;

      mbinx = mbinxhi - mbinxlo + 1;
      mbiny = mbinyhi - mbinylo + 1;
      mbinz = mbinzhi - mbinzlo + 1;

#ifdef LIGGGHTS_DEBUG
  printf("setting insertion bounding box: [%g, %g] x [%g, %g] x [%g, %g]\n", bsubboxlo[0], bsubboxhi[0], bsubboxlo[1], bsubboxhi[1], bsubboxlo[2], bsubboxhi[2]);
  printf("nbinx %d nbiny %d, nbinz %d\n",nbinx,nbiny,nbinz);
  printf("mbinxlo: %d, mbinxhi: %d\n", mbinxlo, mbinxhi);
  printf("mbinylo: %d, mbinyhi: %d\n", mbinylo, mbinyhi);
  printf("mbinzlo: %d, mbinzhi: %d\n", mbinzlo, mbinzhi);
  printf("mbinx %d mbiny %d, mbinz %d\n",mbinx,mbiny,mbinz);
#endif

  }
  else
  {
      mbinxlo = mbinylo = mbinzlo = 0;
      mbinx = nbinx;
      mbiny = nbiny;
      mbinz = nbinz;

#ifdef LIGGGHTS_DEBUG
  printf("setting reduced! insertion bounding box: [%g, %g] x [%g, %g] x [%g, %g]\n", bsubboxlo[0], bsubboxhi[0], bsubboxlo[1], bsubboxhi[1], bsubboxlo[2], bsubboxhi[2]);
  printf("nbinx %d nbiny %d, nbinz %d\n",nbinx,nbiny,nbinz);
  printf("mbinxlo: %d\n", mbinxlo);
  printf("mbinylo: %d\n", mbinylo);
  printf("mbinzlo: %d\n", mbinzlo);
  printf("mbinx %d mbiny %d, mbinz %d\n",mbinx,mbiny,mbinz);
#endif

  }

  // allocate bins
  bigint bbin = ((bigint) mbinx) * ((bigint) mbiny) * ((bigint) mbinz);
  if (bbin > max_small_int) {
    printf("ERROR: Too many neighbor bins\n");
    return false;
  }
  bins.resize(bbin);
  // set cell center
  for(std::vector<Bin>::iterator it = bins.begin(); it != bins.end(); ++it)
  {
    const int ibin = it - bins.begin();
    int itmp = ibin;
    const int iz = ibin/(mbinx*mbiny);
    itmp -= iz*mbinx*mbiny;
    const int iy = itmp/mbinx;
    itmp -= iy*mbinx;
    const int ix = itmp;
    it->center[0] = bsubboxlo[0] + binsizex*((double)(ix+mbinxlo)+0.5);
    it->center[1] = bsubboxlo[1] + binsizey*((double)(iy+mbinylo)+0.5);
    it->center[2] = bsubboxlo[2] + binsizez*((double)(iz+mbinzlo)+0.5);
//#ifdef LIGGGHTS_DEBUG
//    printf("Bin (%d) with indizes: %d %d %d\n",ibin,ix,iy,iz);
//    printf("Center of bin (%d): %g %g %g\n",ibin,it->center[0],it->center[1],it->center[2]);
//#endif
  }

  // generate stencil which will look at all bins 27 bins
  for (int k = -1; k <= 1; k++)
    for (int j = -1; j <= 1; j++)
      for (int i = -1; i <= 1; i++)
        stencil.push_back(k*mbiny*mbinx + j*mbinx + i);

  bbox_set = true;

  return true;
}

/**
 * @brief Compute the closest distance between central bin (0,0,0) and bin (i,j,k)
 * @param i bin coordinate along x-axis
 * @param j bin coordinate along y-axis
 * @param k bin coordinate along z-axis
 * @return closest distance between central bin (0,0,0) and bin (i,j,k)
 */
double RegionNeighborList::bin_distance(int i, int j, int k)
{
  double delx,dely,delz;

  if (i > 0) delx = (i-1)*binsizex;
  else if (i == 0) delx = 0.0;
  else delx = (i+1)*binsizex;

  if (j > 0) dely = (j-1)*binsizey;
  else if (j == 0) dely = 0.0;
  else dely = (j+1)*binsizey;

  if (k > 0) delz = (k-1)*binsizez;
  else if (k == 0) delz = 0.0;
  else delz = (k+1)*binsizez;

  return (delx*delx + dely*dely + delz*delz);
}

/**
 * @brief Calc local bin index (m) of point x
 * @param x point in 3D
 * @return bin index of the given point x
 */
int RegionNeighborList::coord2binLocal(double *x) const
{
  int ix,iy,iz;

  if (x[0] >= bboxhi[0])
    ix = static_cast<int> ((x[0]-bboxhi[0])*bininvx) + nbinx;
  else if (x[0] >= bboxlo[0]) {
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx);
    ix = std::min(ix,nbinx-1);
  } else
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx) - 1;

  if (x[1] >= bboxhi[1])
    iy = static_cast<int> ((x[1]-bboxhi[1])*bininvy) + nbiny;
  else if (x[1] >= bboxlo[1]) {
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy);
    iy = std::min(iy,nbiny-1);
  } else
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy) - 1;

  if (x[2] >= bboxhi[2])
    iz = static_cast<int> ((x[2]-bboxhi[2])*bininvz) + nbinz;
  else if (x[2] >= bboxlo[2]) {
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz);
    iz = std::min(iz,nbinz-1);
  } else
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz) - 1;

  return (iz-mbinzlo)*mbiny*mbinx + (iy-mbinylo)*mbinx + (ix-mbinxlo);
}

/**
 * @brief Calc global bin index (n) of point x
 * @param x point in 3D
 * @return bin index of the given point x
 */
int RegionNeighborList::coord2binGlobal(double *x) const
{
    
    return RegionNeighborList::coord2binLocal(x);
    /*
    int ix,iy,iz;

    if (x[0] >= bboxhi[0])
      ix = static_cast<int> ((x[0]-bboxhi[0])*bininvx) + nbinx;
    else if (x[0] >= bboxlo[0]) {
      ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx);
      ix = std::min(ix,nbinx-1);
    } else
      ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx) - 1;

    if (x[1] >= bboxhi[1])
      iy = static_cast<int> ((x[1]-bboxhi[1])*bininvy) + nbiny;
    else if (x[1] >= bboxlo[1]) {
      iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy);
      iy = std::min(iy,nbiny-1);
    } else
      iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy) - 1;

    if (x[2] >= bboxhi[2])
      iz = static_cast<int> ((x[2]-bboxhi[2])*bininvz) + nbinz;
    else if (x[2] >= bboxlo[2]) {
      iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz);
      iz = std::min(iz,nbinz-1);
    } else
      iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz) - 1;

    return (iz)*nbiny*nbinx + (iy)*nbinx + (ix);*/
}

/**
 * @brief Calc global bin index (n) of point x
 * @param x point in 3D
 * @return bin index of the given point x
 */

bool RegionNeighborList::isInBoundingBox(double *pos)
{
    if(!boundingBoxSet())
        error->one(FLERR,"internal error: call before bbox is set");
    if
    (
        pos[0] >= bboxlo[0] && pos[0] <= bboxhi[0] &&
        pos[1] >= bboxlo[1] && pos[1] <= bboxhi[1] &&
        pos[2] >= bboxlo[2] && pos[2] <= bboxhi[2]
    )   return true;
    return false;
}
