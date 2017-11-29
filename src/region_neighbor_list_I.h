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
    Josef Kerbl (DCS Computing GmbH)

    Copyright 2014-2015 JKU Linz
    Copyright 2015-     DCS Computing GmbH
------------------------------------------------------------------------- */

// include last to ensure correct macros
#include "domain_definitions.h"

static const double SMALL_REGION_NEIGHBOR_LIST = 1.0e-6;
static const double BIG_REGION_NEIGHBOR_LIST = 1.0e20;

/**
 * @brief Default constructor which will create an empty neighbor list
 */
template<bool INTERPOLATE>
RegionNeighborList<INTERPOLATE>::RegionNeighborList(LAMMPS *lmp) :
    IRegionNeighborList(),
    Pointers(lmp),
    ncount(0),
    bbox_set(false)
{
}

/**
 * @brief Determine if the given particle overlaps with any particle in this neighbor list
 * @param x        position of particle to check
 * @param radius   radius of particle to check
 * @return true if particle has an overlap with a particle in this neighbor list, false otherwise
 */
template<bool INTERPOLATE>
bool RegionNeighborList<INTERPOLATE>::hasOverlap(double * x, double radius) const
{
  int ibin = coord2bin(x);

  for(std::vector<int>::const_iterator it = stencil.begin(); it != stencil.end(); ++it)
  {
    const int offset = *it;
    if((ibin+offset < 0) || ((size_t)(ibin+offset) >= bins.size()))
    {
        
        error->one(FLERR,"assertion failed");
    }
    const std::vector<Particle<INTERPOLATE> > & plist = bins[ibin+offset].particles;

    for(typename std::vector<Particle<INTERPOLATE> >::const_iterator pit = plist.begin(); pit != plist.end(); ++pit)
    {
      const Particle<INTERPOLATE> & p = *pit;
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
template<bool INTERPOLATE>
bool RegionNeighborList<INTERPOLATE>::hasOverlap_superquadric(double * x, double radius, double *quaternion, double *shape, double *blockiness) const
{
  int ibin = coord2bin(x);

  for(std::vector<int>::const_iterator it = stencil.begin(); it != stencil.end(); ++it) {
    const int offset = *it;
    if((ibin+offset < 0) || ((size_t)(ibin+offset) >= bins.size()))
    {

        error->one(FLERR,"assertion failed");
    }
    const std::vector<Particle<INTERPOLATE> > & plist = bins[ibin+offset].particles;

    Superquadric particle1(x, quaternion, shape, blockiness);
    for(typename std::vector<Particle<INTERPOLATE> >::const_iterator pit = plist.begin(); pit != plist.end(); ++pit) {
      const Particle<INTERPOLATE> & p = *pit;
      double del[3];
      vectorSubtract3D(x, p.x, del);
      const double rsq = vectorMag3DSquared(del);
      const double radsum = radius + p.radius;
      if(check_obb_flag) {
        double x_copy[3], quaternion_copy[4], shape_copy[3], blockiness_copy[2];
        vectorCopy3D(p.x, x_copy);
        vectorCopy4D(p.quaternion, quaternion_copy);
        vectorCopy3D(p.shape, shape_copy);
        vectorCopy2D(p.blockiness, blockiness_copy);
        Superquadric particle2(x_copy, quaternion_copy, shape_copy, blockiness_copy);
        double contact_point[3];

        if (rsq <= radsum*radsum and (MathExtraLiggghtsNonspherical::capsules_intersect(&particle1, &particle2, contact_point) and
                                      MathExtraLiggghtsNonspherical::obb_intersect(&particle1, &particle2))) {
            return true;
        }
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

template<bool INTERPOLATE>
bool RegionNeighborList<INTERPOLATE>::hasOverlapWith(double * x, double radius, std::vector<int> &overlap_list) const
{
  int ibin = coord2bin(x);

  bool overlap = false;

  for(std::vector<int>::const_iterator it = stencil.begin(); it != stencil.end(); ++it)
  {
    const int offset = *it;
    if((ibin+offset < 0) || ((size_t)(ibin+offset) >= bins.size()))
    {
        
        error->one(FLERR,"assertion failed");
    }
    const std::vector<Particle<INTERPOLATE> > & plist = bins[ibin+offset].particles;

    for(typename std::vector<Particle<INTERPOLATE> >::const_iterator pit = plist.begin(); pit != plist.end(); ++pit)
    {
      const Particle<INTERPOLATE> & p = *pit;
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
template<bool INTERPOLATE>
void RegionNeighborList<INTERPOLATE>::insert(double * x, double radius,int index)
{
  int quadrant;
  double wx,wy,wz;
  int ibin = coord2bin(x,quadrant,wx,wy,wz);
  if((ibin < 0) || ((size_t)(ibin) >= bins.size()))
  {
      
      error->one(FLERR,"assertion failed");
  }

  bins[ibin].particles.push_back(Particle<INTERPOLATE>(index,x, radius,ibin,quadrant,wx,wy,wz));

  ++ncount;
}

#ifdef SUPERQUADRIC_ACTIVE_FLAG
template<bool INTERPOLATE>
void RegionNeighborList<INTERPOLATE>::insert_superquadric(double * x, double radius, double *quaternion, double *shape, double *blockiness, int index) {
  int quadrant;
  double wx,wy,wz;
  int ibin = coord2bin(x,quadrant,wx,wy,wz);
  if((ibin < 0) || ((size_t)(ibin) >= bins.size()))
  {

      error->one(FLERR,"assertion failed");
  }

  bins[ibin].particles.push_back(Particle<INTERPOLATE>(index,x,radius,quaternion, shape, blockiness, ibin,quadrant,wx,wy,wz));
  ++ncount;
}
#endif

/**
 * @brief Reset neighbor list and brings it into initial state
 */
template<bool INTERPOLATE>
void RegionNeighborList<INTERPOLATE>::reset()
{
  // reset all the settings from setBoundaryBox
  bins.clear();
  stencil.clear();
  bbox_set = false;
  nbinx = nbiny = nbinz = 0;
  binsizex = binsizey = binsizez = 0;
  bininvx = bininvy = bininvz = 0;
  mbinxlo = mbinylo = mbinzlo = 0;
  mbinx = mbiny = mbinz = 0;

  // reset particle counter
  ncount = 0;
}

/**
 * @brief Return size of buffer for one bin
 *
 * \sa pushBinToBuffer
 */
template<bool INTERPOLATE>
int RegionNeighborList<INTERPOLATE>::getSizeOne() const
{
    return 4;
}

/**
 * @brief Push bin information to buffer
 *
 * \sa getSizeOne
 */
template<bool INTERPOLATE>
int RegionNeighborList<INTERPOLATE>::pushBinToBuffer(int i, double *buf) const
{
    int m = 0;
    buf[m++] = this->bins[i].center[0];
    buf[m++] = this->bins[i].center[1];
    buf[m++] = this->bins[i].center[2];
    buf[m++] = this->bins[i].id;
    return m;
}

/**
 * @brief Clear bins (remove all inserted particles)
 */
template<bool INTERPOLATE>
void RegionNeighborList<INTERPOLATE>::clear()
{
  for(typename std::vector<Bin<INTERPOLATE> >::iterator it = bins.begin(); it != bins.end(); ++it)
  {
    (*it).particles.clear();
  }
  ncount = 0;
}

/**
 * @brief Returns the number of particles inserted into the neighbor list
 * @return number of particles in neighbor list
 */

template<bool INTERPOLATE>
size_t RegionNeighborList<INTERPOLATE>::count() const
{
  return ncount;
}

/**
 * @brief Set or Update the region bounding box
 *
 * This will update internal data structures to ensure they can handle the new
 * region, which is defined by its bounding box.
 *
 * @param bb        bounding box of region
 * @param maxrad    largest particle radius
 * @param extend    enable extending the box for ghost particles
 * @param failsafe  limit max. number of bins
 * @return true if bounding box was set successfully, false bounding box could not
 * be set and neighbor list is not usable
 */

template<bool INTERPOLATE>
inline void RegionNeighborList<INTERPOLATE>::setBoundingBox_calc_interpolation_stencil(Bin<INTERPOLATE> &it,int ibin,int ix,int iy, int iz) const
{

}

template<>
inline void RegionNeighborList<true>::setBoundingBox_calc_interpolation_stencil(Bin<true> &it,int ibin,int ix,int iy, int iz) const
{
        Stencil stencilTmp;

        it.stencils.clear();

        int stencil_shift_up_quadrant = 0;
        int stencil_shift_down_quadrant = 0;
        int ii,jj,kk;
        for(int iquad = 0; iquad < 8; iquad++)
        {
            
            stencil_shift_up_quadrant = 0;
            stencil_shift_down_quadrant = 0;
            const int qx = (iquad & 1) >> 0; // 0 or 1
            const int qy = (iquad & 2) >> 1; // 0 or 1
            const int qz = (iquad & 4) >> 2; // 0 or 1

            if(0 == ix && 0 == qx)
            {
                stencil_shift_up_quadrant |= 1;
                ii = 1;
            }
            else if(ix == (mbinx-1)  && 1 == qx)
            {
                stencil_shift_down_quadrant |= 1;
                ii = -1;
            }
            else
                ii = 0;

            if(0 == iy && 0 == qy)
            {
                stencil_shift_up_quadrant |= 2;
                jj = 1;
            }
            else if(iy == (mbiny-1) && 1 == qy)
            {
                stencil_shift_down_quadrant |= 2;
                jj = -1;
            }
            else
                jj = 0;

            if(0 == iz && 0 == qz)
            {
                stencil_shift_up_quadrant |= 4;
                kk = 1;
            }
            else if(iz == (mbinz-1) && 1 == qz)
            {
                stencil_shift_down_quadrant |= 4;
                kk = -1;
            }
            else
                kk = 0;

            // generate stencil which will look at 8 relevant bins - self + 7 quadrant neighs
            // self can be located anywhere within these 8 bins
            stencilTmp.clear();
            for (int k = -1+qz; k <= qz; k++)
            {
                for (int j = -1+qy; j <= qy; j++)
                {
                    for (int i = -1+qx; i <= qx; i++)
                    {
                        int isubbin = ((iz+k+kk)*mbiny*mbinx + (iy+j+jj)*mbinx + (ix+i+ii));
                        if( isubbin < 0 || isubbin >= mbins())
                            stencilTmp.push_back(-1);
                        else
                            stencilTmp.push_back(isubbin);
                        
                    }
                }
            }

            it.stencils.push_back(stencilTmp);
            
            it.stencil_shift_up.push_back(stencil_shift_up_quadrant);
            it.stencil_shift_down.push_back(stencil_shift_down_quadrant);

        }
}

template<bool INTERPOLATE>
bool RegionNeighborList<INTERPOLATE>::setBoundingBox(BoundingBox & bb, double maxrad, bool extend, bool failsafe)
{
  double extent[3];
  bb.getExtent(extent);

  if(extent[0] <= 0.0 || extent[1] <= 0.0 || extent[2] <= 0.0) {
    // empty or invalid region
    bins.clear();
    stencil.clear();
    return false;
  }

  bb.getBoxBounds(bboxlo, bboxhi);

  // direct comparison of doubles to detect INF
  if ( bboxlo[0] == -BIG || bboxlo[1] == -BIG || bboxlo[2] == -BIG ||
       bboxhi[0] == BIG || bboxhi[1] == BIG || bboxhi[2] == BIG )
  {
      error->one(FLERR,"'INF' not allowed for definiton of region used for RegionNeighborList.\n"
                       "You may want to use 'EDGE' instead.");
  }

  // testing code
  double binsize_optimal = 4*maxrad;
  double binsizeinv = 1.0/binsize_optimal;

  // test for too many global bins in any dimension due to huge global domain or small maxrad
  const int max_bins = 8000000;
  // we may repeat this calculation in case of failsafe
  bool repeat = false;
  bigint bbin = -1; // final number of bins
  double bsubboxlo[3], bsubboxhi[3]; // box bounds
  do
  {
    // create actual bins
    nbinx = static_cast<int>((extent[0]+binsize_optimal*1e-10)*binsizeinv);
    nbiny = static_cast<int>((extent[1]+binsize_optimal*1e-10)*binsizeinv);
    nbinz = static_cast<int>((extent[2]+binsize_optimal*1e-10)*binsizeinv);

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
    bb.getBoxBounds(bsubboxlo, bsubboxhi);

    // list is extended for ghost atoms or
    // the list covers just exactly the region
    if (extend)
    {
      double coord = bsubboxlo[0] - SMALL_REGION_NEIGHBOR_LIST*extent[0];
      mbinxlo = static_cast<int> ((coord-bboxlo[0])*bininvx);
      if (coord < bboxlo[0]) mbinxlo = mbinxlo - 1;
      coord = bsubboxhi[0] + SMALL_REGION_NEIGHBOR_LIST*extent[0];
      int mbinxhi = static_cast<int> ((coord-bboxlo[0])*bininvx);

      coord = bsubboxlo[1] - SMALL_REGION_NEIGHBOR_LIST*extent[1];
      mbinylo = static_cast<int> ((coord-bboxlo[1])*bininvy);
      if (coord < bboxlo[1]) mbinylo = mbinylo - 1;
      coord = bsubboxhi[1] + SMALL_REGION_NEIGHBOR_LIST*extent[1];
      int mbinyhi = static_cast<int> ((coord-bboxlo[1])*bininvy);

      coord = bsubboxlo[2] - SMALL_REGION_NEIGHBOR_LIST*extent[2];
      mbinzlo = static_cast<int> ((coord-bboxlo[2])*bininvz);
      if (coord < bboxlo[2]) mbinzlo = mbinzlo - 1;
      coord = bsubboxhi[2] + SMALL_REGION_NEIGHBOR_LIST*extent[2];
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

    }
    else
    {
      mbinxlo = mbinylo = mbinzlo = 0;
      mbinx = nbinx;
      mbiny = nbiny;
      mbinz = nbinz;

    }

    // final check of number of bins
    bbin = ((bigint) mbinx) * ((bigint) mbiny) * ((bigint) mbinz);

    if (bbin > max_bins) { // too many bins
      // check for failsafe mode
      // check for repeat - ensure that we run the loop only once!
      if(failsafe && !repeat)
      {
        binsizeinv = 1./ (vectorMax3D(extent) / 100.);
        repeat = true;
      } else {
        printf("ERROR: Too many neighbor bins\n");
        return false;
      }
    }
    else
        repeat = false;
  } while (repeat);

  // allocate bins
  bins.resize(bbin);

  // set cell center and stencils for each bin
  int cId = 0;
  for(typename std::vector<Bin<INTERPOLATE> >::iterator it = bins.begin(); it != bins.end(); ++it)
  {
    const int ibin = it - bins.begin();
    int itmp = ibin;
    const int iz = ibin/(mbinx*mbiny);
    itmp -= iz*mbinx*mbiny;
    const int iy = itmp/mbinx;
    itmp -= iy*mbinx;
    const int ix = itmp;
    it->id = cId++;
    it->center[0] = bsubboxlo[0] + binsizex*((double)(ix+mbinxlo)+0.5);
    it->center[1] = bsubboxlo[1] + binsizey*((double)(iy+mbinylo)+0.5);
    it->center[2] = bsubboxlo[2] + binsizez*((double)(iz+mbinzlo)+0.5);

    setBoundingBox_calc_interpolation_stencil(*it,ibin,ix,iy,iz);

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
template<bool INTERPOLATE>
double RegionNeighborList<INTERPOLATE>::bin_distance(int i, int j, int k)
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

template<bool INTERPOLATE>
inline void RegionNeighborList<INTERPOLATE>::coord2bin_calc_interpolation_weights(double *x,int ibin,int ix,int iy, int iz,int &quadrant,double &wx,double &wy,double &wz) const
{
    quadrant = 0;
    wx = wy = wz = 0.;
}

template<>
inline void RegionNeighborList<true>::coord2bin_calc_interpolation_weights(double *x,int ibin,int ix,int iy, int iz,int &quadrant,double &wx,double &wy,double &wz) const
{
      double dx = (((x[0]-bboxlo[0])*bininvx) - static_cast<int> ((x[0]-bboxlo[0])*bininvx));//*bininvx;
      double dy = (((x[1]-bboxlo[1])*bininvy) - static_cast<int> ((x[1]-bboxlo[1])*bininvy));//*bininvy;
      double dz = (((x[2]-bboxlo[2])*bininvz) - static_cast<int> ((x[2]-bboxlo[2])*bininvz));//*bininvz;
      quadrant = (dx >= 0.5) * 1 + (dy >= 0.5) * 2 + (dz >= 0.5) * 4;

      if(dx >= 0.5)
        wx = dx-0.5;
      else
        wx = 0.5+dx;

      if(dy >= 0.5)
        wy = dy-0.5;
      else
        wy = 0.5+dy;

      if(dz >= 0.5)
        wz = dz-0.5;
      else
        wz = 0.5+dz;

      int stencil_shift_up = bins[ibin].stencil_shift_up[quadrant];
      int stencil_shift_down = bins[ibin].stencil_shift_down[quadrant];
      if( 0 == stencil_shift_up && 0 == stencil_shift_down)
        return;

      if(stencil_shift_down & 1)
        wx = 1.;
      else if(stencil_shift_up & 1)
        wx = 0.;

      if(stencil_shift_down & 2)
        wy = 1.;
      else if(stencil_shift_up & 2)
        wy = 0.;

      if(stencil_shift_down & 4)
        wz = 1.;
      else if(stencil_shift_up & 4)
        wz = 0.;

}

template<bool INTERPOLATE>
int RegionNeighborList<INTERPOLATE>::coord2bin(double *x,int &quadrant,double &wx,double &wy,double &wz) const
{
  int ix,iy,iz;

  /*if(x[0] < bboxlo[0] || x[0] > bboxhi[0] || x[1] < bboxlo[1] || x[1] > bboxhi[1] || x[2] < bboxlo[2] || x[2] > bboxhi[2])
    return -1;*/

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

  int ibin = (iz-mbinzlo)*mbiny*mbinx + (iy-mbinylo)*mbinx + (ix-mbinxlo);

  if(ibin < 0 || ibin >= mbins())
     return -1; // ibin = -1
  coord2bin_calc_interpolation_weights(x,ibin,ix,iy,iz,quadrant,wx,wy,wz);

  return ibin;
}

/**
 * @brief Calc global bin index (n) of point x
 * @param x point in 3D
 * @return bin index of the given point x
 */

template<bool INTERPOLATE>
bool RegionNeighborList<INTERPOLATE>::isInBoundingBox(double *pos) const
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

/**
 * @brief Set bounding box from region
 * @param region   Region used for creating the neighbor list
 * @param maxrad    largest particle radius
 * @param extend    enable extending the box for ghost particles
 * @param failsafe  limit max. number of bins
 * @return BoundingBox the bounding box of the region
 */

template<bool INTERPOLATE>
BoundingBox RegionNeighborList<INTERPOLATE>::setBoundingBoxRegion(const Region &region, double maxrad, bool extend, bool failsafe)
{
    // use region to determine and set bounding box
    BoundingBox bb(region.extent_xlo, region.extent_xhi, region.extent_ylo, region.extent_yhi, region.extent_zlo, region.extent_zhi);
    if (!setBoundingBox(bb, maxrad, extend, failsafe))
        error->one(FLERR,"internal error: could not set bounding box");

    return bb;
}
