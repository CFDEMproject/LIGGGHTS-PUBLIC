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

#ifndef REGION_NEIGHBOR_LIST_DEFINITIONS_H
#define REGION_NEIGHBOR_LIST_DEFINITIONS_H

#include <vector>
#include "vector_liggghts.h"
#include "pointers.h"
#include "bounding_box.h"

namespace LAMMPS_NS {

class Region;

/**
 * @brief A small particle structure
 */

class ParticleBase
{
public:
    int index;
    double x[3];
    double radius;

#ifdef SUPERQUADRIC_ACTIVE_FLAG
  double shape[3];
  double quaternion[4];
  double blockiness[2];
#endif
};

template<bool INTERPOLATION>
class Particle : public ParticleBase
{};

template<>
class Particle<false /*interpolation*/> : public ParticleBase
{
 public:
  Particle(int _i,double * _pos, double _rad, int,int,double,double,double) {
    index = _i;
    LAMMPS_NS::vectorCopy3D(_pos, x);
    radius = _rad;

#ifdef SUPERQUADRIC_ACTIVE_FLAG
    quaternion[0] = 1.0;
    quaternion[1] = quaternion[2] = quaternion[3] = 0.0;
    shape[0] = shape[1] = shape[2] = radius;
    blockiness[0] = blockiness[1] = 2.0;
#endif
  }

#ifdef SUPERQUADRIC_ACTIVE_FLAG
  Particle(int _i,double * pos, double rad, double *quaternion_, double *shape_, double *blockiness_, int,int,double,double,double) {
      index = _i;
      LAMMPS_NS::vectorCopy3D(pos, x);
      radius = rad;
      LAMMPS_NS::vectorCopy4D(quaternion_, quaternion);
      LAMMPS_NS::vectorCopy3D(shape_, shape);
      LAMMPS_NS::vectorCopy2D(blockiness_, blockiness);
    }
#endif
};

template<>
class Particle<true /*interpolation*/> : public ParticleBase
{
 public:
  
  int ibin;
  int quadrant_bitfield;
  double wx, wy, wz;

  Particle(int _i,double * _pos, double _rad,int _ibin,int _quadrant,double _wx = -1.,double _wy = -1.,double _wz = -1.) {
    index = _i;
    LAMMPS_NS::vectorCopy3D(_pos, x);
    radius = _rad;
    ibin = _ibin;
    quadrant_bitfield = _quadrant;
    wx = _wx;
    wy = _wy;
    wz = _wz;
#ifdef SUPERQUADRIC_ACTIVE_FLAG
    quaternion[0] = 1.0;
    quaternion[1] = quaternion[2] = quaternion[3] = 0.0;
    shape[0] = shape[1] = shape[2] = radius;
    blockiness[0] = blockiness[1] = 2.0;
#endif
  }

#ifdef SUPERQUADRIC_ACTIVE_FLAG
  Particle(int _i, double * pos, double rad, double *quaternion_, double *shape_, double *blockiness_, int _ibin,int _quadrant,double _wx = -1.,double _wy = -1.,double _wz = -1.) {
      index = _i;
      LAMMPS_NS::vectorCopy3D(pos, x);
      radius = rad;
      LAMMPS_NS::vectorCopy4D(quaternion_, quaternion);
      LAMMPS_NS::vectorCopy3D(shape_, shape);
      LAMMPS_NS::vectorCopy2D(blockiness_, blockiness);

      ibin = _ibin;
      quadrant_bitfield = _quadrant;
      wx = _wx;
      wy = _wy;
      wz = _wz;
    }
#endif
};

/**
 * @brief typedefs
 */

typedef std::vector<Particle<true> > ParticleListInterpolate;
typedef std::vector<Particle<false> > ParticleListNoInterpolate;

typedef std::vector<int> Stencil;

const static bool interpolate_yes = true;
const static bool interpolate_no  = false;

/**
 * @brief Bin class
 *
 * Defines data structure for binning
 * Variants with and without interpolation.  Interpolation
 * variant stores lists with offsets
 */

class BinBase
{
public:
    int id;
    double center[3];
};

template<bool INTERPOLATION>
class Bin : public BinBase
{};

template<>
class Bin<false> : public BinBase
{
  public:
    ParticleListNoInterpolate particles;
};

template<>
class Bin<true> : public BinBase
{
  public:
    ParticleListInterpolate particles;

    std::vector<Stencil> stencils;
    
    std::vector<int> stencil_shift_down;

    std::vector<int> stencil_shift_up;
};

}

#endif
