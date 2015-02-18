/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
   the producer of the LIGGGHTS® software and the CFDEM®coupling software
   See http://www.cfdem.com/terms-trademark-policy for details.

   LIGGGHTS® is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/*
Contributing authors:
Stefan Amberger (JKU Linz)
Christoph Kloss (DCS Computing GmbH, Linz and JKU Linz)
*/

#ifdef REGION_CLASS

RegionStyle(wedge,RegWedge)

#else

#ifndef domain_region_wedge_h
#define domain_region_wedge_h

#include "region.h"

// command: region id wedge dim c1 c2 radius lo hi angle1 angle2
//   dim ... x or y or z
//   c1 ... if dim == x then y-coord
//          if dim == y then z-coord
//          if dim == z then x-coord
//   c2 ... if dim == x then z-coord
//          if dim == y then x-coord
//          if dim == z then y-coord
//   radius ... radius of cylinder
//   lo ... dim-coord of lower flat face of cylinder of wedge
//   hi ... dim-coord of higher flat face of cylinder of wedge
//   angle1 ... mathematically positive angle of the wedge's starting face
//              if axis == x then starting from y-axis
//              if axis == y then starting from z-axis
//              if axis == z then starting from x-axis
//   dang ... angle between the wedge's ending face

namespace LAMMPS_NS
{
class RegWedge : public Region
{
  friend class DomainWedge;

  public:
    RegWedge(class LAMMPS *, int, char **);
    ~RegWedge();

    int inside(double,double,double);
    int surface_exterior(double *, double);
    int surface_interior(double *, double);

  private:

    // internal data
    char axis;
    double c1,c2;
    double radius;
    double lo,hi;
    double angle1, cosang1, sinang1;
    double dang, cosdang, sindang, cosmdang, sinmdang;
    double angle2, cosang2, sinang2;

    // helper data
    double pi_half;
    double onedivr;

    // normal vectors
    // the normal vectors point outside the wedge and are normalized
    double normal1[2]; // normal vector to 1st plane of section (angle1)
    double normal2[2]; // normal vector to 2nd plane of section (angle2)

    // helper functions
    inline void snormalize2(const double length, const double *v, double *ans);
    void printRegion();
    void printProperty(const char *name, double val);
    void printContacts(double *, int);
};

}

#endif
#endif
