/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Philippe Seil (JKU Linz)
------------------------------------------------------------------------- */

#include "bounding_box.h"
#include "mpi_liggghts.h"

namespace LAMMPS_NS
{

  BoundingBox::BoundingBox()
  : xLo(0.), xHi(0.), yLo(0.), yHi(0.), zLo(0.), zHi(0.), initGiven(false)
  {}
  BoundingBox::BoundingBox(double xLo_, double xHi_, double yLo_, double yHi_, double zLo_, double zHi_)
  : xLo(xLo_), xHi(xHi_), yLo(yLo_), yHi(yHi_), zLo(zLo_), zHi(zHi_), initGiven(true)
  {}

  BoundingBox::~BoundingBox()
  {}

  void BoundingBox::extendToContain(double const *pt)
  {
    if(initGiven){
        if(pt[0] < xLo) xLo = pt[0];
        else if(pt[0] > xHi) xHi = pt[0];

        if(pt[1] < yLo) yLo = pt[1];
        else if(pt[1] > yHi) yHi = pt[1];

        if(pt[2] < zLo) zLo = pt[2];
        else if(pt[2] > zHi) zHi = pt[2];
    } else{
      xLo = pt[0]; xHi = pt[0];
      yLo = pt[1]; yHi = pt[1];
      zLo = pt[2]; zHi = pt[2];
      initGiven = true;
    }
  }

  void BoundingBox::extendToParallel(MPI_Comm comm)
  {
      double limit[6];
      limit[0] = -xLo;
      limit[1] =  xHi;
      limit[2] = -yLo;
      limit[3] =  yHi;
      limit[4] = -zLo;
      limit[5] =  zHi;

     MPI_Max_Vector(limit,6,comm);

      xLo = -limit[0];
      xHi =  limit[1];
      yLo = -limit[2];
      yHi =  limit[3];
      zLo = -limit[4];
      zHi =  limit[5];
  }

  void BoundingBox::reset()
  {
    xLo = 0.; xHi = 0.;
    yLo = 0.; yHi = 0.;
    zLo = 0.; zHi = 0.;
    initGiven = false;
  }

} /* namespace LAMMPS_NS */
