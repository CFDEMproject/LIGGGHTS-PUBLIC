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

#ifndef LMP_BOUNDING_BOX
#define LMP_BOUNDING_BOX

#include "mpi.h"

namespace LAMMPS_NS
{

class BoundingBox
{
  friend class FixNeighlistMesh;

  public:

    BoundingBox();
    BoundingBox(double xLo, double xHi, double yLo, double yHi, double zLo, double zHi);
    virtual ~BoundingBox();

    void extendToContain(double const *pt);
    void extendToParallel(MPI_Comm comm);
    void reset();

    void getBoxBounds(double *lo,double *hi)
    {
        lo[0] = xLo;
        lo[1] = yLo;
        lo[2] = zLo;
        hi[0] = xHi;
        hi[1] = yHi;
        hi[2] = zHi;
    }

    void getBoxBoundsExtendedByDelta(double *lo,double *hi,double delta)
    {
        lo[0] = xLo-delta;
        lo[1] = yLo-delta;
        lo[2] = zLo-delta;
        hi[0] = xHi+delta;
        hi[1] = yHi+delta;
        hi[2] = zHi+delta;
    }

    void shrinkToSubbox(double *sublo,double *subhi)
    {
        if(xLo < sublo[0])
            xLo = sublo[0];
        if(xHi > subhi[0])
            xHi = subhi[0];

        if(yLo < sublo[1])
            yLo = sublo[1];
        if(yHi > subhi[1])
            yHi = subhi[1];

        if(zLo < sublo[2])
            zLo = sublo[2];
        if(xHi > subhi[2])
            zHi = subhi[2];
    }

    bool isInitialized()
    { return initGiven; }

    bool isInside(double *p)
    {
       // check bbox
       // test for >= and < as in Domain class
        if (p[0] >= xLo && p[0] < xHi &&
            p[1] >= yLo && p[1] < yHi &&
            p[2] >= zLo && p[2] < zHi)
            return true;
        return false;
    }

  private:

    double xLo, xHi, yLo, yHi, zLo, zHi;

    bool initGiven;
};

} /* LAMMPS_NS */
#endif /* BOUNDINGBOX_H_ */
