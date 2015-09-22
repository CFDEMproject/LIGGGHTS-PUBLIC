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

    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Philippe Seil (JKU Linz)
    Richard Berger (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2015 JKU Linz
------------------------------------------------------------------------- */

#include "bounding_box.h"
#include <algorithm>

namespace LAMMPS_NS
{

  BoundingBox::BoundingBox()
  : xLo(0.), xHi(0.), yLo(0.), yHi(0.), zLo(0.), zHi(0.), initGiven(false), dirty(true)
  {}
  BoundingBox::BoundingBox(double xLo_, double xHi_, double yLo_, double yHi_, double zLo_, double zHi_)
  : xLo(xLo_), xHi(xHi_), yLo(yLo_), yHi(yHi_), zLo(zLo_), zHi(zHi_), initGiven(true), dirty(true)
  {}

  BoundingBox::~BoundingBox()
  {}

  void BoundingBox::reset()
  {
    xLo = 0.; xHi = 0.;
    yLo = 0.; yHi = 0.;
    zLo = 0.; zHi = 0.;
    initGiven = false;
    dirty = true;
  }

  void BoundingBox::extendByDelta(double delta)
  {
    xLo = xLo-delta;
    yLo = yLo-delta;
    zLo = zLo-delta;
    xHi = xHi+delta;
    yHi = yHi+delta;
    zHi = zHi+delta;
  }

  void BoundingBox::extrude(double length, const double * vec)
  {
    xLo = std::min(xLo, (xLo + length * vec[0]));
    yLo = std::min(yLo, (yLo + length * vec[1]));
    zLo = std::min(zLo, (zLo + length * vec[2]));
    xHi = std::max(xHi, (xHi + length * vec[0]));
    yHi = std::max(yHi, (yHi + length * vec[1]));
    zHi = std::max(zHi, (zHi + length * vec[2]));
  }

} /* namespace LAMMPS_NS */
