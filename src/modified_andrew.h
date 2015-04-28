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
    Stefan Amberger (JKU Linz)
    Christoph Kloss (DCS Computing GmbH, Linz; JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef LMP_MODIFIED_ANDREW_H
#define LMP_MODIFIED_ANDREW_H

#include <vector>
#include "pointers.h"

using namespace std;

namespace MODIFIED_ANDREW_AUX{

struct Point {
  double x, y;

  bool operator <(const Point &p) const {
    return x < p.x || (x == p.x && y < p.y);
  }
};

struct Circle {
  double x, y, r;
};

}

using MODIFIED_ANDREW_AUX::Point;
using MODIFIED_ANDREW_AUX::Circle;

namespace LAMMPS_NS{

class ModifiedAndrew : protected Pointers {

public:
  ModifiedAndrew(LAMMPS *lmp);
  ~ModifiedAndrew();

  double area();

  void add_contact(Circle c);

  inline void clear_contacts()
  { contacts_.clear(); }

private:

  double area(vector<Point> H);
  double area(Point &p, Point &m, Point &q);
  double cross(Point O, Point A, Point B);

  Point mean_point(vector<Point> P);

  vector<Point> convex_hull(vector<Point> P);
  vector<Point> construct_hull_c_all(double *data0, int ndata0);
  int construct_data(vector<Point> hull_c, double *&data);

  // container for contacts: x, y
  vector<Point> contacts_;

  int npoints_per_circle_;
  double **directions_;
};

}

#endif
