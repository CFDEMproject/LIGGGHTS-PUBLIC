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
   Stefan Amberger (JKU Linz)
------------------------------------------------------------------------- */

#ifndef LMP_MODIFIED_ANDREW_H
#define LMP_MODIFIED_ANDREW_H

#include <vector>
#include "pointers.h"

using namespace std;

namespace MODIFIED_ANDREW_AUX{

struct Point {
  double x, y;
};

struct Circle {
  double x, y, r;

  bool operator <(const Circle &p) const {
    return x-r < p.x-p.r || (x-r == p.x-p.r && y-r < p.y-p.r);
  }
};

struct Line {
  double a,b,c;
};

}

using MODIFIED_ANDREW_AUX::Point;
using MODIFIED_ANDREW_AUX::Circle;
using MODIFIED_ANDREW_AUX::Line;

namespace LAMMPS_NS{

class ModifiedAndrew : protected Pointers {

public:
  ModifiedAndrew(LAMMPS *lmp);
  ~ModifiedAndrew();

  double area();

  inline void add_contact(Circle c)
  { contacts_.push_back(c); }

private:

  double area(vector<Point> H);
  double area(Point &p, Point &m, Point &q);
  double turn(const Circle &O, const Circle &A, const Circle &B);

  Line find_right_tangent(const Circle &p1, const Circle &p2);

  Point intersectLines(Line &l, Line &j);
  Point mean_point(vector<Point> P);

  vector<Circle> convex_hull(vector<Circle> P);
  vector<Point> hull_points(vector<Circle> C);
  vector<Circle> construct_hull_c_all(double *data0, int ndata0);
  int construct_data(vector<Circle> hull_c, double *&data);

  void project_contancts();

  // container for contacts: x, y, r
  vector<Circle> contacts_;
};

}

#endif
