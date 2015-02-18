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

/* ----------------------------------------------------------------------
   Contributing authors:
   Stefan Amberger (JKU Linz)
------------------------------------------------------------------------- */

#include "modified_andrew.h"
#include "mpi_liggghts.h"
#include "vector_liggghts.h"
#include "error.h"
#include "comm.h"

#include <algorithm>
#include <cmath>
#include <vector>

using MODIFIED_ANDREW_AUX::Circle;
using MODIFIED_ANDREW_AUX::Point;
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ModifiedAndrew::ModifiedAndrew(LAMMPS *lmp) : Pointers(lmp) {
  npoints_per_circle_ = 6;
  directions_ = new double*[npoints_per_circle_];
  for (int i = 0; i < npoints_per_circle_; ++i)
  {
    directions_[i] = new double[2];
  }
  double angle = (2.0*M_PI)/((double)npoints_per_circle_);
  for (int i = 0; i < npoints_per_circle_; ++i)
  {
    directions_[i][0] = sin(((double)i)*angle);
    directions_[i][1] = cos(((double)i)*angle);
  }

}
ModifiedAndrew::~ModifiedAndrew(){
  for (int i = 0; i < npoints_per_circle_; ++i)
    delete []directions_[i];
  delete []directions_;
}

/* ---------------------------------------------------------------------- */

double ModifiedAndrew::area()
{
    double A;
    vector<Point> hull_c = convex_hull(contacts_);

    // multi-proc case
    if(1 < comm->nprocs) //(1)//
    {
      
        // perform allgather
        double *data = 0,*data0 = 0;
        int ndata = 0, ndata0 = 0;

        if(hull_c.size() > 2)
            hull_c.erase(hull_c.begin());

        // copy from STL vector to C array to safely use MPI
        ndata = construct_data(hull_c,data);
        ndata0 = MPI_Gather0_Vector(data,ndata,data0,world);

        // proc0 calculates convex hull
        if(0 == comm->me) {

            vector<Point> hull_c_allreduced = construct_hull_c_all(data0,ndata0);
            vector<Point> hull_c_global = convex_hull(hull_c_allreduced);

            if (hull_c_global.size() < 3)
                A = 0.;
            else
                A = area(hull_c_global);
        }
        
        if(data) delete []data;
        if(data0) delete []data0;

        // broadcast result
        MPI_Bcast(&A,1,MPI_DOUBLE,0,world);
    }
    // single proc case
    else
    {
        if (hull_c.size() < 3)
           A = 0.;
        else
            A = area(hull_c);
    }

    // clear contact data
    contacts_.clear();

    return A;
}

/* ---------------------------------------------------------------------- */

void ModifiedAndrew::add_contact(Circle c){
  // generate npoints_per_circle_ points for this circle
  Point p;
  for (int i = 0; i < npoints_per_circle_; ++i)
  {
    p.x = c.x+directions_[i][0]*c.r;
    p.y = c.y+directions_[i][1]*c.r;
    contacts_.push_back(p);
  }
}

/* ---------------------------------------------------------------------- */

int ModifiedAndrew::construct_data(vector<Point> hull_c, double *&data)
{
    int size = hull_c.size();
    int datasize = 2*size;
    data = new double[datasize];

    for(int i = 0; i < size; i++)
    {
        data[2*i+0] = hull_c[i].x;
        data[2*i+1] = hull_c[i].y;
    }

    return datasize;
}

/* ---------------------------------------------------------------------- */

vector<Point> ModifiedAndrew::construct_hull_c_all(double *data0, int ndata0)
{
    vector<Point> result;

    Point c;

    for(int i = 0; i < ndata0/2; i++)
    {
        c.x = data0[i*2+0];
        c.y = data0[i*2+1];
        result.push_back(c);
    }

    return result;
}

/* ---------------------------------------------------------------------- */

// area of triangle spanned by three points p,m,q
// based on cross product
double ModifiedAndrew::area(Point &p, Point &m, Point &q){
return 0.5*(-(m.y*p.x) + m.x*p.y +
      m.y*q.x - p.y*q.x - m.x*q.y + p.x*q.y);
}

/* ---------------------------------------------------------------------- */

// area of convex hull
double ModifiedAndrew::area(vector<Point> H){
  double a = 0.0;
  Point m = mean_point(H);

  for (size_t i = 1; i < H.size(); i++){
    a += area(H[i-1],m,H[i]);
  }

  return a;
}

/* ---------------------------------------------------------------------- */

// 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
// Returns a positive value, if OAB makes a counter-clockwise turn,
// negative for clockwise turn, and zero if the points are collinear.
double ModifiedAndrew::cross(Point O, Point A, Point B)
{
  return (A.x - O.x)*(B.y - O.y) - (A.y - O.y)*(B.x - O.x);
}

/* ---------------------------------------------------------------------- */

// find the mean of a vector of points (excluding the first)
Point ModifiedAndrew::mean_point(vector<Point> P){
  Point mean;

  double x_accum = 0.0, y_accum = 0.0;

  for (size_t i = 1; i < P.size(); i++){
    x_accum += P[i].x;
    y_accum += P[i].y;
  }

  mean.x = x_accum/(static_cast<double>(P.size())-1.);
  mean.y = y_accum/(static_cast<double>(P.size())-1.);

  return mean;
}

/* ---------------------------------------------------------------------- */

// Returns a list of points on the convex hull in counter-clockwise order.
// Note: the last circle in the returned list is the same as the first one.
vector<Point> ModifiedAndrew::convex_hull(vector<Point> P)
{
    int n = P.size(), k = 0;
    vector<Point> H(2*n);

    // Sort points lexicographically
    sort(P.begin(), P.end());

    // Build lower hull
    for (int i = 0; i < n; i++) {
        while (k >= 2 && cross(H[k-2], H[k-1], P[i]) <= 0) k--;
        H[k++] = P[i];
    }

    // Build upper hull
    for (int i = n-2, t = k+1; i >= 0; i--) {
        while (k >= t && cross(H[k-2], H[k-1], P[i]) <= 0) k--;
        H[k++] = P[i];
    }

    H.resize(k);
    return H;
}
