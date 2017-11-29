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

    Alexander Podlozhnyuk
    Copyright 2015-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifndef LMP_SUPERQUADRIC_H
#define LMP_SUPERQUADRIC_H

#include <cstddef>
#include "math_extra.h"

struct Superquadric
{ //interface structure for easy coding
  double *center;
  double *quat;
  double *shape;
  double *blockiness;
  double gradient[3];
  double hessian[9];
  double rotation_matrix[9];
  double koef;
  double shape_inv[3];
  bool isEllipsoid;
  bool isCylinder;
  bool useIntBlockiness;

  void local2global(const double *input_coord, double *result);
  void global2local(const double *input_coord, double *result);
  void rotate_local2global(const double *input_coord, double *result);
  void rotate_global2local(const double *input_coord, double *result);

  void apply_homothety(double scale_factor, double *input_point, double *output_point); //shrink or blow a particle
  double shape_function_local(const double *point);  //calculates shape function value at the point, defined in the local(particle based) reference frame
  double shape_function_global(const double *point);  //calculates shape function value at the point, defined in the global reference frame
  void shape_function_gradient_local(const double *input_coord, double *result);  //calculates gradient of the shape function at the point, defined in the local(particle based) reference frame
  void shape_function_gradient_global(const double *input_coord, double *result);  //calculates gradient of the shape function at the point, defined in the global reference frame
  void shape_function_hessian_local(const double *input_coord, double *result);  //calculates hessian of the shape function at the point, defined in the local(particle based) reference frame
  void shape_function_hessian_global(const double *input_coord, double *result);  //calculates hessian of the shape function at the point, defined in the global reference frame
  void map_point(const double *input_coord, double *result);  //
  void map_point_local(double *input_coord, double scale);
  void reference_point(int iphi, int nphi, int itheta, int ntheta, double *point);
  void pre_initial_estimate(const double *input_point, int nphi, int ntheta, int *iphi, int *itheta);
  bool plane_intersection(const double *normal, const double *x0, double *result_point, double *point_of_lowest_potential);  //finds the point of maximal penetration on the particle surface and the point of minimal shape function value on the surface
  double line_intersection(const double *pointA, const double *pointB, double *closestPoint);  //finds a point on edge incl.corners with a minimal value of particle shape function (point of lowest potential)
  bool edge_intersection(const double *pointA, const double *pointB);  //checks superquadric-edge intersection (only yes or no)
  void tensor_quat_rotate(double *tensor, double *result);  //rotate tensor
  double calc_curvature_coefficient(int curvature_radius, const double *input_point);  //calculates the mean curvature radius of a particle surface at a given point
  void shape_function_props_local(const double *input_coord, double *f, double *grad, double *hess);  //calculates shape function value, gradient and hessian at the point defined in local basis by one function
  void shape_function_props_global(const double *input_coord, double *f, double *grad, double *hess);  //calculates shape function value, gradient and hessian at the point defined in global basis by one function

  void set_shape(double a, double b, double c);  //sets particle shape parameters
  void set_blockiness(double n1, double n2);  //sets particle blockiness parameters
  void calc_koef();

  double surface_line_intersection(bool use_alhpa, const double *start_point, const double *normal_vector, double alpha1, double *result);  //calculates the intersection point between a line and particle surface with the Newton's method
  double surface_line_intersection(const int max_num_iters, bool use_alhpa, const double *start_point, const double *direction_vector, double alpha1, double *result);
  double surface_line_intersection1(const double *start_point, const double *direction_vector, double *result);
  void ellipsoid_line_intersection(const double *start_point, const double *direction, double& alpha1, double& alpha2);
  void ellipsoid_line_intersection_local(const double *start_point_local, const double *direction_vector_local, double& alpha1, double& alpha2);
  double point_surface_projection(const int max_num_iters, bool use_alpha_ellipsoid, const double *start_point, double alpha1, double *result);
  double calc_F(double *p, double *p0, double l, double *F);

  Superquadric():
    center(NULL),
    quat(NULL),
    shape(NULL),
    blockiness(NULL)
  {
    isEllipsoid = false;
    isCylinder = false;
    useIntBlockiness = false;
    koef = 1.0;
  }
  Superquadric(double *center_, double *quat_, double *shape_, double *blockiness_);
  void set(double *center_, double *quat_, double *shape_, double *blockiness_);
};

#endif

