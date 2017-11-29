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

#include "superquadric.h"
#include <algorithm>

#ifdef SUPERQUADRIC_ACTIVE_FLAG
#include "math_extra_liggghts_superquadric.h"

#define PHI_INV 0.61803398874989479 //: = 1/phi = 2.0/(1.0+sqrt(5.0))

Superquadric::Superquadric(double *center_, double *quat_, double *shape_, double *blockiness_)
{
  set(center_, quat_, shape_, blockiness_);
}

//rotate from body fixed to global coordinate system
inline void Superquadric::rotate_local2global(const double *input_coord, double *result)
{
  MathExtraLiggghtsNonspherical::matvec(rotation_matrix, input_coord, result);
}

//rotate from global coordinate to body fixed system
inline void Superquadric::rotate_global2local(const double *input_coord, double *result)
{
  MathExtraLiggghtsNonspherical::transpose_matvec(rotation_matrix, input_coord, result);
}

//transform point in local coordinates to global coordinates
void Superquadric::local2global(const double *input_coord, double *result)
{
   double x_temp[3];
   rotate_local2global(input_coord, x_temp);
   MathExtra::add3(x_temp, center, result); //result = x_temp + center
 }

 //transform point in global coordinates to local particle based coordinates
void Superquadric::global2local(const double *input_coord, double *result)
{
   double x_temp[3];
   MathExtra::sub3(input_coord, center, x_temp); //x_temp = x -  center
   rotate_global2local(x_temp, result); //rotate x_temp from global to local
 }

void Superquadric::set(double *center_, double *quat_, double *shape_, double *blockiness_)
{
  center = center_;
  quat = quat_;
  shape = shape_;
  blockiness = blockiness_;
  MathExtraLiggghtsNonspherical::quat_to_mat(quat, rotation_matrix);
  calc_koef();
  if(shape != NULL) {
    shape_inv[0] = 1.0 / shape[0];
    shape_inv[1] = 1.0 / shape[1];
    shape_inv[2] = 1.0 / shape[2];
  }
  isEllipsoid = MathExtraLiggghts::compDouble(blockiness[0], 2.0, 1e-2) and MathExtraLiggghts::compDouble(blockiness[1], 2.0, 1e-2);
  isCylinder = !MathExtraLiggghts::compDouble(blockiness[0], 2.0, 1e-2) and MathExtraLiggghts::compDouble(blockiness[1], 2.0, 1e-2);

  if(isEllipsoid)
    useIntBlockiness = true;
  else {
    if(MathExtraLiggghtsNonspherical::isInteger(blockiness[0]) and MathExtraLiggghtsNonspherical::isInteger(blockiness[1]))
      useIntBlockiness = true;
    else
      useIntBlockiness = false;
  }
}

//value of the particle shape function at x in local reference frame
double Superquadric::shape_function_local(const double *input_coord)
{
  const double n1 = blockiness[0];
  const double n2 = blockiness[1];
  double f;
  if(useIntBlockiness) {
    const int n1_int = floor(blockiness[0] + 1e-2);
    const int n2_int = floor(blockiness[1] + 1e-2);
    if(n1_int == n2_int) {
      f =  MathExtraLiggghtsNonspherical::pow_abs_int(input_coord[0]* shape_inv[0], n2_int) +
           MathExtraLiggghtsNonspherical::pow_abs_int(input_coord[1]* shape_inv[1], n2_int) +
           MathExtraLiggghtsNonspherical::pow_abs_int(input_coord[2]* shape_inv[2], n1_int) - 1.0;
    } else {
      f = (MathExtraLiggghtsNonspherical::pow_abs(
           MathExtraLiggghtsNonspherical::pow_abs_int(input_coord[0]* shape_inv[0], n2_int) +
           MathExtraLiggghtsNonspherical::pow_abs_int(input_coord[1]* shape_inv[1], n2_int), n1/n2) +
           MathExtraLiggghtsNonspherical::pow_abs_int(input_coord[2]* shape_inv[2], n1_int) - 1.0);
    }
  } else {
    f = (MathExtraLiggghtsNonspherical::pow_abs(
         MathExtraLiggghtsNonspherical::pow_abs(input_coord[0]* shape_inv[0], n2) +
         MathExtraLiggghtsNonspherical::pow_abs(input_coord[1]* shape_inv[1], n2), n1/n2) +
         MathExtraLiggghtsNonspherical::pow_abs(input_coord[2]* shape_inv[2], n1) - 1.0);
  }
  return koef*f;
}

//gradient of the shape function, defined in canonical frame
void Superquadric::shape_function_gradient_local(const double *input_coord, double *result)
{
  const double n1 = blockiness[0];
  const double n2 = blockiness[1];

  const double a = shape_inv[0];
  const double b = shape_inv[1];
  const double c = shape_inv[2];

  //straightforward calculation of the shape function 1st derivatives

  const double xa = input_coord[0] * a;
  const double yb = input_coord[1] * b;
  const double zc = input_coord[2] * c;

  double xan21, ybn21, zcn11;
  if(useIntBlockiness) {
    const int n1_int = floor(n1 + 1e-2);
    const int n2_int = floor(n2 + 1e-2);
    xan21 = MathExtraLiggghtsNonspherical::pow_abs_int(xa, n2_int - 1);
    ybn21 = MathExtraLiggghtsNonspherical::pow_abs_int(yb, n2_int - 1);
    zcn11 = MathExtraLiggghtsNonspherical::pow_abs_int(zc, n1_int - 1);
  } else {
    xan21 = MathExtraLiggghtsNonspherical::pow_abs(xa, n2 - 1.0);
    ybn21 = MathExtraLiggghtsNonspherical::pow_abs(yb, n2 - 1.0);
    zcn11 = MathExtraLiggghtsNonspherical::pow_abs(zc, n1 - 1.0);
  }

  if(MathExtraLiggghts::compDouble(n1, n2, 1e-2)) {
    result[0] = n1 * (koef * a) * MathExtraLiggghtsNonspherical::sign(xa) * xan21;
    result[1] = n1 * (koef * b) * MathExtraLiggghtsNonspherical::sign(yb) * ybn21;
    result[2] = n1 * (koef * c) * MathExtraLiggghtsNonspherical::sign(zc) * zcn11;
  } else {
    const double xy_term = xan21 * fabs(xa) + ybn21 * fabs(yb);
    double xy_term_powed1 = MathExtraLiggghtsNonspherical::pow_abs(xy_term, n1 / n2 - 1.0);
    result[0] = n1 * (koef * a) * MathExtraLiggghtsNonspherical::sign(xa) * xan21 * xy_term_powed1;
    result[1] = n1 * (koef * b) * MathExtraLiggghtsNonspherical::sign(yb) * ybn21 * xy_term_powed1;
    result[2] = n1 * (koef * c) * MathExtraLiggghtsNonspherical::sign(zc) * zcn11;
  }
}

//hessian (matrix of the 2nd derivatives) of the shape function in local frame
void Superquadric::shape_function_hessian_local(const double *input_coord, double *result)
{
  const double n1 = blockiness[0];
  const double n2 = blockiness[1];

  const double a = shape_inv[0];
  const double b = shape_inv[1];
  const double c = shape_inv[2];

  const double xa = input_coord[0] * a;
  const double yb = input_coord[1] * b;
  const double zc = input_coord[2] * c;

//straightforward calculation of the 2nd derivatives
  double xan22, ybn22, zcn12;
  if(useIntBlockiness) {
    const int n1_int = floor(n1 + 1e-2);
    const int n2_int = floor(n2 + 1e-2);
    xan22 = MathExtraLiggghtsNonspherical::pow_abs_int(xa, n2_int - 2);
    ybn22 = MathExtraLiggghtsNonspherical::pow_abs_int(yb, n2_int - 2);
    zcn12 = MathExtraLiggghtsNonspherical::pow_abs_int(zc, n1_int - 2);
  } else {
    xan22 = MathExtraLiggghtsNonspherical::pow_abs(xa, n2 - 2.0);
    ybn22 = MathExtraLiggghtsNonspherical::pow_abs(yb, n2 - 2.0);
    zcn12 = MathExtraLiggghtsNonspherical::pow_abs(zc, n1 - 2.0);
  }
  const double xan21 = xan22 * fabs(xa);  // = MathExtraLiggghtsNonspherical::pow_abs(xa, n2 - 1.0)
  const double ybn21 = ybn22 * fabs(yb);  // = MathExtraLiggghtsNonspherical::pow_abs(yb, n2 - 1.0)

  if(MathExtraLiggghts::compDouble(n1, n2, 1e-2)) { //MathExtraLiggghtsNonspherical::pow_abs(xy_term, n1/n2 - 1.0);
    result[0] = (koef * a) * (n1 * (n2 - 1.0) * xan22) * a;
    result[4] = (koef * b) * (n1 * (n2 - 1.0) * ybn22) * b;
    result[8] = (koef * c) * (n1 * (n1 - 1.0) * zcn12) * c;
    result[1] = result[3] = result[2] = result[6] = result[5] = result[7] = 0.0;
  } else {
    const double xy_term = xan21 * fabs(xa) + ybn21 * fabs(yb);  // MathExtraLiggghtsNonspherical::pow_abs(xa, n2) + MathExtraLiggghtsNonspherical::pow_abs(yb, n2)
    double xy_term_powed2;
    if(MathExtraLiggghtsNonspherical::isInteger(n1/n2)) {
      int n1_over_n2_int = floor(n1/n2 + 1e-2);
      xy_term_powed2 = MathExtraLiggghtsNonspherical::pow_abs_int(xy_term, n1_over_n2_int - 2.0);
    } else
      xy_term_powed2 = MathExtraLiggghtsNonspherical::pow_abs(xy_term, n1/n2 - 2.0);
    const double xy_term_powed1 = xy_term_powed2 * xy_term;

    result[0] = (koef * a) * (n1 * (n2 - 1.0) * xan22 * xy_term_powed1 +
        (n1 - n2) * n1 * xan21 * xan21 * xy_term_powed2) * a;
    result[4] = (koef * b) * (n1 * (n2 - 1.0) * ybn22 * xy_term_powed1 +
        (n1 - n2) * n1 * ybn21 * ybn21 * xy_term_powed2) * b;
    result[1] = koef * (n1 - n2) * n1 * MathExtraLiggghtsNonspherical::sign(xa) * MathExtraLiggghtsNonspherical::sign(yb) * (xan21 * a) * (ybn21 * b) * xy_term_powed2;

    result[3] = result[1];
    result[8] = (koef * c) * (n1 * (n1 - 1.0) * zcn12) * c;
    result[2] = result[6] = result[5] = result[7] = 0.0;
  }
}

//compute shape function value, gradient and hessian matrix altogether in one function
void Superquadric::shape_function_props_local(const double *input_coord, double *f, double *grad, double *hess)
{
  const double n1 = blockiness[0];
  const double n2 = blockiness[1];

  const double a = shape_inv[0];
  const double b = shape_inv[1];
  const double c = shape_inv[2];

  const double xa = input_coord[0] * a;
  const double yb = input_coord[1] * b;
  const double zc = input_coord[2] * c;

  double xan22, ybn22, zcn12;
  if(useIntBlockiness) {
    const int n1_int = floor(n1 + 1e-2);
    const int n2_int = floor(n2 + 1e-2);
    xan22 = MathExtraLiggghtsNonspherical::pow_abs_int(xa, n2_int - 2);
    ybn22 = MathExtraLiggghtsNonspherical::pow_abs_int(yb, n2_int - 2);
    zcn12 = MathExtraLiggghtsNonspherical::pow_abs_int(zc, n1_int - 2);
  } else {
    xan22 = MathExtraLiggghtsNonspherical::pow_abs(xa, n2 - 2.0);
    ybn22 = MathExtraLiggghtsNonspherical::pow_abs(yb, n2 - 2.0);
    zcn12 = MathExtraLiggghtsNonspherical::pow_abs(zc, n1 - 2.0);
  }
  const double xan21 = xan22 * fabs(xa);  // = MathExtraLiggghtsNonspherical::pow_abs(xa, n2 - 1.0)
  const double ybn21 = ybn22 * fabs(yb);  // = MathExtraLiggghtsNonspherical::pow_abs(yb, n2 - 1.0)
  const double zcn11 = zcn12 * fabs(zc);  // = MathExtraLiggghtsNonspherical::pow_abs(zc, n1 - 1.0)
  const double xan2 = xan21 * fabs(xa);  // = MathExtraLiggghtsNonspherical::pow_abs(xa, n2)
  const double ybn2 = ybn21 * fabs(yb);  // = MathExtraLiggghtsNonspherical::pow_abs(yb, n2)
  const double zcn1 = zcn11 * fabs(zc);  // = MathExtraLiggghtsNonspherical::pow_abs(zc, n1)

  if(MathExtraLiggghts::compDouble(n1, n2, 1e-2)) { //MathExtraLiggghtsNonspherical::pow_abs(xy_term, n1/n2 - 1.0);
    *f = koef * (xan2 + ybn2 + zcn1 - 1.0);
    grad[0] = koef * a * n1 * MathExtraLiggghtsNonspherical::sign(xa) * xan21;
    grad[1] = koef * b * n1 * MathExtraLiggghtsNonspherical::sign(yb) * ybn21;
    grad[2] = koef * c * n1 * MathExtraLiggghtsNonspherical::sign(zc) * zcn11;

    if(hess != NULL) {
      hess[0] = (koef * a) * (n1 * (n2 - 1.0) * xan22) * a;
      hess[4] = (koef * b) * (n1 * (n2 - 1.0) * ybn22) * b;
      hess[8] = (koef * c) * (n1 * (n1 - 1.0) * zcn12) * c;
      hess[1] = hess[3] = hess[2] = hess[6] = hess[5] = hess[7] = 0.0;
    }
  }
  else {
    *f = koef * (MathExtraLiggghtsNonspherical::pow_abs(xan2 + ybn2, n1/n2) + zcn1 - 1.0);
    const double xy_term = xan21 * fabs(xa) + ybn21 * fabs(yb);  // MathExtraLiggghtsNonspherical::pow_abs(xa, n2) + MathExtraLiggghtsNonspherical::pow_abs(yb, n2)
    double xy_term_powed2;
    if(MathExtraLiggghtsNonspherical::isInteger(n1/n2)) {
      int n1_over_n2_int = floor(n1/n2 + 1e-2);
      xy_term_powed2 = MathExtraLiggghtsNonspherical::pow_abs_int(xy_term, n1_over_n2_int - 2.0);
    } else
      xy_term_powed2 = MathExtraLiggghtsNonspherical::pow_abs(xy_term, n1/n2 - 2.0);
    double xy_term_powed1 = xy_term_powed2 * xy_term;

    grad[0] = koef * a * n1 * MathExtraLiggghtsNonspherical::sign(xa) * xan21 * xy_term_powed1;
    grad[1] = koef * b * n1 * MathExtraLiggghtsNonspherical::sign(yb) * ybn21 * xy_term_powed1;
    grad[2] = koef * c * n1 * MathExtraLiggghtsNonspherical::sign(zc) * zcn11;

    if(hess != NULL) {
      hess[0] = (koef * a) * (n1 * (n2 - 1.0) * xan22 * xy_term_powed1 +
            (n1 - n2) * n1 * xan21 * xan21 * xy_term_powed2) * a;
      hess[4] = (koef * b) * (n1 * (n2 - 1.0) * ybn22 * xy_term_powed1 +
            (n1 - n2) * n1 * ybn21 * ybn21 * xy_term_powed2) * b;
      hess[1] =  koef * (n1 - n2) * n1 * MathExtraLiggghtsNonspherical::sign(xa) * MathExtraLiggghtsNonspherical::sign(yb) * (xan21 * a) * (ybn21 * b) * xy_term_powed2;

      hess[3] = hess[1];
      hess[8] = (koef * c) * (n1 * (n1 - 1.0) * zcn12) * c;
      hess[2] = hess[6] = hess[5] = hess[7] = 0.0;
    }
  }
}

//shape function props (value, gradient and hessian matrix) in global coordinates
void Superquadric::shape_function_props_global(const double *input_coord, double *f, double *grad, double *hess)
{
  double coord_local[3];
  global2local(input_coord, coord_local);  //move to particle based reference frame
  double grad_local[3];
  if(hess != NULL) {  //if we need hessian calculation
    double hess_local[9];
    shape_function_props_local(coord_local, f, grad_local, hess_local);
    tensor_quat_rotate(hess_local, hess);  // transform hessian to global frame
  }
  else  //if we don't need hessian calculation
    shape_function_props_local(coord_local, f, grad_local, NULL);
  rotate_local2global(grad_local, grad);  //transform gradient in global frame
}

//value of the particle shape function at x in global reference frame
double Superquadric::shape_function_global(const double *input_coord)
{
  double coord_local[3];
  global2local(input_coord, coord_local);  //move to particle based (local) reference frame
  return shape_function_local(coord_local);  //calculate shape function value in local reference frame
}

//gradient of the shape function in global frame
void Superquadric::shape_function_gradient_global(const double *input_coord, double *result)
{
  double x_local[3], grad_local[3];
  global2local(input_coord, x_local);  //transform point x to canonical frame
  shape_function_gradient_local(x_local, grad_local);  //calculate gradient in canonical frame
  rotate_local2global(grad_local, result);  //transform gradient in global frame
}

//hessian (matrix of the 2nd derivatives) of the shape function in global frame
void Superquadric::shape_function_hessian_global(const double *input_coord, double result[9])
{
  double hessian_local[9], x_local[3];
  global2local(input_coord, x_local);  //transform point x to canonical frame
  shape_function_hessian_local(x_local, hessian_local);  //calculate hessian in canonical frame
  tensor_quat_rotate(hessian_local, result);  // transform hessian to global frame
}

//map point onto particle surface
void Superquadric::map_point(const double *input_coord, double *result)
{
  double local_point[3];
  global2local(input_coord, local_point);
  double f = shape_function_local(local_point);
  double alpha = 1.0 / MathExtraLiggghtsNonspherical::pow_abs(f/koef + 1.0, 1.0/blockiness[0]);
  for(int i = 0; i < 3; i++)
    local_point[i] *= alpha;
  local2global(local_point, result);
}

//map point onto particle surface
void Superquadric::map_point_local(double *input_coord, double scale)
{
  double f = shape_function_local(input_coord);
  double alpha = 1.0 / MathExtraLiggghtsNonspherical::pow_abs(f/koef + 1.0, 1.0/blockiness[0]);
  for(int i = 0; i < 3; i++)
    input_coord[i] *= alpha*scale;
}

//calculate reference point
//old code, obsolete
void Superquadric::reference_point(int iphi, int nphi, int itheta, int ntheta, double *point)
{
  double phi_start = 1.0e-2;
  double theta_start = 1.0e-2;
  double dphi = (2.0*M_PI - 2.0*phi_start) / static_cast<double>(nphi - 1);
  double dtheta =  (M_PI - 2.0*theta_start) / static_cast<double>(ntheta - 1);
  double phi = dphi*static_cast<double>(iphi) + phi_start;
  double theta = dtheta*static_cast<double>(itheta) + theta_start;
  double r = std::max(std::max(shape[0], shape[1]), shape[2]);
  double circle_point[3], local_point[3];
  circle_point[0] = r * sin(theta)*cos(phi);
  circle_point[1] = r * sin(theta)*sin(phi);
  circle_point[2] = r * cos(theta);
  local2global(circle_point, local_point);
  map_point(local_point, point);
}

//old code, obsolete
void Superquadric::pre_initial_estimate(const double *input_point, int nphi, int ntheta, int *iphi, int *itheta)
{
  double phi_start = 1.0e-2;
  double theta_start = 1.0e-2;
  double dphi = (2.0*M_PI -2.0*phi_start) / static_cast<double>(nphi - 1);
  double dtheta =  (M_PI - 2.0*theta_start) / static_cast<double>(ntheta - 1);
  double phi, theta;
  double local_point[3];
  global2local(input_point, local_point);
  double r = MathExtra::len3(local_point);
  double cos_theta = local_point[2] / r;
  double sin_theta = sqrt(1.0 - cos_theta*cos_theta);
  if(sin_theta < sin(dtheta)) {
    if(cos_theta > 0.0)
      theta = theta_start;
    else
      theta = M_PI - theta_start;
    phi = M_PI;
  } else {
    theta = acos(cos_theta);
    double cos_phi = local_point[0] / sin_theta / r;
    double sin_phi = local_point[1] / sin_theta / r;
    if(sin_phi > 0.0)
      phi = acos(cos_phi);
    else
      phi = 2.0 * M_PI - acos(cos_phi);
  }
  *iphi =   static_cast<int>(round((phi - phi_start)/ dphi + nphi - 1)) % (nphi - 1);
  *itheta = static_cast<int>(round((theta - theta_start) / dtheta + ntheta - 1)) % (ntheta - 1);
}

//shrink or blow particle
void Superquadric::apply_homothety(double scale_factor, double *input_point, double *result_point)
{
  for(int i = 0; i < 3; i++)
    shape[i] *= scale_factor;
  double delta[3];
  if(input_point) {
    MathExtra::sub3(input_point, center, delta);
    MathExtra::scale3(scale_factor, delta);
    MathExtra::add3(center, delta, result_point);
  }
}

//particle-wall intersection check
//finds the point of maximal penetration on the particle surface and the point of minimal shape function value on the surface
bool Superquadric::plane_intersection(const double *normal_global, const double *x0_global, double *point_of_maximal_penetration, double *point_of_lowest_potential) {
  double point_of_maximal_penetration_local[3], x0_local[3], normal_local[3];
  rotate_global2local(normal_global, normal_local);
  global2local(x0_global, x0_local);
  if(LAMMPS_NS::vectorDot3D(normal_local, x0_local) < 0.0) //check that normal vector aligned outwards
    LAMMPS_NS::vectorNegate3D(normal_local);  //make the normal aligned outwards

  const double n1 = blockiness[0];
  const double n2 = blockiness[1];
  const double nx = normal_local[0];
  const double ny = normal_local[1];
  const double nz = normal_local[2];
  const double a = shape[0];
  const double b = shape[1];
  const double c = shape[2];

  double point_of_lowest_potential_local[3]; //point on wall with the lowest value of shape function of a given particle

  const double d = LAMMPS_NS::vectorDot3D(normal_local, x0_local);

  if(MathExtraLiggghts::compDouble(nx, 0.0, 1e-10) and MathExtraLiggghts::compDouble(ny, 0.0, 1e-10)) {
    point_of_maximal_penetration_local[0] = point_of_maximal_penetration_local[1] = 0.0;
    point_of_maximal_penetration_local[2] = 1.0;

    point_of_lowest_potential_local[0] = point_of_lowest_potential_local[1] = 0.0;
    point_of_lowest_potential_local[2] = d / nz;
  } else {

    if(fabs(nx) > fabs(ny)) {
      double alpha = MathExtraLiggghtsNonspherical::pow_abs((ny * b) / (nx * a), 1.0 / (n2 - 1.0));
      double gamma1 = 1.0 + MathExtraLiggghtsNonspherical::pow_abs(alpha, n2);
      double gamma = MathExtraLiggghtsNonspherical::pow_abs(gamma1, n1/n2 - 1.0);
      double beta = MathExtraLiggghtsNonspherical::pow_abs((nz * c) / (nx * a) * gamma, 1.0/(n1 - 1.0));

      point_of_maximal_penetration_local[0] = 1.0 / MathExtraLiggghtsNonspherical::pow_abs(MathExtraLiggghtsNonspherical::pow_abs(gamma1, n1/n2) + MathExtraLiggghtsNonspherical::pow_abs(beta, n1), 1.0 / n1);
      point_of_maximal_penetration_local[1] = alpha * point_of_maximal_penetration_local[0];
      point_of_maximal_penetration_local[2] = beta * point_of_maximal_penetration_local[0];

      point_of_lowest_potential_local[0] = fabs(d) / (fabs(nx)*a + alpha * fabs(ny) * b  + beta * fabs(nz) * c);
      point_of_lowest_potential_local[1] = alpha * point_of_lowest_potential_local[0];
      point_of_lowest_potential_local[2] =  beta * point_of_lowest_potential_local[0];
    } else {
      double alpha = MathExtraLiggghtsNonspherical::pow_abs((nx * a) / (ny * b), 1.0 / (n2 - 1.0));
      double gamma1 = 1.0 + MathExtraLiggghtsNonspherical::pow_abs(alpha, n2);
      double gamma = MathExtraLiggghtsNonspherical::pow_abs(gamma1, n1/n2 - 1.0);
      double beta = MathExtraLiggghtsNonspherical::pow_abs((nz * c) / (ny * b) * gamma, 1.0/(n1 - 1.0));
      point_of_maximal_penetration_local[1] = 1.0 / MathExtraLiggghtsNonspherical::pow_abs(MathExtraLiggghtsNonspherical::pow_abs(gamma1, n1/n2) + MathExtraLiggghtsNonspherical::pow_abs(beta, n1), 1.0 / n1);
      point_of_maximal_penetration_local[0] = alpha * point_of_maximal_penetration_local[1];
      point_of_maximal_penetration_local[2] = beta * point_of_maximal_penetration_local[1];

      point_of_lowest_potential_local[1] = fabs(d) / (fabs(ny) * b + alpha * fabs(nx) * a + beta * fabs(nz) * c);
      point_of_lowest_potential_local[0] = alpha * point_of_lowest_potential_local[1];
      point_of_lowest_potential_local[2] =  beta * point_of_lowest_potential_local[1];
    }
    for(int i = 0; i < 3; i++)
      point_of_lowest_potential_local[i] *= shape[i] * MathExtraLiggghtsNonspherical::sign(normal_local[i]);
  }

  for(int i = 0; i < 3; i++)
    point_of_maximal_penetration_local[i] *= shape[i] * MathExtraLiggghtsNonspherical::sign(normal_local[i]);

  double ff = shape_function_local(point_of_lowest_potential_local);
  local2global(point_of_maximal_penetration_local, point_of_maximal_penetration);
  local2global(point_of_lowest_potential_local, point_of_lowest_potential);
  if(ff < 0.0) //particle and wall intersect
    return true;
  else //particle and wall do not intersect
    return false;
}

//find a point on edge incl.corners with a minimal value of particle shape function (point of lowest potential)
double Superquadric::line_intersection(const double *pointA, const double *pointB, double *closestPoint)
{
  double A[3], B[3];
  double closestPointLocal[3], closestPointLocal1[3], closestPointLocal2[3];
  //using "golden section" search algorithm
  //https://en.wikipedia.org/wiki/Golden_section_search

  global2local(pointA, A);
  global2local(pointB, B);

  double a = 0.0;
  double b = 1.0;
  double eps = fabs(b-a);
  double alpha1 = b - (b - a)*PHI_INV;
  double alpha2 = a + (b - a)*PHI_INV;
  for(int k = 0; k < 3; k++) {
    closestPointLocal1[k] = A[k] +  alpha1*(B[k] - A[k]);
    closestPointLocal2[k] = A[k] +  alpha2*(B[k] - A[k]);
  }
  double f1 = shape_function_local(closestPointLocal1);
  double f2 = shape_function_local(closestPointLocal2);
  while(eps > 1e-8) {
    alpha1 = b - (b - a)*PHI_INV;
    alpha2 = a + (b - a)*PHI_INV;
    for(int k = 0; k < 3; k++) {
      closestPointLocal1[k] = A[k] +  alpha1*(B[k] - A[k]);
      closestPointLocal2[k] = A[k] +  alpha2*(B[k] - A[k]);
    }
    f1 = shape_function_local(closestPointLocal1);
    f2 = shape_function_local(closestPointLocal2);
    if(f1 >= f2) {
      a = alpha1;
    } else
      b = alpha2;
    eps = fabs(b-a);
  }
  for(int k = 0; k < 3; k++)
    closestPointLocal[k] = 0.5*(closestPointLocal1[k] + closestPointLocal2[k]);
  double f = shape_function_local(closestPointLocal);
  double fA = shape_function_local(A);
  double fB = shape_function_local(B);
  if(fA < f) {
    LAMMPS_NS::vectorCopy3D(pointA, closestPoint);
    f = fA;
  }
  else if(fB < f) {
    LAMMPS_NS::vectorCopy3D(pointB, closestPoint);
    f = fB;
  }
  else
    local2global(closestPointLocal, closestPoint);
  return f;
}

//check superquadric-edge intersection (only yes or no)
bool Superquadric::edge_intersection(const double *pointA, const double *pointB)
{
  double A[3], B[3];
  double closestPointLocal[3], closestPointLocal1[3], closestPointLocal2[3];
  //using "golden section" search algorithm
  //https://en.wikipedia.org/wiki/Golden_section_search

  global2local(pointA, A);
  global2local(pointB, B);

  double a = 0.0;
  double b = 1.0;
  double eps = fabs(b-a);
  double alpha1 = b - (b - a)*PHI_INV;
  double alpha2 = a + (b - a)*PHI_INV;
  for(int k = 0; k < 3; k++) {
    closestPointLocal1[k] = A[k] +  alpha1*(B[k] - A[k]);
    closestPointLocal2[k] = A[k] +  alpha2*(B[k] - A[k]);
  }
  double f1 = shape_function_local(closestPointLocal1);
  double f2 = shape_function_local(closestPointLocal2);
  if(f1 < 0.0 or f2 < 0.0)
    return true;
  else {
    while(eps > 1e-4) {
      alpha1 = b - (b - a)*PHI_INV;
      alpha2 = a + (b - a)*PHI_INV;
      for(int k = 0; k < 3; k++) {
        closestPointLocal1[k] = A[k] +  alpha1*(B[k] - A[k]);
        closestPointLocal2[k] = A[k] +  alpha2*(B[k] - A[k]);
      }
      f1 = shape_function_local(closestPointLocal1);
      f2 = shape_function_local(closestPointLocal2);
      if(f1 < 0.0 or f2 < 0.0) {
        return true;
      } else {
        if(f1 >= f2) {
          a = alpha1;
        } else
          b = alpha2;
        eps = fabs(b-a);
      }
    }
    for(int k = 0; k < 3; k++)
      closestPointLocal[k] = 0.5*(closestPointLocal1[k] + closestPointLocal2[k]);
    double f = shape_function_local(closestPointLocal);
    if(f < 0.0)
      return true; //particle and edge intersect
    else
      return false; //particle and edge do not intersect
  }
}

//rotate tensor
void Superquadric::tensor_quat_rotate(double *tensor, double *result)
{
  double temp[9];
  MathExtraLiggghtsNonspherical::times3(rotation_matrix, tensor, temp); // calc temp = A * M
  MathExtraLiggghtsNonspherical::times3_transpose(temp, rotation_matrix, result); //calc result = A * M * A^T
}

//https://en.wikipedia.org/wiki/Mean_curvature
//calculates the mean curvature radius of a particle surface at a given point
double Superquadric::calc_curvature_coefficient(int curvature_radius, const double *input_point)
{
  double local_point[3],hessian[9], F[3];
  global2local(input_point, local_point);
  double f = shape_function_local(local_point);
  double alpha = 1.0 / MathExtraLiggghtsNonspherical::pow_abs(f/koef + 1.0, 1.0/blockiness[0]);
  for(int i = 0; i < 3; i++)
    local_point[i] *= alpha;
  shape_function_gradient_local(local_point, F);
  rotate_local2global(F, gradient);
  shape_function_hessian_local(local_point, hessian);
  double F_mag = LAMMPS_NS::vectorMag3D(F);
  double n[3];
  LAMMPS_NS::vectorScalarMult3D(F, 1.0/F_mag, n);
  double temp[3];
  MathExtraLiggghtsNonspherical::matvec(hessian, n, temp);
  if(curvature_radius == 0) //mean curvature
  return fabs(MathExtra::dot3(n, temp) - (hessian[0] + hessian[4] + hessian[8])) / fabs(2.0 * F_mag);
  else { //gaussian curvature
    double fx = F[0];
  double fy = F[1];
  double fz = F[2];

  double fxx = hessian[0];
  double fxy = hessian[1];
  double fxz = hessian[2];

  double fyy = hessian[4];
  double fyz = hessian[5];

  double fzz = hessian[8];

    double mat[] = {
        fxx, fxy, fxz, fx,
        fxy, fyy, fyz, fy,
        fxz, fyz, fzz, fz,
        fx,  fy,  fz, 0
    };

    double K = -MathExtraLiggghtsNonspherical::determinant_4x4(mat) / (F_mag*F_mag*F_mag*F_mag);
    return sqrt(fabs(K));
  }
}

void Superquadric::set_shape(double a, double b, double c)
{
  shape[0] = a;
  shape[1] = b;
  shape[2] = c;

  shape_inv[0] = 1.0 / shape[0];
  shape_inv[1] = 1.0 / shape[1];
  shape_inv[2] = 1.0 / shape[2];
}

void Superquadric::set_blockiness(double n1, double n2)
{
  blockiness[0] = n1;
  blockiness[1] = n2;
  calc_koef();
  isEllipsoid = MathExtraLiggghts::compDouble(blockiness[0], 2.0, 1e-2) and MathExtraLiggghts::compDouble(blockiness[1], 2.0, 1e-2);
  isCylinder = !MathExtraLiggghts::compDouble(blockiness[0], 2.0, 1e-2) and MathExtraLiggghts::compDouble(blockiness[1], 2.0, 1e-2);
  if(isEllipsoid)
    useIntBlockiness = true;
  else {
    if(MathExtraLiggghtsNonspherical::isInteger(blockiness[0]) and MathExtraLiggghtsNonspherical::isInteger(blockiness[1]))
      useIntBlockiness = true;
    else
      useIntBlockiness = false;
  }
}

void Superquadric::calc_koef() {
  koef = MathExtraLiggghtsNonspherical::pow_abs(0.5, std::max(blockiness[0], blockiness[1])-2.0);
}

//calculates the intersection point between a line and particle surface with the Newton's method
double Superquadric::surface_line_intersection(bool use_alpha_ellipsoid, const double *start_point, const double *direction_vector, double alpha_prev, double *result)
{
  return Superquadric::surface_line_intersection(20, use_alpha_ellipsoid, start_point, direction_vector, alpha_prev, result);
}

//calculates the intersection point between a line and particle surface with the Newton's method
double Superquadric::surface_line_intersection(const int max_num_iters, bool use_alpha_ellipsoid, const double *start_point, const double *direction_vector, double alpha_prev, double *result)
{
  double start_point_local[3], direction_vector_local[3], point[3];
  double grad_local[3];

  global2local(start_point, start_point_local);
  rotate_global2local(direction_vector, direction_vector_local); //rotate direction vector to the local reference frame
  double alpha = alpha_prev;
  double eps0 = 1e-14;
  double alpha_ellipsoid = 0.0;

  if(isEllipsoid or use_alpha_ellipsoid) {
    double alpha1, alpha2;
    ellipsoid_line_intersection_local(start_point_local, direction_vector_local, alpha1, alpha2);

    if(isEllipsoid) { //particle is an ellipsoid/sphere
      if(alpha1 > -eps0 and alpha2 < -eps0)
        alpha = alpha1;
      else if(alpha1 < -eps0 and alpha2 > -eps0)
        alpha = alpha2;
      else if(alpha1 < -eps0 and alpha2 < -eps0)
        alpha = std::max(alpha1, alpha2);
      else
        alpha = std::min(alpha1, alpha2);
    } else { //inscribed ellipsoid
      double f = shape_function_local(start_point_local);

      if(f > eps0)
        alpha_ellipsoid = std::min(alpha1, alpha2);
      else if (f < -eps0)
        alpha_ellipsoid = std::max(alpha1, alpha2);
      else {
        alpha_ellipsoid = alpha = alpha_prev = 0.0;
      }
    }
  }

  if(!isEllipsoid) {
    double f_der, delta;

    double a = 0.5*std::min(std::min(shape[0], shape[1]), shape[2]) / LAMMPS_NS::vectorMag3D(direction_vector);
    const double eps = 1e-16;  //small epsilon

    alpha = MathExtraLiggghtsNonspherical::clamp(alpha, -a, a);
    MathExtraLiggghtsNonspherical::yabx3D(start_point_local, alpha, direction_vector_local, point);
    double f0 = shape_function_local(point);
    if(use_alpha_ellipsoid) {
      double f1;
      double point_ellipsoid[3];
      MathExtraLiggghtsNonspherical::yabx3D(start_point_local, alpha_ellipsoid, direction_vector_local, point_ellipsoid);
      f1 = shape_function_local(point_ellipsoid);
      if(fabs(f1) < fabs(f0)) {
        alpha = alpha_ellipsoid;
        f0 = f1;
        LAMMPS_NS::vectorCopy3D(point_ellipsoid, point);
      }
    }

    if(fabs(f0) < eps) {
      local2global(point, result);
      return alpha;
    }

    for(int i = 0; i < max_num_iters; i++) {
       if(fabs(f0) < eps*koef)
         break;
       shape_function_gradient_local(point, grad_local);
       f_der = LAMMPS_NS::vectorDot3D(grad_local, direction_vector_local);

       if(fabs(f_der)<1e-10) {
         delta = 1e-10;
         double point_[3], grad_local_[3];
         MathExtraLiggghtsNonspherical::yabx3D(point, delta, direction_vector_local, point_);
         shape_function_gradient_local(point_, grad_local_);
         double f_der_ = LAMMPS_NS::vectorDot3D(grad_local_, direction_vector_local);
         if(f_der_ < 0.0)
           delta = -delta;
       }
       else
       delta = MathExtraLiggghtsNonspherical::clamp(-f0 / f_der, -a, a);
       if(fabs(delta) <= eps) {
         alpha += delta;
         MathExtraLiggghtsNonspherical::yabx3D(point, delta, direction_vector_local, point);
         break;
       }
       else {
         while(fabs(delta) > eps) {
           double point_[3];
           MathExtraLiggghtsNonspherical::yabx3D(point, delta, direction_vector_local, point_);
           double f_ = shape_function_local(point_);
           if(fabs(f_) < fabs(f0) or fabs(f_) < eps*koef) {
             f0 = f_;
             alpha += delta;
             LAMMPS_NS::vectorCopy3D(point_, point);
             break;
           } else
             delta *= 0.5;
         }
       }
    }
  }
  result[0] = start_point[0] + alpha*direction_vector[0];
  result[1] = start_point[1] + alpha*direction_vector[1];
  result[2] = start_point[2] + alpha*direction_vector[2];
  return alpha;
}

double Superquadric::calc_F(double *p, double *p0, double l, double *F)
{
  double grad_local[3];
  shape_function_gradient_local(p, grad_local);
  F[0] = shape_function_local(p);
  F[1] = l*grad_local[0] - p[0] + p0[0];
  F[2] = l*grad_local[1] - p[1] + p0[1];
  F[3] = l*grad_local[2] - p[2] + p0[2];

  return LAMMPS_NS::vectorMag4D(F);
}

double Superquadric::point_surface_projection(const int max_num_iters, bool use_alpha_ellipsoid, const double *start_point, double alpha_prev, double *result)
{
  const double tol = 1e-15;
  const double tol2 = 1e-5;
  const double size = std::max(std::max(shape[0],shape[1]),shape[2]);

  double initial_estimate[3];
  double p0[3], p[3];
  double n[3], v[3];
  double grad_local[3];
  double F[4];

  double direction_vector[3];
  shape_function_gradient_global(start_point, direction_vector);
  LAMMPS_NS::vectorNormalize3D(direction_vector);
  if(shape_function_global(start_point) > 0.0)
    LAMMPS_NS::vectorNegate3D(direction_vector); //{-7.7462629689949625e-01, -1.7969401213388331e-01, -6.0635316619524737e-01};

  surface_line_intersection(max_num_iters, use_alpha_ellipsoid, start_point, direction_vector, alpha_prev, initial_estimate);
  global2local(start_point, p0);
  global2local(initial_estimate, p);

  shape_function_gradient_local(p, grad_local);
  LAMMPS_NS::vectorSubtract3D(p, p0, n);
  MathExtra::cross3(n, grad_local, v);
  double sine = sqrt(fabs(LAMMPS_NS::vectorDot3D(v,v) / (LAMMPS_NS::vectorDot3D(grad_local, grad_local) * LAMMPS_NS::vectorDot3D(n, n))));

  double l = LAMMPS_NS::vectorDot3D(grad_local, n) / LAMMPS_NS::vectorDot3D(grad_local, grad_local);

  double hess[9];
  double delta0[] = {0.0, 0.0, 0.0, 0.0};

  double res0 = 1e+99;
  double res1 = res0;
  double res2 = res1;
  bool converged= false;

  int iters = 0;
  double magDelta = -1;
  double f = -1;
  for(int iter = 0; iter < max_num_iters; iter++) {

    shape_function_gradient_local(p, grad_local);
    shape_function_hessian_local(p, hess);

    double J[] = { grad_local[0], grad_local[1], grad_local[2], 0.0,
                   l*hess[0]-1.0, l*hess[1],     l*hess[2],     grad_local[0],
                   l*hess[3],     l*hess[4]-1.0, l*hess[5],     grad_local[1],
                   l*hess[6],     l*hess[7],     l*hess[8]-1.0, grad_local[2] };

    calc_F(p, p0, l, F);
    double delta[4];
    MathExtraLiggghtsNonspherical::GMRES<4,4>(J, F, delta0, delta);
    //printf("iter: %d, point: %e %e %e, mag(delta)=%e\n", iter, p[0],p[1],p[2],LAMMPS_NS::vectorMag4D(delta));
    magDelta = std::max(std::max(fabs(delta[0]),fabs(delta[1])),fabs(delta[2]));
    if(magDelta < tol*size)
      converged = true;

    double p_new[3] = { p[0] - delta[0],
                        p[1] - delta[1],
                        p[2] - delta[2] };
    double l_new = l - delta[3];
    double F_new[4];
    res2 = calc_F(p_new, p0, l_new, F_new);
    if(res2 < res1 || converged) {
      LAMMPS_NS::vectorCopy3D(p_new, p);
      l = l_new;
      res1 = res2;
    } else {

      double a = 0.0;
      double b = 1.0;
      double eps = fabs(b-a);
      double pointa[3], pointb[3];
      double Fa[4], Fb[4];
      double res2a, res2b;
      double la, lb;
      while(eps > 1e-12) {
        double alpha1 = b - (b - a)*PHI_INV;
        double alpha2 = a + (b - a)*PHI_INV;

        MathExtraLiggghtsNonspherical::yabx3D(p, -alpha1, delta, pointa);
        la = l - alpha1*delta[3];
        res2a = calc_F(pointa, p0, la, Fa);

        MathExtraLiggghtsNonspherical::yabx3D(p, -alpha2, delta, pointb);
        lb = l - alpha2*delta[3];
        res2b = calc_F(pointb, p0, lb, Fb);
        if(std::min(res2a, res2b) < res1) {
          if(res2a < res2b) {
            res1 = res2a;
            LAMMPS_NS::vectorCopy4D(Fa, F);
            l = la;
            LAMMPS_NS::vectorCopy3D(pointa, p);
          } else {
            res1 = res2b;
            LAMMPS_NS::vectorCopy4D(Fb, F);
            l = lb;
            LAMMPS_NS::vectorCopy3D(pointb, p);
  }
          eps = 0.0;
        } else  {
          if(res2a > res2b)
            a = alpha1;
          else
            b = alpha2;
          eps = fabs(b-a);
        }
      }
    }

    LAMMPS_NS::vectorSubtract3D(p, p0, n);
    MathExtra::cross3(grad_local, n, v);
    shape_function_gradient_local(p, grad_local);
    sine = sqrt(fabs(LAMMPS_NS::vectorDot3D(v,v) / (LAMMPS_NS::vectorDot3D(grad_local, grad_local) * LAMMPS_NS::vectorDot3D(n, n))));
    f = fabs(shape_function_local(p));
    if(sine < tol2 and f < tol)
      converged = true;

    if(converged)
      break;
    iters ++;
  }

  if(iters >= max_num_iters - 1) {
    double start_point_temp[3];
    local2global(p, start_point_temp);
    shape_function_gradient_global(start_point_temp, direction_vector);
    LAMMPS_NS::vectorNormalize3D(direction_vector);
    f = shape_function_global(start_point_temp);
    if(f > 0.0)
      LAMMPS_NS::vectorNegate3D(direction_vector);
    surface_line_intersection(max_num_iters, false, start_point_temp, direction_vector, alpha_prev, initial_estimate);
    LAMMPS_NS::vectorCopy3D(initial_estimate, result);
    f = fabs(shape_function_global(result));
    LAMMPS_NS::vectorSubtract3D(start_point, result, n);
    MathExtra::cross3(direction_vector, n, v);
    sine = sqrt(fabs(LAMMPS_NS::vectorDot3D(v,v) / (LAMMPS_NS::vectorDot3D(direction_vector, direction_vector) * LAMMPS_NS::vectorDot3D(n, n))));
  }
  else
  local2global(p, result);
  return LAMMPS_NS::pointDistance(result, start_point);
}

//TODO: reduce code duplication
double Superquadric::surface_line_intersection1(const double *start_point, const double *direction_vector, double *result)
{
  double start_point_local[3], direction_vector_local[3], point[3];
  double grad_local[3];

  global2local(start_point, start_point_local);
  LAMMPS_NS::vectorCopy3D(start_point_local, point);
  rotate_global2local(direction_vector, direction_vector_local); //rotate direction vector to the local reference frame

  double f_der, delta;
  double alpha = 0.0;
  const double eps = 1e-16;  //small epsilon

  double f0 = shape_function_local(start_point_local);
  int num_iters = 10;  //max number of Newton iterations

  for(int i = 0; i < num_iters; i++) {
    if(fabs(f0) < eps*koef)
      break;
    shape_function_gradient_local(point, grad_local);
    f_der = LAMMPS_NS::vectorDot3D(grad_local, direction_vector_local);

    delta = -f0 / f_der;
    alpha += delta;
    MathExtraLiggghtsNonspherical::yabx3D(point, delta, direction_vector_local, point);
    if(fabs(delta) <= eps)
      break;
    f0 = shape_function_local(point);
  }
  result[0] = start_point[0] + alpha*direction_vector[0];
  result[1] = start_point[1] + alpha*direction_vector[1];
  result[2] = start_point[2] + alpha*direction_vector[2];
  return alpha;
}

void Superquadric::ellipsoid_line_intersection_local(const double *start_point_local, const double *direction_vector_local, double& alpha1, double& alpha2)
{
  const double x0 = start_point_local[0] * shape_inv[0];
  const double nx = direction_vector_local[0] * shape_inv[0];

  const double y0 = start_point_local[1] * shape_inv[1];
  const double ny = direction_vector_local[1] * shape_inv[1];

  const double z0 = start_point_local[2] * shape_inv[2];
  const double nz = direction_vector_local[2] * shape_inv[2];

  const double a = nx*nx + ny*ny + nz*nz;
  const double b = 2.0*(x0*nx + y0*ny + z0*nz);
  const double c = x0*x0 + y0*y0 + z0*z0 - 1.0;
  const double D = b*b - 4.0*a*c;
  const double two_a_inverted = 0.5 / a;
  if(D < 0) {
    alpha1 = -b * two_a_inverted;
    alpha2 = alpha1;
  } else {
    double sqrtD = sqrt(D);
    alpha1 = (-b + sqrtD) * two_a_inverted;
    alpha2 = (-b - sqrtD) * two_a_inverted;
  }
}

void Superquadric::ellipsoid_line_intersection(const double *start_point, const double *direction_vector, double& alpha1, double& alpha2)
{
  double start_point_local[3], direction_vector_local[3];
  global2local(start_point, start_point_local);
  rotate_global2local(direction_vector, direction_vector_local);

  ellipsoid_line_intersection_local(start_point_local, direction_vector_local, alpha1, alpha2);
}
#endif
