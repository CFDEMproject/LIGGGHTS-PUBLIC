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
    Alexander Podlozhnyuk (DCS Computing GmbH, Linz)

    Copyright 2015-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#include "math_extra_liggghts_nonspherical.h"

#include <map>
#include "math_const.h"

namespace MathExtraLiggghtsNonspherical {

//project point onto the surface defined by normal vector and some point on the surface
double point_wall_projection(const double *normal_vector, const double *point_on_wall, const double *input_point, double *output_point) {
  const double x0 = point_on_wall[0];
  const double y0 = point_on_wall[1];
  const double z0 = point_on_wall[2];

  const double x = input_point[0];
  const double y = input_point[1];
  const double z = input_point[2];

  const double A = normal_vector[0];
  const double B = normal_vector[1];
  const double C = normal_vector[2];

  const double alpha = -(A*(x - x0) + B*(y - y0) + C*(z - z0)) / (A*A + B*B + C*C);
  output_point[0] = x + alpha * A;
  output_point[1] = y + alpha * B;
  output_point[2] = z + alpha * C;
  return fabs((A*(x - x0) + B*(y - y0) + C*(z - z0)) / sqrt(A*A + B*B + C*C));
}

//calculate right hand side term in dynamic euler equation: dw/dt = f(w), w - angular velocity
void dynamic_euler_right_term(const double *omega, const double *torque, const double *inertia, double *result)
{
  const double wx = omega[0];
  const double wy = omega[1];
  const double wz = omega[2];

  const double Ix = inertia[0];
  const double Iy = inertia[1];
  const double Iz = inertia[2];

  result[0] = (wy*wz*(Iy - Iz) + torque[0]) / Ix;
  result[1] = (wz*wx*(Iz - Ix) + torque[1]) / Iy;
  result[2] = (wx*wy*(Ix - Iy) + torque[2]) / Iz;
}

//integrate angular velocity
void integrate_omega(double *quat, double *torque, double *inertia, double *omega, double dt)
{
  double torque_body[3], omega_body[3];

  MathExtraLiggghtsNonspherical::rotate_global2local(quat, omega, omega_body);
  MathExtraLiggghtsNonspherical::rotate_global2local(quat, torque, torque_body);

  double f[3];
  dynamic_euler_right_term(omega_body, torque_body, inertia, f);
  for(int i = 0; i < 3; i++)
    omega_body[i] += dt*f[i];
  MathExtraLiggghtsNonspherical::rotate_local2global(quat, omega_body, omega);
}

//distance between point and line
double point_line_distance(double *pointA, double *pointB, double *inputPoint, double *closestPoint, double *delta)
{
  const double x0 = inputPoint[0];
  const double y0 = inputPoint[1];
  const double z0 = inputPoint[2];

  const double nx = pointB[0] - pointA[0];
  const double ny = pointB[1] - pointA[1];
  const double nz = pointB[2] - pointA[2];

  const double xA = pointA[0];
  const double yA = pointA[1];
  const double zA = pointA[2];

  double alpha = ((x0 - xA)*nx + (y0 - yA)*ny + (z0 - zA)*nz) / (nx*nx + ny*ny + nz*nz);
  if(alpha > 1.0)
    alpha = 1.0;
  if(alpha < 0.0)
    alpha = 0.0;
  closestPoint[0] = xA + alpha*nx;
  closestPoint[1] = yA + alpha*ny;
  closestPoint[2] = zA + alpha*nz;

  LAMMPS_NS::vectorSubtract3D(closestPoint,inputPoint,delta);
  return LAMMPS_NS::vectorMag3D(delta);
}

// integrate quaternion
void integrate_quat(double *quat, double *omega, double dt)
{
  double omega_mag = LAMMPS_NS::vectorMag3D(omega);
  double deltaq[4];
  if(omega_mag > 1e-8) {
    deltaq[0] = cos(0.5*omega_mag*dt);
    deltaq[1] = sin(0.5*omega_mag*dt) * omega[0]/omega_mag;
    deltaq[2] = sin(0.5*omega_mag*dt) * omega[1]/omega_mag,
    deltaq[3] = sin(0.5*omega_mag*dt) * omega[2]/omega_mag;
  } else {
    deltaq[0] = 1.0;
    deltaq[1] = 0.5 * dt * omega[0];
    deltaq[2] = 0.5 * dt * omega[1];
    deltaq[3] = 0.5 * dt * omega[2];
  }
  double quat_temp[4];
  MathExtra::quatquat(deltaq, quat, quat_temp);
  MathExtra::qnormalize(quat_temp);
  LAMMPS_NS::vectorCopy4D(quat_temp, quat);

  /*double wq[4];
  MathExtra::vecquat(omega, quat, wq);
  for(int i = 0; i < 4; i++)
    quat[i] += 0.5*dt * wq[i];
  MathExtra::qnormalize(quat);*/
}

//calculate angular velocity from angular moment
void angmom_to_omega(const double *quat, const double *angmom, const double *inertia, double *omega)
{
  double omega_body[3], angmom_body[3];
  rotate_global2local(quat, angmom, angmom_body);
  for(int k = 0; k < 3; k++)
    omega_body[k] = angmom_body[k] / inertia[k];
  rotate_local2global(quat, omega_body, omega);
}

//calculate angular moment from angular velocity
void omega_to_angmom(const double *quat, const double *omega, const double *inertia, double *angmom)
{
  double omega_body[3], angmom_body[3];
  rotate_global2local(quat, omega, omega_body);
  for(int k = 0; k < 3; k++)
    angmom_body[k] = omega_body[k] * inertia[k];
  rotate_local2global(quat, angmom_body, angmom);
}

void no_squish_rotate(int k, double *p, double *q, double *inertia, double dt)
{
  double phi,c_phi,s_phi,kp[4]={},kq[4]={};

  // apply permuation operator on p and q, get kp and kq

  if (k == 1) {
    kq[0] = -q[1];  kp[0] = -p[1];
    kq[1] =  q[0];  kp[1] =  p[0];
    kq[2] =  q[3];  kp[2] =  p[3];
    kq[3] = -q[2];  kp[3] = -p[2];
  } else if (k == 2) {
    kq[0] = -q[2];  kp[0] = -p[2];
    kq[1] = -q[3];  kp[1] = -p[3];
    kq[2] =  q[0];  kp[2] =  p[0];
    kq[3] =  q[1];  kp[3] =  p[1];
  } else if (k == 3) {
    kq[0] = -q[3];  kp[0] = -p[3];
    kq[1] =  q[2];  kp[1] =  p[2];
    kq[2] = -q[1];  kp[2] = -p[1];
    kq[3] =  q[0];  kp[3] =  p[0];
  }

  // obtain phi, cosines and sines

  phi = p[0]*kq[0] + p[1]*kq[1] + p[2]*kq[2] + p[3]*kq[3];
  phi /= 4.0 * inertia[k-1];
  c_phi = cos(dt * phi);
  s_phi = sin(dt * phi);

  // advance p and q

  p[0] = c_phi*p[0] + s_phi*kp[0];
  p[1] = c_phi*p[1] + s_phi*kp[1];
  p[2] = c_phi*p[2] + s_phi*kp[2];
  p[3] = c_phi*p[3] + s_phi*kp[3];

  q[0] = c_phi*q[0] + s_phi*kq[0];
  q[1] = c_phi*q[1] + s_phi*kq[1];
  q[2] = c_phi*q[2] + s_phi*kq[2];
  q[3] = c_phi*q[3] + s_phi*kq[3];
}

void calc_conjqm(double *quat, double *angmom_body, double *conjqm)
{
  MathExtra::quatvec(quat,angmom_body,conjqm);
  conjqm[0] *= 2.0;
  conjqm[1] *= 2.0;
  conjqm[2] *= 2.0;
  conjqm[3] *= 2.0;
}

double clamp(double val, double min, double max)
{
  if(val < min) return min;
  if(val > max) return max;
  return val;
}

void line_segments_distance(double *P1, double *Q1, double *P2, double *Q2, double *C1, double *C2)
{
  double tol = 1e-14;
  double d1[3], d2[3], r[3];
  LAMMPS_NS::vectorSubtract3D(Q1, P1, d1);
  LAMMPS_NS::vectorSubtract3D(Q2, P2, d2);
  LAMMPS_NS::vectorSubtract3D(P1, P2, r);

  const double a = LAMMPS_NS::vectorDot3D(d1, d1); //always >= 0
  const double e = LAMMPS_NS::vectorDot3D(d2, d2); //always >= 0
  const double f = LAMMPS_NS::vectorDot3D(d2, r);

  double s, t;

  if(a <= tol && e <= tol) { //line segments with zero length (points)
    s = t = 0.0;
    LAMMPS_NS::vectorCopy3D(P1, C1);
    LAMMPS_NS::vectorCopy3D(P2, C2);
    return;
  }
  if(a <= tol) {
    s = 0.0;
    t = clamp(f / e, 0.0, 1.0);
  } else {
    const double c = LAMMPS_NS::vectorDot3D(d1, r);
    if(e <= tol) {
      t = 0.0;
      s = clamp(-c / a, 0.0, 1.0);
    } else {
      const double b = LAMMPS_NS::vectorDot3D(d1, d2);
      const double determinant = a*e - b*b; //always non-negative
      if(determinant >= tol*tol) { //segments non parallel
        s = clamp( (b*f - c*e)/determinant, 0.0, 1.0);
        t = (b*s + f) / e;
        if(t < 0.0) {
          t = 0.0;
          s = clamp(-c / a, 0.0, 1.0);
        } else if(t > 1.0) {
          t = 1.0;
          s = clamp((b-c) / a, 0.0, 1.0);
        }
      } else { //segments are parallel
        double s1 = 0.0;
        double s2 = 1.0;
        double t1 = (f + s1*b) / e;
        double t2 = (f + s2*b) / e;
        if(t1 < 0.0) {
          t1 = 0.0;
          s1 = clamp(-c / a, 0.0, 1.0);
        } else if(t1 > 1.0) {
          t1 = 1.0;
          s1 = clamp((b-c) / a, 0.0, 1.0);
        }
        if(t2 < 0.0) {
          t2 = 0.0;
          s2 = clamp(-c / a, 0.0, 1.0);
        } else if(t2 > 1.0) {
          t2 = 1.0;
          s2 = clamp((b-c) / a, 0.0, 1.0);
        }
        s = 0.5*(s1 + s2);
        t = 0.5*(t1 + t2);
      }
    }
  }
  for(int i = 0; i < 3; i++) {
    C1[i] = P1[i] + d1[i]*s;
    C2[i] = P2[i] + d2[i]*t;
  }
}

//effective curvature radius of a pair of particles at the contact point
double get_effective_radius(SurfacesIntersectData & sidata, double *blockiness_i, double *blockiness_j, double koefi, double koefj, double curvatureLimitFactor, LAMMPS_NS::Error *error)
{
  #ifdef SUPERQUADRIC_ACTIVE_FLAG
  if(blockiness_i == NULL)
    error->one(FLERR,"blockiness_i array in get_effective_radius() is NULL");
  if(blockiness_j == NULL)
    error->one(FLERR,"blockiness_j array in get_effective_radius() is NULL");

  double keff = koefi + koefj;
  //if(blockiness_i[0] == 2.0 and blockiness_i[1] == 2.0 and blockiness_j[0] == 2.0 and blockiness_j[1] == 2.0)
  //  return 1.0 / keff;
  double rmax = sidata.reff * curvatureLimitFactor;
  return (keff < 1.0/rmax)? rmax: 1.0/keff;
  #else
  return 0.;
  #endif
}

double get_effective_radius_wall(SurfacesIntersectData & sidata, double *blockiness_i, double koefi, double curvatureLimitFactor, LAMMPS_NS::Error *error)
{
  #ifdef SUPERQUADRIC_ACTIVE_FLAG

  double keff = koefi;
  if(blockiness_i == NULL)
    error->one(FLERR,"blockiness_i array in get_effective_radius_wall() is NULL");
  //if(blockiness_i[0] == 2.0 and blockiness_i[1] == 2.0)
  //  return 1.0 / keff;
  double rmax = sidata.radi * curvatureLimitFactor;
  return (keff < 1.0/rmax)? rmax: 1.0/keff;
  #else
  return 0.;
  #endif
}

//distance between two points in N-dimensional space
inline double dist(const double *v1, const double *v2, int len)
{
  return sqrt(distsq(v1, v2, len));
}

}
