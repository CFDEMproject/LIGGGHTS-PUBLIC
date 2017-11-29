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

#ifndef LMP_MATH_EXTRA_LIGGGHTS_NONSPHERICAL_H
#define LMP_MATH_EXTRA_LIGGGHTS_NONSPHERICAL_H

#include "pointers.h"
#include <cmath>
#include <stdio.h>
#include <string.h>
#include "vector_liggghts.h"
#include "math_extra.h"
#include "random_park.h"
#include "ctype.h"
#include "contact_interface.h"
#include "math_extra_liggghts.h"

#define TOLERANCE_ORTHO 1e-10

using namespace LIGGGHTS::ContactModels;

namespace MathExtraLiggghtsNonspherical {

  const double scale_koef = 1.0001;

  struct point {
    double x,y,z;
    point():
      x(0), y(0), z(0) {};
    point(double *coords):
      x(coords[0]), y(coords[1]), z(coords[2]) {}
  };

  inline double sign(const double val);
  inline double pow_abs(double a, double b);
  inline double normN(const double *a, const int len);
  inline double normNsq(const double *a, const int len);
  inline void vectorCopyN(const double *from, const int len, double *to);
  inline double dotN(const double *a, const double *b, const int len);
  inline void addN(const double *a, const double *b, const int len, double *result);
  inline void subN(const double *a, const double *b, const int len, double *result);
  inline void mulN(const double *a, const double k, const int len, double *result);
  inline void mulN(const double k, const double *a, const int len, double *result);
  inline void matvecN(const double *A, const double *b, const int len, double *c);
  inline void matmat(const double *A, double *B, const int N, const int K, const int M, double *C);
  inline void zeros(double *a, const int len);
  inline void tensor_quat_rotate(const double *mat, const double *quat, double *result);
  template<int N, int m_start> void GMRES(const double *A, const double *b, const double *x0, double *x);
  inline double distsq(const double *v1, const double *v2, int len);
  inline double dist(const double *v1, const double *v2, int len);
  double point_wall_projection(const double *normal_vector, const double *point_on_wall, const double *input_point, double *output_point);
  inline void rotate_local2global(const double *quat, const double *input_coord, double *result);
  inline void rotate_global2local(const double *quat, const double *input_coord, double *result);
  double point_line_distance(double *pointA, double *pointB, double *input_point, double *output_point, double *delta);
  void integrate_quat(double *quat, double *omega, double dtq);
  void integrate_omega(double *quat, double *torque, double *inertia, double *omega, double dt);

  void angmom_to_omega(const double *quat, const double *angmom, const double *inertia, double *omega);
  void omega_to_angmom(const double *quat, const double *omega, const double *inertia, double *angmom);
  inline void rotate_local2global(const double *mat, const double *input_coord, double *result);
  inline void rotate_global2local(const double *mat, const double *input_coord, double *result);
  void no_squish_rotate(int k, double *p, double *q, double *inertia, double dt);
  void calc_conjqm(double *quat, double *angmom_body, double *conjqm);
  inline void invquat(double *quat);
  double clamp(double val, double min, double max);
  void line_segments_distance(double *P1, double *Q1, double *P2, double *Q2, double *C1, double *C2);
  double get_effective_radius(SurfacesIntersectData & sidata, double *blockiness_i, double *blockiness_j,double koefi, double koefj, double curvatureLimitFactor,LAMMPS_NS::Error *error);
  double get_effective_radius_wall(SurfacesIntersectData & sidata, double *blockiness_i, double koefi, double curvatureLimitFactor,LAMMPS_NS::Error *error);
  inline void yabx3D(double *a, double b, double *x, double *y);
  inline double pow_abs_int(double a, int b);
  inline int round_int(double x);
  inline void matvec(const double *m, const double *v, double *ans);
  inline void transpose_matvec(const double *m, const double *v, double *ans);
  inline void times3(const double *m, const double *m2, double *ans);
  inline void transpose_times3(const double *m, const double *m2, double *ans);
  inline void times3_transpose(const double *m, const double *m2, double *ans);
  inline void quat_to_mat(const double *quat, double *rotation_matrix);

  inline bool isInteger(double x);
  inline void surfacesIntersectNonSpherical(SurfacesIntersectData & sidata, double **x);
};

//rotate tensor by quaternion
inline void MathExtraLiggghtsNonspherical::tensor_quat_rotate(const double *tensor, const double *quat, double *result)
{
  double rotation_matrix[9], temp[9];
  MathExtraLiggghtsNonspherical::quat_to_mat(quat, rotation_matrix); //calc rotation matrix (A) from quaternion
  MathExtraLiggghtsNonspherical::times3(rotation_matrix, tensor, temp); // calc temp = A * M
  MathExtraLiggghtsNonspherical::times3_transpose(temp, rotation_matrix, result); //calc result = A * M * A^T
}

//calculation of 2-norm of the vector
inline double MathExtraLiggghtsNonspherical::normN(const double *a, const int len)
{
  double result = 0.0;
  for(int i = 0; i < len; i++)
    result += a[i]*a[i];
  return ::sqrt(result);
}

inline double MathExtraLiggghtsNonspherical::pow_abs_int(double a, int b)
{
  double result;
  a = fabs(a);
  switch(b) {
  case 0:
    result = 1.0;
    break;
  case 1:
    result = a;
    break;
  case 2:
    result = a*a;
    break;
  case 3:
    result = a*a*a;
    break;
  case 4: {
    double a2 = a*a;
    result = a2*a2;
    break;
  }
  case 5: {
    double a2 = a*a;
    result = a2*a2*a;
    break;
  }
  case 6: {
    double a3 = a*a*a;
    result = a3*a3;
    break;
  }
  case 7: {
    double a2 = a*a;
    double a3 = a2*a;
    result = a2*a2*a3;
    break;
  }
  case 8: {
    double a2 = a*a;
    double a4 = a2*a2;
    result = a4*a4;
    break;
  }
  case 9: {
    double a3 = a*a*a;
    result = a3*a3*a3;
    break;
  }
  default:
    result = a;
    for(int i = 0; i < b - 1; i++)
      result *= a;
    break;
  }
  return result;
}

inline double MathExtraLiggghtsNonspherical::pow_abs(double base, double exponent)
{
  base = fabs(base);
  double tol = 1e-12;
  if(fabs(exponent)<tol or fabs(base - 1.0) < tol)
    return 1.0;
  else if(fabs(base) < tol)
    return 0.0;
  else
  {
    double n = floor(fabs(exponent));
    int n_int = static_cast<int>(n);
    double result = pow_abs_int(base, n_int);
    double n_delta = fabs(exponent) - n;
    if(n_delta > 0.0) {
      if(fabs(2.0*n_delta - 1.0) < 2.0*tol)
        result *= ::sqrt(base);
      else
      {
       result *= ::pow(base, n_delta);
      }
    }
    if(exponent > 0.0)
      return result;
    else
      return 1.0 / result;
  }
}

//GMRES Linear System Solver
//https://math.berkeley.edu/~mgu/MA228A/saad-schultz.pdf
template<int N, int m_start> void MathExtraLiggghtsNonspherical::GMRES(const double *A, const double *b, const double *x0, double *x)
{
  int m = m_start;
  double res[N];
  double Ax0[N];
  matvecN(A, x0, N, Ax0); // calc A*x0
  subN(b, Ax0, N, res); //res = b - A*x0
  vectorCopyN(x0, N, x);
  double beta_sq = normNsq(res, N); //beta = ||res||
  if(beta_sq < 1e-20) {
    return;
  }

  double beta = ::sqrt(beta_sq);

  double y[m_start];
  double v[(m_start+1)*N];
  double h[(m_start+1)*m_start];
  double w[m_start*N];

  zeros(y, m);
  zeros(v, (m+1)*N);
  zeros(h, (m+1)*m);
  zeros(w, m*N);

  mulN(1.0/beta, res, N, v); //v0 = res / beta

  for(int j = 0; j < m; j++) {
    matvecN(A, v + j*N, N, w + j*N);
    for(int i = 0; i < j + 1; i++) {
      for(int ix = 0; ix < N; ix++) {
        h[i*m+j] += w[j*N + ix]*v[i*N + ix];
      }
      for(int ix = 0; ix < N; ix++)
        w[j*N + ix] -= h[i*m+j]*v[i*N + ix];
    }
    h[(j+1)*m + j] = normN(w + j*N, N);

    if(h[(j+1)*m + j] < 1e-12) {
      m = j + 1;
      break;
    } else
      mulN(1.0/h[(j+1)*m + j], w+j*N, N, v+(j+1)*N);
  }

  double hm[(m_start+1)*m_start];
  double g[m_start+1];
  zeros(g, m+1);
  g[0] = beta;

  for(int i = 0; i < m + 1; i++) {
    for(int j = 0; j < m; j++) {
      hm[i*m_start+j] = h[i*m_start+j];
    }
  }

  double row1[m_start];
  double row2[m_start];

  for(int j = 0; j < m; j++) {
    double root_inv = 1.0 / sqrt(hm[j*m_start+j]*hm[j*m_start+j] + hm[(j+1)*m_start+j]*hm[(j+1)*m_start+j]);
    double s =   hm[(j+1)*m_start+j] * root_inv;
    double c =       hm[j*m_start+j] * root_inv;

    for(int i = 0; i < m; i++) {
      row1[i] =  c*hm[j*m_start+i] + s*hm[(j+1)*m_start+i];
      row2[i] = -s*hm[j*m_start+i] + c*hm[(j+1)*m_start+i];
    }
    for(int i = 0; i < m; i++) {
      hm[j*m_start+i] = row1[i];
      hm[(j+1)*m_start+i] = row2[i];
    }
    double g1_temp =  c*g[j] + s*g[j+1];
    double g2_temp = -s*g[j] + c*g[j+1];

    g[j] = g1_temp;
    g[j + 1] = g2_temp;
  }

  for(int j = 0; j < m; j++) {
    int j1 = m - 1 - j;
    double summ = 0.0;
    for(int i = 0; i < j; i++) {
      int i1 = m - 1 - i;
      summ += hm[j1*m_start+i1] * y[i1];
    }
    y[j1] = (g[j1] - summ) / hm[j1*m_start + j1];
    for(int i = 0; i < N; i++)
      x[i] += y[j1] * v[j1*N+i];
  }
}

//signum function
inline double MathExtraLiggghtsNonspherical::sign(const double val)
{
  if(val > 1e-14)
    return 1.0;
  else if (val < -1e-14)
    return -1.0;
  else // -1e-14 <= val <= 1e-14
    return 0.0;
}

//calculation of 2-norm of the vector
inline double MathExtraLiggghtsNonspherical::normNsq(const double *a, const int len)
{
  double result = 0.0;
  for(int i = 0; i < len; i++)
    result += a[i]*a[i];
  return result;
}

//copy vector
inline void MathExtraLiggghtsNonspherical::vectorCopyN(const double *from, const int len, double *to)
{
  for(int i = 0; i < len; i++)
    to[i] = from[i];
}

//vector addition: result = a+b
inline void MathExtraLiggghtsNonspherical::addN(const double *a, const double *b, const int len, double *result)
{
  for(int i = 0; i < len; i++)
    result[i] = a[i] + b[i];
}

//vector substract: result = a-b
inline void MathExtraLiggghtsNonspherical::subN(const double *a, const double *b, const int len, double *result)
{
  for(int i = 0; i < len; i++)
    result[i] = a[i] - b[i];
}

//vector-scalar multiplication
inline void MathExtraLiggghtsNonspherical::mulN(const double *a, const double k, const int len, double *result)
{
  for(int i = 0; i < len; i++)
    result[i] = a[i]*k;
}

//vector-scalar multiplication
inline void MathExtraLiggghtsNonspherical::mulN(const double k, const double *a, const int len, double *result)
{
  for(int i = 0; i < len; i++)
    result[i] = a[i]*k;
}

//dot product of two vectors a and b of length N
inline double MathExtraLiggghtsNonspherical::dotN(const double *a, const double *b, const int len)
{
  double result = 0.0;
  for(int i = 0; i < len; i++)
    result += a[i]*b[i];
  return result;
}

//matrix-vector multiplication, size of A is NxN, size of b is N, (N = len)
inline void MathExtraLiggghtsNonspherical::matvecN(const double *A, const double *b, const int len, double *c)
{
  for(int i = 0; i < len; i++) {
    c[i] = 0.0;
    for(int j = 0; j < len; j++)
    {
      if(A[i*len+j] != 0.0 and b[j] != 0.0)
        c[i] += A[i*len+j]*b[j];
    }
  }
}
//matrix-matrix multiplication: C = A*B, A(NxK), B(KxM), C(NxM)
inline void MathExtraLiggghtsNonspherical::matmat(const double *A, double *B, const int N, const int K, const int M, double *C)
{
  for(int iM = 0; iM < M; iM++) {
    for(int iN = 0; iN < N; iN ++) {
      C[iN * M + iM] = 0.0;
      for(int iK = 0; iK < K; iK++) {
        if(A[iN * K + iK] != 0.0 and B[iK * M + iM] != 0.0)
          C[iN * M + iM] += A[iN * K + iK] * B[iK * M + iM];
      }
    }
  }
}

//fill the first <len> elements of a with zeros
inline void MathExtraLiggghtsNonspherical::zeros(double *a, const int len)
{
  memset(a, 0.0, len*sizeof(double));
}

//squared distance between two points in N-dimensional space
inline double MathExtraLiggghtsNonspherical::distsq(const double *v1, const double *v2, int len)
{
  double result = 0.0;
  for(int i = 0; i < len; i++)
    result += (v1[i]-v2[i])*(v1[i]-v2[i]);
  return result;
}

//rotate vector from local to global reference frame with a given quaternion
inline void MathExtraLiggghtsNonspherical::rotate_local2global(const double *quat, const double *input_coord, double *result) {
  double A[9];
  MathExtraLiggghtsNonspherical::quat_to_mat(quat, A);
  MathExtraLiggghtsNonspherical::matvec(A, input_coord, result);
}

//rotate vector from global to local reference frame with a given quaternion
inline void MathExtraLiggghtsNonspherical::rotate_global2local(const double *quat, const double *input_coord, double *result) {
  double A[9];
  MathExtraLiggghtsNonspherical::quat_to_mat(quat, A);
  MathExtraLiggghtsNonspherical::transpose_matvec(A, input_coord, result);
}

inline void MathExtraLiggghtsNonspherical::invquat(double *quat)
{
  double norm_sq_inv = 1.0 / LAMMPS_NS::vectorMag4D(quat);
  quat[0] *= norm_sq_inv;
  quat[1] *= -norm_sq_inv;
  quat[2] *= -norm_sq_inv;
  quat[3] *= -norm_sq_inv;
}

inline void MathExtraLiggghtsNonspherical::yabx3D(double *a, double b, double *x, double *y)
{
  y[0] = a[0] + b*x[0];
  y[1] = a[1] + b*x[1];
  y[2] = a[2] + b*x[2];
}

inline int MathExtraLiggghtsNonspherical::round_int(double x)
{
  return static_cast<int>(x+0.5);
}

inline void MathExtraLiggghtsNonspherical::matvec(const double *m, const double *v, double *ans)
{
  ans[0] = m[0]*v[0] + m[1]*v[1] + m[2]*v[2];
  ans[1] = m[3]*v[0] + m[4]*v[1] + m[5]*v[2];
  ans[2] = m[6]*v[0] + m[7]*v[1] + m[8]*v[2];
}

inline void MathExtraLiggghtsNonspherical::transpose_matvec(const double *m, const double *v, double *ans)
{
  ans[0] = m[0]*v[0] + m[3]*v[1] + m[6]*v[2];
  ans[1] = m[1]*v[0] + m[4]*v[1] + m[7]*v[2];
  ans[2] = m[2]*v[0] + m[5]*v[1] + m[8]*v[2];
}

inline void MathExtraLiggghtsNonspherical::times3(const double *m, const double *m2, double *ans)
{
  ans[0*3+0] = m[0*3+0]*m2[0*3+0] + m[0*3+1]*m2[1*3+0] + m[0*3+2]*m2[2*3+0];
  ans[0*3+1] = m[0*3+0]*m2[0*3+1] + m[0*3+1]*m2[1*3+1] + m[0*3+2]*m2[2*3+1];
  ans[0*3+2] = m[0*3+0]*m2[0*3+2] + m[0*3+1]*m2[1*3+2] + m[0*3+2]*m2[2*3+2];
  ans[1*3+0] = m[1*3+0]*m2[0*3+0] + m[1*3+1]*m2[1*3+0] + m[1*3+2]*m2[2*3+0];
  ans[1*3+1] = m[1*3+0]*m2[0*3+1] + m[1*3+1]*m2[1*3+1] + m[1*3+2]*m2[2*3+1];
  ans[1*3+2] = m[1*3+0]*m2[0*3+2] + m[1*3+1]*m2[1*3+2] + m[1*3+2]*m2[2*3+2];
  ans[2*3+0] = m[2*3+0]*m2[0*3+0] + m[2*3+1]*m2[1*3+0] + m[2*3+2]*m2[2*3+0];
  ans[2*3+1] = m[2*3+0]*m2[0*3+1] + m[2*3+1]*m2[1*3+1] + m[2*3+2]*m2[2*3+1];
  ans[2*3+2] = m[2*3+0]*m2[0*3+2] + m[2*3+1]*m2[1*3+2] + m[2*3+2]*m2[2*3+2];
}

inline void MathExtraLiggghtsNonspherical::transpose_times3(const double *m, const double *m2, double *ans)
{
  ans[0*3+0] = m[0*3+0]*m2[0*3+0] + m[1*3+0]*m2[1*3+0] + m[2*3+0]*m2[2*3+0];
  ans[0*3+1] = m[0*3+0]*m2[0*3+1] + m[1*3+0]*m2[1*3+1] + m[2*3+0]*m2[2*3+1];
  ans[0*3+2] = m[0*3+0]*m2[0*3+2] + m[1*3+0]*m2[1*3+2] + m[2*3+0]*m2[2*3+2];
  ans[1*3+0] = m[0*3+1]*m2[0*3+0] + m[1*3+1]*m2[1*3+0] + m[2*3+1]*m2[2*3+0];
  ans[1*3+1] = m[0*3+1]*m2[0*3+1] + m[1*3+1]*m2[1*3+1] + m[2*3+1]*m2[2*3+1];
  ans[1*3+2] = m[0*3+1]*m2[0*3+2] + m[1*3+1]*m2[1*3+2] + m[2*3+1]*m2[2*3+2];
  ans[2*3+0] = m[0*3+2]*m2[0*3+0] + m[1*3+2]*m2[1*3+0] + m[2*3+2]*m2[2*3+0];
  ans[2*3+1] = m[0*3+2]*m2[0*3+1] + m[1*3+2]*m2[1*3+1] + m[2*3+2]*m2[2*3+1];
  ans[2*3+2] = m[0*3+2]*m2[0*3+2] + m[1*3+2]*m2[1*3+2] + m[2*3+2]*m2[2*3+2];
}

inline void MathExtraLiggghtsNonspherical::times3_transpose(const double *m, const double *m2, double *ans)
{
  ans[0*3+0] = m[0*3+0]*m2[0*3+0] + m[0*3+1]*m2[0*3+1] + m[0*3+2]*m2[0*3+2];
  ans[0*3+1] = m[0*3+0]*m2[1*3+0] + m[0*3+1]*m2[1*3+1] + m[0*3+2]*m2[1*3+2];
  ans[0*3+2] = m[0*3+0]*m2[2*3+0] + m[0*3+1]*m2[2*3+1] + m[0*3+2]*m2[2*3+2];
  ans[1*3+0] = m[1*3+0]*m2[0*3+0] + m[1*3+1]*m2[0*3+1] + m[1*3+2]*m2[0*3+2];
  ans[1*3+1] = m[1*3+0]*m2[1*3+0] + m[1*3+1]*m2[1*3+1] + m[1*3+2]*m2[1*3+2];
  ans[1*3+2] = m[1*3+0]*m2[2*3+0] + m[1*3+1]*m2[2*3+1] + m[1*3+2]*m2[2*3+2];
  ans[2*3+0] = m[2*3+0]*m2[0*3+0] + m[2*3+1]*m2[0*3+1] + m[2*3+2]*m2[0*3+2];
  ans[2*3+1] = m[2*3+0]*m2[1*3+0] + m[2*3+1]*m2[1*3+1] + m[2*3+2]*m2[1*3+2];
  ans[2*3+2] = m[2*3+0]*m2[2*3+0] + m[2*3+1]*m2[2*3+1] + m[2*3+2]*m2[2*3+2];
}

inline void MathExtraLiggghtsNonspherical::quat_to_mat(const double *quat, double *rotation_matrix) {
  double w2 = quat[0]*quat[0];
  double i2 = quat[1]*quat[1];
  double j2 = quat[2]*quat[2];
  double k2 = quat[3]*quat[3];
  double twoij = 2.0*quat[1]*quat[2];
  double twoik = 2.0*quat[1]*quat[3];
  double twojk = 2.0*quat[2]*quat[3];
  double twoiw = 2.0*quat[1]*quat[0];
  double twojw = 2.0*quat[2]*quat[0];
  double twokw = 2.0*quat[3]*quat[0];

  rotation_matrix[0] = w2+i2-j2-k2;
  rotation_matrix[1] = twoij-twokw;
  rotation_matrix[2] = twojw+twoik;

  rotation_matrix[3] = twoij+twokw;
  rotation_matrix[4] = w2-i2+j2-k2;
  rotation_matrix[5] = twojk-twoiw;

  rotation_matrix[6] = twoik-twojw;
  rotation_matrix[7] = twojk+twoiw;
  rotation_matrix[8] = w2-i2-j2+k2;
}

inline void MathExtraLiggghtsNonspherical::surfacesIntersectNonSpherical(SurfacesIntersectData & sidata, double **x)
{
  #ifdef NONSPHERICAL_ACTIVE_FLAG
  double xci[3], xcj[3];
  double v_rot_i[3], v_rot_j[3];
  double omega_relative[3], v_relative[3], v_rot_relative[3];
  if(sidata.is_wall) {
    LAMMPS_NS::vectorSubtract3D(sidata.contact_point, x[sidata.i], xci);
    LAMMPS_NS::vectorCross3D(sidata.omega_i, xci, v_rot_i);
    LAMMPS_NS::vectorCopy3D(v_rot_i, v_rot_relative);
    LAMMPS_NS::vectorCopy3D(sidata.omega_i, omega_relative);
    sidata.cri = LAMMPS_NS::vectorMag3D(xci);
  } else {
    LAMMPS_NS::vectorSubtract3D(sidata.contact_point, x[sidata.i], xci);
    LAMMPS_NS::vectorSubtract3D(sidata.contact_point, x[sidata.j], xcj);
    LAMMPS_NS::vectorCross3D(sidata.omega_i, xci, v_rot_i);
    LAMMPS_NS::vectorCross3D(sidata.omega_j, xcj, v_rot_j);
    LAMMPS_NS::vectorSubtract3D(v_rot_i, v_rot_j, v_rot_relative);
    LAMMPS_NS::vectorSubtract3D(sidata.omega_i, sidata.omega_j, omega_relative);
    sidata.cri = LAMMPS_NS::vectorMag3D(xci);
    sidata.crj = LAMMPS_NS::vectorMag3D(xcj);
  }
      // relative velocity
  v_relative[0] = sidata.v_i[0] - sidata.v_j[0] + v_rot_relative[0];
  v_relative[1] = sidata.v_i[1] - sidata.v_j[1] + v_rot_relative[1];
  v_relative[2] = sidata.v_i[2] - sidata.v_j[2] + v_rot_relative[2];

  // normal component
  const double vn = LAMMPS_NS::vectorDot3D(v_relative, sidata.en);
  // tangential components
  sidata.vtr1 = v_relative[0] - vn * sidata.en[0];
  sidata.vtr2 = v_relative[1] - vn * sidata.en[1];
  sidata.vtr3 = v_relative[2] - vn * sidata.en[2];

  sidata.wr1 = omega_relative[0];
  sidata.wr2 = omega_relative[1];
  sidata.wr3 = omega_relative[2];

  sidata.vn = vn;
  #endif
}

inline bool MathExtraLiggghtsNonspherical::isInteger(double x)
{
  double x_temp = floor(x + 1e-2);
  return MathExtraLiggghts::compDouble(x, x_temp, 1e-2);
}
#endif
