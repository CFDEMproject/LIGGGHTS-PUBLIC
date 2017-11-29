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

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef LMP_MATH_EXTRA_LIGGGHTS_H
#define LMP_MATH_EXTRA_LIGGGHTS_H

#include "pointers.h"
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include "error.h"
#include "vector_liggghts.h"
#include "math_extra.h"
#include "random_park.h"
#include "ctype.h"

#define TOLERANCE_ORTHO 1e-10

namespace MathExtraLiggghts {

  inline void col_times3(const double m[3][3],const double *v, double *ans);

  inline double mdet(const double m[3][3],LAMMPS_NS::Error *error);

  //cubic root approx
  inline double cbrt_5d(double d);
  inline double cbrta_halleyd(const double a, const double R);
  inline double halley_cbrt1d(double d);

  //exp aproximation
  inline double exp_fast(double x);

  inline double min(const double a, const double b, const double c);
  inline double max(const double a, const double b, const double c);
  inline double min(const double a, const double b, const double c, const double d);
  inline double max(const double a, const double b, const double c, const double d);
  template <typename T>
  inline T min(const T * const input, const unsigned int n, int &which);
  template <typename T>
  inline T max(const T * const input, const unsigned int n, int &which);
  template <typename T>
  inline T abs(const T a);

  inline void matrix_invert_4x4_special(double matrix[4][4]);
  inline void transpose3(const double m[3][3], double ans[3][3]);
  inline int is_inside_tet(double *pos,double invmatrix[4][4]);

  inline void local_coosys_to_cartesian(double * const global, const double * const local, const double * const ex_local, const double * const ey_local, const double * const ez_local);
  inline void cartesian_coosys_to_local(double *local,double *global, double *ex_local, double *ey_local, double *ez_local,LAMMPS_NS::Error *error);
  inline void cartesian_coosys_to_local_orthogonal(double *local,double *global, double *ex_local, double *ey_local, double *ez_local,LAMMPS_NS::Error *error);

  // quaternion operations
  inline bool is_unit_quat(const double *q);
  inline void quat_normalize(double *q);
  inline void qconjugate(const double * const q, double * const qc);
  inline void quat_from_vec(const double *v, double *q);
  inline void vec_from_quat(const double *q, double * const v);
  inline void vec_quat_rotate(const double * const vec, const double * const quat, double *result);
  inline void vec_quat_rotate(double * const vec, const double * const quat);
  inline void vec_quat_rotate(int * const vec, const double * const quat) { UNUSED(vec); UNUSED(quat); }
  inline void vec_quat_rotate(bool * const vec, const double * const quat) { UNUSED(vec); UNUSED(quat); }
  inline void quat_diff(double *q_new, double *q_old, double *q_diff);
  inline void angmom_from_omega(double *w,
                                  double *ex, double *ey, double *ez,
                                  double *idiag, double *m);

  // double comparison
  inline bool compDouble(double const a, double const b, double const prec = 1e-13);

  // calculate barycentrc coordinates of p w.r.t node
  inline void calcBaryTriCoords(double *p, double **edgeVec, double *edgeLen, double *bary);
  inline void calcBaryTriCoords(double *p, double *edgeVec0, double *edgeVec1, double *edgeVec2, double *edgeLen, double *bary);

  inline void random_unit_quat(LAMMPS_NS::RanPark *random,double *quat);

  inline bool is_int(char *str);
  inline void generateComplementBasis(double *uVec, double *vVec, double *direction);

  // template signum function
  template <typename T>
  int sgn(T val) {
      return (T(0) < val) - (val < T(0));
  }

  // prime number test
  inline bool isPrime(int val);
};

/* ----------------------------------------------------------------------
   matrix  times col vector
------------------------------------------------------------------------- */

void MathExtraLiggghts::col_times3(const double m[3][3],const double *v, double *ans)
{
  ans[0] = m[0][0]*v[0]+v[1]*m[1][0]+v[2]*m[2][0];
  ans[1] = v[0]*m[0][1]+m[1][1]*v[1]+v[2]*m[2][1];
  ans[2] = v[0]*m[0][2]+v[1]*m[1][2]+m[2][2]*v[2];
}

/* ----------------------------------------------------------------------
Matrix determinant
------------------------------------------------------------------------- */

double MathExtraLiggghts::mdet(const double m[3][3],LAMMPS_NS::Error *error)
{
    UNUSED(error);
    return ( -m[0][2]*m[1][1]*m[2][0] + m[0][1]*m[1][2]*m[2][0] + m[0][2]*m[1][0]*m[2][1] - m[0][0]*m[1][2]*m[2][1] - m[0][1]*m[1][0]*m[2][2] + m[0][0]*m[1][1]*m[2][2] );

}

/* ----------------------------------------------------------------------
   Cubic root approx.
------------------------------------------------------------------------- */

inline double MathExtraLiggghts::cbrt_5d(double d)
{
   const unsigned int B1 = 715094163;
   double t = 0.0;
   unsigned int* pt = (unsigned int*) &t;
   unsigned int* px = (unsigned int*) &d;
   pt[1]=px[1]/3+B1;
   return t;
}

inline double MathExtraLiggghts::cbrta_halleyd(const double a, const double R)
{
    const double a3 = a*a*a;
    const double b= a * (a3 + R + R) / (a3 + a3 + R);
    return b;
}

// cube root approximation using 1 iteration of Halley's method (double)
inline double MathExtraLiggghts::halley_cbrt1d(double d)
{
    double a = cbrt_5d(d);
    return cbrta_halleyd(a, d);
}

/* ----------------------------------------------------------------------
   exp approx
------------------------------------------------------------------------- */

inline double MathExtraLiggghts::exp_fast(double x)
{
      x = 1.0 + x / 256.0;
      x *= x; x *= x; x *= x; x *= x;
      x *= x; x *= x; x *= x; x *= x;
      return x;
}

/* ----------------------------------------------------------------------
   min max stuff
------------------------------------------------------------------------- */

double MathExtraLiggghts::min(const double a, const double b, const double c)
{
    return std::min(std::min(a,b),c);
}

double MathExtraLiggghts::max(const double a, const double b, const double c)
{
    return std::max(std::max(a,b),c);
}

double MathExtraLiggghts::min(const double a, const double b, const double c, const double d)
{
    return std::min(std::min(a,b),std::min(c,d));
}

double MathExtraLiggghts::max(const double a, const double b, const double c, const double d)
{
    return std::max(std::max(a,b),std::max(c,d));
}

template <typename T>
T MathExtraLiggghts::min(const T * const input, const unsigned int n, int &which)
{
    T min = input[0];
    which = 0;

    for(unsigned int i = 1; i < n; i++)
    {
        if(input[i] < min)
        {
            which = i;
            min = input[i];
        }
    }
    return min;
}

template <typename T>
T MathExtraLiggghts::max(const T * const input, const unsigned int n, int &which)
{
    T max = input[0];
    which = 0;

    for(unsigned int i = 1; i < n; i++)
    {
        if(input[i] > max)
        {
            which = i;
            max = input[i];
        }
    }
    return max;
}

template <typename T>
T MathExtraLiggghts::abs(const T a) { return a < static_cast<T>(0) ? -a : a; }

/*----------------------------------------------------------------------
   inverts a special 4x4 matrix that looks like this

   m11 m12 m13 m14
   m21 m22 m23 m24
   m31 m32 m33 m34
   1   1   1   1
------------------------------------------------------------------------- */
inline void MathExtraLiggghts::matrix_invert_4x4_special(double matrix[4][4])
{
    double fac,invfac,v1x,v2x,v3x,v4x,v1y,v2y,v3y,v4y,v1z,v2z,v3z,v4z;
    v1x = matrix[0][0]; v1y = matrix[1][0]; v1z = matrix[2][0];
    v2x = matrix[0][1]; v2y = matrix[1][1]; v2z = matrix[2][1];
    v3x = matrix[0][2]; v3y = matrix[1][2]; v3z = matrix[2][2];
    v4x = matrix[0][3]; v4y = matrix[1][3]; v4z = matrix[2][3];

    fac = -v1x*v2z*v3y+v1x*v2y*v3z+v2z*v3y*v4x-v2y*v3z*v4x+v1x*v2z*v4y-
        v2z*v3x*v4y-v1x*v3z*v4y+v2x*v3z*v4y+v1z*
        (v2x*v3y-v3y*v4x+v2y*(-v3x+v4x)-v2x*v4y+v3x*v4y)-
        v1x*v2y*v4z+v2y*v3x*v4z+v1x*v3y*v4z-v2x*v3y*v4z+v1y*
        (v2z*v3x-v2x*v3z-v2z*v4x+v3z*v4x+v2x*v4z-v3x*v4z);

    invfac = 1./fac;

    matrix[0][0] = (-v3z*v4y+v2z*(-v3y+v4y)+v2y*(v3z-v4z)+v3y*v4z)*invfac;
    matrix[1][0] = (v1z*(v3y-v4y)+v3z*v4y-v3y*v4z+v1y*(-v3z+v4z))*invfac;
    matrix[2][0] = (-v2z*v4y+v1z*(-v2y+v4y)+v1y*(v2z-v4z)+v2y*v4z)*invfac;
    matrix[3][0] = (v1z*(v2y-v3y)+v2z*v3y-v2y*v3z+v1y*(-v2z+v3z))*invfac;

    matrix[0][1] = (v2z*(v3x-v4x)+v3z*v4x-v3x*v4z+v2x*(-v3z+v4z))*invfac;
    matrix[1][1] = (-v3z*v4x+v1z*(-v3x+v4x)+v1x*(v3z-v4z)+v3x*v4z)*invfac;
    matrix[2][1] = (v1z*(v2x-v4x)+v2z*v4x-v2x*v4z+v1x*(-v2z+v4z))*invfac;
    matrix[3][1] = (-v2z*v3x+v1z*(-v2x+v3x)+v1x*(v2z-v3z)+v2x*v3z)*invfac;

    matrix[0][2] = (-v3y*v4x+v2y*(-v3x+v4x)+v2x*(v3y-v4y)+v3x*v4y)*invfac;
    matrix[1][2] = (v1y*(v3x-v4x)+v3y*v4x-v3x*v4y+v1x*(-v3y+v4y))*invfac;
    matrix[2][2] = (-v2y*v4x+v1y*(-v2x+v4x)+v1x*(v2y-v4y)+v2x*v4y)*invfac;
    matrix[3][2] = (v1y*(v2x-v3x)+v2y*v3x-v2x*v3y+v1x*(-v2y+v3y))*invfac;

    matrix[0][2] = (v2z*v3y*v4x-v2y*v3z*v4x-v2z*v3x*v4y+v2x*v3z*v4y+v2y*v3x*v4z-v2x*v3y*v4z)*invfac;
    matrix[1][2] = (-v1z*v3y*v4x+v1y*v3z*v4x+v1z*v3x*v4y-v1x*v3z*v4y-v1y*v3x*v4z+v1x*v3y*v4z)*invfac;
    matrix[2][2] = (v1z*v2y*v4x-v1y*v2z*v4x-v1z*v2x*v4y+v1x*v2z*v4y+v1y*v2x*v4z-v1x*v2y*v4z)*invfac;
    matrix[3][2] = (-v1z*v2y*v3x+v1y*v2z*v3x+v1z*v2x*v3y-v1x*v2z*v3y-v1y*v2x*v3z+v1x*v2y*v3z)*invfac;
}

/* ----------------------------------------------------------------------
   transpose mat1
------------------------------------------------------------------------- */

inline void MathExtraLiggghts::transpose3(const double m[3][3], double ans[3][3])
{
  ans[0][0] = m[0][0];
  ans[0][1] = m[1][0];
  ans[0][2] = m[2][0];

  ans[1][0] = m[0][1];
  ans[1][1] = m[1][1];
  ans[1][2] = m[2][1];

  ans[2][0] = m[0][2];
  ans[2][1] = m[1][2];
  ans[2][2] = m[2][2];
}

/*----------------------------------------------------------------------
   checks if a point is inside a tetrader, described by an inverse matrix
   This inverse matrix is the the inverse of a special 4x4 matrix that looks like this

   m11 m12 m13 m14
   m21 m22 m23 m24
   m31 m32 m33 m34
   1   1   1   1

   where m11,m21,m31 is vertex 1 etc.
------------------------------------------------------------------------- */

inline int MathExtraLiggghts::is_inside_tet(double *pos,double invmatrix[4][4])
{
    double result[4];

    result[0] = invmatrix[0][0] * pos[0] + invmatrix[0][1] * pos[1] + invmatrix[0][2] * pos[2] + invmatrix[0][3];
    result[1] = invmatrix[1][0] * pos[0] + invmatrix[1][1] * pos[1] + invmatrix[1][2] * pos[2] + invmatrix[1][3];
    result[2] = invmatrix[2][0] * pos[0] + invmatrix[2][1] * pos[1] + invmatrix[2][2] * pos[2] + invmatrix[2][3];
    result[3] = invmatrix[3][0] * pos[0] + invmatrix[3][1] * pos[1] + invmatrix[3][2] * pos[2] + invmatrix[3][3];

    if(max(result[0],result[1],result[2],result[3]) > 1.0) return 0;
    return 1;
}

/*----------------------------------------------------------------------
   transform from local to global coords
------------------------------------------------------------------------- */

void MathExtraLiggghts::local_coosys_to_cartesian(double * const global, const double * const local, const double * const ex_local, const double * const ey_local, const double * const ez_local)
{
    global[0] = local[0]*ex_local[0] + local[1]*ey_local[0] + local[2]*ez_local[0];
    global[1] = local[0]*ex_local[1] + local[1]*ey_local[1] + local[2]*ez_local[1];
    global[2] = local[0]*ex_local[2] + local[1]*ey_local[2] + local[2]*ez_local[2];
}

/*----------------------------------------------------------------------
   transform from global to local coords
------------------------------------------------------------------------- */

void MathExtraLiggghts::cartesian_coosys_to_local(double *local,double *global, double *ex_local, double *ey_local, double *ez_local,LAMMPS_NS::Error *error)
{
  UNUSED(error);
  double M[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
  double Mt[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};

  // set up the matrix
  LAMMPS_NS::vectorCopy3D(ex_local,M[0]);
  LAMMPS_NS::vectorCopy3D(ey_local,M[1]);
  LAMMPS_NS::vectorCopy3D(ez_local,M[2]);
  MathExtraLiggghts::transpose3(M,Mt);

  // solve
  MathExtra::mldivide3(Mt,global,local);
}

/*----------------------------------------------------------------------
   transform from global to local coords
   faster for orthogonal matrix
------------------------------------------------------------------------- */

void MathExtraLiggghts::cartesian_coosys_to_local_orthogonal(double *local,double *global, double *ex_local, double *ey_local, double *ez_local,LAMMPS_NS::Error *error)
{
  // check if orthogonal
  double dot1 = LAMMPS_NS::vectorDot3D(ex_local,ey_local);
  double dot2 = LAMMPS_NS::vectorDot3D(ey_local,ez_local);
  double dot3 = LAMMPS_NS::vectorDot3D(ez_local,ex_local);
  int flag = dot1 > TOLERANCE_ORTHO || dot2 > TOLERANCE_ORTHO || dot3 > TOLERANCE_ORTHO;
  if(flag) error->one(FLERR,"Insufficient accuracy: using MathExtraLiggghts::cartesian_coosys_to_local_orthogonal() for non-orthogonal coo-sys");

  // solve
  local[0] = global[0]*ex_local[0] + global[1]*ex_local[1] + global[2]*ex_local[2];
  local[1] = global[0]*ey_local[0] + global[1]*ey_local[1] + global[2]*ey_local[2];
  local[2] = global[0]*ez_local[0] + global[1]*ez_local[1] + global[2]*ez_local[2];
}

/* ----------------------------------------------------------------------
   conjugate of a quaternion: qc = conjugate of q
   assume q is of unit length
------------------------------------------------------------------------- */

void MathExtraLiggghts::qconjugate(const double * const q, double * const qc)
{
  qc[0] = q[0];
  qc[1] = -q[1];
  qc[2] = -q[2];
  qc[3] = -q[3];
}

/* ----------------------------------------------------------------------
   construct quaternion4 from vector3
------------------------------------------------------------------------- */

void MathExtraLiggghts::quat_from_vec(const double *v, double *q)
{
  q[0] = 0.;
  q[1] = v[0];
  q[2] = v[1];
  q[3] = v[2];
}

/* ----------------------------------------------------------------------
   construct vector3 from quaternion4
------------------------------------------------------------------------- */

void MathExtraLiggghts::vec_from_quat(const double *q, double * const v)
{
  v[0] = q[1];
  v[1] = q[2];
  v[2] = q[3];
}

/*----------------------------------------------------------------------
   rotoate vector by quaternion
------------------------------------------------------------------------- */

void MathExtraLiggghts::vec_quat_rotate(const double * const vec, const double * const quat, double * const result)
{
    double vecQ[4], resultQ[4], quatC[4], temp[4];

    // construct quaternion (0,vec)
    quat_from_vec(vec,vecQ);

    // conjugate initial quaternion
    qconjugate(quat,quatC);

    // rotate by quaternion multiplications
    MathExtra::quatquat(quat,vecQ,temp);
    MathExtra::quatquat(temp,quatC,resultQ);

    // return result
    vec_from_quat(resultQ,result);
}

/*----------------------------------------------------------------------
   rotoate vector by quaternion
------------------------------------------------------------------------- */

void MathExtraLiggghts::vec_quat_rotate(double * const vec, const double * const quat)
{
    double vecQ[4], resultQ[4], quatC[4], temp[4], result[3];

    // construct quaternion (0,vec)
    quat_from_vec(vec,vecQ);

    // conjugate initial quaternion
    qconjugate(quat,quatC);

    // rotate by quaternion multiplications
    MathExtra::quatquat(quat,vecQ,temp);
    MathExtra::quatquat(temp,quatC,resultQ);

    // return result
    vec_from_quat(resultQ,result);
    LAMMPS_NS::vectorCopy3D(result,vec);
}

/* ----------------------------------------------------------------------
   compute angular momentum from omega, both in space frame
   only know Idiag so need to do M = Iw in body frame
   ex,ey,ez are column vectors of rotation matrix P
   wbody = P_transpose wspace
   Mbody = Idiag wbody
   Mspace = P Mbody
------------------------------------------------------------------------- */

inline void MathExtraLiggghts::angmom_from_omega(double *w,
                                  double *ex, double *ey, double *ez,
                                  double *idiag, double *m)
{
  double mbody[3];

  mbody[0] = (w[0]*ex[0] + w[1]*ex[1] + w[2]*ex[2]) * idiag[0];
  mbody[1] = (w[0]*ey[0] + w[1]*ey[1] + w[2]*ey[2]) * idiag[1];
  mbody[2] = (w[0]*ez[0] + w[1]*ez[1] + w[2]*ez[2]) * idiag[2];

  m[0] = mbody[0]*ex[0] + mbody[1]*ey[0] + mbody[2]*ez[0];
  m[1] = mbody[0]*ex[1] + mbody[1]*ey[1] + mbody[2]*ez[1];
  m[2] = mbody[0]*ex[2] + mbody[1]*ey[2] + mbody[2]*ez[2];
}

/* ----------------------------------------------------------------------
   Check if is unit quaternion
------------------------------------------------------------------------- */

inline bool MathExtraLiggghts::is_unit_quat(const double *q)
{
    return MathExtraLiggghts::compDouble(LAMMPS_NS::vectorMag4DSquared(q),1.0,1e-6);
}

/* ----------------------------------------------------------------------
   normalize a quaternion
------------------------------------------------------------------------- */

inline void MathExtraLiggghts::quat_normalize(double *q)
{
  double norm = 1.0 / ::sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
  q[0] *= norm;
  q[1] *= norm;
  q[2] *= norm;
  q[3] *= norm;
}

/* ----------------------------------------------------------------------
   calculate the quaternion that would rotate q_old into q_new
------------------------------------------------------------------------- */

inline void MathExtraLiggghts::quat_diff(double *q_new, double *q_old, double *q_diff)
{
    double q_old_c[4];

    // q_diff = q_old^-1 * q_new
    qconjugate(q_old,q_old_c);
    MathExtra::quatquat(q_old_c,q_new,q_diff);
}

/* -----------------------------------------------------------------------------
 * compare two doubles by using their integer representation
 * source: http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
 -------------------------------------------------------------------------------*/

bool MathExtraLiggghts::compDouble(double const a, double const b, double const prec)
{
  if (a == b)
    return true;

  if (b == 0)
    return a < prec && a > -prec;

  double x = (a-b);//b;

  return x < prec && x > -prec;
}

/* -----------------------------------------------------------------------------
 * calculate barycentric coordinates of a given point (in the triangle plane)
 * should work for any point, at least analytics claim this ...
 * ap is a vector from node[0] to the point
 * edgeVec are assumed to be unit vectors
 * source: http://www.blackpawn.com/texts/pointinpoly/default.html
 * hints on _which_ barycentric coordinates are computed by this method
 * can be found on wikipedia - u_{link} = bary[1] and v_{link} = bary[2]
 -------------------------------------------------------------------------------*/

void MathExtraLiggghts::calcBaryTriCoords(double *ap, double **edgeVec, double *edgeLen, double *bary)
{
  
  double a = LAMMPS_NS::vectorDot3D(ap,edgeVec[0]);
  double b = LAMMPS_NS::vectorDot3D(ap,edgeVec[2]);
  double c = LAMMPS_NS::vectorDot3D(edgeVec[0],edgeVec[2]);
  double oneMinCSqr = 1 - c*c;

  bary[1] = (a - b*c)/(edgeLen[0] * oneMinCSqr);
  bary[2] = (a*c - b)/(edgeLen[2] * oneMinCSqr);
  bary[0] = 1. - bary[1] - bary[2];
}

void MathExtraLiggghts::calcBaryTriCoords(double *ap, double *edgeVec0, double *edgeVec1, double *edgeVec2,
                                           double *edgeLen, double *bary)
{
  UNUSED(edgeVec1);
  
  double a = LAMMPS_NS::vectorDot3D(ap,edgeVec0);
  double b = LAMMPS_NS::vectorDot3D(ap,edgeVec2);
  double c = LAMMPS_NS::vectorDot3D(edgeVec0,edgeVec2);
  double oneMinCSqr = 1 - c*c;

  bary[1] = (a - b*c)/(edgeLen[0] * oneMinCSqr);
  bary[2] = (a*c - b)/(edgeLen[2] * oneMinCSqr);
  bary[0] = 1. - bary[1] - bary[2];
}

/* ----------------------------------------------------------------------
   generate random unit quaternion
   from http://planning.cs.uiuc.edu/node198.html
------------------------------------------------------------------------- */

void MathExtraLiggghts::random_unit_quat(LAMMPS_NS::RanPark *random,double *quat)
{
    double u1 = random->uniform();
    double u2 = random->uniform();
    double u3 = random->uniform();

    double h1 = ::sqrt(1.-u1);
    double h2 = ::sqrt(u1);

    quat[0] = h1 * ::sin(2.*M_PI*u2);
    quat[1] = h1 * ::cos(2.*M_PI*u2);
    quat[2] = h2 * ::sin(2.*M_PI*u3);
    quat[3] = h2 * ::cos(2.*M_PI*u3);
}

/* ----------------------------------------------------------------------
   check if char * string is int
------------------------------------------------------------------------- */

bool MathExtraLiggghts::is_int(char *str)
{
    size_t n = strlen(str);
    for (size_t i = 0; i < n; ++i)
      if (isdigit(str[i]) == 0)
        return false;

    return true;
}

/* ----------------------------------------------------------------------
   generate complement basis
------------------------------------------------------------------------- */
void MathExtraLiggghts::generateComplementBasis(double *uVec, double *vVec, double *direction)
{

    double invLength;

    if ( abs(direction[0]) >= abs(direction[1]) )
    {
        // direction.x or direction.z is the largest magnitude component, swap them
        invLength = 1.0/::sqrt ( direction[0]*direction[0]
                              +direction[2]*direction[2]
                             );
        uVec[0] = -direction[2]*invLength;
        uVec[1] =  0.0;
        uVec[2] =  direction[0]*invLength;

        vVec[0] =  direction[1]*uVec[2];
        vVec[1] =  direction[2]*uVec[0]
                -  direction[0]*uVec[2];
        vVec[2] = -direction[1]*uVec[0];
    }
    else
    {
        // direction.y or direction.z is the largest magnitude component, swap them
        invLength = 1.0/::sqrt ( direction[1]*direction[1]
                              +direction[2]*direction[2]
                             );

        uVec[0] =  0.0;
        uVec[1] =  direction[2]*invLength;
        uVec[2] = -direction[1]*invLength;

        vVec[0] =  direction[1]*uVec[2]
                -  direction[2]*uVec[1];
        vVec[1] = -direction[0]*uVec[2];
        vVec[2] =  direction[0]*uVec[1];
    }
}

/* ----------------------------------------------------------------------
   check if integer is a prime number (primes are 6k+-1)
------------------------------------------------------------------------- */

bool MathExtraLiggghts::isPrime(int val)
{
  if (val < 2)
    return false;
  else if (val == 2)
    return true;
  else if (val == 3)
    return true;
  else if (val % 2 == 0)
    return false;
  else if (val % 3 == 0)
    return false;

  // max range is up to square-root
  int testTo = static_cast<int>(floor(::sqrt(static_cast<double>(val))));
  int test = 5;
  int width = 2;

  while (  test <= testTo )
  {
    if (val % test == 0)
      return false;

    test += width;
    width = 6 - width;
  }

  return true;
}

#endif
