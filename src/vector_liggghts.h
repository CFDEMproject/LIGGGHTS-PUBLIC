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

#ifndef LMP_VECTOR_LIGGGHTS_H
#define LMP_VECTOR_LIGGGHTS_H

#include <cmath>
#include "lammps.h"

namespace LAMMPS_NS {

//================================================
//SOME VERY SIMPLE VECTOR OPERATIONS
//================================================

template<typename T>
inline void vectorConstruct3D(T *v,T x, T y, T z)
{
  v[0] = x;
  v[1] = y;
  v[2] = z;
}

inline void vectorNormalize3D(double *v)
{
    const double norm = ::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    const double invnorm = (norm == 0.) ? 0. : 1./norm;
    v[0] *= invnorm;
    v[1] *= invnorm;
    v[2] *= invnorm;
}

inline double vectorMag3D(const double *v)
{
  return (  ::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])  );
}

inline double vectorMag3DSquared(const double *v)
{
  return (  v[0]*v[0]+v[1]*v[1]+v[2]*v[2]  );
}

inline double vectorMag4D(const double *v)
{
  return (  ::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]+v[3]*v[3])  );
}

inline double pointDistanceSquared(const double *point1, const double *point2)
{
  return
      (point1[0]-point2[0]) * (point1[0]-point2[0]) +
      (point1[1]-point2[1]) * (point1[1]-point2[1]) +
      (point1[2]-point2[2]) * (point1[2]-point2[2]);
}

inline double pointDistance(const double *point1, const double *point2)
{
  return ::sqrt(pointDistanceSquared(point1, point2));
}

inline double vectorMag4DSquared(const double *v)
{
  return (  v[0]*v[0]+v[1]*v[1]+v[2]*v[2]+v[3]*v[3]  );
}

inline double vectorDot3D(const double *v1, const double *v2)
{
  return (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]);
}

inline double vectorDot2D(const double *v1, const double *v2)
{
  return (v1[0]*v2[0]+v1[1]*v2[1]);
}

template<typename T>
inline void vectorCopy2D(const T *from, T *to)
{
  to[0]=from[0];
  to[1]=from[1];
}

inline void vectorFlip3D(double *v)
{
  v[0]=-v[0];
  v[1]=-v[1];
  v[2]=-v[2];
}

template<typename T>
inline void vectorCopyN(const T *from, T *to, int N)
{
    for(int i = 0; i < N; ++i)
       to[i] = from[i];
}

template<typename T>
inline void vectorCopy3D(const T *from, T *to)
{
  to[0]=from[0];
  to[1]=from[1];
  to[2]=from[2];
}

inline void vectorRoundN(double *vec, int N)
{
    for(int i = 0; i < N; ++i)
       vec[i] = static_cast<double>(round(vec[i]));
}

inline void vectorAbs3D(double *v)
{
    if(v[0] < 0) v[0] = -v[0];
    if(v[1] < 0) v[1] = -v[1];
    if(v[2] < 0) v[2] = -v[2];
}

template<typename T>
inline T vectorMin3D(const T *v,int &dim)
{
    if(v[0] < v[1] && v[0] < v[2])
    {
        dim = 0;
        return v[0];
    }

    if(v[1] < v[2])
    {
        dim = 1;
        return v[1];
    }

    dim = 2;
    return v[2];
}

template<typename T>
inline T vectorMin3D(const T *v)
{
    if(v[0] < v[1] && v[0] < v[2])
    return v[0];

    if(v[1] < v[2])
    return v[1];

    return v[2];
}

template<typename T>
inline T vectorMax3D(const T *v)
{
    if(v[0] > v[1] && v[0] > v[2])
    return v[0];

    if(v[1] > v[2])
    return v[1];

    return v[2];
}

template<typename T>
inline T vectorMaxN(const T *v, const int n)
{
    T max = v[0];
    for (int i=1;i<n;++i)
        max = max > v[i] ? max : v[i];
    return max;
}

template<typename T>
inline T vectorMinN(const T *v, const int n)
{
    T min = v[0];
    for (int i=1;i<n;++i)
        min = min < v[i] ? min : v[i];
    return min;
}

inline void vectorComponentMin3D(double const*v1,double const*v2,double *min)
{
    if(v1[0] > v2[0])
        min[0] = v2[0];
    else
        min[0] = v1[0];

    if(v1[1] > v2[1])
        min[1] = v2[1];
    else
        min[1] = v1[1];

    if(v1[2] > v2[2])
        min[2] = v2[2];
    else
        min[2] = v1[2];
}

inline void vectorComponentMax3D(double const*v1,double const*v2,double *max)
{
    if(v1[0] > v2[0])
        max[0] = v1[0];
    else
        max[0] = v2[0];

    if(v1[1] > v2[1])
        max[1] = v1[1];
    else
        max[1] = v2[1];

    if(v1[2] > v2[2])
        max[2] = v1[2];
    else
        max[2] = v2[2];
}

inline void vectorCopy4D(const double *from, double *to)
{
  to[0]=from[0];
  to[1]=from[1];
  to[2]=from[2];
  to[3]=from[3];
}

inline void vectorScalarMultN(int n,double *v, double s)
{
    for(int i = 0; i < n; ++i)
        v[i] = s*v[i];
}

inline void vectorScalarMultN(int n,int *v, double s)
{
    for(int i = 0; i < n; ++i)
        v[i] = static_cast<int>(static_cast<double>(s)*v[i]);
}

inline void vectorScalarMult3D(double *v, double s)
{
  v[0]=s*v[0];
  v[1]=s*v[1];
  v[2]=s*v[2];
}

inline void vectorScalarMult3D(const double *v, double s, double *result)
{
  result[0]=s*v[0];
  result[1]=s*v[1];
  result[2]=s*v[2];
}

inline void vectorScalarDiv3D(double *v, double s)
{
  const double sinv = 1./s;
  v[0]=sinv*v[0];
  v[1]=sinv*v[1];
  v[2]=sinv*v[2];
}

inline void vectorComponentDiv3D(const double *nom, const double *denom,double *result)
{
  result[0]=nom[0]/denom[0];
  result[1]=nom[1]/denom[1];
  result[2]=nom[2]/denom[2];
}

inline void vectorScalarAdd3D(double *v, double s)
{
  v[0]+=s;
  v[1]+=s;
  v[2]+=s;
}

inline void vectorScalarSubtract3D(double *v, double s)
{
  v[0]-=s;
  v[1]-=s;
  v[2]-=s;
}

inline void vectorNegate3D(const double * const v, double * const result)
{
  result[0]=-v[0];
  result[1]=-v[1];
  result[2]=-v[2];
}

inline void vectorNegate3D(double * const v)
{
  v[0]=-v[0];
  v[1]=-v[1];
  v[2]=-v[2];
}

inline void vectorScalarDiv3D(double *v, double s, double *result)
{
  const double sinv = 1./s;
  result[0]=sinv*v[0];
  result[1]=sinv*v[1];
  result[2]=sinv*v[2];
}

inline void vectorAdd3D(const double *v1, const double *v2, double *result)
{
  result[0]=v1[0]+v2[0];
  result[1]=v1[1]+v2[1];
  result[2]=v1[2]+v2[2];
}

inline void vectorAdd4D(double *v1, const double *v2)
{
  v1[0]+=v2[0];
  v1[1]+=v2[1];
  v1[2]+=v2[2];
  v1[3]+=v2[3];
}

template<typename T>
inline void vectorAddN(T *v1, const T *v2, int n)
{
  for(int i = 0; i < n; ++i)
    v1[i] += v2[i];
}

inline void vectorAddMultiple3D(const double *v1, double v2factor, const double *v2, double *result)
{
  result[0]=v1[0]+v2factor*v2[0];
  result[1]=v1[1]+v2factor*v2[1];
  result[2]=v1[2]+v2factor*v2[2];
}

inline void vectorSubtract4D(const double *v1,const double *v2, double *result)
{
  result[0]=v1[0]-v2[0];
  result[1]=v1[1]-v2[1];
  result[2]=v1[2]-v2[2];
  result[3]=v1[3]-v2[3];
}

inline void vectorSubtract3D(const double *v1,const double *v2, double *result)
{
  result[0]=v1[0]-v2[0];
  result[1]=v1[1]-v2[1];
  result[2]=v1[2]-v2[2];
}

inline void vectorSubtract2D(const double *v1,const double *v2, double *result)
{
  result[0]=v1[0]-v2[0];
  result[1]=v1[1]-v2[1];
}

template<typename T>
inline void vectorMultiN(const T *v1, const T *v2, T *result, const int n)
{
  for(int i = 0; i < n; ++i)
    result[i] = v1[i] * v2[i];
}

template<typename T>
inline void vectorComponentDivN(const T *v1, const T *v2, T *result, const int n)
{
  for(int i = 0; i < n; ++i)
    result[i] = v1[i] / v2[i];
}

inline void vectorCross3D(const double *v1,const double *v2, double *result)
{
  result[0]=v1[1]*v2[2]-v1[2]*v2[1];
  result[1]=v1[2]*v2[0]-v1[0]*v2[2];
  result[2]=v1[0]*v2[1]-v1[1]*v2[0];
}

inline double vectorCrossMag3D(const double *v1,const double *v2)
{
  double res[3];
  res[0]=v1[1]*v2[2]-v1[2]*v2[1];
  res[1]=v1[2]*v2[0]-v1[0]*v2[2];
  res[2]=v1[0]*v2[1]-v1[1]*v2[0];
  return vectorMag3D(res);
}

inline double triangleArea(const double * const v1, const double * const v2, const double * const v3, const double * const n)
{
  // formula: |((v1-v3) x (v2 - v3)).n|
  return fabs(0.5*(
    ((v1[1]-v3[1])*(v2[2]-v3[2])-(v1[2]-v3[2])*(v2[1]-v3[1]))*n[0] +
    ((v1[2]-v3[2])*(v2[0]-v3[0])-(v1[0]-v3[0])*(v2[2]-v3[2]))*n[1] +
    ((v1[0]-v3[0])*(v2[1]-v3[1])-(v1[1]-v3[1])*(v2[0]-v3[0]))*n[2]
  ));
}

template<typename T>
inline void vectorProject3D(const T *v, const T *on, T *result)
{
  T norm[3] = {on[0], on[1], on[2]};
  vectorNormalize3D(norm);
  vectorScalarMult3D(norm,vectorDot3D(v,norm),result);
}

template<typename T>
inline void vectorZeroize3D(T *v)
{
  v[0]=0;
  v[1]=0;
  v[2]=0;
}

template<typename T>
inline void vectorZeroize4D(T *v)
{
  v[0]=0;
  v[1]=0;
  v[2]=0;
  v[3]=0;
}

template<typename T>
inline void vectorZeroizeN(T *v,const int n)
{
  for(int i = 0; i < n; ++i)
     v[i]=0;
}

template<typename T>
inline void vectorInitialize3D(T *v,T init)
{
  v[0]=init;
  v[1]=init;
  v[2]=init;
}

template<typename T>
inline void vectorInitializeN(T *v,const int n,const T init)
{
  for(int i = 0; i < n; ++i)
     v[i]=init;
}

template<typename T>
inline T vectorSumN(const T * const v, const unsigned int n)
{
  T sum = 0;
  for (unsigned int i = 0; i < n; ++i)
    sum += v[i];
  return sum;
}

inline void quatIdentity4D(double *q)
{
  q[0]=1.;
  q[1]=0.;
  q[2]=0.;
  q[3]=0.;
}

inline void quatNormalize4D(double *q)
{
    const double norm = vectorMag4D(q);
    const double invnorm = (norm == 0.) ? 0. : 1./norm;
    q[0] *= invnorm;
    q[1] *= invnorm;
    q[2] *= invnorm;
    q[3] *= invnorm;
}

inline bool isIdentityQuat4D(const double *q)
{
    return
    (
        q[0] == 1. &&
        q[1] == 0. &&
        q[2] == 0. &&
        q[3] == 0.
    );
}

inline void quatMult4D(const double * const q, const double * const p, double * const result)
{
    result[0] = q[0]*p[0] - q[1]*p[1] - q[2]*p[2] - q[3]*p[3];
    result[1] = q[0]*p[1] + q[1]*p[0] + q[2]*p[3] - q[3]*p[2];
    result[2] = q[0]*p[2] - q[1]*p[3] + q[2]*p[0] + q[3]*p[1];
    result[3] = q[0]*p[3] + q[1]*p[2] - q[2]*p[1] + q[3]*p[0];
}

inline void quatMult4D(double * const q, const double * const p)
{
    double tmp[4];
    quatMult4D(q, p, tmp);
    vectorCopy4D(tmp, q);
}

inline void quatInverse4D(double *q, double *result)
{
    double invNorm = 1.0/vectorMag4DSquared(q);
    result[0] =  q[0]*invNorm;
    result[1] = -q[1]*invNorm;
    result[2] = -q[2]*invNorm;
    result[3] = -q[3]*invNorm;
}

inline void phiToQuat(const double phi, const double * const axis, double * const q)
{
    q[0] = ::cos(phi*0.5);
    const double sinPhi = ::sin(phi*0.5);
    q[1] = axis[0]*sinPhi;
    q[2] = axis[1]*sinPhi;
    q[3] = axis[2]*sinPhi;
}

inline void normalize_bary(double *v)
{
  const double mag = v[0]+v[1]+v[2];
  v[0]/=mag;
  v[1]/=mag;
  v[2]/=mag;
}

inline void vectorToBuf3D(double const*vec,double *buf,int &m)
{
  buf[m++] = vec[0];
  buf[m++] = vec[1];
  buf[m++] = vec[2];
}

inline void bufToVector3D(double *vec,double const*buf,int &m)
{
  vec[0] = buf[m++];
  vec[1] = buf[m++];
  vec[2] = buf[m++];
}

inline void vectorToBuf4D(double const*vec,double *buf,int &m)
{
  buf[m++] = vec[0];
  buf[m++] = vec[1];
  buf[m++] = vec[2];
  buf[m++] = vec[3];
}

inline void bufToVector4D(double *vec,double const*buf,int &m)
{
  vec[0] = buf[m++];
  vec[1] = buf[m++];
  vec[2] = buf[m++];
  vec[3] = buf[m++];
}

inline void vectorToBufN(double const*vec,double *buf,int &m, int nvalues)
{
  for(int i = 0; i < nvalues; ++i)
    buf[m++] = vec[i];
}

inline void bufToVectorN(double *vec,double const*buf,int &m, int nvalues)
{
  for(int i = 0; i < nvalues; ++i)
    vec[i] = buf[m++];
}

inline void printVec2D(FILE *out, const char *name, double const *vec)
{
    fprintf(out," vector %s: %e %e\n",name,vec[0],vec[1]);
}

inline void printVec3D(FILE *out, const char *name, double const*vec)
{
    fprintf(out," vector %s: %e %e %e\n",name,vec[0],vec[1],vec[2]);
}

inline void printVec3D(FILE *out, const char *name, int const*vec)
{
    fprintf(out," vector %s: %d %d %d\n",name,vec[0],vec[1],vec[2]);
}

inline void printVec4D(FILE *out, const char *name, double const*vec)
{
    fprintf(out," vector %s: %e %e %e %e\n",name,vec[0],vec[1],vec[2],vec[3]);
}

inline void printVecN(FILE *out, const char *name, double const *vec, int n)
{
    if(name) fprintf(out," vector %s:",name);
    for(int i = 0; i < n; ++i)
        fprintf(out,"%f ",vec[i]);
    fprintf(out,"\n");
}

inline void printVecN(FILE *out, const char *name, int const *vec, int n)
{
    if(name) fprintf(out," vector %s:",name);
    for(int i = 0; i < n; ++i)
        fprintf(out,"%d ",vec[i]);
    fprintf(out,"\n");
}

inline void printMat33(FILE *out, const char *name, double const* const* mat)
{
    fprintf(out," matrix %s: %f %f %f\n",name,mat[0][0],mat[0][1],mat[0][2]);
    fprintf(out,"        %s: %f %f %f\n",name,mat[1][0],mat[1][1],mat[1][2]);
    fprintf(out,"        %s: %f %f %f\n",name,mat[2][0],mat[2][1],mat[2][2]);
}

}

#endif
