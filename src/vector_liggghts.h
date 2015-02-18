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

#ifndef LMP_VECTOR_LIGGGHTS_H
#define LMP_VECTOR_LIGGGHTS_H

#include "math.h"
#include "lammps.h"

namespace LAMMPS_NS {

//================================================
//SOME VERY SIMPLE VECTOR OPERATIONS
//================================================

inline void vectorConstruct3D(double *v,double x, double y, double z)
{
  v[0] = x;
  v[1] = y;
  v[2] = z;
}

inline void vectorConstruct3D(int *v,int x, int y, int z)
{
  v[0] = x;
  v[1] = y;
  v[2] = z;
}

inline void vectorNormalize3D(double *v)
{
    double norm = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    double invnorm = (norm == 0.) ? 0. : 1./norm;
    v[0] *= invnorm;
    v[1] *= invnorm;
    v[2] *= invnorm;
}

inline double vectorMag3D(const double *v)
{
  return (  sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])  );
}

inline double vectorMag3DSquared(const double *v)
{
  return (  v[0]*v[0]+v[1]*v[1]+v[2]*v[2]  );
}

inline double vectorMag4D(const double *v)
{
  return (  sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]+v[3]*v[3])  );
}

inline double pointDistance(const double *point1, const double *point2)
{
  return
  (
     sqrt
     (
          (point1[0]-point2[0]) * (point1[0]-point2[0]) +
          (point1[1]-point2[1]) * (point1[1]-point2[1]) +
          (point1[2]-point2[2]) * (point1[2]-point2[2])
     )
  );
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

inline void vectorCopy2D(const double *from, double *to)
{
  to[0]=from[0];
  to[1]=from[1];
}

inline void vectorCopy3D(const double *from, double *to)
{
  to[0]=from[0];
  to[1]=from[1];
  to[2]=from[2];
}

inline void vectorFlip3D(double *v)
{
  v[0]=-v[0];
  v[1]=-v[1];
  v[2]=-v[2];
}

inline void vectorCopyN(const double *from, double *to, int N)
{
    for(int i = 0; i < N; i++)
       to[i] = from[i];
}

inline void vectorCopy3D(const int *from, int *to)
{
  to[0]=from[0];
  to[1]=from[1];
  to[2]=from[2];
}

inline void vectorCopy3D(const bool *from, bool *to)
{
  to[0]=from[0];
  to[1]=from[1];
  to[2]=from[2];
}

inline void vectorAbs3D(double *v)
{
    if(v[0] < 0) v[0] = -v[0];
    if(v[1] < 0) v[1] = -v[1];
    if(v[2] < 0) v[2] = -v[2];
}

inline double vectorMin3D(double *v,int &dim)
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

inline double vectorMin3D(double *v)
{
    if(v[0] < v[1] && v[0] < v[2])
    return v[0];

    if(v[1] < v[2])
    return v[1];

    return v[2];
}

inline int vectorMin3D(int *v)
{
    if(v[0] < v[1] && v[0] < v[2])
    return v[0];

    if(v[1] < v[2])
    return v[1];

    return v[2];
}

inline double vectorMax3D(double *v)
{
    if(v[0] > v[1] && v[0] > v[2])
    return v[0];

    if(v[1] > v[2])
    return v[1];

    return v[2];
}

inline int vectorMax3D(int *v)
{
    if(v[0] > v[1] && v[0] > v[2])
    return v[0];

    if(v[1] > v[2])
    return v[1];

    return v[2];
}

inline void vectorComponentMin3D(double *v1,double *v2,double *min)
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

inline void vectorComponentMax3D(double *v1,double *v2,double *max)
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
    for(int i = 0; i < n; i++)
        v[i] = s*v[i];
}

inline void vectorScalarMult3D(double *v, double s)
{
  v[0]=s*v[0];
  v[1]=s*v[1];
  v[2]=s*v[2];
}

inline void vectorScalarMult3D(double *v, double s, double *result)
{
  result[0]=s*v[0];
  result[1]=s*v[1];
  result[2]=s*v[2];
}

inline void vectorScalarDiv3D(double *v, double s)
{
  double sinv = 1./s;
  v[0]=sinv*v[0];
  v[1]=sinv*v[1];
  v[2]=sinv*v[2];
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

inline void vectorNegate3D(double *v, double *result)
{
  result[0]=-v[0];
  result[1]=-v[1];
  result[2]=-v[2];
}

inline void vectorNegate3D(double *v)
{
  v[0]=-v[0];
  v[1]=-v[1];
  v[2]=-v[2];
}

inline void vectorScalarDiv3D(double *v, double s, double *result)
{
  double sinv = 1./s;
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

inline void vectorAddN(double *v1, const double *v2, int n)
{
  for(int i = 0; i < n; i++)
    v1[i] += v2[i];
}

inline void vectorAddN(int *v1, const int *v2, int n)
{
  for(int i = 0; i < n; i++)
    v1[i] += v2[i];
}

inline void vectorAddMultiple3D(const double *v1, double v2factor, const double *v2, double *result)
{
  result[0]=v1[0]+v2factor*v2[0];
  result[1]=v1[1]+v2factor*v2[1];
  result[2]=v1[2]+v2factor*v2[2];
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

inline void vectorZeroize3D(double *v)
{
  v[0]=0.;
  v[1]=0.;
  v[2]=0.;
}

inline void vectorZeroize3D(int *v)
{
  v[0]=0;
  v[1]=0;
  v[2]=0;
}

inline void vectorZeroize4D(double *v)
{
  v[0]=0.;
  v[1]=0.;
  v[2]=0.;
  v[3]=0.;
}

inline void vectorZeroizeN(double *v,int n)
{
  for(int i = 0; i < n; i++)
     v[i]=0.;
}

inline void vectorZeroizeN(int *v,int n)
{
  for(int i = 0; i < n; i++)
     v[i]=0;
}

inline void vectorInitialize3D(double *v,double init)
{
  v[0]=init;
  v[1]=init;
  v[2]=init;
}

inline void vectorInitializeN(int *v,int n,int init)
{
  for(int i = 0; i < n; i++)
     v[i]=init;
}

inline double vectorSumN(double *v,int n)
{
  double sum = 0.;
  for(int i = 0; i < n; i++)
     sum+=v[i];
  return sum;
}

inline void quatUnitize4D(double *q)
{
  q[0]=1.;
  q[1]=0.;
  q[2]=0.;
  q[3]=0.;
}

inline bool isUnitQuat4D(double *q)
{
    return
    (
        q[0] == 1. &&
        q[1] == 0. &&
        q[2] == 0. &&
        q[3] == 0.
    );
}

inline void normalize_bary(double *v)
{
  double mag = v[0]+v[1]+v[2];
  v[0]/=mag;
  v[1]/=mag;
  v[2]/=mag;
}

inline void vectorToBuf3D(double *vec,double *buf,int &m)
{
  buf[m++] = vec[0];
  buf[m++] = vec[1];
  buf[m++] = vec[2];
}

inline void bufToVector3D(double *vec,double *buf,int &m)
{
  vec[0] = buf[m++];
  vec[1] = buf[m++];
  vec[2] = buf[m++];
}

inline void vectorToBuf4D(double *vec,double *buf,int &m)
{
  buf[m++] = vec[0];
  buf[m++] = vec[1];
  buf[m++] = vec[2];
  buf[m++] = vec[3];
}

inline void bufToVector4D(double *vec,double *buf,int &m)
{
  vec[0] = buf[m++];
  vec[1] = buf[m++];
  vec[2] = buf[m++];
  vec[3] = buf[m++];
}

inline void vectorToBufN(double *vec,double *buf,int &m, int nvalues)
{
  for(int i = 0; i < nvalues; i++)
    buf[m++] = vec[i];
}

inline void bufToVectorN(double *vec,double *buf,int &m, int nvalues)
{
  for(int i = 0; i < nvalues; i++)
    vec[i] = buf[m++];
}

inline void printVec2D(FILE *out, const char *name, double *vec)
{
    fprintf(out," vector %s: %e %e\n",name,vec[0],vec[1]);
}

inline void printVec3D(FILE *out, const char *name, double *vec)
{
    fprintf(out," vector %s: %e %e %e\n",name,vec[0],vec[1],vec[2]);
}

inline void printVec3D(FILE *out, const char *name, int *vec)
{
    fprintf(out," vector %s: %d %d %d\n",name,vec[0],vec[1],vec[2]);
}

inline void printVec4D(FILE *out, const char *name, double *vec)
{
    fprintf(out," vector %s: %e %e %e %e\n",name,vec[0],vec[1],vec[2],vec[3]);
}

inline void printVecN(FILE *out, const char *name, double *vec, int n)
{
    fprintf(out," vector %s:",name);
    for(int i = 0; i < n; i++)
        fprintf(out,"%f ",vec[i]);
    fprintf(out,"\n");
}

inline void printVecN(FILE *out, const char *name, int *vec, int n)
{
    fprintf(out," vector %s:",name);
    for(int i = 0; i < n; i++)
        fprintf(out,"%d ",vec[i]);
    fprintf(out,"\n");
}

inline void printMat33(FILE *out, const char *name, double **mat)
{
    fprintf(out," matrix %s: %f %f %f\n",name,mat[0][0],mat[0][1],mat[0][2]);
    fprintf(out,"        %s: %f %f %f\n",name,mat[1][0],mat[1][1],mat[1][2]);
    fprintf(out,"        %s: %f %f %f\n",name,mat[2][0],mat[2][1],mat[2][2]);
}

}

#endif
