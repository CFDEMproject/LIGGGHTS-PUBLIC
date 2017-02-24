//#**************************************************************
//#
//# filename:             linalg3d.cpp
//#
//# author:               Gerstmayr Johannes
//#
//# generated:            15.05.97
//# description:          Classes for linear and nonlinear algebra which is
//#												thought to be used in finite element or similar applications
//#												There are 2D and 3D and arbitrary size Vectors (Vector2D, Vector3D, Vector),
//#												arbitrary size matrices (Matrix), and a nonlinear system (NumNLSys)
//#												and a nonlinear solver (NumNLSolver)
//# remarks:							Indizes run from 1 to n except in Vector3D/2D
//#
//# Copyright (c) 2003-2013 Johannes Gerstmayr, Linz Center of Mechatronics GmbH, Austrian
//# Center of Competence in Mechatronics GmbH, Institute of Technical Mechanics at the 
//# Johannes Kepler Universitaet Linz, Austria. All rights reserved.
//#
//# This file is part of HotInt.
//# HotInt is free software: you can redistribute it and/or modify it under the terms of 
//# the HOTINT license. See folder 'licenses' for more details.
//#
//# bug reports are welcome!!!
//# WWW:		www.hotint.org
//# email:	bug_reports@hotint.org or support@hotint.org
//#**************************************************************

#include "mbs_interface.h"

extern UserOutputInterface * global_uo;

double DistToLine(const Vector3D& lp1, const Vector3D& lp2, const Vector3D& p)
{
  Vector3D vn(lp1-lp2);
  Vector3D v1(lp1-p);

  double vnl = vn.Norm();

	if (vnl == 0)
	{
		return (lp1-p).Norm();
	}
	else
	{
		return (vn.Cross(v1)).Norm() / vnl;
	}
};

double DistToLine(const Vector3D& lp, const Vector3D& lv, const Vector3D& p, Vector3D& pp)
{
  Vector3D vlp = p-lp;

  double num = lv*vlp;
  double den = lv*lv;

	pp = lp + (num/den)*lv;

	return sqrt(vlp*vlp - num * num /den);
};

double DistToLine(const Vector3D& lp1, const Vector3D& lp2, const Vector3D& p, Vector3D& pp);

double DistToLine(const Vector2D& lp1, const Vector2D& lp2, const Vector2D& p)
{
	Vector3D lp13d = lp1.MakeV3D();
	Vector3D lp23d = lp2.MakeV3D();
	Vector3D p3d = p.MakeV3D();
	return DistToLine(lp13d, lp23d, p3d);
}

//intersect two planes (plane normals pn1, pn2; plane points pp1, pp2): result=Line (represented by point lp and line direction lv)
//compute the intersection of two planes (if an intersection return 1, if planes are equal return 2, otherwise return 0)
int PlaneIntersection(const Vector3D& pn1, const Vector3D& pp1, 
											 const Vector3D& pn2, const Vector3D& pp2, 
											 Vector3D& lp, Vector3D& lv)
{
  lp = Vector3D(0.); //line point default
  lv = pn1.Cross(pn2); //line vector (must be perpendicular to both plane normals!


	if (lv.Norm() <= 1e-14) //planes are nearly parallel
	{
		lv = pp1 - pp2;
		if (pn1*lv <= 1e-16) return 1; //planes coincide
		else return 0; //planes are parallel
	}

	//find largest component of line vector
	int maxind=1; //x-coordinate of line vector is has largest value

	if (fabs(lv.Y()) > fabs(lv.X()) && fabs(lv.Y()) > fabs(lv.Z())) {maxind = 2;} //y-coordinate is largest
	else if (fabs(lv.Z()) > fabs(lv.X()) && fabs(lv.Z()) > fabs(lv.Y())) {maxind = 3;} //z-coordinate is largest

	double d1 = -(pn1*pp1);
	double d2 = -(pn2*pp2);
	switch(maxind)
	{
	case 1: 
		lp.X() = 0; 
		lp.Y() = (d2*pn1.Z() - d1*pn2.Z()) / lv.X();
		lp.Z() = (d1*pn2.Y() - d2*pn1.Y()) / lv.X();
		break;
	case 2: 
		lp.X() = (d1*pn2.Z() - d2*pn1.Z()) / lv.Y();
		lp.Y() = 0; 
		lp.Z() = (d2*pn1.X() - d1*pn2.X()) / lv.Y();
		break;
	case 3: 
		lp.X() = (d2*pn1.Y() - d1*pn2.Y()) / lv.Z();
		lp.Y() = (d1*pn2.X() - d2*pn1.X()) / lv.Z();
		lp.Z() = 0; 
		break;
	}
	return 1;
}



double Deflection(const Vector3D& lp1, const Vector3D& lp2, const Vector3D& p, const Vector3D& pn)
{
/* how to use this sensor: (D.R. 01.2011)

  lp1                                         lp2
  x-------------------------------------------x
              |
			  |
			  |
              x p

The normal distance (in this sketch vertical) of the point p to the line (lp1-lp2) is returned.
The vector pn is just needed for defining the sign.

*/
  Vector3D vn(lp1-lp2);
  Vector3D v1(lp1-p);
	double sign = 1;
	//pn defines normal to deflection (just for sign)
	if ((vn.Cross(v1))*pn < 0) sign = -1;

  double vnl = vn.Norm();

	if (vnl == 0)
	{
		return (lp1-p).Norm();
	}
	else
	{
		return sign*(vn.Cross(v1)).Norm() / vnl;
	}
};

//get deflection in 2D, distance to line segment, right side is positive
double RighthandMinDist2D(const Vector2D& lp1, const Vector2D& lp2, const Vector2D& p, Vector2D& pp)
{
  Vector2D v = lp2-lp1;
  Vector2D vlp = p-lp1;
	Vector2D n(v.Y(),-v.X());

	double sgn = 1;
	if (n*vlp < 0) sgn = -1;

  double num = v*vlp;
  double den = v*v;

	pp = lp1;
  if (num <= 0) 
	{
		//pp = lp1;
		// AP: Why sign, if num <= 0 then p is not between lp1 and lp2, the distance should always be positive
    return /*sgn**/Dist(lp1, p);
	}

  if (num >= den) 
	{
		pp = lp2;
		// AP: Why sign, if num >= den then p is not between lp1 and lp2, the distance should always be positive
    return /*sgn**/Dist(lp2, p);
	}
  
	if (den > 0)
	{
		//return DistToLine(lp1,lp2,p);
		pp = lp1 + (num/den)*v;
		//pp = lp1 + ((v*vlp)/Sqr(v.Norm()))*v;
		// AP I dont understand this formula, doesnt give correct results
		//return sgn*sqrt(fabs(vlp*vlp - num * num /den)); //for security, if (vlp*vlp - num * num /den) becomes slightly < 0
		return sgn*Dist(pp, p);
	}
	else
    return vlp.Norm();
}


double Dist(const Vector3D& p1, const Vector3D& p2)
{
	return sqrt(Sqr(p2.X()-p1.X())+Sqr(p2.Y()-p1.Y())+Sqr(p2.Z()-p1.Z()));
}

double Dist(const Vector2D& p1, const Vector2D& p2)
{
	return sqrt(Sqr(p2.X()-p1.X())+Sqr(p2.Y()-p1.Y()));
}

double Dist2(const Vector3D& p1, const Vector3D& p2)
{
	return Sqr(p2.X()-p1.X())+Sqr(p2.Y()-p1.Y())+Sqr(p2.Z()-p1.Z());
}

double Dist2(const Vector2D& p1, const Vector2D& p2)
{
	return Sqr(p2.X()-p1.X())+Sqr(p2.Y()-p1.Y());
}

double NormalizedVectorAngle(const Vector3D& v1, const Vector3D& v2)
{
	double x = v1*v2;

	if (fabs(x) < 0.707106781186547524) //if less than sqrt(0.5), then arc-cosine is more accurate
	{
		return acos(x); //returns values in [0, pi], see C++ reference of math.h at http://www.cplusplus.com/reference/clibrary/cmath/acos/
	}
	else //sine is more accurate
	{
		double y = (v1.Cross(v2)).Norm();
		if (y > 1) y = 1;

		double phi = asin(y);
		if (x < 0) phi = MY_PI-phi;

		return phi;
	}
}

double VectorAngle(Vector3D v1, Vector3D v2)
{
	v1.Normalize();
	v2.Normalize();

	return NormalizedVectorAngle(v1,v2);
}


double NormalizedVectorAngle(const Vector2D& v1, const Vector2D& v2)
{
	double x = v1*v2;

	if (fabs(x) < 0.707106781186547524) //if less than sqrt(0.5), then arc-cosine is more accurate
	{
		return acos(x); //returns values in [0, pi], see C++ reference of math.h at http://www.cplusplus.com/reference/clibrary/cmath/acos/
	}
	else //arc-sine is more accurate
	{
		double y = fabs(v1.Cross(v2));
		if (y > 1) y = 1;

		double phi = asin(y);
		if (x < 0) phi = MY_PI-phi;

		return phi;
	}
}

double VectorAngle(Vector2D v1, Vector2D v2)
{
	v1.Normalize();
	v2.Normalize();

	return NormalizedVectorAngle(v1,v2);
}

// returns parameters k and d of linear equation y = k*x+d, defined by 2 points pi=(xi,yi)
Vector2D GetLinearEquationParameters(double x1,double y1,double x2,double y2)	//$ DR 2012-04-13
{
	if(x1==x2) {assert(0);}	// x1 and x2 must not be equal!
	double k = (y2-y1)/(x2-x1);
	double d = y1-k*x1;
	return Vector2D(k,d);
}

// returns polar angle [0,2pi)
double PolarAngle(double real, double imaginary)
{
  if(imaginary == 0. && real >= 0.) return 0.;
	else if(imaginary == 0. && real >= 0.) return MY_PI;
	else return atan2(imaginary,real);
}

// returns polar angle [0,2pi)
double PolarAngle(Vector2D& v)
{
  return PolarAngle(v.X(), v.Y());
}

//Minimal Distance of p to plane (p0,n0), normalized normal, projected point in plane pp:
double DistToPlane(const Vector3D& p0, const Vector3D& n0, const Vector3D& p, Vector3D& pp)
{
  pp = ((p0-p)*n0)*n0;
	double l = pp.Norm();
	pp+=p;
	return l;
}

//Minimal Distance of p to plane(HNF)
double DistToPlane(const Vector3D& p, const Vector3D& nplane, double cplane)
{ 
	double d = (nplane*p - cplane) / nplane.Norm();
	return d;
} 

//Minimal Distance of p to plane(HNF)
double DistToPlane(const Vector2D& p, const Vector2D& nplane, double cplane)
{ 
	double d = (nplane*p - cplane) / nplane.Norm();
	return d;
} 

//Minimal Distance of p to plane(HNF)
double DistToPlane(const Vector3D& p, const Vector3D& nplane, const Vector3D& pplane)
{ 
	double d = ((pplane-p) * nplane) / nplane.Norm();
	return d;
} 

//Minimal Distance of p to plane(HNF)
double DistToPlane(const Vector2D& p, const Vector2D& nplane, const Vector2D& pplane)
{ 
	double d = ((pplane-p) * nplane) / nplane.Norm();
	return d;
} 


//Mirror Point/Vector at plane(HNF)
Vector3D MirrorAtPlane3D(const Vector3D& p, const Vector3D& nplane, double cplane) 
{
  double d = DistToPlane(p,nplane,cplane);
	return p - 2.0 * d * nplane;
}

//distance from line, including boundary:
double MinDistLP(const Vector3D& lp1, const Vector3D& lp2, const Vector3D& p)
{
  Vector3D v = lp2-lp1;
  Vector3D vlp = p-lp1;

  double num = v*vlp;
  double den = v*v;

  if (num <= 0) 
    return Dist(lp1, p);

  if (num >= den) 
    return Dist(lp2, p);
  
  if (den > 0)
    {
      return sqrt(vlp*vlp - num * num /den);
    }
  else
    return vlp.Norm();
}

//distance from line, including boundary:
double MinDistLP(const Vector3D& lp1, const Vector3D& lp2, const Vector3D& p, Vector3D& pp)
{
  Vector3D v = lp2-lp1;
  Vector3D vlp = p-lp1;

  double num = v*vlp;
  double den = v*v;

	pp = lp1;
  if (num <= 0) 
	{
    return Dist(lp1, p);
	}

  if (num >= den) 
	{
		pp = lp2;
    return Dist(lp2, p);
	}
  
	if (den > 0)
	{
		pp += (num/den)*v;
		return sqrt(vlp*vlp - num * num /den);
	}
	else
		return vlp.Norm();
}

void LocalCoordinates (const Vector3D & e1, const Vector3D & e2,
		  const Vector3D & v, double & lam1, double & lam2)
{
  double m11 = e1 * e1;
  double m12 = e1 * e2;
  double m22 = e2 * e2;
  double rs1 = v * e1;
  double rs2 = v * e2;
  
  double det = m11 * m22 - m12 * m12;
  lam1 = (rs1 * m22 - rs2 * m12)/det;
  lam2 = (m11 * rs2 - m12 * rs1)/det;
}

//minimum distance between point p and triangle
double MinDistTP (const Vector3D & tp1, const Vector3D & tp2, 
		   const Vector3D & tp3, const Vector3D & p)
{
	double lam1, lam2;
	double res;

	LocalCoordinates (tp2-tp1, tp3-tp1,	p - tp1, lam1, lam2);

	int in1 = lam1 >= 0;
	int in2 = lam2 >= 0;
	int in3 = lam1+lam2 <= 1;

	if (in1 && in2 && in3)
	{
		Vector3D pp = tp1 + lam1 * (tp2 - tp1) + lam2 *  (tp3 - tp1);
		res = Dist (p, pp);
	}
	else
	{
		res = Dist (tp1, p);
		if (!in1)
		{
			double hv = MinDistLP (tp1, tp3, p);
			if (hv < res) res = hv; 
		}
		if (!in2)
		{
			double hv = MinDistLP (tp1, tp2, p);
			if (hv < res) res = hv; 
		}
		if (!in3)
		{
			double hv = MinDistLP (tp2, tp3, p);
			if (hv < res) res = hv; 
		}
	}
	return res;
}

//minimum distance between point p and triangle
double MinDistTP (const Vector3D & tp1, const Vector3D & tp2, 
		   const Vector3D & tp3, const Vector3D & p, Vector3D& pp)
{
	double lam1, lam2;
	double res;

	LocalCoordinates (tp2-tp1, tp3-tp1,	p - tp1, lam1, lam2);

	int in1 = lam1 >= 0;
	int in2 = lam2 >= 0;
	int in3 = lam1+lam2 <= 1;

	if (in1 && in2 && in3)
	{
		pp = tp1 + lam1 * (tp2 - tp1) + lam2 *  (tp3 - tp1);
		res = Dist (p, pp);
	}
	else
	{
		res = Dist (tp1, p);
		pp = tp1;
		if (!in1)
		{
			double hv = MinDistLP(tp1, tp3, p, pp);
			if (hv < res) res = hv; 
		}
		if (!in2)
		{
			double hv = MinDistLP (tp1, tp2, p, pp);
			if (hv < res) res = hv; 
		}
		if (!in3)
		{
			double hv = MinDistLP (tp2, tp3, p, pp);
			if (hv < res) res = hv; 
		}
	}
	return res;
}

//Minimal Distance to Triangle:
double DistToTrig(const Vector3D& p1, const Vector3D& p2, const Vector3D& p3, const Vector3D& p)
{
	return MinDistTP(p1,p2,p3,p);
}

//Minimal Distance to Quad:
double DistToQuad(const Vector3D& p1, const Vector3D& p2, const Vector3D& p3, const Vector3D& p4, const Vector3D& p)
{
	return Minimum(MinDistTP(p1,p2,p3,p),MinDistTP(p1,p3,p4,p));

	/*double min = Minimum(MinDistTP(p1,p2,p3,p),MinDistTP(p2,p3,p4,p));
	min = Minimum(min, MinDistTP(p3,p4,p1,p));
	return Minimum(min,MinDistTP(p4,p1,p2,p));*/
}
int IntersectLineWithPlane(const Vector3D& pplane, const Vector3D& nplane, const Vector3D& vline, Vector3D& pline)
{
	//line equation: pline + lambda*vline = pplane
	//condition: nplane * (pplane - pline - lambda*vline) = 0

  double nv = nplane * vline;

	if (fabs(nv) <= 1e-14) //no projection possible
	{
		return 0;
	}

	double lambda = (nplane*(pplane - pline))/nv;
	pline = pline + lambda*vline;

	return 1;
}

void ProjectInPlane(const Vector3D& p1, const Vector3D& p2, const Vector3D& p3, Vector3D& pp)
{
	Vector3D nt;
	Normal3D(p1,p2,p3,nt);

	double c = - (p1.X()*nt.X() + p1.Y()*nt.Y() + p1.Z()*nt.Z());

	double nfact = -(pp.X()*nt.X() + pp.Y()*nt.Y() + pp.Z()*nt.Z() + c);

	pp = pp + (nfact) * nt;
}

void ProjectInPlane(const Vector3D& p1, const Vector3D& n, Vector3D& pp)
{
	double c = - (p1.X()*n.X() + p1.Y()*n.Y() + p1.Z()*n.Z());
	double nfact = -(pp.X()*n.X() + pp.Y()*n.Y() + pp.Z()*n.Z() + c);

	pp = pp + (nfact) * n;
}


double NormSym(const Matrix3D& m)
{
  return sqrt(Sqr(m(1,1))+Sqr(m(2,2))+Sqr(m(3,3))+
    2*Sqr(m(2,3))+2*Sqr(m(1,3))+2*Sqr(m(1,2)));
}

Matrix3D Deviator(const Matrix3D& m)
{
  Matrix3D dev(m);
  double tr = dev.Trace();
  tr /= 3.;
  dev(1,1) -= tr;
  dev(2,2) -= tr;
  dev(3,3) -= tr;
  return dev;
}



double DoubleProd3D(const Matrix3D& m, const Matrix3D& n)
{
	//return (m.GetTp()*n).Trace();

	double dp=0;
	for (int i=1; i <= 3; i++)
	{  
		for (int j=1; j <= 3; j++)
		{  
			dp+=m(i,j)*n(i,j);
		}
	}
	return dp;


}

double Area(const Vector3D& p1, const Vector3D& p2, const Vector3D& p3)
{
	Vector3D v1 = p2-p1;
	Vector3D v2 = p3-p1;

	return 0.5*(v1.Cross(v2)).Norm();
}

double Area(const Vector3D& p1, const Vector3D& p2, const Vector3D& p3, const Vector3D& p4)
{
	Vector3D v1 = p2-p1;
	Vector3D v2 = p4-p1;
	Vector3D v3 = p4-p3;
	Vector3D v4 = p2-p3;

	return 0.5*((v1.Cross(v2)).Norm()+(v3.Cross(v4)).Norm());
}


//int IsQuadConvex2D(Vector3D& p1, Vector3D& p2, Vector3D& p3, Vector3D& p4)
//{
//	return 0;	
//}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Matrix3D RotMatrix1(double p)
{
	return Matrix3D(1.,0.,0.,
		0.,cos(p),-sin(p),
		0.,sin(p), cos(p));
}

Matrix3D RotMatrix2(double p)
{
	return Matrix3D(cos(p),0.,sin(p),
		0.,1.,0.,
		-sin(p),0., cos(p));
}

Matrix3D RotMatrix3(double p)
{
	return Matrix3D(cos(p),-sin(p),0.,
		sin(p), cos(p),0.,
		0.,0.,1.);
}

double GetRotation(int rotaxis, const Matrix3D& rotmat)
{
	double s1,s2,c1,c2;
	if (rotaxis == 1) //rotation around x-axis
	{
		c1 = rotmat(2,2);
		s1 = rotmat(3,2);
		c2 = rotmat(3,3);
		s2 =-rotmat(2,3);
	}
	else if (rotaxis == 2) //rotation around y-axis
	{
		c1 = rotmat(1,1);
		s1 = rotmat(3,1);
		c2 = rotmat(3,3);
		s2 =-rotmat(1,3);
	}
	else if (rotaxis == 3) //rotation around z-axis
	{
		c1 = rotmat(1,1);
		s1 = rotmat(2,1);
		c2 = rotmat(2,2);
		s2 =-rotmat(1,2);
	}

	return atan2(s1+s2,c1+c2);
}

Vector3D BaseVector3D(int i)
{
	switch(i)
	{
	case 1: return Vector3D(1.,0.,0.); break;
	case 2: return Vector3D(0.,1.,0.); break;
	case 3: return Vector3D(0.,0.,1.); break;
	default: return Vector3D(1.,1.,1.);
	}
}

Vector2D BaseVector2D(int i)
{
	switch(i)
	{
	case 1: return Vector2D(1.,0.); break;
	case 2: return Vector2D(0.,1.); break;
	default: return Vector2D(1.,1.);
	}
}


//-----------------------------------------------------------------
//--------------    CLASS VECTOR2D   ------------------------------
//-----------------------------------------------------------------

Vector2D operator* (const Matrix& m, const Vector2D& v)
{
	Vector2D v2;
	for (int i=1; i <= 2; i++)
	{
		for (int j=1; j <= 2; j++)
		{
			v2[i-1]+=v[j-1]*m(i,j);
		}
	}
	return v2;
}

//Vector2D::Vector2D( Vector3D& veci)
//{
//	vec[0] = veci[0]; vec[1] = veci[1];
//}

Vector3D Vector2D::MakeV3D() const
{
	return Vector3D(vec[0],vec[1],0.0);
}

Matrix GetRot(double phi)
{
	return Matrix(cos(phi),-sin(phi),
		sin(phi), cos(phi));
}
void GetRot(double phi, Matrix& m)
{
	m(1,1) = cos(phi);
	m(1,2) =-sin(phi);
	m(2,1) = sin(phi);
	m(2,2) = cos(phi);
}

void GetRotT(double phi, Matrix& m)
{
	m(1,1) = cos(phi);
	m(1,2) = sin(phi);
	m(2,1) =-sin(phi);
	m(2,2) = cos(phi);
}

void Rotate(Vector2D& v, double phi) 
{
	v = GetRot(phi)*v;
}


//-----------------------------------------------------------------
// CLASS VECTOR3D   CLASS VECTOR3D   CLASS VECTOR3D   CLASS VECTOR3D
//-----------------------------------------------------------------

Vector3D operator* (const Matrix& m, const Vector3D& v)
{
	Vector3D v2;
	for (int i=1; i <= 3; i++)
	{
		for (int j=1; j <= 3; j++)
		{
			v2[i-1]+=v[j-1]*m(i,j);
		}
	}
	return v2;
}

//Vector3D::Vector3D(const Vector2D& avec)
//{
//	vec[0] = avec.X(); vec[1] = avec.Y(); vec[2] = 0.0;
//}

Vector2D Vector3D::MakeV2D() const
{
	return Vector2D(vec[0],vec[1]);
}

void Normal3D(const Vector3D& p1,const Vector3D& p2,const Vector3D& p3, Vector3D& n)
{
	Vector3D v1 = p2-p1;
	Vector3D v2 = p3-p1;

	n = v1.Cross(v2);
	n.Normalize();
}

double Norm3D(double x, double y, double z)
{
	return sqrt(x*x+y*y+z*z);
}

Matrix GetRot1(double p)
{
	return Matrix(1.,0.,0.,
		0.,cos(p),-sin(p),
		0.,sin(p), cos(p));
}
Matrix GetRot2(double p)
{
	return Matrix(cos(p),0.,-sin(p),
		0.,1.,0.,
		sin(p),0., cos(p));
}
Matrix GetRot3(double p)
{
	return Matrix(cos(p),-sin(p),0.,
		sin(p), cos(p),0.,
		0.,0.,1.);
}

Matrix GetRot(double p1,double p2,double p3)
{
	return GetRot1(p3)*GetRot2(p2)*GetRot3(p1);
}





//-----------------------------------------------------------------
// CLASS MATRIX2D   CLASS MATRIX2D   CLASS MATRIX2D   CLASS MATRIX2D
//-----------------------------------------------------------------


int Matrix2D::GetInverse(Matrix2D &m2)
{
#ifndef __QUICKMATH
	release_assert(m2.cols == 2);
	release_assert(m2.rows == 2);
#endif
	double det = mat[0]*mat[3] - mat[1]*mat[2];
	if (det == 0) { return 0; }
	det = 1./det;
	m2.mat[0] =  det*mat[3];  m2.mat[1] = -det*mat[1];
	m2.mat[2] = -det*mat[2];  m2.mat[3] =  det*mat[0];

	return 1;
}

void Mult(const Matrix2D& m, const Vector2D& v, Vector& res)
{
	int cols = m.Getcols(); //-->must be 2!!!
	int rows = m.Getrows();
	if (res.Length() < rows) res.SetLen(rows);

	for (int i=0; i < rows; i++)
	{
		res.vec[i] = m.mat[i*cols]*v.vec[0]+m.mat[i*cols+1]*v.vec[1];
	}
}

Matrix2D operator+ (const Matrix2D& m1, const Matrix2D& m2)
{
#ifndef __QUICKMATH
	release_assert(m1.cols==m2.cols && m1.rows==m2.rows);
#endif
	Matrix2D m;
	m.cols=m1.cols;
	m.rows=m1.rows;
	for (int i=0; i < m1.rows*m1.cols; i++) m.mat[i] = m1.mat[i]+m2.mat[i];
	return m;
}
Matrix2D operator- (const Matrix2D& m1, const Matrix2D& m2)
{
#ifndef __QUICKMATH
	release_assert(m1.cols==m2.cols && m1.rows==m2.rows);
#endif
	Matrix2D m;
	m.cols=m1.cols;
	m.rows=m1.rows;
	for (int i=0; i < m1.rows*m1.cols; i++) m.mat[i] = m1.mat[i]-m2.mat[i];
	return m;
}
Matrix2D operator* (const Matrix2D& m1, const double& val)
{
	Matrix2D m;
	m.cols=m1.cols;
	m.rows=m1.rows;
	for (int i=0; i < m1.rows*m1.cols; i++) m.mat[i] = val*m1.mat[i];
	return m;
}
Matrix2D operator* (const double& val, const Matrix2D& m1)
{
	Matrix2D m;
	m.cols=m1.cols;
	m.rows=m1.rows;
	for (int i=0; i < m1.rows*m1.cols; i++) m.mat[i] = val*m1.mat[i];
	return m;
}
Matrix2D operator* (const Matrix2D& m1, const Matrix2D& m2)
{
#ifndef __QUICKMATH
	release_assert(m1.cols==m2.rows);
#endif
	Matrix2D m;
	m.rows = m1.rows;
	m.cols = m2.cols;
	for (int i=0;i<m.rows;i++)
	{
		for (int j=0;j<m.cols;j++)
		{
			m.mat[i*m.cols+j] = 0;
			for (int k=0;k<m1.cols;k++)
			{
				m.mat[i*m.cols+j]+=m1.mat[i*m1.cols+k]*m2.mat[k*m2.cols+j];
			}
		}
	}
	return m;
}
Vector2D operator* (const Matrix2D& m, const Vector2D& v)
{
#ifndef __QUICKMATH
	release_assert(m.cols==2 && m.rows==2);
#endif
	Vector2D res;
	for (int i=0; i < 2; i++)
	{
		res.vec[i] = m.mat[i*2]*v.vec[0]+m.mat[i*2+1]*v.vec[1];
	}
	return res;
}

Vector2D operator* (const Vector2D& v, const Matrix2D& m)
{
#ifndef __QUICKMATH
	release_assert(m.cols==2 && m.rows==2);
#endif
	Vector2D res;
	for (int i=0; i < 2; i++)
	{
		res.vec[i] = m.mat[i]*v.vec[0]+m.mat[2+i]*v.vec[1];
	}
	return res;
}

Vector2D operator* (const Matrix2D& m, const Vector& v)
{
#ifndef __QUICKMATH
	release_assert(m.cols==v.l && m.rows==2);
#endif
	Vector2D res;
	for (int i=0; i < 2; i++)
	{
		res.vec[i] = 0;
		for (int j=0; j < v.l; j++)
		{
			res.vec[i] += m.Get0(i,j)*v.vec[j];
		}
	}
	return res;
}

//transform strain type vector to matrix (by PG)
void StrainVectorToMatrix2D(Matrix2D& m, const Vector& v)
{
#ifndef __QUICKMATH
	release_assert(v.Length() >= 3 && m.Getcols() >= 2 && m.Getrows() >= 2);
#endif
	//  0 1
	//  2 3
	m.mat[0] = v.vec[0];
	m.mat[3] = v.vec[1];
	m.mat[1] = m.mat[2] = 0.5*v.vec[2];
}

void StrainVectorToMatrix2D(Matrix2D& m, const Vector3D& v)
{
	//  0 1
	//  2 3
	m.mat[0] = v.vec[0];
	m.mat[3] = v.vec[1];
	m.mat[1] = m.mat[2] = 0.5*v.vec[2];
}

//transform stress type vector to matrix (by PG)
void StressVectorToMatrix2D(Matrix2D& m, const Vector& v)
{
#ifndef __QUICKMATH
	release_assert(v.Length() >= 3 && m.Getcols() >= 2 && m.Getrows() >= 2);
#endif
	//  0 1
	//  2 3
	m.mat[0] = v.vec[0];
	m.mat[3] = v.vec[1];
	m.mat[1] = m.mat[2] = v.vec[2];
}

//transform matrix to strain type vector (by PG)
void Matrix2DToStrainVector(Vector& v, const Matrix2D& m)
{
#ifndef __QUICKMATH
	release_assert(v.Length() >= 3 && m.Getcols() >= 2 && m.Getrows() >= 2);
#endif
	//  0 1
	//  2 3
	v.vec[0] = m.mat[0];
	v.vec[1] = m.mat[3];
	v.vec[2] = 2.*m.mat[1];
}

//transform matrix to strain type vector (by PG)
void Matrix2DToStressVector(Vector& v, const Matrix2D& m)
{
#ifndef __QUICKMATH
	release_assert(v.Length() >= 3 && m.Getcols() >= 2 && m.Getrows() >= 2);
#endif
	//  0 1
	//  2 3
	v.vec[0] = m.mat[0];
	v.vec[1] = m.mat[3];
	v.vec[2] = m.mat[1];
}





//-----------------------------------------------------------------
// CLASS MATRIX3D   CLASS MATRIX3D   CLASS MATRIX3D   CLASS MATRIX3D
//-----------------------------------------------------------------
Matrix3D::Matrix3D(const Matrix& mat)
	{		
		rows = mat.Getrows();
		cols = mat.Getcols();

		if (rows < 3 || rows > 4 || cols < 3 || cols > 4) {assert(0 && "Linalg3D::Matrix3D(const Matrix& mat)");}

		for (int i=0; i<rows; i++)
		{
			for (int j=0; j<cols; j++)
			{
				Get0(i,j) = mat.Get0(i,j);
			}
		}
	}
//returns 1, if sucsess
int Matrix3D::GetInverse(Matrix3D& m2)
{
	if (cols == 3)
	{
		m2.cols=3; m2.rows=3;
		double det = mat[0]*mat[4]*mat[8]-mat[0]*mat[5]*mat[7]-mat[3]*mat[1]*mat[8]+mat[3]*mat[2]*mat[7]+mat[6]*mat[1]*mat[5]-mat[6]*mat[2]*mat[4];
		if (det == 0) return 0;
		det=1./det;
		m2.mat[0] = det*(mat[4]*mat[8]-mat[5]*mat[7]);
		m2.mat[1] = det*(-mat[1]*mat[8]+mat[2]*mat[7]);
		m2.mat[2] = det*(mat[1]*mat[5]-mat[2]*mat[4]);
		m2.mat[3] = det*(-mat[3]*mat[8]+mat[5]*mat[6]);
		m2.mat[4] = det*(mat[0]*mat[8]-mat[2]*mat[6]);
		m2.mat[5] = det*(-mat[0]*mat[5]+mat[2]*mat[3]);
		m2.mat[6] = det*(mat[3]*mat[7]-mat[4]*mat[6]);
		m2.mat[7] = det*(-mat[0]*mat[7]+mat[1]*mat[6]);
		m2.mat[8] = det*(mat[0]*mat[4]-mat[1]*mat[3]);
		return 1;
	}
	if (cols == 2)
	{
		m2.cols=2; m2.rows=2;
		double det=mat[0]*mat[3]-mat[1]*mat[2];
		if (det==0) {return 0;}
		det=1./det;
		m2.mat[0]= det*mat[3]; m2.mat[1]=-det*mat[1];
		m2.mat[2]=-det*mat[2]; m2.mat[3]= det*mat[0];
		return 1;
	}
	return 0;
}


Matrix3D operator+ (const Matrix3D& m1, const Matrix3D& m2)
{
#ifndef __QUICKMATH
	release_assert(m1.cols==m2.cols && m1.rows==m2.rows);
#endif
	Matrix3D m;
	m.cols=m1.cols;
	m.rows=m1.rows;
	for (int i=0; i < m1.rows*m1.cols; i++) m.mat[i] = m1.mat[i]+m2.mat[i];
	return m;
}
Matrix3D operator- (const Matrix3D& m1, const Matrix3D& m2)
{
#ifndef __QUICKMATH
	release_assert(m1.cols==m2.cols && m1.rows==m2.rows);
#endif
	Matrix3D m;
	m.cols=m1.cols;
	m.rows=m1.rows;
	for (int i=0; i < m1.rows*m1.cols; i++) m.mat[i] = m1.mat[i]-m2.mat[i];
	return m;
}
Matrix3D operator* (const Matrix3D& m1, const double& val)
{
	Matrix3D m;
	m.cols=m1.cols;
	m.rows=m1.rows;
	for (int i=0; i < m1.rows*m1.cols; i++) m.mat[i] = val*m1.mat[i];
	return m;
}
Matrix3D operator* (const double& val, const Matrix3D& m1)
{
	Matrix3D m;
	m.cols=m1.cols;
	m.rows=m1.rows;
	for (int i=0; i < m1.rows*m1.cols; i++) m.mat[i] = val*m1.mat[i];
	return m;
}

Matrix3D operator* (const Matrix3D& m1, const Matrix3D& m2)
{
#ifndef __QUICKMATH
	release_assert(m1.cols==m2.rows);
#endif
	Matrix3D m;
	m.rows = m1.rows;
	m.cols = m2.cols;
	for (int i=0;i<m.rows;i++)
	{
		for (int j=0;j<m.cols;j++)
		{
			m.mat[i*m.cols+j] = 0;
			for (int k=0;k<m1.cols;k++)
			{
				m.mat[i*m.cols+j]+=m1.mat[i*m1.cols+k]*m2.mat[k*m2.cols+j];
			}
		}
	}
	return m;
}

Vector3D operator* (const Matrix3D& m, const Vector3D& v)
{
#ifndef __QUICKMATH
	release_assert(m.cols==3 && m.rows==3);
#endif
	Vector3D res;
	for (int i=0; i < 3; i++)
	{
		res.vec[i] = m.mat[i*3]*v.vec[0]+m.mat[i*3+1]*v.vec[1]+m.mat[i*3+2]*v.vec[2];
	}
	return res;
}

Vector3D operator* (const Vector3D& v, const Matrix3D& m)
{
#ifndef __QUICKMATH
	release_assert(m.cols==3 && m.rows==3);
#endif
	Vector3D res;
	for (int i=0; i < 3; i++)
	{
		res.vec[i] = m.mat[i]*v.vec[0]+m.mat[3+i]*v.vec[1]+m.mat[6+i]*v.vec[2];
	}
	return res;
}

Vector3D operator* (const Matrix3D& m, const Vector& v)
{
#ifndef __QUICKMATH
	release_assert(m.cols==v.l && m.rows==3);
#endif
	Vector3D res;
	for (int i=0; i < 3; i++)
	{
		res.vec[i] = 0;
		for (int j=0; j < v.l; j++)
		{
			res.vec[i] += m.Get0(i,j)*v.vec[j];
		}
	}
	return res;
}

#define __DEPREC_MSG "**this is a deprecated function**"
__declspec(deprecated(__DEPREC_MSG)) Vector2D operator* (const Matrix3D& m, const Vector2D& v)
{
#ifndef __QUICKMATH
	release_assert(m.cols * m.rows >= 4);
#endif
	Vector2D res;
	for (int i=0; i < 2; i++)
	{
		res.vec[i] = m.mat[i*2]*v.vec[0]+m.mat[i*2+1]*v.vec[1];
	}
	return res;
}

void Mult(const Matrix3D& m, const Vector3D& v, Vector& res)
{
	int cols = m.Getcols(); //-->must be 3!!!
	int rows = m.Getrows();
	if (res.Length() < rows) res.SetLen(rows);

	for (int i=0; i < rows; i++)
	{
		res.vec[i] = m.mat[i*cols]*v.vec[0]+m.mat[i*cols+1]*v.vec[1]+m.mat[i*cols+2]*v.vec[2];
	}
}

void Mult(const Matrix3D& m, const Vector3D& v, Vector3D& res)
{
	int cols = 3; //-->must be 3!!!
	int rows = 3;

	for (int i=0; i < 3; i++)
	{
		res.vec[i] = m.mat[i*3]*v.vec[0]+m.mat[i*3+1]*v.vec[1]+m.mat[i*3+2]*v.vec[2];
	}
}

//transform strain type vector to matrix (by PG)
//takes the first three entries of v as the diagonal entries of m
//and the second three entries as the shear components
void StrainVectorToMatrix3D(Matrix3D& m, const Vector& v)
{
#ifndef __QUICKMATH
	release_assert(v.Length() >= 6 && m.Getcols() >= 3 && m.Getrows() >= 3);
#endif
	//  0 1 2
	//  3 4 5
	//  6 7 8
	m.mat[0] = v.vec[0];
	m.mat[4] = v.vec[1];
	m.mat[8] = v.vec[2];
	m.mat[5] = m.mat[7] = 0.5*v.vec[3];
	m.mat[2] = m.mat[6] = 0.5*v.vec[4];
	m.mat[1] = m.mat[3] = 0.5*v.vec[5];
}

//transform stress type vector to matrix (by PG)
//takes the first three entries of v as the diagonal entries of m
//and the second three entries as the shear components
void StressVectorToMatrix3D(Matrix3D& m, const Vector& v)
{
#ifndef __QUICKMATH
	release_assert(v.Length() >= 6 && m.Getcols() >= 3 && m.Getrows() >= 3);
#endif
	//  0 1 2
	//  3 4 5
	//  6 7 8
	m.mat[0] = v.vec[0];
	m.mat[4] = v.vec[1];
	m.mat[8] = v.vec[2];
	m.mat[5] = m.mat[7] = v.vec[3];
	m.mat[2] = m.mat[6] = v.vec[4];
	m.mat[1] = m.mat[3] = v.vec[5];
}

//transform matrix to strain type vector (by PG)
//takes the diagonal entries of m as the first three entries of v 
//and the shear components as the  second three entries
void Matrix3DToStrainVector(Vector& v, const Matrix3D& m)
{
#ifndef __QUICKMATH
	release_assert(v.Length() >= 6 && m.Getcols() >= 3 && m.Getrows() >= 3);
#endif
	//  0 1 2
	//  3 4 5
	//  6 7 8
	v.vec[0] = m.mat[0];
	v.vec[1] = m.mat[4];
	v.vec[2] = m.mat[8];
	v.vec[3] = 2.*m.mat[7];
	v.vec[4] = 2.*m.mat[6];
	v.vec[5] = 2.*m.mat[3];
}

//transform matrix to stress type vector (by PG)
//takes the diagonal entries of m as the first three entries of v 
//and the shear components as the  second three entries
void Matrix3DToStressVector(Vector& v, const Matrix3D& m)
{
#ifndef __QUICKMATH
	release_assert(v.Length() >= 6 && m.Getcols() >= 3 && m.Getrows() >= 3);
#endif
	//  0 1 2
	//  3 4 5
	//  6 7 8
	v.vec[0] = m.mat[0];
	v.vec[1] = m.mat[4];
	v.vec[2] = m.mat[8];
	v.vec[3] = m.mat[7];
	v.vec[4] = m.mat[6];
	v.vec[5] = m.mat[3];
}


//-----------------------------------------------------------------
// CLASS MATRIXXD   CLASS MATRIXXD   CLASS MATRIXXD   CLASS MATRIXXD
//-----------------------------------------------------------------

void Mult(const MatrixXD& m, const Vector& v, Vector& res)
{
	int cols = m.Getcols();
	int rows = m.Getrows();
	if (res.Length() < rows) res.SetLen(rows);

	for (int i=0; i < rows; i++)
	{
		res.vec[i] = 0;
		for (int j=0; j < cols; j++)
		{
			res.vec[i] += m.mat[i*cols+j]*v.vec[j];
		}
	}
}

mystr& Vector3D::MakeString(mystr intro)
{ // always 3 double numbers
	mystr* buffer = new mystr(intro);
	*buffer +=("[3] = { "); 
	*buffer += mystr(X()); *buffer += ", ";
	*buffer += mystr(Y()); *buffer += ", ";
	*buffer += mystr(Z()); *buffer += mystr(" };");
	return *buffer;
}

// writes out MatrixXD to a mystr intro+"[][]= { , , }"
mystr& MatrixXD::MakeString(mystr intro)
{
	mystr* buffer = new mystr(intro);
	*buffer +=("["); *buffer += mystr(rows); *buffer +=("]");
	*buffer +=("["); *buffer += mystr(cols); *buffer +=("]");
	*buffer +=(" = { ");
	for (int i=0; i < rows; i++)
	{
		for (int j=0; j < cols; j++)
		{
			if( i || j ) *buffer += ", ";
			*buffer += mystr(Get0(i,j));
		}
	}
	*buffer += mystr(" };");
	return *buffer;
}

//writes out MatrixXD with constant width and factor normalized
ostream& operator<<(ostream& os, const MatrixXD& m)
{
	for (int i=0; i < m.Getrows(); i++)
	{
		os << "[";
		for (int j=0; j < m.Getcols() - 1; j++)
		{
			os << m.Get0(i,j) << " ";
		}
		os << m.Get0(i,m.Getcols() - 1) << "]\n";
	}
	return os;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//SEARCHTREE
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//return 6 indices for box: minx, maxx, miny, maxy, minz, maxz
void SearchTree::GetBoxIndizes(const Box3D& b, int* ind) const
	{
		ind[0] = IndX(b.PMin().X());
		ind[1] = IndX(b.PMax().X());
		ind[2] = IndY(b.PMin().Y());
		ind[3] = IndY(b.PMax().Y());
		ind[4] = IndZ(b.PMin().Z());
		ind[5] = IndZ(b.PMax().Z());
	}

	//return items in box defined by 6 indices: minx, maxx, miny, maxy, minz, maxz
	void SearchTree::GetItemsInBox(int* ind, TArray<int>& items) const
	{
		if(!data) return;
		items.SetLen(0);
		int ld;
		for (int ix=ind[0]; ix <= ind[1]; ix++)
		{
			for (int iy=ind[2]; iy <= ind[3]; iy++)
			{
				for (int iz=ind[4]; iz <= ind[5]; iz++)
				{
					ld = data[Index(ix,iy,iz)].Length();
					for (int i = 1; i <= ld; i++)
						items.Add((data[Index(ix,iy,iz)])(i));
				}
			}
		}
	}

	//add items in box defined by 6 indices: minx, maxx, miny, maxy, minz, maxz
	//does not reset items list
	void SearchTree::AddItemsInBox(int* ind, TArray<int>& items) const
	{
		if(!data) return;
		int ld;
		for (int ix=ind[0]; ix <= ind[1]; ix++)
		{
			for (int iy=ind[2]; iy <= ind[3]; iy++)
			{
				for (int iz=ind[4]; iz <= ind[5]; iz++)
				{
					ld = data[Index(ix,iy,iz)].Length();
					for (int i = 1; i <= ld; i++)
						items.Add((data[Index(ix,iy,iz)])(i));
				}
			}
		}
	}

	//get only items in box
	void SearchTree::GetItemsInBox(const Box3D& b, TArray<int>& items) const
	{
		items.SetLen(0);
		int ind[6];
		ind[0] = IndX(b.PMin().X());
		ind[1] = IndX(b.PMax().X());
		ind[2] = IndY(b.PMin().Y());
		ind[3] = IndY(b.PMax().Y());
		ind[4] = IndZ(b.PMin().Z());
		ind[5] = IndZ(b.PMax().Z());
		int ix, iy, iz, i;
		IVector* id;

		for (ix=ind[0]; ix <= ind[1]; ix++)
		{
			for (iy=ind[2]; iy <= ind[3]; iy++)
			{
				for (iz=ind[4]; iz <= ind[5]; iz++)
				{
					id = &data[Index(ix,iy,iz)];
					for (i = 1; i <= id->Length(); i++)
						items.Add(id->Get(i));
				}
			}
		}
		

		//AddItemsInBox(b, items);
	}

	//get items in box, do not reset items list
	void SearchTree::AddItemsInBox(const Box3D& b, TArray<int>& items) const
	{
		int ind[6];
		ind[0] = IndX(b.PMin().X());
		ind[1] = IndX(b.PMax().X());
		ind[2] = IndY(b.PMin().Y());
		ind[3] = IndY(b.PMax().Y());
		ind[4] = IndZ(b.PMin().Z());
		ind[5] = IndZ(b.PMax().Z());
		int ix, iy, iz, i;
		IVector* id;

		for (ix=ind[0]; ix <= ind[1]; ix++)
		{
			for (iy=ind[2]; iy <= ind[3]; iy++)
			{
				for (iz=ind[4]; iz <= ind[5]; iz++)
				{
					id = &data[Index(ix,iy,iz)];
					for (i = 1; i <= id->Length(); i++)
						items.Add(id->Get(i));
				}
			}
		}
	/*
		for (int ix=IndX(b.PMin().X()); ix <= IndX(b.PMax().X()); ix++)
		{
			for (int iy=IndY(b.PMin().Y()); iy <= IndY(b.PMax().Y()); iy++)
			{
				for (int iz=IndZ(b.PMin().Z()); iz <= IndZ(b.PMax().Z()); iz++)
				{
					for (int i = 1; i <= data[Index(ix,iy,iz)].Length(); i++)
						items.Add((data[Index(ix,iy,iz)])(i));
				}
			}
		}
	*/
	}

	void SearchTree::AddItem(const Box3D& b, int identifier)
	{
		for (int ix=IndX(b.PMin().X()); ix <= IndX(b.PMax().X()); ix++)
		{
			for (int iy=IndY(b.PMin().Y()); iy <= IndY(b.PMax().Y()); iy++)
			{
				for (int iz=IndZ(b.PMin().Z()); iz <= IndZ(b.PMax().Z()); iz++)
				{
					data[Index(ix,iy,iz)].Add(identifier);
				}
			}
		}
	}

	//return i-index for a double i-value in Box
	int SearchTree::IndX(double x) const 
	{ 
		int i = (int)(((x-box.PMin().X())*(double)sx)/box.SizeX());
		if (i < 0) i = 0;
		if (i >= sx) i = sx-1;
		return i;
	}
	//return i-index for a double i-value in Box
	int SearchTree::IndY(double x) const 
	{ 
		int i = (int)(((x-box.PMin().Y())*(double)sy)/box.SizeY());
		if (i < 0) i = 0;
		if (i >= sy) i = sy-1;
		return i;
	}
	//return i-index for a double i-value in Box
	int SearchTree::IndZ(double x) const 
	{ 
		int i = (int)(((x-box.PMin().Z())*(double)sz)/box.SizeZ());
		if (i < 0) i = 0;
		if (i >= sz) i = sz-1;
		return i;
	}

	int SearchTree::Index(int x, int y, int z) const
	{
		return x+y*sx+z*(sx*sy);
	}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//ROTATION MATRICES
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



Matrix3D ComputeRotMatrixFromQuaternions(double b0, double b1, double b2, double b3)
{
	Matrix3D rot;

	rot.Get0(0,0)= -2.0*b3*b3-2.0*b2*b2+1.0;
	rot.Get0(0,1)= -2.0*b3*b0+2.0*b2*b1;
	rot.Get0(0,2)= 2.0*b3*b1+2.0*b2*b0;
	rot.Get0(1,0)= 2.0*b3*b0+2.0*b2*b1;
	rot.Get0(1,1)= -2.0*b3*b3-2.0*b1*b1+1.0;
	rot.Get0(1,2)= 2.0*b3*b2-2.0*b1*b0;
	rot.Get0(2,0)= -2.0*b2*b0+2.0*b3*b1;
	rot.Get0(2,1)= 2.0*b3*b2+2.0*b1*b0;
	rot.Get0(2,2)= -2.0*b2*b2-2.0*b1*b1+1.0;

	return rot;
}



Matrix3D ComputeRotMatrixEulerParam(const double& beta0, const double& beta1,
																	 const double& beta2, const double& beta3)
{
	//shall be copy of function from linalg3d
	Matrix3D rot;

	rot.Get0(0,0)= -2.0*beta3*beta3-2.0*beta2*beta2+1.0;
	rot.Get0(0,1)= -2.0*beta3*beta0+2.0*beta2*beta1;
	rot.Get0(0,2)= 2.0*beta3*beta1+2.0*beta2*beta0;
	rot.Get0(1,0)= 2.0*beta3*beta0+2.0*beta2*beta1;
	rot.Get0(1,1)= -2.0*beta3*beta3-2.0*beta1*beta1+1.0;
	rot.Get0(1,2)= 2.0*beta3*beta2-2.0*beta1*beta0;
	rot.Get0(2,0)= -2.0*beta2*beta0+2.0*beta3*beta1;
	rot.Get0(2,1)= 2.0*beta3*beta2+2.0*beta1*beta0;
	rot.Get0(2,2)= -2.0*beta2*beta2-2.0*beta1*beta1+1.0;
	return rot;
}

Matrix3D ComputeRotMatrixPEulerParam(const double& beta0, const double& beta1,
																		const double& beta2, const double& beta3, const double& betap0, const double& betap1,
																		const double& betap2, const double& betap3)
{
	//shall be copy of function from linalg3d
	Matrix3D rotp;
	rotp.Get0(0,0)= -4.0*beta3*betap3-4.0*beta2*betap2;
	rotp.Get0(0,1)= -2.0*betap3*beta0-2.0*beta3*betap0+2.0*betap2*beta1+2.0*beta2*betap1;
	rotp.Get0(0,2)= 2.0*betap3*beta1+2.0*beta3*betap1+2.0*betap2*beta0+2.0*beta2*betap0;
	rotp.Get0(1,0)= 2.0*betap3*beta0+2.0*beta3*betap0+2.0*betap2*beta1+2.0*beta2*betap1;
	rotp.Get0(1,1)= -4.0*beta3*betap3-4.0*beta1*betap1;
	rotp.Get0(1,2)= 2.0*betap3*beta2+2.0*beta3*betap2-2.0*betap1*beta0-2.0*beta1*betap0;
	rotp.Get0(2,0)= -2.0*betap2*beta0-2.0*beta2*betap0+2.0*betap3*beta1+2.0*beta3*betap1;
	rotp.Get0(2,1)= 2.0*betap3*beta2+2.0*beta3*betap2+2.0*betap1*beta0+2.0*beta1*betap0;
	rotp.Get0(2,2)= -4.0*beta2*betap2-4.0*beta1*betap1;
	return rotp;
}

Matrix3D ComputeRotMatrixEuler(const double phi1, const double phi2, const double phi3)
{
	//A(phi1):=<<cos(phi1),sin(phi1),0>|<-sin(phi1),cos(phi1),0>|<0,0,1>>;#AI1
	//A(phi2):=<<1,0,0>|<0,cos(phi2),sin(phi2)>|<0,-sin(phi2),cos(phi2)>>;#A12
	//A(phi3):=<<cos(phi3),sin(phi3),0>|<-sin(phi3),cos(phi3),0>|<0,0,1>>;#A23
	//Z-X-Z rotations
	Matrix3D rot; //rot = A(phi1)*A(phi2)*A(phi3); //Acc. to Shabana, Comput. Dynamics, 1994, p.363
	double sphi1=sin(phi1); double cphi1=cos(phi1);
	double sphi2=sin(phi2); double cphi2=cos(phi2);
	double sphi3=sin(phi3); double cphi3=cos(phi3);

	rot.Get0(0,0) = cphi1*cphi3-sphi1*cphi2*sphi3;
	rot.Get0(0,1) = -cphi1*sphi3-sphi1*cphi2*cphi3;
	rot.Get0(0,2) = sphi1*sphi2;
	rot.Get0(1,0) = sphi1*cphi3+cphi1*cphi2*sphi3;
	rot.Get0(1,1) = -sphi1*sphi3+cphi1*cphi2*cphi3;
	rot.Get0(1,2) = -cphi1*sphi2;
	rot.Get0(2,0) = sphi2*sphi3;
	rot.Get0(2,1) = sphi2*cphi3;
	rot.Get0(2,2) = cphi2;

	return rot;
}

Matrix3D ComputeRotMatrixWithKardanAngles(const double phi1, const double phi2, const double phi3) //this is the HOTINT Kardan angle definition from local to global coordinates
{
	//positive rotations around x,y and z-axis
	return RotMatrix3(phi3)*RotMatrix2(phi2)*RotMatrix1(phi1);
}

void RotMatToKardanAngles(const Matrix3D& A, Vector3D& phi) 
{
	//this function gives all angles between -pi to +pi such that the original rotation matrix can be recomputed
	//in cases of small angles about x or z, the functions ComputeRotMatrixWithKardanAngles and RotMatToKardanAngles lead to the same angles in the end
	//
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Kardan angles (rotation of global coordinate system about Z (gamma), then around local Y-axis about beta, than around X about alpha):
	// A rotates a vector in body fixed coordinates into global coordinate system
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// A_Kardan: see http://www-mechanik.uni-duisburg.de/Lehre/Lehrunterlagen/AdvancedDynamics/EulerBryantAngles.pdf
	//  		        [cos(gamma) cos(beta),  -sin(gamma) cos(alpha) + cos(gamma) sin(beta) sin(alpha),  sin(gamma) sin(alpha) + cos(gamma) sin(beta) cos(alpha)], 
	//  A_Kardan =  [sin(gamma) cos(beta),   cos(gamma) cos(alpha) + sin(gamma) sin(beta) sin(alpha), -cos(gamma) sin(alpha) + sin(gamma) sin(beta) cos(alpha)], 
  //              [-sin(beta)          ,   cos(beta) sin(alpha)                                   ,  cos(beta) cos(alpha)]
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	// atan2 is used instead of asin --> angles are correct in all 4 quadrants
	phi.Y() = asin(-A(3,1)); //beta € [-pi/2,+pi/2]
		
	if(fabs(fabs(A(3,1))-1.)> 16e-16)  // cos(beta)!=0
	{		
		phi.Z()=atan2(A(2,1), A(1,1));
		phi.X()=atan2(A(3,2), A(3,3));

		//determine correct sign:
		double sinphiZ = sin(phi.Z());
		if (fabs(sinphiZ) > 0.5)
		{
			phi.Y() = atan2(-A(3,1), A(2,1)/sinphiZ);
		}
		else
		{
			phi.Y() = atan2(-A(3,1), A(1,1)/cos(phi.Z()));
		}
	}
	else
	{
		//phi.Y() is approximately pi/2 or -pi/2, but uniquely defined between -pi to pi

		//cos(beta)==0, sin(beta) €{-1,1}, beta €{-pi/2,+pi/2}
		phi.X() = 0; //we do not need twice rotation around x
		
		//the rotation matrix in this case is approximately:
		//A=[[0, sin(beta)*sin(alpha), sin(beta)*cos(alpha)], 
		//   [0, cos(alpha), -sin(alpha)],
		//   [-sin(beta), 0, 0]]		
		//double fact = 1; 
		//if (-A(3,1) < 0) fact = -fact; //change sign, because in this case, phi.X() has opposite definition
	 // 
		//phi.X() = fact*atan2(A(1,2), A(1,3));
		phi.Z() = atan2(-A(1,2),A(2,2));
	}
}

void QuaternionsToKardanAngles(double beta0, double beta1, double beta2, double beta3, Vector3D& phi)
{  		
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Kardan angles (rotation of global coordinate system about Z (gamma), then around local Y-axis about beta, than around X about alpha):
	// A rotates a vector in body fixed coordinates into global coordinate system
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// A_Kardan: see http://www-mechanik.uni-duisburg.de/Lehre/Lehrunterlagen/AdvancedDynamics/EulerBryantAngles.pdf
	//  		        [cos(gamma) cos(beta),  -sin(gamma) cos(alpha) + cos(gamma) sin(beta) sin(alpha),  sin(gamma) sin(alpha) + cos(gamma) sin(beta) cos(alpha)], 
	//  A_Kardan =  [sin(gamma) cos(beta),   cos(gamma) cos(alpha) + sin(gamma) sin(beta) sin(alpha), -cos(gamma) sin(alpha) + sin(gamma) sin(beta) cos(alpha)], 
  //              [-sin(beta)          ,   cos(beta) sin(alpha)                                   ,  cos(beta) cos(alpha)]
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// A_EulerParameter: see Shabana, 1998, p. 32, equ. 2.14 a
	//                    [beta1^2+beta0^2-beta3^2-beta2^2,	2*beta1*beta2-2*beta0*beta3    , 2*beta1*beta3+2*beta0*beta2],
 	// A_EulerParameter = [2*beta1*beta2+2*beta0*beta3    , beta2^2-beta3^2+beta0^2-beta1^2, 2*beta2*beta3-2*beta0*beta1], 
	//                    [2*beta1*beta3-2*beta0*beta2, 2*beta2*beta3+2*beta0*beta1, beta3^2-beta2^2-beta1^2+beta0^2]
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//TODO:MSax: please replace and check
	//Matrix3D A = ComputeRotMatrixFromQuaternions(beta0, beta1, beta2, beta3);
	//RotMatToKardanAngles(A, phi);

	// atan2 is used instead of asin --> angles are correct in all 4 quadrants
	double A31 = 2.0*(beta1*beta3-beta0*beta2);
	phi.Y() = asin(-A31); //beta € [-pi/2,+pi/2]
		
	if(fabs(fabs(A31)-1.)> 16e-16)  // cos(beta)!=0
	{		
		double A11 = 1.0-2.0*(Sqr(beta3)+Sqr(beta2));
		double A21 = 2.0*(beta1*beta2+beta0*beta3);
		phi.Z()=atan2(A21, A11);

		double A32 = 2.0*beta3*beta2+2.0*beta1*beta0;
		double A33 = 1.0-2.0*(Sqr(beta2)+Sqr(beta1));
		phi.X()=atan2(A32, A33);

		//JG 2013-01-11: please check MSax
		//determine correct sign:
		double sinphiZ = sin(phi.Z());
		if (fabs(sinphiZ) > 0.5)
		{
			phi.Y() = atan2(-A31, A21/sinphiZ);
		}
		else
		{
			phi.Y() = atan2(-A31, A11/cos(phi.Z()));
		}


	}
	else
	{
		//phi.Y() is approximately pi/2 or -pi/2, but uniquely defined between -pi to pi

		//cos(beta)==0, sin(beta) €{-1,1}, beta €{-pi/2,+pi/2}
		phi.Z() = 0; //we do not need twice rotation around z 
		
		//A=[[0, sin(beta)*sin(alpha), sin(beta)*cos(alpha)], 
		//   [0, cos(alpha), -sin(alpha)],
		//   [-sin(beta), 0, 0]]		
	  
		double A12 = 2.0*(beta1*beta2-beta0*beta3);
		double A13 = 2.0*(beta3*beta1+beta2*beta0);
		phi.X() = atan2(A12, A13);
		for(int i=1;i<=20;i++)
		{
			(*global_uo) << "Warning: angles phi are not tested yet in this special case of beta € {-pi/2, pi/2}.\n";
		}
	}
}


Matrix3D ComputeRotMatrixEulerP(const double& phi1, const double& phi2, const double& phi3, 
																				 const double& phi1p, const double& phi2p, const double& phi3p)
{
	Matrix3D rotp;
	double sphi1=sin(phi1); double cphi1=cos(phi1);
	double sphi2=sin(phi2); double cphi2=cos(phi2);
	double sphi3=sin(phi3); double cphi3=cos(phi3);

	rotp.Get0(0,0)= -sphi1*phi1p*cphi3-cphi1*sphi3*phi3p-cphi1*phi1p*cphi2*sphi3+sphi1*sphi2*phi2p*sphi3-sphi1*cphi2*cphi3*phi3p;
	rotp.Get0(0,1)= sphi1*phi1p*sphi3-cphi1*cphi3*phi3p-cphi1*phi1p*cphi2*cphi3+sphi1*sphi2*phi2p*cphi3+sphi1*cphi2*sphi3*phi3p;
	rotp.Get0(0,2)= cphi1*phi1p*sphi2+sphi1*cphi2*phi2p;
	rotp.Get0(1,0)= cphi1*phi1p*cphi3-sphi1*sphi3*phi3p-sphi1*phi1p*cphi2*sphi3-cphi1*sphi2*phi2p*sphi3+cphi1*cphi2*cphi3*phi3p;
	rotp.Get0(1,1)= -cphi1*phi1p*sphi3-sphi1*cphi3*phi3p-sphi1*phi1p*cphi2*cphi3-cphi1*sphi2*phi2p*cphi3-cphi1*cphi2*sphi3*phi3p;
	rotp.Get0(1,2)= sphi1*phi1p*sphi2-cphi1*cphi2*phi2p;
	rotp.Get0(2,0)= cphi2*phi2p*sphi3+sphi2*cphi3*phi3p;
	rotp.Get0(2,1)= cphi2*phi2p*cphi3-sphi2*sphi3*phi3p;
	rotp.Get0(2,2)= -sphi2*phi2p;

	return rotp;
}

void QuaternionsToEulerAngles(double beta0, double beta1, double beta2, double beta3, Vector3D& phi)
{
	//compute directly the Euler angles from the quaternions:

	//Rotation Matrix -> Euler angles (phi-Z,theta-X,psi-Z):
	//phi = arctan(A13, -A23);
	//theta = arccos(A33);
	//psi = arctan(A31, A32);
	// -->if theta approx 0
	// psi = 0; -->no twice rotation for Z
	// theta = arccos(A33); -->leave as is
	// phi = arctan(-A12, A11);
	//Quaternions -> Euler angles (phi-Z,theta-X,psi-Z):
	//phi   = arctan(q1*q3+q2*q4, -(q2*q3-q1*q4)); //problems e.g. if q1=q2=q4=0 --> arctan(0,0)
	//theta = arccos(-q1*q1-q2*q2, q3*q3+q4*q4);
	//psi   = arctan(q1*q3-q2*q4, q2*q3+q1*q4);

	double A33 = -2.0*beta2*beta2-2.0*beta1*beta1+1.0;
	phi.Y() = acos(A33);
	if (fabs(phi.Y()) > 16e-16)
	{
		double A13 = 2.0*beta3*beta1+2.0*beta2*beta0;
		double A23 = 2.0*beta3*beta2-2.0*beta1*beta0;
		double A31 = -2.0*beta2*beta0+2.0*beta3*beta1;
		double A32 = 2.0*beta3*beta2+2.0*beta1*beta0;
		phi.X() = atan2(A13, -A23);
		phi.Z() = atan2(A31, A32);

		//determine correct sign:
		double sinpsi = sin(phi.Z());
		if (fabs(sinpsi) > 0.5)
		{
			phi.Y() = atan2(A13/sinpsi, A33);
		}
		else
		{
			phi.Y() = atan2(A23/cos(phi.Z()), A33);
		}
	} 
	else
	{
		double A11 = -2.0*beta3*beta3-2.0*beta2*beta2+1.0;
		double A12 = -2.0*beta3*beta0+2.0*beta2*beta1;
		phi.Z() = 0;
		phi.X() = atan2(-A12, A11);
	}
}

void RotMatToEulerAngles(const Matrix3D& A, Vector3D& phi)
{
	//Rotation Matrix -> Euler angles (phi-Z,theta-X,psi-Z):

	double A33 = A(3,3);
	phi.Y() = acos(A33);
	if (fabs(phi.Y()) > 16e-16)
	{
		double A13 = A(1,3);
		double A23 = A(2,3);
		double A31 = A(3,1);
		double A32 = A(3,2);


		phi.X() = atan2(A13, -A23);
		phi.Z() = atan2(A31, A32);

		//determine correct sign:
		double sinpsi = sin(phi.Z());
		if (fabs(sinpsi) > 0.5)
		{
			phi.Y() = atan2(A31/sinpsi, A33);
		}
		else
		{
			phi.Y() = atan2(A32/cos(phi.Z()), A33);
		}
	} 
	else
	{
		double A11 = A(1,1);
		double A12 = A(1,2);
		phi.Z() = 0;
		phi.X() = atan2(-A12, A11);
	}
}


void QuaternionsPToAngularVelocities(double beta0, double beta1, double beta2, double beta3, 
																							double beta0p, double beta1p, double beta2p, double beta3p, 
																			 Vector3D& phip)
{
	//only for initialization, not for computation!!!
	beta0 *= 2; beta1 *= 2; beta2 *= 2; beta3 *= 2;
	//the betas are already multiplied with 2, compared to Shabana
	Matrix G(3,4);	//G=2E
	G(1,1) = -beta1; G(1,2) =  beta0; G(1,3) = -beta3; G(1,4) =  beta2;
	G(2,1) = -beta2; G(2,2) =  beta3; G(2,3) =  beta0; G(2,4) = -beta1;
	G(3,1) = -beta3; G(3,2) = -beta2; G(3,3) =  beta1; G(3,4) =  beta0;

	Vector q(beta0p, beta1p, beta2p, beta3p); 
	Vector f = G*q;
	phip.X() = f(1);
	phip.Y() = f(2);
	phip.Z() = f(3);
}

void RotMatToQuaternions(const Matrix3D& A, 
																			 double& beta0, double& beta1, double& beta2, double& beta3)
{

	double trace = A(1,1) + A(2,2) + A(3,3) + 1.0;
	double M_EPSILON = 1e-15;

	if( fabs(trace) > M_EPSILON ) 
	{
		double s = 0.5 / sqrt(fabs(trace));
		beta0 = 0.25 / s;
		beta1 = ( A(3,2) - A(2,3) ) * s;
		beta2 = ( A(1,3) - A(3,1) ) * s;
		beta3 = ( A(2,1) - A(1,2) ) * s;
	} 
	else 
	{
		if ( A(1,1) > A(2,2) && A(1,1) > A(3,3) ) {
			double s = 2.0 * sqrt(fabs(1.0 + A(1,1) - A(2,2) - A(3,3)));
			beta1 = 0.25 * s;
			beta2 = (A(1,2) + A(2,1) ) / s;
			beta3 = (A(1,3) + A(3,1) ) / s;
			beta0 = (A(2,3) - A(3,2) ) / s;

		} else if (A(2,2) > A(3,3)) {
			double s = 2.0 * sqrt(fabs(1.0 + A(2,2) - A(1,1) - A(3,3)));
			beta1 = (A(1,2) + A(2,1) ) / s;
			beta2 = 0.25 * s;
			beta3 = (A(2,3) + A(3,2) ) / s;
			beta0 = (A(1,3) - A(3,1) ) / s;
		} else {
			double s = 2.0 * sqrt(fabs(1.0 + A(3,3) - A(1,1) - A(2,2) ));
			beta1 = (A(1,3) + A(3,1) ) / s;
			beta2 = (A(2,3) + A(3,2) ) / s;
			beta3 = 0.25 * s;
			beta0 = (A(1,2) - A(2,1) ) / s;
		}
	}

}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// functions for DH-Transformation
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// indices Aij similar to Bremer H. - definitions
// computes the 4x4 denavit hartenberg transformation matrix from given dh-parameters
Matrix DH_ij(double a, double alpha, double d, double theta)
{
	///////////////////////////////////////////////////////////////////////////////////////////
	//                      4x4 matrix transformation	                                       //
	//     (see Kinematic Analysis of Robot Manipulators, Crane&Duffy,p.5ff)	  		         //
	//     given: point P	                                                                   //
	//		         vectors: rj_P=rj_j0->P																										 //
	//                              ri_P = [DH_ij]. rj_P , ri_P=ri_i0->P                     //
	//                                                                                       //
	///////////////////////////////////////////////////////////////////////////////////////////
	//                                                                                       //
	// vectors to a point P                                                                  //
	// rI_P =  rI_1 + AI1*r1_P                                                               //
	//      =  DH_I1 * r1_P                                                                  //
	//      =  DH_I1*DH_12*r2_P                         = DH_I2 * r2_P                       //
	//      =  ...																																					 //
	//                                -----                                                  //
	//      =  DH_I1*DH_12 * ... * ( |DH_ij| ) * rj_P   = DH_Ij * rj_P                       //
	//                                -----                                                  //
	//  ri_P = DH_ij.rj_P = [  Aij       | ri_j ]																		     		 //
	//                      [------------------ ] . rj_P																		 //
	//                      [ 0   0   0  |   1  ]																						 //
	//																																											 //
	///////////////////////////////////////////////////////////////////////////////////////////
	// 
	//    given: r*_P=r*_*o->P, * ... any coordinate system
	//    
	//-------------------------------------------------------------------------------------------------------------------------------------
	//      I = diagonalMatrix(1,1,1); e1_vec,e3_vec ... Einheitsvektoren
	//        ri_jo  +   rj_P                   r1_P=r1_1o->P
	//      o------->o------->P      =        o--------------->P
	//      1o       jo                      1o
	//          
	//              
	//      ri_jo=[a,0,0]^t=a*e1_vec
	//   
	//    r1_P = T1j . rj_P = [  I    |a*e1_vec]   |rj_P|
	//                        [--------------- ] . |----| = rj_P + a*e1*1	  ... coordinate system 1 has same basis vectors like c.s. j (parallel)
	//                        [ 0 0 0 |  1     ]   |  1 |                       it's origin is translated by a along the xj-axis
	//-------------------------------------------------------------------------------------------------------------------------------------
	//
	//    r2_P = R21 . r1_P = [  Ax(-alpha) |  0 ]   |r1_P|
	//                        [----------------- ] . |----| = A(-alpha).r1_P ... coordinate system 2 is a rotated coord. sys. around the x-axis
	//                        [ 0   0   0  |  1  ]   |  1 |                
	//------------------------------------------------------------------------------------------------------------------------------------- 
	//    r3_P = T32 . r2_P = ...                           = r2_P + d*e3    ...        translation along zi
	//-------------------------------------------------------------------------------------------------------------------------------------
	//    ri_P = Ri3 . r3_P = [  Az(-theta) |  0  ]  |r3_P|
	//                        [-------------------] .|----| = A(-theta).r3_P ...        rotation around zi
	//                        [ 0   0   0  |  1   ]  |  1 |                 
	//-------------------------------------------------------------------------------------------------------------------------------------
	Matrix Ri3(4,4);		 
	Ri3(1,1) = cos(theta); Ri3(1,2) =-sin(theta); Ri3(1,3) = 0;          Ri3(1,4) = 0; 
	Ri3(2,1) = sin(theta); Ri3(2,2) = cos(theta); Ri3(2,3) = 0;          Ri3(2,4) = 0; 
	Ri3(3,1) = 0;          Ri3(3,2) = 0;          Ri3(3,3) = 1;          Ri3(3,4) = 0; 
	Ri3(4,1) = 0;          Ri3(4,2) = 0;          Ri3(4,3) = 0;          Ri3(4,4) = 1; 
	Matrix T32(4,4);
	T32(1,1) = 1;          T32(1,2) = 0;          T32(1,3) = 0;          T32(1,4) = 0;
	T32(2,1) = 0;          T32(2,2) = 1;          T32(2,3) = 0;          T32(2,4) = 0; 
	T32(3,1) = 0;          T32(3,2) = 0;          T32(3,3) = 1;          T32(3,4) = d; 
	T32(4,1) = 0;          T32(4,2) = 0;          T32(4,3) = 0;          T32(4,4) = 1; 
	Matrix R21(4,4);
	R21(1,1) = 1;          R21(1,2) = 0;					R21(1,3) = 0;          R21(1,4) = 0; 
	R21(2,1) = 0;          R21(2,2) = cos(alpha); R21(2,3) =-sin(alpha); R21(2,4) = 0; 
	R21(3,1) = 0;          R21(3,2) = sin(alpha); R21(3,3) = cos(alpha); R21(3,4) = 0;
	R21(4,1) = 0;          R21(4,2) = 0;          R21(4,3) = 0;          R21(4,4) = 1; 
	Matrix T1j(4,4);	
	T1j(1,1) = 1; T1j(1,2) = 0; T1j(1,3) = 0; T1j(1,4) = a;
	T1j(2,1) = 0; T1j(2,2) = 1; T1j(2,3) = 0; T1j(2,4) = 0;
	T1j(3,1) = 0; T1j(3,2) = 0; T1j(3,3) = 1; T1j(3,4) = 0;
	T1j(4,1) = 0; T1j(4,2) = 0; T1j(4,3) = 0; T1j(4,4) = 1;

	Matrix DH_ij(4,4);
	DH_ij=Ri3*T32*R21*T1j;
	// component of DH-matrix :          DH_ij=Ri3*T32*R21*T1j
	// DH_ij = [cos(theta), -sin(theta)*cos(alpha), sin(theta)*sin(alpha) , cos(theta)*a]
	//         [sin(theta), cos(theta)*cos(alpha),  -cos(theta)*sin(alpha), sin(theta)*a]
	//         [    0     ,       sin(alpha),            cos(alpha)       ,      d      ] 
	//         [    0     ,           0,                     0,                  1      ]  
	return DH_ij;
}

void GetAij(const Matrix* T4_ij,Matrix3D& Aij) //T4_ij...4x4 matrix, Aij...rotation matrix
{
	for(int i=1;i<=3;i++)
	{	
		for(int j=1;j<=3;j++)	
		{		
			Aij(i,j)=(*T4_ij)(i,j);		
 		}	
	}  // get Aij from 4x4 matrix
}

void Get_ri_jo(const Matrix* T4_ij,Vector3D& ri_jo)    //T_ij...4x4 matrix, ri_jo...from io to jo; components of vector in i-system
{
	for(int i=1 ;i<=3;i++) 
	{
		ri_jo(i) = (*T4_ij)(i,4); // last column
	}
}

void GetT4x4_ij(const Matrix3D& Aij,const Vector3D& ri_jo, Matrix* T4_ij)
{
	//		ri_P = DH_ij.rj_P = [  Aij       | ri_j0]		[ rj_P ]																					 
	//										    [------------------ ] . [------]																				 
	//										    [ 0   0   0  |   1  ]		[   1  ]																					  
	for(int i=1;i<=3;i++)
	{	
		for(int j=1;j<=3;j++)	
		{		
			(*T4_ij)(i,j) = Aij(i,j);		
		}	
	}  // copy Aij into 4x4-matrix

	for(int i=1 ;i<=3;i++) 
	{
		(*T4_ij)(4,i) = 0.;      // last row
		(*T4_ij)(i,4) = ri_jo(i); // last column
	}

	(*T4_ij)(4,4) = 1.;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// compute intersection of two lines in vector paramater representation 
// d(L1,L2) = |vector between points of minimal distance|
// L1 = pt1+s*dir1,   L2 = pt2+t*dir2
// with a = dir1*dir1, b = dir1*dir2, c = dir2*dir2, d = dir1*(pt1-pt2), e = dir2*(pt1-pt2)
// s_mindist = (be-cd)/(ac-bb), t_mindist = (ae-bd)/(ac-bb)
// d(L1,L2) = |(pt1-pt2)+s*dir1-t*dir2|
// 
// returns 1 on success 
// returns 2 for identical lines
// retirns 0 for parallel or skew lines(pti is not changed)
int IntersectLines3D( Vector3D pt1, Vector3D dir1, Vector3D pt2, Vector3D dir2, Vector3D& pti, Vector2D& lambda_mu  )
{
	int rv = 0;
	double eps = 1e-10;
  pti = pt1;
	double s = 0.,t = 0., mindist = 0.;

	Vector3D w0 = pt1-pt2;
	double a = dir1*dir1;
	double b = dir1*dir2;
	double c = dir2*dir2;
	double d = dir1*w0;
	double e = dir2*w0;

  double denom = a*c-b*b;
	if( abs(denom) < eps ) 
	{
// parallel or identical lines
		t = d/b; 
		Vector3D wc = pt1-pt2+s*dir1-t*dir2;
		mindist = wc.Norm();
		if( mindist < eps )
		{
			rv = 2; // identical 
			pti = pt1;
			lambda_mu.X() = 0;
			lambda_mu.Y() = t;		
		}
		else
		{
			rv = 0; // parallel
			pti = pt1;
			lambda_mu.X() = 0;
			lambda_mu.Y() = t;		
		}
	}
	else
	{
// intersecting or skew lines
		s = (b*e-c*d)/denom;
		t = (a*e-b*d)/denom;
		Vector3D wc = pt1-pt2+s*dir1-t*dir2;
		mindist = wc.Norm();
		if( mindist < eps )
		{
		  rv = 1; // intersection
			pti = pt1+s*dir1;
			lambda_mu.X() = s;
			lambda_mu.Y() = t;		
		}
		else
		{
			rv = 0; // skew 
			pti = pt1+s*dir1;
			lambda_mu.X() = s;
			lambda_mu.Y() = t;		
		}
	}

  return rv;
}



// compute intersection of two lines in vector paramater representation
// returns 1 on success 
// returns 2 for identical lines
// retirns 0 for parallel or skew lines(pti is not changed)
int IntersectLines2D( Vector2D pt1, Vector2D dir1, Vector2D pt2, Vector2D dir2, Vector2D& pti, Vector2D& lambda_mu )
{
	int rv = 0;
	double eps = 1e-10;
  pti = pt1;
	double s = 0.,t = 0., mindist = 0.;

	Vector2D w0 = pt1-pt2;
	double a = dir1*dir1;
	double b = dir1*dir2;
	double c = dir2*dir2;
	double d = dir1*w0;
	double e = dir2*w0;

  double denom = a*c-b*b;
	if( abs(denom) < eps ) 
	{
// parallel or identical lines
		t = d/b; 
		Vector2D wc = pt1-pt2+s*dir1-t*dir2;
		mindist = wc.Norm();
		if( mindist < eps )
		{
			rv = 2; // identical 
			pti = pt1;
			lambda_mu.X() = 0;
			lambda_mu.Y() = t;		
		}
		else
		{
			rv = 0; // parallel
			pti = pt1;
			lambda_mu.X() = 0;
			lambda_mu.Y() = t;		
		}
	}
	else
	{
// intersecting or skew lines
		s = (b*e-c*d)/denom;
		t = (a*e-b*d)/denom;
		Vector2D wc = pt1-pt2+s*dir1-t*dir2;
		mindist = wc.Norm();
		if( mindist < eps )
		{
		  rv = 1; // intersection
			pti = pt1+s*dir1;
			lambda_mu.X() = s;
			lambda_mu.Y() = t;		
		}
		else
		{
			assert(0); //there is nw skew in 2D
			rv = 2; // skew 
		}
	}

  return rv;
}


//dont use this function!, use KArdan functions instead!
//this is not the HOTINT Kardan definition:
Matrix3D ComputeRot1Rot2Rot3Matrix(const double phi1, const double phi2, const double phi3)
{
	//positive rotations around x,y and z-axis
	//note that rotation matrix for rot around y has different sign in 'sin' term!
	//see linalg3d.cpp definitions of RotMatrix1(phi), etc.
  //A(phi1):=<<1.,0.,0.>|<0.,cos(phi1),sin(phi1)>|<0.,-sin(phi1), cos(phi1)>>;#=AI1
  //A(phi2):=<<cos(phi2),0.,-sin(phi2)>|<0.,1.,0.>|<sin(phi2),0., cos(phi2)>>;#=A12
  //A(phi3):=<<cos(phi3),sin(phi3),0.>|<-sin(phi3), cos(phi3),0.>|<0.,0.,1.>>;#=A23
	return RotMatrix1(phi1)*RotMatrix2(phi2)*RotMatrix3(phi3);
}

//dont use this function!, use KArdan functions instead!
void QuaternionsToRot1Rot2Rot3Angles(double beta0, double beta1, double beta2, double beta3, Vector3D& phi)
{
	//-------------------------------------------------------------------------------
	//Rotation Matrix -> angles: phi1-X,phi2-Y,phi3-Z:
	//-------------------------------------------------------------------------------
  //A(phi1):=<<1.,0.,0.>|<0.,cos(phi1),sin(phi1)>|<0.,-sin(phi1), cos(phi1)>>;#=AI1
  //A(phi2):=<<cos(phi2),0.,-sin(phi2)>|<0.,1.,0.>|<sin(phi2),0., cos(phi2)>>;#=A12
  //A(phi3):=<<cos(phi3),sin(phi3),0.>|<-sin(phi3), cos(phi3),0.>|<0.,0.,1.>>;#=A23
	//
	//A(phi)=[cos(phi2)*cos(phi3), -cos(phi2)*sin(phi3), sin(phi2)],
	//	     [sin(phi1)*sin(phi2)*cos(phi3)+cos(phi1)*sin(phi3), -sin(phi1)*sin(phi2)*sin(phi3)+cos(phi1)*cos(phi3), -sin(phi1)*cos(phi2)], 
	//			 [-cos(phi1)*sin(phi2)*cos(phi3)+sin(phi1)*sin(phi3), cos(phi1)*sin(phi2)*sin(phi3)+sin(phi1)*cos(phi3), cos(phi1)*cos(phi2)]])
	//-------------------------------------------------------------------------------
	//A(beta)...Shabana, 1998, p. 32, equ. 2.14 a	(Euler parameter)
	//-------------------------------------------------------------------------------
	//phi2 = arcsin(A13)
	//phi2 = arcsin(Rot[1,3]); 
	//
	//if cos(phi2) != 0 //phi2 != Pi/2
	//phi3 = arcsin(-Rot[1,2]/cos(phi2)); 
	//phi1 = arcsin(-Rot[2,3]/cos(phi2)); 
  //
	//else
	//phi3 = 0; //we do not need twice rotation around z
	//phi1 = arctan(Rot[3,2],Rot[2,2]);
  //-------------------------------------------------------------------------------

	
	double A13 = 2.0*beta3*beta1+2.0*beta2*beta0;
	phi.Y() = asin(A13);

	if (fabs(fabs(A13)-1.) > 16e-16)
	{
		double A23 = 2.0*beta3*beta2-2.0*beta1*beta0;
		double A12 = -2.0*beta3*beta0+2.0*beta2*beta1;
		phi.Z() = asin(-A12/cos(phi.Y())); 
		phi.X() = asin(-A23/cos(phi.Y())); 
	} 
	else
	{
		double A22 = -2.0*beta3*beta3-2.0*beta1*beta1+1.0;
		double A32 = 2.0*beta3*beta2+2.0*beta1*beta0;
		phi.Z() = 0;
		phi.X() = atan2(A32, A22);
	}
}



