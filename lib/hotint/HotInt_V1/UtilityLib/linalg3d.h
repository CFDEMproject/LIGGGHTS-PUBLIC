//#**************************************************************
//#
//# filename:             linalg3d.h
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

void Normal3D(const Vector3D& p1,const Vector3D& p2,const Vector3D& p3, Vector3D& n); //normalized normal

//distance from point p to infinite line (line consists of two points lp1 and lp2)
double DistToLine(const Vector3D& lp1, const Vector3D& lp2, const Vector3D& p);
double DistToLine(const Vector2D& lp1, const Vector2D& lp2, const Vector2D& p);

//distance from point p to line defined by linepoint lp and line vector lv; pp is nearest (projected) point on line:
//lv must have finite length!
double DistToLine(const Vector3D& lp, const Vector3D& lv, const Vector3D& p, Vector3D& pp);

//intersect two planes (plane normals pn1, pn2; plane points pp1, pp2): result=Line (represented by point lp and line direction lv)
//compute the intersection of two planes (if an intersection return 1, if planes are equal return 2, otherwise return 0)
int PlaneIntersection(const Vector3D& pn1, const Vector3D& pp1, 
											 const Vector3D& pn2, const Vector3D& pp2, 
											 Vector3D& lp, Vector3D& lv);
//{
//P=[0 0 0];
//N=cross(N1,N2);
//
//%  test if the two planes are parallel
//if norm(N) < 10^-7                % Plane 1 and Plane 2 are near parallel
//    V=A1-A2;
//        if (dot(N1,V) == 0)         
//            check=1;                    % Plane 1 and Plane 2 coincide
//           return
//        else 
//            check=0;                   %Plane 1 and Plane 2 are disjoint
//            return
//        end
//end
//            
// check=2;
// 
//% Plane 1 and Plane 2 intersect in a line
//%first determine max abs coordinate of cross product
//maxc=find(abs(N)==max(abs(N)));
//
//
//%next, to get a point on the intersection line and
//%zero the max coord, and solve for the other two
//      
//d1 = -dot(N1,A1);   %the constants in the Plane 1 equations
//d2 = -dot(N2, A2);  %the constants in the Plane 2 equations
//
//switch maxc
//    case 1                   % intersect with x=0
//        P(1)= 0;
//        P(2) = (d2*N1(3) - d1*N2(3))/ N(1);
//        P(3) = (d1*N2(2) - d2*N1(2))/ N(1);
//    case 2                    %intersect with y=0
//        P(1) = (d1*N2(3) - d2*N1(3))/ N(2);
//        P(2) = 0;
//        P(3) = (d2*N1(1) - d1*N2(1))/ N(2);
//    case 3                    %intersect with z=0
//        P(1) = (d2*N1(2) - d1*N2(2))/ N(3);
//        P(2) = (d1*N2(1) - d2*N1(1))/ N(3);
//        P(3) = 0;
//end
//
//}

//distance from line, including boundary:
double MinDistLP(const Vector3D& lp1, const Vector3D& lp2, const Vector3D& p);

//distance from line, including boundary, get projected point:
double MinDistLP(const Vector3D& lp1, const Vector3D& lp2, const Vector3D& p, Vector3D& pp);

double Deflection(const Vector3D& lp1, const Vector3D& lp2, const Vector3D& p, const Vector3D& pn);

//get deflection in 2D, distance to line segment, right side is positive
double RighthandMinDist2D(const Vector2D& lp1, const Vector2D& lp2, const Vector2D& p, Vector2D& pp);

//Distance of 2 points:
double Dist(const Vector3D& p1, const Vector3D& p2);
double Dist(const Vector2D& p1, const Vector2D& p2);

//squared Distance:
double Dist2(const Vector3D& p1, const Vector3D& p2);
double Dist2(const Vector2D& p1, const Vector2D& p2);

double NormalizedVectorAngle(const Vector3D& v1, const Vector3D& v2);
double VectorAngle(Vector3D v1, Vector3D v2);

double NormalizedVectorAngle(const Vector2D& v1, const Vector2D& v2);
double VectorAngle(Vector2D v1, Vector2D v2);
double PolarAngle(double real, double imaginary);
double PolarAngle(Vector2D& v);

// returns parameters k and d of linear equation y = k*x+d, defined by 2 points pi=(xi,yi)
Vector2D GetLinearEquationParameters(double x1,double y1,double x2,double y2);	//$ DR 2012-04-13

//Minimal Distance of p to plane (p0,n0), normalized normal, projected point in plane pp:
double DistToPlane(const Vector3D& p0, const Vector3D& n0, const Vector3D& p, Vector3D& pp);

//Minimal Distance of p to plane(HNF)
double DistToPlane(const Vector3D& p, const Vector3D& nplane, double cplane);
double DistToPlane(const Vector2D& p, const Vector2D& nplane, double cplane);
double DistToPlane(const Vector3D& p, const Vector3D& nplane, const Vector3D& n0);
double DistToPlane(const Vector2D& p, const Vector2D& nplane, const Vector2D& n0);

//Mirror Point/Vactor at plane(HNF)
Vector3D MirrorAtPlane3D(const Vector3D& p, const Vector3D& nplane, double cplane);

//Minimal Distance to Triangle:
double DistToTrig(const Vector3D& p1, const Vector3D& p2, const Vector3D& p3, 
									const Vector3D& p);
//Minimal Distance to Quad:
double DistToQuad(const Vector3D& p1, const Vector3D& p2, const Vector3D& p3, const Vector3D& p4, 
									const Vector3D& p);

//compute local coordinates lam1, lam2 of a vector v in a triangle:
void LocalCoordinates (const Vector3D & e1, const Vector3D & e2,
		  const Vector3D & v, double & lam1, double & lam2);

//minimum distance between point p and triangle
double MinDistTP(const Vector3D & tp1, const Vector3D & tp2, 
		   const Vector3D & tp3, const Vector3D & p);

//minimum distance between point p and triangle
double MinDistTP(const Vector3D& tp1, const Vector3D& tp2, 
		   const Vector3D& tp3, const Vector3D& p, Vector3D& pp);

void ProjectInPlane(const Vector3D& p1, const Vector3D& p2, const Vector3D& p3, Vector3D& pp);
void ProjectInPlane(const Vector3D& p1, const Vector3D& n, Vector3D& pp);

//return intersection point; rv=1 if success, rv=0 if no success, means line is normal to plane normal; pline is projected in direction of vline in plane
int IntersectLineWithPlane(const Vector3D& pplane, const Vector3D& nplane, const Vector3D& vline, Vector3D& pline);

//int PointInside(const Vector3D& p1, const Vector3D& p2, const Vector3D& p3, const Vector3D& pp);
//double MinDistTP2(const Vector3D& p1, const Vector3D& p2, const Vector3D& p3, const Vector3D& p3d);

double NormSym(const Matrix3D& m);
Matrix3D Deviator(const Matrix3D& m);
double DoubleProd3D(const Matrix3D& m, const Matrix3D& n);

double Area(const Vector3D& p1, const Vector3D& p2, const Vector3D& p3);
double Area(const Vector3D& p1, const Vector3D& p2, const Vector3D& p3, const Vector3D& p4);


////// compute intersection of two lines in vector paramater representation
////// returns 1 on success , -1 on skew lines(pti is not changed)
////// returns 0 when identical(returns pti == pt1), -1 when parallel(pti is not changed)
////int IntersectionPoint( Vector3D pt1, Vector3D dir1, Vector3D pt2, Vector3D dir2, Vector3D& pti );

// compute intersection of two lines in vector paramater representation
// returns 1 on success 
// returns 2 for identical lines
// retirns 0 for parallel or skew lines(pti is not changed)
int IntersectLines3D( Vector3D pt1, Vector3D dir1, Vector3D pt2, Vector3D dir2, Vector3D& pti, Vector2D& lambda_mu );

// compute intersection of two lines in vector paramater representation
// returns 1 on success 
// returns 2 for identical lines
// retirns 0 for parallel or skew lines(pti is not changed)
int IntersectLines2D( Vector2D pt1, Vector2D dir1, Vector2D pt2, Vector2D dir2, Vector2D& pti, Vector2D& lambda_mu );

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Rotation parameters:

//compute rotation matrix from Quaternions (euler parameters); b0 ..b3 are the euler parameters
Matrix3D ComputeRotMatrixFromQuaternions(double b0, double b1, double b2, double b3);


Matrix3D ComputeRotMatrixEulerParam(const double& beta0, const double& beta1,
																	 const double& beta2, const double& beta3);

Matrix3D ComputeRotMatrixPEulerParam(const double& beta0, const double& beta1,
																		const double& beta2, const double& beta3, const double& betap0, const double& betap1,
																		const double& betap2, const double& betap3);

Matrix3D ComputeRotMatrixEuler(const double phi1, const double phi2, const double phi3);

Matrix3D ComputeRotMatrixWithKardanAngles(const double phi1, const double phi2, const double phi3); //this is the HOTINT Kardan angle definition from local to global coordinates

void RotMatToKardanAngles(const Matrix3D& A, Vector3D& phi); //this function gives all angles between -pi to +pi

void QuaternionsToKardanAngles(double beta0, double beta1, double beta2, double beta3, Vector3D& phi);//RL

Matrix3D ComputeRotMatrixEulerP(const double& phi1, const double& phi2, const double& phi3, 
																				 const double& phi1p, const double& phi2p, const double& phi3p);

void QuaternionsToEulerAngles(double beta0, double beta1, double beta2, double beta3, Vector3D& phi);

void RotMatToEulerAngles(const Matrix3D& A, Vector3D& phi);

void QuaternionsPToAngularVelocities(double beta0, double beta1, double beta2, double beta3, 
																							double beta0p, double beta1p, double beta2p, double beta3p, Vector3D& phip);

void RotMatToQuaternions(const Matrix3D& A, double& beta0, double& beta1, double& beta2, double& beta3);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Matrix DH_ij(double a, double alpha, double d, double theta);              // compute denavit hartenberg transformaion matrix from dh-parameters
void GetAij(const Matrix* T4_ij,Matrix3D& Aij);                            // compute rotation matrix from 4x4 transformation matrix
void Get_ri_jo(const Matrix* T4_ij,Vector3D& ri_jo);                       // compute displacement vector from 4x4 transformation matrix
void GetT4x4_ij(const Matrix3D& Aij,const Vector3D& ri_jo, Matrix* T4_ij); // compute 4x4 transformation matrix from rotation matrix and displacement vector


//+++++++++++++++++++++++++++++++++++++++++

//compute rotation matrices with respect to axis 1, 2 and 3
Matrix3D RotMatrix1(double p);
Matrix3D RotMatrix2(double p);
Matrix3D RotMatrix3(double p);

double GetRotation(int rotaxis, const Matrix3D& rotmat);

class Vector3D BaseVector3D(int i);
class Vector2D BaseVector2D(int i);
//-----------------------------------------------------------------
//--------------    CLASS VECTOR 2D   -----------------------------
//-----------------------------------------------------------------
//Generation of a vector with user-defined length
//Operations: +,-,* (also matrices), <<
//initialized with 0 !!!
class Vector2D
{
public:

	//Constructors, Destructor
	Vector2D() { vec[0]=0; vec[1]=0;};
	Vector2D(double x) { vec[0]=x; vec[1]=x; };
	Vector2D(double x, double y)
	{
		vec[0]=x; vec[1]=y;
	}
	Vector2D(const Vector2D& veci)
	{ 
		vec[0] = veci[0]; vec[1] = veci[1];
	}
//	Vector2D( Vector3D& veci);
  Vector3D MakeV3D() const;
	//virtual ~Vector2D() { };

	const double* GetVecPtr() const {return &vec[0];}
	double* GetVecPtr() {return &vec[0];}

	//Assignment-operator
	Vector2D& operator= (const Vector2D& veci)
	{ 
		if (&veci == this) return *this;

		vec[0] = veci[0];
		vec[1] = veci[1];
		return *this;
	}

	Vector2D& operator= (const double& val)
	{ 
		vec[0] = val;
		vec[1] = val;
		return *this;
	}
	//Logical operator
	int operator== (const Vector2D& veci)
	{
		return (vec[0] == veci[0] && 
			vec[1] == veci[1]);
	}

	double& X() {return vec[0];}
	double& Y() {return vec[1];}

	const double& X() const {return vec[0];}
	const double& Y() const {return vec[1];}

	//Referencing access-operator
	double& operator[](int elem)
	{
#ifndef __QUICKMATH
		release_assert((elem>=0) && (elem<2));
#endif
		return vec[elem];
	};
	//Referencing access-operator
	const double& operator[](int elem) const
	{
#ifndef __QUICKMATH
		release_assert((elem>=0) && (elem<2));
#endif
		return vec[elem];
	};
	//Referencing access-operator
	double& operator()(int elem)
	{
#ifndef __QUICKMATH
		release_assert((elem>0) && (elem<=2));
#endif
		return vec[elem-1];
	};
	//Referencing access-operator
	const double& operator()(int elem) const
	{
#ifndef __QUICKMATH
		release_assert((elem>0) && (elem<=2));
#endif
		return vec[elem-1];
	};
	void Set(const double& x, const double& y)
	{
		vec[0]=x; vec[1]=y; 
	}
	void Get(double& x, double& y) //$JG2012-01: this is for easier access to vector3d data
	{
		x=vec[0]; y=vec[1]; 
	}

	//Returns the length of a vector
	int Length() const {return 2; };
	int GetLen() const {return 2; };

	double Cross(const Vector2D& v) const
	{
		return vec[0]*v.vec[1]-vec[1]*v.vec[0];
	}
	//Returns the quadratic norm of a vector
	double Norm() const {return sqrt(vec[0]*vec[0]+vec[1]*vec[1]);}
	//Returns the quadratic norm of a vector squared
	double Norm2() const {return vec[0]*vec[0]+vec[1]*vec[1];}

	//norms the vector
	void Normalize() 
	{
		double l=Norm();
		if (l>1e-16) {vec[0]/=l;vec[1]/=l;}
	}

	//Returns the maximum-norm of a vector
	double MaxNorm() const {return Maximum(vec[0],vec[1]);}

	//Returns the minimum-norm of a vector
	double MinNorm() const {return Minimum(vec[0],vec[1]);};

	void Scale(double x, double y)
	{
		vec[0]/=x; vec[1]/=y;
	}

	//Arithmetic operations with one parameter
	Vector2D operator+= (const Vector2D& v) 
	{
		vec[0]+=v.vec[0];
		vec[1]+=v.vec[1];
		return *this;
	}
	Vector2D operator-= (const Vector2D& v) 
	{
		vec[0]-=v.vec[0];
		vec[1]-=v.vec[1];
		return *this;
	}

	Vector2D operator*= (const double& v)
	{
		vec[0]*=v;
		vec[1]*=v;
		return *this;
	}

	Vector2D operator/= (const double& v)
	{
		vec[0]/=v;
		vec[1]/=v;
		return *this;
	}

	//Arithmetic operations with two parameters
	friend Vector2D operator+ (const Vector2D& v1, const Vector2D& v2)
	{
		return Vector2D(v1.vec[0]+v2.vec[0],
			v1.vec[1]+v2.vec[1]);
	}
	friend Vector2D operator- (const Vector2D& v1, const Vector2D& v2)
	{
		return Vector2D(v1.vec[0]-v2.vec[0],
			v1.vec[1]-v2.vec[1]);
	}
	friend Vector2D operator* (const Vector2D& v, const double& val)
	{
		return Vector2D(v.vec[0]*val,
			v.vec[1]*val);
	}

	friend Vector2D operator* (const double& val, const Vector2D& v)
	{
		return Vector2D(v.vec[0]*val,
			v.vec[1]*val);
	}
	friend double operator* (const Vector2D& v1, const Vector2D& v2)
	{
		return v1.vec[0]*v2.vec[0]+v1.vec[1]*v2.vec[1];
	}

	friend Vector2D operator* (const Matrix& m, const Vector2D& v);
	friend Vector2D operator* (const Matrix3D& m, const Vector2D& v);
	friend Vector2D operator* (const Matrix2D& m, const Vector2D& v);
	friend Vector2D operator* (const Vector2D& v, const Matrix2D& m); //res=v^T*m =m^T*v 
	friend Vector2D operator* (const Matrix2D& m, const Vector& v);
	friend void Mult(const Matrix2D& m, const Vector2D& v, Vector& res);
	friend void Mult(const Matrix& m, const Vector2D& v, Vector& res); //computes res=m*v

	//Output parameter
	friend ostream& operator<<(ostream& os, const Vector2D& v) 
	{
		os << "[" << v.X() << ", " << v.Y() << "]";
		return os;
	}

private:
	double vec[2];

};


//-----------------------------------------------------------------
//--------------    CLASS VECTOR 3D   -----------------------------
//-----------------------------------------------------------------
//Generation of a vector with user-defined length
//Operations: +,-,* (also matrices), <<
//initialized with 0 !!!
class Vector3D
{
public:

	//Constructors, Destructor
	Vector3D() { vec[0]=0; vec[1]=0; vec[2]=0;};
	Vector3D(double x) { vec[0]=x; vec[1]=x; vec[2]=x;};
	Vector3D(double x, double y, double z)
	{
		vec[0]=x; vec[1]=y; vec[2]=z;
	}
	Vector3D(const Vector3D& avec)
	{
		vec[0] = avec.X(); vec[1] = avec.Y(); vec[2] = avec.Z();
	}
	
//	Vector3D(const Vector2D& avec);
	Vector2D MakeV2D() const;

//	virtual ~Vector3D() { };

	//Assignment-operator
	Vector3D& operator= (const Vector3D& veci)
	{ 
		if (&veci == this) return *this;

		vec[0] = veci[0];
		vec[1] = veci[1];
		vec[2] = veci[2];
		return *this;
	}

	Vector3D& operator= (const double& val)
	{ 
		vec[0] = val;
		vec[1] = val;
		vec[2] = val;
		return *this;
	}
	//Logical operator
	int operator== (const Vector3D& veci)
	{
		return (vec[0] == veci[0] && 
			vec[1] == veci[1] &&
			vec[2] == veci[2]);
	}

	int operator== (const Vector3D& veci) const
	{
		return (vec[0] == veci[0] && 
			vec[1] == veci[1] &&
			vec[2] == veci[2]);
	}

	double& X() {return vec[0];}
	double& Y() {return vec[1];}
	double& Z() {return vec[2];}

	const double& X() const {return vec[0];}
	const double& Y() const {return vec[1];}
	const double& Z() const {return vec[2];}

	//Referencing access-operator
	double& operator[](int elem)
	{
#ifndef __QUICKMATH
		release_assert((elem>=0) && (elem<3));
#endif
		return vec[elem];
	};
	//Referencing access-operator
	const double& operator[](int elem) const
	{
#ifndef __QUICKMATH
		release_assert((elem>=0) && (elem<3));
#endif
		return vec[elem];
	};
	//Referencing access-operator
	double& operator()(int elem)
	{
#ifndef __QUICKMATH
		release_assert((elem>0) && (elem<=3));
#endif
		return vec[elem-1];
	};
	//Referencing access-operator
	const double& operator()(int elem) const
	{
#ifndef __QUICKMATH
		release_assert((elem>0) && (elem<=3));
#endif
		return vec[elem-1];
	};

	const double* GetVecPtr() const {return &vec[0];}
	double* GetVecPtr() {return &vec[0];}
	//Returns the length of a vector
	int Length() const {return 3; };
	int GetLen() const {return 3; };
	void Set(const double& x, const double& y, const double& z)
	{
		vec[0]=x; vec[1]=y; vec[2]=z; 
	}
	void Get(double& x, double& y, double& z) //$JG2012-01: this is for easier access to vector3d data
	{
		x=vec[0]; y=vec[1]; z=vec[2]; 
	}

	//Returns the quadratic norm of a vector
	double Norm() const {return sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);}
	//Returns the quadratic norm of a vector squared
	double Norm2() const {return vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2];}

	//norms the vector
	void Normalize() 
	{
		double l=Norm();
		if (l!=0.) {vec[0]/=l;vec[1]/=l;vec[2]/=l;}
	}
	void Scale(double x, double y, double z)
	{
		vec[0]/=x; vec[1]/=y; vec[2]/=z; 
	}
	//normalizes *this and gets some orthogonal vectors n1 and n2
 	void SetNormalBasis(Vector3D& n1, Vector3D& n2)
	{
		Normalize();
		Vector3D nx;
		if (fabs(vec[0]) > 0.5 && fabs(vec[1]) < 0.1 && fabs(vec[2]) < 0.1) nx.Set(0.,1.,0.);
		else nx.Set(1.,0.,0.);

		double h = nx*(*this);
		n1 = nx-h*(*this);
		n1.Normalize();
		n2=this->Cross(n1);
	}

	//Project n into normal plane of *this
	void GramSchmidt(Vector3D& n)
	{
		double h = n*(*this)/((*this)*(*this));
		n -= h*(*this);
		n.Normalize();
	}

	//Project n into normal plane of *this
	//PG: added const-version of this routine
	void GramSchmidt(Vector3D& n) const
	{
		double h = n*(*this)/((*this)*(*this));
		n -= h*(*this);
		n.Normalize();
	}

	//Returns the maximum-norm of a vector
	double MaxNorm() const {return Maximum(vec[0],vec[1],vec[2]);}

	//Returns the minimum-norm of a vector
	double MinNorm() const {return Minimum(vec[0],vec[1],vec[2]);};

	//Returns the normalized planar perpendicular vector of a vector
	Vector3D Cross(const Vector3D& v) const
	{
		return Vector3D(vec[1]*v.vec[2]-vec[2]*v.vec[1], 
			vec[2]*v.vec[0]-vec[0]*v.vec[2], 
			vec[0]*v.vec[1]-vec[1]*v.vec[0]);
	}

	mystr& MakeString(mystr intro = mystr(""));

	//Arithmetic operations with one parameter
	Vector3D operator+= (const Vector3D& v) 
	{
		vec[0]+=v.vec[0];
		vec[1]+=v.vec[1];
		vec[2]+=v.vec[2];
		return *this;
	}
	Vector3D operator-= (const Vector3D& v) 
	{
		vec[0]-=v.vec[0];
		vec[1]-=v.vec[1];
		vec[2]-=v.vec[2];
		return *this;
	}

	Vector3D operator*= (const double& v)
	{
		vec[0]*=v;
		vec[1]*=v;
		vec[2]*=v;
		return *this;
	}

	Vector3D operator/= (const double& v)
	{
		vec[0]/=v;
		vec[1]/=v;
		vec[2]/=v;
		return *this;
	}

	//Arithmetic operations with two parameters
	friend Vector3D operator+ (const Vector3D& v1, const Vector3D& v2)
	{
		return Vector3D(v1.vec[0]+v2.vec[0],
			v1.vec[1]+v2.vec[1],
			v1.vec[2]+v2.vec[2]);
	}
	friend Vector3D operator- (const Vector3D& v1, const Vector3D& v2)
	{
		return Vector3D(v1.vec[0]-v2.vec[0],
			v1.vec[1]-v2.vec[1],
			v1.vec[2]-v2.vec[2]);
	}
	friend Vector3D operator* (const Vector3D& v, const double& val)
	{
		return Vector3D(v.vec[0]*val,
			v.vec[1]*val,
			v.vec[2]*val);
	}

	friend Vector3D operator* (const double& val, const Vector3D& v)
	{
		return Vector3D(v.vec[0]*val,
			v.vec[1]*val,
			v.vec[2]*val);
	}
	friend double operator* (const Vector3D& v1, const Vector3D& v2)
	{
		return v1.vec[0]*v2.vec[0]+v1.vec[1]*v2.vec[1]+v1.vec[2]*v2.vec[2];
	}

	friend Vector3D operator* (const Matrix& m, const Vector3D& v);
	friend Vector3D operator* (const Matrix3D& m, const Vector3D& v);
	friend Vector3D operator* (const Vector3D& v, const Matrix3D& m); //res=v^T*m =m^T*v 
	friend Vector3D operator* (const Matrix3D& m, const Vector& v);
	friend Vector3D operator* (const Matrix3D& m, const Vector& v);

	friend void Mult(const Matrix3D& m, const Vector3D& v, Vector& res);
	friend void Mult(const Matrix& m, const Vector3D& v, Vector& res); //computes res=m*v
	friend void Mult(const MatrixXD& m, const Vector& v, Vector3D& res);
	friend void Mult(const Matrix3D& m, const Vector3D& v, Vector3D& res);
	friend void Mult(const Matrix& m, const Vector& v, Vector3D& res); //computes res=m*v
	friend void MultTp(const Matrix& m, const Vector& v, Vector3D& res); //computes res=m*v

	friend void StrainVectorToMatrix2D(Matrix2D& m, const Vector3D& v);

	//Output parameter
	friend ostream& operator<<(ostream& os, const Vector3D& v) 
	{
		os << "[" << v.X() << ", " << v.Y() << ", " << v.Z() << "]";
		return os;
	}

private:
	double vec[3];

};


//-----------------------------------------------------------------
//CLASS MatrixXD  CLASS MatrixXD  CLASS MatrixXD  CLASS MatrixXD  
//-----------------------------------------------------------------
//Generation of a Matrix with user-defined size
//Operations: +,-,*,<<
class MatrixXD
{
protected:
	int rows;
	int cols;
	double mat[16];
public:

	//Constructors,destructor
	MatrixXD()
	{
		rows = 0; cols = 0;
	}

	//fill Diagonal matrix with x
	MatrixXD(double x)
	{
		rows = 0; cols = 0;
	}

	// copy constructor
	MatrixXD(const MatrixXD& mati)
	{
		cols = mati.cols; rows = mati.rows;

		for (int i=0; i < cols*rows; i++)
		{
			mat[i] = mati.mat[i];
		}
	}

	virtual ~MatrixXD() {};

	virtual int Getcols() const {return cols;}
	virtual int Getrows() const {return rows;}
	virtual double* GetMatPtr() {return mat;}
	
	//Assignment-operator
	virtual MatrixXD& operator= (const MatrixXD& mati)
	{
		rows=mati.rows;
		cols=mati.cols;
		for (int i=0; i < rows*cols; i++)
		{
			mat[i] = mati.mat[i];
		}
		return *this;
	}

	virtual void SetAll(double x) 
	{
		for (int i=0; i < rows*cols; i++)
		{
			mat[i] = x;
		}
	}

	virtual void SetSize(int r, int c) {rows = r; cols = c;}

	//fill Diagonal matrix with x -- pure virtual: overwrite in derived class
	virtual void SetDiag(double x) = 0;

	//returnsum of absolute values
	virtual double AbsNorm() const 
	{
		double val = 0.;

		for (int i=0; i < rows*cols; i++)
		{
			val += fabs(mat[i]);
		}

		return val;
	}

	//Returns Determinate of matrix
	virtual double Det() const
	{
		if (cols == 3)
		{
			return mat[0]*mat[4]*mat[8]-mat[0]*mat[5]*mat[7]-mat[3]*mat[1]*mat[8]
			+mat[3]*mat[2]*mat[7]+mat[6]*mat[1]*mat[5]-mat[6]*mat[2]*mat[4];
		}
		else
		{
			return mat[0]*mat[3]-mat[1]*mat[2];
		}
	}
	
	virtual double Trace() const
	{
		double sum=0;
		for (int i=0; i < min(rows,cols); i++)
		{
			sum += Get0(i,i);  
		}
		return sum;
	}
	
	// synonim of DoubleContract
	virtual double InnerProduct(const MatrixXD& rhs) const
	{
#ifndef __QUICKMATH
		release_assert(rows == rhs.rows);
		release_assert(cols == rhs.rows);
#endif
		double val = 0.;
		for(int i = 0; i < rows*cols; i++)
		{
			val += mat[i]*rhs.mat[i];
		}
		return val;
	}

	//transpose me
	virtual void TpYs()
	{
#ifndef __QUICKMATH
		release_assert(rows == cols);
#endif
		for (int i=1; i<rows; i++)
		{
			for (int j=0; j<i; j++) 
			{
				swap(mat[j*rows+i],mat[i*cols+j]);
			}
		}
	}

	virtual void MakeSym() //make matrix symmetric: (A+A^T)/2
	{
#ifndef __QUICKMATH
		release_assert(rows==cols);
#endif
		for (int i=0; i < cols; i++)
		{
			for (int j=0; j < i; j++)
			{
				Get0(i,j) += Get0(j,i);
				Get0(i,j) *= 0.5;
				Get0(j,i) = Get0(i,j);
			}
		}
	}

	// synonim of InnerProduct
	virtual double DoubleContract(const MatrixXD& m) const
	{
		return InnerProduct(m);
	}

	//Referencing access-operator on element using row- and column-values, 0-based!!!
	virtual double& Get0(int row, int col)
	{
		return mat[row*cols+col];
	};

	//Referencing constant access-operator on element using row- and column-values, 0-based!!!
	virtual const double& Get0(int row, int col) const
	{
		return mat[row*cols+col];
	};


	//Referencing access-operator on element using row- and column-values, 1-based!!!
	virtual double& operator()(int row, int col)
	{
		return mat[(row-1)*cols+col-1];
	};

	//Referencing constant access-operator on element using row- and column-values, 1-based!!!
	virtual const double& operator()(int row, int col) const
	{
		return mat[(row-1)*cols+col-1];
	};


	virtual void Transpose()
	{
#ifndef __QUICKMATH
		release_assert(rows==cols);
#endif
		double help;
		int idx1;
		int idx2;
		
		for(int i=0; i < rows; i++)
		{
			for(int j=0; j < i; j++)
			{
				idx1 = i*cols + j;
				idx2 = j*cols + i;
		
				help = mat[idx1];
				mat[idx1] = mat[idx2];
				mat[idx2] = help;
			}
		}

	}

	mystr& MakeString(mystr intro = mystr(""));

	friend double InnerProduct(const MatrixXD& m1, const MatrixXD& m2) { return m1.InnerProduct(m2); }
	friend void Mult(const MatrixXD& m, const Vector& v, Vector& res);
	friend void Mult(const MatrixXD& m1, const Matrix& m2, Matrix& res);
	friend void Mult(const Matrix& m1, const MatrixXD& m2, Matrix& res);

	//writes out MatrixXD with constant width and factor normalized
	friend ostream& operator<<(ostream& os, const MatrixXD& m);
};



//-----------------------------------------------------------------
//CLASS Matrix2D  CLASS Matrix2D  CLASS Matrix2D  CLASS Matrix2D  
//-----------------------------------------------------------------

class Matrix2D : public MatrixXD
{
public:
	
	Matrix2D()
	{
		rows = 2; cols = 2;
		
		mat[0] = 0.; mat[1] = 0.;
		mat[2] = 0.; mat[3] = 0.;
	}

	Matrix2D(double x1)
	{
		rows = 2; cols = 2;
		
		mat[0] = x1; mat[1] = 0.;
		mat[2] = 0.; mat[3] = x1;
	}
	
	Matrix2D(double x1, double x2)
	{
		rows=2; cols=2;
		
		mat[0] = x1; mat[1] = 0.;
		mat[2] = 0.; mat[3] = x2;
	}
	
	Matrix2D(
		double x11, double x12,
		double x21, double x22)
	{
		rows=2; cols=2;
		
		mat[0] = x11; mat[1] = x12;
		mat[2] = x21; mat[3] = x22;
	}
	
	void SetDiag(double x)
	{
		rows = 2; cols = 2;

		mat[0]=x; mat[1]=0; 
		mat[2]=0; mat[3]=x; 
	}

	//set Matrix as Skew matrix (to rotate a Vector2D around an angle phi)
	void SetSkew(double phi)
	{
		rows = 2; cols = 2;
	
		double cos_phi = cos(phi);
		double sin_phi = sin(phi);

		mat[0] = cos_phi; mat[1] = -sin_phi;
		mat[2] = sin_phi; mat[3] = cos_phi;
	}
	
	void Set(const Vector2D& r1, const Vector2D& r2)
	{
		mat[0] = r1.X(); mat[1] = r2.X();
		mat[2] = r1.Y(); mat[3] = r2.Y();
	}

	//only for 2x2, fast!!!
	void GetATA2(MatrixXD& m)
	{
		for (int i=0;i<2;i++)
		{
			for (int j=i;j<2;j++)
			{
				m.Get0(i,j) = 0.5*(mat[i]*mat[j] + mat[2+i]*mat[2+j]);
			}
		}
		// 0 1
		// 2 3
		m.Get0(1,0) = m.Get0(0,1);
	}

	Matrix2D GetTp() const
	{
		Matrix2D mt;
		mt.SetSize(cols, rows);

		for (int i=0; i<rows; i++)
		{
			for (int j=0; j<cols; j++) 
			{ 
				mt.Get0(i,j)=Get0(j,i);
			}
		}
		return mt;
	}

	int GetInverse(Matrix2D& m2);
	int Invert()
	{
		Matrix2D m2;
		int rv = GetInverse(m2);
		
		if (rv)
		{
			*this = m2;
		}

		return rv;
	}

	//Arithmetic operations with 1 parameter
	Matrix2D& operator+= (const Matrix2D& m)
	{
#ifndef __QUICKMATH
		release_assert(rows == m.Getrows());
		release_assert(cols == m.Getcols());
#endif
		for (int i=0; i < rows; i++)
		{
			for (int j=0; j < cols; j++)
			{
				Get0(i,j) += m.Get0(i,j);
			}
		}
		return *this;
	}

	Matrix2D& operator-= (const Matrix2D& m)
	{
#ifndef __QUICKMATH
		release_assert(rows == m.Getrows());
		release_assert(cols == m.Getcols());
#endif
		for (int i=0; i < rows; i++)
		{
			for (int j=0; j < cols; j++)
			{
				Get0(i,j) -= m.Get0(i,j);
			}
		}
		return *this;
	}

	Matrix2D& operator*= (const double& val)
	{
		for (int i=0; i < rows; i++)
		{
			for (int j=0; j < cols; j++)
			{
				Get0(i,j) *= val;
			}
		}
		return *this;
	}

	//get a orthogonal matrix from vectors v1 and v2
	void SetOrthogonalBasis(const Vector2D& v1, const Vector2D& v2)
	{
		Vector2D r1 = v1;
		Vector2D r2 = v2;
		r1.Normalize();

		double h = r2*r1;
		r2 -= h*r1;
		r2.Normalize();

		Set(r1,r2);
	}

	//Arithmetic operations with 2 parameters
	friend Matrix2D operator+ (const Matrix2D& m1, const Matrix2D& m2);
	friend Matrix2D operator- (const Matrix2D& m1, const Matrix2D& m2);
	friend Matrix2D operator* (const Matrix2D& m1, const double& val);
	friend Matrix2D operator* (const double& val, const Matrix2D& m1);
	friend Matrix2D operator* (const Matrix2D& m1, const Matrix2D& m2);
	friend Vector2D operator* (const Matrix2D& m, const Vector2D& v);
	friend Vector2D operator* (const Vector2D& v, const Matrix2D& m); //res=v^T*m =m^T*v 
	friend Vector2D operator* (const Matrix2D& m, const Vector& v);

	friend void Mult(const Matrix2D& m, const Vector2D& v, Vector& res);

	//transform strain and stress vectors to matrices and vice versa (by PG)
	friend void StrainVectorToMatrix2D(Matrix2D& m, const Vector& v);
	friend void StrainVectorToMatrix2D(Matrix2D& m, const Vector3D& v);
	friend void StressVectorToMatrix2D(Matrix2D& m, const Vector& v);
	friend void Matrix2DToStrainVector(Vector& v, const Matrix2D& m);
	friend void Matrix2DToStressVector(Vector& v, const Matrix2D& m);
};




//-----------------------------------------------------------------
//CLASS Matrix3D  CLASS Matrix3D  CLASS Matrix3D  CLASS Matrix3D  
//-----------------------------------------------------------------

class Matrix3D : public MatrixXD
{
public:
	
	Matrix3D()
	{
		rows = 3; cols = 3;
		
		mat[0] = 0.; mat[1] = 0.; mat[2] = 0.;
		mat[3] = 0.; mat[4] = 0.; mat[5] = 0.;
		mat[6] = 0.; mat[7] = 0.; mat[8] = 0.;
	}

	Matrix3D(double x1)
	{
		rows = 3; cols = 3;
		
		mat[0] = x1; mat[1] = 0.; mat[2] = 0.;
		mat[3] = 0.; mat[4] = x1; mat[5] = 0.;
		mat[6] = 0.; mat[7] = 0.; mat[8] = x1;
	}
	
	Matrix3D(double x1, double x2, double x3)
	{
		rows = 3; cols = 3;
		
		mat[0] = x1; mat[1] = 0.; mat[2] = 0.;
		mat[3] = 0.; mat[4] = x2; mat[5] = 0.;
		mat[6] = 0.; mat[7] = 0.; mat[8] = x3;
	}
	
	Matrix3D(
		double x11, double x12, double x13,
		double x21, double x22, double x23,
		double x31, double x32, double x33)
	{	
		rows = 3; cols = 3;
		
		mat[0] = x11; mat[1] = x12; mat[2] = x13; 
		mat[3] = x21; mat[4] = x22; mat[5] = x23; 
		mat[6] = x31; mat[7] = x32; mat[8] = x33; 
	}

	Matrix3D(
		double x11, double x12, double x13, double x14,
		double x21, double x22, double x23, double x24,
		double x31, double x32, double x33, double x34)
	{		
		rows=3; cols=4;

		mat[0]=x11; mat[1]=x12; mat[2]=x13; mat[3]=x14; 
		mat[4]=x21; mat[5]=x22; mat[6]=x23; mat[7]=x24; 
		mat[8]=x31; mat[9]=x32; mat[10]=x33;mat[11]=x34; 
	}
	
	//only possible for 
	Matrix3D(const Matrix& mat);

	void SetDiag(double x)
	{
		rows = 3; cols = 3;

		mat[0]=x; mat[1]=0; mat[2]=0; 
		mat[3]=0; mat[4]=x; mat[5]=0; 
		mat[6]=0; mat[7]=0; mat[8]=x; 
	}

	//set Matrix as Skew matrix
	void SetSkew(double x1, double x2, double x3)
	{
		rows = 3; cols = 3;

		mat[0]=  0; mat[1]=-x3; mat[2]= x2; 
		mat[3]= x3; mat[4]=  0; mat[5]=-x1; 
		mat[6]=-x2; mat[7]= x1; mat[8]= 0; 
	}

	//set Matrix as Skew matrix
	void SetSkew(const Vector3D& x)
	{
		rows = 3; cols = 3;

		mat[0]=     0; mat[1]=-x.Z(); mat[2]= x.Y(); 
		mat[3]= x.Z(); mat[4]=     0; mat[5]=-x.X(); 
		mat[6]=-x.Y(); mat[7]= x.X(); mat[8]=     0; 
	}

	void Set43(
		double x11, double x12, double x13, 
		double x21, double x22, double x23,
		double x31, double x32, double x33,
		double x41, double x42, double x43)
	{
		rows = 4; cols = 3;

		mat[0]=x11; mat[1]=x12; mat[2]=x13; 
		mat[3]=x21; mat[4]=x22; mat[5]=x23;
		mat[6]=x31; mat[7]=x32; mat[8]=x33;
		mat[9]=x41;mat[10]=x42;mat[11]=x43; 
	}

	//PG: is this redundant since introduction of class Matrix2D (?)
	void Set22(
		double x11, double x12, 
		double x21, double x22)
	{
		rows = 2; cols = 2;

		mat[0]=x11; mat[1]=x12;  
		mat[2]=x21; mat[3]=x22; 
	}

	void Set(const Vector3D& r1, const Vector3D& r2, const Vector3D& r3)
	{
		mat[0] = r1.X(); mat[1] = r2.X(); mat[2] = r3.X();
		mat[3] = r1.Y(); mat[4] = r2.Y(); mat[5] = r3.Y();
		mat[6] = r1.Z(); mat[7] = r2.Z(); mat[8] = r3.Z();
	}

	double Mises() const
	{
		Matrix3D s;
		if (cols == 3)
		{
			s = *this - (1./3.)*Trace()*Matrix3D(1.);
		}
		else if (cols == 2)
		{
			s = *this;
			double tr = (1./3.)*Trace();
			s(1,1) -= tr; s(2,2) -= tr;
		}
		return sqrt(3./2.*(s.InnerProduct(s)));
	}

	//only for 3x3, fast!!!
	void GetATA2(MatrixXD& m)
	{
		for (int i=0; i < 3; i++)
		{
			for (int j=i; j < 3; j++)
			{
				m.Get0(i,j) = 0.5*(mat[i]*mat[j] + mat[3+i]*mat[3+j] + mat[2*3+i]*mat[2*3+j]);
			}
		}
		// 0 1 2
		// 3 4 5
		// 6 7 8
		m.Get0(1,0) = m.Get0(0,1);
		m.Get0(2,1) = m.Get0(1,2);
		m.Get0(2,0) = m.Get0(0,2);
	}

	Matrix3D GetTp() const
	{
		Matrix3D mt;
		mt.SetSize(cols, rows);

		for (int i=0; i<rows; i++)
		{
			for (int j=0; j<cols; j++) 
			{ 
				mt.Get0(j,i)=Get0(i,j);
			}
		}
		return mt;
	}

	int GetInverse(Matrix3D& m2);
	int Invert()
	{
		Matrix3D m2;
		int rv = GetInverse(m2);
		
		if (rv)
		{
			*this = m2;
		}

		return rv;
	}

	//Arithmetic operations with 1 parameter
	Matrix3D& operator+= (const Matrix3D& m)
	{
#ifndef __QUICKMATH
		release_assert(rows == m.Getrows());
		release_assert(cols == m.Getcols());
#endif
		for (int i=0; i < rows; i++)
		{
			for (int j=0; j < cols; j++)
			{
				Get0(i,j) += m.Get0(i,j);
			}
		}
		return *this;
	}

	Matrix3D& operator-= (const Matrix3D& m)
	{
#ifndef __QUICKMATH
		release_assert(rows == m.Getrows());
		release_assert(cols == m.Getcols());
#endif
		for (int i=0; i < rows; i++)
		{
			for (int j=0; j < cols; j++)
			{
				Get0(i,j) -= m.Get0(i,j);
			}
		}
		return *this;
	}

	Matrix3D& operator*= (const double& val)
	{
		for (int i=0; i < rows; i++)
		{
			for (int j=0; j < cols; j++)
			{
				Get0(i,j) *= val;
			}
		}
		return *this;
	}

	//get a orthogonal matrix from vectors v1 and v2
	void SetOrthogonalBasis(const Vector3D& v1, const Vector3D& v2)
	{
		Vector3D r1 = v1;
		Vector3D r2 = v2;
		r1.Normalize();

		double h = r2*r1;
		r2 -= h*r1;
		r2.Normalize();

		Set(r1,r2,r1.Cross(r2));
	}


	//Arithmetic operations with 2 parameters
	friend Matrix3D operator+ (const Matrix3D& m1, const Matrix3D& m2);
	friend Matrix3D operator- (const Matrix3D& m1, const Matrix3D& m2);
	friend Matrix3D operator* (const Matrix3D& m1, const double& val);
	friend Matrix3D operator* (const double& val, const Matrix3D& m1);
	friend Matrix3D operator* (const Matrix3D& m1, const Matrix3D& m2);
	friend Vector3D operator* (const Matrix3D& m, const Vector3D& v);
	friend Vector2D operator* (const Matrix3D& m, const Vector2D& v);	//this routine is redundant (since introduction of class Matrix2D)
	friend Vector3D operator* (const Vector3D& v, const Matrix3D& m); //res=v^T*m =m^T*v 
	friend Vector3D operator* (const Matrix3D& m, const Vector& v);

	friend void Mult(const Matrix3D& m, const Vector3D& v, Vector& res);
	friend void Mult(const Matrix3D& m, const Vector3D& v, Vector3D& res);
	
	//transform strain and stress vectors to matrices and vice versa (by PG)
	friend void StrainVectorToMatrix3D(Matrix3D& m, const Vector& v);
	friend void StressVectorToMatrix3D(Matrix3D& m, const Vector& v);
	friend void Matrix3DToStrainVector(Vector& v, const Matrix3D& m);
	friend void Matrix3DToStressVector(Vector& v, const Matrix3D& m);
};














//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++    BOX3D     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Box3D
{
public:
	Box3D()
	{
		//set empty box:
		pmin = Vector3D(1.,1.,1.);
		pmax = Vector3D(-1.,-1.,-1.);
	}
	Box3D(const Vector3D& p1, const Vector3D& p2)
	{
		pmin.X() = Minimum(p1.X(),p2.X());    
		pmin.Y() = Minimum(p1.Y(),p2.Y());    
		pmin.Z() = Minimum(p1.Z(),p2.Z());    

		pmax.X() = Maximum(p1.X(),p2.X());    
		pmax.Y() = Maximum(p1.Y(),p2.Y());    
		pmax.Z() = Maximum(p1.Z(),p2.Z());    
	}
	Box3D(const Box3D& b)
	{
		pmin = b.pmin;
		pmax = b.pmax;
	}
	Box3D(const Vector3D& c, double r)
	{
		pmin = c;
		pmax = c;
		Increase(r);
	}
	int Empty() const 
	{
		if (pmin.X() == 1 && pmax.X() == -1) {return 1;}
		return 0;
	}
	void Clear()
	{
		//set empty box:
		pmin = Vector3D(1.,1.,1.);
		pmax = Vector3D(-1.,-1.,-1.);
	}
	void Add(const Vector3D & p)
	{
		if (Empty())
		{
			pmin = p;
			pmax = p;
		}
		else
		{
			pmin.X() = Minimum(pmin.X(),p.X());    
			pmin.Y() = Minimum(pmin.Y(),p.Y());    
			pmin.Z() = Minimum(pmin.Z(),p.Z());    

			pmax.X() = Maximum(pmax.X(),p.X());    
			pmax.Y() = Maximum(pmax.Y(),p.Y());    
			pmax.Z() = Maximum(pmax.Z(),p.Z()); 
		}
	}
	void Add(const Box3D& b)
	{
		if (b.Empty()) return;

		if (Empty())
		{
			pmin = b.PMin();
			pmax = b.PMax();
		}
		else
		{
			pmin.X() = Minimum(pmin.X(),b.pmin.X());    
			pmin.Y() = Minimum(pmin.Y(),b.pmin.Y());    
			pmin.Z() = Minimum(pmin.Z(),b.pmin.Z());    

			pmax.X() = Maximum(pmax.X(),b.pmax.X());    
			pmax.Y() = Maximum(pmax.Y(),b.pmax.Y());    
			pmax.Z() = Maximum(pmax.Z(),b.pmax.Z()); 
		}
	}
	const Vector3D& PMin() const {return pmin;}
	const Vector3D& PMax() const {return pmax;}
	const double SizeX() const {return pmax[0]-pmin[0];}
	const double SizeY() const {return pmax[1]-pmin[1];}
	const double SizeZ() const {return pmax[2]-pmin[2];}
	Vector3D Center() const {return 0.5*(pmin+pmax);}
	double Radius() const {return 0.5*(pmax-pmin).Norm();}

	void Increase(double x, double y, double z)
	{
		pmin.X() -= x;
		pmin.Y() -= y;
		pmin.Z() -= z;
		pmax.X() += x;
		pmax.Y() += y;
		pmax.Z() += z;
	}
	void Increase(double x) 
	{
		Increase(x,x,x);
	}
  void InflateFactor(double x)
	{
		Vector3D pmi = PMin();
		Vector3D pma = PMax();
		Vector3D pc = Center();
		pmin = pc + ( pmi-pc ) * x;
		pmax = pc + ( pma-pc ) * x;
	}

	int Intersect(const Box3D& b) const 
	{
		if ( pmin[0] > b.pmax[0] || pmax[0] < b.pmin[0]
			|| pmin[1] > b.pmax[1] || pmax[1] < b.pmin[1]
			|| pmin[2] > b.pmax[2] || pmax[2] < b.pmin[2])
			return 0;

		return 1;
	}

	/// return 1 if point p in closure
	int IsIn (const Vector3D & p) const
	{
		if ( pmin[0] <= p[0] && pmax[0] >= p[0]
			&& pmin[1] <= p[1] && pmax[1] >= p[1]
			&& pmin[2] <= p[2] && pmax[2] >= p[2])
			return 1;

		return 0;
	}

private:
	Vector3D pmin, pmax;
};



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Zero-based!!!!
//generate a searchtree with a physical size of Box3D b
//all items must fit into this box and should be equally distributed
//then use AddItems to add items with a bounding box and an identifier
//GetItemsInBox() gives you all identifiers which have a bounding box within the specified box
class SearchTree
{
public:
	SearchTree():data(0), sx(0), sy(0), sz(0) {};
	SearchTree(int sizex, int sizey, int sizez, Box3D b)
	{
		box = b;
		if (box.SizeX()*box.SizeY()*box.SizeZ() == 0) box.Increase(1e-4);

		sx = sizex;
		sy = sizey;
		sz = sizez;

		data = new IVector[sx*sy*sz]();
		for (int i = 0; i < sx*sy*sz; i++)
		{
			data[i].Init();
		}
	}
	SearchTree(const SearchTree& tree)
	{
		box = tree.box;

		sx = tree.sx;
		sy = tree.sy;
		sz = tree.sz;

		data = new IVector[sx*sy*sz]();
		for (int i = 0; i < sx*sy*sz; i++)
		{
			data[i].Init();
			data[i] = tree.data[i];
		}
	}

	SearchTree& operator=(const SearchTree& tree)
	{
		if (&tree == this) return *this;

		if (data)
		{
			for (int i = 0; i < sx*sy*sz; i++)
			{
				data[i].Flush();
			}
			delete [] data;
		}

		box = tree.box;

		sx = tree.sx;
		sy = tree.sy;
		sz = tree.sz;

		data = new IVector[sx*sy*sz]();
		for (int i = 0; i < sx*sy*sz; i++)
		{
			data[i] = tree.data[i];
		}
		return *this;
	}
	
	virtual ~SearchTree()
	{
		if (data)
		{
			for (int i = 0; i < sx*sy*sz; i++)
			{
				data[i].Flush();
			}
			delete [] data;
		}
	}

	void ResetSearchTree(int sizex, int sizey, int sizez, Box3D b)
	{
		ClearItems(); //empty all items
		box = b;
		if (box.SizeX()*box.SizeY()*box.SizeZ() == 0) box.Increase(1e-4);

		if (sx != sizex || sy != sizey || sz != sizez)
		{
			if (data)
			{
				for (int i = 0; i < sx*sy*sz; i++)
				{
					data[i].Flush();
				}
				delete [] data;
			}

			sx = sizex;
			sy = sizey;
			sz = sizez;

			data = new IVector[sx*sy*sz]();
			for (int i = 0; i < sx*sy*sz; i++)
			{
				data[i].Init();
			}
		}
	}

	//return 6 indices for box: minx, maxx, miny, maxy, minz, maxz
	void GetBoxIndizes(const Box3D& b, int* ind) const;

	//return items in box defined by 6 indices: minx, maxx, miny, maxy, minz, maxz
	void GetItemsInBox(int* ind, TArray<int>& items) const;

	//add items in box defined by 6 indices: minx, maxx, miny, maxy, minz, maxz
	//does not reset items list
	void AddItemsInBox(int* ind, TArray<int>& items) const;

	//get only items in box
	void GetItemsInBox(const Box3D& b, TArray<int>& items) const;

	//get items in box, do not reset items list
	void AddItemsInBox(const Box3D& b, TArray<int>& items) const;

	void AddItem(const Box3D& b, int identifier);

	//return i-index for a double i-value in Box
	int IndX(double x) const;

	//return i-index for a double i-value in Box
	int IndY(double x) const;

	//return i-index for a double i-value in Box
	int IndZ(double x) const;

	int Index(int x, int y, int z) const;

	void ClearItems()
	{
		for (int ix=0; ix < sx; ix++)
		{
			for (int iy=0; iy < sy; iy++)
			{
				for (int iz=0; iz < sz; iz++)
				{
					data[Index(ix,iy,iz)].SetLen(0);
				}
			}
		}
	}

private:
	int sx, sy, sz;
	Box3D box;
	IVector* data;

};


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++    BOX2D     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Box2D
{
public:
	Box2D()
	{
		//set empty box:
		pmin = Vector2D(1.,1.);
		pmax = Vector2D(-1.,-1.);
	}
	Box2D(const Vector2D& p1, const Vector2D& p2)
	{
		pmin.X() = Minimum(p1.X(),p2.X());    
		pmin.Y() = Minimum(p1.Y(),p2.Y());    

		pmax.X() = Maximum(p1.X(),p2.X());    
		pmax.Y() = Maximum(p1.Y(),p2.Y());    
	}
	Box2D(const Box2D& b)
	{
		pmin = b.pmin;
		pmax = b.pmax;
	}
	Box2D(const Vector2D& c, double r)
	{
		pmin = c;
		pmax = c;
		Increase(r);
	}
	int Empty() const 
	{
		if (pmin.X() > pmax.X()) {return 1;}
		return 0;
	}
	void Clear()
	{
		pmin = Vector2D(1.,1.);
		pmax = Vector2D(-1.,-1.);
	}

	void Add(const Vector2D & p)
	{
		if (Empty())
		{
			pmin = p;
			pmax = p;
		}
		else
		{
			pmin.X() = Minimum(pmin.X(),p.X());    
			pmin.Y() = Minimum(pmin.Y(),p.Y());    

			pmax.X() = Maximum(pmax.X(),p.X());    
			pmax.Y() = Maximum(pmax.Y(),p.Y());    
		}
	}
	void Add(const Box2D& b)
	{
		if (Empty())
		{
			pmin = b.PMin();
			pmax = b.PMax();
		}
		else
		{
			pmin.X() = Minimum(pmin.X(),b.pmin.X());    
			pmin.Y() = Minimum(pmin.Y(),b.pmin.Y());    

			pmax.X() = Maximum(pmax.X(),b.pmax.X());    
			pmax.Y() = Maximum(pmax.Y(),b.pmax.Y());    
		}
	}
	const Vector2D& PMin() const {return pmin;}
	const Vector2D& PMax() const {return pmax;}
	const double SizeX() const {return pmax[0]-pmin[0];}
	const double SizeY() const {return pmax[1]-pmin[1];}
	Vector2D Center() const {return 0.5*(pmin+pmax);}
	double Radius() const {return (pmax-pmin).Norm();}
	void Increase(double x)
	{
		Increase(x,x);
	}
	void Increase(double x, double y)
	{
		pmin.X() -= x;
		pmin.Y() -= y;
		pmax.X() += x;
		pmax.Y() += y;
	}
  void InflateFactor(double x)
	{
		Vector2D pmi = PMin();
		Vector2D pma = PMax();
		Vector2D pc = Center();
		pmin = pc + ( pmi-pc ) * x;
		pmax = pc + ( pma-pc ) * x;
	}
	int Intersect(const Box2D& b) const 
	{
		if ( pmin[0] > b.pmax[0] || pmax[0] < b.pmin[0]
			|| pmin[1] > b.pmax[1] || pmax[1] < b.pmin[1])	return 0;

		return 1;
	}

	/// return 1 if point p in closure
	int IsIn (const Vector2D & p) const
	{
		if ( pmin[0] <= p[0] && pmax[0] >= p[0]
			&& pmin[1] <= p[1] && pmax[1] >= p[1])
			return 1;

		return 0;
	}

private:
	Vector2D pmin, pmax;
};



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++    SEARCHTREE 2D    ++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Zero-based!!!!
//generate a searchtree with a physical size of Box2D b
//all items must fit into this box and should be equally distributed
//then use AddItems to add items with a bounding box and an identifier
//GetItemsInBox() gives you all identifiers which have a bounding box within the specified box
class SearchTree2D
{
public:
	SearchTree2D():data(0), sx(0), sy(0) {};
	SearchTree2D(int sizex, int sizey, Box2D b)
	{
		box = b;
		if (box.SizeX()*box.SizeY() == 0) box.Increase(1e-4);

		sx = sizex;
		sy = sizey;

		data = new IVector[sx*sy]();
		for (int i = 0; i < sx*sy; i++)
		{
			data[i].Init();
		}
	}
	SearchTree2D(const SearchTree2D& tree)
	{
		box = tree.box;

		sx = tree.sx;
		sy = tree.sy;

		data = new IVector[sx*sy]();
		for (int i = 0; i < sx*sy; i++)
		{
			data[i].Init();
			data[i] = tree.data[i];
		}
	}

	SearchTree2D& operator=(const SearchTree2D& tree)
	{
		if (&tree == this) return *this;

		if (data)
		{
			for (int i = 0; i < sx*sy; i++)
			{
				data[i].Flush();
			}
			delete [] data;
		}

		box = tree.box;

		sx = tree.sx;
		sy = tree.sy;

		data = new IVector[sx*sy]();
		for (int i = 0; i < sx*sy; i++)
		{
			data[i] = tree.data[i];
		}
		return *this;
	}

	virtual ~SearchTree2D()
	{
		if (data)
		{
			for (int i = 0; i < sx*sy; i++)
			{
				data[i].Flush();
			}
			delete [] data;
		}
	}

	void ResetSearchTree(int sizex, int sizey, Box2D b)
	{
		ClearItems(); //empty all items
		box = b;
		if (box.SizeX()*box.SizeY() == 0) box.Increase(1e-4);

		if (sx != sizex || sy != sizey)
		{
			if (data)
			{
				for (int i = 0; i < sx*sy; i++)
				{
					data[i].Flush();
				}
				delete [] data;
			}

			sx = sizex;
			sy = sizey;

			data = new IVector[sx*sy]();
			for (int i = 0; i < sx*sy; i++)
			{
				data[i].Init();
			}
		}
	}

	void GetItemsInBox(const Box2D& b, TArray<int>& items) const
	{
		items.SetLen(0);
		for (int ix=IndX(b.PMin().X()); ix <= IndX(b.PMax().X()); ix++)
		{
			for (int iy=IndY(b.PMin().Y()); iy <= IndY(b.PMax().Y()); iy++)
			{
				for (int i = 1; i <= data[Index(ix,iy)].Length(); i++)
					items.Add((data[Index(ix,iy)])(i));
			}
		}
	}

	void AddItem(const Box2D& b, int identifier)
	{
		for (int ix=IndX(b.PMin().X()); ix <= IndX(b.PMax().X()); ix++)
		{
			for (int iy=IndY(b.PMin().Y()); iy <= IndY(b.PMax().Y()); iy++)
			{
				data[Index(ix,iy)].Add(identifier);
			}
		}
	}
	//return i-index for a double i-value in Box
	int IndX(double x) const 
	{ 
		int i = (int)(((x-box.PMin().X())*sx)/box.SizeX());
		if (i < 0) i = 0;
		if (i >= sx) i = sx-1;
		return i;
	}
	//return i-index for a double i-value in Box
	int IndY(double x) const 
	{ 
		int i = (int)(((x-box.PMin().Y())*sy)/box.SizeY());
		if (i < 0) i = 0;
		if (i >= sy) i = sy-1;
		return i;
	}

	int Index(int x, int y) const
	{
		return x+y*sx;
	}

	void ClearItems()
	{
		for (int ix=0; ix < sx; ix++)
		{
			for (int iy=0; iy < sy; iy++)
			{
				data[Index(ix,iy)].SetLen(0);
			}
		}
	}

private:
	int sx, sy;
	Box2D box;
	IVector* data;

};

//dont use this function!, use KArdan functions instead!
Matrix3D ComputeRot1Rot2Rot3Matrix(const double phi1, const double phi2, const double phi3); //this is not the same as HOTINT Kardan definition !!!!
void QuaternionsToRot1Rot2Rot3Angles(double beta0, double beta1, double beta2, double beta3, Vector3D& phi);

