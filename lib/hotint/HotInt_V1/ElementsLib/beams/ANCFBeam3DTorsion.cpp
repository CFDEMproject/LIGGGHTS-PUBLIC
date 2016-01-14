//#**************************************************************
//# filename:             ANCFBeam3DTorsion.cpp
//#
//# author:               PG & KN
//#
//# generated:						2012
//# description:          3D ANCF beam element with torsion, without shear deformation (BE beam theory)
//# comments:							
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
//#***************************************************************************************
 
#include "body3d.h"
#include "femathhelperfunctions.h"
#include "material.h"
#include "ANCFBeam3Dtorsion.h"
#include "Node.h"
#include "graphicsconstants.h"
#include "elementdataaccess.h"
#include "solversettings_auto.h"


//#define OLD_SHAPEFUNCTIONS

double ANCFBeam3DTorsion::CalculateElementLength() const
{
	//
	// calculation of the length:
	// --------------------------
	//
	// length = \int_{-1}^{1} 1/2 |r'_0(x)| dx,   where x is in [-1,1]   (which is the interval of Gauss points for numerical integration)
	//
	//
	// representation of r'_0(x):
	// ----------------------------
	//
	// by definition, the beam element axis r_0(x) (r_0 is the axis in reference configuration, and x is in [-1,1])
	// is a cubic spline defined by the nodes a and b, i.e., the nodal positions r_0^a and r_0^b, and the nodal slopes r'_0^a and r'_0^b, respectively.
	// let us fix the i-th component (in space), then the i-th component of the beam axis is defined by 
	// (r_0(x))_i = a_0 + a_1 x_0 + a_2 x_0^2 + a_3 x_0^3, where the coefficients a_0,..,a_3 are yet unknown.
	// the coefficients a_0,..,a_3 are found by matching the boundary conditions, 
	// (r_0(-1))_i=(r_0^a)_i, (r_0(1))_i=(r_0^b)_i, (r'_0(-1))_i=(r'_0^a)_i, (r'_0(1))_i=(r'_0^b)_i
	// which leads to solving the linear system
	// [ (r_0^a)_i  ]   [ 1 -1  1 -1 ] [ a_0 ]
	// [ (r'_0^a)_i ] = [ 0  1 -2  3 ] [ a_1 ]
	// [ (r_0^b)_i  ]   [ 1  1  1  1 ] [ a_2 ]
	// [ (r'_0^b)_i ]   [ 0  1  2  3 ] [ a_3 ]
	// or, in short:
	// q0_i = \AMat \aVec.
	// hence, the coefficients read
	// \aVec = \AInv q0_i, where \AInv = \AMat^{-1},
	// and we end up with the representation as an inner product
	// (r'_0(x_0))_i = [  0  1  2(x_0)  3(x_0)^2] \AInv q0_i              (*)
	// 

	double length = 0;

	Vector3D r0a(q0(1),q0(2),q0(3));
	Vector3D r0b(q0(1+6),q0(2+6),q0(3+6));
	Vector3D r0a_prime(q0(4),q0(5),q0(6));
	Vector3D r0b_prime(q0(4+6),q0(5+6),q0(6+6));
	r0a_prime.Normalize();
	r0b_prime.Normalize();
	
	// calculate some geometrical quantities
	Vector3D vector_AB(r0b-r0a);                 // vector from node a to node b
	double distance_AB = vector_AB.Norm();       // distance between node a and node b
	
	// "full circle" beam elements or elements of length 0 are not allowed
	assert(distance_AB > 0);                     
	
	double alpha = acos((r0a_prime*vector_AB)/distance_AB);    // acos returns value between [0,pi]; alpha is angle between vector_AB and slope in node a (which in an arc is equal to the angle between vector_AB and slope in node b) (r0a_prime is already normalized!)
	double sin_alpha = sin(alpha);   // is non-negative, since alpha is non-negative
	
	if (alpha < 1e-14) // if this is a straight beam element
	{
		// keep distance_AB as it is
	}
	else if (alpha > MY_PI - 1e-14)  // if this is an almost full circle
	{
		assert (0 && "full circle beam elements elements are not allowed");
	}
	else // general curved beam element
	{
		distance_AB *= alpha/sin_alpha;    // arc_length = radius * theta,  where theta is the angle between the two sekants of the arc; there holds  radius = distance_AB / (2*sin(alpha)), and theta = 2*alpha
	}

	//multiply slopes at nodes with arclength and divide by two, since x is in [-1,1]
	r0a_prime *= distance_AB/2;
	r0b_prime *= distance_AB/2;
	
	// Matrix containing nodal positions and nodal slope vectors of arc-length divided by two (size of interval [-1,1])
	Matrix Q(4,3);
	
	// fill columns of Q
	for (int j=1; j<=3; j++)
	{
		Q(1,j)=r0a(j);
		Q(2,j)=r0a_prime(j);
		Q(3,j)=r0b(j);
		Q(4,j)=r0b_prime(j);
	}

	// first row (would be nr. 0) can be neglected due to multiplication with 0 in representation of (r'_0(x_0))_i   (see *)
	Matrix AInv(3,4);
	AInv(1,1) = -3;
	AInv(1,2) = -1;
	AInv(1,3) = 3;
	AInv(1,4) = -1;
	AInv(2,1) = 0;
	AInv(2,2) = -1;
	AInv(2,3) = 0;
	AInv(2,4) = 1;
	AInv(3,1) = 1;
	AInv(3,2) = 1;
	AInv(3,3) = -1;
	AInv(3,4) = 1;

	static Vector w,x; //integration points and weights
	GetIntegrationRule(x,w,intorder_axial_strain); 
	
	//loop over all integration points
	for (int j=1; j<=x.GetLen(); j++)  
	{
		
		Vector r0_prime_at_xj = Vector(1.,2.*x(j),3.*Sqr(x(j)))*AInv*Q;
		//mbs->UO() << "r'r0_prime(" << x(j) << ") = " << r0_prime_at_xj << "\n";
		length += w(j)*r0_prime_at_xj.GetNorm();
	}	
	length /= 4.; // factor of AInv: 1/4

	return length;
}

// standard set function
void ANCFBeam3DTorsion::SetANCFBeam3DTorsion(const Vector& xc1, const Vector& xc2, const double& theta1, const double& theta2,
		const Vector3D& director1, const Vector3D& director2, int n1i, int n2i, int materialnumi,
		const Vector3D& coli, const int do_update_directorsi, const int kinematic_computation_modei)
{ 
	n1=n1i; n2=n2i;

	x_init.SetLen(2*SOS()); //2*SOS -> add initial conditions for positions and velocities
	x_init.SetAll(0);
	
	materialnum = materialnumi;
	
	col = coli;
	
	//set vector of generalized coordinates of reference configuration
	q0.SetLen(SOS());
	//position and slope of first node
	for(int i=1;i<=SOSPos()/NNodes();i++)
	{
		q0(i)=xc1(i);
	}
	//position and slope of second node
	for(int i=SOSPos()/NNodes()+1;i<=SOSPos();i++)
	{
		q0(i)=xc2(i-SOSPos()/2);
	}
	//torsional angle at first and second node (measured w.r.t. the director: projection of the director into cross section and normalization -> e30, torsion of e30 around e1=r'/|r'| by the angle theta -> e3 of local frame, finally: e2 = e3 x e1)
	q0(13) = theta1;
	q0(14) = theta2;

	size.X() = CalculateElementLength();
	Material mat = GetMaterial();
	if (mat.IsMaterialOfBeamWithRectangularCrossSection())
	{
		size.Y() = GetMaterial().GetBeamThicknessY();
		size.Z() = GetMaterial().GetBeamThicknessZ();
	}
	else
	{
		assert (GetMaterial().IsMaterialOfBeamWithRectangularCrossSection());   //$ PG 2012-12-18: (as yet) a rectangular cross section (which is specified by the material) is required in ANCFBeam3DTorsion
	}
	mass = size.X()*size.Y()*size.Z()*GetMaterial().Density();

	do_update_directors = do_update_directorsi;

	//kinematic_computation_mode =   0 ... exact terms, 5th order gaussian integration (slow);
	//                               1 ... exact terms, low order (1st order lobatto) integration (fast);
	//                               2 ... approximate mass matrix (torsional terms approximated), no quadratic velocity vector (fastest) --- see also paper by dmitrochenko
	kinematic_computation_mode = kinematic_computation_modei;
	assert(kinematic_computation_mode >= 0 && kinematic_computation_mode <= 2);

	// set directors
	Vector xdatainit(director1.X(), director1.Y(), director1.Z(), director2.X(), director2.Y(), director2.Z());
	SetDataInit(xdatainit);
}

// standard set function with oriented nodes (default set function for script language)
void ANCFBeam3DTorsion::SetANCFBeam3DTorsion(int n1nr, int n2nr, int matnr, const Vector3D& coli, const int do_update_directorsi, const int kinematic_computation_modei)
{
	ANCFNodeS1rot1_3D& n1 = (ANCFNodeS1rot1_3D&)mbs->GetNode(n1nr);
	ANCFNodeS1rot1_3D& n2 = (ANCFNodeS1rot1_3D&)mbs->GetNode(n2nr);
	Matrix3D f1 = n1.GetLocalFrame();
	Matrix3D f2 = n2.GetLocalFrame();
	Vector3D n1_pos = n1.Pos();
	Vector3D n2_pos = n2.Pos();
	Vector3D n1_e1(f1(1,1),f1(2,1),f1(3,1));
	Vector3D n2_e1(f2(1,1),f2(2,1),f2(3,1));
	Vector3D n1_e3(f1(1,3),f1(2,3),f1(3,3));   // direction of director1, when torsional angle theta1 = 0
	Vector3D n2_e3(f2(1,3),f2(2,3),f2(3,3));   // direction of director2, when torsional angle theta2 = 0

	Vector xc1(6);
	Vector xc2(6);

	for (int i=1; i<=Dim(); i++)
	{
		xc1(i) = n1_pos(i);
		xc1(i+Dim()) = n1_e1(i);
		xc2(i) = n2_pos(i);
		xc2(i+Dim()) = n2_e1(i);
	}
	
	// call standard set function
	SetANCFBeam3DTorsion(xc1, xc2, 0, 0, n1_e3, n2_e3, n1.NodeNum(), n2.NodeNum(), matnr, coli, do_update_directorsi, kinematic_computation_modei);
}


// alternative set function with explicit element size (si)
void ANCFBeam3DTorsion::SetANCFBeam3DTorsion(const Vector& xc1, const Vector& xc2, const double& theta1, const double& theta2,
		const Vector3D& director1, const Vector3D& director2, int n1i, int n2i, int materialnumi, const Vector3D& si,
		const Vector3D& coli, const int do_update_directorsi, const int kinematic_computation_modei)
{ 
	// call standard set function
	SetANCFBeam3DTorsion(xc1, xc2, theta1, theta2, director1, director2, n1i, n2i, materialnumi, coli, do_update_directorsi, kinematic_computation_modei);

	// reset size
	size = si;
}

// alternative set function, which determines directors automatically, caution: angles theta1 and theta2 are measured w.r.t. these directors
void ANCFBeam3DTorsion::SetANCFBeam3DTorsion(const Vector& xc1, const Vector& xc2, const double& theta1, const double& theta2,
		int n1i, int n2i, int materialnumi, const Vector3D& si,	const Vector3D& coli, const int do_update_directorsi, const int kinematic_computation_modei)
{
	// determine directors
	Vector3D director;
	Vector3D slopevector1(xc1(4), xc1(5), xc1(6));
	Vector3D slopevector2(xc2(4), xc2(5), xc2(6));
	director = DetermineStandardDirector(slopevector1, slopevector2);

	// call standard set function
	SetANCFBeam3DTorsion(xc1, xc2, theta1, theta2, director, director, n1i, n2i, materialnumi, si, coli, do_update_directorsi, kinematic_computation_modei);
}

// determine one standard director (for the whole element) by the knowledge of the axial slopes in the nodes
// if slopes are significantly distinct, then the director is chosen to be the cross product of both
// else
//   if axial direction is not identical with 3rd kartesian direction then 3rd kartesian direction is chosen for director
//   else 2nd kartesian direction is chosen for director
Vector3D ANCFBeam3DTorsion::DetermineStandardDirector(Vector3D slopevector1, Vector3D slopevector2) const
{
	Vector3D director;

	slopevector1.Normalize();
	slopevector2.Normalize();

	double minimal_distance = 1e-8;
	if ((slopevector1 - slopevector2).Norm2() < minimal_distance)   // slopevectors are almost identical
	{
		director.Set(0., 0., 1.);
		double sqrt_half = sqrt(0.5);
		if (((director - slopevector1).Norm2() < sqrt_half) || ((director + slopevector1).Norm2() < sqrt_half))
		{
			director.Set(0., 1., 0.);
		}
	}
	else   // slopevectors are significantly distinct
	{
		director = slopevector1.Cross(slopevector2);
		director.Normalize();
	}

	return director;
}


void ANCFBeam3DTorsion::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Body3D::GetElementData(edc);

	//ElementData ed;
	//Vector3D director1(GetDirector1()), director2(GetDirector2());
	//ed.SetVector3D(director1.X(), director1.Y(), director1.Z(), "director at node 1"); ed.SetToolTipText("director at node 1"); edc.Add(ed);
	//ed.SetVector3D(director2.X(), director2.Y(), director2.Z(), "director at node 2"); ed.SetToolTipText("director at node 2"); edc.Add(ed);
}

//user has changed data ==> write new data and parameters to element
int ANCFBeam3DTorsion::SetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	int rv = Body3D::SetElementData(edc);

	SetANCFBeam3DTorsion(n1, n2, materialnum, col, do_update_directors, kinematic_computation_mode);

	//Vector3D director1, director2;
	//GetElemDataVector3D(GetMBS(), edc, "director at node 1", director1); 
	//GetElemDataVector3D(GetMBS(), edc, "director at node 2", director2); 
	//SetDirector1(director1);
	//SetDirector2(director2);

	return rv;
}

void ANCFBeam3DTorsion::SetDirector1(const Vector3D& director)
// save director1 in XData
{
	XData(1) = director.X();
	XData(2) = director.Y();
	XData(3) = director.Z();
}

void ANCFBeam3DTorsion::SetDirector2(const Vector3D& director)
// save director2 in XData
{
	XData(4) = director.X();
	XData(5) = director.Y();
	XData(6) = director.Z();
}

Vector3D ANCFBeam3DTorsion::GetDirector1(int flagd) const
// return director1 from XData
{
	if (flagd)
	{
		return Vector3D(XDataD(1), XDataD(2), XDataD(3));
	}
	return Vector3D(XData(1), XData(2), XData(3));
}

Vector3D ANCFBeam3DTorsion::GetDirector2(int flagd) const
// return director2 from XData
{
	if (flagd)
	{
		return Vector3D(XDataD(4), XDataD(5), XDataD(6));
	}
	return Vector3D(XData(4), XData(5), XData(6));
}

int ANCFBeam3DTorsion::UpdateDirectors(const Vector3D& r1dx, const Vector3D& r2dx, const double theta1, const double theta2)
{
	
	//test, if directori is collinear with ridx
	Vector3D dir1 = GetDirector1();
	Vector3D dir2 = GetDirector2();
	double abs_cos1 = fabs(r1dx*dir1) / (r1dx.Norm() * dir1.Norm());
	double abs_cos2 = fabs(r2dx*dir2) / (r2dx.Norm() * dir2.Norm());
	double critical_value = 0.98;     //has to be just a little smaller than one
	if (abs_cos1 > critical_value || abs_cos2 > critical_value)
	{
		return 1;
	}

	Vector3D director_projection(dir1);
	r1dx.GramSchmidt(director_projection);  //project director1 into normal plane of r1dx -> director_projection; GramSchmidt returns normalized director_projection
	SetDirector1(director_projection);   //update director1
	
	// same for second director
	director_projection = dir2;
	r2dx.GramSchmidt(director_projection);
	SetDirector2(director_projection);

	return 0;
}

void ANCFBeam3DTorsion::Initialize()
{
	Body3D::Initialize();
}

void ANCFBeam3DTorsion::LinkToElements()
{
	//Vector of degrees of freedom: [ r1 , r1' , r2 , r2' , phi1 , phi2 ]
	//= 14 DOFs
	//1=left node, 2=right node
	//ri=position vector, ri'=slope vector along the beamcenterline, phi=rotation (torsion)

	if (SOSowned() == 0)
	{
		//UO() << "Link nodes to elements";
		const Node& node1 = GetMBS()->GetNode(n1);
		const Node& node2 = GetMBS()->GetNode(n2);

		//nodes must have dimension 7: r, r', phi (rotation of cross section)
		const int nodedim = 7;
		assert(node1.SOS() == nodedim);
		assert(node2.SOS() == nodedim);

		for (int i=1; i <= nodedim-1; i++)
		{
			AddLTG(node1.Get(i)); //LocalToGlobal-Element-list - left node (6 DOFs = pos. vector and slope vector)
		}
		for (int i=1; i <= nodedim-1; i++)
		{
			AddLTG(node2.Get(i));//right node(6 DOFs = pos. vector and slope vector)
		}
		//add rotation (torsion) degrees of freedom
		AddLTG(node1.Get(nodedim));//1 DOF = rotation of cross section in left node
		AddLTG(node2.Get(nodedim));//1 DOF = rotation of cross section in right node

		//same for velocities:
		for (int i=1; i <= nodedim-1; i++)
		{
			AddLTG(node1.Get(i+nodedim));
		}
		for (int i=1; i <= nodedim-1; i++)
		{
			AddLTG(node2.Get(i+nodedim));
		}
		AddLTG(node1.Get(nodedim*2));
		AddLTG(node2.Get(nodedim*2));
	}
}

//---------------------------------------------------
//Shapefunctions: between -lx/2 and lx/2
//---------------------------------------------------

//---------------------------------------------------
//interpolation of positions r_i and derivatives r_i'
//---------------------------------------------------

//GetShapesPos(sf,X) computes Shapefunctions at position X and saves values in vector sf
void ANCFBeam3DTorsion::GetShapesPos(Vector& sf, double xi) const
{ 
	sf.SetLen(NSPos());  //NSPos=2*nnodes, length(sf)=4
	double lx = size.X();

#ifndef OLD_SHAPEFUNCTIONS // xi e [-1,1]
	double xi2 = xi*xi;
	double xi3 = xi2*xi;
	sf(1) = (2. - 3.*xi + xi3)/4.;
	sf(2) = (1. - xi - xi2 + xi3)*lx/8.;
	sf(3) = (2. + 3.*xi - xi3)/4.;
	sf(4) = (-1. - xi + xi2 + xi3)*lx/8.;
#else // xi e [-L/2,L/2]
	double s=2*xi/lx;//-1<=s<=+1
	double s2 = s*s;
	double s3 = Cub(s);

	sf(1)=(s3 - 3.*s + 2.)/4.;
	sf(2)=(s3 - s2 - s + 1.)/2./lx;
	sf(3)=(-s3 + 3.*s + 2.)/4.;
	sf(4)=(s3 + s2 - s - 1.)/2./lx;
#endif
}

//GetSFPos(i,X) computes i-th Shapefunction at X
double ANCFBeam3DTorsion::GetSFPos(int i, double xi) const
{ 
	double lx = size.X();

#ifndef OLD_SHAPEFUNCTIONS
	//xi between -1 and 1
	double xi2 = xi*xi;
	double xi3 = xi2*xi;
	switch(i)
	{
	case 1: return (2. - 3.*xi + xi3)/4.;
	case 2: return (1. - xi - xi2 + xi3)*lx/8.;
	case 3: return (2. + 3.*xi - xi3)/4.;
	case 4: return (-1. - xi + xi2 + xi3)*lx/8.;
	default: mbs->UO() << "only 4 Shapefunctions for position\n"; return 0.;
	}
#else
	//xi between -length/2 and length/2
	double s=2*xi/lx;
	double s2 = s*s;
	double s3 = Cub(s);

	switch(i)
	{
	case 1: return (s3 - 3.*s + 2.)/4.;
	case 2: return (s3 - s2 - s + 1.)/2./lx;
	case 3: return (-s3 + 3.*s + 2.)/4.;
	case 4: return (s3 + s2 - s - 1.)/2./lx;
	default: mbs->UO()<<"only 4 Shapefunctions for position\n"; return 0.;
	}
#endif

	return 0.;
}

//1st derivative dS/dx
void ANCFBeam3DTorsion::GetShapesPosx(Vector& sf, double xi) const
{
	sf.SetLen(NSPos());
	double lx = size.X();

#ifndef OLD_SHAPEFUNCTIONS 
	//xi between -1 and 1
	double f = 2./lx;
	double xi2 = xi*xi;
	sf(1) = f*(-3. + 3.*xi2)/4.;
	sf(2) = f*(-1. - 2.*xi + 3.*xi2)*lx/8.;
	sf(3) = f*(3. - 3.*xi2)/4.;
	sf(4) = f*(-1. + 2.*xi + 3.*xi2)*lx/8.;
#else
	//xi between -length/2 and length/2
	double s=2*xi/lx;
	double s2 = s*s;
	double lx2 = lx*lx;

	sf(1)=1.5*(s2 - 1.)/lx;
	sf(2)=(3.*s2 - 2.*s -1.)/lx2;
	sf(3)=1.5*(-s2 + 1.)/lx;
	sf(4)=(3.*s2 + 2.*s - 1.)/lx2;
#endif
}

double ANCFBeam3DTorsion::GetSFPosx(int i, double xi) const
{
	double lx = size.X();

#ifndef OLD_SHAPEFUNCTIONS 
	//xi between -1/2 and 1/2
	double xi2 = xi*xi;
	double f = 2./lx;
	switch(i)
	{
	case 1: return f*(- 3. + 3.*xi2)/4.;
	case 2: return f*(-1. - 2.*xi + 3.*xi2)*lx/8.;
	case 3: return f*(3. - 3.*xi2)/4.;
	case 4: return f*(-1. + 2.*xi + 3.*xi2)*lx/8.;
	default: mbs->UO() << "only 4 Shapefunctions for position\n"; return 0.;
	}
#else
	//xi between -length/2 and length/2
	double s=2*xi/lx;
	double s2 = s*s;
	double lx2 = lx*lx;

	switch(i)
	{
	case 1: return 1.5*(s2 - 1.)/lx;
	case 2: return (3.*s2 - 2.*s -1.)/(lx2);
	case 3: return 1.5*(-s2 + 1.)/lx;
	case 4: return(3.*s2 + 2.*s - 1.)/(lx2);
	default: mbs->UO()<<"only 4 Shapefunctions for position\n"; return 0.;
	}
#endif

	return 0.;
}

//2nd derivatives d^2S/dx^2
void ANCFBeam3DTorsion::GetShapesPosxx(Vector& sf, double xi) const
{
	sf.SetLen(NSPos());  
	double lx = size.X();

#ifndef OLD_SHAPEFUNCTIONS 
	//xi between -1/2 and 1/2
	double f = 4./Sqr(lx);
	sf(1) = f*(6.*xi)/4.;
	sf(2) = f*(-2. + 6.*xi)*lx/8.;
	sf(3) = f*(-6.*xi)/4.;
	sf(4) = f*(2. + 6.*xi)*lx/8.;
#else
	//xi between -length/2 and length/2
	double s=2*xi/lx;
	double lx2 = lx*lx;
	double lx3 = lx2*lx;

	sf(1)=6.*s/lx2;
	sf(2)=(12.*s - 4.)/lx3;
	sf(3)=-6.*s/lx2;
	sf(4)=(12.*s + 4.)/lx3;
#endif
}

double ANCFBeam3DTorsion::GetSFPosxx(int i, double xi) const
{
	double lx = size.X();

#ifndef OLD_SHAPEFUNCTIONS 
	//xi between -1/2 and 1/2
	double f = 4./Sqr(lx);
	switch(i)
	{
	case 1: return f*(6.*xi)/4.;
	case 2: return f*(-2. + 6.*xi)*lx/8.;
	case 3: return f*(-6.*xi)/4.;
	case 4: return f*(2. + 6.*xi)*lx/8.;
	default: mbs->UO()<<"only 4 Shapefunctions for position\n"; return 0.;
	}
#else
	//xi between -length/2 and length/2
	double s=2*xi/lx;
	double lx2 = lx*lx;
	double lx3 = lx2*lx;
	switch(i)
	{

	case 1: return 6.*s/lx2;
	case 2: return (12.*s - 4.)/lx3;
	case 3: return -6.*s/lx2;
	case 4: return (12.*s + 4.)/lx3;
	default: mbs->UO()<<"only 4 Shapefunctions for position\n"; return 0.;
	}
#endif

	return 0.;
}

//3rd derivatives d^3S/dx^3
void ANCFBeam3DTorsion::GetShapesPosxxx(Vector& sf, double xi) const
{
	sf.SetLen(NSPos());
	double lx = size.X();

#ifndef OLD_SHAPEFUNCTIONS 
	//xi between -1/2 and 1/2
	double f = 8./Cub(lx);
	sf(1) = f*(6.)/4.;
	sf(2) = f*(6.)*lx/8.;
	sf(3) = f*(-6.)/4.;
	sf(4) = f*(6.)*lx/8.;
#else
	//xi between -length/2 and length/2
	double lx3 = Cub(lx);
	double lx4 = lx3*lx;

	sf(1)=12./lx3;
	sf(2)=24./lx4;
	sf(3)=-12./lx3;
	sf(4)=24./lx4;
#endif
}

double ANCFBeam3DTorsion::GetSFPosxxx(int i, double xi) const
{
	double lx = size.X();

#ifndef OLD_SHAPEFUNCTIONS 
	//xi between -1/2 and 1/2
	double f = 8./Cub(lx);
	switch(i)
	{
	case 1: return f*(6.)/4.;
	case 2: return f*(6.)*lx/8.;
	case 3: return f*(-6.)/4.;
	case 4: return f*(6.)*lx/8.;
	default: mbs->UO()<<"only 4 Shapefunctions for position\n"; return 0.;
	}
#else
	//xi between -length/2 and length/2
	double lx3 = Cub(lx);
	double lx4 = lx3*lx;

	switch(i)
	{
	case 1: return 12./lx3;
	case 2: return 24./lx4;
	case 3: return -12./lx3;
	case 4: return 24./lx4;
	default: mbs->UO()<<"only 4 Shapefunctions for position\n"; return 0.;
	}
#endif

	return 0.;
}

//---------------------------------------------------
//interpolation of angle theta, -1/2<\xi<+L/2
//---------------------------------------------------

void ANCFBeam3DTorsion::GetShapesRot(Vector& sf, double xi) const
{
	double lx = size.X();

	sf.SetLen(NSRot());  //NSRot=2
	sf(1)= 0.5-xi/lx;
	sf(2)= 0.5+xi/lx;
}

double ANCFBeam3DTorsion::GetSFRot(int i, double xi) const
{
	double lx = size.X();

	switch(i)
	{
	case 1: return 0.5-xi/lx;
	case 2: return 0.5+xi/lx;
	default: mbs->UO()<<"only 2 Shapefunctions for rotation\n"; return 0.;
	}
	return 0.;
}

// 1st derivative
void ANCFBeam3DTorsion::GetShapesRotx(Vector& sf, double xi) const
{
	double lx = size.X();

	sf.SetLen(NSRot());  //NSRot=2
	sf(1)=-1./lx;
	sf(2)=1./lx;
}

double ANCFBeam3DTorsion::GetSFRotx(int i, double xi) const
{
	double lx = size.X();

	switch(i)
	{
	case 1: return -1./lx;
	case 2: return 1./lx;
	default: mbs->UO()<<"only 2 Shapefunctions for rotation\n"; return 0.;
	}
	return 0.;
}


//----------------------------------------------
//Position and Displacement for Drawing
//----------------------------------------------
Vector3D ANCFBeam3DTorsion::GetRefPosD() const
{
	return GetPosD(Vector3D(0.,0.,0.));
}

Vector3D ANCFBeam3DTorsion::GetPosD(const Vector3D& p_loc) const
{
	return GetPos(p_loc(1),1) + GetRotMatrix(p_loc(1), 1)*Vector3D(0., p_loc.Y(), p_loc.Z());
}

//Get actual position of relative point p_loc in range [-lx/2..lx/2, etc.]
Vector3D ANCFBeam3DTorsion::GetPos(const Vector3D& p_loc) const
{
	return GetPos(p_loc(1),0) + GetRotMatrix(p_loc(1), 0)*Vector3D(0., p_loc.Y(), p_loc.Z());
}

Vector3D ANCFBeam3DTorsion::GetDisplacement(const Vector3D& p0, int flagD) const
{
	// p0 is in [-lx/2., lx/2.] x [-ly/2., ly/2.] x [-lz/2., lz/2.]
	// compute u(x,y,z) = r(x) - r0(x) + (A(x)-A0(x))(0,y,z)^T
	Vector3D u(0.,0.,0.);
	Vector3D ploc = p0;
	ploc.Scale(size.X()/2.,size.Y()/2.,size.Z()/2.);   // ploc is in [-1, 1]^3

	Matrix3D A = GetRotMatrix(p0(1), flagD);
	Matrix3D A0 = GetInitRotMatrix3D(p0.X());
	for (int i = 1; i <= Dim(); i++)
	{
		//A(i,i) -= 1.;
		u(i) += (A(i,2)-A0(i,2))*p0.Y() + (A(i,3)-A0(i,3))*p0.Z();
		for (int j = 1; j <= NSPos(); j++)
		{
#ifndef OLD_SHAPEFUNCTIONS
			if (flagD)
				u(i) += GetSFPos(j, ploc(1))*(XGD((j-1)*Dim()+i));
			else
				u(i) += GetSFPos(j, ploc(1))*(XG((j-1)*Dim()+i));
#else
			if (flagD)
				u(i) += GetSFPos(j, p0(1))*(XGD((j-1)*Dim()+i));
			else
				u(i) += GetSFPos(j, p0(1))*(XG((j-1)*Dim()+i));
#endif
		}
	}
	return u;
}

Vector3D ANCFBeam3DTorsion::GetVel(double xi, int flagD) const
{
	//status: -l/2<=xi<=l/2
#ifndef OLD_SHAPEFUNCTIONS 
	xi *= 2./size.X();   //status: -1<=xi<=1
#endif
	Vector3D v(0.,0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NSPos(); j++)
		{
			if(flagD==0)
				v(i) += GetSFPos(j,xi)*XGP((j-1)*Dim()+i);
			else
				v(i) += GetSFPos(j,xi)*XGPD((j-1)*Dim()+i);
		}
	}
	return v;
}

Vector3D ANCFBeam3DTorsion::GetVel(const Vector3D& p_loc) const
{
	return GetVel(p_loc(1), 0) + GetRotMatrixP(p_loc(1), 0)*Vector3D(0., p_loc.Y(), p_loc.Z());
}

Vector3D ANCFBeam3DTorsion::GetVelD(const Vector3D& p_loc) const
{
	return GetVel(p_loc(1), 1) + GetRotMatrixP(p_loc(1), 1)*Vector3D(0., p_loc.Y(), p_loc.Z());
}

Vector3D ANCFBeam3DTorsion::GetNodePos(int node_idx) const //returns position of i-th node
{
	switch(node_idx)
	{
	case 1:	return GetPos(Vector3D(-size.X()*.5, 0, 0));
	case 2:	return GetPos(Vector3D(+size.X()*.5, 0, 0));
	}
	assert(0);
	return Vector3D(0.);
}

Vector3D ANCFBeam3DTorsion::GetNodePosD(int node_idx) const //returns position of i-th node (draw mode)
{
	switch(node_idx)
	{
	case 1:	return GetPosD(Vector3D(-size.X()*.5, 0, 0));
	case 2:	return GetPosD(Vector3D(+size.X()*.5, 0, 0));
	}
	assert(0);
	return Vector3D(0.);
}

Vector3D ANCFBeam3DTorsion::GetDOFPosD(int idof) const //returns position of i-th DOF
{
	if (idof > SOS() && idof <= 2*SOS()) idof -= SOS();
	if (idof > 0)
	{
		if (idof < 7 || idof == 13)          //left node
		{
			return GetPosD(Vector3D(-0.5*size.X(),0.,0.));
		}
		else if (idof < 13 || idof == 14)    //right node
		{
			return GetPosD(Vector3D(0.5*size.X(),0.,0.));
		}
	}
	assert(0);
	return Vector3D();
}

Vector3D ANCFBeam3DTorsion::GetDOFDirD(int idof) const //returns direction of action of i-th DOF
{
	if (idof > SOS() && idof <= 2*SOS()) idof -= SOS();

	Matrix3D rot;
	if (idof < 7)
		rot = Matrix3D(1.);//GetInitRotMatrix3D(-size.X()*.5);
	else
		rot = Matrix3D(1.);//GetInitRotMatrix3D(size.X()*.5);

	switch(idof)
	{
	case 1: case 7: return Vector3D(rot(1,1),rot(2,1),rot(3,1)); break;
	case 2: case 8: return Vector3D(rot(1,2),rot(2,2),rot(3,2)); break;
	case 3: case 9: return Vector3D(rot(1,3),rot(2,3),rot(3,3)); break;
	}
	return Vector3D(0.,0.,0.);
}


//-------------------------------------------------------
//GetPos, shapefunctions are defined between -lx/2 and lx/2
//-------------------------------------------------------

//GetPos berechnet Positionsvektor r im deformierten Element
//GetPos computes position r in deformed element
//xg = generalized coordinates (= our vector of unknowns q)
Vector3D ANCFBeam3DTorsion::GetPos(double xi, const Vector& xg) const
{
#ifndef OLD_SHAPEFUNCTIONS
	xi *= 2./size.X();
#endif
	Vector3D p(0.,0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NSPos(); j++)
		{
			p(i) += GetSFPos(j,xi)*(xg((j-1)*Dim()+i)+q0((j-1)*Dim()+i));  //r=u+r_0=Sq+Sq_0 (q_0=initial values, x_init=0-Vector -> no displacement in the beginning))
		}
	}
	return p;
}


//GetPos von oben, aber falls im Aufruf als zweite Eingabe kein Vektor, sondern ein Integer steht,
//wird diese Fkt für Grafik aufgerufen; globale Variable XG wird für Rechnung verwendet anstatt Eingabe xg
//same GetPos function as above, but with the second parameter (flag D) the function for the visualization is activated
//Aufruf mit xg braucht weniger Zeit, weil es sich xg nicht holen muss
Vector3D ANCFBeam3DTorsion::GetPos(double xi, int flagD) const
{
	//status: -l/2<=xi<=l/2
#ifndef OLD_SHAPEFUNCTIONS 
	xi *= 2./size.X();   //status: -1<=xi<=1
#endif
	Vector3D p(0.,0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NSPos(); j++)
		{
			if(flagD==0)
				p(i) += GetSFPos(j,xi)*(XG((j-1)*Dim()+i)+q0((j-1)*Dim()+i));
			else
				p(i) += GetSFPos(j,xi)*(XGD((j-1)*Dim()+i)+q0((j-1)*Dim()+i));
		}
	}
	return p;
}


//Ableitungen von GetPos entsprechen unseren Richtungsvektoren
//GetPosx (Ableitung in Richtung der Balkenachse: dr/dx=r,x)
//GetPosx is the derivative with respect to the direction of the beam centerline
Vector3D ANCFBeam3DTorsion::GetPosx(double xi, const Vector& xg) const
{
#ifndef OLD_SHAPEFUNCTIONS
	xi *= 2./size.X();
#endif
	Vector3D p(0.,0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NSPos(); j++)
		{
			p(i) += GetSFPosx(j,xi)*(xg((j-1)*Dim()+i)+q0((j-1)*Dim()+i));  //r,x=S,x*u+S,x*q_0
		}
	}
	return p;
}


Vector3D ANCFBeam3DTorsion::GetPosx(double xi, int flagD) const
{
#ifndef OLD_SHAPEFUNCTIONS
	xi *= 2./size.X();
#endif
	Vector3D p(0.,0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NSPos(); j++)
		{
			if(flagD==0)
				p(i) += GetSFPosx(j,xi)*(XG((j-1)*Dim()+i)+q0((j-1)*Dim()+i));
			else
				p(i) += GetSFPosx(j,xi)*(XGD((j-1)*Dim()+i)+q0((j-1)*Dim()+i));
		}
	}
	return p;
}


//second derivatives of position vector r with respect to x
Vector3D ANCFBeam3DTorsion::GetPosxx(double xi, const Vector& xg) const
{
#ifndef OLD_SHAPEFUNCTIONS
	xi *= 2./size.X();
#endif
	Vector3D p(0.,0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NSPos(); j++)
		{
			p(i) += GetSFPosxx(j,xi)*(xg((j-1)*Dim()+i)+q0((j-1)*Dim()+i));
		}
	}
	return p;
}


Vector3D ANCFBeam3DTorsion::GetPosxx(double xi, int flagD) const
{
#ifndef OLD_SHAPEFUNCTIONS
	xi *= 2./size.X();
#endif
	Vector3D p(0.,0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NSPos(); j++)
		{
			if(flagD==0)
				p(i) += GetSFPosxx(j,xi)*(XG((j-1)*Dim()+i)+q0((j-1)*Dim()+i));
			else
				p(i) += GetSFPosxx(j,xi)*(XGD((j-1)*Dim()+i)+q0((j-1)*Dim()+i));
		}
	}
	return p;
}


//third derivatives of position vector r with respect to x
Vector3D ANCFBeam3DTorsion::GetPosxxx(double xi, const Vector& xg) const
{
#ifndef OLD_SHAPEFUNCTIONS
	xi *= 2./size.X();
#endif
	Vector3D p(0.,0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NSPos(); j++)
		{
			p(i) += GetSFPosxxx(j,xi)*(xg((j-1)*Dim()+i)+q0((j-1)*Dim()+i));
		}
	}
	return p;
}


Vector3D ANCFBeam3DTorsion::GetPosxxx(double xi, int flagD) const
{
#ifndef OLD_SHAPEFUNCTIONS
	xi *= 2./size.X();
#endif
	Vector3D p(0.,0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NSPos(); j++)
		{
			if(flagD==0)
				p(i) += GetSFPosxxx(j,xi)*(XG((j-1)*Dim()+i)+q0((j-1)*Dim()+i));
			else
				p(i) += GetSFPosxxx(j,xi)*(XGD((j-1)*Dim()+i)+q0((j-1)*Dim()+i));
		}
	}
	return p;
}

Vector3D ANCFBeam3DTorsion::GetPosxP(double xi, int flagD) const
{
#ifndef OLD_SHAPEFUNCTIONS
	xi *= 2./size.X();
#endif
	Vector3D p(0.,0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NSPos(); j++)
		{
			if(flagD==0)
				p(i) += GetSFPosx(j,xi)*XGP((j-1)*Dim()+i);
			else
				p(i) += GetSFPosx(j,xi)*XGPD((j-1)*Dim()+i);
		}
	}
	return p;
}

//---------------------------------------------
//GetRot
//---------------------------------------------

double ANCFBeam3DTorsion::GetRot(double xi, const Vector& xg) const
{
	double a=0;
	for (int j = 1; j <= NSRot(); j++)
	{
		a += GetSFRot(j,xi)*(xg(SOSPos()+j)+q0(SOSPos()+j));  //aus xg und q0 werden nur die letzten beide Einträge verwendet
	}
	return a;
}

double ANCFBeam3DTorsion::GetRot(double xi, int flagD) const
{
	double a=0;
	for (int j = 1; j <= NSRot(); j++)
	{
		if(flagD==0)
			a += GetSFRot(j,xi)*(XG(SOSPos()+j)+q0(SOSPos()+j));
		else
			a += GetSFRot(j,xi)*(XGD(SOSPos()+j)+q0(SOSPos()+j));
	}
	return a;
}

double ANCFBeam3DTorsion::GetRotx(double xi, const Vector& xg) const
{
	double a=0;
	for (int j = 1; j <= NSRot(); j++)
	{
		a += GetSFRotx(j,xi)*(xg(SOSPos()+j)+q0(SOSPos()+j));
	}
	return a;
}

double ANCFBeam3DTorsion::GetRotx(double xi, int flagD) const
{
	double a=0;
	for (int j = 1; j <= NSRot(); j++)
	{
		if(flagD==0)
			a += GetSFRotx(j,xi)*(XG(SOSPos()+j)+q0(SOSPos()+j));
		else
			a += GetSFRotx(j,xi)*(XGD(SOSPos()+j)+q0(SOSPos()+j));
	}
	return a;
}

double ANCFBeam3DTorsion::GetRotP(double xi, int flagD) const
{
	double a=0;
	for (int j = 1; j <= NSRot(); j++)
	{
		if(flagD==0)
			a += GetSFRot(j,xi)*XGP(SOSPos()+j);
		else
			a += GetSFRot(j,xi)*XGPD(SOSPos()+j);
	}
	return a;
}


//---------------------------------------------
//GetRotMatrix
//---------------------------------------------
Matrix3D ANCFBeam3DTorsion::GetRotMatrix(double xi, int flagD) const
{
	// xi in [-lx/2,lx/2]

	Vector3D rx,ry,rz,d;
	double Ang,c,s;

	Ang = GetRot(xi, flagD); //-lx/2<=xi<=+lx/2
	c = cos(Ang);
	s = sin(Ang);

	rx = GetPosx(xi, flagD); //-lx/2<=xi<=+lx/2
	rx.Normalize();

	//d = Vector3D(0.,0.,1.); //first trial, later put this vector into XData and update every step
	d = GetDirector(xi, flagD);
	Vector3D e30 = d;

	rx.GramSchmidt(e30); //d is projected into plane normal to rx (GramSchmidt() already returns a normalized Vector3D)

	Vector3D e20 = e30.Cross(rx);
	//e20.Normalize();   ... is redundant, since |rx|=1, |e30| = 1, rx ortho to e30
	
	ry = e20*c + e30*s;
	rz = e30*c - e20*s;
	
	Matrix3D A(rx.X(),ry.X(),rz.X(),
						 rx.Y(),ry.Y(),rz.Y(),
						 rx.Z(),ry.Z(),rz.Z());
	return A;
}

Matrix3D ANCFBeam3DTorsion::GetRotMatrixP(double xi, int flagD) const
{
	// xi in [-lx/2,lx/2]

	double Ang = GetRot(xi, flagD);
	double AngP = GetRotP(xi, flagD);

	double c = cos(Ang);
	double s = sin(Ang);
	double cP = -s*AngP;
	double sP = c*AngP;

	Vector3D e1 = GetPosx(xi, flagD);		// non normalized e1
	double length_rx = e1.Norm();
	e1 *= (1./length_rx);		// e1 now normalized

	Vector3D e1P = GetPosxP(xi, flagD);
	Vector3D rxP(e1P);		// store for later use (see return)
	e1P *= (1./length_rx);
	e1P -= (e1P*e1)*e1;		// time derivative of normalized e1
	
	Vector3D d = GetDirector(xi, flagD);

	Vector3D e30 = d - (e1*d)*e1;		// GramSchmidt() w/o Normalize()
	Vector3D e30P = -(e1P*d)*e1 - (e1*d)*e1P;		// time derivative of non-normalized e30:   Director d is constant during time step
	double length_e30hat = e30.Norm();
	e30 *= (1./length_e30hat);		// e30 now normalized
	
	e30P *= (1./length_e30hat);
	e30P -= (e30P*e30)*e30;		// time derivative of normalized e30
	
	Vector3D e20 = e30.Cross(e1);		// e20 has unit length, since e30 and e1 are of unit length and orthogonal
	Vector3D e20P = e30P.Cross(e1) + e30.Cross(e1P);
	
	Vector3D ryP = e20P*c + e30P*s + e20*cP + e30*sP;
	Vector3D rzP = e30P*c - e20P*s + e30*cP - e20*sP;
	
	Matrix3D AP(rxP.X(),ryP.X(),rzP.X(),
						 rxP.Y(),ryP.Y(),rzP.Y(),
						 rxP.Z(),ryP.Z(),rzP.Z());
	return AP;
}

Matrix3D ANCFBeam3DTorsion::GetInitRotMatrix3D(double xi) const
{
	// xi in [-lx/2,lx/2]

	Vector3D rx,ry,rz,d;
	double Ang,c,s;

	Ang = GetInitRot3D(xi); //-lx/2<=xi<=+lx/2
	c = cos(Ang);
	s = sin(Ang);

	rx = GetInitPosx3D(xi); //-lx/2<=xi<=+lx/2
	rx.Normalize();

	//d = Vector3D(0.,0.,1.); //first trial, later put this vector into XData and update every step
	d = GetInitDirector(xi);
	Vector3D e30 = d;

	rx.GramSchmidt(e30); //d is projected into plane normal to rx (GramSchmidt() already returns a normalized Vector3D)

	Vector3D e20 = e30.Cross(rx);
	//e20.Normalize();   ... is redundant, since |rx|=1, |e30| = 1, rx ortho to e30
	
	ry = e20*c + e30*s;
	rz = e30*c - e20*s;

	Matrix3D A(rx.X(),ry.X(),rz.X(),
						 rx.Y(),ry.Y(),rz.Y(),
						 rx.Z(),ry.Z(),rz.Z());
	return A;
}

void ANCFBeam3DTorsion::GetdPosdqT(const Vector3D& ploc, Matrix& d)
{
	if (ploc.Y() != 0 || ploc.Z() != 0) {mbs->UO() << "ERROR: ANCFBeam3DTorsion::GetdPosdqT only implemented for y=z=0!!!!!\n";}

	//r = Spos*q; dr/dq = Spos 
	d.SetSize(SOS(),3);
	d.FillWithZeros();

	double xi = ploc.X();
#ifndef OLD_SHAPEFUNCTIONS
	xi *= 2./size.X();
#endif

	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NSPos(); j++)
		{
			d((j-1)*Dim()+i, i) = GetSFPos(j,xi);
			//p(i) += GetSFPos(j,xi)*(xg((j-1)*Dim()+i)+q0((j-1)*Dim()+i));  //r=u+r_0=Sq+Sq_0 (q_0=initial values, x_init=0-Vector -> no displacement in the beginning))
		}
	}
}

void ANCFBeam3DTorsion::GetIntDuDq(Matrix& dudq)
{
	//integrate du/dq over beam volume   => (crossection area) * (integration over axis)
	//u(q) = r(q) - r0  ==>  du/dq = dr/dq

	dudq.SetSize(NSPos()*Dim()+NSRot(), Dim());
	dudq.SetAll(0);

	GetIntegrationRule(x1, w1, 5);		//Interpolation positions: third order 3*2-1 = 5
	double fact = size.X()*0.5*size.Y()*size.Z();				//crosssection area times determinant for line integral

	for (int i1=1; i1<=x1.GetLen(); i1++)
	{
		double scaled_weight = w1(i1) * fact;			//scale weight
		double xi = x1(i1);

#ifdef OLD_SHAPEFUNCTIONS
		xi *= size.X()*0.5;
#endif
		for (int j=1; j<=NSPos(); j++)
		{
			double sf_value = GetSFPos(j, xi);			//value returned by i'th (positional) shape function at point xi
			double block_j = Dim()*(j-1);						//addresses (j-1)st block in dudq
			for (int k=1; k<=Dim(); k++)
			{
				dudq(k+block_j, k) += scaled_weight*sf_value;
			}
		}
	}
}


void ANCFBeam3DTorsion::GetdRotvdqT(const Vector3D& vloc, const Vector3D& ploc, Matrix& drotvdqT)
{
	GetdRotvdqT(vloc, ploc, drotvdqT, 0);
}

//// method calculates gradient of (A(ploc)*vloc) with respect to generalized coordinates (gc), where
//// A(ploc) denotes the Rotation Matrix at the axial position ploc,
//// vloc denotes an input vector, 
//// the result is returned through Matrix& drotvdqT (the T (transpose) means, that drotvdqT has
//// as many rows as partial derivatives and 3 columns (dimension of vector w=A(ploc)*vloc)
////
//// IMPORTANT: if (!use_transposed), then everything is done as explained above
////            if (use_transposed), then vloc is multiplied with transposed of A
//void ANCFBeam3DTorsion::GetdRotvdqT(const Vector3D& vloc, const Vector3D& ploc, Matrix& drotvdqT, int use_transposed)
//{
//	// ploc is in [-lx/2,lx/2] x [-ly/2,ly/2] x [-lz/2,lz/2]
//	// drot/dq = drot/dr . dr/dq = dRot/dr . S, since r=Sq
//
//	double xi=ploc.X();
//	
//	Vector3D rdx = GetPosx(xi);
//	Vector3D d = GetDirector(xi);
//	double angle = GetRot(xi);
//
//	// rot = (rot11 rot12 rot13 rot21 rot22 rot23 rot31 rot32 rot33)
//	// r = (r1', r2', r3', d1, d2, d3, angle)   //notice, (d1,d2,d3)/dq = (0,0,0)
//	// r_pos = (r1', r2', r3')
//	// r_rot = (angle)
//	Matrix drotdrpos(Dim()*Dim(),3);
//
//	//drotdrpos = 
//	//     [dA11dr1 dA11dr2 dA11dr3]
//	//     [dA12dr1 dA12dr2 dA12dr3]
//	//     [dA13dr1 dA13dr2 dA13dr3]
//	//     [dA21dr1 dA21dr2 dA21dr3]
//	//     [dA22dr1 dA22dr2 dA22dr3]
//	//     [dA23dr1 dA23dr2 dA23dr3]
//	//     [dA31dr1 dA31dr2 dA31dr3]
//	//     [dA32dr1 dA32dr2 dA32dr3]
//	//     [dA33dr1 dA33dr2 dA33dr3]
//	GetdRotdrpos(drotdrpos, rdx.X(), rdx.Y(), rdx.Z(), d.X(), d.Y(), d.Z(), angle);
//
//	Matrix drotvdrpos(Dim(),3);    
//	for (int i=1; i<=Dim(); i++)
//	{
//		for (int j=1; j<=3; j++)
//		{
//			if (!use_transposed)
//			{
//				drotvdrpos(i,j) += vloc.X()*drotdrpos(i*Dim()-2,j) + vloc.Y()*drotdrpos(i*Dim()-1,j) + vloc.Z()*drotdrpos(i*Dim(),j);
//			}
//			else
//			{
//			  drotvdrpos(i,j) += vloc.X()*drotdrpos(i,j) + vloc.Y()*drotdrpos(i+Dim(),j) + vloc.Z()*drotdrpos(i+2*Dim(),j);
//			}
//		}
//	}
//	
//	Matrix drotdrrot;
//	drotdrrot.SetSize(Dim()*Dim(),1);
//	GetdRotdrrot(drotdrrot, rdx.X(), rdx.Y(), rdx.Z(), d.X(), d.Y(), d.Z(), angle);
//
//	Matrix drotvdrrot(Dim(),1);
//	for (int i=1; i<=Dim(); i++)
//	{
//		for (int j=1; j<=1; j++)
//		{
//			if (!use_transposed)
//			{
//				drotvdrrot(i,j) += vloc.X()*drotdrrot(i*Dim()-2,j) + vloc.Y()*drotdrrot(i*Dim()-1,j) + vloc.Z()*drotdrrot(i*Dim(),j);
//			}
//			else
//			{
//			  drotvdrrot(i,j) += vloc.X()*drotdrrot(i,j) + vloc.Y()*drotdrrot(i+Dim(),j) + vloc.Z()*drotdrrot(i+2*Dim(),j);
//			}
//		}
//	}
//
//	drotvdqT.SetSize(SOS(),Dim());
//	drotvdqT.SetAll(0.);
//
//	for (int h=1; h<=Dim(); h++)
//	{
//		for (int i=1; i<=Dim(); i++)
//		{
//			for (int j=1; j<=NSPos(); j++)
//			{
//#ifndef OLD_SHAPEFUNCTIONS
//				drotvdqT(i+Dim()*(j-1), h) += drotvdrpos(h,i)*GetSFPosx(j,xi*2./size.X());
//#else
//				drotvdqT(i+Dim()*(j-1), h) += drotvdrpos(h,i)*GetSFPosx(j,xi);
//#endif
//			}
//		}
//		
//		for (int j=1; j<=NSRot(); j++)
//		{
//			drotvdqT(SOSPos()+j, h) += drotvdrrot(h,1)*GetSFRot(j,xi);
//		}
//	}
//}
// method calculates gradient of (A(ploc)*vloc) with respect to generalized coordinates (gc), where
// A(ploc) denotes the Rotation Matrix at the axial position ploc,
// vloc denotes an input vector, 
// the result is returned through Matrix& drotvdqT (the T (transpose) means, that drotvdqT has
// as many rows as partial derivatives and 3 columns (dimension of vector w=A(ploc)*vloc)
//
// IMPORTANT: if (!use_transposed), then everything is done as explained above
//            if (use_transposed), then vloc is multiplied with transposed of A

void ANCFBeam3DTorsion::GetdRotvdqT(const Vector3D& vloc, const Vector3D& ploc, Matrix& drotvdqT, int use_transposed)
{
	// ploc is in [-lx/2,lx/2] x [-ly/2,ly/2] x [-lz/2,lz/2]
	// drot/dq = drot/dr . dr/dq = dRot/dr . S, since r=Sq

	double xi=ploc.X();
	
	Vector3D rdx = GetPosx(xi);
	Vector3D d = GetDirector(xi);
	double angle = GetRot(xi);

	// rot = (rot11 rot12 rot13 rot21 rot22 rot23 rot31 rot32 rot33)
	// r = (r1', r2', r3', d1, d2, d3, angle)   //notice, (d1,d2,d3)/dq = (0,0,0)
	// r_pos = (r1', r2', r3')
	// r_rot = (angle)
	Matrix drotdrpos(Dim()*Dim(),3);
	Matrix drotdrposT(Dim()*Dim(),3);

	//drotdrpos = 
	//     [dA11dr1 dA11dr2 dA11dr3]
	//     [dA12dr1 dA12dr2 dA12dr3]
	//     [dA13dr1 dA13dr2 dA13dr3]
	//     [dA21dr1 dA21dr2 dA21dr3]
	//     [dA22dr1 dA22dr2 dA22dr3]
	//     [dA23dr1 dA23dr2 dA23dr3]
	//     [dA31dr1 dA31dr2 dA31dr3]
	//     [dA32dr1 dA32dr2 dA32dr3]
	//     [dA33dr1 dA33dr2 dA33dr3]
	GetdRotdrpos(drotdrpos, rdx.X(), rdx.Y(), rdx.Z(), d.X(), d.Y(), d.Z(), angle);
	GetdRotdrposT(drotdrposT, rdx.X(), rdx.Y(), rdx.Z(), d.X(), d.Y(), d.Z(), angle);

	Matrix drotvdrpos(Dim(),3);    
	for (int i=1; i<=Dim(); i++)
	{
		for (int j=1; j<=3; j++)
		{
			if (!use_transposed)
			{
				drotvdrpos(i,j) += vloc.X()*drotdrpos(i*Dim()-2,j) + vloc.Y()*drotdrpos(i*Dim()-1,j) + vloc.Z()*drotdrpos(i*Dim(),j);
			}
			else
			{
			  drotvdrpos(i,j) += vloc.X()*drotdrposT(i*Dim()-2,j) + vloc.Y()*drotdrposT(i*Dim()-1,j) + vloc.Z()*drotdrposT(i*Dim(),j);
			}
		}
	}
	
	Matrix drotdrrot;
	Matrix drotdrrotT;
	drotdrrot.SetSize(Dim()*Dim(),1);
	drotdrrotT.SetSize(Dim()*Dim(),1);
	
	GetdRotdrrot(drotdrrot, rdx.X(), rdx.Y(), rdx.Z(), d.X(), d.Y(), d.Z(), angle);
	GetdRotdrrotT(drotdrrotT, rdx.X(), rdx.Y(), rdx.Z(), d.X(), d.Y(), d.Z(), angle);

	Matrix drotvdrrot(Dim(),1);
	for (int i=1; i<=Dim(); i++)
	{
		for (int j=1; j<=1; j++)
		{
			if (!use_transposed)
			{
				drotvdrrot(i,j) += vloc.X()*drotdrrot(i*Dim()-2,j) + vloc.Y()*drotdrrot(i*Dim()-1,j) + vloc.Z()*drotdrrot(i*Dim(),j);
			}
			else
			{
			  drotvdrrot(i,j) += vloc.X()*drotdrrotT(i*Dim()-2,j) + vloc.Y()*drotdrrotT(i*Dim()-1,j) + vloc.Z()*drotdrrotT(i*Dim(),j);
			}
		}
	}

	drotvdqT.SetSize(SOS(),Dim());
	drotvdqT.SetAll(0.);

	for (int h=1; h<=Dim(); h++)
	{
		for (int i=1; i<=Dim(); i++)
		{
			for (int j=1; j<=NSPos(); j++)
			{
#ifndef OLD_SHAPEFUNCTIONS
				drotvdqT(i+Dim()*(j-1), h) += drotvdrpos(h,i)*GetSFPosx(j,xi*2./size.X());
#else
				drotvdqT(i+Dim()*(j-1), h) += drotvdrpos(h,i)*GetSFPosx(j,xi);
#endif
			}
		}
		
		for (int j=1; j<=NSRot(); j++)
		{
			drotvdqT(SOSPos()+j, h) += drotvdrrot(h,1)*GetSFRot(j,xi);
		}
	}
}

void ANCFBeam3DTorsion::GetdRotdqT(const Vector3D& ploc, Matrix& d)
{
	//ploc in [-lx/2,lx/2] x [-ly/2,ly/2] x [-lz/2,lz/2]
	////////////////////////////////////////////////////
	
	int approach = 3;     // 1..JG, 2..YV, 3..KN,PG

	if (approach == 1)
	{
		//JG's approach (probably just for special cases)
		//here, the routine GetdRotvdqT has to sum up (e_j)_n v_j (instead of (e_j)_n v_n)
		//UO() << "dRotdq\n";

		Vector3D p0=ploc;
		//p0.Scale(0.5*size.X(),0.5*size.Y(),0.5*size.Z());

		d.SetSize(SOS(),3);
		d.FillWithZeros(); //needed???

		Matrix dAvdq(SOS(),3);
		double x = ploc.X();

		Matrix3D AT = GetRotMatrix(ploc).GetTp();

		Vector3D v;
		v = Vector3D(AT(1,1),AT(2,1),AT(3,1));
		GetdRotvdqT(v, ploc, dAvdq, 0);
		for (int j=1; j <= SOS(); j++)
		{
			d(j,3) = dAvdq(j,2);
			d(j,2) =-dAvdq(j,3);
		}
		v = Vector3D(AT(1,2),AT(2,2),AT(3,2));
		GetdRotvdqT(v, ploc, dAvdq, 0);
		for (int j=1; j <= SOS(); j++)
		{
			d(j,1) = dAvdq(j,3);
		}
	}
	else if (approach == 2)
	{
	////////////////////////////////////////////////////////////
	// YV's approach
	// let e_k denote the k-th base vector of the actual system
	// \delta \theta = 0.5 * (e_k) \times \delta (e_k)
	// (\grad \theta)_{l,m} = \eps_{ijm} (e_k)_i (\partial (e_k)_j / \partial q_l)
	//here, the routine GetdRotvdqT has to sum up (e_j)_n v_n (instead of (e_j)_n v_j)
	
		d.SetSize(SOS(),Dim()); d.SetAll(0.);
		Matrix drotvdqT;
		drotvdqT.SetSize(SOS(),Dim());
		Vector3D e;
		Vector3D epse;  //  Levi-Civita-Symbol multiplyed with base vector e
		Matrix3D AT = GetRotMatrix(ploc).GetTp();
		for (int m=1; m<=Dim(); m++) {
			for (int k=1; k<=Dim(); k++) {
				e = Vector3D(AT(k,1),AT(k,2),AT(k,3));
				for (int j=1; j<=Dim(); j++) {
					epse(j) = 0.;
					for (int i=1; i<=Dim(); i++) {
						epse(j) += (i-j)*(j-m)*(m-i)*e(i);
					}
				}
				GetdRotvdqT(0.25*epse, ploc, drotvdqT, 1);
				for (int l=1; l<=SOS(); l++) {
					d(l,m) += drotvdqT(l,k);
				}
			}
		}
	}
	else
	{
		////////////////////////////////////////////////////////////
		//KN&PG's approach from \cite{book_hodges} 
		//let e_i denote the i-th base vector of the actual system (i-th column of the rotation matrix)
		//\delta \theta = \sum_{i=1}^3  e_i (\delta e_j . e_k)    where j = i%3+1, k = j%3+1
		//(\grad_q \theta^T)_{l,m} = (\partial q_l/ \partial (e_j)_n) (e_k)_n (e_i)_m
		//here, the routine GetdRotvdqT has to sum up (e_j)_n v_n (instead of (e_j)_n v_j)
		//note: ploc is in [-lx/2,lx/2] x [-ly/2,ly/2] x [-lz/2,lz/2]
		d.SetSize(SOS(),Dim());
		d.FillWithZeros();

		Matrix3D A = GetRotMatrix(ploc);
		//mbs->UO() << "RotMatrix = \n" << A << "\n";

		Matrix drotvdqT(SOS(),Dim());
		for (int i=1; i<=Dim(); i++) 
		{
			int j = i%3 + 1;
			int k = j%3 + 1;
			Vector3D ei(A(1,i), A(2,i), A(3,i));
			GetdRotvdqT(Vector3D(A(1,k), A(2,k), A(3,k)), ploc, drotvdqT, 1);
			for (int l=1; l<=SOS(); l++)
				for (int m=1; m<=Dim(); m++)
					d(l,m) += drotvdqT(l,j)*ei(m);
		}
	}

}


Vector3D ANCFBeam3DTorsion::GetInitDirector(const double& xi) const
{
	Vector d = GetDataInit();
	Vector3D director1(d(1), d(2), d(3));
	Vector3D director2(d(4), d(5), d(6));
	return GetSFRot(1,xi)*director1 + GetSFRot(2,xi)*director2;
}

Vector3D ANCFBeam3DTorsion::GetDirector(const double& xi, int flagD) const
{
	return GetSFRot(1,xi)*GetDirector1(flagD) + GetSFRot(2,xi)*GetDirector2(flagD);
}

Vector3D ANCFBeam3DTorsion::GetDirectorx(const double& xi, int flagD) const
{
	return GetSFRotx(1,xi)*GetDirector1(flagD) + GetSFRotx(2,xi)*GetDirector2(flagD);
}

//------------------------------------------------------------
//computes position vector in the reference element: r_0 
//------------------------------------------------------------
Vector3D ANCFBeam3DTorsion::GetInitPos3D(double xi) const
{
#ifndef OLD_SHAPEFUNCTIONS
	xi *= 2./size.X();
#endif
	Vector3D p(0.,0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NSPos(); j++)
		{
			p(i) += GetSFPos(j,xi)*(q0((j-1)*Dim()+i));  //r_0=Sq_0; (q_0 includes initial values)
		}
	}
	return p;
}

Vector3D ANCFBeam3DTorsion::GetInitPosx3D(double xi) const
{
#ifndef OLD_SHAPEFUNCTIONS
	xi *= 2./size.X();
#endif
	Vector3D p(0.,0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NSPos(); j++)
		{
			p(i) += GetSFPosx(j,xi)*(q0((j-1)*Dim()+i));  //rx_0=Sx q_0; (q_0 includes initial values)
		}
	}
	return p;
}

double ANCFBeam3DTorsion::GetInitRot3D(double xi) const
{
	double a=0;
	for (int j = 1; j <= NSRot(); j++)
	{
		a += GetSFRot(j,xi)*q0(SOSPos()+j);
	}
	return a;
}

double ANCFBeam3DTorsion::GetInitRotx3D(double xi) const
{
	double a=0;
	for (int j = 1; j <= NSRot(); j++)
	{
		a += GetSFRotx(j,xi)*q0(SOSPos()+j);
	}
	return a;
}

//-------
//strains
//-------

void ANCFBeam3DTorsion::GetGamma1(const double& xi, double& Gamma1, int flagD) const
{
	// Gamma1 = (|r'|-|r0'|)/(|r0'|)
	Vector3D rdx = GetPosx(xi,flagD);
	Gamma1 = rdx.Norm();
	rdx = GetInitPosx3D(xi);
	double rdxnorm = rdx.Norm();
	Gamma1 = (Gamma1 - rdxnorm)/rdxnorm;
}

void ANCFBeam3DTorsion::GetDeltaGamma1(const double& xi, Vector& DeltaGamma1) const
{
	// delta gamma1 = (|r0'| |r'|)^{-1} r' . S'
	double fac = 1.;
	Vector3D rdx = GetInitPosx3D(xi);
	fac /= rdx.Norm();
	rdx = GetPosx(xi);
	fac /= rdx.Norm();

	for (int i=1; i<=Dim(); i++)
	{
		for (int j=1; j<=NSPos(); j++)
		{
#ifndef OLD_SHAPEFUNCTIONS
			DeltaGamma1(3*(j-1)+i) = fac*rdx(i)*GetSFPosx(j,xi*2./size.X());
#else	
			DeltaGamma1(3*(j-1)+i) = fac*rdx(i)*GetSFPosx(j,xi);
#endif
		}
	}

	DeltaGamma1(13) = DeltaGamma1(14) = 0.;
}

void ANCFBeam3DTorsion::GetKappa1(const double& xi, double& Kappa1, int flagD) const
{
	Vector3D rx = GetPosx(xi, flagD);
	Vector3D rxx = GetPosxx(xi, flagD);
	double ang = GetRot(xi, flagD) - GetInitRot3D(xi);
	double angx = GetRotx(xi, flagD) - GetInitRotx3D(xi);
	Vector3D d = GetDirector(xi, flagD);
	Vector3D dx = GetDirectorx(xi, flagD);

	Kappa1 = GetKappa1(rx.X(), rx.Y(), rx.Z(), rxx.X(), rxx.Y(), rxx.Z(), ang, angx, d.X(), d.Y(),d.Z(), dx.X(), dx.Y(), dx.Z());
}

void ANCFBeam3DTorsion::GetKappa2(const double& xi, double& Kappa2, int flagD) const
{
	Vector3D rx = GetPosx(xi, flagD);
	Vector3D rxx = GetPosxx(xi, flagD);
	double ang = GetRot(xi, flagD) - GetInitRot3D(xi);
	double angx = GetRotx(xi, flagD) - GetInitRotx3D(xi);
	Vector3D d = GetDirector(xi, flagD);
	Vector3D dx = GetDirectorx(xi, flagD);

	Kappa2 = GetKappa2(rx.X(), rx.Y(), rx.Z(), rxx.X(), rxx.Y(), rxx.Z(), ang, angx, d.X(), d.Y(),d.Z(), dx.X(), dx.Y(), dx.Z());
}

void ANCFBeam3DTorsion::GetKappa3(const double& xi, double& Kappa3, int flagD) const
{
	Vector3D rx = GetPosx(xi, flagD);
	Vector3D rxx = GetPosxx(xi, flagD);
	double ang = GetRot(xi, flagD) - GetInitRot3D(xi);
	double angx = (GetRotx(xi, flagD) - GetInitRotx3D(xi));
	Vector3D d = GetDirector(xi, flagD);
	Vector3D dx = GetDirectorx(xi, flagD);

	Kappa3 = GetKappa3(rx.X(), rx.Y(), rx.Z(), rxx.X(), rxx.Y(), rxx.Z(), ang, angx, d.X(), d.Y(),d.Z(), dx.X(), dx.Y(), dx.Z());
}

void ANCFBeam3DTorsion::GetDeltaKappa1(const double& xi, Vector& DeltaKappa1) const   //derivative of Kappa1 wrt q at xi
{
	Vector3D rx = GetPosx(xi);
	Vector3D rxx = GetPosxx(xi);
	double ang = GetRot(xi) - GetInitRot3D(xi);
	double angx = GetRotx(xi) - GetInitRotx3D(xi);
	Vector3D d = GetDirector(xi);
	Vector3D dx = GetDirectorx(xi);

	// delta kappa1 = [ d kappa1 / d (r',r'',ang,ang') ] . [ d (r',r'',ang,ang') / dq ]
	Vector3D kappa_pos;
	double kappa_angle;

	kappa_pos = GetDKappa1DPosx(rx.X(), rx.Y(), rx.Z(), rxx.X(), rxx.Y(), rxx.Z(), ang, angx, d.X(), d.Y(),d.Z(), dx.X(), dx.Y(), dx.Z());
	DeltaKappa1.SetAll(0.);
	for (int i=1; i<=Dim(); i++)
	{
		for (int j=1; j<=NSPos(); j++)
		{
#ifndef OLD_SHAPEFUNCTIONS
			DeltaKappa1(i+Dim()*(j-1)) += kappa_pos(i)*GetSFPosx(j,xi*2./size.X());
#else
			DeltaKappa1(i+Dim()*(j-1)) += kappa_pos(i)*GetSFPosx(j,xi);
#endif
		}
	}
	
	kappa_pos = GetDKappa1DPosxx(rx.X(), rx.Y(), rx.Z(), rxx.X(), rxx.Y(), rxx.Z(), ang, angx, d.X(), d.Y(),d.Z(), dx.X(), dx.Y(), dx.Z());
	for (int i=1; i<=Dim(); i++)
	{
		for (int j=1; j<=NSPos(); j++)
		{
#ifndef OLD_SHAPEFUNCTIONS
			DeltaKappa1(i+Dim()*(j-1)) += kappa_pos(i)*GetSFPosxx(j,xi*2./size.X());
#else
			DeltaKappa1(i+Dim()*(j-1)) += kappa_pos(i)*GetSFPosxx(j,xi);
#endif
		}
	}

	kappa_angle = GetDKappa1DAng(rx.X(), rx.Y(), rx.Z(), rxx.X(), rxx.Y(), rxx.Z(), ang, angx, d.X(), d.Y(),d.Z(), dx.X(), dx.Y(), dx.Z());
	for (int j=1; j<=NSRot(); j++)
	{
		DeltaKappa1(Dim()*NSPos()+j) += kappa_angle*GetSFRot(j,xi);
	}

	kappa_angle = GetDKappa1DAngx(rx.X(), rx.Y(), rx.Z(), rxx.X(), rxx.Y(), rxx.Z(), ang, angx, d.X(), d.Y(),d.Z(), dx.X(), dx.Y(), dx.Z());
	for (int j=1; j<=NSRot(); j++)
	{
		DeltaKappa1(Dim()*NSPos()+j) += kappa_angle*GetSFRotx(j,xi);
	}
}

void ANCFBeam3DTorsion::GetDeltaKappa2(const double& xi, Vector& DeltaKappa2) const   //derivative of Kappa2 wrt q at xi
{
	Vector3D rx = GetPosx(xi);
	Vector3D rxx = GetPosxx(xi);
	double ang = GetRot(xi) - GetInitRot3D(xi);
	double angx = GetRotx(xi) - GetInitRotx3D(xi);
	Vector3D d = GetDirector(xi);
	Vector3D dx = GetDirectorx(xi);

	// delta Kappa2 = [ d Kappa2 / d (r',r'',ang,ang') ] . [ d (r',r'',ang,ang') / dq ]
	Vector3D kappa_pos;
	double kappa_angle;

	kappa_pos = GetDKappa2DPosx(rx.X(), rx.Y(), rx.Z(), rxx.X(), rxx.Y(), rxx.Z(), ang, angx, d.X(), d.Y(),d.Z(), dx.X(), dx.Y(), dx.Z());
	DeltaKappa2.SetAll(0.);
	for (int i=1; i<=Dim(); i++)
	{
		for (int j=1; j<=NSPos(); j++)
		{
#ifndef OLD_SHAPEFUNCTIONS
			DeltaKappa2(i+Dim()*(j-1)) += kappa_pos(i)*GetSFPosx(j,xi*2./size.X());
#else
			DeltaKappa2(i+Dim()*(j-1)) += kappa_pos(i)*GetSFPosx(j,xi);
#endif
		}
	}
	
	kappa_pos = GetDKappa2DPosxx(rx.X(), rx.Y(), rx.Z(), rxx.X(), rxx.Y(), rxx.Z(), ang, angx, d.X(), d.Y(),d.Z(), dx.X(), dx.Y(), dx.Z());
	for (int i=1; i<=Dim(); i++)
	{
		for (int j=1; j<=NSPos(); j++)
		{
#ifndef OLD_SHAPEFUNCTIONS
			DeltaKappa2(i+Dim()*(j-1)) += kappa_pos(i)*GetSFPosxx(j,xi*2./size.X());
#else
			DeltaKappa2(i+Dim()*(j-1)) += kappa_pos(i)*GetSFPosxx(j,xi);
#endif
		}
	}

	kappa_angle = GetDKappa2DAng(rx.X(), rx.Y(), rx.Z(), rxx.X(), rxx.Y(), rxx.Z(), ang, angx, d.X(), d.Y(),d.Z(), dx.X(), dx.Y(), dx.Z());
	for (int j=1; j<=NSRot(); j++)
	{
		DeltaKappa2(Dim()*NSPos()+j) += kappa_angle*GetSFRot(j,xi);
	}

	kappa_angle = GetDKappa2DAngx(rx.X(), rx.Y(), rx.Z(), rxx.X(), rxx.Y(), rxx.Z(), ang, angx, d.X(), d.Y(),d.Z(), dx.X(), dx.Y(), dx.Z());
	for (int j=1; j<=NSRot(); j++)
	{
		DeltaKappa2(Dim()*NSPos()+j) += kappa_angle*GetSFRotx(j,xi);
	}
}

void ANCFBeam3DTorsion::GetDeltaKappa3(const double& xi, Vector& DeltaKappa3) const   //derivative of Kappa3 wrt q at xi
{
	Vector3D rx = GetPosx(xi);
	Vector3D rxx = GetPosxx(xi);
	double ang = GetRot(xi) - GetInitRot3D(xi);
	double angx = (GetRotx(xi) - GetInitRotx3D(xi));
	Vector3D d = GetDirector(xi);
	Vector3D dx = GetDirectorx(xi);

	// delta Kappa3 = [ d Kappa3 / d (r',r'',ang,ang') ] . [ d (r',r'',ang,ang') / dq ]
	Vector3D kappa_pos;
	double kappa_angle;

	kappa_pos = GetDKappa3DPosx(rx.X(), rx.Y(), rx.Z(), rxx.X(), rxx.Y(), rxx.Z(), ang, angx, d.X(), d.Y(),d.Z(), dx.X(), dx.Y(), dx.Z());
	DeltaKappa3.SetAll(0.);
	for (int i=1; i<=Dim(); i++)
	{
		for (int j=1; j<=NSPos(); j++)
		{
#ifndef OLD_SHAPEFUNCTIONS
			DeltaKappa3(i+Dim()*(j-1)) += kappa_pos(i)*GetSFPosx(j,xi*2./size.X());
#else
			DeltaKappa3(i+Dim()*(j-1)) += kappa_pos(i)*GetSFPosx(j,xi);
#endif
		}
	}
	
	kappa_pos = GetDKappa3DPosxx(rx.X(), rx.Y(), rx.Z(), rxx.X(), rxx.Y(), rxx.Z(), ang, angx, d.X(), d.Y(),d.Z(), dx.X(), dx.Y(), dx.Z());
	for (int i=1; i<=Dim(); i++)
	{
		for (int j=1; j<=NSPos(); j++)
		{
#ifndef OLD_SHAPEFUNCTIONS
			DeltaKappa3(i+Dim()*(j-1)) += kappa_pos(i)*GetSFPosxx(j,xi*2./size.X());
#else
			DeltaKappa3(i+Dim()*(j-1)) += kappa_pos(i)*GetSFPosxx(j,xi);
#endif
		}
	}

	kappa_angle = GetDKappa3DAng(rx.X(), rx.Y(), rx.Z(), rxx.X(), rxx.Y(), rxx.Z(), ang, angx, d.X(), d.Y(),d.Z(), dx.X(), dx.Y(), dx.Z());
	for (int j=1; j<=NSRot(); j++)
	{
		DeltaKappa3(Dim()*NSPos()+j) += kappa_angle*GetSFRot(j,xi);
	}

	kappa_angle = GetDKappa3DAngx(rx.X(), rx.Y(), rx.Z(), rxx.X(), rxx.Y(), rxx.Z(), ang, angx, d.X(), d.Y(),d.Z(), dx.X(), dx.Y(), dx.Z());
	for (int j=1; j<=NSRot(); j++)
	{
		DeltaKappa3(Dim()*NSPos()+j) += kappa_angle*GetSFRotx(j,xi);
	}
}

//------------------------------------------------------------------------------
//variational formulation: Delta U (Beam3D-Paper: equation (15))
//we calculate residual
//------------------------------------------------------------------------------
void ANCFBeam3DTorsion::EvalF2(Vector& f, double t)     // dynamic test example: 74.0118% of CPU Time
{

	Body3D::EvalF2(f,t);

	static Vector fadd;//element residual
	fadd.SetLen(SOS());
	fadd.SetAll(0.);

	double lx = size.X();

	{ // dynamic test example: 62.2240% of CPU Time

		double beamGJ = GetMaterial().BeamGJkx();
		double beamEIy = GetMaterial().BeamEIy();
		double beamEIz = GetMaterial().BeamEIz();
		double beamEA = GetMaterial().BeamEA();

		static Vector wT,xT; //integration points and weights
		GetIntegrationRule(xT,wT,intorder_curvature); 

		for (int i1=1; i1<=xT.GetLen(); i1++)  //i1-loop over all integration points
		{
			double x = xT(i1);		//x in unit configuration [-1 , 1]
			double xi = x*lx*0.5;	//xi in scaled configuration [-lx/2 , lx/2]
			double det = 0.5*lx;	//determinant of the jacobian (from element transformation [-1 , 1] -> [-lx/2 , lx/2])

			double fact_kappa1= beamGJ*wT(i1)*det;		// beamGJ  = torsional stiffness
			double fact_kappa2= beamEIy*wT(i1)*det;		// beamEIy = bending stiffness (in e2=y direction)
			double fact_kappa3= beamEIz*wT(i1)*det;		// beamEIz = bending stiffness (in e3=z direction)


			ConstVector<14> DeltaKappa1; GetDeltaKappa1(xi,DeltaKappa1);
			double Kappa1; GetKappa1(xi,Kappa1);

			ConstVector<14> DeltaKappa2; GetDeltaKappa2(xi,DeltaKappa2);
			double Kappa2; GetKappa2(xi,Kappa2);

			ConstVector<14> DeltaKappa3; GetDeltaKappa3(xi,DeltaKappa3);
			double Kappa3; GetKappa3(xi,Kappa3);

			fadd +=
				(fact_kappa1*Kappa1)*DeltaKappa1 +    // around x axis (beam axis)
				(fact_kappa2*Kappa2)*DeltaKappa2 +		// around y axis
				(fact_kappa3*Kappa3)*DeltaKappa3; 		// around z-axis
		}

		GetIntegrationRule(xT,wT,intorder_axial_strain); 
		for (int i1=1; i1<=xT.GetLen(); i1++)  //i1-loop over all integration points
		{
			double x = xT(i1);		//x in unit configuration [-1 , 1]
			double xi = x*lx*0.5;	//xi in scaled configuration [-lx/2 , lx/2]
			double det = 0.5*lx;	//determinant of the jacobian (from element transformation [-1 , 1] -> [-lx/2 , lx/2])

			double fact_gamma1= beamEA*wT(i1)*det;		// beamEA  = axial stiffness

			ConstVector<14> DeltaGamma1; GetDeltaGamma1(xi,DeltaGamma1);
			double Gamma1; GetGamma1(xi,Gamma1);

			fadd +=
				(fact_gamma1*Gamma1)*DeltaGamma1;
		}

		//Residual:
		//in ANCFBeam3DTorsion::EvalF2:
		//from Element::EvalF2 -> Loads, Constraints (=f)
		//additional -> fadd = Ku
		//f-Ku (right hand side of second order DE) -> return value
		f -= fadd;  //f=f-fadd; Ku=f -> Residual=f-Ku, Ku=fadd
	}

	{ // dynamic test example:  10.3733% of CPU Time
		//// add quadratic velocity vector
		if (!mbs->GetSolSet().dostaticcomputation && !GetMBS()->IsJacobianComputation() && kinematic_computation_mode != 2)
		{
			EvalQVV(fadd,t);
			f -= fadd;
		}
	}
}

void ANCFBeam3DTorsion::EvalQVV(Vector& f, double t)
{
	f.SetAll(0.);

	double det = 0.5*size.X();
	ConstVector<14> qvv(14, 1); //nofill

	if (kinematic_computation_mode == 0)
	{			
		Vector x1, w1;
		GetIntegrationRule(x1,w1,intorder_mass);
		
		for (int i1=1; i1<=x1.GetLen(); i1++)
		{
			double xi = x1(i1)*det;
			GetQuadraticVelocityVector(xi, qvv); //get quadratic velocity vector at integration point
			f += (det*w1(i1))*qvv;
		}
	}
	else if (kinematic_computation_mode == 1)
	{
		for (int n=1; n<=2; n++)   // 2nd order Lobatto integration (IPs are at node n=1 and n=2)
		{
			GetQuadraticVelocityVectorAtNode(qvv, n);   // determinant multiplied inside!
			f += qvv;
		}
	}
}

void ANCFBeam3DTorsion::EvalM(Matrix& m, double t) 
{
	if (kinematic_computation_mode == 2)   // 2 ... approximate mass matrix (torsional terms approximated), no quadratic velocity vector (fastest)
	{
		if (massmatrix.Getcols() == SOS())
		{
			m = massmatrix;
			return;
		}
		else
		{
			m.SetAll(0);
			// space dimension, = 3
			int dim = Dim();
			// number of shape functions
			int nspos = NSPos();
			int nsrot = NSRot();
			// storage vector for shape functions
			Vector SVpos(nspos);
			Vector SVrot(nsrot);
			// shape function matrix
			// (r1 )   ( Spos1   0     0   Spos2   0      0  Spos3   0     0  Spos4   0     0     0     0   )
			// (r2 ) = (  0    Spos1   0     0   Spos2    0    0   Spos3   0    0   Spos4   0     0     0   )   q^T
			// (r3 )   (  0      0   Spos1   0     0   Spos2   0     0   Spos3  0     0   Spos4   0     0   )
			// (ang)   (  0      0     0     0     0      0    0     0     0    0     0     0   ///////////////////Srot1 Srot2 )

			Matrix HL(SOS(),dim+1);

			// integration rule
			Vector x1, w1;
			GetIntegrationRule(x1,w1,intorder_mass);
			
			// material parameters
			double rhoA = GetMaterial().BeamRhoA();
			double rhoIx = GetMaterial().BeamRhoIx();

			for (int i1=1; i1<=x1.GetLen(); i1++)
			{
				double x0 = x1(i1);    //in [-1,1]
				double xi = x0*size.X()*0.5; //in [-lx/2, lx/2]
#ifndef OLD_SHAPEFUNCTIONS 
				GetShapesPos(SVpos, x0);
#else
				GetShapesPos(SVpos, xi);
#endif
				GetShapesRot(SVrot, xi);

				// fill in matrix HL
				for (int i=1; i<=nspos; i++)
				{
					for (int j=1; j<=dim; j++)
					{
						HL((i-1)*dim+j,j) = SVpos(i);
					}
				}
				for (int i=1; i<=nsrot; i++)
				{
					HL(nspos*dim+i,4) = SVrot(i);
				}

				// jacobi determinant
				Vector3D r0x = GetInitPosx3D(xi); //this should be the initial and not the current configuration!
				double rx0n = r0x.Norm();
				rx0n = 1; //
				double det = 0.5*size.X()*rx0n;
				// add to mass matrix
				m += det * w1(i1) * rhoA * (HL*HL.GetTp());
			}

			for (int i=1; i<=nsrot; i++)
			{
				for (int j=1; j<=nsrot; j++)
				{
					m(dim*nspos+i,dim*nspos+j) *= rhoIx/rhoA;
				}
			}

			massmatrix = m;
		}
	}
	else if (kinematic_computation_mode == 0 || kinematic_computation_mode == 1) // kinematic_computation_mode:  0 ... exact terms, 5th order gaussian integration (slow); 
	{
		m.SetAll(0);

		ConstMatrix<14*14> m_at_ip(14,14,1); //storage for mass matrix at integration point  (nofill)
		double det = 0.5*size.X();
		Vector x1, w1;
		GetIntegrationRule(x1,w1,intorder_mass);
		
		for (int i1=1; i1<=x1.GetLen(); i1++)
		{
			double xi = x1(i1)*det;
			GetM(xi, m_at_ip);  //get mass matrix at integration point
			m += (det*w1(i1))*m_at_ip;
		}
	}
	//else  // commented out: too inaccurate approximation, above exact evaluation of GetM (at IP) is not CPU-time consuming compared to QVV     kinematic_computation_mode: 1 ... exact terms, low order (1st order lobatto) integration (fast); 
	//{
	//	m.SetAll(0);
	//	double det = 0.5*size.X();
	//	double det_rho_A = det*GetMaterial().BeamRhoA();
	//	{
	//		// evaluation at left node
	//		TArray<double> de;
	//		Getdeidposx_and_deidtheta(de,XG(4)+q0(4),XG(5)+q0(5),XG(6)+q0(6),XG(13)+q0(13),XData(1),XData(2),XData(3));
	//		Matrix3D de2dposx(de(1),de(4),de(7),de(2),de(5),de(8),de(3),de(6),de(9));
	//		Matrix3D de3dposx(de(13),de(16),de(19),de(14),de(17),de(20),de(15),de(18),de(21));
	//		Vector3D de2dtheta(de(10),de(11),de(12));
	//		Vector3D de3dtheta(de(22),de(23),de(24));

	//		Matrix3D m2 = GetMaterial().BeamRhoIz()*de2dposx.GetTp()*de2dposx + GetMaterial().BeamRhoIy()*de3dposx.GetTp()*de3dposx;
	//		for (int i=1; i<=3 ; i++)
	//		{
	//			m(i,i) = det_rho_A*1.;
	//			for (int j=1; j<=3 ; j++)
	//			{
	//				m(3+i,3+j) = det*m2(i,j);
	//			}
	//		}
	//		m(13,13) = det*(GetMaterial().BeamRhoIz()*de2dtheta*de2dtheta + GetMaterial().BeamRhoIy()*de3dtheta*de3dtheta);
	//	}
	//	{
	//		// evaluation at right node
	//		TArray<double> de;
	//		Getdeidposx_and_deidtheta(de,XG(10)+q0(10),XG(11)+q0(11),XG(12)+q0(12),XG(14)+q0(14),XData(4),XData(5),XData(6));
	//		Matrix3D de2dposx(de(1),de(4),de(7),de(2),de(5),de(8),de(3),de(6),de(9));
	//		Matrix3D de3dposx(de(13),de(16),de(19),de(14),de(17),de(20),de(15),de(18),de(21));
	//		Vector3D de2dtheta(de(10),de(11),de(12));
	//		Vector3D de3dtheta(de(22),de(23),de(24));

	//		Matrix3D m2 = GetMaterial().BeamRhoIz()*de2dposx.GetTp()*de2dposx + GetMaterial().BeamRhoIy()*de3dposx.GetTp()*de3dposx;
	//		for (int i=1; i<=3 ; i++)
	//		{
	//			m(6+i,6+i) = det_rho_A*1.;
	//			for (int j=1; j<=3 ; j++)
	//			{
	//				m(9+i,9+j) = det*m2(i,j);
	//			}
	//		}
	//		m(14,14) = det*(GetMaterial().BeamRhoIz()*de2dtheta*de2dtheta + GetMaterial().BeamRhoIy()*de3dtheta*de3dtheta);
	//	}
	//}
};


double ANCFBeam3DTorsion::PostNewtonStep(double t)
{	
	if (do_update_directors)
	{
		Vector3D r1dx, r2dx;
		double theta1, theta2;

		r1dx.X() = q0(4) + XG(4);
		r1dx.Y() = q0(5) + XG(5);
		r1dx.Z() = q0(6) + XG(6);
		r2dx.X() = q0(10) + XG(10);
		r2dx.Y() = q0(11) + XG(11);
		r2dx.Z() = q0(12) + XG(12);

		theta1 = q0(13) + XG(13);
		theta2 = q0(14) + XG(14);

		int err = UpdateDirectors(r1dx, r2dx, theta1, theta2); // returns 1 in case update was not possible (old director colinear with slope), and 0 else
		
		// return error for nonlinear iteration,
		// if the sum of these return value (over all elements) is larger than mbs->Discontinuous(), 
		// then the Discontinuous iteration step is repeated.
		// If this step is repeated a several times, then the time/load increment is reduced by a factor,
		// see e.g. SolverOptions.Static.load_inc_up or SolverOptions.Static.load_inc_down 
		return (mbs->DiscontinuousAccuracy() + 1.)*err; 
	}

	return 0.;
}



//for visualization
Vector3D ANCFBeam3DTorsion::GetPos3D0D(const Vector3D& p_loc) const 
{
	Vector3D plocscaled;
	plocscaled(1)=p_loc(1)*size.X()/2.;
	plocscaled(2)=p_loc(2)*size.Y()/2.;
	plocscaled(3)=p_loc(3)*size.Z()/2.;
	return GetPos(plocscaled(1),1) + GetRotMatrix(plocscaled(1), 1)*Vector3D(0., plocscaled.Y(), plocscaled.Z());
};

Vector3D ANCFBeam3DTorsion::GetPos3D0D(const Vector3D& p_loc, double defscale) const 
{
	Vector3D plocscaled;
	plocscaled(1)=p_loc(1)*size.X()/2.;
	plocscaled(2)=p_loc(2)*size.Y()/2.;
	plocscaled(3)=p_loc(3)*size.Z()/2.;
	return GetPos(plocscaled(1),1) + GetRotMatrix(plocscaled(1), 1)*Vector3D(0., plocscaled.Y(), plocscaled.Z()); 
};






static bool test_output_written = false;
void ANCFBeam3DTorsion::StartTimeStep()
{
	//int test_output_at_TIit = 1;

	//if(mbs->GetTIit() == test_output_at_TIit)
	//{
	//	if (GetNode(1).NodeNum() == 2 || GetNode(1).NodeNum() == 8 || GetNode(1).NodeNum() == 14)
	//	{
	//		ConstMatrix<14*14> m(14,14,1);
	//		ConstVector<14> qvv(14,1);
	//		EvalM(m,GetMBS()->GetTime());
	//		EvalQVV(qvv,GetMBS()->GetTime());

	//		mbs->UO() << "\nmode " << kinematic_computation_mode << ":\n";
	//		mbs->UO() << "Mass Matrix in mode " << kinematic_computation_mode << ":\n";
	//		mbs->UO() << m << "\n";
	//		mbs->UO() << "Quadratic Velocity vector:\n";
	//		mbs->UO() << qvv << "\n";
	//		mbs->UO() << "Node Position:\n";
	//		mbs->UO() << GetNode(1).Pos() << "\n";
	//	}
	//}

	//		if (kinematic_computation_mode == 0)
	//		{
	//			for (int n=1; n<=2; n++)
	//			{
	//				double xi;
	//				if (n==1)
	//				{
	//					mbs->UO() << "left node (mode 0) -------------------------- :\n";
	//					xi = -size.X()*0.5;
	//				}
	//				else
	//				{
	//					mbs->UO() << "right node (mode 0) -------------------------- :\n";
	//					xi = size.X()*0.5;
	//				}
	//				
	//				ConstMatrix<3*3> de2dposx(3,3,1);
	//				ConstMatrix<3*3> de3dposx(3,3,1);
	//				ConstVector<3> de3dtheta(3,1);
	//				ConstVector<3> de2dtheta(3,1);
	//				//ConstMatrix<3*3*3> dde2dposxdposx(3,3*3,1);
	//				//ConstMatrix<3*3*3> dde3dposxdposx(3,3*3,1);
	//				ConstMatrix<3*3> dde2dthetadposx(3,3,1);
	//				ConstMatrix<3*3> dde3dthetadposx(3,3,1);
	//				ConstVector<3> dde2dthetadtheta(3,1);
	//				ConstVector<3> dde3dthetadtheta(3,1);

	//				
	//				Getdeidposx(xi, 2, de2dposx);
	//				Getdeidposx(xi, 3, de3dposx);
	//				Getdeidrot(xi, 2, de2dtheta);
	//				Getdeidrot(xi, 3, de3dtheta);
	//				//Getddeidposxdposx(xi, 2, dde2dposxdposx);
	//				//Getddeidposxdposx(xi, 3, dde3dposxdposx);
	//				Getddeidposxdrot(xi,2,dde2dthetadposx);
	//				Getddeidposxdrot(xi,3,dde3dthetadposx);
	//				Getddeidrotdrot(xi,2,dde2dthetadtheta);
	//				Getddeidrotdrot(xi,3,dde3dthetadtheta);
	//				
	//				Matrix3D Dde2dposx;
	//				Matrix3D Dde3dposx;
	//				Vector3D Dde2dtheta;
	//				Vector3D Dde3dtheta;
	//				Vector3D Ddde2dthetadtheta;
	//				Vector3D Ddde3dthetadtheta;
	//				Matrix3D Ddde2dthetadposx;
	//				Matrix3D Ddde3dthetadposx;
	//				for (int i=1; i<=3; i++)
	//				{
	//					Ddde2dthetadtheta(i)=dde2dthetadtheta(i);
	//					Ddde3dthetadtheta(i)=dde3dthetadtheta(i);

	//					Dde2dtheta(i) = de2dtheta(i);
	//					Dde3dtheta(i) = de3dtheta(i);

	//					for (int j=1; j<=3; j++)
	//					{
	//						Ddde2dthetadposx(i,j)=dde2dthetadposx(i,j);
	//						Ddde3dthetadposx(i,j)=dde3dthetadposx(i,j);

	//						Dde2dposx(i,j) = de2dposx(i,j);
	//						Dde3dposx(i,j) = de3dposx(i,j);							
	//					}
	//				}

	//				//Matrix3D dde2dposxdposx1;
	//				//Matrix3D dde2dposxdposx2;
	//				//Matrix3D dde2dposxdposx3;
	//				//Matrix3D dde3dposxdposx1;
	//				//Matrix3D dde3dposxdposx2;
	//				//Matrix3D dde3dposxdposx3;

	//				for (int i=1; i<=3; i++)
	//				{
	//					for (int j=1; j<=3; j++)
	//					{
	//						//dde2dposxdposx1(i,j)=dde2dposxdposx(i,j);
	//						//dde2dposxdposx2(i,j)=dde2dposxdposx(i,j+3);
	//						//dde2dposxdposx3(i,j)=dde2dposxdposx(i,j+6);
	//						//dde3dposxdposx1(i,j)=dde3dposxdposx(i,j);
	//						//dde3dposxdposx2(i,j)=dde3dposxdposx(i,j+3);
	//						//dde3dposxdposx3(i,j)=dde3dposxdposx(i,j+6);

	//					}
	//				}
	//				
	//				mbs->UO() << "de2dposx =\n" << Dde2dposx;
	//				mbs->UO() << "de3dposx =\n" << Dde3dposx;
	//				mbs->UO() << "de2dtheta =\n" << Dde2dtheta << "\n";
	//				mbs->UO() << "de3dtheta =\n" << Dde3dtheta << "\n";
	//				//mbs->UO() << "dde2dposxdposx1 =\n" << dde2dposxdposx1;
	//				//mbs->UO() << "dde2dposxdposx2 =\n" << dde2dposxdposx2;
	//				//mbs->UO() << "dde2dposxdposx3 =\n" << dde2dposxdposx3;
	//				//mbs->UO() << "dde3dposxdposx1 =\n" << dde3dposxdposx1;
	//				//mbs->UO() << "dde3dposxdposx2 =\n" << dde3dposxdposx2;
	//				//mbs->UO() << "dde3dposxdposx3 =\n" << dde3dposxdposx3;
	//				mbs->UO() << "dde2dthetadposx =\n" << Ddde2dthetadposx;
	//				mbs->UO() << "dde3dthetadposx =\n" << Ddde3dthetadposx;
	//				mbs->UO() << "dde2dthetadtheta =\n" << Ddde2dthetadtheta << "\n";
	//				mbs->UO() << "dde3dthetadtheta =\n" << Ddde3dthetadtheta << "\n";

	//				/*for (int k=1; k<=14; k++)
	//				{
	//					ConstMatrix<9*14> dLdq
	//				}*/
	//				ConstMatrix<9*14> L0(9,14,1);
	//				GetL0(xi, L0);
	//				mbs->UO() << "L0 =\n" << L0;
	//			
	//				ConstMatrix<9*14> dL0dqk(9,14,1);
	//				GetdL0dqk(xi, 12+n, dL0dqk);
	//				mbs->UO() << "dL0dqk =\n" << dL0dqk;

	//				ConstVector<14> qdot(14);
	//				for (int k=1; k<=14; k++)
	//				{
	//					qdot(k)=XGP(k);
	//				}
	//				double Iz = GetMaterial().BeamRhoIz();
	//				double Iy = GetMaterial().BeamRhoIy();
	//				ConstMatrix<9*9> D(9,9);
	//				for (int i=1; i<=9; i++)
	//				{
	//					D(i,i) = (i<4)?GetMaterial().BeamRhoA():((i<7)?Iz:Iy);
	//				}
	//				ConstVector<9> L0_qdot(9); L0_qdot=D*L0*qdot;
	//				ConstVector<9> dL0dqk_qdot(9); dL0dqk_qdot=dL0dqk*qdot;				
	//				double L0_qdot_dL0dqk_qdot = 0;
	//				for (int a=1; a<=9; a++) L0_qdot_dL0dqk_qdot += L0_qdot(a)*dL0dqk_qdot(a);
	//				double det=0.5*size.X();
	//				double qvv12n = det*0.5*L0_qdot_dL0dqk_qdot;

	//				mbs->UO() << "qdot =\n" << qdot << "\n";
	//				mbs->UO() << "L0_qdot =\n" << L0_qdot << "\n";
	//				mbs->UO() << "Iz = " << Iz << ",  Iy = " << Iy << "\n";
	//				mbs->UO() << "D =\n" << D << "\n";
	//				mbs->UO() << "dL0dqk_qdot =\n" << dL0dqk_qdot << "\n";
	//				mbs->UO() << "L0_qdot_D_dL0dqk_qdot = " << L0_qdot_dL0dqk_qdot << "\n";
	//				mbs->UO() << "det = " << det << "\n";
	//				mbs->UO() << "qvv(" << 12+n << ") = " << qvv12n << "\n";
	//			}
	//		}
	//		else if (kinematic_computation_mode == 1)
	//		{
	//			for (int n=1; n<=2; n++)
	//			{
	//				if (n==1)
	//				{
	//					mbs->UO() << "left node (mode 1) -------------------------- :\n";
	//				}
	//				else 
	//				{
	//					mbs->UO() << "right node (mode 1) -------------------------- :\n";
	//				}
	//				
	//				int offset_ang = n-1;
	//				int offset_dir = 3*offset_ang;
	//				int offset_pos = 6*offset_ang;

	//				double Iy = GetMaterial().BeamRhoIy();
	//				double Iz = GetMaterial().BeamRhoIz();

	//				// evaluation at left node
	//				double drdx1 = XG(4+offset_pos)+q0(4+offset_pos);
	//				double drdx2 = XG(5+offset_pos)+q0(5+offset_pos);
	//				double drdx3 = XG(6+offset_pos)+q0(6+offset_pos);
	//				double theta = XG(13+offset_ang)+q0(13+offset_ang);
	//				double d1 = XData(1+offset_dir);
	//				double d2 = XData(2+offset_dir);
	//				double d3 = XData(3+offset_dir);
	//				Vector3D qdotpos(XGP(4+offset_pos),XGP(5+offset_pos),XGP(6+offset_pos));
	//				double qdottheta(XGP(13+offset_ang));

	//				TArray<double> de, de11, de21, de31, de22, de32, de33, de4, de44;

	//				Getdeidposx_and_deidtheta(de,drdx1,drdx2,drdx3,theta,d1,d2,d3);
	//				//Getddeidposx1dposx1(de11,drdx1,drdx2,drdx3,theta,d1,d2,d3);
	//				//Getddeidposx2dposx1(de21,drdx1,drdx2,drdx3,theta,d1,d2,d3);
	//				//Getddeidposx3dposx1(de31,drdx1,drdx2,drdx3,theta,d1,d2,d3);
	//				//Getddeidposx2dposx2(de22,drdx1,drdx2,drdx3,theta,d1,d2,d3);
	//				//Getddeidposx3dposx2(de32,drdx1,drdx2,drdx3,theta,d1,d2,d3);
	//				//Getddeidposx3dposx3(de33,drdx1,drdx2,drdx3,theta,d1,d2,d3);
	//				Getddeidthetadposx(de4,drdx1,drdx2,drdx3,theta,d1,d2,d3);
	//				Getddeidthetadtheta(de44,drdx1,drdx2,drdx3,theta,d1,d2,d3);

	//				Matrix3D de2dposx(de(1),de(4),de(7),de(2),de(5),de(8),de(3),de(6),de(9));
	//				Matrix3D de3dposx(de(13),de(16),de(19),de(14),de(17),de(20),de(15),de(18),de(21));
	//				Vector3D de2dtheta(de(10),de(11),de(12));
	//				Vector3D de3dtheta(de(22),de(23),de(24));
	//				//Matrix3D dde2dposxdposx1(de11(1),de21(1),de31(1),de11(2),de21(2),de31(2),de11(3),de21(3),de31(3));
	//				//Matrix3D dde3dposxdposx1(de11(4),de21(4),de31(4),de11(5),de21(5),de31(5),de11(6),de21(6),de31(6));
	//				//Matrix3D dde2dposxdposx2(de21(1),de22(1),de32(1),de21(2),de22(2),de32(2),de21(3),de22(3),de32(3));
	//				//Matrix3D dde3dposxdposx2(de21(4),de22(4),de32(4),de21(5),de22(5),de32(5),de21(6),de22(6),de32(6));
	//				//Matrix3D dde2dposxdposx3(de31(1),de32(1),de33(1),de31(2),de32(2),de33(2),de31(3),de32(3),de33(3));
	//				//Matrix3D dde3dposxdposx3(de31(4),de32(4),de33(4),de31(5),de32(5),de33(5),de31(6),de32(6),de33(6));
	//				Matrix3D dde2dthetadposx(de4(1),de4(4),de4(7),de4(2),de4(5),de4(8),de4(3),de4(6),de4(9));
	//				Matrix3D dde3dthetadposx(de4(10),de4(13),de4(16),de4(11),de4(14),de4(17),de4(12),de4(15),de4(18));
	//				Vector3D dde2dthetadposx1(de4(1),de4(2),de4(3));
	//				Vector3D dde3dthetadposx1(de4(10),de4(11),de4(12));
	//				Vector3D dde2dthetadposx2(de4(4),de4(5),de4(6));
	//				Vector3D dde3dthetadposx2(de4(13),de4(14),de4(15));
	//				Vector3D dde2dthetadposx3(de4(7),de4(8),de4(9));
	//				Vector3D dde3dthetadposx3(de4(16),de4(17),de4(18));
	//				Vector3D dde2dthetadtheta(de44(1),de44(2),de44(3));
	//				Vector3D dde3dthetadtheta(de44(4),de44(5),de44(6));

	//				mbs->UO() << "de2dposx =\n" << de2dposx;
	//				mbs->UO() << "de3dposx =\n" << de3dposx;
	//				mbs->UO() << "de2dtheta =\n" << de2dtheta << "\n";
	//				mbs->UO() << "de3dtheta =\n" << de3dtheta << "\n";
	//				//mbs->UO() << "dde2dposxdposx1 =\n" << dde2dposxdposx1;
	//				//mbs->UO() << "dde2dposxdposx2 =\n" << dde2dposxdposx2;
	//				//mbs->UO() << "dde2dposxdposx3 =\n" << dde2dposxdposx3;
	//				//mbs->UO() << "dde3dposxdposx1 =\n" << dde3dposxdposx1;
	//				//mbs->UO() << "dde3dposxdposx2 =\n" << dde3dposxdposx2;
	//				//mbs->UO() << "dde3dposxdposx3 =\n" << dde3dposxdposx3;
	//				mbs->UO() << "dde2dthetadposx =\n" << dde2dthetadposx;
	//				mbs->UO() << "dde3dthetadposx =\n" << dde3dthetadposx;
	//				mbs->UO() << "dde2dthetadtheta =\n" << dde2dthetadtheta << "\n";
	//				mbs->UO() << "dde3dthetadtheta =\n" << dde3dthetadtheta << "\n";

	//				Vector3D qdot_de2dposx = de2dposx*(Iz*qdotpos);
	//				Vector3D qdot_de2dtheta = (Iz*qdottheta)*de2dtheta;
	//				Vector3D qdot_de3dposx = de3dposx*(Iy*qdotpos);
	//				Vector3D qdot_de3dtheta = (Iy*qdottheta)*de3dtheta;

	//				Vector3D dde2dthetadtheta_qdot = dde2dthetadtheta*qdottheta;
	//				Vector3D dde3dthetadtheta_qdot = dde3dthetadtheta*qdottheta;
	//				Vector3D dde2dposxdtheta_qdot = dde2dthetadposx*qdotpos;
	//				Vector3D dde3dposxdtheta_qdot = dde3dthetadposx*qdotpos;


	//				mbs->UO() << "qdotpos = " << qdotpos << "\n";
	//				mbs->UO() << "qdottheta = " << qdottheta << "\n";
	//				mbs->UO() << "dde2dposxdtheta_qdot =\n" << dde2dposxdtheta_qdot << "\n";
	//				mbs->UO() << "dde3dposxdtheta_qdot =\n" << dde3dposxdtheta_qdot << "\n";
	//				mbs->UO() << "dde2dthetadtheta_qdot =\n" << dde2dthetadtheta_qdot << "\n";
	//				mbs->UO() << "dde3dthetadtheta_qdot =\n" << dde3dthetadtheta_qdot << "\n";
	//				mbs->UO() << "dL02dq_qdot =\n" << dde2dposxdtheta_qdot+dde2dthetadtheta_qdot << "\n";
	//				mbs->UO() << "dL03dq_qdot =\n" << dde3dposxdtheta_qdot+dde3dthetadtheta_qdot << "\n";
	//				mbs->UO() << "Iz = " << Iz << ",  Iy = " << Iy << "\n";
	//				mbs->UO() << "qdot_de2dposx = de2dposx*(Iz*qdotpos) =\n" << qdot_de2dposx << "\n";
	//				mbs->UO() << "qdot_de2dtheta = (Iz*qdottheta)*de2dtheta =\n" << qdot_de2dtheta << "\n";
	//				mbs->UO() << "L02_qdot =\n" << qdot_de2dposx+qdot_de2dtheta << "\n";
	//				mbs->UO() << "qdot_de3dposx = de3dposx*(Iy*qdotpos) =\n" << qdot_de3dposx << "\n";
	//				mbs->UO() << "qdot_de3dtheta = (Iy*qdottheta)*de3dtheta =\n" << qdot_de3dtheta << "\n";
	//				mbs->UO() << "L03_qdot =\n" << qdot_de3dposx+qdot_de3dtheta << "\n";

	//				mbs->UO() << "fac = " << 0.25*size.X() << "\n";
	//				

	//				//double ftheta = (qdot_de2dposx+qdot_de2dtheta)*(dde2dposxdtheta_qdot+dde2dthetadtheta_qdot) + (qdot_de3dposx+qdot_de3dtheta)*(dde3dposxdtheta_qdot+dde3dthetadtheta_qdot);

	//			}
	//		}
	//	}
	//}
	//		
	//if(mbs->GetTIit() > test_output_at_TIit)
	//{
	//	GetMBS()->StopCalculation();
	//}


//	int test_output_at_TIit = mbs->GetModelDataContainer()->TreeGetInt("ANCFBeam3DTorsionTest.test_output_at_TIit", 0);
//
//	if (test_output_at_TIit < 0)
//	{
//		return;
//	}
//	
//	if(mbs->GetTIit() > test_output_at_TIit)
//	{
//		GetMBS()->GetSolSet().endtime = GetMBS()->GetTime();
//	}
//
//	if(mbs->GetTIit() == test_output_at_TIit && !test_output_written)
//	{
//		double xi = 0.38729833;
//		double x0 = xi*2./size.X(); 
//
//		Vector xg(SOS()); for (int i=1; i<=SOS(); i++) xg(i) = XG(i);
//		//--------------------------------------------------------------
//		ConstMatrix<9*14> dL0dq1(9,14,1); GetdL0dqk(xi, 1, dL0dq1);
//		ConstMatrix<3*3*3> dde2dposxdposx(3,3*3,1); Getddeidposxdposx(xi, 2, dde2dposxdposx);
//		ConstMatrix<3*3*3> dde3dposxdposx(3,3*3,1); Getddeidposxdposx(xi, 3, dde3dposxdposx);
//		ConstMatrix<3*3> dde2dposxdrot(3,3,1); Getddeidposxdrot(xi, 2, dde2dposxdrot);
//		ConstMatrix<3*3> dde3dposxdrot(3,3,1); Getddeidposxdrot(xi, 3, dde3dposxdrot);
//		ConstVector<3> dde2drotdrot(3,1); Getddeidrotdrot(xi, 2, dde2drotdrot);
//		ConstVector<3> dde3drotdrot(3,1); Getddeidrotdrot(xi, 3, dde3drotdrot);
//		ConstMatrix<3*3*3> dde20dposxdposx(3,3*3,1); Getddei0dposxdposx(xi, 2, dde20dposxdposx);
//		ConstMatrix<3*3*3> dde30dposxdposx(3,3*3,1); Getddei0dposxdposx(xi, 3, dde30dposxdposx);
//		ConstMatrix<3*3*3> dde20hatdposxdposx(3,3*3,1); Getddei0hatdposxdposx(xi, 2, dde20hatdposxdposx);
//		ConstMatrix<3*3*3> dde30hatdposxdposx(3,3*3,1); Getddei0hatdposxdposx(xi, 3, dde30hatdposxdposx);
//		ConstMatrix<3*3> de20hatdposx(3,3,1); Getdei0hatdposx(xi, 2, de20hatdposx);
//		ConstMatrix<3*3> de30hatdposx(3,3,1); Getdei0hatdposx(xi, 3, de30hatdposx);
//		ConstVector<3> e20hat(3,1); Getei0hat(xi, 2, e20hat);	
//		ConstVector<3> e30hat(3,1); Getei0hat(xi, 3, e30hat);	
//		ConstMatrix<3*3*3> dde30hatde1de1(3,3*3,1); Getdde30hatde1de1(xi, dde30hatde1de1);
//		ConstMatrix<3*3*3> dde1dposxdposx(3,3*3,1); Getdde1dposxdposx(xi, dde1dposxdposx);
//		ConstMatrix<3*3> de1dposx(3,3,1);	Getde1dposx(xi, de1dposx);
//		ConstMatrix<3*3> de30hatde1(3,3,1); Getde30hatde1(xi, de30hatde1);
//		//--------------------------------------------------------------
//		ConstMatrix<9*14> L0(9,14,1); GetL0(xi, L0);
//		ConstMatrix<3*14> de2dq(3, 14, 1); Getdeidq(xi, 2, de2dq);
//		ConstMatrix<3*14> de3dq(3, 14, 1); Getdeidq(xi, 3, de3dq);
//		ConstMatrix<3*3> de2dposx(3,3,1); Getdeidposx(xi, 2, de2dposx);
//		ConstMatrix<3*3> de3dposx(3,3,1); Getdeidposx(xi, 3, de3dposx);
//		ConstMatrix<3*3> de20dposx(3,3,1); Getdei0dposx(xi, 2, de20dposx);
//		ConstMatrix<3*3> de30dposx(3,3,1); Getdei0dposx(xi, 3, de30dposx);
//		ConstVector<3> de2drot(3,1); Getdeidrot(xi, 2, de2drot);
//		ConstVector<3> de3drot(3,1); Getdeidrot(xi, 3, de3drot);
//		Vector3D posx = GetPosx(xi);
//		double theta = GetRot(xi);
//		ConstVector<4> sfposx(4,1); for(int j=1; j<=4; j++) sfposx(j) = GetSFPosx(j,x0);
//		ConstVector<2> sfrot(2,1); for(int j=1; j<=2; j++) sfrot(j) = GetSFRot(j,xi);
//
//		mbs->UO(UO_LVL_0) << "--------------------------------------------------------------------------------------------------------\n";
//		mbs->UO(UO_LVL_0) << "Test Output at step " << mbs->GetTIit() << "(" << mbs->GetTime() << " s):\n";
//		mbs->UO(UO_LVL_0) << "xi=" << xi << ", director = " << GetDirector(xi) << "\n";
//		mbs->UO(UO_LVL_0) << "XG()=\n" << xg << "\n" << "q0 = \n" << q0 << "\n" << "XG()+q0=\n" << xg+q0 << "\n";
//		mbs->UO(UO_LVL_0) << "--------------------------------------------------------------------------------------------------------\n";
//		mbs->UO(UO_LVL_0) << "dL0dq1=\n" << dL0dq1;
//		mbs->UO(UO_LVL_0) << "dde2dposxdposx=\n" << dde2dposxdposx;
//		mbs->UO(UO_LVL_0) << "dde3dposxdposx=\n" << dde3dposxdposx;
//		mbs->UO(UO_LVL_0) << "dde2dposxdrot=\n" << dde2dposxdrot;
//		mbs->UO(UO_LVL_0) << "dde3dposxdrot=\n" << dde3dposxdrot;
//		mbs->UO(UO_LVL_0) << "dde2drotdrot=\n" << dde2drotdrot << "\n";
//		mbs->UO(UO_LVL_0) << "dde3drotdrot=\n" << dde3drotdrot << "\n";
//		mbs->UO(UO_LVL_0) << "dde20dposxdposx=\n" << dde20dposxdposx;
//		mbs->UO(UO_LVL_0) << "dde30dposxdposx=\n" << dde30dposxdposx;
//		mbs->UO(UO_LVL_0) << "dde20hatdposxdposx=\n" << dde20hatdposxdposx;
//		mbs->UO(UO_LVL_0) << "dde30hatdposxdposx=\n" << dde30hatdposxdposx;
//		mbs->UO(UO_LVL_0) << "de20hatdposx=\n" << de20hatdposx;
//		mbs->UO(UO_LVL_0) << "de30hatdposx=\n" << de30hatdposx;
//		mbs->UO(UO_LVL_0) << "e20hat=\n" << e20hat << "\n";
//		mbs->UO(UO_LVL_0) << "e30hat=\n" << e30hat << "\n";
//		mbs->UO(UO_LVL_0) << "dde30hatde1de1=\n" << dde30hatde1de1;
//		mbs->UO(UO_LVL_0) << "dde1dposxdposx=\n" << dde1dposxdposx;
//		mbs->UO(UO_LVL_0) << "de1dposx=\n" << de1dposx;
//		mbs->UO(UO_LVL_0) << "de30hatde1=\n" << de30hatde1;
//		mbs->UO(UO_LVL_0) << "--------------------------------------------------------------------------------------------------------\n";
//		mbs->UO(UO_LVL_0) << "L0=\n" << L0;
//		mbs->UO(UO_LVL_0) << "de2dq=\n" << de2dq.PrintToMathematica() << "\n";
//		mbs->UO(UO_LVL_0) << "de3dq=\n" << de3dq.PrintToMathematica() << "\n";
//		mbs->UO(UO_LVL_0) << "de2dposx=\n" << de2dposx.PrintToMathematica() << "\n";
//		mbs->UO(UO_LVL_0) << "de3dposx=\n" << de3dposx.PrintToMathematica() << "\n";
//		mbs->UO(UO_LVL_0) << "de20dposx=\n" << de20dposx.PrintToMathematica() << "\n";
//		mbs->UO(UO_LVL_0) << "de30dposx=\n" << de30dposx.PrintToMathematica() << "\n";
//		mbs->UO(UO_LVL_0) << "Cos(theta)=" << Cos(theta) << "\n";
//		mbs->UO(UO_LVL_0) << "Sin(theta)=" << Sin(theta) << "\n";
//		mbs->UO(UO_LVL_0) << "de2drot=" << de2drot << "\n";
//		mbs->UO(UO_LVL_0) << "de3drot=" << de3drot << "\n";
//		mbs->UO(UO_LVL_0) << "posx=" << posx << "\n";
//		mbs->UO(UO_LVL_0) << "theta=" << theta << "\n";
//		mbs->UO(UO_LVL_0) << "sfposx=" << sfposx << "\n";
//		mbs->UO(UO_LVL_0) << "sfrot=" << sfrot << "\n";
//		mbs->UO(UO_LVL_0) << "--------------------------------------------------------------------------------------------------------\n";
//
//		test_output_written = true;
//	}
}


// standard version:
void ANCFBeam3DTorsion::GetQuadraticVelocityVector(const double& xi, Vector& qvv) const //eqs.(47)-(48), xi in [-lx/2, lx/2]
{
	qvv.SetAll(0);
	ConstMatrix<14*14> dMdqk(14, 14, 1); //nofill
	ConstVector<14> qdot(14, 1); //nofill

	for (int k=1; k<=14; k++)
	{
		qdot(k)=XGP(k);
	}
	for (int k=1; k<=14; k++)
	{
		GetdMdqk(xi, k, dMdqk);
		qvv(k) += 0.5*(qdot*dMdqk*qdot);
		//qvv(k) -= qdot*dMdqk*qdot;    //eq.(47)    // both lines are equal!   and summation wrong...  right would be:  (48) - 0.5 * (47)       now, replace by    qvv = -0.5 (47)      which is faster than (48)!
		//qvv += qdot(k)*(qdot*dMdqk);  //eq.(48)
	}
}

// used for 2nd order lobatto integration (IPs are in nodes n=1 and n=2 (local indexing)):
void ANCFBeam3DTorsion::GetQuadraticVelocityVectorAtNode(Vector& f, int n) const 
{
	// for speedup of dynamic computation we assume 2nd order integration rule by Lobatto (integration points only at nodes),
	// then:  1.) dde2dposxdposx, dde3dposxdposx, dde2dposxdtheta, dde3dposxdtheta, dde2dthetadtheta, dde3dthetadtheta are calculated at nodes only, s.t. generalized coordinates may be used as input directly (instead of time consuming shape function evaluations).      2.) dL0dq(x0=-1 or x0=1) = (0, dde2dqdq, dde3dqdq) is sparse.

	f.SetAll(0.);

	int offset_ang = n-1;
	int offset_dir = 3*offset_ang;
	int offset_pos = 6*offset_ang;

	double Iy = GetMaterial().BeamRhoIy();
	double Iz = GetMaterial().BeamRhoIz();

	// evaluation at left node
	double drdx1 = XG(4+offset_pos)+q0(4+offset_pos);
	double drdx2 = XG(5+offset_pos)+q0(5+offset_pos);
	double drdx3 = XG(6+offset_pos)+q0(6+offset_pos);
	double theta = XG(13+offset_ang)+q0(13+offset_ang);
	double d1 = XData(1+offset_dir);
	double d2 = XData(2+offset_dir);
	double d3 = XData(3+offset_dir);
	Vector3D qdotpos(XGP(4+offset_pos),XGP(5+offset_pos),XGP(6+offset_pos));
	double qdottheta(XGP(13+offset_ang));

	TArray<double> de, de11, de21, de31, de22, de32, de33, de4; //de44;

	{ // this part takes 3/4 CPU time of the whole function GetQuadraticVelocityVectorAtNode, thus a compound calculation (not split into 6 functions, where almost the same calculations are done) would be advisable.
		// however, the kinematic terms in total just need the third of CPU-time compared to the calculations on the work of interior forces... 
		Getdeidposx_and_deidtheta(de,drdx1,drdx2,drdx3,theta,d1,d2,d3);
		Getddeidposx1dposx1(de11,drdx1,drdx2,drdx3,theta,d1,d2,d3);
		Getddeidposx2dposx1(de21,drdx1,drdx2,drdx3,theta,d1,d2,d3);
		Getddeidposx3dposx1(de31,drdx1,drdx2,drdx3,theta,d1,d2,d3);
		Getddeidposx2dposx2(de22,drdx1,drdx2,drdx3,theta,d1,d2,d3);
		Getddeidposx3dposx2(de32,drdx1,drdx2,drdx3,theta,d1,d2,d3);
		Getddeidposx3dposx3(de33,drdx1,drdx2,drdx3,theta,d1,d2,d3);
		Getddeidthetadposx(de4,drdx1,drdx2,drdx3,theta,d1,d2,d3);
		//Getddeidthetadtheta(de44,drdx1,drdx2,drdx3,theta,d1,d2,d3);
	}

	Matrix3D de2dposx(de(1),de(4),de(7),de(2),de(5),de(8),de(3),de(6),de(9));
	Matrix3D de3dposx(de(13),de(16),de(19),de(14),de(17),de(20),de(15),de(18),de(21));
	Vector3D de2dtheta(de(10),de(11),de(12));
	Vector3D de3dtheta(de(22),de(23),de(24));
	Matrix3D dde2dposxdposx1(de11(1),de21(1),de31(1),de11(2),de21(2),de31(2),de11(3),de21(3),de31(3));
	Matrix3D dde3dposxdposx1(de11(4),de21(4),de31(4),de11(5),de21(5),de31(5),de11(6),de21(6),de31(6));
	Matrix3D dde2dposxdposx2(de21(1),de22(1),de32(1),de21(2),de22(2),de32(2),de21(3),de22(3),de32(3));
	Matrix3D dde3dposxdposx2(de21(4),de22(4),de32(4),de21(5),de22(5),de32(5),de21(6),de22(6),de32(6));
	Matrix3D dde2dposxdposx3(de31(1),de32(1),de33(1),de31(2),de32(2),de33(2),de31(3),de32(3),de33(3));
	Matrix3D dde3dposxdposx3(de31(4),de32(4),de33(4),de31(5),de32(5),de33(5),de31(6),de32(6),de33(6));
	//Matrix3D dde2dthetadposx(de4(1),de4(4),de4(7),de4(2),de4(5),de4(8),de4(3),de4(6),de4(9));
	//Matrix3D dde3dthetadposx(de4(10),de4(13),de4(16),de4(11),de4(14),de4(17),de4(12),de4(15),de4(18));
	Vector3D dde2dthetadposx1(de4(1),de4(2),de4(3));
	Vector3D dde3dthetadposx1(de4(10),de4(11),de4(12));
	Vector3D dde2dthetadposx2(de4(4),de4(5),de4(6));
	Vector3D dde3dthetadposx2(de4(13),de4(14),de4(15));
	Vector3D dde2dthetadposx3(de4(7),de4(8),de4(9));
	Vector3D dde3dthetadposx3(de4(16),de4(17),de4(18));
	//Vector3D dde2dthetadtheta(de44(1),de44(2),de44(3));
	//Vector3D dde3dthetadtheta(de44(4),de44(5),de44(6));

	Vector3D qdot_de2dposx = de2dposx*(Iz*qdotpos);
	Vector3D qdot_de2dtheta = (Iz*qdottheta)*de2dtheta;

	Vector3D qdot_de3dposx = de3dposx*(Iy*qdotpos);
	Vector3D qdot_de3dtheta = (Iy*qdottheta)*de3dtheta;

	Vector3D dde2dposxdposx1_qdot = dde2dposxdposx1*qdotpos;
	Vector3D dde2dposxdposx2_qdot = dde2dposxdposx2*qdotpos;
	Vector3D dde2dposxdposx3_qdot = dde2dposxdposx3*qdotpos;
	Vector3D dde2dthetadposx1_qdot = dde2dthetadposx1*qdottheta;
	Vector3D dde2dthetadposx2_qdot = dde2dthetadposx2*qdottheta;
	Vector3D dde2dthetadposx3_qdot = dde2dthetadposx3*qdottheta;
	//Vector3D dde2dposxdtheta_qdot = dde2dthetadposx*qdotpos;
	//Vector3D dde2dthetadtheta_qdot = dde2dthetadtheta*qdottheta;

	Vector3D dde3dposxdposx1_qdot = dde3dposxdposx1*qdotpos;
	Vector3D dde3dposxdposx2_qdot = dde3dposxdposx2*qdotpos;
	Vector3D dde3dposxdposx3_qdot = dde3dposxdposx3*qdotpos;
	Vector3D dde3dthetadposx1_qdot = dde3dthetadposx1*qdottheta;
	Vector3D dde3dthetadposx2_qdot = dde3dthetadposx2*qdottheta;
	Vector3D dde3dthetadposx3_qdot = dde3dthetadposx3*qdottheta;
	//Vector3D dde3dposxdtheta_qdot = dde3dthetadposx*qdotpos;
	//Vector3D dde3dthetadtheta_qdot = dde3dthetadtheta*qdottheta;

	Vector3D fposx(
		(qdot_de2dposx+qdot_de2dtheta)*(dde2dposxdposx1_qdot+dde2dthetadposx1_qdot) + (qdot_de3dposx+qdot_de3dtheta)*(dde3dposxdposx1_qdot+dde3dthetadposx1_qdot),
		(qdot_de2dposx+qdot_de2dtheta)*(dde2dposxdposx2_qdot+dde2dthetadposx2_qdot) + (qdot_de3dposx+qdot_de3dtheta)*(dde3dposxdposx2_qdot+dde3dthetadposx2_qdot),
		(qdot_de2dposx+qdot_de2dtheta)*(dde2dposxdposx3_qdot+dde2dthetadposx3_qdot) + (qdot_de3dposx+qdot_de3dtheta)*(dde3dposxdposx3_qdot+dde3dthetadposx3_qdot));
	//double ftheta = (qdot_de2dposx+qdot_de2dtheta)*(dde2dposxdtheta_qdot+dde2dthetadtheta_qdot) + (qdot_de3dposx+qdot_de3dtheta)*(dde3dposxdtheta_qdot+dde3dthetadtheta_qdot);   // equals zero!!!!

	double fac = 0.25*size.X();   // = 0.5*det,  det = 0.5*lx;
	f(4+offset_pos) = fac*fposx(1);
	f(5+offset_pos) = fac*fposx(2);
	f(6+offset_pos) = fac*fposx(3);
	//f(13+offset_ang) = fac*ftheta;    // always 0, since d(deidt^2)/dq13=0 at xi=-0.5*lx   and d(deidt^2)/dq14=0 at xi=0.5*lx
}

void ANCFBeam3DTorsion::GetM(const double& xi, Matrix& m) const //eq.(49), xi in [-lx/2, lx/2]
{
	m.SetAll(0);
	ConstMatrix<9*14> L0(9, 14, 1); //nofill
	GetL0(xi, L0);

	double sqrt_rhoA = sqrt(GetMaterial().BeamRhoA());
	double sqrt_rhoIy = sqrt(GetMaterial().BeamRhoIy());
	double sqrt_rhoIz = sqrt(GetMaterial().BeamRhoIz());

	L0.MulRow(1, sqrt_rhoA);
	L0.MulRow(2, sqrt_rhoA);
	L0.MulRow(3, sqrt_rhoA);
	L0.MulRow(4, sqrt_rhoIz);
	L0.MulRow(5, sqrt_rhoIz);
	L0.MulRow(6, sqrt_rhoIz);
	L0.MulRow(7, sqrt_rhoIy);
	L0.MulRow(8, sqrt_rhoIy);
	L0.MulRow(9, sqrt_rhoIy);

	MultTp(L0,L0,m);
}

void ANCFBeam3DTorsion::GetdMdqk(const double& xi, const int k, Matrix& dMdqk) const //eq.(51), xi in [-lx/2, lx/2], k in {1,..,14}
{
	dMdqk.SetAll(0);

	ConstMatrix<9*14> L0(9, 14, 1); //nofill
	ConstMatrix<9*14> dL0dqk(9, 14, 1); //nofill
	GetL0(xi, L0);
	GetdL0dqk(xi, k, dL0dqk);

	L0.MulRow(1, GetMaterial().BeamRhoA());
	L0.MulRow(2, GetMaterial().BeamRhoA());
	L0.MulRow(3, GetMaterial().BeamRhoA());
	L0.MulRow(4, GetMaterial().BeamRhoIz());
	L0.MulRow(5, GetMaterial().BeamRhoIz());
	L0.MulRow(6, GetMaterial().BeamRhoIz());
	L0.MulRow(7, GetMaterial().BeamRhoIy());
	L0.MulRow(8, GetMaterial().BeamRhoIy());
	L0.MulRow(9, GetMaterial().BeamRhoIy());

	MultTp(L0, dL0dqk, dMdqk);
}

void ANCFBeam3DTorsion::GetL0(const double& xi, Matrix& L0) const //eq.(42), xi in [-lx/2, lx/2]
{
	double x0 = xi*2./size.X();

	if(0)
	{
		L0.SetAll(0);

		ConstMatrix<3*14> deidq(3, 14, 1); //nofill

		Getdeidq(xi, 2, deidq);
		L0.SetSubmatrix(deidq, 3+1, 1); 

		Getdeidq(xi, 3, deidq);
		L0.SetSubmatrix(deidq, 2*3+1, 1);

		int offset = 0;
		for (int i = 1; i<=NSPos(); i++)
		{
			double sf_val = GetSFPos(i, x0);
			for (int j = 1; j<=Dim(); j++)
			{
				L0(j, j+offset) = sf_val;
			}
			offset += Dim();
		}
	}
	else   //faster
	{
		TArray<double> de;
		Vector3D d = GetDirector(xi);
		Vector3D drdx = GetPosx(xi);
		double theta = GetRot(xi);
		Getdeidposx_and_deidtheta(de, drdx(1), drdx(2), drdx(3), theta, d(1), d(2), d(3));

		Matrix3D de2dposx(de(1),de(4),de(7),de(2),de(5),de(8),de(3),de(6),de(9));
		Matrix3D de3dposx(de(13),de(16),de(19),de(14),de(17),de(20),de(15),de(18),de(21));
		Vector3D de2dtheta(de(10),de(11),de(12));
		Vector3D de3dtheta(de(22),de(23),de(24));

		Vector sf(4);
		GetShapesPos(sf, x0);
		Vector sfx(4);
		GetShapesPosx(sfx, x0);
		Vector sfrot(2);
		GetShapesRot(sfrot, xi);

		for (int i=1; i<=3; i++)
		{
			int ip3 = i+3;
			int ip6 = i+6;
			for (int j=1; j<=14; j++)
			{
				L0(i,j)=0.;
			}
			for (int j=1; j<=4; j++)
			{
				L0(i,i+(j-1)*3)=sf(j);
			}
			for (int j=1; j<=3; j++)
			{
				for (int k=1; k<=4; k++)
				{
					L0(ip3,j+(k-1)*3)=de2dposx(i,j)*sfx(k);
					L0(ip6,j+(k-1)*3)=de3dposx(i,j)*sfx(k);
				}
			}
			for (int k=1; k<=2; k++)
			{
				L0(ip3,k+12)=de2dtheta(i)*sfrot(k);
				L0(ip6,k+12)=de3dtheta(i)*sfrot(k);
			}
		}
	}
}

void ANCFBeam3DTorsion::GetdL0dqk(const double& xi, const int k, Matrix& dL0dqk) const //eq.(53), xi in [-lx/2, lx/2], k in {1,..,14}
{
	dL0dqk.SetAll(0);

	ConstMatrix<3*14> ddeidqdqk(3, 14, 1); //nofill
	
	Getddeidqdqk(xi, 2, k, ddeidqdqk);
	dL0dqk.SetSubmatrix(ddeidqdqk, 3+1, 1); 
	
	Getddeidqdqk(xi, 3, k, ddeidqdqk);
	dL0dqk.SetSubmatrix(ddeidqdqk, 2*3+1, 1);
}

void ANCFBeam3DTorsion::Getddeidqdqk(const double& xi, const int i, const int k, Matrix& ddeidqdqk) const //eq.(81), xi in [-lx/2, lx/2], i in {2,3}, k in {1,..,14}
{
	ddeidqdqk.SetAll(0);

	// dd ei / dq dqk = ( [ dd ei / dd (r',ang) ] . [ d (r',ang) / dqk ] ) . [ d (r',ang) / dq ] + [ d ei / d (r',ang) ] . [ dd (r',ang) / dqdqk ]
	//                = ( [ dd ei / dd (r',ang) ] . [ d (r',ang) / dqk ] ) . [ d (r',ang) / dq ] + 0
	
	ConstMatrix<3*3*3> ddeidposxdposx(3,3*3,1); //nofill
	ConstMatrix<3*3> ddeidposxdrot(3,3,1); //nofill
	ConstVector<3> ddeidrotdrot(3,1); //nofill

	Getddeidposxdposx(xi, i, ddeidposxdposx);
	Getddeidposxdrot(xi, i, ddeidposxdrot);
	Getddeidrotdrot(xi, i, ddeidrotdrot);

	ConstMatrix<3*3> ddeidposxdposx_times_dposxdqk_plus_ddeidrotdposx_times_drotdqk(3,3,1); //nofill
	ConstVector<3> ddeidrotdrot_times_drotdqk_plus_ddeidposxdrot_times_dposxdqk(3,1); //nofill

	double x0 = xi*2./size.X();

	if (k >= 1 && k <= 12)
	{
		int posx_dim = (k-1)%Dim() + 1;
		int sf_idx = (k-1)/Dim()+1;
		for (int d=1; d<=Dim(); d++)
		{
			ddeidrotdrot_times_drotdqk_plus_ddeidposxdrot_times_dposxdqk(d) = ddeidposxdrot(d,posx_dim)*GetSFPosx(sf_idx,x0);
			for (int e=1; e<=Dim(); e++)
			{
				ddeidposxdposx_times_dposxdqk_plus_ddeidrotdposx_times_drotdqk(d,e) = ddeidposxdposx(d,(posx_dim-1)*Dim()+e)*GetSFPosx(sf_idx,x0);
			}
		}
	}
	else if (k==13 || k==14)
	{
		int sf_idx = k-NSPos()*Dim();
		for (int d=1; d<=Dim(); d++)
		{
			ddeidrotdrot_times_drotdqk_plus_ddeidposxdrot_times_dposxdqk(d) = ddeidrotdrot(d)*GetSFRot(sf_idx,xi);
			for (int e=1; e<=Dim(); e++)
			{
				ddeidposxdposx_times_dposxdqk_plus_ddeidrotdposx_times_drotdqk(d,e) = ddeidposxdrot(d,e)*GetSFRot(sf_idx,xi);
			}
		}
	}
	else
	{
		assert(0);
	}

	for (int d=1; d<=Dim(); d++)
	{
		for (int j=1; j<=NSPos(); j++)
		{
			for (int e=1; e<=Dim(); e++)
			{
				ddeidqdqk(d,e+Dim()*(j-1)) = ddeidposxdposx_times_dposxdqk_plus_ddeidrotdposx_times_drotdqk(d,e)*GetSFPosx(j,x0);
			}
		}
		for (int j=1; j<=NSRot(); j++)
		{
			ddeidqdqk(d,Dim()*NSPos()+j) = ddeidrotdrot_times_drotdqk_plus_ddeidposxdrot_times_dposxdqk(d)*GetSFRot(j,xi);
		}
	}
}

void ANCFBeam3DTorsion::Getddeidposxdposx(const double& xi, const int i, Matrix& ddeidposxdposx) const //eq.(79)+(91), xi in [-lx/2, lx/2], i in {2,3}
{
	ConstMatrix<3*3*3> ddei0dposxdposx(3,3*3,1); //nofill
	double theta = GetRot(xi); 

	if (i==2)
	{
		Getddei0dposxdposx(xi, 2, ddei0dposxdposx);
		ddeidposxdposx = Cos(theta)*ddei0dposxdposx;
		Getddei0dposxdposx(xi, 3, ddei0dposxdposx);
		ddeidposxdposx += Sin(theta)*ddei0dposxdposx;
	}
	else if (i == 3)
	{
		Getddei0dposxdposx(xi, 2, ddei0dposxdposx);
		ddeidposxdposx = -Sin(theta)*ddei0dposxdposx;
		Getddei0dposxdposx(xi, 3, ddei0dposxdposx);
		ddeidposxdposx += Cos(theta)*ddei0dposxdposx;
	}
	else
	{
		assert(0);
	}
}

void ANCFBeam3DTorsion::Getddei0hatdposxdposx(const double& xi, const int i, Matrix& ddei0hatdposxdposx) const //eq.(80)-(84), xi in [-lx/2, lx/2], i in {2,3}
{
	if ( i<2 || i>3 )
	{
		assert(0);
	}

	if (i==2)
	{
		ddei0hatdposxdposx.SetAll(0.);
	}
	else if (i==3)
	{ 
		// dde30hatdposxdposx  = ((ddei0hatde1de1 de1dposx) de1dposx) + dei0hatde1 dde1dposxdposx      where e1 = posx/|posx|

		// v.. vector, depending on vector x, and vector x depending on vector y
		// chain rule:
		// dv(j)dy(c1) = dv(j)dx(s1) dx(s1)dy(c1)      (j: row index, c1: first column index, s1: first summation index)
		// ddv(j)dy(c1)dy(c2) = ddv(j)dx(s1)dx(s2) dx(s2)dy(c2) dx(s1)dy(c1) + dv(j)dx(s1) ddx(s1)dy(c1)dy(c2)      (c2: second column index, s2: second summation index)
	
		IVector offset(Dim());    //helper
		for(int j=1; j<=Dim(); j++)
		{
			offset(j) = (j-1)*Dim();
		}

		ConstMatrix<27> dde30hatde1de1(3,9,1); //nofill
		Getdde30hatde1de1(xi, dde30hatde1de1);
		ConstMatrix<9> de1dposx(3,3,1); //nofill
		Getde1dposx(xi, de1dposx);
		ConstMatrix<9> de30hatde1(3,3,1); //nofill
		Getde30hatde1(xi, de30hatde1);
		ConstMatrix<27> dde1dposxdposx(3,9,1); //nofill
		Getdde1dposxdposx(xi, dde1dposxdposx);

		ddei0hatdposxdposx.SetAll(0.);
		for (int j=1; j<=Dim(); j++)  //row index
		{
			for (int c1=1; c1<=Dim(); c1++)  // first column index
			{
				for (int c2=1; c2<=c1; c2++)  // second column index
				{
					for (int s1=1; s1<=Dim(); s1++)  // first summation index
					{
						ddei0hatdposxdposx(j, offset(c1) + c2) += de30hatde1(j, s1)*dde1dposxdposx(s1, offset(c1) + c2);

						for (int s2=1; s2<=Dim(); s2++)  // second summation index
						{
							ddei0hatdposxdposx(j, offset(c1) + c2) += dde30hatde1de1(j, s1 + offset(s2))*de1dposx(s2, c2)*de1dposx(s1, c1);
						}
					}
					
					if (c2 < c1)
					{
						ddei0hatdposxdposx(j, offset(c2) + c1) = ddei0hatdposxdposx(j, offset(c1) + c2); //symmetric components w.r.t. differentiation variables
					}
				}
			}
		}
	}
}



void ANCFBeam3DTorsion::Getddei0dposxdposx(const double& xi, const int i, Matrix& ddei0dposxdposx) const //eq.(57)-(58)+(79)+(91), xi in [-lx/2, lx/2], i in {2,3}
{
	ConstVector<3> ei0hat(3,1); //nofill
	Getei0hat(xi, i, ei0hat);	
	ConstMatrix<3*3> dei0hatdposx(3,3,1); //nofill
	Getdei0hatdposx(xi, i, dei0hatdposx);
	ConstMatrix<3*3*3> ddei0hatdposxdposx(3,3*3,1); //nofill
	Getddei0hatdposxdposx(xi, i, ddei0hatdposxdposx);

	ConstVector<3> dei0hatdposxj(3,1); //nofill
	ConstVector<3> dei0hatdposxk(3,1); //nofill
	ConstVector<3> ddei0hatdposxjdposxk(3,1); //nofill
	ConstVector<3> ddei0dposxjdposxk(3,1); //nofill

	for (int j=1; j<=3; j++)
	{
		int offset = (j-1)*Dim();
		dei0hatdposx.GetColVec(j, dei0hatdposxj);
		for (int k=1; k<=j; k++)
		{
			ddei0hatdposxdposx.GetColVec(offset+k, ddei0hatdposxjdposxk);

			if (j==k)
			{
				GetddFOverAbsFdxdx(ei0hat, dei0hatdposxj, ddei0hatdposxjdposxk, ddei0dposxjdposxk);
			} 
			else
			{
				dei0hatdposx.GetColVec(k, dei0hatdposxk);
				GetddFOverAbsFdxdy(ei0hat, dei0hatdposxj, dei0hatdposxk, ddei0hatdposxjdposxk, ddei0dposxjdposxk);
			}

			ddei0dposxdposx.SetColVec(ddei0dposxjdposxk, offset+k);
			if (k<j)
			{
				ddei0dposxdposx.SetColVec(ddei0dposxjdposxk, (k-1)*Dim()+j);  //symmetric components
			}
		}
	}
}

void ANCFBeam3DTorsion::Getddeidposxdrot(const double& xi, const int i, Matrix& ddeidposxdrot) const //eq.(87)+(92), xi in [-lx/2, lx/2], i in {2,3}
{
	ConstMatrix<3*3> dei0dposx(3,3,1); //nofill
	double theta = GetRot(xi); 

	if (i==2)
	{
		Getdei0dposx(xi, 2, dei0dposx);
		ddeidposxdrot = -Sin(theta)*dei0dposx;
		Getdei0dposx(xi, 3, dei0dposx);
		ddeidposxdrot += Cos(theta)*dei0dposx;
	}
	else if (i == 3)
	{
		Getdei0dposx(xi, 2, dei0dposx);
		ddeidposxdrot = -Cos(theta)*dei0dposx;
		Getdei0dposx(xi, 3, dei0dposx);
		ddeidposxdrot += -Sin(theta)*dei0dposx;
	}
	else
	{
		assert(0);
	}
}

void ANCFBeam3DTorsion::Getddeidrotdrot(const double& xi, const int i, Vector& ddeidrotdrot) const //eq.(88)+(93), xi in [-lx/2, lx/2], i in {2,3}
{
	ConstVector<3> ei0(3,1); //nofill
	double theta = GetRot(xi); 

	if (i==2)
	{
		Getei0(xi, 2, ei0);
		ddeidrotdrot = -Cos(theta)*ei0;
		Getei0(xi, 3, ei0);
		ddeidrotdrot += -Sin(theta)*ei0;
	}
	else if (i == 3)
	{
		Getei0(xi, 2, ei0);
		ddeidrotdrot = Sin(theta)*ei0;
		Getei0(xi, 3, ei0);
		ddeidrotdrot += -Cos(theta)*ei0;
	}
	else
	{
		assert(0);
	}
}

void ANCFBeam3DTorsion::Getdeidq(const double& xi, const int i, Matrix& deidq) const //eq.(66), xi in [-lx/2, lx/2]
{
	deidq.SetAll(0);

	// delta deidq = [ d ei / d (r',ang) ] . [ d (r',ang) / dq ]
	ConstMatrix<3*3> deidposx(3,3,1); //nofill
	ConstVector<3> deidrot(3,1); //nofill

	Getdeidposx(xi, i, deidposx);
	Getdeidrot(xi, i, deidrot);

	double x0 = xi*2./size.X();

	for (int d=1; d<=Dim(); d++)
	{
		for (int k=1; k<=Dim(); k++)
		{
			for (int j=1; j<=NSPos(); j++)
			{
				deidq(d,k+Dim()*(j-1)) = deidposx(d,k)*GetSFPosx(j,x0);
			}
		}
		for (int k=1; k<=NSRot(); k++)
		{
			deidq(d,Dim()*NSPos()+k) = deidrot(d)*GetSFRot(k,xi);
		}
	}
}

void ANCFBeam3DTorsion::Getdeidposx(const double& xi, const int i, Matrix& deidposx) const //eq.(67)+(79), xi in [-lx/2, lx/2], i in {2,3}
{
	ConstMatrix<3*3> dei0dposx(3,3,1); //nofill
	double theta = GetRot(xi);

	if (i==2)
	{
		Getdei0dposx(xi, 2, dei0dposx);
		deidposx = Cos(theta)*dei0dposx;
		Getdei0dposx(xi, 3, dei0dposx);
		deidposx += Sin(theta)*dei0dposx;
	}
	else if (i==3)
	{
		Getdei0dposx(xi, 2, dei0dposx);
		deidposx = -Sin(theta)*dei0dposx;
		Getdei0dposx(xi, 3, dei0dposx);
		deidposx += Cos(theta)*dei0dposx;
	}
	else
	{
		assert(0);
	}
}

void ANCFBeam3DTorsion::Getdei0dposx(const double& xi, const int i, Matrix& dei0dposx) const //eq.(69)-(76) + (55)-(56), xi in [-lx/2, lx/2], i in {2,3}
{
	ConstVector<3> ei0hat(3,1); //nofill
	Getei0hat(xi, i, ei0hat);	
	ConstMatrix<3*3> dei0hatdposx(3,3,1); //nofill
	Getdei0hatdposx(xi, i, dei0hatdposx);

	ConstVector<3> dei0hatdposxj(3,1);//nofill
	ConstVector<3> dei0dposxj(3,1);//nofill

	for (int j=1; j<=3; j++)
	{
		
		dei0hatdposx.GetColVec(j, dei0hatdposxj);
		GetdFOverAbsFdx(ei0hat, dei0hatdposxj, dei0dposxj);
		dei0dposx.SetColVec(dei0dposxj, j);
	}
}

void ANCFBeam3DTorsion::Getdei0hatdposx(const double& xi, const int i, Matrix& dei0hatdposx) const //eq.(69)-(76), xi in [-lx/2, lx/2], i in {2,3}
{	
	if (i==2)
	{
		Vector3D d = GetDirector(xi);
		// dei0hatdposx = skew(d)
		dei0hatdposx(1,1) = 0;
		dei0hatdposx(2,2) = 0;
		dei0hatdposx(3,3) = 0;
		dei0hatdposx(1,2) = -d(3);
		dei0hatdposx(2,1) = d(3);
		dei0hatdposx(1,3) = d(2);
		dei0hatdposx(3,1) = -d(2);
		dei0hatdposx(2,3) = -d(1);
		dei0hatdposx(3,2) = d(1);
	}
	else if (i==3)
	{
		// dei0hatdposx = de30hatde1.de1dposx
		ConstMatrix<3*3> de30hatde1(3,3,1); //nofill
		Getde30hatde1(xi, de30hatde1);
		ConstMatrix<3*3> de1dposx(3,3,1); //nofill
		Getde1dposx(xi, de1dposx);

		dei0hatdposx = de30hatde1*de1dposx;
	}
	else
	{
		assert(0);
	}
}

void ANCFBeam3DTorsion::Getde1dposx(const double& xi, Matrix& de1dposx) const //xi in [-lx/2, lx/2]
{
	//de1dposx = (1/|posx|) I - (1/|posx|³) posx posx^T           where e1=posx/|posx|
	Vector3D posx = GetPosx(xi);
	double posx_abs_square = posx*posx;
	double posx_abs = sqrt(posx_abs_square);
	double posx_abs_cubic = posx_abs_square*posx_abs;
	Vector3D d = GetDirector(xi);

	double diag = 1./posx_abs;
	for (int j=1; j<=3; j++)
	{
		for (int k=1; k<=3; k++)
		{
			de1dposx(j,k) = -posx(j)*posx(k)/posx_abs_cubic;
		}
		de1dposx(j,j) += diag;
	}
}

void ANCFBeam3DTorsion::Getdde1dposxdposx(const double& xi, Matrix& dde1dposxdposx) const //xi in [-lx/2, lx/2], k in {1,2,3}, dde1dposxdposx has 3 rows (according to dde1) and 3*3 columns (according to dposxdposx)
{
	//dde1dposxdposx = (dde1dposxdposxk)_k\in{1,2,3}, dde1dposxdposxk = posxk/|posx|^3 (e1 e1^T - I) - 1/|posx|(de1dposxk e1^T + e1 de1dposxk^T),  where e1 = posx/|posx| and I is identity
	Vector3D posx = GetPosx(xi);
	double posx_abs_square = posx*posx;
	double posx_abs = sqrt(posx_abs_square);
	double posx_abs_cubic = posx_abs_square*posx_abs;

	Vector3D e1((1./posx_abs)*posx);
	ConstMatrix<3*3> de1dposx(3,3);
	Getde1dposx(xi, de1dposx);

	for (int k=1; k<=3; k++)
	{
		double fact = posx(k)/posx_abs_cubic;
		for (int i=1; i<=3; i++)
		{
			for (int j=1; j<=k; j++)
			{
				dde1dposxdposx(i,j+(k-1)*3) = fact*(e1(i)*e1(j) - (double)(i==j))	- (1./posx_abs)*(de1dposx(i,k)*e1(j) + e1(i)*de1dposx(j,k));
				if (j < k)
				{
					dde1dposxdposx(i,k+(j-1)*3) = dde1dposxdposx(i,j+(k-1)*3); // symmetric components
				}
			}
		}
	}
}

void ANCFBeam3DTorsion::Getde30hatde1(const double& xi, Matrix& de30hatde1) const //xi in [-lx/2, lx/2]
{
	// de30hatde1 = - e1 d^T - (e1^T d) I        where e1=posx/|posx|
	Vector3D posx = GetPosx(xi);
	double posx_abs_square = posx*posx;
	double posx_abs = sqrt(posx_abs_square);
	Vector3D d = GetDirector(xi);

	Vector3D e1 = (1/posx_abs)*posx;
	double diag = -(e1*d);
	for (int j=1; j<=3; j++)
	{
		for (int k=1; k<=3; k++)
		{
			de30hatde1(j,k) = -e1(j)*d(k);
		}
		de30hatde1(j,j) += diag;
	}
}

void ANCFBeam3DTorsion::Getdde30hatde1de1(const double& xi, Matrix& dde30hatde1de1) const //xi in [-lx/2, lx/2], dde30hatde1de1 has 3 rows (according to de30hat) and 3*3 columns (according to de1de1)
{
	// dde30hatde1de1 = (dde30hatde1de1k)_k\in{1,2,3}, dde30hatde1de1k = - Ek d^T - d_k I, where E1 = (1,0,0)^T, E2 = (0,1,0)^T, E3 = (0,0,1)^T,   e1=posx/|posx|

	dde30hatde1de1.SetAll(0.);
	Vector3D d = GetDirector(xi);

	for(int k=1; k<=3; k++)
	{
		int offset = Dim()*(k-1);
		for (int j=1; j<=3; j++)
		{
			dde30hatde1de1(k,j+offset) -= d(j);
			dde30hatde1de1(j,j+offset) -= d(k);
		}
	}
}

void ANCFBeam3DTorsion::Getdeidrot(const double& xi, const int i, Vector& deidrot) const //eq.(68)+(80), xi in [-lx/2, lx/2], i in {2,3}
{
	ConstVector<3> ei0(3,1); //nofill
	double theta = GetRot(xi);

	if (i==2)
	{
		Getei0(xi, 2, ei0);
		deidrot = -Sin(theta)*ei0;
		Getei0(xi, 3, ei0);
		deidrot += Cos(theta)*ei0;
	}
	else if (i==3)
	{
		Getei0(xi, 2, ei0);
		deidrot = -Cos(theta)*ei0;
		Getei0(xi, 3, ei0);
		deidrot += -Sin(theta)*ei0;
	}
	else
	{
		assert(0);
	}
}

void ANCFBeam3DTorsion::Getei0(const double& xi, const int i, Vector& ei0) const //eq.(59)-(62) + (54), xi in [-lx/2, lx/2], i in {2,3}
{	
	ConstVector<3> ei0hat(3,1); //nofill
	Getei0hat(xi, i, ei0hat);	
	GetFOverAbsF(ei0hat, ei0);
}

void ANCFBeam3DTorsion::Getei0hat(const double& xi, const int i, Vector& ei0hat) const //eq.(59)-(62), xi in [-lx/2, lx/2], i in {2,3}
{
	Vector3D d = GetDirector(xi);
	Vector3D posx = GetPosx(xi);
	Vector3D ei0hat_3d;
	
	if (i==2)
	{
		ei0hat_3d = d.Cross(posx);
	}
	else if (i==3)
	{
		// Gramm Schmitt Projektion
		ei0hat_3d = d-(d*posx)/(posx*posx)*posx;
	}
	else
	{
		assert(0);
	}

	ei0hat = ei0hat_3d;
}

void ANCFBeam3DTorsion::GetFOverAbsF(const Vector& f, Vector& f_over_abs_f) const //eq.(54)  computes (f/|f|)', where f maps a scalar to a vector   //move to linalg.h/cpp
{
	f_over_abs_f = (1/sqrt(f*f))*f;
}

void ANCFBeam3DTorsion::GetdFOverAbsFdx(const Vector& f, const Vector& dfdx, Vector& d_f_over_abs_f_prime_dx) const //eq.(55)-(56)  computes (f/|f|)', where f maps a scalar to a vector   //move to linalg.h/cpp
{
	double f_abs_square = f*f;
	double f_abs = sqrt(f_abs_square);
	d_f_over_abs_f_prime_dx = (-(f*dfdx)/(f_abs_square*f_abs))*f + (1/f_abs)*dfdx;
}

void ANCFBeam3DTorsion::GetddFOverAbsFdxdx(const Vector& f, const Vector& dfdx, const Vector& ddfdxdx, Vector& dd_f_over_abs_f_dxdx) const //eq.(57)-(58)  computes (f/|f|)'', where f maps a scalar to a vector   //move to linalg.h/cpp
{
	double f_f = f*f;
	double f_abs_to_the_fifth = pow(f_f,5./2.);
	double f_dfdx = f*dfdx;
	double dfdx_dfdx = dfdx*dfdx;
	double f_ddfdxdx = f*ddfdxdx;
	
	dd_f_over_abs_f_dxdx =
		((3*f_dfdx*f_dfdx - (dfdx_dfdx + f_ddfdxdx)*f_f)/f_abs_to_the_fifth) * f
		+ (-2.*f_f*f_dfdx/f_abs_to_the_fifth) * dfdx
		+ (f_f*f_f/f_abs_to_the_fifth) * ddfdxdx;
}

void ANCFBeam3DTorsion::GetddFOverAbsFdxdy(const Vector& f, const Vector& dfdx, const Vector& dfdy, const Vector& ddfdxdy, Vector& dd_f_over_abs_f_dxdy) const //eq.(57)-(58)  computes (f/|f|)'', where f maps a scalar to a vector   //move to linalg.h/cpp
{
	double f_f = f*f;
	double f_abs_to_the_fifth = pow(f_f,5./2.);
	double f_dfdx = f*dfdx;
	double f_dfdy = f*dfdy;
	double dfdx_dfdy = dfdx*dfdy;
	double f_ddfdxdy = f*ddfdxdy;
	
	dd_f_over_abs_f_dxdy =
		((3*f_dfdx*f_dfdy - (dfdx_dfdy + f_ddfdxdy)*f_f)/f_abs_to_the_fifth) * f
		+ (-f_f*f_dfdy/f_abs_to_the_fifth) * dfdx
		+ (-f_f*f_dfdx/f_abs_to_the_fifth) * dfdy
		+ (f_f*f_f/f_abs_to_the_fifth) * ddfdxdy;
}







//kopiert von ANCFBeam3D und ComputeStress-Aufrufe auskommentiert
void ANCFBeam3DTorsion::DrawElement()
{
	mbs->SetColor(col);

  //Vector3D p1 = GetPos3D0D(Vector3D(-1,0.,0.));
  //Vector3D p2 = GetPos3D0D(Vector3D(1,0.,0.));

	//UO() << "p1=" << p1 << "\n";
	//UO() << "p2=" << p2 << "\n";


	double lx1 = 1; double ly1 = 1*GetMBS()->GetMagnifyYZ(); double lz1 = 1*GetMBS()->GetMagnifyYZ();
	double def_scale = GetMBS()->GetDOption(105); //deformation scaling

	int linemode = 1; //0=no lines, 1=outline+color, 2=outline, 3=elementline+color, 4=elementline
	if (GetMBS()->GetIOption(110) && !GetMBS()->GetIOption(111))
	{
		linemode = 2;
	}
	if (!GetMBS()->GetIOption(110) && GetMBS()->GetIOption(111))
	{
		linemode = 0;
	}
	if (!GetMBS()->GetIOption(110) && !GetMBS()->GetIOption(111))
	{
		linemode = 4;
	}

	int colormode = 0;
	if (GetMBS()->GetActualPostProcessingFieldVariable() != NULL) colormode = 1;
	//if (GetMBS()->GetDrawMode() != 0) colormode = 1;//Yuri
	else
	{
		colormode = 0; 
	}

	//GetMBS()->UO(0) << "colormode = " << colormode << ", linemode = " << linemode << "\n";

	//int comp1, comp2, type;
	//GetMBS()->GetDrawTypeComponents(type, comp1, comp2);

	if (colormode)
	{
		double modeval = 0;
		int xgset = 0;

		double tilex = GetMBS()->GetIOption(137);
		double tiley = GetMBS()->GetIOption(138);
		TArray<Vector3D> points((int)(tilex+1)*(int)(tiley+1));
		TArray<double> vals((int)(tilex+1)*(int)(tiley+1));
		double v=0;

		for (int side=1; side <= 6; side++)
		{
			points.SetLen(0); vals.SetLen(0);
			Vector3D p0, vx, vy;
			int tileyn = (int)tiley;
			int tilexn = (int)tilex;

			if (side == 1)
			{ //bottom
				p0 = Vector3D(-lx1,-ly1,-lz1);
				vx = Vector3D(2.*lx1/tilexn,0,0);
				vy = Vector3D(0,0,2.*lz1/tileyn);
			}
			else if (side == 2)
			{ //top
				p0 = Vector3D(-lx1, ly1, lz1);
				vx = Vector3D(2.*lx1/tilexn,0,0);
				vy = Vector3D(0,0,-2.*lz1/tileyn);
			}
			else if (side == 3)
			{ //front
				p0 = Vector3D(-lx1,-ly1, lz1);
				vx = Vector3D(2.*lx1/tilexn,0,0);
				vy = Vector3D(0,2.*ly1/tileyn,0);
			}
			else if (side == 4)
			{ //back
				p0 = Vector3D(-lx1, ly1,-lz1);
				vx = Vector3D(2.*lx1/tilexn,0,0);
				vy = Vector3D(0,-2.*ly1/tileyn,0);
			}
			else if (side == 5)
			{ //left
				p0 = Vector3D(-lx1, ly1,-lz1);
				vx = Vector3D(0,-2.*ly1/tilexn,0);
				vy = Vector3D(0,0,2.*lz1/tileyn);
			}
			else if (side == 6)
			{ //right
				p0 = Vector3D( lx1,-ly1,-lz1);
				vx = Vector3D(0,2.*ly1/tilexn,0);
				vy = Vector3D(0,0,2.*lz1/tileyn);
			}

			for (double iy = 0; iy <= tileyn+1e-10; iy++)
			{
				for (double ix = 0; ix <= tilexn+1e-10; ix++)
				{
					Vector3D ploc = (p0+ix*vx+iy*vy);//-1<=p0<=+1
					Vector3D ploc_scaled(ploc(1)*size.X()*0.5,ploc(2)*size.Y()*0.5,ploc(3)*size.Z()*0.5);//-l/2<=ploc_scaled<=+l/2
					Vector3D pg = GetPos3D0D(ploc, def_scale);
					//Vector3D p(pg);//p=pg
					points.Add(pg);
						if (colormode)
							v = GetFieldVariableValue(*GetMBS()->GetActualPostProcessingFieldVariable(), ploc_scaled, true);
						vals.Add(v);

					//if (colormode) ComputeStressD(ploc,type,comp1,comp2,v);//old_KN
					//vals.Add(v);
				}
			}
			mbs->DrawColorQuads(points,vals,(int)tilexn+1,(int)tileyn+1,colormode,linemode);
			//mbs->UO()<<"points"<<points(1)<<"\n";
		}
	}
	else
	{
		int drawgrid = 0;
		double thickness = 1;
		//double tiling = 20;
		double tiling = GetMBS()->GetIOption(136);
		mbs->SetDrawlines(0);
		mbs->SetDrawlinesH(0);
		Vector3D p1,p2,p3,p4,p5,p6,p7,p8;
		for (double i = 0; i < tiling; i++)
		{
			double l1 = -lx1+2*lx1*i/tiling;
			double l2 = -lx1+2*lx1*(i+1)/tiling;
			if (i == 0)
			{
				p8 = Vector3D(GetPos3D0D(Vector3D(l1,-ly1,-lz1)));
				p7 = Vector3D(GetPos3D0D(Vector3D(l1,-ly1, lz1)));
				p4 = Vector3D(GetPos3D0D(Vector3D(l1, ly1,-lz1))); 
				p3 = Vector3D(GetPos3D0D(Vector3D(l1, ly1, lz1)));
				if (drawgrid)
				{
					mbs->MyDrawLine(p8,p7,thickness);
					mbs->MyDrawLine(p7,p3,thickness);
					mbs->MyDrawLine(p4,p3,thickness);
					mbs->MyDrawLine(p4,p8,thickness);
				}
			}
			else
			{
				p8 = p6;
				p7 = p5;
				p4 = p2;
				p3 = p1;
			}
			p6 = Vector3D(GetPos3D0D(Vector3D(l2,-ly1,-lz1)));
			p5 = Vector3D(GetPos3D0D(Vector3D(l2,-ly1, lz1)));
			p2 = Vector3D(GetPos3D0D(Vector3D(l2, ly1,-lz1)));
			p1 = Vector3D(GetPos3D0D(Vector3D(l2, ly1, lz1)));
			int dout = 0;
			if (i == 0 || i == tiling-1) dout = 1;
			if (drawgrid)
			{
				mbs->MyDrawLine(p6,p5,thickness);
				mbs->MyDrawLine(p5,p1,thickness);
				mbs->MyDrawLine(p2,p1,thickness);
				mbs->MyDrawLine(p2,p6,thickness);
				mbs->MyDrawLine(p6,p8,thickness);
				mbs->MyDrawLine(p5,p7,thickness);
				mbs->MyDrawLine(p2,p4,thickness);
				mbs->MyDrawLine(p1,p3,thickness);
			}
			else
			{
				mbs->DrawHex(p1,p2,p3,p4,p5,p6,p7,p8,dout);
			}
		}

		mbs->SetDrawlines(1);

		// draw directors
		Vector3D position, director;
		double drawlength = max(size.Y(),size.Z())*2.;
		double linethickness = max(2.,0.04*mbs->GetDOption(120));
		Vector3D color(.0, .5, .0);
		
		position = GetPos3D0D(Vector3D(-1, 0., 0.));   //left director
		director = GetDirector1(1);
		director.Normalize();
		director *= drawlength;
		mbs->MyDrawLine(position, position+director, linethickness, color);

		position = GetPos3D0D(Vector3D(1, 0., 0.));    //right director
		director = GetDirector2(1);
		director.Normalize();
		director *= drawlength;
		mbs->MyDrawLine(position, position+director, linethickness, color);
	}
};


void ANCFBeam3DTorsion::GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables)
{
	// add all the field variables of the parent class
	Body3D::GetAvailableFieldVariables(variables);
	//FVT_position,FVT_displacement, FVT_velocity: implemented in Element
	// possibility to remove entries does not exist yet, just do not use the Function of the parent class

	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_axial_extension);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_force_axial);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_curvature, FieldVariableDescriptor::FVCI_z);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_moment, FieldVariableDescriptor::FVCI_z);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_torsion);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_moment_torsional);
}

double ANCFBeam3DTorsion::GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & ploc, bool flagD)
{
	// ploc is local position in [-lx/2., +lx/2.] x [-ly/2., +ly/2.] x [-lz/2., +lz/2.]

	double xi = ploc.X();   

	double kappa1; GetKappa1(xi,kappa1,flagD);
	double kappa2; GetKappa2(xi,kappa2,flagD);
	double kappa3; GetKappa3(xi,kappa3,flagD);
	double gamma1; GetGamma1(xi,gamma1,flagD);

	double beamGJ = GetMaterial().BeamGJkx();
	double beamEIy = GetMaterial().BeamEIy();
	double beamEIz = GetMaterial().BeamEIz();
	double beamEA = GetMaterial().BeamEA();

	switch(fvd.VariableType())
	{
	case FieldVariableDescriptor::FVT_displacement: return fvd.GetComponent(GetDisplacement(ploc, flagD));
	case FieldVariableDescriptor::FVT_beam_axial_extension: return gamma1;
	case FieldVariableDescriptor::FVT_beam_curvature: return fvd.GetComponent(Vector3D(kappa1, kappa2, kappa3));
	case FieldVariableDescriptor::FVT_beam_torsion: return kappa1;
	case FieldVariableDescriptor::FVT_beam_force_axial: return beamEA*gamma1;
	case FieldVariableDescriptor::FVT_beam_moment: return fvd.GetComponent(Vector3D(beamGJ*kappa1, beamEIy*kappa2, beamEIz*kappa3));
	case FieldVariableDescriptor::FVT_beam_moment_torsional: return beamGJ*kappa1;
	default: return Body3D::GetFieldVariableValue(fvd, ploc, flagD);
	}

	return FIELD_VARIABLE_NO_VALUE;
}


