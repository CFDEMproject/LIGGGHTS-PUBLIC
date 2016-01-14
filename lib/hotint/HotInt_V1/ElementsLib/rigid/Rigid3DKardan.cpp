//#**************************************************************
//#
//# filename:             Rigid3DKardan.cpp
//#
//# author:               Ludwig Rafael, Gerstmayr Johannes, Yury Vetyukov
//#
//# generated:					  22.04.2010
//# description:          3D Element Library - rigid body with kardan angles
//# remarks:						  
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
 
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++         Rigid3DKARDAN      ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "element.h"
#include "body3d.h"
#include "femathhelperfunctions.h"
#include "material.h"
#include "node.h"
#include "rigid3d.h"
#include "rigid3dkardan.h"
#include "myfile.h"
#include "elementdataaccess.h"

#define DEFAULT_ROTATIONS_SEQUENCE xyz

Rigid3DKardan::RotationsSequence Rigid3DKardan::defaultRotationsSequence = DEFAULT_ROTATIONS_SEQUENCE;

Rigid3DKardan::Rigid3DKardan(MBS* mbsi):Rigid3D(mbsi), rotationsSequence(defaultRotationsSequence)
{ 
	ElementDefaultConstructorInitialization(); //$ DR 2012-12: ElementDefaultConstructorInitialization added
}

Rigid3DKardan::Rigid3DKardan(MBS* mbsi, const Vector& x0, const Vector& phi0, double rhoi, double Vi, const Vector3D& Ip,
														 const Vector3D& si, const Vector3D& coli):
Rigid3D(mbsi, x0, phi0, rhoi, Vi, Ip, si, coli),
rotationsSequence(defaultRotationsSequence)
{
	Vector3D xp(x0(1),x0(2),x0(3));
	Vector3D vp(x0(4),x0(5),x0(6));
	Vector3D phi(phi0(1),phi0(2),phi0(3));
	Vector3D phip(phi0(4),phi0(5),phi0(6));
	ComputeInitialConditions(xp, vp, phi, phip, x_init);

};

Rigid3DKardan::Rigid3DKardan(MBS* mbsi, const Vector& x0, const Vector& phi0, double rhoi, double Vi, const Matrix3D& Ip,
														 const Vector3D& si, const Vector3D& coli):
Rigid3D(mbsi, x0, phi0, rhoi, Vi, Ip, si, coli),
rotationsSequence(defaultRotationsSequence)
{	
	Vector3D xp(x0(1),x0(2),x0(3));
	Vector3D vp(x0(4),x0(5),x0(6));
	Vector3D phi(phi0(1),phi0(2),phi0(3));
	Vector3D phip(phi0(4),phi0(5),phi0(6));
	ComputeInitialConditions(xp, vp, phi, phip, x_init);
	elementname = GetElementSpec();
};

Rigid3DKardan::Rigid3DKardan(MBS* mbsi, const Vector& x0, const Vector& phi0, double rhoi, double Vi, const Matrix3D& Ip,
														 const Vector3D& si, const Vector3D& coli, RotationsSequence rotSeq):
Rigid3D(mbsi, x0, phi0, rhoi, Vi, Ip, si, coli),
rotationsSequence(rotSeq)
{	
	Vector3D xp(x0(1),x0(2),x0(3));
	Vector3D vp(x0(4),x0(5),x0(6));
	Vector3D phi(phi0(1),phi0(2),phi0(3));
	Vector3D phip(phi0(4),phi0(5),phi0(6));
	ComputeInitialConditions(xp, vp, phi, phip, x_init);
	elementname = GetElementSpec();
};

Rigid3DKardan::Rigid3DKardan(MBS* mbsi, const Vector& x0, const Vector& phi0, double rhoi,
														 const Vector3D& si, const Vector3D& coli):
Rigid3D(mbsi, x0,phi0, rhoi, si, coli),
rotationsSequence(defaultRotationsSequence)
{
	Vector3D xp(x0(1),x0(2),x0(3));
	Vector3D vp(x0(4),x0(5),x0(6));
	Vector3D phi(phi0(1),phi0(2),phi0(3));
	Vector3D phip(phi0(4),phi0(5),phi0(6));
	ComputeInitialConditions(xp, vp, phi, phip, x_init);
	elementname = GetElementSpec();

	//////// HACK AP
	//double xinit_double[12] = {0.9771035473772125, 0.8119080368041125, 0.8149451665816472, 0.4345721941229590, 0.4196667118043127, 0.4654793137759562, 0.6468043502831272, 0.1137455299827921, 0.8447456427769386, 0.7828729837346822, -0.05939125082671293, -0.4017941539650536};
	//for (int i=1; i<=12; i++)
	//	x_init(i) = xinit_double[i-1];

};

//$ DR 2012-12: ElementDefaultConstructorInitialization added
void Rigid3DKardan::ElementDefaultConstructorInitialization()
{
	ComputeInitialConditions(Vector3D(0.),Vector3D(0.),Vector3D(0.), Vector3D(0.), x_init);
	elementname = GetElementSpec();
}

void Rigid3DKardan::ComputeInitialConditions(const Vector3D& xp, const Vector3D& vp, 
																						 const Vector3D& phi, const Vector3D& phip, Vector& xinit)
{
	xinit.SetLen(SS()); 
	xinit(1) = xp.X();
	xinit(2) = xp.Y();
	xinit(3) = xp.Z();
	xinit(7) = vp.X();
	xinit(8) = vp.Y();
	xinit(9) = vp.Z();

	Matrix3D A = ComputeRotMatrixEuler(phi.X(),phi.Y(),phi.Z());

	Vector3D theta;
	switch(rotationsSequence)
	{
	case xyz:
		theta(2) = -asin(A(3,1));
		theta(3) = atan2(A(2,1)/*/cos(theta(2))*/,A(1,1)/*/cos(theta(2))*/);
		theta(1) = atan2(A(3,2)/*/cos(theta(2))*/,A(3,3)/*/cos(theta(2))*/);
		break;
	case zxy:
		theta(2) = -asin(A(2,3));
		theta(3) = atan2(A(1,3)/*/cos(theta(2))*/,A(3,3)/*/cos(theta(2))*/);
		theta(1) = atan2(A(2,1)/*/cos(theta(2))*/,A(2,2)/*/cos(theta(2))*/);
		break;
	case zxz:
		theta(2) = acos(A(3,3));
		theta(3) = atan2(A(1,3)/*/sin(theta(2))*/,-A(2,3)/*/sin(theta(2))*/);
		theta(1) = atan2(A(3,1)/*/sin(theta(2))*/,A(3,2)/*/sin(theta(2))*/);
		break;
	default: assert(0);
	}
	xinit(4) = theta.X();
	xinit(5) = theta.Y();
	xinit(6) = theta.Z();	

	for(int i=NRotParam()+7; i<=2*NRotParam()+6;i++)
	{
		xinit(i) = 0.;
	}

	Matrix3D G3 = GetG_Kardan(xinit(4), xinit(5), xinit(6));
	Matrix G(3,3);
	for(int i=1;i<=3;i++)
	{
		for(int j=1;j<=3;j++)
		{
			G(i,j) = G3(i,j); // Solve is only defined for Matrix
		}
	}

	Vector dthetadt(3);
	int rv = G.Solve(Vector(phip(1), phip(2), phip(3)), dthetadt);
	//UO() << "ComputeInitialConditions: theta = " << theta << ", dthetadt = " << dthetadt << "\n";
	if (!rv) {GetMBS()->UO() << "Rigid3DKardan:Initialization: could not determine initial Kardan parameter velocities due to singularity!!!\n";}
	else
	{									 
		xinit(10) = dthetadt(1);
		xinit(11) = dthetadt(2);
		xinit(12) = dthetadt(3);		
	}
}


void Rigid3DKardan::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Body3D::GetElementData(edc);
	ElementData ed;

	//$ DR 2013-01-30 removed this option, see change log 379
	//ed.SetInt(1,"inertia_mode",1,3); ed.SetToolTipText("1.. directly entered mass (m) and inertia values (I) are used. 2..m and I are computed from density and body dimensions. 3..m and I are computed from density and GeomElement.");  edc.TreeAdd("Physics",ed);

	Vector3D phi;
	// Rotation matrix --> Quaternions --> Euler Angles
	Matrix3D A;
	if(ltg.Length())	{	A = GetRotMatrix();	}		
	else {	A = Matrix3D(1.);}	//$ DR 2012-12-13: necessary if GetElementData is called before assemble, ltg is not initialized at this moment

	double b0, b1, b2, b3;
	RotMatToQuaternions(A, b0,b1,b2,b3);
	QuaternionsToKardanAngles(b0,b1,b2,b3, phi);

	ed.SetVector3D(phi.X(), phi.Y(), phi.Z(), "initial_rotation"); 
	ed.SetToolTipText(mystr("3 consecutive rotations (global rotation axes): [rot3_X, rot2_Y, rot1_Z] in rad")); edc.TreeAdd("Initialization",ed);

	Vector3D phip;
	Matrix3D G = GetG(Vector(x_init(4), x_init(5), x_init(6)));
	phip = G*Vector3D(x_init(10), x_init(11), x_init(12)); //omega = G beta_p

	ed.SetVector3D(phip.X(), phip.Y(), phip.Z(), "initial_angular_velocity"); 
	ed.SetToolTipText(mystr("Angular velocity vector in global coordinates: [ang_X, ang_Y, ang_Z] in rad/s")); edc.TreeAdd("Initialization",ed);
}

int Rigid3DKardan::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Body3D::SetElementData(edc);

	Vector3D xp, vp, phi, phip;

	xp=Vector3D(x_init(1),x_init(2),x_init(3)); //$ DR 2012-07: is now in Get/SetElementDataAuto]
	vp=Vector3D(x_init(NRotParam()+4),x_init(NRotParam()+5),x_init(NRotParam()+6)); //$ DR 2012-07: is now in Get/SetElementDataAuto]

	GetElemDataVector3D(GetMBS(), edc, "Initialization.initial_rotation", phi, 0);

	GetElemDataVector3D(GetMBS(), edc, "Initialization.initial_angular_velocity", phip, 0);

	Matrix3D rot = ComputeRotMatrixWithKardanAngles(phi.X(),phi.Y(),phi.Z());
	RotMatToEulerAngles(rot, phi);

	ComputeInitialConditions(xp, vp, phi, phip, x_init);

	//$ DR 2013-01-30 removed this option, see change log 379
	//SetElementDataMassAndInertia(edc);

	return rv;
}


void Rigid3DKardan::CopyFrom(const Element& e)
{
	Rigid3D::CopyFrom(e);
	const Rigid3DKardan& ce = (const Rigid3DKardan&)e;
	rotationsSequence = ce.rotationsSequence;
}

Element* Rigid3DKardan::GetCopy()
{
	Element* ec = new Rigid3DKardan(*this);
	return ec;
}

void Rigid3DKardan::EvalF2(Vector& f, double t) 
{
	Body3D::EvalF2(f,t);
	//compare Shabana 1998, p. 158, eq. 3.156/157 and p.159, eq. 3.164
	//quadratic velocity vector (Shabana Computational Dynamics 1994, p. 396 or p. 414 in newer version)
	//right hand side of equation 
	//	Term1:			  GbarP^T*Iphi*(GbarP*beta-Gbar*betaP)
	//  Term2:      - Gbar^T*Iphi*GbarP*betaP   ... 3.148 (zero in case of Euler parameter)                 
	//---------------------------------------------------------------------
	//	Term1:			  GbarP^T*Iphi*(GbarP*beta-Gbar*betaP)
	//---------------------------------------------------------------------
	//center of mass must be at origin(0,0,0) !!!!!!!!!!!!!!!!!!!!

#ifdef OLD_VERSION		// old version

	//theta
	ConstVector<4> beta_cv(NRotParam());
	Vector3D beta;
	beta_cv.LinkWith(beta.GetVecPtr(), NRotParam()); 
	GetBeta(beta_cv);

	// thetaP
	ConstVector<4> betap_cv(NRotParam());
	Vector3D betap;
	betap_cv.LinkWith(betap.GetVecPtr(), NRotParam()); 
	GetBetaP(betap_cv);
 
	// Term1
	Matrix3D GbarpTp = GetGbarpT();
	Matrix3D Gbarp = GetGbarp();
	Matrix3D Gbar = GetGbar();
	//alternative:
	Vector3D omegabar = -1.0*Gbar*betap; //                        -Gbar*betaP
	omegabar += Gbarp * beta;          // tp: JG            (GbarP*beta-Gbar*betaP)
	Vector3D temp = Iphi*omegabar;       //        Iphi*(GbarP*beta-Gbar*betaP)
	Mult(GbarpTp,temp,betap_cv);         //GbarP^T*Iphi*(GbarP*beta-Gbar*betaP)
	
	for(int i = 1;i<=NRotParam();i++)
	{
		//add GbarP^T*Iphi*(GbarP*beta-Gbar*betaP) to right hand side of dynamic equation
		f(i+3) -= betap(i); //betap is now result of the above evaluation! ///JG - added
	}
	////---------- begin test ------------------
	//Vector3D theta(XG(4),XG(5),XG(6));
	//Vector3D thetaP(XG(4+6),XG(5+6),XG(6+6));
	//Vector3D t1 = GetGbarpT()*Iphi*(GetGbarp()*theta-GetGbar()*thetaP);
	//Vector3D t2 = (-1.)*GetGbarT()*Iphi*GetGbarp()*thetaP;
	//UO() << "!!! Term1 = " <<	 betap << " = " << t1 << "\n";
	////---------- end test ------------------

	//---------------------------------------------------------------------
	//  Term2:      - Gbar^T*Iphi*Gbarp*thetaP   ... 3.148 (zero in case of Euler parameter)
	//---------------------------------------------------------------------
	GetBetaP(betap_cv);

	//result=vector^T*matrix=matrix^T*vector, if vector is first element of multiplication 
	temp = Gbarp * betap;				//tp: JG //tmp is now Gbarp*betap;	 old: Matrix3D Gbarp = GetGbarp();		
	betap	= Iphi * temp;          //betap is now Iphi*GbarP*omegabarP     old: betap = Iphi*tmp;          

	temp = Gbar.GetTp() * betap;  //tp:JG  		  //tmp is now Gbar.GetTp()*Iphi*GbarP*omegabarP	 old: Matrix GbarT = GetGbarT();Mult(GbarT, betap, tmp);

	f(4) -= -temp(1); //temp is now result of the above evaluation! //JG
	f(5) -= -temp(2);
	f(6) -= -temp(3);
	//UO() << "!!! Term2 = " <<	 -1.*temp << " = " << t2 << "\n";

#else		// new version - YV

	//center of mass must be at origin(0,0,0) !!!!!!!!!!!!!!!!!!!!
	// WHY?

	// we form the right-hand side of the equation
	/*
	// why doesn't it work this way?
	Vector beta;
	GetBeta(beta);
	Vector betaDot;
	GetBetaP(betaDot);
	*/
	//theta
	ConstVector<4> beta_cv(NRotParam());
	Vector3D beta;
	beta_cv.LinkWith(beta.GetVecPtr(), NRotParam()); 
	GetBeta(beta_cv);

	// thetaP
	ConstVector<4> betap_cv(NRotParam());
	Vector3D betaDot;
	betap_cv.LinkWith(betaDot.GetVecPtr(), NRotParam()); 
	GetBetaP(betap_cv);

	// first what follows from the term d/dt(dT/d(betaDot))
	// the matrix m_{\theta\theta} is on the left-hand side,
	// the right hand side follows as
	Matrix3D m = GetGbarT() * Iphi * GetGbarp(); 
	Matrix3D m1 = m;     //(PG) why not m.MakeSym()?
	m.Transpose();
	m1 += m;
	Vector3D rhs1 = m1 * betaDot;

	// and now what follows from the term dT/dq
	Vector3D v = Iphi * GetGbar() * betaDot;	// part of the product before dGbar/dbeta_i
	Vector3D rhs2(
		- ((GetDGbarDBeta1() * betaDot) * v),
		- ((GetDGbarDBeta2() * betaDot) * v),
		0
		);
	for(int i = 1; i <= 3; i++)
	{
		double linearTerm = rhs1(i) + rhs2(i);
		f(i + 3) -= linearTerm;
	}	

#endif	// old/new version
};

Matrix3D Rigid3DKardan::GetG_Kardan(double beta0, double beta1, double beta2)	const
{
	double beta = beta1;
	double gamma = beta2;
	double cbeta = cos(beta);
	double sbeta = sin(beta);
	double cgamma = cos(gamma);
	double sgamma = sin(gamma);

	switch(rotationsSequence)
	{
	case xyz: return Matrix3D(cbeta*cgamma,-sgamma,0,cbeta*sgamma,cgamma,0,-sbeta,0,1);
	case zxy: return Matrix3D(cbeta*sgamma,cgamma,0,-sbeta,0,1,cbeta*cgamma,-sgamma,0);
	case zxz: return Matrix3D(sbeta*sgamma,cgamma,0,-(cgamma*sbeta),sgamma,0,cbeta,0,1);
	}
	assert(0);
	return Matrix3D(0,0,0,0,0,0,0,0,0);
}

Matrix3D Rigid3DKardan::GetGT_Kardan(double beta0, double beta1, double beta2) const
{
	//Matrix3D G = GetG_Kardan(beta0,beta1,beta2);
	//G.Transpose();
	//return G;
	return (GetG_Kardan(beta0,beta1,beta2)).GetTp();
}

Matrix3D Rigid3DKardan::GetG() const
{
	ConstVector<4> beta(NRotParam());
	GetBeta(beta);

	return GetG_Kardan(beta(1),beta(2),beta(3));	
}

Matrix3D Rigid3DKardan::GetG(Vector& beta) const
{
	return GetG_Kardan(beta(1),beta(2),beta(3));	
}

Matrix3D Rigid3DKardan::GetGT() const
{
	ConstVector<4> beta(NRotParam());
	GetBeta(beta);

	return GetGT_Kardan(beta(1),beta(2),beta(3));	
}

Matrix3D Rigid3DKardan::GetGbar() const
{
	ConstVector<4> b(NRotParam());
	GetBeta(b);

	double alpha = b(1);
	double beta = b(2);
	double calpha = cos(alpha);
	double salpha = sin(alpha);
	double cbeta = cos(beta);
	double sbeta = sin(beta);

	switch(rotationsSequence)
	{
	case xyz: return Matrix3D(1,0,-sbeta,0,calpha,cbeta*salpha,0,-salpha,calpha*cbeta);
	case zxy: return Matrix3D(0,calpha,cbeta*salpha,0,-salpha,calpha*cbeta,1,0,-sbeta);
	case zxz: return Matrix3D(0,calpha,salpha*sbeta,0,-salpha,calpha*sbeta,1,0,cbeta);
	}
	assert(0);
	return Matrix3D(0,0,0,0,0,0,0,0,0);
}

Matrix3D Rigid3DKardan::GetGbarT() const
{
	//Matrix3D G = GetGbar();
	//G.Transpose();
	//return G;
	return GetGbar().GetTp();
}

Matrix3D Rigid3DKardan::GetDGbarDBeta1() const
{
	ConstVector<4> b(NRotParam());
	GetBeta(b);

	double alpha = b(1);
	double beta = b(2);
	double calpha = cos(alpha);
	double salpha = sin(alpha);
	double cbeta = cos(beta);
	double sbeta = sin(beta);

	switch(rotationsSequence)
	{
	case xyz: return Matrix3D(0,0,0,0,-salpha,calpha*cbeta,0,-calpha,-(cbeta*salpha));
	case zxy: return Matrix3D(0,-salpha,calpha*cbeta,0,-calpha,-(cbeta*salpha),0,0,0);
	case zxz: return Matrix3D(0,-salpha,calpha*sbeta,0,-calpha,-(salpha*sbeta),0,0,0);
	}
	assert(0);
	return Matrix3D(0,0,0,0,0,0,0,0,0);
}

Matrix3D Rigid3DKardan::GetDGbarDBeta2() const
{
	ConstVector<4> b(NRotParam());
	GetBeta(b);

	double alpha = b(1);
	double beta = b(2);
	double calpha = cos(alpha);
	double salpha = sin(alpha);
	double cbeta = cos(beta);
	double sbeta = sin(beta);

	switch(rotationsSequence)
	{
	case xyz: return Matrix3D(0,0,-cbeta,0,0,-(salpha*sbeta),0,0,-(calpha*sbeta));
	case zxy: return Matrix3D(0,0,-(salpha*sbeta),0,0,-(calpha*sbeta),0,0,-cbeta);
	case zxz: return Matrix3D(0,0,cbeta*salpha,0,0,calpha*cbeta,0,0,-sbeta);
	}
	assert(0);
	return Matrix3D(0,0,0,0,0,0,0,0,0);
}

Matrix3D Rigid3DKardan::GetGbarp() const
{
	ConstVector<4> betap(NRotParam());
	GetBetaP(betap);
	return betap(1) * GetDGbarDBeta1() + betap(2) * GetDGbarDBeta2();		// GetDGbarDBeta3 === 0
}

Matrix3D Rigid3DKardan::GetGbarpT() const
{
	//Matrix3D m = GetGbarp();
	//m.Transpose();
	return GetGbarp().GetTp();
}

Matrix3D Rigid3DKardan::GetDOmegaDTheta(Vector& betap) const
{
	Matrix3D DGbarDbeta1 = GetDGbarDBeta1();
	Matrix3D DGbarDbeta2 = GetDGbarDBeta2();
	Vector3D DGbarDbeta1_betap = DGbarDbeta1*betap;
	Vector3D DGbarDbeta2_betap = DGbarDbeta2*betap;
	// d GbarT / d theta thetaP = d omega / d theta
	return Matrix3D (DGbarDbeta1_betap(1), DGbarDbeta2_betap(1), 0,
		                     DGbarDbeta1_betap(2), DGbarDbeta2_betap(2), 0,
												 DGbarDbeta1_betap(3), DGbarDbeta2_betap(3), 0);
}

Vector3D Rigid3DKardan::GetPosD(const Vector3D& p_loc) const
{
	if (GetMBS()->GetIOption(151) && IsRigid() && GetMBS()->GetDOption(105) != 1.) //for deformation scaling
	{
		double fact = GetMBS()->GetDOption(105);

		ConstVector<3> beta(3);
		ConstVector<3> betainit(3);
		GetBetaD(beta);
		GetBetaInitD(betainit);

		beta *= fact;
		betainit *= 1-fact;
	
		beta += betainit;

		Matrix3D A = ComputeRotMatrix(beta);

		return A*p_loc+GetRefPosD();
	}
	else
	{
		return GetRotMatrixD()*p_loc+GetRefPosD();
	}
};


Matrix3D Rigid3DKardan::ComputeRotMatrix(const Vector& b) const
{
	double alpha = b(1);
	double beta = b(2);
	double gamma = b(3);
	double calpha = cos(alpha);
	double salpha = sin(alpha);
	double cbeta = cos(beta);
	double sbeta = sin(beta);
	double cgamma = cos(gamma);
	double sgamma = sin(gamma);

	switch(rotationsSequence)
	{
	case xyz: return Matrix3D(cbeta*cgamma,cgamma*salpha*sbeta - calpha*sgamma,
		 calpha*cgamma*sbeta + salpha*sgamma,cbeta*sgamma,
		 calpha*cgamma + salpha*sbeta*sgamma,
		 -(cgamma*salpha) + calpha*sbeta*sgamma,-sbeta,cbeta*salpha,
		 calpha*cbeta);
	case zxy: return Matrix3D(calpha*cgamma + salpha*sbeta*sgamma,
		 -(cgamma*salpha) + calpha*sbeta*sgamma,cbeta*sgamma,
		 cbeta*salpha,calpha*cbeta,-sbeta,
		 cgamma*salpha*sbeta - calpha*sgamma,
		 calpha*cgamma*sbeta + salpha*sgamma,cbeta*cgamma);
	case zxz: return Matrix3D(calpha*cgamma - cbeta*salpha*sgamma,
		 -(cgamma*salpha) - calpha*cbeta*sgamma,sbeta*sgamma,
		 cbeta*cgamma*salpha + calpha*sgamma,
		 calpha*cbeta*cgamma - salpha*sgamma,-(cgamma*sbeta),
		 salpha*sbeta,calpha*sbeta,cbeta);
	}
	assert(0);
	return Matrix3D(0,0,0,0,0,0,0,0,0);
}

Matrix3D Rigid3DKardan::ComputeRotMatrixP(const Vector& b, const Vector& dbdt) const
{
	double alpha = b(1);
	double beta = b(2);
	double gamma = b(3);

	double dalphadt = dbdt(1);
	double dbetadt = dbdt(2);
	double dgammadt = dbdt(3);

	double calpha = cos(alpha);
	double salpha = sin(alpha);
	double cbeta = cos(beta);
	double sbeta = sin(beta);
	double cgamma = cos(gamma);
	double sgamma = sin(gamma);

	switch(rotationsSequence)
	{
	case xyz: return Matrix3D(-(cgamma*dbetadt*sbeta) - cbeta*dgammadt*sgamma,
		 cbeta*cgamma*dbetadt*salpha + 
			dalphadt*(calpha*cgamma*sbeta + salpha*sgamma) - 
			dgammadt*(calpha*cgamma + salpha*sbeta*sgamma),
		 calpha*cbeta*cgamma*dbetadt + 
			dalphadt*(-(cgamma*salpha*sbeta) + calpha*sgamma) + 
			dgammadt*(cgamma*salpha - calpha*sbeta*sgamma),
		 cbeta*cgamma*dgammadt - dbetadt*sbeta*sgamma,
		 cbeta*dbetadt*salpha*sgamma + 
			dgammadt*(cgamma*salpha*sbeta - calpha*sgamma) + 
			dalphadt*(-(cgamma*salpha) + calpha*sbeta*sgamma),
		 calpha*cbeta*dbetadt*sgamma + 
			dgammadt*(calpha*cgamma*sbeta + salpha*sgamma) - 
			dalphadt*(calpha*cgamma + salpha*sbeta*sgamma),
		 -(cbeta*dbetadt),calpha*cbeta*dalphadt - dbetadt*salpha*sbeta,
		 -(cbeta*dalphadt*salpha) - calpha*dbetadt*sbeta);
	case zxy: return Matrix3D(cbeta*dbetadt*salpha*sgamma + 
			dgammadt*(cgamma*salpha*sbeta - calpha*sgamma) + 
			dalphadt*(-(cgamma*salpha) + calpha*sbeta*sgamma),
		 calpha*cbeta*dbetadt*sgamma + 
			dgammadt*(calpha*cgamma*sbeta + salpha*sgamma) - 
			dalphadt*(calpha*cgamma + salpha*sbeta*sgamma),
		 cbeta*cgamma*dgammadt - dbetadt*sbeta*sgamma,
		 calpha*cbeta*dalphadt - dbetadt*salpha*sbeta,
		 -(cbeta*dalphadt*salpha) - calpha*dbetadt*sbeta,
		 -(cbeta*dbetadt),cbeta*cgamma*dbetadt*salpha + 
			dalphadt*(calpha*cgamma*sbeta + salpha*sgamma) - 
			dgammadt*(calpha*cgamma + salpha*sbeta*sgamma),
		 calpha*cbeta*cgamma*dbetadt + 
			dalphadt*(-(cgamma*salpha*sbeta) + calpha*sgamma) + 
			dgammadt*(cgamma*salpha - calpha*sbeta*sgamma),
		 -(cgamma*dbetadt*sbeta) - cbeta*dgammadt*sgamma);
		case zxz: return Matrix3D(dbetadt*salpha*sbeta*sgamma - 
			dgammadt*(cbeta*cgamma*salpha + calpha*sgamma) - 
			dalphadt*(cgamma*salpha + calpha*cbeta*sgamma),
		 calpha*dbetadt*sbeta*sgamma + 
			dgammadt*(-(calpha*cbeta*cgamma) + salpha*sgamma) + 
			dalphadt*(-(calpha*cgamma) + cbeta*salpha*sgamma),
		 cgamma*dgammadt*sbeta + cbeta*dbetadt*sgamma,
		 -(cgamma*dbetadt*salpha*sbeta) + 
			dalphadt*(calpha*cbeta*cgamma - salpha*sgamma) + 
			dgammadt*(calpha*cgamma - cbeta*salpha*sgamma),
		 -(calpha*cgamma*dbetadt*sbeta) - 
			dalphadt*(cbeta*cgamma*salpha + calpha*sgamma) - 
			dgammadt*(cgamma*salpha + calpha*cbeta*sgamma),
		 -(cbeta*cgamma*dbetadt) + dgammadt*sbeta*sgamma,
		 cbeta*dbetadt*salpha + calpha*dalphadt*sbeta,
		 calpha*cbeta*dbetadt - dalphadt*salpha*sbeta,-(dbetadt*sbeta));
	}
	assert(0);
	return Matrix3D(0,0,0,0,0,0,0,0,0);
}

// creates a matlab/simulink .c file for this element
void Rigid3DKardan::MatlabExportF2()//(Vector& f, double t)
{
// subroutine still incomplete, missing (M)inv ...
// mex file already works
// nametags
	mystr name_iphi =            "Iphiphi";
  mystr name_mass =            "Mass"; 
	mystr name_beta =            "Theta";
	mystr name_betap =           "ThetaP";

	mystr name_gbar =            "Gbar";
	mystr name_gbart =           "GbarT";
	mystr name_gbarp =           "GbarP";
	mystr name_dgbardbeta1 =     "dGbar_dTheta_1";
	mystr name_dgbardbeta2 =     "dGbar_dTheta_2";
	mystr name_m =               "M";
	mystr name_rhs1 =            "Rhs1";
	mystr name_rhs2 =            "Rhs2";
	mystr name_linearterm =      "LinearTerm";
	mystr name_lineartermfull =  "LT6";
	mystr name_massmatrix =      "MassMatrix";
	mystr name_q =               "Q";

	mystr name_dummymatrix =     "temp_mat";
	mystr name_dummyvector =     "temp_vec";
	mystr name_dummyvector2 =    "temp_vec2";
	mystr name_dummydouble =     "temp_double";

	mystr rotseq;
	switch(rotationsSequence)
	{
	case xyz: rotseq = "xyz"; break;
	case zxy: rotseq = "zxy"; break;
	case zxz: rotseq = "zxz"; break;
	}
	mystr buffer("");
	ConstVector<3> vec;

	int standalonec = 1;
  int singlestepmex = 1;
	int sfunction = 1;
	int evals = 1;
	int ansimath = 1;
/*
if(standalonec)
{
	mystr filename(mystr("Rigid3DKardan")+mystr("_test.c"));    
	CMFile file_export(filename,TFMwrite);
	ConstVector<4> beta(3);
	ConstVector<4> betap(3);
	GetBeta(beta);
	GetBetaP(betap);
//header
	file_export.RWSoleStr(mystr("//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"));
	file_export.RWSoleStr(mystr("//+ single ansi c file for better debugging"));
	file_export.RWSoleStr(mystr("//+ automatically created by Rigid3DKardan::MatlabExportF2"));
	file_export.RWSoleStr(mystr(""));
  file_export.RWSoleStr(mystr("#include \"Rigid3DKardan")+mystr(".h\""));
	file_export.RWSoleStr(mystr(""));
	file_export.RWSoleStr(mystr("void main( )"));
	file_export.RWSoleStr(mystr("{"));
	file_export.RWSoleStr(mystr("  double ")+name_lineartermfull+mystr("[6] = {0.0};"));
	file_export.RWSoleStr(mystr("  double* ")+name_linearterm+mystr(" = &")+name_lineartermfull+mystr("[3];"));
	file_export.RWSoleStr(mystr("  double ")+name_massmatrix+mystr("[6][6];"));
  file_export.RWSoleStr(mystr("  double* ")+name_m+mystr(" = (double*) ")+name_massmatrix+mystr(";"));
	file_export.RWSoleStr(mystr("  double ")+name_q+mystr("[6]};"));

	file_export.RWSoleStr(beta.MakeString(mystr("  double ")+name_beta));
	file_export.RWSoleStr(betap.MakeString(mystr("  double ")+name_betap));
	file_export.RWSoleStr(mystr("  EvalF2( ")+name_linearterm+mystr(", ")+name_beta+mystr(", ")+name_betap+mystr(");")); 
	file_export.RWSoleStr(mystr("  EvalM( ")+name_massmatrix+mystr(", ")+name_beta+mystr(", ")+name_betap+mystr(");")); 
	file_export.RWSoleStr(mystr("  MatrixInvert( ")+name_massmatrix+mystr(", ")+name_massmatrix+mystr(", 6, 6);"));
	file_export.RWSoleStr(mystr("}"));
}
*/
if(singlestepmex)
{
	mystr filename(mystr("Rigid3DKardan")/*+mystr(this->GetElnum())*/+mystr(".c"));
	CMFile file_export(filename,TFMwrite);
//header
	file_export.RWSoleStr(mystr("//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"));
	file_export.RWSoleStr(mystr("//+ single step MEX-wrapping: use [LinTerm,Minv] = fctn (Theta,Theta_dot)"));
	file_export.RWSoleStr(mystr("//+ automatically created by Rigid3DKardan::MatlabExportF2"));
	file_export.RWSoleStr(mystr(""));
	file_export.RWSoleStr(mystr("#include \"mex.h\""));
  file_export.RWSoleStr(mystr("#include \"Rigid3DKardan")/*+mystr(this->GetElnum())*/+mystr(".h\""));
	file_export.RWSoleStr(mystr(""));
	file_export.RWSoleStr(mystr("void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )"));
	file_export.RWSoleStr(mystr("{"));
	file_export.RWSoleStr(mystr("// c-function := void F2(double* )")+name_linearterm+mystr(", double* ")+name_beta+mystr(", double* ")+name_betap+mystr("), result is double* ")); 
	file_export.RWSoleStr(mystr("  double *")+name_beta+mystr(", *")+name_betap+mystr(", *")+name_linearterm+mystr(", *")+name_massmatrix+mystr(";"));
	file_export.RWSoleStr(mystr("  ")+name_beta+mystr(" = mxGetPr(prhs[0]);"));
	file_export.RWSoleStr(mystr("  ")+name_betap+mystr(" = mxGetPr(prhs[1]);"));
	file_export.RWSoleStr(mystr("  plhs[0] = mxCreateDoubleMatrix(1, 3, mxREAL);"));  
	file_export.RWSoleStr(mystr("  plhs[1] = mxCreateDoubleMatrix(6, 6, mxREAL);"));  
	file_export.RWSoleStr(mystr("  ")+name_linearterm+mystr(" = mxGetPr(plhs[0]);"));
	file_export.RWSoleStr(mystr("  EvalF2( ")+name_linearterm+mystr(", ")+name_beta+mystr(", ")+name_betap+mystr(");")); 
	file_export.RWSoleStr(mystr("  ")+name_massmatrix+mystr(" = mxGetPr(plhs[1]);"));
	file_export.RWSoleStr(mystr("  EvalM( ")+name_massmatrix+mystr(", ")+name_beta+mystr(", ")+name_betap+mystr(");")); 
	file_export.RWSoleStr(mystr("  MatrixInvert( ")+name_massmatrix+mystr(", ")+name_massmatrix+mystr(", 6, 6);"));
	file_export.RWSoleStr(mystr("}"));  
}
if(sfunction)
{
	mystr filename(mystr("sRigid3DKardan")/*+mystr(this->GetElnum())*/+mystr(".c"));
	CMFile file_export(filename,TFMwrite);
//header
	file_export.RWSoleStr(mystr("//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"));
	file_export.RWSoleStr(mystr("//+ Simulink S-Function template:  "));
	file_export.RWSoleStr(mystr("//+ specify initial conditions as ic = [,,,,,,,,,,,] 12D vector"));
	file_export.RWSoleStr(mystr("//+ Iphi and Mr taken from .h-file"));
	file_export.RWSoleStr(mystr("//+ automatically created by Rigid3DKardan::MatlabExportF2"));
	file_export.RWSoleStr(mystr(""));
  file_export.RWSoleStr(mystr("#define S_FUNCTION_NAME sRigid3DKardan")/*+mystr(this->GetElnum())*/);
  file_export.RWSoleStr(mystr("#define S_FUNCTION_LEVEL 2"));
  file_export.RWSoleStr(mystr("#include \"simstruc.h\""));
  file_export.RWSoleStr(mystr("#include \"Rigid3DKardan")/*+mystr(this->GetElnum())*/+mystr(".h\""));
  file_export.RWSoleStr(mystr(""));
  file_export.RWSoleStr(mystr("/* defines for simulation sizes */"));
  file_export.RWSoleStr(mystr("#define NUMSFCNPARAMS    1"));
  file_export.RWSoleStr(mystr("#define NUMCONTSTATES    12"));
  file_export.RWSoleStr(mystr("#define NUMDISCSTATES    0"));
  file_export.RWSoleStr(mystr("#define NUMINPUTPORTS    1"));
  file_export.RWSoleStr(mystr("#define INPUTPORT0WIDTH  6"));
  file_export.RWSoleStr(mystr("#define FEEDTHROUGHPORT0 1"));
  file_export.RWSoleStr(mystr("#define NUMOUTPUTPORTS   1"));
  file_export.RWSoleStr(mystr("#define OUTPUTPORT0WIDTH 12"));
  file_export.RWSoleStr(mystr(""));
  file_export.RWSoleStr(mystr("static void mdlInitializeSizes(SimStruct *S)"));
  file_export.RWSoleStr(mystr("{"));
  file_export.RWSoleStr(mystr("  ssSetNumSFcnParams(S,NUMSFCNPARAMS); /* Number of expected parameters */"));
  file_export.RWSoleStr(mystr("  if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) return;"));
  file_export.RWSoleStr(mystr("  ssSetNumContStates(S,NUMCONTSTATES);"));
  file_export.RWSoleStr(mystr("  ssSetNumDiscStates(S,NUMDISCSTATES);"));
  file_export.RWSoleStr(mystr("  if (!ssSetNumInputPorts(S,NUMINPUTPORTS)) return;"));
  file_export.RWSoleStr(mystr("  ssSetInputPortWidth(S,0,INPUTPORT0WIDTH);"));
  file_export.RWSoleStr(mystr("  ssSetInputPortDirectFeedThrough(S,0,FEEDTHROUGHPORT0);"));
  file_export.RWSoleStr(mystr("  if (!ssSetNumOutputPorts(S,NUMOUTPUTPORTS)) return;"));
  file_export.RWSoleStr(mystr("  ssSetOutputPortWidth(S,0,OUTPUTPORT0WIDTH);"));
  file_export.RWSoleStr(mystr("  ssSetNumSampleTimes(S,1);"));
  file_export.RWSoleStr(mystr("  ssSetNumRWork(S,0);"));
  file_export.RWSoleStr(mystr("  ssSetNumIWork(S,0);"));
  file_export.RWSoleStr(mystr("  ssSetNumPWork(S,0);"));
  file_export.RWSoleStr(mystr("  ssSetNumModes(S,0);"));
  file_export.RWSoleStr(mystr("  ssSetNumNonsampledZCs(S,0);"));
  file_export.RWSoleStr(mystr("  ssSetOptions(S,0);"));
  file_export.RWSoleStr(mystr("}\n"));
	file_export.RWSoleStr(mystr("static void mdlInitializeSampleTimes(SimStruct *S)"));
  file_export.RWSoleStr(mystr("{"));
  file_export.RWSoleStr(mystr("  ssSetSampleTime(S,0,CONTINUOUS_SAMPLE_TIME);"));
  file_export.RWSoleStr(mystr("  ssSetOffsetTime(S,0,0.0);"));
  file_export.RWSoleStr(mystr("}\n"));
  file_export.RWSoleStr(mystr("#define MDL_INITIALIZE_CONDITIONS"));
  file_export.RWSoleStr(mystr("static void mdlInitializeConditions(SimStruct *S)"));
  file_export.RWSoleStr(mystr("{"));
  file_export.RWSoleStr(mystr("//\"external\" input parameters must be defined in \"S-Function Parameters\" and in Matbab itself !"));
  file_export.RWSoleStr(mystr("  real_T *x0 = ssGetContStates(S);"));
  file_export.RWSoleStr(mystr("  real_T *ic = mxGetPr(ssGetSFcnParam(S,0));"));
  file_export.RWSoleStr(mystr("  int_T i;"));
  file_export.RWSoleStr(mystr("  for(i=0; i<NUMCONTSTATES; i++) { x0[i]=ic[i]; }"));
  file_export.RWSoleStr(mystr("}\n"));
  file_export.RWSoleStr(mystr("static void calcOutputs(real_T *y,real_T *x, SimStruct *S)"));
  file_export.RWSoleStr(mystr("{"));
  file_export.RWSoleStr(mystr("  int_T i;"));
  file_export.RWSoleStr(mystr("  for(i=0; i<NUMCONTSTATES; i++) { y[i]=x[i]; }"));
  file_export.RWSoleStr(mystr("// add derived outputs"));
  file_export.RWSoleStr(mystr("}\n"));
  file_export.RWSoleStr(mystr("static void mdlOutputs(SimStruct *S,int_T tid)"));
  file_export.RWSoleStr(mystr("{"));
  file_export.RWSoleStr(mystr("  real_T *y = ssGetOutputPortRealSignal(S,0);"));
  file_export.RWSoleStr(mystr("  real_T *x = ssGetContStates(S);"));
  file_export.RWSoleStr(mystr("  calcOutputs(y,x,S);"));
  file_export.RWSoleStr(mystr("}\n"));
  file_export.RWSoleStr(mystr("#define MDL_DERIVATIVES"));
  file_export.RWSoleStr(mystr("static void calcDerivatives(real_T *dx, real_T *x, real_T *u, SimStruct *S)"));
  file_export.RWSoleStr(mystr("{"));
  file_export.RWSoleStr(mystr("  real_T Minv[36];"));
  file_export.RWSoleStr(mystr("  real_T LinVec6[6] = {0.0};"));
  file_export.RWSoleStr(mystr("  real_T Q[6];"));
  file_export.RWSoleStr(mystr("  real_T *LinVec3 = &LinVec6[3];"));
  file_export.RWSoleStr(mystr("  real_T *Theta = &x[3];"));
  file_export.RWSoleStr(mystr("  real_T *ThetaP = &x[9];"));
  file_export.RWSoleStr(mystr("  int_T i;"));
  file_export.RWSoleStr(mystr(""));
  file_export.RWSoleStr(mystr("  EvalF2(LinVec3,Theta,ThetaP);"));
  file_export.RWSoleStr(mystr("  EvalM(Minv,Theta,ThetaP);"));
  file_export.RWSoleStr(mystr("  MatrixInvert(Minv,Minv,6,6);"));
  file_export.RWSoleStr(mystr("  MatrixMultiplyVector(Q,Minv,6,6,LinVec6,6);"));
  file_export.RWSoleStr(mystr(""));
  file_export.RWSoleStr(mystr("  for(i=0; i< 6; i++) {dx[ i] = x[ i+6];}"));
  file_export.RWSoleStr(mystr("  for(i=6; i<12; i++) {dx[ i] = Q[ i-6] + u[ i-6];}"));
	file_export.RWSoleStr(mystr("}\n"));
	file_export.RWSoleStr(mystr("static void mdlDerivatives(SimStruct *S)"));
	file_export.RWSoleStr(mystr("{"));
	file_export.RWSoleStr(mystr("  real_T *dx = ssGetdX(S);"));
	file_export.RWSoleStr(mystr("  real_T *x  = ssGetContStates(S);"));
	file_export.RWSoleStr(mystr("  InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0);"));
	file_export.RWSoleStr(mystr("  real_T u[INPUTPORT0WIDTH];"));
	file_export.RWSoleStr(mystr("  int_T i;"));
	file_export.RWSoleStr(mystr("  for(i=0; i<INPUTPORT0WIDTH; i++) { u[i]=*uPtrs[i]; }"));
	file_export.RWSoleStr(mystr("  calcDerivatives(dx,x,u,S);"));
	file_export.RWSoleStr(mystr("}\n"));
	file_export.RWSoleStr(mystr("static void mdlTerminate(SimStruct *S)"));
	file_export.RWSoleStr(mystr("{"));
	file_export.RWSoleStr(mystr("}\n"));
	file_export.RWSoleStr(mystr("#ifdef MATLAB_MEX_FILE   /* Compile as a MEX-file? */"));
	file_export.RWSoleStr(mystr("  #include \"simulink.c\"  /* MEX-file interface mechanism */"));
	file_export.RWSoleStr(mystr("#else"));
	file_export.RWSoleStr(mystr("  #include \"cg_sfun.h\"   /* Code generation registration */"));
	file_export.RWSoleStr(mystr("#endif"));
}
if(evals)
{
	mystr filename(mystr("Rigid3DKardan")/*+mystr(this->GetElnum())*/+mystr(".h"));
	CMFile file_export(filename,TFMwrite);
//header
	file_export.RWSoleStr(mystr("//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"));
	file_export.RWSoleStr(mystr("//+ evaluation functions, element specific"));
	file_export.RWSoleStr(mystr("//+ automatically created by Rigid3DKardan::MatlabExportF2"));
	file_export.RWSoleStr(mystr(""));
	file_export.RWSoleStr(mystr("#include \"stdlib.h\""));
	file_export.RWSoleStr(mystr("#include \"math.h\""));
	file_export.RWSoleStr(mystr("#include \"ansimath.h\""));
	file_export.RWSoleStr(mystr(""));
// constants
  file_export.RWSoleStr(mystr("//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"));
  file_export.RWSoleStr(mystr("// constants - element specific"));
	file_export.RWSoleStr(Iphi.MakeString("  double "+name_iphi));
	file_export.RWSoleStr(mystr("  double ")+name_mass+mystr(" = ")+mystr(this->GetMass())+mystr(";"));
	file_export.RWSoleStr(mystr(""));
//calculation suroutines
  file_export.RWSoleStr(mystr("//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"));
	file_export.RWSoleStr(mystr("//+ general calculations - rotationsequence specific"));
//Gbar
	file_export.RWSoleStr(mystr("int Calc")+name_gbar+mystr("(void* i")+name_gbar+mystr(", double* ")+name_beta+mystr(") // ")+mystr(rotseq));
	file_export.RWSoleStr(mystr("{")); 
	file_export.RWSoleStr(mystr("  double* ")+name_gbar+mystr(" = (double*) i")+name_gbar+mystr(";")); 
	switch(rotationsSequence)
	{
	case xyz:
	file_export.RWSoleStr(mystr("  ")+name_gbar+mystr("[0] = 1;"));
	file_export.RWSoleStr(mystr("  ")+name_gbar+mystr("[1] = 0;"));         
	file_export.RWSoleStr(mystr("  ")+name_gbar+mystr("[2] =-sin(")+name_beta+mystr("[1]);"));	
	file_export.RWSoleStr(mystr("  ")+name_gbar+mystr("[3] = 0;"));
	file_export.RWSoleStr(mystr("  ")+name_gbar+mystr("[4] = cos(")+name_beta+mystr("[0]);"));	
	file_export.RWSoleStr(mystr("  ")+name_gbar+mystr("[5] = cos(")+name_beta+mystr("[1])*sin(")+name_beta+mystr("[0]);"));	
	file_export.RWSoleStr(mystr("  ")+name_gbar+mystr("[6] = 0;"));
	file_export.RWSoleStr(mystr("  ")+name_gbar+mystr("[7] = -sin(")+name_beta+mystr("[0]);"));	
	file_export.RWSoleStr(mystr("  ")+name_gbar+mystr("[8] = cos(")+name_beta+mystr("[0])*cos(")+name_beta+mystr("[1]);"));	
	break;
	case zxy:
	file_export.RWSoleStr(mystr("  ")+name_gbar+mystr("[0] = 0;"));
	file_export.RWSoleStr(mystr("  ")+name_gbar+mystr("[1] = cos(")+name_beta+mystr("[0]);"));         
	file_export.RWSoleStr(mystr("  ")+name_gbar+mystr("[2] = cos(")+name_beta+mystr("[1])*sin(")+name_beta+mystr("[0]);"));     	
	file_export.RWSoleStr(mystr("  ")+name_gbar+mystr("[3] = 0;"));
	file_export.RWSoleStr(mystr("  ")+name_gbar+mystr("[4] =-sin(")+name_beta+mystr("[0]);")); 	
	file_export.RWSoleStr(mystr("  ")+name_gbar+mystr("[5] = cos(")+name_beta+mystr("[1])*cos(")+name_beta+mystr("[0]);"));	
	file_export.RWSoleStr(mystr("  ")+name_gbar+mystr("[6] = 1;"));
	file_export.RWSoleStr(mystr("  ")+name_gbar+mystr("[7] = 0;"));
	file_export.RWSoleStr(mystr("  ")+name_gbar+mystr("[8] =-sin(")+name_beta+mystr("[1]);"));	
	break;
	case zxz:
	file_export.RWSoleStr(mystr("  ")+name_gbar+mystr("[0] = 0;"));
	file_export.RWSoleStr(mystr("  ")+name_gbar+mystr("[1] = cos(")+name_beta+mystr("[0]);"));         
	file_export.RWSoleStr(mystr("  ")+name_gbar+mystr("[2] = sin(")+name_beta+mystr("[1])*sin(")+name_beta+mystr("[0]);"));     	
	file_export.RWSoleStr(mystr("  ")+name_gbar+mystr("[3] = 0;"));
	file_export.RWSoleStr(mystr("  ")+name_gbar+mystr("[4] =-sin(")+name_beta+mystr("[0]);")); 	
	file_export.RWSoleStr(mystr("  ")+name_gbar+mystr("[5] = sin(")+name_beta+mystr("[1])*cos(")+name_beta+mystr("[0]);"));	
	file_export.RWSoleStr(mystr("  ")+name_gbar+mystr("[6] = 1;"));
	file_export.RWSoleStr(mystr("  ")+name_gbar+mystr("[7] = 0;"));
	file_export.RWSoleStr(mystr("  ")+name_gbar+mystr("[8] = cos(")+name_beta+mystr("[1]);"));	
	break;
	}
	file_export.RWSoleStr(mystr("  return 0;"));
  file_export.RWSoleStr(mystr("}\n"));  
//DGbarDBeta1()
	file_export.RWSoleStr(mystr("int Calc")+name_dgbardbeta1+mystr("(void* i")+name_dgbardbeta1+mystr(", double* ")+name_beta+mystr(") // ")+mystr(rotseq));
	file_export.RWSoleStr(mystr("{"));  
	file_export.RWSoleStr(mystr("  double* ")+name_dgbardbeta1+mystr(" = (double*) i")+name_dgbardbeta1+mystr(";")); 
	switch(rotationsSequence)
	{
	case xyz:
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta1+mystr("[0] = 0;"));
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta1+mystr("[1] = 0;"));         
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta1+mystr("[2] = 0;"));	
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta1+mystr("[3] = 0;"));
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta1+mystr("[4] =-sin(")+name_beta+mystr("[0]);"));	
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta1+mystr("[5] = cos(")+name_beta+mystr("[0])*cos(")+name_beta+mystr("[1]);"));	
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta1+mystr("[6] = 0;"));
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta1+mystr("[7] =-cos(")+name_beta+mystr("[0]);"));	
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta1+mystr("[8] =-cos(")+name_beta+mystr("[1])*sin(")+name_beta+mystr("[0]);"));	
	break;
	case zxy:
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta1+mystr("[0] = 0;"));
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta1+mystr("[1] =-sin(")+name_beta+mystr("[0]);"));         
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta1+mystr("[2] = cos(")+name_beta+mystr("[0])*cos(")+name_beta+mystr("[1]);"));
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta1+mystr("[3] = 0;"));
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta1+mystr("[4] =-cos(")+name_beta+mystr("[0]);")); 
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta1+mystr("[5] =-cos(")+name_beta+mystr("[1])*sin(")+name_beta+mystr("[0]);"));	
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta1+mystr("[6] = 0;"));
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta1+mystr("[7] = 0;"));	
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta1+mystr("[8] = 0;"));
	break;
	case zxz:
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta1+mystr("[0] = 0;"));
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta1+mystr("[1] =-sin(")+name_beta+mystr("[0]);"));         
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta1+mystr("[2] = cos(")+name_beta+mystr("[0])*sin(")+name_beta+mystr("[1]);"));
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta1+mystr("[3] = 0;"));
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta1+mystr("[4] =-cos(")+name_beta+mystr("[0]);")); 
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta1+mystr("[5] =-sin(")+name_beta+mystr("[1])*sin(")+name_beta+mystr("[0]);"));	
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta1+mystr("[6] = 0;"));
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta1+mystr("[7] = 0;"));	
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta1+mystr("[8] = 0;"));
	break;
	}
	file_export.RWSoleStr(mystr("  return 0;"));
  file_export.RWSoleStr(mystr("}\n"));  
//DGbarDBeta2()
	file_export.RWSoleStr(mystr("int Calc")+name_dgbardbeta2+mystr("(void* i")+name_dgbardbeta2+mystr(", double* ")+name_beta+mystr(") // ")+mystr(rotseq));
	file_export.RWSoleStr(mystr("{"));  
	file_export.RWSoleStr(mystr("  double* ")+name_dgbardbeta2+mystr(" = (double*) i")+name_dgbardbeta2+mystr(";")); 
	switch(rotationsSequence)
	{
	case xyz:
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta2+mystr("[0] = 0;"));
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta2+mystr("[1] = 0;"));         
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta2+mystr("[2] =-cos(")+name_beta+mystr("[1]);"));	
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta2+mystr("[3] = 0;"));
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta2+mystr("[4] = 0;"));	
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta2+mystr("[5] =-sin(")+name_beta+mystr("[0])*sin(")+name_beta+mystr("[1]);"));	
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta2+mystr("[6] = 0;"));
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta2+mystr("[7] = 0;"));	
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta2+mystr("[8] =-cos(")+name_beta+mystr("[0])*sin(")+name_beta+mystr("[1]);"));	
	break;
	case zxy:
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta2+mystr("[0] = 0;"));
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta2+mystr("[1] = 0;"));         
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta2+mystr("[2] =-sin(")+name_beta+mystr("[0])*sin(")+name_beta+mystr("[1]);"));	
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta2+mystr("[3] = 0;"));
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta2+mystr("[4] = 0;"));	
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta2+mystr("[5] =-cos(")+name_beta+mystr("[0])*sin(")+name_beta+mystr("[1]);"));	
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta2+mystr("[6] = 0;"));
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta2+mystr("[7] = 0;"));	
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta2+mystr("[8] =-cos(")+name_beta+mystr("[1]);"));		
	break;
	case zxz:
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta2+mystr("[0] = 0;"));
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta2+mystr("[1] = 0;"));         
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta2+mystr("[2] = sin(")+name_beta+mystr("[0])*cos(")+name_beta+mystr("[1]);"));	
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta2+mystr("[3] = 0;"));
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta2+mystr("[4] = 0;"));	
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta2+mystr("[5] = cos(")+name_beta+mystr("[0])*cos(")+name_beta+mystr("[1]);"));	
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta2+mystr("[6] = 0;"));
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta2+mystr("[7] = 0;"));	
	file_export.RWSoleStr(mystr("  ")+name_dgbardbeta2+mystr("[8] =-sin(")+name_beta+mystr("[1]);"));		
	break;
	}
	file_export.RWSoleStr(mystr("  return 0;"));
  file_export.RWSoleStr(mystr("}\n")); 
//evaluation functions
	file_export.RWSoleStr(mystr("//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"));
	file_export.RWSoleStr(mystr("// main c - functions"));
//EvalF2
	file_export.RWSoleStr(mystr("void EvalF2(double* ")+name_linearterm+mystr(", double* ")+name_beta+mystr(", double* ")+name_betap+mystr(")")); 
	file_export.RWSoleStr(mystr("{"));  
	file_export.RWSoleStr(mystr("// variables"));
	file_export.RWSoleStr(mystr("  double ")+name_gbar+mystr("[3][3];"));
	file_export.RWSoleStr(mystr("  double ")+name_gbart+mystr("[3][3];"));
	file_export.RWSoleStr(mystr("  double ")+name_gbarp+mystr("[3][3];"));
	file_export.RWSoleStr(mystr("  double ")+name_dgbardbeta1+mystr("[3][3];"));
	file_export.RWSoleStr(mystr("  double ")+name_dgbardbeta2+mystr("[3][3];"));
	file_export.RWSoleStr(mystr("  double ")+name_m+mystr("[3][3];"));
	file_export.RWSoleStr(mystr("  double ")+name_rhs1+mystr("[3];"));
	file_export.RWSoleStr(mystr("  double ")+name_rhs2+mystr("[3];"));
	file_export.RWSoleStr(mystr("  double ")+name_dummymatrix+mystr("[3][3];"));
	file_export.RWSoleStr(mystr("  double ")+name_dummyvector+mystr("[3];"));
	file_export.RWSoleStr(mystr("  double ")+name_dummyvector2+mystr("[3];"));
	file_export.RWSoleStr(mystr("  double ")+name_dummydouble+mystr(";"));
	file_export.RWSoleStr(mystr(""));  
	file_export.RWSoleStr(mystr("// calculations:"));
//rh1
	file_export.RWSoleStr(mystr("// calculation of term ")+name_rhs1+mystr(" := [(")+name_gbart+mystr(" * ")+name_iphi+mystr(" * ")+name_gbarp+mystr(")^T + (..)] * ")+name_betap);       
	file_export.RWSoleStr(mystr("  Calc")+name_gbar+mystr("( ")+name_gbar+mystr(", ")+name_beta+mystr(");"));
	file_export.RWSoleStr(mystr("  MatrixTranspose( ")+name_gbart+mystr(", ")+name_gbar+mystr(", 3 ,3);"));
	file_export.RWSoleStr(mystr("  Calc")+name_dgbardbeta1+mystr("( ")+name_dgbardbeta1+mystr(", ")+name_beta+mystr(");"));
	file_export.RWSoleStr(mystr("  Calc")+name_dgbardbeta2+mystr("( ")+name_dgbardbeta2+mystr(", ")+name_beta+mystr(");"));
	file_export.RWSoleStr(mystr("  MatrixMultiplyScalar( ")+name_dgbardbeta1+mystr(", ")+name_dgbardbeta1+mystr(", 3, 3, ")+name_betap+mystr("[0]);"));
	file_export.RWSoleStr(mystr("  MatrixMultiplyScalar( ")+name_dgbardbeta2+mystr(", ")+name_dgbardbeta2+mystr(", 3, 3, ")+name_betap+mystr("[1]);"));
	file_export.RWSoleStr(mystr("  MatrixAddMatrix( ")+name_gbarp+mystr(", ")+name_dgbardbeta1+mystr(", 3, 3, ")+name_dgbardbeta2+mystr(", 3, 3);"));
	file_export.RWSoleStr(mystr("  MatrixMultiplyMatrix( ")+name_dummymatrix+mystr(", ")+name_iphi+mystr(", 3, 3, ")+name_gbarp+mystr(", 3, 3);"));
	file_export.RWSoleStr(mystr("  MatrixMultiplyMatrix( ")+name_m+mystr(", ")+name_gbart+mystr(", 3, 3, ")+name_dummymatrix+mystr(", 3, 3);"));
	file_export.RWSoleStr(mystr("  MatrixTranspose( ")+name_dummymatrix+mystr(", ")+name_m+mystr(", 3 ,3);"));
	file_export.RWSoleStr(mystr("  MatrixAddMatrix( ")+name_m+mystr(", ")+name_m+mystr(", 3, 3, ")+name_dummymatrix+mystr(", 3, 3);"));
	file_export.RWSoleStr(mystr(""));  
	file_export.RWSoleStr(mystr("  //MatrixMultiplyMatrix( ")+name_rhs1+mystr(", ")+name_m+mystr(", 3, 3, ")+name_betap+mystr(", 1, 3);"));
	file_export.RWSoleStr(mystr("  MatrixMultiplyVector( ")+name_rhs1+mystr(", ")+name_m+mystr(", 3, 3, ")+name_betap+mystr(", 3);"));
	file_export.RWSoleStr(mystr(""));  
//rh2
	file_export.RWSoleStr(mystr("// calculation of term ")+name_rhs2+mystr(" := ( (")+name_dgbardbeta1+mystr(" * ")+name_betap+mystr(") * v ), (")+name_dgbardbeta2+mystr(" * ")+
		name_betap+mystr(") * v ), 0 ) with v :=")+name_iphi+mystr(" * ")+name_gbar+mystr(" * ")+name_betap); 
	file_export.RWSoleStr(mystr("  MatrixMultiplyMatrix( ")+name_dummymatrix+mystr(", ")+name_iphi+mystr(", 3, 3, ")+name_gbar+mystr(", 3, 3);"));
	file_export.RWSoleStr(mystr("  //MatrixMultiplyMatrix( ")+name_dummyvector+mystr(", ")+name_dummymatrix+mystr(", 3, 3, ")+name_betap+mystr(", 1, 3);"));
	file_export.RWSoleStr(mystr("  MatrixMultiplyVector( ")+name_dummyvector+mystr(", ")+name_dummymatrix+mystr(", 3, 3, ")+name_betap+mystr(", 3);"));
	file_export.RWSoleStr(mystr("  Calc")+name_dgbardbeta1+mystr("( ")+name_dgbardbeta1+mystr(", ")+name_beta+mystr(");"));
	file_export.RWSoleStr(mystr("  Calc")+name_dgbardbeta2+mystr("( ")+name_dgbardbeta2+mystr(", ")+name_beta+mystr(");"));
	file_export.RWSoleStr(mystr("  //MatrixMultiplyMatrix( ")+name_dummyvector2+mystr(", ")+name_dgbardbeta1+mystr(", 3, 3, ")+name_betap+mystr(", 1, 3);"));
	file_export.RWSoleStr(mystr("  MatrixMultiplyVector( ")+name_dummyvector2+mystr(", ")+name_dgbardbeta1+mystr(", 3, 3, ")+name_betap+mystr(", 3);"));
	file_export.RWSoleStr(mystr("  //MatrixMultiplyMatrix( &")+name_dummydouble+mystr(", ")+name_dummyvector2+mystr(", 3, 1, ")+name_dummyvector+mystr(", 1, 3);"));
	file_export.RWSoleStr(mystr("  VectorMultiplyVector( &")+name_dummydouble+mystr(", ")+name_dummyvector2+mystr(", 3, ")+name_dummyvector+mystr(", 3);"));
 	file_export.RWSoleStr(mystr("  ")+name_rhs2+mystr("[0] = -")+name_dummydouble+mystr(";"));  
	file_export.RWSoleStr(mystr(""));  
	file_export.RWSoleStr(mystr("  //MatrixMultiplyMatrix( ")+name_dummyvector2+mystr(", ")+name_dgbardbeta2+mystr(", 3, 3, ")+name_betap+mystr(", 1, 3);"));
	file_export.RWSoleStr(mystr("  MatrixMultiplyVector( ")+name_dummyvector2+mystr(", ")+name_dgbardbeta2+mystr(", 3, 3, ")+name_betap+mystr(", 3);"));
	file_export.RWSoleStr(mystr("  //MatrixMultiplyMatrix( &")+name_dummydouble+mystr(", ")+name_dummyvector2+mystr(", 3, 1, ")+name_dummyvector+mystr(", 1, 3);"));
	file_export.RWSoleStr(mystr("  VectorMultiplyVector( &")+name_dummydouble+mystr(", ")+name_dummyvector2+mystr(", 3, ")+name_dummyvector+mystr(", 3);"));
	file_export.RWSoleStr(mystr("  ")+name_rhs2+mystr("[1] = -")+name_dummydouble+mystr(";"));  
	file_export.RWSoleStr(mystr(""));  
	file_export.RWSoleStr(mystr("  ")+name_rhs2+mystr("[2] = 0.0;"));  
	file_export.RWSoleStr(mystr(""));  
	file_export.RWSoleStr(mystr("  //MatrixAddMatrix( ")+name_linearterm+mystr(", ")+name_rhs1+mystr(", 1, 3, ")+name_rhs2+mystr(", 1, 3);"));
	file_export.RWSoleStr(mystr("  VectorAddVector( ")+name_linearterm+mystr(", ")+name_rhs1+mystr(", 3, ")+name_rhs2+mystr(", 3);"));
	file_export.RWSoleStr(mystr("}\n"));  
//EvalM
	file_export.RWSoleStr(mystr("void EvalM(double* ")+name_massmatrix+mystr(", double* ")+name_beta+mystr(", double* ")+name_betap+mystr(")")); 
	file_export.RWSoleStr(mystr("{"));
	file_export.RWSoleStr(mystr("// variables"));
	file_export.RWSoleStr(mystr("  double ")+name_gbar+mystr("[3][3];"));
	file_export.RWSoleStr(mystr("  double ")+name_gbart+mystr("[3][3];"));
	file_export.RWSoleStr(mystr("  double ")+name_m+mystr("[3][3];"));
	file_export.RWSoleStr(mystr("  double ")+name_dummymatrix+mystr("[3][3];"));
	file_export.RWSoleStr(mystr("  int i;"));
	file_export.RWSoleStr(mystr(""));  
	file_export.RWSoleStr(mystr("// calculations:"));
//M
	file_export.RWSoleStr(mystr("// calculation of term ")+name_m+mystr(" := Matrix[(diag(m))_3x3, (null)_3x3, (null)_3x3, (")+name_gbart+mystr(" * ")+name_iphi+mystr(" * ")+name_gbar+mystr(")_3x3]_6x6"));       
	file_export.RWSoleStr(mystr("  Calc")+name_gbar+mystr("( ")+name_gbar+mystr(", ")+name_beta+mystr(");"));
	file_export.RWSoleStr(mystr("  MatrixTranspose( ")+name_gbart+mystr(", ")+name_gbar+mystr(", 3 ,3);"));
	file_export.RWSoleStr(mystr("  MatrixMultiplyMatrix( ")+name_dummymatrix+mystr(", ")+name_iphi+mystr(", 3, 3, ")+name_gbar+mystr(", 3, 3);"));
	file_export.RWSoleStr(mystr("  MatrixMultiplyMatrix( ")+name_m+mystr(", ")+name_gbart+mystr(", 3, 3, ")+name_dummymatrix+mystr(", 3, 3);"));
	file_export.RWSoleStr(mystr(""));  
//diag elements: 0,7,14... Mphiphi 21-23,27-29,33-35
	file_export.RWSoleStr(mystr("  for (i= 1; i<21; i++) {")+name_massmatrix+mystr("[i] = 0.0;}"));
	file_export.RWSoleStr(mystr("  for (i=24; i<27; i++) {")+name_massmatrix+mystr("[i] = 0.0;}"));
	file_export.RWSoleStr(mystr("  for (i=30; i<33; i++) {")+name_massmatrix+mystr("[i] = 0.0;}"));
	file_export.RWSoleStr(mystr("  for (i= 0; i< 3; i++) {")+name_massmatrix+mystr("[i*7] = ")+name_mass+mystr(";}"));
	file_export.RWSoleStr(mystr("  for (i= 0; i< 3; i++) {")+name_massmatrix+mystr("[i+21] = ")+name_m+mystr("[i][0]; ")+name_massmatrix+mystr("[i+27] = ")+name_m+mystr("[i][1]; ")+name_massmatrix+mystr("[i+33] = ")+name_m+mystr("[i][2]; }"));
	file_export.RWSoleStr(mystr("}\n")); 
	file_export.RWSoleStr(mystr("// end of main c - function"));
	file_export.RWSoleStr(mystr("//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"));
	file_export.RWSoleStr(mystr("\n"));  
}

if(ansimath)
{	
	mystr filename(mystr("ansimath.h"));
	CMFile file_export(filename,TFMwrite);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ ansimath.h
	file_export.RWSoleStr(mystr("//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"));
	file_export.RWSoleStr(mystr("//+ ansimath.h - matrix and vector operations"));
	file_export.RWSoleStr(mystr("//+ automatically created by Rigid3DKardan::MatlabExportF2"));
	file_export.RWSoleStr(mystr("//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"));
	file_export.RWSoleStr(mystr(""));  
	file_export.RWSoleStr(mystr("// calculates maximum of two double numbers"));
	file_export.RWSoleStr(mystr("double Maximum(double a, double b) {return (a > b) ? a : b;}"));
	file_export.RWSoleStr(mystr(""));  
	file_export.RWSoleStr(mystr("// swaps rows of a double** matrix (swaps row pointers)"));
	file_export.RWSoleStr(mystr("void SwapRowsPtr(double** mat, int r1,int r2) {double* d = mat[r1]; mat[r1] = mat[r2]; mat[r2] = d;}"));
	file_export.RWSoleStr(mystr(""));  
	file_export.RWSoleStr(mystr("// multiplies a row of a double** matrix"));
	file_export.RWSoleStr(mystr("void MultRowPtr(double** mat, int r1, int cols, double fact) {int i; for(i=0;i<cols;i++) mat[r1][i]=mat[r1][i]*fact;}"));
	file_export.RWSoleStr(mystr(""));  
	file_export.RWSoleStr(mystr("// adds multiple of a row to another row of a double** matrix"));
	file_export.RWSoleStr(mystr("void AddRowVecPtr(double** mat, int fromRow, int toRow, double Fact, int fromCol, int toCol)"));
	file_export.RWSoleStr(mystr("  {int i; for (i=fromCol;i<=toCol;i++)	mat[toRow][i] = mat[toRow][i] + Fact*mat[fromRow][i];}"));
	file_export.RWSoleStr(mystr(""));  

	file_export.RWSoleStr(mystr("// transposes a double** matrix, result==source allowed"));  
	file_export.RWSoleStr(mystr("int MatrixTranspose(void* iresult, void* isource, int scol, int srow)"));
	file_export.RWSoleStr(mystr("{"));  
	file_export.RWSoleStr(mystr("  int i,j;"));
	file_export.RWSoleStr(mystr("  double * result = (double*) iresult;"));
	file_export.RWSoleStr(mystr("  double * source = (double*) isource;"));
	file_export.RWSoleStr(mystr("  double dummy;"));
	file_export.RWSoleStr(mystr("  if (result != source)"));
	file_export.RWSoleStr(mystr("  {"));  
	file_export.RWSoleStr(mystr("    for(i=0;i<srow;i++)"));  
	file_export.RWSoleStr(mystr("    {"));  
	file_export.RWSoleStr(mystr("      for(j=0;j<scol;j++)"));  
	file_export.RWSoleStr(mystr("      {"));  
	file_export.RWSoleStr(mystr("         result[j*srow+i] = source[i*scol+j];"));  
	file_export.RWSoleStr(mystr("      }"));  
	file_export.RWSoleStr(mystr("    }"));  
	file_export.RWSoleStr(mystr("  }"));  
	file_export.RWSoleStr(mystr("  else"));  
	file_export.RWSoleStr(mystr("  {"));  
	file_export.RWSoleStr(mystr("    for(i=0;i<srow;i++)"));  
	file_export.RWSoleStr(mystr("    {"));  
	file_export.RWSoleStr(mystr("      for(j=i+1;j<scol;j++)"));  
	file_export.RWSoleStr(mystr("      {"));  
	file_export.RWSoleStr(mystr("         dummy = source[i*scol+j];"));  
	file_export.RWSoleStr(mystr("         result[j*srow+i] = source[i*scol+j];"));  
	file_export.RWSoleStr(mystr("         source[i*scol+j] = dummy;"));  
	file_export.RWSoleStr(mystr("      }"));  
	file_export.RWSoleStr(mystr("    }"));  
	file_export.RWSoleStr(mystr("  }"));  
	file_export.RWSoleStr(mystr("  return 1;"));  
	file_export.RWSoleStr(mystr("}\n"));  

	file_export.RWSoleStr(mystr("// adds two double** matrices, result==source allowed"));
	file_export.RWSoleStr(mystr("int MatrixAddMatrix(void* iresult, void* isource1, int s1col, int s1row, void* isource2, int s2col, int s2row)"));
	file_export.RWSoleStr(mystr("{"));  
	file_export.RWSoleStr(mystr("  int i;"));  
	file_export.RWSoleStr(mystr("  double * result = (double*) iresult;"));
	file_export.RWSoleStr(mystr("  double * source1 = (double*) isource1;"));
	file_export.RWSoleStr(mystr("  double * source2 = (double*) isource2;"));
	file_export.RWSoleStr(mystr("  if (s1col != s2col ) return 0;"));  
	file_export.RWSoleStr(mystr("  if (s1row != s2row ) return 0;"));  
	file_export.RWSoleStr(mystr("  for(i=0;i<(s1col*s1row);i++) { result[i] = source1[i]+source2[i]; }")); 
	file_export.RWSoleStr(mystr("  return 1;"));  
	file_export.RWSoleStr(mystr("}\n"));  
	
	file_export.RWSoleStr(mystr("// adds two double* vectors, result==source allowed"));  
	file_export.RWSoleStr(mystr("int VectorAddVector(double* result, double* source1, int s1len, double* source2, int s2len)"));
	file_export.RWSoleStr(mystr("{"));  
	file_export.RWSoleStr(mystr("  int i;"));  
	file_export.RWSoleStr(mystr("  if (s1len != s2len ) return 0;"));  
	file_export.RWSoleStr(mystr("  for(i=0;i<(s1len);i++) { result[i] = source1[i]+source2[i]; }")); 
	file_export.RWSoleStr(mystr("  return 1;"));  
	file_export.RWSoleStr(mystr("}\n"));  

	file_export.RWSoleStr(mystr("// multiplies two double** matrices, checks matching condition, result==source NOT ALLOWED (to be changed...)"));  
  file_export.RWSoleStr(mystr("int MatrixMultiplyMatrix(void* iresult, void* isource1, int s1col, int s1row, void* isource2, int s2col, int s2row)"));
	file_export.RWSoleStr(mystr("{"));  
	file_export.RWSoleStr(mystr("  int i,j,k;"));  
	file_export.RWSoleStr(mystr("  double * result = (double*) iresult;"));
	file_export.RWSoleStr(mystr("  double * source1 = (double*) isource1;"));
	file_export.RWSoleStr(mystr("  double * source2 = (double*) isource2;"));
	file_export.RWSoleStr(mystr("  if ((result == source1) || (result ==source2)) return 0;"));  
	file_export.RWSoleStr(mystr("  if (s1col != s2row ) return 0;"));  
	file_export.RWSoleStr(mystr("  for(i=0;i<s1row;i++)"));  
	file_export.RWSoleStr(mystr("  {"));  
	file_export.RWSoleStr(mystr("    for(j=0;j<s2col;j++)"));  
	file_export.RWSoleStr(mystr("    {"));  
	file_export.RWSoleStr(mystr("      result[i*s2col+j] = 0.0;"));  
	file_export.RWSoleStr(mystr("      for(k=0;k<s1col;k++)"));  
	file_export.RWSoleStr(mystr("      {"));  
	file_export.RWSoleStr(mystr("        result[i*s2col+j] += source1[i*s1col+k]*source2[k*s2col+j];"));  
	file_export.RWSoleStr(mystr("      }"));  
	file_export.RWSoleStr(mystr("    }"));  
	file_export.RWSoleStr(mystr("  }"));  
	file_export.RWSoleStr(mystr("  return 1;"));  
	file_export.RWSoleStr(mystr("}\n"));  

	file_export.RWSoleStr(mystr("// multiplies a double** matrix with a double* vector (|), checks matching condition, result==source NOT ALLOWED (to be changed...)"));  
	file_export.RWSoleStr(mystr("int MatrixMultiplyVector(double* result, void* isource1, int s1col, int s1row, double* source2, int s2row)"));
	file_export.RWSoleStr(mystr("{"));  
	file_export.RWSoleStr(mystr("  int i,k;"));  
	file_export.RWSoleStr(mystr("  double * source1 = (double*) isource1;"));
	file_export.RWSoleStr(mystr("  if (result == source2) return 0;"));  
	file_export.RWSoleStr(mystr("  if (s1col != s2row ) return 0;"));  
	file_export.RWSoleStr(mystr("  for(i=0;i<s1row;i++) //(j==0)"));  
	file_export.RWSoleStr(mystr("  {"));  
	file_export.RWSoleStr(mystr("    result[i] = 0.0;"));  
	file_export.RWSoleStr(mystr("    for(k=0;k<s1col;k++)"));  
	file_export.RWSoleStr(mystr("    {"));  
	file_export.RWSoleStr(mystr("      result[i] += source1[i*s1col+k]*source2[k];"));  
	file_export.RWSoleStr(mystr("    }"));  
	file_export.RWSoleStr(mystr("  }"));  
	file_export.RWSoleStr(mystr("  return 1;"));  
	file_export.RWSoleStr(mystr("}\n")); 

	file_export.RWSoleStr(mystr("// multiplies a double* vector (-) with a double** matrix, checks matching condition, result==source NOT ALLOWED (to be changed...)"));  
	file_export.RWSoleStr(mystr("int VectorMultiplyMatrix(double* result, double* source1, int s1col, void* isource2, int s2col, int s2row)"));
	file_export.RWSoleStr(mystr("{"));  
	file_export.RWSoleStr(mystr("  int j,k;"));  
	file_export.RWSoleStr(mystr("  double * source2 = (double*) isource2;"));
	file_export.RWSoleStr(mystr("  if (result == source1) return 0;"));  
	file_export.RWSoleStr(mystr("  if (s1col != s2row ) return 0;"));  
	file_export.RWSoleStr(mystr("  for(j=0;j<s2col;j++) //(i==0)"));  
	file_export.RWSoleStr(mystr("  {"));  
	file_export.RWSoleStr(mystr("    result[j] = 0.0;"));  
	file_export.RWSoleStr(mystr("    for(k=0;k<s1col;k++)"));  
	file_export.RWSoleStr(mystr("    {"));  
	file_export.RWSoleStr(mystr("      result[j] += source1[k]*source2[k*s2col+j];"));  
	file_export.RWSoleStr(mystr("    }"));  
	file_export.RWSoleStr(mystr("  }"));  
	file_export.RWSoleStr(mystr("  return 1;"));  
	file_export.RWSoleStr(mystr("}\n"));  

	file_export.RWSoleStr(mystr("// multiplies two double* vectos, checks matching condition (scalar product)"));  
	file_export.RWSoleStr(mystr("int VectorMultiplyVector(double* result, double* source1, int s1len, double* source2, int s2len)"));
	file_export.RWSoleStr(mystr("{"));  
	file_export.RWSoleStr(mystr("  int i;"));  
	file_export.RWSoleStr(mystr("  if ((result == source1) || (result ==source2)) return 0;"));  
	file_export.RWSoleStr(mystr("  if (s1len != s2len ) return 0;"));  
	file_export.RWSoleStr(mystr("  *result = 0.0;"));  
	file_export.RWSoleStr(mystr("  for(i=0;i<s1len;i++)"));  
	file_export.RWSoleStr(mystr("  {"));  
	file_export.RWSoleStr(mystr("    *result += source1[i]*source2[i];"));  
	file_export.RWSoleStr(mystr("  }"));  
	file_export.RWSoleStr(mystr("  return 1;"));  
	file_export.RWSoleStr(mystr("}\n"));  

	file_export.RWSoleStr(mystr("// multiplies a double** matrix with a scalar, result==source allowed"));  
	file_export.RWSoleStr(mystr("int MatrixMultiplyScalar(void* iresult, void* isource, int scol, int srow, double factor)"));
	file_export.RWSoleStr(mystr("{"));  
	file_export.RWSoleStr(mystr("  int i;"));  
	file_export.RWSoleStr(mystr("  double * result = (double*) iresult;"));
	file_export.RWSoleStr(mystr("  double * source = (double*) isource;"));
	file_export.RWSoleStr(mystr("  for(i=0;i<(scol*srow);i++) { result[i] = source[i]*factor; }")); 
	file_export.RWSoleStr(mystr("  return 1;"));  
	file_export.RWSoleStr(mystr("}\n"));  

	file_export.RWSoleStr(mystr("// multiplies a double* vector with a scalar, result==source allowed"));  
	file_export.RWSoleStr(mystr("int VectorMultiplyScalar(double* result, double* source, int slen, double factor)"));
	file_export.RWSoleStr(mystr("{"));  
	file_export.RWSoleStr(mystr("  int i;"));  
	file_export.RWSoleStr(mystr("  for(i=0;i<(slen);i++) { result[i] = source[i]*factor; }")); 
	file_export.RWSoleStr(mystr("  return 1;"));  
	file_export.RWSoleStr(mystr("}\n"));  

	file_export.RWSoleStr(mystr("// inverts a (regular) matrix, result==source allowed"));  
	file_export.RWSoleStr(mystr("int MatrixInvert(void* result, void* source, int scol, int srow)"));
  file_export.RWSoleStr(mystr("{ // copied from linalg.cpp Matrix::Invert(), 1-based -> 0 based."));
  file_export.RWSoleStr(mystr("  int i,j,k,pivotpos;"));
	file_export.RWSoleStr(mystr("  double mij;"));
  file_export.RWSoleStr(mystr("  int n=srow;"));
	file_export.RWSoleStr(mystr("  int maxj=0;"));
	file_export.RWSoleStr(mystr("  double * r = (double*) result;"));
	file_export.RWSoleStr(mystr("  double * s = (double*) source;"));
	file_export.RWSoleStr(mystr("  double** a = (double**) malloc (scol*sizeof(double*)); // copy source to 2d array a"));
	file_export.RWSoleStr(mystr("  double ** m = (double**) malloc (srow*sizeof(double*)); // unity matrix"));
	file_export.RWSoleStr(mystr("  for (i=0; i<srow; i++)"));
	file_export.RWSoleStr(mystr("  {"));
	file_export.RWSoleStr(mystr("    a[i] = (double*) malloc (srow*(sizeof(double)));"));
	file_export.RWSoleStr(mystr("    for (j=0; j<scol; j++) {a[i][j] = s[i*scol+j];}"));
	file_export.RWSoleStr(mystr("  }"));
	file_export.RWSoleStr(mystr("  for (i=0; i<srow; i++)"));
	file_export.RWSoleStr(mystr("  {"));
  file_export.RWSoleStr(mystr("    m[i] = (double*) malloc (srow*(sizeof(double)));"));
	file_export.RWSoleStr(mystr("    for (j=0; j<srow; j++) {if(i==j) m[i][j] = 1.0; else m[i][j] = 0.0;}"));
	file_export.RWSoleStr(mystr("  }"));
	file_export.RWSoleStr(mystr("// rewrite from Matrix::Invert(): unreferenced calls are calls on \"a\" ++ zerobased!"));
	file_export.RWSoleStr(mystr("  for (j=0; j < n; j++)"));
	file_export.RWSoleStr(mystr("  {"));
	file_export.RWSoleStr(mystr("    double pivot=fabs(a[j][j]);"));
	file_export.RWSoleStr(mystr("    pivotpos=j;"));
	file_export.RWSoleStr(mystr("    for (k=j+1; k < n; k++)"));
	file_export.RWSoleStr(mystr("    {"));
	file_export.RWSoleStr(mystr("      if (fabs(a[k][j]) > pivot) {pivotpos=k; pivot=fabs(a[k][j]);}"));
	file_export.RWSoleStr(mystr("    }"));
	file_export.RWSoleStr(mystr("    if (pivot == 0.0) return 0;"));
	file_export.RWSoleStr(mystr("      maxj=(int)Maximum(pivotpos,maxj);"));
	file_export.RWSoleStr(mystr("    SwapRowsPtr(m,pivotpos,j);//m.SwapRows(pivotpos,j);"));
	file_export.RWSoleStr(mystr("    SwapRowsPtr(a,pivotpos,j);//SwapRows(pivotpos,j);"));
	file_export.RWSoleStr(mystr("    MultRowPtr(m,j,srow,(1.0/a[j][j])); // cols in m: srow! //m.MulRow(j,1/Get(j,j));"));
	file_export.RWSoleStr(mystr("    MultRowPtr(a,j,scol,(1.0/a[j][j]));//MulRow(j,1/Get(j,j));"));
	file_export.RWSoleStr(mystr("    for (i=j+1; i <n; i++)"));
	file_export.RWSoleStr(mystr("    {"));
	file_export.RWSoleStr(mystr("      mij = a[i][j];"));
	file_export.RWSoleStr(mystr("      if (mij !=0)"));
	file_export.RWSoleStr(mystr("      {"));
	file_export.RWSoleStr(mystr("        AddRowVecPtr(a,j,i,-mij,j,n-1);  //AddRowVec(j,i,-mij,j,n); //j..n"));
	file_export.RWSoleStr(mystr("        AddRowVecPtr(m,j,i,-mij,0,n-1);  //m.AddRowVec(j,i,-mij,1,maxj); //1..j"));
	file_export.RWSoleStr(mystr("      }"));
	file_export.RWSoleStr(mystr("    }"));
	file_export.RWSoleStr(mystr("  }"));
	file_export.RWSoleStr(mystr("  for (j=n-1; j > 0; j--) //backsubstitution, unfortunately this takes most of the time!!!"));
	file_export.RWSoleStr(mystr("  {"));
	file_export.RWSoleStr(mystr("    for (i=0; i < j; i++)"));
	file_export.RWSoleStr(mystr("    {"));
	file_export.RWSoleStr(mystr("      mij = a[i][j];"));
	file_export.RWSoleStr(mystr("      if (mij!=0) {AddRowVecPtr(m,j,i,-mij,0,n-1);} // m.AddRowVec(j+1,i,-mij); //1..n"));
	file_export.RWSoleStr(mystr("    }"));
	file_export.RWSoleStr(mystr("  }"));
	file_export.RWSoleStr(mystr("  for(i=0;i<srow;i++)"));
	file_export.RWSoleStr(mystr("  {"));
	file_export.RWSoleStr(mystr("    for(j=0;j<scol;j++) {r[i*scol+j] = m[i][j];}"));
  file_export.RWSoleStr(mystr("  }"));
	file_export.RWSoleStr(mystr("  return 1;"));
  file_export.RWSoleStr(mystr("}\n"));  
}
}
