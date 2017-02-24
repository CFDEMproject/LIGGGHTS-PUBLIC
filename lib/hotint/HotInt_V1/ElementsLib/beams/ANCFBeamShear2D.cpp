//#**************************************************************
//# filename:             ANCFBeamShear2D.cpp
//#
//# author:               Karin & Astrid
//#
//# generated:						2009 / 2010
//# description:          2D ANCF beam element with shear deformation
//# comments:							formulation is based on the following paper: 
//#												Gerstmayr, Matikainen, Mikkola: A geometrically exact beam element based on the ANCF
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
 
#include "rigid2d.h"
#include "femathhelperfunctions.h"
#include "material.h"
#include "ANCFBeamShear2D.h"
#include "Node.h"
//#include "graphicsconstants.h"
//#include "elementdataaccess.h"
//#include "solversettings_auto.h"

//SetANCFBeamShear-Fkt. for the linear element:
void ANCFBeamShear2D::SetANCFBeamShear2D(const Vector& xc1, const Vector& xc2, int n1i, int n2i, int materialnumi,
																				 const Vector3D& si, const Vector3D& coli)
{ 
	nnodes=2;
	n1=n1i; n2=n2i; sos2=0;
	size = si;
	n3=-1;

	//x_init = xc1.Append(xc2);
	//xg = xc1.Append(xc2);
	xg.SetLen(SOS());
	xg.SetAll(0.);
	x_init.SetLen(SOS());
	x_init.SetAll(0.);
	x_init = x_init.Append(Vector(SOS())); //Velocity initial conditions can also be transformed by Tinv!!!

	lx = size.X(); ly = size.Y(); lz = size.Z(); //lx is parameter to identify reference configuration
	materialnum = materialnumi;
	mass = lx*ly*lz*GetMaterial().Density(); //use lx, because this is the parameter of the real geometry ...
	//lx=1;ly=1;lz=1;
	col = coli;
	BuildDSMatrices();

	q0 = xc1.Append(xc2);
};

//SetANCFBeamShear-Fkt. for quadratic element:
void ANCFBeamShear2D::SetANCFBeamShear2D(const Vector& xc1, const Vector& xc2, const Vector& xc3, int n1i, int n2i, int n3i, int materialnumi,
																				 const Vector3D& si, const Vector3D& coli)
{
	nnodes=3;
	n1=n1i; n2=n2i;	n3=n3i;
	sos2=0;
	size = si;

	//x_init = xc1.Append(xc2);
	//xg = xc1.Append(xc2);
	xg.SetLen(SOS());
	xg.SetAll(0.);
	x_init.SetLen(SOS());
	x_init.SetAll(0.);
	x_init = x_init.Append(Vector(SOS())); //Velocity initial conditions can also be transformed by Tinv!!!

	lx = size.X(); ly = size.Y(); lz = size.Z(); //lx is parameter to identify reference configuration

	//size.X() = ComputeCurvedLength(x_init);			 //size.X() is true curve length, but other curve length is better for sliding joint!
	//UO() << "linit=" << lx << ", ltrue=" << size.X() << "\n";
	materialnum = materialnumi;

	mass = lx*ly*lz*GetMaterial().Density(); //use lx, because this is the parameter of the real geometry ...
	//lx=1;ly=1;lz=1;

	col = coli;
	//UO() << "Cable: rho=" << rho << ", E=" << Em << ", size=" << size 	<< ", n1=" << n1 << ", n2=" << n2 << "\n";

	BuildDSMatrices();

	//q0 is generated from vector xc1(4),... :q0=xc1.Append(xc2), q0.Append(xc3)
	q0.SetLen(SOS());
	q0.SetAll(0.);
	int i;

	//should be realized with append. ...
	for(i=1;i<=4;i++)
	{
		q0(i) = xc1(i);
	}
	for(i=5;i<=8;i++)
	{
		q0(i) = xc2(i-4);
	}
	for(i=9;i<=12;i++)
	{
		q0(i) = xc3(i-8);
	}
	//if (GetBeamEIy() == 0)  // if Beam parameters were not initialized in Material..
	//{
	//	GetMaterial().BeamEIy() = -1;
	//	GetMaterial().BeamEA() = 0;
	//	GetMaterial().BeamRhoA() = 0;
	//	GetMaterial().BeamA() = ly*lz;
	//}
};


void ANCFBeamShear2D::LinkToElements()
{
	if (SOSowned() == 0)
	{
		//UO() << "Link nodes to elements in Cable\n";
		const Node& node1 = GetMBS()->GetNode(n1);
		const Node& node2 = GetMBS()->GetNode(n2);
		for (int i=1; i <= node1.SOS(); i++)
		{
			AddLTG(node1.Get(i)); //LocalToGlobal-Element-list
		}
		for (int i=1; i <= node2.SOS(); i++)
		{
			AddLTG(node2.Get(i));
		}
		if(nnodes == 3)
		{
			const Node& node3 = GetMBS()->GetNode(n3);
			for (int i=1; i <= node3.SOS(); i++)
			{
				AddLTG(node3.Get(i));
			}
		}
		for (int i=1; i <= node1.SOS(); i++)
		{
			AddLTG(node1.Get(i+node1.SOS()));
		}
		for (int i=1; i <= node2.SOS(); i++)
		{
			AddLTG(node2.Get(i+node2.SOS()));
		}
		if(nnodes == 3)
		{
			const Node& node3 = GetMBS()->GetNode(n3);
			for (int i=1; i <= node3.SOS(); i++)
			{
				AddLTG(node3.Get(i+node3.SOS()));
			}
		}
	}
}

void ANCFBeamShear2D::BuildDSMatrices()
{
	//orderx = 9; //max 9x5, otherwise array grad too small!!!!
	//GetIntegrationRule(x1,w1,orderx);
};


//----------------
//Shapefunctions
//
//defined on scaled rectangular element: ploc=(\xi,\eta)
//-L/2<\xi<+L/2, -H/2<\eta<+H/2
//----------------

//GetShapes(sf,ploc) computes Shapefunctions at position ploc and saves values in vector sf
void ANCFBeamShear2D::GetShapes(Vector& sf, const Vector2D& ploc) const  //ploc=(x,y)
{ 
	if(nnodes == 2)
	{
		sf.SetLen(NS());  //NS=2*nnodes, vector sf has 4 entries
		sf(1)=1./lx*(lx/2.-ploc(1));
		sf(2)=ploc(2)/lx*(lx/2.-ploc(1));
		sf(3)=1./lx*(lx/2.+ploc(1));
		sf(4)=ploc(2)/lx*(lx/2.+ploc(1));
	}
	else if(nnodes == 3)
	{
		sf.SetLen(NS());  //NS=2*nnodes, vector sf has 6 entries
		sf(1)=-2./(lx*lx)*ploc(1)*(lx/2.-ploc(1));
		sf(2)=-2./(lx*lx)*ploc(1)*ploc(2)*(lx/2.-ploc(1));
		sf(3)=2./(lx*lx)*ploc(1)*(lx/2.+ploc(1));
		sf(4)=2./(lx*lx)*ploc(1)*ploc(2)*(lx/2.+ploc(1));
		sf(5)=-4./(lx*lx)*(ploc(1)-lx/2.)*(ploc(1)+lx/2.);
		sf(6)=-4./(lx*lx)*ploc(2)*(ploc(1)-lx/2.)*(ploc(1)+lx/2.);
	}
}

//GetSF(i,ploc) computes i-th Shapefunctions at ploc=(x,y)
double ANCFBeamShear2D::GetSF(int i, const Vector2D& ploc) const
{ 
	if(nnodes == 2)
	{
		switch(i)
		{
		case 1: return 1./lx*(lx/2.-ploc(1));
		case 2: return ploc(2)/lx*(lx/2.-ploc(1));
		case 3: return 1./lx*(lx/2.+ploc(1));
		case 4: return ploc(2)/lx*(lx/2.+ploc(1));
		default: mbs->UO()<<"only 4 Shapefunctions\n"; return 0.; //in case i>4 
		}
	}
	else if(nnodes == 3)
	{
		switch(i)
		{
		case 1: return -2./(lx*lx)*ploc(1)*(lx/2.-ploc(1));
		case 2: return -2./(lx*lx)*ploc(1)*ploc(2)*(lx/2.-ploc(1));
		case 3: return 2./(lx*lx)*ploc(1)*(lx/2.+ploc(1));
		case 4: return 2./(lx*lx)*ploc(1)*ploc(2)*(lx/2.+ploc(1));
		case 5: return -4./(lx*lx)*(ploc(1)-lx/2.)*(ploc(1)+lx/2.);
		case 6: return -4./(lx*lx)*ploc(2)*(ploc(1)-lx/2.)*(ploc(1)+lx/2.);
		default: mbs->UO()<<"only 6 Shapefunctions\n"; return 0.; //fin case i>6
		}
	}
	return 0.;
}
//computes derivatives dS/dx at ploc
void ANCFBeamShear2D::GetShapesx(Vector& sfx, const Vector2D& ploc) const
{
	if(nnodes == 2)
	{
		sfx.SetLen(NS());
		sfx(1)=-1./lx;
		sfx(2)=-ploc(2)/lx;
		sfx(3)=1./lx;
		sfx(4)=ploc(2)/lx;
	}
	else if(nnodes == 3)
	{
		sfx.SetLen(NS());
		sfx(1)=-1./lx+4.*ploc(1)/(lx*lx);
		sfx(2)=-ploc(2)/lx+4.*ploc(1)*ploc(2)/(lx*lx);
		sfx(3)=1./lx+4.*ploc(1)/(lx*lx);
		sfx(4)=ploc(2)/lx+4.*ploc(1)*ploc(2)/(lx*lx);
		sfx(5)=-8.*ploc(1)/(lx*lx);
		sfx(6)=-8.*ploc(1)*ploc(2)/(lx*lx);
	}
}

double ANCFBeamShear2D::GetSFx(int i, const Vector2D& ploc) const
{
	if(nnodes == 2)
	{
		switch(i)
		{
		case 1: return -1./lx;
		case 2: return -ploc(2)/lx;
		case 3: return 1./lx;
		case 4: return ploc(2)/lx;
		default: mbs->UO()<<"only 4 Shapefunctions\n"; return 0.;
		}
	}
	else if(nnodes == 3)
	{
		switch(i)
		{
		case 1: return -1./lx+4.*ploc(1)/(lx*lx);
		case 2: return -ploc(2)/lx+4.*ploc(1)*ploc(2)/(lx*lx);
		case 3: return 1./lx+4.*ploc(1)/(lx*lx);
		case 4: return ploc(2)/lx+4.*ploc(1)*ploc(2)/(lx*lx);
		case 5: return -8.*ploc(1)/(lx*lx);
		case 6: return -8.*ploc(1)*ploc(2)/(lx*lx);
		default: mbs->UO()<<"only 6 Shapefunctions\n"; return 0.;
		}
	}
	return 0.;
}

//2nd derivatives d^2S/dx^2 are 0, GetShapesxx creates vector sfxx with zeros
void ANCFBeamShear2D::GetShapesxx(Vector& sfxx, const Vector2D& ploc) const
{
	if(nnodes == 2)
	{
		sfxx.SetLen(NS());  //NS=4
		sfxx.SetAll(0);     //all derivatives are 0
	}
	else if(nnodes == 3)
	{
		sfxx.SetLen(NS());  //NS=6
		sfxx(1)=4./(lx*lx);
		sfxx(2)=4.*ploc(2)/(lx*lx);
		sfxx(3)=4./(lx*lx);
		sfxx(4)=4.*ploc(2)/(lx*lx);
		sfxx(5)=-8./(lx*lx);
		sfxx(6)=-8.*ploc(2)/(lx*lx);
	}
}

double ANCFBeamShear2D::GetSFxx(int i, const Vector2D& ploc) const
{
	if(nnodes == 2)
	{
		return 0;
	}
	else if(nnodes == 3)
	{
		switch(i)
		{
		case 1: return 4./(lx*lx);
		case 2: return 4.*ploc(2)/(lx*lx);
		case 3: return 4./(lx*lx);
		case 4: return 4.*ploc(2)/(lx*lx);
		case 5: return -8./(lx*lx);
		case 6: return -8.*ploc(2)/(lx*lx);
		default: mbs->UO()<<"only 6 Shapefunctions\n"; return 0.;
		}
	}
	return 0.;
}

//computes derivatives dS/dy at ploc
void ANCFBeamShear2D::GetShapesy(Vector& sfy, const Vector2D& ploc) const
{
	if(nnodes == 2)
	{
		sfy.SetLen(NS());
		sfy(1)=0;
		sfy(2)=1./lx*(lx/2.-ploc(1));
		sfy(3)=0;
		sfy(4)=1./lx*(lx/2.+ploc(1));
	}
	else if(nnodes == 3)
	{
		sfy.SetLen(NS());
		sfy(1)=0;
		sfy(2)=-ploc(1)/lx+2*ploc(1)*ploc(1)/(lx*lx);
		sfy(3)=0;
		sfy(4)=ploc(1)/lx+2*ploc(1)*ploc(1)/(lx*lx);
		sfy(5)=0;
		sfy(6)=1.-4.*ploc(1)*ploc(1)/(lx*lx);
	}
}

double ANCFBeamShear2D::GetSFy(int i, const Vector2D& ploc) const
{
	if(nnodes == 2)
	{
		switch(i)
		{
		case 1: return 0;
		case 2: return 1./lx*(lx/2.-ploc(1));
		case 3: return 0;
		case 4: return 1./lx*(lx/2.+ploc(1));
		default: mbs->UO()<<"only 4 Shapefunctions\n"; return 0.;
		}
	}
	else if(nnodes == 3)
	{
		switch(i)
		{
		case 1: return 0;
		case 2: return -ploc(1)/lx+2*ploc(1)*ploc(1)/(lx*lx);
		case 3: return 0;
		case 4: return ploc(1)/lx+2*ploc(1)*ploc(1)/(lx*lx);
		case 5: return 0;
		case 6: return 1.-4.*ploc(1)*ploc(1)/(lx*lx);
		default: mbs->UO()<<"only 6 Shapefunctions\n"; return 0.;
		}
	}
	return 0.;
}

//computes mixed derivatives d^2S/dxdy at ploc
void ANCFBeamShear2D::GetShapesxy(Vector& sfy, const Vector2D& ploc) const
{
	if(nnodes == 2)
	{
		sfy.SetLen(NS());
		sfy(1)=0;
		sfy(2)=-1./lx;
		sfy(3)=0;
		sfy(4)=1./lx;
	}
	else if(nnodes == 3)
	{
		sfy.SetLen(NS());
		sfy(1)=0;
		sfy(2)=-1./lx+4.*ploc(1)/(lx*lx);
		sfy(3)=0;
		sfy(4)=1./lx+4.*ploc(1)/(lx*lx);
		sfy(5)=0;
		sfy(6)=-8.*ploc(1)/(lx*lx);
	}
}

double ANCFBeamShear2D::GetSFxy(int i, const Vector2D& ploc) const
{
	if(nnodes == 2)
	{
		switch(i)
		{
		case 1: return 0;
		case 2: return -1./lx;
		case 3: return 0;
		case 4: return 1./lx;
		default: mbs->UO()<<"only 4 Shapefunctions\n"; return 0.;
		}
	}
	else if(nnodes == 3)
	{
		switch(i)
		{
		case 1: return 0;
		case 2: return -1./lx+4.*ploc(1)/(lx*lx);
		case 3: return 0;
		case 4: return 1./lx+4.*ploc(1)/(lx*lx);
		case 5: return 0;
		case 6: return -8.*ploc(1)/(lx*lx);
		default: mbs->UO()<<"only 6 Shapefunctions\n"; return 0.;
		}
	}
	return 0.;
}

void ANCFBeamShear2D::Gradu(const Vector2D& ploc, const Vector& u, Matrix3D& gradu) const
{

	Matrix3D jac, jacinv;
	GetJacobi(jac,ploc,q0);

	jac.GetInverse(jacinv);
	jacinv = jacinv.GetTp();

	gradu.SetSize(2,2);
	gradu.SetAll(0);

	int dim = Dim();
	int l;
	for (int j = 1; j <= dim; j++) 
	{
		for (int i = 1; i <= NS(); i++)
		{
			l = (i-1)*dim+j;
				gradu(j,1) += GetSFx(i,ploc)*u(l);
				gradu(j,2) += GetSFy(i,ploc)*u(l);
		}
	}
}

//Kappa
void ANCFBeamShear2D::GetDeltaKappa(const double& x, const Vector& xg, Vector& dkappa, double& kappa) const
{

	int dim = Dim();
	int ns = NS();

	static Vector SVx;
	static Vector SVxx;
	SVx.SetLen(NS());
	SVxx.SetLen(NS());

	Vector2D ploc(x,0.);
	GetShapesx(SVx,ploc);
	GetShapesxx(SVxx,ploc);

	Vector2D rx(0.,0.);
	rx = GetPosx2D(ploc);

	Vector2D rxx(0.,0.);
	rxx = GetPos2D(ploc);

	double rxn = rx.Norm();

	if (rxn == 0) 
	{
		dkappa *= 0; kappa=0; return;
	}  
	double rxcrxx = rx.Cross(rxx);
	double f = rxcrxx;

	double g = Sqr(rxn);

	if (rxn == 0) 
	{
		kappa=0; //f/0
	}
	else 
	{
		kappa = f/g; 
	}

	double g2inv = 1./Sqr(g);
	double fn = f*g2inv;
	fn *= 2; //==cable element kappa mode = 2

	double gn = g*g2inv;

	double t1;
	double df, dg;

	for (int i=1; i <= dim; i++)
	{
		for (int j=1; j <= ns; j++)
		{
			//dr,x/de x r,xx + r,x x dr,xx/de: 
			//delta r,x x d^2(r)/dx^2 + ... =>paper 
			switch (i) {
					case 1: {	t1 =  SVx(j)*rxx.Y()-SVxx(j)*rx.Y(); break; }
					case 2: {	t1 = -SVx(j)*rxx.X()+SVxx(j)*rx.X(); break; }
					default: ;
			}

			dg = (rx(i)*SVx(j)); //normed
			df = t1;
			dkappa((j-1)*dim+i) = df*gn-fn*dg; //normed
		}
	}
}

//----------------
//Velocity
//----------------

//flagD = 0 für Berechnungen, flagD = 1 for Visualization

Vector2D ANCFBeamShear2D::GetVel2D(const Vector2D& p_loc, int flagD) const
{
	Vector2D p(0.,0.);
	for (int i = 1; i <= Dim(); i++) //i=1,2
	{
		for (int j = 1; j <= NS(); j++)//j=1,...,4
		{	
			if(flagD==0)
				p(i) += GetSF(j,p_loc)*XGP((j-1)*Dim()+i);  //XGP=actual solution vector

			else
				p(i) += GetSF(j,p_loc)*XGPD((j-1)*Dim()+i);
		}
	}
	return p;
};

//----------------
//Displacement
//----------------
//only for visualization
Vector2D ANCFBeamShear2D::GetDisplacement2DD(const Vector2D& p_loc) const
{
	Vector2D p(0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += GetSF(j,p_loc)*(XGD((j-1)*Dim()+i));
		}
	}
	return p;
};

//----------------
//GetPos
//----------------

//GetPos2D berechnet Vektor Positionsvektor r im deformierten Element
//GetPos2D computes r in deformed element
//xg = generalized coordinates (= our unknowns q)
Vector2D ANCFBeamShear2D::GetPos2D(const Vector2D& p_loc, const Vector& xg) const
{
	Vector2D p(0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += GetSF(j,p_loc)*(xg((j-1)*Dim()+i)+q0((j-1)*Dim()+i));  //r=u+r_0=Sq+Sq_0 (q_0=initial values, x_init=0-Vector -> no displacement in the beginning))
		}
	}
	return p;
};

//GetPos2D von oben, aber falls im Aufruf als zweite Eingabe kein Vektor, sondern ein Integer steht,
//wird diese Fkt für Grafik aufgerufen; globale Variable XG wird für Rechnung verwendet anstatt Eingabe xg

//same GetPos2D function as above, but with the second parameter flag D the function for the visualization is activated
Vector2D ANCFBeamShear2D::GetPos2D(const Vector2D& p_loc, int flagD) const
{
	Vector2D p(0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			if(flagD==0)
				p(i) += GetSF(j,p_loc)*(XG((j-1)*Dim()+i)+q0((j-1)*Dim()+i));
			else
				p(i) += GetSF(j,p_loc)*(XGD((j-1)*Dim()+i)+q0((j-1)*Dim()+i));
		}
	}
	return p;
};

//Ableitungen von GetPos entsprechen unseren Richtungsvektoren
//GetPosx2D (Ableitung in Richtung der Balkenachse: dr/dx=r,x)
//GetPosx2D is the derivative with respect to the direction of the beam centerline
Vector2D ANCFBeamShear2D::GetPosx2D(const Vector2D& p_loc, const Vector& xg) const
{
	Vector2D p(0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += GetSFx(j,p_loc)*(xg((j-1)*Dim()+i)+q0((j-1)*Dim()+i));  //r,x=S,x*u+S,x*q_0
		}
	}
	return p;
}

Vector2D ANCFBeamShear2D::GetPosx2D(const Vector2D& p_loc, int flagD) const
{
	Vector2D p(0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			if(flagD==0)
				p(i) += GetSFx(j,p_loc)*(XG((j-1)*Dim()+i)+q0((j-1)*Dim()+i));
			else
				p(i) += GetSFx(j,p_loc)*(XGD((j-1)*Dim()+i)+q0((j-1)*Dim()+i));
		}
	}
	return p;
}

//GetPosy2D (Ableitung in Richtung des Querschnitts: dr/dy=r,y)
//GetPosy2D is the derivative with respect to the direction of the cross section
Vector2D ANCFBeamShear2D::GetPosy2D(const Vector2D& p_loc, const Vector& xg) const
{
	Vector2D p(0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += GetSFy(j,p_loc)*(xg((j-1)*Dim()+i)+q0((j-1)*Dim()+i)); //r,y=S,y*u+S,y*q_0
		}
	}
	return p;
}

Vector2D ANCFBeamShear2D::GetPosy2D(const Vector2D& p_loc, int flagD) const
{
	Vector2D p(0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			if(flagD==0)
				p(i) += GetSFy(j,p_loc)*(XG((j-1)*Dim()+i)+q0((j-1)*Dim()+i));
			else
				p(i) += GetSFy(j,p_loc)*(XGD((j-1)*Dim()+i)+q0((j-1)*Dim()+i));
		}
	}
	return p;
}

//not used at the moment:
//second derivatives of position vector r with respect to x

//Vector2D ANCFBeamShear2D::GetPosxx2D(const Vector2D& p_loc, const Vector& xg) const
//	{
//		Vector2D p(0.,0.);
//		for (int i = 1; i <= Dim(); i++)
//		{
//			for (int j = 1; j <= NS(); j++)
//			{
//				p(i) += GetSFxx(j,p_loc)*(xg((j-1)*Dim()+i)+q0((j-1)*Dim()+i));
//			}
//		}
//		return p;
//	}
//
//Vector2D ANCFBeamShear2D::GetPosxx2D(const Vector2D& p_loc, int flagD) const
//	{
//		Vector2D p(0.,0.);
//		for (int i = 1; i <= Dim(); i++)
//		{
//			for (int j = 1; j <= NS(); j++)
//			{
//				if(flagD==0)
//				p(i) += GetSFxx(j,p_loc)*(XG((j-1)*Dim()+i)+q0((j-1)*Dim()+i));
//				else
//				p(i) += GetSFxx(j,p_loc)*(XGD((j-1)*Dim()+i)+q0((j-1)*Dim()+i));
//			}
//		}
//		return p;
//	}


//GetPosxy2D: d2r/dydx=r,xy
//is used in DeltaThetax
Vector2D ANCFBeamShear2D::GetPosxy2D(const Vector2D& p_loc, const Vector& xg) const
{
	Vector2D p(0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += GetSFxy(j,p_loc)*(xg((j-1)*Dim()+i)+q0((j-1)*Dim()+i));
		}
	}
	return p;
}

Vector2D ANCFBeamShear2D::GetPosxy2D(const Vector2D& p_loc, int flagD) const
{
	Vector2D p(0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			if(flagD==0)
				p(i) += GetSFxy(j,p_loc)*(XG((j-1)*Dim()+i)+q0((j-1)*Dim()+i));
			else
				p(i) += GetSFxy(j,p_loc)*(XGD((j-1)*Dim()+i)+q0((j-1)*Dim()+i));
		}
	}
	return p;
}


//computes position vector in the reference element: r_0 
Vector2D ANCFBeamShear2D::GetInitPos2D(const Vector2D& p_loc) const
{
	Vector2D p(0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += GetSF(j,p_loc)*(q0((j-1)*Dim()+i));  //r_0=Sq_0; (q_0 includes initial values)
		}
	}
	return p;
};

//for visualization
Vector2D ANCFBeamShear2D::GetPos2D0D(const Vector2D& p_loc) const 
{
	Vector2D plocscaled;
	plocscaled(1)=p_loc(1)*lx/2.;
	plocscaled(2)=p_loc(2)*ly/2.;
	return GetPos2D(plocscaled,1); 
};

Vector2D ANCFBeamShear2D::GetPos2D0D(const Vector2D& p_loc, double defscale) const 
{
	Vector2D plocscaled;
	plocscaled(1)=p_loc(1)*lx/2.;
	plocscaled(2)=p_loc(2)*ly/2.;
	return GetInitPos2D(plocscaled) + defscale*GetDisplacement2DD(plocscaled); 
};


//----------------
//directions T1, T2
//----------------
//T1 steht normal auf Querschnittsrichtung T2
//T1 is perpendicular to the direction of the cross section T2
Vector2D ANCFBeamShear2D::GetT12D(const double& xbar) const
{
	Vector2D T1(0.,0.), T2(0.,0.);
	Vector2D ploc(xbar,0.);
	T2=GetPosy2D(ploc);  //r,y
	T2.Normalize();      //T1 has length 1
	T1(1)=T2(2);         //Rechts-Kippregel: vector T1 is perpendicular to T2
	T1(2)=-T2(1);
	return T1;
}

//direction of cross section: T2=r,y/Norm(r,y)
Vector2D ANCFBeamShear2D::GetT22D(const double& xbar) const
{
	Vector2D T2(0.,0.);
	Vector2D ploc(xbar,0.);
	T2=GetPosy2D(ploc);  //r,y
	T2.Normalize();      //T2 has length 1
	return T2;
}


//----------------
//Gamma, Theta'
//----------------
//Gamma1,Gamma2 (and corresponding Delta-Variations)
double ANCFBeamShear2D::GetGamma12D(const double& xbar) const
{
	Vector2D ploc(xbar,0.); 
	return GetT12D(xbar)*GetPosx2D(ploc)-1;  //Gamma1 = T1.r,x-1  (r,x=direction of beam axis, in Paper denoted by r_0
}

double ANCFBeamShear2D::GetGamma22D(const double& xbar) const
{
	Vector2D ploc(xbar,0.);
	return GetT22D(xbar)*GetPosx2D(ploc);    //Gamma2 = T2.r,x
}

void ANCFBeamShear2D::GetDeltaGamma22D(const double& xbar, Vector& DeltaGamma2) const
{
	//Paper-equation (56)
	//DeltaGamma2 = r,y.S,x/Norm + r,x.S,y/Norm - (r,x.r,y).r,y.S,y/Norm^3

	DeltaGamma2.SetLen(SOS());
	DeltaGamma2.SetAll(0);

	Vector2D ry, rx, ploc(xbar,0.);
	rx=GetPosx2D(ploc);
	ry=GetPosy2D(ploc);

	double ryNorm, ry3;
	ryNorm=ry.Norm();
	ry3=ryNorm*ryNorm*ryNorm;

	if(ryNorm<1e-10)
		return;  //if denominator gets 0

	for (int i = 1; i <= Dim(); i++)  //i=1,2
	{
		for (int j = 1; j <= NS(); j++)//j=1,...,4
		{
			int k=(j-1)*Dim()+i;
			DeltaGamma2(k)=ry(i)*GetSFx(j,ploc)/ryNorm + rx(i)*GetSFy(j,ploc)/ryNorm - (rx*ry)*ry(i)*GetSFy(j,ploc)/ry3;
		}
	}
}

void ANCFBeamShear2D::GetDeltaGamma12D(const double& xbar, Vector& DeltaGamma1) const
{
	//Paper-equation (57)
	//DeltaGamma1 = \hat{r,y}.S,x/Norm - \hat{r,x}.S,y/Norm + (\hat{r,x}.r,y).r,y.S,y/Norm^3

	DeltaGamma1.SetLen(SOS());
	DeltaGamma1.SetAll(0);

	Vector2D ry, rx, ploc(xbar,0.), rxhat, ryhat;
	rx=GetPosx2D(ploc);
	ry=GetPosy2D(ploc);
	rxhat(1)=rx(2);
	rxhat(2)=-rx(1);
	ryhat(1)=ry(2);
	ryhat(2)=-ry(1);

	double ryNorm, ry3;
	ryNorm=ry.Norm();
	ry3=ryNorm*ryNorm*ryNorm;

	if(ryNorm<1e-10)
		return;  //if denominator gets 0

	for (int i = 1; i <= Dim(); i++)  //i=1,2
	{
		for (int j = 1; j <= NS(); j++)//j=1,...,4
		{
			int k=(j-1)*Dim()+i;
			DeltaGamma1(k)=ryhat(i)*GetSFx(j,ploc)/ryNorm - rxhat(i)*GetSFy(j,ploc)/ryNorm + (rxhat*ry)*ry(i)*GetSFy(j,ploc)/ry3;
		}
	}
}

//angle (for curvature): Theta'
double ANCFBeamShear2D::GetThetax2D(const double& xbar) const
{
	//Paper-equation (50)
	//Theta'= (r,y Cross r,xy)/Norm^2

	Vector2D ploc(xbar,0.);
	double f,g;
	f=GetPosy2D(ploc).Cross(GetPosxy2D(ploc));
	g=sqr(GetPosy2D(ploc).Norm());
	if(fabs(g)<1e-10)
		return 0;  //if denominator gets 0
	else
		return f/g;
}

void ANCFBeamShear2D::GetDeltaThetax2D(const double& xbar, Vector& DeltaThetax) const
{
	//Paper-equations (54) und (55)

	DeltaThetax.SetLen(SOS());
	DeltaThetax.SetAll(0);

	double f,g;
	Vector2D ploc(xbar,0);
	f=GetPosy2D(ploc).Cross(GetPosxy2D(ploc));
	g=sqr(GetPosy2D(ploc).Norm());
	if(fabs(g)<1e-10)
		return;  //if denominator gets 0

	Vector Deltaf, Deltag;
	Deltaf.SetLen(SOS());
	Deltaf.SetAll(0);
	Deltag.SetLen(SOS());
	Deltag.SetAll(0);

	Vector2D ry, rxy;
	ry=GetPosy2D(ploc);
	rxy=GetPosxy2D(ploc);

	//Deltaf(k)=(r,y Cross S,xy) - (r,xy Cross S,y)
	//Deltag(k)=2*r,y*S,y

	for (int i = 1; i <= Dim(); i++)  //i=1,2
	{
		for (int j = 1; j <= NS(); j++)//j=1,...,4
		{
			int k=(j-1)*Dim()+i; //k = coordinate in Deltaf-vector
			if(i==1)
				Deltaf(k)=-ry(2)*GetSFxy(j,ploc)+rxy(2)*GetSFy(j,ploc);  
			else
				Deltaf(k)=ry(1)*GetSFxy(j,ploc)-rxy(1)*GetSFy(j,ploc);

			Deltag(k)=2.*ry(i)*GetSFy(j,ploc);
		}
	}

	DeltaThetax=(g*Deltaf-f*Deltag);
	DeltaThetax*=1./sqr(g);
}

//is used for \Pi^{thickness}:
//Eyy = transverse strain
double ANCFBeamShear2D::GetEyy2D(const double& xbar) const
{
	Vector2D ploc(xbar,0);
	Vector2D ry;
	ry=GetPosy2D(ploc);
	double Eyy;
	Eyy = 0.5*(ry*ry)-0.5;  //equation (51), Fehler: -1/2 statt -1, (bei Ableitung DeltaEyy wirkt sich dieser Faktor nicht aus)
	return Eyy;
}

//DeltaEyy
void ANCFBeamShear2D::GetDeltaEyy2D(const double& xbar, Vector& DeltaEyy) const
{
	Vector2D ploc(xbar,0);
	Vector2D ry;
	ry=GetPosy2D(ploc);

	for (int i = 1; i <= Dim(); i++)  //i=1,2
	{
		for (int j = 1; j <= NS(); j++)//j=1,...,4
		{
			int k=(j-1)*Dim()+i;
			DeltaEyy(k)=ry(i)*GetSFy(j,ploc);   //DeltaEyy = 0.5*2*r,y*Delta(r,y)=r,y.S,y
		}
	}
}

double ANCFBeamShear2D::GetAngle2D(const Vector2D& ploc) const
{
	mbs->UO() << "called GetAngle2D";
	//static Vector xg(SOS()); 
	//xg.SetLen(SOS());

	//for (int i = 1; i <= SOS(); i++)
	//	xg(i) = XG(i);  //e=e1..e8

	//Vector2D rx = GetPosx2D(ploc.X(),xg);
	//double ang = atan2(rx.Y(),rx.X());

	return 0.;//ang;
}
//
//double ANCFBeamShear2D::GetAngle2DP(const Vector2D& ploc) const
//{
//	static Vector xg(SOS()); 
//	xg.SetLen(SOS());
//
//	for (int i = 1; i <= SOS(); i++)
//		xg(i) = XG(i);  //e=e1..e8
//
//	double p0=ploc.X()/(0.5*lx);
//	static Vector SVx;
//
//	GetShapesx(SVx, p0);
//	Vector2D rx(0.,0.);
//	for (int i = 1; i <= Dim(); i++)
//	{
//		for (int j = 1; j <= NS(); j++)
//		{
//			rx(i) += SVx(j)*xg((j-1)*Dim()+i);
//		}
//	}
//
//	double rxn = Sqr(rx.Norm());
//	if (rxn == 0) rxn = 1;
//
//	double angp = 0;
//	for (int j = 1; j <= NS(); j++)
//	{
//		angp += XGP((j-1)*Dim()+1) * (-rx.Y()*SVx(j))/rxn;
//		angp += XGP((j-1)*Dim()+2) * (SVx(j)*rx.X())/rxn;
//	}
//
//	return angp;
//}

void ANCFBeamShear2D::GetdPosdqT(const Vector2D& ploc, Matrix& d)
{
	d.SetSize(NS()*Dim(),Dim());//12 x 2   (eigentlich:   6*(2x2))
	d.FillWithZeros();

	for (int i = 1; i <= NS(); i++)
	{
		d((i-1)*Dim()+1,1)=GetSF(i,ploc);
		d((i-1)*Dim()+2,2)=GetSF(i,ploc);
	}
}

//H = int(S,dV,V)
void ANCFBeamShear2D::GetH(Matrix& H) 
{
	if (Hmatrix.Getrows() == SOS())
	{
		H = Hmatrix;
		return;
	}
	else
	{
		int dim = Dim();
		int ns = NS();

		H.SetSize(ns*dim,dim);
		H.SetAll(0);
		Matrix3D jac;

		ConstVector<1> x2, w2;
		GetIntegrationRule(x1,w1,3); //3x1x1 !!!!!
		GetIntegrationRule(x2,w2,1);

		for (int i1=1; i1<=x1.GetLen(); i1++)
		{
			for (int i2=1; i2<=x2.GetLen(); i2++)
			{
				Vector2D ploc(x1(i1)*lx*0.5,x2(i2)*ly*0.5);
				GetJacobi(jac,ploc,q0);
				double jacdet = jac.Det() *0.25*lx*ly;
				//mbs->UO() << "jacobian = " << jac << "\njacdet = " << jacdet << "\n";
				double fact = fabs (jacdet) * w1(i1)*w2(i2) * lz;

				for (int i=0; i<ns; i++)
				{
					for (int j=1; j<=dim; j++)
					{
						H(i*dim+j,j)+=fact*GetSF(i+1, ploc);
					}
				}
			}
		}
		Hmatrix = H;
	}
}


//----------------
//Massmatrix M
//----------------
//M = int(rho*((S)^T).S, dV,V)
void ANCFBeamShear2D::EvalM(Matrix& m, double t)
{
	if (massmatrix.Getcols() == SOS())
	{
		m = massmatrix;
		return;
	}
	else
	{
		int dim = Dim();
		int ns = NS();

		Vector x1,w1;
		int order=4;

		GetIntegrationRule(x1,w1,order);

		//we have a general cross section; instead of Em() we use Rho()
		double rhoI_w = Rho()*lz*Cub(ly)/12.; 
		double A = ly*lz;										   
		double rhoA = Rho()*A;                 

		double Shape_i, Shape_j;
		double Shapey_i, Shapey_j;

		for (int i1=1; i1<=x1.GetLen(); i1++)//sum over integration points
		{
			double xi = x1(i1);    //unit element
			double x = xi*lx*0.5;  //scaled rect. reference element: \xi
			Vector2D ploc(x,0);

			Vector2D rx0 = GetPosx2D(ploc);
			double rxn = rx0.Norm(); //curvature in reference element
			double det = 0.5*lx*rxn; //element transformation

			for (int i_dim=1; i_dim<=dim; i_dim++)
			{
				for(int i_ns=1; i_ns<=ns; i_ns++)
				{
					int i=(i_ns-1)*dim+i_dim;
					Shape_i=GetSF(i_ns,ploc);
					Shapey_i=GetSFy(i_ns,ploc);

					for (int j_dim=1; j_dim<=dim; j_dim++)
					{
						for(int j_ns=1; j_ns<=ns; j_ns++)
						{ 
							int j=(j_ns-1)*dim+j_dim;
							Shape_j=GetSF(j_ns,ploc);
							Shapey_j=GetSFy(j_ns,ploc);
							if(j_dim==i_dim)//i not equal j->Si*Sj=0 
							{
								m(i,j)+=w1(i1)*det*rhoA*Shape_i*Shape_j+w1(i1)*det*rhoI_w*Shapey_i*Shapey_j;
							}
						}
					}
				}
			}

		}
		massmatrix = m;
		//UO() << "m-invert=" << m.Invert2() << "\n";
	}  
};
//last entry: number of rows, number of columns -> 0

//variational formulation: Delta\Pi (equation (52))
//we calculate residual
void ANCFBeamShear2D::EvalF2(Vector& f, double t)
{
	Body2D::EvalF2(f,t);
	TMStartTimer(22);  //CPU timing
	int sos = SOS();
	ConstVector<12> fadd;//element residual

	int ns = NS();
	int dim = Dim();

	fadd.SetLen(SOS());
	fadd.SetAll(0);

	ConstVector<12> temp;
	temp.SetLen(SOS());
	temp.SetAll(0);

	SV.SetLen(NS());

	//choose formulation:
	//Paper: A geom. exact beam element based on the ANCF (Gerstmayr, Matikainen, Mikkola)
	//int BeamFormulation = 0;//Finnland Paper: Simo&Vu-Quoc
	//int BeamFormulation = 1;//continuum mechanics based formulation (locking)
	//int BeamFormulation = 2;//continuum mechanics based formulation (locking compensation)

	double EI_w = Em()*lz*Cub(ly)/12.;   //Flächenträgheitsmoment: moment of inertia of area
	double A = ly*lz;										 //Querschnittsfläche: cross section area
	double EA = Em()*A;                  //E-Modul*area
	//Shear-Modul: G = E/(2+2nu)
	double GAks = A*ks*Em()/(2.+2.*Nu());
	//double GAks = A*Em()/(2.+2.*Nu());//GA (without ks)
	SetComputeCoordinates();

	if(BeamFormulation == 0)//Finnland Paper
	{

		//if (IsBeamParameters())
		//{
		//	EI_w = GetBeamEIy();
		//	EA = GetBeamEA();
		//	GAks = GetMaterial().BeamGAky();
		//}
		//mbs->UO()<<"EI="<<EI_w<<", EA="<<EA<<", GAks="<<GAks<<"\n";
		//mbs->UO()<<"Nu="<<Nu()<<"\n";

		//Integration order:
		//for Gauß Integration:
		int	order_Gamma, order_Theta;
		//for Lobatto Integration:
		int order_Thickness;

		//int_order
		if(nnodes == 2)
		{
			if(RI == 0)
			{
				order_Gamma = 4;
				order_Theta = 4;
			}
			else if(RI ==1)
			{
				order_Gamma = 1;
				order_Theta = 1;
			}
			order_Thickness =2;

		}
		else if(nnodes == 3)
		{
			if(RI == 0)
			{
				order_Gamma = 6;
				order_Theta = 6;
			}
			else if(RI == 1)
			{
				order_Gamma = 3;
				order_Theta = 3;
			}
			order_Thickness = 3;  //this gives order 4 in Lobatto
		}
		//Lobatto: 
		// order=1 and order=2 -> both cases 2 Integration points on the borders
		// order=3 -> 3 integration points, 2 on the borders, 1 on the mid point

		static Vector xG,xT,wG,wT,xE,wE; //xE,wE for Lobatto Integration

		//Integration points and weights:
		GetIntegrationRule(xG,wG,order_Gamma); 
		GetIntegrationRule(xT,wT,order_Theta);    
		GetIntegrationRuleLobatto(xE,wE, order_Thickness);  
		//mbs->UO()<<"Integrationspunkte Lobatto:"<<xE<<"\n";//How many Lobatto-points do we have?
		SetComputeCoordinates();  //values from global XG -> xg

		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//Integration
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		static Vector x0;
		x0.SetLen(SOS());
		x0.SetAll(0.);

		//Simo-Vu-Quoc-Strain:

		//Theta*DeltaTheta (bending stiffness):
		for (int i1=1; i1<=xT.GetLen(); i1++)  //i1-loop over all integration points
		{
			double xi = xT(i1);    //unit element
			double x = xi*lx*0.5;  //scaled rect. reference element: xi
			Vector2D ploc(x,0);

			Vector2D rx0 = GetPosx2D(ploc, x0);
			double rxn = rx0.Norm(); //Vorkrümmung im Referenzelement: curvature in reference element
			double det = 0.5*lx*rxn; //element transformation

			double fact_EI= EI_w*wT(i1)*det; //Det from element transformation

			GetDeltaThetax2D(x,temp);  //get DeltaTheta'-Werte and write in temp-vector

			temp *= fact_EI*GetThetax2D(x); //temp=EI*Theta'*DeltaTheta'
			//mbs->UO()<<"temp="<<temp<<"\n";
			fadd += temp;
		}

		//Gamma1*DeltaGamma1 (axial stiffness):
		for (int i1=1; i1<=xG.GetLen(); i1++)
		{
			double xi = xG(i1);    //unit element
			double x = xi*lx*0.5;  //scaled rect. reference element: xi
			Vector2D ploc(x,0);

			Vector2D rx0 = GetPosx2D(ploc, x0);
			double rxn = rx0.Norm(); //Vorkrümmung im Referenzelement: curvature in reference element
			double det = 0.5*lx*rxn; //element transformation

			double fact_EA= EA*wG(i1)*det; //Det from element transformation

			GetDeltaGamma12D(x,temp);  //get DelatGamma1-Werte and write in temp-vector

			temp *= fact_EA*GetGamma12D(x); //temp=EA*Gamma1*DeltaGamma1
			fadd += temp;
		}

		//Gamma2*DeltaGamma2 (shear stiffness):
		for (int i1=1; i1<=xG.GetLen(); i1++)
		{
			double xi = xG(i1);    //unit element
			double x = xi*lx*0.5;  //scaled rect. reference element: xi
			Vector2D ploc(x,0);

			Vector2D rx0 = GetPosx2D(ploc, x0);
			double rxn = rx0.Norm(); //Vorkrümmung im Referenzelement: curvature in reference element
			double det = 0.5*lx*rxn; //element transformation

			double fact_GAks= GAks*wG(i1)*det; //Det from element transformation
			//double fact_GAks= 1./(1./GAks+(lx*lx)/(12*EI_w))*wG(i1)*det; //Substitution from Felippa-document

			GetDeltaGamma22D(x,temp);  //get DeltaGamma2-Werte and write in temp-vector

			temp *= fact_GAks*GetGamma22D(x); //temp=ks*GA*Gamma2*DeltaGamma2
			fadd += temp;//fadd=fadd+temp
		}

		//thickness-Strain:

		//Eyy*DeltaEyy:
		for (int i1=1; i1<=xE.GetLen(); i1++)
		{
			double xi = xE(i1);    //unit element
			double x = xi*lx*0.5;  //scaled rect. reference element: xi
			Vector2D ploc(x,0);

			Vector2D rx0 = GetPosx2D(ploc, x0);
			double rxn = rx0.Norm();
			double det = 0.5*lx*rxn;
			//double update, Eyy, DeltaEyy;
			double fact_Eyy= EA*wE(i1)*det;

			GetDeltaEyy2D(x,temp);  //get DeltaEyy-Werte and write in temp-vector

			temp *= fact_Eyy*GetEyy2D(x); //temp=EA*Eyy*DeltaEyy
			fadd += temp;

			//linearized thickness-strain:
			//linearisiertes Eyy:

			//if(i1==1)
			//{
			//Eyy=GetSFy(2,ploc)*XG(4);  //XG = Vektor der FG
			//DeltaEyy=GetSFy(2,ploc);
			//update = fact_Eyy*Eyy*DeltaEyy; //temp=EA*Eyy(lin)*DeltaEyy(lin)
			//fadd(4) += update;
			//}
			//else
			//{
			//Eyy=GetSFy(4,ploc)*XG(8);
			//DeltaEyy=GetSFy(4,ploc);
			//update = fact_Eyy*Eyy*DeltaEyy; //temp=EA*Eyy(lin)*DeltaEyy(lin)
			//fadd(8) += update;
			//}

		}
	}
	//= 1;//continuum mechanics based formulation (locking)
	//= 2;//continuum mechanics based formulation (locking compensation)
	else if(BeamFormulation == 1 || BeamFormulation == 2)
	{

		//Integration Order
		int order_Poisson;
		int order_Thickness;

		if(nnodes == 2)
		{
			if(RI == 0)
			{
				order_Poisson = 4;
			}
			else if(RI ==1)
			{
				order_Poisson = 1;
			}
			order_Thickness =2;

		}
		else if(nnodes == 3)
		{
			if(RI == 0)
			{
				order_Poisson = 6;
			}
			else if(RI == 1)
			{
				order_Poisson = 3;
			}
			order_Thickness = 3;
		}
	//mbs->UO()<<"order="<<order_Poisson<<"\n";
	//mbs->UO()<<"order="<<order_Thickness<<"\n";

	//UO() << "standard ANCForig\n";

		Matrix3D Id;
		Id.Set22(1.,0.,0.,1.);

		Matrix3D strain, piola1, F;
		F.SetSize(2,2);

		temp.SetLen(SOS());
		fadd.SetLen(SOS());
		fadd.SetAll(0);

		int poissoncorrection = 0; //1==reduced integration of poisson part
		if (BeamFormulation == 2) poissoncorrection = 1;
		
		static Vector x1,x2,w1,w2;

		for (int kk=1; kk <= 1+poissoncorrection; kk++)
		{
			Matrix3D Dm;

			double ks = 10.*(1.+Nu())/(12.+11.*Nu()); //for rectangular cross-section only!

			if (poissoncorrection)
			{
				if (BeamFormulation == 2)
				{
					Dm.SetAll(0);
					if (kk == 1)//D^0 (equ. (25))
					{
						//poisson ratio zero:
						Dm(1,1) = Em();
						Dm(2,2) = Em();
						Dm(3,3) = ks*Em()/(2.*(1.+Nu()));//=ks*G

						//GetDMatrix(Dm,nu,Em);
						GetIntegrationRule(x1,w1,order_Poisson);//x_i ... integration points, w_i ... integration weights
						GetIntegrationRule(x2,w2,order_Thickness);
					}
					else if (kk == 2)//D^v (equ. (26))
					{
						Dm(1,1) = Em()/(1.-Sqr(Nu()))-Em();//=Em*nu^2/(1-nu^2)
						Dm(1,2) = Nu()*Em()/(1.-Sqr(Nu()));
						Dm(2,1) = Nu()*Em()/(1.-Sqr(Nu()));
						Dm(2,2) = Em()/(1.-Sqr(Nu()))-Em();

						GetIntegrationRule(x1,w1,order_Poisson);//x_i ... integration points, w_i ... integration weights
						GetIntegrationRule(x2,w2,1);//only 1 int. point along cross section
					}
				}
			}
			else
			{
				GetIntegrationRule(x1,w1,order_Poisson);//x_i ... integration points, w_i ... integration weights
				GetIntegrationRule(x2,w2,order_Thickness);
				double nu2 = Nu();

				double f = Em()/(1.-Sqr(nu2));

				Dm(1,1)=f;
				Dm(1,2)=nu2*f;
				Dm(1,3)=0;
				
				Dm(2,1)=nu2*f;
				Dm(2,2)=f;
				Dm(2,3)=0;
				
				Dm(3,1)=0;
				Dm(3,2)=0;
				Dm(3,3)=.5*f*(1.-nu2);
			}

			int kx1 = x2.Length();//number of integration points

			for (int i1=1; i1<=x1.GetLen(); i1++)
			{
				for (int i2=1; i2<=x2.GetLen(); i2++)
				{
					//TMStartTimer(20);
					int i,j,k;
					int ind = (i1-1)*kx1+(i2-1);
					Vector2D p(x1(i1)*0.5*lx,x2(i2)*0.5*ly);

					// compute F
					//F.SetAll(0);
					//F(1,1)=GetPosx2D(p).X();
					//F(2,1)=GetPosx2D(p).Y();
					//F(1,2)=GetPosy2D(p).X();
					//F(2,2)=GetPosy2D(p).Y();

					Matrix3D jac;
					GetJacobi(jac,p,q0);
					double jacdet = jac.Det();
					
					Matrix3D jacinv;
					jac.GetInverse(jacinv);
					jacinv = jacinv.GetTp();

					static Matrix grad0, grad;//am Einheitselement
					grad0.SetSize(Dim(),NS());
					for(i=1; i<=NS(); i++)
					{
						grad0(1,i)=GetSFx(i,p);
						grad0(2,i)=GetSFy(i,p);
					}
					Mult(jacinv, grad0, grad);
					Gradu(p,xg,F);//Verschiebungsgrad
					F(1,1) += 1;//Positionsgrad
					F(2,2) += 1;

					int linear=0;
					if (linear)
					{
						F -= Id;
						strain = 0.5*(F+F.GetTp());
						Vector3D s3(strain(1,1), strain(2,2), 2*strain(1,2));
						Vector3D stress = Dm*s3;
						piola1.Set22(stress.X(),stress.Z(),stress.Z(),stress.Y()); //linearized!
					}
					else if (BeamFormulation == 1 || BeamFormulation == 2)
					{
						//Green-Lagrange strain tensor
						//strain = 0.5 * (F.GetTp() * F - I);

						//F.GetATA2(strain); //does not exist in 2D
						strain = 0.5 * (F.GetTp() * F);
						strain(1,1) -= 0.5; strain(2,2) -= 0.5;

						//this formula does not work for plane stress/strain:
						//piola1 = F * ((2*mu) * strain + (la * strain.Trace())*Id );
						Vector3D s3(strain(1,1), strain(2,2), 2.*strain(1,2));
						Vector3D stress;
						stress = Dm * s3;

						piola1.Set22(stress.X(),stress.Z(),stress.Z(),stress.Y()); //linearized!

						piola1 = F * piola1;//2.PK
					}

					for (int j=1; j <= dim; j++)
					{
						for (int i = 0; i < ns; i++)
						{
							//temp(dim*i+j) = grad[ind](1, i+1)*piola1(j,1) + grad[ind](2, i+1)*piola1(j,2);
							temp(dim*i+j) = grad(1, i+1)*piola1(j,1) + grad(2, i+1)*piola1(j,2);
						}
					}

					//fadd.MultAdd(fabs (jacdet[ind]) * lz * w1(i1)*w2(i2),temp);
					fadd.MultAdd(fabs(jacdet)*0.5*lx*0.5*ly * lz * w1(i1)*w2(i2),temp);
				}
			}
		}	

	}
	//mbs->UO()<<"fadd="<<fadd<<"\n";
	//Residual:
	f -= fadd;  //f=f-fadd, (Ku=f -> Residual=f-Ku, Ku=fadd)

	//mbs->UO()<<"general.coord.: "<<XG(1)<<","<<XG(4)<<"\n";

	TMStopTimer(22);

	//if (GetMassDamping() != 0)
	//{
	//	// +++++++ damping: +++++++
	//	for (int i = 1; i <= SOS(); i++)
	//		xg(i) = XGP(i);

	//	if (massmatrix.Getcols() == SOS())
	//	{
	//		Mult(massmatrix,xg,temp);
	//	}
	//	else
	//	{
	//		Matrix dmat;
	//		dmat.SetSize(SOS(),SOS());
	//		EvalM(dmat,t);
	//		massmatrix = dmat;
	//		Mult(dmat,xg,temp);
	//	}
	//	//double k=1;
	//	//if (GetMBS()->GetTime() < 0.4) k=500;
	//	temp *= GetMassDamping();
	//	f -= temp;
	//}

};



//++++++++++++++++++++++++ Plasticity ++++++++++++++++++++++

double ANCFBeamShear2D::PostNewtonStep(double t)		// changes plastic strains, returns yieldfunction(sigma)/Em
{
	return 0;
}

void ANCFBeamShear2D::PostprocessingStep()
{
}


//--------------------------------------------------------
//for Displacement/Stress/Strain/Stress Resultants-Visualization:
//---------------------------------------------------------

void ANCFBeamShear2D::GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables)
{
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_displacement,
		FieldVariableDescriptor::FVCI_z);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_stress,
		FieldVariableDescriptor::FVCI_z, true);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_stress_mises);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_total_strain,
		FieldVariableDescriptor::FVCI_y, true);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_force_axial);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_force_transversal);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_moment_bending);
}

double ANCFBeamShear2D::GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector2D & local_position, bool flagD)
{
	if(!flagD)
		return FIELD_VARIABLE_NO_VALUE;

	//Displacements:
	//ploc_xi_0 is from -1 .. +1
	Vector2D ploc;
	ploc(1)=local_position(1)*lx*0.5;  //ploc_xi_0  -> ploc_\xi
	ploc(2)=local_position(2)*ly*0.5;  //ploc_eta_0 -> ploc_\eta

	xgd.SetLen(SOS());
	GetDrawCoordinates(xgd);  //xgd=XGD

	if(fvd.VariableType() == FieldVariableDescriptor::FVT_displacement)
		return fvd.GetComponent(GetDisplacement2DD(ploc));			// displacements

	//Stress and Strain:

	double Thetax = GetThetax2D(ploc(1));  //später: event. hier Vorkrümmung abziehen: -Vorkrümmung;
	double Gamma1 = GetGamma12D(ploc(1));
	double Gamma2 = GetGamma22D(ploc(1));

	Matrix3D strain(0);
	Vector2D rx, ry;
	rx=GetPosx2D(ploc, xgd);
	ry=GetPosy2D(ploc, xgd);
	strain(1,1) = 0.5*(rx(1)*rx(1)+rx(2)*rx(2)-1);
	strain(2,1) = 0.5*(rx(1)*ry(1)+rx(2)*ry(2));
	strain(1,2) = strain(2,1);
	strain(2,2) = 0.5*(ry(1)*ry(1)+ry(2)*ry(2)-1);

	Matrix3D stress(0);
	double fact=Em()/(1.+Nu())/(1.-2.*Nu());
	stress(1,1)=fact*((1-Nu())*strain(1,1)+Nu()*strain(2,2));
	stress(2,2)=fact*(Nu()*strain(1,1)+(1-Nu())*strain(2,2));
	stress(1,2)=fact*(1-2*Nu())*strain(1,2);
	stress(2,1)=stress(1,2);
	stress(3,3)=fact*(Nu()*strain(1,1)+Nu()*strain(2,2));

	double A = ly*lz;
	double EI_w = Em()*lz*Cub(ly)/12.;
	double EA = Em()*A; 
	// for Gamma2:
	// G=E/(2+2nu)
	double GAks = ks*Em()/(2+2*Nu());

	switch(fvd.VariableType())
	{
	case FieldVariableDescriptor::FVT_stress_mises:							return stress.Mises();
	case FieldVariableDescriptor::FVT_stress:										return fvd.GetComponent(stress);
	case FieldVariableDescriptor::FVT_total_strain:							return fvd.GetComponent(strain);
	case FieldVariableDescriptor::FVT_beam_force_axial:					return EA * Gamma1;
	case FieldVariableDescriptor::FVT_beam_force_transversal:		return GAks * Gamma2;
	case FieldVariableDescriptor::FVT_beam_moment_bending:			return EI_w * Thetax;
	}

	return FIELD_VARIABLE_NO_VALUE;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++Änderungen bis hier++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ANCFBeamShear2D::DrawElement()
{
	mbs->SetColor(col);

	// draw elements at z = 0.1
	Vector3D offset(0,0,0.);

	double lx1 = lx; double ly1 = ly*GetMBS()->GetMagnifyYZ();

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
	if (GetMBS()->GetActualPostProcessingFieldVariable() != NULL)
		colormode = 1;

	double def_scale = GetMBS()->GetDOption(105); // deformation scaling

	int drawlinemode = GetMBS()->GetIOption(117); // Draw Flat Elements, in this case, function is plotted normal to element axis

	if ( (1||colormode) && !drawlinemode) // plot coloured elements
	{
		lx1 /= lx;
		ly1 /= ly;

		double modeval = 0;
		int xgset = 0;

		double tilex = GetMBS()->GetIOption(137);
		double tiley = GetMBS()->GetIOption(138);

		TArray<Vector3D> points((int)(tilex+1)*(int)(tiley+1));
		TArray<double> vals((int)(tilex+1)*(int)(tiley+1));
		double v=0;

		points.SetLen(0); vals.SetLen(0);
		Vector2D p0, p0_mag, vx, vy, vy_mag;
		int tileyn = (int)tiley;
		int tilexn = (int)tilex;

		p0 = Vector2D(-lx1,-1.);
		vx = Vector2D(2.*lx1/tilexn,0);
		vy = Vector2D(0,2./tileyn);

		p0_mag = Vector2D(-lx1,-ly1);
		vy_mag = Vector2D(0,2.*ly1/tileyn);

		for (double iy = 0; iy <= tileyn+1e-10; iy++)
		{
			for (double ix = 0; ix <= tilexn+1e-10; ix++)
			{
				Vector2D ploc = (p0+ix*vx+iy*vy);
				Vector2D ploc_mag = (p0_mag+ix*vx+iy*vy_mag);
				Vector2D pg = GetPos2D0D(ploc_mag, def_scale);

				Vector3D p(ToP3D(pg)+offset);
				points.Add(p);
				if (colormode)
					v = GetFieldVariableValue(*GetMBS()->GetActualPostProcessingFieldVariable(), ploc, true);
				vals.Add(v);
			}
		}
		mbs->DrawColorQuads(points,vals,(int)tilexn+1,(int)tileyn+1,colormode,linemode);
	}
	else if (drawlinemode) // plot function on center line over beam axis
	{
		double scalex = lx1/lx;
		double scaley = ly1/ly;

		int tilex = GetMBS()->GetIOption(137);

		TArray<Vector3D> points(tilex+1), normals(tilex+1), drawpoints(tilex+1);
		TArray<double> vals(tilex+1);
		double v=0;

		points.SetLen(0); vals.SetLen(0); normals.SetLen(0); drawpoints.SetLen(0);
		Vector2D p0, vx, vy;

		p0 = Vector2D(-scalex,0);
		vx = Vector2D(2.*scalex/tilex,0);
		vy = Vector2D(0.,0.);

		int cnt = 0;
		for (double ix = 0.; ix <= tilex+1e-10; ix++)
		{
			// point on the middle line
			Vector2D ploc = (p0+ix*vx);
			Vector2D pg = GetPos2D0D(ploc, def_scale);

			Vector3D p(ToP3D(pg)+offset);
			points.Add(p);

			// normal direction (numerical differentiation for the moment)
			Vector2D ploc_eps = ploc + 1e-5*vx;
			Vector2D pg_eps = GetPos2D0D(ploc_eps, def_scale);

			Vector3D normal = Vector3D(-pg_eps.Y()+pg.Y(), pg_eps.X()-pg.X(), 0.);
			double scale = ly1/normal.Norm();
			normal *= scale;
			normals.Add(normal);

			// computation of value
			if (colormode)
				v = GetFieldVariableValue(*GetMBS()->GetActualPostProcessingFieldVariable(), ploc, true);
			vals.Add(v);

			MBS* constmbs = const_cast<MBS*>(mbs);
			constmbs->UpdateFEMinMaxCol(v);
		}

		for (int cnt=1; cnt<=tilex+1; cnt++)
		{
			// scale normal according to minimum/maximum
			double scale = Maximum(fabs(GetMBS()->GetTImincol()), fabs(GetMBS()->GetTImaxcol()));
			if (scale==0) scale = 1.;
			// point on the function line
			drawpoints.Add(points(cnt)+(vals(cnt)/scale)*normals(cnt));

			// draw arrow for value of function, in normal direction to middle line, colored matching legend
			double diff = mbs->GetFEmaxcol()-mbs->GetFEmincol();
			Vector3D fecolor = mbs->FEColor( (vals(cnt)-mbs->GetFEmincol())/diff );
			mbs->MyDrawArrow(points(cnt), drawpoints(cnt), fecolor);
			if (cnt>1)
			{
				// red line for plotting function value
				mbs->MyDrawLine(drawpoints(cnt-1), drawpoints(cnt), ly1*1e2, col);
				// middle line of element
				mbs->MyDrawLine(points(cnt-1), points(cnt), ly1*5e1);
				// upper and lower line of element, unscaled
				mbs->MyDrawLine(points(cnt-1)+(0.5/scaley)*normals(cnt-1), points(cnt)+(0.5/scaley)*normals(cnt), ly1*5e1);
				mbs->MyDrawLine(points(cnt-1)-(0.5/scaley)*normals(cnt-1), points(cnt)-(0.5/scaley)*normals(cnt), ly1*5e1);
			}

		}

		TArray<Vector3D> nodepoints(4);
		nodepoints(1) = ToP3D(GetPos2D0D(Vector2D(-scalex, 1.), def_scale)) + offset;
		nodepoints(2) = ToP3D(GetPos2D0D(Vector2D(-scalex,-1.), def_scale)) + offset;
		nodepoints(3) = ToP3D(GetPos2D0D(Vector2D( scalex,-1.), def_scale)) + offset;
		nodepoints(4) = ToP3D(GetPos2D0D(Vector2D( scalex, 1.), def_scale)) + offset;
		//mbs->DrawPolygon(nodepoints, 1, ly*ly1*1e-1);
		// left and right line of element
		mbs->MyDrawLine(nodepoints(1), nodepoints(2), ly1*4e1);
		mbs->MyDrawLine(nodepoints(3), nodepoints(4), ly1*4e1);
	}
	else
	{
		double tiling = 24;

		Vector2D p1,p2,p3,p4;
		mbs->SetColor(col);

		for (double i = 0; i < tiling; i++)
		{
			double l1 = -lx1*0.5+lx1*i/tiling;
			double l2 = -lx1*0.5+lx1*(i+1)/tiling;
			if (i == 0)
			{
				p1 = GetPos2D0D(Vector2D(l1/(0.5*lx), 1.), def_scale);
				p2 = GetPos2D0D(Vector2D(l1/(0.5*lx),-1.), def_scale);

			}
			else
			{
				p1 = p4;
				p2 = p3;
			}
			p3 = GetPos2DD(Vector2D(l2,-ly*0.5));
			p4 = GetPos2DD(Vector2D(l2, ly*0.5));
			mbs->DrawQuad(ToP3D(p4),ToP3D(p3),ToP3D(p2),ToP3D(p1));
		}
	}
};


