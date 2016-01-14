//#**************************************************************
//#
//# filename:             ANCFBeam2D.cpp
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						23. October 2007
//# description:          ANCF beam formulation according to Omar and Shabana
//#                       Modifications according to geometrically exact beam theory Reissner / Simo and VuQuoc
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
 
//#include "stdafx.h"

#include "ANCFBeam2D.h"
#include "elementDataAccess.h"
#include "Material.h"
#include "node.h"
#include "femathHelperFunctions.h"
#include "graphicsConstants.h"

//Element shares nodes with other elements, n1 and n2 are nodenumbers; element sets initial conditions for nodes
ANCFBeam2D::ANCFBeam2D(MBS* mbsi, const Vector& xc1, const Vector& xc2, const Vector& vc1, const Vector& vc2, 
											 int n1i, int n2i, double rhoi, double Emi, double nui,
											 const Vector3D& si, const Vector3D& coli, int beammodel):
Body2D(mbsi), massmatrix(), Hmatrix(), SV(), DS(), x1(), x2(), w1(), w2(),
T1(), T2(), T(), xg(), xgd(), e0()
{
	elasticforce_beam = beammodel;
	n1=n1i; n2=n2i; sos2=0;
	size = si;
	lx = si.X(); ly = si.Y(); lz = si.Z();

	mass = lx*ly*lz*rhoi;
	//lx=1;ly=1;lz=1;

	x_init = xc1.Append(xc2);
	xg = xc1.Append(xc2);
	Vector v_init = vc1.Append(vc2);

	nu = nui;
	Em = Emi;
	rho = rhoi;
	col = coli;

	double EA, EI, GAks, rhoA, rhoI;
	SetRectangularBeamParameters(EA, EI, GAks, rhoA, rhoI);
	SetBeamParameters(EA, EI, GAks, rhoA, rhoI);

	InitConstructor();
	BuildDSMatrices();
	Matrix Tinv = T;
	if (!Tinv.Invert2()) {UO() << "ERROR: T matrix not invertible!!!\n";}
	e0 = x_init;
	x_init = Tinv * x_init;
	v_init = Tinv * v_init;
	x_init = x_init.Append(v_init); //Velocity initial conditions can also be transformed by Tinv!!!

	SetInitialCurvatures(); //must be done after x_init!!!

	elementname = GetElementSpec();
};



void ANCFBeam2D::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{

	Body2D::GetElementData(edc);

	ElementData ed;

	ed.SetDouble(Em, "Youngs_modulus"); edc.Add(ed);
	ed.SetDouble(nu, "Poisson_ratio"); edc.Add(ed);
	ed.SetVector3D(lx, ly, lz, "Beam_dimensions"); edc.Add(ed);

	ed.SetInt(n1, "Node_number1", 1, GetMBS()->NNodes()); edc.Add(ed);
	ed.SetInt(n2, "Node_number2", 1, GetMBS()->NNodes()); edc.Add(ed);

	ed.SetVector2D(concentratedmass1, concentratedmass2, "Node_masses"); edc.Add(ed);

	ed.SetBool(elasticforce_beam!=0, "Elastic_line_model"); edc.Add(ed);
/*
	ed.SetDouble(beamEA, "Axial_stiffness_EA"); edc.Add(ed);
	ed.SetDouble(beamEI, "Bending_stiffness_EI"); edc.Add(ed);
	ed.SetDouble(beamGAks, "Shear_stiffness_GA"); edc.Add(ed);
	ed.SetDouble(beamRhoA, "Cross_section_inertia_RhoA"); edc.Add(ed);
	ed.SetDouble(beamRhoI, "Moment_of_inertia_RhoI"); edc.Add(ed);
*/
	Matrix Tinv = T;
	Vector xinit = x_init.SubVector(1, SOS());
	Vector vinit = x_init.SubVector(SOS()+1, 2*SOS());
	xinit = T*xinit;
	vinit = T*vinit;

	ed.SetVector2D(xinit( 1), xinit( 2), "Node1_r"); edc.Add(ed);
	ed.SetVector2D(xinit( 3), xinit( 4), "Node1_rx"); edc.Add(ed);
	ed.SetVector2D(xinit( 5), xinit( 6), "Node1_ry"); edc.Add(ed);

	ed.SetVector2D(xinit( 1+6), xinit( 2+6), "Node2_r"); edc.Add(ed);
	ed.SetVector2D(xinit( 3+6), xinit( 4+6), "Node2_rx"); edc.Add(ed);
	ed.SetVector2D(xinit( 5+6), xinit( 6+6), "Node2_ry"); edc.Add(ed);

	ed.SetVector2D(vinit( 1), vinit( 2), "Node1_v"); edc.Add(ed);
	ed.SetVector2D(vinit( 3), vinit( 4), "Node1_vx"); edc.Add(ed);
	ed.SetVector2D(vinit( 5), vinit( 6), "Node1_vy"); edc.Add(ed);

	ed.SetVector2D(vinit( 1+6), vinit( 2+6), "Node2_v"); edc.Add(ed);
	ed.SetVector2D(vinit( 3+6), vinit( 4+6), "Node2_vx"); edc.Add(ed);
	ed.SetVector2D(vinit( 5+6), vinit( 6+6), "Node2_vy"); edc.Add(ed);

}

int ANCFBeam2D::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Body2D::SetElementData(edc);

	InitConstructor();
	double EA, EI, GAks, rhoA, rhoI;
	SetRectangularBeamParameters(EA, EI, GAks, rhoA, rhoI);
	SetBeamParameters(EA, EI, GAks, rhoA, rhoI);

	sos2 = 0;
	Em = 0;
	elasticforce_beam = 0;

	GetElemDataDouble(GetMBS(), edc, "Youngs_modulus", Em, 1);
	GetElemDataDouble(GetMBS(), edc, "Poisson_ratio", nu, 1);

	GetElemDataVector3D(GetMBS(), edc, "Beam_dimensions", lx, ly, lz, 1); 

	GetElemDataInt(GetMBS(), edc, "Node_number1", n1, 1);
	GetElemDataInt(GetMBS(), edc, "Node_number2", n2, 1);


	if (n1 < 1 || n1 > GetMBS()->NNodes() || n2 < 1 || n2 > GetMBS()->NNodes())
	{
		GetMBS()->EDCError("Invalid node numbers in ANCF beam element. Element has not been created!");
		return 0;
	}

	GetElemDataVector2D(GetMBS(), edc, "Node_masses", concentratedmass1, concentratedmass2, 0);

	GetElemDataBool(GetMBS(), edc, "Elastic_line_model", elasticforce_beam, 0);
	if (elasticforce_beam) elasticforce_beam = 1;

	/*
	GetElemDataDouble(GetMBS(), edc, "Axial_stiffness_EA", beamEA, 1);
	GetElemDataDouble(GetMBS(), edc, "Bending_stiffness_EI", beamEI, 1);
	GetElemDataDouble(GetMBS(), edc, "Shear_stiffness_GA", beamGAks, 1);
	GetElemDataDouble(GetMBS(), edc, "Cross_section_inertia_RhoA", beamRhoA, 1);
	GetElemDataDouble(GetMBS(), edc, "Moment_of_inertia_RhoI", beamRhoI, 1);
	*/

	Vector xinit(SOS());
	Vector vinit(SOS());

	GetElemDataVector2D(GetMBS(), edc, "Node1_r" , xinit( 1), xinit( 2), 1);
	GetElemDataVector2D(GetMBS(), edc, "Node1_rx", xinit( 3), xinit( 4), 1);
	GetElemDataVector2D(GetMBS(), edc, "Node1_ry", xinit( 5), xinit( 6), 1);

	GetElemDataVector2D(GetMBS(), edc, "Node2_r" , xinit( 1+6), xinit( 2+6), 1);
	GetElemDataVector2D(GetMBS(), edc, "Node2_rx", xinit( 3+6), xinit( 4+6), 1);
	GetElemDataVector2D(GetMBS(), edc, "Node2_ry", xinit( 5+6), xinit( 6+6), 1);

	GetElemDataVector2D(GetMBS(), edc, "Node1_v" , vinit( 1), vinit( 2), 1);
	GetElemDataVector2D(GetMBS(), edc, "Node1_vx", vinit( 3), vinit( 4), 1);
	GetElemDataVector2D(GetMBS(), edc, "Node1_vy", vinit( 5), vinit( 6), 1);

	GetElemDataVector2D(GetMBS(), edc, "Node2_v" , vinit( 1+6), vinit( 2+6), 1);
	GetElemDataVector2D(GetMBS(), edc, "Node2_vx", vinit( 3+6), vinit( 4+6), 1);
	GetElemDataVector2D(GetMBS(), edc, "Node2_vy", vinit( 5+6), vinit( 6+6), 1);

	size = Vector3D(lx,ly,lz);
	mass = GetMaterial().BeamRhoA()*lx;


	BuildDSMatrices();
	Matrix Tinv = T;
	if (!Tinv.Invert2()) {UO() << "ERROR: T matrix not invertible!!!\n";}
	e0 = xinit;
	xinit = Tinv * xinit;
	vinit = Tinv * vinit;
	x_init = xinit.Append(vinit); //Velocity initial conditions can also be transformed by Tinv!!!

	SetInitialCurvatures(); //must be done after x_init!!!	

	return rv;
}

int ANCFBeam2D::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = Element::CheckConsistency(errorstr);
	if (rv) return rv;


	if (n1 && n2)
	{
		if (GetMBS()->GetNode(n1).SOS() != NodeSize()) 
		{
			GetMBS()->GetNode(n1).SOS() = NodeSize();
			errorstr += mystr("ANCFbeam ")+ mystr(GetOwnNum()) + mystr(": Node does not have ")+mystr(NodeSize())+mystr(" coordinates; the size has been corrected");
		}
		if (GetMBS()->GetNode(n2).SOS() != NodeSize()) 
		{
			GetMBS()->GetNode(n2).SOS() = NodeSize();
			errorstr += mystr("ANCFbeam ")+ mystr(GetOwnNum()) + mystr(": Node does not have ")+mystr(NodeSize())+mystr(" coordinates; the size has been corrected");
		}
	}

	return rv;
}



void ANCFBeam2D::LinkToElements()
{
	if (SOSowned() == 0)
	{
		const Node& node1 = GetMBS()->GetNode(n1);
		const Node& node2 = GetMBS()->GetNode(n2);
		for (int i=1; i <= node1.SOS(); i++)
		{
			AddLTG(node1.Get(i));
		}
		for (int i=1; i <= node2.SOS(); i++)
		{
			AddLTG(node2.Get(i));
		}
		for (int i=1; i <= node1.SOS(); i++)
		{
			AddLTG(node1.Get(i+node1.SOS()));
		}
		for (int i=1; i <= node2.SOS(); i++)
		{
			AddLTG(node2.Get(i+node2.SOS()));
		}
	}
}

void ANCFBeam2D::BuildDSMatrices() 
{

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//for original ANCF:
	orderx = 9; //max 9x5, sonst array grad zu klein!!!!
	orderyz = 5; //9x5 und 7x3 ergibt fast das gleiche!!!

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//for Arend Schwab ANCF
	//for all integration points:
	orderWl = 5; //5; 5 and 6 coincide with 8 up to 10 digits, 7=8
	orderWs = 2; //2; more than 2 gives locking, 1 is rather inaccurate
	orderWsHR=4; //4; for Hellinger Reissner, 4 is exact
	orderWb = 3; //3; 3=4=6, 2 gives shear locking

	factstiffWl = 1; //reduce stiffness of thickness modes, standard==1

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	GetIntegrationRule(x1,w1,orderx);
	GetIntegrationRule(x2,w2,orderyz);
	Matrix3D jac, jacinv;
	int kx1 = x2.Length();

	//UO() << "Initialize Element\n";

	for (int i1=1; i1<=x1.GetLen(); i1++)
	{
		for (int i2=1; i2<=x2.GetLen(); i2++)
		{
			int ind = (i1-1)*kx1+(i2-1);
			//UO() << "ind=" << ind << "\n";

			Vector2D p(x1(i1),x2(i2));

			GetDSMatrix0(p,DS);

			GetJacobi(jac,p,DS,x_init);
			jacdet[ind] = jac.Det();
			jac.GetInverse(jacinv);
			jacinv = jacinv.GetTp();

			grad[ind].Init();
			grad[ind].SetSize(Dim(),NS());
			Mult(jacinv, DS, grad[ind]);
		}
	}
	//Build T-matrices:
	T1.Init();
	T2.Init();
	T1.SetSize(NodeSize(),NodeSize());
	T1.SetAll(0);
	T2=T1;
	if (SOS()!=8)
	{
		Vector2D p1(-1,0);
		Vector2D p2( 1,0);

		GetDSMatrix0(p1,DS);
		GetJacobi(jac1,p1,DS,x_init);
		jac1 = jac1.GetTp();
		GetDSMatrix0(p2,DS);
		GetJacobi(jac2,p2,DS,x_init);
		jac2 = jac2.GetTp();

		//UO() << "jac1A=\n" << jac1 << "\n";
		//UO() << "jac2A=\n" << jac2 << "\n";

		for (int i=1; i <= Dim(); i++)
		{
			for (int j=1; j <= Dim(); j++)
			{
				jac1(i,j) = x_init((i-1)*Dim()+j+Dim());
				jac2(i,j) = x_init((i-1)*Dim()+j+NodeSize()+Dim());
			}
		}

		T1(1,1) = 1; T1(2,2) = 1;
		T2 = T1;
		for (int i=1; i <= Dim(); i++)
		{
			for (int j=1; j <= Dim(); j++)
			{
				T1(i*Dim()+1,j*Dim()+1) = jac1(i,j);
				T1(i*Dim()+2,j*Dim()+2) = jac1(i,j);
				T2(i*Dim()+1,j*Dim()+1) = jac2(i,j);
				T2(i*Dim()+2,j*Dim()+2) = jac2(i,j);
			}
		}
	}
	else
	{
		T1.SetDiagMatrix(1);
		T2.SetDiagMatrix(1);
	}


	T.Init();
	T.SetSize(SOS(),SOS());
	T.SetAll(0);
	T.SetSubmatrix(T1,1,1);
	T.SetSubmatrix(T2,1+T1.Getrows(),1+T1.Getcols());
};

int writtenbp = 0;
void ANCFBeam2D::SetBeamParameters(double EA, double EI, double GAks, double rhoA, double rhoI)
{
	if (!GetMaterialNum())
	{
		Material m;
		AddMaterial(m); //materialnum is internally set!
	}
	GetMaterial().BeamEA() = EA;
	GetMaterial().BeamEIz() = EI;
	GetMaterial().BeamGAky() = GAks;
	GetMaterial().BeamRhoA() = rhoA;
	GetMaterial().BeamRhoIz() = rhoI;

	if (!writtenbp)
	{
		UO() << "**** Beam parameters rhoA and rhoI not yet implemented in ANCFBeam2D!\n";
		writtenbp = 1;
	}

	SetInitialCurvatures();
}

void ANCFBeam2D::SetRectangularBeamParameters(double& beamEA, double& beamEI, double& beamGAks, double& beamRhoA, double& beamRhoI) //set beam parameters for rectangular cross-section, Schwab-Meijaard model
{

	double A = ly*lz;
	double Gm = Em / (2.*(1+nu));
	double GA = Gm*A;
	double ks = 10*(1+nu)/(12.+11.*nu); //shear correction factors (square): 10*(1+nu)/(12.+11.*nu) 

	beamEA = Em * A;
	beamEI = Em * lz*Cub(ly)/12.;
	beamGAks = GA*ks;
	beamRhoA = rho*A;
	beamRhoI = rho*lz*Cub(ly)/12.;
}



void ANCFBeam2D::ApplyT(Vector& v) const
{
	if (SOS() == 8) return; //T does not work in linear element

	static Vector x;
	x = v;
	for (int i=1; i<=Dim(); i++)
	{
		for (int j=1; j<=Dim(); j++)
		{
			v((i-1)*Dim()+j+Dim()) =            jac1(i,1)*x(0+j+Dim())            + jac1(i,2)*x(Dim()+j+Dim());
			v((i-1)*Dim()+j+Dim()+NodeSize()) = jac2(i,1)*x(0+j+Dim()+NodeSize()) + jac2(i,2)*x(Dim()+j+Dim()+NodeSize());
		}
	}
	//global_uo << "T*v2=" << v << "\n";
}

void ANCFBeam2D::ApplyTD(Vector& v) const
{
	if (SOS() == 8) return; //T does not work in linear element

	static Vector x;
	x = v;
	for (int i=1; i<=Dim(); i++)
	{
		for (int j=1; j<=Dim(); j++)
		{
			v((i-1)*Dim()+j+Dim())            = jac1(i,1)*x(0+j+Dim()) + jac1(i,2)*x(Dim()+j+Dim());
			v((i-1)*Dim()+j+Dim()+NodeSize()) = jac2(i,1)*x(0+j+Dim()+NodeSize()) + jac2(i,2)*x(Dim()+j+Dim()+NodeSize());
		}
	}
}

void ANCFBeam2D::ApplyTtp(Vector& v) const
{
	if (SOS() == 8) return; //T does not work in linear element

	//v=T.GetTp()*v; return;
	static Vector x;
	x = v;
	for (int i=1; i<=Dim(); i++)
	{
		for (int j=1; j<=Dim(); j++)
		{
			v((i-1)*Dim()+j+Dim())            = jac1(1,i)*x(0+j+Dim())            + jac1(2,i)*x(Dim()+j+Dim());
			v((i-1)*Dim()+j+Dim()+NodeSize()) = jac2(1,i)*x(0+j+Dim()+NodeSize()) + jac2(2,i)*x(Dim()+j+Dim()+NodeSize());
		}
	}
}

void ANCFBeam2D::ApplyTtp(Matrix& m) const
{
	if (SOS() == 8) return; //T does not work in linear element

	//v=T.GetTp()*v; return;
	static Vector x;
	int sns = SOS();
	x.SetLen(sns);

	for (int k = 1; k <= m.Getcols(); k++)
	{
		for (int i = 1+Dim(); i <= sns; i++)
		{
			x(i) = m(i,k);
		}
		for (int i=1; i<=Dim(); i++)
		{
			for (int j=1; j<=Dim(); j++)
			{
				m((i-1)*Dim()+j+Dim(),k) = jac1(1,i)*x(0+j+Dim()) + jac1(2,i)*x(Dim()+j+Dim());
				m((i-1)*Dim()+j+Dim()+NodeSize(),k) = jac2(1,i)*x(0+j+Dim()+NodeSize()) + jac2(2,i)*x(Dim()+j+Dim()+NodeSize());
			}
		}
	}
}

void ANCFBeam2D::GetS0(Vector& sf, const Vector2D& ploc) const
{
	double xb = ploc.X();
	double yb = ploc.Y();
	sf.SetLen(6);
	double xb2 = xb*xb;
	double xb3 = xb2*xb;
	sf(1) = 1.0/2.0-3.0/4.0*xb+xb3/4.0;
	sf(2) = lx*(1.0-xb-xb2+xb3)/8.0;
	sf(3) = -yb*ly*(-1.0+xb)/4.0;
	sf(4) = 1.0/2.0+3.0/4.0*xb-xb3/4.0;
	sf(5) = lx*(1.0+xb)*(1.0+xb)*(-1.0+xb)/8.0;
	sf(6) = (1.0+xb)*yb*ly/4.0;
}

void ANCFBeam2D::GetDSMatrix0(const Vector2D& ploc, Matrix& sf) const
{
	double xb = ploc.X();
	double yb = ploc.Y();
	sf.SetSize(2,6);

	sf(1,1) = -3.0/4.0+3.0/4.0*xb*xb;
	sf(1,2) = lx*(-1.0-2.0*xb+3.0*xb*xb)/8.0;
	sf(1,3) = -yb*ly/4.0;
	sf(1,4) = 3.0/4.0-3.0/4.0*xb*xb;
	sf(1,5) = lx*(1.0+xb)*(-1.0+xb)/4.0+lx*(1.0+xb)*(1.0+xb)/8.0;
	sf(1,6) = yb*ly/4.0;
	sf(2,1) = 0.0;
	sf(2,2) = 0.0;
	sf(2,3) = -ly*(-1.0+xb)/4.0;
	sf(2,4) = 0.0;
	sf(2,5) = 0.0;
	sf(2,6) = (1.0+xb)*ly/4.0;

}

void ANCFBeam2D::GetS0x(Vector& sfx, const Vector2D& ploc) const
{
	//at y=0 !!!!!!!!!!!!!!!!!!!!!
	double xb = ploc.X();
	sfx.SetLen(NS());
	double f = 2./lx;

	sfx(1) = f*(-3.0/4.0+3.0/4.0*xb*xb);
	sfx(2) = f*(lx*(-1.0-2.0*xb+3.0*xb*xb)/8.0);
	sfx(3) = 0.;
	sfx(4) = f*(3.0/4.0-3.0/4.0*xb*xb);
	sfx(5) = f*(lx*(1.0+xb)*(-1.0+xb)/4.0+lx*(1.0+xb)*(1.0+xb)/8.0);
	sfx(6) = 0.;

}

void ANCFBeam2D::GetS0xlin(Vector& sfx, const Vector2D& ploc) const
{
	//leads to small error in large deformation problems!!!
	//at y=0 !!!!!!!!!!!!!!!!!!!!!
	double xb = ploc.X();
	sfx.SetLen(NS());
	double f = 2./lx;

	sfx(1) = 0;
	sfx(2) = f*(lx*(1.0-1.0*xb)/4.0);
	sfx(3) = 0.;
	sfx(4) = 0;
	sfx(5) = f*(lx*(1.0+1.0*xb)/4.0);
	sfx(6) = 0.;

}

void ANCFBeam2D::GetS0xlin2(Vector& sfx, const Vector2D& ploc) const
{
	//leads to locking
	//at y=0 !!!!!!!!!!!!!!!!!!!!!
	double xb = ploc.X();
	sfx.SetLen(NS());
	double f = 2./lx;

	double k = 0.5;

	sfx(1) = f*(-k);
	sfx(2) = f*(-0.5*k*lx*xb);
	sfx(3) = 0.;
	sfx(4) = f*(k);
	sfx(5) = f*(0.5*k*lx*xb);
	sfx(6) = 0.;
}

void ANCFBeam2D::GetS0y(Vector& sfx, const Vector2D& ploc) const
{
	double xb = ploc.X();
	sfx.SetLen(NS());
	double f = 2./ly;
	sfx(1) = 0.0;
	sfx(2) = 0.0;
	sfx(3) = -f*ly*(-1.0+xb)/4.0;
	sfx(4) = 0.0;
	sfx(5) = 0.0;
	sfx(6) = f*(1.0+xb)*ly/4.0;
}

void ANCFBeam2D::GetS0yx(Vector& sfx, const Vector2D& ploc) const
{
	double xb = ploc.X();
	sfx.SetLen(NS());
	double f = 4./(ly*lx);
	sfx(1) = 0.0;
	sfx(2) = 0.0;
	sfx(3) = -f*ly/4.0;
	sfx(4) = 0.0;
	sfx(5) = 0.0;
	sfx(6) = f*ly/4.0;
}

void ANCFBeam2D::GetS0xx(Vector& sfxx, const Vector2D& ploc) const
{
	double xb = ploc.X();
	sfxx.SetLen(NS());
	double f = 4./Sqr(lx);
	sfxx(1) = f*3.0/2.0*xb;
	sfxx(2) = f*(-1.0+3.0*xb)*lx/4.0;
	sfxx(3) = 0.0;
	sfxx(4) = -f*3.0/2.0*xb;
	sfxx(5) = f*(1.0+3.0*xb)*lx/4.0;
	sfxx(6) = 0.0;

}



void ANCFBeam2D::GetRot(Matrix3D& rot, const Vector& xg) const
{
	Vector2D r1(xg(1+NodeSize())-xg(1),xg(2+NodeSize())-xg(2)); //r_B - r_A

	r1.Normalize();
	Vector2D r2(-r1.Y(),r1.X());

	rot.Set22(r1.X(),r2.X(),r1.Y(),r2.Y());
}

Vector2D ANCFBeam2D::GetDisplacement2D(const Vector2D& p_loc, const Vector& xg, int flagD) const
{
	Vector2D p0 = p_loc;
	p0.Scale(0.5*lx,0.5*ly);
	
	ConstVector<ANCFBeam2D_nshapes> SV;
	GetS0(SV, p0);

	Vector2D p(0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += SV(j)*(xg((j-1)*Dim()+i) - e0((j-1)*Dim()+i));
		}
	}
	return p;
}

Vector2D ANCFBeam2D::GetPos2D(const Vector2D& p_loc) const
{
	static Vector xg;
	xg.SetLen(SOS());
	GetCoordinates(xg);
	ApplyT(xg);
	Vector2D p0=p_loc;
	p0.Scale(0.5*lx,0.5*ly);
	static Vector SV;
	GetS0(SV, p0);
	Vector2D p(0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += SV(j)*xg((j-1)*Dim()+i);
		}
	}
	return p;
};

Vector2D ANCFBeam2D::GetVel2D(const Vector2D& p_loc) const
{
	static Vector xg;
	xg.SetLen(SOS());
	GetCoordinatesP(xg);
	ApplyT(xg);
	Vector2D p0=p_loc;
	p0.Scale(0.5*lx,0.5*ly);
	static Vector SV;
	GetS0(SV, p0);
	Vector2D p(0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += SV(j)*xg((j-1)*Dim()+i);
		}
	}
	return p;
};

Vector2D ANCFBeam2D::GetPos2DD(const Vector2D& p_loc) const
{
	static Vector xgd;
	xgd.SetLen(SOS());
	GetDrawCoordinates(xgd);
	xgd = T*xgd;
	Vector2D p0=p_loc;
	p0.Scale(0.5*lx,0.5*ly);
	static Vector SV;
	GetS0(SV, p0);
	Vector2D p(0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += SV(j)*xgd((j-1)*Dim()+i);
		}
	}
	return p;
};

//in reference element coordinates (-1..1)
Vector2D ANCFBeam2D::GetPos2D0D(const Vector2D& p_loc) const
{
	static Vector xgd;
	xgd.SetLen(SOS());
	GetDrawCoordinates(xgd);
	//xgd = T*xgd;
	ApplyTD(xgd);
	static Vector SV;
	GetS0(SV, p_loc);
	Vector2D p(0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += SV(j)*xgd((j-1)*Dim()+i);
		}
	}
	if (GetMBS()->GetDOption(105) != 1.)
	{
		Vector2D p0 = GetPos2D0Dinit(p_loc);
		return (p - p0)*GetMBS()->GetDOption(105) + p0;
	}

	return p;
};

//in reference element coordinates (-1..1)
Vector2D ANCFBeam2D::GetPos2D0Dinit(const Vector2D& p_loc) const
{
	static Vector xgd;
	xgd.SetLen(SOS());
	xgd = x_init;
	ApplyTD(xgd);
	static Vector SV;
	GetS0(SV, p_loc);
	Vector2D p(0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += SV(j)*xgd((j-1)*Dim()+i);
		}
	}
	return p;
};

Vector2D ANCFBeam2D::GetVel2DD(const Vector2D& p_loc) const
{
	static Vector xgd;
	xgd.SetLen(SOS());
	GetDrawCoordinatesP(xgd);
	xgd = T*xgd;
	Vector2D p0=p_loc;
	p0.Scale(0.5*lx,0.5*ly);
	static Vector SV;
	GetS0(SV, p0);
	Vector2D p(0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += SV(j)*xgd((j-1)*Dim()+i);
		}
	}
	return p;
};


Vector2D ANCFBeam2D::GetPosx2D(const Vector2D& p_loc) const
{
	static Vector xg;
	xg.SetLen(SOS());
	GetCoordinates(xg);
	ApplyT(xg);
	Vector2D p0=p_loc;
	p0.Scale(0.5*lx,0.5*ly);
	static Vector SV;
	GetS0x(SV, p0);
	Vector2D p(0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += SV(j)*xg((j-1)*Dim()+i);
		}
	}
	return p;
};

Vector2D ANCFBeam2D::GetPosxx2D(const Vector2D& p_loc) const
{
	static Vector xg;
	xg.SetLen(SOS());
	GetCoordinates(xg);
	ApplyT(xg);
	Vector2D p0=p_loc;
	p0.Scale(0.5*lx,0.5*ly);
	static Vector SV;
	GetS0xx(SV, p0);
	Vector2D p(0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += SV(j)*xg((j-1)*Dim()+i);
		}
	}
	return p;
};

void ANCFBeam2D::GetH(Matrix& H) 
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

		DS.SetSize(dim,ns);
		SV.SetLen(ns);

		GetIntegrationRule(x1,w1,3); //3x1x1 !!!!!
		GetIntegrationRule(x2,w2,1);

		for (int i1=1; i1<=x1.GetLen(); i1++)
		{
			for (int i2=1; i2<=x2.GetLen(); i2++)
			{
				Vector2D p(x1(i1),x2(i2));
				GetS0(SV,p);

				GetDSMatrix0(p,DS);
				GetJacobi(jac,p,DS, e0);
				double jacdet = jac.Det();
				double fact = fabs (jacdet) * w1(i1)*w2(i2) * lz;

				for (int i=0; i<ns; i++)
				{
					for (int j=1; j<=dim; j++)
					{
						H(i*dim+j,j)+=fact*SV(i+1);
					}
				}
			}
		}
		ApplyTtp(H);
		Hmatrix = H;
	}
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double ANCFBeam2D::GetKappa(const double& x, const Vector& xg) const
{
	Vector2D rx = GetPosx2D(Vector2D(x*0.5*lx,0.)); 
	Vector2D rxx = GetPosxx2D(Vector2D(x*0.5*lx,0.)); 
	double rxn = rx.Norm();

	if (rxn == 0) return 0;
	else
	{
		return (rx.Cross(rxx))/Sqr(rxn);
	}
}


//x:-1...+1
void ANCFBeam2D::GetDeltaKappa(const double& x, const Vector& xg, Vector& dkappa, double& kappa) const
{

	int dim = Dim();
	int ns = NS();

	static Vector SVx;
	static Vector SVxx;
	SVx.SetLen(NS());
	SVxx.SetLen(NS());

	Vector2D ploc(x,0.);
	GetS0x(SVx,ploc);
	GetS0xx(SVxx,ploc);

	Vector2D rx(0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= ns; j++)
		{
			rx(i) += SVx(j)*xg((j-1)*dim+i);
		}
	}

	Vector2D rxx(0.,0.);
	for (int i = 1; i <= dim; i++)
	{
		for (int j = 1; j <= ns; j++)
		{
			rxx(i) += SVxx(j)*xg((j-1)*dim+i);
		}
	}

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

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ANCFBeam2D::EvalM(Matrix& m, double t) 
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

		SV.SetLen(ns);
		Matrix3D jac;

		Matrix HL(ns*dim,dim);
		DS.SetSize(dim,ns);

		//Vector x1,x2,w1,w2;

		GetIntegrationRule(x1,w1,6); //optimal: 6x3x3, 6x1x1 geht auch!!!!
		GetIntegrationRule(x2,w2,3);

		for (int i1=1; i1<=x1.GetLen(); i1++)
		{
			for (int i2=1; i2<=x2.GetLen(); i2++)
			{
				Vector2D p(x1(i1),x2(i2));
				GetS0(SV,p);

				for (int i=0; i<ns; i++)
				{
					for (int j=1; j<=dim; j++)
					{
						HL(i*dim+j,j)=SV(i+1);
					}
				}

				GetDSMatrix0(p,DS);
				//UO() << "DS=\n" << DS << "\n";
				GetJacobi(jac,p,DS,e0);
				//UO() << "jac=\n" << jac << "\n";

				double jacdet = jac.Det();
				//UO() << "jacdet=\n" << jacdet << ", rho=" << rho << "\n";
				m += fabs (jacdet) * rho * lz * w1(i1)*w2(i2) * (HL*HL.GetTp());
			}
		}
		m(1,1) += concentratedmass1;
		m(2,2) += concentratedmass1;
		m(1+NodeSize(),1+NodeSize()) += concentratedmass2;
		m(2+NodeSize(),2+NodeSize()) += concentratedmass2;

		massmatrix = T.GetTp()*m*T;
		m = massmatrix;
	}
	//Matrix minv = m;
	//UO() << "Mass matrix invertable=" << minv.Invert2() << "\n";
	//UO() << "m=" << m << "\n";
};


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//compute curvature, stretch, shear and torsion for initial values in order to make precurved beams possible
void ANCFBeam2D::SetInitialCurvatures()
{
	if (elasticforce_beam != 1) return;

	int sos = SOS();
	//SetComputeCoordinates();
	//--> here: set xg = x_init;
	xg = x_init; 

	int ns = NS();
	int dim = Dim();

	ApplyT(xg);


	Vector Sx; Sx.SetLen(NS()); Sx.SetAll(0);  //rx
	Vector Sy; Sy.SetLen(NS()); Sy.SetAll(0);   //ry
	Vector Sxx; Sxx.SetLen(NS()); Sxx.SetAll(0); //rxx
	Vector Syx; Syx.SetLen(NS()); Syx.SetAll(0); //ryx


	Vector2D rx, ry, rxx, ryx;

	//needs to be computed for bending and shear!!!
}

void ANCFBeam2D::GetDMatrix(Matrix3D& D, double nu, double em) const
{
	//double nu = 1./2.*la/(la+mu);
	//double em = mu*(3.*la+2.*mu)/(la+mu);

	//for plane stress:
	D.SetSize(3,3);
	double nu2 = nu;

	double f = em/(1.-Sqr(nu2));

	D(1,1)=f;
	D(1,2)=nu2*f;
	D(1,3)=0;
	
	D(2,1)=nu2*f;
	D(2,2)=f;
	D(2,3)=0;
	
	D(3,1)=0;
	D(3,2)=0;

	//double ks = 10.*(1.+nu)/(12.+11.*nu); //for rectangular cross-section only!
	D(3,3)=.5*f*(1.-nu2);//*ks;

	//D(3,3) = Em/(2.*(1.+nu));
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ANCFBeam2D::EvalF2(Vector& f, double t) 
{
	Body2D::EvalF2(f,t);
	TMStartTimer(22);

	int sos = SOS();
	SetComputeCoordinates();
	static Vector fadd;

	int ns = NS();
	int dim = Dim();
	//double u;

	//0 is Original fully parametrized ANCF, continuum mechanics formulation
	//1 Geometrically exact approach
	//2 Geometrically exact approach + locking compensation
	//3 is Original fully parametrized ANCF + locking compensation, continuum mechanics formulation
	//4 is Original fully parametrized ANCF + locking compensation, continuum mechanics formulation + hyperelastic material

	if (elasticforce_beam == 1 || elasticforce_beam == 2)
	{
		//this approach uses the beam parameters!
		
		ApplyT(xg); // u = T*u

		temp.SetLen(SOS());
		fadd.SetLen(SOS());
		fadd.SetAll(0);


		static Vector Sx; Sx.SetLen(NS()); Sx.SetAll(0);  //rx
		static Vector Sy; Sy.SetLen(NS()); Sy.SetAll(0);   //ry
		static Vector Syx; Syx.SetLen(NS()); Syx.SetAll(0); //ryx
		static Vector Sxlin; Sxlin.SetLen(NS()); Sxlin.SetAll(0);  //rx
		Vector2D rx, ry, ryx, rxlin;

		double intK = 0;
		double intthetax = 0;
		double intGam2W1 = 0;
		double intGam2W2 = 0;
		Vector2D intrxW1 = 0;
		Vector2D intrxW2 = 0;

		GetIntegrationRule(x1,w1,orderx); //x_i ... integration points, w_i ... integration weights
		for (int i1=1; i1<=x1.GetLen(); i1++)
		{
			double x = x1(i1);

			Vector2D p0(x,0.);
			GetS0x(Sx, p0);
			GetS0y(Sy, p0);
			GetS0yx(Syx, p0);

			rx = 0;
			ry = 0;
			ryx = 0;
			for (int i = 1; i <= dim; i++)
				for (int j = 1; j <= ns; j++)
					rx(i) += Sx(j)*xg((j-1)*dim+i);

			for (int i = 1; i <= dim; i++)
				for (int j = 1; j <= ns; j++)
					ry(i) += Sy(j)*xg((j-1)*dim+i);

			for (int i = 1; i <= dim; i++)
				for (int j = 1; j <= ns; j++)
					ryx(i) += Syx(j)*xg((j-1)*dim+i);


			Vector2D t2 = ry;
			t2.Normalize();
			Vector2D t1(t2.Y(),-t2.X()); //90° CW !!!!!!!!!!!!!!!!!

			double ryn2 = Sqr(ry.X()) + Sqr(ry.Y());
			if (ryn2 == 0) ryn2 = 1; //this case should not be possible!
			double ryn = sqrt(ryn2);
			double ryn3 = Cub(ryn);


			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			//bending terms (theta_x):
			double g = ryn2;
			double ginv2 = 1./(g*g);

			double f = ry.Cross(ryx);
			double thetax = f/g;

			double df, dg;
			double fact =  GetMaterial().BeamEIz()  * thetax * 0.5*lx * w1(i1);


			if (elasticforce_beam == 2)
			{
				//locking compensation, only integrate bending part with high order
				double factkappa = GetMaterial().BeamEIz() * 0.5*lx * w1(i1);

				double kappa;
				GetDeltaKappa(p0.X(),xg,temp, kappa);
				temp *= kappa*factkappa;
				fadd += temp;

				intK += kappa * 0.5 * w1(i1);
				intthetax += thetax * 0.5 * w1(i1);
			}
			else
			{
				for (int i=1; i <= dim; i++)
				{
					for (int j=1; j <= ns; j++) 
					{
						//dr,y/de x r,yx + r,y x dr,yx/de: 

						switch (i) {
					case 1:
						{
							df =  Sy(j)*ryx.Y()-Syx(j)*ry.Y(); break; 
						}
					case 2:
						{
							df = -Sy(j)*ryx.X()+Syx(j)*ry.X(); break; 
						}
					default: ;
						}

						dg = 2.*(ry(i)*Sy(j));
						temp((j-1)*dim+i) = ginv2*(df*g-f*dg); //normed
					}
				}
				fadd.MultAdd(fact,temp); //not yet for precurved elements!!!
			}

			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			//axial deformation:

			//UO() << "t1=" << t1 << ", rx=" << rx << "\n";
			double Gam1 = t1*rx - 1.;
			fact = GetMaterial().BeamEA()  * Gam1 * 0.5*lx * w1(i1);
			//double dt1;
			//double hatf;
			Vector2D ryhat(ry.Y(),-ry.X());  //90° CW !!!!!!!!!!!!!!!!!
			Vector2D rxhat(rx.Y(),-rx.X());  //90° CW !!!!!!!!!!!!!!!!!

			for (int i=1; i <= dim; i++)
			{
				//rx*ryhat = -rxhat*ry; //!!!!!!!!!
				for (int j=1; j <= ns; j++)
				{
					//dt1 =  - ryhat(i)/ryn3*(ry(i)*Sy(j));
					temp((j-1)*dim+i) = t1(i)*Sx(j) - rxhat(i)/ryn*Sy(j) - (rx*ryhat)/ryn3*(ry(i)*Sy(j));
				}
			}
			fadd.MultAdd(fact,temp); //not yet for precurved elements!!!


			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			//shear deformation:

			double Gam2 = t2*rx;
			fact = GetMaterial().BeamGAky()  * Gam2 * 0.5*lx * w1(i1);
			//double dt2;

			intGam2W1 += 0.5*w1(i1) * 0.5*(1-x) * Gam2; //for Hellinger Reissner of shear strain
			intGam2W2 += 0.5*w1(i1) * 0.5*(1+x) * Gam2; //for Hellinger Reissner of shear strain
			intrxW1 += 0.5*w1(i1) * 0.5*(1-x) * rx; //for Hellinger Reissner of shear strain
			intrxW2 += 0.5*w1(i1) * 0.5*(1+x) * rx; //for Hellinger Reissner of shear strain

			if (elasticforce_beam != 2)
			{
				for (int i=1; i <= dim; i++)
				{
					for (int j=1; j <= ns; j++)
					{
						//dt2 =  - ry(i)/ryn3*(ry(i)*Sy(j));
						temp((j-1)*dim+i) = t2(i)*Sx(j) + 1./ryn*Sy(j)*rx(i) - (rx*ry)/ryn3*(ry(i)*Sy(j));
					}
				}

				fadd.MultAdd(fact,temp); //not yet for precurved elements!!!
			}
		}

		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//cross section deformation:

		//just integration equals sum of both nodes ==> Lobatto order 2
		for (int i1=1; i1<=2; i1++)
		{
			//TMStartTimer(20);

			double x = -1./sqrt(3.);
			if (i1 == 2) x = 1./sqrt(3.);

			//the following is equal to ry1=Vector(xg(5),xg(6)), ry2=Vector(xg(5+6),xg(6+6))
			Vector2D p0(x,0.);
			GetS0y(Sy, p0);
			ry = 0;
			for (int i = 1; i <= dim; i++)
				for (int j = 1; j <= ns; j++)
					ry(i) += Sy(j)*xg((j-1)*dim+i);


			double epsyy = 0.5*(ry*ry - 1.);
			double fact = Em * lz * ly * 0.5*lx * epsyy;
			//double fact = 1./Sqr(lx/ly)*Em * lz * 0.5*ly * epsyy; //just use cross section -> less stiff !!!

			for (int i=1; i <= dim; i++)
			{
				for (int j=1; j <= ns; j++)
				{
					temp((j-1)*dim+i) = ry(i)*Sy(j);
				}
			}
			fadd.MultAdd(fact,temp); //not yet for precurved elements!!!
		}

		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//locking compensation:

		if (elasticforce_beam == 2)
		{
			for (int kk=1; kk <= 2; kk++)
			{
				//UO() << "in\n";

				//does not work ==> apply Hellinger Reissner!!!
				//reduced integration: order-3 shear deformation + order-1 bending works for 1 element
				if (kk == 1)
					//GetIntegrationRuleLobatto(x1,w1,2); //shear deformation Lobatto2 works without HR
					GetIntegrationRule(x1,w1,5); //shear deformation ### use at least 5 in order to work without HR!
				else
					GetIntegrationRule(x1,w1,5); //use 5, does not influence eigenfrequ, but slightly influences large deformation!!!!

				for (int i1=1; i1<=x1.GetLen(); i1++)
				{
					double x = x1(i1);
					//TMStartTimer(20);
					Vector2D p0(x,0.);
					GetS0x(Sx, p0);
					GetS0y(Sy, p0);
					GetS0yx(Syx, p0);
					GetS0xlin(Sxlin, p0);

					rx = 0;
					ry = 0;
					ryx = 0;
					rxlin = 0; 
					for (int i = 1; i <= dim; i++)
						for (int j = 1; j <= ns; j++)
							rx(i) += Sx(j)*xg((j-1)*dim+i);

					for (int i = 1; i <= dim; i++)
						for (int j = 1; j <= ns; j++)
							rxlin(i) += Sxlin(j)*xg((j-1)*dim+i);

					for (int i = 1; i <= dim; i++)
						for (int j = 1; j <= ns; j++)
							ry(i) += Sy(j)*xg((j-1)*dim+i);

					for (int i = 1; i <= dim; i++)
						for (int j = 1; j <= ns; j++)
							ryx(i) += Syx(j)*xg((j-1)*dim+i);


					Vector2D t2 = ry;
					t2.Normalize();
					Vector2D t1(t2.Y(),-t2.X()); //90° CW !!!!!!!!!!!!!!!!!

					double ryn2 = Sqr(ry.X()) + Sqr(ry.Y());
					if (ryn2 == 0) ryn2 = 1; //this case should not be possible!
					double ryn = sqrt(ryn2);
					double ryn3 = Cub(ryn);


					//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					//shear deformation:

					if (kk == 1) 
					{
						double Gam2 = t2*rx;
						int HR = 2; //HR does not yet give better results?
						double factshear, fact1, fact2;
						Vector2D factshearrx;

						if (HR==0)
						{
						  factshear = GetMaterial().BeamGAky() * Gam2 * 0.5*lx * w1(i1); //without HellingerReissner
						}
						else if (HR == 1)
						{
							//++++++++++++++++++++++++++++++++++
							//does not work ....
							//HR: assume only linear rx: 
							Vector2D Ws(0.5*(1.-x), 0.5*(1.+x));
							Vector2D intWs1(intrxW1(1), intrxW2(1));
							Vector2D intWs2(intrxW1(2), intrxW2(2));

							//if (!GetMBS()->IsJacobianComputation() && GetMBS()->GetTime() > 0.999) 
							//	UO() << "rxlin" << rxlin << ",intrx=[" << intrxW1 << "," << intrxW2 << "]\n";

							Matrix3D Ss; //stiffness matrix for HR
							Ss.Set22(4.,-2.,-2.,4.);

							factshearrx(1) = (intWs1*(Ss*Ws));
							factshearrx(2) = (intWs2*(Ss*Ws));

							factshear = 1;
							fact1 = (t2*factshearrx) * GetMaterial().BeamGAky() * 0.5*lx * w1(i1);
						  fact2 = GetMaterial().BeamGAky() * Gam2 * 0.5*lx * w1(i1); //without HellingerReissner
						}
						else if (HR == 2)
						{
							//Vector2D rxlin2 = (0.5*(1.-x)*intrxW1 + 0.5*(1.+x)*intrxW2);
							//if (!GetMBS()->IsJacobianComputation() && GetMBS()->GetTime() > 0.9999) UO() << "rxlin" << rxlin << ",rx=" << rx << "\n";
							Gam2 = t2*rxlin; //works, but slightly wrong results ...
						  factshear = GetMaterial().BeamGAky() * Gam2 * 0.5*lx * w1(i1); //without HellingerReissner
						}

						for (int i=1; i <= dim; i++)
						{
							for (int j=1; j <= ns; j++)
							{
								//delta Gam2:
								if (HR == 0)
									temp((j-1)*dim+i) = t2(i)*Sx(j) + 1./ryn*Sy(j)*rx(i) - (rx*ry)/ryn3*(ry(i)*Sy(j));
								else if (HR == 1)
								{
									temp((j-1)*dim+i) = fact1*t2(i)*Sx(j) + fact1*(1./ryn*Sy(j)*rx(i) - (rx*ry)/ryn3*(ry(i)*Sy(j)));
								}
								else if (HR == 2)
									temp((j-1)*dim+i) = t2(i)*Sxlin(j) + 1./ryn*Sy(j)*rxlin(i) - (rxlin*ry)/ryn3*(ry(i)*Sy(j));
							}
						}
						fadd.MultAdd(factshear,temp); //not yet for precurved elements!!!
					}


					//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					//bending deformation:

					if (kk == 2)
					{
						double g = ryn2;
						double ginv2 = 1./(g*g);

						double f = ry.Cross(ryx);
						double thetax = f/g;
						double df, dg;

						//locking compensation, subtract low order bending part!!!

						double kappa;
						GetDeltaKappa(p0.X(),xg,temp, kappa);

						double factkappa = GetMaterial().BeamEIz() * (-intK) * 0.5*lx * w1(i1);
						temp *= factkappa;
						fadd += temp;

						//change of rotation of the cross section:
						double fact =  GetMaterial().BeamEIz() * (intthetax) * 0.5*lx * w1(i1);

						for (int i=1; i <= dim; i++)
						{
							for (int j=1; j <= ns; j++) 
							{
								//dr,y/de x r,yx + r,y x dr,yx/de: 

								switch (i) 
								{
								case 1:
									{
										df =  Sy(j)*ryx.Y()-Syx(j)*ry.Y(); break; 
									}
								case 2:
									{
										df = -Sy(j)*ryx.X()+Syx(j)*ry.X(); break; 
									}
								default: ;
								}

								dg = 2.*(ry(i)*Sy(j));
								temp((j-1)*dim+i) = ginv2*(df*g-f*dg); //normed
							}
						}
						fadd.MultAdd(fact,temp); //not yet for precurved elements!!!

					}
				}
			}
		}


	}
	else
	{
		//UO() << "standard ANCForig\n";
		static Vector u;
		u.SetLen(sos);
		int useu = 1;
		Matrix3D Id;
		Id.Set22(1.,0.,0.,1.);

		if (useu)
		{
			for (int i=1; i <= sos; i++)
				u(i) = xg(i) - x_init(i);
		}
		else
		{
			for (int i=1; i <= sos; i++)
				u(i) = xg(i);
		}
		ApplyT(u); // u = T*u
		//UO() << "u=" << u << "\n";

		//double la=Em * nu / ((1.+nu)*(1.-2.*nu));
		//double mu=Em / 2. / (1.+nu);

		Matrix3D strain, piola1, F;
		F.SetSize(2,2);

		temp.SetLen(SOS());
		fadd.SetLen(SOS());
		fadd.SetAll(0);


		int poissoncorrection = 0; //1==reduced integration of poisson part
		if (elasticforce_beam == 3) poissoncorrection = 1;
		//if (elasticforce_beam == 4) ==> hyperelastic material;

		for (int kk=1; kk <= 1+poissoncorrection; kk++)
		{
			Matrix3D Dm;

			double ks = 10.*(1.+nu)/(12.+11.*nu); //for rectangular cross-section only!

			if (poissoncorrection)
			{
				if (elasticforce_beam == 3)
				{
					Dm.SetAll(0);
					if (kk == 1)
					{
						//poisson ratio zero:
						Dm(1,1) = Em;
						Dm(2,2) = Em;
						Dm(3,3) = ks*Em/(2.*(1.+nu));

						//GetDMatrix(Dm,nu,Em);
						GetIntegrationRule(x1,w1,orderx);   //x_i ... integration points, w_i ... integration weights
						GetIntegrationRule(x2,w2,orderyz);
						if (SOS() == 8)
						{
							GetIntegrationRule(x1,w1,1);   //x_i ... integration points, w_i ... integration weights
							GetIntegrationRule(x2,w2,3);
							Dm(2,2) = 0;
						}
					}
					else if (kk == 2)
					{
						Dm(1,1) = Em/(1.-Sqr(nu))-Em;
						Dm(1,2) = nu*Em/(1.-Sqr(nu));
						Dm(2,1) = nu*Em/(1.-Sqr(nu));
						Dm(2,2) = Em/(1.-Sqr(nu))-Em;

						GetIntegrationRule(x1,w1,orderx);        //x_i ... integration points, w_i ... integration weights
						GetIntegrationRule(x2,w2,1);
						if (SOS() == 8)
						{
							GetIntegrationRule(x1,w1,3);   //x_i ... integration points, w_i ... integration weights
							GetIntegrationRule(x2,w2,1);
							Dm(2,2) += Em;
						}

					}
				}
				else if (elasticforce_beam == 3)
				{
					Dm.SetAll(0);
					if (kk == 1)
					{
						//poisson ratio zero:
						Dm(1,1) = Em;
						Dm(2,2) = Em;
						Dm(3,3) = ks*Em/(2.*(1.+nu));

						//GetDMatrix(Dm,nu,Em);
						GetIntegrationRule(x1,w1,orderx);   //x_i ... integration points, w_i ... integration weights
						GetIntegrationRule(x2,w2,orderyz);
						if (SOS() == 8)
						{
							GetIntegrationRule(x1,w1,1);   //x_i ... integration points, w_i ... integration weights
							GetIntegrationRule(x2,w2,3);
							Dm(2,2) = 0;
						}
					}
					else if (kk == 2)
					{
						Dm(1,1) = Em/(1.-Sqr(nu))-Em;
						Dm(1,2) = nu*Em/(1.-Sqr(nu));
						Dm(2,1) = nu*Em/(1.-Sqr(nu));
						Dm(2,2) = Em/(1.-Sqr(nu))-Em;

						GetIntegrationRule(x1,w1,orderx);        //x_i ... integration points, w_i ... integration weights
						GetIntegrationRule(x2,w2,1);
						if (SOS() == 8)
						{
							GetIntegrationRule(x1,w1,3);   //x_i ... integration points, w_i ... integration weights
							GetIntegrationRule(x2,w2,1);
							Dm(2,2) += Em;
						}

					}
				}
			}
			else
			{
			  GetDMatrix(Dm, nu, Em);
			}

			int kx1 = x2.Length();

			for (int i1=1; i1<=x1.GetLen(); i1++)
			{
				for (int i2=1; i2<=x2.GetLen(); i2++)
				{
					//TMStartTimer(20);
					int ind = (i1-1)*kx1+(i2-1);
					Vector2D p(x1(i1),x2(i2));

					// compute F 
					F.SetAll(0);
					/*
					int l;
					const Matrix& agrad = grad[ind]; //\bar S^D
					for (j = 1; j <= dim; j++) 
					{
						for (i = 1; i <= ns; i++)
						{
							l = (i-1)*dim+j;
							//u = xg(l) - x_init(l);
							for (k = 1; k <= dim; k++)
							{
								F(j,k) += grad[ind](k,i)*u(l);
							}
						}
						if (useu) F(j,j) += 1;
					}*/
					Gradu(p, u, F);

					GetDSMatrix0(p,DS);
					Matrix3D jac;

					GetJacobi(jac,p,DS,x_init);
					double jacdet = jac.Det();
					Matrix3D jacinv;
					jac.GetInverse(jacinv);
					jacinv = jacinv.GetTp();

					static Matrix grad;
					grad.SetSize(Dim(),NS());
					Mult(jacinv, DS, grad);

					if (useu)
					{
						F(1,1) += 1;
						F(2,2) += 1;
					}

					//jacinv = jacinv.GetTp();
					//G = G * jacinv; //also possible!!!
					//F = Matrix3D(1) + G;

					int linear=0;
					if (linear)
					{
						F -= Id;
						strain = 0.5*(F+F.GetTp());
						Vector3D s3(strain(1,1), strain(2,2), 2*strain(1,2));
						Vector3D stress = Dm*s3;
						piola1.Set22(stress.X(),stress.Z(),stress.Z(),stress.Y()); //linearized!
					}
					else if (elasticforce_beam == 4)
					{
						// Green-Lagrange strain tensor
						//strain = 0.5 * (F.GetTp() * F - I);

						//F.GetATA2(strain); //does not exist in 2D
						strain = 0.5 * (F.GetTp() * F);
						strain(1,1) -= 0.5; strain(2,2) -= 0.5;


						Matrix3D strain3D(strain(1,1),strain(1,2),0.,strain(2,1),0.,0.,0.,0.,0.);
						Matrix3D F3D(F(1,1),F(1,2),0.,F(2,1),F(2,2),0.,0.,0.,1.);

						Matrix3D stress3D;

						//replaced Matrix3D inel_strain by Vector inel_variables (by PG)
						Vector inel_variables(0);
						GetMaterial().ComputeStressFromStrain(strain3D, stress3D, inel_variables, F3D);

						//this formula does not work for plane stress/strain:
						//piola1 = F * ((2*mu) * strain + (la * strain.Trace())*Id );

						Matrix3D s2 = stress3D;
						//override stress22:
						stress3D(1,1) = Em*strain(1,1);
						stress3D(1,2) = Em*0.5*strain(1,2);
						stress3D(2,1) = Em*0.5*strain(2,1);
						stress3D(2,2) = Em*strain(2,2);
						
						stress3D = s2;
						stress3D(2,2) = Em*strain(2,2);
						double ks = 10.*(1.+nu)/(12.+11.*nu); //for rectangular cross-section only!
						stress3D(1,2) *= ks;
						stress3D(2,1) *= ks;

						piola1.Set22(stress3D(1,1),stress3D(1,2),stress3D(2,1),stress3D(2,2)); //linearized!

						piola1 = F * piola1;
					}
					else if (elasticforce_beam != 4)
					{
						// Green-Lagrange strain tensor
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
					fadd.MultAdd(fabs (jacdet) * lz * w1(i1)*w2(i2),temp);
				}
			}
		}
	}

	ApplyTtp(fadd);

	f -= fadd;
	TMStopTimer(22);

	if (GetMassDamping() != 0)
	{
		// +++++++ damping: +++++++
		for (int i = 1; i <= SOS(); i++)
			xg(i) = XGP(i);

		if (massmatrix.Getcols() == SOS())
		{
			Mult(massmatrix,xg,temp);
		}
		else
		{
			Matrix dmat;
			EvalM(dmat,t);
			Mult(dmat,xg,temp);
		}
		temp *= GetMassDamping();
		f -= temp;
	}
}; 

void ANCFBeam2D::GraduD(const Vector2D& ploc, const Vector& u, Matrix3D& gradu) const
{
	static Matrix DS;
	GetDSMatrix0(ploc,DS);

	Matrix3D jac, jacinv;
	GetJacobi(jac,ploc,DS,e0);

	jac.GetInverse(jacinv);
	jacinv = jacinv.GetTp();

	static Matrix grad;
	grad.SetSize(Dim(), NS());
	Mult(jacinv, DS, grad);

	gradu.SetSize(2,2);
	gradu.SetAll(0);

	int dim = Dim();
	int l;
	for (int j = 1; j <= dim; j++) 
	{
		for (int i = 1; i <= NS(); i++)
		{
			l = (i-1)*dim+j;
			for (int k = 1; k <= dim; k++)
			{
				gradu(j,k) += grad(k,i)*u(l);
			}
		}
	}
}

int linearstrainwarned = 0;

void ANCFBeam2D::GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables)
{
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_displacement,
		FieldVariableDescriptor::FVCI_y);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_stress,
		FieldVariableDescriptor::FVCI_y, true);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_total_strain,
		FieldVariableDescriptor::FVCI_y, true);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_stress_mises);
}

double ANCFBeam2D::GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector2D & local_position, bool flagD)
{
	// local_position is from -1 .. +1
	const Vector2D & ploc = local_position;
	Vector2D p(ploc.X()*lx*0.5, ploc.Y()*ly*0.5);

	Matrix3D F;
	Vector2D u(0.,0.);
	if (flagD == 0)
	{
		xg.SetLen(SOS());
		GetCoordinates(xg);

		xg = T*xg;
		//transform xgd to displacements, actually not necessary ...
		u = GetDisplacement2D(p, xg, 0);

		xg -= e0;
		Gradu(ploc,xg,F);
	}
	else
	{
		xgd.SetLen(SOS());
		GetDrawCoordinates(xgd);

		xgd = T*xgd;
		u = GetDisplacement2D(p, xgd, 1);
		//transform xgd to displacements, actually not necessary ...
		xgd -= e0;
		GraduD(ploc,xgd,F);
	}

	if (fvd.VariableType() == FieldVariableDescriptor::FVT_displacement)
		return fvd.GetComponent(u);

	if(
		fvd.VariableType() == FieldVariableDescriptor::FVT_stress ||
		fvd.VariableType() == FieldVariableDescriptor::FVT_total_strain ||
		fvd.VariableType() == FieldVariableDescriptor::FVT_stress_mises
		)
	{
		if (elasticforce_beam != 1 || 1)
		{
			F(1,1) += 1;
			F(2,2) += 1;
			Matrix3D Dm;
			GetDMatrix(Dm, nu, Em);

			Matrix3D stress;

			int linear = 0;
			if (!linear)
			{
				//Matrix3D strain = 0.5*(F.GetTp()*F);
				Matrix3D strain = 0.5*(F.GetTp()*F); //linear
				strain(1,1) -= 0.5; //linear
				strain(2,2) -= 0.5;
				Vector3D s3(strain(1,1), strain(2,2), 2.*strain(1,2));
				if (fvd.VariableType() != FieldVariableDescriptor::FVT_total_strain)
				{
					s3 = Dm*s3;
					stress(1,1) = s3(1);
					stress(2,2) = s3(2);
					stress(1,2) = s3(3);
					stress(2,1) = s3(3);
				}
				else
				{
					stress = strain;
				}

			}
			else
			{
				if (!linearstrainwarned)
				{
					linearstrainwarned = 1;
					UO() << "***************************\n";
					UO() << "Linear strain computation used in plane ANCF element\n";
					UO() << "***************************\n";
				}
				//linear:
				Matrix3D strain = 0.5*(F.GetTp()+F); //linear
				strain(1,1) -= 0.5*2; //linear
				strain(2,2) -= 0.5*2;
				Vector3D s3(strain(1,1), strain(2,2), 2.*strain(1,2));
				if (fvd.VariableType() != FieldVariableDescriptor::FVT_total_strain)
				{
					s3 = Dm*s3;
					stress(1,1) = s3(1);
					stress(2,2) = s3(2);
					stress(1,2) = s3(3);
					stress(2,1) = s3(3);
				}
				else
				{
					stress = strain;
				}
			}

			if(fvd.VariableType() == FieldVariableDescriptor::FVT_stress_mises)
			{
				return stress.Mises();
			}
			else
			{
				return fvd.GetComponent(stress);
			}
		}
		else
		{
			//only for rectangular cross-section!:
			xgd.SetLen(SOS());
			GetDrawCoordinates(xgd);

			xgd = T*xgd;

			double x0 = ploc.X();
			double yb = 0.5*ploc.Y()*ly;

			double Gm = Em / (2.*(1+nu));
			double A = ly*lz;

			double ks = 10.*(1+nu)/(12.+11.*nu); //shear correction factors (square): 10*(1+nu)/(12.+11.*nu) 

			static Vector Sx; Sx.SetLen(NS());   //rx
			static Vector Sy; Sy.SetLen(NS());   //ry
			static Vector Sxx; Sxx.SetLen(NS()); //rxx
			static Vector Syx; Syx.SetLen(NS()); //ryx
			Vector2D rx,ry,rxx,ryx;

			Vector2D rx1,ry1,ryx1;
			Vector2D rx2,ry2,ryx2;

			Matrix3D stress(0);

			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			//insert new computation of stresses according to Reissner!!!

			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			stress(2,1) = stress(1,2);
			stress(3,1) = stress(1,3);

			if(fvd.VariableType() == FieldVariableDescriptor::FVT_stress_mises)
			{
				return stress.Mises();
			}
			else
			{
				return fvd.GetComponent(stress);
			}
		}
	}
	return FIELD_VARIABLE_NO_VALUE;
}

void ANCFBeam2D::DrawElement() 
{
	mbs->SetColor(col);

	double lx1 = 1; double ly1 = 1*GetMBS()->GetMagnifyYZ();

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

	double oldlinethick = GetMBS()->GetDOption(102);


	int colormode = 0;
	if (GetMBS()->GetActualPostProcessingFieldVariable() != NULL) colormode = 1;

	if (colormode || (linemode != 0))
	{
		//for (int kk=1; kk <= 1+GetMBS()->GetIOption(116); kk++) //draw original shape ...
		if (1)
		{
			const FieldVariableDescriptor * fvd = GetMBS()->GetActualPostProcessingFieldVariable();

			double modeval = 0;
			int xgset = 0;

			double tilex = GetMBS()->GetIOption(137);
			double tiley = GetMBS()->GetIOption(138);

			TArray<Vector3D> points((int)(tilex+1)*(int)(tiley+1));
			TArray<double> vals((int)(tilex+1)*(int)(tiley+1));
			double v=0;

			points.SetLen(0); vals.SetLen(0);
			Vector2D p0, vx, vy;
			int tileyn = (int)tiley;
			int tilexn = (int)tilex;

			p0 = Vector2D(-lx1,-ly1);
			vx = Vector2D(2.*lx1/tilexn,0);
			vy = Vector2D(0,2.*ly1/tileyn);

			for (double iy = 0; iy <= tileyn+1e-10; iy++)
			{
				for (double ix = 0; ix <= tilexn+1e-10; ix++)
				{
					Vector2D ploc = (p0+ix*vx+iy*vy);
					Vector2D pg;
					//if (kk==1)
					//{
					//	pg = GetPos2D0D(ploc);
					//	if (GetMBS()->GetIOption(116) && 0)
					//	{
					//		Vector2D pgi = GetPos2D0Dinit(ploc);
					//		pg += fact*(pg-pgi);
					//	}

					//	GetMBS()->GetDOption(102) = oldlinethick;
					//	//global_uo << "p=" << pg << ", ploc=" << ploc << "\n";
					//}
					//else
					//{
					//	pg = GetPos2D0Dinit(ploc);
					//	GetMBS()->GetDOption(102) = 1;
					//}
					pg = GetPos2D0D(ploc);

					Vector3D p(ToP3D(pg));
					points.Add(p);
					if (colormode)
						v = GetFieldVariableValue(*fvd, ploc, true);
					vals.Add(v);
				}
			}
			mbs->DrawColorQuads(points,vals,(int)tilexn+1,(int)tileyn+1,colormode,linemode);
		}
		GetMBS()->GetDOption(102) = oldlinethick;
	}
	else
	{
		int drawgrid = 0;
		double thickness = 1;
		//double tiling = 20;
		double tiling = GetMBS()->GetIOption(136);
		mbs->SetDrawlines(0);
		mbs->SetDrawlinesH(0);

		Vector3D p1,p2,p3,p4;
		for (double i = 0; i < tiling; i++)
		{
			double l1 = -lx1+2*lx1*i/tiling;
			double l2 = -lx1+2*lx1*(i+1)/tiling;

			p1 = Vector3D(ToP3D(GetPos2D0D(Vector2D(l1,-ly1))));
			p2 = Vector3D(ToP3D(GetPos2D0D(Vector2D(l2,-ly1))));
			p3 = Vector3D(ToP3D(GetPos2D0D(Vector2D(l2, ly1))));
			p4 = Vector3D(ToP3D(GetPos2D0D(Vector2D(l1, ly1))));

			if (drawgrid)
			{
				mbs->MyDrawLine(p1,p2,thickness);
				mbs->MyDrawLine(p2,p3,thickness);
				mbs->MyDrawLine(p3,p4,thickness);
				mbs->MyDrawLine(p4,p1,thickness);
			}
			else
			{
				mbs->DrawQuad(p1,p2,p3,p4);
			}
		}
		mbs->SetDrawlines(0);
	}

	if (concentratedmass1 != 0) 
	{
		Vector3D mp1(ToP3D(GetPos2D0D(Vector2D(-lx1,0)))); 
		mbs->SetColor(colred);
		//double r = pow(3.*concentratedmass1/(4.*MY_PI*GetRho()),1./3.); 
		double r = ly/10.;	//$ DR 2013-02-04 deleted rho from class element, do not use it here!
		//double r = ly*GetMBS()->GetMagnifyYZ()*4.;
		mbs->DrawSphere(mp1,r,10);
	}
	if (concentratedmass2 != 0) 
	{
		Vector3D mp1(ToP3D(GetPos2D0D(Vector2D(lx1,0)))); 
		mbs->SetColor(colred);
		//double r = pow(3.*concentratedmass2/(4.*MY_PI*GetRho()),1./3.); 
		double r = ly/10.;	//$ DR 2013-02-04 deleted rho from class element, do not use it here!
		//double r = ly*GetMBS()->GetMagnifyYZ()*4.;
		mbs->DrawSphere(mp1,r,12);
	}
};




//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//Element shares nodes with other elements, n1 and n2 are nodenumbers; element sets initial conditions for nodes
ANCFBeam2Dlin::ANCFBeam2Dlin(MBS* mbsi, const Vector& xc1, const Vector& xc2, const Vector& vc1, const Vector& vc2, 
											 int n1i, int n2i, double rhoi, double Emi, double nui,
											 const Vector3D& si, const Vector3D& coli, int beammodel): ANCFBeam2D(mbsi)
{
	
	elasticforce_beam = beammodel;
	n1=n1i; n2=n2i; sos2=0;
	size = si;
	lx = si.X(); ly = si.Y(); lz = si.Z();

	mass = lx*ly*lz*rhoi;
	//lx=1;ly=1;lz=1;

	x_init = xc1.Append(xc2);
	xg = xc1.Append(xc2);
	Vector v_init = vc1.Append(vc2);

	nu = nui;
	Em = Emi;
	rho = rhoi;
	col = coli;

	double EA, EI, GAks, rhoA, rhoI;
	SetRectangularBeamParameters(EA, EI, GAks, rhoA, rhoI);
	SetBeamParameters(EA, EI, GAks, rhoA, rhoI);

	InitConstructor();
	BuildDSMatrices();

	e0 = x_init;
	x_init = x_init.Append(v_init); //Velocity initial conditions can also be transformed by Tinv!!!

	SetInitialCurvatures(); //must be done after x_init!!!

	elementname = GetElementSpec();
};




void ANCFBeam2Dlin::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{

	Body2D::GetElementData(edc);

	ElementData ed;

	ed.SetDouble(Em, "Youngs_modulus"); edc.Add(ed);
	ed.SetDouble(nu, "Poisson_ratio"); edc.Add(ed);
	ed.SetVector3D(lx, ly, lz, "Beam_dimensions"); edc.Add(ed);

	ed.SetInt(n1, "Node_number1", 1, GetMBS()->NNodes()); edc.Add(ed);
	ed.SetInt(n2, "Node_number2", 1, GetMBS()->NNodes()); edc.Add(ed);

	ed.SetVector2D(concentratedmass1, concentratedmass2, "Node_masses"); edc.Add(ed);

	ed.SetBool(elasticforce_beam!=0, "Elastic_line_model"); edc.Add(ed);
/*
	ed.SetDouble(beamEA, "Axial_stiffness_EA"); edc.Add(ed);
	ed.SetDouble(beamEI, "Bending_stiffness_EI"); edc.Add(ed);
	ed.SetDouble(beamGAks, "Shear_stiffness_GA"); edc.Add(ed);
	ed.SetDouble(beamRhoA, "Cross_section_inertia_RhoA"); edc.Add(ed);
	ed.SetDouble(beamRhoI, "Moment_of_inertia_RhoI"); edc.Add(ed);
*/

	Vector xinit = x_init.SubVector(1, SOS());
	Vector vinit = x_init.SubVector(SOS()+1, 2*SOS());

	ed.SetVector2D(xinit( 1), xinit( 2), "Node1_r"); edc.Add(ed);
	ed.SetVector2D(xinit( 3), xinit( 4), "Node1_ry"); edc.Add(ed);

	ed.SetVector2D(xinit( 1+4), xinit( 2+4), "Node2_r"); edc.Add(ed);
	ed.SetVector2D(xinit( 3+4), xinit( 4+4), "Node2_ry"); edc.Add(ed);

	ed.SetVector2D(vinit( 1), vinit( 2), "Node1_v"); edc.Add(ed);
	ed.SetVector2D(vinit( 3), vinit( 4), "Node1_vy"); edc.Add(ed);

	ed.SetVector2D(vinit( 1+4), vinit( 2+4), "Node2_v"); edc.Add(ed);
	ed.SetVector2D(vinit( 3+4), vinit( 4+4), "Node2_vy"); edc.Add(ed);

}

int ANCFBeam2Dlin::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Body2D::SetElementData(edc);

	InitConstructor();

	double EA, EI, GAks, rhoA, rhoI;
	SetRectangularBeamParameters(EA, EI, GAks, rhoA, rhoI);
	SetBeamParameters(EA, EI, GAks, rhoA, rhoI);

	sos2 = 0;
	Em = 0;
	elasticforce_beam = 0;

	GetElemDataDouble(GetMBS(), edc, "Youngs_modulus", Em, 1);
	GetElemDataDouble(GetMBS(), edc, "Poisson_ratio", nu, 1);

	GetElemDataVector3D(GetMBS(), edc, "Beam_dimensions", lx, ly, lz, 1); 

	GetElemDataInt(GetMBS(), edc, "Node_number1", n1, 1);
	GetElemDataInt(GetMBS(), edc, "Node_number2", n2, 1);


	if (n1 < 1 || n1 > GetMBS()->NNodes() || n2 < 1 || n2 > GetMBS()->NNodes())
	{
		GetMBS()->EDCError("Invalid node numbers in ANCF beam element. Element has not been created!");
		return 0;
	}

	GetElemDataVector2D(GetMBS(), edc, "Node_masses", concentratedmass1, concentratedmass2, 0);

	GetElemDataBool(GetMBS(), edc, "Elastic_line_model", elasticforce_beam, 0);
	if (elasticforce_beam) elasticforce_beam = 1;
/*
	GetElemDataDouble(GetMBS(), edc, "Axial_stiffness_EA", beamEA, 1);
	GetElemDataDouble(GetMBS(), edc, "Bending_stiffness_EI", beamEI, 1);
	GetElemDataDouble(GetMBS(), edc, "Shear_stiffness_GA", beamGAks, 1);
	GetElemDataDouble(GetMBS(), edc, "Cross_section_inertia_RhoA", beamRhoA, 1);
	GetElemDataDouble(GetMBS(), edc, "Moment_of_inertia_RhoI", beamRhoI, 1);
*/
	Vector xinit(SOS());
	Vector vinit(SOS());

	GetElemDataVector2D(GetMBS(), edc, "Node1_r" , xinit( 1), xinit( 2), 1);
	GetElemDataVector2D(GetMBS(), edc, "Node1_ry", xinit( 3), xinit( 4), 1);

	GetElemDataVector2D(GetMBS(), edc, "Node2_r" , xinit( 1+4), xinit( 2+4), 1);
	GetElemDataVector2D(GetMBS(), edc, "Node2_ry", xinit( 3+4), xinit( 4+4), 1);

	GetElemDataVector2D(GetMBS(), edc, "Node1_v" , vinit( 1), vinit( 2), 1);
	GetElemDataVector2D(GetMBS(), edc, "Node1_vy", vinit( 3), vinit( 4), 1);

	GetElemDataVector2D(GetMBS(), edc, "Node2_v" , vinit( 1+4), vinit( 2+4), 1);
	GetElemDataVector2D(GetMBS(), edc, "Node2_vy", vinit( 3+4), vinit( 4+4), 1);

	size = Vector3D(lx,ly,lz);
	mass = GetMaterial().BeamRhoA()*lx;


	BuildDSMatrices();
	e0 = xinit;
	x_init = xinit.Append(vinit); //Velocity initial conditions can also be transformed by Tinv!!!

	SetInitialCurvatures(); //must be done after x_init!!!	

	return rv;
}



void ANCFBeam2Dlin::GetS0(Vector& sf, const Vector2D& ploc) const
{
	double xb = ploc.X();
	double yb = ploc.Y();
	sf.SetLen(NS());

	sf(1) = 1.0/2.0-xb/2.0;
	sf(2) = -(-1.0+xb)*yb*ly/4.0;
	sf(3) = 1.0/2.0+xb/2.0;
	sf(4) = (1.0+xb)*yb*ly/4.0;
}

void ANCFBeam2Dlin::GetDSMatrix0(const Vector2D& ploc, Matrix& sf) const
{
	double xb = ploc.X();
	double yb = ploc.Y();
	sf.SetSize(2,NS());

	sf(1,1) = -1.0/2.0;
	sf(1,2) = -yb*ly/4.0;
	sf(1,3) = 1.0/2.0;
	sf(1,4) = yb*ly/4.0;
	sf(2,1) = 0.0;
	sf(2,2) = -(-1.0+xb)*ly/4.0;
	sf(2,3) = 0.0;
	sf(2,4) = (1.0+xb)*ly/4.0;
}

void ANCFBeam2Dlin::GetS0x(Vector& sfx, const Vector2D& ploc) const
{
	//at y=0 !!!!!!!!!!!!!!!!!!!!!
	//double xb = ploc.X();
	sfx.SetLen(NS());
	double f = 2./lx;

	sfx(1) = f*(-1.0/2.0);
	sfx(2) = 0;
	sfx(3) = f*(1.0/2.0);
	sfx(4) = 0;
}

void ANCFBeam2Dlin::GetS0y(Vector& sfx, const Vector2D& ploc) const
{
	double xb = ploc.X();
	sfx.SetLen(NS());
	double f = 2./ly;

	sfx(1) = 0.0;
	sfx(2) = f*(-(-1.0+xb)*ly/4.0);
	sfx(3) = 0.0;
	sfx(4) = f*((1.0+xb)*ly/4.0);
}

void ANCFBeam2Dlin::GetS0yx(Vector& sfx, const Vector2D& ploc) const
{
	double xb = ploc.X();
	sfx.SetLen(NS());
	double f = 4./(ly*lx);

	sfx(1) = 0.0;
	sfx(2) = f*(-ly/4.0);
	sfx(3) = 0.0;
	sfx(4) = f*(ly/4.0);
}

void ANCFBeam2Dlin::GetS0xx(Vector& sfxx, const Vector2D& ploc) const
{
	double xb = ploc.X();
	sfxx.SetLen(NS());

	sfxx(1) = 0;
	sfxx(2) = 0;
	sfxx(3) = 0;
	sfxx(4) = 0;
}

