//#**************************************************************
//#
//# filename:             ANCFBeam3D.cpp
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						9. November 2004
//# description:          ANCF beam formulation according to Yakoub and Shabana
//#                       Elastic line approach according to Schwab and Meijaard
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
 
#include "ANCFBeam3D.h"
#include "elementDataAccess.h"
#include "Material.h"
#include "node.h"
#include "femathHelperFunctions.h"
#include "graphicsConstants.h"

const int usetangentframe = 0;

ANCFBeam3D::ANCFBeam3D(MBS* mbsi, const Vector& xc1, const Vector& xc2, double rhoi, double Emi, double nui,
											 const Vector3D& si, const Vector3D& coli, int beammodel):
Body3D(mbsi),massmatrix(), Hmatrix(), SV(), DS(), x1(), x2(), x3(), w1(), w2(), w3(),
T1(), T2(), T(), xg(), xgd(), e0()
{
	elasticforce_beam = beammodel;
	n1=0; n2=0; sos2=SOS();
	size = si;
	lx = si.X(); ly = si.Y(); lz = si.Z();

	mass = lx*ly*lz*rhoi;
	//lx=1;ly=1;lz=1;

	x_init = xc1.Append(xc2);
	xg = xc1.Append(xc2);

	nu = nui;
	Em = Emi;
	//rho = rhoi; //$ DR 2013-02-04 deleted rho from class element, do not use it here!
	col = coli;

	SetRectangularBeamParameters();
	InitConstructor();
	BuildDSMatrices();
	Matrix Tinv = T;
	if (!Tinv.Invert2()) {UO() << "ERROR: T matrix not invertible!!!\n";}
	e0 = x_init;
	x_init = Tinv * x_init;
	x_init = x_init.Append(Vector(SOS())); //Velocity initial conditions can also be transformed by Tinv!!!

	SetInitialCurvatures(); //must be done after x_init!!!

	elementname = GetElementSpec();
};

//Element shares nodes with other elements, n1 and n2 are nodenumbers; element sets initial conditions for nodes
ANCFBeam3D::ANCFBeam3D(MBS* mbsi, const Vector& xc1, const Vector& xc2, int n1i, int n2i, double rhoi, double Emi, double nui,
											 const Vector3D& si, const Vector3D& coli, int beammodel):
Body3D(mbsi),massmatrix(), Hmatrix(), SV(), DS(), x1(), x2(), x3(), w1(), w2(), w3(),
T1(), T2(), T(), xg(), xgd(), e0()
{
	elasticforce_beam = beammodel;
	n1=n1i; n2=n2i; sos2=0;
	size = si;
	lx = si.X(); ly = si.Y(); lz = si.Z();

	//global_uo << "beam(" << n1i << "," << n2i << ")\n";

	mass = lx*ly*lz*rhoi;
	//lx=1;ly=1;lz=1;

	x_init = xc1.Append(xc2);
	xg = xc1.Append(xc2);

	nu = nui;
	Em = Emi;
	//rho = rhoi; //$ DR 2013-02-04 deleted rho from class element, do not use it here!
	col = coli;
	SetRectangularBeamParameters();

	InitConstructor();
	BuildDSMatrices();
	Matrix Tinv = T;
	if (!Tinv.Invert2()) {UO() << "ERROR: T matrix not invertible!!!\n";}
	e0 = x_init;
	x_init = Tinv * x_init;
	x_init = x_init.Append(Vector(SOS())); //Velocity initial conditions can also be transformed by Tinv!!!

	SetInitialCurvatures(); //must be done after x_init!!!	

	elementname = GetElementSpec();
};

//Element shares nodes with other elements, n1 and n2 are nodenumbers; element sets initial conditions for nodes
ANCFBeam3D::ANCFBeam3D(MBS* mbsi, const Vector& xc1, const Vector& xc2, const Vector& vc1, const Vector& vc2, int n1i, int n2i, double rhoi, double Emi, double nui,
											 const Vector3D& si, const Vector3D& coli, int beammodel):
Body3D(mbsi),massmatrix(), Hmatrix(), SV(), DS(), x1(), x2(), x3(), w1(), w2(), w3(),
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
	SetRectangularBeamParameters();

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



void ANCFBeam3D::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{

	Body3D::GetElementData(edc);

	ElementData ed;

	ed.SetDouble(Em, "Youngs_modulus"); edc.Add(ed);
	ed.SetDouble(nu, "Poisson_ratio"); edc.Add(ed);
	ed.SetVector3D(lx, ly, lz, "Beam_dimensions"); edc.Add(ed);

	ed.SetInt(n1, "Node_number1", 1, GetMBS()->NNodes()); edc.Add(ed);
	ed.SetInt(n2, "Node_number2", 1, GetMBS()->NNodes()); edc.Add(ed);

	ed.SetVector2D(concentratedmass1, concentratedmass2, "Node_masses"); edc.Add(ed);

	ed.SetBool(elasticforce_beam!=0, "Elastic_line_model"); edc.Add(ed);

	ed.SetDouble(beamEA, "Axial_stiffness"); edc.Add(ed);
	ed.SetVector2D(beamEIy, beamEIz, "Bending_stiffness"); edc.Add(ed);
	ed.SetVector2D(beamGAky, beamGAkz, "Shear_stiffness"); edc.Add(ed);
	ed.SetDouble(beamGJkx, "Torsional_stiffness"); edc.Add(ed);

	Matrix Tinv = T;
	Vector xinit = x_init.SubVector(1, SOS());
	Vector vinit = x_init.SubVector(SOS()+1, 2*SOS());
	xinit = T*xinit;
	vinit = T*vinit;

	ed.SetVector3D(xinit( 1), xinit( 2), xinit( 3), "Node1_r"); edc.Add(ed);
	ed.SetVector3D(xinit( 4), xinit( 5), xinit( 6), "Node1_rx"); edc.Add(ed);
	ed.SetVector3D(xinit( 7), xinit( 8), xinit( 9), "Node1_ry"); edc.Add(ed);
	ed.SetVector3D(xinit(10), xinit(11), xinit(12), "Node1_rz"); edc.Add(ed);

	ed.SetVector3D(xinit( 1+12), xinit( 2+12), xinit( 3+12), "Node2_r"); edc.Add(ed);
	ed.SetVector3D(xinit( 4+12), xinit( 5+12), xinit( 6+12), "Node2_rx"); edc.Add(ed);
	ed.SetVector3D(xinit( 7+12), xinit( 8+12), xinit( 9+12), "Node2_ry"); edc.Add(ed);
	ed.SetVector3D(xinit(10+12), xinit(11+12), xinit(12+12), "Node2_rz"); edc.Add(ed);

	ed.SetVector3D(vinit( 1), vinit( 2), vinit( 3), "Node1_v"); edc.Add(ed);
	ed.SetVector3D(vinit( 4), vinit( 5), vinit( 6), "Node1_vx"); edc.Add(ed);
	ed.SetVector3D(vinit( 7), vinit( 8), vinit( 9), "Node1_vy"); edc.Add(ed);
	ed.SetVector3D(vinit(10), vinit(11), vinit(12), "Node1_vz"); edc.Add(ed);

	ed.SetVector3D(vinit( 1+12), vinit( 2+12), vinit( 3+12), "Node2_v"); edc.Add(ed);
	ed.SetVector3D(vinit( 4+12), vinit( 5+12), vinit( 6+12), "Node2_vx"); edc.Add(ed);
	ed.SetVector3D(vinit( 7+12), vinit( 8+12), vinit( 9+12), "Node2_vy"); edc.Add(ed);
	ed.SetVector3D(vinit(10+12), vinit(11+12), vinit(12+12), "Node2_vz"); edc.Add(ed);

}

int ANCFBeam3D::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Body3D::SetElementData(edc);

	InitConstructor();
	SetRectangularBeamParameters();

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
	/*
	if (GetMBS()->GetNode(n1).SOS() != 12) 
	{
	GetMBS()->GetNode(n1).SOS() = 12;
	GetMBS()->EDCError("ANCFbeam Node does not have 12 coordinates; the size has been corrected");
	}
	if (GetMBS()->GetNode(n2).SOS() != 12) 
	{
	GetMBS()->GetNode(n2).SOS() = 12;
	GetMBS()->EDCError("ANCFbeam Node does not have 12 coordinates; the size has been corrected");
	}*/

	GetElemDataVector2D(GetMBS(), edc, "Node_masses", concentratedmass1, concentratedmass2, 0);

	GetElemDataBool(GetMBS(), edc, "Elastic_line_model", elasticforce_beam, 0);
	if (elasticforce_beam) elasticforce_beam = 3;

	GetElemDataDouble(GetMBS(), edc, "Axial_stiffness", beamEA, 1);
	GetElemDataVector2D(GetMBS(), edc, "Bending_stiffness", beamEIy, beamEIz, 1);
	GetElemDataVector2D(GetMBS(), edc, "Shear_stiffness", beamGAky, beamGAkz, 1);
	GetElemDataDouble(GetMBS(), edc, "Torsional_stiffness", beamGJkx, 1);

	Vector xinit(SOS());
	Vector vinit(SOS());

	GetElemDataVector3D(GetMBS(), edc, "Node1_r" , xinit( 1), xinit( 2), xinit( 3), 1); 
	GetElemDataVector3D(GetMBS(), edc, "Node1_rx", xinit( 4), xinit( 5), xinit( 6), 1); 
	GetElemDataVector3D(GetMBS(), edc, "Node1_ry", xinit( 7), xinit( 8), xinit( 9), 1); 
	GetElemDataVector3D(GetMBS(), edc, "Node1_rz", xinit(10), xinit(11), xinit(12), 1); 

	GetElemDataVector3D(GetMBS(), edc, "Node2_r" , xinit( 1+12), xinit( 2+12), xinit( 3+12), 1); 
	GetElemDataVector3D(GetMBS(), edc, "Node2_rx", xinit( 4+12), xinit( 5+12), xinit( 6+12), 1); 
	GetElemDataVector3D(GetMBS(), edc, "Node2_ry", xinit( 7+12), xinit( 8+12), xinit( 9+12), 1); 
	GetElemDataVector3D(GetMBS(), edc, "Node2_rz", xinit(10+12), xinit(11+12), xinit(12+12), 1); 

	GetElemDataVector3D(GetMBS(), edc, "Node1_v" , vinit( 1), vinit( 2), vinit( 3), 1); 
	GetElemDataVector3D(GetMBS(), edc, "Node1_vx", vinit( 4), vinit( 5), vinit( 6), 1); 
	GetElemDataVector3D(GetMBS(), edc, "Node1_vy", vinit( 7), vinit( 8), vinit( 9), 1); 
	GetElemDataVector3D(GetMBS(), edc, "Node1_vz", vinit(10), vinit(11), vinit(12), 1); 

	GetElemDataVector3D(GetMBS(), edc, "Node2_v" , vinit( 1+12), vinit( 2+12), vinit( 3+12), 1); 
	GetElemDataVector3D(GetMBS(), edc, "Node2_vx", vinit( 4+12), vinit( 5+12), vinit( 6+12), 1); 
	GetElemDataVector3D(GetMBS(), edc, "Node2_vy", vinit( 7+12), vinit( 8+12), vinit( 9+12), 1); 
	GetElemDataVector3D(GetMBS(), edc, "Node2_vz", vinit(10+12), vinit(11+12), vinit(12+12), 1); 

	size = Vector3D(lx,ly,lz);
	mass = lx*ly*lz*rho;


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


int ANCFBeam3D::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = Element::CheckConsistency(errorstr);
	if (rv) return rv;


	if (n1 && n2)
	{
		if (GetMBS()->GetNode(n1).SOS() != 12) 
		{
			GetMBS()->GetNode(n1).SOS() = 12;
			errorstr += mystr("ANCFbeam ")+ mystr(GetOwnNum()) + mystr(": Node does not have 12 coordinates; the size has been corrected");
		}
		if (GetMBS()->GetNode(n2).SOS() != 12) 
		{
			GetMBS()->GetNode(n2).SOS() = 12;
			errorstr += mystr("ANCFbeam ")+ mystr(GetOwnNum()) + mystr(": Node does not have 12 coordinates; the size has been corrected");
		}
	}

	return rv;
}



void ANCFBeam3D::LinkToElements()
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

void ANCFBeam3D::BuildDSMatrices() 
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
	orderWt = 2; //2; 1=2=4, 1 does not work in L-shape!!!
	orderWb = 6; //3; 3=4=6, 2 gives shear locking

	factstiffWl = 1; //reduce stiffness of thickness modes, standard==1

	GetIntegrationRuleLobatto(x1,w1,orderWl); int nipWl = x1.Length();
	//GetIntegrationRuleLobatto(x1,w1,orderWs); int nipWs = x1.Length();
	GetIntegrationRuleLobatto(x1,w1,orderWsHR); int nipWsHR = x1.Length();
	GetIntegrationRule(x1,w1,orderWt); int nipWt = x1.Length();
	GetIntegrationRule(x1,w1,orderWb); int nipWb = x1.Length();

	Wlepsx0.SetLen(nipWl);
	Wlepsy0.SetLen(nipWl);
	Wlepsz0.SetLen(nipWl);
	Wlgamyz0.SetLen(nipWl);
	WsHRk10.SetLen(nipWsHR);
	WsHRk20.SetLen(nipWsHR);
	Wtkapx0.SetLen(nipWt);
	Wbkapy0.SetLen(nipWb);
	Wbkapz0.SetLen(nipWb);

	Wlepsx0.SetAll(0);
	Wlepsy0.SetAll(0);
	Wlepsz0.SetAll(0);
	Wlgamyz0.SetAll(0);
	WsHRk10.SetAll(0);
	WsHRk20.SetAll(0);
	Wtkapx0.SetAll(0);
	Wbkapy0.SetAll(0);
	Wbkapz0.SetAll(0);

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	GetIntegrationRule(x1,w1,orderx);
	GetIntegrationRule(x2,w2,orderyz);
	GetIntegrationRule(x3,w3,orderyz);
	Matrix3D jac, jacinv;
	int kx1 = x2.Length()*x3.Length();
	int kx2 = x3.Length();

	//UO() << "Initialize Element\n";

	for (int i1=1; i1<=x1.GetLen(); i1++)
	{
		for (int i2=1; i2<=x2.GetLen(); i2++)
		{
			for (int i3=1; i3<=x3.GetLen(); i3++)
			{
				int ind = (i1-1)*kx1+(i2-1)*kx2+(i3-1);
				//UO() << "ind=" << ind << "\n";

				Vector3D p(x1(i1),x2(i2),x3(i3));

				GetDSMatrix0(p,DS);

				GetJacobi(jac,p,DS,x_init);
				jacdet[ind] = jac.Det();
				jac.GetInverse(jacinv);
				jacinv = jacinv.GetTp();

				grad[ind].Init();
				grad[ind].SetSize(3,NS());
				Mult(jacinv, DS, grad[ind]);
			}
		}
	}
	//Build T-matrices:
	T1.Init();
	T2.Init();
	T1.SetSize(Dim()*NS()/2,Dim()*NS()/2);
	T1.SetAll(0);
	T2=T1;
	Vector3D p1(-1,0,0);
	Vector3D p2( 1,0,0);

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
			jac1(i,j) = x_init((i-1)*3+j+3);
			jac2(i,j) = x_init((i-1)*3+j+12+3);
		}
	}

	T1(1,1) = 1; T1(2,2) = 1; T1(3,3) = 1;
	T2 = T1;
	for (int i=1; i <= Dim(); i++)
	{
		for (int j=1; j <= Dim(); j++)
		{
			T1(i*3+1,j*3+1) = jac1(i,j);
			T1(i*3+2,j*3+2) = jac1(i,j);
			T1(i*3+3,j*3+3) = jac1(i,j);
			T2(i*3+1,j*3+1) = jac2(i,j);
			T2(i*3+2,j*3+2) = jac2(i,j);
			T2(i*3+3,j*3+3) = jac2(i,j);
		}
	}


	T.Init();
	T.SetSize(SOS(),SOS());
	T.SetAll(0);
	T.SetSubmatrix(T1,1,1);
	T.SetSubmatrix(T2,1+T1.Getrows(),1+T1.Getcols());

	//UO() << "T=" << T << "\n";
};

void ANCFBeam3D::SetBeamParameters(double EA, double EIw, double EIv, double GJ, double GA)
{
	//Em and nu already set
	//old function that assumes rectangular cross-section
	double ky = 10*(1+nu)/(12.+11.*nu); //shear correction factors (square): 10*(1+nu)/(12.+11.*nu) 
	double kz = 10*(1+nu)/(12.+11.*nu); //shear correction factors (square): 10*(1+nu)/(12.+11.*nu)
	double kx = 0.8436; //shear correction factor torsion  : 0.8436 (square)

	double A = ly*lz;
	double Gm = Em / (2.*(1+nu));

	beamEA	 = EA; 
	beamEIy  = EIw;
	beamEIz  = EIv;
	beamGJkx = GJ*kx;
	beamGAky = GA*ky;
	beamGAkz = GA*kz;

	SetInitialCurvatures();
}

void ANCFBeam3D::SetBeamParameters2(double beamEAi, double beamEIyi, double beamEIzi, double beamGJkxi, 
																		double beamGAkyi, double beamGAkzi)
{
	beamEA	 = beamEAi; 
	beamEIy  = beamEIyi;
	beamEIz  = beamEIzi;
	beamGJkx = beamGJkxi;
	beamGAky = beamGAkyi;
	beamGAkz = beamGAkzi;

	SetInitialCurvatures();
}

void ANCFBeam3D::SetRectangularBeamParameters() //set beam parameters for rectangular cross-section, Schwab-Meijaard model
{
	double A = ly*lz;
	double Gm = Em / (2.*(1+nu));
	double GA = Gm*A;

	beamEIy = Em * ly*Cub(lz)/12.;
	beamEIz = Em * lz*Cub(ly)/12.;
	beamEA = Em * A;

	double GIT = Gm * ly*lz*(Sqr(ly) + Sqr(lz))/12.; //Schwab/Meijaard
	//GIT = Gm * 0.140 * Sqr(ly)*Sqr(lz); //acc. to Mang or Gieck for ly==lz

	double ky = 10*(1+nu)/(12.+11.*nu); //shear correction factors (square): 10*(1+nu)/(12.+11.*nu) 
	double kz = 10*(1+nu)/(12.+11.*nu); //shear correction factors (square): 10*(1+nu)/(12.+11.*nu)
	double kx = 0.8436; //shear correction factor torsion  : 0.8436 (square)

	beamGJkx = GIT * kx;
	beamGAky = GA*ky;
	beamGAkz = GA*kz;

	//UO() << "GAky=" << beamGAky << ", GAkz=" << beamGAkz << "\n";
}



void ANCFBeam3D::ApplyT(Vector& v) const
{
	static Vector x;
	x = v;
	for (int i=1; i<=3; i++)
	{
		for (int j=1; j<=3; j++)
		{
			v((i-1)*3+j+3) =    jac1(i,1)*x(0+j+3) +    jac1(i,2)*x(3+j+3) +    jac1(i,3)*x(6+j+3);
			v((i-1)*3+j+3+12) = jac2(i,1)*x(0+j+3+12) + jac2(i,2)*x(3+j+3+12) + jac2(i,3)*x(6+j+3+12);
		}
	}
	//global_uo << "T*v2=" << v << "\n";
}

void ANCFBeam3D::ApplyTD(Vector& v) const
{
	static Vector x;
	x = v;
	for (int i=1; i<=3; i++)
	{
		for (int j=1; j<=3; j++)
		{
			v((i-1)*3+j+3) = jac1(i,1)*x(0+j+3) + jac1(i,2)*x(3+j+3) + jac1(i,3)*x(6+j+3);
			v((i-1)*3+j+3+12) = jac2(i,1)*x(0+j+3+12) + jac2(i,2)*x(3+j+3+12) + jac2(i,3)*x(6+j+3+12);
		}
	}
}

void ANCFBeam3D::ApplyTtp(Vector& v) const
{
	//v=T.GetTp()*v; return;
	static Vector x;
	x = v;
	for (int i=1; i<=3; i++)
	{
		for (int j=1; j<=3; j++)
		{
			v((i-1)*3+j+3) = jac1(1,i)*x(0+j+3) + jac1(2,i)*x(3+j+3) + jac1(3,i)*x(6+j+3);
			v((i-1)*3+j+3+12) = jac2(1,i)*x(0+j+3+12) + jac2(2,i)*x(3+j+3+12) + jac2(3,i)*x(6+j+3+12);
		}
	}
}

void ANCFBeam3D::ApplyTtp(Matrix& m) const
{
	//v=T.GetTp()*v; return;
	static Vector x;
	int sns = SOS();
	x.SetLen(sns);

	for (int k = 1; k <= m.Getcols(); k++)
	{
		for (int i = 1+3; i <= sns; i++)
		{
			x(i) = m(i,k);
		}
		for (int i=1; i<=3; i++)
		{
			for (int j=1; j<=3; j++)
			{
				m((i-1)*3+j+3,k) = jac1(1,i)*x(0+j+3) + jac1(2,i)*x(3+j+3) + jac1(3,i)*x(6+j+3);
				m((i-1)*3+j+3+12,k) = jac2(1,i)*x(0+j+3+12) + jac2(2,i)*x(3+j+3+12) + jac2(3,i)*x(6+j+3+12);
			}
		}
	}
}

void ANCFBeam3D::GetS0(Vector& sf, const Vector3D& ploc) const
{
	double xb = ploc.X();
	double yb = ploc.Y();
	double zb = ploc.Z();
	sf.SetLen(8);
	double xb2 = xb*xb;
	double xb3 = xb2*xb;
	sf(1) = 1.0/2.0-3.0/4.0*xb+xb3/4.0;
	sf(2) = lx*(1.0-xb-xb2+xb3)/8.0;
	sf(3) = -yb*ly*(-1.0+xb)/4.0;
	sf(4) = -zb*lz*(-1.0+xb)/4.0;
	sf(5) = 1.0/2.0+3.0/4.0*xb-xb3/4.0;
	sf(6) = lx*(1.0+xb)*(1.0+xb)*(-1.0+xb)/8.0;
	sf(7) = (1.0+xb)*yb*ly/4.0;
	sf(8) = (1.0+xb)*zb*lz/4.0;
}

void ANCFBeam3D::GetDSMatrix0(const Vector3D& ploc, Matrix& sf) const
{
	double xb = ploc.X();
	double yb = ploc.Y();
	double zb = ploc.Z();
	sf.SetSize(3,8);

	sf(1,1) = -3.0/4.0+3.0/4.0*xb*xb;
	sf(1,2) = lx*(-1.0-2.0*xb+3.0*xb*xb)/8.0;
	sf(1,3) = -yb*ly/4.0;
	sf(1,4) = -zb*lz/4.0;
	sf(1,5) = 3.0/4.0-3.0/4.0*xb*xb;
	sf(1,6) = lx*(1.0+xb)*(-1.0+xb)/4.0+lx*(1.0+xb)*(1.0+xb)/8.0;
	sf(1,7) = yb*ly/4.0;
	sf(1,8) = zb*lz/4.0;
	sf(2,1) = 0.0;
	sf(2,2) = 0.0;
	sf(2,3) = -ly*(-1.0+xb)/4.0;
	sf(2,4) = 0.0;
	sf(2,5) = 0.0;
	sf(2,6) = 0.0;
	sf(2,7) = (1.0+xb)*ly/4.0;
	sf(2,8) = 0.0;
	sf(3,1) = 0.0;
	sf(3,2) = 0.0;
	sf(3,3) = 0.0;
	sf(3,4) = -lz*(-1.0+xb)/4.0;
	sf(3,5) = 0.0;
	sf(3,6) = 0.0;
	sf(3,7) = 0.0;
	sf(3,8) = (1.0+xb)*lz/4.0;

}

void ANCFBeam3D::GetS0x(Vector& sfx, const Vector3D& ploc) const
{
	//at y=z=0 !!!!!!!!!!!!!!!!!!!!!
	double xb = ploc.X();
	sfx.SetLen(NS());
	double f = 2./lx;

	sfx(1) = f*(-3.0/4.0+3.0/4.0*xb*xb);
	sfx(2) = f*(lx*(-1.0-2.0*xb+3.0*xb*xb)/8.0);
	sfx(3) = 0.;
	sfx(4) = 0.;
	sfx(5) = f*(3.0/4.0-3.0/4.0*xb*xb);
	sfx(6) = f*(lx*(1.0+xb)*(-1.0+xb)/4.0+lx*(1.0+xb)*(1.0+xb)/8.0);
	sfx(7) = 0.;
	sfx(8) = 0.;

}

void ANCFBeam3D::GetS0y(Vector& sfx, const Vector3D& ploc) const
{
	double xb = ploc.X();
	sfx.SetLen(NS());
	double f = 2./ly;
	sfx(1) = 0.0;
	sfx(2) = 0.0;
	sfx(3) = -f*ly*(-1.0+xb)/4.0;
	sfx(4) = 0.0;
	sfx(5) = 0.0;
	sfx(6) = 0.0;
	sfx(7) = f*(1.0+xb)*ly/4.0;
	sfx(8) = 0.0;
}

void ANCFBeam3D::GetS0z(Vector& sfx, const Vector3D& ploc) const
{
	double xb = ploc.X();
	sfx.SetLen(NS());
	double f = 2./lz;

	sfx(1) = 0.0;
	sfx(2) = 0.0;
	sfx(3) = 0.0;
	sfx(4) = -f*lz*(-1.0+xb)/4.0;
	sfx(5) = 0.0;
	sfx(6) = 0.0;
	sfx(7) = 0.0;
	sfx(8) = f*(1.0+xb)*lz/4.0;
}

void ANCFBeam3D::GetS0yx(Vector& sfx, const Vector3D& ploc) const
{
	double xb = ploc.X();
	sfx.SetLen(NS());
	double f = 4./(ly*lx);
	sfx(1) = 0.0;
	sfx(2) = 0.0;
	sfx(3) = -f*ly/4.0;
	sfx(4) = 0.0;
	sfx(5) = 0.0;
	sfx(6) = 0.0;
	sfx(7) = f*ly/4.0;
	sfx(8) = 0.0;
}

void ANCFBeam3D::GetS0zx(Vector& sfx, const Vector3D& ploc) const
{
	double xb = ploc.X();
	sfx.SetLen(NS());
	double f = 4./(lz*lx);
	sfx(1) = 0.0;
	sfx(2) = 0.0;
	sfx(3) = 0.0;
	sfx(4) = -f*lz/4.0;
	sfx(5) = 0.0;
	sfx(6) = 0.0;
	sfx(7) = 0.0;
	sfx(8) = f*lz/4.0;
}

void ANCFBeam3D::GetS0xx(Vector& sfxx, const Vector3D& ploc) const
{
	double xb = ploc.X();
	sfxx.SetLen(NS());
	double f = 4./Sqr(lx);
	sfxx(1) = f*3.0/2.0*xb;
	sfxx(2) = f*(-1.0+3.0*xb)*lx/4.0;
	sfxx(3) = 0.0;
	sfxx(4) = 0.0;
	sfxx(5) = -f*3.0/2.0*xb;
	sfxx(6) = f*(1.0+3.0*xb)*lx/4.0;
	sfxx(7) = 0.0;
	sfxx(8) = 0.0;

}




void ANCFBeam3D::GetRot(Matrix3D& rot, const Vector& xg) const
{

	Vector3D r1(xg(13)-xg(1),xg(14)-xg(2),xg(15)-xg(3)); //r_B - r_A
	Vector3D r2(0.5*(xg(12+7)+xg(7)),0.5*(xg(12+8)+xg(8)),0.5*(xg(12+9)+xg(9))); //0.5*(r,y_B + r,y_A)
	Vector3D r3;

	r1.Normalize();
	double h = r2*r1;
	r2 -= h*r1;
	r2.Normalize();
	r3 = r1.Cross(r2);

	rot.Set(r1,r2,r3);
}

Vector3D ANCFBeam3D::GetPos(const Vector3D& p_loc) const
{
	static Vector xg;
	xg.SetLen(SOS());
	GetCoordinates(xg);
	ApplyT(xg);
	Vector3D p0=p_loc;
	p0.Scale(0.5*lx,0.5*ly,0.5*lz);
	static Vector SV;
	GetS0(SV, p0);
	Vector3D p(0.,0.,0.);
	for (int i = 1; i <= 3; i++)
	{
		for (int j = 1; j <= 8; j++)
		{
			p(i) += SV(j)*xg((j-1)*3+i);
		}
	}
	return p;
};

Vector3D ANCFBeam3D::GetVel(const Vector3D& p_loc) const
{
	static Vector xg;
	xg.SetLen(SOS());
	GetCoordinatesP(xg);
	ApplyT(xg);
	Vector3D p0=p_loc;
	p0.Scale(0.5*lx,0.5*ly,0.5*lz);
	static Vector SV;
	GetS0(SV, p0);
	Vector3D p(0.,0.,0.);
	for (int i = 1; i <= 3; i++)
	{
		for (int j = 1; j <= 8; j++)
		{
			p(i) += SV(j)*xg((j-1)*3+i);
		}
	}
	return p;
};

Vector3D ANCFBeam3D::GetPosD(const Vector3D& p_loc) const
{
	static Vector xgd;
	xgd.SetLen(SOS());
	GetDrawCoordinates(xgd);
	xgd = T*xgd;
	Vector3D p0=p_loc;
	p0.Scale(0.5*lx,0.5*ly,0.5*lz);
	static Vector SV;
	GetS0(SV, p0);
	Vector3D p(0.,0.,0.);
	for (int i = 1; i <= 3; i++)
	{
		for (int j = 1; j <= 8; j++)
		{
			p(i) += SV(j)*xgd((j-1)*3+i);
		}
	}
	return p;
};

//in reference element coordinates (-1..1)
Vector3D ANCFBeam3D::GetPos0D(const Vector3D& p_loc, double def_scale) const
{
	static Vector xgd;
	xgd.SetLen(SOS());
	GetDrawCoordinates(xgd);
	//xgd = T*xgd;
	ApplyTD(xgd);

	int i;
	for (i=1; i <=SOS(); i++)
	{
		xgd(i) = def_scale*(xgd(i)-x_init(i))+x_init(i);
	}

	static Vector SV;
	GetS0(SV, p_loc);
	Vector3D p(0.,0.,0.);
	for (i = 1; i <= 3; i++)
	{
		for (int j = 1; j <= 8; j++)
		{
			p(i) += SV(j)*xgd((j-1)*3+i);
		}
	}
	return p;
};

Vector3D ANCFBeam3D::GetVelD(const Vector3D& p_loc) const
{
	static Vector xgd;
	xgd.SetLen(SOS());
	GetDrawCoordinatesP(xgd);
	xgd = T*xgd;
	Vector3D p0=p_loc;
	p0.Scale(0.5*lx,0.5*ly,0.5*lz);
	static Vector SV;
	GetS0(SV, p0);
	Vector3D p(0.,0.,0.);
	for (int i = 1; i <= 3; i++)
	{
		for (int j = 1; j <= 8; j++)
		{
			p(i) += SV(j)*xgd((j-1)*3+i);
		}
	}
	return p;
};


Vector3D ANCFBeam3D::GetPosx(const Vector3D& p_loc) const
{
	static Vector xg;
	xg.SetLen(SOS());
	GetCoordinates(xg);
	ApplyT(xg);
	Vector3D p0=p_loc;
	p0.Scale(0.5*lx,0.5*ly,0.5*lz);
	static Vector SV;
	GetS0x(SV, p0);
	Vector3D p(0.,0.,0.);
	for (int i = 1; i <= 3; i++)
	{
		for (int j = 1; j <= 8; j++)
		{
			p(i) += SV(j)*xg((j-1)*3+i);
		}
	}
	return p;
};

Vector3D ANCFBeam3D::GetDisplacement(const Vector3D& p_loc) const
{
	static Vector xg;
	xg.SetLen(SOS());
	GetCoordinates(xg);

	for (int i=1; i <= SOS(); i++)
		xg(i) -= x_init(i);

	ApplyT(xg);
	Vector3D p0=p_loc;
	p0.Scale(0.5*lx,0.5*ly,0.5*lz);
	static Vector SV;
	GetS0(SV, p0);
	Vector3D p(0.,0.,0.);
	for (int i = 1; i <= 3; i++)
	{
		for (int j = 1; j <= 8; j++)
		{
			p(i) += SV(j)*xg((j-1)*3+i);
		}
	}
	return p;
};

Vector3D ANCFBeam3D::GetDisplacementD(const Vector3D& p_loc) const
{
	static Vector xgd;
	xgd.SetLen(SOS());
	GetDrawCoordinates(xgd);

	for (int i=1; i <= SOS(); i++)
		xgd(i) -= x_init(i);

	xgd = T*xgd;
	Vector3D p0=p_loc;
	p0.Scale(0.5*lx,0.5*ly,0.5*lz);
	static Vector SV;
	GetS0(SV, p0);
	Vector3D p(0.,0.,0.);
	for (int i = 1; i <= 3; i++)
	{
		for (int j = 1; j <= 8; j++)
		{
			p(i) += SV(j)*xgd((j-1)*3+i);
		}
	}
	return p;
};



void ANCFBeam3D::GetH(Matrix& H) 
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

		DS.SetSize(3,ns);
		SV.SetLen(8);

		GetIntegrationRule(x1,w1,3); //3x1x1 !!!!!
		GetIntegrationRule(x2,w2,1);
		GetIntegrationRule(x3,w3,1);

		for (int i1=1; i1<=x1.GetLen(); i1++)
		{
			for (int i2=1; i2<=x2.GetLen(); i2++)
			{
				for (int i3=1; i3<=x3.GetLen(); i3++)
				{
					Vector3D p(x1(i1),x2(i2),x3(i3));
					GetS0(SV,p);

					GetDSMatrix0(p,DS);
					GetJacobi(jac,p,DS, e0);
					double jacdet = jac.Det();
					double fact = fabs (jacdet) * w1(i1)*w2(i2)*w3(i3);

					for (int i=0; i<ns; i++)
					{
						for (int j=1; j<=dim; j++)
						{
							H(i*dim+j,j)+=fact*SV(i+1);
						}
					}
				}
			}
		}
		ApplyTtp(H);
		//H=T.GetTp()*H;
		Hmatrix = H;
	}
}

void ANCFBeam3D::EvalM(Matrix& m, double t) 
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
		DS.SetSize(3,ns);

		Vector x1,x2,x3,w1,w2,w3;

		GetIntegrationRule(x1,w1,6); //optimal: 6x3x3, 6x1x1 geht auch!!!!
		GetIntegrationRule(x2,w2,3);
		GetIntegrationRule(x3,w3,3);
		//UO() << "intrule x1=" << x1 << ", w1=" << w1 << "\n";
		//UO() << "intrule x2=" << x2 << ", w2=" << w2 << "\n";

		for (int i1=1; i1<=x1.GetLen(); i1++)
		{
			for (int i2=1; i2<=x2.GetLen(); i2++)
			{
				for (int i3=1; i3<=x3.GetLen(); i3++)
				{
					Vector3D p(x1(i1),x2(i2),x3(i3));
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
					m += fabs (jacdet) * rho * w1(i1)*w2(i2)*w3(i3) * (HL*HL.GetTp());
				}
			}
		}
		m(1,1) += concentratedmass1;
		m(2,2) += concentratedmass1;
		m(3,3) += concentratedmass1;
		m(1+12,1+12) += concentratedmass2;
		m(2+12,2+12) += concentratedmass2;
		m(3+12,3+12) += concentratedmass2;

		massmatrix = T.GetTp()*m*T;
		m = massmatrix;
	}
	//Matrix minv = m;
	//UO() << "Mass matrix invertable=" << minv.Invert2() << "\n";
	//UO() << "m=" << m << "\n";
};

double ANCFBeam3D::GetTau(const Vector& xg) const
{
	//simplified torsion:
	Vector3D ryA(xg(7),xg(8),xg(9));
	Vector3D ryB(xg(12+7),xg(12+8),xg(12+9));
	Vector3D rA(xg(1),xg(2),xg(3));
	Vector3D rB(xg(12+1),xg(12+2),xg(12+3));
	Vector3D rAB = rB-rA;

	//Project into normal plane:
	rAB.GramSchmidt(ryA);
	rAB.GramSchmidt(ryB);

	//torsion
	double ryAryB = (ryA.Cross(ryB)).Norm();
	double val = ryAryB/(ryA.Norm()*ryB.Norm());
	double tau;
	if (fabs(val) >= 1) tau = val;
	//else tau = acos(ryAryB/(ryA.Norm()*ryB.Norm()));
	else tau = asin(val);
	return tau;
}

void ANCFBeam3D::GetDeltaTau(const Vector& xg, Vector& deltatau) const
{
	static Vector xgdiff;
	xgdiff = xg;

	deltatau.SetLen(SOS());
	double tau0 = GetTau(xgdiff);
	double eps = 1e-7;
	for (int i = 1; i <= SOS(); i++)
	{
		xgdiff(i) += eps;
		deltatau(i) = 1./eps*(GetTau(xgdiff) - tau0);
		xgdiff(i) -= eps;
	}
}


void ANCFBeam3D::GetDeltaKappa(const double& x, const Vector& xg, Vector& dkappa, double& kappa) const
{
	int dim = Dim();
	int ns = NS();

	static Vector SVx;
	static Vector SVxx;
	SVx.SetLen(NS());
	SVxx.SetLen(NS());

	GetS0x(SVx,x);
	GetS0xx(SVxx,x);

	Vector3D rx(0.,0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			rx(i) += SVx(j)*xg((j-1)*3+i);
		}
	}

	Vector3D rxx(0.,0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			rxx(i) += SVxx(j)*xg((j-1)*3+i);
		}
	}

	//Vector3D rx = GetPosx0(x,xg);
	//Vector3D rxx = GetPosxx0(x,xg);
	double rxn = rx.Norm();

	if (rxn == 0) {dkappa *= 0; return;}
	Vector3D rxcrxx = rx.Cross(rxx);
	double f = (rxcrxx).Norm();
	double g = Cub(rxn);

	if (rxn == 0) kappa=0;
	else kappa = f/g;

	double g2inv = 1./Sqr(g);
	double fn = f*g2inv*(3*rxn);
	double gn = g*g2inv;
	if (f != 0) {gn/=f;}


	Vector3D t1;
	double df, dg;

	for (int i=1; i <= dim; i++)
	{
		for (int j=1; j <= ns; j++)
		{
			//dr,x/de x r,xx + r,x x dr,xx/de:
			switch (i) {
					case 1: 
						t1.X() = 0;     
						t1.Y() =-(SVx(j)*rxx.Z()-SVxx(j)*rx.Z()); 
						t1.Z() =  SVx(j)*rxx.Y()-SVxx(j)*rx.Y(); break; //dy=dz=0
					case 2: 
						t1.X() =  SVx(j)*rxx.Z()-SVxx(j)*rx.Z(); 
						t1.Y() =  0;     
						t1.Z() =-(SVx(j)*rxx.X()-SVxx(j)*rx.X()); break; //dx=dz=0
					case 3: 
						t1.X() =-(SVx(j)*rxx.Y()-SVxx(j)*rx.Y()); 
						t1.Y() =  SVx(j)*rxx.X()-SVxx(j)*rx.X(); 
						t1.Z() =  0; break; //dx=dy=0
						//case 1: t1.X() = 0;     t1.Y() =-dx*rz; t1.Z() = dx*ry; break; //dy=dz=0
						//case 2: t1.X() = dy*rz; t1.Y() = 0;     t1.Z() =-dy*rx; break; //dx=dz=0
						//case 3: t1.X() =-dz*ry; t1.Y() = dz*rx; t1.Z() = 0; break; //dx=dy=0
					default: ;
			}
			//dg = (3*rxn)*(rx(i)*SVx(j)); 
			dg = (rx(i)*SVx(j)); //normed

			//df = (rxcrxx*t1);
			df = (rxcrxx*t1); //normed

			//dkappa((j-1)*dim+i) = (df*g-f*dg)/Sqr(g);
			dkappa((j-1)*dim+i) = df*gn-fn*dg; //normed
		}
	}
	//oldkappa = dkappa;
}




int ANCFBeam3D::FastStiffnessMatrix() const 
{
	return 0; //1==Position derivatives done analytically, damping added with numerical differentiation
} 


void ANCFBeam3D::StiffnessMatrix(Matrix& m) //fill in sos x sos components, m might be larger
{
	//size of m given!!!!!!! DO NOT CHANGE size, m contains derivatives w.r.t. velocity coordinates, etc.
	int sos= SOS();
	int dim = Dim();
	int ns = NS();

	TMStartTimer(24);

	if (elasticforce_beam == 4 || elasticforce_beam == 0)
	{

		static Vector u;
		u.SetLen(sos);
		int useu = 1;
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

		double la=Em * nu / ((1.+nu)*(1.-2.*nu));
		double mu=Em / 2. / (1.+nu);
		double ks = 10.*(1+nu)/(12+11*nu); //shear correction factor

		Matrix3D strain, piola1, F, stress;
		Matrix3D jac, jacinv;
		double jacdet;

		static Matrix Dmat;
		static Matrix Dmat2;
		Dmat.SetSize(6,6);
		Dmat2.SetSize(6,6);
		ConstVector<6> strainvec(6);
		ConstVector<6> stressvec(6);

		for (int i=1; i <= sos; i++)
		{
			for (int j=1; j <= sos; j++)
			{
				m(i,j) = 0;
			}
		}

		GetDMatrix(Dmat, nu, Em);

		int kkend = 1;
		if (elasticforce_beam == 4)
		{
			kkend = 2;
			Dmat(4,4) *= ks;
			Dmat(5,5) *= ks;
			Dmat(6,6) *= ks;

			Dmat2.SetAll(0);
			Dmat2(1,1) = Em;
			Dmat2(2,2) = Em;
			Dmat2(3,3) = Em;
			Dmat2(4,4) = Dmat(4,4);
			Dmat2(5,5) = Dmat(5,5);
			Dmat2(6,6) = Dmat(6,6);

			Dmat -= Dmat2;
		}

		int reduceorder = 2; //reduce order of stiffness matrix
		Matrix3D A2;
		Matrix3D A1;
		Matrix3D dStress;
		Matrix3D A3;

		for (int kk=1; kk <= kkend; kk++)
		{
			//kk=1: just poisson effect at element line
			//kk=2: no poisson effect, whole body

			if (kk == 2 || kkend == 1)
			{
				GetIntegrationRule(x1,w1,orderx-reduceorder); //x_i ... integration points, w_i ... integration weights
				GetIntegrationRule(x2,w2,orderyz-reduceorder);
				GetIntegrationRule(x3,w3,orderyz-reduceorder);
			}
			else
			{
				GetIntegrationRule(x1,w1,orderx-reduceorder); //x_i ... integration points, w_i ... integration weights
				GetIntegrationRule(x2,w2,1);
				GetIntegrationRule(x3,w3,1);
			}

			int kx1 = x2.Length()*x3.Length();
			int kx2 = x3.Length();

			for (int ip1=1; ip1<=x1.GetLen(); ip1++)
			{
				for (int ip2=1; ip2<=x2.GetLen(); ip2++)
				{
					for (int ip3=1; ip3<=x3.GetLen(); ip3++)
					{
						//TMStartTimer(20);
						int i,j,k;
						int ind = (ip1-1)*kx1+(ip2-1)*kx2+(ip3-1);
						Vector3D p(x1(ip1),x2(ip2),x3(ip3));

						// compute F 
						F.SetAll(0);
						int l;

						static Matrix agrad; //\bar S^D
						GetDSMatrix0(p,DS);
						GetJacobi(jac,p,DS,x_init);
						jacdet = jac.Det();
						jac.GetInverse(jacinv);
						jacinv = jacinv.GetTp();

						agrad.SetSize(3,NS());
						Mult(jacinv, DS, agrad);

						for (j = 1; j <= dim; j++) 
						{
							for (i = 1; i <= ns; i++)
							{
								l = (i-1)*dim+j;
								//u = xg(l) - x_init(l);
								for (k = 1; k <= dim; k++)
								{
									F(j,k) += agrad(k,i)*u(l);
									//G(j,k) += DS(k,i)*u((i-1)*dim+j);
								}
							}
							if (useu) F(j,j) += 1;
						}

						// Green-Lagrange strain tensor
						//strain = 0.5 * (F.GetTp() * F - I);
						F.GetATA2(strain);
						strain(1,1) -= 0.5; strain(2,2) -= 0.5; strain(3,3) -= 0.5;

						strainvec(1) = strain(1,1);
						strainvec(2) = strain(2,2);
						strainvec(3) = strain(3,3);
						strainvec(4) = 2.*strain(2,3);
						strainvec(5) = 2.*strain(3,1);
						strainvec(6) = 2.*strain(1,2);

						if (kk == 1)
							Mult(Dmat,strainvec,stressvec);
						else
							Mult(Dmat2,strainvec,stressvec);

						//piola1 = F * ((2*mu) * strain + Matrix3D(la * strain.Trace()));
						stress(1,1) = stressvec(1);
						stress(2,2) = stressvec(2);
						stress(3,3) = stressvec(3);
						stress(2,3) = stressvec(4);
						stress(3,1) = stressvec(5);
						stress(1,2) = stressvec(6);
						stress(3,2) = stressvec(4);
						stress(1,3) = stressvec(5);
						stress(2,1) = stressvec(6);

						//loop for all components of stiffness matrix:
						for (int i1=1; i1 <= ns; i1++)
						{
							for (int j1=1; j1 <= dim; j1++)
							{
								for (int i2=1; i2 <= ns; i2++)
								{
									for (int j2=1; j2 <= dim; j2++)
									{
										double dKij = 0;
										//dKij = dS/dqj : F^T dF/dqi + S : dF^T/dqj dF/dqi
										//     = 0.5*(F^T dF/dqj + dF/dqj^T F) : D : F^T dF/dqi + S : dF^T/dqj dF/dqi
										//     = 0.5*(A1         + A1^T      ) : D : A2         + S : A3

										//A2(r,s) = F(r,t)^T dF(t,s)/dq(i) = F(r,t)^T dF(t,s)/dq((i1-1)*3+j1) = F(j1,r)*agrad(s,i1)
										//A1(r,s) = F(r,t)^T dF(t,s)/dq(j) = F(r,t)^T dF(t,s)/dq((i2-1)*3+j2) = F(j2,r)*agrad(s,i2)
										for (int r=1; r <= dim; r++)
											for (int s=1; s <= dim; s++)
											{
												A2(r,s) = F(j1, r) * agrad(s,i1);
												A1(r,s) = F(j2, r) * agrad(s,i2);
											}

											A1.MakeSym();

											//dStress = D:A1
											/*
											double laA1tr = la * A1.Trace();
											A1 *= 2*mu;
											A1(1,1) += laA1tr;
											A1(2,2) += laA1tr;
											A1(3,3) += laA1tr;
											dKij += A2.DoubleContract(A1);
											*/

											strainvec(1) = A1(1,1);
											strainvec(2) = A1(2,2);
											strainvec(3) = A1(3,3);
											strainvec(4) = 2.*A1(2,3);
											strainvec(5) = 2.*A1(3,1);
											strainvec(6) = 2.*A1(1,2);

											if (kk == 1)
											{
												//Mult(Dmat,strainvec,stressvec);

												stressvec(1) = Dmat(1,1)*strainvec(1) + Dmat(1,2)*strainvec(2) + Dmat(1,3)*strainvec(3);
												stressvec(2) = Dmat(2,1)*strainvec(1) + Dmat(2,2)*strainvec(2) + Dmat(2,3)*strainvec(3);
												stressvec(3) = Dmat(3,1)*strainvec(1) + Dmat(3,2)*strainvec(2) + Dmat(3,3)*strainvec(3);
												stressvec(4) = Dmat(4,4)*strainvec(4);
												stressvec(5) = Dmat(5,5)*strainvec(5);
												stressvec(6) = Dmat(6,6)*strainvec(6);

											}
											else
											{
												//Mult(Dmat2,strainvec,stressvec);

												stressvec(1) = Dmat2(1,1)*strainvec(1);
												stressvec(2) = Dmat2(2,2)*strainvec(2);
												stressvec(3) = Dmat2(3,3)*strainvec(3);
												stressvec(4) = Dmat2(4,4)*strainvec(4);
												stressvec(5) = Dmat2(5,5)*strainvec(5);
												stressvec(6) = Dmat2(6,6)*strainvec(6);
											}

											//piola1 = F * ((2*mu) * strain + Matrix3D(la * strain.Trace()));
											dStress(1,1) = stressvec(1);
											dStress(2,2) = stressvec(2);
											dStress(3,3) = stressvec(3);
											dStress(2,3) = stressvec(4);
											dStress(3,1) = stressvec(5);
											dStress(1,2) = stressvec(6);
											dStress(3,2) = stressvec(4);
											dStress(1,3) = stressvec(5);
											dStress(2,1) = stressvec(6);
											//dKij += A2.DoubleContract(dStress);
											//* YV 17.11.2010: the previous line was commented out as the function Matrix3D::DoubleContract
											// used to have a strange meaning, and now it is unclear what should appear there:
											// either dStress : dStress or A2 : dStress.
											GetMBS()->UO() << "Warning: reached a point in ANCFBeam3D::StiffnessMatrix() in which a decision is yet to be made, the results may be incorrect.";




											//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
											//geometric stiffening:
											//A3(r,s) = dF(k,r)/dqj dF(k,s)/dqi = {DS(r,i1) * DS(s,i2) if (j1==j2), 0 else}
											if (j1 == j2)
											{
												//compute this once for all j1 and j2 !!!!
												/*
												for (int r=1; r <= dim; r++)
												for (int s=1; s <= dim; s++)
												A3(r,s) =  agrad(r,i1) * agrad(s,i2);
												dKij += stress.DoubleContract(A3);
												*/
												for (int r=1; r <= dim; r++)
													for (int s=1; s <= dim; s++)
													{
														dKij += agrad(r,i1) * agrad(s,i2) * stress(r,s);
													}
											}

											m((i1-1)*3+j1,(i2-1)*3+j2) -= dKij * fabs (jacdet) * w1(ip1)*w2(ip2)*w3(ip3); //negative, because K is put on right-hand-side
									}
								}
							}
						}
					}
				}
			}
		}
	}


	TMStopTimer(24);


	if (0)
	{
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//Compute initial stiffness matrix K0 at first step!!!!!
		//u = u_0 + R*(X+u_f^*) - X
		//R is rotation from initial position to actual position, R=Ortho([r_x r_y r_z])
		//K_T = R*K0*R^T
		if (K0.Getrows() != sos) 
		{
			K0.SetSize(sos,sos);
			static Vector locjacf20; 
			locjacf20.SetLen(sos);
			locjacf20.SetAll(0);
			static Vector locjacf21; 
			locjacf21.SetLen(sos);

			double numdiffepsi = GetMBS()->NumSolver().NumDiffepsi();
			double eps;

			double t = GetMBS()->GetTime();
			EvalF2(locjacf20,t);
			double storex;


			for (int i = 1; i <= sos; i++)
			{
				eps = numdiffepsi*Maximum(1e-2,fabs(XG(i)));

				storex = XG(i);
				XG(i) += 2*eps;
				locjacf21.SetAll(0);
				EvalF2(locjacf21,t);
				XG(i) = storex;

				for (int j=1; j<=sos;j++)
				{
					K0(j,i)= (locjacf21(j)-locjacf20(j))/(2.*eps);
				}
			}
			//UO() << "K0=" << K0 << "\n";
		}
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		SetComputeCoordinates();
		ApplyT(xg); //xg is transformed version of coordinates

		//compute R 3x3

		//???try if two different rotation matrices at each node are better???
		int mode = 2; //1= single rotation matrix, 2=two rotation matrices

		Matrix3D R1, R1T, R2, R2T;
		if (mode == 1)
		{
			//average of 2 gradients:
			Vector3D x1(0.5*(xg(4)), 0.5*(xg(5)), 0.5*(xg(6)));
			Vector3D y1(0.5*(xg(7)), 0.5*(xg(8)), 0.5*(xg(9)));
			Vector3D z1(0.5*(xg(10)),0.5*(xg(11)),0.5*(xg(12)));
			//Vector3D x1(0.5*(xg(4)+ xg(4+12)), 0.5*(xg(5)+ xg(5+12)), 0.5*(xg(6)+ xg(6+12)));
			//Vector3D y1(0.5*(xg(7)+ xg(7+12)), 0.5*(xg(8)+ xg(8+12)), 0.5*(xg(9)+ xg(9+12)));
			//Vector3D z1(0.5*(xg(10)+xg(10+12)),0.5*(xg(11)+xg(11+12)),0.5*(xg(12)+xg(12+12)));

			x1.Normalize();
			Vector3D y1a = y1;
			Vector3D z1a = z1;
			x1.GramSchmidt(y1a);
			x1.GramSchmidt(z1);
			y1 = -1*x1.Cross(z1);
			y1 = 0.5*(y1+y1a);
			x1.GramSchmidt(y1);
			z1 = x1.Cross(y1);

			//??? Rotationsmatrix stimmt nicht!!!

			R1 = Matrix3D(x1.X(),y1.X(),z1.X(),
				x1.Y(),y1.Y(),z1.Y(),
				x1.Z(),y1.Z(),z1.Z());
			R1T = R1.GetTp();
		}
		else
		{
			//first node:
			Vector3D x1(0.5*(xg(4)), 0.5*(xg(5)), 0.5*(xg(6)));
			Vector3D y1(0.5*(xg(7)), 0.5*(xg(8)), 0.5*(xg(9)));
			Vector3D z1(0.5*(xg(10)),0.5*(xg(11)),0.5*(xg(12)));

			x1.Normalize();
			Vector3D y1a = y1;
			Vector3D z1a = z1;
			x1.GramSchmidt(y1a);
			x1.GramSchmidt(z1);
			y1 = -1.*x1.Cross(z1);
			y1 = 0.5*(y1+y1a);
			x1.GramSchmidt(y1);
			z1 = x1.Cross(y1);

			//??? Rotationsmatrix stimmt nicht!!!

			R1 = Matrix3D(x1.X(),y1.X(),z1.X(),
				x1.Y(),y1.Y(),z1.Y(),
				x1.Z(),y1.Z(),z1.Z());
			R1T = R1.GetTp();

			//second node:
			x1 = Vector3D(0.5*(xg(4+12)), 0.5*(xg(5+12)), 0.5*(xg(6+12)));
			y1 = Vector3D(0.5*(xg(7+12)), 0.5*(xg(8+12)), 0.5*(xg(9+12)));
			z1 = Vector3D(0.5*(xg(10+12)),0.5*(xg(11+12)),0.5*(xg(12+12)));
			x1.Normalize();
			y1a = y1;
			z1a = z1;
			x1.GramSchmidt(y1a);
			x1.GramSchmidt(z1);
			y1 = -1.*x1.Cross(z1);
			y1 = 0.5*(y1+y1a);
			x1.GramSchmidt(y1);
			z1 = x1.Cross(y1);

			//??? Rotationsmatrix stimmt nicht!!!

			R2 = Matrix3D(x1.X(),y1.X(),z1.X(),
				x1.Y(),y1.Y(),z1.Y(),
				x1.Z(),y1.Z(),z1.Z());
			R2T = R2.GetTp();
		}

		//UO() << "R=" << R << "\n";
		//UO() << "R1=" << R1 << "\n";
		//UO() << "RT=" << RT << "\n";

		//write K0 in m
		//transform = R*K0(3x3)*R^T
		for (int i=0; i<sos/3; i++)
		{
			for (int j=0; j<sos/3; j++)
			{
				int io = i*3;
				int jo = j*3;

				Matrix3D Ksub;
				for (int i0=1; i0<=3; i0++)
				{
					for (int j0=1; j0<=3; j0++)
					{
						Ksub(i0,j0) = K0(io+i0, jo+j0);
					}
				}
				if (mode == 1)
				{
					Ksub = (R1*Ksub)*R1T;
				}
				else
				{
					//Ksub = (R2*Ksub)*R2T;

					if (i < sos/6 && j < sos/6)
					{
						Ksub = (R1*Ksub)*R1T;
					}
					else if (i >=sos/6 && j < sos/6)
					{
						Ksub = (R2*Ksub)*R1T;
					}
					else if (i < sos/6 && j >=sos/6)
					{
						Ksub = (R1*Ksub)*R2T;
					}
					else if (i >=sos/6 && j >=sos/6)
					{
						Ksub = (R2*Ksub)*R2T;
					}
				}
				for (int i0=1; i0<=3; i0++)
				{
					for (int j0=1; j0<=3; j0++)
					{
						m(io+i0, jo+j0) = Ksub(i0,j0);
					}
				}
			}
		}
		//UO() << "RKR=" << m << "\n";
	}
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//compute curvature, stretch, shear and torsion for initial values in order to make precurved beams possible
void ANCFBeam3D::SetInitialCurvatures()
{
	if (elasticforce_beam != 3) return;

	int sos = SOS();
	//SetComputeCoordinates(); //only for EvalF2
	//--> here: set xg = x_init;
	xg = x_init; 

	int ns = NS();
	int dim = 3;

	ApplyT(xg);


	Vector Sx; Sx.SetLen(NS()); Sx.SetAll(0);  //rx
	Vector Sy; Sy.SetLen(NS()); Sy.SetAll(0);   //ry
	Vector Sz;	Sz.SetLen(NS()); Sz.SetAll(0);   //rz
	Vector Sxx; Sxx.SetLen(NS()); Sxx.SetAll(0); //rxx
	Vector Syx; Syx.SetLen(NS()); Syx.SetAll(0); //ryx
	Vector Szx; Szx.SetLen(NS()); Szx.SetAll(0); //rzx


	Vector3D rx, ry, rz, rxx, ryx, rzx;
	//integration of elongation and cross-section change:
	GetIntegrationRuleLobatto(x1,w1,orderWl);
	for (int i1=1; i1<=x1.GetLen(); i1++)
	{
		double x = x1(i1); 
		double x2 = x*x;

		//rx:
		//GetS0x(Sx,Vector3D(x,0,0));

		Sx(1) = (-1.5+1.5*x2)/lx;
		Sx(2) = (-0.25-0.5*x+0.75*x2);
		Sx(5) = -Sx(1);//(1.5-1.5*x2)/lx;
		Sx(6) = Sx(2)+x;//(-0.25+0.5*x+0.75*x2);

		rx.X()  = Sx(1)*xg(3*0+1);
		rx.Y()  = Sx(1)*xg(3*0+2);
		rx.Z()  = Sx(1)*xg(3*0+3);					
		rx.X() += Sx(2)*xg(3*1+1);
		rx.Y() += Sx(2)*xg(3*1+2);
		rx.Z() += Sx(2)*xg(3*1+3);					
		rx.X() += Sx(5)*xg(3*4+1);
		rx.Y() += Sx(5)*xg(3*4+2);
		rx.Z() += Sx(5)*xg(3*4+3);					
		rx.X() += Sx(6)*xg(3*5+1);
		rx.Y() += Sx(6)*xg(3*5+2);
		rx.Z() += Sx(6)*xg(3*5+3);					

		//GetS0y(Sy,Vector3D(x,0,0));
		Sy(3) = -0.5*(-1.0+x);
		Sy(7) = 0.5*(1.0+x);
		ry.X()  = Sy(3)*xg(3*2+1);
		ry.Y()  = Sy(3)*xg(3*2+2);
		ry.Z()  = Sy(3)*xg(3*2+3);					
		ry.X() += Sy(7)*xg(3*6+1);
		ry.Y() += Sy(7)*xg(3*6+2);
		ry.Z() += Sy(7)*xg(3*6+3);					

		//GetS0z(Sz,Vector3D(x,0,0));
		Sz(4) = -0.5*(-1.0+x);
		Sz(8) = 0.5*(1.0+x);
		rz.X()  = Sz(4)*xg(3*3+1);
		rz.Y()  = Sz(4)*xg(3*3+2);
		rz.Z()  = Sz(4)*xg(3*3+3);					
		rz.X() += Sz(8)*xg(3*7+1);
		rz.Y() += Sz(8)*xg(3*7+2);
		rz.Z() += Sz(8)*xg(3*7+3);					


		//integration factors:
		double fact = lx * 0.5 * w1(i1);

		Wlepsx0(i1) = 0.5*(rx*rx-1.);
		Wlepsy0(i1) = 0.5*(ry*ry-1.)*factstiffWl;
		Wlepsz0(i1) = 0.5*(rz*rz-1.)*factstiffWl;
		Wlgamyz0(i1)= ry*rz*factstiffWl ;

	}


	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//Hellinger Reissner principle:
	//note that xi=0..1 in the paper of Schwab and x = -1..+1 

	//1. compute Wxz1, Wxz2, Wxy1, Wxy2
	//2. compute delta Wx... times Wx... and sum up in temp
	//UO() << "Hellinger Reisner\n";

	double Wxz1f = 0; 
	double Wxz2f = 0;
	double Wxy1f = 0;
	double Wxy2f = 0;
	//shear deformation:
	GetIntegrationRuleLobatto(x1,w1,orderWsHR); 
	for (int i1=1; i1<=x1.GetLen(); i1++)
	{
		double x = x1(i1);

		double x2 = x*x;
		//rx:
		//GetS0x(Sx,Vector3D(x,0,0));
		Sx(1) = (-1.5+1.5*x2)/lx;
		Sx(2) = (-0.25-0.5*x+0.75*x2);
		Sx(5) = -Sx(1);//(1.5-1.5*x2)/lx;
		Sx(6) = Sx(2)+x;//(-0.25+0.5*x+0.75*x2);
		rx.X()  = Sx(1)*xg(3*0+1);
		rx.Y()  = Sx(1)*xg(3*0+2);
		rx.Z()  = Sx(1)*xg(3*0+3);					
		rx.X() += Sx(2)*xg(3*1+1);
		rx.Y() += Sx(2)*xg(3*1+2);
		rx.Z() += Sx(2)*xg(3*1+3);					
		rx.X() += Sx(5)*xg(3*4+1);
		rx.Y() += Sx(5)*xg(3*4+2);
		rx.Z() += Sx(5)*xg(3*4+3);					
		rx.X() += Sx(6)*xg(3*5+1);
		rx.Y() += Sx(6)*xg(3*5+2);
		rx.Z() += Sx(6)*xg(3*5+3);					

		//GetS0y(Sy,Vector3D(x,0,0));
		Sy(3) = -0.5*(-1.0+x);
		Sy(7) = 0.5*(1.0+x);
		ry.X()  = Sy(3)*xg(3*2+1);
		ry.Y()  = Sy(3)*xg(3*2+2);
		ry.Z()  = Sy(3)*xg(3*2+3);					
		ry.X() += Sy(7)*xg(3*6+1);
		ry.Y() += Sy(7)*xg(3*6+2);
		ry.Z() += Sy(7)*xg(3*6+3);					

		//GetS0z(Sz,Vector3D(x,0,0));
		Sz(4) = -0.5*(-1.0+x);
		Sz(8) = 0.5*(1.0+x);
		rz.X()  = Sz(4)*xg(3*3+1);
		rz.Y()  = Sz(4)*xg(3*3+2);
		rz.Z()  = Sz(4)*xg(3*3+3);					
		rz.X() += Sz(8)*xg(3*7+1);
		rz.Y() += Sz(8)*xg(3*7+2);
		rz.Z() += Sz(8)*xg(3*7+3);					


		//integration factors:
		double fact = lx * 0.5 * w1(i1);

		double gamxy= fact*rx*ry;
		double gamxz= fact*rx*rz;
		//UO() << "rx=" << rx << "\n";
		//UO() << "ry=" << ry << "\n";
		//UO() << "gamxy=" << gamxy << "\n";
		//UO() << "gamxz=" << gamxz << "\n";

		Wxy1f += 0.5*(1-x) * gamxy;
		Wxy2f += 0.5*(1+x) * gamxy;
		Wxz1f += 0.5*(1-x) * gamxz;
		Wxz2f += 0.5*(1+x) * gamxz;

	}

	GetIntegrationRuleLobatto(x1,w1,orderWsHR);
	//now compute variations:
	for (int i1=1; i1<=x1.GetLen(); i1++)
	{
		double x = x1(i1);

		double x2 = x*x;
		//rx:
		//GetS0x(Sx,Vector3D(x,0,0));
		Sx(1) = (-1.5+1.5*x2)/lx;
		Sx(2) = (-0.25-0.5*x+0.75*x2);
		Sx(5) = -Sx(1);//(1.5-1.5*x2)/lx;
		Sx(6) = Sx(2)+x;//(-0.25+0.5*x+0.75*x2);
		rx.X()  = Sx(1)*xg(3*0+1);
		rx.Y()  = Sx(1)*xg(3*0+2);
		rx.Z()  = Sx(1)*xg(3*0+3);					
		rx.X() += Sx(2)*xg(3*1+1);
		rx.Y() += Sx(2)*xg(3*1+2);
		rx.Z() += Sx(2)*xg(3*1+3);					
		rx.X() += Sx(5)*xg(3*4+1);
		rx.Y() += Sx(5)*xg(3*4+2);
		rx.Z() += Sx(5)*xg(3*4+3);					
		rx.X() += Sx(6)*xg(3*5+1);
		rx.Y() += Sx(6)*xg(3*5+2);
		rx.Z() += Sx(6)*xg(3*5+3);					

		//GetS0y(Sy,Vector3D(x,0,0));
		Sy(3) = -0.5*(-1.0+x);
		Sy(7) = 0.5*(1.0+x);
		ry.X()  = Sy(3)*xg(3*2+1);
		ry.Y()  = Sy(3)*xg(3*2+2);
		ry.Z()  = Sy(3)*xg(3*2+3);					
		ry.X() += Sy(7)*xg(3*6+1);
		ry.Y() += Sy(7)*xg(3*6+2);
		ry.Z() += Sy(7)*xg(3*6+3);					

		//GetS0z(Sz,Vector3D(x,0,0));
		Sz(4) = -0.5*(-1.0+x);
		Sz(8) = 0.5*(1.0+x);
		rz.X()  = Sz(4)*xg(3*3+1);
		rz.Y()  = Sz(4)*xg(3*3+2);
		rz.Z()  = Sz(4)*xg(3*3+3);					
		rz.X() += Sz(8)*xg(3*7+1);
		rz.Y() += Sz(8)*xg(3*7+2);
		rz.Z() += Sz(8)*xg(3*7+3);					


		//integration factors:
		double fact = 0.5 * w1(i1);

		//factor for delta gamxy components:
		WsHRk10(i1) = beamGAky * fact*(0.5*(1-x)*(4.*Wxy1f - 2.*Wxy2f) + 0.5*(1+x)*(-2.*Wxy1f + 4.*Wxy2f));
		//factor for delta gamxz components:
		WsHRk20(i1) = beamGAkz * fact*(0.5*(1-x)*(4.*Wxz1f - 2.*Wxz2f) + 0.5*(1+x)*(-2.*Wxz1f + 4.*Wxz2f));

	}


	//torsion:
	GetIntegrationRule(x1,w1,orderWt);
	temp.SetAll(0);
	for (int i1=1; i1<=x1.GetLen(); i1++)
	{
		double x = x1(i1);
		//double x2 = x*x;
		//GetS0y(Sy,Vector3D(x,0,0));
		Sy(3) = -0.5*(-1.0+x);
		Sy(7) = 0.5*(1.0+x);
		ry.X()  = Sy(3)*xg(3*2+1);
		ry.Y()  = Sy(3)*xg(3*2+2);
		ry.Z()  = Sy(3)*xg(3*2+3);					
		ry.X() += Sy(7)*xg(3*6+1);
		ry.Y() += Sy(7)*xg(3*6+2);
		ry.Z() += Sy(7)*xg(3*6+3);					

		//GetS0z(Sz,Vector3D(x,0,0));
		Sz(4) = -0.5*(-1.0+x);
		Sz(8) = 0.5*(1.0+x);
		rz.X()  = Sz(4)*xg(3*3+1);
		rz.Y()  = Sz(4)*xg(3*3+2);
		rz.Z()  = Sz(4)*xg(3*3+3);					
		rz.X() += Sz(8)*xg(3*7+1);
		rz.Y() += Sz(8)*xg(3*7+2);
		rz.Z() += Sz(8)*xg(3*7+3);					

		//ryx:
		//GetS0yx(Syx,Vector3D(x,0,0));
		Syx(3) = -1./lx;
		Syx(7) = 1./lx;
		ryx.X()  = Syx(3)*xg(3*2+1);
		ryx.Y()  = Syx(3)*xg(3*2+2);
		ryx.Z()  = Syx(3)*xg(3*2+3);					
		ryx.X() += Syx(7)*xg(3*6+1);
		ryx.Y() += Syx(7)*xg(3*6+2);
		ryx.Z() += Syx(7)*xg(3*6+3);					

		//rzx:
		//GetS0zx(Szx,Vector3D(x,0,0));
		Szx(4) = -1./lx;
		Szx(8) = 1./lx;
		rzx.X()  = Szx(4)*xg(3*3+1);
		rzx.Y()  = Szx(4)*xg(3*3+2);
		rzx.Z()  = Szx(4)*xg(3*3+3);					
		rzx.X() += Szx(8)*xg(3*7+1);
		rzx.Y() += Szx(8)*xg(3*7+2);
		rzx.Z() += Szx(8)*xg(3*7+3);					


		//integration factors:
		double fact = beamGJkx * lx * 0.5 * w1(i1);

		Wtkapx0(i1) = 0.5 * (rz*ryx - ry*rzx);
	}

	//bending:
	GetIntegrationRule(x1,w1,orderWb);
	for (int i1=1; i1<=x1.GetLen(); i1++)
	{
		double x = x1(i1);

		double x2 = x*x;
		//GetS0y(Sy,Vector3D(x,0,0));
		Sy(3) = -0.5*(-1.0+x);
		Sy(7) = 0.5*(1.0+x);
		ry.X()  = Sy(3)*xg(3*2+1);
		ry.Y()  = Sy(3)*xg(3*2+2);
		ry.Z()  = Sy(3)*xg(3*2+3);					
		ry.X() += Sy(7)*xg(3*6+1);
		ry.Y() += Sy(7)*xg(3*6+2);
		ry.Z() += Sy(7)*xg(3*6+3);					

		//GetS0z(Sz,Vector3D(x,0,0));
		Sz(4) = -0.5*(-1.0+x);
		Sz(8) = 0.5*(1.0+x);
		rz.X()  = Sz(4)*xg(3*3+1);
		rz.Y()  = Sz(4)*xg(3*3+2);
		rz.Z()  = Sz(4)*xg(3*3+3);					
		rz.X() += Sz(8)*xg(3*7+1);
		rz.Y() += Sz(8)*xg(3*7+2);
		rz.Z() += Sz(8)*xg(3*7+3);					

		//rxx:
		//GetS0xx(Sxx,Vector3D(x,0,0));
		Sxx(1) =  4./Sqr(lx)*3.0/2.0*x;
		Sxx(2) =  1./lx*(-1.0+3.0*x);
		Sxx(5) = -Sxx(1);
		Sxx(6) =  1./lx*(1.0+3.0*x);
		rxx.X()  = Sxx(1)*xg(3*(1-1)+1);
		rxx.Y()  = Sxx(1)*xg(3*(1-1)+2);
		rxx.Z()  = Sxx(1)*xg(3*(1-1)+3);					
		rxx.X() += Sxx(2)*xg(3*(2-1)+1);
		rxx.Y() += Sxx(2)*xg(3*(2-1)+2);
		rxx.Z() += Sxx(2)*xg(3*(2-1)+3);					
		rxx.X() += Sxx(5)*xg(3*(5-1)+1);
		rxx.Y() += Sxx(5)*xg(3*(5-1)+2);
		rxx.Z() += Sxx(5)*xg(3*(5-1)+3);					
		rxx.X() += Sxx(6)*xg(3*(6-1)+1);
		rxx.Y() += Sxx(6)*xg(3*(6-1)+2);
		rxx.Z() += Sxx(6)*xg(3*(6-1)+3);					


		//integration factors:
		double facty = beamEIy * lx * 0.5 * w1(i1);
		double factz = beamEIz * lx * 0.5 * w1(i1);

		Wbkapy0(i1) =  rz*rxx;
		Wbkapz0(i1) =-(ry*rxx);
	}


}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ANCFBeam3D::EvalF2(Vector& f, double t) 
{
	Body3D::EvalF2(f,t);

	int sos = SOS();
	SetComputeCoordinates();
	static Vector fadd;

	int ns = NS();
	int dim = 3;
	//double u;

	if (!((elasticforce_beam == 4 || elasticforce_beam == 0) && GetMBS()->IsJacobianComputation() && FastStiffnessMatrix() != 0 && GetMBS()->NumSolver().UseSparseSolver()))
	{
		TMStartTimer(22);

		//0 is Original ANCF
		//3 Arend Schwab && Jaap Meijaard
		//4 Poisson locking eliminated in original ANCF
		//5 Corotational

		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		if (elasticforce_beam == 5)
		{
			int includedeltaterms = 1;
			//elastic line linearized
			mystatic Vector temp2;
			temp2.SetLen(SOS());

			temp.SetLen(SOS());
			fadd.SetLen(SOS());
			fadd.SetAll(0);
			//optimized version:

			ApplyT(xg);

			Matrix3D A, AT;
			Vector3D pref;
			ComputeCorotationalFrame(pref, A);
			AT = A.GetTp();


			mystatic Vector Sx; Sx.SetLen(NS()); Sx.SetAll(0);  //rx
			mystatic Vector Sy; Sy.SetLen(NS()); Sy.SetAll(0);   //ry
			mystatic Vector Sz;	Sz.SetLen(NS()); Sz.SetAll(0);   //rz
			mystatic Vector Sxx; Sxx.SetLen(NS()); Sxx.SetAll(0); //rxx
			mystatic Vector Syx; Syx.SetLen(NS()); Syx.SetAll(0); //ryx
			mystatic Vector Szx; Szx.SetLen(NS()); Szx.SetAll(0); //rzx

			Vector3D rx, ry, rz, rxx, ryx, rzx;
			double /*epsx,*/ epsy, epsz, gamyz, gamxy, gamxz;
			double depsy, depsz;

			Vector3D tx[2]; //rx vectors at both ends
			Vector3D ty[2]; //ry vectors at both ends
			Vector3D tz[2]; //rz vectors at both ends
			
			Vector3D txlin[2]; //rx linearized
			Vector3D tylin[2]; //rx linearized
			Vector3D tzlin[2]; //rx linearized
			txlin[0] = 0.;
			txlin[1] = 0.;
			tylin[0] = 0.;
			tylin[1] = 0.;
			tzlin[0] = 0.;
			tzlin[1] = 0.;

			ConstMatrix<72> temp3(3,SOS());
			ConstVector<8> txlin_delta[2];
			txlin_delta[0].SetLen(8);
			txlin_delta[1].SetLen(8);
			txlin_delta[0].SetAll(0.);
			txlin_delta[1].SetAll(0.);


			Vector3D rxx_int(0.);
			ConstVector<8> rxx_int_delta(8);
			rxx_int_delta.SetAll(0.);

			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			//apply penalty terms for thickness deformation and cross section distorsion at both sides
			//shear deformation included with reduced integration (Lobatto)
			double ix;
			int ind = 0;
			temp.SetAll(0.);
			for (ix = -1; ix <= 1; ix += 2)
			{
				double x = ix;
				double x2 = x*x;

				//GetS0y(Sy,Vector3D(x,0,0));
				Sy(3) = -0.5*(-1.0+x);
				Sy(7) = 0.5*(1.0+x);
				ry.X()  = Sy(3)*xg(3*2+1);
				ry.Y()  = Sy(3)*xg(3*2+2);
				ry.Z()  = Sy(3)*xg(3*2+3);					
				ry.X() += Sy(7)*xg(3*6+1);
				ry.Y() += Sy(7)*xg(3*6+2);
				ry.Z() += Sy(7)*xg(3*6+3);					

				//GetS0z(Sz,Vector3D(x,0,0));
				Sz(4) = -0.5*(-1.0+x);
				Sz(8) = 0.5*(1.0+x);
				rz.X()  = Sz(4)*xg(3*3+1);
				rz.Y()  = Sz(4)*xg(3*3+2);
				rz.Z()  = Sz(4)*xg(3*3+3);					
				rz.X() += Sz(8)*xg(3*7+1);
				rz.Y() += Sz(8)*xg(3*7+2);
				rz.Z() += Sz(8)*xg(3*7+3);					

				Sx(1) = (-1.5+1.5*x2)/lx;
				Sx(2) = (-0.25-0.5*x+0.75*x2);
				Sx(5) = -Sx(1);//(1.5-1.5*x2)/lx;
				Sx(6) = Sx(2)+x;//(-0.25+0.5*x+0.75*x2);

				rx.X()  = Sx(1)*xg(3*0+1);
				rx.Y()  = Sx(1)*xg(3*0+2);
				rx.Z()  = Sx(1)*xg(3*0+3);					
				rx.X() += Sx(2)*xg(3*1+1);
				rx.Y() += Sx(2)*xg(3*1+2);
				rx.Z() += Sx(2)*xg(3*1+3);					
				rx.X() += Sx(5)*xg(3*4+1);
				rx.Y() += Sx(5)*xg(3*4+2);
				rx.Z() += Sx(5)*xg(3*4+3);					
				rx.X() += Sx(6)*xg(3*5+1);
				rx.Y() += Sx(6)*xg(3*5+2);
				rx.Z() += Sx(6)*xg(3*5+3);					


				tx[ind  ] = rx;
				ty[ind  ] = ry;
				tz[ind++] = rz;

				epsy = 0.5*(ry*ry-1.);
				epsz = 0.5*(rz*rz-1.);
				gamyz= ry*rz;
				gamxy= rx*ry;
				gamxz= rx*rz;
				double GAcs = beamEA; //use beam cross-section to define this factor!
				double EAcs = GAcs;			//stiffness factor for cross-section deformation and axial elongation

				double fact = 0.5*lx;
				//UO() << "EAcs=" << EAcs << "\n"; 

				for (int i=1; i <= dim; i++)
				{
					//j=3,7: Sx(j)=Sz(j)=depsx=depsz=0
					depsy = Sy(3)*ry(i);
					temp((3-1)*dim+i) = epsy*depsy*EAcs + gamyz * GAcs * Sy(3)*rz(i);
					depsy = Sy(7)*ry(i);
					temp((7-1)*dim+i) = epsy*depsy*EAcs + gamyz * GAcs * Sy(7)*rz(i);

					//j=4,8: Sx(j)=Sy(j)=depsx=depsy=0
					depsz = Sz(4)*rz(i);
					temp((4-1)*dim+i) = epsz*depsz*EAcs + gamyz * GAcs * (Sz(4)*ry(i));
					depsz = Sz(8)*rz(i);
					temp((8-1)*dim+i) = epsz*depsz*EAcs + gamyz * GAcs * (Sz(8)*ry(i));

				}
				fadd += fact*temp;
			}

			Vector3D rxs[2];
			//Vector3D rxslin[2];
			Vector3D rys[2];
			Vector3D rzs[2];
			double gammay[2]; //shear angle, rotation around y
			double gammaz[2]; //shear angle, rotation around z

			const int uselin = 0;
			const int altmode = 0; //use average gamma + relative gamma

			for (int i=0; i<2; i++) 
			{
				rxs[i] = AT * tx[i];

				rys[i] = AT * ty[i];
				rzs[i] = AT * tz[i];

				gammaz[i] =( rxs[i].Y() + rys[i].X()); //in fact would be -1*
				gammay[i] =( rxs[i].Z() + rzs[i].X());
			}

			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			//integration of elongation and bending:

			GetIntegrationRule(x1,w1,orderWl);
			for (int i1=1; i1<=x1.GetLen(); i1++)
			{
				double x = x1(i1); 
				double x2 = x*x;

				//rx:
				//GetS0x(Sx,Vector3D(x,0,0));

				Sx(1) = (-1.5+1.5*x2)/lx;
				Sx(2) = (-0.25-0.5*x+0.75*x2);
				Sx(5) = -Sx(1);//(1.5-1.5*x2)/lx;
				Sx(6) = Sx(2)+x;//(-0.25+0.5*x+0.75*x2);

				rx.X()  = Sx(1)*xg(3*0+1);
				rx.Y()  = Sx(1)*xg(3*0+2);
				rx.Z()  = Sx(1)*xg(3*0+3);					
				rx.X() += Sx(2)*xg(3*1+1);
				rx.Y() += Sx(2)*xg(3*1+2);
				rx.Z() += Sx(2)*xg(3*1+3);					
				rx.X() += Sx(5)*xg(3*4+1);
				rx.Y() += Sx(5)*xg(3*4+2);
				rx.Z() += Sx(5)*xg(3*4+3);					
				rx.X() += Sx(6)*xg(3*5+1);
				rx.Y() += Sx(6)*xg(3*5+2);
				rx.Z() += Sx(6)*xg(3*5+3);					

				//rxx:
				//GetS0xx(Sxx,Vector3D(x,0,0));
				Sxx(1) =  4./Sqr(lx)*3.0/2.0*x;
				Sxx(2) =  1./lx*(-1.0+3.0*x);
				Sxx(5) = -Sxx(1);
				Sxx(6) =  1./lx*(1.0+3.0*x);
				rxx.X()  = Sxx(1)*xg(3*(1-1)+1);
				rxx.Y()  = Sxx(1)*xg(3*(1-1)+2);
				rxx.Z()  = Sxx(1)*xg(3*(1-1)+3);					
				rxx.X() += Sxx(2)*xg(3*(2-1)+1);
				rxx.Y() += Sxx(2)*xg(3*(2-1)+2);
				rxx.Z() += Sxx(2)*xg(3*(2-1)+3);					
				rxx.X() += Sxx(5)*xg(3*(5-1)+1);
				rxx.Y() += Sxx(5)*xg(3*(5-1)+2);
				rxx.Z() += Sxx(5)*xg(3*(5-1)+3);					
				rxx.X() += Sxx(6)*xg(3*(6-1)+1);
				rxx.Y() += Sxx(6)*xg(3*(6-1)+2);
				rxx.Z() += Sxx(6)*xg(3*(6-1)+3);					

				//for shear deformation:
				txlin[0] += rx*(1.-x)*0.5 * w1(i1);
				txlin[1] += rx*(1.+x)*0.5 * w1(i1);
				
				for (int i = 1; i <= NS(); i++)
				{
					txlin_delta[0](i) += Sx(i)*(1.-x)*0.5 * w1(i1);
					txlin_delta[1](i) += Sx(i)*(1.+x)*0.5 * w1(i1);
				}

				rxx_int += rxx * 0.5 * w1(i1) * lx;

				for (int i = 1; i <= NS(); i++)
				{
					rxx_int_delta(i) += Sxx(i)* 0.5 * w1(i1) * lx;
				}


				//GetS0y(Sy,Vector3D(x,0,0));
				Sy(3) = -0.5*(-1.0+x);
				Sy(7) = 0.5*(1.0+x);
				ry.X()  = Sy(3)*xg(3*2+1);
				ry.Y()  = Sy(3)*xg(3*2+2);
				ry.Z()  = Sy(3)*xg(3*2+3);					
				ry.X() += Sy(7)*xg(3*6+1);
				ry.Y() += Sy(7)*xg(3*6+2);
				ry.Z() += Sy(7)*xg(3*6+3);					

				//GetS0z(Sz,Vector3D(x,0,0));
				Sz(4) = -0.5*(-1.0+x);
				Sz(8) = 0.5*(1.0+x);
				rz.X()  = Sz(4)*xg(3*3+1);
				rz.Y()  = Sz(4)*xg(3*3+2);
				rz.Z()  = Sz(4)*xg(3*3+3);					
				rz.X() += Sz(8)*xg(3*7+1);
				rz.Y() += Sz(8)*xg(3*7+2);
				rz.Z() += Sz(8)*xg(3*7+3);					

				tylin[0] += ry*(1.-x)*0.5 * w1(i1);
				tylin[1] += ry*(1.+x)*0.5 * w1(i1);
				tzlin[0] += rz*(1.-x)*0.5 * w1(i1);
				tzlin[1] += rz*(1.+x)*0.5 * w1(i1);





				Vector3D rxxs = AT*rxx;
				double wxxs = rxxs(3)-(gammay[1]-gammay[0])/lx;
				double vxxs = rxxs(2)+(gammaz[1]-gammaz[0])/lx;

				Vector3D rxs = AT*rx;
				double uxs = rxs(1)-1.;

				//integration factors:
				double facty = beamEIy * lx * 0.5 * w1(i1);
				double factz = beamEIz * lx * 0.5 * w1(i1);

				//integration factors:
				double factx = lx * 0.5 * w1(i1) * beamEA;
				//epsx = 0.5*(rx*rx-1.);

				temp.SetAll(0.);
				for (int i=1; i <= dim; i++)
				{
					//j=1,2,5,6:
					//temp((1-1)*dim+i) += factx*Sx(1)*rx(i)*(epsx);
					//temp((2-1)*dim+i) += factx*Sx(2)*rx(i)*(epsx);
					//temp((5-1)*dim+i) += factx*Sx(5)*rx(i)*(epsx);
					//temp((6-1)*dim+i) += factx*Sx(6)*rx(i)*(epsx);
					temp((1-1)*dim+i) += factx*uxs*Sx(1)*AT(1,i);
					temp((2-1)*dim+i) += factx*uxs*Sx(2)*AT(1,i);
					temp((5-1)*dim+i) += factx*uxs*Sx(5)*AT(1,i);
					temp((6-1)*dim+i) += factx*uxs*Sx(6)*AT(1,i);


					//j=1,2,5,6:
					temp((1-1)*dim+i) += facty*wxxs*Sxx(1)*AT(3,i);
					temp((2-1)*dim+i) += facty*wxxs*Sxx(2)*AT(3,i);
					temp((5-1)*dim+i) += facty*wxxs*Sxx(5)*AT(3,i);
					temp((6-1)*dim+i) += facty*wxxs*Sxx(6)*AT(3,i);

					temp((1-1)*dim+i) += factz*vxxs*Sxx(1)*AT(2,i);
					temp((2-1)*dim+i) += factz*vxxs*Sxx(2)*AT(2,i);
					temp((5-1)*dim+i) += factz*vxxs*Sxx(5)*AT(2,i);
					temp((6-1)*dim+i) += factz*vxxs*Sxx(6)*AT(2,i);

					for (int k=0; k <=1; k++)
					{
						double sign = 1;
						if (k==1) sign = -1;
						//gammaz[i] =( rxs[i].Y() + rys[i].X()); //in fact would be -1*
						//gammay[i] =( rxs[i].Z() + rzs[i].X());
						//double wxxs = rxxs(3)-(gammay[1]-gammay[0])/lx;
						//double vxxs = rxxs(2)+(gammaz[1]-gammaz[0])/lx;

						//delta rx:
						temp(3+i+12*k) += -1.*factz*vxxs*AT(2,i)/lx*sign; //y-component of AT
						temp(3+i+12*k) +=     facty*wxxs*AT(3,i)/lx*sign; //z-component of AT

						//delta ry:
						temp(6+i+12*k)+= -1.*factz*vxxs*AT(1,i)/lx*sign; //x-component of AT

						//delta rz:
						temp(9+i+12*k)+=     facty*wxxs*AT(1,i)/lx*sign; //x-component of AT
					}

				}

				//include terms with delta AT
				if (includedeltaterms)
				{
					ComputeCorotationalFrameDATvdq(rx,temp3);
					for (int i=1; i <= SOS(); i++)
					{
						temp(i) += temp3(1,i)*uxs*factx;
					}

					ComputeCorotationalFrameDATvdq(rxx,temp3);
					for (int i=1; i <= SOS(); i++)
					{
						temp(i) += temp3(2,i)*vxxs*factz;
						temp(i) += temp3(3,i)*wxxs*facty;
					}
				}

				fadd += temp;
			}

			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			//Torsion:
			//linearized angle of tz[0]=rz(1) and tz[1]=rz(2)
			//theta = [AT*(tz[0]-tz[1])].Y
			Vector3D rz1rz2 = tz[0] - tz[1];
			Vector3D ATrz1rz2 = AT * rz1rz2;

			double theta = ATrz1rz2.Y();
			//delta AT * rz1rz2:
			if (includedeltaterms)
			{
				ComputeCorotationalFrameDATvdq(rz1rz2, temp3);
				for (int i=1; i <= SOS(); i++)
				{
					temp(i) = temp3(2,i); //take y-component of difference of rz2 and rz1
				}
			}
			else
			{
				temp.SetAll(0.);
			}
			//AT * delta rz1rz2 = AT * (delta[q10,q11,q12] - delta[q22,q23,q24]):

			for (int i = 1; i <= Dim(); i++)
			{
				temp(9+i)    += AT(2,i); //y-component of AT
				temp(12+9+i) +=-AT(2,i); //y-component of AT
			}

			double fact_theta = beamGJkx / lx * theta;
			temp *= fact_theta;
			fadd += temp;

			
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			//shear deformation:

			for (int i=0; i<2; i++) 
			{
				/*if (GetMBS()->GetTime() > 0.9) 
				{
					UO() << "x=" << x_init(1) << ", ";
					UO() << "rzlin=" << tzlin[i] << ",rz=" << tz[i]
				 << "rx=" << tx[i] << ",rx=" << tx[i]<< ", gammay=" << gammay[i] << "\n";
				  //UO() << "beamGAky=" << beamGAky << ", beamGAkz=" << beamGAkz << "\n";
				}*/

				//if (GetMBS()->GetTime() > 0.02 && i == 0) 
				//	UO() << "gammay=" << gammay[i] << ", rxz=" << rxs[i].Z() << ", rzx=" << rzs[i].X() << "\n";
				//UO() << "beamGAky=" << beamGAky << ", beamGAkz=" << beamGAkz << "\n";
				//UO() << "beamEIy=" << beamEIy << ", beamEIz=" << beamEIz << "\n";

				double facty = 0.5 * lx * beamGAky * gammay[i];
				double factz = 0.5 * lx * beamGAkz * gammaz[i];
				if (altmode)
				{
					facty = 0.5*lx * beamGAky * (gammay[1] + gammay[2]);
					factz = 0.5*lx * beamGAkz * (gammaz[1] + gammaz[2]);
				}

				temp.SetAll(0.); //not necessary!


				//compute terms with delta AT
				if (includedeltaterms && 0)
				{
					ComputeCorotationalFrameDATvdq(tx[i], temp3);
					for (int j=1; j <= SOS(); j++)
					{
						temp(j)+= factz*temp3(2,j);
						temp(j)+= facty*temp3(3,j);
					}

					ComputeCorotationalFrameDATvdq(ty[i], temp3);
					for (int j=1; j <= SOS(); j++)
					{
						temp(j)+= factz*temp3(1,j);
					}

					ComputeCorotationalFrameDATvdq(tz[i], temp3);
					for (int j=1; j <= SOS(); j++)
					{
						temp(j)+= facty*temp3(1,j);
					}
				}

				//compute terms with delta rx, delta ry, delta rz
				for (int j = 1; j <= Dim(); j++)
				{
					//delta rx:
					temp(3+j+12*i) += factz*AT(2,j); //y-component of AT
					temp(3+j+12*i) += facty*AT(3,j); //z-component of AT

					//delta ry:
					temp(6+j+12*i)+= factz*AT(1,j); //x-component of AT

					//delta rz:
					temp(9+j+12*i)+= facty*AT(1,j); //x-component of AT
				}
				fadd += temp;
			}
		}

		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		if (elasticforce_beam == 3)
		{
			//elastic line approach from Arend Schwab
			//Version of elastic forces in Hiro paper/Sound and Vibration
			temp.SetLen(SOS());
			fadd.SetLen(SOS());
			fadd.SetAll(0);
			//optimized version:

			ApplyT(xg);


			static Vector Sx; Sx.SetLen(NS()); Sx.SetAll(0);  //rx
			static Vector Sy; Sy.SetLen(NS()); Sy.SetAll(0);   //ry
			static Vector Sz;	Sz.SetLen(NS()); Sz.SetAll(0);   //rz
			static Vector Sxx; Sxx.SetLen(NS()); Sxx.SetAll(0); //rxx
			static Vector Syx; Syx.SetLen(NS()); Syx.SetAll(0); //ryx
			static Vector Szx; Szx.SetLen(NS()); Szx.SetAll(0); //rzx

			double epsx;
			double epsy;
			double epsz;
			double gamyz;
			double depsy;
			double depsz;

			Vector3D rx, ry, rz, rxx, ryx, rzx;
			//integration of elongation and cross-section change:
			GetIntegrationRuleLobatto(x1,w1,orderWl);
			for (int i1=1; i1<=x1.GetLen(); i1++)
			{
				double x = x1(i1); 
				double x2 = x*x;

				//rx:
				//GetS0x(Sx,Vector3D(x,0,0));

				Sx(1) = (-1.5+1.5*x2)/lx;
				Sx(2) = (-0.25-0.5*x+0.75*x2);
				Sx(5) = -Sx(1);//(1.5-1.5*x2)/lx;
				Sx(6) = Sx(2)+x;//(-0.25+0.5*x+0.75*x2);

				rx.X()  = Sx(1)*xg(3*0+1);
				rx.Y()  = Sx(1)*xg(3*0+2);
				rx.Z()  = Sx(1)*xg(3*0+3);					
				rx.X() += Sx(2)*xg(3*1+1);
				rx.Y() += Sx(2)*xg(3*1+2);
				rx.Z() += Sx(2)*xg(3*1+3);					
				rx.X() += Sx(5)*xg(3*4+1);
				rx.Y() += Sx(5)*xg(3*4+2);
				rx.Z() += Sx(5)*xg(3*4+3);					
				rx.X() += Sx(6)*xg(3*5+1);
				rx.Y() += Sx(6)*xg(3*5+2);
				rx.Z() += Sx(6)*xg(3*5+3);					

				//GetS0y(Sy,Vector3D(x,0,0));
				Sy(3) = -0.5*(-1.0+x);
				Sy(7) = 0.5*(1.0+x);
				ry.X()  = Sy(3)*xg(3*2+1);
				ry.Y()  = Sy(3)*xg(3*2+2);
				ry.Z()  = Sy(3)*xg(3*2+3);					
				ry.X() += Sy(7)*xg(3*6+1);
				ry.Y() += Sy(7)*xg(3*6+2);
				ry.Z() += Sy(7)*xg(3*6+3);					

				//GetS0z(Sz,Vector3D(x,0,0));
				Sz(4) = -0.5*(-1.0+x);
				Sz(8) = 0.5*(1.0+x);
				rz.X()  = Sz(4)*xg(3*3+1);
				rz.Y()  = Sz(4)*xg(3*3+2);
				rz.Z()  = Sz(4)*xg(3*3+3);					
				rz.X() += Sz(8)*xg(3*7+1);
				rz.Y() += Sz(8)*xg(3*7+2);
				rz.Z() += Sz(8)*xg(3*7+3);					


				//integration factors:
				double fact = lx * 0.5 * w1(i1);

				epsx = 0.5*(rx*rx-1.) - Wlepsx0(i1);
				epsy = 0.5*(ry*ry-1.)*factstiffWl - Wlepsy0(i1);
				epsz = 0.5*(rz*rz-1.)*factstiffWl - Wlepsz0(i1);
				gamyz= ry*rz*factstiffWl - Wlgamyz0(i1);

				double GAcs = beamEA / (2.*(1+nu)); //use beam cross-section to define this factor!
				double gf = 2.*GAcs/(1.-2.*nu);			//stiffness factor for cross-section deformation and axial elongation


				//delta eps:
				double sumepsi = epsx+epsy+epsz;
				for (int i=1; i <= dim; i++)
				{
					//j=1,2,5,6: Sy(j)=Sz(j)=depsy=depsz=0
					temp((1-1)*dim+i) = gf*Sx(1)*rx(i)*(sumepsi*nu+epsx*(1.-2.*nu));
					temp((2-1)*dim+i) = gf*Sx(2)*rx(i)*(sumepsi*nu+epsx*(1.-2.*nu));
					temp((5-1)*dim+i) = gf*Sx(5)*rx(i)*(sumepsi*nu+epsx*(1.-2.*nu));
					temp((6-1)*dim+i) = gf*Sx(6)*rx(i)*(sumepsi*nu+epsx*(1.-2.*nu));

					//j=3,7: Sx(j)=Sz(j)=depsx=depsz=0
					depsy = Sy(3)*ry(i);
					temp((3-1)*dim+i) = sumepsi*nu*gf*depsy+epsy*depsy*(1.-2.*nu)*gf+gamyz * GAcs * Sy(3)*rz(i);
					depsy = Sy(7)*ry(i);
					temp((7-1)*dim+i) = sumepsi*nu*gf*depsy+epsy*depsy*(1.-2.*nu)*gf+gamyz * GAcs * Sy(7)*rz(i);

					//j=4,8: Sx(j)=Sy(j)=depsx=depsy=0
					depsz = Sz(4)*rz(i);
					temp((4-1)*dim+i) = sumepsi * nu*gf*depsz+epsz*depsz*(1.-2.*nu)*gf+gamyz * GAcs * (Sz(4)*ry(i));
					depsz = Sz(8)*rz(i);
					temp((8-1)*dim+i) = sumepsi * nu*gf*depsz+epsz*depsz*(1.-2.*nu)*gf+gamyz * GAcs * (Sz(8)*ry(i));

				}
				temp *= fact;
				fadd += temp;
			}


			//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			//Hellinger Reissner principle:
			//note that xi=0..1 in the paper of Schwab and x = -1..+1 

			//1. compute Wxz1, Wxz2, Wxy1, Wxy2
			//2. compute delta Wx... times Wx... and sum up in temp
			//UO() << "Hellinger Reisner\n";

			double Wxz1f = 0; 
			double Wxz2f = 0;
			double Wxy1f = 0;
			double Wxy2f = 0;
			//shear deformation:
			GetIntegrationRuleLobatto(x1,w1,orderWsHR); 
			for (int i1=1; i1<=x1.GetLen(); i1++)
			{
				double x = x1(i1);

				double x2 = x*x;
				//rx:
				//GetS0x(Sx,Vector3D(x,0,0));
				Sx(1) = (-1.5+1.5*x2)/lx;
				Sx(2) = (-0.25-0.5*x+0.75*x2);
				Sx(5) = -Sx(1);//(1.5-1.5*x2)/lx;
				Sx(6) = Sx(2)+x;//(-0.25+0.5*x+0.75*x2);
				rx.X()  = Sx(1)*xg(3*0+1);
				rx.Y()  = Sx(1)*xg(3*0+2);
				rx.Z()  = Sx(1)*xg(3*0+3);					
				rx.X() += Sx(2)*xg(3*1+1);
				rx.Y() += Sx(2)*xg(3*1+2);
				rx.Z() += Sx(2)*xg(3*1+3);					
				rx.X() += Sx(5)*xg(3*4+1);
				rx.Y() += Sx(5)*xg(3*4+2);
				rx.Z() += Sx(5)*xg(3*4+3);					
				rx.X() += Sx(6)*xg(3*5+1);
				rx.Y() += Sx(6)*xg(3*5+2);
				rx.Z() += Sx(6)*xg(3*5+3);					

				//GetS0y(Sy,Vector3D(x,0,0));
				Sy(3) = -0.5*(-1.0+x);
				Sy(7) = 0.5*(1.0+x);
				ry.X()  = Sy(3)*xg(3*2+1);
				ry.Y()  = Sy(3)*xg(3*2+2);
				ry.Z()  = Sy(3)*xg(3*2+3);					
				ry.X() += Sy(7)*xg(3*6+1);
				ry.Y() += Sy(7)*xg(3*6+2);
				ry.Z() += Sy(7)*xg(3*6+3);					

				//GetS0z(Sz,Vector3D(x,0,0));
				Sz(4) = -0.5*(-1.0+x);
				Sz(8) = 0.5*(1.0+x);
				rz.X()  = Sz(4)*xg(3*3+1);
				rz.Y()  = Sz(4)*xg(3*3+2);
				rz.Z()  = Sz(4)*xg(3*3+3);					
				rz.X() += Sz(8)*xg(3*7+1);
				rz.Y() += Sz(8)*xg(3*7+2);
				rz.Z() += Sz(8)*xg(3*7+3);					


				//integration factors:
				double fact = lx * 0.5 * w1(i1);

				double gamxy= fact*rx*ry;
				double gamxz= fact*rx*rz;

				Wxy1f += 0.5*(1-x) * gamxy;
				Wxy2f += 0.5*(1+x) * gamxy;
				Wxz1f += 0.5*(1-x) * gamxz;
				Wxz2f += 0.5*(1+x) * gamxz;

			}

			GetIntegrationRuleLobatto(x1,w1,orderWsHR);
			//now compute variations:
			for (int i1=1; i1<=x1.GetLen(); i1++)
			{
				double x = x1(i1);

				double x2 = x*x;
				//rx:
				//GetS0x(Sx,Vector3D(x,0,0));
				Sx(1) = (-1.5+1.5*x2)/lx;
				Sx(2) = (-0.25-0.5*x+0.75*x2);
				Sx(5) = -Sx(1);//(1.5-1.5*x2)/lx;
				Sx(6) = Sx(2)+x;//(-0.25+0.5*x+0.75*x2);
				rx.X()  = Sx(1)*xg(3*0+1);
				rx.Y()  = Sx(1)*xg(3*0+2);
				rx.Z()  = Sx(1)*xg(3*0+3);					
				rx.X() += Sx(2)*xg(3*1+1);
				rx.Y() += Sx(2)*xg(3*1+2);
				rx.Z() += Sx(2)*xg(3*1+3);					
				rx.X() += Sx(5)*xg(3*4+1);
				rx.Y() += Sx(5)*xg(3*4+2);
				rx.Z() += Sx(5)*xg(3*4+3);					
				rx.X() += Sx(6)*xg(3*5+1);
				rx.Y() += Sx(6)*xg(3*5+2);
				rx.Z() += Sx(6)*xg(3*5+3);					

				//GetS0y(Sy,Vector3D(x,0,0));
				Sy(3) = -0.5*(-1.0+x);
				Sy(7) = 0.5*(1.0+x);
				ry.X()  = Sy(3)*xg(3*2+1);
				ry.Y()  = Sy(3)*xg(3*2+2);
				ry.Z()  = Sy(3)*xg(3*2+3);					
				ry.X() += Sy(7)*xg(3*6+1);
				ry.Y() += Sy(7)*xg(3*6+2);
				ry.Z() += Sy(7)*xg(3*6+3);					

				//GetS0z(Sz,Vector3D(x,0,0));
				Sz(4) = -0.5*(-1.0+x);
				Sz(8) = 0.5*(1.0+x);
				rz.X()  = Sz(4)*xg(3*3+1);
				rz.Y()  = Sz(4)*xg(3*3+2);
				rz.Z()  = Sz(4)*xg(3*3+3);					
				rz.X() += Sz(8)*xg(3*7+1);
				rz.Y() += Sz(8)*xg(3*7+2);
				rz.Z() += Sz(8)*xg(3*7+3);					


				//integration factors:
				double fact = 0.5 * w1(i1); // length cancels out with length omitted in Fy and Fz, see A.Schwab

				//factor for delta gamxy components:
				double k1 = beamGAky * fact*(0.5*(1-x)*(4.*Wxy1f - 2.*Wxy2f) + 0.5*(1+x)*(-2.*Wxy1f + 4.*Wxy2f)) - WsHRk10(i1);
				//factor for delta gamxz components:
				double k2 = beamGAkz * fact*(0.5*(1-x)*(4.*Wxz1f - 2.*Wxz2f) + 0.5*(1+x)*(-2.*Wxz1f + 4.*Wxz2f)) - WsHRk20(i1);

				//delta eps:
				for (int i=1; i <= dim; i++)
				{
					//j=1,2,5,6: Sy=Sz=0;
					temp((1-1)*dim+i) = (k1 * (Sx(1)*ry(i)) + k2 * (Sx(1)*rz(i)));
					temp((2-1)*dim+i) = (k1 * (Sx(2)*ry(i)) + k2 * (Sx(2)*rz(i)));
					temp((5-1)*dim+i) = (k1 * (Sx(5)*ry(i)) + k2 * (Sx(5)*rz(i)));
					temp((6-1)*dim+i) = (k1 * (Sx(6)*ry(i)) + k2 * (Sx(6)*rz(i)));
					//j=3,7: Sx=Sz=0;
					temp((3-1)*dim+i) = k1 * Sy(3)*rx(i);
					temp((7-1)*dim+i) = k1 * Sy(7)*rx(i);
					//j=4,8: Sx=Sy=0;
					temp((4-1)*dim+i) = k2 * Sz(4)*rx(i);
					temp((8-1)*dim+i) = k2 * Sz(8)*rx(i);
				}
				fadd += temp;
			}


			//torsion:
			GetIntegrationRule(x1,w1,orderWt);
			temp.SetAll(0);
			for (int i1=1; i1<=x1.GetLen(); i1++)
			{
				double x = x1(i1);
				//double x2 = x*x;
				//GetS0y(Sy,Vector3D(x,0,0));
				Sy(3) = -0.5*(-1.0+x);
				Sy(7) = 0.5*(1.0+x);
				ry.X()  = Sy(3)*xg(3*2+1);
				ry.Y()  = Sy(3)*xg(3*2+2);
				ry.Z()  = Sy(3)*xg(3*2+3);					
				ry.X() += Sy(7)*xg(3*6+1);
				ry.Y() += Sy(7)*xg(3*6+2);
				ry.Z() += Sy(7)*xg(3*6+3);					

				//GetS0z(Sz,Vector3D(x,0,0));
				Sz(4) = -0.5*(-1.0+x);
				Sz(8) = 0.5*(1.0+x);
				rz.X()  = Sz(4)*xg(3*3+1);
				rz.Y()  = Sz(4)*xg(3*3+2);
				rz.Z()  = Sz(4)*xg(3*3+3);					
				rz.X() += Sz(8)*xg(3*7+1);
				rz.Y() += Sz(8)*xg(3*7+2);
				rz.Z() += Sz(8)*xg(3*7+3);					

				//ryx:
				//GetS0yx(Syx,Vector3D(x,0,0));
				Syx(3) = -1./lx;
				Syx(7) = 1./lx;
				ryx.X()  = Syx(3)*xg(3*2+1);
				ryx.Y()  = Syx(3)*xg(3*2+2);
				ryx.Z()  = Syx(3)*xg(3*2+3);					
				ryx.X() += Syx(7)*xg(3*6+1);
				ryx.Y() += Syx(7)*xg(3*6+2);
				ryx.Z() += Syx(7)*xg(3*6+3);					

				//rzx:
				//GetS0zx(Szx,Vector3D(x,0,0));
				Szx(4) = -1./lx;
				Szx(8) = 1./lx;
				rzx.X()  = Szx(4)*xg(3*3+1);
				rzx.Y()  = Szx(4)*xg(3*3+2);
				rzx.Z()  = Szx(4)*xg(3*3+3);					
				rzx.X() += Szx(8)*xg(3*7+1);
				rzx.Y() += Szx(8)*xg(3*7+2);
				rzx.Z() += Szx(8)*xg(3*7+3);					


				//integration factors:
				double fact = beamGJkx * lx * 0.5 * w1(i1);

				double kapx= 0.5 * (rz*ryx - ry*rzx) - Wtkapx0(i1);

				//delta eps:
				for (int i=1; i <= dim; i++)
				{
					//j=3,7: Sx=Sz=0;
					temp((3-1)*dim+i) = 0.5*(Syx(3)*rz(i) - Sy(3)*rzx(i));
					temp((7-1)*dim+i) = 0.5*(Syx(7)*rz(i) - Sy(7)*rzx(i));
					//j=4,8: Sx=Sy=0;
					temp((4-1)*dim+i) = 0.5*(Sz(4)*ryx(i) - Szx(4)*ry(i));
					temp((8-1)*dim+i) = 0.5*(Sz(8)*ryx(i) - Szx(8)*ry(i));

					/*
					for (int j=1; j <= ns; j++)
					{
					temp((j-1)*dim+i) = 0.5*(Sz(j)*ryx(i)+Syx(j)*rz(i) - (Sy(j)*rzx(i)+Szx(j)*ry(i)));
					}*/
				}
				temp *= fact*kapx;
				fadd += temp;
			}

			//bending:
			GetIntegrationRule(x1,w1,orderWb);
			for (int i1=1; i1<=x1.GetLen(); i1++)
			{
				double x = x1(i1);

				double x2 = x*x;
				//GetS0y(Sy,Vector3D(x,0,0));
				Sy(3) = -0.5*(-1.0+x);
				Sy(7) = 0.5*(1.0+x);
				ry.X()  = Sy(3)*xg(3*2+1);
				ry.Y()  = Sy(3)*xg(3*2+2);
				ry.Z()  = Sy(3)*xg(3*2+3);					
				ry.X() += Sy(7)*xg(3*6+1);
				ry.Y() += Sy(7)*xg(3*6+2);
				ry.Z() += Sy(7)*xg(3*6+3);					

				//GetS0z(Sz,Vector3D(x,0,0));
				Sz(4) = -0.5*(-1.0+x);
				Sz(8) = 0.5*(1.0+x);
				rz.X()  = Sz(4)*xg(3*3+1);
				rz.Y()  = Sz(4)*xg(3*3+2);
				rz.Z()  = Sz(4)*xg(3*3+3);					
				rz.X() += Sz(8)*xg(3*7+1);
				rz.Y() += Sz(8)*xg(3*7+2);
				rz.Z() += Sz(8)*xg(3*7+3);					

				//rxx:
				//GetS0xx(Sxx,Vector3D(x,0,0));
				Sxx(1) =  4./Sqr(lx)*3.0/2.0*x;
				Sxx(2) =  1./lx*(-1.0+3.0*x);
				Sxx(5) = -Sxx(1);
				Sxx(6) =  1./lx*(1.0+3.0*x);
				rxx.X()  = Sxx(1)*xg(3*(1-1)+1);
				rxx.Y()  = Sxx(1)*xg(3*(1-1)+2);
				rxx.Z()  = Sxx(1)*xg(3*(1-1)+3);					
				rxx.X() += Sxx(2)*xg(3*(2-1)+1);
				rxx.Y() += Sxx(2)*xg(3*(2-1)+2);
				rxx.Z() += Sxx(2)*xg(3*(2-1)+3);					
				rxx.X() += Sxx(5)*xg(3*(5-1)+1);
				rxx.Y() += Sxx(5)*xg(3*(5-1)+2);
				rxx.Z() += Sxx(5)*xg(3*(5-1)+3);					
				rxx.X() += Sxx(6)*xg(3*(6-1)+1);
				rxx.Y() += Sxx(6)*xg(3*(6-1)+2);
				rxx.Z() += Sxx(6)*xg(3*(6-1)+3);					


				//integration factors:
				double facty = beamEIy * lx * 0.5 * w1(i1);
				double factz = beamEIz * lx * 0.5 * w1(i1);

				double kapy=  rz*rxx  - Wbkapy0(i1);
				double kapz=-(ry*rxx) - Wbkapz0(i1);

				//delta eps:
				for (int i=1; i <= dim; i++)
				{
					//j=1,2,5,6: Sy=Sz=0;
					temp((1-1)*dim+i) = facty*kapy*(Sxx(1)*rz(i)) + factz*kapz*(-Sxx(1)*ry(i));
					temp((2-1)*dim+i) = facty*kapy*(Sxx(2)*rz(i)) + factz*kapz*(-Sxx(2)*ry(i));
					temp((5-1)*dim+i) = facty*kapy*(Sxx(5)*rz(i)) + factz*kapz*(-Sxx(5)*ry(i));
					temp((6-1)*dim+i) = facty*kapy*(Sxx(6)*rz(i)) + factz*kapz*(-Sxx(6)*ry(i));
					//j=3,7: Sx=Sz=0;
					temp((3-1)*dim+i) = factz*kapz*(-Sy(3)*rxx(i));
					temp((7-1)*dim+i) = factz*kapz*(-Sy(7)*rxx(i));
					//j=4,8: Sx=Sy=0;
					temp((4-1)*dim+i) = facty*kapy*(Sz(4)*rxx(i));
					temp((8-1)*dim+i) = facty*kapy*(Sz(8)*rxx(i));


					/*				for (int j=1; j <= ns; j++)
					{
					temp((j-1)*dim+i) = facty*kapy*(Sz(j)*rxx(i)+Sxx(j)*rz(i)) + factz*kapz*(-Sy(j)*rxx(i)-Sxx(j)*ry(i));
					}*/
				}
				//temp *= 1;
				fadd += temp;
			}

		}
		else if (elasticforce_beam < 3)
		{
			//UO() << "standard ANCForig\n";
			static Vector u;
			u.SetLen(sos);
			int useu = 1;
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

			double la=Em * nu / ((1.+nu)*(1.-2.*nu));
			double mu=Em / 2. / (1.+nu);

			Matrix3D strain, piola1, F;

			temp.SetLen(SOS());
			fadd.SetLen(SOS());
			fadd.SetAll(0);

			GetIntegrationRule(x1,w1,orderx); //x_i ... integration points, w_i ... integration weights
			GetIntegrationRule(x2,w2,orderyz);
			GetIntegrationRule(x3,w3,orderyz);
			int kx1 = x2.Length()*x3.Length();
			int kx2 = x3.Length();

			for (int i1=1; i1<=x1.GetLen(); i1++)
			{
				for (int i2=1; i2<=x2.GetLen(); i2++)
				{
					for (int i3=1; i3<=x3.GetLen(); i3++)
					{
						//TMStartTimer(20);
						int i,j,k;
						int ind = (i1-1)*kx1+(i2-1)*kx2+(i3-1);

						// compute F 
						F.SetAll(0);
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
									//G(j,k) += DS(k,i)*u((i-1)*dim+j);
								}
							}
							if (useu) F(j,j) += 1;
						}

						//jacinv = jacinv.GetTp();
						//G = G * jacinv; //also possible!!!
						//F = Matrix3D(1) + G;

						int linear=0;
						if (linear)
						{
							F += Matrix3D(-1.);
							strain = 0.5*(F+F.GetTp());
							piola1 = ((2*mu) * strain + Matrix3D(la * strain.Trace()));
						}
						else
						{
							// Green-Lagrange strain tensor
							//strain = 0.5 * (F.GetTp() * F - I);
							F.GetATA2(strain);
							strain(1,1) -= 0.5; strain(2,2) -= 0.5; strain(3,3) -= 0.5;

							// Matrix3D Eplast = ....
							// strain -= Eplast; //strain = strain - Eplast;

							// Matrix3D piola2 = (2*mu) * strain + Matrix3D(la * strain.Trace());
							// piola1 = F*piola2;
							piola1 = F * ((2*mu) * strain + Matrix3D(la * strain.Trace()));
						}

						for (int j=1; j <= 3; j++)
						{
							for (int i = 0; i < ns; i++)
							{
								temp(3*i+j) = grad[ind](1, i+1)*piola1(j,1)
									+ grad[ind](2, i+1)*piola1(j,2)
									+ grad[ind](3, i+1)*piola1(j,3);
							}
						}

						//temp *= fabs (jacdet[ind]) * w1(i1)*w2(i2)*w3(i3);
						//fadd += temp;
						fadd.MultAdd(fabs (jacdet[ind]) * w1(i1)*w2(i2)*w3(i3),temp);
						//f -= fabs (jacdet[ind]) * w1(i1)*w2(i2)*w3(i3) * (B*piola1v);
						//TMStopTimer(21);
					}
				}
			}
		}
		else if (elasticforce_beam == 4)
		{
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			//Poisson locking eliminated

			static Vector u;
			u.SetLen(sos);
			int useu = 1;
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
			//ApplyT(u); // u = T*u
			//UO() << "u=" << u << "\n";

			double la=Em * nu / ((1.+nu)*(1.-2.*nu));
			double mu=Em / 2. / (1.+nu);
			//double ksT = 1./1.248424513011470; //shear correction for h=0.02, b=0.01
			//double ksT = 3.51176216306567e-009/1.4857142857142857e-8;//10.*(1+nu)/(12+11*nu); //shear correction for h=0.02, b=0.005
			double ksT = (beamGJkx) / (mu * ly*lz*(Sqr(ly) + Sqr(lz))/12.);//10.*(1+nu)/(12+11*nu); //shear correction for h=0.02, b=0.005
			double ks = 1;//10.*(1+nu)/(12+11*nu); //shear correction factor GA


			Matrix3D strain, piola1, F;
			Matrix3D jac, jacinv;
			double jacdet;

			static Matrix Dmat;
			static Matrix Dmat2;
			Dmat.SetSize(6,6);
			Dmat2.SetSize(6,6);
			static Vector strainvec;
			strainvec.SetLen(6);
			static Vector stressvec;
			stressvec.SetLen(6);

			GetDMatrix(Dmat, nu, Em);
			Dmat(4,4) *= 1; //no shear correction for Syz necessary
			Dmat(5,5) *= ksT;
			Dmat(6,6) *= ksT; 

			Dmat2.SetAll(0);
			Dmat2(1,1) = Em;
			Dmat2(2,2) = Em;
			Dmat2(3,3) = Em;
			Dmat2(4,4) = Dmat(4,4);
			Dmat2(5,5) = Dmat(5,5);
			Dmat2(6,6) = Dmat(6,6);

			Dmat -= Dmat2;
			Dmat(5,5) = ks-ksT; //at elastic line
			Dmat(6,6) = ks-ksT; //at elastic line

			//GetDMatrix(Dmat2,nu,Em);

			temp.SetLen(SOS());
			fadd.SetLen(SOS());
			fadd.SetAll(0);

			for (int kk=1; kk <= 2; kk++)
			{
				//kk=1: just poisson effect at element line
				//kk=2: no poisson effect, whole body

				if (kk == 2)
				{
					GetIntegrationRule(x1,w1,orderx); //x_i ... integration points, w_i ... integration weights
					GetIntegrationRule(x2,w2,orderyz);
					GetIntegrationRule(x3,w3,orderyz);
				}
				else
				{
					GetIntegrationRule(x1,w1,orderx); //x_i ... integration points, w_i ... integration weights
					GetIntegrationRule(x2,w2,1);
					GetIntegrationRule(x3,w3,1);
				}

				int kx1 = x2.Length()*x3.Length();
				int kx2 = x3.Length();

				for (int i1=1; i1<=x1.GetLen(); i1++)
				{
					for (int i2=1; i2<=x2.GetLen(); i2++)
					{
						for (int i3=1; i3<=x3.GetLen(); i3++)
						{
							//TMStartTimer(20);
							int i,j,k;
							int ind = (i1-1)*kx1+(i2-1)*kx2+(i3-1);
							Vector3D p(x1(i1),x2(i2),x3(i3));

							// compute F 
							F.SetAll(0);
							int l;
							static Matrix agrad; //\bar S^D

							GetDSMatrix0(p,DS);

							GetJacobi(jac,p,DS,x_init);
							jacdet = jac.Det();
							jac.GetInverse(jacinv);
							jacinv = jacinv.GetTp();

							agrad.SetSize(3,NS());
							Mult(jacinv, DS, agrad);

							for (j = 1; j <= dim; j++) 
							{
								for (i = 1; i <= ns; i++)
								{
									l = (i-1)*dim+j;
									//u = xg(l) - x_init(l);
									for (k = 1; k <= dim; k++)
									{
										F(j,k) += agrad(k,i)*u(l);
										//G(j,k) += DS(k,i)*u((i-1)*dim+j);
									}
								}
								if (useu) F(j,j) += 1;
							}

							// Green-Lagrange strain tensor
							//strain = 0.5 * (F.GetTp() * F - I);
							F.GetATA2(strain);
							strain(1,1) -= 0.5; strain(2,2) -= 0.5; strain(3,3) -= 0.5;

							strainvec(1) = strain(1,1);
							strainvec(2) = strain(2,2);
							strainvec(3) = strain(3,3);
							strainvec(4) = 2.*strain(2,3);
							strainvec(5) = 2.*strain(3,1);
							strainvec(6) = 2.*strain(1,2);

							if (kk == 1)
								Mult(Dmat,strainvec,stressvec);
							else
								Mult(Dmat2,strainvec,stressvec);

							//piola1 = F * ((2*mu) * strain + Matrix3D(la * strain.Trace()));
							piola1(1,1) = stressvec(1);
							piola1(2,2) = stressvec(2);
							piola1(3,3) = stressvec(3);
							piola1(2,3) = stressvec(4);
							piola1(3,1) = stressvec(5);
							piola1(1,2) = stressvec(6);
							piola1(3,2) = stressvec(4);
							piola1(1,3) = stressvec(5);
							piola1(2,1) = stressvec(6);

							piola1 = F*piola1;

							for (int j=1; j <= 3; j++)
							{
								for (int i = 0; i < ns; i++)
								{
									temp(3*i+j) = agrad(1, i+1)*piola1(j,1)
										+ agrad(2, i+1)*piola1(j,2)
										+ agrad(3, i+1)*piola1(j,3);
								}
							}

							//temp *= fabs (jacdet[ind]) * w1(i1)*w2(i2)*w3(i3);
							//fadd += temp;
							fadd.MultAdd(fabs (jacdet) * w1(i1)*w2(i2)*w3(i3),temp);
							//f -= fabs (jacdet[ind]) * w1(i1)*w2(i2)*w3(i3) * (B*piola1v);
							//TMStopTimer(21);
						}
					}
				}
			}
		}



		if (!(elasticforce_beam == 4)) ApplyTtp(fadd);

		f -= fadd;
		TMStopTimer(22);

	}

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

void ANCFBeam3D::GraduD(const Vector3D& ploc, const Vector& u, Matrix3D& gradu) const
{
	static Matrix DS;
	GetDSMatrix0(ploc,DS);

	Matrix3D jac, jacinv;
	GetJacobi(jac,ploc,DS,e0);

	jac.GetInverse(jacinv);
	jacinv = jacinv.GetTp();

	static Matrix grad;
	grad.SetSize(3,NS());
	Mult(jacinv, DS, grad);

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

void ANCFBeam3D::GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables)
{
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_displacement,
		FieldVariableDescriptor::FVCI_z);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_total_strain,
		FieldVariableDescriptor::FVCI_z, true);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_stress,
		FieldVariableDescriptor::FVCI_z, true);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_stress_mises);
}

double ANCFBeam3D::GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & local_position, bool flagD)
{
	if (!flagD)
		return FIELD_VARIABLE_NO_VALUE; //not implemented for computation!
	if (fvd.VariableType() == FieldVariableDescriptor::FVT_displacement) //show displacements
	{
		Vector3D u;
		Vector3D p0 = local_position;
		p0.Scale(2./lx,2./ly,2./lz);
		if(flagD) 
			u = GetDisplacementD(p0);
		else
			u = GetDisplacement(p0);

		return fvd.GetComponent(u);
	}

	if (elasticforce_beam != 3 || fvd.VariableType() == FieldVariableDescriptor::FVT_total_strain)
	{
		xgd.SetLen(SOS());
		GetDrawCoordinates(xgd);

		xgd = T*xgd;
		//transform xgd to displacements, actually not necessary ...
		xgd -= e0;
		Matrix3D F;
		GraduD(local_position,xgd,F);
		F += Matrix3D(1);
		Matrix3D strain = 0.5*(F.GetTp()*F-Matrix3D(1));
		//F.GetATA2(strain);
		//strain-=Matrix3D(0.5);
		//F -= Matrix3D(1);

		double la=Em * nu / ((1.+nu)*(1.-2.*nu));
		double mu=Em / 2. / (1.+nu);
		Matrix3D stress = (2*mu) * strain + Matrix3D(la * strain.Trace());

		if(fvd.VariableType() == FieldVariableDescriptor::FVT_stress_mises)
		{
			return stress.Mises();
		}
		else if (fvd.VariableType() == FieldVariableDescriptor::FVT_stress)
		{
			return fvd.GetComponent(stress);
		}
		else if (fvd.VariableType() == FieldVariableDescriptor::FVT_total_strain)
		{
			return fvd.GetComponent(strain);
		}
	}
	else
	{
		//only for rectangular cross-section!:
		xgd.SetLen(SOS());
		GetDrawCoordinates(xgd);

		xgd = T*xgd;

		double x0 = local_position.X();
		double yb = 0.5*local_position.Y()*ly;
		double zb = 0.5*local_position.Z()*lz;
		//global_uo << "z=" << local_position.Z() << ", zb=" << zb << "\n";

		double Gm = Em / (2.*(1+nu));
		double A = ly*lz;

		double ky = 10.*(1+nu)/(12.+11.*nu); //shear correction factors (square): 10*(1+nu)/(12.+11.*nu) 
		double kz = 10.*(1+nu)/(12.+11.*nu); //shear correction factors (square): 10*(1+nu)/(12.+11.*nu)
		double kx = 0.8436; //shear correction factor torsion  : 0.8436 (square)
		//kx=ky=kz=1;

		static Vector Sx; Sx.SetLen(NS());   //rx
		static Vector Sy; Sy.SetLen(NS());   //ry
		static Vector Sz;	Sz.SetLen(NS());   //rz
		static Vector Sxx; Sxx.SetLen(NS()); //rxx
		static Vector Syx; Syx.SetLen(NS()); //ryx
		static Vector Szx; Szx.SetLen(NS()); //rzx
		Vector3D rx,ry,rz,rxx,ryx,rzx;

		Vector3D rx1,ry1,rz1,ryx1,rzx1;
		Vector3D rx2,ry2,rz2,ryx2,rzx2;


		Matrix3D stress(0);

		//shear stresses
		//need integration of coefficitents Wxz and Wxy:
		double Wxz1f = 0; 
		double Wxz2f = 0;
		double Wxy1f = 0;
		double Wxy2f = 0;

		static Vector x1d; //not allowed to take variables from computation!!!
		static Vector w1d;

		GetIntegrationRuleLobatto(x1d,w1d,2); 
		for (int i1=1; i1<=x1d.GetLen(); i1++)
		{
			double x = x1d(i1);

			double x2 = x*x;
			//rx:
			//GetS0x(Sx,Vector3D(x,0,0));
			Sx(1) = (-1.5+1.5*x2)/lx;
			Sx(2) = (-0.25-0.5*x+0.75*x2);
			Sx(5) = -Sx(1);//(1.5-1.5*x2)/lx;
			Sx(6) = Sx(2)+x;//(-0.25+0.5*x+0.75*x2);
			rx.X()  = Sx(1)*xgd(3*0+1);
			rx.Y()  = Sx(1)*xgd(3*0+2);
			rx.Z()  = Sx(1)*xgd(3*0+3);					
			rx.X() += Sx(2)*xgd(3*1+1);
			rx.Y() += Sx(2)*xgd(3*1+2);
			rx.Z() += Sx(2)*xgd(3*1+3);					
			rx.X() += Sx(5)*xgd(3*4+1);
			rx.Y() += Sx(5)*xgd(3*4+2);
			rx.Z() += Sx(5)*xgd(3*4+3);					
			rx.X() += Sx(6)*xgd(3*5+1);
			rx.Y() += Sx(6)*xgd(3*5+2);
			rx.Z() += Sx(6)*xgd(3*5+3);					

			//GetS0y(Sy,Vector3D(x,0,0));
			Sy(3) = -0.5*(-1.0+x);
			Sy(7) = 0.5*(1.0+x);
			ry.X()  = Sy(3)*xgd(3*2+1);
			ry.Y()  = Sy(3)*xgd(3*2+2);
			ry.Z()  = Sy(3)*xgd(3*2+3);					
			ry.X() += Sy(7)*xgd(3*6+1);
			ry.Y() += Sy(7)*xgd(3*6+2);
			ry.Z() += Sy(7)*xgd(3*6+3);					

			//GetS0z(Sz,Vector3D(x,0,0));
			Sz(4) = -0.5*(-1.0+x);
			Sz(8) = 0.5*(1.0+x);
			rz.X()  = Sz(4)*xgd(3*3+1);
			rz.Y()  = Sz(4)*xgd(3*3+2);
			rz.Z()  = Sz(4)*xgd(3*3+3);					
			rz.X() += Sz(8)*xgd(3*7+1);
			rz.Y() += Sz(8)*xgd(3*7+2);
			rz.Z() += Sz(8)*xgd(3*7+3);					

			//integration factors:
			double fact = A * lx * 0.5 * w1d(i1);

			double gamxy= fact*rx*ry;
			double gamxz= fact*rx*rz;

			Wxy1f += 0.5*(1-x) * gamxy;
			Wxy2f += 0.5*(1+x) * gamxy;
			Wxz1f += 0.5*(1-x) * gamxz;
			Wxz2f += 0.5*(1+x) * gamxz;

			GetS0yx(Syx,Vector3D(x,0,0));
			GetS0zx(Szx,Vector3D(x,0,0));
			ryx.Set(0,0,0);
			rzx.Set(0,0,0);

			for (int i = 1; i <= NS(); i++)
			{
				ryx.X() += Syx(i)*xgd(3*(i-1)+1);
				ryx.Y() += Syx(i)*xgd(3*(i-1)+2);
				ryx.Z() += Syx(i)*xgd(3*(i-1)+3);					
				rzx.X() += Szx(i)*xgd(3*(i-1)+1);
				rzx.Y() += Szx(i)*xgd(3*(i-1)+2);
				rzx.Z() += Szx(i)*xgd(3*(i-1)+3);
			}


			if (x == -1) {rx1 = rx; ry1 = ry; rz1 = rz; ryx1 = ryx; rzx1 = rzx;}
			if (x == 1)  {rx2 = rx; ry2 = ry; rz2 = rz; ryx2 = ryx; rzx2 = rzx;}
		}

		GetS0x(Sx,Vector3D(x0,0,0));
		GetS0y(Sy,Vector3D(x0,0,0));
		GetS0z(Sz,Vector3D(x0,0,0));
		GetS0xx(Sxx,Vector3D(x0,0,0));
		GetS0yx(Syx,Vector3D(x0,0,0));
		GetS0zx(Szx,Vector3D(x0,0,0));

		rx.Set(0,0,0);
		ry.Set(0,0,0);
		rz.Set(0,0,0);
		rxx.Set(0,0,0);
		ryx.Set(0,0,0);
		rzx.Set(0,0,0);

		for (int i = 1; i <= NS(); i++)
		{
			rx.X() += Sx(i)*xgd(3*(i-1)+1);
			rx.Y() += Sx(i)*xgd(3*(i-1)+2);
			rx.Z() += Sx(i)*xgd(3*(i-1)+3);					
			ry.X() += Sy(i)*xgd(3*(i-1)+1);
			ry.Y() += Sy(i)*xgd(3*(i-1)+2);
			ry.Z() += Sy(i)*xgd(3*(i-1)+3);					
			rz.X() += Sz(i)*xgd(3*(i-1)+1);
			rz.Y() += Sz(i)*xgd(3*(i-1)+2);
			rz.Z() += Sz(i)*xgd(3*(i-1)+3);					
			rxx.X() += Sxx(i)*xgd(3*(i-1)+1);
			rxx.Y() += Sxx(i)*xgd(3*(i-1)+2);
			rxx.Z() += Sxx(i)*xgd(3*(i-1)+3);					
			ryx.X() += Syx(i)*xgd(3*(i-1)+1);
			ryx.Y() += Syx(i)*xgd(3*(i-1)+2);
			ryx.Z() += Syx(i)*xgd(3*(i-1)+3);					
			rzx.X() += Szx(i)*xgd(3*(i-1)+1);
			rzx.Y() += Szx(i)*xgd(3*(i-1)+2);
			rzx.Z() += Szx(i)*xgd(3*(i-1)+3);					
		}

		//bending curvatures:
		double kapy=  (rz*rxx);
		double kapz=-(ry*rxx);

		//(mid strain+curvatures*heigth)*Em:
		stress(1,1) = Em*(0.5*(rx*rx-1)+kapz*yb-kapy*zb); //checked sign by comparison!!!

		//shear due to shear forces:

		Vector2D Nf(0.5*(1-x0),0.5*(1+x0)); //shape function for shear ... linear!!!

		double Sxym; 
		double Sxzm; 

		//Sxym = Gm*ky/(A*lx)*(Vector2D(4.*Wxy1f-2.*Wxy2f,-2.*Wxy1f+4.*Wxy2f)*Nf);
		//Sxzm = Gm*kz/(A*lx)*(Vector2D(4.*Wxz1f-2.*Wxz2f,-2.*Wxz1f+4.*Wxz2f)*Nf);


		//alternatively for shear???
		//take shear directly, interpolate linearly ...

		Sxym = Gm*ky*(Vector2D(rx1*ry1,rx2*ry2)*Nf);
		Sxzm = Gm*kz*(Vector2D(rx1*rz1,rx2*rz2)*Nf);

		//linearly interpolate shear parameters computed in Hellinger Reissner?
		//Sxym = Gm*ky/(A*lx)*(Vector2D(2.*Wxy1f, 2.*Wxy2f)*Nf);
		//Sxzm = Gm*kz/(A*lx)*(Vector2D(2.*Wxz1f, 2.*Wxz2f)*Nf);


		stress(1,2) = Sxym*1.5*(Sqr(ly)-4.*Sqr(yb))/(lz*Cub(ly))*A; //put formula for shear for different cross-sections here!!!
		stress(1,3) =-Sxzm*1.5*(Sqr(lz)-4.*Sqr(zb))/(ly*Cub(lz))*A; //Vorzeichen nicht gleich mit General.cpp


		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//Torsion by St. Venant's theory:

		double kapx1	= 0.5 * (rz1*ryx1 - ry1*rzx1);
		double kapx2	= 0.5 * (rz2*ryx2 - ry2*rzx2);

		double Fxy, Fxz, theta;
		//linear interpolation of kapx:
		Fxy = Sxym;
		Fxz = Sxzm;
		theta = (Vector2D(kapx1,kapx2)*Nf); //twist

		double sum12 = 0;
		double sum13 = 0;
		double a = 0.5*ly;
		double b = 0.5*lz;

		double nn = 10; //terms to add from summation formula ... 2 gives 5 digits for maximum strain, but only 1 digit for edge strain
		//10 gives about 2.7% error in edge strain

		for (double n=0; n <= nn; n++)
		{
			sum12 += 16.*pow(-1,n)*a/Sqr(2*n+1)/Sqr(MY_PI)/cosh(0.5*(2*n+1)*MY_PI/a*b)*cos(0.5*(2*n+1)*MY_PI/a*yb)*sinh(0.5*(2*n+1)*MY_PI/a*zb);
			sum13 += 16.*pow(-1,n)*a/Sqr(2*n+1)/Sqr(MY_PI)/cosh(0.5*(2*n+1)*MY_PI/a*b)*sin(0.5*(2*n+1)*MY_PI/a*yb)*cosh(0.5*(2*n+1)*MY_PI/a*zb);
		}

		stress(1,2) +=-Gm*theta*sum12;
		stress(1,3) += Gm*theta*(2.*yb-sum13);

		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		stress(2,1) = stress(1,2);
		stress(3,1) = stress(1,3);

		if(fvd.VariableType() == FieldVariableDescriptor::FVT_stress_mises)
		{
			return stress.Mises();
		}
		else if (fvd.VariableType() == FieldVariableDescriptor::FVT_stress)
		{
			return fvd.GetComponent(stress);
		}
	}

	return FIELD_VARIABLE_NO_VALUE;
}

void ANCFBeam3D::ComputeCorotationalFrame(Vector3D& pref, Matrix3D& Aref)
{

	xg.SetLen(SOS());
	GetCoordinates(xg);
	ApplyT(xg);

	pref = 0.5*Vector3D(xg(13)+xg(1), xg(14)+xg(2), xg(15)+xg(3));

	/*
	Vector3D x1(xg(13)-xg(1), xg(14)-xg(2), xg(15)-xg(3));
	Vector3D y1(xg(7), xg(8), xg(9));
	Vector3D z1(xg(10),xg(11),xg(12));

	x1.Normalize();
	Vector3D y1a = y1;
	Vector3D z1a = z1;
	x1.GramSchmidt(y1);
	y1.Normalize();
	z1 = x1.Cross(y1);

	Aref = Matrix3D(
	x1.X(),y1.X(),z1.X(),
	x1.Y(),y1.Y(),z1.Y(),
	x1.Z(),y1.Z(),z1.Z());
	*/


	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	double q1  = xg(1);  double q2  = xg(2);  double q3  = xg(3);
	double q7  = xg(7);  double q8  = xg(8);  double q9  = xg(9);
	double q13 = xg(13); double q14 = xg(14); double q15 = xg(15);
	double q19 = xg(19); double q20 = xg(20); double q21 = xg(21);

	double t4 =  ((double) q13 * (double) q13);
	double t7 =  ((double) q1 * (double) q1);
	double t8 =  ((double) q14 * (double) q14);
	double t11 =  ((double) q2 * (double) q2);
	double t12 =  ((double) q15 * (double) q15);
	double t15 =  ((double) q3 * (double) q3);
	double t17 = sqrt((double) (t4 - 2 * q13 * q1 + t7 + t8 - 2 * q14 * q2 + t11 + t12 - 2 * q15 * q3 + t15));
	double t18 = 0.1e1 / t17;
	double t19 = q13 - q1;
	double t20 = t18 * (double) t19;
	double t23 = t19 * t19;
	double t24 = q14 - q2;
	double t25 = t24 * t24;
	double t26 = q15 - q3;
	double t27 = t26 * t26;
	double t37 = 0.1e1 / (double) (t23 + t25 + t27) * ((double) ((q7 + q19) * t19) / 0.2e1 + (double) ((q8 + q20) * t24) / 0.2e1 + (double) ((q9 + q21) * t26) / 0.2e1);
	double t39 = (double) q7 / 0.2e1 + (double) q19 / 0.2e1 - t37 * (double) t19;
	double t40 = t39 * t39;
	double t44 = (double) q8 / 0.2e1 + (double) q20 / 0.2e1 - t37 * (double) t24;
	double t45 = t44 * t44;
	double t49 = (double) q9 / 0.2e1 + (double) q21 / 0.2e1 - t37 * (double) t26;
	double t50 = t49 * t49;
	double t52 = sqrt(t40 + t45 + t50);
	double t53 = 0.1e1 / t52;
	double t54 = t53 * t39;
	double t55 = t18 * (double) t24;
	double t56 = t53 * t49;
	double t58 = t18 * (double) t26;
	double t59 = t53 * t44;
	Aref(1,1) = t20;
	Aref(1,2) = t54;
	Aref(1,3) = t55 * t56 - t58 * t59;
	Aref(2,1) = t55;
	Aref(2,2) = t59;
	Aref(2,3) = t58 * t54 - t20 * t56;
	Aref(3,1) = t58;
	Aref(3,2) = t56;
	Aref(3,3) = t20 * t59 - t55 * t54;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

}

void ANCFBeam3D::ComputeCorotationalFrameD(Vector3D& pref, Matrix3D& Aref)
{

	xgd.SetLen(SOS());
	GetDrawCoordinates(xgd);
	xgd = T*xgd;

	pref = 0.5*Vector3D(xgd(13)+xgd(1), xgd(14)+xgd(2), xgd(15)+xgd(3));

	/*
	Vector3D x1(xgd(13)-xgd(1), xgd(14)-xgd(2), xgd(15)-xgd(3));
	Vector3D y1(xgd(7), xgd(8), xgd(9));
	Vector3D z1(xgd(10),xgd(11),xgd(12));

	x1.Normalize();
	Vector3D y1a = y1;
	Vector3D z1a = z1;
	x1.GramSchmidt(y1);
	y1.Normalize();
	z1 = x1.Cross(y1);

	Aref = Matrix3D(
	x1.X(),y1.X(),z1.X(),
	x1.Y(),y1.Y(),z1.Y(),
	x1.Z(),y1.Z(),z1.Z());
	*/


	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	double q1  = xgd(1);  double q2  = xgd(2);  double q3  = xgd(3);
	double q7  = xgd(7);  double q8  = xgd(8);  double q9  = xgd(9);
	double q13 = xgd(13); double q14 = xgd(14); double q15 = xgd(15);
	double q19 = xgd(19); double q20 = xgd(20); double q21 = xgd(21);


	double t4 =  ((double) q13 * (double) q13);
	double t7 =  ((double) q1 * (double) q1);
	double t8 =  ((double) q14 * (double) q14);
	double t11 =  ((double) q2 * (double) q2);
	double t12 =  ((double) q15 * (double) q15);
	double t15 =  ((double) q3 * (double) q3);
	double t17 = sqrt((double) (t4 - 2 * q13 * q1 + t7 + t8 - 2 * q14 * q2 + t11 + t12 - 2 * q15 * q3 + t15));
	double t18 = 0.1e1 / t17;
	double t19 = q13 - q1;
	double t20 = t18 * (double) t19;
	double t23 = t19 * t19;
	double t24 = q14 - q2;
	double t25 = t24 * t24;
	double t26 = q15 - q3;
	double t27 = t26 * t26;
	double t37 = 0.1e1 / (double) (t23 + t25 + t27) * ((double) ((q7 + q19) * t19) / 0.2e1 + (double) ((q8 + q20) * t24) / 0.2e1 + (double) ((q9 + q21) * t26) / 0.2e1);
	double t39 = (double) q7 / 0.2e1 + (double) q19 / 0.2e1 - t37 * (double) t19;
	double t40 = t39 * t39;
	double t44 = (double) q8 / 0.2e1 + (double) q20 / 0.2e1 - t37 * (double) t24;
	double t45 = t44 * t44;
	double t49 = (double) q9 / 0.2e1 + (double) q21 / 0.2e1 - t37 * (double) t26;
	double t50 = t49 * t49;
	double t52 = sqrt(t40 + t45 + t50);
	double t53 = 0.1e1 / t52;
	double t54 = t53 * t39;
	double t55 = t18 * (double) t24;
	double t56 = t53 * t49;
	double t58 = t18 * (double) t26;
	double t59 = t53 * t44;
	Aref(1,1) = t20;
	Aref(1,2) = t54;
	Aref(1,3) = t55 * t56 - t58 * t59;
	Aref(2,1) = t55;
	Aref(2,2) = t59;
	Aref(2,3) = t58 * t54 - t20 * t56;
	Aref(3,1) = t58;
	Aref(3,2) = t56;
	Aref(3,3) = t20 * t59 - t55 * t54;


}

void ANCFBeam3D::ComputeCorotationalFrameDAvdq(Vector3D& v, Matrix& dAvdq)
{

}


void ANCFBeam3D::ComputeCorotationalFrameDATvdq(Vector3D& v, Matrix& dATvdq)
{
	xg.SetLen(SOS());
	GetCoordinates(xg);
	ApplyT(xg);

	dATvdq.SetSize(3,24);
	double v1 = v.X();
	double v2 = v.Y();
	double v3 = v.Z();
	double q1  = xg(1);  double q2  = xg(2);  double q3  = xg(3);
	double q7  = xg(7);  double q8  = xg(8);  double q9  = xg(9);
	double q13 = xg(13); double q14 = xg(14); double q15 = xg(15);
	double q19 = xg(19); double q20 = xg(20); double q21 = xg(21);

	double t4 =  ( q13 *  q13);
	double t7 =  ( q1 *  q1);
	double t8 =  ( q14 *  q14);
	double t11 =  ( q2 *  q2);
	double t12 =  ( q15 *  q15);
	double t15 =  ( q3 *  q3);
	double t16 = t4 - 2 * q13 * q1 + t7 + t8 - 2 * q14 * q2 + t11 + t12 - 2 * q15 * q3 + t15;
	double t17 = sqrt( t16);
	double t19 = 0.1e1 / t17 /  t16;
	double t20 = q13 - q1;
	double t21 = t19 *  t20;
	double t25 = 0.1e1 / t17;
	double t26 = t25 * v1;
	double t27 = q14 - q2;
	double t28 = t19 *  t27;
	double t32 = q15 - q3;
	double t33 = t19 *  t32;
	double t44 = t25 * v2;
	double t58 = t25 * v3;
	double t92 = t20 * t20;
	double t93 = t27 * t27;
	double t94 = t32 * t32;
	double t95 = t92 + t93 + t94;
	double t96 = 1 / t95;
	double t97 =  (q7 + q19);
	double t99 =  (q8 + q20);
	double t101 =  (q9 + q21);
	double t103 =  (t97 * t20) / 0.2e1 +  (t99 * t27) / 0.2e1 +  (t101 * t32) / 0.2e1;
	double t104 =  t96 * t103;
	double t106 = q7 / 0.2e1 + q19 / 0.2e1 - t104 *  t20;
	double t107 = t106 * t106;
	double t111 = q8 / 0.2e1 + q20 / 0.2e1 - t104 *  t27;
	double t112 = t111 * t111;
	double t116 = q9 / 0.2e1 + q21 / 0.2e1 - t104 *  t32;
	double t117 = t116 * t116;
	double t118 = t107 + t112 + t117;
	double t119 = sqrt(t118);
	double t121 = 0.1e1 / t119 / t118;
	double t122 = t121 * t106;
	double t123 = t95 * t95;
	double t125 = 0.1e1 /  t123 * t103;
	double t128 = - (t96 * t97) / 0.2e1;
	double t130 = -0.2e1 * t125 *  (t20 * t20) - t128 *  t20 + t104;
	double t133 = -0.2e1 * t125 *  t27 *  t20;
	double t135 = t133 - t128 *  t27;
	double t138 = -0.2e1 * t125 *  t32 *  t20;
	double t140 = t138 - t128 *  t32;
	double t142 = t106 * t130 + t111 * t135 + t116 * t140;
	double t146 = 0.1e1 / t119;
	double t147 = t146 * t130;
	double t149 = t121 * t111;
	double t153 = t146 * t135;
	double t155 = t121 * t116;
	double t159 = t146 * t140;
	double t162 = - (t96 * t99) / 0.2e1;
	double t164 = t133 - t162 *  t20;
	double t169 = -0.2e1 * t125 *  (t27 * t27) - t162 *  t27 + t104;
	double t172 = -0.2e1 * t125 *  t32 *  t27;
	double t174 = t172 - t162 *  t32;
	double t176 = t106 * t164 + t111 * t169 + t116 * t174;
	double t180 = t146 * t164;
	double t185 = t146 * t169;
	double t190 = t146 * t174;
	double t193 = - (t96 * t101) / 0.2e1;
	double t195 = t138 - t193 *  t20;
	double t198 = t172 - t193 *  t27;
	double t203 = -0.2e1 * t125 *  (t32 * t32) - t193 *  t32 + t104;
	double t205 = t106 * t195 + t111 * t198 + t116 * t203;
	double t209 = t146 * t195;
	double t214 = t146 * t198;
	double t219 = t146 * t203;
	double t222 =  (t96 * t20) / 0.2e1;
	double t224 = 0.1e1 / 0.2e1 - t222 *  t20;
	double t226 = t111 *  t96;
	double t227 =  (t27 * t20) / 0.2e1;
	double t229 = t116 *  t96;
	double t230 =  (t32 * t20) / 0.2e1;
	double t232 = t106 * t224 - t226 * t227 - t229 * t230;
	double t236 = t146 * t224;
	double t241 = t146 *  t96;
	double t249 = -t122 * v1 * t232 + t236 * v1 - t149 * v2 * t232 - t241 * t227 * v2 - t155 * v3 * t232 - t241 * t230 * v3;
	double t250 = t106 *  t96;
	double t252 =  (t96 * t27) / 0.2e1;
	double t254 = 0.1e1 / 0.2e1 - t252 *  t27;
	double t256 =  (t32 * t27) / 0.2e1;
	double t258 = -t250 * t227 + t111 * t254 - t229 * t256;
	double t267 = t146 * t254;
	double t274 = -t122 * v1 * t258 - t241 * t227 * v1 - t149 * v2 * t258 + t267 * v2 - t155 * v3 * t258 - t241 * t256 * v3;
	double t277 =  (t96 * t32) / 0.2e1;
	double t279 = 0.1e1 / 0.2e1 - t277 *  t32;
	double t281 = -t250 * t230 - t226 * t256 + t116 * t279;
	double t295 = t146 * t279;
	double t297 = -t122 * v1 * t281 - t241 * t230 * v1 - t149 * v2 * t281 - t241 * t256 * v2 - t155 * v3 * t281 + t295 * v3;
	double t300 =  (t96 * t97) / 0.2e1;
	double t302 = 0.2e1 * t125 *  (t20 * t20) - t300 *  t20 - t104;
	double t305 = 0.2e1 * t125 *  t27 *  t20;
	double t307 = t305 - t300 *  t27;
	double t310 = 0.2e1 * t125 *  t32 *  t20;
	double t312 = t310 - t300 *  t32;
	double t314 = t106 * t302 + t111 * t307 + t116 * t312;
	double t318 = t146 * t302;
	double t323 = t146 * t307;
	double t328 = t146 * t312;
	double t331 =  (t96 * t99) / 0.2e1;
	double t333 = t305 - t331 *  t20;
	double t338 = 0.2e1 * t125 *  (t27 * t27) - t331 *  t27 - t104;
	double t341 = 0.2e1 * t125 *  t32 *  t27;
	double t343 = t341 - t331 *  t32;
	double t345 = t106 * t333 + t111 * t338 + t116 * t343;
	double t349 = t146 * t333;
	double t354 = t146 * t338;
	double t359 = t146 * t343;
	double t362 =  (t96 * t101) / 0.2e1;
	double t364 = t310 - t362 *  t20;
	double t367 = t341 - t362 *  t27;
	double t372 = 0.2e1 * t125 *  (t32 * t32) - t362 *  t32 - t104;
	double t374 = t106 * t364 + t111 * t367 + t116 * t372;
	double t378 = t146 * t364;
	double t383 = t146 * t367;
	double t388 = t146 * t372;
	double t391 = t146 * t116;
	double t392 = -0.2e1 * t391 *  t20;
	double t395 = t25 *  t27;
	double t396 = 0.2e1 * t155 * t142;
	double t400 = t146 * t111;
	double t401 = -0.2e1 * t400 *  t20;
	double t404 = t25 *  t32;
	double t405 = 0.2e1 * t149 * t142;
	double t411 = t146 * t106;
	double t412 = -0.2e1 * t411 *  t20;
	double t415 = 0.2e1 * t122 * t142;
	double t421 = t25 * t146;
	double t422 = t421 * t116;
	double t423 = t25 *  t20;
	double t431 = t421 * t111;
	double t443 = -0.2e1 * t391 *  t27;
	double t446 = 0.2e1 * t155 * t176;
	double t450 = -0.2e1 * t400 *  t27;
	double t453 = 0.2e1 * t149 * t176;
	double t459 = -0.2e1 * t411 *  t27;
	double t462 = 0.2e1 * t122 * t176;
	double t480 = t421 * t106;
	double t487 = -0.2e1 * t391 *  t32;
	double t490 = 0.2e1 * t155 * t205;
	double t494 = -0.2e1 * t400 *  t32;
	double t497 = 0.2e1 * t149 * t205;
	double t503 = -0.2e1 * t411 *  t32;
	double t506 = 0.2e1 * t122 * t205;
	double t530 = 0.2e1 * t155 * t232;
	double t532 = 0.2e1 * t149 * t232;
	double t536 = 0.2e1 * t122 * t232;
	double t542 = t423 * t146;
	double t556 = (-t395 * t530 + t404 * t532) * v1 / 0.2e1 + (-t404 * t536 / 0.2e1 + t404 * t236 + t423 * t530 / 0.2e1 + t542 * t222 *  t32) * v2 + (-t423 * t532 / 0.2e1 - t542 * t222 *  t27 + t395 * t536 / 0.2e1 - t395 * t236) * v3;
	double t557 = 0.2e1 * t155 * t258;
	double t560 = t395 * t146;
	double t563 = 0.2e1 * t149 * t258;
	double t569 = 0.2e1 * t122 * t258;
	double t583 = (-t395 * t557 / 0.2e1 - t560 * t252 *  t32 + t404 * t563 / 0.2e1 - t404 * t267) * v1 + (-t404 * t569 + t423 * t557) * v2 / 0.2e1 + (-t423 * t563 / 0.2e1 + t423 * t267 + t395 * t569 / 0.2e1 + t560 * t252 *  t20) * v3;
	double t584 = 0.2e1 * t155 * t281;
	double t588 = 0.2e1 * t149 * t281;
	double t591 = t404 * t146;
	double t596 = 0.2e1 * t122 * t281;
	double t610 = (-t395 * t584 / 0.2e1 + t395 * t295 + t404 * t588 / 0.2e1 + t591 * t277 *  t27) * v1 + (-t404 * t596 / 0.2e1 - t591 * t277 *  t20 + t423 * t584 / 0.2e1 - t423 * t295) * v2 + (-t423 * t588 + t395 * t596) * v3 / 0.2e1;
	double t611 = 0.2e1 * t391 *  t20;
	double t614 = 0.2e1 * t155 * t314;
	double t618 = 0.2e1 * t400 *  t20;
	double t621 = 0.2e1 * t149 * t314;
	double t627 = 0.2e1 * t411 *  t20;
	double t630 = 0.2e1 * t122 * t314;
	double t654 = 0.2e1 * t391 *  t27;
	double t657 = 0.2e1 * t155 * t345;
	double t661 = 0.2e1 * t400 *  t27;
	double t664 = 0.2e1 * t149 * t345;
	double t670 = 0.2e1 * t411 *  t27;
	double t673 = 0.2e1 * t122 * t345;
	double t697 = 0.2e1 * t391 *  t32;
	double t700 = 0.2e1 * t155 * t374;
	double t704 = 0.2e1 * t400 *  t32;
	double t707 = 0.2e1 * t149 * t374;
	double t713 = 0.2e1 * t411 *  t32;
	double t716 = 0.2e1 * t122 * t374;
	dATvdq(0+1,1+0 ) = t21 * v1 *  t20 - t26 + t28 * v2 *  t20 + t33 * v3 *  t20;
	dATvdq(0+1,1+1 ) = t21 * v1 *  t27 + t28 * v2 *  t27 - t44 + t33 * v3 *  t27;
	dATvdq(0+1,1+2 ) = t21 * v1 *  t32 + t28 * v2 *  t32 + t33 * v3 *  t32 - t58;
	dATvdq(0+1,1+3 ) = 0.0e0;
	dATvdq(0+1,1+4 ) = 0.0e0;
	dATvdq(0+1,1+5 ) = 0.0e0;
	dATvdq(0+1,1+6 ) = 0.0e0;
	dATvdq(0+1,1+7 ) = 0.0e0;
	dATvdq(0+1,1+8 ) = 0.0e0;
	dATvdq(0+1,1+9 ) = 0.0e0;
	dATvdq(0+1,1+10) = 0.0e0;
	dATvdq(0+1,1+11) = 0.0e0;
	dATvdq(0+1,1+12) = -t21 * v1 *  t20 + t26 - t28 * v2 *  t20 - t33 * v3 *  t20;
	dATvdq(0+1,1+13) = -t21 * v1 *  t27 - t28 * v2 *  t27 + t44 - t33 * v3 *  t27;
	dATvdq(0+1,1+14) = -t21 * v1 *  t32 - t28 * v2 *  t32 - t33 * v3 *  t32 + t58;
	dATvdq(0+1,1+15) = 0.0e0;
	dATvdq(0+1,1+16) = 0.0e0;
	dATvdq(0+1,1+17) = 0.0e0;
	dATvdq(0+1,1+18) = 0.0e0;
	dATvdq(0+1,1+19) = 0.0e0;
	dATvdq(0+1,1+20) = 0.0e0;
	dATvdq(0+1,1+21) = 0.0e0;
	dATvdq(0+1,1+22) = 0.0e0;
	dATvdq(0+1,1+23) = 0.0e0;
	dATvdq(1+1,1+0 ) = -t122 * v1 * t142 + t147 * v1 - t149 * v2 * t142 + t153 * v2 - t155 * v3 * t142 + t159 * v3;
	dATvdq(1+1,1+1 ) = -t122 * v1 * t176 + t180 * v1 - t149 * v2 * t176 + t185 * v2 - t155 * v3 * t176 + t190 * v3;
	dATvdq(1+1,1+2 ) = -t122 * v1 * t205 + t209 * v1 - t149 * v2 * t205 + t214 * v2 - t155 * v3 * t205 + t219 * v3;
	dATvdq(1+1,1+3 ) = 0.0e0;
	dATvdq(1+1,1+4 ) = 0.0e0;
	dATvdq(1+1,1+5 ) = 0.0e0;
	dATvdq(1+1,1+6 ) = t249;
	dATvdq(1+1,1+7 ) = t274;
	dATvdq(1+1,1+8 ) = t297;
	dATvdq(1+1,1+9 ) = 0.0e0;
	dATvdq(1+1,1+10) = 0.0e0;
	dATvdq(1+1,1+11) = 0.0e0;
	dATvdq(1+1,1+12) = -t122 * v1 * t314 + t318 * v1 - t149 * v2 * t314 + t323 * v2 - t155 * v3 * t314 + t328 * v3;
	dATvdq(1+1,1+13) = -t122 * v1 * t345 + t349 * v1 - t149 * v2 * t345 + t354 * v2 - t155 * v3 * t345 + t359 * v3;
	dATvdq(1+1,1+14) = -t122 * v1 * t374 + t378 * v1 - t149 * v2 * t374 + t383 * v2 - t155 * v3 * t374 + t388 * v3;
	dATvdq(1+1,1+15) = 0.0e0;
	dATvdq(1+1,1+16) = 0.0e0;
	dATvdq(1+1,1+17) = 0.0e0;
	dATvdq(1+1,1+18) = t249;
	dATvdq(1+1,1+19) = t274;
	dATvdq(1+1,1+20) = t297;
	dATvdq(1+1,1+21) = 0.0e0;
	dATvdq(1+1,1+22) = 0.0e0;
	dATvdq(1+1,1+23) = 0.0e0;
	dATvdq(2+1,1+0 ) = (-t28 * t392 / 0.2e1 - t395 * t396 / 0.2e1 + t395 * t159 + t33 * t401 / 0.2e1 + t404 * t405 / 0.2e1 - t404 * t153) * v1 + (-t33 * t412 / 0.2e1 - t404 * t415 / 0.2e1 + t404 * t147 + t21 * t392 / 0.2e1 + t422 + t423 * t396 / 0.2e1 - t423 * t159) * v2 + (-t21 * t401 / 0.2e1 - t431 - t423 * t405 / 0.2e1 + t423 * t153 + t28 * t412 / 0.2e1 + t395 * t415 / 0.2e1 - t395 * t147) * v3;
	dATvdq(2+1,1+1 ) = (-t28 * t443 / 0.2e1 - t422 - t395 * t446 / 0.2e1 + t395 * t190 + t33 * t450 / 0.2e1 + t404 * t453 / 0.2e1 - t404 * t185) * v1 + (-t33 * t459 / 0.2e1 - t404 * t462 / 0.2e1 + t404 * t180 + t21 * t443 / 0.2e1 + t423 * t446 / 0.2e1 - t423 * t190) * v2 + (-t21 * t450 / 0.2e1 - t423 * t453 / 0.2e1 + t423 * t185 + t28 * t459 / 0.2e1 + t480 + t395 * t462 / 0.2e1 - t395 * t180) * v3;
	dATvdq(2+1,1+2 ) = (-t28 * t487 / 0.2e1 - t395 * t490 / 0.2e1 + t395 * t219 + t33 * t494 / 0.2e1 + t431 + t404 * t497 / 0.2e1 - t404 * t214) * v1 + (-t33 * t503 / 0.2e1 - t480 - t404 * t506 / 0.2e1 + t404 * t209 + t21 * t487 / 0.2e1 + t423 * t490 / 0.2e1 - t423 * t219) * v2 + (-t21 * t494 / 0.2e1 - t423 * t497 / 0.2e1 + t423 * t214 + t28 * t503 / 0.2e1 + t395 * t506 / 0.2e1 - t395 * t209) * v3;
	dATvdq(2+1,1+3 ) = 0.0e0;
	dATvdq(2+1,1+4 ) = 0.0e0;
	dATvdq(2+1,1+5 ) = 0.0e0;
	dATvdq(2+1,1+6 ) = t556;
	dATvdq(2+1,1+7 ) = t583;
	dATvdq(2+1,1+8 ) = t610;
	dATvdq(2+1,1+9 ) = 0.0e0;
	dATvdq(2+1,1+10) = 0.0e0;
	dATvdq(2+1,1+11) = 0.0e0;
	dATvdq(2+1,1+12) = (-t28 * t611 / 0.2e1 - t395 * t614 / 0.2e1 + t395 * t328 + t33 * t618 / 0.2e1 + t404 * t621 / 0.2e1 - t404 * t323) * v1 + (-t33 * t627 / 0.2e1 - t404 * t630 / 0.2e1 + t404 * t318 + t21 * t611 / 0.2e1 - t422 + t423 * t614 / 0.2e1 - t423 * t328) * v2 + (-t21 * t618 / 0.2e1 + t431 - t423 * t621 / 0.2e1 + t423 * t323 + t28 * t627 / 0.2e1 + t395 * t630 / 0.2e1 - t395 * t318) * v3;
	dATvdq(2+1,1+13) = (-t28 * t654 / 0.2e1 + t422 - t395 * t657 / 0.2e1 + t395 * t359 + t33 * t661 / 0.2e1 + t404 * t664 / 0.2e1 - t404 * t354) * v1 + (-t33 * t670 / 0.2e1 - t404 * t673 / 0.2e1 + t404 * t349 + t21 * t654 / 0.2e1 + t423 * t657 / 0.2e1 - t423 * t359) * v2 + (-t21 * t661 / 0.2e1 - t423 * t664 / 0.2e1 + t423 * t354 + t28 * t670 / 0.2e1 - t480 + t395 * t673 / 0.2e1 - t395 * t349) * v3;
	dATvdq(2+1,1+14) = (-t28 * t697 / 0.2e1 - t395 * t700 / 0.2e1 + t395 * t388 + t33 * t704 / 0.2e1 - t431 + t404 * t707 / 0.2e1 - t404 * t383) * v1 + (-t33 * t713 / 0.2e1 + t480 - t404 * t716 / 0.2e1 + t404 * t378 + t21 * t697 / 0.2e1 + t423 * t700 / 0.2e1 - t423 * t388) * v2 + (-t21 * t704 / 0.2e1 - t423 * t707 / 0.2e1 + t423 * t383 + t28 * t713 / 0.2e1 + t395 * t716 / 0.2e1 - t395 * t378) * v3;
	dATvdq(2+1,1+15) = 0.0e0;
	dATvdq(2+1,1+16) = 0.0e0;
	dATvdq(2+1,1+17) = 0.0e0;
	dATvdq(2+1,1+18) = t556;
	dATvdq(2+1,1+19) = t583;
	dATvdq(2+1,1+20) = t610;
	dATvdq(2+1,1+21) = 0.0e0;
	dATvdq(2+1,1+22) = 0.0e0;
	dATvdq(2+1,1+23) = 0.0e0;
}


void ANCFBeam3D::DrawElement() 
{
	mbs->SetColor(col);

	//mbs->uout << "ANCFBeam3D: m=" << mass << ", I=" << Iphi << ", Size=" << size << "\n";
	//GetMBS()->UO() << "r_posD=" << GetPosD(Vector3D(1,1,1)) << "\n";
	//char str[100];
	//Vector3D p0 = GetPosD(Vector3D(-1,0,0));
	//sprintf(str,"    posx=%6.3g,posy=%6.3g,posz=%6.3g  ", p0.X(), p0.Y(), p0.Z());
	//GetMBS()->GetRC()->PrintTextStruct(5,-1,str);

	//double lx = -0.5*size.X(); double ly = 0.5*size.Y(); double lz = 0.5*size.Z();
	double lx1 = 1; double ly1 = 1*GetMBS()->GetMagnifyYZ(); double lz1 = 1*GetMBS()->GetMagnifyYZ();

	double def_scale = GetMBS()->GetDOption(105); //deformation scaling


	if(0) //draw underlying rigid body configuration:
	{
		Vector3D pmid;
		Matrix3D Rot;
		ComputeCorotationalFrameD(pmid, Rot);
		Vector3D pmid0(x_init(13)+x_init(1),x_init(14)+x_init(2),x_init(15)+x_init(3));
		pmid0*=0.5;

		//Matrix3D A;
		//A = GetRotMatrixD(Vector3D(0.,0.,0.));


		Vector3D p[10];
		double l1 = -0.5*lx;
		double l2 =  0.5*lx;
		double thickness = 1;
		double ly1 = 0.5*ly;
		double lz1 = 0.5*lz;

		p[8] = Vector3D(l1,-ly1,-lz1);
		p[7] = Vector3D(l1,-ly1, lz1);
		p[4] = Vector3D(l1, ly1,-lz1); 
		p[3] = Vector3D(l1, ly1, lz1);
		p[6] = Vector3D(l2,-ly1,-lz1);
		p[5] = Vector3D(l2,-ly1, lz1);
		p[2] = Vector3D(l2, ly1,-lz1);
		p[1] = Vector3D(l2, ly1, lz1);

		for (int i=1; i <= 8; i++) p[i] = def_scale*(pmid-pmid0)+pmid0 + (def_scale*(Rot-Matrix3D(1.))+Matrix3D(1.))*p[i];

		mbs->MyDrawLine(p[8],p[7],thickness);
		mbs->MyDrawLine(p[7],p[3],thickness);
		mbs->MyDrawLine(p[4],p[3],thickness);
		mbs->MyDrawLine(p[4],p[8],thickness);
		mbs->MyDrawLine(p[6],p[5],thickness);
		mbs->MyDrawLine(p[5],p[1],thickness);
		mbs->MyDrawLine(p[2],p[1],thickness);
		mbs->MyDrawLine(p[2],p[6],thickness);
		mbs->MyDrawLine(p[6],p[8],thickness);
		mbs->MyDrawLine(p[5],p[7],thickness);
		mbs->MyDrawLine(p[2],p[4],thickness);
		mbs->MyDrawLine(p[1],p[3],thickness);

	}

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
	else
	{
		colormode = 0; 
	}

	if (colormode)
	{

		//draw modes:
		//for (int ii = 1; ii <= 2; ii++)
		{
			double modeval = 0;
			int xgset = 0;

			double tilex = GetMBS()->GetIOption(137);
			double tiley = GetMBS()->GetIOption(138);
			//if (!colormode) {tiley = 1; tilex = 20;}
			TArray<Vector3D> points((int)(tilex+1)*(int)(tiley+1));
			TArray<double> vals((int)(tilex+1)*(int)(tiley+1));
			double v=0;

			for (int side=1; side <= 6; side++)
			{
				points.SetLen(0); vals.SetLen(0);
				Vector3D p0, vx, vy;
				int tileyn = (int)tiley;
				int tilexn = (int)tilex;

				/*if (side >= 5) 
				{
				tilexn = 2;
				tileyn = 2;
				}*/

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
						Vector3D ploc = (p0+ix*vx+iy*vy);
						Vector3D pg = GetPos0D(ploc, def_scale );
						Vector3D p(pg);
						points.Add(p);
						if (colormode)
							v = GetFieldVariableValue(*GetMBS()->GetActualPostProcessingFieldVariable(), ploc, true);
						vals.Add(v);
					}
				}
				mbs->DrawColorQuads(points,vals,(int)tilexn+1,(int)tileyn+1,colormode,linemode);
			}

			/*
			if (ii == 2) 
			{
			int nmode = GetMBS()->GetDrawResolution()+1;
			if (nmode >= 1 && nmode <= 12)
			{
			XGD((nmode)) -= modeval;
			}
			}*/

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
		//lx1*=0.9;ly1*=0.9;lz1*=0.9;
		Vector3D p1,p2,p3,p4,p5,p6,p7,p8;
		for (double i = 0; i < tiling; i++)
		{
			double l1 = -lx1+2*lx1*i/tiling;
			double l2 = -lx1+2*lx1*(i+1)/tiling;
			if (i == 0)
			{
				p8 = Vector3D(GetPos0D(Vector3D(l1,-ly1,-lz1)));
				p7 = Vector3D(GetPos0D(Vector3D(l1,-ly1, lz1)));
				p4 = Vector3D(GetPos0D(Vector3D(l1, ly1,-lz1))); 
				p3 = Vector3D(GetPos0D(Vector3D(l1, ly1, lz1)));
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
			p6 = Vector3D(GetPos0D(Vector3D(l2,-ly1,-lz1)));
			p5 = Vector3D(GetPos0D(Vector3D(l2,-ly1, lz1)));
			p2 = Vector3D(GetPos0D(Vector3D(l2, ly1,-lz1)));
			p1 = Vector3D(GetPos0D(Vector3D(l2, ly1, lz1)));
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
		mbs->SetDrawlines(0);
	}

	if (concentratedmass1 != 0) 
	{
		Vector3D mp1(GetPos0D(Vector3D(-lx1,0,0))); 
		mbs->SetColor(colred);
		//double r = pow(3.*concentratedmass1/(4.*MY_PI*GetRho()),1./3.); 
		double r = ly/10.;	//$ DR 2013-02-04 deleted rho from class element, do not use it here!
		//double r = ly*GetMBS()->GetMagnifyYZ()*4.;
		mbs->DrawSphere(mp1,r,12);
	}
	if (concentratedmass2 != 0) 
	{
		Vector3D mp1(GetPos0D(Vector3D(lx1,0,0))); 
		mbs->SetColor(colred);
		//double r = pow(3.*concentratedmass2/(4.*MY_PI*GetRho()),1./3.); 
		double r = ly/10.;	//$ DR 2013-02-04 deleted rho from class element, do not use it here!
		//double r = ly*GetMBS()->GetMagnifyYZ()*4.;
		mbs->DrawSphere(mp1,r,12);
	}
};


void ANCFBeam3D::GetdPosdqT(const Vector3D& ploc, Matrix& d)
{
	//p = S(p.x,p.y,p.z)*q; d(p)/dq 
	Vector3D p0=ploc;
	p0.Scale(0.5*lx,0.5*ly,0.5*lz);
	static Vector SV;
	GetS0(SV, p0);
	d.SetSize(NS()*Dim(),3);
	d.FillWithZeros();
	for (int i = 1; i <= NS(); i++)
	{
		d((i-1)*3+1,1)=SV(i);
		d((i-1)*3+2,2)=SV(i);
		d((i-1)*3+3,3)=SV(i);
	}
	ApplyTtp(d);
}

void ANCFBeam3D::GetdRotvdqT(const Vector3D& vloc, const Vector3D& ploc, Matrix& d)
{
	d.SetSize(NS()*Dim(),3);
	/*

	static Matrix d2;
	d2.SetSize(NS()*Dim(),3);
	double diffpar = 1e-8;
	GetdPosdqT(ploc+diffpar*vloc, d2);
	d = d2;
	GetdPosdqT(ploc,d2);
	d -= d2;
	d *= 1./diffpar;
	*/

	Vector3D p0(ploc);
	p0.Scale(0.5*lx,0.5*ly,0.5*lz);

	static Vector xg;
	xg.SetLen(SOS());
	GetCoordinates(xg);
	ApplyT(xg);

	if (usetangentframe) GetTangentFramedRotvdqT(vloc, p0.X(), xg, d);
	else GetCrosssectionFramedRotvdqT(vloc, p0.X(), xg, d);

	ApplyTtp(d); //****????
}

void ANCFBeam3D::GetdRotdqT(const Vector3D& ploc, Matrix& d)
{
	//UO() << "dRotdq\n";

	Vector3D p0=ploc;
	p0.Scale(0.5*lx,0.5*ly,0.5*lz);

	d.SetSize(NS()*Dim(),3);
	d.FillWithZeros(); //needed???

	if (0) //old model: inconsistent with GetRotMatrix / TangentFrame !!!
	{
		static Vector Sy; Sy.SetLen(NS()); 
		static Vector Sz;	Sz.SetLen(NS()); 

		static Vector xg;
		xg.SetLen(SOS());
		GetCoordinates(xg);
		//ApplyT(xg);

		Vector3D ry, rz;

		GetS0y(Sy,p0);
		GetS0z(Sz,p0);

		for (int i = 1; i <= NS(); i++)
		{
			ry.X() += Sy(i)*xg(3*(i-1)+1);
			ry.Y() += Sy(i)*xg(3*(i-1)+2);
			ry.Z() += Sy(i)*xg(3*(i-1)+3);					
			rz.X() += Sz(i)*xg(3*(i-1)+1);
			rz.Y() += Sz(i)*xg(3*(i-1)+2);
			rz.Z() += Sz(i)*xg(3*(i-1)+3);					
		}

		double dyx = Sqr(ry.X())+Sqr(ry.Y());
		double dzx = Sqr(rz.X())+Sqr(rz.Z());
		double dyz = Sqr(ry.Y())+Sqr(ry.Z());
		double dzy = Sqr(rz.Y())+Sqr(rz.Z());

		for (int i = 1; i <= NS(); i++)
		{
			//from delta gamma_x
			d((i-1)*3+1,1) = 0; //only delta r.X() terms!
			d((i-1)*3+2,1) =-(ry.Z()*Sy(i))/dyz - (rz.Z()*Sz(i))/dzy; //only delta r.Y() terms!
			d((i-1)*3+3,1) = (ry.Y()*Sy(i))/dyz + (rz.Y()*Sz(i))/dzy; //only delta r.Z() terms!

			//from delta gamma_y
			d((i-1)*3+1,2) = (rz.Z()*Sz(i))/dzx; //only delta r.X() terms!
			d((i-1)*3+2,2) = 0; //only delta r.Y() terms!
			d((i-1)*3+3,2) =-(rz.X()*Sz(i))/dzx; //only delta r.Z() terms!

			//from delta gamma_z
			d((i-1)*3+1,3) =-(ry.Y()*Sy(i))/dyx; //only delta r.X() terms!
			d((i-1)*3+2,3) = (ry.X()*Sy(i))/dyx; //only delta r.Y() terms!
			d((i-1)*3+3,3) = 0; //only delta r.Z() terms!

		}
		//ApplyTtp(d);
	}
	else
	{
		//alternatively:
		ConstMatrix<72> dAvdq(24,3);
		double x = ploc.X();

		Matrix3D AT = GetRotMatrix(ploc).GetTp();

		Vector3D v;
		v = Vector3D(AT(1,1),AT(2,1),AT(3,1));
		GetdRotvdqT(v, ploc, dAvdq);
		for (int j=1; j <= SOS(); j++)
		{
			d(j,3) = dAvdq(j,2);
			d(j,2) =-dAvdq(j,3);
		}
		v = Vector3D(AT(1,2),AT(2,2),AT(3,2));
		GetdRotvdqT(v, ploc, dAvdq);
		for (int j=1; j <= SOS(); j++)
		{
			d(j,1) = dAvdq(j,3);
		}
		//UO() << "d_new=" << d << "\n";
	}
}

Matrix3D ANCFBeam3D::GetRotMatrix(const Vector3D& ploc) const 
{
	Vector3D p0(ploc);
	p0.Scale(0.5*lx,0.5*ly,0.5*lz);
	/*
	//Compute Gradient ...
	static Vector u;
	u.SetLen(SOS());
	GetCoordinates(u);

	ApplyT(u);
	u -= e0;
	Matrix3D G;
	Gradu(p0, u, G);
	G += Matrix3D(1);

	//use cross-section frame, use average of ry and rz for rotation
	Vector3D ry(G(1,2),G(2,2),G(3,2));
	Vector3D rz(G(1,3),G(2,3),G(3,3));
	ry.Normalize();
	rz.Normalize();
	Vector3D ry2=ry;
	Vector3D rz2=rz;
	ry.GramSchmidt(rz2);
	rz = 0.5*(rz+rz2);
	rz.GramSchmidt(ry);
	Vector3D rx = ry.Cross(rz);

	Matrix3D A1(rx.X(),ry.X(),rz.X(),
	rx.Y(),ry.Y(),rz.Y(),
	rx.Z(),ry.Z(),rz.Z());

	return A1;
	*/
	Matrix3D A;

	static Vector xg;
	xg.SetLen(SOS());
	GetCoordinates(xg);
	ApplyT(xg);

	double x = p0.X();


	return ComputeTangentFrame(x, xg);

}

Matrix3D ANCFBeam3D::GetRotMatrixP(const Vector3D& ploc) const 
{
	/*
	//Compute Gradient ...
	static Vector u;
	u.SetLen(SOS());
	GetCoordinatesP(u);

	ApplyT(u);
	Vector3D p0(ploc);
	p0.Scale(0.5*lx,0.5*ly,0.5*lz);
	Matrix3D G;
	Gradu(p0, u, G);
	return G;
	*/

	Vector3D p0(ploc);
	p0.Scale(0.5*lx,0.5*ly,0.5*lz);

	static Vector xg;
	xg.SetLen(SOS());
	GetCoordinates(xg);
	ApplyT(xg);

	static Vector xgp;
	xgp.SetLen(SOS());
	GetCoordinatesP(xgp);
	ApplyT(xgp);

	double x = p0.X();

	return ComputeTangentFrameP(x, xg, xgp);//JG: 2012 falsch, hier sollte cross section frame 
}	

Matrix3D ANCFBeam3D::GetRotMatrixD(const Vector3D& ploc) const 
{	
	Vector3D p0(ploc);
	p0.Scale(0.5*lx,0.5*ly,0.5*lz);
	/*
	//Compute Gradient ...
	static Vector u;
	u.SetLen(SOS());
	GetDrawCoordinates(u);

	ApplyTD(u);
	u -= e0;
	Matrix3D G;
	GraduD(p0, u, G);
	G += Matrix3D(1);

	//use cross-section frame, use average of ry and rz for rotation
	Vector3D ry(G(1,2),G(2,2),G(3,2));
	Vector3D rz(G(1,3),G(2,3),G(3,3));
	ry.Normalize();
	rz.Normalize();

	Vector3D ry2=ry;
	Vector3D rz2=rz;
	ry.GramSchmidt(rz2);
	rz = 0.5*(rz+rz2);
	rz.GramSchmidt(ry);
	Vector3D rx = ry.Cross(rz);

	Matrix3D A1(rx.X(),ry.X(),rz.X(),
	rx.Y(),ry.Y(),rz.Y(),
	rx.Z(),ry.Z(),rz.Z());
	*/

	//the following is the tangent frame
	//Reason: The rotation of the cross-section is only interpolated linearly
	Matrix3D A;

	static Vector xgd;
	xgd.SetLen(SOS());
	GetCoordinates(xgd);
	ApplyT(xgd);

	double q1  = xgd(1);  double q2  = xgd(2);  double q3  = xgd(3);
	double q4  = xgd(4);  double q5  = xgd(5);  double q6  = xgd(6);
	double q7  = xgd(7);  double q8  = xgd(8);  double q9  = xgd(9);
	double q10 = xgd(10); double q11 = xgd(11); double q12 = xgd(12);
	double q13 = xgd(13); double q14 = xgd(14); double q15 = xgd(15);
	double q16 = xgd(16); double q17 = xgd(17); double q18 = xgd(18);
	double q19 = xgd(19); double q20 = xgd(20); double q21 = xgd(21);
	double q22 = xgd(22); double q23 = xgd(23); double q24 = xgd(24);

	double x = p0.X();

	double t4 = 1 /  lx;
	double t5 =  (x * x);
	double t6 = -1 + t5;
	double t7 = 0.3e1 / 0.4e1 * (double) t4 * (double) t6;
	double t12 = -0.1e1 / 0.4e1 - x / 0.2e1 + 0.3e1 / 0.4e1 * (double) t5;
	double t14 = -0.3e1 / 0.4e1 * (double) t4 * (double) t6;
	double t17 = 0.1e1 + x;
	double t19 = -0.1e1 + x;
	double t22 = t17 * t17;
	double t26 = (double) t4 * (lx * t17 * t19 / 0.4e1 + lx * t22 / 0.8e1);
	double t29 = 0.2e1 * t7 * q1 + t12 * q4 + 0.2e1 * t14 * q13 + 0.2e1 * t26 * q16;
	double t30 = t29 * t29;
	double t38 = 0.2e1 * t7 * q2 + t12 * q5 + 0.2e1 * t14 * q14 + 0.2e1 * t26 * q17;
	double t39 = t38 * t38;
	double t47 = 0.2e1 * t7 * q3 + t12 * q6 + 0.2e1 * t14 * q15 + 0.2e1 * t26 * q18;
	double t48 = t47 * t47;
	double t49 = t30 + t39 + t48;
	double t50 = sqrt(t49);
	double t51 = 0.1e1 / t50;
	double t52 = t51 * t29;
	double t53 = -t19 * q7 / 0.2e1;
	double t54 = t17 * q19 / 0.2e1;
	double t55 = 0.1e1 / t49;
	double t58 = -t19 * q8 / 0.2e1;
	double t59 = t17 * q20 / 0.2e1;
	double t62 = -t19 * q9 / 0.2e1;
	double t63 = t17 * q21 / 0.2e1;
	double t67 = t55 * ((t53 + t54) * t29 + (t58 + t59) * t38 + (t62 + t63) * t47);
	double t69 = t53 + t54 - t67 * t29;
	double t70 = t69 * t69;
	double t72 = t58 + t59 - t67 * t38;
	double t73 = t72 * t72;
	double t75 = t62 + t63 - t67 * t47;
	double t76 = t75 * t75;
	double t78 = sqrt(t70 + t73 + t76);
	double t79 = 0.1e1 / t78;
	double t81 = t51 * t38;
	double t82 = -t19 * q10 / 0.2e1;
	double t83 = t17 * q22 / 0.2e1;
	double t86 = -t19 * q11 / 0.2e1;
	double t87 = t17 * q23 / 0.2e1;
	double t90 = -t19 * q12 / 0.2e1;
	double t91 = t17 * q24 / 0.2e1;
	double t95 = t55 * ((t82 + t83) * t29 + (t86 + t87) * t38 + (t90 + t91) * t47);
	double t97 = t82 + t83 - t95 * t29;
	double t98 = t97 * t97;
	double t100 = t86 + t87 - t95 * t38;
	double t101 = t100 * t100;
	double t103 = t90 + t91 - t95 * t47;
	double t104 = t103 * t103;
	double t106 = sqrt(t98 + t101 + t104);
	double t107 = 0.1e1 / t106;
	double t108 = t107 * t103;
	double t110 = t51 * t47;
	double t111 = t107 * t100;
	double t113 = t79 * t69 - t81 * t108 + t110 * t111;
	double t116 = t107 * t97;
	double t118 = t79 * t75 - t52 * t111 + t81 * t116;
	double t123 = t79 * t72 - t110 * t116 + t52 * t108;
	A(0+1,1+0) = t52;
	A(0+1,1+1) = t113 / 0.2e1;
	A(0+1,1+2) = t81 * t118 / 0.2e1 - t110 * t123 / 0.2e1;
	A(1+1,1+0) = t81;
	A(1+1,1+1) = t123 / 0.2e1;
	A(1+1,1+2) = t110 * t113 / 0.2e1 - t52 * t118 / 0.2e1;
	A(2+1,1+0) = t110;
	A(2+1,1+1) = t118 / 0.2e1;
	A(2+1,1+2) = t52 * t123 / 0.2e1 - t81 * t113 / 0.2e1;

	//if (GetMBS()->GetTime() > 0.05)
	//	global_uo << "A=" << A << ", A1=" << A1 << "\n";

	return A;

}

