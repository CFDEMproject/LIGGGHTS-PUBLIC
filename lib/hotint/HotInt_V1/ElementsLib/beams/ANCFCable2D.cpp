//#**************************************************************
//#
//# filename:             ANCFCable2D.cpp
//#
//# author:               Gerstmayr Johannes & Rafael Ludwig
//#
//# generated:						24. April 2006
//# description:          Driver and model for timeintegration
//#                       2D ANCF - element
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
 
#include "rigid2d.h"
#include "femathhelperfunctions.h"
#include "material.h"
#include "ANCFCable2D.h"
#include "Node.h"
#include "graphicsconstants.h"
#include "elementdataaccess.h"
//#include "solversettings_auto.h"
#include "sensors.h"


int kappamode = 4-2; //1=konventional 3rd power in denominator of curvature, 2=quadratic=RIGHT, 
//3=only rxn in denominator, 4=GreenPiola
//kappa=|abs(cross(r,x ,r,xx))|/abs(rx)^3=:kappa_num/kappa_den
int epsmode = 2; //1==Green, 2=Biot=Dmitrochenko
int usekappafabs = 0; //use absolute value of curvature --> leads to problems for precurved case


	// -------------------------------------------
	// CONSTRUCTOR, COPY-ROUTINES...

	//nodal coordinates of first and second point, no global nodes given
void ANCFCable2D::SetANCFCable2D(const Vector& xc1, const Vector& xc2, double rhoi, double Emi,
												 const Vector3D& si, const Vector3D& coli, int kappanodalI)
{
	ElementDefaultConstructorInitialization();

	//plasticstrains_height = 0;
	//plasticstrains_width = 0;
	//n1=0; n2=0; 
	sos2=SOS();
	//concentratedmass1 = 0;
	//concentratedmass2 = 0;

	kappanodal = kappanodalI;
	size = si;

	x_init = xc1.Append(xc2);
	//xg = xc1.Append(xc2);

	//lx = size.X(); ly = size.Y(); lz = size.Z();

	mass = GetLx()*GetLy()*GetLz()*rhoi;
	GetMBS()->UO() << "mass=" << mass << "\n";
	//lx=1;ly=1;lz=1;

	//        rho, Youngs modulus, poisson ratio = 0, plane = 1, planestress = 0
	//rho = rhoi; //$ DR 2013-02-04 deleted rho from class element, do not use it here!
	Material mat(GetMBS(), rhoi, Emi, 0., 1, 0);
	AddMaterial(mat);

	col = coli;

	BuildDSMatrices();
	e0 = x_init;
	x_init = x_init.Append(Vector(SOS())); //Velocity initial conditions can also be transformed by Tinv!!!

	//UO() << "xinit=" << x_init << "\n";

	GetMaterial().BeamEIy() = -1;
	GetMaterial().BeamEA() = 0;
	GetMaterial().BeamRhoA() = 0;
	GetMaterial().BeamA() = GetLy()*GetLz();
};

	// Element shares nodes n1 and n2 with other elements; element sets initial conditions for nodes
	// a new material with rhoi and Emi is added
void ANCFCable2D::SetANCFCable2D(const Vector& xc1, const Vector& xc2, int n1i, int n2i, double rhoi, double Emi,
												 const Vector3D& si, const Vector3D& coli, int kappanodalI)
{
	ElementDefaultConstructorInitialization();

	//plasticstrains_height = 0;
	//plasticstrains_width = 0;
	//concentratedmass1 = 0;
	//concentratedmass2 = 0;

	kappanodal = kappanodalI;
	n1=n1i; n2=n2i; //sos2=0;
	size = si;


	x_init = xc1.Append(xc2);
	//xg = xc1.Append(xc2);

	//lx = size.X(); ly = size.Y(); lz = size.Z(); //lx is parameter to identify reference configuration

	//size.X() = ComputeCurvedLength(x_init);			 //size.X() is true curve length, but other curve length is better for sliding joint!
	//UO() << "linit=" << lx << ", ltrue=" << size.X() << "\n";

	mass = GetLx()*GetLy()*GetLz()*rhoi; //use lx, because this is the parameter of the real geometry ...
	//lx=1;ly=1;lz=1;

	//        rho, Youngs modulus, poisson ratio = 0, plane = 1, planestress = 0
	Material mat(GetMBS(), rhoi, Emi, 0., 1, 0);
	AddMaterial(mat);
	//rho = rhoi; //$ DR 2013-02-04 deleted rho from class element, do not use it here!
	col = coli;

	//UO() << "Cable: rho=" << rho << ", E=" << Em << ", size=" << size 	<< ", n1=" << n1 << ", n2=" << n2 << "\n";

	BuildDSMatrices();
	e0 = x_init;
	x_init = x_init.Append(Vector(SOS())); //Velocity initial conditions can also be transformed by Tinv!!!

	//UO() << "xinit=" << x_init << "\n";

	GetMaterial().BeamEIy() = -1;
	GetMaterial().BeamEA() = 0;
	GetMaterial().BeamRhoA() = 0;
	GetMaterial().BeamA() = GetLy()*GetLz();
};

	// Element shares nodes n1 and n2 with other elements; element sets initial conditions and initial velocities for nodes
	// a new material with rhoi and Emi is added
void ANCFCable2D::SetANCFCable2D(const Vector& xc1, const Vector& xc2, const Vector& vc1, const Vector& vc2,
												 int n1i, int n2i, double rhoi, double Emi,
												 const Vector3D& si, const Vector3D& coli, int kappanodalI)
{
	ElementDefaultConstructorInitialization();
	//plasticstrains_height = 0;
	//plasticstrains_width = 0;
	//concentratedmass1 = 0;
	//concentratedmass2 = 0;

	kappanodal = kappanodalI;
	n1=n1i; n2=n2i; //sos2=0;
	size = si;

	x_init = xc1.Append(xc2);
	//xg = xc1.Append(xc2);

	//lx = size.X(); ly = size.Y(); lz = size.Z();

	mass = GetLx()*GetLy()*GetLz()*rhoi;
	//lx=1;ly=1;lz=1;

	Vector v_init = vc1.Append(vc2);

	//        rho, Youngs modulus, poisson ratio = 0, plane = 1, planestress = 0
	Material mat(GetMBS(), rhoi, Emi, 0., 1, 0);
	AddMaterial(mat);
	//rho = rhoi; //$ DR 2013-02-04 deleted rho from class element, do not use it here!
	col = coli;

	BuildDSMatrices();
	e0 = x_init;
	x_init = x_init.Append(v_init); //Velocity initial conditions can also be transformed by Tinv!!!

	//UO() << "xinit=" << x_init << "\n";

	GetMaterial().BeamEIy() = -1;
	GetMaterial().BeamEA() = 0;
	GetMaterial().BeamRhoA() = 0;
	GetMaterial().BeamA() = GetLy()*GetLz();
};

	// Element shares nodes n1 and n2 with other elements; element sets initial conditions for nodes
	// material is given by number
void ANCFCable2D::SetANCFCable2D(const Vector& xc1, const Vector& xc2, int n1i, int n2i, int materialnumi,
		const Vector3D& si, const Vector3D& coli, int kappanodalI)
{
	ElementDefaultConstructorInitialization();

	//plasticstrains_height = 0;
	//plasticstrains_width = 0;
	//concentratedmass1 = 0;
	//concentratedmass2 = 0;

	kappanodal = kappanodalI;
	n1=n1i; n2=n2i; //sos2=0;
	size = si;

	x_init = xc1.Append(xc2);
	//xg = xc1.Append(xc2);

	//lx = size.X(); ly = size.Y(); lz = size.Z(); 

	//size.X() = ComputeCurvedLength(x_init);			 //size.X() is true curve length, but other curve length is better for sliding joint!

	materialnum = materialnumi;

	mass = GetLx()*GetLy()*GetLz()*GetMaterial().Density();

	col = coli;

	BuildDSMatrices();
	e0 = x_init;
	x_init = x_init.Append(Vector(SOS())); //Velocity initial conditions can also be transformed by Tinv!!!


	if (GetBeamEIy() == 0)  // if Beam parameters were not initialized in Material..
	{
		GetMaterial().BeamEIy() = -1;
		GetMaterial().BeamEA() = 0;
		GetMaterial().BeamRhoA() = 0;
		GetMaterial().BeamA() = GetLy()*GetLz();
	}

	if (GetMaterial().IsInelasticMaterial())
	{
		Vector datainit(DataS()); //automatically initialized with zeros
		SetDataInit(datainit); //initialize data variables with zero = initial inelastic strain == zero
	}

};







	// ------------------------------------------------
	// INITIALIZING

void ANCFCable2D::LinkToElements()
{
	if (SOSowned() == 0)
	{
		//UO() << "Link nodes to elements in Cable\n";
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


void ANCFCable2D::BuildDSMatrices()
{
	//orderx = 9; //max 9x5, sonst array grad zu klein!!!!
	//GetIntegrationRule(x1,w1,orderx);
};






	// --------------------------------------------------------
	// Computational routines:

//H = int(S,dV,V)
void ANCFCable2D::GetH(Matrix& H)
{
	if (Hmatrix.Getrows() == SOS())
	{
		H = Hmatrix;
		return;
	}
	else
	{
		// space dimension = 2
		int dim = Dim();
		// number of shape functions
		int ns = NS();
		// material parameter
		double rhoA = GetLy()*GetLz() * Rho();
		// area
		double A = GetLy()*GetLz();
		if (IsBeamParameters() && Rho() != 0) 
		{
			rhoA = GetRhoA();
			A = rhoA/Rho(); //this is the correct cross-section for body loads
		}

		H.SetSize(ns*dim,dim);
		H.SetAll(0);
	ConstVector<ANCFCable2D_MaxNShapes> SV(ns);
		ConstVector<ANCFCable2D_MaxIP> x1, w1;
		GetIntegrationRule(x1,w1,3+2*kappanodal); //3x1x1 !!!!!

		for (int i1=1; i1<=x1.GetLen(); i1++)
		{
			double p=x1(i1);
			GetS0(SV,p);

			// jacobi determinant
			Vector2D rx0 = GetPosx2D(p*GetLx()*0.5, x_init);
			double rxn = rx0.Norm();
			double det = 0.5*GetLx()*rxn;

			double fact = A * det * w1(i1);

			// integration
			for (int i=0; i<ns; i++)
			{
				for (int j=1; j<=dim; j++)
				{
					H(i*dim+j,j)+=fact*SV(i+1);
				}
			}
		}
		Hmatrix = H;
	}
}

void ANCFCable2D::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Body2D::GetElementData(edc);
}

int ANCFCable2D::SetElementData(ElementDataContainer& edc) //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
{
	int rv = Body2D::SetElementData(edc);
	// Element shares nodes n1 and n2 with other elements; element sets initial conditions for nodes
	// material is given by number

	//==> bodies with materials can only be added if at least one material exists
	if (GetMaterialNum() == 0 || GetMaterialNum() > GetMBS()->NMaterials()) {return 0;}

	mass = GetLx()*GetLy()*GetLz()*GetMaterial().Density();

	BuildDSMatrices();
	e0 = x_init.SubVector(1,SOS());
	//x_init = x_init.Append(Vector(SOS())); //Velocity initial conditions can also be transformed by Tinv!!!

	if (GetMaterial().IsInelasticMaterial())
	{
		Vector datainit(DataS()); //automatically initialized with zeros
		SetDataInit(datainit); //initialize data variables with zero = initial inelastic strain == zero
	}

	return rv;
}

//void ANCFCable2D::GetElementDataAuto(ElementDataContainer& edc) 		//fill in all element data
//{
//}
//
//int ANCFCable2D::SetElementDataAuto(const ElementDataContainer& edc) //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
//{
//	return 1;
//}


// Mass matrix
// M = int(rho*((S)^T).S, dV,V)
void ANCFCable2D::EvalM(Matrix& m, double t)
{
	if (massmatrix.Getcols() == SOS())
	{
		m = massmatrix;
		return;
	}
	else
	{
		m.SetAll(0);
		// space dimension, = 2
		int dim = Dim();
		// number of shape functions
		int ns = NS();
		// storage vector for shape functions
	ConstVector<ANCFCable2D_MaxNShapes> SV(ns);
		// shape function matrix
		// ( S1  0 S2  0 S3  0 ...
		// (  0 S1  0 S2  0 S3 ...
		Matrix HL(SOS(),dim);

		// integration rule
		ConstVector<ANCFCable2D_MaxIP> x1, w1;
		GetIntegrationRule(x1,w1,6+4*kappanodal);

		// material parameters
		double rhoA = GetLy()*GetLz() * Rho();
		if (IsBeamParameters()) rhoA = GetRhoA();

		for (int i1=1; i1<=x1.GetLen(); i1++)
		{
			double p = x1(i1);
			GetS0(SV,p);

			// fill in matrix HL
			for (int i=0; i<ns; i++)
			{
				for (int j=1; j<=dim; j++)
				{
					HL(i*dim+j,j)=SV(i+1);
				}
			}

			// jacobi determinant
			Vector2D rx0 = GetPosx2D(p*GetLx()*0.5, x_init);
			double rxn = rx0.Norm();
			double det = 0.5*GetLx()*rxn;

			// add to mass matrix
			m += det * rhoA * w1(i1) * (HL*HL.GetTp());
		}

		// add concentrated masses in nodal end points to mass matrix
		m(1,1) += concentratedmass1;
		m(2,2) += concentratedmass1;
		m(1+4+2*kappanodal,1+4+2*kappanodal) += concentratedmass2;
		m(2+4+2*kappanodal,2+4+2*kappanodal) += concentratedmass2;
		massmatrix = m;
	}  
		//UO() << "Mass_Cable2D =" << m << "\n";
};

// external forces
void ANCFCable2D::EvalF2(Vector& f, double t)
{
	Body2D::EvalF2(f,t);
	TMStartTimer(22);

	// residual vector
	ConstVector<ANCFCable2D_MaxDOF> fadd(SOS());
	fadd.SetAll(0);
	// temporary storage vector
	ConstVector<ANCFCable2D_MaxDOF> temp(SOS());
	// actual coordinates
	ConstVector<ANCFCable2D_MaxDOF> xg(SOS());

	// number of shape functions
	int ns = NS();
	// dimension, = 2
	int dim = Dim();
	// size of vector f
	int sos = SOS();

	// storage vector for shape functions
	ConstVector<ANCFCable2D_MaxNShapes> SV(ns);
	SV.SetLen(NS());

	// material parameters
	double EI_w = Em() * GetLz()*Cub(GetLy())/12.;
	double EI4 = Em() * GetLz() * To5(GetLy())/80.; //for Green Piola approach
	double EA = Em() * GetLy() * GetLz(); 
	double A = GetLy()*GetLz();//for Hellinger Reissner
	if (IsBeamParameters())
	{
		EI_w = GetBeamEIy();
		EA = GetBeamEA();
	}

	// integration order
	int order_kappa = 5+4*kappanodal; //regular: 5
	int order_eps   = 9+4*kappanodal; //regular: 9
	ConstVector<ANCFCable2D_MaxIP> xk, wk, xe, we;

	// integration rules
	GetIntegrationRule(xk,wk,order_kappa); 
	GetIntegrationRule(xe,we,order_eps);   


	// initial curvature:
	// compute kappa0 and eps0 as curvature/strain in initial state
	ConstVector<ANCFCable2D_MaxIP> kappa0, eps0;
	kappa0.SetLen(xk.Length());
	kappa0.SetAll(0);
	eps0.SetLen(xe.Length());
	eps0.SetAll(0);

	int usedmitrochenko = 0;  //for comparison with Euler elastica
	if (epsmode == 2) usedmitrochenko = 1;

	for (int i1=1; i1<=xk.GetLen(); i1++)
	{
		double x = xk(i1);
		kappa0(i1) = GetKappa(x,x_init);
	}

	for (int i1=1; i1<=xe.GetLen(); i1++)
	{
		double x = xe(i1);
		eps0(i1) = GetEpsAxial(x,x_init);
	}

	//SetComputeCoordinates();
	GetCoordinates(xg);

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//kappa*delta kappa:
	int hellingerreissnerkappa = 0;
	int hellingerreissnereps = 0;

	for (int i1=1; i1<=xk.GetLen(); i1++)
	{
		double x = xk(i1);

		// compute jacobi determinant
		Vector2D rx0 = GetPosx2D(x*GetLx()*0.5, x_init);
		double rxn = rx0.Norm();
		double det = 0.5*GetLx()*rxn; 

		// bending moments:
		double factkappa = EI_w*wk(i1)*det; 
		double kappa_plast = GetPlasticKappa(x*0.5*GetLx());

		// set temp to delta_kappa
		double kappa;
		GetDeltaKappa(x,xg,temp, kappa);

		if (!UseTangentStiffness() || !GetMaterial().IsInelasticMaterial() )
		{
		////// add EI kappa * delta_kappa to vector f
		temp *= (kappa - kappa0(i1) - kappa_plast)*factkappa;
		}
		else // TangentStiffness  // compute the integral int_{-H/2}^{H/2} y*stress(strain) * W dy
		{
			//////////// integrate int_(-h/2)^(h/2) y*sigma dy
			if (fabs(kappa0(i1))>1e-10) mbs->UO() << "warning ANCFCable2D: Tangent stiffness matrix not correct for precurved elements!\n";
			Matrix plasticstrain_mat, internalvar_mat;
			GetPlasticStrainMatrix(plasticstrain_mat);
			GetInternalVariableMatrix(internalvar_mat);
			// Vector plasticstrains contains plastic strains in gridpoints over height interpolated at integration point
			ConstVector<ANCFCable2D_maxppy> plasticstrains(plasticstrains_height);
			ConstVector<ANCFCable2D_maxppy> internalvars(plasticstrains_height);
			plasticstrain_mat.InterpolateX(x, plasticstrains);
			internalvar_mat.InterpolateX(x, internalvars);

			double epsax = GetEpsAxial(x,xg);
			double bendingmom = 0;
			double tangmod = GetMaterial().TangentModule();

			double delta_y = GetLy()/(plasticstrains_height-1);
			double factor = 0.5*GetLz()*delta_y;
			// integration loop - trapezoidal rule
			for (int ii=1; ii<=plasticstrains_height; ii++)
			{
				Vector2D p(x*GetLx()*0.5, GetLy()*0.5-(ii-1)*delta_y);
				double strain = epsax - p.Y()*kappa - plasticstrains(ii);
				double stress = Em()*(strain);
				double yield = GetMaterial().YieldStress()+tangmod*internalvars(ii);
				if (stress > yield) 
				{
					stress = yield + tangmod*(strain-yield/Em());
				}
				if (stress < -yield)
				{
					stress = -yield + tangmod*(strain+yield/Em());
				}
				bendingmom -= factor*stress*p.Y();
				if (ii>1 && ii<plasticstrains_height)
				{
					bendingmom -= factor*stress*p.Y();
				}
			}

			temp *= bendingmom*wk(i1)*det;
		}
		fadd += temp;
	}

	//eps*delta eps:
	for (int i1=1; i1<=xe.GetLen(); i1++)
	{
		double x = xe(i1);

		// compute jacobi determinant
		Vector2D rx0 = GetPosx2D(x*GetLx()*0.5, x_init);
		double rxn = rx0.Norm();
		double det = 0.5*GetLx()*rxn; 

		// axial strains:
		double facteps = EA*we(i1)*det;

		// plastic axial strain
		double eps_plast = GetPlasticAxStrain(x*0.5*GetLx());

		int optimized = 1;
		if (UseTangentStiffness()) optimized = 0;

		if (!optimized)
		{
			// set temp to delta_eps
			GetDeltaEpsAxial(x,xg,temp);

			if (UseTangentStiffness()) // compute the integral int_{-H/2}^{H/2} stress(strain) * W dy
			{
				if (fabs(eps0(i1))>1e-10) mbs->UO() << "warning ANCFCable2D: Tangent stiffness matrix not correct for precurved elements!\n";
				Matrix plasticstrain_mat, internalvar_mat;
				GetPlasticStrainMatrix(plasticstrain_mat);
				GetInternalVariableMatrix(internalvar_mat);
				// Vector plasticstrains contains plastic strains in gridpoints over height interpolated at x-coord x0
				ConstVector<ANCFCable2D_maxppy> plasticstrains(plasticstrains_height);
				ConstVector<ANCFCable2D_maxppy> internalvars(plasticstrains_height);
				plasticstrain_mat.InterpolateX(x, plasticstrains);
				internalvar_mat.InterpolateX(x, internalvars);

				double kappa = GetKappa(x, xg);
				double epsax = GetEpsAxial(x,xg);
				double axforce = 0;
				double tangmod = GetMaterial().TangentModule();
				double delta_y = GetLy()/(plasticstrains_height-1);
				double factor = 0.5*GetLz()*delta_y;
				// integration loop - trapezoidal rule
				for (int ii=1; ii<=plasticstrains_height; ii++)
				{
					Vector2D p(x*GetLx()*0.5, GetLy()*0.5-(ii-1)*delta_y);
					double strain = epsax - p.Y()*kappa-plasticstrains(ii);
					double stress = Em()*(strain);
					double yield = GetMaterial().YieldStress() + tangmod*internalvars(ii);
					if (stress > yield) stress = yield + tangmod*(strain-yield/Em());
					if (stress < -yield) stress = -yield + tangmod*(strain+yield/Em());
					axforce += factor*stress;
					if (ii>1 && ii<plasticstrains_height)
					{
						axforce += factor*stress;
					}
				}


				temp *= axforce*we(i1)*det;
			}
			else // No TangentStiffness
			{
				// add EA (eps-epsplast-epsinit) * delta_eps to vector f
				temp *= (GetEpsAxial(x,xg)-eps_plast)*facteps;
			}
			fadd += temp;

		}
		else
		{
			GetS0x(SV,x);
			//compute rx:
			Vector2D rx(0.,0.);
			for (int i = 1; i <= Dim(); i++)
			{
				for (int j = 1; j <= NS(); j++)
				{
					rx(i) += SV(j)*xg((j-1)*Dim()+i);
				}
			}

			//compute deltaeps:
			int dim = Dim();
			int ns = NS();

			double fact = 1;
			double eps;
			double kappa; //for formulation based on Green strain!

			//1==Green, 2=Biot=Dmitrochenko, 3=linearized strain
			if (epsmode == 2)
			{
				eps = (sqrt(rx*rx)-1.0) - eps0(i1) - eps_plast;
				fact = 1./(sqrt(rx*rx)); //additional part for variation!
			}
			else if (epsmode == 1)
			{
				if (kappamode == 1 || kappamode == 2)
				{ 
					eps = 0.5*(rx*rx-1.0) - eps0(i1) - eps_plast; //Green strain
				}
				else if (kappamode == 3)
				{
					kappamode = 2; //get K, not K^GP, this is a little hack !!!
					GetDeltaKappa(x,xg,temp, kappa);
					kappamode = 3; //switch back to old kappa mode

					eps = 0.5*(rx*rx-1.0) + 0.5*(kappa*kappa)*EI_w/EA - eps0(i1) - eps_plast;
				}
				else if (kappamode == 4)
				{
					kappamode = 2; //get K, not K^GP, this is a little hack !!!
					GetDeltaKappa(x,xg,temp, kappa);
					kappamode = 4; //switch back to old kappa mode

					eps = 0.5*(rx*rx-1.0)*rx.Norm() + 1.5*EI_w/EA*Sqr(kappa)*rx.Norm() - eps0(i1) - eps_plast;
				}

			}

			//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			if (0) //this is the geometrical correction of the strain for the original curvature with 3rd power in denominator!
			{
				double kappa=GetKappa(x,xg);
				eps += (kappa*kappa)*EI_w/EA; //correction term for correct normal forces and stretch
			}
			if (0) //this is the geometrical correction of the strain for the curvature times |rx|^2! (Studienblätter)
			{
				double kappa=GetKappa(x,xg);
				eps -= (kappa*kappa)*EI_w/EA; //correction term for correct normal forces and stretch
			}
			//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			rx *= eps * facteps * fact;

			for (int i=1; i <= dim; i++)
			{
				for (int j=1; j <= ns; j++)
				{
					fadd((j-1)*dim+i) += SV(j)*rx(i);
				}
			}
			if (kappamode == 3 && epsmode == 1)
			{
				temp *= 1.*(EI_w/EA*facteps*kappa*(eps - 0.5*(kappa*kappa)*EI_w/EA + 0.5*EI4/EI_w*kappa*kappa ));
				fadd += temp;
			}

			//compute epsaxial*deltaeps*fact:
			//temp *= 0.5*(rx*rx-1.0)*facteps;
		}
	}


	// add local residual vector to global vector f
	f -= fadd;

	TMStopTimer(22);

		// +++++++ damping: +++++++
	if (GetMassDamping() != 0)
	{
		for (int i = 1; i <= SOS(); i++)
			xg(i) = XGP(i);

		if (massmatrix.Getcols() == SOS())
		{
			Mult(massmatrix,xg,temp);
		}
		else
		{
			Matrix dmat;
			dmat.SetSize(SOS(),SOS());
			EvalM(dmat,t);
			massmatrix = dmat;
			Mult(dmat,xg,temp);
		}
		temp *= GetMassDamping();
		f -= temp;
	}
};



	// -------------------------------------------------------
	// Nonlinear iteration for plasticity

// changes plastic strains, returns yieldfunction(sigma)/Em
double ANCFCable2D::PostNewtonStep(double t)		
{

	if (!GetMaterial().IsInelasticMaterial()) return 0;

	Material& mat = GetMaterial();

	TMStartTimer(24);
	// size of a plastic cell
	double hx = 2./(plasticstrains_width-1);
	double hy = 2./(plasticstrains_height-1);

	// get current plastic information
	Matrix plasticstrains, internalvar;
	GetPlasticStrainMatrix(plasticstrains);
	if (mat.IsHardeningMaterial())
		GetInternalVariableMatrix(internalvar);

	ConstVector<ANCFCable2D_MaxDOF> xgen(SOS());
	GetCoordinates(xgen);

	// get plastic information from last approved iterative step
	const Matrix plasticstrains_old, internalvar_old;
#ifdef AND_INCLUDES
	// AND: Transport of plastic strains is included
	GetPlasticStrainMatrixLastStepTransported(plasticstrains_old);
	if (mat.IsHardeningMaterial())
		GetInternalVariableMatrixLastStepTransported(internalvar_old);
#else // no AND_INCLUDES
	// no AND: no transport available, get plastic strains from last approved iterative step
	GetPlasticStrainMatrixLastStep(plasticstrains_old);
	if (mat.IsHardeningMaterial())
		GetInternalVariableMatrixLastStep(internalvar_old);
#endif

	// loop over all plastic integration points
	double nlerror = 0;

	TMStopTimer(24);
	TMStartTimer(25);
	double loadfact_init = GetInitLoadfact();
	for (int j=1; j<=plasticstrains_width; j++)
	{
		double x0 = (j-1.)*hx-1.;	
		double eps_axial = GetEpsAxial(x0, xgen);
		double kappa = GetKappa(x0,xgen);
		double eps_axial_init = GetEpsAxial(x0, x_init);
		double kappa_init = GetKappa(x0,x_init);
		for (int i=1; i<=plasticstrains_height; i++)
		{
			// local position
			double y0 = 1.-(i-1.)*hy;
			double yb = 0.5*y0*GetLy();
			
			// total strain
			// slow version:
			//Vector2D ploc(x0*0.5*GetLx(), yb);
			//double strain = GetEps(ploc) - GetEpsInit(ploc);
			// faster:
			double strain = (eps_axial-yb*kappa) - loadfact_init*(eps_axial_init-yb*kappa_init);

			double stress;
			double plasticstrain_old = plasticstrains(i,j);
			double error = 0;
			// call material routine

			if (mat.IsHardeningMaterial())
				error += mat.DoNonlinStepBeam2D(stress, strain, plasticstrains(i,j), internalvar(i,j), plasticstrains_old(i,j), internalvar_old(i,j), UseTangentStiffness());
			else
			{
				double dummy = 0;
				error += mat.DoNonlinStepBeam2D(stress, strain, plasticstrains(i,j), dummy, plasticstrains_old(i,j), dummy, UseTangentStiffness());
			}
			nlerror += error;
		}
	}

	nlerror *= GetLx()/(plasticstrains_height*plasticstrains_width);
#ifdef COMPILE_AND
	// HACK AP:
	// this factor 2 is to regain the same error level as in an earlier implementation, 
	// such that AND-input files return same results for same accuracy
	nlerror *= 0.5;
#endif
	TMStopTimer(25);
	if (UseTangentStiffness()) return 0;
	return nlerror;
}

void ANCFCable2D::PostprocessingStep()
{
	if (!GetMaterial().IsInelasticMaterial()) return;
}


	// ------------------------------------------------------
	// VISUALIZATION

const char * str_plastic_curvature = "plastic curvature";
const char * str_plastic_ax_strain = "plastic axial strain";
const char * str_internal_variable = "internal variable";

void ANCFCable2D::GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables)
{
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_displacement,
		FieldVariableDescriptor::FVCI_y, false);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_stress,
		FieldVariableDescriptor::FVCI_y, true, true);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_stress_mises);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_total_strain);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_axial_extension);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_curvature);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_force_axial);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_moment_bending);
	if(IsMaterialAvailable() && GetMaterial().IsInelasticMaterial())
	{
		FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_inelastic_strain);
		// problem specific variables
		variables.Add(FieldVariableDescriptor(str_plastic_curvature));
		variables.Add(FieldVariableDescriptor(str_plastic_ax_strain));
		if(GetMaterial().IsHardeningMaterial())
			variables.Add(FieldVariableDescriptor(str_internal_variable));
	}
}

double ANCFCable2D::GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector2D & local_position, bool flagD)
{
	//ploc is from -1 .. +1

	//only for rectangular cross-section!:
	ConstVector<ANCFCable2D_MaxDOF> xgd(SOS());

	if (flagD)
		GetDrawCoordinates(xgd);
	else
		GetCoordinates(xgd);

	// displacements
	if(fvd.VariableType() == FieldVariableDescriptor::FVT_displacement)
		return fvd.GetComponent(GetDisplacement2D(local_position.X()*GetLx()*0.5,flagD));

	// x0 in [-1,1]
	double x0 = local_position.X();
	// AP: yb has to be in [-Ly/2, Ly/2] ! see definition of ploc_scal below!
	double yb = local_position.Y()*0.5*GetLy();

	////  COMPUTE INITIAL CURVATURE
	double kappa_init = 0, eps_init = 0;
	double kappa, eps0;

	if (flagD)
	{
		kappa_init = GetKappaD(x0,x_init);
		eps_init = GetEpsAxialD(x0,x_init);
		kappa = GetKappaD(x0,xgd) - kappa_init;
		eps0 = GetEpsAxialD(x0, xgd) - eps_init; 
	}
	else
	{
		kappa_init = GetKappa(x0,x_init);
		eps_init = GetEpsAxial(x0,x_init);
		kappa = GetKappa(x0,xgd) - kappa_init;
		eps0 = GetEpsAxial(x0, xgd) - eps_init; 
	}


	Vector2D ploc_scal(x0*0.5*GetLx(), yb);
	double kappa_plast = GetPlasticKappa(x0*0.5*GetLx(),flagD);
	double eps_plast = GetPlasticAxStrain(x0*0.5*GetLx(),flagD);
	double eps_plast_loc = GetPlasticStrain(ploc_scal,flagD);

	Matrix plasticstrains, internalvars;
	if (IsMaterialAvailable() && GetMaterial().IsInelasticMaterial())
		GetPlasticStrainMatrix(plasticstrains,1);
	if (IsMaterialAvailable() && GetMaterial().IsInelasticMaterial()&&GetMaterial().IsHardeningMaterial())
		GetInternalVariableMatrix(internalvars,	1);

	double eps = eps0 - kappa*yb;



	Matrix3D stress(0);
	Matrix3D strain(0);
	stress(1,1) = (eps-eps_plast_loc)*Em();
	strain(1,1) = eps;

	double A = GetLy()*GetLz();
	double EI_w = Em() * GetLz()*Cub(GetLy())/12.;
	double EA = Em() * A; 
	if (IsBeamParameters())
	{
		EI_w = GetBeamEIy();
		EA = GetBeamEA();
	}

	switch(fvd.VariableType())
	{
	case FieldVariableDescriptor::FVT_stress: return fvd.GetComponent(stress);
	case FieldVariableDescriptor::FVT_stress_mises: return stress.Mises();
	case FieldVariableDescriptor::FVT_total_strain: return strain(1,1);
	case FieldVariableDescriptor::FVT_beam_axial_extension: return eps0;
	case FieldVariableDescriptor::FVT_beam_curvature: return kappa;
	case FieldVariableDescriptor::FVT_beam_force_axial: return EA * (eps0-eps_plast);
	case FieldVariableDescriptor::FVT_beam_moment_bending: return -EI_w * (kappa-kappa_plast);		//not true for prescribed case!
	}

	if(IsMaterialAvailable() && GetMaterial().IsInelasticMaterial())
	{
		switch(fvd.VariableType())
		{
		case FieldVariableDescriptor::FVT_inelastic_strain: return eps_plast_loc;		// inelastic strain  eps_xx^p
		//case FieldVariableDescriptor::FVT_plastic_strain: return plasticstrains.Interpolate(local_position);	// plastic strain matrix
		case FieldVariableDescriptor::FVT_problem_specific:
			{
				if(fvd.GetTextualIdentifierWithoutComponents() == str_plastic_curvature)		// plastic curvature kappa^p
					return kappa_plast;
				if(fvd.GetTextualIdentifierWithoutComponents() == str_plastic_ax_strain)		// plastic axial strain eps_ax^p
					return eps_plast;
				if(fvd.GetTextualIdentifierWithoutComponents() == str_internal_variable && GetMaterial().IsHardeningMaterial())		// internal variable matrix
					return internalvars.Interpolate(local_position);
			}
		}
	}

	return FIELD_VARIABLE_NO_VALUE;
}



void ANCFCable2D::DrawElement()
{
	mbs->SetColor(col);

	//// draw elements in front of origin, at z = 0.1
	//Vector3D offset(0,0,0.1);
	//Vector3D offset(0, 0, 0.01);
	Vector3D offset(0.,0.,0.);
	//$!AD: [$%&@"§$& ARGGGGGHHHH] 

	// length and magnified thickness of element
	double lx1 = GetLx(); double ly1 = GetLy()*GetMBS()->GetMagnifyYZ();

	// get draw mode from Visualization window
	//0=no lines, 1=outline+color, 2=outline, 3=elementline+color, 4=elementline
	int linemode = 1; 
	if (GetMBS()->GetIOption(110) && !GetMBS()->GetIOption(111))
	{
		linemode = 2;
	}
	if (!GetMBS()->GetIOption(110) && GetMBS()->GetIOption(111))
	{
		linemode = 0;
	}

	int colormode = 0;
	if (GetMBS()->GetActualPostProcessingFieldVariable() != NULL)
		colormode = 1;

	double def_scale = GetMBS()->GetDOption(105); //deformation scaling

	int drawlinemode = GetMBS()->GetIOption(117); // Draw Flat Elements, in this case, function is plotted normal to element axis

	if ( (1||colormode) && !drawlinemode) // plot coloured elements
	{
		lx1 /= GetLx();
		ly1 /= GetLy();

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
		double scalex = lx1/GetLx();
		double scaley = ly1/GetLy();

		Vector3D col_line(0.8,0,0);

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

			// normal direction (numerical differetiation for the moment)
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
				mbs->MyDrawLine(drawpoints(cnt-1), drawpoints(cnt), ly1*1e2, col_line);
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
		//mbs->DrawPolygon(nodepoints, 1, GetLy()*ly1*1e-1);
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
				p1 = GetPos2D0D(Vector2D(l1/(0.5*GetLx()), 1.), def_scale);
				p2 = GetPos2D0D(Vector2D(l1/(0.5*GetLx()),-1.), def_scale);

				if (concentratedmass1 != 0)
				{
					mbs->SetColor(colred);
					double r = pow(3.*concentratedmass1/(4.*MY_PI*Rho()),1./3.);
					//double r = GetLy()*GetMBS()->GetMagnifyYZ()*1;
					mbs->DrawSphere(Vector3D(p1.X(),p1.Y(),0.),r,12); // ToP3d!!!!!!!!
					//mbs->DrawSphere(p1,ly1*2,8);
					mbs->SetColor(col);
				}
			}
			else
			{
				p1 = p4;
				p2 = p3;
			}
			p3 = GetPos2DD(Vector2D(l2,-GetLy()*0.5));
			p4 = GetPos2DD(Vector2D(l2, GetLy()*0.5));
			mbs->DrawQuad(ToP3D(p4),ToP3D(p3),ToP3D(p2),ToP3D(p1));
		}

		if (concentratedmass2 != 0)
		{
			mbs->SetColor(colred);
			double r = pow(3.*concentratedmass2/(4.*MY_PI*Rho()),1./3.);
			mbs->DrawSphere(Vector3D(p2.X(),p2.Y(),0),r,12);
		}
	}
};


	// Shape function routines ------------------------------------
void ANCFCable2D::GetS0(Vector& sf, const double& ploc) const
{
	double xb = ploc;
	sf.SetLen(NS()); //4

	if (!kappanodal)
	{
		double xb2 = xb*xb;
		double xb3 = xb2*xb;
		sf(1) = 1.0/2.0-3.0/4.0*xb+xb3/4.0;
		sf(2) = (1.0-xb-xb2+xb3)*GetLx()/8.0;
		sf(3) = 1.0/2.0+3.0/4.0*xb-xb3/4.0;
		sf(4) = Sqr(1.0+xb)*(-1.0+xb)*GetLx()/8.0;
	}
	else
	{
		sf(1) = 1.0/2.0-15.0/16.0*xb+5.0/8.0*xb*xb*xb-3.0/16.0*xb*xb*xb*xb*xb;
		sf(2) = -(-5.0+7.0*xb+6.0*xb*xb-10.0*xb*xb*xb-xb*xb*xb*xb+3.0*xb*xb*xb*xb*xb)*GetLx()/32.0;
		sf(3) = -pow(1.0+xb,2.0)*(-1.0+3.0*xb-3.0*xb*xb+xb*xb*xb)*GetLx()*GetLx()/64.0;
		sf(4) = 1.0/2.0+15.0/16.0*xb-5.0/8.0*xb*xb*xb+3.0/16.0*xb*xb*xb*xb*xb;
		sf(5) = -pow(1.0+xb,3.0)*(5.0-8.0*xb+3.0*xb*xb)*GetLx()/32.0;
		sf(6) = pow(1.0+xb,3.0)*(-2.0*xb+1.0+xb*xb)*GetLx()*GetLx()/64.0;

		/*
		sf(1) = 1.0/2.0-15.0/16.0*xb+5.0/8.0*xb*xb*xb-3.0/16.0*xb*xb*xb*xb*xb;
		sf(2) = -(-5.0+7.0*xb+6.0*xb*xb-10.0*xb*xb*xb-xb*xb*xb*xb+3.0*xb*xb*xb*xb*xb)*GetLx()/32.0;
		sf(3) = -Sqr(1.0+xb)*(-1.0+3.0*xb-3.0*xb*xb+xb*xb*xb)*GetLx()*GetLx()/64.0;
		sf(4) = 1.0/2.0+15.0/16.0*xb-5.0/8.0*xb*xb*xb+3.0/16.0*xb*xb*xb*xb*xb;
		sf(5) = -Cub(1.0+xb)*(5.0-8.0*xb+3.0*xb*xb)*GetLx()/32.0;
		sf(6) = Cub(1.0+xb)*(-2.0*xb+1.0+xb*xb)*GetLx()*GetLx()/64.0;
		*/
	}
}
// Get Shape Function i
double ANCFCable2D::GetS0(int i, const double& ploc) const
{
	double xb = ploc;

	if (!kappanodal)
	{
		double xb2 = xb*xb;
		double xb3 = xb2*xb;
		switch (i)
		{
		case 1:
			return 1.0/2.0-3.0/4.0*xb+xb3/4.0;
		case 2:
			return (1.0-xb-xb2+xb3)*GetLx()/8.0;
		case 3:
			return 1.0/2.0+3.0/4.0*xb-xb3/4.0;
		case 4:
			return Sqr(1.0+xb)*(-1.0+xb)*GetLx()/8.0;
		default:
			mbs->UO() << "ANCFCable2D::GetS0 Error: Shape function " << i << " not available!\n";
			return 0;
		}
	}
	else
	{
		switch (i)
		{
		case 1:
			return 1.0/2.0-15.0/16.0*xb+5.0/8.0*xb*xb*xb-3.0/16.0*xb*xb*xb*xb*xb;
		case 2:
			return -(-5.0+7.0*xb+6.0*xb*xb-10.0*xb*xb*xb-xb*xb*xb*xb+3.0*xb*xb*xb*xb*xb)*GetLx()/32.0;
		case 3:
			return -pow(1.0+xb,2.0)*(-1.0+3.0*xb-3.0*xb*xb+xb*xb*xb)*GetLx()*GetLx()/64.0;
		case 4:
			return 1.0/2.0+15.0/16.0*xb-5.0/8.0*xb*xb*xb+3.0/16.0*xb*xb*xb*xb*xb;
		case 5:
			return -pow(1.0+xb,3.0)*(5.0-8.0*xb+3.0*xb*xb)*GetLx()/32.0;
		case 6:
			return pow(1.0+xb,3.0)*(-2.0*xb+1.0+xb*xb)*GetLx()*GetLx()/64.0;
		default:
			mbs->UO() << "ANCFCable2D::GetS0 Error: Shape function " << i << " not available!\n";
			return 0;
		}
	}
	return 0;
}

// Get Shape Function Vector dS/dxi
void ANCFCable2D::GetS0x(Vector& sfx, const double& ploc) const
{
	double xb = ploc;
	sfx.SetLen(NS());
	double f = 2./GetLx();

	if (!kappanodal)
	{
		sfx(1) = f*(-3.0/4.0+3.0/4.0*xb*xb);
		sfx(2) = f*((-1.0-2.0*xb+3.0*xb*xb)*GetLx()/8.0);
		sfx(3) = f*(3.0/4.0-3.0/4.0*xb*xb);
		sfx(4) = f*((1.0+xb)*(-1.0+xb)*GetLx()/4.0+Sqr(1.0+xb)*GetLx()/8.0);
	}
	else
	{
		//f included:
		sfx(1) = -0.125*(15.0-30.0*xb*xb+15.0*xb*xb*xb*xb)/GetLx();
		sfx(2) = -0.4375-0.75*xb+0.1875E1*xb*xb+0.25*xb*xb*xb-0.9375*xb*xb*xb*xb;
		sfx(3) = -0.625E-1*GetLx()*(1.0+xb)*(-1.0+3.0*xb-3.0*xb*xb+xb*xb*xb)-0.3125E-1*GetLx()*pow(1.0+xb,2.0)*(3.0-6.0*xb+3.0*xb*xb);
		sfx(4) = 0.375*pow(1.0+xb,2.0)*(8.0-9.0*xb+3.0*xb*xb)/GetLx()+0.125*pow(1.0+xb,3.0)*(-9.0+6.0*xb)/GetLx();
		sfx(5) = -0.4375+0.75*xb+0.1875E1*xb*xb-0.25*xb*xb*xb-0.9375*xb*xb*xb*xb;
		sfx(6) = 0.9375E-1*GetLx()*pow(1.0+xb,2.0)*(-2.0*xb+1.0+xb*xb)+0.3125E-1*GetLx()*pow(1.0+xb,3.0)*(-2.0+2.0*xb);

		/*
		sfx(1) = f*(-15.0/16.0+15.0/8.0*xb*xb-15.0/16.0*xb*xb*xb*xb);
		sfx(2) = f*(-(7.0+12.0*xb-30.0*xb*xb-4.0*xb*xb*xb+15.0*xb*xb*xb*xb)*GetLx()/32.0);
		sfx(3) = f*(-(1.0+xb)*(-1.0+3.0*xb-3.0*xb*xb+xb*xb*xb)*GetLx()*GetLx()/32.0-pow(1.0+xb,2.0)*(3.0-6.0*xb+3.0*xb*xb)*GetLx()*GetLx()/64.0);
		sfx(4) = f*(15.0/16.0-15.0/8.0*xb*xb+15.0/16.0*xb*xb*xb*xb);
		sfx(5) = f*(-3.0/32.0*pow(1.0+xb,2.0)*(5.0-8.0*xb+3.0*xb*xb)*GetLx()-pow(1.0+xb,3.0)*(-8.0+6.0*xb)*GetLx()/32.0);
		sfx(6) = f*(3.0/64.0*pow(1.0+xb,2.0)*(-2.0*xb+1.0+xb*xb)*GetLx()*GetLx()+pow(1.0+xb,3.0)*(-2.0+2.0*xb)*GetLx()*GetLx()/64.0);
		*/
	}
}

// Get Shape Function  dS(i)/dxi
double ANCFCable2D::GetS0x(int i, const double& ploc) const
{
	double xb = ploc;
	double f = 2./GetLx();

	if (!kappanodal)
	{
		switch (i)
		{
		case 1:
			return f*(-3.0/4.0+3.0/4.0*xb*xb);
		case 2:
			return f*((-1.0-2.0*xb+3.0*xb*xb)*GetLx()/8.0);
		case 3:
			return f*(3.0/4.0-3.0/4.0*xb*xb);
		case 4:
			return f*((1.0+xb)*(-1.0+xb)*GetLx()/4.0+Sqr(1.0+xb)*GetLx()/8.0);
		default:
			mbs->UO() << "ANCFCable2D::GetS0x Error: Shape function " << i << " not available!\n";
			return 0;
		}
	}

	else
	{
		//f included:
		switch (i)
		{
		case 1:
			return -0.125*(15.0-30.0*xb*xb+15.0*xb*xb*xb*xb)/GetLx();
		case 2:
			return -0.4375-0.75*xb+0.1875E1*xb*xb+0.25*xb*xb*xb-0.9375*xb*xb*xb*xb;
		case 3:
			return -0.625E-1*GetLx()*(1.0+xb)*(-1.0+3.0*xb-3.0*xb*xb+xb*xb*xb)-0.3125E-1*GetLx()*pow(1.0+xb,2.0)*(3.0-6.0*xb+3.0*xb*xb);
		case 4:
			return 0.375*pow(1.0+xb,2.0)*(8.0-9.0*xb+3.0*xb*xb)/GetLx()+0.125*pow(1.0+xb,3.0)*(-9.0+6.0*xb)/GetLx();
		case 5:
			return -0.4375+0.75*xb+0.1875E1*xb*xb-0.25*xb*xb*xb-0.9375*xb*xb*xb*xb;
		case 6:
			return 0.9375E-1*GetLx()*pow(1.0+xb,2.0)*(-2.0*xb+1.0+xb*xb)+0.3125E-1*GetLx()*pow(1.0+xb,3.0)*(-2.0+2.0*xb);
		default:
			mbs->UO() << "ANCFCable2D::GetS0x Error: Shape function " << i << " not available!\n";
			return 0;
		}
	}
	return 0;
}


//d^2S/dxi^2
void ANCFCable2D::GetS0xx(Vector& sfxx, const double& ploc) const
{
	double xb = ploc;
	sfxx.SetLen(NS());
	double f = 4./Sqr(GetLx());

	if (!kappanodal)
	{
		sfxx(1) = f*(3.0/2.0*xb);
		sfxx(2) = f*((-2.0+6.0*xb)*GetLx()/8.0);
		sfxx(3) = f*(-3.0/2.0*xb);
		sfxx(4) = f*((-1.0+xb)*GetLx()/4.0+(1.0+xb)*GetLx()/2.0);
	}
	else
	{
		//f included:
		sfxx(1) = -0.25*(-60.0*xb+60.0*xb*xb*xb)/(GetLx()*GetLx());
		sfxx(2) = -0.125*(12.0-60.0*xb-12.0*xb*xb+60.0*xb*xb*xb)/GetLx();
		sfxx(3) = -0.25+0.75*xb+0.75*xb*xb-0.125E1*xb*xb*xb;
		sfxx(4) = 0.15E1*(1.0+xb)*(8.0-9.0*xb+3.0*xb*xb)/(GetLx()*GetLx())+0.15E1*pow(1.0+xb,2.0)*(-9.0+6.0*xb)/(GetLx()*GetLx())+0.15E1*pow(1.0+xb,3.0)/(GetLx()*GetLx());
		sfxx(5) = -0.75*(1.0+xb)*(5.0-8.0*xb+3.0*xb*xb)/GetLx()-0.75*pow(1.0+xb,2.0)*(-8.0+6.0*xb)/GetLx()-0.75*pow(1.0+xb,3.0)/GetLx();
		sfxx(6) = -0.25-0.75*xb+0.75*xb*xb+0.125E1*xb*xb*xb;

		/*
		sfxx(1) = f*(15.0/4.0*xb-15.0/4.0*xb*xb*xb);
		sfxx(2) = f*(-(12.0-60.0*xb-12.0*xb*xb+60.0*xb*xb*xb)*GetLx()/32.0);
		sfxx(3) = f*(-(-1.0+3.0*xb-3.0*xb*xb+xb*xb*xb)*GetLx()*GetLx()/32.0-(1.0+xb)*(3.0-6.0*xb+3.0*xb*xb)*GetLx()*GetLx()/16.0-pow(1.0+xb,2.0)*(-6.0+6.0*xb)*GetLx()*GetLx()/64.0);
		sfxx(4) = f*(-15.0/4.0*xb+15.0/4.0*xb*xb*xb);
		sfxx(5) = f*(-3.0/16.0*(1.0+xb)*(5.0-8.0*xb+3.0*xb*xb)*GetLx()-3.0/16.0*pow(1.0+xb,2.0)*(-8.0+6.0*xb)*GetLx()-3.0/16.0*pow(1.0+xb,3.0)*GetLx());
		sfxx(6) = f*(3.0/32.0*(1.0+xb)*(-2.0*xb+1.0+xb*xb)*GetLx()*GetLx()+3.0/32.0*pow(1.0+xb,2.0)*(-2.0+2.0*xb)*GetLx()*GetLx()+pow(1.0+xb,3.0)*GetLx()*GetLx()/32.0);
		*/
	}
}

//Get Shape Function d^2S(i)/dxi^2
double ANCFCable2D::GetS0xx(int i, const double& ploc) const
{
	double xb = ploc;
	double f = 4./Sqr(GetLx());

	if (!kappanodal)
	{
		switch (i)
		{
		case 1:
			return f*(3.0/2.0*xb);
		case 2:
			return f*((-2.0+6.0*xb)*GetLx()/8.0);
		case 3:
			return f*(-3.0/2.0*xb);
		case 4:
			return f*((-1.0+xb)*GetLx()/4.0+(1.0+xb)*GetLx()/2.0);
		default:
			mbs->UO() << "ANCFCable2D::GetS0xx Error: Shape function " << i << " not available!\n";
			return 0;
		}
	}
	else
	{
		//f included:
		switch (i)
		{
		case 1:
			return -0.25*(-60.0*xb+60.0*xb*xb*xb)/(GetLx()*GetLx());
		case 2:
			return -0.125*(12.0-60.0*xb-12.0*xb*xb+60.0*xb*xb*xb)/GetLx();
		case 3:
			return -0.25+0.75*xb+0.75*xb*xb-0.125E1*xb*xb*xb;
		case 4:
			return 0.15E1*(1.0+xb)*(8.0-9.0*xb+3.0*xb*xb)/(GetLx()*GetLx())+0.15E1*pow(1.0+xb,2.0)*(-9.0+6.0*xb)/(GetLx()*GetLx())+0.15E1*pow(1.0+xb,3.0)/(GetLx()*GetLx());
		case 5:
			return -0.75*(1.0+xb)*(5.0-8.0*xb+3.0*xb*xb)/GetLx()-0.75*pow(1.0+xb,2.0)*(-8.0+6.0*xb)/GetLx()-0.75*pow(1.0+xb,3.0)/GetLx();
		case 6:
			return -0.25-0.75*xb+0.75*xb*xb+0.125E1*xb*xb*xb;
		default:
			mbs->UO() << "ANCFCable2D::GetS0xx Error: Shape function " << i << " not available!\n";
			return 0;
		}
	}
}




	// -----------------------------------------------
	// Routines for computing positions/velocities...

// Get Shape Function Vector
//x->r_p(x), ploc € [-0.5*GetLx() .. 0.5*GetLx()]
Vector2D ANCFCable2D::GetVel2D(const double& p_loc) const
{
	double p0=p_loc/(0.5*GetLx());
	Vector2D p(0.,0.);
	int k = 0;
	for (int j = 1; j <= NS(); j++)
	{
		double sf = GetS0(j,p0);
		for (int i = 1; i <= Dim(); i++)
		{
		k++;
			p(i) += sf*XGP(k);
		}
	}
	return p;
};

Vector2D ANCFCable2D::GetVel2DD(const double& p_loc) const
{
	double p0=p_loc/(0.5*GetLx());
	Vector2D p(0.,0.);
	int k = 0;
	for (int j = 1; j <= NS(); j++)
	{
		double sf = GetS0(j,p0);
		for (int i = 1; i <= Dim(); i++)
		{
		k++;
			p(i) += sf*XGPD(k);
		}
	}
	return p;
};

//x->rD(x),p0 €[-1,1], ploc € [-GetLx()/2,GetLx()/2]
Vector2D ANCFCable2D::GetPos2DD(const double& p_loc) const
{
	//static Vector xgd;
	//xgd.SetLen(SOS());
	//GetDrawCoordinates(xgd);
	double p0=p_loc/(0.5*GetLx());
	////mbs->UO() << "p0=" << p0 << "\n";
	//static Vector SV;
	//GetS0(SV, p0);
	Vector2D p(0.,0.);
	// mbs->UO() << "p0 = ploc/0.5*GetLx() = " << p0 << "= 1/2*" << p_loc << "/" << GetLx() << "\n";
	int k = 0;
	for (int j = 1; j <= NS(); j++)
	{
		double sf = GetS0(j,p0);
		for (int i = 1; i <= Dim(); i++)
		{
			k++;
			p(i) += sf*XGD(k);
		}
	}
	//mbs->UO() << "p=" << p << "\n";
	return p;
};

//x->rD(x),p0 €[-1,1], ploc € [-GetLx()/2,GetLx()/2]
Vector2D ANCFCable2D::GetDisplacement2D(const double& p_loc, bool flagD) const
{
	double p0=p_loc/(0.5*GetLx());
	Vector2D p(0.,0.);
	int k = 0;
	for (int j = 1; j <= NS(); j++)
	{
		double sf = GetS0(j,p0);
		for (int i = 1; i <= Dim(); i++)
		{
			k++;
			if (flagD)
				p(i) += sf*(XGD(k) - x_init(k));
			else
				p(i) += sf*(XG(k) - x_init(k));
		}
	}
	return p;
};

//xi->rD(xi) in reference element coordinates xi(==p_loc) € [-1..1]
Vector2D ANCFCable2D::GetPos2D0D(const double& p_loc) const
{
	//static Vector xgd;
	//xgd.SetLen(SOS()); //8
	//GetDrawCoordinates(xgd);
	double p0=p_loc;
	//static Vector SV;
	//GetS0(SV, p0);
	Vector2D p(0.,0.);
	int k = 0;
	for (int j = 1; j <= NS(); j++)
	{	
		double sf = GetS0(j,p0);
		for (int i = 1; i <= Dim(); i++)
		{
			k++;
			p(i) += sf*XGD(k);
		}
	}
	//mbs->UO() << "p = " << p << ", ploc = " << p_loc << ", xgd = " << xgd << "\n";
	return p;
};



void ANCFCable2D::GetdAngle2Ddx(const Vector2D& ploc, double& dphidx)
{
	ConstVector<ANCFCable2D_MaxDOF> xg(SOS());
	for (int i = 1; i <= SOS(); i++)
		xg(i) = XG(i);  //e=e1..e8

	Vector2D rx = GetPosx2D(ploc.X(),xg);
	Vector2D rxx = GetPosxx2D(ploc.X(),xg);

	dphidx = (rxx.Y()*rx.X() - rx.Y()*rxx.X())/(Sqr(rx.X()) + Sqr(rx.Y()));
}


//for rotational constraint
void ANCFCable2D::GetdAngle2DdqT(const Vector2D& ploc, Matrix& d)
{
	ConstVector<ANCFCable2D_MaxDOF> xg(SOS());
	xg.SetLen(SOS());

	for (int i = 1; i <= SOS(); i++)
		xg(i) = XG(i);  //e=e1..e8

	d.SetSize(SOS(),1);
	d.FillWithZeros();

	double p0=ploc.X()/(0.5*GetLx());
	ConstVector<ANCFCable2D_MaxNShapes> SVx(NS());

	GetS0x(SVx, p0);
	Vector2D rx(0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			rx(i) += SVx(j)*xg((j-1)*Dim()+i);
		}
	}
	double rxn = Sqr(rx.Norm());
	if (rxn == 0) rxn = 1;

	for (int j = 1; j <= NS(); j++)
	{
		d((j-1)*Dim()+1, 1) = (-rx.Y()*SVx(j))/rxn;
		d((j-1)*Dim()+2, 1) = (SVx(j)*rx.X())/rxn;
	}

}

double ANCFCable2D::GetAngle2D(const Vector2D& ploc) const
{
	ConstVector<ANCFCable2D_MaxDOF> xg(SOS());

	for (int i = 1; i <= SOS(); i++)
		xg(i) = XG(i);  //e=e1..e8

	Vector2D rx = GetPosx2D(ploc.X(),xg);
	double ang = atan2(rx.Y(),rx.X());

	return ang;
}

double ANCFCable2D::GetAngle2DP(const Vector2D& ploc) const
{
	ConstVector<ANCFCable2D_MaxDOF> xg(SOS());

	for (int i = 1; i <= SOS(); i++)
		xg(i) = XG(i);  //e=e1..e8

	double p0=ploc.X()/(0.5*GetLx());

	ConstVector<ANCFCable2D_MaxNShapes> SVx(NS());
	GetS0x(SVx, p0);
	Vector2D rx(0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			rx(i) += SVx(j)*xg((j-1)*Dim()+i);
		}
	}

	double rxn = Sqr(rx.Norm());
	if (rxn == 0) rxn = 1;

	double angp = 0;
	for (int j = 1; j <= NS(); j++)
	{
		angp += XGP((j-1)*Dim()+1) * (-rx.Y()*SVx(j))/rxn;
		angp += XGP((j-1)*Dim()+2) * (SVx(j)*rx.X())/rxn;
	}

	return angp;
}


int ANCFCable2D::GetKappaMode() const {return kappamode;}

//does not work because of transformation:
double ANCFCable2D::ComputeCurvedLength(const Vector& xg) const
{
	//compute int_-1^+1 s_i(xi)*e0_i;
	/*
	double l1 = 0;

	double n = 10000;
	for (double i = 1; i <= n;i++)
	{
	double x1 = (i-1)/n;
	double x2 = (i)/n;

	Vector2D p1 = GetPos2D((x1-0.5)*GetLx(), xg);
	Vector2D p2 = GetPos2D((x2-0.5)*GetLx(), xg);
	l1+=(p1-p2).Norm();
	}
	*/
	Vector x1,w1;
	double l = 0;

	GetIntegrationRule(x1,w1,17); //

	for (int i1=1; i1<=x1.GetLen(); i1++)
	{
		double p = x1(i1);

		Vector2D rx = GetLx()*GetPosx2D(p*GetLx()*0.5, xg);
		double rxn = rx.Norm();

		l += 0.5*w1(i1)*rxn;
	}

	//mbs->UO() << "l1=" << l1 << ", l=" << l << "\n";
	//char str[32];
	//sprintf(str, "l=%.17g\n", l);
	//mbs->UO() << str;


	return l;
}





double ANCFCable2D::GetKappa(const double& x, const Vector& xg) const
{
	Vector2D rx = GetPosx2D(x*0.5*GetLx(),xg); //xg == e, x is element of -1,1, (=xi)
	Vector2D rxx = GetPosxx2D(x*0.5*GetLx(),xg); // XG(i)
	double rxn = rx.Norm();

	if (rxn == 0) return 0;
	else
	{
		double k;
		if (kappamode == 1)
			k = (rx.Cross(rxx))/Cub(rxn);
		else if (kappamode == 2 || kappamode == 4)
			k = (rx.Cross(rxx))/Sqr(rxn);
		else if (kappamode == 3)
			k = (rx.Cross(rxx))/(rxn);
		/* //not used in computation
		else if (kappamode == 4)
		{
		double h = GetSize().Y();
		double w = GetSize().Z();

		double I2 = w * Cub(h)/12.;
		double I4 = w * To5(h)/80.;

		k = (rx.Cross(rxx)) - 0.5*I4/I2*Cub((rx.Cross(rxx))/Sqr(rxn)); //GreenPiola - see paper Appendix
		}*/

		if (usekappafabs) return fabs(k);
		else return k;
	}
}

double ANCFCable2D::GetCurvature2D(const Vector2D& p_loc) const
{
	double x = p_loc.X()/(0.5*GetLx());
	ConstVector<ANCFCable2D_MaxDOF> xg(Dim()*NS());
	for (int i=1; i <= xg.Length(); i++) xg(i) = XG(i);

	return GetKappa(x, xg);
}

double ANCFCable2D::GetCurvature2DP(const Vector2D& p_loc) const
{
	double x = p_loc.X()/(0.5*GetLx());

	//Pos-DOF:
	ConstVector<ANCFCable2D_MaxDOF> xg(Dim()*NS());
	for (int i=1; i <= xg.Length(); i++) xg(i) = XG(i);

	//Vel-DOF:
	ConstVector<ANCFCable2D_MaxDOF> xgp(Dim()*NS());
	for (int i=1; i <= xgp.Length(); i++) xgp(i) = XGP(i);

	Vector2D rx = GetPosx2D(x*0.5*GetLx(),xg); //xg == e, x is element of -1,1, (=xi)
	Vector2D rxx = GetPosxx2D(x*0.5*GetLx(),xg); // XG(i)

	Vector2D rxp = GetPosx2D(x*0.5*GetLx(),xgp); //xg == e, x is element of -1,1, (=xi)
	Vector2D rxxp = GetPosxx2D(x*0.5*GetLx(),xgp); // XG(i)

	double rxn = rx.Norm();

	double k = 0;
	if (rxn == 0) 
	{
		return 0;
	}
	else
	{
		k = (rxp.Cross(rxx) + rx.Cross(rxxp))/(rx*rx) + (rx.Cross(rxx))*(-2.)*pow(rx*rx,-3.)*2.*(rxp*rx);
	}

	return -k;
}


double ANCFCable2D::GetKappaD(const double& x, const Vector& xg) const
{
	Vector2D rx = GetPosx2D(x*0.5*GetLx(),xg); //xg == e, x is element of -1,1, (=xi)
	Vector2D rxx = GetPosxx2D(x*0.5*GetLx(),xg); // XG(i)
	double rxn = rx.Norm();

	if (rxn == 0) return 0;
	else
	{
		double k;
		if (kappamode == 1)
			k = (rx.Cross(rxx))/Cub(rxn);
		else if (kappamode == 2 || kappamode == 4)
			k = (rx.Cross(rxx))/Sqr(rxn);
		else if (kappamode == 3)
			k = (rx.Cross(rxx))/(rxn);

		if (usekappafabs) return fabs(k);
		else return k;
	}
}

//epsilonxx=1/2*(r_x_vec.r_x_vec-1);
double ANCFCable2D::GetEpsAxial(const double& x, const Vector& xg) const
{
	Vector2D rx = GetPosx2D(x*0.5*GetLx(),xg);

	//1==Green, 2=Biot=Dmitrochenko, 3=linearized strain
	switch (epsmode)
	{
	case 1:
		return 0.5*(rx*rx-1.0);
	case 2:
		return (sqrt(rx*rx)-1.0);
	default:
		mbs->UO() << "ANCFCable2D::GetEpsAxial epsmode = " << epsmode << " not defined!\n";
		return 0;
	}
} 
double ANCFCable2D::GetEpsAxialD(const double& x, const Vector& xg) const
{
	Vector2D rx = GetPosx2D(x*0.5*GetLx(),xg);

	//1==Green, 2=Biot=Dmitrochenko, 3=linearized strain
	switch (epsmode)
	{
	case 1:
		return 0.5*(rx*rx-1.0);
	case 2:
		return (sqrt(rx*rx)-1.0);
	default:
		mbs->UO() << "ANCFCable2D::GetEpsAxialD epsmode = " << epsmode << " not defined!\n";
		return 0;
	}
} 
double ANCFCable2D::GetEps(const Vector2D& ploc) const
{
	double x = ploc.X();
	double x0 = x*2./GetLx();
	ConstVector<ANCFCable2D_MaxDOF> xgen(SOS());
	GetCoordinates(xgen);

	double eps = GetEpsAxial(x0, xgen);

	double yb = ploc.Y();
	double kappa = GetKappa(x0,xgen);
	return eps - yb*kappa;
} 

double ANCFCable2D::GetEpsInit(const Vector2D& ploc, int flagD) const
{
	double x = ploc.X();
	double x0 = x*2./GetLx();

	double eps;
	if (flagD) eps = GetEpsAxialD(x0, x_init);
	else eps = GetEpsAxial(x0, x_init);
	double yb = ploc.Y();
	double kappa;
	if (flagD) kappa = GetKappaD(x0,x_init);
	else kappa = GetKappa(x0, x_init);

	return eps - yb*kappa;
} 

double ANCFCable2D::GetEpsD(const Vector2D& ploc) const
{
	double x = ploc.X();
	double x0 = x*2./GetLx();
	ConstVector<ANCFCable2D_MaxDOF> xgD(SOS());
	GetDrawCoordinates(xgD);

	double eps = GetEpsAxialD(x0, xgD);

	double yb = ploc.Y();
	double kappa = GetKappaD(x0,xgD);
	return eps - yb*kappa;
} 


//x:-1...+1
double deps2 = 1e-9;
void ANCFCable2D::GetDeltaKappa(const double& x, const Vector& xg, const Vector& SVx, const Vector& SVxx, Vector& dkappa, double& kappa) const
{

	if (1)
	{
		int dim = Dim();
		int ns = NS();

		Vector2D rx(0.,0.);
		Vector2D rxx(0.,0.);
		for (int i = 1; i <= dim; i++)
		{
			for (int j = 1; j <= ns; j++)
			{
				rx(i) += SVx(j)*xg((j-1)*dim+i);
				rxx(i) += SVxx(j)*xg((j-1)*dim+i);
			}
		}

		double rxn = rx.Norm();

		if (rxn == 0) 
		{
			dkappa.SetAll(0); kappa=0; return;
		}  
		double rxcrxx = rx.Cross(rxx);
		double f;

		if (usekappafabs) f = fabs(rxcrxx);
		else f = rxcrxx;

		double g;
		if (kappamode == 1)
			g = Cub(rxn);
		else if (kappamode == 2 || kappamode == 4)
			g = Sqr(rxn);
		else if (kappamode == 3)
			g = (rxn);

			kappa = f/g; 
			if (kappamode == 4)
			{ //according to GreenPiola
				double h = GetSize().Y();
				double w = GetSize().Z();

				double I2 = w * Cub(h)/12.;
				double I4 = w * To5(h)/80.;

				kappa = rxcrxx - 0.5*I4/I2*Cub(rxcrxx/Sqr(rxn)); //GreenPiola - see paper Appendix
			}

		double g2inv = 1./Sqr(g);
		double fn = f*g2inv;
		if (kappamode == 1) fn *= (3.*rxn); //third power of rx
		else if (kappamode == 2 || kappamode == 4) fn *= 2; //only second power of rx
		else if (kappamode == 3) fn *= 1./rxn; //only first power of rx

		double gn = g*g2inv;
		if (f != 0 && usekappafabs) 
		{
			gn=gn/f;
		}

		double t1;
		double df, dg;

		for (int i=1; i <= dim; i++)
		{
			for (int j=1; j <= ns; j++)
			{
				//dr,x/de x r,xx + r,x x dr,xx/de: 
				//delta r,x x d^2(r)/dx^2 + ... =>paper 

				switch (i) {
					case 1:
						{
							t1 =  SVx(j)*rxx.Y()-SVxx(j)*rx.Y(); break; //dy=dz=0
						}
					case 2:
						{
							t1 = -SVx(j)*rxx.X()+SVxx(j)*rx.X(); break; //dx=dz=0
						}
						//correct:
						//					case 1:
						//						t1 = +SVx(j)*rxx.Y()-SVxx(j)*rx.Y(); break; //dy=dz=0
						//					case 2:
						//						t1 = -SVx(j)*rxx.X()+SVxx(j)*rx.X(); break; //dx=dz=0

					default: ;
				}

				dg = (rx(i)*SVx(j)); //normed

				if (usekappafabs) df = (rxcrxx*t1); //normed
				else df = t1;

				//dkappa((j-1)*dim+i) = (t1*g*Sgn(rxcrxx)-f*dg*(3.*rxn))/(g*g);
				dkappa((j-1)*dim+i) = df*gn-fn*dg; //normed
			}
		}
		//oldkappa = dkappa;
	}
	else
	{
		ConstVector<ANCFCable2D_MaxDOF> xgdiff(SOS());
		xgdiff = xg;
		double val0;
		val0 = GetKappa(x,xgdiff);
		kappa = val0;

		for (int i=1; i <= SOS(); i++)
		{
			xgdiff(i) += deps2;
			dkappa(i) = GetKappa(x,xgdiff);
			dkappa(i) -= val0;
			dkappa(i) *= 1./deps2;
			xgdiff(i) -= deps2;
		}

		/*
		Vector diff = dkappa - oldkappa;
		if (GetMBS()->GetTime() > 0.99 && diff.GetNorm() > 1e-8)
		{
		mbs->UO() << "kappaanal=" << oldkappa << "\n";
		mbs->UO() << "kappanum=" << dkappa << "\n";
		}
		dkappa = oldkappa;
		*/
	}

}
void ANCFCable2D::GetDeltaDeltaKappa(const double& x, const Vector& xg, const Vector& SVx, const Vector& SVxx, 
																		 double& kappa, Vector& dkappa, Matrix& ddkappa) const
{
	if (1)
	{
		kappa = 0; dkappa.SetAll(0); ddkappa.SetAll(0.); 
		int dim = Dim();
		int ns = NS();

		Vector2D rx(0.,0.), rxx(0.,0.);

		for (int i = 1; i <= dim; i++)
		{
			for (int j = 1; j <= ns; j++)
			{
				rx(i) += SVx(j)*xg((j-1)*dim+i);
				rxx(i) += SVxx(j)*xg((j-1)*dim+i);
			}
		}
		double rxn = rx.Norm();

		// kappa = f/g
		double f, g;
		if (rxn == 0) 
		{
			return;
		}  
		double rxcrxx = rx.Cross(rxx);

		f = rxcrxx;
		g = Sqr(rxn);
		kappa = f/g;

		// dkappa = (df g - dg f)/g^2
		ConstVector<ANCFCable2D_MaxDOF> df(SOS()), dg(SOS());
		df.SetAll(0);
		dg.SetAll(0);

		// ddkappa = [ (ddf_ij g + df_i dg_j - df_j dg_i - f ddg_ij) g^2 - (df_i g - dg_i f) 2 dg_j g ] / g^4
		ConstMatrix<ANCFCable2D_MaxDOF*ANCFCable2D_MaxDOF> ddf, ddg;
		ddf.SetSize(SOS(), SOS());
		ddg.SetSize(SOS(), SOS());
		ddf.SetAll(0.);
		ddg.SetAll(0.);

		double matrixentry;
		for (int di=1; di <= ns; di++)
		{
			// Shape = (S(di), 0)
			df((di-1)*dim+1) = SVx(di)*rxx(2) - SVxx(di)*rx(2);
			dg((di-1)*dim+1) = 2* rx(1)*SVx(di);
			// Shape = (0, S(di))
			df((di-1)*dim+2) = SVxx(di)*rx(1) - SVx(di)*rxx(1);
			dg((di-1)*dim+2) = 2* rx(2)*SVx(di);
			for (int dj=1; dj <= di; dj++)
			{
				// i==1, j==2
				matrixentry = SVx(di)*SVxx(dj) - SVxx(di)*SVx(dj);
				ddf((di-1)*dim+1, (dj-1)*dim+2) = matrixentry;
				ddf((dj-1)*dim+2, (di-1)*dim+1) = matrixentry;
				// i==2, j==1)
				matrixentry = SVxx(di)*SVx(dj) - SVx(di)*SVxx(dj);
				ddf((di-1)*dim+2, (dj-1)*dim+1) = matrixentry;
				ddf((dj-1)*dim+1, (di-1)*dim+2) = matrixentry;
				//i == j == 1
				matrixentry = 2*SVx(di)*SVx(dj);
				ddg((di-1)*dim+1, (dj-1)*dim+1) = matrixentry;
				ddg((dj-1)*dim+1, (di-1)*dim+1) = matrixentry;
				// i == j == 2
				ddg((di-1)*dim+2, (dj-1)*dim+2) = matrixentry;
				ddg((dj-1)*dim+2, (di-1)*dim+2) = matrixentry;


			}
		}

		for (int ii=1; ii<=SOS(); ii++)
		{
			dkappa(ii) = (df(ii)*g - dg(ii)*f)/(g*g);
			double g4 = g*g*g*g;
			for (int jj=1; jj<=ii; jj++)
			{
				matrixentry = ( (ddf(ii,jj)*g+df(ii)*dg(jj)-df(jj)*dg(ii)-ddg(ii,jj)*f)*g*g - (df(ii)*g-dg(ii)*f)*2*dg(jj)*g) / g4;
				ddkappa(jj,ii) = matrixentry;
				ddkappa(ii,jj) = matrixentry;
			}
		}

		//mbs->UO() << "ddkappa = " << ddkappa << "\n";

	}
	else
	{
		GetDeltaKappa(x, xg, dkappa, kappa);
		ConstVector<ANCFCable2D_MaxDOF> xgdiff(SOS()), dkappa2(SOS());
		xgdiff = xg;
		dkappa.SetLen(SOS());
		dkappa.SetAll(0.);
		dkappa2.SetLen(SOS());
		dkappa2.SetAll(0.);

		double dummy;
		double deps2=1e-9;
		for (int i=1; i <= SOS(); i++)
		{
			xgdiff(i) += deps2;
			GetDeltaKappa(x,xgdiff, dkappa2, dummy);
			dkappa2 -= dkappa;
			dkappa2 *= 1./deps2;
			for (int j=1; j<=SOS(); j++)
				ddkappa(i,j) = dkappa2(j);
			xgdiff(i) -= deps2;
		}
		//mbs->UO() << "ddkappa2=" << ddkappa << "\n";
	}
}

//+
void ANCFCable2D::GetDeltaEpsAxial(const double& x, const Vector& xg, const Vector& SVx, Vector& depsaxial) const
{
	if (1)
	{
		int dim = Dim();
		int ns = NS();
		Vector2D rx(0.,0.);
		for (int i = 1; i <= dim; i++)
		{
			for (int j = 1; j <= ns; j++)
			{
				rx(i) += SVx(j)*xg((j-1)*dim+i);
			}
		}
		double rxn = rx.Norm();
		for (int i=1; i <= dim; i++)
		{
			for (int j=1; j <= ns; j++)
			{
				if (epsmode==1)
					depsaxial((j-1)*dim+i) = SVx(j)*rx(i);
				else if (epsmode==2)
					depsaxial((j-1)*dim+i) = 1. / rxn * SVx(j)*rx(i);
			}
		}
		//mbs->UO() << "deps1=" << depsaxial << ", rx=" << rx << "\n";
	}
	else
	{
		ConstVector<ANCFCable2D_MaxDOF> xgdiff(SOS());
		xgdiff = xg;
		double val0;
		val0 = GetEpsAxial(x,xgdiff);

		for (int i=1; i <= SOS(); i++)
		{
			xgdiff(i) += 0.1*deps2;
			depsaxial(i) = GetEpsAxial(x,xgdiff);
			depsaxial(i) -= val0;
			depsaxial(i) *= 1./deps2/0.1;
			xgdiff(i) -= deps2*0.1;
		}

		//mbs->UO() << "deps2=" << depsaxial << "\n";
	}
}

void ANCFCable2D::GetDeltaDeltaEpsAxial(const double& x, const Vector& xg, const Vector& SVx, 
																				double& eps, Vector& deltaeps, Matrix& ddepsaxial) const
{
	if (1)
	{
		int dim = Dim();
		int ns = NS();
		Vector2D rx(0.,0.);
		for (int i = 1; i <= dim; i++)
		{
			for (int j = 1; j <= ns; j++)
			{
				rx(i) += SVx(j)*xg((j-1)*dim+i);
			}
		}

		eps = 0;
		deltaeps.SetAll(0.);
		ddepsaxial.SetAll(0.);

		double rxrx = rx*rx;
		double rxn = sqrt(rxrx);

		// eps
		if (epsmode == 1) eps = 0.5*(rxrx-1.0);
		if (epsmode == 2) eps = (rxn-1.0);

		double matrixentry;
		for (int i=1; i <= dim; i++)
		{
			for (int di=1; di <= ns; di++)
			{
				// depsax
				if (epsmode==1)
					deltaeps((di-1)*dim+i) = SVx(di)*rx(i);
				else if (epsmode==2)
					deltaeps((di-1)*dim+i) = 1. / rxn * SVx(di)*rx(i);

				for (int j=1; j <= dim; j++)
				{
					for (int dj=1; dj <= di; dj++)
					{
						// ddepsax(i,j) = S_i ^T S_j
						if (epsmode==1 && i==j)
						{
							matrixentry = SVx(di)*SVx(dj);
							ddepsaxial((di-1)*dim+i, (dj-1)*dim+j) = matrixentry;
							ddepsaxial((dj-1)*dim+i, (di-1)*dim+j) = matrixentry;
						}
						// ddepsax(i,j) = - 1/|r|^3 (s_i^T r) (S_j^T r) + 1/|r| S_i^T S_j
						else if (epsmode==2 && i!=j)
						{
							matrixentry =  -1. / (rxrx*rxn) * SVx(di)*rx(i) * SVx(dj)*rx(j);
							ddepsaxial((di-1)*dim+i, (dj-1)*dim+j) = matrixentry;
							ddepsaxial((dj-1)*dim+j, (di-1)*dim+i) = matrixentry;
						}
						if (epsmode==2 && i==j)
						{
							matrixentry =  -1. / (rxrx*rxn) * SVx(di)*rx(i) * SVx(dj)*rx(j) + 1./rxn * SVx(di)*SVx(dj);
							ddepsaxial((di-1)*dim+i, (dj-1)*dim+j) = matrixentry;
							ddepsaxial((dj-1)*dim+j, (di-1)*dim+i) = matrixentry;
						}
					}
				}
			}
		}
		//mbs->UO() << "deps1=" << depsaxial << ", rx=" << rx << "\n";
	}
	else
	{
		eps = GetEpsAxial(x,xg);
		GetDeltaEpsAxial(x,xg,deltaeps);

		ConstVector<ANCFCable2D_MaxDOF> xgdiff(SOS()), depsax(SOS()), depsax2(SOS());
		xgdiff = xg;
		depsax.SetAll(0.);
		depsax2.SetAll(0.);

		GetDeltaEpsAxial(x, xg, depsax);
		double deps2=1e-9;
		for (int i=1; i <= SOS(); i++)
		{
			xgdiff(i) += deps2;
			GetDeltaEpsAxial(x,xgdiff, depsax2);
			depsax2 -= depsax;
			depsax2 *= 1./deps2;
			for (int j=1; j<=SOS(); j++)
				ddepsaxial(i,j) = depsax2(j);
			xgdiff(i) -= deps2;
		}
		//mbs->UO() << "deps2=" << depsaxial << "\n";
	}
}

//++++++++++++++++++++++++ Plasticity ++++++++++++++++++++++
// x0 € [-GetLx()/2, GetLx()/2]
double ANCFCable2D::GetPlasticAxStrain(double x0, int flagD) const    // integral mean value of plastic strains due to plastic cells
{
	if (!GetMaterial().IsInelasticMaterial()) return 0.;

	double eps_ax = 0;
	double hy = GetLy()/(plasticstrains_height-1.);
	// Matrix contains plastic strain matrix
	Matrix plasticstrain_mat;
	GetPlasticStrainMatrix(plasticstrain_mat, flagD);
	// Vector plasticstrains contains plastic strains in gridpoints over height interpolated at x-coord x0
	ConstVector<ANCFCable2D_maxppy> plasticstrains(plasticstrains_height);
	plasticstrain_mat.InterpolateX(x0/GetLx()*2., plasticstrains);

	// Integrate plasticstrain over height - trapezoidal rule
	for (int cell=1; cell<plasticstrains_height; cell++)
	{
		eps_ax += 0.5*(plasticstrains(plasticstrains_height-(cell-1))+plasticstrains(plasticstrains_height-cell));
	}
	return eps_ax/(plasticstrains_height-1.);
}

// ip is integration point number, € 1 ... plasticstrains_width
double ANCFCable2D::GetPlasticAxStrainIP(int ip, int flagD) const    // integral mean value of plastic strains due to plastic cells
{
	if (!GetMaterial().IsInelasticMaterial()) return 0.;

	Matrix plasticstrains;
	GetPlasticStrainMatrix(plasticstrains);
	double eps_ax = 0;
	for (int cell=1; cell<plasticstrains_height; cell++)
	{
		eps_ax += plasticstrains(cell, ip);
	}
	return eps_ax/(plasticstrains_height-1.);
}

double ANCFCable2D::GetPlasticKappa(double x0, int flagD) const				// negative first moment of plastic strains due to plastic cells
{
	if (!GetMaterial().IsInelasticMaterial()) return 0.;
	double kappa = 0;
	double hy = GetLy()/(plasticstrains_height-1.);
	// Matrix contains plastic strain matrix
	Matrix plasticstrain_mat;
	GetPlasticStrainMatrix(plasticstrain_mat, flagD);
	// Vector plasticstrains contains plastic strains in gridpoints over height interpolated at x-coord x0
	ConstVector<ANCFCable2D_maxppy> plasticstrains(plasticstrains_height);
	plasticstrain_mat.InterpolateX(x0/GetLx()*2., plasticstrains);

	// Integrate plasticstrain*y over height - exact formula in each plastic cell
	for (int cell=1; cell<plasticstrains_height; cell++)
	{
		double y1 = hy*(cell-1)-0.5*GetLy();
		double y2 = hy*(cell)-0.5*GetLy();
		double eps1 = plasticstrains(plasticstrains_height-(cell-1));
		double eps2 = plasticstrains(plasticstrains_height-cell);
		kappa += 2./Cub(GetLy()) * ( eps2*(2*y2*y2-y1*y1) + (eps1-eps2)*y1*y2 + eps1*(y2*y2-2*y1*y1) );
	}
	return -kappa;
}

double ANCFCable2D::GetPlasticStrain(Vector2D& ploc, int flagD) const  // inter/extrapolated value of plastic strain at local coord. ploc due to plastic cells
{
	if (!GetMaterial().IsInelasticMaterial()) return 0.;
	Vector2D ploc_unit(ploc.X()/GetLx()*2, 2*ploc.Y()/GetLy());
	const Matrix plasticstrains;
	GetPlasticStrainMatrix(plasticstrains, flagD);
	return plasticstrains.Interpolate(ploc_unit);
}


void ANCFCable2D::SetPlasticStrainMatrix( const Matrix& initstrains, int height, int width)   // initialize plastic strain matrix
{
	plasticstrains_height = height;
	plasticstrains_width = width;
	Vector datainit(DataS()); //automatically initialized with zeros, set to correct values from initstrains

	double hy = 2./(height-1.), hx = 2./(width-1.);
	for (int i=1; i<=height; i++)
		for (int j=1; j<=width; j++)
		{
			Vector2D ploc( (j-1)*hx-1, 1-(i-1)*hy );
			datainit( (i-1)*(width) + j ) = initstrains.Interpolate(ploc);
		}

	SetDataInit(datainit); //initialize data variables with initial inelastic strain

}
void ANCFCable2D::SetPlasticStrainMatrix(int height, int width)   // initialize plastic strain matrix
{
	plasticstrains_height = height;
	plasticstrains_width = width;
	Vector datainit(DataS()); //automatically initialized with zeros
	SetDataInit(datainit); //initialize data variables with zero = initial inelastic strain == zero
}

void ANCFCable2D::GetPlasticStrainMatrix(Matrix& mat, int flagD)
{
	if (GetMaterial().IsInelasticMaterial() && !flagD)
	{
		mat.LinkWith(plasticstrains_height, plasticstrains_width, &XData(1));
	}
	else if (GetMaterial().IsInelasticMaterial() && flagD)
	{
		mat.LinkWith(plasticstrains_height, plasticstrains_width, &XDataD(1));
	}
}

void ANCFCable2D::GetPlasticStrainMatrix(const Matrix& mat, int flagD) const
{
	if (GetMaterial().IsInelasticMaterial() && !flagD)
	{
		ANCFCable2D& constthis = const_cast<ANCFCable2D&>(*this);
		Matrix& constmat = const_cast<Matrix&>(mat);
		constmat.LinkWith(plasticstrains_height, plasticstrains_width, &(constthis.XData(1)) );
	}
	else if (GetMaterial().IsInelasticMaterial() && flagD)
	{
		ANCFCable2D& constthis = const_cast<ANCFCable2D&>(*this);
		Matrix& constmat = const_cast<Matrix&>(mat);
		constmat.LinkWith(plasticstrains_height, plasticstrains_width, &(constthis.XDataD(1)) );
	}
}

void ANCFCable2D::GetPlasticStrainMatrixLastStep(const Matrix& mat) const
{
	ANCFCable2D& constthis = const_cast<ANCFCable2D&>(*this);
	Matrix& constmat = const_cast<Matrix&>(mat);
	double* constptr = const_cast<double*>(& GetMBS()->GetLastDataVector()(LTGdata(1)) );
	constmat.LinkWith(plasticstrains_height, plasticstrains_width, constptr );
}

#ifdef AND_INCLUDES
// Transport scheme only available for AND
void ANCFCable2D::GetPlasticStrainMatrixLastStepTransported(const Matrix& mat) const
{
	ANCFCable2D& constthis = const_cast<ANCFCable2D&>(*this);
	Matrix& constmat = const_cast<Matrix&>(mat);
	// last step data vector, AND: last step + transport
	double* constptr = const_cast<double*>(& GetMBS()->GetTempDataVector()(LTGdata(1)) );
	constmat.LinkWith(plasticstrains_height, plasticstrains_width, constptr );
}
#endif // AND_INCLUDES

// matrix containing internal variable xi for hardening, xi = int_0^t |plasticstrainP| ds
void ANCFCable2D::GetInternalVariableMatrix( const Matrix& mat, int flagD) const
{
	if (GetMaterial().IsInelasticMaterial() && !flagD)
	{
		ANCFCable2D& constthis = const_cast<ANCFCable2D&>(*this);
		Matrix& constmat = const_cast<Matrix&>(mat);
		constmat.LinkWith(plasticstrains_height, plasticstrains_width, &(constthis.XData(plasticstrains_height*plasticstrains_width+1)) );
	}
	else if (GetMaterial().IsInelasticMaterial() && flagD)
	{
		ANCFCable2D& constthis = const_cast<ANCFCable2D&>(*this);
		Matrix& constmat = const_cast<Matrix&>(mat);
		constmat.LinkWith(plasticstrains_height, plasticstrains_width, &(constthis.XDataD(plasticstrains_height*plasticstrains_width+1)) );
	}
}
void ANCFCable2D::GetInternalVariableMatrix( Matrix& mat, int flagD)
{
	if (GetMaterial().IsInelasticMaterial() && !flagD)
	{
		mat.LinkWith(plasticstrains_height, plasticstrains_width, &XData(plasticstrains_height*plasticstrains_width+1));
	}
	else if (GetMaterial().IsInelasticMaterial() && flagD)
	{
		mat.LinkWith(plasticstrains_height, plasticstrains_width, &XDataD(plasticstrains_height*plasticstrains_width+1));
	}
}
void ANCFCable2D::GetInternalVariableMatrixLastStep(const Matrix& mat) const
{
	ANCFCable2D& constthis = const_cast<ANCFCable2D&>(*this);
	Matrix& constmat = const_cast<Matrix&>(mat);
	double* constptr = const_cast<double*>(& GetMBS()->GetLastDataVector()(LTGdata(plasticstrains_height*plasticstrains_width+1)) );
	constmat.LinkWith(plasticstrains_height, plasticstrains_width, constptr );
}

#ifdef AND_INCLUDES
// Transport scheme only available for AND
void ANCFCable2D::GetInternalVariableMatrixLastStepTransported(const Matrix& mat) const
{
	ANCFCable2D& constthis = const_cast<ANCFCable2D&>(*this);
	Matrix& constmat = const_cast<Matrix&>(mat);
	// last step data vector, AND: last step + transport
	double* constptr = const_cast<double*>(& GetMBS()->GetTempDataVector()(LTGdata(plasticstrains_height*plasticstrains_width+1)) );
	constmat.LinkWith(plasticstrains_height, plasticstrains_width, constptr );
}
#endif // AND_INCLUDES


// inter/extrapolated value of plastic strain at local coord. ploc due to plastic cells
double ANCFCable2D::GetInternalVariable(Vector2D& ploc, int flagD) const
{
	if (!GetMaterial().IsInelasticMaterial()) return 0.;
	Vector2D ploc_unit(ploc.X()/GetLx()*2, 2*ploc.Y()/GetLy());
	const Matrix internalvar;
	GetInternalVariableMatrix(internalvar, flagD);
	return internalvar.Interpolate(ploc_unit);
}



// Energy
double ANCFCable2D::GetKineticEnergy()
{
	//xg.SetLen(SOS());
	//for (int i = 1; i <= SOS(); i++)
	//	xg(i) = XGP(i);
	ConstVector<ANCFCable2D_MaxDOF> xgp(SOS()), temp(SOS());
	GetCoordinatesP(xgp);

	if (massmatrix.Getcols() != SOS()) return 0;

	Mult(massmatrix,xgp,temp);

	double EK = 0.5 * temp * xgp;

	return EK;
}

double ANCFCable2D::GetPotentialEnergy()
{
	//0.5 * int_L EI*kappa^2 dx + 
	//0.5 * int_L EA*eps^2 dx

	double E_pot = Body2D::GetPotentialEnergy();

	int sos = SOS();
	int ns = NS();
	int dim = Dim();


	double EI_w = Em() * GetLz()*Cub(GetLy())/12.;
	double EA = Em() * GetLy() * GetLz(); 
	double A = GetLy()*GetLz();//for Hellinger Reissner
	if (IsBeamParameters())
	{
		EI_w = GetBeamEIy();
		EA = GetBeamEA();
	}

	int order_kappa = 5+4*kappanodal; //regular: 5
	int order_eps   = 9+4*kappanodal; //regular: 9
	ConstVector<ANCFCable2D_MaxIP> xk, wk, xe, we;

	GetIntegrationRule(xk,wk,order_kappa); 
	GetIntegrationRule(xe,we,order_eps);   

	//compute kappa0 and eps0
	ConstVector<ANCFCable2D_MaxIP> kappa0, eps0;
	kappa0.SetLen(xk.Length());
	eps0.SetLen(xe.Length());

	int usedmitrochenko = 0;  //for comparison with Euler elastica
	if (epsmode == 2) usedmitrochenko = 1;

	ConstVector<ANCFCable2D_MaxDOF> xg(SOS());
	GetCoordinates(xg);
	//xg.SetLen(sos);

	for (int i1=1; i1<=xk.GetLen(); i1++)
	{
		double x = xk(i1);
		kappa0(i1) = GetKappa(x,x_init);
		//UO() << "K=" << kappa0(i1) << "\n";
	}

	for (int i1=1; i1<=xe.GetLen(); i1++)
	{
		double x = xe(i1);
		//eps0(i1) = GetEpsAxial(x,xg);

		//1==Green, 2=Biot=Dmitrochenko
		Vector2D rx = GetPosx2D(x*0.5*GetLx(),x_init);
		if (epsmode == 2)
		{
			eps0(i1) = (sqrt(rx*rx)-1.0);
		}
		else if (epsmode == 1)
		{
			eps0(i1) = 0.5*(rx*rx-1.0);
		}
	}

	//SetComputeCoordinates();

	//KAPPA:   0.5 * int_L EI*kappa^2 dx
	for (int i1=1; i1<=xk.GetLen(); i1++)
	{
		double x = xk(i1);

		Vector2D rx0 = GetPosx2D(x*GetLx()*0.5, x_init);
		double rxn = rx0.Norm();
		double det = 0.5*GetLx()*rxn; //is approximately 0.5*GetLx(); in straight elements exactly 0.5*GetLx()!!!

		double factkappa = EI_w*wk(i1)*det; //?!!!!Achtung: Determinante der Elementransformation
		double kappa = GetKappa(x, xg);

		E_pot += 0.5 * factkappa * Sqr(kappa - kappa0(i1));
	}


	//EPS:     0.5 * int_L EA*eps^2 dx
	for (int i1=1; i1<=xe.GetLen(); i1++)
	{
		double x = xe(i1);

		Vector2D rx0 = GetPosx2D(x*GetLx()*0.5, x_init);
		double rxn = rx0.Norm();

		double det = 0.5*GetLx()*rxn; //is approximately 0.5*GetLx(); in straight elements exactly 0.5*GetLx()!!!
		double facteps = EA*we(i1)*det;

		double eps = GetEpsAxial(x,xg);
		E_pot += 0.5 * facteps * Sqr(eps - eps0(i1));
	}

	return E_pot;
}

 void ANCFCable2D::GetIntDkappaDq2D(Vector& dudq)
 {
 
 	if (kappamode != 2) UO() << "Error: PiezoCable2D, kappmode is wrong!\n";
 
 	int order_kappa = 5+4*kappanodal; //regular: 5
	ConstVector<ANCFCable2D_MaxIP> xk, wk;
 	
	ConstVector<ANCFCable2D_MaxDOF> temp(SOS()), xg(SOS());
 	temp.SetLen(SOS());
 	temp.SetAll(0.);
 	dudq.SetAll(0.);
 
 	GetIntegrationRule(xk,wk,order_kappa); 
 
	GetCoordinates(xg);
 	//SetComputeCoordinates();
 
 	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 	//kappa*delta kappa:
 
 	for (int i1=1; i1<=xk.GetLen(); i1++)
 	{
 		double x = xk(i1);
 
 		Vector2D rx0 = GetPosx2D(x*GetLx()*0.5, x_init);
 		double rxn = rx0.Norm();
 		double det = 0.5*GetLx()*rxn; //is approximately 0.5*GetLx(); in straight elements exactly 0.5*GetLx()!!!
 
 		//bending moments:
 		double factkappa = wk(i1)*det; //?!!!!Achtung: Determinante der Elementransformation
 
 		double kappa;
 		GetDeltaKappa(x,xg,temp, kappa);
 
 		temp *= factkappa;
 		dudq += temp;
 	}
 
 }






//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double ANCFPiezoCable2D::EvaluatePiezoMoment(double x_loc)
{
	//return 0;

	double EI_w = Em() * GetLz()*Cub(GetLy())/12.;
	double EA = Em() * GetLy() * GetLz(); 
	double A = GetLy()*GetLz();//for Hellinger Reissner
	if (IsBeamParameters())
	{
		EI_w = GetBeamEIy();
		EA = GetBeamEA();
	}

	if (GetPar(3) == 1) //one beam
	{
		double ang = GetMBS()->GetSensor(sensors[0]).GetCurrentValueWithSensorProcessing(mbs->GetTime());
		double angvel = GetMBS()->GetSensor(sensors[1]).GetCurrentValueWithSensorProcessing(mbs->GetTime());
		double angacc = GetMBS()->GetSensor(sensors[2]).GetCurrentValueWithSensorProcessing(mbs->GetTime());

		double xacc = 0;
		double yacc = 0;
		double g = 9.81;
		//if (!GetMBS()->IsJacobianComputation()) UO() << "sens: x=" << ang << ", v=" << angvel << ", a=" << angacc << "\n";

		//return 0;
		double rhoA = Rho()*GetLy()*GetLz();
		double L = GetPar(1);
		double x = (x_loc+1.)*0.5*GetLx() + (double)(GetPar(2)-1)*GetLx();
		double xi = x/L;
		double S1 = Sqr(1.-xi);
		double S2 = 0.5*Sqr(1.-xi)*(2.+xi);

		double Ma = -(0.5*rhoA*Sqr(L)*( (xacc+g)*sin(ang)+yacc*cos(ang) )*S1 - 1./3.*rhoA*Cub(L)*angacc*S2);

		return Ma;
	}
	else if (GetPar(3) == 2) //two-arm
	{
		double phi1    = GetMBS()->GetSensor(sensors[0]).GetCurrentValueWithSensorProcessing(mbs->GetTime()); //angle at support of beam1
		double phi1_t  = GetMBS()->GetSensor(sensors[1]).GetCurrentValueWithSensorProcessing(mbs->GetTime()); //angular velocity at x=0 of beam1
		double phi1_tt = GetMBS()->GetSensor(sensors[2]).GetCurrentValueWithSensorProcessing(mbs->GetTime()); //angular acceleration at x=0 of beam1

		double phi2    = -(GetMBS()->GetSensor(sensors[3]).GetCurrentValueWithSensorProcessing(mbs->GetTime()) - GetMBS()->GetSensor(sensors[6]).GetCurrentValueWithSensorProcessing(mbs->GetTime())); //angle at support of beam2
		double phi2_t  = -(GetMBS()->GetSensor(sensors[4]).GetCurrentValueWithSensorProcessing(mbs->GetTime()) - GetMBS()->GetSensor(sensors[7]).GetCurrentValueWithSensorProcessing(mbs->GetTime())); //angular velocity at x=0 of beam2
		double phi2_tt = -(GetMBS()->GetSensor(sensors[5]).GetCurrentValueWithSensorProcessing(mbs->GetTime()) - GetMBS()->GetSensor(sensors[8]).GetCurrentValueWithSensorProcessing(mbs->GetTime())); //angular acceleration at x=0 of beam2

		double m3 = GetPar(4); //Endmass
		double A1 = GetLy()*GetLz();
		double A2 = GetLy()*GetLz();
		double L1 = GetPar(1);
		double L2 = GetPar(1);

		double x = (x_loc+1.)*0.5*GetLx() + (double)(GetPar(2)-1)*GetLx();

		double Ma = 0;
		if (GetPar(5) == 1)
		{ //first arm
			double xi = x/L1;
			//++++++++++++++++++++++++++++
double MapleGenVar3 = L1;
      double MapleGenVar5 = 1.0-xi;
      double MapleGenVar7 = -sin(phi2)*(-m3*(-cos(phi1+phi2)*L1*cos(phi1)*phi1_t*
phi1_t-cos(phi1+phi2)*L1*sin(phi1)*phi1_tt-sin(phi1+phi2)*L1*sin(phi1)*phi1_t*
phi1_t+sin(phi1+phi2)*L1*cos(phi1)*phi1_tt-L2*phi1_t*phi1_t-2.0*L2*phi1_t*
phi2_t-L2*phi2_t*phi2_t)+L2*(-Rho()*A2*(-L2*phi1_t*phi1_t-2.0*L2*phi1_t*phi2_t-L2
*phi2_t*phi2_t)/2.0-Rho()*A2*(-cos(phi1+phi2)*L1*cos(phi1)*phi1_t*phi1_t-cos(phi1
+phi2)*L1*sin(phi1)*phi1_tt-sin(phi1+phi2)*L1*sin(phi1)*phi1_t*phi1_t+sin(phi1+
phi2)*L1*cos(phi1)*phi1_tt)));
      double MapleGenVar8 = cos(phi2)*(-m3*(-sin(phi1+phi2)*L1*cos(phi1)*phi1_t*phi1_t
-sin(phi1+phi2)*L1*sin(phi1)*phi1_tt-L2*phi1_tt-L2*phi2_tt+cos(phi1+phi2)*L1*
sin(phi1)*phi1_t*phi1_t-cos(phi1+phi2)*L1*cos(phi1)*phi1_tt)+L2*(-Rho()*A2*(-L2*
phi1_tt-L2*phi2_tt)/2.0-Rho()*A2*(-sin(phi1+phi2)*L1*cos(phi1)*phi1_t*phi1_t-sin(
phi1+phi2)*L1*sin(phi1)*phi1_tt+cos(phi1+phi2)*L1*sin(phi1)*phi1_t*phi1_t-cos(
phi1+phi2)*L1*cos(phi1)*phi1_tt)));
      double MapleGenVar6 = MapleGenVar7+MapleGenVar8;
      double MapleGenVar4 = MapleGenVar5*MapleGenVar6;
      double MapleGenVar2 = MapleGenVar3*MapleGenVar4;
      MapleGenVar3 = L1*L1*(Rho()*A1*L1*phi1_tt*(1.0-xi*xi*xi)/3.0-Rho()*A1*xi*L1*
phi1_tt*(1.0-xi*xi)/2.0);
      double MapleGenVar1 = MapleGenVar2+MapleGenVar3;
      Ma = MapleGenVar1-L2*m3*(-sin(phi1+phi2)*L1*cos(phi1)*phi1_t*phi1_t-sin(
phi1+phi2)*L1*sin(phi1)*phi1_tt-L2*phi1_tt-L2*phi2_tt+cos(phi1+phi2)*L1*sin(
phi1)*phi1_t*phi1_t-cos(phi1+phi2)*L1*cos(phi1)*phi1_tt)+L2*L2*(-Rho()*A2*(-L2*
phi1_tt-L2*phi2_tt)/3.0-Rho()*A2*(-sin(phi1+phi2)*L1*cos(phi1)*phi1_t*phi1_t-sin(
phi1+phi2)*L1*sin(phi1)*phi1_tt+cos(phi1+phi2)*L1*sin(phi1)*phi1_t*phi1_t-cos(
phi1+phi2)*L1*cos(phi1)*phi1_tt)/2.0);
		}
		else
		{
			double xi = x/L2;
			//++++++++++++++++++++++++++++
double MapleGenVar1 = -L2*(1.0-xi)*m3*(-sin(phi1+phi2)*L1*cos(phi1)*phi1_t*
phi1_t-sin(phi1+phi2)*L1*sin(phi1)*phi1_tt-L2*phi1_tt-L2*phi2_tt+cos(phi1+phi2)
*L1*sin(phi1)*phi1_t*phi1_t-cos(phi1+phi2)*L1*cos(phi1)*phi1_tt);
      double MapleGenVar2 = L2*L2*(-Rho()*A2*(-L2*phi1_tt-L2*phi2_tt)*(1.0-xi*xi*xi)/3.0
+(xi*Rho()*A2*(-L2*phi1_tt-L2*phi2_tt)-Rho()*A2*(-sin(phi1+phi2)*L1*cos(phi1)*
phi1_t*phi1_t-sin(phi1+phi2)*L1*sin(phi1)*phi1_tt+cos(phi1+phi2)*L1*sin(phi1)*
phi1_t*phi1_t-cos(phi1+phi2)*L1*cos(phi1)*phi1_tt))*(1.0-xi*xi)/2.0+xi*Rho()*A2*(
-sin(phi1+phi2)*L1*cos(phi1)*phi1_t*phi1_t-sin(phi1+phi2)*L1*sin(phi1)*phi1_tt+
cos(phi1+phi2)*L1*sin(phi1)*phi1_t*phi1_t-cos(phi1+phi2)*L1*cos(phi1)*phi1_tt)*
(1.0-xi));
      Ma = MapleGenVar1+MapleGenVar2;
		}
		//if (!GetMBS()->IsJacobianComputation() && GetMBS()->GetTime() < 0.001) UO() << "Ma=" << Ma << "\n";
		//if (!GetMBS()->IsJacobianComputation() && GetMBS()->GetTime() <= 1) UO() << "phi1=" << phi1 << ",phi2=" << phi2 << "\n";
		return Ma;
	}
}


//externe Kräfte
void ANCFPiezoCable2D::EvalF2(Vector& f, double t)
{
	Body2D::EvalF2(f,t);
	TMStartTimer(22);
	int sos = SOS();
	ConstVector<ANCFCable2D_MaxDOF> fadd(SOS());
	fadd.SetLen(SOS());
	fadd.SetAll(0);

	int ns = NS();
	int dim = Dim();


	ConstVector<ANCFCable2D_MaxDOF> temp(SOS());
	ConstVector<ANCFCable2D_MaxDOF> xg(SOS());
	ConstVector<ANCFCable2D_MaxNShapes> SV(NS());

	if (kappamode != 2) UO() << "Error: PiezoCable2D, kappmode is wrong!\n";

	double EI_w = Em() * GetLz()*Cub(GetLy())/12.;
	double EA = Em() * GetLy() * GetLz(); 
	double A = GetLy()*GetLz();//for Hellinger Reissner
	if (IsBeamParameters())
	{
		EI_w = GetBeamEIy();
		EA = GetBeamEA();
	}

	int order_kappa = 5+4*kappanodal; //regular: 5
	int order_eps   = 9+4*kappanodal; //regular: 9
	ConstVector<ANCFCable2D_MaxIP> xk, wk, xe, we;

	GetIntegrationRule(xk,wk,order_kappa); 
	GetIntegrationRule(xe,we,order_eps);   

	//GetIntegrationRuleLobatto(xk,wk,4);   
	//GetIntegrationRuleLobatto(xe,we,4);   

	//compute kappa0 and eps0
	ConstVector<ANCFCable2D_MaxIP> kappa0, eps0;
	kappa0.SetLen(xk.Length());
	eps0.SetLen(xe.Length());

	int usedmitrochenko = 0;  //for comparison with Euler elastica
	if (epsmode == 2) usedmitrochenko = 1;

	for (int i1=1; i1<=xk.GetLen(); i1++)
	{
		double x = xk(i1);
		kappa0(i1) = GetKappa(x,x_init);
		//UO() << "K=" << kappa0(i1) << "\n";
	}

	for (int i1=1; i1<=xe.GetLen(); i1++)
	{
		double x = xe(i1);
		//eps0(i1) = GetEpsAxial(x,xg_init);

		//1==Green, 2=Biot=Dmitrochenko, 3=linearized strain
		Vector2D rx = GetPosx2D(x*0.5*GetLx(),x_init);
		if (epsmode == 2)
		{
			eps0(i1) = (sqrt(rx*rx)-1.0);
		}
		else if (epsmode == 1)
		{
			eps0(i1) = 0.5*(rx*rx-1.0);
		}
	}

	//SetComputeCoordinates();
	GetCoordinates(xg);

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//kappa*delta kappa:

	for (int i1=1; i1<=xk.GetLen(); i1++)
	{
		double x = xk(i1);

		Vector2D rx0 = GetPosx2D(x*GetLx()*0.5, x_init);
		double rxn = rx0.Norm();
		double det = 0.5*GetLx()*rxn; //is approximately 0.5*GetLx(); in straight elements exactly 0.5*GetLx()!!!

		//bending moments:
		double factkappa = EI_w*wk(i1)*det; //?!!!!Achtung: Determinante der Elementransformation

		if (GetPar(6) == 0)
		{
			double kappa;
			GetDeltaKappa(x,xg,temp, kappa);

			double m0 = EvaluatePiezoMoment(x);

			temp *= (kappa - kappa0(i1) - m0/EI_w)*factkappa;
			fadd += temp;
		}
		else
		{
			//use discrete patches:
			double kappa;
			GetDeltaKappa(x,xg,temp, kappa);

			double m0 = 0;
			if (GetPar(7))
			{
				//estimated feed-forward moment for discrete patches (evaluated only at one position!)
				double mfact = GetPar(10); //total length / piezo actuated length
				//+++++++++++++++++++++++++++++++++++++++++++++++++++++
				//feed-forward control:
				m0 = EvaluatePiezoMoment(GetPar(9)); //amplification because of reduced areas of patches!

				//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
				//feed back control:

				//measure curvature:
				double curv = 0;
				curv += GetMBS()->GetSensor(sensors[8+(int)GetPar(8)*4-3]).GetCurrentValueWithSensorProcessing(mbs->GetTime());
				curv += GetMBS()->GetSensor(sensors[8+(int)GetPar(8)*4-2]).GetCurrentValueWithSensorProcessing(mbs->GetTime());
				curv += GetMBS()->GetSensor(sensors[8+(int)GetPar(8)*4-1]).GetCurrentValueWithSensorProcessing(mbs->GetTime());
				curv += GetMBS()->GetSensor(sensors[8+(int)GetPar(8)*4-0]).GetCurrentValueWithSensorProcessing(mbs->GetTime());
				curv *= 0.25;
				double P = GetPar(11);

				//
				double curvp = 0;
				curvp +=GetMBS()->GetSensor(sensors[8+(int)GetPar(13)*4-3]).GetCurrentValueWithSensorProcessing(mbs->GetTime());
				curvp +=GetMBS()->GetSensor(sensors[8+(int)GetPar(13)*4-2]).GetCurrentValueWithSensorProcessing(mbs->GetTime());
				curvp +=GetMBS()->GetSensor(sensors[8+(int)GetPar(13)*4-1]).GetCurrentValueWithSensorProcessing(mbs->GetTime());
				curvp +=GetMBS()->GetSensor(sensors[8+(int)GetPar(13)*4-0]).GetCurrentValueWithSensorProcessing(mbs->GetTime());
				curvp *= 0.25;
				double D = GetPar(12);

				curv -= (mfact-1.)*m0/EI_w;
				m0 *= mfact; //correction of smaller area of piezo compared to theoretical beam values

				double factp = fabs(m0*200);
				if (fabs(factp) <= 2) factp = 1+0.5*factp;
				//factp = 1;

				m0 += P*EI_w*curv+D/factp*EI_w*curvp; //add a moment which causes opposite curvature than existing curvature
				/*if (fabs(curvp) > 1e-7 || GetMBS()->GetTime() > 0.5) 
				{
					int kk=1;
				}*/

				//+++++++++++++++++++++++++++++++++++++++++++++++++++++
			}

			temp *= (kappa - kappa0(i1) - m0/EI_w)*factkappa;
			fadd += temp;
		}
	}

	//SetComputeCoordinates();
	GetCoordinates(xg);

	//eps*delta eps:
	for (int i1=1; i1<=xe.GetLen(); i1++)
	{
		double x = xe(i1);

		Vector2D rx0 = GetPosx2D(x*GetLx()*0.5, x_init);
		double rxn = rx0.Norm();
		double det = 0.5*GetLx()*rxn; //is approximately 0.5*GetLx(); in straight elements exactly 0.5*GetLx()!!!

		//axial strains:
		double facteps = EA*we(i1)*det;

		int optimized = 1;
		if (optimized)
		{
			GetS0x(SV,x);
			//compute rx:
			Vector2D rx(0.,0.);
			for (int i = 1; i <= Dim(); i++)
			{
				for (int j = 1; j <= NS(); j++)
				{
					rx(i) += SV(j)*xg((j-1)*Dim()+i);
				}
			}

			//compute deltaeps:
			int dim = Dim();
			int ns = NS();

			double fact = 1;
			double eps;
			double kappa; //for formulation based on Green strain!

			//1==Green, 2=Biot=Dmitrochenko, 3=linearized strain
			if (epsmode == 2)
			{
				eps = (sqrt(rx*rx)-1.0) - eps0(i1);
				fact = 1./(sqrt(rx*rx)); //additional part for variation!
			}


			rx *= eps * facteps * fact;
			for (int i=1; i <= dim; i++)
			{
				for (int j=1; j <= ns; j++)
				{
					fadd((j-1)*dim+i) += SV(j)*rx(i);
				}
			}

			//compute epsaxial*deltaeps*fact:
			//temp *= 0.5*(rx*rx-1.0)*facteps;
		}
		else
		{
			GetDeltaEpsAxial(x,xg,temp);
			temp *= GetEpsAxial(x,xg)*facteps;
			fadd += temp;
		}
	}

	if (numConstraintElem)
	{
		fadd(loaddof) += mbs->GetElement(numConstraintElem).XG(1) * gcload;
	}

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
			dmat.SetSize(SOS(),SOS());
			EvalM(dmat,t);
			massmatrix = dmat;
			Mult(dmat,xg,temp);
		}
		//double k=1;
		//if (GetMBS()->GetTime() < 0.4) k=500;
		temp *= GetMassDamping();
		f -= temp;
	}

};

void ANCFCable2D::SetArcLengthParameter(int elem, double load, int ldof)
{
	numConstraintElem = elem;
	gcload = load;
	loaddof = ldof;
}