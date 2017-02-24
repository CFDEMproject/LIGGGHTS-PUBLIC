//#**************************************************************
//# filename:             ANCFBeamShearFE2D.cpp
//#
//# author:               Gerstmayr, Vetyukov
//#
//# generated:						
//# description:          
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
 

#include "element.h"
#include "material.h"
#include "Node.h"
#include "body2d.h"
#include "ANCFBeamShearFE2D.h"
#include "femathhelperfunctions.h"




// Set-Routine for linear element
void ANCFBeamShearFE2DLinear::SetANCFBeamShearFE2DLinear(int bodyind, const Vector& xc1, const Vector& xc2, int n1, int n2, 
																						 int material_num,	const Vector3D& size,	const Vector3D& color)
{
	// reuse Set-Routine from FiniteElementGenericBeam2D
	SetFiniteElementGenericBeam2D(bodyind, xc1, xc2, n1, n2, material_num, size, color);

}

void ANCFBeamShearFE2DLinear::SetANCFBeamShearFE2DLinear(int n1, int n2,
																						 int material_num,	const Vector3D& asize,	const Vector3D& color)
{
	this->nodes.SetLen(2);
	this->nodes(1) = n1;
	this->nodes(2) = n2;

	SetGeometricNonlinearityStatus(GNS_NonlinearLargeStrain);		// YV

	this->elementname = GetElementSpec();
	SetMaterialNum(material_num);
	col = color;

	this->size = asize;

	q0.SetLen(2*SOS());
	q0.SetAll(0);
	x_init = Vector(2*SOS());
	int count=0;
	for (int nn=1; nn<=NNodes(); nn++)
	{
		const ANCFNodeS2_2D &node = dynamic_cast<const ANCFNodeS2_2D&> (GetNode(nn));
		q0(++count) = node.GetRefPos2D()(1);
		q0(++count) = node.GetRefPos2D()(2);
		q0(++count) = node.GetRefSlope2()(1);
		q0(++count) = node.GetRefSlope2()(2);
		for (int i=1; i<=node.SOS(); i++)
		{
			x_init((nn-1)*node.SOS()+i) = node.X_Init()(i);
			x_init(SOS()+(nn-1)*node.SOS()+i) = node.X_Init()(node.SOS()+i);
		}
	}

	xgReferenceState.SetLen(SOS());
	xgReferenceState.SetAll(0);
	xgInit.SetXGProvider(this, &xgReferenceState);
	xgCompute.SetXGProvider(this, &xgReferenceState);
	xgDraw.SetXGProvider(this, &xgReferenceState);

}


// Set-Routine for quadratic element
void ANCFBeamShearFE2DQuadratic::SetANCFBeamShearFE2DQuadratic(int bodyind,const Vector& xc1, const Vector& xc2, const Vector& xc3, int n1, int n2, int n3,
																						 int material_num,	const Vector3D& asize,	const Vector3D& color)
{
	TArray<int>node_list(3); 		// generate list of node ids
	node_list(1) = n1;
	node_list(2) = n2;
	node_list(3) = n3;
	this->nodes = node_list;

	SetGeometricNonlinearityStatus(GNS_NonlinearLargeStrain);		// YV

	this->bodyind = bodyind;
	this->elementname = GetElementSpec();
	SetMaterialNum(material_num);
	col = color;

	//this->lx = size.X();
	//this->ly = size.Y();
	//this->lz = size.Z();
	this->size = asize;

	//q0 = (xc1.Append(xc2)).Append(Vector(SOS()));
	q0.SetLen(2*SOS());
	for (int i=1; i<=DOFPerNode(); i++)
	{
		q0(i) = xc1(i);
		q0(i+DOFPerNode()) = xc2(i);
		q0(i+DOFPerNode()*2) = xc3(i);
		q0(SOS()+i) = 0;
		q0(SOS()+i+DOFPerNode()) = 0;
		q0(SOS()+i+DOFPerNode()*2) = 0;
	}
	x_init = Vector(2*FlexDOF());
	x_init.SetAll(0);

	xgReferenceState.SetLen(FlexDOF());
	xgReferenceState.SetAll(0);
	xgInit.SetXGProvider(this, &xgReferenceState);
	xgCompute.SetXGProvider(this, &xgReferenceState);
	xgDraw.SetXGProvider(this, &xgReferenceState);

}

// Set-Routine for quadratic element using data from nodes
void ANCFBeamShearFE2DQuadratic::SetANCFBeamShearFE2DQuadratic(int n1, int n2, int n3,
																						 int material_num,	const Vector3D& asize,	const Vector3D& color)
{
	this->nodes.SetLen(3);
	this->nodes(1) = n1;
	this->nodes(2) = n2;
	this->nodes(3) = n3;

	SetGeometricNonlinearityStatus(GNS_NonlinearLargeStrain);		// YV

	this->elementname = GetElementSpec();
	SetMaterialNum(material_num);
	col = color;

	this->size = asize;

	q0.SetLen(2*SOS());
	q0.SetAll(0);
	x_init = Vector(2*SOS());
	int count=0;
	for (int nn=1; nn<=NNodes(); nn++)
	{
		const ANCFNodeS2_2D &node = dynamic_cast<const ANCFNodeS2_2D&> (GetNode(nn));
		q0(++count) = node.GetRefPos2D()(1);
		q0(++count) = node.GetRefPos2D()(2);
		q0(++count) = node.GetRefSlope2()(1);
		q0(++count) = node.GetRefSlope2()(2);
		for (int i=1; i<=node.SOS(); i++)
		{
			x_init((nn-1)*node.SOS()+i) = node.X_Init()(i);
			x_init(SOS()+(nn-1)*node.SOS()+i) = node.X_Init()(node.SOS()+i);
		}
	}

	xgReferenceState.SetLen(SOS());
	xgReferenceState.SetAll(0);
	xgInit.SetXGProvider(this, &xgReferenceState);
	xgCompute.SetXGProvider(this, &xgReferenceState);
	xgDraw.SetXGProvider(this, &xgReferenceState);

}

void ANCFBeamShearFE2DGeneric::SetMaterialNum(int mnum)
{
	materialnum = mnum;
	if (GetMaterial().IsInelasticMaterial() && GetMaterial().PoissonRatio() != 0)
	{
		mbs->UO() << "Attention: ANCFBeamShearFE2DGeneric: when using inelastic materials, the Poisson ratio has to be zero!!!!!\n\
								 Otherwise the solution is not consistent\n";
	}
}

void ANCFBeamShearFE2DLinear::GetElementData(ElementDataContainer& edc)
{
	Body2D::GetElementData(edc);
}
int ANCFBeamShearFE2DLinear::SetElementData(ElementDataContainer& edc)
{
	int rv = Body2D::SetElementData(edc);
	SetANCFBeamShearFE2DLinear(nodes(1), nodes(2), materialnum, size, col);
	return rv;
}
void ANCFBeamShearFE2DQuadratic::GetElementData(ElementDataContainer& edc)
{
	Body2D::GetElementData(edc);
}
int ANCFBeamShearFE2DQuadratic::SetElementData(ElementDataContainer& edc)
{
	int rv = Body2D::SetElementData(edc);
	SetANCFBeamShearFE2DQuadratic(nodes(1), nodes(2), nodes(3), materialnum, size, col);
	return rv;
}

void ANCFBeamShearFE2DLinear::DefineIntegrationRulesStructuralMech(IntegrationRule& ruleGamma, IntegrationRule& ruleTheta, IntegrationRule& ruleThickness)
{
	int	order_Gamma, order_Theta;
	//for Lobatto Integration:
	int order_Thickness;

	if(!UseReducedIntegration())
	{
		order_Gamma = 4;
		order_Theta = 4;
	}
	else
	{
		order_Gamma = 1;
		order_Theta = 1;
	}
	order_Thickness =2;

	IntegrationRule::DefineIntegrationRuleLine(ruleTheta, order_Theta);
	IntegrationRule::DefineIntegrationRuleLine(ruleGamma, order_Gamma);

	// we use a ready function, which provides a one-dimensional integration rule
	ConstVector<MAX_IP> x, w;
	GetIntegrationRuleLobatto(x, w, order_Thickness);

	ruleThickness.integrationPoints.SetLen(x.Length());
	for (int i = 1; i <= x.GetLen(); i++)
	{
		ruleThickness.integrationPoints(i).x = x(i);
		ruleThickness.integrationPoints(i).y = 0.;
		ruleThickness.integrationPoints(i).z = 0.;
		ruleThickness.integrationPoints(i).weight = w(i);
	}
}

void ANCFBeamShearFE2DQuadratic::DefineIntegrationRulesStructuralMech(IntegrationRule& ruleGamma, IntegrationRule& ruleTheta, IntegrationRule& ruleThickness)
{
	int	order_Gamma, order_Theta;
	//for Lobatto Integration:
	int order_Thickness;

	if(!UseReducedIntegration())
	{
		order_Gamma = 6;
		order_Theta = 6;
	}
	else
	{
		order_Gamma = 3;
		order_Theta = 3;
		}
		order_Thickness = 3;  //this gives order 4 in Lobatto

	IntegrationRule::DefineIntegrationRuleLine(ruleTheta, order_Theta);
	IntegrationRule::DefineIntegrationRuleLine(ruleGamma, order_Gamma);

	// we use a ready function, which provides a one-dimensional integration rule
	ConstVector<MAX_IP> x, w;
	GetIntegrationRuleLobatto(x, w, order_Thickness);

	ruleThickness.integrationPoints.SetLen(x.Length());
	for (int i = 1; i <= x.GetLen(); i++)
	{
		ruleThickness.integrationPoints(i).x = x(i);
		ruleThickness.integrationPoints(i).y = 0.;
		ruleThickness.integrationPoints(i).z = 0.;
		ruleThickness.integrationPoints(i).weight = w(i);
	}
}

// define IntegrationRule for continuum element:
// use_red_integration_x = 1 if reduced integration along axis is required
// use_red_integration_y = 1 if reduced integration in thickness direction to avoid poisson locking
void ANCFBeamShearFE2DGeneric::DefineIntegrationRule(IntegrationRule& integrationRule, int use_red_integration_y)
{
	//assert(integrationRule.settings.elementType == TFE_Beam2D);

	int ruleOrder_x = 0;
	int ruleOrder_y = 0;

	if (integrationRule.settings.integratedValueType == IntegrationRule::IVT_Load)
	{
		ruleOrder_x = NNodes();
		ruleOrder_y = 1;
	}
	else
	{
		if (!UseReducedIntegration())
		{
			ruleOrder_x = 2*NNodes();
		}
		else // reduced integration in x component
		{
			if (NNodes()==2)
			{
				ruleOrder_x = 1;
			}
			else if (NNodes()==3)
			{
				ruleOrder_x = 3;
			}
		}
		if (!use_red_integration_y)
		{
			ruleOrder_y = 2;
		}
		else // reduced integration in y component
		{
			ruleOrder_y = 1;
		}
	}

	if (!GetMaterial().IsInelasticMaterial()) 
	{
		IntegrationRule::DefineIntegrationRuleSquareAnisotropic(integrationRule, ruleOrder_x, ruleOrder_y);
	}
	else
	{
		// integration on plastic grid points ... define integration rule by hand
		Vector xy, wy;
		IntegrationRule rule_x;
		IntegrationRule::DefineIntegrationRuleLine(rule_x, ruleOrder_x);
		//GetIntegrationRule(xy, wy, ruleOrder_y);
		xy.SetLen(plasticip_y);
		wy.SetLen(plasticip_y);
		xy(1) = -1.; xy(plasticip_y) = 1.;
		wy(1) = 1./(plasticip_y-1.);
		wy(plasticip_y) = 1./(plasticip_y-1.);
		for (int i=2; i<plasticip_y; i++)
		{
			xy(i) = -1. + (i-1)*2./(plasticip_y-1.);
			wy(i) = 2./(plasticip_y-1.);
		}

		// and now we create a 2D integration rule
		int n_x = rule_x.integrationPoints.Length();
		integrationRule.integrationPoints.SetLen(n_x * xy.GetLen());
		int cnt = 1;
		for (int i1 = 1; i1 <= n_x; i1++)
		{
			for (int i2 = 1; i2 <= xy.GetLen(); i2++)
			{
				integrationRule.integrationPoints(cnt).x = rule_x.integrationPoints(i1).x;
				integrationRule.integrationPoints(cnt).y = xy(i2);
				integrationRule.integrationPoints(cnt).z = 0;
				integrationRule.integrationPoints(cnt).weight = rule_x.integrationPoints(i1).weight * wy(i2);
				cnt++;
			}
		}
	}
}

#pragma region FE-matrices

// if use_red_integration==0, return the standard elasticity matrix
// if use_red_integration==1, return the locking-free part for poisson-case=1, the poisson-part for poisson-case=2
void ANCFBeamShearFE2DGeneric::GetElasticityMatrix(Matrix3D& Dm, int use_red_integration, int poisson_case)
{
	double ks = ShearCorrectionFactor(); //10.*(1.+Nu())/(12.+11.*Nu()); //for rectangular cross-section only!

	if (use_red_integration==0) // standard elasticity matrix
	{
		if (!GetMaterial().IsPlaneStrain()) // plane stress
		{
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
			Dm(3,3)=0.5*f*(1.-nu2);
		}
		else // plane strain
		{
			Dm.SetAll(0.);
			double la = Em() * Nu() / ((1.+Nu())*(1.-2.*Nu()));
			double mu = Em() / 2. / (1.+Nu());
			Dm(1,1) = 2*mu+la;
			Dm(2,2) = 2*mu+la;
			Dm(1,2) = Dm(2,1) = la;
			Dm(3,3) = mu;
		}
	}
	else // reduced integration for second part of matrix
	{
		Dm.SetAll(0);
		if (poisson_case == 1)//D^0 (equ. (25))
		{
			if (!GetMaterial().IsPlaneStrain()) // plane stress
			{
				//poisson ratio zero:
				Dm(1,1) = Em();
				Dm(2,2) = Em();
				Dm(3,3) = ks*Em()/(2.*(1.+Nu()));//=ks*G
			}
			else // plane strain case
			{
				double mu = Em() / 2. / (1.+Nu());
				Dm(1,1) = Em();
				Dm(2,2) = Em();
				Dm(3,3) = ks*mu;
			}
		}
		else if (poisson_case == 2)//D^v (equ. (26))
		{
			if (!GetMaterial().IsPlaneStrain())
			{
				Dm(1,1) = Em()/(1.-Sqr(Nu()))-Em();//=Em*nu^2/(1-nu^2)
				Dm(1,2) = Nu()*Em()/(1.-Sqr(Nu()));
				Dm(2,1) = Nu()*Em()/(1.-Sqr(Nu()));
				Dm(2,2) = Em()/(1.-Sqr(Nu()))-Em();
			}
			else
			{
				double la = Em() * Nu() / ((1.+Nu())*(1.-2.*Nu()));
				double mu = Em() / 2. / (1.+Nu());
				Dm(1,1) = 2*mu+la-Em();
				Dm(2,2) = 2*mu+la-Em();
				Dm(1,2) = Dm(2,1) = la;
				//Dm(3,3) = mu;
			}
		}
	}

}


void ANCFBeamShearFE2DGeneric::ComputeStrainMatrix2D(Vector2D ploc, Matrix2D& strain, const XGProvider& xg) const
{
	if (UseContinuumMechanicsFormulation())
	{
		Matrix3D jac;
		GetJacobi(jac,ploc,q0);

		Matrix3D jacinv;
		jac.GetInverse(jacinv);
		jacinv = jacinv.GetTp();

		ConstMatrix<MAX_NS*DIM> grad0, grad;// grad0: am Einheitselement, grad: am physischen Element, undeformiert
		grad0.SetSize(Dim(),NS());
		for(int i=1; i<=NS(); i++)
		{
			grad0(1,i)=GetS0x(ploc, i);
			grad0(2,i)=GetS0y(ploc, i);
		}

		Mult(jacinv, grad0, grad);

		Matrix2D F;
		// compute position gradient F = deformation gradient + identity
		F.SetAll(0.);
		for (int i=1; i<=NS(); i++)
		{
			F(1,1) += grad(1,i)*xg.XGdispl(2*i-1);
			F(2,1) += grad(1,i)*xg.XGdispl(2*i);
			F(1,2) += grad(2,i)*xg.XGdispl(2*i-1);
			F(2,2) += grad(2,i)*xg.XGdispl(2*i);
		}
		F(1,1) += 1; 
		F(2,2) += 1;
		////Green-Lagrange strain tensor
		////strain = 0.5 * (F.GetTp() * F - I);
		for (int i=1; i<=2; i++)
		{
			for (int j=1; j<=2; j++)
			{
				strain(i,j) = 0.5*(F(1,i)*F(1,j)+F(2,i)*F(2,j));
			}
			strain(i,i) -= 0.5;
		}
	}
	else
	{
		//double factor = 3./2. * (Sqr(size(2))-Sqr(ploc.Y()*size(2)))/(Sqr(size(2)));
		strain(1,1) = GetGamma12D(ploc.X(),xg)-ploc.Y()*size(2)*0.5*GetThetax2D(ploc.X(),xg);
		strain(2,2) = (GetStretchy(ploc.X(), xg)-1.); //Eyy;
		strain(1,2) = strain(2,1) = 0.5*GetGamma22D(ploc.X(),xg);
	}
}


void ANCFBeamShearFE2DGeneric::EvalF2(Vector& f, double t)
{
	Body2D::EvalF2(f,t);
	TMStartTimer(22);  //CPU timing: 22 = EvalF2 CMS

	if (UseContinuumMechanicsFormulation())
	{
		EvalF2ContMech(f, t);
	}
	else
	{
		EvalF2StructuralMech(f, t);
	}
	TMStopTimer(22);

}

void ANCFBeamShearFE2DGeneric::EvalF2StructuralMech(Vector& f, double t)
{
	int sos = SOS();
	ConstVector<12> fadd;//element residual

	int ns = NS();
	int dim = Dim();

	fadd.SetLen(SOS());
	fadd.SetAll(0);

	ConstVector<12> temp;
	temp.SetLen(SOS());
	temp.SetAll(0);
	ConstVector<12> temp2;
	temp2.SetLen(SOS());
	temp2.SetAll(0);
	ConstVector<6> shape;
	shape.SetLen(NS());

	XGProviderCached<ANCFBeamShearFE2DGenericmaxDOF> xg(xgCompute, XGLength());

	//double EI_w = Em()*GetBeamEIy();   //Flächenträgheitsmoment: moment of inertia of area
	double EI_w = Em() * size(3)*Cub(size(2))/12.;
	double A = GetLy()*GetLz();										 //Querschnittsfläche: cross section area
	double EA = Em()*A; // GetBeamEA();                  //E-Modul*area
	//Shear-Modul: G = E/(2+2nu)
	//double GAks = GetBeamGAky();
	double ks = ShearCorrectionFactor(); //10.*(1.+Nu())/(12.+11.*Nu()); //for rectangular cross-section only!
	double GAks = A*ks*Em()/(2+2*Nu());

	IntegrationRule rule_gamma, rule_theta, rule_thickness;
	DefineIntegrationRulesStructuralMech(rule_gamma, rule_theta, rule_thickness);

	//Theta*DeltaTheta (bending stiffness):
	for (IntegrationPointsIterator ip(&rule_theta); !ip.IsEnd();++ip)  //loop over all integration points
	{
		Vector2D ploc(ip.Point2D());

		Vector2D rx0 = GetPosx2D(ploc, xgInit);
		double det = rx0.Norm(); 

		double fact_EI= EI_w*ip.Weight()*det; //Det from element transformation

		double kappa;// = GetThetax2D(ploc.X(),xg);
		GetDeltaThetax2D(ploc.X(),kappa,temp,xg);  //get DeltaTheta'-Werte and write in temp-vector
		double kappa_plast = GetPlasticKappa(ploc.X());
		temp *= fact_EI*(kappa-kappa_plast); //temp=EI*Theta'*DeltaTheta'
		//mbs->UO()<<"temp="<<temp<<"\n";
		fadd += temp;
	}

	//Gamma1*DeltaGamma1 (axial stiffness):
	//Gamma2*DeltaGamma2 (shear stiffness):
	for (IntegrationPointsIterator ip(&rule_gamma); !ip.IsEnd();++ip)  //loop over all integration points
	{
		Vector2D ploc(ip.Point2D());

		Vector2D rx0 = GetPosx2D(ploc, xgInit);
		double det = rx0.Norm(); 

		double fact_EA= EA*ip.Weight()*det; //Det from element transformation

		double epsax;// = GetGamma12D(ploc.X(),xg);
		double epssh;// = GetGamma22D(ploc.X(), xg);
		GetDeltaGamma1Gamma22D(ploc.X(),epsax,epssh,temp,temp2, xg);  //get DelatGamma1-Werte and write in temp-vector
		
		double epsax_plast = GetPlasticAxialStrain(ploc.X());
		double epssh_plast = 2*GetPlasticShearStrain(ploc.X());

		temp *= fact_EA*(epsax-epsax_plast); //temp=EA*Gamma1*DeltaGamma1
		fadd += temp;


		double fact_GAks= GAks*ip.Weight()*det; //Det from element transformation
		temp2 *= fact_GAks*(epssh-epssh_plast); //temp=ks*GA*Gamma2*DeltaGamma2
		fadd += temp2;//fadd=fadd+temp
	}


	//thickness-Strain:
	//Eyy*DeltaEyy:
	for (IntegrationPointsIterator ip(&rule_thickness); !ip.IsEnd();++ip)  //loop over all integration points
	{
		Vector2D ploc(ip.Point2D());

		Vector2D rx0 = GetPosx2D(ploc, xgInit);
		double det = rx0.Norm(); 

		//double update, Eyy, DeltaEyy;
		double fact_Eyy= EA*ip.Weight()*det;

		double Eyy_plast = GetPlasticThicknessStrain(ploc.X());
		double Eyy; // = GetEyy(ploc.X(), xg);
		// Eyy = stretchy - 1.
		GetDeltaStretchy(ploc.X(), Eyy, temp, xg);
		Eyy -= 1;

		temp *= fact_Eyy*(Eyy-Eyy_plast); //temp=EA*Eyy*DeltaEyy
		fadd += temp;

	}
	
	f -= fadd;
}




void ANCFBeamShearFE2DGeneric::EvalF2ContMech(Vector& f, double t)
{
	int sos = SOS();
	ConstVector<12> fadd;//element residual

	int ns = NS();
	int dim = Dim();

	fadd.SetLen(SOS());
	fadd.SetAll(0);

	ConstVector<12> temp;
	temp.SetLen(SOS());
	temp.SetAll(0);
	ConstVector<6> shape;
	shape.SetLen(NS());

	XGProviderCached<ANCFBeamShearFE2DGenericmaxDOF> xg(xgCompute, XGLength());

	//double EI_w = GetBeamEIy();   //Flächenträgheitsmoment: moment of inertia of area
	double A = GetLy()*GetLz();										 //Querschnittsfläche: cross section area
	//double EA = GetBeamEA();                  //E-Modul*area
	//Shear-Modul: G = E/(2+2nu)
	//double GAks = GetBeamGAky();

	int reduced_integration = UseReducedIntegration();
	int poissoncorrection = UseReducedIntegrationPoisson();

	Matrix3D strain, piola1, F;
	F.SetSize(2,2);

	IntegrationRule integrationRule;
	// kk=1: integrate lockingfree part
	// kk=2: integrate poisson part with reduced order
	for (int kk=1; kk <= 1+poissoncorrection; kk++)
	{
		// get integration rule: for y-direction, no reduced integration for kk=1, reduced integration for kk=2
		DefineIntegrationRule(integrationRule, kk-1);	

		Matrix3D Dm;

		GetElasticityMatrix(Dm,reduced_integration, kk);

		for (IntegrationPointsIterator ip(&integrationRule); !ip.IsEnd(); ++ip)
		{
			//int ind = (i1-1)*kx1+(i2-1);
			Vector2D ploc(ip.Point2D());

			Matrix3D jac;
			GetJacobi(jac,ploc,q0);
			double jacdet = jac.Det();

			Matrix3D jacinv;
			jac.GetInverse(jacinv);
			jacinv = jacinv.GetTp();

			ConstMatrix<MAX_NS*DIM> grad0, grad;// grad0: am Einheitselement, grad: am physischen Element, undeformiert
			grad0.SetSize(Dim(),NS());
			for(int i=1; i<=NS(); i++)
			{
				grad0(1,i)=GetS0x(ploc, i);
				grad0(2,i)=GetS0y(ploc, i);
			}
			Mult(jacinv, grad0, grad);
			// compute position gradient F = deformation gradient + identity
			F.SetAll(0.);
			for (int i=1; i<=NS(); i++)
			{
				F(1,1) += grad(1,i)*xg.XGdispl(2*i-1);
				F(2,1) += grad(1,i)*xg.XGdispl(2*i);
				F(1,2) += grad(2,i)*xg.XGdispl(2*i-1);
				F(2,2) += grad(2,i)*xg.XGdispl(2*i);
			}
			F(1,1) += 1; 
			F(2,2) += 1;

			//Green-Lagrange strain tensor
			//strain = 0.5 * (F.GetTp() * F - I);
			strain = 0.5 * (F.GetTp() * F);
			strain(1,1) -= 0.5; strain(2,2) -= 0.5;	


			Vector3D s3(strain(1,1), strain(2,2), 2.*strain(1,2));
			Vector3D stress;
			stress = Dm * s3;
			// subtract plastic strain
			// if reduced integration for the poisson term is used, the plastic strain
			// is considered only in the first loop, where the non-deviatoric part of the elasticity matrix is integrated
			if (kk==1 && GetMaterial().IsInelasticMaterial())
			{
				ConstVector<MAX_PLASTIC_QUANTITIES> plasticstrain;
				DataToPlasticStrain(plasticstrain, ploc);
				Vector3D trial_s3(s3);
				trial_s3(2) = s3(2) - GetPlasticThicknessStrain(ploc.X()) + plasticstrain(2);
				stress = Dm*trial_s3;

				if (GetMaterial().IsPlaneStrain())
				{
					double mu = Em() / 2. / (1.+Nu());
					stress(1) -= 2*mu*plasticstrain(1);
					stress(2) -= 2*mu*plasticstrain(2);
					stress(3) -= 2*mu*plasticstrain(3);
				}
				else // plane stress
				{
					Vector3D plasticeps(plasticstrain(1),plasticstrain(2),2*plasticstrain(3));
					stress -= Dm*plasticeps;
				}
			}
			piola1.Set22(stress.X(),stress.Z(),stress.Z(),stress.Y()); 

			piola1 = F * piola1;//2.PK

			for (int j=1; j <= dim; j++)
			{
				for (int i = 0; i < ns; i++)
				{
					temp(dim*i+j) = grad(1, i+1)*piola1(j,1) + grad(2, i+1)*piola1(j,2);
				}
			}

			fadd.MultAdd(fabs(jacdet) * GetLz() * ip.Weight(),temp);
		}

	}
	f -= fadd;  

}
void ANCFBeamShearFE2DGeneric::EvalM(Matrix& m, double t)
{
	if (massmatrix.Getcols() == SOS()) // if constant massmatrix is already computed
	{
		m = massmatrix;
		return;
	}
	else // constant massmatrix has to be computed
	{
		int dim = Dim();
		int ns = NS();

		Matrix3D jac;

		IntegrationRule integrationRule;
		IntegrationRule::DefineIntegrationRuleSquareAnisotropic(integrationRule, 6, 2);
		for (IntegrationPointsIterator ip(&integrationRule); !ip.IsEnd(); ++ip)
		{
			Vector2D ploc(ip.Point().X(), ip.Point().Y());
			GetJacobi(jac,ploc,q0);
			double jacdet = jac.Det();
			double fact = fabs (jacdet) * ip.Weight() * GetLz();

			for (int i=0; i<ns; i++)
			{
				double shape_i = GetS0(ploc, i+1);
				for (int j=0; j<ns; j++)
				{
					double shape_j = GetS0(ploc,j+1);
					for (int idim=1; idim<=dim; idim++)
					{
						m(i*dim+idim,j*dim+idim)+=fact*shape_i*shape_j;
					}
				}
			}

		}
		massmatrix = m;
	}  // end computation massmatrix
	return;
}




void ANCFBeamShearFE2DGeneric::GetH(Matrix& H)
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

		IntegrationRule integrationRule;
		IntegrationRule::DefineIntegrationRuleSquareAnisotropic(integrationRule, 3, 1);

		for (IntegrationPointsIterator ip(&integrationRule); !ip.IsEnd(); ++ip)
		{
			Vector2D ploc(ip.Point2D());
			GetJacobi(jac,ploc,q0);
			double jacdet = jac.Det();
			//mbs->UO() << "jacobian = " << jac << "\njacdet = " << jacdet << "\n";
			double fact = fabs (jacdet) * ip.Weight() * GetLz();

			for (int i=0; i<ns; i++)
			{
				for (int j=1; j<=dim; j++)
				{
					H(i*dim+j,j)+=fact*GetS0(ploc,i+1);
				}
			}

		}
		Hmatrix = H;
	}
}

#pragma endregion

void ANCFBeamShearFE2DGeneric::GetdPosdqT(const Vector2D& ploc, Matrix& dpdqi)
{
	dpdqi.SetSize(SOS(),Dim());//12 x 2   (eigentlich:   6*(2x2))
	dpdqi.FillWithZeros();

	for (int i = 1; i <= NS(); i++)
	{
		double s = GetS0(ploc, i);
		dpdqi((i-1)*Dim()+1,1) = s;
		dpdqi((i-1)*Dim()+2,2) = s;
	}
}

// get jacobian [du/dx, du/dy] 
// continuum mechanics formulation: for u corresponding to dof-vector xg0, x, y in reference coordinates, works for pre-curved elements
// reissner formulation: simple scaling ((lx/2, 0), (0, ly/2)), in this formulation pre-curvature is not possible up to now anyway
void ANCFBeamShearFE2DGeneric::GetJacobi(Matrix3D& jac, const Vector2D& ploc, const Vector& xg0) const
{
	jac.SetSize(2,2);
	jac.SetAll(0.);
	if (UseContinuumMechanicsFormulation())
	{
		int ns = NS();
		for (int i = 1; i <= Dim(); i++)
		{ 
			for (int k=1; k <= ns; k++)
			{ 
				jac(i,1) += GetS0x(ploc,k)*xg0((k-1)*Dim()+i);
				jac(i,2) += GetS0y(ploc,k)*xg0((k-1)*Dim()+i);
			}
		}
	}
	else
	{
		jac(1,1) = size(1)*0.5;
		jac(2,2) = size(2)*0.5;
	}
}

#pragma region DOF dir/pos
Vector3D ANCFBeamShearFE2DGeneric::GetDOFDirD(int idof) const
{
	if (idof%4 == 1) return Vector3D(1.,0.,0.);
	else if (idof%4 == 2) return Vector3D(0.,1.,0.);
	else if (idof%4 == 3) return Vector3D(1.,0.,0.);
	else if (idof%4 == 0) return ToP3D(Vector3D(0.,1.,0.));
	return Vector3D(0.,0.,0.);
}

Vector3D ANCFBeamShearFE2DGeneric::GetDOFPosD(int idof) const
{
	if (idof <=4) return GetPosD(Vector3D(-1.,0.,0.));
	else if (idof <= 8) return GetPosD(Vector3D(1.,0.,0.));
	else return Vector3D(0.,0.,0.);
}
#pragma endregion

Vector2D ANCFBeamShearFE2DGeneric::GetRefPos2D(const Vector2D& p_loc) const
{
	Vector2D p(0.,0.);
	for (int i = 1; i <= NS(); i++)
	{
		double s = GetS0(p_loc, i);
		p(1) += s * q0(2*i-1);
		p(2) += s * q0(2*i  );
	}
	return p;
}
Vector2D ANCFBeamShearFE2DGeneric::GetPos2D(const Vector2D& p_loc)	const
{
	return GetPos2D(p_loc, xgCompute);
}
Vector2D ANCFBeamShearFE2DGeneric::GetPos2DD(const Vector2D& p_loc)	const
{
	return GetPos2D(p_loc, xgDraw);
}
Vector2D ANCFBeamShearFE2DGeneric::GetPos2D(const Vector2D& p_loc, const XGProvider& xg)	const
{
	Vector2D p(0.,0.);
	for (int i = 1; i <= NS(); i++)
	{
		double s = GetS0(p_loc, i);
		p(1) += s * (q0(2*i-1)+xg.XGdispl(2*i-1));
		p(2) += s * (q0(2*i  )+xg.XGdispl(2*i  ));
	}
	return p;
}



Vector2D ANCFBeamShearFE2DGeneric::GetPos2DD(const Vector2D& p_loc, double def_scale) const
{
	return GetPos2D(p_loc, def_scale, xgDraw);
}
Vector2D ANCFBeamShearFE2DGeneric::GetPos2D(const Vector2D& p_loc, double def_scale, const XGProvider& xg) const
{
	Vector2D p(0.,0.);
	for (int i = 1; i <= NS(); i++)
	{
		double s = GetS0(p_loc, i);
		p(1) += s * (q0(2*i-1)+def_scale*xg.XGdispl(2*i-1));
		p(2) += s * (q0(2*i  )+def_scale*xg.XGdispl(2*i  ));
	}
	return p;
}


Vector2D ANCFBeamShearFE2DGeneric::GetDisplacement2D(const Vector2D& p_loc) const
{
	return GetDisplacement2D(p_loc, xgCompute);
}
Vector2D ANCFBeamShearFE2DGeneric::GetDisplacement2D(const Vector2D& p_loc, const XGProvider& xg) const
{
	Vector2D p(0.,0.);
	for (int i = 1; i <= NS(); i++)
	{
		double s = GetS0(p_loc, i);
		p(1) += s * xg.XGdispl(2*i-1);
		p(2) += s * xg.XGdispl(2*i  );
	}
	return p;
}
Vector2D ANCFBeamShearFE2DGeneric::GetDisplacement2DD(const Vector2D& p_loc) const
{
	return GetDisplacement2D(p_loc, xgDraw);
}
Vector2D ANCFBeamShearFE2DGeneric::GetVel2D(const Vector2D& p_loc) const
{
	return GetVel2D(p_loc, xgCompute);
}
Vector2D ANCFBeamShearFE2DGeneric::GetVel2DD(const Vector2D& p_loc) const
{
	return GetVel2D(p_loc, xgDraw);
}
Vector2D ANCFBeamShearFE2DGeneric::GetVel2D(const Vector2D& p_loc, const XGProvider& xg) const
{
	Vector2D p(0.,0.);
	for (int i = 1; i <= NS(); i++)
	{
		double s = GetS0(p_loc, i);
		p(1) += s * xg.XGP(2*i-1);
		p(2) += s * xg.XGP(2*i  );
	}
	return p;
}
Vector2D ANCFBeamShearFE2DGeneric::GetRefPosy2D(const Vector2D& p_loc) const
{
	Vector2D p(0.,0.);
	for (int i = 1; i <= NS(); i++)
	{
		double s = GetS0y(p_loc, i);
		p(1) += s * q0(2*i-1);
		p(2) += s * q0(2*i  );
	}
	return p;
}
Vector2D ANCFBeamShearFE2DGeneric::GetPosy2D(const Vector2D& p_loc, const XGProvider& xg) const
{
	Vector2D p(0.,0.);
	for (int i = 1; i <= NS(); i++)
	{
		double s = GetS0y(p_loc, i);
		p(1) += s * (q0(2*i-1)+xg.XGdispl(2*i-1));
		p(2) += s * (q0(2*i)  +xg.XGdispl(2*i  ));
	}
	return p;
}
Vector2D ANCFBeamShearFE2DGeneric::GetPosx2D(const Vector2D& p_loc, const XGProvider& xg) const
{
	Vector2D p(0.,0.);
	for (int i = 1; i <= NS(); i++)
	{
		double s = GetS0x(p_loc, i);
		p(1) += s * (q0(2*i-1)+xg.XGdispl(2*i-1));
		p(2) += s * (q0(2*i)  +xg.XGdispl(2*i  ));
	}
	return p;
}
Vector2D ANCFBeamShearFE2DGeneric::GetVely2D(const Vector2D& p_loc, const XGProvider& xg) const
{
	Vector2D p(0.,0.);
	for (int i = 1; i <= NS(); i++)
	{
		double s = GetS0y(p_loc, i);
		p(1) += s * xg.XGP(2*i-1);
		p(2) += s * xg.XGP(2*i  );
	}
	return p;
}


#pragma region inplane unit vector
Vector2D ANCFBeamShearFE2DGeneric::GetInplaneUnitVector2D(const double& p_loc) const
{
	Vector2D ry = GetPosy2D(Vector2D(p_loc, 0),xgCompute);
	ry.Normalize();
	return ry;
}
Vector2D ANCFBeamShearFE2DGeneric::GetInplaneUnitVector2DD(const double& p_loc) const
{
	Vector2D ry = GetPosy2D(Vector2D(p_loc, 0), xgDraw);
	ry.Normalize();
	return ry;
}
Vector2D ANCFBeamShearFE2DGeneric::GetRefInplaneUnitVector2D(const double& p_loc) const
{
	Vector2D ry0 = GetRefPosy2D(Vector2D(p_loc, 0));
	ry0.Normalize();
	return ry0;
}
Vector2D ANCFBeamShearFE2DGeneric::GetInplaneUnitVectorP2D(const double& p_loc) const
{
	Vector2D ry = GetPosy2D(Vector2D(p_loc, 0), xgCompute);
	Vector2D vy = GetVely2D(Vector2D(p_loc, 0), xgCompute);

	double norm_ry = ry.Norm();
	double norm_ry_p = 1./norm_ry*(ry(1)*vy(1) + ry(2)*vy(2));
	return 1./norm_ry * vy + norm_ry_p/(norm_ry*norm_ry)*ry;
}
#pragma endregion

#pragma region shapefunction routines
double ANCFBeamShearFE2DLinear::GetS0(const Vector2D& ploc, int i) const
{
	switch(i)
	{
	case 1: return 0.5*(1.-ploc(1));
	case 2: return ploc(2)*0.5*(1.-ploc(1));
	case 3: return 0.5*(1.+ploc(1));
	case 4: return ploc(2)*0.5*(1.+ploc(1));
	default: 
		mbs->UO()<<"ANCFBeamShearFE2DLinear::GetS0:Error: Shape function " << i << " requested, " << NS() << " functions available\n"; 
		return 0;
	}
}
double ANCFBeamShearFE2DQuadratic::GetS0(const Vector2D& ploc, int i) const
{
	switch(i)
	{
	case 1: return -0.5*ploc(1)*(1.-ploc(1));
	case 2: return -0.5*ploc(1)*ploc(2)*(1.-ploc(1));
	case 3: return 0.5*ploc(1)*(1.+ploc(1));
	case 4: return 0.5*ploc(1)*ploc(2)*(1.+ploc(1));
	case 5: return -(ploc(1)-1.)*(ploc(1)+1.);
	case 6: return -ploc(2)*(ploc(1)-1.)*(ploc(1)+1.);
	default: 
		mbs->UO()<<"ANCFBeamShearFE2DQuadratic::GetS0:Error: Shape function " << i << " requested, " << NS() << " functions available\n"; 
		return 0;
	}
}

double ANCFBeamShearFE2DLinear::GetS0x(const Vector2D& ploc, int i) const
{
	switch(i)
	{
	case 1: return -0.5;
	case 2: return -ploc(2)*0.5;
	case 3: return 0.5;
	case 4: return ploc(2)*0.5;
	default: 
		mbs->UO()<<"ANCFBeamShearFE2DLinear::GetS0x:Error: Shape function " << i << " requested, " << NS() << " functions available\n"; 
		return 0;
	}
}
double ANCFBeamShearFE2DQuadratic::GetS0x(const Vector2D& ploc, int i) const
{
	switch(i)
	{
	case 1: return -0.5+ploc(1);
	case 2: return ploc(2)*(-0.5+ploc(1));
	case 3: return 0.5+ploc(1);
	case 4: return ploc(2)*(0.5+ploc(1));
	case 5: return -2*ploc(1);
	case 6: return -2*ploc(1)*ploc(2);
	default: 
		mbs->UO()<<"ANCFBeamShearFE2DQuadratic::GetS0x:Error: Shape function " << i << " requested, " << NS() << " functions available\n"; 
		return 0;
	}
}
double ANCFBeamShearFE2DLinear::GetS0xx(const Vector2D& ploc, int i) const
{
	return 0;
}
double ANCFBeamShearFE2DQuadratic::GetS0xx(const Vector2D& ploc, int i) const
{
	switch(i)
	{
	case 1: return 1.;
	case 2: return ploc(2);
	case 3: return 1.;
	case 4: return ploc(2);
	case 5: return -2;
	case 6: return -2*ploc(2);
	default: 
		mbs->UO()<<"ANCFBeamShearFE2DQuadratic::GetS0xx:Error: Shape function " << i << " requested, " << NS() << " functions available\n"; 
		return 0;
	}
}
double ANCFBeamShearFE2DLinear::GetS0y(const Vector2D& ploc, int i) const
{
	switch(i)
	{
	case 1: return 0;
	case 2: return 0.5*(1.-ploc(1));
	case 3: return 0;
	case 4: return 0.5*(1.+ploc(1));
	default: 
		mbs->UO()<<"ANCFBeamShearFE2DLinear::GetS0y:Error: Shape function " << i << " requested, " << NS() << " functions available\n"; 
		return 0;
	}
}
double ANCFBeamShearFE2DQuadratic::GetS0y(const Vector2D& ploc, int i) const
{
	switch(i)
	{
	case 1: return 0;
	case 2: return -0.5*ploc(1)*(1.-ploc(1));
	case 3: return 0;
	case 4: return 0.5*ploc(1)*(1.+ploc(1));
	case 5: return 0;
	case 6: return -(ploc(1)-1.)*(ploc(1)+1.);
	default: 
		mbs->UO()<<"ANCFBeamShearFE2DQuadratic::GetS0y:Error: Shape function " << i << " requested, " << NS() << " functions available\n"; 
		return 0;
	}
}
double ANCFBeamShearFE2DLinear::GetS0xy(const Vector2D& ploc, int i) const
{
	switch(i)
	{
	case 1: return 0;
	case 2: return -0.5;
	case 3: return 0;
	case 4: return 0.5;
	default: 
		mbs->UO()<<"ANCFBeamShearFE2DLinear::GetS0xy:Error: Shape function " << i << " requested, " << NS() << " functions available\n"; 
		return 0;
	}
}
double ANCFBeamShearFE2DQuadratic::GetS0xy(const Vector2D& ploc, int i) const
{
	switch(i)
	{
	case 1: return 0;
	case 2: return (-0.5+ploc(1));
	case 3: return 0;
	case 4: return (0.5+ploc(1));
	case 5: return 0;
	case 6: return -2*ploc(1);
	default: 
		mbs->UO()<<"ANCFBeamShearFE2DQuadratic::GetS0xy:Error: Shape function " << i << " requested, " << NS() << " functions available\n"; 
		return 0;
	}
}

#pragma endregion

const char * str_plastic_curvature_ = "plastic curvature";
const char * str_plastic_ax_strain_ = "plastic axial strain";
const char * str_plastic_shear_strain_ = "plastic shear strain";
const char * str_plastic_thick_strain_ = "plastic thickness strain";
void ANCFBeamShearFE2DGeneric::GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables)
{
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_displacement,
		FieldVariableDescriptor::FVCI_y);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_stress,
		FieldVariableDescriptor::FVCI_y, true, true);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_stress_mises);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_total_strain,
		FieldVariableDescriptor::FVCI_y, true, true);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_axial_extension);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_curvature);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_force_axial);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_moment_bending);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_shear);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_force);
	if(GetMaterial().IsInelasticMaterial())
	{
		FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_inelastic_strain,
			FieldVariableDescriptor::FVCI_y, true, true);
		FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_hardening_parameter);
		FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_yield_function);
		variables.Add(FieldVariableDescriptor(str_plastic_curvature_));
		variables.Add(FieldVariableDescriptor(str_plastic_ax_strain_));
		variables.Add(FieldVariableDescriptor(str_plastic_shear_strain_));
		variables.Add(FieldVariableDescriptor(str_plastic_thick_strain_));
	}
}

double ANCFBeamShearFE2DGeneric::GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector2D& p_loc, bool flagD)
{
	// displacements
	if(fvd.VariableType() == FieldVariableDescriptor::FVT_displacement)
		return fvd.GetComponent(GetDisplacement2DD(p_loc));

	if (fvd.VariableType() == FieldVariableDescriptor::FVT_beam_axial_extension)
		return GetGamma12D(p_loc.X(), xgDraw);

	if (fvd.VariableType() == FieldVariableDescriptor::FVT_beam_force_axial)
		return Em()*size(2)*size(3)*(GetGamma12D(p_loc.X(), xgDraw)-GetPlasticAxialStrain(p_loc.X(),flagD));

	if (fvd.VariableType() == FieldVariableDescriptor::FVT_beam_curvature)
		return GetThetax2D(p_loc.X(), xgDraw);

	if (fvd.VariableType() == FieldVariableDescriptor::FVT_beam_moment_bending)
		return Em()*size(3)*Cub(size(2))/12.*(GetThetax2D(p_loc.X(), xgDraw)-GetPlasticKappa(p_loc.X(),flagD));

	if (fvd.VariableType() == FieldVariableDescriptor::FVT_beam_shear)
		return GetGamma22D(p_loc.X(), xgDraw);

	double ks = ShearCorrectionFactor(); //10.*(1.+Nu())/(12.+11.*Nu()); //for rectangular cross-section only!
	double GAks = size(2)*size(3)*ks*Em()/(2+2*Nu());

	if (fvd.VariableType() == FieldVariableDescriptor::FVT_beam_force)
		return GAks*(GetGamma22D(p_loc.X(), xgDraw)-2*GetPlasticShearStrain(p_loc.X(),flagD));

	if(fvd.VariableType() == FieldVariableDescriptor::FVT_problem_specific)
	{
		if(fvd.GetTextualIdentifierWithoutComponents() == str_plastic_curvature_)		// plastic curvature kappa^p
			return GetPlasticKappa(p_loc.X(), flagD);
		if(fvd.GetTextualIdentifierWithoutComponents() == str_plastic_ax_strain_)		// plastic axial strain eps_ax^p
			return GetPlasticAxialStrain(p_loc.X(), flagD);
		if(fvd.GetTextualIdentifierWithoutComponents() == str_plastic_shear_strain_)		// plastic curvature kappa^p
			return GetPlasticShearStrain(p_loc.X(), flagD);
		if(fvd.GetTextualIdentifierWithoutComponents() == str_plastic_thick_strain_)		// plastic axial strain eps_ax^p
			return GetPlasticThicknessStrain(p_loc.X(), flagD);
	}


	// compute strain and stress
	//Matrix3D Dm;
	//GetElasticityMatrix(Dm,0,1);
	
	int draw_stress_with_constant_interpolation = 1;  // constant interpolation of stress components & mises along axis direction
	Vector2D ploc(p_loc);
	if (flagD && draw_stress_with_constant_interpolation &&
			(
			fvd.VariableType() == FieldVariableDescriptor::FVT_stress_mises ||
			fvd.VariableType() == FieldVariableDescriptor::FVT_stress ||
			fvd.VariableType() == FieldVariableDescriptor::FVT_total_strain
			)
		)
	{
		ploc.X()=0.;
	}

	Matrix2D strain;
	ComputeStrainMatrix2D(ploc, strain, xgDraw);


	Matrix3D Dmat;
	Vector3D strain3(strain(1,1), strain(2,2), 2*strain(1,2));
	GetElasticityMatrix(Dmat, 0, 0);


	// generate deviator, 
	// in case of plane strain, sigma_zz = nu(sigma_xx+sigma_yy), sigma_xx + sigma_yy + sigma_zz = (1+nu)*(sigma_xx+sigma_yy)
	// in case of plane stress, sigma_zz = 0
	double dval;
	if (GetMaterial().IsPlaneStrain()) { dval = (Nu()+1.)/3.;}
	else { dval = 1./3.; }
	Matrix3D deviator(1.);
	deviator(1,1) -= dval;
	deviator(1,2) -= dval;
	deviator(2,1) -= dval;
	deviator(2,2) -= dval;

	Vector3D devstress;
	Vector3D stress;
	// interpolate plastic variables from grid
	ConstVector<MAX_PLASTIC_QUANTITIES> plastic_vars(NPlasticParameters());
	Matrix2D plasticstrainmatrix;
	if (GetMaterial().IsInelasticMaterial())
	{
		DataToPlasticVariablesD(plastic_vars, ploc);
		plasticstrainmatrix(1,1) = plastic_vars(1);
		plasticstrainmatrix(2,2) = plastic_vars(2);
		plasticstrainmatrix(1,2) = plasticstrainmatrix(2,1) = plastic_vars(3);
		double mu = Em() / 2. / (1.+Nu());
			// quadratic interpolation of shear stresses over height
		double factor = 3./2. * (Sqr(size(2))-Sqr(ploc.Y()*size(2)))/(Sqr(size(2)));
		Vector3D trialstrain(strain3);
		trialstrain(2) += - GetPlasticThicknessStrain(ploc.X(), flagD) + plastic_vars(2);
		trialstrain(3) = factor*(strain3(3) - 2*GetPlasticShearStrain(ploc.X(), flagD)) + 2*plastic_vars(3);
		// inel_strain has e_xy in third component, then stress = 2mu*inel_strain
		Vector3D inel_strain(plastic_vars(1), plastic_vars(2), plastic_vars(3));
		if (GetMaterial().IsPlaneStrain())
		{
			devstress = deviator*Dmat*trialstrain - 2.*mu*inel_strain;
		}
		else // plane stress  
		{
			devstress = deviator*Dmat*(trialstrain - inel_strain);
		}
		stress = Dmat*(trialstrain-inel_strain);
	}
	else // elastic material
	{
		stress = Dmat*strain3;
		devstress = deviator*stress;
	}
	// plot the first piola-kirchhoff stresses
	//Matrix3D piola1;
	//piola1.Set22(stress.X(),stress.Z(),stress.Z(),stress.Y()); 
	Matrix3D piola1;
	piola1.Set22(stress(1), stress(3), stress(3), stress(2)); 

	double mises;
	// norm dev(sigma)^2 = devsigma_xx^2 + devsigma_yy^2 + (devsigma_xx+devsigma_yy)^2
	Matrix3D norm_mat(2.);
	norm_mat(1,2) = 1.;
	norm_mat(2,1) = 1.;
	mises = devstress*(norm_mat*devstress);
	mises = sqrt(mises);
	switch(fvd.VariableType())
	{
	case FieldVariableDescriptor::FVT_stress_mises:							return mises;
	case FieldVariableDescriptor::FVT_stress:										return fvd.GetComponent(piola1);
	case FieldVariableDescriptor::FVT_total_strain:							return fvd.GetComponent(strain);
	case FieldVariableDescriptor::FVT_inelastic_strain:		      return fvd.GetComponent(plasticstrainmatrix);
	case FieldVariableDescriptor::FVT_hardening_parameter:      return plastic_vars(NPlasticParameters()-1);
	case FieldVariableDescriptor::FVT_yield_function:           return plastic_vars(NPlasticParameters());
	default: return FIELD_VARIABLE_NO_VALUE;
	}
};

//angle (for curvature): Theta', xbar in (-1,1)
double ANCFBeamShearFE2DGeneric::GetThetax2D(double x0, const XGProvider& xg) const
{

	//Paper-equation (50)
	//Theta'= (r,y Cross r,xy)/Norm^2

	// ry, rxy are derivatives with respect to scaled coordinates
	// x in (-L/2, L/2), y in (-H/2, H/2)
	Vector2D ploc(x0,0.);
	double scalex = 0.5*size(1), scaley = 0.5*size(2);
	Vector2D ry(0), rxy(0);
	int ii=1;
	for (int i=1; i<=NS(); i++)
	{
		double Sy = GetS0y(ploc, i)/scaley;
		double Sxy = GetS0xy(ploc, i)/(scalex*scaley);
		for (int j=1; j<=Dim(); j++)
		{
			ry(j) += Sy*(xg.XGdispl(ii)+q0(ii));
			rxy(j) += Sxy*(xg.XGdispl(ii)+q0(ii));
			ii++;
		}
	}
	double f,g;
	f=ry.Cross(rxy);
	g=ry.Norm2();
	if(fabs(g)<1e-10)
		return 0;  //if denominator gets 0
	else
		return f/g;
}

void ANCFBeamShearFE2DGeneric::GetDeltaThetax2D(double x0, double& kappa, Vector& DeltaThetax, const XGProvider& xg) const
{
	//Paper-equations (54) und (55)

	Vector2D ploc(x0,0.);
	DeltaThetax.SetLen(SOS());
	DeltaThetax.SetAll(0);

	// ry, rxy are derivatives with respect to scaled coordinates
	// x in (-L/2, L/2), y in (-H/2, H/2)
	double scalex = 0.5*size(1), scaley = 0.5*size(2);
	Vector2D ry(0), rxy(0);
	int ii=1;
	ConstVector<MAX_NS> Sy(NS()), Sxy(NS());
	for (int i=1; i<=NS(); i++)
	{
		Sy(i) = GetS0y(ploc, i)/scaley;
		Sxy(i) =  GetS0xy(ploc, i)/(scaley*scalex);
		for (int j=1; j<=Dim(); j++)
		{
			double xgd = xg.XGdispl(ii);
			ry(j) += Sy(i)*(xg.XGdispl(ii)+q0(ii));
			rxy(j) += Sxy(i)*(xg.XGdispl(ii)+q0(ii));
			ii++;
		}
	}

	double f,g;
	f=ry.Cross(rxy);
	g=ry.Norm2();
	if(fabs(g)<1e-10)
	{
		return;  //if denominator gets 0
		kappa = 0;
	}

	kappa = f/g;


	//DeltaThetax = 1/g^2 (g*Deltaf - f*Deltag)
	//Deltaf(k)=(r,y Cross S,xy) - (r,xy Cross S,y)
	//Deltag(k)=2*r,y*S,y

	int k=1;
	for (int j = 1; j <= NS(); j++)//j=1,...,4
	{
		// i==1
		double deltaf_k = -ry(2)*Sxy(j)+rxy(2)*Sy(j);
		double deltag_k = (2.*ry(1)*Sy(j));
		DeltaThetax(k) = g*deltaf_k-f*deltag_k;
		k++;
		// i==2
		deltaf_k = (ry(1)*Sxy(j) -rxy(1)*Sy(j));
		deltag_k = 2.*ry(2)*Sy(j);
		DeltaThetax(k) = g*deltaf_k - f*deltag_k;
		k++;
	}

	DeltaThetax*=1./sqr(g);
}


//Gamma1, Gamma2 (and corresponding delta)
double ANCFBeamShearFE2DGeneric:: GetGamma12D(double x0, const XGProvider& xg) const
{

	Vector2D ploc(x0,0.);
	double scalex = 0.5*size(1);
	// compute r_x and unit vector t2 = r_y / ||r_y||
	Vector2D ty(0), rx(0);
	int ii=1;
	for (int i=1; i<=NS(); i++)
	{
		for (int j=1; j<=Dim(); j++)
		{
			ty(j) += GetS0y(ploc, i)*(xg.XGdispl(ii)+q0(ii));
			rx(j) += GetS0x(ploc, i)*(xg.XGdispl(ii)+q0(ii))/scalex;
			ii++;
		}
	}
	ty.Normalize();
	// inner product t1*rx, t1 = rotated t2
	double t1rx = ty.Y()*rx.X() - ty.X()*rx.Y();
	return t1rx - 1.;
}

double ANCFBeamShearFE2DGeneric:: GetGamma22D(double x0, const XGProvider& xg) const
{

	Vector2D ploc(x0,0.);
	double scalex = 0.5*size(1);
	// compute r_x and unit vector t2 = r_y / ||r_y||
	Vector2D ty(0), rx(0);
	int ii=1;
	for (int i=1; i<=NS(); i++)
	{
		double sy = GetS0y(ploc, i);
		double sx = GetS0x(ploc, i);
		for (int j=1; j<=Dim(); j++)
		{
			ty(j) += sy*(xg.XGdispl(ii)+q0(ii));
			rx(j) += sx*(xg.XGdispl(ii)+q0(ii))/scalex;
			ii++;
		}
	}
	ty.Normalize();
	// inner product t2*rx
	return ty.X()*rx.X() + ty.Y()*rx.Y();
}


void ANCFBeamShearFE2DGeneric:: GetDeltaGamma1Gamma22D(double x0, double& Gamma1, double& Gamma2, Vector& DeltaGamma1, Vector& DeltaGamma2, const XGProvider& xg) const
{
	//Paper-equation (57)
	Vector2D ploc(x0,0.);

	DeltaGamma1.SetLen(SOS());
	DeltaGamma1.SetAll(0);

	double scalex = 0.5*size(1), scaley = 0.5*size(2);
	// compute r_x and unit vector t2 = r_y / ||r_y||
	Vector2D ry(0), rx(0);
	// rotated vectors
	Vector2D rxhat, ryhat;

	// Vectors storing shape function values, already scaled to scaled element gradients
	ConstVector<MAX_NS> Sy(NS()), Sx(NS());

	int ii=1;
	for (int i=1; i<=NS(); i++)
	{
		Sy(i) = GetS0y(ploc, i)/scaley;
		Sx(i) = GetS0x(ploc, i)/scalex;
		for (int j=1; j<=Dim(); j++)
		{
			ry(j) += Sy(i)*(xg.XGdispl(ii)+q0(ii));
			rx(j) += Sx(i)*(xg.XGdispl(ii)+q0(ii));
			ii++;
		}
	}
	rxhat(1)=rx(2);
	rxhat(2)=-rx(1);
	ryhat(1)=ry(2);
	ryhat(2)=-ry(1);

	double ryNorm, ry3;
	ryNorm=ry.Norm();
	ry3=ryNorm*ryNorm*ryNorm;

	// compute Gamma1
	double t1rx = (ry.Y()*rx.X() - ry.X()*rx.Y())/ryNorm;
	Gamma1 = t1rx - 1.;

	// compute Gamma2
	Gamma2 = (rx.X()*ry.X()+rx.Y()*ry.Y())/ryNorm;


	if(ryNorm<1e-10)
		return;  //if denominator gets 0

	int k=1;
	for (int j = 1; j <= NS(); j++)//j=1,...,4
	{
		for (int i = 1; i <= Dim(); i++)  //i=1,2
		{
			//int k=(j-1)*Dim()+i;
			DeltaGamma1(k)=ryhat(i)*Sx(j)/(ryNorm) - rxhat(i)*Sy(j)/(ryNorm) + (rxhat*ry)*ry(i)*Sy(j)/(ry3);
			DeltaGamma2(k)=ry(i)*Sx(j)/(ryNorm) + rx(i)*Sy(j)/(ryNorm) - (rx*ry)*ry(i)*Sy(j)/(ry3);
			k++;
		}
	}
}


void ANCFBeamShearFE2DGeneric:: GetDeltaEyy(double x0, double& Eyy, Vector& DeltaEyy, const XGProvider& xg) const
{
	Vector2D ploc(x0,0.);

	DeltaEyy.SetLen(SOS());
	DeltaEyy.SetAll(0);

	double scaley = 0.5*size(2);
	// compute r_y
	Vector2D ry(0);
	int ii=1;
	// Vector storing shape function values, already scaled to scaled element gradients
	ConstVector<MAX_NS> Sy(NS());
	for (int i=1; i<=NS(); i++)
	{
		Sy(i) = GetS0y(ploc, i)/scaley;
		for (int j=1; j<=Dim(); j++)
		{
			ry(j) += Sy(i)*(xg.XGdispl(ii)+q0(ii));
			ii++;
		}
	}

	double ryNorm = ry.Norm();

	Eyy = 0.5*(ry.Norm2()-1);


	if(ryNorm<1e-10)
		return;  //if denominator gets 0

	int k=1;
	for (int j = 1; j <= NS(); j++)//j=1,...,4
	{
		for (int i = 1; i <= Dim(); i++)  //i=1,2
		{
			//int k=(j-1)*Dim()+i;

			DeltaEyy(k)=ry(i)*Sy(j);

			k++;
		}
	}

}

double ANCFBeamShearFE2DGeneric:: GetEyy(double x0, const XGProvider& xg) const
{
	Vector2D ploc(x0,0.);
	double scaley = 0.5*size(2);
	// compute r_y
	Vector2D ry(0);
	int ii=1;
	for (int i=1; i<=NS(); i++)
	{
		for (int j=1; j<=Dim(); j++)
		{
			ry(j) += GetS0y(ploc, i)*(xg.XGdispl(ii)+q0(ii))/scaley;
			ii++;
		}
	}

	return 0.5*(ry.Norm2()-1);

}

double ANCFBeamShearFE2DGeneric:: GetStretchy(double x0, const XGProvider& xg) const
{
	Vector2D ploc(x0,0.);
	double scaley = 0.5*size(2);
	// compute r_y
	Vector2D ry(0);
	int ii=1;
	for (int i=1; i<=NS(); i++)
	{
		for (int j=1; j<=Dim(); j++)
		{
			ry(j) += GetS0y(ploc, i)*(xg.XGdispl(ii)+q0(ii))/scaley;
			ii++;
		}
	}

	return (ry.Norm());

}

void ANCFBeamShearFE2DGeneric:: GetDeltaStretchy(double x0, double& stretchy, Vector& DeltaStretchy, const XGProvider& xg) const
{
	Vector2D ploc(x0,0.);

	DeltaStretchy.SetLen(SOS());
	DeltaStretchy.SetAll(0);

	double scaley = 0.5*size(2);
	// compute r_y
	Vector2D ry(0);
	int ii=1;
	// Vector storing shape function values, already scaled to scaled element gradients
	ConstVector<MAX_NS> Sy(NS());
	for (int i=1; i<=NS(); i++)
	{
		Sy(i) = GetS0y(ploc, i)/scaley;
		for (int j=1; j<=Dim(); j++)
		{
			ry(j) += Sy(i)*(xg.XGdispl(ii)+q0(ii));
			ii++;
		}
	}

	double ryNorm = ry.Norm();

	stretchy = ryNorm;


	if(ryNorm<1e-10)
		return;  //if denominator gets 0

	int k=1;
	for (int j = 1; j <= NS(); j++)//j=1,...,4
	{
		for (int i = 1; i <= Dim(); i++)  //i=1,2
		{
			//int k=(j-1)*Dim()+i;
			DeltaStretchy(k)= ry(i)*Sy(j)/ryNorm;
			k++;
		}
	}

}

double ANCFBeamShearFE2DGeneric::PostNewtonStep(double t)
{
	// no PostNewtonStep for elastic materials
	if (!GetMaterial().IsInelasticMaterial())
	{
		return 0;
	}
	// Plasticity... FixedPointIteration procedure

	//// in continuum mechanics formulation, PostNewtonStep is inherited from FiniteElementGenericBeam2D
	//// in Structural mechanics formulation, a faster algorithm requiring fewer computations of the strain tensor is given below
	////   it should give exactly the same results as the inherited PostNewtonStep
	//if (UseContinuumMechanicsFormulation())
	//{
	//	return FiniteElementGenericBeam2D::PostNewtonStep(t);
	//}


	double error = 0;
	// loop over all grid points
	double hx = 2./(plasticip_x-1.);
	double hy = 2./(plasticip_y-1.);
	XGProviderCached<FE2DmaxDOF> xg(xgCompute, SOS());

	Matrix3D D, Dinv;
	GetElasticityMatrix(D);
	D.GetInverse(Dinv);

	double gamma1, gamma2, Eyy, thetax;
	for (int j=1; j<=plasticip_x; j++)
	{
		// grid points are ordered as in matrix
		double x0 = -1.+(j-1)*hx;
		// precompute strain resultants at position x0 for structural mechanics formulation
		if (!UseContinuumMechanicsFormulation())
		{
			gamma1 = GetGamma12D(x0, xg);
			gamma2 = GetGamma22D(x0, xg);
			Eyy = GetStretchy(x0, xg)-1.;
			thetax = GetThetax2D(x0, xg);
		}
		// precompute plastic strain resultants at position x0
		double Eyy_plast = GetPlasticThicknessStrain(x0);
		double Exy_plast = GetPlasticShearStrain(x0);

		for (int i=1; i<=plasticip_y; i++)
		{
			double y0 = 1.-(i-1)*hy;
			Vector2D ploc(x0, y0);
			// quadratic interpolation of shear stresses over height
			double factor = 3./2. * (Sqr(size(2))-Sqr(ploc.Y()*size(2)))/(Sqr(size(2)));
			//double factor = 1;
			// strain matrix
			Matrix2D strain;
			//structural mechanics: faster definition of strain tensor by pre-computed gamma1, gamma2 and thetax
			if (!UseContinuumMechanicsFormulation())
			{
				strain(1,1) = gamma1-ploc.Y()*size(2)*0.5*thetax;
				strain(2,2) = Eyy; //Eyy;
				strain(1,2) = strain(2,1) = 0.5*gamma2;
			}
			else	// continuum mechanics: compute strain tensor
			{
				ComputeStrainMatrix2D(ploc, strain, xg);
			}
			// write plasticVariables from XData in struct;
			// values used as starting point for this iteration - may be from last converged time step or from last iteration
			// see 10 lines below
			ConstVector<MAX_PLASTIC_QUANTITIES> plasticvars(NPlasticParameters());
			// actual values to be stored during PostNewtonStep
			ConstVector<MAX_PLASTIC_QUANTITIES> plasticvars_store(NPlasticParameters());
			DataToPlasticVariables(plasticvars_store, j, i);
			// values from last converged time step
			ConstVector<MAX_PLASTIC_QUANTITIES> plasticvars_last(NPlasticParameters());
			DataToPlasticVariablesLastTimeStep(plasticvars_last, j, i);

			// choose here which actual values to use
			plasticvars = plasticvars_last;
			//plasticvars = plasticvars_store;

			// compute trial strain which mimicks constant stress (yy) over height and quadratic distribution of shear stress sigma_xy
			Matrix2D trial_strain(strain);
			trial_strain(2,2) = strain(2,2) - Eyy_plast + plasticvars(2);
			trial_strain(1,2) = factor*(strain(1,2) - Exy_plast) + plasticvars(3);
			trial_strain(2,1) = trial_strain(1,2);

			if (GetMaterial().IsPlaneStrain())
			{
				GetMaterial().ComputeInelasticVariablesFromStrainPlaneStrain2D(plasticvars, trial_strain);
			}
			else // plane stress
			{
				GetMaterial().ComputeInelasticVariablesFromStrainPlaneStress2D(plasticvars, trial_strain);
			}

			// add to error the difference between the stored plasticvars_store and the new plasticvars (only strain components
			// scale by length of element and by number of integration points
			double delta_epsplast_laststep = 0;
			double delta_epsplast_thisstep = 0;
			for (int ii=1; ii<=NPlasticParameters()-2; ii++) 
			{
				delta_epsplast_laststep += fabs(plasticvars(ii)-plasticvars_last(ii)); 
				delta_epsplast_thisstep += fabs(plasticvars(ii)-plasticvars_store(ii)); 
			}
			// if delta_epsplast > 0, then yield function has to be zero (compatibility condition)
			if (plasticvars(NPlasticParameters()) < - GetMBS()->DiscontinuousAccuracy() && delta_epsplast_laststep > GetMBS()->DiscontinuousAccuracy()/Em())
			{
				error += size(1)*delta_epsplast_laststep	/(plasticip_x*plasticip_y);
				plasticvars = plasticvars_last;
			}
			else // otherwise, add error quantity and use plastic variables returned in vector
			{
				error += size(1)*delta_epsplast_thisstep	/(plasticip_x*plasticip_y);
			}

			// map vector back to plasticvars and then to XData
			PlasticVariablesToData(plasticvars, j, i);


		}
	}

	return error;
}