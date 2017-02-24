//#**************************************************************
//# filename:             BeamShear2D.cpp
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
 

#include "MBS_includes_element.h"
#include "BeamShear2D.h"

BeamShear2D::~BeamShear2D(void) {}

void BeamShear2D::SetBeamShear2D(int bodyind,
		const Vector& xc1, 
		const Vector& xc2, 
		int n1, 
		int n2, 
		int material_num,
		const Vector3D& size,
		const Vector3D& color)
{
	TArray<int>node_list(2); node_list.Add(n1); node_list.Add(n2);		// generate list of node ids
	
	SetGeometricNonlinearityStatus(GNS_NonlinearSmallStrain);		// YV

	this->bodyind = bodyind;
	this->nodes = node_list;
	this->elementname = GetElementSpec();

	SetMaterialNum(material_num);
	col = color;

	this->lx = size.X();
	this->ly = size.Y();
	this->lz = size.Z();
	this->size = size;

	this->geometricNonlinearityStatus = GNS_NonlinearLargeStrain;

	q0 = Vector(6);
	q0(1) = xc1(1); q0(2) = xc1(2); q0(3) = xc2(1); q0(4) = xc2(2);

	x_init = Vector(2*FlexDOF());
	x_init.SetAll(0);
};

void BeamShear2D::LinkToElements()
{
	// dofs are sorted as follows: u,w,u,w,theta,theta
	// then, the getpos, ... functions from generic element can be used
	for (int i = 1; i <= NNodes(); i++)
	{
		const Node& node = this->GetNode(i);

		for (int j = 1; j <= node.SOS()-1; j++)
		{
			this->AddLTG(node.Get(j));
		}
	}
	// angle theta
	this->AddLTG(this->GetNode(1).Get(3));
	this->AddLTG(this->GetNode(2).Get(3));

	for (int i = 1; i <= NNodes(); i++)
	{
		const Node& node = this->GetNode(i);

		for (int j = 1; j <= node.SOS()-1; j++)
		{
			this->AddLTG(node.Get(j+node.SOS()));
		}
	}
	// anglular velocity
	this->AddLTG(this->GetNode(1).Get(3+this->GetNode(1).SOS()));
	this->AddLTG(this->GetNode(2).Get(3+this->GetNode(2).SOS()));
};

Vector3D BeamShear2D::GetDOFDirD(int idof) const
{
	if (idof == 1 || idof == 3) return ToP3D(Vector2D(1.,0.));
	else if (idof == 2 || idof == 4) return ToP3D(Vector2D(0.,1.));
	else if (idof == 5 || idof == 6) return ToP3D(Vector3D(0.,0.,2.));
	return Vector3D(0.,0.,0.);
};

// this could be moved to generic
Vector3D BeamShear2D::GetDOFPosD(int idof) const //returns postion of i-th DOF
{
	if (idof <= 2 || idof == 5)
		return ToP3D(GetPos2DD(Vector2D(-1.,0.)));
	else
		return ToP3D(GetPos2DD(Vector2D(1.,0.)));
};

void BeamShear2D::EvalF2(Vector& f, double t)
{
	Body2D::EvalF2(f, t);
	TMStartTimer(22);

	ConstVector<BeamShear2DmaxDOF> fadd(BeamShear2DmaxDOF);	// temporary storage vector
	ConstVector<BeamShear2DmaxDOF> delta_eps(BeamShear2DmaxDOF), delta_kappa(BeamShear2DmaxDOF), delta_gamma(BeamShear2DmaxDOF);
	
	double EI = GetBeamEIy();
	double EA = GetBeamEA();
	double ksGA = GetBeamGAky();

	for (IntegrationPointsIterator ip(integrationRuleStiffness); !ip.IsEnd(); ++ip)
	{
		double x = ip.Point2D().X();

		GetDeltaEpsAxial(x, delta_eps);
		GetDeltaKappa(x, delta_kappa);
		GetDeltaGamma(x, delta_gamma);

		double kappa = GetKappa(x);
		double eps = GetEpsAxial(x);
		double gamma = GetGamma(x);

		fadd += ((EI*kappa)*delta_kappa + (ksGA*gamma)*delta_gamma + (EA*eps)*delta_eps) * (ip.Weight() * lx*0.5);
	}
	f -= fadd;
};

void BeamShear2D::EvalM(Matrix& m, double t)
{
	if (this->massmatrix.Getrows() != 0)	// mass matrix already computed
	{
		m = massmatrix;
		return;
	}
	massmatrix.SetSize(SOS(), SOS());
	massmatrix.SetAll(0.);

	double rhoA = GetBeamRhoA();
	double rhoIz = GetBeamRhoIz();
	
	for (IntegrationPointsIterator ip(integrationRuleMass); !ip.IsEnd(); ++ip)
	{
		double x = ip.Point2D().X();
		for (int i = 1; i <= NS(); i++)
		{
			for (int j = 1; j <= NS(); j++)
			{
				double d = GetS0(x,i)*GetS0(x,j);
				massmatrix(2*i-1, 2*j-1) += rhoA*d;
				massmatrix(2*i, 2*j) += rhoA*d;
				massmatrix(4+i, 4+j) += rhoIz*d;
			}
		}
		massmatrix *= (0.5*lx*ip.Weight());
	}
	m = massmatrix;
};

Vector2D BeamShear2D::GetDirector1(const double& p_loc) const
{
	double theta = GetTheta(p_loc);
	return Vector2D(cos(theta), sin(theta));
};
Vector2D BeamShear2D::GetDirector1D(const double& p_loc) const
{
	double theta = GetThetaD(p_loc);
	return Vector2D(cos(theta), sin(theta));
};

Vector2D BeamShear2D::GetDirector2(const double& p_loc) const
{
	double theta = GetTheta(p_loc);
	return Vector2D(-sin(theta), cos(theta));
};
Vector2D BeamShear2D::GetDirector2D(const double& p_loc) const
{
	double theta = GetThetaD(p_loc);
	return Vector2D(-sin(theta), cos(theta));
};
Vector2D BeamShear2D::GetRefDirector2(const double& p_loc) const
{
	double theta = GetRefTheta(p_loc);
	return Vector2D(-sin(theta), cos(theta));
};

double BeamShear2D::GetEpsAxial(const double& p_loc) const
{
	Vector2D rx = GetPosx2D(p_loc);
	Vector2D t1 = GetDirector1(p_loc);

	return rx(1)*t1(1) + rx(2)*t1(2) - 1.;
};

double BeamShear2D::GetEpsAxialD(const double& p_loc) const
{
	Vector2D rx = GetPosx2DD(p_loc);
	Vector2D t1 = GetDirector1D(p_loc);

	return rx(1)*t1(1) + rx(2)*t1(2) - 1.;
};

double BeamShear2D::GetGamma(const double& p_loc) const
{
	Vector2D rx = GetPosx2D(p_loc);
	Vector2D t2 = GetDirector2(p_loc);

	return rx(1)*t2(1) + rx(2)*t2(2);
};

double BeamShear2D::GetGammaD(const double& p_loc) const
{
	Vector2D rx = GetPosx2DD(p_loc);
	Vector2D t2 = GetDirector2D(p_loc);

	return rx(1)*t2(1) + rx(2)*t2(2);
};

double BeamShear2D::GetTheta(const double& p_loc) const
{
	double theta = 0.;
	for (int i = 1; i <= NS(); i++)
	{
		double s = GetS0(p_loc, i); 
		theta += XG(4+i)*s;
	}
	return theta;
};
double BeamShear2D::GetThetaD(const double& p_loc) const
{
	double theta = 0.;
	for (int i = 1; i <= NS(); i++)
	{
		double s = GetS0(p_loc, i); 
		theta += XGD(4+i)*s;
	}
	return theta;
};
double BeamShear2D::GetRefTheta(const double& p_loc) const
{
	double theta = 0.;
	for (int i = 1; i <= NS(); i++)
	{
		double s = GetS0(p_loc, i); 
		theta += x_init(4+i)*s;
	}
	return theta;
};

double BeamShear2D::GetKappa(const double& p_loc) const
{
	double kappa = 0.;
	for (int i = 1; i <= NS(); i++)
	{
		double sx = GetS0x(p_loc, i) * 2./lx; 
		kappa += XG(4+i)*sx;
	}
	return kappa;
};

double BeamShear2D::GetKappaD(const double& p_loc) const
{
	double kappa = 0.;
	for (int i = 1; i <= NS(); i++)
	{
		double sx = GetS0x(p_loc, i) * 2./lx;  
		kappa += XGD(4+i)*sx;
	}
	return kappa;
};

void BeamShear2D::GetDeltaEpsAxial(const double& p_loc, Vector& delta_eps)
{
	double theta = GetTheta(p_loc);
	double cos_theta = cos(theta);
	double sin_theta = sin(theta);
	Vector2D rx = GetPosx2D(p_loc);

	for (int i = 1; i <= NS(); i++)
	{
		double sx = GetS0x(p_loc, i) * 2./lx;
		double s  = GetS0(p_loc, i);

		delta_eps(2*i-1) = cos_theta * sx;
		delta_eps(2*i  ) = sin_theta * sx;
		delta_eps(i+4  ) = (-sin_theta*rx(1) + rx(2)*cos_theta) * s;
	}
};

void BeamShear2D::GetDeltaKappa(const double& p_loc, Vector& delta_kappa)
{
	delta_kappa(5) = GetS0x(p_loc, 1) * 2./lx;
	delta_kappa(6) = GetS0x(p_loc, 2) * 2./lx;
};

void BeamShear2D::GetDeltaGamma(const double& p_loc, Vector& delta_gamma)
{
	double theta = GetTheta(p_loc);
	double cos_theta = cos(theta);
	double sin_theta = sin(theta);
	Vector2D rx = GetPosx2D(p_loc);

	for (int i = 1; i <= NS(); i++)
	{
		double sx = GetS0x(p_loc, i) * 2./lx;
		double s  = GetS0(p_loc, i);

		delta_gamma(2*i-1) = -sin_theta * sx;
		delta_gamma(2*i  ) =  cos_theta * sx;
		delta_gamma(i+4  ) = (-cos_theta * rx(1) - rx(2)*sin_theta) * s;
	}
};


#pragma region shapefunction routines
double BeamShear2D::GetS0(const Vector2D& p_loc, int shape) const
{
	double x = p_loc.X();
	switch (shape)
	{
		case 1: return 0.5*(1. - x);
		case 2: return 0.5*(1. + x);	
		default: mbs->UO() << "BeamShear2D::GetS0 Error: Shape function " << shape << " not available!\n"; return 0;
	}
};

double BeamShear2D::GetS0x(const Vector2D& p_loc, int shape) const
{
	double x = p_loc.X();
	switch (shape)
	{
	case 1:	return -0.5;
	case 2:	return  0.5;	
	default: mbs->UO() << "BeamShear2D::GetS0x Error: Shape function " << shape << " not available!\n";	return 0;
	}
};

double BeamShear2D::GetS0xx(const Vector2D& p_loc, int shape) const
{
	return 0.;
};
#pragma endregion

void BeamShear2D::DefineIntegrationRule(IntegrationRule& integrationRule)
{
	assert(integrationRule.settings.elementType == TFE_Beam2D);

	int ruleOrder = 1;

	if (integrationRule.settings.integratedValueType == IntegrationRule::IVT_Mass)
		ruleOrder = 2;
	else
		ruleOrder = 1;

	IntegrationRule::DefineIntegrationRuleLine(integrationRule, ruleOrder);
};