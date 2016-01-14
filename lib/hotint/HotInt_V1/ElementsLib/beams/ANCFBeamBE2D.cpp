//#**************************************************************
//# filename:             ANCFBeamBE2D.cpp
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
#include "ANCFBeamBE2D.h"
#include "femathhelperfunctions.h"

ANCFBeamBE2D::~ANCFBeamBE2D(void) {};

void ANCFBeamBE2D::SetANCFBeamBE2D(int bodyind, const Vector& xc1, const Vector& xc2, int n1, int n2, int material_num, const Vector3D& size, const Vector3D& color)
{
	SetFiniteElementGenericBeam2D(bodyind, xc1, xc2, n1, n2, material_num, size, color);
};

// get direction of of action of i-th DOF
// AH: this is for drawing cones that represent the coordinate constraints
Vector3D ANCFBeamBE2D::GetDOFDirD(int idof) const
{
	if (idof > 6) 
		idof -= 6;

	if (idof == 1) return ToP3D(Vector2D(1.,0.));
	else if (idof == 2) return ToP3D(Vector2D(0.,1.));
	else if (idof == 3) return Vector3D(0.,0.,0.);
	else if (idof == 4) return ToP3D(Vector3D(0.,0.,2.));
	return Vector3D(0.,0.,0.);
};

Vector3D ANCFBeamBE2D::GetDOFPosD(int idof) const //returns postion of i-th DOF
{
	if (idof <= 4)
		return ToP3D(GetPos2DD(Vector2D(-1.,0.)));
	else
		return ToP3D(GetPos2DD(Vector2D(1.,0.)));
};

void ANCFBeamBE2D::EvalF2(Vector& f, double t)
{
	Body2D::EvalF2(f, t);
	TMStartTimer(22);

	// do the computation
	/*ConstVector<ANCFBeamBE2DmaxDOF> fadd(2*NS());	// temporary storage vector
	ConstVector<ANCFBeamBE2DmaxDOF> delta_eps(2*NS()), delta_kappa(2*NS());
	
	double EI = GetBeamEIy();
	double EA = GetBeamEA();
	for (IntegrationPointsIterator ip(integrationRuleStiffness); !ip.IsEnd(); ++ip)
	{
		double x = ip.Point2D().X();
		double eps = GetEpsAxial(x);
		double kappa = GetKappa(x);
		GetDeltaEpsAxial(x, delta_eps); 
		GetDeltaKappa(x, delta_kappa);

		fadd += ((EA*eps)*delta_eps + (EI*kappa)*delta_kappa) * (ip.Weight() * lx*0.5);
	}	
	f -= fadd;*/

	ConstVector<ANCFBeamBE2DmaxDOF> fadd(2*NS());	// temporary storage vector
	ConstVector<ANCFBeamBE2DmaxDOF> delta_eps(2*NS()), delta_kappa(2*NS());
	
	double EI = GetBeamEIy();
	double EA = GetBeamEA();

	for (IntegrationPointsIterator ip(integrationRuleStiffness); !ip.IsEnd(); ++ip)
	{
		double x = ip.Point2D().X();

		//Vector2D rx = GetPosx2D(x);
		//Vector2D rxx = GetPosxx2D(x);

		Vector2D rx, rxx;
		for (int i = 1; i <= NS(); i++)
		{
			//double sx = ip.GetVectorData(1)->Get(i) * 2./lx;
			//double sxx = ip.GetVectorData(2)->Get(i) * 4./(lx*lx);
			double sx = GetS0x(x, i) * 2./GetLx();
			double sxx = GetS0xx(x, i) * 4./(GetLx()*GetLx());
			double d1 = q0(2*i-1) + XG(2*i-1);
			double d2 = q0(2*i  ) + XG(2*i  );

			rx(1) += sx * d1;
			rx(2) += sx * d2;
			rxx(1) += sxx * d1;
			rxx(2) += sxx * d2;
		}

		double rxn2 = rx.Norm2();

		double eps = rx.Norm() - 1.0;
		double kappa = (rx.Cross(rxx))/rxn2;

		double d = 1./rx.Norm();

		double kdux = -2.*rx(1)*kappa / rxn2 + rxx(2) / rxn2;
		double kduxx = -rx(2) / rxn2;
		double kdwx = -2.*rx(2)*kappa / rxn2 - rxx(1) / rxn2;
		double kdwxx = rx(1) / rxn2;

		for (int j = 1; j <= NS(); j++)
		{
			double sx = GetS0x(x, j) * 2./GetLx();
			double sxx = GetS0xx(x, j) * 4./(GetLx()*GetLx());
			//double sx = ip.GetVectorData(1)->Get(j) * 2./lx;
			//double sxx = ip.GetVectorData(2)->Get(j) * 4./(lx*lx);

			delta_eps(2*j-1) = d * rx(1) * sx;
			delta_eps(2*j  ) = d * rx(2) * sx;

			delta_kappa(2*j-1) = (kdux * sx + kduxx * sxx);
			delta_kappa(2*j  ) = (kdwx * sx + kdwxx * sxx);
		}

		fadd += ((EA*eps)*delta_eps + (EI*kappa)*delta_kappa) * (ip.Weight() * GetLx()*0.5);
	}	
	f -= fadd;

	TMStopTimer(22);
};

void ANCFBeamBE2D::EvalM(Matrix& m, double t)
{
	if (this->massmatrix.Getrows() != 0)	// mass matrix already computed
	{
		m = massmatrix;
		return;
	}
	
	Matrix H(Dim(), SOS());
	massmatrix.SetSize(SOS(), SOS());
	massmatrix.SetAll(0.);

	double rhoA = GetBeamRhoA();
	
	for (IntegrationPointsIterator ip(integrationRuleMass); !ip.IsEnd(); ++ip)
	{
		for (int i = 1; i <= NS(); i++)
		{
			H(1, 2*i-1) = GetS0(ip.Point2D().X(), i);
			H(2, 2*i) = H(1, 2*i-1);
		}
		H = (0.5*GetLx()*rhoA*ip.Weight()) * (H.GetTp() * H);
		massmatrix += H;
	}
	m = massmatrix;
};

void ANCFBeamBE2D::GetH(Matrix& H)
{
	if (Hmatrix.Getrows() == SOS())
	{
		H = Hmatrix;
		return;
	}
	else
	{
		double A = this->GetMaterial().BeamRhoA() / this->GetMaterial().Density();

		H.SetSize(SOS(), Dim());
		H.SetAll(0);

		for (IntegrationPointsIterator ip(integrationRuleLoad); !ip.IsEnd(); ++ip)
		{
			double x = ip.Point2D().X();

			// jacobi determinant
			//Vector2D rx0 = GetRefPosx2D(x);
			double det = 0.5*GetLx(); //*rx0.Norm();

			double d = A * det * ip.Weight();

			for (int i = 1; i <= NS(); i++)
			{
				double s = GetS0(x, i);
				H(2*i-1, 1) += d * s;
				H(2*i,   2) += d * s; //H(2*i-1, 1);
			}	
		}
		Hmatrix = H;
	}
};

void ANCFBeamBE2D::GetdPosdqT(const Vector2D& p_loc, Matrix& dpdqi)
{
	//p = S(p.x,p.y,p.z)*q; d(p)/dq
	dpdqi.SetSize(SOS(),Dim());
	dpdqi.FillWithZeros();
	//d = S + ...
	for (int i = 1; i <= NS(); i++)
	{
		double s = GetS0(p_loc.X(), i);
		dpdqi((i-1)*Dim()+1,1) = s;
		dpdqi((i-1)*Dim()+2,2) = s;
	}
	if (p_loc.Y() != 0)
	{
		double y = p_loc.Y();
		Vector2D rx = GetPosx2D(p_loc.X());
		Vector2D n(-rx.X(), rx.Y());
		n /= rx.Norm();

		for (int i = 1; i <= NS(); i++)
		{
			double sx = GetS0x(p_loc.X(), i) * 2./GetLx();
			//y/|n|*dn/dq
			dpdqi((i-1)*Dim()+1,2) +=  y*sx;
			dpdqi((i-1)*Dim()+2,1) += -y*sx;

			//y*n/|n|*(r_x1*S_x1 + r_x2*S_x2)
			dpdqi((i-1)*Dim()+1,1) +=  y*n.X()*(rx.X()*sx);
			dpdqi((i-1)*Dim()+1,2) +=  y*n.Y()*(rx.X()*sx);
			dpdqi((i-1)*Dim()+2,1) +=  y*n.X()*(rx.Y()*sx);
			dpdqi((i-1)*Dim()+2,2) +=  y*n.Y()*(rx.Y()*sx);
		}
	}
};

Vector2D ANCFBeamBE2D::GetInplaneUnitVector2D(const double& p_loc) const
{
	Vector2D rx = GetPosx2D(p_loc);
	rx.Normalize();
	return Vector2D(-rx(2), rx(1));
};

Vector2D ANCFBeamBE2D::GetInplaneUnitVector2DD(const double& p_loc) const
{
	Vector2D rx = GetPosx2DD(p_loc);
	rx.Normalize();
	return Vector2D(-rx(2), rx(1));
};

Vector2D ANCFBeamBE2D::GetInplaneUnitVectorP2D(const double& p_loc) const
{
	Vector2D rx = GetPosx2D(p_loc);
	Vector2D vx = GetVelx2D(p_loc);
	
	double rxn = rx.Norm();
	double rxnp = 1./rxn*(rx(1)*vx(1) + rx(2)*vx(2));
	return 1./rxn * (Vector2D(-vx(2),vx(1)) - rxnp/rxn * Vector2D(-rx(2),rx(1)));
};

Vector2D ANCFBeamBE2D::GetRefInplaneUnitVector2D(const double& p_loc) const
{
	Vector2D rx0 = GetRefPosx2D(p_loc);
	rx0.Normalize();
	return Vector2D(-rx0(2), rx0(1));
};

double ANCFBeamBE2D::GetEpsAxial(const double& p_loc) const
{
	Vector2D rx = GetPosx2D(p_loc);
	return rx.Norm() - 1.0;
};

double ANCFBeamBE2D::GetEpsAxialD(const double& p_loc) const
{
	Vector2D rx = GetPosx2DD(p_loc);
	return rx.Norm() - 1.0;
};

double ANCFBeamBE2D::GetKappa(const double& p_loc) const
{
	Vector2D rx = GetPosx2D(p_loc);
	Vector2D rxx = GetPosxx2D(p_loc);
	return (rx.Cross(rxx))/rx.Norm2();
};

double ANCFBeamBE2D::GetKappaD(const double& p_loc) const
{
	Vector2D rx = GetPosx2DD(p_loc);
	Vector2D rxx = GetPosxx2DD(p_loc);
	return (rx.Cross(rxx))/rx.Norm2();
}

void ANCFBeamBE2D::GetDeltaEpsAxial(const double& p_loc, Vector& delta_eps)
{
	Vector2D rx = GetPosx2D(p_loc);
	double d = 1./rx.Norm();
	for (int j = 1; j <= NS(); j++)
	{
		double sx = GetS0x(p_loc, j) * 2./GetLx();
		delta_eps(2*j-1) = d * rx(1) * sx;
		delta_eps(2*j  ) = d * rx(2) * sx;
	}
};

void ANCFBeamBE2D::GetDeltaKappa(const double& p_loc, Vector& delta_kappa)
{
	Vector2D rx = GetPosx2D(p_loc);
	Vector2D rxx = GetPosxx2D(p_loc);
	double rxn2 = rx.Norm2();
	double kappa = (rx.Cross(rxx))/rxn2;

	//AH: compute deltakappa
	double kdux = -2.*rx(1)*kappa / (rxn2*rxn2) + rxx(2) / rxn2;
	double kduxx = -rx(2) / rxn2;
	double kdwx = -2.*rx(2)*kappa / (rxn2*rxn2) - rxx(1) / rxn2;
	double kdwxx = rx(1) / rxn2;

	for (int j = 1; j <= NS(); j++)
	{
		double sx = GetS0x(p_loc, j) * 2./GetLx();
		double sxx = GetS0xx(p_loc, j) * 4./(GetLx()*GetLx());
		delta_kappa(2*j-1) = (kdux * sx + kduxx * sxx);
		delta_kappa(2*j  ) = (kdwx * sx + kdwxx * sxx);
	}
};

#pragma region shapefunction routines
double ANCFBeamBE2D::GetS0(const Vector2D& p_loc, int shape) const
{
	double x = p_loc.X();
	switch (shape)
	{
		case 1: return 0.125*( 4.	- 7.5*x				+ 5.*x*x*x					- 1.5*x*x*x*x*x);
		case 2: return 0.125*( 2.5	- 3.5*x - 3.*x*x	+ 5.*x*x*x	+ 0.5*x*x*x*x	- 1.5*x*x*x*x*x);		
		case 3: return 0.125*( 0.5	- 0.5*x	- 1.*x*x	+ 1.*x*x*x	+ 0.5*x*x*x*x	- 0.5*x*x*x*x*x);		
		case 4: return 0.125*( 4.	+ 7.5*x				- 5.*x*x*x					+ 1.5*x*x*x*x*x);				
		case 5: return 0.125*(-2.5	- 3.5*x	+ 3.*x*x	+ 5.*x*x*x	- 0.5*x*x*x*x	- 1.5*x*x*x*x*x);	
		case 6: return 0.125*( 0.5	+ 0.5*x	- 1.*x*x	- 1.*x*x*x	+ 0.5*x*x*x*x	+ 0.5*x*x*x*x*x);
		default: mbs->UO() << "ANCFBeamBE2D::GetS0 Error: Shape function " << shape << " not available!\n"; return 0;
	}
};

double ANCFBeamBE2D::GetS0x(const Vector2D& p_loc, int shape) const
{
	double x = p_loc.X();
	switch (shape)
	{
	case 1:	return 0.125*(-7.5			+ 15.*x*x				- 7.5*x*x*x*x);
	case 2:	return 0.125*(-3.5	- 6.*x	+ 15.*x*x	+ 2.*x*x*x	- 7.5*x*x*x*x);	
	case 3:	return 0.125*(-0.5	- 2.*x	+  3.*x*x	+ 2.*x*x*x	- 2.5*x*x*x*x);
	case 4:	return 0.125*( 7.5			- 15.*x*x				+ 7.5*x*x*x*x);
	case 5:	return 0.125*(-3.5	+ 6.*x	+ 15.*x*x	- 2.*x*x*x	- 7.5*x*x*x*x);
	case 6:	return 0.125*( 0.5	- 2.*x	-  3.*x*x	+ 2.*x*x*x	+ 2.5*x*x*x*x);
	default: mbs->UO() << "ANCFBeamBE2D::GetS0x Error: Shape function " << shape << " not available!\n";	return 0;
	}
};

double ANCFBeamBE2D::GetS0xx(const Vector2D& p_loc, int shape) const
{
	double x = p_loc.X();
	switch (shape)
	{
	case 1: return 0.25*(		+ 15.*x				- 15.*x*x*x);
	case 2:	return 0.25*(-3.	+ 15.*x	+ 3.*x*x	- 15.*x*x*x);
	case 3:	return 0.25*(-1.	+  3.*x	+ 3.*x*x	-  5.*x*x*x);
	case 4:	return 0.25*(		- 15.*x				+ 15.*x*x*x);
	case 5:	return 0.25*( 3.	+ 15.*x	- 3.*x*x	- 15.*x*x*x);
	case 6:	return 0.25*(-1.	-  3.*x	+ 3.*x*x	+  5.*x*x*x);
	default: mbs->UO() << "ANCFBeamBE2D::GetS0xx Error: Shape function " << shape << " not available!\n";	return 0;
	}
};
#pragma endregion

void ANCFBeamBE2D::DefineIntegrationRule(IntegrationRule& integrationRule)
{
	assert(integrationRule.settings.elementType == TFE_ThinBeam2D);

	int ruleOrder = 0;

	if (integrationRule.settings.integratedValueType == IntegrationRule::IVT_Load)
		ruleOrder = 5;
	else
		ruleOrder = 10;

	IntegrationRule::DefineIntegrationRuleLine(integrationRule, ruleOrder);
};