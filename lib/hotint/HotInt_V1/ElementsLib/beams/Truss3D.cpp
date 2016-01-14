//#**************************************************************
//#
//# filename:             Truss3D.cpp
//#
//# author:               Michael Stangl
//#
//# generated:			  13.01.2012
//#
//# description:          3D Element Library
//#                       
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
 
#include "element.h"
#include "body3d.h"
#include "ANCFCable3d.h"
#include "Truss3D.h"
#include "femathhelperfunctions.h"
#include "graphicsconstants.h"


void Truss3D::GetS0(Vector& sf, const double& ploc) const
{
	sf.SetLen(NS());
	sf(1) = 0.5*(1.-ploc);
	sf(2) = 0.5*(1.+ploc);
}


void Truss3D::GetS0x(Vector& sfx, const double& ploc) const
{
	double xb = ploc;
	sfx.SetLen(NS());
	double f = 2./lx;
	sfx(1) = f*(-0.5);
	sfx(2) = f*0.5;

}

void Truss3D::GetS0xx(Vector& sfxx, const double& ploc) const
{
	double xb = ploc;
	sfxx.SetLen(NS());
	sfxx(1) = 0.;
	sfxx(2) = 0.;
}

void Truss3D::EvalM(Matrix& m, double t) 
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

		Matrix HL(SOS(),dim);

		Vector x1,w1;

		GetIntegrationRule(x1,w1,6); //optimal: 6x3x3, 6x1x1 geht auch!!!!

		for (int i1=1; i1<=x1.GetLen(); i1++)
		{
			double p = x1(i1);
			GetS0(SV,p);

			for (int i=0; i<ns; i++)
			{
				for (int j=1; j<=dim; j++)
				{
					HL(i*dim+j,j)=SV(i+1);
				}
			}

			m += 0.5*lx *ly*lz * rho * w1(i1) * (HL*HL.GetTp());
		}
		
		m(1,1) += concentratedmass1;
		m(2,2) += concentratedmass1;
		m(3,3) += concentratedmass1;
		m(1+3,1+3) += concentratedmass2;
		m(2+3,2+3) += concentratedmass2;
		m(3+3,3+3) += concentratedmass2;
		massmatrix = m;
		//UO() << "m=" << m << "\n";
		//UO() << "m-invert=" << m.Invert2() << "\n";
	}
};

double Truss3D::GetKappa(const double& x, const Vector& xg) const
{
	assert(0 && "Truss3D: GetKappa called!");return 0.;
}

void Truss3D::GetDeltaKappa(const double& x, const Vector& xg, Vector& dkappa, double& kappa) const
{
	assert(0 && "Truss3D: GetDeltaKappa called!");return;
}

void Truss3D::EvalF2(Vector& f, double t) 
{
	Body3D::EvalF2(f,t);
	TMStartTimer(22);
	int sos = SOS();
	static Vector fadd;

	int ns = NS();
	int dim = 3;

	fadd.SetLen(SOS());
	fadd.SetAll(0);

	temp.SetLen(SOS());

	SV.SetLen(NS());

	double EI_w = Em * lz*Cub(ly)/12.;
	double EA = Em * ly * lz;

	int order_eps = 1;   
	static Vector xk,xe,wk,we;

	GetIntegrationRule(xe,we,order_eps);   

	static Vector eps0;
	eps0.SetLen(xe.Length());

	for (int i=1; i <= sos; i++)
		xg(i) = x_init(i);

	for (int i1=1; i1<=xe.GetLen(); i1++)
	{
		double x = xe(i1);
		eps0(i1) = GetEpsAxial(x,xg);
	}

	SetComputeCoordinates();

	//eps*delta eps:
	for (int i1=1; i1<=xe.GetLen(); i1++)
	{
		double x = xe(i1);

		//axial strains:
		double facteps = EA*we(i1)*lx*0.5;

		int optimized = 1;
		if (optimized)
		{
			GetS0x(SV,x);
			//compute rx:
			Vector3D rx(0.,0.,0.);
			for (int i = 1; i <= Dim(); i++)
			{
				for (int j = 1; j <= NS(); j++)
				{
					rx(i) += SV(j)*xg((j-1)*3+i);
				}
			}

			//compute deltaeps:
			int dim = Dim();
			int ns = NS();

			int usedmitrochenko = 1; 
			double fact = 1;
			double eps;

			if (usedmitrochenko)
			{
				eps = (sqrt(rx*rx)-1.0) - eps0(i1);
				fact = 1./(sqrt(rx*rx));
			}
			else
			{
				eps = 0.5*(rx*rx-1.0) - eps0(i1);
			}

			rx *= fact * eps * facteps;
			for (int i=1; i <= dim; i++)
			{
				for (int j=1; j <= ns; j++)
				{
					fadd((j-1)*dim+i) += SV(j)*rx(i);
				}
			}

		}
		else
		{
			GetDeltaEpsAxial(x,xg,temp);
			temp *= GetEpsAxial(x,xg)*facteps;
			fadd += temp;
		}
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

//---------------------------------------------------------------
//for Displacement/Stress/Strain/Stress Resultants-Visualization:
//---------------------------------------------------------------
double Truss3D::GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & local_position, bool flagD)
{	

	if (!flagD)	return FIELD_VARIABLE_NO_VALUE;

	Vector3D p0 = local_position;//-l/2 <= p0 <= +l/2

	xgd.SetLen(SOS());
	GetDrawCoordinates(xgd);  //xgd=XGD

	double eps_init = 0;
	double eps = 0;

	eps_init = GetEpsAxial(p0.X(),x_init);
	eps = GetEpsAxial(p0.X(), xgd) - eps_init; 

	double stress;
	stress = eps*Em;

	switch(fvd.VariableType())
	{
	case FieldVariableDescriptor::FVT_stress: return stress;
	case FieldVariableDescriptor::FVT_beam_axial_extension: return eps;
	case FieldVariableDescriptor::FVT_beam_force_axial: return stress*ly*lz;
	}

	return FIELD_VARIABLE_NO_VALUE;

}


void Truss3D::GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables)
{
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_stress,
		FieldVariableDescriptor::FVCI_z, true);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_axial_extension,
		FieldVariableDescriptor::FVCI_z, true);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_force_axial,
		FieldVariableDescriptor::FVCI_z, true);
}

void Truss3D::DrawElement()
{
	mbs->SetColor(col);

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
	else
	{
		colormode = 0; 
	}

	if (colormode)
	{

		double lx1 = 1; double ly1 = ly*GetMBS()->GetMagnifyYZ();

		double tiling = 2.*pow(2.,GetMBS()->GetDrawResolution());

		Vector3D p1,p2,p1dir,p2dir;
		for (double i = 0; i < tiling; i++)
		{
			double l1 = -lx1+2*lx1*i/tiling;
			double l2 = -lx1+2*lx1*(i+1)/tiling;
			
			if (i == 0)
			{
				p1 = GetPos0D(l1);
				p1dir = GetPosx0D(l1);

				mbs->SetColor(col);
				if (concentratedmass1 != 0) 
				{
					mbs->SetColor(colred);
					//double r = pow(3.*concentratedmass1/(4.*MY_PI*GetRho()),1./3.); 
					double r = ly/10.;	//$ DR 2013-02-04 deleted rho from class element, do not use it here!
					mbs->DrawSphere(p1,r,12);
					mbs->SetColor(col);
				}
			}
			else
			{
				p1 = p2;
				p1dir = p2dir;
			}
			p2 = GetPos0D(l2);
			p2dir = GetPosx0D(l2);

			Vector3D tt = 0.1*(p2-p1);

			Vector3D ploc(.5*(l1+l2));
			double val = GetFieldVariableValue(*GetMBS()->GetActualPostProcessingFieldVariable(), ploc, true);
			mbs->DrawColorZyl(val,p1-tt,p2+tt,0.5*ly1,16);
		}

		if (concentratedmass2 != 0) 
		{
			mbs->SetColor(colred);
			//double r = pow(3.*concentratedmass2/(4.*MY_PI*GetRho()),1./3.); 
			double r = ly/10.;	//$ DR 2013-02-04 deleted rho from class element, do not use it here!
			//double r = ly*GetMBS()->GetMagnifyYZ()*1;
			mbs->DrawSphere(p2,r,12);
			//mbs->DrawSphere(p2,ly1*2,8);
		}

	}
	else
	{
		ANCFCable3D::DrawElement();
	}
};

