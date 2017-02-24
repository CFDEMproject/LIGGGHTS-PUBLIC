//#**************************************************************
//#
//# filename:             ANCFAxMovBeam2D.cpp
//#
//# author:               Astrid Sinwel, Gerstmayr Johannes
//#
//# generated:						August 4, 2009
//# description:          axially moving 2D ANCF - element
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
#include "constraint.h"
#include "material.h"
#include "ancfcable2d.h"
#include "ancfaxmovbeam2d.h"
#include "femathhelperfunctions.h"
#include "graphicsconstants.h"


extern int epsmode, kappamode, usekappafabs;

// init dets, kappa0, eps0
void ANCFAxMovBeam2D::Initialize()
{
	int ns = NS();
	order_kappa = 5; //regular: 5
	order_eps = 9;   //regular: 9
	ConstVector<ANCFCable2D_MaxIP> xk,xe,wk,we;
	GetIntegrationRule(xk,wk,order_kappa); //3 leads to 0.01%error, 5 leads to 1e-5 relative error
	GetIntegrationRule(xe,we,order_eps);   //7 leads to 0.1% error

	kappa0.SetLen(xk.Length());
	kappa0.SetAll(0.);
	eps0.SetLen(xe.Length());
	eps0.SetAll(0.);
	det_intrule_kappa.SetLen(xk.Length());
	det_intrule_eps.SetLen(xe.Length());

	Vector SV(0);
	for (int i1=1; i1<= xk.Length(); i1++)
	{
		double x = xk(i1);
		kappa0(i1) = GetKappa(x,x_init);
		Vector2D rx0 = GetPosx2D(x*GetLx()*0.5, x_init);
		double rxn = rx0.Norm();
		det_intrule_kappa(i1) = 0.5*GetLx()*rxn; //is approximately 0.5*GetLx(); in straight elements exactly 0.5*GetLx()!!!
	}

	for (int i1=1; i1<= xe.Length(); i1++)
	{
		double x = xe(i1);
		eps0(i1) = GetEpsAxial(x,x_init);
		Vector2D rx0 = GetPosx2D(x*GetLx()*0.5, x_init);
		double rxn = rx0.Norm();
		det_intrule_eps(i1) = 0.5*GetLx()*rxn; //is approximately 0.5*GetLx(); in straight elements exactly 0.5*GetLx()!!!
	}
	ANCFCable2D::Initialize();
}

//M = int(rho*((S)^T).S, dV,V)
void ANCFAxMovBeam2D::EvalM(Matrix& m, double t)
{
	//double dummy;
	if (massmatrix.Getcols() == SOS())
	{
		m = massmatrix;
		return;
	}
	else
	{
		m.SetAll(0);

		int dim = Dim();
		int ns = NS();

		ConstVector<ANCFCable2D_MaxNShapes> SV;

		SV.SetLen(ns);

		ConstMatrix<2*ANCFCable2D_MaxDOF> HL(SOS(),dim);

		static Vector x1,w1;

		GetIntegrationRule(x1,w1,6); //optimal: 6x3x3, 6x1x1 geht auch!!!!

		double rhoA = GetLy()*GetLz() * Rho();
		if (IsBeamParameters()) rhoA = GetRhoA();

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

			Vector2D rx0 = GetPosx2D(p*GetLx()*0.5, x_init);
			double rxn = rx0.Norm();
			double det = 0.5*GetLx()*rxn; //is approximately 0.5*GetLx(); in straight elements exactly 0.5*GetLx()!!!

			m += det * rhoA * w1(i1) * (HL*HL.GetTp());
		}

		m(1,1) += concentratedmass1;
		m(2,2) += concentratedmass1;
		m(1+4,1+4) += concentratedmass2;
		m(2+4,2+4) += concentratedmass2;
		massmatrix = m;
		//UO() << "m-invert=" << m.Invert2() << "\n";
	} 

	//UO() << "gesmasse=" << dummy << "\n";
};

//externe Kräfte
void ANCFAxMovBeam2D::EvalF2(Vector& f, double t)
{
	Body2D::EvalF2(f,t);
	TMStartTimer(22);
	int sos = SOS();
	static Vector fadd;

	int ns = NS();
	int dim = Dim();

	ConstVector<ANCFCable2D_MaxDOF> temp(SOS()), xg(SOS());
	ConstVector<ANCFCable2D_MaxDOF> xgp(SOS());
	ConstVector<ANCFCable2D_MaxNShapes> SV, SVx, SVxx;

	fadd.SetLen(SOS());
	fadd.SetAll(0);

	temp.SetLen(SOS());

	SV.SetLen(NS());
	SVx.SetLen(NS());
	SVxx.SetLen(NS());

	//
	//  Beam equations
	//

	double EI_w = Em() * GetLz()*Cub(GetLy())/12.;
	double EA = Em() * GetLy() * GetLz();
	if (IsBeamParameters())
	{
		EI_w = GetBeamEIy();
		EA = GetBeamEA();
	}
	double axialvelocity = AxialVelocity();

	ConstVector<ANCFCable2D_MaxIP> xk,xe,wk,we;

	GetIntegrationRule(xk,wk,order_kappa); //3 leads to 0.01%error, 5 leads to 1e-5 relative error

	GetIntegrationRule(xe,we,order_eps);   //7 leads to 0.1% error

	double loadfact_initcurv = GetInitLoadfact();


	//SetComputeCoordinates();
	//for (int i=1; i<=SOS(); i++)
	//	xgp(i) = XGP(i);
	GetCoordinates(xg);
	GetCoordinatesP(xgp);

	//kappa*delta kappa:
	for (int i1=1; i1<=xk.GetLen(); i1++)
	{
		double x = xk(i1);

		//bending moments:
		// factor EI * weight * det
		double factkappa = EI_w*wk(i1)*det_intrule_kappa(i1); 
		double kappa_plast = GetPlasticKappa(x*0.5*GetLx());
		double kappa;
		GetDeltaKappa(x,xg,temp, kappa);

		temp *= (kappa - loadfact_initcurv*kappa0(i1) - kappa_plast)*factkappa;
		fadd += temp;
	}

	//eps*delta eps + kinetic volume terms
	double rhoA = GetLy()*GetLz() * Rho();
	if (IsBeamParameters()) rhoA = GetRhoA();
	int optimized = 1;
	for (int i1=1; i1<=xe.GetLen(); i1++)
	{
		// x in (-1,1)
		double x = xe(i1);
		//axial strains:
		double facteps = EA*we(i1)*det_intrule_eps(i1);

		if (optimized)
		{
			//compute deltaeps:
			int dim = Dim();
			int ns = NS();

			GetS0xx(SVxx,x); // SVxx = d2S/dx2
			GetS0x(SVx,x);  // SVx = dS/dx
			GetS0(SV,x);  // SV = S shape


			//compute rx:
			Vector2D rx(0.,0.);
			Vector2D rt(0.,0.); //dr/dt
			Vector2D rxt(0.,0.); //drx/dt
			Vector2D rxx(0.,0.); //d2r/dx2
			for (int j = 1; j <= NS(); j++)
			{
				int k = (j-1)*Dim();
				rx(1) += SVx(j)*xg(k+1);
				rx(2) += SVx(j)*xg(k+2);
				rxx(1) += SVxx(j)*xg(k+1);
				rxx(2) += SVxx(j)*xg(k+2);
#ifndef COMPILE_AND
				// for strain damping
				rt(1) += SV(j)*XGP(k+1);
				rxt(1) += SVx(j)*XGP(k+1);
				rt(2) += SV(j)*XGP(k+2);
				rxt(2) += SVx(j)*XGP(k+2);
#endif
			}


			double fact = 1;
			double eps;
			double eps_plast = GetPlasticAxStrain(x*0.5*GetLx());

			double rxnorm = rx.Norm();
			// strain - initial strain - plastic strain
			eps = (rxnorm-1.0) - loadfact_initcurv*eps0(i1) - eps_plast;
			fact = 1./rxnorm * (eps * facteps);
#ifndef COMPILE_AND
			// for strain damping
			double epsp = rxt*rx / rx.Norm() + AxialVelocity(GetMBS()->GetTime())* rxx*rx / rx.Norm();
			double factepsp = GetEpsDamping()*we(i1)*det_intrule_eps(i1);
			fact += 1./rxnorm * (epsp*factepsp);
#endif

			double rhoAv_weightdet = rhoA * axialvelocity * we(i1)* det_intrule_eps(i1);
			Vector2D rtnew = rt+axialvelocity*rx;

			for (int i=1; i <= dim; i++)
			{
				for (int j=1; j <= ns; j++)
				{
					// EA eps * deltaeps
					fadd((j-1)*dim+i) += SVx(j)*rx(i)*fact;
					//  dT/dqi
					fadd((j-1)*dim+i) += -rhoAv_weightdet * ( rtnew(i)*SVx(j) );	

#ifndef COMPILE_AND
					// d/dt(dT/dqp)
					fadd((j-1)*dim+i) += rhoAv_weightdet * ( rxt(i)*SV(j) ); 
#endif
				}
			}




		}
		else // not optimized
		{
			GetDeltaEpsAxial(x,xg,temp);
			temp *= GetEpsAxial(x,xg)*facteps;
			fadd += temp;
		}
	}
	//
	//  Kinetic terms
	//

	// volume terms -- in optimized case already done before!
	if (!optimized)
	{
		GetIntegrationRule(xe,we,order_eps);   //7 leads to 0.1% error
		for (int i1=1; i1<=xe.GetLen(); i1++)
		{
			double x = xe(i1);


			//axial strains:
			double rhoAv_weightdet = rhoA * axialvelocity * we(i1)*det_intrule_eps(i1);

			GetS0x(SVx,x);  // SVx = dS/dx
			GetS0(SV,x);  // SV = S
			//compute rx:
			Vector2D rx(0.,0.); //dr/dx
			Vector2D rt(0.,0.); //dr/dt
			Vector2D rxt(0.,0.); //drx/dt
			Vector2D tau(0.,0.); // tangent vector tau

			for (int i = 1; i <= dim; i++)
			{
				for (int j = 1; j <= ns; j++)
				{
					rx(i) += SVx(j)*xg((j-1)*dim+i);
					rt(i) += SV(j)*XGP((j-1)*dim+i);
					rxt(i) += SVx(j)*XGP((j-1)*dim+i);
				}
			}

			Vector2D rtnew = rt+axialvelocity*rx;
			for (int i=1; i <= dim; i++)
			{
				for (int j=1; j <= ns; j++)
				{
					//  dT/dqi
					fadd((j-1)*dim+i) += -rhoAv_weightdet * ( rtnew(i)*SVx(j) );	

					// d/dt(dT/dqp)
					fadd((j-1)*dim+i) += rhoAv_weightdet * ( rxt(i)*SV(j) ); 

				}
			}

		}
	}

	////global_uo << "Volume terms: f = " << fadd << "\n";
	//f -= fadd;
	//fadd.SetAll(0.);

	//surface forces:
	for (int ii=1; ii <= 2; ii++)
	{
		double x;
		double sign=1;

		if (ii==1)
		{
			sign=-1;
			x=-1;
		}
		else
		{
			sign=1;
			x=1;
		}

		GetS0x(SVx,x); //  dS/dx
		GetS0(SV,x); //  S

		//compute rx:
		Vector2D rx(0.,0.); //dr/dx
		Vector2D rt(0.,0.); //dr/dx
		for (int i = 1; i <= dim; i++)
		{
			for (int j = 1; j <= ns; j++)
			{
				rx(i) += SVx(j)*xg((j-1)*dim+i);
#ifndef COMPILE_AND
				rt(i) += SV(j)*XGP((j-1)*dim+i);
#endif
			}
		}

		double rhoA = GetLy()*GetLz() * Rho();
		if (IsBeamParameters()) rhoA = GetRhoA();
		double rhoAv = sign*rhoA*axialvelocity;

		for (int i=1; i <= dim; i++)
		{
			for (int j=1; j <= ns; j++)
			{
				// Surface Term   rhoA v ( v rx + rt )^T S
				fadd((j-1)*dim+i) += rhoAv* ( axialvelocity*rx(i) + rt(i) )*SV(j);
			}
		}
	}

	//global_uo << "Surface terms: f = " << fadd << "\n";

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
		double k=GetMassDamping();
		if (GetMBS()->GetTime() < 0.1) k=10;
		temp *= k;
		f -= temp;
	}

};

// AND: use stiffness matrix, otherwise dont use it (does not contain damping terms!!!)
#ifdef COMPILE_AND
int ANCFAxMovBeam2D::FastStiffnessMatrix() const {return 2;}
#else
int ANCFAxMovBeam2D::FastStiffnessMatrix() const {return 0;}
#endif

void ANCFAxMovBeam2D::StiffnessMatrix(Matrix& m)
{
	TMStartTimer(23);
	int ns = NS();
	int dim = Dim();

	ConstVector<ANCFCable2D_MaxNShapes> SV, SVx, SVxx;
	SV.SetLen(NS());
	SVx.SetLen(NS());
	SVxx.SetLen(NS());

	ConstVector<ANCFCable2D_MaxDOF> xg;

	for (int i=1; i<=SOS(); i++)
		for (int j=1; j<=2*SOS(); j++)
			m(i,j) = 0;

	//
	//  Beam equations
	//

	double EI_w = Em() * GetLz()*Cub(GetLy())/12.;
	double ET_Iw = GetMaterial().TangentModule() * GetLz()*Cub(GetLy())/12.;
	double EA = Em() * GetLy() * GetLz();
	double ET_A = GetMaterial().TangentModule() *GetLy()*GetLz();
	if (IsBeamParameters())
	{
		EI_w = GetBeamEIy();
		EA = GetBeamEA();
	}

	double axialvelocity = AxialVelocity();

	static Vector dkappa;
	dkappa.SetLen(SOS());
	dkappa.SetAll(0.);
	static Matrix ddkappa;
	ddkappa.SetSize(SOS(), SOS());
	ddkappa.SetAll(0.);

	static Vector xk,xe,wk,we;

	GetIntegrationRule(xk,wk,order_kappa); //3 leads to 0.01%error, 5 leads to 1e-5 relative error

	GetIntegrationRule(xe,we,order_eps);   //7 leads to 0.1% error

	double loadfact_initcurv = GetInitLoadfact();

	//SetComputeCoordinates();
	GetCoordinates(xg);

	//kappa*delta kappa:
	for (int i1=1; i1<=xk.GetLen(); i1++)
	{
		double x = xk(i1);
		GetS0x(SVx,x);  // SVx = dS/dx
		GetS0xx(SVxx,x);  // SVxx = d2S/dx2


		//bending moments:
		double factkappa = EI_w*wk(i1)*det_intrule_kappa(i1);

		double kappa, kappaplast;
		kappaplast = GetPlasticKappa(x*GetLx()*0.5);
		GetDeltaDeltaKappa(x,xg,SVx,SVxx,kappa, dkappa, ddkappa);
		kappa -= kappaplast+loadfact_initcurv*kappa0(i1);


		//matrix(i,j) = factor*(delta_k(j))*delta_k(i) + factor*deltadelta_k(i,j)*(kappa-kappaplast)
		for (int i=1; i<=SOS(); i++)
		{
			for (int j=1; j<i; j++)
			{
				double matrixentry = factkappa*(dkappa(j)*dkappa(i) + ddkappa(i,j)*(kappa));
				m(i,j) -= matrixentry;
				m(j,i) -= matrixentry;
			}
			m(i,i) -= factkappa*(dkappa(i)*dkappa(i) + ddkappa(i,i)*(kappa));
		}
	}


	//eps*delta eps + Kinetic volume terms
	double rhoA = GetLy()*GetLz() * Rho();
	if (IsBeamParameters()) rhoA = GetRhoA();
	for (int i1=1; i1<=xe.GetLen(); i1++)
	{
		double x = xe(i1);

		GetS0x(SVx,x);  // SVx = dS/dx

		// Internal forces
		// deps : deps + ddeps(eps - epsplast - eps0)
		double epsax;
		GetDeltaDeltaEpsAxial(x,xg, SVx, epsax, dkappa, ddkappa);
		double epsaxplast = GetPlasticAxStrain(x*GetLx()*0.5);
		epsax -= epsaxplast + loadfact_initcurv*eps0(i1);
		//axial strains:
		double facteps = EA*we(i1)*det_intrule_eps(i1);


		for (int i=1; i<=SOS(); i++)
		{
			for (int j=1; j<i; j++)
			{
				double matrixentry = facteps*(dkappa(j)*dkappa(i) + ddkappa(j,i)*epsax);
				m(i,j) -= matrixentry;
				m(j,i) -= matrixentry;
			}
			m(i,i) -= facteps*(dkappa(i)*dkappa(i) + ddkappa(i,i)*epsax);
		}


		//
		//  Kinetic volume terms
		//
		double rhoAv_weightdet = rhoA * axialvelocity * we(i1)*GetLx()*0.5;

		// only terms for i==j
		// use symmetry of matrix!!
		// -rho A v^2*dr'*dr'
		for (int di=1; di <= ns; di++)
		{
			// Entry (di, di)
			double matrixentry = -rhoAv_weightdet * axialvelocity*SVx(di)*SVx(di) ;	
			m((di-1)*dim+1,(di-1)*dim+1) -= matrixentry;
			m((di-1)*dim+2,(di-1)*dim+2) -= matrixentry;

			for (int dj=1; dj < di; dj++)
			{
				matrixentry = -rhoAv_weightdet * axialvelocity*SVx(dj)*SVx(di) ;	
				//  dT/dqi
				// Entry (di, dj)
				m((di-1)*dim+1,(dj-1)*dim+1) -= matrixentry;
				m((di-1)*dim+2,(dj-1)*dim+2) -= matrixentry;

				// Entry (dj, di)
				m((dj-1)*dim+1,(di-1)*dim+1) -= matrixentry;
				m((dj-1)*dim+2,(di-1)*dim+2) -= matrixentry;
			}
		}
	}


	//Kinetic surface forces:
	for (int ii=1; ii <= 2; ii++)
	{
		double x;
		double sign=1;

		if (ii==1)
		{
			sign=-1;
			x=-1;
		}
		else
		{
			sign=1;
			x=1;
		}

		GetS0x(SVx,x); //  dS/dx
		GetS0(SV,x); //  S

		double rhoA = GetLy()*GetLz() * Rho();
		if (IsBeamParameters()) rhoA = GetRhoA();
		double sign_rhoAv2 = sign*rhoA*axialvelocity*axialvelocity;

		for (int di=1; di <= ns; di++)
		{
			for (int dj=1; dj <= ns; dj++)
			{
				double matrix_entry = sign_rhoAv2*SVx(dj)*SV(di);
				for (int i=1; i <= dim; i++)
				{
					// Surface Term   rhoA v ( v drx )^T dr
					m((di-1)*dim+i,(dj-1)*dim+i) -= matrix_entry;
				}
			}

		}
	}


	if (loads.Length() > 0 || constraints.Length() > 0)
	{
	static Vector load(SOS()), load2(SOS());
	load.SetLen(SOS());
	load2.SetLen(SOS());
	load.SetAll(0);
	load2.SetAll(0);

	for (int il=1; il <= loads.Length(); il++)
	{
		//loads(il)->AddElementLoad(load,GetMBS()->GetTime());		//$ DR 2012-10: loads moved from element to mbs, old code
		GetMBS()->GetLoad(loads(il)).AddElementLoad(load,GetMBS()->GetTime(),this);		//$ DR 2012-10: loads moved from element to mbs
	}
	for (int il=1; il <= constraints.Length(); il++)
	{
		constraints(il)->AddElementCqTLambda(GetMBS()->GetTime(),constraintindices(il),load);
	}

	double deps2=1e-10;
	for (int i=1; i <= SOS(); i++)
	{
		XG(i) += deps2;
		load2.SetAll(0);
		for (int il=1; il <= loads.Length(); il++)
		{
			//loads(il)->AddElementLoad(load2,GetMBS()->GetTime());			//$ DR 2012-10: loads moved from element to mbs, old code
			GetMBS()->GetLoad(loads(il)).AddElementLoad(load2,GetMBS()->GetTime(),this);		//$ DR 2012-10: loads moved from element to mbs
		}
		for (int il=1; il <= constraints.Length(); il++)
		{
			constraints(il)->AddElementCqTLambda(GetMBS()->GetTime(),constraintindices(il),load2);
		}
		load2 -= load;
		load2 *= 1./deps2;

		for (int j=1; j<=SOS(); j++)
			m(j,i) += load2(j);
		XG(i) -= deps2;
	}
	}

	TMStopTimer(23);
}

//
void ANCFAxMovBeam2D::DrawElement()
{
	ANCFCable2D::DrawElement();
};


Vector2D ANCFAxMovBeam2D::GetVel2D(const Vector2D& p_loc) const
{
	Vector2D rx, rt;
	rx = GetPosx2D(p_loc);
	rt = ANCFCable2D::GetVel2D(p_loc);
	return AxialVelocity(GetMBS()->GetTime())*rx + rt;
}

Vector2D ANCFAxMovBeam2D::GetVel2DD(const Vector2D& p_loc) const
{
	Vector2D rx, rt;
	rx = GetPosx2DD(p_loc);
	rt = ANCFCable2D::GetVel2DD(p_loc);
	return AxialVelocity(GetMBS()->GetDrawTime())*rx + rt;
}

Vector2D ANCFAxMovBeam2D::GetGeneralizedVel(const Vector2D& ploc, Vector& xgen, double time) const
{
	double p0=ploc.X()/(0.5*GetLx());// ploc€[-GetLx()/2,GetLx()/2]=>p0€[-1,1]
	ConstVector<ANCFCable2D_MaxNShapes> SVx(NS()), SV(NS());
	GetS0(SV, p0);
	GetS0x(SVx, p0);
	Vector2D rx(0.,0.), rt(0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			rx(i) += SVx(j)*xgen((j-1)*Dim()+i);  //p = S.e, NS() = 4, Dim() = 2
			rt(i) += SV(j)*XGP((j-1)*Dim()+i);  //p = S.e, NS() = 4, Dim() = 2
		}
	}
	return AxialVelocity(time)*rx;
}

Vector2D ANCFAxMovBeam2D::GetGeneralizedVelD(const Vector2D& ploc, Vector& xgen, double time) const
{
	double p0=ploc.X()/(0.5*GetLx());// ploc€[-GetLx()/2,GetLx()/2]=>p0€[-1,1]
	ConstVector<ANCFCable2D_MaxNShapes> SVx(NS()), SV(NS());
	GetS0x(SVx, p0);
	GetS0(SV, p0);
	Vector2D rx(0.,0.), rt(0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			rx(i) += SVx(j)*xgen((j-1)*Dim()+i);  //p = S.e, NS() = 4, Dim() = 2
			rt(i) += SV(j)*XGP((j-1)*Dim()+i);  //p = S.e, NS() = 4, Dim() = 2
		}
	}
	return AxialVelocity(time)*rx;
}

double ANCFAxMovBeam2D::GetEpsInit(const Vector2D& ploc, int flagD) const
{
	if (!UseInitialCurvature()) return 0;

	double x = ploc.X();
	double x0 = x*2./GetLx();

	double eps;
	if (flagD) eps = GetEpsAxialD(x0, x_init);
	else eps = GetEpsAxial(x0, x_init);
	//Vector2D rx;
	//if (flagD) rx = GetPosx2DD(x,xg_init);
	//else rx = GetPosx2D(x, xg_init);

	////1==Green, 2=Biot=Dmitrochenko, 3=linearized strain
	//double eps;
	//if (epsmode == 1) eps = 0.5*(rx*rx-1.0);
	//if (epsmode == 2) eps = (sqrt(rx*rx)-1.0);

	double yb = ploc.Y();
	double kappa;
	if (flagD) kappa = GetKappaD(x0,x_init);
	else kappa = GetKappa(x0,x_init);

	return GetInitLoadfact()*(eps - yb*kappa);
} 


void ANCFAxMovBeam2D::ComputeStress(const Vector2D& ploc, int type, int c1, int c2, double& val, bool flagD)
{
	if (type == 7 )   // acceleration
	{
		val = 0;
		return;
	}
	//ploc is from -1 .. +1

	//only for rectangular cross-section!:
	ConstVector<ANCFCable2D_MaxDOF> xgd;
	xgd.SetLen(SOS());

	if (flagD)
		GetDrawCoordinates(xgd);
	else
		GetCoordinates(xgd);

	if (type == 5) //show displacements
	{ 
		Vector2D u = GetDisplacement2D(ploc.X()*GetLx()*0.5, flagD);
		if (c1 == 1) val = u.X();
		else if (c1 == 2) val = u.Y();
		else if (c1 == 3) val = u.Norm();
		return;
	}

	if (type == 6 )   // velocity, actual velocity due to axial velocity v_0 * r' is computed
	{
		Vector2D p_loc_scal(ploc.X()*0.5*GetLx(), ploc.Y()*0.5*GetLy());
		Vector2D vel;
		if (flagD)
			vel = GetVel2DD(p_loc_scal);
		else
			vel = GetVel2D(p_loc_scal);
		if (c1 == 1) val = vel.X();
		else if (c1==2) val = vel.Y();
		else if (c1 == 3) val = vel.Norm();
		return;
	}


	double x0 = ploc.X();

	double yb = 0.5*ploc.Y()*GetLy();

	////  COMPUTE INITIAL CURVATURE
	double kappa0d = 0, eps0d = 0;
	double time;
	if (flagD)
		time = GetMBS()->GetDrawTime();
	else 
		time = GetMBS()->GetTime();

	if (UseInitialCurvature() )
	{
		if (flagD)
			kappa0d = GetKappaD(x0,x_init);
		else
			kappa0d = GetKappa(x0,x_init);
		kappa0d *= GetInitLoadfact();

		if (flagD)
			eps0d = GetEpsAxialD(x0,x_init);
		else
			eps0d = GetEpsAxial(x0,x_init);
		eps0d *= GetInitLoadfact();
	}

	//// INITIAL CURVATURE END 

	Vector2D ploc_scal(x0*0.5*GetLx(), yb);
	double kappa_plast = GetPlasticKappa(x0*0.5*GetLx(), flagD);
	double eps_plast_centerline = GetPlasticAxStrain(x0*0.5*GetLx(), flagD);
	double eps_plast_loc = GetPlasticStrain(ploc_scal, flagD);

	double kappa, eps_centerline;
	if (flagD)
	{
		kappa = GetKappaD(x0,xgd)-kappa0d;
		eps_centerline = GetEpsAxialD(x0,xgd) - eps0d; // (sqrt(rx*rx)-1.0) - eps0;
	}
	else
	{
		kappa = GetKappa(x0,xgd)-kappa0d;
		eps_centerline = GetEpsAxial(x0,xgd) - eps0d; // (sqrt(rx*rx)-1.0) - eps0;
	}


	double eps = eps_centerline - kappa*yb;

	Matrix3D stress(0);
	Matrix3D strain(0);
	stress(1,1) = (eps - eps_plast_loc)*Em();
	strain(1,1) = eps;

	double A = GetLy()*GetLz();
	double EI_w = Em() * GetLz()*Cub(GetLy())/12.;
	double EA = Em() * A; 
	if (IsBeamParameters())
	{
		EI_w = GetBeamEIy();
		EA = GetBeamEA();
	}

	Matrix plasticstrains, internalvars;
	if (GetMaterial().IsInelasticMaterial())
		GetPlasticStrainMatrix(plasticstrains, 1);
	if (GetMaterial().IsInelasticMaterial()&&GetMaterial().IsHardeningMaterial())
		GetInternalVariableMatrix(internalvars,	1);

	//if (internalvars(1,1))
	//global_uo << "internalvars.." << internalvars << "\n";

	if (type == 1 && c1+c2 == 0) //von Mises
	{
		val = stress.Mises();
	}
	else if (type == 2 && c1 == 1 && c2 == 1) //stress
	{
		val = EA * (eps_centerline-eps_plast_centerline);
		// HACK
		//val = stress(c1,c2); 
	}
	else if (type == 3 && c1 ==1 && c2 == 1) //strain
	{
		val = strain(c1,c2);
	}
	else if (type == 3 && c1 ==2 && c2 ==2) //strain
	{
		val = eps_centerline;
	}
	else if (type == 3 && c1 ==3 && c2 ==3) //strain
	{
		val = kappa;
	}
	else if (type == 8 && c1 == 1) // Nx, Normalkraft
	{
		val = EA * (eps_centerline-eps_plast_centerline);
	}
	else if (type == 9 && c1 == 3) // Mz, Biegemoment
	{
		val = -EI_w * (kappa-kappa_plast); 
	}
	else if (type == 4 && c1 == 2 && c2==2 && GetMaterial().IsInelasticMaterial() ) // inelastic strain  eps_xx^p
	{
		val = eps_plast_centerline; 
	}
	else if (type == 4 && c1 == 3 && c2==3 && GetMaterial().IsInelasticMaterial() ) // plastic curvature kappa^p
	{
		val = kappa_plast; 
	}
	else if (type == 4 && c1 == 1 && c2==1 &&GetMaterial().IsInelasticMaterial() ) // plastic  strain matrix
	{
		val = plasticstrains.Interpolate(ploc);
	}
	else if (type == 4 && ((c1==1&&c2==2)||(c1==2&&c2==1)) && GetMaterial().IsInelasticMaterial() && GetMaterial().IsHardeningMaterial() ) // internal variable matrix
	{
		val = internalvars.Interpolate(ploc);
	}
	else val = 0;

}



	// ------------------------------------------------------
	// VISUALIZATION

const char * str_plastic_curvature_amb = "plastic curvature";
const char * str_plastic_ax_strain_amb = "plastic axial strain";
const char * str_internal_variable_amb = "internal variable";

void ANCFAxMovBeam2D::GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables)
{
	ANCFCable2D::GetAvailableFieldVariables(variables);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_velocity,
		FieldVariableDescriptor::FVCI_y, true);
}

double ANCFAxMovBeam2D::GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector2D & local_position, bool flagD)
{
	//ploc is from -1 .. +1
	Vector2D ploc(local_position);

	//only for rectangular cross-section!:
	ConstVector<ANCFCable2D_MaxDOF> xgd;
	xgd.SetLen(SOS());

	if (flagD)
		GetDrawCoordinates(xgd);
	else
		GetCoordinates(xgd);

	// displacements
	if(fvd.VariableType() == FieldVariableDescriptor::FVT_displacement)
		return fvd.GetComponent(GetDisplacement2D(local_position.X()*GetLx()*0.5,flagD));

	// velocity, actual velocity due to axial velocity v_0 * r' is computed
	if(fvd.VariableType() == FieldVariableDescriptor::FVT_velocity)
	{
		Vector2D p_loc_scal(ploc.X()*0.5*GetLx(), ploc.Y()*0.5*GetLy());
		Vector2D vel;
		if (flagD)
			vel = GetVel2DD(p_loc_scal);
		else
			vel = GetVel2D(p_loc_scal);
		return fvd.GetComponent(vel);
	}


	double x0 = ploc.X();

	double yb = 0.5*ploc.Y()*GetLy();

	////  COMPUTE INITIAL CURVATURE
	double kappa0d = 0, eps0d = 0;
	double time;
	if (flagD)
		time = GetMBS()->GetDrawTime();
	else 
		time = GetMBS()->GetTime();

	if (UseInitialCurvature() )
	{
		if (flagD)
			kappa0d = GetKappaD(x0,x_init);
		else
			kappa0d = GetKappa(x0,x_init);
		kappa0d *= GetInitLoadfact();

		if (flagD)
			eps0d = GetEpsAxialD(x0,x_init);
		else
			eps0d = GetEpsAxial(x0,x_init);
		eps0d *= GetInitLoadfact();
	}

	//// INITIAL CURVATURE END 

	Vector2D ploc_scal(x0*0.5*GetLx(), yb);
	double kappa_plast = GetPlasticKappa(x0*0.5*GetLx(), flagD);
	double eps_plast_centerline = GetPlasticAxStrain(x0*0.5*GetLx(), flagD);
	double eps_plast_loc = GetPlasticStrain(ploc_scal, flagD);

	double kappa, eps_centerline;
	if (flagD)
	{
		kappa = GetKappaD(x0,xgd)-kappa0d;
		eps_centerline = GetEpsAxialD(x0,xgd) - eps0d; // (sqrt(rx*rx)-1.0) - eps0;
	}
	else
	{
		kappa = GetKappa(x0,xgd)-kappa0d;
		eps_centerline = GetEpsAxial(x0,xgd) - eps0d; // (sqrt(rx*rx)-1.0) - eps0;
	}


	double eps = eps_centerline - kappa*yb;

	Matrix3D stress(0);
	Matrix3D strain(0);
	stress(1,1) = (eps - eps_plast_loc)*Em();
	strain(1,1) = eps;

	double A = GetLy()*GetLz();
	double EI_w = Em() * GetLz()*Cub(GetLy())/12.;
	double EA = Em() * A; 
	if (IsBeamParameters())
	{
		EI_w = GetBeamEIy();
		EA = GetBeamEA();
	}

	Matrix plasticstrains, internalvars;
	if (GetMaterial().IsInelasticMaterial())
		GetPlasticStrainMatrix(plasticstrains, 1);
	if (GetMaterial().IsInelasticMaterial()&&GetMaterial().IsHardeningMaterial())
		GetInternalVariableMatrix(internalvars,	1);






	switch(fvd.VariableType())
	{
	case FieldVariableDescriptor::FVT_stress: return fvd.GetComponent(stress);
	case FieldVariableDescriptor::FVT_stress_mises: return stress.Mises();
	case FieldVariableDescriptor::FVT_total_strain: return strain(1,1);
	case FieldVariableDescriptor::FVT_beam_axial_extension: return eps_centerline;
	case FieldVariableDescriptor::FVT_beam_curvature: return kappa;
	case FieldVariableDescriptor::FVT_beam_force_axial: return EA * (eps_centerline-eps_plast_centerline);
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
				if(fvd.GetTextualIdentifierWithoutComponents() == str_plastic_curvature_amb)		// plastic curvature kappa^p
					return kappa_plast;
				if(fvd.GetTextualIdentifierWithoutComponents() == str_plastic_ax_strain_amb)		// plastic axial strain eps_ax^p
					return eps_plast_centerline;
				if(fvd.GetTextualIdentifierWithoutComponents() == str_internal_variable_amb && GetMaterial().IsHardeningMaterial())		// internal variable matrix
					return internalvars.Interpolate(local_position);
			}
		}
	}

	return FIELD_VARIABLE_NO_VALUE;
}

