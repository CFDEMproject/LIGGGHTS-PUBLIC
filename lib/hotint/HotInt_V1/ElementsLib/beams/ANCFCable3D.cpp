//#**************************************************************
//#
//# filename:             ANCFCable3D.cpp
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						9. November 2004
//# description:          Driver and model for timeintegration
//#                       Model of a rigid arm and a hydraulic zylinder
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
 
#include "body3d.h"
#include "femathhelperfunctions.h"
#include "material.h"
#include "ANCFCable3D.h"
#include "Node.h"
#include "graphicsconstants.h"
#include "elementdataaccess.h"
//#include "solversettings_auto.h"

void ANCFCable3D::LinkToElements()
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
void ANCFCable3D::BuildDSMatrices() 
{
	//orderx = 9; //max 9x5, sonst array grad zu klein!!!!
	//GetIntegrationRule(x1,w1,orderx);

};

void ANCFCable3D::GetS0(Vector& sf, const double& ploc) const
{
	double xb = ploc;
	sf.SetLen(NS());
	double xb2 = xb*xb;
	double xb3 = xb2*xb;
	sf(1) = 1.0/2.0-3.0/4.0*xb+xb3/4.0;
	sf(2) = (1.0-xb-xb2+xb3)*lx/8.0;
	sf(3) = 1.0/2.0+3.0/4.0*xb-xb3/4.0;
	sf(4) = Sqr(1.0+xb)*(-1.0+xb)*lx/8.0;
}


void ANCFCable3D::GetS0x(Vector& sfx, const double& ploc) const
{
	double xb = ploc;
	sfx.SetLen(NS());
	double f = 2./lx;
	sfx(1) = f*(-3.0/4.0+3.0/4.0*xb*xb);
	sfx(2) = f*((-1.0-2.0*xb+3.0*xb*xb)*lx/8.0);
	sfx(3) = f*(3.0/4.0-3.0/4.0*xb*xb);
	sfx(4) = f*((1.0+xb)*(-1.0+xb)*lx/4.0+pow(1.0+xb,2.0)*lx/8.0);
}

void ANCFCable3D::GetS0xx(Vector& sfxx, const double& ploc) const
{
	double xb = ploc;
	sfxx.SetLen(NS());
	double f = 4./Sqr(lx);
	sfxx(1) = f*3.0/2.0*xb;
	sfxx(2) = f*(-2.0+6.0*xb)*lx/8.0;
	sfxx(3) = -f*3.0/2.0*xb;
	sfxx(4) = f*((-1.0+xb)*lx/4.0+(1.0+xb)*lx/2.0);
}

Vector3D ANCFCable3D::GetPos(const double& p_loc) const
{
	static Vector xg;
	xg.SetLen(SOS());
	GetCoordinates(xg);
	double p0=p_loc/(0.5*lx);
	static Vector SV;
	GetS0(SV, p0);
	Vector3D p(0.,0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += SV(j)*xg((j-1)*3+i);
		}
	}
	return p;
};

Vector3D ANCFCable3D::GetVel(const double& p_loc) const
{
	static Vector xg;
	xg.SetLen(SOS());
	GetCoordinatesP(xg);
	double p0=p_loc/(0.5*lx);
	static Vector SV;
	GetS0(SV, p0);
	Vector3D p(0.,0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += SV(j)*xg((j-1)*3+i);
		}
	}
	return p;
};

Vector3D ANCFCable3D::GetPosD(const double& p_loc) const
{
	static Vector xgd;
	xgd.SetLen(SOS());
	GetDrawCoordinates(xgd);
	double p0=p_loc/(0.5*lx);
	static Vector SV;
	GetS0(SV, p0);
	Vector3D p(0.,0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += SV(j)*xgd((j-1)*3+i);
		}
	}
	return p;
};

//in reference element coordinates (-1..1)
Vector3D ANCFCable3D::GetPos0D(const double& p_loc) const
{
	static Vector xgd;
	xgd.SetLen(SOS());
	GetDrawCoordinates(xgd);
	double p0=p_loc;
	static Vector SV;
	GetS0(SV, p0);
	Vector3D p(0.,0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += SV(j)*xgd((j-1)*3+i);
		}
	}
	return p;
};

Vector3D ANCFCable3D::GetVelD(const double& p_loc) const
{
	static Vector xgd;
	xgd.SetLen(SOS());
	GetDrawCoordinatesP(xgd);
	double p0=p_loc/(0.5*lx);
	static Vector SV;
	GetS0(SV, p0);
	Vector3D p(0.,0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += SV(j)*xgd((j-1)*3+i);
		}
	}
	return p;
};


void ANCFCable3D::GetH(Matrix& H) 
{
	UO() << "****************************************************" << "\n";
	UO() << "**********  ANCFCable3D::GetIntDuDq -> GetH ********" << "\n";
	UO() << "****************************************************" << "\n";
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

		SV.SetLen(ns);

		GetIntegrationRule(x1,w1,3); //3x1x1 !!!!!

		for (int i1=1; i1<=x1.GetLen(); i1++)
		{
			double p=x1(i1);
			GetS0(SV,p);

			double fact = ly*lz * lx*0.5 * w1(i1);
			//UO() << "scaled_weight = " << fact << "\n";
			//UO() << "xi = " << p << "\n";

			for (int i=0; i<ns; i++)
			{
				for (int j=1; j<=dim; j++)
				{
					H(i*dim+j,j)+=fact*SV(i+1);
				}
			}
		}
		//H=T.GetTp()*H;
		Hmatrix = H;
	}
}

void ANCFCable3D::EvalM(Matrix& m, double t) 
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
		//UO() << "********** Integration **********" << "\n";
		//UO() << "integration points = " << x1 << "\n";
		//UO() << "weights = " << w1 << "\n";

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
			//UO() << "********** scaled integration points **********" << "\n";
			//UO() << "scaled integration point = " << x1(i1) << "\n";
			//UO() << "********** H Matrix **********" << "\n";
			//UO() << "HL_pos = " << HL << "\n";

			m += 0.5*lx *ly*lz * rho * w1(i1) * (HL*HL.GetTp());
		}
		
		m(1,1) += concentratedmass1;
		m(2,2) += concentratedmass1;
		m(3,3) += concentratedmass1;
		m(1+6,1+6) += concentratedmass2;
		m(2+6,2+6) += concentratedmass2;
		m(3+6,3+6) += concentratedmass2;
		massmatrix = m;
		//UO() << "m_Cable=" << m << "\n";
		//UO() << "LTG=" << ltg << "\n";
		//UO() << "m-invert=" << m.Invert2() << "\n";
	}
};

double ANCFCable3D::GetKappa(const double& x, const Vector& xg) const
{
	Vector3D rx = GetPosx0(x,xg);
	Vector3D rxx = GetPosxx0(x,xg);
	double rxn = rx.Norm();

	if (rxn == 0) return 0;
	else
	{
		return (rx.Cross(rxx)).Norm()/Cub(rxn);
	}
}

double ANCFCable3D::GetEpsAxial(const double& x, const Vector& xg) const
{
	Vector3D rx = GetPosx0(x,xg);
	return 0.5*(rx*rx-1.0);
	//return sqrt(rx*rx)-1; //Dmitrochenko
}

double deps = 1e-8;
void ANCFCable3D::GetDeltaKappa(const double& x, const Vector& xg, Vector& dkappa, double& kappa) const
{
	//Vector oldkappa;
	if (1)
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

		//global_uo << "ERROR: Cable element needs to be corrected for DeltaKappa!!!!\n";
		//global_uo << "       Compare bending with Cable2D with respect to all 3 directions!\n";

		for (int i=1; i <= dim; i++)
		{
			for (int j=1; j <= ns; j++)
			{
				//dr,x/de x r,xx + r,x x dr,xx/de:
				switch (i) {
					case 1: 
						t1.X() = 0;     
						t1.Y() = -SVx(j)*rxx.Z()+SVxx(j)*rx.Z(); 
						t1.Z() =  SVx(j)*rxx.Y()-SVxx(j)*rx.Y(); break; //dy=dz=0
					case 2: 
						t1.X() =  SVx(j)*rxx.Z()-SVxx(j)*rx.Z(); 
						t1.Y() =  0;     
						t1.Z() = -SVx(j)*rxx.X()+SVxx(j)*rx.X(); break; //dx=dz=0
					case 3: 
						t1.X() = -SVx(j)*rxx.Y()+SVxx(j)*rx.Y(); 
						t1.Y() = +SVx(j)*rxx.X()-SVxx(j)*rx.X(); 
						t1.Z() =  0; break; //dx=dy=0
/* //original, in papers with Shabana, coincides with actual version:
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
*/
					default: ;
				}
				dg = (rx(i)*SVx(j)); //normed

				df = (rxcrxx*t1); //normed

				dkappa((j-1)*dim+i) = df*gn-fn*dg; //normed
			}
		}
		//oldkappa = dkappa;
	}
	else
	{
		static Vector xgdiff;
		xgdiff = xg;
		double val0;
		val0 = GetKappa(x,xgdiff);

		for (int i=1; i <= SOS(); i++)
		{
			xgdiff(i) += deps;
			dkappa(i) = GetKappa(x,xgdiff);
			dkappa(i) -= val0;
			dkappa(i) *= 1./deps;
			xgdiff(i) -= deps;
		}
		//global_uo << "kappaanal=" << oldkappa << "\n";
		//global_uo << "kappanum=" << dkappa << "\n";
	}
}
void ANCFCable3D::GetDeltaEpsAxial(const double& x, const Vector& xg, Vector& depsaxial) const
{
	if (1)
	{
		static Vector SVx;
		SVx.SetLen(NS());

		Vector3D rx = GetPosx0(x,xg);
		GetS0x(SVx,x);

		int dim = Dim();
		int ns = NS();
		for (int i=1; i <= dim; i++)
		{
			for (int j=1; j <= ns; j++)
			{
				depsaxial((j-1)*dim+i) = SVx(j)*rx(i);
			}
		}
		//global_uo << "deps1=" << depsaxial << ", rx=" << rx << "\n";
	}
	else
	{
		static Vector xgdiff;
		xgdiff = xg;
		double val0;
		val0 = GetEpsAxial(x,xgdiff);

		for (int i=1; i <= SOS(); i++)
		{
			xgdiff(i) += deps;
			depsaxial(i) = GetEpsAxial(x,xgdiff);
			depsaxial(i) -= val0;
			depsaxial(i) *= 1./deps;
			xgdiff(i) -= deps;
		}
		
		//global_uo << "deps2=" << depsaxial << "\n";
	}
}


void ANCFCable3D::EvalF2(Vector& f, double t) 
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

	int order_kappa = 5; //regular: 5
	int order_eps = 9;   //regular: 9
	static Vector xk,xe,wk,we;

	GetIntegrationRule(xk,wk,order_kappa); //3 leads to 0.01%error, 5 leads to 1e-5 relative error
	GetIntegrationRule(xe,we,order_eps);   //7 leads to 0.1% error

	//compute kappa0 and eps0
	static Vector kappa0;
	kappa0.SetLen(xk.Length());

	static Vector eps0;
	eps0.SetLen(xe.Length());

	for (int i=1; i <= sos; i++)
		xg(i) = x_init(i);

	for (int i1=1; i1<=xk.GetLen(); i1++)
	{
		double x = xk(i1);
		kappa0(i1) = GetKappa(x,xg);
	}

	for (int i1=1; i1<=xe.GetLen(); i1++)
	{
		double x = xe(i1);
		eps0(i1) = GetEpsAxial(x,xg);
	}

	SetComputeCoordinates();

	//kappa*delta kappa:
	for (int i1=1; i1<=xk.GetLen(); i1++)
	{
		double x = xk(i1);

		//bending moments:
		double factkappa = EI_w*wk(i1)*lx*0.5;

		if (1)
		{
			double kappa;
			GetDeltaKappa(x,xg,temp, kappa);
			//UO() << "delta_k=" << temp << "\n";
			temp *= (kappa - kappa0(i1))*factkappa;
			fadd += temp;
		}
	}

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

			//compute epsaxial*deltaeps*fact:
			//temp *= 0.5*(rx*rx-1.0)*facteps;
		}
		else
		{
			GetDeltaEpsAxial(x,xg,temp);
			//UO() << "delta_eps=" << temp << "\n";
			temp *= GetEpsAxial(x,xg)*facteps;
			fadd += temp;
		}
	}

	f -= fadd;
	//UO() << "********** ANCFCable3D::EvalF2 **********" << "\n";
	//UO() << "f_before_mass_damping = " << f << "\n";

	TMStopTimer(22);

	if (GetMassDamping() != 0)
	{
		UO() << "!!!!DAMPING!!!\n";
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
	//UO() << "********** ANCFCable3D::EvalF2 **********" << "\n";
	//UO() << "f_after_mass_damping = " << f << "\n";

}; 


void ANCFCable3D::DrawElement() 
{
	mbs->SetColor(col);

	double lx1 = 1; double ly1 = ly*GetMBS()->GetMagnifyYZ();

	double tiling = 2.*pow(2.,GetMBS()->GetDrawResolution());

	//lx1*=0.9;ly1*=0.9;lz1*=0.9;
	Vector3D p1,p2,p1dir,p2dir;
	for (double i = 0; i < tiling; i++)
	{
		double l1 = -lx1+2*lx1*i/tiling;
		double l2 = -lx1+2*lx1*(i+1)/tiling;
		if (i == 0)
		{
			p1 = GetPos0D(l1);
			p1dir = GetPosx0D(l1);

			//mbs->SetColor(colred);
			//mbs->DrawSphere(p1,0.65*ly1,4);
			mbs->SetColor(col);
			if (concentratedmass1 != 0) 
			{
				mbs->SetColor(colred);
				//double r = pow(3.*concentratedmass1/(4.*MY_PI*GetRho()),1./3.); 
				double r = ly/10.;	//$ DR 2013-02-04 deleted rho from class element, do not use it here!
				//double r = ly*GetMBS()->GetMagnifyYZ()*1;
				mbs->DrawSphere(p1,r,12);
				//mbs->DrawSphere(p1,ly1*2,8);
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

		//mbs->DrawZyl(p1,p2,p1dir,p2dir,0.5*ly1,0,0,16); //does not work well due to drawing problems ...
		Vector3D tt = 0.1*(p2-p1);
		mbs->DrawZyl(p1-tt,p2+tt,0.5*ly1,16);
	}

	if (use_director)
	{
		// draw local frame (e1 red, e2 green, e3 blue)
		Vector3D position, director;
		double drawlength = max(ly,lz)*2.;
		double linethickness = max(2.,0.04*mbs->GetDOption(120));
			
		position = GetPosD(Vector3D());   //in mid of element
		Matrix3D rotmat = GetRotMatrixD(0);
		for (int i=1; i<=3; i++)
		{
			Vector3D ei;
			ei(1) = rotmat(1,i);
			ei(2) = rotmat(2,i);
			ei(3) = rotmat(3,i);
			ei *= drawlength;
			
			Vector3D color;
			color(i) = 0.5;
			
			mbs->MyDrawLine(position, position+ei, linethickness, color);
		}
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
};


