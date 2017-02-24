//#**************************************************************
//#
//# filename:             FiniteElement3DFFRF.cpp
//#
//# author:               Gerstmayr Johannes, Sinwell Astrid, YV
//#
//# generated:						October 2010
//# description:          functionality of 3D finite elements related to FFRF & CMS
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
 
#include "finiteElement3D.h"
#include "Material.h"
#include "node.h"
#include "solversettings_auto.h"
#include "femathhelperfunctions.h"
#include "rigid3d.h"
#include "rigid3dkardan.h"
#include "femeshinterface.h"
#include "referenceframe3d.h"
#include "basecmselement.h"
#include "finiteelement3dffrf.h"

void FiniteElement3DFFRF::SetFFRFElement(int CMSelementI)
{ 
	if (CMSelementI==0)
	{
		SetFFRF(NULL);
		return;
	}
	// set FFRF: init ffrfdata, set type TFFRF
	SetFFRF(1); 
	ffrfdata->CMSelnum = CMSelementI;
	// set type TCMS in case of CMS
	if (ReferenceFrame().IsCMS())
		AddType(TCMS);
	if (ReferenceFrame().IsGCMS())
		AddType(TGCMS);
	
	// in GCMS-case, no change of x_init
	if (IsGCMS()) return;

	// new Init-Vector x_init: keep flexible init values, add reference frame init values at the right positions
	// reset length of x_init if dofs were added due to reference frame, redistribute x_init values
	if (x_init.Length() != SOS())
	{
		Vector x_init_save(x_init);
		x_init = Vector(2*SOS());
		// flexible init values
		for (int i=1; i<=FlexDOF(); i++)
		{
			x_init(i) = x_init_save(i);
			x_init(SOS()+i) = x_init_save(FlexDOF()+i);
		}
	}
	// reference frame init values
	for (int i=1; i<=FFRFDim(); i++)
	{
		x_init(FlexDOF()+i) = ReferenceFrame().GetXInit()(i);
		x_init(SOS()+FlexDOF()+i) = ReferenceFrame().GetXInit()(FFRFDim()+i);
	}
}

void FiniteElement3DFFRF::CopyFrom(const Element& e)
{
	FiniteElement3D::CopyFrom(e);
	const FiniteElement3DFFRF& ce = (const FiniteElement3DFFRF&)e;
	ffrfdata = new FFRFData(*ce.ffrfdata);
}

void FiniteElement3DFFRF::Initialize() 
{
	FiniteElement3D::Initialize();

	//stored functions:
	// are only computed in non-reduced FFRF case, in reduced CMS case, these precomputed matrices are not needed
	if (!IsCMS() )
	{
		//precompute mass matrix at time 0:
		ConstMatrix<FEmaxDOF*FEmaxDOF> tmp(SOS(),SOS());
		EvalMff(tmp, 0);

		GetI1(ffrfdata->I1S);
		GetSbar(ffrfdata->SbarS);

		//Matrix3D ikl;
		for (int k=1; k <= 3; k++)
		{
			for (int l=1; l <= 3; l++)
			{
				GetSbarkl(k, l, ffrfdata->SbarklS[k-1][l-1]);

				ffrfdata->IklS[k-1][l-1] = GetIkl(k, l);
				//ikl(k,l) = GetIkl(k, l);

				GetIbarkl(k, l, ffrfdata->IbarklS[k-1][l-1]);
			}
		}

		GetSbarkl(2,3,ffrfdata->SbarTilde[0]);
		GetSbarkl(3,2,tmp);
		ffrfdata->SbarTilde[0] -= tmp; //Sbar23-Sbar32

		GetSbarkl(3,1,ffrfdata->SbarTilde[1]);
		GetSbarkl(1,3,tmp);
		ffrfdata->SbarTilde[1] -= tmp; //Sbar31-Sbar13

		GetSbarkl(1,2,ffrfdata->SbarTilde[2]);
		GetSbarkl(2,1,tmp);
		ffrfdata->SbarTilde[2] -= tmp; //Sbar12-Sbar21
	}
}

void FiniteElement3DFFRF::SetFFRF(int isffrf)
{ 
	if (isffrf) 
		SetGeometricNonlinearityStatus(GNS_Linear);

	if (isffrf && ffrfdata) // if already ffrf, do nothing
		return;
	else if (isffrf && !ffrfdata) // create new ffrf-data
		ffrfdata = new FFRFData(GetMBS());
	else // delete old ffrfdata
	{ 
		type = TMBSElement(type & (~TCMS));   // set type to NOT TCMS
		type = TMBSElement(type & (~TGCMS));   // set type to NOT TGCMS
		delete ffrfdata;
		ffrfdata = NULL;
	}
}

void FiniteElement3DFFRF::LinkToElements()
{
	// LTG map for reduced CMS element
	if (IsCMS())
	{
		int Nnode = FlexDOF();

		LTGreset();

		//Position:
		for (int j = 1; j <= NNodes(); j++)
		{
			const Node& node = ReferenceFrame().GetNode(NodeNum(j));

			for (int i=1; i <= node.SOS(); i++)
			{
				int l=node.Get(i);
				AddLTG(node.Get(i));
			}
		}

		//+++++++++++++++++++++++++++++
		//velocities:
		for (int j = 1; j <= NNodes(); j++)
		{
			//AH: IMHO WRONG!!!! const Node& node = ReferenceFrame().GetNode(j);
			const Node& node = ReferenceFrame().GetNode(NodeNum(j));
			for (int i=1; i <= node.SOS(); i++)
			{
				AddLTG(node.Get(i+node.SOS()));
			}
		}
	}
	// LTG map for standard element or Floating Frame element
	else
	{
		if (SOSowned() != SOS())
		{

			int Nnonode = SOSowned();
			int Nnode = SOS()-Nnonode;
			TArray<int> storeltg(Nnonode);
			for (int i=1; i <= Nnonode*2; i++)
			{
				storeltg.Add(LTG(i));
			}
			LTGreset();

			// Node positions
			for (int j = 1; j <= NNodes(); j++)
			{
				const Node& node = GetNode(j);
				//Position:
				for (int i=1; i <= node.SOS(); i++)
				{
					AddLTG(node.Get(i));
				}
			}
			//reference frame, if FloatingFrame element, otherwise FFRFDim==0
			for (int i = 1; i <= FFRFDim(); i++)
			{
				AddLTG(ReferenceFrame().LTG(i));
			}
			// remaining, if exist
			for (int i=1; i <= Nnonode; i++)
			{
				AddLTG(storeltg(i));
			}

			// Node velocities
			for (int j = 1; j <= NNodes(); j++)
			{
				const Node& node = GetNode(j);
				//velocity:
				for (int i=1; i <= node.SOS(); i++)
				{
					AddLTG(node.Get(i+node.SOS()));
				}
			}
			//reference frame, if FloatingFrame element, otherwise FFRFDim==0
			for (int i = 1; i <= FFRFDim(); i++)
			{
				AddLTG(ReferenceFrame().LTG(i+FFRFDim()));
			}
			// remaining, if exist
			for (int i=1; i <= Nnonode; i++)
			{
				AddLTG(storeltg(i+Nnonode));
			}
		}
	}
}

void FiniteElement3DFFRF::GetI1(Vector& I1) //Shabana p. 206 ff.
{
	Vector3D I1local;
	Matrix3D jac;
	for (IntegrationPointsIterator ip(integrationRuleMass); !ip.IsEnd(); ++ip)
	{
		GetJacobi(jac, ip.Point());
		double jacdet = jac.Det();
		Vector3D pp = GetRefConfPos(ip.Point());
		I1local += fabs(jacdet) * Rho() * ip.Weight() * pp;
	}
	I1 = I1local;
}

double FiniteElement3DFFRF::GetIkl(int k, int l)
{
	double ikl = 0;

	Matrix3D jac;
	for (IntegrationPointsIterator ip(integrationRuleMass); !ip.IsEnd(); ++ip)
	{
		GetJacobi(jac, ip.Point());
		double jacdet = jac.Det();
		Vector3D pp = GetRefConfPos(ip.Point());
		ikl += fabs(jacdet) * Rho() * ip.Weight() * pp(k) * pp(l) ;
	}
	
	return ikl;
}

void FiniteElement3DFFRF::GetIbarkl(int k, int l, Vector& I1)
{
	int dim = Dim();
	int ns = NS();

	I1.SetLen(dim * ns);
	I1.SetAll(0);

	Matrix3D jac;
	for (IntegrationPointsIterator ip(integrationRuleMass); !ip.IsEnd(); ++ip)
	{
		GetJacobi(jac, ip.Point());
		double jacdet = jac.Det();
		Vector3D pp = GetRefConfPos(ip.Point());
		for (int i = 1; i <= ns; i++)
		{
			I1(dim*(i-1)+l) += fabs(jacdet) * Rho() * pp(k) * GetS0(ip.Point(),i) * ip.Weight();
		}
	}
}

void FiniteElement3DFFRF::GetSbar(Matrix& Sbar)
{
	int dim = Dim(); 
	int ns = NS();

	Sbar.SetSize(dim,dim*ns);
	Sbar.SetAll(0);

	Matrix3D jac;
	for (IntegrationPointsIterator ip(integrationRuleMass); !ip.IsEnd(); ++ip)
	{
		GetJacobi(jac, ip.Point());
		double jacdet = jac.Det();
		for (int i = 1; i <= ns; i++)
		{
			double Sbarval =  fabs(jacdet) * Rho() * GetS0(ip.Point(),i) * ip.Weight();
			Sbar(1,dim*(i-1)+1) += Sbarval;
			Sbar(2,dim*(i-1)+2) += Sbarval;
			Sbar(3,dim*(i-1)+3) += Sbarval;
		}
	}
}

void FiniteElement3DFFRF::GetSbarkl(int k, int l, Matrix& Sbar)
{
	int dim = Dim(); 
	int ns = NS();

	Sbar.SetSize(dim*ns,dim*ns);
	Sbar.SetAll(0);

	Matrix3D jac;
	for (IntegrationPointsIterator ip(integrationRuleMass); !ip.IsEnd(); ++ip)
	{
		GetJacobi(jac, ip.Point());
		double jacdet = jac.Det();
		for (int i = 1; i<= ns; i++)
		{
			for (int j = 1; j<= ns; j++)
			{
				Sbar(dim*(i-1)+k, dim*(j-1)+l) +=
					fabs(jacdet) * Rho() * GetS0(ip.Point(),i) * GetS0(ip.Point(),j) * ip.Weight();
			}
		}	
	}
}

// FFRF: compute integrals int_V \rho ubar_k * ubar_l dV
double FiniteElement3DFFRF::GetIntRhoUkUl(int k, int l, const Vector& xgloc)
{
	//Compute: int_V \rho ubar_k * ubar_l dV = Ikl + Ibarlk*qf + Ibarkl*qf + qf*Sbarkl*qf
	int dim = Dim();
	int ns = NS();
	int len = ffrfdata->IbarklS[l-1][k-1].Length();

	double v = GetIkl(k, l);   // Ikl

	// temp = IbarklS * xg
	ConstVector<FEmaxDOF> temp;
	temp.SetLen(dim*ns);

	for (int i=1; i<=ns; i++)
	{
		v += ffrfdata->IbarklS[l-1][k-1](dim*(i-1)+k)*xgloc(dim*(i-1)+k);  // Ibarlk*qf
		v += ffrfdata->IbarklS[k-1][l-1](dim*(i-1)+l)*xgloc(dim*(i-1)+l);  // Ibarkl*qf
		for (int j=1; j<=ns; j++)
		{
			v += xgloc(dim*(i-1)+k)*ffrfdata->SbarklS[k-1][l-1](dim*(i-1)+k,dim*(j-1)+l)*xgloc(dim*(j-1)+l);  // qf*Sbarkl*qf
		}
	}

	return v;
}

// FFRF: compute matrix of integrals int_V \rho ubar_k * ubar_l dV
void FiniteElement3DFFRF::GetIntRhoUkUlMat(Matrix3D& mat, const Vector& xgloc)
{
	//Compute: int_V \rho ubar_k * ubar_l dV = Ikl + Ibarlk*qf + Ibarkl*qf + qf*Sbarkl*qf
	int dim = Dim();
	int ns = NS();
	int len = ffrfdata->IbarklS[0][0].Length();


	for (int k=1; k<=3; k++)
		for (int l=1; l<=k; l++)
		{
			mat(k,l) = GetIkl(k, l);   // Ikl
			for (int i=1; i<=ns; i++)
			{
				int ii = (i-1)*dim;
				mat(k,l) += ffrfdata->IbarklS[l-1][k-1](ii+k)*xgloc(ii+k) + ffrfdata->IbarklS[k-1][l-1](ii+l)*xgloc(ii+l);   // Ibarlk*qf
				for (int j=1; j<=ns; j++)
				{
					mat(k,l) += xgloc(ii+k)*ffrfdata->SbarklS[k-1][l-1](ii+k,dim*(j-1)+l)*xgloc(dim*(j-1)+l);  // qf*Sbarkl*qf
				}
			}
		}

		// symmetric part
		mat(1,2) = mat(2,1);
		mat(1,3) = mat(3,1);
		mat(2,3) = mat(3,2);

}

// FFRF: compute integrals d/dt( int_V \rho ubar_k * ubar_l dV )
double FiniteElement3DFFRF::GetIntRhoUkUlP(int k, int l, const Vector& xgloc, const Vector& xglocp)
{
	// Compute: d/dt {int_V \rho ubar_k * ubar_l dV} =  0  + Ibarlk*qf_dot + Ibarkl*qf_dot + qf_dot*Sbarkl*qf + qf*Sbarkl*qf_dot
	double v = 0;
	int len = ffrfdata->IbarklS[l-1][k-1].Length();
	int dim = Dim();
	int ns = NS();

	for (int i=1; i<=ns; i++)
	{
		v += ffrfdata->IbarklS[l-1][k-1](dim*(i-1)+k)*xglocp(dim*(i-1)+k);  // Ibarlk*qf_dot
		v += ffrfdata->IbarklS[k-1][l-1](dim*(i-1)+l)*xglocp(dim*(i-1)+l);  // Ibarkl*qf_dot
		for (int j=1; j<=ns; j++)
		{
			v += xglocp(dim*(i-1)+k)*ffrfdata->SbarklS[k-1][l-1](dim*(i-1)+k,dim*(j-1)+l)*xgloc(dim*(j-1)+l);  // qfdot * Sbarkl * qf
			v += xgloc(dim*(i-1)+k)*ffrfdata->SbarklS[k-1][l-1](dim*(i-1)+k,dim*(j-1)+l)*xglocp(dim*(j-1)+l);  // qf * Sbarkl * qfdot
		}
	}
	return v;
}

// (negative) quadratic velocity vector -Q is added to fadd
void FiniteElement3DFFRF::AddQuadraticVelocityVector(Vector& fadd, double t)
{
	// GCMS: no FFRF-quadratic velocity vector is used
	if (IsGCMS()) return;

	TMStartTimer(17);
	TMStartTimer(24);

	int off = FlexDOF(); //size of M_ff
	Matrix3D Gbar(ReferenceFrame().GetGbar());
	Matrix3D GbarT(ReferenceFrame().GetGbarT());
	Matrix3D GbarPT(ReferenceFrame().GetGbarpT());
	Matrix3D A(ReferenceFrame().GetRotMatrix());

	int ffrf_rot = FFRFRotationDim();
	int ffrf_trans = FFRFTranslationDim();
	int flexdof = FlexDOF();
	int ns = NS();
	int dim = Dim();

	ConstVector<FEmaxDOF> xgloc(SOS());
	for (int i=1; i<=SOS(); i++)
		xgloc(i) = XG(i);
	ConstVector<FEmaxDOF> xglocp(SOS());
	for (int i=1; i<=SOS(); i++)
		xglocp(i) = XGP(i);


	Vector3D omega_bar;
	// omega_bar
	omega_bar = ReferenceFrame().GetAngularVelLocal();
	Matrix3D omega_bartilde, omega_bartilde2;
	omega_bartilde.SetSkew(omega_bar);
	omega_bartilde2 = omega_bartilde*omega_bartilde;

	TMStartTimer(27);
	// Sbar_t =  I1 + Sbar*qf
	Vector3D Sbar_t(ffrfdata->I1S(1),ffrfdata->I1S(2),ffrfdata->I1S(3));
	for (int i=1; i <= 3; i++)
	{
		for (int j=1; j <= ns; j++)
		{
			Sbar_t(i) += ffrfdata->SbarS(i,dim*(j-1)+i)*xgloc(dim*(j-1)+i);
		}
	}

	//TMStartTimer(28);
	//Matrix3D Ibar_theta_theta
	Matrix3D Ibar_theta_theta;
	GetIbarThetaTheta(Ibar_theta_theta, xgloc);

	////Matrix3D Ibar_theta_thetaP
	Matrix3D Ibar_theta_thetaP;
	GetIbarThetaThetaP(Ibar_theta_thetaP, xgloc, xglocp);

	//TMStopTimer(28);
	TMStopTimer(27);

	TMStartTimer(29);
	//Ibar_theta_f = int_V \rho[q_f*Sbar_tilde23, q_f*Sbar_tilde31, q_f*Sbar_tilde12]dV + int_V \rho [x2*S3-x3*S2, x3*S1-x1*S3, x1*S2-x2*S1] dV
	//             = int_V \rho[q_f*Sbar_tilde23, q_f*Sbar_tilde31, q_f*Sbar_tilde12]dV + [Ibarkl23-Ibarkl32, Ibarkl31-Ibarkl13, Ibarkl12-Ibarkl21]
	ConstMatrix<FEmaxDOF*3> Ibar_theta_f;
	GetIbarThetaF(Ibar_theta_f, xgloc);

	TMStopTimer(29);

	TMStopTimer(24);
	TMStartTimer(25);

	// Quadratic velocity vector component for flexible dofs:
	// Q_f
	//get integration points:
	if (1)
	{
		for (IntegrationPointsIterator ip(integrationRuleStiffness); !ip.IsEnd(); ++ip)
		{
			double weightdet = fabs(integrationPointStiffnessMatrixLocalData(ip.GetIndex())->jacdet) * Rho() * ip.Weight();

			Vector3D omega_bartilde2_ubar, omega_bartilde_ufbardot;
			omega_bartilde2_ubar = omega_bartilde2 * GetPosRel(ip.Point());
			omega_bartilde_ufbardot = omega_bartilde * GetVelRel(ip.Point());
			//omega_tilde_ufbardot = omega_tilde2*GetVelRel(p);
			//fadd(i) += Rho() * Shapevec(i) * ( (omega_tilde2) * ubar - 2*omega_tilde*ufbardot ) * weight * det;
			for (int i = 1; i <= ns; i++)
			{
				double sv = GetS0(ip.Point(), i);
				for (int j=1; j <= dim; j++)
				{
					fadd(dim*(i-1)+j) += weightdet * ( sv * (omega_bartilde2_ubar(j) -  2*omega_bartilde_ufbardot(j)) );
				}
			}

		}
	}

	// Quadratic velocity vector component for translation dofs:
	// Q_R
	if (1)
	{
		Vector3D Qr(0.,0.,0.), help(0.,0.,0.);
		help = omega_bartilde2*Sbar_t;
		Vector3D Sbarqfdot(0.,0.,0.);
		for (int i=1; i<=FlexDOF(); i++)
			for (int j=1; j<=3; j++)
				Sbarqfdot(j) += ffrfdata->SbarS(j,i)*XGP(i);

		help += 2*omega_bartilde*Sbarqfdot;
		Qr = -1*A*help;

		for (int i=1; i<=3; i++)
			fadd(FlexDOF()+i) -= Qr(i);
	}


	// Quadratic velocity vector component for rotation dofs:
	// Q_theta = -2GbarP^T Ibar_thth omega - 2 GbarP^T Ibar_thetaf qfP  - Gbar^T Ibar_ththP omega = -2GbarP^T*AA - 2GbarP^T*BB - Gbar^T*CC
		ConstVector<4> Qtheta(4);
		Qtheta.SetAll(0.);
		Vector3D AA, BB(0.,0.,0.), CC;
		AA = Ibar_theta_theta*omega_bar;
		//B = Ibar_theta_f * qfP;
		for(int i=1; i<=FlexDOF(); i++)
		{
			BB(1) += Ibar_theta_f(1,i)*XGP(i);
			BB(2) += Ibar_theta_f(2,i)*XGP(i);
			BB(3) += Ibar_theta_f(3,i)*XGP(i);
		}
		CC = Ibar_theta_thetaP * omega_bar;

		for (int i=1; i<=ffrf_rot; i++)
			for (int j=1; j<=ffrf_trans; j++)
			{
				Qtheta(i) -= 2*GbarPT(i,j)*(AA(j)+BB(j)) + GbarT(i,j) * CC(j);
			}	

			for (int i=1; i<=ffrf_rot; i++)
				fadd(FlexDOF()+ffrf_trans+i) -= Qtheta(i);


	TMStopTimer(25);
	TMStopTimer(17);


}

//Ibar_theta_f = int_V \rho[q_f*Sbar_tilde23, q_f*Sbar_tilde31, q_f*Sbar_tilde12]dV + int_V \rho [x2*S3-x3*S2, x3*S1-x1*S3, x1*S2-x2*S1] dV
//             = int_V \rho[q_f*Sbar_tilde23, q_f*Sbar_tilde31, q_f*Sbar_tilde12]dV + [Ibarkl23-Ibarkl32, Ibarkl31-Ibarkl13, Ibarkl12-Ibarkl21]
void FiniteElement3DFFRF::GetIbarThetaF(Matrix& Ibar_theta_f, const Vector& xgloc)
{
	int flexdof = FlexDOF();
	Ibar_theta_f.SetSize(Dim(), flexdof);
	Ibar_theta_f.SetAll(0);

	// compute temp = Ibar_theta_f 
	// could be implemented more efficiently using that several components of IbarklS, SbarTilde are zero
	for (int i=1; i <= flexdof; i++)
	{
		// first step: temp = [I-term] skript p. 82
		Ibar_theta_f(1,i) = ffrfdata->IbarklS[1][2](i)-ffrfdata->IbarklS[2][1](i);
		Ibar_theta_f(2,i) = ffrfdata->IbarklS[2][0](i)-ffrfdata->IbarklS[0][2](i);
		Ibar_theta_f(3,i) = ffrfdata->IbarklS[0][1](i)-ffrfdata->IbarklS[1][0](i);

		// second step: add q * N-term from skript p. 82 (5.61)
		for (int j=1; j <= flexdof; j++)
		{
			for (int k=1; k <= Dim(); k++)
			{
				Ibar_theta_f(k,i) += xgloc(j)*ffrfdata->SbarTilde[k-1](j,i);
			}
		}
	}
}

void FiniteElement3DFFRF::GetIbarThetaTheta(Matrix3D& Ibar_theta_theta, const Vector& xgloc)
{
	double IntRhoU1U1 = GetIntRhoUkUl(1,1, xgloc);
	double IntRhoU2U2 = GetIntRhoUkUl(2,2, xgloc);
	double IntRhoU3U3 = GetIntRhoUkUl(3,3, xgloc);
	double IntRhoU1U2 = GetIntRhoUkUl(1,2, xgloc);
	double IntRhoU1U3 = GetIntRhoUkUl(1,3, xgloc);
	double IntRhoU2U3 = GetIntRhoUkUl(2,3, xgloc);

	Ibar_theta_theta =
	Matrix3D(  //Shabana page 208, Eq.(5.69)
		IntRhoU2U2+IntRhoU3U3,-IntRhoU1U2,          -IntRhoU1U3,
		-IntRhoU1U2,          IntRhoU1U1+IntRhoU3U3,-IntRhoU2U3,
		-IntRhoU1U3,          -IntRhoU2U3,          IntRhoU1U1+IntRhoU2U2);

}

void FiniteElement3DFFRF::GetIbarThetaThetaP(Matrix3D& Ibar_theta_thetaP, const Vector& xgloc, const Vector& xglocp)
{
	double IntRhoU1U1P = GetIntRhoUkUlP(1,1, xgloc, xglocp);
	double IntRhoU2U2P = GetIntRhoUkUlP(2,2, xgloc, xglocp);
	double IntRhoU3U3P = GetIntRhoUkUlP(3,3, xgloc, xglocp);
	double IntRhoU1U2P = GetIntRhoUkUlP(1,2, xgloc, xglocp);
	double IntRhoU1U3P = GetIntRhoUkUlP(1,3, xgloc, xglocp);
	double IntRhoU2U3P = GetIntRhoUkUlP(2,3, xgloc, xglocp);

	Ibar_theta_thetaP = 
		Matrix3D(  //Shabana page 208, Eq.(5.69)
		IntRhoU2U2P+IntRhoU3U3P,-IntRhoU1U2P,          -IntRhoU1U3P,
		-IntRhoU1U2P,          IntRhoU1U1P+IntRhoU3U3P,-IntRhoU2U3P,
		-IntRhoU1U3P,          -IntRhoU2U3P,          IntRhoU1U1P+IntRhoU2U2P);

}

const double& FiniteElement3DFFRF::GetXact(int i) const 
{
	if (IsCMS())
	{
		const Rigid3D& framerigid = ReferenceFrame();
		const BaseCMSElement<Rigid3D>* cmsrigid = dynamic_cast<const BaseCMSElement<Rigid3D>*>(&framerigid);
		const BaseCMSElement<Rigid3DKardan>* cmsrigidkardan = dynamic_cast<const BaseCMSElement<Rigid3DKardan>*>(&framerigid);
		if (cmsrigid)
		{
			return cmsrigid->GetXactFull(i);
		}
		else
		{
			return cmsrigidkardan->GetXactFull(i);
		}
	}
	else return mbs->GetXact(i);
}

const double& FiniteElement3DFFRF::GetDrawValue(int iloc) const 
{
	if (IsCMS())
	{
		const Rigid3D& framerigid = ReferenceFrame();
		const BaseCMSElement<Rigid3D>* cmsrigid = dynamic_cast<const BaseCMSElement<Rigid3D>*>(&framerigid);
		const BaseCMSElement<Rigid3DKardan>* cmsrigidkardan = dynamic_cast<const BaseCMSElement<Rigid3DKardan>*>(&framerigid);
		if (cmsrigid)
		{
			return cmsrigid->GetDrawValueFull(iloc);
		}
		else
		{
			return cmsrigidkardan->GetDrawValueFull(iloc);
		}
	}
	else return mbs->GetDrawValue(iloc);
}
Vector3D FiniteElement3DFFRF::GetNodePos(int i) const
{
	return GetPos(GetNodeLocPos(i));
}

Vector3D FiniteElement3DFFRF::GetNodePosD(int i) const 
{
	// GetPosD(locpos, 0) -> 0 means no magnification
	return GetPosD(GetNodeLocPos(i), 0);
}

Vector3D FiniteElement3DFFRF::GetNodeVel(int i) const
{
	return GetVel(GetNodeLocPos(i));
}

Vector3D FiniteElement3DFFRF::GetPos(const Vector3D& p_loc) const 
{ 
	if (IsGCMS())
	{
		return GetPosRel(p_loc);
	}
	else
	{
		return RefFramePos() + RefFrameRot() * GetPosRel(p_loc);
	}
}
Vector3D FiniteElement3DFFRF::GetVel(const Vector3D& p_loc) const 
{
	if (IsGCMS())
	{
		return GetVelRel(p_loc);	
	}
	else
	{
		return RefFrameVel() + RefFrameRotP() * GetPosRel(p_loc) + RefFrameRot() * GetVelRel(p_loc);
	}
}
Vector3D FiniteElement3DFFRF::GetPosD(const Vector3D& p_loc, int use_magnification) const
{ 
	if (IsGCMS())
	{
		// AP: GetPosRelD yields the position vector relative to the absolute coordinate system
		//  to obtain magnification only of the flexible and not of the rigid body deformation
		//  take:
		//  (1-factof)*(RefFrameTranslation()+RefFrameRotation()*InitPosition()) + factor*ActualPosition()
		// Faster computation possible in case of use_magnification = 0, then return GetPosRelD(p_loc,0)
		double factor = GetMBS()->GetDOption(105);
		if (!use_magnification) {factor = 1;}
		return (1-factor)*( RefFramePosD()-RefFramePosInit() +RefFrameRotD()*GetRefConfPos(p_loc) )
			+ factor*GetPosRelD(p_loc,0);
	}
	else
	{
		return RefFramePosD() + RefFrameRotD() * GetPosRelD(p_loc, use_magnification);
	}
}
Vector3D FiniteElement3DFFRF::GetVelD(const Vector3D& p_loc) const 
{ 
	if (IsGCMS())
	{
		return GetVelRelD(p_loc); 
	}
	else
	{
		return RefFrameVelD() + RefFrameRotPD() * GetPosRelD(p_loc, 0) + RefFrameRotD() * GetVelRelD(p_loc); 
	}
}
Vector3D FiniteElement3DFFRF::GetDisplacement(const Vector3D& p_loc) const 
{ 
	if (IsGCMS())
	{
		return GetPosRel(p_loc) - GetRefConfPos(p_loc);
	}
	else
	{
		//AP: added  " -GetRefConfPos(p_loc) "
		return RefFramePos() + RefFrameRot() * GetPosRel(p_loc) - (ReferenceFrame().GetRefPosInit()+GetRefConfPos(p_loc)); 
	}
}
Vector3D FiniteElement3DFFRF::GetDisplacementD(const Vector3D& p_loc) const 
{ 
	if (IsGCMS())
	{
		return GetPosRelD(p_loc, 0) - GetRefConfPos(p_loc);
	}
	else
	{
		return RefFramePosD() + RefFrameRotD() * (GetPosRelD(p_loc, 0)) - (ReferenceFrame().GetRefPosInit()+GetRefConfPos(p_loc)); 
	}
}

void FiniteElement3DFFRF::GetdPosdqT(const Vector3D& ploc, Matrix& d)
{
	if (IsGCMS()) // GCMS: dpos/dq as in FiniteElement3D
	{
		FiniteElement3D::GetdPosdqT(ploc, d);
	}
	else // CMS or FFRF: dp_ref/dq + dA/dq p_rel + A dp_rel/dq
	{
	d.SetSize(SOS(),Dim());
	d.FillWithZeros();
	Matrix3D A = ReferenceFrame().GetRotMatrix();
	Vector3D prel;
	prel = GetPosRel(ploc);
	int j=1;
	for (int i = 1; i <= NS(); i++)
	{
		double sv =	GetS0(ploc, i);
		Matrix3D svA(A);
		svA *= sv;

		d(j,1) = svA(1,1); 
		d(j,2) = svA(2,1);
		d(j,3) = svA(3,1);
		j++;
		d(j,1) = svA(1,2);
		d(j,2) = svA(2,2);
		d(j,3) = svA(3,2);
		j++;
		d(j,1) = svA(1,3);
		d(j,2) = svA(2,3);
		d(j,3) = svA(3,3);
		j++;
	}

	//d_pref/dq:
	d(FlexDOF()+1,1) = 1;
	d(FlexDOF()+2,2) = 1;
	d(FlexDOF()+3,3) = 1;

	//d_A/dq*prel:
	Matrix3D G;
	ReferenceFrame().GetH3T(prel, G);
	for (int i=1; i<=FFRFRotationDim(); i++)
		for (int j=1; j<=3; j++)
			d(FlexDOF()+FFRFTranslationDim()+i,j) = G(i,j);
	}
	//UO() << "dpdq=" << d << "\n";
}

void FiniteElement3DFFRF::GetIntDuDq(Matrix& H)
{
	// GCMS: absolute coordinates -> IntDuDq is the same as for standard FiniteElement3D 
	if (IsGCMS())
	{
		FiniteElement3D::GetIntDuDq(H);
		return;
	}

	// FFRF/CMS: IntDuDq contains terms from floating frame
	int dim = Dim();
	int ns = NS();
	int sos = SOS();

	H.SetSize(sos,dim);
	H.SetAll(0);
	ConstMatrix<(FEmaxDOF+FFRFsize)*3> dudq;
	Matrix3D jac;
	// flexible part Hf
	for (IntegrationPointsIterator ip(integrationRuleLoad); !ip.IsEnd(); ++ip)
	{
		GetJacobi(jac, ip.Point());
		double jacdet = jac.Det();
		double fact = fabs (jacdet) * ip.Weight();
		GetdPosdqT(ip.Point(), dudq);
		dudq *= fact;

		H += dudq;
	}
	Hmatrix = H;
}

void FiniteElement3DFFRF::EvalM(Matrix& m, double t) 
{
	// Initialize mass matrix
	m.SetSize(SOS(), SOS());
	m.SetAll(0.);

	// GCMS: mass matrix is constant matrix M_ff
	if (IsGCMS())
	{
		EvalMff(m, t);
		return;
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//for FFRF mode!

	//sorting of mass matrix:
	// ( M_ff     M_rf     M_thetaf    )
	// ( M_rf     M_rr     M_thetar    )
	// ( M_thetaf M_thetar M_thetatheta)

	//order of DOF in LTG: 0 .. FlexDOF() x_ref y_ref z_ref theta_0 theta_1 theta_2 theta_3
	TMStartTimer(23);

	int off = FlexDOF(); //size of M_ff
	// matrix Gbar, omega = Gbar theta
	Matrix3D Gbar(ReferenceFrame().GetGbar());
	Matrix3D GbarT(ReferenceFrame().GetGbarT());
	// rotation matrix
	Matrix3D A(ReferenceFrame().GetRotMatrix());

	int ffrf_rot = FFRFRotationDim();
	int ffrf_trans = FFRFTranslationDim();
	int ns = NS();
	int dim = Dim();
	int flexdof = FlexDOF();
	int sos = SOS();

	ConstVector<FEmaxDOF> xgloc(sos), xglocp(sos);
	for (int i=1; i<=sos; i++)
	{
		xgloc(i) = XG(i);
		xglocp(i) = XGP(i);
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//m_ff

	// fill only FlexibleDof entries of massmatrix, other entries remain unchanged
	EvalMff(m, t);


	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//mRR = UnitMatrix(3,3) * mass
	m(1+off,1+off) = GetMass();
	m(2+off,2+off) = GetMass();
	m(3+off,3+off) = GetMass();

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//m_Rtheta = -A * SbarTilde_t * Gbar
	//SbarTilde_t = Skew(Sbar_t)
	//Sbar_t = int_V \rho*[ubar_0 + ubar_f] dV = I1 + Sbar*qf

	// Sbar_t =  I1 + Sbar*qf
	Vector3D Sbar_t(ffrfdata->I1S(1),ffrfdata->I1S(2),ffrfdata->I1S(3));
	for (int i=1; i <= 3; i++)
	{
		for (int j=1; j <= ns; j++)
		{
			Sbar_t(i) += ffrfdata->SbarS(i,dim*(j-1)+i)*xgloc(dim*(j-1)+i);
		}
	}

	// SbarTilde_t = Skew(Sbar_t)
	Matrix3D SbarTilde_t;
	SbarTilde_t.SetSkew(Sbar_t);

	// mass matrix component
	Matrix3D mRtheta, help;
	help = SbarTilde_t * Gbar;
	mRtheta = -1.*(A * help);

	for (int i=1; i <= ffrf_trans; i++) //position coordinates
	{
		for (int j=1; j <= ffrf_rot; j++) //rotational coordinates
		{
			m(i+off, ffrf_trans+j+off) = mRtheta(i, j);
			m(ffrf_trans+j+off, i+off) = mRtheta(i, j);
		}
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//m_Rf = A*Sbar //very important term
	for (int j=1; j <= FlexDOF(); j++)
	{
		Vector3D v(ffrfdata->SbarS(1,j),ffrfdata->SbarS(2,j),ffrfdata->SbarS(3,j));
		v = A*v;
		m(j, 1+off) = v(1);
		m(j, 2+off) = v(2);
		m(j, 3+off) = v(3);
		m(1+off, j) = v(1);
		m(2+off, j) = v(2);
		m(3+off, j) = v(3);
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//m_theta_theta = GbarT * Ibar_theta_theta * Gbar, Eq. (5.67)

	Matrix3D Ibar_theta_theta;
	GetIbarThetaTheta(Ibar_theta_theta, xgloc);


	Matrix3D help2 = Ibar_theta_theta*Gbar;
	Matrix3D m_theta_theta = GbarT*(help2);

	for (int i=1; i <= ffrf_rot; i++)
	{
		for (int j=1; j <= ffrf_rot; j++)
		{
			m(i+off+ffrf_trans, j+off+ffrf_trans) = m_theta_theta(i, j);
		}
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//m_theta_f = -Gbar^T*Ibar_theta_f
	//Ibar_theta_f = int_V \rho[q_f*Sbar_tilde23, q_f*Sbar_tilde31, q_f*Sbar_tilde12]dV + int_V \rho [x2*S3-x3*S2, x3*S1-x1*S3, x1*S2-x2*S1] dV
	//             = int_V \rho[q_f*Sbar_tilde23, q_f*Sbar_tilde31, q_f*Sbar_tilde12]dV + [Ibarkl23-Ibarkl32, Ibarkl31-Ibarkl13, Ibarkl12-Ibarkl21]
	ConstMatrix<3*FEmaxDOF> IbarThetaF;
	GetIbarThetaF(IbarThetaF, xgloc);


	for (int i=1; i <= flexdof; i++)
	{
		Matrix3D v;
		v.SetSize(3,1);
		v(1,1) = IbarThetaF(1,i);
		v(2,1) = IbarThetaF(2,i);
		v(3,1) = IbarThetaF(3,i);

		//v = GbarT*v;

		for (int j=1; j <= ffrf_rot; j++)
			for (int k=1; k<=Dim(); k++)
			{
				m(i, off+ffrf_trans+j) += GbarT(j,k)*v(k,1); // = v(j,1); // 
				m(off+ffrf_trans+j, i) += GbarT(j,k)*v(k,1); // = v(j,1); // 
			}
	}

	TMStopTimer(23);
}

int FiniteElement3DFFRF::AddBodyNode(Node & n)
{
	if(IsCMS())
		return ReferenceFrame().AddNode(&n);
	else
		return mbs->AddBodyNode(&n);
}


