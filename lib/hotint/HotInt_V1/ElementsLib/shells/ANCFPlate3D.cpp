//#**************************************************************
//#
//# filename:             ANCFPlate3D.cpp
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
 
#include "body3D.h"
#include "femathhelperfunctions.h"
//#include "material.h"
#include "ANCFPlate3D.h"
#include "Node.h"
//#include "graphicsconstants.h"
//#include "elementdata.h"
//#include "elementdataaccess.h"
//#include "solversettings_auto.h"


void ANCFPlate3D::LinkToElements()
{
	if (SOSowned() != SOS())
	{
		const Node& node1 = GetMBS()->GetNode(node[0]);
		const Node& node2 = GetMBS()->GetNode(node[1]);
		const Node& node3 = GetMBS()->GetNode(node[2]);
		const Node& node4 = GetMBS()->GetNode(node[3]);

		int Nnonode = SOSowned();
		int Nnode = SOS()-Nnonode;
		TArray<int> storeltg(Nnonode);
		for (int i=1; i <= Nnonode*2; i++)
		{
			storeltg.Add(LTG(i));
		}
		LTGreset();

		for (int i=1; i <= node1.SOS(); i++)
		{
			AddLTG(node1.Get(i));
		}
		for (int i=1; i <= node2.SOS(); i++)
		{
			AddLTG(node2.Get(i));
		}
		for (int i=1; i <= node3.SOS(); i++)
		{
			AddLTG(node3.Get(i));
		}
		for (int i=1; i <= node4.SOS(); i++)
		{
			AddLTG(node4.Get(i));
		}
		for (int i=1; i <= Nnonode; i++)
		{
			AddLTG(storeltg(i));
		}

		for (int i=1; i <= node1.SOS(); i++)
		{
			AddLTG(node1.Get(i+node1.SOS()));
		}
		for (int i=1; i <= node2.SOS(); i++)
		{
			AddLTG(node2.Get(i+node2.SOS()));
		}
		for (int i=1; i <= node3.SOS(); i++)
		{
			AddLTG(node3.Get(i+node3.SOS()));
		}
		for (int i=1; i <= node4.SOS(); i++)
		{
			AddLTG(node4.Get(i+node4.SOS()));
		}
		for (int i=1; i <= Nnonode; i++)
		{
			AddLTG(storeltg(i+Nnonode));
		}
		/*
		UO() << "ltg=[";
		for (int i=1; i <= LTGlength(); i++)
		{
			UO() << LTG(i) << " ";
		}
		UO() << "]\n";
		*/		
	}
}


void ANCFPlate3D::BuildDSMatrices() 
{
	orderxy = 9; //9x5 und 7x3 ergibt fast das gleiche????
	orderz = 5; 


	GetIntegrationRule(x1,w1,orderxy);
	GetIntegrationRule(x2,w2,orderxy);
	GetIntegrationRule(x3,w3,orderz);

	Matrix3D jac0, jacinv;
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

				GetJacobi(jac0,p,DS,x_init);
				jacdet[ind] = jac0.Det();
				jac0.GetInverse(jacinv);
				jacinv = jacinv.GetTp();

				grad[ind].Init();
				grad[ind].SetSize(Dim(),NS());
				Mult(jacinv, DS, grad[ind]);
			}
		}
	}
	//Build T-matrices:

	Vector3D p[4];
	p[0] = Vector3D(-1,-1,0);
	p[1] = Vector3D( 1,-1,0);
	p[2] = Vector3D( 1, 1,0);
	p[3] = Vector3D(-1, 1,0);


	for (int i=0; i < NNodes(); i++)
	{
		Ti[i].Init();
		Ti[i].SetSize(NodeSize(),NodeSize());
		Ti[i].SetAll(0);

		/*
		GetDSMatrix0(p[i],DS);
		GetJacobi(jac[i],p[i],DS,x_init);
		jac[i] = jac[i].GetTp();
		*/
	}

	for (int k = 0; k < NNodes(); k++)
	{
		for (int i=1; i <= Dim(); i++)
		{
			for (int j=1; j <= Dim(); j++)
			{
				(jac[k])(i,j) = x_init((i-1)*Dim()+j+Dim()+k*NodeSize());
			}
		}
		Ti[k](1,1) = 1; Ti[k](2,2) = 1; Ti[k](3,3) = 1;
	}

	T.Init();
	T.SetSize(SOS(),SOS());
	T.SetAll(0);
	for (int i = NNodes()*NodeSize()+1; i <= SOS(); i++)
	{
		T(i,i) = 1;
	}

	for (int k = 0; k < NNodes(); k++)
	{
		for (int i=1; i <= Dim(); i++)
		{
			for (int j=1; j <= Dim(); j++)
			{
				Ti[k](i*Dim()+1,j*Dim()+1) = jac[k](i,j);
				Ti[k](i*Dim()+2,j*Dim()+2) = jac[k](i,j);
				Ti[k](i*Dim()+3,j*Dim()+3) = jac[k](i,j);
			}
		}
		T.SetSubmatrix(Ti[k],1+k*NodeSize(),1+k*NodeSize());
	}

/*
	UO() << "T1=" << Ti[0] << "\n";
	UO() << "T2=" << Ti[1] << "\n";
	UO() << "T3=" << Ti[2] << "\n";
	UO() << "T4=" << Ti[3] << "\n";

	UO() << "T=" << T << "\n";
	UO() << "jac1=" << jac[0] << "\n";
	UO() << "jac2=" << jac[1] << "\n";
	UO() << "jac3=" << jac[2] << "\n";
	UO() << "jac4=" << jac[3] << "\n";
*/
};

void ANCFPlate3D::ApplyT(Vector& v) const
{

	static Vector x;
	x = v;
	for (int k = 0; k < NNodes(); k++)
	{
		for (int i=1; i<=3; i++)
		{
			for (int j=1; j<=3; j++)
			{
				v((i-1)*3+j+3+k*12) = jac[k](i,1)*x(0+j+3+k*12) + jac[k](i,2)*x(3+j+3+k*12) + jac[k](i,3)*x(6+j+3+k*12);
			}
		}
	}

	//global_uo << "v1=" << v << "\n";
	//v = T*x;
	//global_uo << "v2=" << v << "\n";

}

void ANCFPlate3D::ApplyTtp(Vector& v) const 
{
	//v = T.GetTp()*v; return;
	static Vector x;
	x = v;
	for (int k = 0; k < NNodes(); k++)
	{
		for (int i=1; i<=3; i++)
		{
			for (int j=1; j<=3; j++)
			{
				v((i-1)*3+j+3+k*12) = jac[k](1,i)*x(0+j+3+k*12) + jac[k](2,i)*x(3+j+3+k*12) + jac[k](3,i)*x(6+j+3+k*12);
			}
		}
	}
}

void ANCFPlate3D::ApplyTtp(Matrix& m) const
{
	//m = T.GetTp()*m; return;

	static Vector x;
	int sns = SOS();
	x.SetLen(sns);

	for (int l = 0; l < NNodes(); l++)
	{
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
					m((i-1)*3+j+3+l*12,k) = jac[l](1,i)*x(0+j+3+l*12) + jac[l](2,i)*x(3+j+3+l*12) + jac[l](3,i)*x(6+j+3+l*12);
				}
			}
		}
	}
}

int ANCFPlate3D::NS() const 
{
	return 16;
}
//double t0;

void ANCFPlate3D::GetS0(Vector& sf, const Vector3D& ploc) const
{
	double xb = ploc.X();
	double yb = ploc.Y();
	double zb = ploc.Z();
	sf.SetLen(NS());

	double xb2 = xb*xb;
	double xb3 = xb2*xb;
	double yb2 = yb*yb;
	double yb3 = yb2*yb;

	sf(1) = -xb*yb3/8.0-xb3*yb/8.0+1.0/4.0+yb3/8.0+xb*yb/2.0-3.0/8.0*xb-3.0/8.0*yb+xb3/8.0;
	sf(2) = -lx*xb3*yb/16.0+lx/16.0-lx*xb/16.0-lx*xb2/16.0+lx*xb3/16.0+lx*xb2*yb/16.0-lx*yb/16.0+lx*xb*yb/16.0;
	sf(3) = -ly*yb2/16.0+ly*yb3/16.0-ly*xb/16.0-ly*yb/16.0+ly*xb*yb2/16.0-ly*xb*yb3/16.0+ly*xb*yb/16.0+ly/16.0;
	sf(4) = lz*(-xb-yb+1.0+xb*yb)*zb/8.0;
	sf(5) = xb*yb3/8.0+xb3*yb/8.0+1.0/4.0+yb3/8.0-xb*yb/2.0+3.0/8.0*xb-3.0/8.0*yb-xb3/8.0;
	sf(6) = -lx/16.0-lx*xb/16.0+lx*xb2/16.0+lx*xb3/16.0+lx*yb/16.0+lx*xb*yb/16.0-lx*xb2*yb/16.0-lx*xb3*yb/16.0;
	sf(7) = ly/16.0-ly*yb/16.0+ly*xb/16.0-ly*xb*yb/16.0-ly*yb2/16.0-ly*xb*yb*yb/16.0+ly*yb3/16.0+ly*xb*yb3/16.0;
	sf(8) = -lz*(1.0+xb)*(-1.0+yb)*zb/8.0;
	sf(9) = -xb*yb3/8.0-xb3*yb/8.0+1.0/4.0-yb3/8.0+xb*yb/2.0+3.0/8.0*xb+3.0/8.0*yb-xb3/8.0;
	sf(10) = -lx/16.0-lx*yb/16.0-lx*xb/16.0-lx*xb*yb/16.0+lx*xb2/16.0+lx*xb2*yb/16.0+lx*xb3/16.0+lx*xb3*yb/16.0;
	sf(11) = -ly/16.0-ly*yb/16.0+ly*yb2/16.0-ly*xb/16.0-ly*xb*yb/16.0+ly*xb*yb*yb/16.0+ly*yb3/16.0+ly*xb*yb3/16.0;
	sf(12) = lz*(1.0+xb)*(1.0+yb)*zb/8.0;
	sf(13) = xb*yb3/8.0+xb3*yb/8.0+1.0/4.0-yb3/8.0-xb*yb/2.0-3.0/8.0*xb+3.0/8.0*yb+xb3/8.0;
	sf(14) = lx/16.0+lx*yb/16.0-lx*xb/16.0-lx*xb*yb/16.0-lx*xb2/16.0-lx*xb2*yb/16.0+lx*xb3/16.0+lx*xb3*yb/16.0;
	sf(15) = -ly/16.0-ly*yb/16.0+ly*yb2/16.0+ly*yb3/16.0+ly*xb/16.0+ly*xb*yb/16.0-ly*xb*yb2/16.0-ly*xb*yb3/16.0;
	sf(16) = -lz*(1.0+yb)*(-1.0+xb)*zb/8.0;
}

void ANCFPlate3D::GetDSMatrix0(const Vector3D& ploc, Matrix& sf) const
{
	double xb = ploc.X();
	double yb = ploc.Y();
	double zb = ploc.Z();
	sf.SetSize(Dim(),NS());

	double xb2 = xb*xb;
	double xb3 = xb2*xb;
	double yb2 = yb*yb;
	double yb3 = yb2*yb;

	sf(1,1) = -yb3/8.0-3.0/8.0*xb2*yb+yb/2.0-3.0/8.0+3.0/8.0*xb2;
	sf(1,2) = -3.0/16.0*lx*xb2*yb-lx/16.0-lx*xb/8.0+3.0/16.0*lx*xb2+lx*xb*yb/8.0+lx*yb/16.0;
	sf(1,3) = -ly/16.0+ly*yb2/16.0-ly*yb3/16.0+ly*yb/16.0;
	sf(1,4) = lz*(-1.0+yb)*zb/8.0;
	sf(1,5) = yb3/8.0+3.0/8.0*xb2*yb-yb/2.0+3.0/8.0-3.0/8.0*xb2;
	sf(1,6) = -lx/16.0+lx*xb/8.0+3.0/16.0*lx*xb2+lx*yb/16.0-lx*xb*yb/8.0-3.0/16.0*lx*xb2*yb;
	sf(1,7) = ly/16.0-ly*yb/16.0-ly*yb2/16.0+ly*yb3/16.0;
	sf(1,8) = -lz*(-1.0+yb)*zb/8.0;
	sf(1,9) = -yb3/8.0-3.0/8.0*xb2*yb+yb/2.0+3.0/8.0-3.0/8.0*xb2;
	sf(1,10) = -lx/16.0-lx*yb/16.0+lx*xb/8.0+lx*xb*yb/8.0+3.0/16.0*lx*xb2+3.0/16.0*lx*xb2*yb;
	sf(1,11) = -ly/16.0-ly*yb/16.0+ly*yb2/16.0+ly*yb3/16.0;
	sf(1,12) = lz*(1.0+yb)*zb/8.0;
	sf(1,13) = yb3/8.0+3.0/8.0*xb2*yb-yb/2.0-3.0/8.0+3.0/8.0*xb2;
	sf(1,14) = -lx/16.0-lx*yb/16.0-lx*xb/8.0-lx*xb*yb/8.0+3.0/16.0*lx*xb2+3.0/16.0*lx*xb2*yb;
	sf(1,15) = ly/16.0+ly*yb/16.0-ly*yb2/16.0-ly*yb3/16.0;
	sf(1,16) = -lz*(1.0+yb)*zb/8.0;
	sf(2,1) = -3.0/8.0*xb*yb2-xb3/8.0+3.0/8.0*yb2+xb/2.0-3.0/8.0;
	sf(2,2) = -lx*xb3/16.0+lx*xb2/16.0-lx/16.0+lx*xb/16.0;
	sf(2,3) = -ly*yb/8.0+3.0/16.0*ly*yb2-ly/16.0+ly*xb*yb/8.0-3.0/16.0*ly*xb*yb*yb+ly*xb/16.0;
	sf(2,4) = lz*(-1.0+xb)*zb/8.0;
	sf(2,5) = 3.0/8.0*xb*yb2+xb3/8.0+3.0/8.0*yb2-xb/2.0-3.0/8.0;
	sf(2,6) = lx/16.0+lx*xb/16.0-lx*xb2/16.0-lx*xb3/16.0;
	sf(2,7) = -ly/16.0-ly*xb/16.0-ly*yb/8.0-ly*xb*yb/8.0+3.0/16.0*ly*yb2+3.0/16.0*ly*xb*yb2;
	sf(2,8) = -lz*(1.0+xb)*zb/8.0;
	sf(2,9) = -3.0/8.0*xb*yb2-xb3/8.0-3.0/8.0*yb2+xb/2.0+3.0/8.0;
	sf(2,10) = -lx/16.0-lx*xb/16.0+lx*xb2/16.0+lx*xb3/16.0;
	sf(2,11) = -ly/16.0+ly*yb/8.0-ly*xb/16.0+ly*xb*yb/8.0+3.0/16.0*ly*yb2+3.0/16.0*ly*xb*yb2;
	sf(2,12) = lz*(1.0+xb)*zb/8.0;
	sf(2,13) = 3.0/8.0*xb*yb2+xb3/8.0-3.0/8.0*yb2-xb/2.0+3.0/8.0;
	sf(2,14) = lx/16.0-lx*xb/16.0-lx*xb2/16.0+lx*xb3/16.0;
	sf(2,15) = -ly/16.0+ly*yb/8.0+3.0/16.0*ly*yb2+ly*xb/16.0-ly*xb*yb/8.0-3.0/16.0*ly*xb*yb2;
	sf(2,16) = -lz*(-1.0+xb)*zb/8.0;
	sf(3,1) = 0.0;
	sf(3,2) = 0.0;
	sf(3,3) = 0.0;
	sf(3,4) = lz*(-xb-yb+1.0+xb*yb)/8.0;
	sf(3,5) = 0.0;
	sf(3,6) = 0.0;
	sf(3,7) = 0.0;
	sf(3,8) = -lz*(1.0+xb)*(-1.0+yb)/8.0;
	sf(3,9) = 0.0;
	sf(3,10) = 0.0;
	sf(3,11) = 0.0;
	sf(3,12) = lz*(1.0+xb)*(1.0+yb)/8.0;
	sf(3,13) = 0.0;
	sf(3,14) = 0.0;
	sf(3,15) = 0.0;
	sf(3,16) = -lz*(1.0+yb)*(-1.0+xb)/8.0;

}


void ANCFPlate3D::GetS0x(Vector& sfx, const Vector3D& ploc) const
{
	double xb = ploc.X();
	double yb = ploc.Y();
	double zb = ploc.Z();
	sfx.SetLen(NS());
	double f = 2./lx;

	double xb2 = xb*xb;
	double yb2 = yb*yb;
	double yb3 = yb*yb2;

	sfx(1) = f*(-yb3/8.0-3.0/8.0*xb2*yb+yb/2.0-3.0/8.0+3.0/8.0*xb2);
	sfx(2) = f*(-3.0/16.0*lx*xb2*yb-lx/16.0-lx*xb/8.0+3.0/16.0*lx*xb2+lx*xb*yb/8.0+lx*yb/16.0);
	sfx(3) = f*(-ly/16.0+ly*yb2/16.0-ly*yb3/16.0+ly*yb/16.0);
	sfx(4) = f*(lz*(-1.0+yb)*zb/8.0);
	sfx(5) = f*(yb3/8.0+3.0/8.0*xb2*yb-yb/2.0+3.0/8.0-3.0/8.0*xb2);
	sfx(6) = f*(-lx/16.0+lx*xb/8.0+3.0/16.0*lx*xb2+lx*yb/16.0-lx*xb*yb/8.0-3.0/16.0*lx*xb2*yb);
	sfx(7) = f*(ly/16.0-ly*yb/16.0-ly*yb2/16.0+ly*yb3/16.0);
	sfx(8) = f*(-lz*(-1.0+yb)*zb/8.0);
	sfx(9) = f*(-yb3/8.0-3.0/8.0*xb2*yb+yb/2.0+3.0/8.0-3.0/8.0*xb2);
	sfx(10) = f*(-lx/16.0-lx*yb/16.0+lx*xb/8.0+lx*xb*yb/8.0+3.0/16.0*lx*xb2+3.0/16.0*lx*xb2*yb);
	sfx(11) = f*(-ly/16.0-ly*yb/16.0+ly*yb2/16.0+ly*yb3/16.0);
	sfx(12) = f*(lz*(1.0+yb)*zb/8.0);
	sfx(13) = f*(yb3/8.0+3.0/8.0*xb2*yb-yb/2.0-3.0/8.0+3.0/8.0*xb2);
	sfx(14) = f*(-lx/16.0-lx*yb/16.0-lx*xb/8.0-lx*xb*yb/8.0+3.0/16.0*lx*xb2+3.0/16.0*lx*xb2*yb);
	sfx(15) = f*(ly/16.0+ly*yb/16.0-ly*yb2/16.0-ly*yb3/16.0);
	sfx(16) = f*(-lz*(1.0+yb)*zb/8.0);

}

//from -1 to +1
Vector3D ANCFPlate3D::GetPosx(const Vector3D& p_loc, const Vector& xg) const
{
	static Vector SVx;
	GetS0x(SVx, p_loc);
	Vector3D p(0.,0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += SVx(j)*xg((j-1)*Dim()+i);
		}
	}
	return p;
};

void ANCFPlate3D::GetRot(Matrix3D& rot, const Vector& xg) const
{

	int n2 = 3;
	int nn2 = (n2-1)*12;
	Vector3D r1(xg(nn2+1)-xg(1),xg(nn2+2)-xg(2),xg(nn2+3)-xg(3)); //r_B - r_A
	Vector3D r2(0.5*(xg(nn2+7)+xg(7)),0.5*(xg(nn2+8)+xg(8)),0.5*(xg(nn2+9)+xg(9))); //0.5*(r,y_B + r,y_A)
	Vector3D r3;

	r1.Normalize();
	double h = r2*r1;
	r2 -= h*r1;
	r2.Normalize();
	r3 = r1.Cross(r2);

	rot.Set(r1,r2,r3);
}


Vector3D ANCFPlate3D::GetPos(const Vector3D& p_loc) const
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

	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += SV(j)*xg((j-1)*Dim()+i);
		}
	}
	return p;
};

//-1 .. +1, provide xg
Vector3D ANCFPlate3D::GetPos(const Vector3D& p_loc, const Vector& xg) const
{
	static Vector SV;
	GetS0(SV, p_loc);
	Vector3D p(0.,0.,0.);

	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += SV(j)*xg((j-1)*Dim()+i);
		}
	}
	return p;
};

Vector3D ANCFPlate3D::GetVel(const Vector3D& p_loc) const
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
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += SV(j)*xg((j-1)*Dim()+i);
		}
	}
	return p;
};

Vector3D ANCFPlate3D::GetPosD(const Vector3D& p_loc) const
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
Vector3D ANCFPlate3D::GetPos0D(const Vector3D& p_loc) const
{
	static Vector xgd;
	xgd.SetLen(SOS());
	GetDrawCoordinates(xgd);
	xgd = T*xgd;
	static Vector SV;
	GetS0(SV, p_loc);
	Vector3D p(0.,0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += SV(j)*xgd((j-1)*Dim()+i);
		}
	}
	return p;
};

Vector3D ANCFPlate3D::GetVelD(const Vector3D& p_loc) const
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
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += SV(j)*xgd((j-1)*Dim()+i);
		}
	}
	return p;
};


void ANCFPlate3D::GetH(Matrix& H) 
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

		GetIntegrationRule(x1,w1,5);
		GetIntegrationRule(x2,w2,5);
		GetIntegrationRule(x3,w3,2);

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

void ANCFPlate3D::EvalM(Matrix& m, double t) 
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

		Vector x1,x2,x3,w1,w2,w3;


		GetIntegrationRule(x1,w1,5); //optimal: 6x3x3, 6x1x1 geht auch!!!!
		GetIntegrationRule(x2,w2,5);
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
					m += fabs (jacdet) * rho * w1(i1)*w2(i2)*w3(i3) * (HL*HL.GetTp());
				}
			}
		}

		for (int k=0; k < NNodes(); k++)
		{
			m(1+k*12,1+k*12) += concentratedmass[k];
			m(2+k*12,2+k*12) += concentratedmass[k];
			m(3+k*12,3+k*12) += concentratedmass[k];
		}
		m = T.GetTp()*m*T;
		massmatrix = m;

		//UO() << "M=\n";
		//PrintMatrix01(massmatrix);
		//UO() << "massmatrix=" << m << "\n";
	}
	//Matrix minv = m;
	//UO() << "Mass matrix invertable=" << minv.Invert2() << "\n";
	//UO() << "m=" << m << "\n";
};


void ANCFPlate3D::EvalF2(Vector& f, double t) 
{
	Body3D::EvalF2(f,t);
	TMStartTimer(22);

	int sos = SOS();
	SetComputeCoordinates();
	static Vector u;
	static Vector fadd;
	u.SetLen(sos);

	//UO() << "xg=" << xg << "\n";

	for (int i=1; i <= sos; i++)
		u(i) = xg(i) - x_init(i);

	ApplyT(u);
	//UO() << "u=" << u << "\n";

	int ns = NS();
	int dim = Dim();
	//double u;

	double la=Em * nu / ((1.+nu)*(1.-2.*nu));
	double mu=Em / 2. / (1.+nu);

	Matrix3D strain, piola1, F;

	temp.SetLen(SOS());
	fadd.SetLen(SOS());
	fadd.SetAll(0);
	//UO() << "orderx=" << orderx << ", orderyz=" << orderyz << "\n";


	GetIntegrationRule(x1,w1,orderxy);
	GetIntegrationRule(x2,w2,orderxy);
	GetIntegrationRule(x3,w3,orderz);

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
				const Matrix& agrad = grad[ind];
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
					F(j,j) += 1;
				}

				//jacinv = jacinv.GetTp();
				//G = G * jacinv; //also possible!!!

				//F = Matrix3D(1) + G;

				// Green-Lagrange strain tensor
				//strain = 0.5 * (F.GetTp() * F);
				F.GetATA2(strain);
				strain(1,1) -= 0.5; strain(2,2) -= 0.5; strain(3,3) -= 0.5;

				//strain = 0.5 * (G.GetTp() + G + G.GetTp()* G);
				//TMStopTimer(20);
				//TMStartTimer(21);
				piola1 = F * ((2*mu) * strain + Matrix3D(la * strain.Trace()));

				//linear:
				//Matrix3D G=F-Matrix3D(1);
				//strain = 0.5*(G+G.GetTp());
				//piola1 = ((2*mu) * strain + Matrix3D(la * strain.Trace()));

				for (int j=1; j <= Dim(); j++)
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
			//Mult(massmatrix,xg,temp);
			for (int i = 1; i <= SOS(); i++)
				temp(i) = xg(i)*massmatrix(i,i);
		}
		else
		{
			Matrix dmat;
			EvalM(dmat,t);
			//Mult(dmat,xg,temp);
			for (int i = 1; i <= SOS(); i++)
				temp(i)=xg(i)*dmat(i,i);
		}
		temp *= GetMassDamping();
		f -= temp;
	}


}

void ANCFPlate3D::GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables)
{
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_stress,
		FieldVariableDescriptor::FVCI_z, true);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_stress_mises);
}

double ANCFPlate3D::GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & local_position, bool flagD)
{
	if(!flagD)
		return FIELD_VARIABLE_NO_VALUE;

	static Vector u1;
	u1.SetLen(SOS());
	GetDrawCoordinates(u1);
	static Vector u;
	u.SetLen(SOS());

	//UO() << "ploc=" << ploc << "\n";

	//u = T*u;
	Mult(T,u1,u);
	u -= e0;
	Matrix3D F;
	GraduD(local_position,u,F);
	F += Matrix3D(1);
	Matrix3D strain = 0.5*(F.GetTp()*F-Matrix3D(1));
	//F.GetATA2(strain);
	//strain-=Matrix3D(0.5);
	//F -= Matrix3D(1);

	double la=Em * nu / ((1.+nu)*(1.-2.*nu));
	double mu=Em / 2. / (1.+nu);
	Matrix3D stress = (2*mu) * strain + Matrix3D(la * strain.Trace());

	if(fvd.VariableType() == FieldVariableDescriptor::FVT_stress_mises)
		return stress.Mises();
	if(fvd.VariableType() == FieldVariableDescriptor::FVT_stress)
		return fvd.GetComponent(stress);

	//UO() << " Mises=" << val << ", c1=" << c1 << "\n";
	//val = strain(1,1);

	//UO() << "ploc=" << ploc << ", strain=" << stress;
	return FIELD_VARIABLE_NO_VALUE;
}

void ANCFPlate3D::DrawElement() 
{
	mbs->SetColor(col);
	double lx1 = 1; double ly1 = 1; double lz1 = 1*GetMBS()->GetMagnifyYZ();

	int chladni = GetMBS()->GetIOption(116); //for chladni figures in eigenmodes
	int drawflat = GetMBS()->GetIOption(117); //draw only element midplane

	int linemode = 1; //0=no lines, 1=outline+color, 2=outline, 3=elementline+color, 4=elementline
	if (GetMBS()->GetIOption(110) && !GetMBS()->GetIOption(111))
	{
		linemode = 2;
	} 
	else if (!GetMBS()->GetIOption(110) && GetMBS()->GetIOption(111))
	{
		linemode = 0;
	}
	else if (!GetMBS()->GetIOption(110) && !GetMBS()->GetIOption(111))
	{
		linemode = 4;
	}

	static Vector xgd;

	//double lx = -0.5*size.X(); double ly = 0.5*size.Y(); double lz = 0.5*size.Z();

	int colormode = 0;
	if (GetMBS()->GetActualPostProcessingFieldVariable() != NULL)
		colormode = 1;

	//if (colormode)
	{
		double sh = 0.999;
		Vector3D p1(GetPos0D(Vector3D(-lx1,-ly1*sh,-lz1*sh))); 
		Vector3D p2(GetPos0D(Vector3D(-lx1, ly1*sh,-lz1*sh)));
		Vector3D p3(GetPos0D(Vector3D(-lx1, ly1*sh, lz1*sh)));
		Vector3D p4(GetPos0D(Vector3D(-lx1,-ly1*sh, lz1*sh)));

		Vector3D p5=(GetPos0D(Vector3D(lx1,-ly1*sh,-lz1*sh))); 
		Vector3D p6=(GetPos0D(Vector3D(lx1, ly1*sh,-lz1*sh)));
		Vector3D p7=(GetPos0D(Vector3D(lx1, ly1*sh, lz1*sh)));
		Vector3D p8=(GetPos0D(Vector3D(lx1,-ly1*sh, lz1*sh)));
		if (linemode != 2)
		{
			//mbs->DrawQuad(p1,p2,p3,p4);
			//mbs->DrawQuad(p8,p7,p6,p5);
		}

		double tilex = pow(2.,GetMBS()->GetDrawResolution());
		/*
		int ne = GetMBS()->NE();
		if (ne <= 16) {tilex = 16;}
		if (ne < 4) {tilex = 24;}
		if (ne >= 80) {tilex = 6;}

		if (GetMBS()->GetDrawResolution() == 0) {tilex = Maximum(1., tilex/6.); }
		else if (GetMBS()->GetDrawResolution() == 1) {tilex = Maximum(1., tilex/3.); }
		else if (GetMBS()->GetDrawResolution() == 3) {tilex = tilex*2; }
		else if (GetMBS()->GetDrawResolution() == 4) {tilex = tilex*4; }
		*/

		double tiley = tilex;

		TArray<Vector3D> points((int)(tilex+1)*(int)(tiley+1));
		TArray<double> vals((int)(tilex+1)*(int)(tiley+1));
		double v=0;

		int starts = 1;
		int ends = 6;
 		if (chladni || drawflat) {starts=3; ends=3;}

		for (int side=starts; side <= ends; side++)
		{
			points.SetLen(0); vals.SetLen(0);
			Vector3D p0, vx, vy;
			int tileyn = (int)tiley;
			if (side <= 2 || side >= 5) tileyn = 2;

			if (side == 1)
			{ //bottom
				p0 = Vector3D(-lx1,-ly1,-lz1);
				vx = Vector3D(2.*lx1/tilex,0,0);
				vy = Vector3D(0,0,2.*lz1/tileyn);
			}
			else if (side == 2)
			{ //top
				p0 = Vector3D(-lx1, ly1, lz1);
				vx = Vector3D(2.*lx1/tilex,0,0);
				vy = Vector3D(0,0,-2.*lz1/tileyn);
			}
			else if (side == 3)
			{ //front
				p0 = Vector3D(-lx1,-ly1, lz1);
				vx = Vector3D(2.*lx1/tilex,0,0);
				vy = Vector3D(0,2.*ly1/tileyn,0);

				if (chladni || drawflat) p0 = Vector3D(-lx1,-ly1, 0);
			}
			else if (side == 4)
			{ //back
				p0 = Vector3D(-lx1, ly1,-lz1);
				vx = Vector3D(2.*lx1/tilex,0,0);
				vy = Vector3D(0,-2.*ly1/tileyn,0);
			}
			else if (side == 5)
			{ //left
				p0 = Vector3D(-lx1, ly1,-lz1);
				vx = Vector3D(0,-2.*ly1/tilex,0);
				vy = Vector3D(0,0,2.*lz1/tileyn);
			}
			else if (side == 6)
			{ //right
				p0 = Vector3D( lx1,-ly1,-lz1);
				vx = Vector3D(0,2.*ly1/tilex,0);
				vy = Vector3D(0,0,2.*lz1/tileyn);
			}

			for (double iy = 0; iy <= tileyn+1e-10; iy++)
			{
				for (double ix = 0; ix <= tilex+1e-10; ix++)
				{
					if (!chladni || ! colormode)
					{
						Vector3D ploc = (p0+ix*vx+iy*vy);
						Vector3D pg = GetPos0D(ploc);
						Vector3D p(pg);
						points.Add(p);

						//if (side == 3) normals.Add(GetNormalNormed0D(ploc, xgd));
						//if (side == 4) normals.Add(-1.*GetNormalNormed0D(ploc, xgd));

						if (colormode)
							v = GetFieldVariableValue(*GetMBS()->GetActualPostProcessingFieldVariable(), ploc, true);

						vals.Add(v);
					}
					else
					{
						//draw chladni figures:
						Vector3D ploc = (p0+ix*vx+iy*vy);
						Vector3D pg = GetPos0D(ploc);
						//Vector3D pginit = GetPos0D(ploc,e0);
						Vector3D pginit = GetPos0D(ploc);
						points.Add(pginit);

						//normals.Add(GetNormalNormed0D(ploc, e0));

						v = pg.Z();
						vals.Add(v);

					}


				}
			}
			mbs->DrawColorQuads(points,vals,(int)tilex+1,(int)tileyn+1,colormode,linemode);
		}

	}

	/*
	if (concentratedmass1 != 0) 
	{
	Vector3D mp1(GetPos0D(Vector3D(-lx1,0,0))); 
	mbs->SetColor(colred);
	double r = pow(3.*concentratedmass1/(4.*MY_PI*GetRho()),1./3.); 
	//double r = ly*GetMBS()->GetMagnifyYZ()*4.;
	mbs->DrawSphere(mp1,r,12);
	}*/
};


