//#**************************************************************
//#
//# filename:             Plate2D.cpp
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						April 2006
//# description:          2D-plane stress plate element for floating frame of reference or absolute coordinates
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
#include "body2d.h"
#include "femathhelperfunctions.h"
#include "material.h"
#include "node.h"
#include "referenceframe2d.h"
#include "plate2d.h"
#include "elementdataaccess.h"
#include "constraint.h"

void Plate2D::SetPlate2D(int bodyindi, const TArray<int> &nodelist, int matnr, double thickness, const Vector3D &coli, int CMSelememti)
{
	bodyind = bodyindi;
	lz = thickness;

	TArray<Vector2D> p(NNodes());
	for(int i=1; i<=NNodes(); i++)
	{
		node[i-1] = nodelist(i);
		p(i) = Vector2D(GetNode(i).Pos().X(),GetNode(i).Pos().Y());
	}
	mass = lz * GetMBS()->GetMaterial(matnr).Density() * Area(ToP3D(p(1)),ToP3D(p(2)),ToP3D(p(3)),ToP3D(p(4)));

	xg = Vector(NNodes()*Dim()); //for intial conditions??? don't know why?
	xgd = xg; //for initial drawing???

	x_p0 = Vector(NNodes()*Dim()); //store nodal positions for drawing and constraints
	for(int i=1; i<=NNodes(); i++)
	{
		x_p0(2*i-1) = p(i).X(); 
		x_p0(2*i) = p(i).Y();
	}
	
	x_init = Vector(2*NNodes()*Dim());
//	for (int i=1; i<=NNodes()*Dim(); i++) x_init(i) = 0.0; //initial nodal displacements
//	for (int i=NNodes()*Dim()+1; i<=2*NNodes()*Dim(); i++) x_init(i) = 0.0; //initial nodal velocities
// vgl. FiniteElement3D::SetFiniteElement3D(...)	
	for (int i=1; i<=NNodes(); i++)
	{
		Node node(GetNode(i));
// {x1,y1,x2,y2,... vx1,...} -> twice the size
		if (node.X_Init().GetLen())
		{
			x_init.Copy(node.X_Init(), 1,       (i-1)*Dim()+1,          Dim());
			x_init.Copy(node.X_Init(), 1+Dim(), (i-1+NNodes())*Dim()+1, Dim());
		}
	}

	materialnum = matnr;
	col = coli;

	BuildDSMatrices();
	elementname = GetElementSpec();
}


void Plate2D::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Body2D::GetElementData(edc);

	ElementData ed;

	IVector v;
	for (int i=1; i <= NNodes(); i++)
	{
    v.Add(NodeNum(i));
	}
	SetElemDataIVector(edc, v, "Node_list");

	ed.SetInt(bodyind, "Body_index"); edc.Add(ed);

	ed.SetDouble(lz, "Thickness"); edc.Add(ed);

}

int Plate2D::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Body2D::SetElementData(edc);

	ElementData ed;

	IVector nodelist(NNodes());
	GetElemDataIVector(GetMBS(), edc, "Node_list", nodelist, 1);

	GetElemDataInt(GetMBS(), edc, "Body_index", bodyind);

	GetElemDataDouble(GetMBS(), edc, "Thickness", lz);

	for (int i=1; i <= NNodes(); i++)
	{
		node[i-1] = nodelist(i);
	}



	return rv;
}



void Plate2D::LinkToElements()
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

		for (int j = 1; j <= NNodes(); j++)
		{
			const Node& node = GetMBS()->GetNode(NodeNum(j));
			//Position:
			for (int i=1; i <= node.SOS(); i++)
			{
				AddLTG(node.Get(i));
			}
		}
		for (int i=1; i <= Nnonode; i++)
		{
			AddLTG(storeltg(i));
		}

		for (int j = 1; j <= NNodes(); j++)
		{
			const Node& node = GetMBS()->GetNode(NodeNum(j));
			//velocity:
			for (int i=1; i <= node.SOS(); i++)
			{
				AddLTG(node.Get(i+node.SOS()));
			}
		}
		for (int i=1; i <= Nnonode; i++)
		{
			AddLTG(storeltg(i+Nnonode));
		}
	}
}

void Plate2D::BuildDSMatrices() 
{
	if (!IsTrig())
	{
		GetIntegrationRule(x1,w1,orderxy);
		GetIntegrationRule(x2,w2,orderxy);
	}
	else
	{
		GetIntegrationRuleTrig(x1,x2,w1,orderxy);
		w2.SetLen(1); //only dummy
	}

	Matrix3D jac0, jacinv;
	jac0.SetSize(2,2); jacinv.SetSize(2,2);
	int kx1 = w2.Length();

	for (int i1=1; i1<=x1.GetLen(); i1++)
	{
		for (int i2=1; i2<=w2.GetLen(); i2++)
		{
			int ind = (i1-1)*kx1+(i2-1);

			Vector2D p;
			if (!IsTrig())
				p = Vector2D(x1(i1),x2(i2));
			else
				p = Vector2D(x1(i1),x2(i1)); //i2 always 1


			GetJacobi(jac0,p);

			jacdet[ind] = jac0.Det();
			jac0.GetInverse(jacinv);
			jacinv.TpYs();

			grad[ind].Init();
			grad[ind].SetSize(Dim(),NS());

			GetDSMatrix0(p,DS);
			Mult(jacinv, DS, grad[ind]);
			//global_uo << "Jacinv*DS=" << grad[ind] << "\n";
		}
	}
};


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//floating frame of reference formulation: for FFRF elements

//fill in sos x sos components, m might be larger
//compute stiffness matrix
void Plate2D::StiffnessMatrix(Matrix& m) 
{
	int dim = Dim(); 
	int ns = NS();

	double nu = Nu();   // parameters from material

	Matrix3D jac;
	jac.SetSize(dim,dim);
	Matrix3D jacinv;
	jacinv.SetSize(dim,dim);

	m.SetSize(ns*dim, ns*dim);
	m.SetAll(0);

	Matrix grad(dim,ns);
	Matrix B((dim-1)*3,dim*ns);

	//$!PG 2011-3-15:[
	//$!PG 2011-3-15: replace GetMaterialTensor(Matrix&) by ComputeElasticityMatrix(Matrix&)
	//Matrix C(3,3);
	//GetMaterial().GetMaterialTensor(C);
	Matrix C;
	GetMaterial().ComputeElasticityMatrix(C);
	//$!PG 2011-3-15:]

	GetIntegrationRule(x1,w1,orderxyM);
	GetIntegrationRule(x2,w2,orderxyM);

	for (int i1=1; i1<=x1.GetLen(); i1++)
	{
		for (int i2=1; i2<=x2.GetLen(); i2++)
		{
			Vector2D p(x1(i1),x2(i2));

			GetJacobi(jac,p);
			double jacdet = jac.Det();

			jac.GetInverse(jacinv);
			jacinv.TpYs();

			GetDSMatrix0(p,DS);
			Mult(jacinv, DS, grad);
			//grad contains derivatives in compressed form:
			//grad(1,1..ns) = da/dx     where a = du/dqi or dv/dqi
			//grad(2,1..ns) = da/dy

			for (int i = 1; i <= ns; i++)
			{
				B(1,2*i-1) = grad(1,i);
				B(2,2*i)   = grad(2,i);
				B(3,2*i-1) = grad(2,i); //compute shearing strain, not green shear strain!
				B(3,2*i)   = grad(1,i); //compute shearing strain, not green shear strain!
			}

			double fact = fabs (jacdet) * w1(i1)*w2(i2) * lz;
			m += fact * (B.GetTp()*C*B);
		}
	}

	//UO() << "Klin=" << (-1.)*m << "\n";
} 

//Shabana p. 209-211, inertia terms, 2D
void Plate2D::GetI1(Vector& I1) 
{
	//3% computime
	//compute 'moment of mass' of element:
	int dim = Dim(); 
	int ns = NS();
	I1.SetLen(dim);

	I1.SetAll(0);

	Matrix3D jac;
	jac.SetSize(dim,dim);

	GetIntegrationRule(x1,w1,orderxyH); //linear, but jac.Det() is nonlinear ...
	GetIntegrationRule(x2,w2,orderxyH);

	for (int i1=1; i1<=x1.GetLen(); i1++)
	{
		for (int i2=1; i2<=x2.GetLen(); i2++)
		{
			Vector2D p(x1(i1),x2(i2));

			GetJacobi(jac,p);
			double jacdet = jac.Det();

			Vector2D pp = GetPos2Drel0(p);

			I1(1) += fabs (jacdet) * Rho() * w1(i1)*w2(i2) * pp.X() * lz;
			I1(2) += fabs (jacdet) * Rho() * w1(i1)*w2(i2) * pp.Y() * lz;
		}
	}

}

double Plate2D::GetIkl(int k, int l) 
{
	//5% computime
	//compute 'mass moment of inertia' of element:
	double Ikl = 0;

	int dim = Dim(); 
	int ns = NS();

	Matrix3D jac;
	jac.SetSize(dim,dim);

	GetIntegrationRule(x1,w1,orderxyM); //quadratic, but jac.Det() is nonlinear ...
	GetIntegrationRule(x2,w2,orderxyM);

	for (int i1=1; i1<=x1.GetLen(); i1++)
	{
		for (int i2=1; i2<=x2.GetLen(); i2++)
		{
			Vector2D p(x1(i1),x2(i2));

			GetJacobi(jac,p);
			double jacdet = jac.Det();
			Vector2D pp = GetPos2Drel0(p);

			Ikl += fabs (jacdet) * Rho() * w1(i1)*w2(i2) * pp(k) * pp(l) * lz;
		}
	}
	return Ikl;
}


void Plate2D::GetIbarkl(int k, int l, Vector& I1) 
{
	//18% computime
	//compute int_V (rho*x_k * S_l) d_V of element:
	int dim = Dim(); 
	int ns = NS();

	I1.SetLen(2*ns);
	I1.SetAll(0);

	Matrix3D jac;
	jac.SetSize(dim,dim);

	GetIntegrationRule(x1,w1,orderxyM); 
	GetIntegrationRule(x2,w2,orderxyM);

	for (int i1=1; i1<=x1.GetLen(); i1++)
	{
		for (int i2=1; i2<=x2.GetLen(); i2++)
		{
			Vector2D p(x1(i1),x2(i2));

			GetJacobi(jac,p);
			double jacdet = jac.Det();
			Vector2D pp = GetPos2Drel0(p);

			for (int i = 1; i<= ns; i++)
			{
				I1(2*(i-1)+l) += fabs (jacdet) * Rho() * w1(i1)*w2(i2) * pp(k) * GetS0(p,i) * lz;
			}
		}
	}
}

//Sbar is a matrix ..., equal to GetH
void Plate2D::GetSbar(Matrix& Sbar) 
{
	//6.2% computime

	//compute int_V (rho*x_k * S_l) d_V of element:
	int dim = Dim(); 
	int ns = NS();

	Sbar.SetSize(dim,dim*ns);
	Sbar.SetAll(0);

	Matrix3D jac;
	jac.SetSize(dim,dim);

	GetIntegrationRule(x1,w1,orderxyH);
	GetIntegrationRule(x2,w2,orderxyH);

	for (int i1=1; i1<=x1.GetLen(); i1++)
	{
		for (int i2=1; i2<=x2.GetLen(); i2++)
		{
			Vector2D p(x1(i1),x2(i2));

			GetJacobi(jac,p);
			double jacdet = jac.Det();

			for (int i = 1; i<= ns; i++)
			{
				Sbar(1,dim*(i-1)+1) += fabs (jacdet) * Rho() * w1(i1)*w2(i2) *  GetS0(p,i) * lz;
				Sbar(2,dim*(i-1)+2) += fabs (jacdet) * Rho() * w1(i1)*w2(i2) *  GetS0(p,i) * lz;
			}
		}
	}
}

void Plate2D::GetSbarkl(int k, int l, Matrix& Sbar) 
{
	//46% computime
	//compute int_V (rho*x_k * S_l) d_V of element:
	int dim = Dim(); 
	int ns = NS();

	Sbar.SetSize(2*ns, 2*ns);
	Sbar.SetAll(0);

	Matrix3D jac;
	jac.SetSize(dim,dim);

	GetIntegrationRule(x1,w1,orderxyM);
	GetIntegrationRule(x2,w2,orderxyM);

	for (int i1=1; i1<=x1.GetLen(); i1++)
	{
		for (int i2=1; i2<=x2.GetLen(); i2++)
		{
			Vector2D p(x1(i1),x2(i2));

			GetJacobi(jac,p);
			double jacdet = jac.Det();

			for (int i = 1; i<= ns; i++)
			{
				for (int j = 1; j<= ns; j++)
				{
					Sbar(2*(i-1)+k, 2*(j-1)+l) += fabs (jacdet) * Rho() * w1(i1)*w2(i2) *  GetS0(p,i) * GetS0(p,j) * lz;
					//original Fehler bei 2*(i-1)+l !!!! Überprüfen!!!
					//Sbar(2*(i-1)+k, 2*(i-1)+l) += fabs (jacdet) * Rho() * w1(i1)*w2(i2) *  GetS0(p,i) * GetS0(p,j) * lz;
				}
			}
		}
	}
	//UO() << "+++++++++++++++++++++++++++++++++\n";
	//UO() << "+++++++++++++++++++++++++++++++++\n";
	//UO() << "+++++++++++++++++++++++++++++++++\n";
	//UO() << "check GetSbarkl in Plate2D !!!!!!\n";
	//UO() << "+++++++++++++++++++++++++++++++++\n";
	//UO() << "+++++++++++++++++++++++++++++++++\n";
	//UO() << "+++++++++++++++++++++++++++++++++\n";
}


//-1..+1 based!!!
Vector2D Plate2D::GetPos2D(const Vector2D& p_loc) const
{
	Vector2D p;

	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += GetS0(p_loc,j)*(XG((j-1)*Dim()+i) + x_p0((j-1)*Dim()+i));
		}
	}
	return p;
};

//-1..+1 based!!!
Vector2D Plate2D::GetDisplacement2D(const Vector2D& p_loc) const
{
	Vector2D p;

	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			double s = GetS0(p_loc,j);
			p(i) += GetS0(p_loc,j)*xg((j-1)*Dim()+i);
		}
	}
	return p;
};

//-1..+1 based!!!
Vector2D Plate2D::GetVel2D(const Vector2D& p_loc) const
{
	Vector2D p;

	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += GetS0(p_loc,j)*(XGP((j-1)*Dim()+i));
		}
	}
	return p;
};

//-1..+1 based!!!
Vector2D Plate2D::GetPos2DD(const Vector2D& p_loc) const
{
	return GetPos2DD(p_loc, 1);
};

//-1..+1 based!!!
Vector2D Plate2D::GetPos2DD(const Vector2D& p_loc, int use_magnification) const
{
	Vector2D p;
	double factor = GetMBS()->GetDOption(105);
	if (!use_magnification) factor = 1;

	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			double s = GetS0(p_loc,j);
			p(i) += GetS0(p_loc,j)*(factor*xgd((j-1)*Dim()+i) + x_p0((j-1)*Dim()+i));
		}
	}
	return p;
};

//-1..+1 based!!!
Vector2D Plate2D::GetDisplacement2DD(const Vector2D& p_loc) const
{
	Vector2D p;

	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			double s = GetS0(p_loc,j);
			p(i) += GetS0(p_loc,j)*xgd((j-1)*Dim()+i);
		}
	}
	return p;
};

//-1..+1 based!!!
Vector2D Plate2D::GetVel2DD(const Vector2D& p_loc) const
{
	Vector2D p;

	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += GetS0(p_loc,j)*(XGPD((j-1)*Dim()+i));
		}
	}
	return p;
};

//-1..+1 based!!!
//relative position!, UNDEFORMED
Vector2D Plate2D::GetPos2Drel0(const Vector2D& p_loc) const
{
	Vector2D p;

	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += GetS0(p_loc,j)*(x_p0((j-1)*Dim()+i));
		}
	}
	return p;
};

//-1..+1 based!!!
//relative position!
Vector2D Plate2D::GetPos2Drel(const Vector2D& p_loc) const
{
	Vector2D p;

	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += GetS0(p_loc,j)*(XG((j-1)*Dim()+i) + x_p0((j-1)*Dim()+i));
		}
	}
	return p;
};

//relative velocity!
Vector2D Plate2D::GetVel2Drel(const Vector2D& p_loc) const
{
	Vector2D p;

	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += GetS0(p_loc,j)*(XGP((j-1)*Dim()+i));
		}
	}
	return p;
};

//relative position!
Vector2D Plate2D::GetPos2DrelD(const Vector2D& p_loc, int use_magnification) const
{
	double factor = GetMBS()->GetDOption(105);
	if (!use_magnification) factor = 1;

	Vector2D p;

	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += GetS0(p_loc,j)*(xgd((j-1)*Dim()+i)*factor + x_p0((j-1)*Dim()+i));
		}
	}
	return p;
};

//relative velocity!
Vector2D Plate2D::GetVel2DrelD(const Vector2D& p_loc) const
{
	Vector2D p;

	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += GetS0(p_loc,j)*(XGPD((j-1)*Dim()+i));
		}
	}
	return p;
};


Vector2D Plate2D::GetRefPos2D() const 
{
	double fact = 1./NNodes();
	Vector2D p(0.,0.);
	for (int j=1; j <= Dim(); j++)
	{
		for (int i=1; i <= NNodes(); i++)
		{
			p(j) += fact*(XG((i-1)*2+j)+x_p0((i-1)*2+j));
		}
	}
	return p;
};

Vector2D Plate2D::GetRefPos2DD() const 
{
	double fact = 1./NNodes();
	Vector2D p(0.,0.);
	for (int j=1; j <= Dim(); j++)
	{
		for (int i=1; i <= NNodes(); i++)
		{
			p(j) += fact*(XGD((i-1)*2+j)+x_p0((i-1)*2+j));
		}
	}
	return p;
};




void Plate2D::GetH(Matrix& H) 
{
	if (Hmatrix.Getrows() == FlexDOF())
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
		jac.SetSize(dim,dim);

		DS.SetSize(dim,ns);
		SV.SetLen(ns);

		GetIntegrationRule(x1,w1,orderxyH);
		GetIntegrationRule(x2,w2,orderxyH);

		for (int i1=1; i1<=x1.GetLen(); i1++)
		{
			for (int i2=1; i2<=x2.GetLen(); i2++)
			{
				Vector2D p(x1(i1),x2(i2));
				GetS0(SV,p);

				GetDSMatrix0(p,DS);
				GetJacobi(jac,p);
				//GetJacobi(jac,p,DS);
				double jacdet = jac.Det();
				double fact = fabs (jacdet) * w1(i1)*w2(i2)*lz;

				for (int i=0; i<ns; i++)
				{
					for (int j=1; j<=dim; j++)
					{
						H(i*dim+j,j)+=fact*SV(i+1);
					}
				}
			}
		}
		Hmatrix = H;
	}
}

void Plate2D::EvalM(Matrix& m, double t) 
{
	if (massmatrix.Getcols() == FlexDOF())
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
		jac.SetSize(dim,dim);

		Matrix HL(ns*dim,dim);
		DS.SetSize(dim,ns);

		GetIntegrationRule(x1,w1,orderxyM); //optimal: 6x3x3, 6x1x1 geht auch!!!!
		GetIntegrationRule(x2,w2,orderxyM);

		for (int i1=1; i1<=x1.GetLen(); i1++)
		{
			for (int i2=1; i2<=x2.GetLen(); i2++)
			{
				Vector2D p(x1(i1),x2(i2));
				GetS0(SV,p);

				for (int i=0; i<ns; i++)
				{
					for (int j=1; j<=dim; j++)
					{
						HL(i*dim+j,j)=SV(i+1);
					}
				}

				//GetDSMatrix0(p,DS);
				//GetJacobi(jac,p,DS);
				GetJacobi(jac,p);

				//global_uo << "jac=" << jac << "\n";

				double jacdet = jac.Det();
				m += fabs (jacdet) * Rho() * w1(i1)*w2(i2)*lz * (HL*HL.GetTp());
			}
		}

		/*
		for (int k=0; k < NNodes(); k++)
		{
		m(1+k*dim,1+k*dim) += concentratedmass[k];
		m(2+k*dim,2+k*dim) += concentratedmass[k];
		}
		*/
		//global_uo << "massmatrix=" << m << "\n";

		massmatrix = m;

		/*
		double m1 = 0;
		for (int i = 1; i <= ns*dim; i++)
		{
		for (int j = 1; j <= ns*dim; j++)
		{
		m1 += m(i,j);
		}
		}
		m1*=0.5; //Kinetic energy T=1/2*v^t*M*v;
		global_uo << "mass2=" << m1 << "\n";
		*/

	}
};

int useG=0;




void Plate2D::EvalF2(Vector& f, double t) 
{
	Body2D::EvalF2(f,t);
	TMStartTimer(22);

	int sos = FlexDOF();
	SetComputeCoordinates();
	static Vector fadd;


	int ns = NS();
	int dim = Dim();
	//double u;

	//$ PG 2011-3-15:[
	//$ PG 2011-3-15: commented, since Matrix3D C is never used throughout this method
	//Matrix3D C;
	//GetMaterial().GetMaterialTensor2D(C);
	//$ PG 2011-3-15:]

	Matrix3D piola1;
	// Vector3D eps, sigma;
	piola1.SetSize(2,2);

	temp.SetLen(SOS());
	fadd.SetLen(SOS());
	fadd.SetAll(0);


	GetIntegrationRule(x1,w1,orderxy);
	GetIntegrationRule(x2,w2,orderxy);

	int kx1 = x2.Length();

			//int_A (P : d F / d q) dA
			//P = F * S
			//S = D : E
			//E = 1/2*(F^T F-I)

	for (int i1=1; i1<=x1.GetLen(); i1++)
	{
		for (int i2=1; i2<=x2.GetLen(); i2++)
		{
			int ind = (i1-1)*kx1+(i2-1);

			FirstPiolaTensor(ind, xg, piola1);

			for (int j=1; j <= Dim(); j++)
			{
				for (int i = 0; i < ns; i++)
				{
					temp(Dim()*i+j) = grad[ind](1, i+1)*piola1(j,1)
						+ grad[ind](2, i+1)*piola1(j,2);
				}
			}

			fadd.MultAdd(fabs (jacdet[ind]) * w1(i1)*w2(i2)*lz,temp);
		} 
	}



	f -= fadd;
	TMStopTimer(22);

	if (GetMassDamping() != 0)
	{
		//UO() << "Warning: damping only includes diagonal of mass matrix!!!\n";
		// +++++++ damping: +++++++
		for (int i = 1; i <= SOS(); i++)
			xg(i) = XGP(i);

		if (massmatrix.Getcols() == SOS())
		{
			Mult(massmatrix,xg,temp);
			/*
			for (int i = 1; i <= SOS(); i++)
			temp(i) = xg(i)*massmatrix(i,i);
			*/
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
}; 

// Jacobian for linear mechanics
void Plate2D::JacobianF2(double t, Matrix& m, IVector& colref)
{
	
	int sos = SOS();

	static Vector locjacf20; //for local jacobians
	static Vector locjacf21; //for local jacobians 
	locjacf20.SetLen(sos);
	locjacf21.SetLen(sos);


	double numdiffepsi = GetMBS()->NumSolver().NumDiffepsi();
	double eps;

	//only non-symmetric mode!!!
	locjacf20.SetAll(0);
	EvalF2(locjacf20,t);


	colref.SetLen(0);

	int usesparseK = UseSparseK();
	if (!usesparseK)
	{
		const	TArray<int>& ltg = GetLTGArray();
		for (int i = 1; i <= SS(); i++)
		{
// (AD) changed () to .Get()
			colref.Add(ltg.Get(i));
//			colref.Add(ltg(i));
		}
	}

	//use constraints_nodouble:
	//special for contact where one element is linked several times to the same contact constraint 
	for (int i=1; i <= constraints_nodouble.Length(); i++)
	{
		const Constraint& c = *(constraints_nodouble(i));
		for (int j=1; j <= c.SS(); j++)
		{
			colref.Add(c.LTG(j));
		}
	}

	m.SetSize(sos,colref.Length());
	m.SetAll(0.);
	Matrix tempm(sos, sos);
	tempm.SetAll(0.);

	int begcnt = 1;
	if (FastStiffnessMatrix() == 1)
	{
		// sos dofs are computed by stiffness matrix routine, remaining dofs by finite differencing
		begcnt = sos+1;
		StiffnessMatrix(tempm);
	}
	else
	{
		begcnt = 1;
	} 

	double storex;
	for (int i = 1; i < begcnt; i++)
		for(int j = 1; j <= sos; j++)
		{
			m(j,i) -= tempm(j,i);
		}	

	for (int i = begcnt; i <= colref.Length(); i++)
	{
		eps = numdiffepsi*Maximum(1e-2,fabs(XGG(colref(i))));

		if (GetMBS()->NumSolver().SymmetricJacobian())
		{
			storex = XGG(colref(i));
			XGG(colref(i)) += eps;

			locjacf21.SetAll(0);
			EvalF2(locjacf21,t);

			XGG(colref(i)) -= 2*eps;
			locjacf20.SetAll(0);
			EvalF2(locjacf20,t);
			XGG(colref(i)) = storex;
		}
		else
		{
			storex = XGG(colref(i));
			XGG(colref(i)) += 2*eps;
			locjacf21.SetAll(0);
			EvalF2(locjacf21,t);
			XGG(colref(i)) = storex;
		}

		//UO() << "f1=" << locjacf21 << "\n";

		for (int j=1; j<=sos;j++)
		{
			m(j,i) += (locjacf21(j)-locjacf20(j))/(2.*eps);
		}
	}

	 //global_uo << "aha, stiffness\n" << m << "\n";
}

double Plate2D::PostNewtonStep(double t)
{
	return 0;  
}


void Plate2D::PostprocessingStep()
{
	//
}

char * str_inelastic_strain_error_output = "inelastic strain error output";

void Plate2D::GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables)
{
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_displacement,
		FieldVariableDescriptor::FVCI_y);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_stress,
		FieldVariableDescriptor::FVCI_y, true);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_stress_mises);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_total_strain,
		FieldVariableDescriptor::FVCI_y, true);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_inelastic_strain,
		FieldVariableDescriptor::FVCI_y, true);
	FieldVariableDescriptor::AddTypeIntoArray(variables, str_inelastic_strain_error_output,
		FieldVariableDescriptor::FVCI_z, true);		//*YV: seems to be symmetric
}

double Plate2D::GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector2D & local_position, bool flagD)
{
	double nu = Nu();   // material parameters

	if(fvd.VariableType() == FieldVariableDescriptor::FVT_displacement)
	{
		Vector2D u;
		if(flagD) 
			u = GetDisplacement2DD(local_position);
		else
			u = GetDisplacement2D(local_position);
		return fvd.GetComponent(u);
	}

	Matrix3D F,G;
	F.SetSize(2,2);
	G.SetSize(2,2);
	Matrix3D strain;
	strain.SetSize(2,2);
	Matrix3D stress;
	stress.SetSize(2,2);

	double la=Em() * nu / ((1.+nu)*(1.-2.*nu));
	double mu=Em() / 2. / (1.+nu);

	Vector3D sigma;
	double k;
	ConstMatrix<36> C(3,3);
	GetMaterial().ComputeElasticityMatrix(C);

	//if(Is_Planestress())	// plane stress	
	//{
	//	k = Em()/(1.-nu*nu);
	//	Matrix3D C_pse(k,k*nu,0,k*nu,k,0,0,0,k*(1.-nu)/2.); //Bathe p. 229, mit Gleitung!!!
	//	C=C_pse;
	//}
	//else //if(!Is_Planestress())		// plane strain
	//{
	//	k = Em()/((1.+nu)*(1.-2*nu));
	//	Matrix3D C_psa(k*(1-nu),k*nu,0,k*nu,k*(1-nu),0,0,0,k*(1.-2*nu)/2.);
	//	C=C_psa;
	//}

	if(flagD) 
	{
		if (NNodes() > 4) GraduD(local_position,F); //small deformations, stresses may look smoother
		else GraduD(Vector2D(0,0),F);     //large deformation, element jacobian not diagonal
	}
	else
	{
		SetComputeCoordinates();
		if (NNodes() > 4) Gradu(local_position,xg,F); //small deformations, stresses may look smoother
		else Gradu(Vector2D(0,0),xg,F);     //large deformation, element jacobian not diagonal
	}

	//global_uo << "xgd=" << xgd << "\n";
	//global_uo << "gradu=" << F << "\n";

	if(!useG)
	{			
		F(1,1) += 1;
		F(2,2) += 1;

		strain = 0.5 * (F.GetTp() * F);
		strain(1,1) -= 0.5; strain(2,2) -= 0.5;
	}
	else
	{
		G=F;
		strain = 0.5 * (G.GetTp() + G + G.GetTp() * G);
	}


	Vector3D eps(strain(1,1), strain(2,2), 2.*strain(1,2));
	sigma = C*eps;

	Vector3D inel_strain(0.,0.,0.);
	Vector3D inel_strain_error(0.,0.,0.);

	stress(1,1) = sigma(1);
	stress(1,2) = sigma(3);
	stress(2,1) = sigma(3);
	stress(2,2) = sigma(2);

	switch(fvd.VariableType())
	{
	case FieldVariableDescriptor::FVT_stress_mises: return stress.Mises();
	case FieldVariableDescriptor::FVT_stress: return fvd.GetComponent(stress);
	case FieldVariableDescriptor::FVT_total_strain: return fvd.GetComponent(strain);
	case FieldVariableDescriptor::FVT_inelastic_strain:
		{
			Matrix3D s;
			s.Set22(inel_strain(1),0.5*inel_strain(3),0.5*inel_strain(3),inel_strain(2));
			return fvd.GetComponent(s);
		}
	case FieldVariableDescriptor::FVT_problem_specific:
		{
			if(fvd.GetTextualIdentifierWithoutComponents() == str_inelastic_strain_error_output)
			{
				Matrix3D se;
				Vector3D r1(0.,0.,0.);
				Vector3D r2(0.,0.,0.);
				Vector3D r3(0.,0.,0.);
				r1.X()=inel_strain(1); r1.Y()=0.5*inel_strain(3); r1.Z()=inel_strain_error(1);
				r2.X()=0.5*inel_strain(3); r2.Y()=inel_strain(2); r2.Z()=inel_strain_error(2);
				r3.X()=inel_strain_error(1); r3.Y()=inel_strain_error(2); r3.Z()=inel_strain_error(3);
				se.Set(r1,r2,r3); 
				return fvd.GetComponent(se);
			}
		}
	}

	return FIELD_VARIABLE_NO_VALUE;
}



void Plate2D::DrawElementPreProc() 
{
	//$MaSch 2012-11-27: moved here from Plate2D::DrawElement() to ensure that SetDrawCoordinates is called also if AltShape is used
	SetDrawCoordinates();

	if (!GetMBS()->GetIOption(118)) return;

	FieldVariableDescriptor * fvd = GetMBS()->GetActualPostProcessingFieldVariable();
	if(fvd != NULL)
	{
		//works only for 4-noded quadrilateral element:
		//nodal stress values
		double v[4] = {
			GetFieldVariableValue(*fvd,Vector2D( 1, 1), true),
			GetFieldVariableValue(*fvd,Vector2D(-1, 1), true),
			GetFieldVariableValue(*fvd,Vector2D(-1,-1), true),
			GetFieldVariableValue(*fvd,Vector2D( 1,-1), true)
		};
		for (int i=1; i <= 4; i++)
		{
			GetNode(i).GetDrawTemp() += v[i-1];
			GetNode(i).GetDrawTempCnt()++;
		}
	}
}

void Plate2D::DrawElement() 
{
	mbs->SetColor(col);

	//$MaSch 2012-11-27: moved this to Plate2D::DrawElementPreProc(), since - if AltShape is used - this is never called and, thus, xgd stays constant at its initial value
	//SetDrawCoordinates();

	int linemode = 1; //0=no lines, 1=outline+color, 2=outline, 3=elementline+color, 4=elementline
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

	int res = GetMBS()->GetDrawResolution()-1;

	if (colormode == 0 && NS() <= 4)
	{
		res = 0;
	}
	if (res < 0) res = 0;

	double tilex = pow(2.,res);

	double tiley = tilex;
	int tileyn = (int)tiley;

	double v=0;

	TArray<Vector3D> points((int)(tilex+1)*(int)(tiley+1));
	TArray<double> vals((int)(tilex+1)*(int)(tiley+1));
	points.SetLen(0); vals.SetLen(0);

	if (GetMBS()->GetIOption(118) && NNodes() == 4)
	{
		//add rectangular grid of points and values
		points.Add(ToP3D(GetPos2DD(Vector2D(-1,-1), 1))); //3
		points.Add(ToP3D(GetPos2DD(Vector2D( 1,-1), 1))); //4
		points.Add(ToP3D(GetPos2DD(Vector2D(-1, 1), 1))); //2
		points.Add(ToP3D(GetPos2DD(Vector2D( 1, 1), 1))); //1
		vals.Add(GetNode(3).GetDrawTemp());
		vals.Add(GetNode(4).GetDrawTemp());
		vals.Add(GetNode(2).GetDrawTemp());
		vals.Add(GetNode(1).GetDrawTemp());

		mbs->DrawColorQuads(points,vals,2,2,colormode,linemode);
	}
	else
	{
		Vector2D p0, vx, vy;

		double shrink = GetMBS()->GetDOption(106);//0.995;
		p0 = Vector2D(-1*shrink,-1*shrink);
		vx = Vector2D(2./tilex*shrink,0);
		vy = Vector2D(0,2./tileyn*shrink);

		for (double iy = 0; iy <= tileyn+1e-8; iy++)
		{
			for (double ix = 0; ix <= tilex+1e-8; ix++)
			{
				Vector2D ploc = (p0+ix*vx+iy*vy);
				Vector3D pg = ToP3D(GetPos2DD(ploc, 1));
				points.Add(pg);
				if (colormode)
					v = GetFieldVariableValue(*GetMBS()->GetActualPostProcessingFieldVariable(), ploc, true);
				vals.Add(v);
			}
		}
		mbs->DrawColorQuads(points,vals,(int)tilex+1,(int)tileyn+1,colormode,linemode);
	}

};



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//nodal coordinates of points in 2D; optional rotation matrix should be included to couple with 3D motion!
//initial velocities!!!!!
Plate2Dlin::Plate2Dlin(MBS* mbsi, int bodyindi, const Vector2D& xc1, const Vector2D& xc2, const Vector2D& xc3, const Vector2D& xc4, 
											 const Vector2D& vc1, const Vector2D& vc2, const Vector2D& vc3, const Vector2D& vc4, 
											 double rhoi, double Emi, double nui, double thickness, const Vector3D& coli, int* nodenums):
Plate2D(mbsi)
{
	bodyind = bodyindi;
	Vector2D p1(xc1(1),xc1(2));
	Vector2D p2(xc2(1),xc2(2));
	Vector2D p3(xc3(1),xc3(2));
	Vector2D p4(xc4(1),xc4(2));

	lz = thickness;

	//Transformation must be set before call to ToP3D
	Node n1(2, bodyind, Vector3D(p1.X(),p1.Y(),0.));
	Node n2(2, bodyind, Vector3D(p2.X(),p2.Y(),0.));
	Node n3(2, bodyind, Vector3D(p3.X(),p3.Y(),0.));
	Node n4(2, bodyind, Vector3D(p4.X(),p4.Y(),0.));

	if (nodenums == 0)
	{
		node[0] = mbsi->AddBodyNode(&n1); //check if node exists in body!!!
		node[1] = mbsi->AddBodyNode(&n2); //add this function to mbs.h, add bodyindex to class Node
		node[2] = mbsi->AddBodyNode(&n3);
		node[3] = mbsi->AddBodyNode(&n4);
	}
	else
	{
		node[0] = nodenums[0];
		node[1] = nodenums[1];
		node[2] = nodenums[2];
		node[3] = nodenums[3];
	}

	mass = lz * rhoi * Area(ToP3D(p1),ToP3D(p2),ToP3D(p3),ToP3D(p4));

	xg = Vector(NNodes()*Dim()); //for intial conditions??? don't know why?
	xgd = xg; //for initial drawing???

	x_p0 = Vector(NNodes()*Dim()); //store nodal positions for drawing and constraints
	x_p0(1) = xc1.X(); x_p0(2) = xc1.Y();
	x_p0(3) = xc2.X(); x_p0(4) = xc2.Y();
	x_p0(5) = xc3.X(); x_p0(6) = xc3.Y();
	x_p0(7) = xc4.X(); x_p0(8) = xc4.Y();

	x_init = Vector(2*NNodes()*Dim());
	for (int i=1; i<=NNodes()*Dim(); i++) x_init(i) = 0; //initial nodal displacements

	x_init( 9) = vc1.X(); x_init(10) = vc1.Y(); //intial nodal velocities
	x_init(11) = vc2.X(); x_init(12) = vc2.Y();
	x_init(13) = vc3.X(); x_init(14) = vc3.Y();
	x_init(15) = vc4.X(); x_init(16) = vc4.Y();

	Material mat(GetMBS(), rhoi, Emi, nui, 1, 0);
	AddMaterial(mat);

	//nu = nui;
	//Em = Emi;
	//rho = rhoi;
	col = coli;
	//concentratedmass[0] = concentratedmass[1] = concentratedmass[2] = concentratedmass[3] = 0;

	orderxy = Nip()*2-1; //2x2 integration points
	orderxyM = 2;
	orderxyH = 2;

	BuildDSMatrices();
	elementname = GetElementSpec();
};

Plate2Dlin::Plate2Dlin(MBS* mbsi, int bodyindi, const Vector2D& xc1, const Vector2D& xc2, const Vector2D& xc3, const Vector2D& xc4, 
											 const Vector2D& vc1, const Vector2D& vc2, const Vector2D& vc3, const Vector2D& vc4, 
											 int matnr, double thickness, const Vector3D& coli, int* nodenums):
Plate2D(mbsi)
{
	bodyind = bodyindi;
	Vector2D p1(xc1(1),xc1(2));
	Vector2D p2(xc2(1),xc2(2));
	Vector2D p3(xc3(1),xc3(2));
	Vector2D p4(xc4(1),xc4(2));

	lz = thickness;

	//Transformation must be set before call to ToP3D
	Node n1(2, bodyind, Vector3D(p1.X(),p1.Y(),0.));
	Node n2(2, bodyind, Vector3D(p2.X(),p2.Y(),0.));
	Node n3(2, bodyind, Vector3D(p3.X(),p3.Y(),0.));
	Node n4(2, bodyind, Vector3D(p4.X(),p4.Y(),0.));

	if (nodenums == 0)
	{
		node[0] = mbsi->AddBodyNode(&n1); //check if node exists in body!!!
		node[1] = mbsi->AddBodyNode(&n2); //add this function to mbs.h, add bodyindex to class Node
		node[2] = mbsi->AddBodyNode(&n3);
		node[3] = mbsi->AddBodyNode(&n4);
	}
	else
	{
		node[0] = nodenums[0];
		node[1] = nodenums[1];
		node[2] = nodenums[2];
		node[3] = nodenums[3];
	}

	mass = lz * GetMBS()->GetMaterial(matnr).Density() * Area(ToP3D(p1),ToP3D(p2),ToP3D(p3),ToP3D(p4));

	xg = Vector(NNodes()*Dim()); //for intial conditions??? don't know why?
	xgd = xg; //for initial drawing???

	x_p0 = Vector(NNodes()*Dim()); //store nodal positions for drawing and constraints
	x_p0(1) = xc1.X(); x_p0(2) = xc1.Y();
	x_p0(3) = xc2.X(); x_p0(4) = xc2.Y();
	x_p0(5) = xc3.X(); x_p0(6) = xc3.Y();
	x_p0(7) = xc4.X(); x_p0(8) = xc4.Y();

	x_init = Vector(2*NNodes()*Dim());
	for (int i=1; i<=NNodes()*Dim(); i++) x_init(i) = 0.0; //initial nodal displacements

	x_init( 9) = vc1.X(); x_init(10) = vc1.Y(); //intial nodal velocities
	x_init(11) = vc2.X(); x_init(12) = vc2.Y();
	x_init(13) = vc3.X(); x_init(14) = vc3.Y();
	x_init(15) = vc4.X(); x_init(16) = vc4.Y();

	materialnum = matnr;

	//nu = nui;
	//Em = Emi;
	//rho = rhoi;
	col = coli;
	//concentratedmass[0] = concentratedmass[1] = concentratedmass[2] = concentratedmass[3] = 0;

	orderxy = Nip()*2-1; //2x2 integration points
	orderxyM = 2;
	orderxyH = 2;

	BuildDSMatrices();
	elementname = GetElementSpec();

};


void Plate2Dlin::SetPlate2Dlin(int bodyindi, const TArray<int> &nodelist, int matnr, double thickness, const Vector3D &coli, int CMSelememti)
{
	orderxy = Nip()*2-1; //2x2 integration points
	orderxyM = 2;
	orderxyH = 2;

	Plate2D::SetPlate2D(bodyindi, nodelist, matnr, thickness, coli, CMSelememti);
	elementname = GetElementSpec();
}

void Plate2Dlin::GetS0(Vector& sf, const Vector2D& ploc) const
{
	//Nodes always sorted counterclock-wise in x-y (r-s) plane
	double r = ploc.X();
	double s = ploc.Y();
	sf.SetLen(NS());

	//order of nodes:
	//
	//  node 2 +------+ node 1
	//         |      |
	//         |      |
	//         |      |
	//  node 3 +------+ node 4
	//

	sf(1) = 0.25*(1+r)*(1+s);  //Node 1
	sf(2) = 0.25*(1-r)*(1+s);  //Node 2
	sf(3) = 0.25*(1-r)*(1-s);  //Node 3
	sf(4) = 0.25*(1+r)*(1-s);  //Node 4
}

void Plate2Dlin::GetDSMatrix0(const Vector2D& ploc, Matrix& sf) const
{
	double r = ploc.X();
	double s = ploc.Y();
	sf.SetSize(Dim(),NS());

	//D_sf/D_r :
	sf(1,1) = 0.25*(1+s);  //Node 1
	sf(1,2) =-0.25*(1+s);  //Node 2
	sf(1,3) =-0.25*(1-s);  //Node 3
	sf(1,4) = 0.25*(1-s);  //Node 4

	//D_sf/D_s :
	sf(2,1) = 0.25*(1+r);  //Node 1
	sf(2,2) = 0.25*(1-r);  //Node 2
	sf(2,3) =-0.25*(1-r);  //Node 3
	sf(2,4) =-0.25*(1+r);  //Node 4

}

double Plate2Dlin::GetS0(const Vector2D& ploc, int i) const
{
	//Nodes always sorted counterclock-wise in x-y (r-s) plane
	double r = ploc.X();
	double s = ploc.Y();

	switch(i)
	{
	case 1: return 0.25*(1+r)*(1+s);  //Node 1
	case 2: return 0.25*(1-r)*(1+s);  //Node 2
	case 3: return 0.25*(1-r)*(1-s);  //Node 3
	case 4: return 0.25*(1+r)*(1-s);  //Node 4
	default: return 0;
	}
}

//D S_i / Dx_j
double Plate2Dlin::GetDS0(const Vector2D& ploc, int shape, int dxj) const
{
	double r = ploc.X();
	double s = ploc.Y();

	if (dxj == 1)
	{
		switch(shape)
		{
			//D_sf/D_r :
		case 1: return  0.25*(1+s);  //Node 1
		case 2: return -0.25*(1+s);  //Node 2
		case 3: return -0.25*(1-s);  //Node 3
		case 4: return  0.25*(1-s);  //Node 4
		default: return 0;
		}
	}
	else
	{
		switch(shape)
		{
			//D_sf/D_s :
		case 1: return  0.25*(1+r);  //Node 1
		case 2: return  0.25*(1-r);  //Node 2
		case 3: return -0.25*(1-r);  //Node 3
		case 4: return -0.25*(1+r);  //Node 4
		default: return 0;
		}
	}

}






//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//nodal coordinates of points in 2D; optional rotation matrix should be included to couple with 3D motion!
//initial velocities!!!!!
Plate2Dquad::Plate2Dquad(MBS* mbsi, int bodyindi, const Vector2D& xc1, const Vector2D& xc2, const Vector2D& xc3, const Vector2D& xc4, 
												 const Vector2D& vc1, const Vector2D& vc2, const Vector2D& vc3, const Vector2D& vc4, 
												 double rhoi, double Emi, double nui, double thickness, const Vector3D& coli):
Plate2D(mbsi)
{
	//order of nodes:
	//
	//          N5
	//   N2 +---+---+ N1
	//      |       |
	//   N6 +   N9  + N8
	//      |       |
	//   N3 +---+---+ N4
	//          N7
	//


	bodyind = bodyindi;
	lz = thickness;

	Vector2D p[10];
	p[1]=Vector2D(xc1(1),xc1(2));
	p[2]=Vector2D(xc2(1),xc2(2));
	p[3]=Vector2D(xc3(1),xc3(2));
	p[4]=Vector2D(xc4(1),xc4(2));

	p[5]=Vector2D(0.5*(p[1]+p[2]));
	p[6]=Vector2D(0.5*(p[2]+p[3]));
	p[7]=Vector2D(0.5*(p[3]+p[4]));
	p[8]=Vector2D(0.5*(p[4]+p[1]));
	p[9]=Vector2D(0.25*(p[1]+p[2]+p[3]+p[4]));

	for (int i=1; i <= 9; i++)
	{
		//global_uo << "Node " << i << "=" << p[i] << "\n";
	}


	//Transformation must be set before call to ToP3D
	Node n[10];
	for (int i=1; i <= NNodes(); i++)
	{
		n[i] = Node(2, bodyind, Vector3D(p[i].X(),p[i].Y(),0.));
		NodeNum(i) = mbsi->AddBodyNode(&n[i]); //check if node exists in body!!!
	}

	mass = lz * rhoi * Area(ToP3D(p[1]),ToP3D(p[2]),ToP3D(p[3]),ToP3D(p[4])); //only approximative

	xg = Vector(NNodes()*Dim()); //for intial conditions??? don't know why?
	xgd = xg; //for initial drawing???

	x_p0 = Vector(NNodes()*Dim()); //store nodal positions for drawing and constraints
	for (int i=1; i <= NNodes(); i++)
	{
		x_p0((i-1)*2+1) = p[i].X(); 
		x_p0((i-1)*2+2) = p[i].Y();
	}

	int nd = NNodes()*Dim();
	x_init = Vector(2*NNodes()*Dim());
	for (int i=1; i<=NNodes()*Dim(); i++) x_init(i) = 0; //initial nodal displacements

	x_init(nd+1) = vc1.X(); x_init(nd+2) = vc1.Y(); //intial nodal velocities
	x_init(nd+3) = vc2.X(); x_init(nd+4) = vc2.Y();
	x_init(nd+5) = vc3.X(); x_init(nd+6) = vc3.Y();
	x_init(nd+7) = vc4.X(); x_init(nd+8) = vc4.Y();

	x_init(nd+ 9) = 0.5*(vc1.X()+vc2.X()); x_init(nd+10) = 0.5*(vc1.Y()+vc2.Y()); //average velocities
	x_init(nd+11) = 0.5*(vc2.X()+vc3.X()); x_init(nd+12) = 0.5*(vc2.Y()+vc3.Y());
	x_init(nd+13) = 0.5*(vc3.X()+vc4.X()); x_init(nd+14) = 0.5*(vc3.Y()+vc4.Y());
	x_init(nd+15) = 0.5*(vc4.X()+vc1.X()); x_init(nd+16) = 0.5*(vc4.Y()+vc1.Y());
	x_init(nd+17) = 0.25*(vc1.X()+vc2.X()+vc3.X()+vc4.X()); 
	x_init(nd+18) = 0.25*(vc1.Y()+vc2.Y()+vc3.Y()+vc4.Y());

	Material mat(GetMBS(), rhoi, Emi, nui, 1, 0);
	AddMaterial(mat);

	//nu = nui;
	//Em = Emi;
	//rho = rhoi;
	col = coli;
	//concentratedmass[0] = concentratedmass[1] = concentratedmass[2] = concentratedmass[3] = 0;

	orderxy = Nip()*2-1; //2x2 integration points
	orderxyM = 4;
	orderxyH = 2;
	BuildDSMatrices();

	elementname = GetElementSpec();
};


Plate2Dquad::Plate2Dquad(MBS* mbsi, int bodyindi, const Vector2D& xc1, const Vector2D& xc2, const Vector2D& xc3, const Vector2D& xc4, 
												 const Vector2D& xc5, const Vector2D& xc6, const Vector2D& xc7, const Vector2D& xc8, const Vector2D& xc9, 
												 const Vector2D& vc1, const Vector2D& vc2, const Vector2D& vc3, const Vector2D& vc4, 
												 const Vector2D& vc5, const Vector2D& vc6, const Vector2D& vc7, const Vector2D& vc8, const Vector2D& vc9, 
												 double rhoi, double Emi, double nui, double thickness, const Vector3D& coli): Plate2D(mbsi)
{
	//order of nodes:
	//
	//          N5
	//   N2 +---+---+ N1
	//      |       |
	//   N6 +   N9  + N8
	//      |       |
	//   N3 +---+---+ N4
	//          N7
	//


	bodyind = bodyindi;
	lz = thickness;

	Vector2D p[10];
	p[1]=xc1;	p[2]=xc2;
	p[3]=xc3;	p[4]=xc4;
	p[5]=xc5;	p[6]=xc6;
	p[7]=xc7;	p[8]=xc8;
	p[9]=xc9;


	//Transformation must be set before call to ToP3D
	Node n[10];
	for (int i=1; i <= NNodes(); i++)
	{
		n[i] = Node(2, bodyind, Vector3D(p[i].X(),p[i].Y(),0.));
		NodeNum(i) = mbsi->AddBodyNode(&n[i]); //check if node exists in body!!!
	}

	mass = lz * rhoi * Area(ToP3D(p[1]),ToP3D(p[2]),ToP3D(p[3]),ToP3D(p[4])); //only approximative

	xg = Vector(NNodes()*Dim()); //for intial conditions??? don't know why?
	xgd = xg; //for initial drawing???

	x_p0 = Vector(NNodes()*Dim()); //store nodal positions for drawing and constraints
	for (int i=1; i <= NNodes(); i++)
	{
		x_p0((i-1)*2+1) = p[i].X(); 
		x_p0((i-1)*2+2) = p[i].Y();
	}

	int nd = NNodes()*Dim();
	x_init = Vector(2*NNodes()*Dim());
	for (int i=1; i<=NNodes()*Dim(); i++) x_init(i) = 0; //initial nodal displacements

	x_init(nd+1) = vc1.X(); x_init(nd+2) = vc1.Y(); //intial nodal velocities
	x_init(nd+3) = vc2.X(); x_init(nd+4) = vc2.Y();
	x_init(nd+5) = vc3.X(); x_init(nd+6) = vc3.Y();
	x_init(nd+7) = vc4.X(); x_init(nd+8) = vc4.Y();

	x_init(nd+ 9) = vc5.X(); x_init(nd+10) = vc5.Y();
	x_init(nd+11) = vc6.X(); x_init(nd+12) = vc6.Y();
	x_init(nd+13) = vc7.X(); x_init(nd+14) = vc7.Y();
	x_init(nd+15) = vc8.X(); x_init(nd+16) = vc8.Y();
	x_init(nd+17) = vc9.X(); x_init(nd+18) = vc9.Y();

	Material mat(GetMBS(), rhoi, Emi, nui, 1, 0);
	AddMaterial(mat);
	//nu = nui;
	//Em = Emi;
	//rho = rhoi;
	col = coli;
	//concentratedmass[0] = concentratedmass[1] = concentratedmass[2] = concentratedmass[3] = 0;

	orderxy = 7;
	orderxyM = 4;
	orderxyH = 2;
	BuildDSMatrices();

	elementname = GetElementSpec();
}

void Plate2Dquad::SetPlate2Dquad(int bodyindi, const TArray<int> &nodelist, int matnr, double thickness, const Vector3D &coli, int CMSelememti)
{
	orderxy = 7;
	orderxyM = 4;
	orderxyH = 2;

	// (AD) make a 9th node here ? NO hack in FEMesh::AddElement and FEMesh::GetPtr...
	Plate2D::SetPlate2D(bodyindi, nodelist, matnr, thickness, coli, CMSelememti);

	elementname = GetElementSpec();
}

void Plate2Dquad::GetS0(Vector& sf, const Vector2D& ploc) const
{
	//Nodes always sorted counterclock-wise in x-y (r-s) plane
	double r = ploc.X();
	double s = ploc.Y();
	sf.SetLen(NS());

	double r2 = r*r;
	double s2 = s*s;

	sf(1) =  r*s/4.0+r2*s/4.0+r2*s2/4.0+s2*r/4.0;
	sf(2) =  -r*s/4.0+r2*s/4.0+r2*s2/4.0-s2*r/4.0;
	sf(3) =  r*s/4.0-s2*r/4.0+r2*s2/4.0-r2*s/4.0;
	sf(4) =  -r*s/4.0-r2*s/4.0+r2*s2/4.0+s2*r/4.0;
	sf(5) =  s/2.0-r2*s/2.0+s2/2.0-r2*s2/2.0;
	sf(6) =  -r/2.0+s2*r/2.0+r2/2.0-r2*s2/2.0;
	sf(7) =  -s/2.0+r2*s/2.0+s2/2.0-r2*s2/2.0;
	sf(8) =  r/2.0-s2*r/2.0+r2/2.0-r2*s2/2.0;
	sf(9) =  (-1.0+r2)*(-1.0+s2);

}

void Plate2Dquad::GetDSMatrix0(const Vector2D& ploc, Matrix& sf) const
{
	double r = ploc.X();
	double s = ploc.Y();
	sf.SetSize(Dim(),NS());

	double r2 = r*r;
	double s2 = s*s;

	//D_sf/D_r:
	sf(1,1) =  s/4.0+r*s/2.0+s2*r/2.0+s2/4.0;
	sf(1,2) =  -s/4.0+r*s/2.0+s2*r/2.0-s2/4.0;
	sf(1,3) =  s/4.0-s2/4.0+s2*r/2.0-r*s/2.0;
	sf(1,4) =  -s/4.0-r*s/2.0+s2*r/2.0+s2/4.0;
	sf(1,5) =  -r*s-s2*r;
	sf(1,6) =  -1.0/2.0+s2/2.0+r-s2*r;
	sf(1,7) =  r*s-s2*r;
	sf(1,8) =  1.0/2.0-s2/2.0+r-s2*r;
	sf(1,9) =  2.0*r*(-1.0+s2);

	//D_sf/D_s:
	sf(2,1) =  r/4.0+r2/4.0+r2*s/2.0+r*s/2.0;
	sf(2,2) =  -r/4.0+r2/4.0+r2*s/2.0-r*s/2.0;
	sf(2,3) =  r/4.0-r*s/2.0+r2*s/2.0-r2/4.0;
	sf(2,4) =  -r/4.0-r2/4.0+r2*s/2.0+r*s/2.0;
	sf(2,5) =  1.0/2.0-r2/2.0+s-r2*s;
	sf(2,6) =  r*s-r2*s;
	sf(2,7) =  -1.0/2.0+r2/2.0+s-r2*s;
	sf(2,8) =  -r*s-r2*s;
	sf(2,9) =  2.0*(-1.0+r2)*s;

}

double Plate2Dquad::GetS0(const Vector2D& ploc, int i) const
{
	//Nodes always sorted counterclock-wise in x-y (r-s) plane
	double r = ploc.X();
	double s = ploc.Y();

	double r2 = r*r;
	double s2 = s*s;

	switch(i)
	{
	case 1: return  r*s/4.0+r2*s/4.0+r2*s2/4.0+s2*r/4.0;
	case 2: return  -r*s/4.0+r2*s/4.0+r2*s2/4.0-s2*r/4.0;
	case 3: return  r*s/4.0-s2*r/4.0+r2*s2/4.0-r2*s/4.0;
	case 4: return  -r*s/4.0-r2*s/4.0+r2*s2/4.0+s2*r/4.0;
	case 5: return  s/2.0-r2*s/2.0+s2/2.0-r2*s2/2.0;
	case 6: return  -r/2.0+s2*r/2.0+r2/2.0-r2*s2/2.0;
	case 7: return  -s/2.0+r2*s/2.0+s2/2.0-r2*s2/2.0;
	case 8: return  r/2.0-s2*r/2.0+r2/2.0-r2*s2/2.0;
	case 9: return  (-1.0+r2)*(-1.0+s2);
	default: return 0;
	}
}

//D S_i / Dx_j
double Plate2Dquad::GetDS0(const Vector2D& ploc, int shape, int dxj) const
{
	double r = ploc.X();
	double s = ploc.Y();

	double r2 = r*r;
	double s2 = s*s;

	if (dxj == 1)
	{
		switch(shape)
		{//D_sf/D_r:
		case 1: return  s/4.0+r*s/2.0+s2*r/2.0+s2/4.0;
		case 2: return  -s/4.0+r*s/2.0+s2*r/2.0-s2/4.0;
		case 3: return  s/4.0-s2/4.0+s2*r/2.0-r*s/2.0;
		case 4: return  -s/4.0-r*s/2.0+s2*r/2.0+s2/4.0;
		case 5: return  -r*s-s2*r;
		case 6: return  -1.0/2.0+s2/2.0+r-s2*r;
		case 7: return  r*s-s2*r;
		case 8: return  1.0/2.0-s2/2.0+r-s2*r;
		case 9: return  2.0*r*(-1.0+s2);
		default: return 0;
		}
	}
	else
	{
		switch(shape)
		{//D_sf/D_s:
		case 1: return  r/4.0+r2/4.0+r2*s/2.0+r*s/2.0;
		case 2: return  -r/4.0+r2/4.0+r2*s/2.0-r*s/2.0;
		case 3: return  r/4.0-r*s/2.0+r2*s/2.0-r2/4.0;
		case 4: return  -r/4.0-r2/4.0+r2*s/2.0+r*s/2.0;
		case 5: return  1.0/2.0-r2/2.0+s-r2*s;
		case 6: return  r*s-r2*s;
		case 7: return  -1.0/2.0+r2/2.0+s-r2*s;
		case 8: return  -r*s-r2*s;
		case 9: return  2.0*(-1.0+r2)*s;
		default: return 0;
		}
	}


}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//nodal coordinates of points in 2D; optional rotation matrix should be included to couple with 3D motion!
//initial velocities!!!!!
Trig2Dlin::Trig2Dlin(MBS* mbsi, int bodyindi, const Vector2D& xc1, const Vector2D& xc2, const Vector2D& xc3, 
										 const Vector2D& vc1, const Vector2D& vc2, const Vector2D& vc3, 
										 double rhoi, double Emi, double nui, double thickness, const Vector3D& coli):
Plate2D(mbsi)
{
	bodyind = bodyindi;
	Vector2D p1(xc1(1),xc1(2));
	Vector2D p2(xc2(1),xc2(2));
	Vector2D p3(xc3(1),xc3(2));

	lz = thickness;

	//Transformation must be set before call to ToP3D
	Node n1(2, bodyind, Vector3D(p1.X(),p1.Y(),0.));
	Node n2(2, bodyind, Vector3D(p2.X(),p2.Y(),0.));
	Node n3(2, bodyind, Vector3D(p3.X(),p3.Y(),0.));
	node[0] = mbsi->AddBodyNode(&n1); //check if node exists in body!!!
	node[1] = mbsi->AddBodyNode(&n2); //add this function to mbs.h, add bodyindex to class Node
	node[2] = mbsi->AddBodyNode(&n3);

	mass = lz * rhoi * Area(ToP3D(p1),ToP3D(p2),ToP3D(p3));

	xg = Vector(NNodes()*Dim()); //for intial conditions??? don't know why?
	xgd = xg; //for initial drawing???

	x_p0 = Vector(NNodes()*Dim()); //store nodal positions for drawing and constraints
	x_p0(1) = xc1.X(); x_p0(2) = xc1.Y();
	x_p0(3) = xc2.X(); x_p0(4) = xc2.Y();
	x_p0(5) = xc3.X(); x_p0(6) = xc3.Y();

	x_init = Vector(2*NNodes()*Dim());
	for (int i=1; i<=NNodes()*Dim(); i++) x_init(i) = 0; //initial nodal displacements

	x_init( 7) = vc1.X(); x_init( 8) = vc1.Y(); //intial nodal velocities
	x_init( 9) = vc2.X(); x_init(10) = vc2.Y();
	x_init(11) = vc3.X(); x_init(12) = vc3.Y();

	Material mat(GetMBS(), rhoi, Emi, nui, 1, 0);
	AddMaterial(mat);
	//nu = nui;
	//Em = Emi;
	//rho = rhoi;
	col = coli;
	//concentratedmass[0] = concentratedmass[1] = concentratedmass[2] = concentratedmass[3] = 0;

	orderxy = 2; //all derivatives of the shape functions are linear, should be actually enought to set = 1
	orderxyM = 2;
	orderxyH = 2;

	BuildDSMatrices();
};

void Trig2Dlin::SetTrig2Dlin(int bodyindi, const TArray<int> &nodelist, int matnr, double thickness, const Vector3D &coli, int CMSelememti)
{
	orderxy = 2; //all derivatives of the shape functions are linear, should be actually enought to set = 1
	orderxyM = 2;
	orderxyH = 2;

	Plate2D::SetPlate2D(bodyindi, nodelist, matnr, thickness, coli, CMSelememti);

	BuildDSMatrices();
}

void Trig2Dlin::GetS0(Vector& sf, const Vector2D& ploc) const
{
	//Nodes always sorted counterclock-wise in x-y (r-s) plane
	double r = ploc.X();
	double s = ploc.Y();
	double t = 1 - r - s;
	sf.SetLen(NS());

	//order of nodes:
	//
	//  node 2 +
	//         |\ 
	//         | \
	//         |  \ 
	//  node 3 +---+ node 1
	//

	sf(1) = r;  //Node 1
	sf(2) = s;  //Node 2
	sf(3) = t;  //Node 3
}

void Trig2Dlin::GetDSMatrix0(const Vector2D& ploc, Matrix& sf) const
{
	double r = ploc.X();
	double s = ploc.Y();
	sf.SetSize(Dim(),NS());

	//D_sf/D_r:
	sf(1,1) = 1;  //Node 1
	sf(1,2) = 0;  //Node 2
	sf(1,3) =-1;  //Node 3

	//D_sf/D_s:
	sf(2,1) = 0;  //Node 1
	sf(2,2) = 1;  //Node 2
	sf(2,3) =-1;  //Node 3

}

double Trig2Dlin::GetS0(const Vector2D& ploc, int i) const
{
	//Nodes always sorted counterclock-wise in x-y (r-s) plane
	double r = ploc.X();
	double s = ploc.Y();

	switch(i)
	{
	case 1: return r;  //Node 1
	case 2: return s;  //Node 2
	case 3: return 1 - r - s;  //Node 3
	default: return 0;
	}
}

//D S_i / Dx_j
double Trig2Dlin::GetDS0(const Vector2D& ploc, int shape, int dxj) const
{
	double r = ploc.X();
	double s = ploc.Y();

	if (dxj == 1)
	{
		switch(shape)
		{
			//D_sf/D_r :
		case 1: return  1;  //Node 1
		case 2: return  0;  //Node 2
		case 3: return -1;  //Node 3
		default: return 0;
		}
	}
	else
	{
		switch(shape)
		{
			//D_sf/D_s :
		case 1: return  0;  //Node 1
		case 2: return  1;  //Node 2
		case 3: return -1;  //Node 3
		default: return 0;
		}
	}

}


void Trig2Dlin::GetH(Matrix& H) 
{
	if (Hmatrix.Getrows() == FlexDOF())
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
		jac.SetSize(dim,dim);

		DS.SetSize(dim,ns);
		SV.SetLen(ns);

		if (IsTrig())
		{
			GetIntegrationRuleTrig(x1,x2,w1,orderxy);
			w2.SetLen(1); //only dummy
		}
		else
		{
			GetIntegrationRule(x1,w1,orderxy);
			GetIntegrationRule(x2,w2,orderxy);
		}

		for (int i1=1; i1<=x1.GetLen(); i1++)
		{
			for (int i2=1; i2<=w2.GetLen(); i2++)
			{
				Vector2D p;
				if (!IsTrig())
					p = Vector2D(x1(i1),x2(i2));
				else
					p = Vector2D(x1(i1),x2(i1)); //i2 always 1

				GetS0(SV,p);

				GetDSMatrix0(p,DS);
				GetJacobi(jac,p);
				//GetJacobi(jac,p,DS);
				double jacdet = jac.Det();
				double fact = fabs (jacdet) * w1(i1)*lz; 

				for (int i=0; i<ns; i++)
				{
					for (int j=1; j<=dim; j++)
					{
						H(i*dim+j,j)+=fact*SV(i+1);
					}
				}
			}
		}
		Hmatrix = H;
	}
}

void Trig2Dlin::EvalM(Matrix& m, double t) 
{
	if (massmatrix.Getcols() == FlexDOF())
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
		jac.SetSize(dim,dim);

		Matrix HL(ns*dim,dim);
		DS.SetSize(dim,ns);

		if (IsTrig())
		{
			GetIntegrationRuleTrig(x1,x2,w1,orderxy);
			w2.SetLen(1); //only dummy
		}
		else
		{
			GetIntegrationRule(x1,w1,orderxy);
			GetIntegrationRule(x2,w2,orderxy);
		}

		for (int i1=1; i1<=x1.GetLen(); i1++)
		{
			for (int i2=1; i2<=w2.GetLen(); i2++)
			{
				Vector2D p;
				if (!IsTrig())
					p = Vector2D(x1(i1),x2(i2));
				else
					p = Vector2D(x1(i1),x2(i1)); //i2 always 1

				GetS0(SV,p);

				for (int i=0; i<ns; i++)
				{
					for (int j=1; j<=dim; j++)
					{
						HL(i*dim+j,j)=SV(i+1);
					}
				}

				GetJacobi(jac,p);

				double jacdet = jac.Det();
				m += fabs (jacdet) * Rho() * w1(i1) * lz * (HL*HL.GetTp());
			}
		}

		massmatrix = m;

	}
};


void Trig2Dlin::EvalF2(Vector& f, double t) 
{
	Body2D::EvalF2(f,t);
	TMStartTimer(22);

	int sos = FlexDOF();
	SetComputeCoordinates();
	static Vector fadd;


	int ns = NS();
	int dim = Dim();
	//double u;
	double nu = Nu();

	double la=Em() * nu / ((1.+nu)*(1.-2.*nu));
	double mu=Em() / 2. / (1.+nu);

	Vector3D sigma;
	double k = Em()/(1.-nu*nu);
	Matrix3D C(k,k*nu,0,k*nu,k,0,0,0,k*(1.-nu)/2.); //Bathe p. 229, mit Gleitung!!!

	Matrix3D strain, piola1, F;
	strain.SetSize(2,2);
	piola1.SetSize(2,2);
	F.SetSize(2,2);

	temp.SetLen(SOS());
	fadd.SetLen(SOS());
	fadd.SetAll(0);


	if (IsTrig())
	{
		GetIntegrationRuleTrig(x1,x2,w1,orderxy);
		w2.SetLen(1); //only dummy
	}
	else
	{
		GetIntegrationRule(x1,w1,orderxy);
		GetIntegrationRule(x2,w2,orderxy);
	}

	int kx1 = w2.Length();
	for (int i1=1; i1<=x1.GetLen(); i1++)
	{
		for (int i2=1; i2<=w2.GetLen(); i2++)
		{
			Vector2D p;
			if (!IsTrig())
				p = Vector2D(x1(i1),x2(i2));
			else
				p = Vector2D(x1(i1),x2(i1)); //i2 always 1


			//int_A (P : d F / d q) dA
			//P = F * S
			//S = D : E
			//E = 1/2*(F^T F-I)

			int i,j,k;
			int ind = (i1-1)*kx1+(i2-1);

			// compute F 
			F.SetAll(0);
			int l;

			for (j = 1; j <= dim; j++) 
			{
				for (i = 1; i <= ns; i++)
				{
					l = (i-1)*dim+j;
					for (k = 1; k <= dim; k++)
					{
						F(j,k) += grad[ind](k,i)*xg(l);
					}
				}
				F(j,j) += 1;
			}

			// Green-Lagrange strain tensor
			strain = 0.5 * (F.GetTp() * F);
			strain(1,1) -= 0.5; strain(2,2) -= 0.5;

			Vector3D eps(strain(1,1), strain(2,2), 2.*strain(1,2));
			sigma = C*eps;

			piola1(1,1) = sigma(1);
			piola1(2,2) = sigma(2);
			piola1(1,2) = sigma(3);
			piola1(2,1) = sigma(3);

			piola1 = F * (piola1);

			//linear:
			//Matrix3D G=F-Matrix3D(1);
			//strain = 0.5*(G+G.GetTp());
			//piola1 = ((2*mu) * strain + Matrix3D(la * strain.Trace()));

			for (int j=1; j <= Dim(); j++)
			{
				for (int i = 0; i < ns; i++)
				{
					temp(2*i+j) = grad[ind](1, i+1)*piola1(j,1)
						+ grad[ind](2, i+1)*piola1(j,2);
				}
			}

			fadd.MultAdd(fabs (jacdet[ind]) * w1(i1) * lz,temp); 
		} 
	}

	f -= fadd;
	TMStopTimer(22);

	if (GetMassDamping() != 0)
	{
		//UO() << "Warning: damping only includes diagonal of mass matrix!!!\n";
		// +++++++ damping: +++++++
		for (int i = 1; i <= SOS(); i++)
			xg(i) = XGP(i);

		if (massmatrix.Getcols() == SOS())
		{
			Mult(massmatrix,xg,temp);
			/*
			for (int i = 1; i <= SOS(); i++)
			temp(i) = xg(i)*massmatrix(i,i);
			*/
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
}; 

void Trig2Dlin::DrawElement() 
{
	mbs->SetColor(col);
	//return;

	SetDrawCoordinates();

	int linemode = 1; //0=no lines, 1=outline+color, 2=outline, 3=elementline+color, 4=elementline
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

	int res = GetMBS()->GetDrawResolution()-1;

	if (colormode == 0 && NS() == 3) 
	{
		//linear element shall not be interpolated without stresses
		res = 0;
	}
	if (res < 0) res = 0;

	double tilex = pow(2.,res);

	double tiley = tilex;

	TArray<Vector3D> points((int)(tilex+1)*(int)(tiley+1));
	TArray<double> vals((int)(tilex+1)*(int)(tiley+1));
	double v=0;

	points.SetLen(0); vals.SetLen(0);
	Vector2D p0, vx, vy;
	int tileyn = (int)tiley;

	double shrink = GetMBS()->GetDOption(106);//0.995;
	p0 = Vector2D(0.5*(1.-shrink),0.5*(1.-shrink));
	vx = Vector2D(1./tilex*shrink,0);
	vy = Vector2D(0,1./tiley*shrink);

	//UO() << "vx=" << vx << "\n";
	//UO() << "vy=" << vy << "\n";

	for (double iy = 0; iy <= tileyn+1e-10; iy++)
	{
		for (double ix = 0; ix <= tilex+1e-10; ix++)
		{
			Vector2D ploc = (p0+ix*vx*(1.-iy/tiley)+iy*vy);
			//UO() << "ploc=" << ploc << "\n";
			Vector3D pg = ToP3D(GetPos2DD(ploc));
			//UO() << "pg=" << pg << "\n";
			points.Add(pg);
			if (colormode)
				v = GetFieldVariableValue(*GetMBS()->GetActualPostProcessingFieldVariable(), ploc, true);
			vals.Add(v);
		}
	}

	mbs->DrawColorQuads(points,vals,(int)tilex+1,(int)tileyn+1,colormode,linemode);

};





//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//nodal coordinates of points in 2D; optional rotation matrix should be included to couple with 3D motion!
//initial velocities!!!!!
Trig2Dquad::Trig2Dquad(MBS* mbsi, int bodyindi, Vector2D* xc, Vector2D* vc, int interpolate_points,
		double rhoi, double Emi, double nui, double thickness, const Vector3D& coli): Trig2Dlin(mbsi)
{
	bodyind = bodyindi;
	Vector2D p[6];
	Vector2D v[6];

	if (interpolate_points)
	{
		for (int i = 0; i < 3; i++)
		{
			p[i] = xc[i];
			p[i+3] = 0.5*(xc[i]+xc[(i+1)%3]);
			v[i] = vc[i];
			v[i+3] = 0.5*(vc[i]+vc[(i+1)%3]);
		}
	}
	else
	{
		for (int i = 0; i < 6; i++)
		{
			p[i] = xc[i];
			v[i] = vc[i];
		}
	}

	InitializeElement(bodyindi, p, v, rhoi, Emi, nui, thickness, coli);
};

void Trig2Dquad::InitializeElement(int bodyindi, Vector2D* p, Vector2D* v,
		double rhoi, double Emi, double nui, double thickness, const Vector3D& coli)
{
	lz = thickness;

	//Transformation must be set before call to ToP3D
	for (int i=1; i <= NS(); i++)
	{
		Node n1(2, bodyind, Vector3D(p[i-1].X(), p[i-1].Y(),0.));
		node[i-1] = mbs->AddBodyNode(&n1); //check if node exists in body!!!
	}

	mass = lz * rhoi * Area(ToP3D(p[0]),ToP3D(p[1]),ToP3D(p[2]));

	xg = Vector(NNodes()*Dim()); //for intial conditions??? don't know why?
	xgd = xg; //for initial drawing???

	x_p0 = Vector(NNodes()*Dim()); //store nodal positions for drawing and constraints
	for (int i=1; i <= NS(); i++)
	{
		x_p0(i*2-1) = p[i-1].X();
		x_p0(i*2-0) = p[i-1].Y();
	}

	x_init = Vector(2*NNodes()*Dim());
	for (int i=1; i<=NNodes()*Dim(); i++) x_init(i) = 0; //initial nodal displacements

	for (int i=1; i <= NS(); i++)
	{
		x_init(i*2-1+2*NS()) = v[i-1].X();
		x_init(i*2-0+2*NS()) = v[i-1].Y();
	}

	Material mat(GetMBS(), rhoi, Emi, nui, 1, 0);
	AddMaterial(mat);

	//nu = nui;
	//Em = Emi;
	//rho = rhoi;
	col = coli;
	//concentratedmass[0] = concentratedmass[1] = concentratedmass[2] = concentratedmass[3] = 0;

	orderxy = 4; //all derivatives of the shape functions are linear, should be actually enought to set = 4
	orderxyM = 4;
	orderxyH = 3;

	BuildDSMatrices();
}

void Trig2Dquad::SetTrig2Dquad(int bodyindi, const TArray<int> &nodelist, int matnr, double thickness, const Vector3D &coli, int CMSelememti)
{
	orderxy = 4; //all derivatives of the shape functions are linear, should be actually enought to set = 1
	orderxyM = 4;
	orderxyH = 3;

	Plate2D::SetPlate2D(bodyindi, nodelist, matnr, thickness, coli, CMSelememti);

	BuildDSMatrices();
}

void Trig2Dquad::GetS0(Vector& sf, const Vector2D& ploc) const
{
	//Nodes always sorted counterclock-wise in x-y (r-s) plane
	double r = ploc.X();
	double s = ploc.Y();
	double t = 1 - r - s;
	sf.SetLen(NS());

	//order of nodes:
	//
	//  node 2 +
	//         |\ 
	//  node 5 + + node 4
	//         |  \ 
	//  node 3 +-+-+ node 1
	//         node 6

	sf(1) = r*(2.*r-1.);  //Node 1
	sf(2) = s*(2.*s-1.);  //Node 2
	sf(3) = t*(2.*t-1.);  //Node 3
	sf(4) = 4.*r*s;      //Node 4
	sf(5) = 4.*s*t;      //Node 5
	sf(6) = 4.*r*t;      //Node 6
}

void Trig2Dquad::GetDSMatrix0(const Vector2D& ploc, Matrix& sf) const
{
	double r = ploc.X();
	double s = ploc.Y();
	sf.SetSize(Dim(),NS());

	//D_sf/D_r:
	sf(1,1) = 4.0*r-1.0;
	sf(1,2) = 0.0;
	sf(1,3) = -3.0+4.0*r+4.0*s;
	sf(1,4) = 4.0*s;
	sf(1,5) = -4.0*s;
	sf(1,6) = 4.0-8.0*r-4.0*s;

	//D_sf/D_s:
	sf(2,1) = 0.0;
	sf(2,2) = 4.0*s-1.0;
	sf(2,3) = -3.0+4.0*r+4.0*s;
	sf(2,4) = 4.0*r;
	sf(2,5) = 4.0-4.0*r-8.0*s;
	sf(2,6) = -4.0*r;

}

double Trig2Dquad::GetS0(const Vector2D& ploc, int i) const
{
	//Nodes always sorted counterclock-wise in x-y (r-s) plane
	double r = ploc.X();
	double s = ploc.Y();
	double t = 1 - r - s;

	switch(i)
	{
	case 1: return r*(2.*r-1.);  //Node 1
	case 2: return s*(2.*s-1.);  //Node 2
	case 3: return t*(2.*t-1.);  //Node 3
	case 4: return 4.*r*s;       //Node 4
	case 5: return 4.*s*t;       //Node 5
	case 6: return 4.*r*t;       //Node 6
	default: return 0;
	}
}

//D S_i / Dx_j
double Trig2Dquad::GetDS0(const Vector2D& ploc, int shape, int dxj) const
{
	double r = ploc.X();
	double s = ploc.Y();

	if (dxj == 1)
	{
		switch(shape)
		{
			//D_sf/D_r :
		case 1: return 4.0*r-1.0;
		case 2: return 0.0;
		case 3: return -3.0+4.0*r+4.0*s;
		case 4: return 4.0*s;
		case 5: return -4.0*s;
		case 6: return 4.0-8.0*r-4.0*s;
		default: return 0;
		}
	}
	else
	{
		switch(shape)
		{
			//D_sf/D_s :
		case 1: return 0.0;
		case 2: return 4.0*s-1.0;
		case 3: return -3.0+4.0*r+4.0*s;
		case 4: return 4.0*r;
		case 5: return 4.0-4.0*r-8.0*s;
		case 6: return -4.0*r;
		default: return 0;
		}
	}

}






//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//nine-node plate element:
//nodal coordinates of points in 2D; optional rotation matrix should be included to couple with 3D motion!
//initial velocities!!!!!
Plate2DquadFFRF::Plate2DquadFFRF(MBS* mbsi, int FFRFindex, int bodyindi, const Vector2D& xc1, const Vector2D& xc2, const Vector2D& xc3, const Vector2D& xc4, 
																 const Vector2D& vc1, const Vector2D& vc2, const Vector2D& vc3, const Vector2D& vc4, 
																 double rhoi, double Emi, double nui, double thickness, 
																 const Vector3D& coli, int isCMSi):Plate2Dquad(mbsi), 
																 Sbar_tilde(), Sbar_tildeSM(), K(),	SbarS(), SbarSM(), Ibar11S(), Ibar12S(), Ibar21S(), Ibar22S(), I1S()
{
	//Add FFRF reference
	FFRFind = FFRFindex;
	if (isCMSi) AddType(TCMS);

	bodyind = bodyindi;
	lz = thickness;

	Vector2D p[10];
	p[1]=Vector2D(xc1(1),xc1(2));
	p[2]=Vector2D(xc2(1),xc2(2));
	p[3]=Vector2D(xc3(1),xc3(2));
	p[4]=Vector2D(xc4(1),xc4(2));

	p[5]=Vector2D(0.5*(p[1]+p[2]));
	p[6]=Vector2D(0.5*(p[2]+p[3]));
	p[7]=Vector2D(0.5*(p[3]+p[4]));
	p[8]=Vector2D(0.5*(p[4]+p[1]));
	p[9]=Vector2D(0.25*(p[1]+p[2]+p[3]+p[4]));

	//Transformation must be set before call to ToP3D
	Node n[10];
	for (int i=1; i <= NNodes(); i++)
	{
		n[i] = Node(2, bodyind, Vector3D(p[i].X(),p[i].Y(),0.));
		if (!IsCMS())
		{
			NodeNum(i) = mbsi->AddBodyNode(&n[i]); //check if node exists in body!!!
		}
		else
		{
			n[i].SetAuxNode();
			NodeNum(i) = ReferenceFrame().AddNode(&n[i]); //check if node exists in body!!!
			//global_uo << "n" << NodeNum(i) << "=" << n[i].Pos() << "\n";
		}

	}

	int nnn = ReferenceFrame().NNodes();

	mass = lz * rhoi * Area(ToP3D(p[1]),ToP3D(p[2]),ToP3D(p[3]),ToP3D(p[4])); //only approximative

	xg = Vector(NNodes()*Dim()+FFRFDim()); //for intial conditions??? don't know why?
	xgd = xg; //for initial drawing???

	x_p0 = Vector(NNodes()*Dim()); //store nodal positions for drawing and constraints
	for (int i=1; i <= NNodes(); i++)
	{
		x_p0((i-1)*2+1) = p[i].X(); 
		x_p0((i-1)*2+2) = p[i].Y();
	}

	if (!IsCMS())
	{
		int nd = NNodes()*Dim()+FFRFDim();
		x_init = Vector(2*(NNodes()*Dim()+FFRFDim()));
		for (int i=1; i<=NNodes()*Dim(); i++) x_init(i) = 0; //initial nodal displacements

		//initial velocities:
		x_init(nd+1) = vc1.X(); x_init(nd+2) = vc1.Y(); //intial nodal velocities
		x_init(nd+3) = vc2.X(); x_init(nd+4) = vc2.Y();
		x_init(nd+5) = vc3.X(); x_init(nd+6) = vc3.Y();
		x_init(nd+7) = vc4.X(); x_init(nd+8) = vc4.Y();

		x_init(nd+ 9) = 0.5*(vc1.X()+vc2.X()); x_init(nd+10) = 0.5*(vc1.Y()+vc2.Y()); //average velocities
		x_init(nd+11) = 0.5*(vc2.X()+vc3.X()); x_init(nd+12) = 0.5*(vc2.Y()+vc3.Y());
		x_init(nd+13) = 0.5*(vc3.X()+vc4.X()); x_init(nd+14) = 0.5*(vc3.Y()+vc4.Y());
		x_init(nd+15) = 0.5*(vc4.X()+vc1.X()); x_init(nd+16) = 0.5*(vc4.Y()+vc1.Y());
		x_init(nd+17) = 0.25*(vc1.X()+vc2.X()+vc3.X()+vc4.X()); 
		x_init(nd+18) = 0.25*(vc1.Y()+vc2.Y()+vc3.Y()+vc4.Y());

		//reference frame:
		x_init(NS()*Dim()+1) = ReferenceFrame().GetXInit()(1);
		x_init(NS()*Dim()+2) = ReferenceFrame().GetXInit()(2);
		x_init(NS()*Dim()+3) = ReferenceFrame().GetXInit()(3);

		x_init(2*NS()*Dim()+3+1) = ReferenceFrame().GetXInit()(4);
		x_init(2*NS()*Dim()+3+2) = ReferenceFrame().GetXInit()(5);
		x_init(2*NS()*Dim()+3+3) = ReferenceFrame().GetXInit()(6);
	}
	else
	{
		x_init.SetLen(0);
	}
	//UO() << "x_init=" << x_init << "\n";

	Material mat(GetMBS(), rhoi, Emi, nui, 1, 0);
	AddMaterial(mat);
	//nu = nui;
	//Em = Emi;
	//rho = rhoi;
	col = coli;
	//concentratedmass[0] = concentratedmass[1] = concentratedmass[2] = concentratedmass[3] = 0;

	orderxy = 7; //only for nonlinear
	orderxyM = 4;
	orderxyH = 2;

	if (!ReferenceFrame().IsACRS())
	{
		BuildDSMatrices();
	}
}

Plate2DquadFFRF::Plate2DquadFFRF(MBS* mbsi, int FFRFindex, int bodyindi, const Vector2D& xc1, const Vector2D& xc2, const Vector2D& xc3, const Vector2D& xc4, 
																 const Vector2D& xc5, const Vector2D& xc6, const Vector2D& xc7, const Vector2D& xc8, const Vector2D& xc9, 
																 const Vector2D& vc1, const Vector2D& vc2, const Vector2D& vc3, const Vector2D& vc4, 
																 const Vector2D& vc5, const Vector2D& vc6, const Vector2D& vc7, const Vector2D& vc8, const Vector2D& vc9, 
																 double rhoi, double Emi, double nui, double thickness, const Vector3D& coli, int isCMSi):Plate2Dquad(mbsi), 
																 Sbar_tilde(), Sbar_tildeSM(), K(),	SbarS(), SbarSM(), Ibar11S(), Ibar12S(), Ibar21S(), Ibar22S(), I1S()
{
	//Add FFRF reference
	FFRFind = FFRFindex;
	if (isCMSi) AddType(TCMS);

	bodyind = bodyindi;
	lz = thickness;

	Vector2D p[10];
	p[1]=xc1; p[2]=xc2; p[3]=xc3; p[4]=xc4; 
	p[5]=xc5; p[6]=xc6; p[7]=xc7; p[8]=xc8; p[9]=xc9;

	//Transformation must be set before call to ToP3D
	Node n[10];
	for (int i=1; i <= NNodes(); i++)
	{
		n[i] = Node(2, bodyind, Vector3D(p[i].X(),p[i].Y(),0.));
		if (!IsCMS())
		{
			NodeNum(i) = mbsi->AddBodyNode(&n[i]); //check if node exists in body!!!
		}
		else
		{
			n[i].SetAuxNode();
			NodeNum(i) = ReferenceFrame().AddNode(&n[i]); //check if node exists in body!!!
			//global_uo << "n" << NodeNum(i) << "=" << n[i].Pos() << "\n";
		}

	}

	mass = lz * rhoi * Area(ToP3D(p[1]),ToP3D(p[2]),ToP3D(p[3]),ToP3D(p[4])); //only approximative

	xg = Vector(NNodes()*Dim()+FFRFDim()); //for intial conditions??? don't know why?
	xgd = xg; //for initial drawing???

	x_p0 = Vector(NNodes()*Dim()); //store nodal positions for drawing and constraints
	for (int i=1; i <= NNodes(); i++)
	{
		x_p0((i-1)*2+1) = p[i].X(); 
		x_p0((i-1)*2+2) = p[i].Y();
	}

	if (!IsCMS())
	{
		int nd = NNodes()*Dim()+FFRFDim();
		x_init = Vector(2*(NNodes()*Dim()+FFRFDim()));
		for (int i=1; i<=NNodes()*Dim(); i++) x_init(i) = 0; //initial nodal displacements

		//initial velocities:
		x_init(nd+1) = vc1.X(); x_init(nd+2) = vc1.Y(); //intial nodal velocities
		x_init(nd+3) = vc2.X(); x_init(nd+4) = vc2.Y();
		x_init(nd+5) = vc3.X(); x_init(nd+6) = vc3.Y();
		x_init(nd+7) = vc4.X(); x_init(nd+8) = vc4.Y();

		x_init(nd+ 9) = vc5.X(); x_init(nd+10) = vc5.Y();
		x_init(nd+11) = vc6.X(); x_init(nd+12) = vc6.Y();
		x_init(nd+13) = vc7.X(); x_init(nd+14) = vc7.Y();
		x_init(nd+15) = vc8.X(); x_init(nd+16) = vc8.Y();
		x_init(nd+17) = vc9.X(); x_init(nd+18) = vc9.Y();

		//reference frame:
		x_init(NS()*Dim()+1) = ReferenceFrame().GetXInit()(1);
		x_init(NS()*Dim()+2) = ReferenceFrame().GetXInit()(2);
		x_init(NS()*Dim()+3) = ReferenceFrame().GetXInit()(3);

		x_init(2*NS()*Dim()+3+1) = ReferenceFrame().GetXInit()(4);
		x_init(2*NS()*Dim()+3+2) = ReferenceFrame().GetXInit()(5);
		x_init(2*NS()*Dim()+3+3) = ReferenceFrame().GetXInit()(6);
	}
	else
	{
		x_init.SetLen(0);
	}
	//UO() << "x_init=" << x_init << "\n";

	Material mat(GetMBS(), rhoi, Emi, nui, 1, 0);
	AddMaterial(mat);
	//nu = nui;
	//Em = Emi;
	//rho = rhoi;
	col = coli;
	//concentratedmass[0] = concentratedmass[1] = concentratedmass[2] = concentratedmass[3] = 0;

	orderxy = 7; //only for nonlinear
	orderxyM = 4;
	orderxyH = 2;

	if (!ReferenceFrame().IsACRS())
	{
		BuildDSMatrices();
	}

}


void Plate2DquadFFRF::Initialize() 
{
	Plate2Dquad::Initialize();

	StiffnessMatrix(K);

	//+++++++++++++++++++++++++++++++++++++++++
	//compute mass matrix --> stored in "massmatrix"
	Matrix tmp;
	EvalMff(tmp,0);

	if (!ReferenceFrame().IsACRS())
	{
		//stored functions:
		GetSbarkl(1,2,Sbar_tilde);
		GetSbarkl(2,1,tmp);
		Sbar_tilde -= tmp; //Sbar12-Sbar21


		GetSbar(SbarS);
		GetIbarkl(1,1,Ibar11S);
		GetIbarkl(1,2,Ibar12S);
		GetIbarkl(2,1,Ibar21S);
		GetIbarkl(2,2,Ibar22S);
		GetI1(I1S);

		I1122S = (GetIkl(1,1)+GetIkl(2,2));

		Sbar_tildeSM.CopyFrom(Sbar_tilde);
		SbarSM.CopyFrom(SbarS);
	}
	//UO() << "SbarSM=" << SbarSM << "\n";
	//UO() << "Sbar~SM=" << Sbar_tildeSM << "\n";
}

void Plate2DquadFFRF::LinkToElements()
{
	//order of DOF in LTG: 0 .. Dim()*NS() Xref Yref phi

	if (IsCMS())
	{
		if (IsType(TCMSflag))
		{
			//do only once during DoModalAnalysis

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
				const Node& node = ReferenceFrame().GetNode(j);
				for (int i=1; i <= node.SOS(); i++)
				{
					AddLTG(node.Get(i+node.SOS()));
				}
			}
		}
	}
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

			//Position:
			for (int j = 1; j <= NNodes(); j++)
			{
				const Node& node = GetMBS()->GetNode(NodeNum(j));
				for (int i=1; i <= node.SOS(); i++)
				{
					AddLTG(node.Get(i));
				}
			}
			//reference frame:
			for (int i = 1; i <= ReferenceFrame().SOS(); i++)
			{
				//UO() << "RF-ltg" << i << "=" << ReferenceFrame().LTG(i) << "\n";
				AddLTG(ReferenceFrame().LTG(i));
			}
			//remaining if exists:
			for (int i=1; i <= Nnonode; i++)
			{
				AddLTG(storeltg(i));
			}

			//+++++++++++++++++++++++++++++
			//velocities:
			for (int j = 1; j <= NNodes(); j++)
			{
				const Node& node = GetMBS()->GetNode(NodeNum(j));
				for (int i=1; i <= node.SOS(); i++)
				{
					AddLTG(node.Get(i+node.SOS()));
				}
			}
			//reference frame
			for (int i = 1; i <= ReferenceFrame().SOS(); i++)
			{
				AddLTG(ReferenceFrame().LTG(i+ReferenceFrame().SOS()));
			}
			//remaining if exists:
			for (int i=1; i <= Nnonode; i++)
			{
				AddLTG(storeltg(i+Nnonode));
			}
		}
	}
}

Vector2D Plate2DquadFFRF::GetPos2D(const Vector2D& p_loc) const
{
	if (ReferenceFrame().IsACRS())
		return GetPos2Drel(p_loc);
	else
		return ReferenceFrame().GetRefPos2D()+ReferenceFrame().GetRotMatrix2D()*GetPos2Drel(p_loc);
}

//-1..+1 based!!!
Vector2D Plate2DquadFFRF::GetVel2D(const Vector2D& p_loc) const
{
	if (ReferenceFrame().IsACRS())
		return GetVel2Drel(p_loc);
	else
		return ReferenceFrame().GetRefVel2D()+ReferenceFrame().GetRotMatrix2D()*GetVel2Drel(p_loc)+
		ReferenceFrame().GetRotMatrix2DP()*GetPos2Drel(p_loc);
}

//-1..+1 based!!!
Vector2D Plate2DquadFFRF::GetPos2DD(const Vector2D& p_loc, int use_magnification) const
{
	if (ReferenceFrame().IsACRS())
		return GetPos2DrelD(p_loc, use_magnification);
	else
		return ReferenceFrame().GetRefPos2DD()+ReferenceFrame().GetRotMatrix2DD()*GetPos2DrelD(p_loc, use_magnification);
}

//-1..+1 based!!!
Vector2D Plate2DquadFFRF::GetVel2DD(const Vector2D& p_loc) const
{
	if (ReferenceFrame().IsACRS())
		return GetVel2DrelD(p_loc);
	else
		return ReferenceFrame().GetRefVel2DD()+ReferenceFrame().GetRotMatrix2DD()*GetVel2DrelD(p_loc)+
		ReferenceFrame().GetRotMatrix2DPD()*GetPos2DrelD(p_loc, 0);
}

Vector2D Plate2DquadFFRF::GetNodePos2D(int i) const
{
	if (ReferenceFrame().IsACRS())
		return Vector2D(XG((i-1)*2+1)+x_p0((i-1)*2+1), XG((i-1)*2+2)+x_p0((i-1)*2+2));
	else
		return ReferenceFrame().GetRefPos2D()+ReferenceFrame().GetRotMatrix2D()*Vector2D(XG((i-1)*2+1)+x_p0((i-1)*2+1), XG((i-1)*2+2)+x_p0((i-1)*2+2));
}

Vector2D Plate2DquadFFRF::GetNodePos2DD(int i) const 
{
	return GetPos2DD(GetNodeLocPos2D(i), 1); //always use magnification
/*
	if (ReferenceFrame().IsACRS())
		return Vector2D(XGD((i-1)*2+1)+x_p0((i-1)*2+1), XGD((i-1)*2+2)+x_p0((i-1)*2+2));
	else
		return ReferenceFrame().GetRefPos2DD()+ReferenceFrame().GetRotMatrix2DD()*Vector2D(XGD((i-1)*2+1)+x_p0((i-1)*2+1), XGD((i-1)*2+2)+x_p0((i-1)*2+2));
*/
}

Vector2D Plate2DquadFFRF::GetNodeVel2D(int i) const
{
	if (ReferenceFrame().IsACRS())
		return Vector2D(XGP((i-1)*2+1), XGP((i-1)*2+2));
	else
		return ReferenceFrame().GetRefVel2DD()+ReferenceFrame().GetRotMatrix2DD()*Vector2D(XGP((i-1)*2+1), XGP((i-1)*2+2))+
		ReferenceFrame().GetRotMatrix2DPD()*Vector2D(XG((i-1)*2+1)+x_p0((i-1)*2+1), XG((i-1)*2+2)+x_p0((i-1)*2+2));
}


//insert all 9 entries of mass matrix
void Plate2DquadFFRF::EvalM(Matrix& m, double t)
{
	//UO() << "Mff=" << massmatrix << "\n";
	//UO() << "S~=" << Sbar_tilde << "\n";

	static Vector temp;
	static Vector temp2;
	//static Matrix mtemp;

	int off = Dim()*NS(); //offset where rigid body entries start!

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//mRR = Unit(2,2) * mass
	m(1+off,1+off) = GetMass();
	m(2+off,2+off) = GetMass();

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//mRtheta = A_theta*[I1+Sbar*qf]
	Vector2D I1;
	I1(1) = I1S(1);
	I1(2) = I1S(2);

	Vector2D pt(0.,0.);

	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS()*Dim(); j++)
		{
			pt(i) += SbarS(i,j)*XG(j);
		}
	}

	I1 += pt;
	Matrix3D ADphi = ReferenceFrame().GetRotMatrixDphi2D();
	I1 = ADphi*I1;
	m(1+off,3+off) = I1(1);
	m(2+off,3+off) = I1(2);
	m(3+off,1+off) = I1(1);
	m(3+off,2+off) = I1(2);


	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//m_Rf: = A*Sbar //very important term

	Matrix3D A = ReferenceFrame().GetRotMatrix2D();

	//if (!GetMBS()->IsJacobianComputation())
	{

		for (int j = 1; j <= NS()*Dim(); j++)
		{
			Vector2D v(SbarS(1,j),SbarS(2,j)); //stimmt mit rho*GetH() überein!!!
			v = A*v;
			m(j,1+off) = v(1);
			m(j,2+off) = v(2);
			m(1+off,j) = v(1);
			m(2+off,j) = v(2);
		}
	}


	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//m_theta_theta_rr: I_11+I_22
	m(3+off,3+off) = I1122S; //(GetIkl(1,1)+GetIkl(2,2));

	//m_theta_theta_rf: 2*(Ibar_11+Ibar_22)*q_f
	double I11, I22;
	I11 = 0;
	for (int i=1; i <= NS()*Dim(); i++)
	{
		I11 += Ibar11S(i)*XG(i);
	}
	I22 = 0;
	for (int i=1; i <= NS()*Dim(); i++)
	{
		I22 += Ibar22S(i)*XG(i);
	}
	m(3+off,3+off) += 2.*(I11+I22);

	//m_theta_theta_ff: q_f^T*m_ff*q_f
	double mthth = 0;
	for (int i=1; i <= NS()*Dim(); i++) xg(i) = XG(i);
	//Mult(massmatrix,xg,temp);

	for (int i=1; i <= NS(); i+=2) //only for isoparametric elements!!!!
	{
		for (int j=1; j <= NS(); j+=2)
		{
			mthth += xg(i)*massmatrix(i,j)*xg(j);
			mthth += xg(i+1)*massmatrix(i+1,j+1)*xg(j+1);
		}
	}
	m(3+off,3+off) += mthth;

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//m_ff
	for (int i=1; i <= NS()*Dim(); i++)
	{
		for (int j=1; j <= NS()*Dim(); j++)
		{
			m(i,j) = massmatrix(i,j);
		}
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//m_theta_f:
	temp.SetLen(Dim()*NS());
	temp = Ibar12S;
	temp -= Ibar21S; //-I_21+I_12

	for (int i=1; i <= Dim()*NS(); i++) // qf^T*(S12-S21)
	{
		for (int j=1; j <= Dim()*NS(); j++)
		{
			temp(i) += XG(j)*Sbar_tilde(j,i);
		}
	}

	//if (!GetMBS()->IsJacobianComputation())
	{
		for (int i=1; i <= Dim()*NS(); i++) 
		{
			m(i,off+3) = temp(i);
			m(off+3,i) = temp(i);
		}
	}


}


void Plate2DquadFFRF::EvalMff(Matrix& m, double t) 
{
	if (massmatrix.Getcols() == NS()*Dim())
	{
		m = massmatrix;
		return;
	}
	else
	{
		int dim = Dim(); 
		int ns = NS();

		m.SetSize(NS()*Dim(),NS()*Dim());
		m.SetAll(0);

		SV.SetLen(ns);
		Matrix3D jac;
		jac.SetSize(dim,dim);

		Matrix HL(ns*dim,dim);
		DS.SetSize(dim,ns);

		GetIntegrationRule(x1,w1,orderxyM); //optimal: 6x3x3, 6x1x1 geht auch!!!!
		GetIntegrationRule(x2,w2,orderxyM);

		for (int i1=1; i1<=x1.GetLen(); i1++)
		{
			for (int i2=1; i2<=x2.GetLen(); i2++)
			{
				Vector2D p(x1(i1),x2(i2));
				GetS0(SV,p);

				for (int i=0; i<ns; i++)
				{
					for (int j=1; j<=dim; j++)
					{
						HL(i*dim+j,j)=SV(i+1);
					}
				}

				GetJacobi(jac,p);

				double jacdet = jac.Det();
				m += fabs (jacdet) * Rho() * w1(i1)*w2(i2)*lz * (HL*HL.GetTp());
			}
		}

		massmatrix = m;
		//UO() << "m=" << m << "\n";

	}
};


//insert quadratic velocity vector
void Plate2DquadFFRF::EvalF2(Vector& f, double t)
{
	Body2D::EvalF2(f,t);
	//UO() << "load-f=" << f << "\n";

	//UO() << "PlateFFRF EvalF2\n";

	double nu = Nu();

	TMStartTimer(22);

	int evalfnonlin = 0;

	if (evalfnonlin)
	{
		//nonlinear:
		SetComputeCoordinates();
		static Vector fadd;

		int ns = NS();
		int dim = Dim();
		int sos = dim*ns;

		double la=Em() * nu / ((1.+nu)*(1.-2.*nu));
		double mu=Em() / 2. / (1.+nu);

		Matrix3D strain, piola1, F;
		strain.SetSize(2,2);
		piola1.SetSize(2,2);
		F.SetSize(2,2);

		temp.SetLen(SOS());
		fadd.SetLen(SOS());
		fadd.SetAll(0);

		Vector3D sigma;
		double k = Em()/(1.-nu*nu);
		Matrix3D C(k,k*nu,0,k*nu,k,0,0,0,k*(1.-nu)/2.); //Bathe p. 229, mit Gleitung!!!

		GetIntegrationRule(x1,w1,orderxy);
		GetIntegrationRule(x2,w2,orderxy);

		int kx1 = x2.Length();

		for (int i1=1; i1<=x1.GetLen(); i1++)
		{
			for (int i2=1; i2<=x2.GetLen(); i2++)
			{
				//TMStartTimer(20);
				int i,j,k;
				int ind = (i1-1)*kx1+(i2-1);

				// compute F 
				F.SetAll(0);
				int l;
				for (j = 1; j <= dim; j++) 
				{
					for (i = 1; i <= ns; i++)
					{
						l = (i-1)*dim+j;
						for (k = 1; k <= dim; k++)
						{
							F(j,k) += grad[ind](k,i)*xg(l);
						}
					}
					F(j,j) += 1;
				}
				// Green-Lagrange strain tensor
				strain = 0.5 * (F.GetTp() * F);
				strain(1,1) -= 0.5; strain(2,2) -= 0.5;

				Vector3D eps(strain(1,1), strain(2,2), 2.*strain(1,2));
				sigma = C*eps;

				piola1(1,1) = sigma(1);
				piola1(2,2) = sigma(2);
				piola1(1,2) = sigma(3);
				piola1(2,1) = sigma(3);

				piola1 = F * (piola1);

				for (int j=1; j <= Dim(); j++)
				{
					for (int i = 0; i < ns; i++)
					{
						temp(2*i+j) = grad[ind](1, i+1)*piola1(j,1)
							+ grad[ind](2, i+1)*piola1(j,2);
					}
				}

				fadd.MultAdd(fabs (jacdet[ind]) * w1(i1)*w2(i2)*lz,temp);
			} 
		}

		for (int i=1; i <= NS()*Dim(); i++)
		{
			f(i) -= fadd(i);
		}

	}
	else
	{

		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//add linear stiffness terms (negative on right hand side!!!!)
		for (int i=1; i <= NS()*Dim(); i++)
		{
			for (int j=1; j <= NS()*Dim(); j++)
			{
				f(i) -= K(i,j)*XG(j);
			}
		}
	}
	TMStopTimer(22);

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//quadratic velocity vector, Shabana p.229:

	static Vector temp; 
	static Vector Qvf;
	Qvf.SetLen(NS()*Dim());
	Qvf.SetAll(0);
	double Qvtheta;
	Vector2D QvR(0.,0.);

	Matrix3D A = ReferenceFrame().GetRotMatrix2D();
	Matrix3D Atheta = ReferenceFrame().GetRotMatrixDphi2D();

	double theta = ReferenceFrame().GetAngle2D();
	double thetap = ReferenceFrame().GetAngle2DP();


	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//Qv_R:
	Vector2D I1;
	I1(1) = I1S(1);
	I1(2) = I1S(2);

	Vector2D pt(0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS()*Dim(); j++)
		{
			pt(i) += SbarS(i,j)*XG(j);
		}
	}

	I1 += pt;
	QvR = Sqr(thetap)*(A*I1);

	pt=Vector2D(0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS()*Dim(); j++)
		{
			pt(i) += SbarS(i,j)*XGP(j);
		}
	}
	QvR += (-2.*thetap)*(Atheta*pt);


	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//Qv_theta:
	temp.SetLen(Dim()*NS());
	temp.SetAll(0);

	double mthth = 0;

	xg.SetLen(Dim()*NS());
	for (int i=1; i <= NS()*Dim(); i++) 
	{
		xg(i) = XG(i);
	}
	Mult(massmatrix,xg,temp);

	/*
	for (int i=1; i <= NS()*Dim(); i++)
	{
	for (int j=1; j <= NS()*Dim(); j++)
	{
	temp(i) += massmatrix(i,j)*xg(j); //m_ff*q_f
	}
	}*/
	temp += Ibar11S;
	temp += Ibar22S; // Ibar_0 = Ibar_11+Ibar_22

	for (int i=1; i <= NS()*Dim(); i++)
	{
		mthth += XGP(i)*temp(i);
	}
	Qvtheta = -2.*thetap * mthth; //does not influence up very much


	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//Qv_f: //temp is still m_ff*q_f+Ibar_0

	//does not influence number of jacobians!!!! maybe error in FFRF?
	/* if (GetMBS()->IsJacobianComputation())
	{
	temp.SetAll(0);
	temp += Ibar11S;
	temp += Ibar22S; // Ibar_0 = Ibar_11+Ibar_22
	}*/


	for (int i=1; i <= NS()*Dim(); i++)
	{
		Qvf(i) = Sqr(thetap)*temp(i);
	}

	/*
	for (int i=1; i <= NS()*Dim(); i++)
	{
	for (int j=1; j <= NS()*Dim(); j++)
	{
	Qvf(i) += 2.*thetap*Sbar_tilde(i,j)*XGP(j);
	}
	}*/

	//sparse version:
	for (int i=1; i <= NS()*Dim(); i++) xg(i) = XGP(i);

	Mult(Sbar_tildeSM,xg,temp);
	Qvf += (2.*thetap)*temp;

	//if (!GetMBS()->IsJacobianComputation())
	{
		//fill in velocity vector terms (positive on right hand side):
		for (int i=1; i <= NS()*Dim(); i++)
		{
			f(i) += Qvf(i);
		}

		//if (!GetMBS()->IsJacobianComputation())
		f(NS()*Dim()+1) += QvR(1);
		f(NS()*Dim()+2) += QvR(2);
		f(NS()*Dim()+3) += Qvtheta;
	}


	//mass damping ... ?
	if (this->GetMassDamping() != 0) UO() << "Error: Mass damping in FFRF not possible!!!\n";

}

//->only for volumeloads (gravity ...)
void Plate2DquadFFRF::GetIntDuDq(Matrix& dudq)
{
	//UO() << "Not yet implemented\n";
	dudq.FillWithZeros();

	int off = NS()*Dim();
	double rho = Rho();

	//same as mRR/rho
	dudq(1+off,1) = GetMass()/Rho();
	dudq(2+off,2) = GetMass()/Rho();

	//same asmRtheta/Rho() = A_theta*[I1+Sbar*qf]/rho
	Vector2D I1;
	I1(1) = I1S(1);
	I1(2) = I1S(2);

	Vector2D pt(0.,0.);

	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS()*Dim(); j++)
		{
			pt(i) += SbarS(i,j)*XG(j);
		}
	}

	I1 += pt;
	Matrix3D ADphi = ReferenceFrame().GetRotMatrixDphi2D();
	I1 = ADphi*I1;
	dudq(3+off,1) = I1(1)/Rho();
	dudq(3+off,2) = I1(2)/Rho();

	//same as m_Rf = A*Sbar:

	Matrix3D A = ReferenceFrame().GetRotMatrix2D();

	for (int j = 1; j <= NS()*Dim(); j++)
	{
		Vector2D v(SbarS(1,j),SbarS(2,j)); //stimmt mit rho*GetH() überein!!!
		v = A*v;
		dudq(j,1) = v(1)/Rho();
		dudq(j,2) = v(2)/Rho();
	}

}

void Plate2DquadFFRF::GetdPosdqT(const Vector2D& ploc, Matrix& d)
{
	//UO() << "Not yet implemented\n";

	Vector2D prel = GetPos2Drel(ploc);
	//double phi = ReferenceFrame().GetAngle2D();
	Matrix3D A = ReferenceFrame().GetRotMatrix2D();

	//dprel/dq
	d.SetSize(NS()*Dim()+FFRFDim(),Dim());
	d.FillWithZeros();

	for (int i = 1; i <= NS(); i++)
	{
		Vector2D tmp = A*Vector2D(GetS0(ploc,i),0.);
		d((i-1)*Dim()+1,1)=tmp.X();
		d((i-1)*Dim()+1,2)=tmp.Y();

		tmp = A*Vector2D(0.,GetS0(ploc,i));
		d((i-1)*Dim()+2,1)=tmp.X();
		d((i-1)*Dim()+2,2)=tmp.Y();
	}

	//d_pref/dq:
	d(NS()*Dim()+1,1) = 1;
	d(NS()*Dim()+2,2) = 1;

	//d_A/dq*prel:
	Vector2D rotdphi_prel = ReferenceFrame().GetRotMatrixDphi2D()*prel;
	//UO() << "GetdPos:rotdphi=" << rotdphi_prel << "\n";
	d(NS()*Dim()+3,1) = rotdphi_prel.X();
	d(NS()*Dim()+3,2) = rotdphi_prel.Y();

	//UO() << "dpdq=" << d << "\n";
}


