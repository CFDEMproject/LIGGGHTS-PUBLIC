//#**************************************************************
//#
//# filename:             FiniteElement2D.cpp
//#
//# project:              
//#
//# author:               Gerstmayr Johannes, Aigner Larissa, PG, YV
//#
//# generated:						November, 2010
//# description:          general 2D FiniteElement
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
 
#include "rigid2d.h"
#include "node.h"
#include "finiteElement2D.h"
#include "finiteelement3d.h"			// for ffrfsize
#include "Material.h"
#include "solversettings_auto.h"
#include "femathhelperfunctions.h"

// this helper class will simplify browsing through all integration points for stiffness
// with the direct access to the associated data units
class IntPointsStiffnessIterator : public IntegrationPointsIterator
{
	FiniteElement2D * fe;
public:
	IntPointsStiffnessIterator(FiniteElement2D * fe_) :
			IntegrationPointsIterator(fe_->integrationRuleStiffness),
			fe(fe_)
			{
			}
	FiniteElement2D::IntegrationPointStiffnessMatrixLocalData & Data()
	{
		return *fe->integrationPointStiffnessMatrixLocalData(GetIndex());
	}
};

void FiniteElement2D::SetFiniteElement2D(int bodyindi, const TArray<int>& nodelist,
																				 int material_num, double thickness, const Vector3D& coli)
{
	FiniteElementGeneric<Body2D>::SetFiniteElementGeneric(bodyindi,nodelist,material_num,coli);
	assert(GetMaterial().IsPlanarMaterial());
	this->thickness = thickness;
};

void FiniteElement2D::CopyFrom(const Element& e)
{
	FiniteElementGeneric<Body2D>::CopyFrom(e);
	const FiniteElement2D& ce = (const FiniteElement2D&)e;
	thickness = ce.thickness;
}

void FiniteElement2D::Initialize() 
{
	Body2D::Initialize();
	ComputeMass();
}

//compute element transformation matrix (jacobian)
void FiniteElement2D::GetJacobi(Matrix2D& jac, const Vector2D& ploc) const
{
	jac.SetSize(2,2);
	int ns = NS();
	int dim = Dim();
	jac.SetAll(0.);
	ConstMatrix<FE2DmaxDOF*FE2DmaxDOF> dsMatrix(FlexDOF(),FlexDOF());
	GetDSMatrix0(ploc,dsMatrix);
	// double aa[2][2] = {	{1,1}, {2,2} };		// *YV - commented out
	for (int j = 1; j <= dim; j++)
	{
		for (int k=1; k <= ns; k++)
		{ 
			Vector2D p0 = GetNodeRefPos2D(k);
			double sum = 0;
			for (int i = 1; i <= dim; i++)
			{ 
				jac(i,j) += dsMatrix(j,k)*p0(i);
			}
		}
	}
}

//compute DSMAtrix and apply element jacobian inverse; also compte determinant of jacobian
double FiniteElement2D::GetJacInvDS(const Vector2D& ploc, Matrix& jacinvDS) const
{
	//jacinvDS should have already Dim() x NS()
	jacinvDS.SetSize(Dim(), NS());

	Matrix2D jac0, jacinv;
	jacinv.SetSize(Dim(),Dim());

	GetJacobi(jac0,ploc);

	jac0.GetInverse(jacinv);
	jacinv.TpYs();

	GetDSMatrix0(ploc,jacinvDS);
	for (int i=1; i <= NS(); i++)
	{
		Vector2D v(jacinvDS(1,i),jacinvDS(2,i));
		v = jacinv*v;
		jacinvDS(1,i) = v(1);
		jacinvDS(2,i) = v(2);
	}
	//Mult(jacinv, DS, grad(i));

	return jac0.Det();
}

void FiniteElement2D::Gradu(const Vector& u, const Matrix& jacinvDS, Matrix2D& gradu) const
{
	gradu.SetAll(0);
	int dim = Dim();
	int l;
	for (int j = 1; j <= dim; j++) 
	{
		for (int i = 1; i <= NS(); i++)
		{
			l = (i-1)*dim+j;
			for (int k = 1; k <= dim; k++)
			{
				gradu(j,k) += jacinvDS(k,i)*u(l);
			}
		}
	}
}

void FiniteElement2D::Gradu(const Vector2D& ploc, const Vector& u, Matrix2D& gradu) const
{
	ConstMatrix<FE2DmaxDOF> jacinvDS(Dim(),NS());

	GetJacInvDS(ploc, jacinvDS);
	Gradu(u, jacinvDS, gradu);
}

//ploc -1 ... +1
void FiniteElement2D::GraduD(const Vector2D& ploc, Matrix2D& gradu) const
{
	ConstVector<FE2DmaxDOF> xgd(SOS(),1);		//*YV,JG
	
	for (int i=1; i<=SOS(); i++)
	{
		xgd(i) = XGD(i);
	}
	Gradu(ploc, xgd, gradu);
}

//-1..+1 based!!!
//FFRF: deformed position relative to floating frame; else: absolute deformed position
Vector2D FiniteElement2D::GetPos2DRel(const Vector2D& p_loc) const
{
	Vector2D p; //initialized with zeros

	for (int j = 1; j <= NS(); j++)
	{
		double s0 = GetS0(p_loc,j);
		Vector2D p0 = GetNodeRefPos2D(j);
		for (int i = 1; i <= Dim(); i++)
		{
			p(i) += s0*(XG((j-1)*Dim()+i) + p0(i));
		}
	}
	return p;
};

//-1..+1 based!!!
Vector2D FiniteElement2D::GetDisplacement2DRel(const Vector2D& p_loc) const
{
	Vector2D p; //initialized with zeros

	for (int j = 1; j <= NS(); j++)
	{
		double s0 = GetS0(p_loc,j);
		for (int i = 1; i <= Dim(); i++)
		{
			p(i) += s0*XG((j-1)*Dim()+i);
		}
	}
	return p;
};

//-1..+1 based!!!
Vector2D FiniteElement2D::GetVel2DRel(const Vector2D& p_loc) const
{
	Vector2D p; //initialized with zeros

	for (int j = 1; j <= NS(); j++)
	{
		double s0 = GetS0(p_loc,j);
		for (int i = 1; i <= Dim(); i++)
		{
			p(i) += s0*XGP((j-1)*Dim()+i);
		}
	}
	return p;
};

//-1..+1 based!!!
Vector2D FiniteElement2D::GetPos2DRelD(const Vector2D& p_loc, int use_magnification) const
{
	Vector2D p; //initialized with zeros
	double factor = GetMBS()->GetDOption(105);
	if (!use_magnification) factor = 1;

	for (int j = 1; j <= NS(); j++)
	{
		double s0 = GetS0(p_loc,j);
		Vector2D p0 = GetNodeRefPos2D(j);
		for (int i = 1; i <= Dim(); i++)
		{
			p(i) += s0*(XGD((j-1)*Dim()+i)*factor + p0(i));
		}
	}
	return p;
};

//-1..+1 based!!!
Vector2D FiniteElement2D::GetDisplacement2DRelD(const Vector2D& p_loc) const
{
	Vector2D p; //initialized with zeros

	for (int j = 1; j <= NS(); j++)
	{
		double s0 = GetS0(p_loc,j);
		for (int i = 1; i <= Dim(); i++)
		{
			p(i) += s0*XGD((j-1)*Dim()+i);
		}
	}
	return p;
};

//-1..+1 based!!!
Vector2D FiniteElement2D::GetVel2DRelD(const Vector2D& p_loc) const
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
//relative position for FFRF elements, in reference configuration
Vector2D FiniteElement2D::GetRefConfPos2D(const Vector2D& p_loc) const
{
	Vector2D p;

	for (int j = 1; j <= NS(); j++)
	{
		Vector2D p0 = GetNodeRefPos2D(j);
		for (int i = 1; i <= Dim(); i++)
		{
			p(i) += GetS0(p_loc,j)*p0(i);
		}
	}
	return p;
};

void FiniteElement2D::GetNodedPosdqT(int node, Matrix& dpdqi)
{
	dpdqi.SetSize(SOS(),Dim());
	dpdqi.FillWithZeros();

	dpdqi((node-1)*2+1,1) = 1; //dnx/dnx
	dpdqi((node-1)*2+2,2) = 1; //dny/dny
}

Vector2D FiniteElement2D::GetNodePos2D(int i) const
{
	return Vector2D(XG((i-1)*2+1), XG((i-1)*2+2)) + GetNodeRefPos2D(i);
}

Vector2D FiniteElement2D::GetNodePos2DD(int i) const 
{
	return Vector2D(XGD((i-1)*2+1), XGD((i-1)*2+2)) + GetNodeRefPos2D(i);
}

Vector2D FiniteElement2D::GetNodeVel2D(int i) const
{
	return Vector2D(XGP((i-1)*2+1), XGP((i-1)*2+2));
}

void FiniteElement2D::ComputeMass()
{
	mass = 0;
	Matrix2D jac;
	for (IntegrationPointsIterator ip(integrationRuleMass); !ip.IsEnd(); ++ip)
	{
		GetJacobi(jac, ip.Point2D());
		double jacdet = jac.Det();
		mass += fabs(jacdet) * Rho() * ip.Weight() * thickness;
	}
}


// fill only FlexibleDof entries of massmatrix, other entries remain unchanged
void FiniteElement2D::EvalMff(Matrix& m, double t)
{
	// Ioption 20 -> use elementwise precomputed mass matrix
	if (mbs->GetSolSet().store_FE_matrices && massmatrix.Getcols() == FlexDOF())
	{
		m.SetSubmatrix(massmatrix,1,1);
		return;
	}
	else
	{
		//computed only once and stored afterwards
		int dim = Dim(); 
		int ns = NS();
		int sos = SOS();
		int flexdof = FlexDOF();
		m.SetAll(0);

		ConstVector<FE2DmaxDOF> SF(ns);
		ConstMatrix<FE2DmaxDOF*FE2DmaxDOF> SST(flexdof,flexdof);
		Matrix2D jac;
		for (IntegrationPointsIterator ip(integrationRuleMass); !ip.IsEnd(); ++ip)
		{
			GetJacobi(jac, ip.Point2D());
			double jacdet = jac.Det();
			// compute SF = vector containing shape functions
			for (int i=1; i<=ns; i++)
			{
				SF(i)=GetS0(ip.Point2D(), i);
			}

			// compute S * S^T  where S is the shape function matrix 
			//     [ SF(1) 0 0 SF(2) 0 0 ...
			// S = [ 0 SF(1) 0 0 SF(2) 0 ..

			SST.SetAll(0);
			for (int i=0; i<ns; i++)
				for (int j=0; j<ns; j++)
					for (int di=1; di<=dim; di++)
							SST(i*dim+di,j*dim+di) += SF(i+1)*SF(j+1);

			m.AddSubmatrix(SST,1,1,1,1,flexdof,flexdof,
				fabs(jacdet) * GetMaterial().Density() * ip.Weight() * thickness);
		}

		// fill only FlexibleDof entries of massmatrix, other entries remain unchanged
		if (mbs->GetSolSet().store_FE_matrices) // if mass matrix is stored
		{
			massmatrix.SetSize(flexdof, flexdof);
			for (int i=1; i<=flexdof; i++)
				for (int j=1; j<=flexdof; j++)
					massmatrix(i,j) = m(i,j);
		}
	}
}

void FiniteElement2D::GetH(Matrix& H) 
{
	int ns = NS();
	int sos = SOS();
	int dim = Dim();
	if ( Hmatrix.Getrows() == ns*dim)
	{
		H = Hmatrix;
		return;
	}
	else
	{
		H.SetSize(ns*dim,dim);
		H.SetAll(0);
		ConstMatrix<(FE2DmaxDOF+FFRFsize)*2> dudq;
		Matrix2D jac;
		// flexible part Hf
		for (IntegrationPointsIterator ip(integrationRuleLoad); !ip.IsEnd(); ++ip)
		{
			GetJacobi(jac, ip.Point2D());
			double jacdet = jac.Det();
			double fact = fabs(jacdet) * ip.Weight();
			for (int i=0; i<ns; i++)
			{
				for (int j=1; j<=dim; j++)
				{
					H(i*dim+j,j) += fact * GetS0(ip.Point2D(), i+1) * thickness;
				}
			}
		}
		Hmatrix = H;
	}
}


void FiniteElement2D::GetdPosdqT(const Vector2D& ploc, Matrix& d)
{
	d.SetSize(SOS(),Dim());
	d.FillWithZeros();

	for (int i = 1; i <= NS(); i++)
	{
		double sv =	GetS0(ploc, i);
		d((i-1)*Dim()+1,1) = sv;
		d((i-1)*Dim()+2,2) = sv;
	}
}

//fill in sos x sos components, m might be larger
//compute stiffness matrix
void FiniteElement2D::StiffnessMatrix(Matrix& m) 
{
	//linear stiffness matrix: for linear problems or component mode synthesis:
	int dim = Dim(); 
	int ns = NS();

	if (
		mbs->GetSolSet().store_FE_matrices &&
		GetGeometricNonlinearityStatus() == GNS_Linear &&
		stiffnessmatrix.Getcols() == FlexDOF()
		)
	{
		m.SetSubmatrix(stiffnessmatrix,1,1);
		return;
	}

	m.SetAll(0);

	// B = deps/dq in the linear case, in the nonlinear case B = dE/dq
	ConstMatrix<6*FE2DmaxDOF> B((dim-1)*3,dim*ns); //filled with zeros
	// intermediate quantities
	ConstMatrix<FE2DmaxDOF*FE2DmaxDOF> helpmat(dim*ns,dim*ns), BDB(dim*ns,dim*ns); 

	// Elasticity Matrix

	ConstMatrix<36> C;
	GetMaterial().ComputeElasticityMatrix(C);

	// needed only in case of geometric nonlinear:
	Matrix2D F, sigma, strain;  // deformation gradient F, second piola kirchhoff tensor, nonlinear strain tensor 1/2(F^T F - I)
	// dB_nl/dq^T sigma (B_nl is the nonlinear part of B in case of geometric nonlinear problem, see Skript p. 38
	ConstMatrix<FE2DmaxDOF*FE2DmaxDOF> dB_nldq_sigma; // Contains dB_nldq_sigma(3*(alpha-1)+j, 3*(beta-1)+k) = dB_nl(i, 3*(alpha-1)+j)/dq(3*(beta-1)+k)*stress(i)

	ConstVector<FE2DmaxDOF> xgCache;
	// AP: xgCache is only needed for geometric nonlinear case
	//         do not set here, otherwise CMS-Element does not work
	//SetComputeCoordinates(xgCache);

	// Integration routine
	for (IntPointsStiffnessIterator ip(this); !ip.IsEnd(); ++ip)
	{
		//grad contains derivatives in compressed form:
		//grad(1,1..ns) = da/dx     where a = du/dqi or dv/dqi
		//grad(2,1..ns) = da/dy
		Matrix dShapei_dxj;
		GetGrad(ip, dShapei_dxj);

		// the computation is implemented only for the geometrically Linear case
		assert(GetGeometricNonlinearityStatus() == GNS_Linear);
		// Geometrically linearized B matrix - as in the old Plate2D
		for (int alpha = 1; alpha <= ns; alpha++)
		{
			B(1,2*alpha-1) = dShapei_dxj(1,alpha);
			B(2,2*alpha) = dShapei_dxj(2,alpha);
			B(3,2*alpha-1) = dShapei_dxj(2,alpha);
			B(3,2*alpha) = dShapei_dxj(1,alpha);
		}

		// compute B^T D B
		Mult(C, B, helpmat);
		MultTp(B, helpmat, BDB);
		
		// Stiffness matrix linear:    K = B^T D B,

		// integration weight factor
		double fact = -fabs(ip.Data().jacdet) * ip.Weight() * thickness;

		// Add part  B^T D B 
		// AddSubmatrix(mat, offset rows submat, offset cols submat, offset rows mat, offset cols mat, # rows, # cols, factor)
		m.AddSubmatrix(BDB, 1, 1, 1, 1, ns*dim, ns*dim, fact);
	}

	if (mbs->GetSolSet().store_FE_matrices && GetGeometricNonlinearityStatus() == GNS_Linear) // if stiffness matrix is stored
	{
		stiffnessmatrix.SetSize(FlexDOF(), FlexDOF());
		for (int i = 1; i <= FlexDOF(); i++)
			for (int j = 1; j <= FlexDOF(); j++)
				stiffnessmatrix(i,j) = m(i,j);
	}
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//COMPUTE RESIDUAL       COMPUTE RESIDUAL       COMPUTE RESIDUAL       COMPUTE RESIDUAL       COMPUTE RESIDUAL       COMPUTE RESIDUAL
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void FiniteElement2D::EvalF2(Vector& f, double t) 
{
	//add loads, constraint forces and other things in parent function:
	Body2D::EvalF2(f,t);
	TMStartTimer(22);

	int sos = SOS();
	int ns = NS();
	int dim = Dim();

	ConstVector<FE2DmaxDOF> fadd;
	fadd.SetLen(sos);
	fadd.SetAll(0);

	assert(GetGeometricNonlinearityStatus() == GNS_Linear);

	ConstVector<FE2DmaxDOF> xgCache;
	SetComputeCoordinates(xgCache);

	//double u;

	Matrix2D strain, stress, gradu, F;

	ConstVector<FE2DmaxDOF> temp;
	temp.SetLen(sos);
	temp.SetAll(0.);

	//compute deformation gradient, compute strain, compute stress and compute elastic forces:
	for (IntPointsStiffnessIterator ip(this); !ip.IsEnd(); ++ip)
	{
		//int_V (P : d grad u / d q) dV
		//P = D : E
		//E = 1/2*(grad u + grad u^T)

		Matrix agrad;
		GetGrad(ip, agrad);

		Gradu(xgCache, agrad, gradu);
		F = gradu;
		F(1,1) += 1.; F(2,2) += 1.;

		// Green-Lagrange strain tensor
		strain = 0.5 * (gradu.GetTp() + gradu);

		// potentially subtract initial strain
		Matrix2D initial_strain;
		if(HasInitialStrains())
		{
			GetInitialStrain(initial_strain, (const Vector2D&)(ip.Point()), false);
		}

		Vector inelastic_variables(0);
		// get old inelastic strain (from last nonlin step, or from last time step, depending on the specific solution method)
		if (IsInelasticMaterial())
		{
			Material::InelasticitySolutionMethod solmeth = GetMaterial().GetInelasticitySolutionMethod();
			switch (solmeth)
			{
			case Material::ISM_ConsistentTangentStiffness:
			case Material::ISM_ReturnMapping:
				DataLastStepToInelasticVariables(inelastic_variables, ip.GetIndex());  // read old inelastic strain from last computed time/load step
				break;
			case Material::ISM_FixedPoint:
			case Material::ISM_Default:
				DataToInelasticVariables(inelastic_variables, ip.GetIndex());  // read old inelastic strain from last nonlinear (discontinuous) step
				break;
			default:
				mbs->UO(UO_LVL_warn) << "WARNING: in FiniteElement2D::EvalF2(double t) not implemented for inelasticity-solution-method " << solmeth << "\n";
			}
		}

		GetMaterial().ComputeStressFromStrain2D(strain, stress, inelastic_variables, initial_strain);

		for (int j=1; j <= Dim(); j++)
		{
			for (int i = 1; i <= ns; i++)
			{
				int k = 2*(i-1)+j;
				temp(k) = agrad(1, i)*stress(j,1)
					+ agrad(2, i)*stress(j,2);
				fadd(k) += fabs(ip.Data().jacdet) * ip.Weight() * temp(k) * thickness;
			}
		}
	}

	f -= fadd;
	TMStopTimer(22);

	// +++++++ damping: +++++++
	if (GetMaterial().DampingM() != 0)
	{
		// set velocities into a temporary vector, modified by YV:
		ConstVector<FE2DmaxDOF> velocities;
		velocities.SetLen(sos);
		for (int i = 1; i <= sos; i++)
			velocities(i) = XGP(i);

		ConstVector<FE2DmaxDOF> temp;
		temp.SetLen(sos);
		temp.SetAll(0.);

		if (mbs->GetSolSet().store_FE_matrices) // mass matrix is stored
		{
			if (massmatrix.Getcols() != sos) {UO() << "ERROR: finite element: mass matrix not built!\n"; return;}
			Mult(massmatrix,velocities,temp);
		}
		else // mass matrix is not stored
		{
			ConstMatrix<FE2DmaxDOF*FE2DmaxDOF> m(FlexDOF(),FlexDOF());
			EvalMff(m,t);
			Mult(m,velocities,temp);
		}
		temp *= GetMaterial().DampingM();
		f -= temp;
	}
};

//EK 2012-07-10 integrate the function given by the values of the shape functions in f1 and f2
double FiniteElement2D::L2InnerProduct(Vector& f1, Vector &f2)
{
	int sos = SOS();
	int ns = NS();
	int dim = Dim();
	double res = 0;

	assert(GetGeometricNonlinearityStatus() == GNS_Linear);

	double tmp;
	ConstVector<FE2DmaxDOF> shape;
	shape.SetLen(ns);

	for (IntPointsStiffnessIterator ip(this); !ip.IsEnd(); ++ip)
	{	
		for (int i = 1; i <= ns; ++i)
			shape(i) = GetS0(ip.Point2D(), i);

		for (int j=1; j <= Dim(); j++)
		{
			for (int i = 1; i <= ns; i++)
			{
				for (int k = 1; k <= ns; k++)
				{
					//int k = 2*(i-1)+j;
					//temp(k) = agrad(1, i)*stress(j,1)
					//	+ agrad(2, i)*stress(j,2);
					tmp = shape(i)*f1((i-1)*Dim() + j)*shape(k)*f2((k-1)*Dim() + j);
					res += fabs(ip.Data().jacdet) * ip.Weight() * tmp * thickness;
				}
			}
		}
	}
	return res;
};


void FiniteElement2D::GetIntDuDqFCentrifugal(Matrix& H, const Vector3D& omega, const Vector3D& r0)
{	
	int dim = Dim();
	int ns = NS();
	int sos = SOS();

	H.SetSize(sos,dim);
	H.SetAll(0);
	Matrix2D jac;
	for (IntegrationPointsIterator ip(integrationRuleLoad); !ip.IsEnd(); ++ip)
	{
		GetJacobi(jac, ip.Point2D());
		double jacdet = jac.Det();
		double fact = fabs (jacdet) * ip.Weight();
		Vector2D r_rel = r0.MakeV2D() - GetPos2D(ip.Point2D());
//		Vector2D force = omega.Cross(omega.Cross(r_rel));
		Vector2D force = r_rel*omega.Z();

		for (int i = 0; i < ns; i++)
		{
			double combinedfactor = fact*GetS0(ip.Point2D(), i+1)*Rho();
			H(i*dim+1,1) += combinedfactor * force.X();
			H(i*dim+2,1) += combinedfactor * force.Y();
			//H(i*dim+3,1) += combinedfactor * force.Z();
		}
	}
}

double FiniteElement2D::PostNewtonStep(double t)
{
	if (!IsInelasticMaterial()) return 0;
	Material::InelasticitySolutionMethod solmeth = GetMaterial().GetInelasticitySolutionMethod();

	//----------------------
	// Variable declaration and Parameters
	double err = 0.; // error for nonlinear iteration
	double fe_volume = 0;

	int sos = FlexDOF();
	ConstVector<FE2DmaxDOF> xgCache;
	SetComputeCoordinates(xgCache);

	int ns = NS();
	int dim = Dim();

	Matrix2D F,jac;
	ConstVector<MAXInelasticVariablesCount> inelastic_variables(GetMaterial().GetInelasticVariablesCount());
	ConstVector<MAXInelasticVariablesCount> inelastic_variables_last_step(GetMaterial().GetInelasticVariablesCount());
	ConstVector<MAXInelasticVariablesCount> inelastic_variables_update(GetMaterial().GetInelasticVariablesCount());

	//----------------------
	// loop over integration points
	for (IntPointsStiffnessIterator ip(this); !ip.IsEnd(); ++ip)
	{
		// Jacobi-matrix and Jacobi-determinant for integration point p
		GetJacobi(jac,ip.Point2D());
		double jacdet = jac.Det();

		// compute F = grad u +I
		Matrix agrad;
		GetGrad(ip, agrad);
		Gradu(xgCache, agrad, F);
		F(1,1) += 1; F(2,2) += 1;
		
		// compute strain from F
		Matrix2D strain;
		if ( GetGeometricNonlinearityStatus() == GNS_Linear || IsFFRF())
		{
			strain = 0.5*(F.GetTp() + F) - Matrix2D(1.);
		}
		else
		{
			strain = 0.5*(F.GetTp()*F - Matrix2D(1.));
		}

		// potentially subtract initial strain
		Matrix2D initial_strain;
		if(HasInitialStrains())
		{
			GetInitialStrain(initial_strain, (const Vector2D&)(ip.Point()), false);
		}

		// get old inelastic strain (from last time step)
		DataLastStepToInelasticVariables(inelastic_variables_last_step, ip.GetIndex());  // read old inelastic strain from last computed time/load step
		
		// and last computed update of inelastic strain (either from last nonlin step (return mapping) or also from last time step (consistent tangent stiffness matrix))
		switch (solmeth)
		{
		case Material::ISM_ConsistentTangentStiffness:
		case Material::ISM_FixedPoint:
			DataLastStepToInelasticVariables(inelastic_variables, ip.GetIndex());  // read old inelastic strain from last computed time/load step
			break;
		case Material::ISM_ReturnMapping:
		case Material::ISM_Default:
			DataToInelasticVariables(inelastic_variables, ip.GetIndex());  // read old inelastic strain from last nonlinear (discontinuous) step
			break;
		default:
			mbs->UO(UO_LVL_warn) << "WARNING: in FiniteElement2D::PostNewtonStep(double t) not implemented for inelasticity-solution-method " << solmeth << "\n";
		}

		// compute actual inelastic strain as function from old inelastic strain and elastic strain (minus initial strain)
		err += fabs(jacdet) * ip.Weight() * GetMaterial().ComputeInelasticVariablesUpdateFromStrain2D(inelastic_variables_update, inelastic_variables, strain - initial_strain);
		fe_volume += fabs(jacdet) * ip.Weight();

		inelastic_variables = inelastic_variables_last_step + inelastic_variables_update;
		InelasticVariablesToData(inelastic_variables, ip.GetIndex());
	}

	return err/fe_volume;
}

void FiniteElement2D::PostprocessingStep()
{
}

void FiniteElement2D::DrawElementPreProc() 
{
	if (!GetMBS()->GetIOption(118)) return;

	FieldVariableDescriptor * fvd = GetMBS()->GetActualPostProcessingFieldVariable();
	if(fvd != NULL)
	{
		for (int i=1; i <= NNodes(); i++)
		{
			GetNode(i).GetDrawTemp() += GetFieldVariableValue(*fvd, (1./sqrt(2.))*GetNodeLocPos2D(i), true);
			GetNode(i).GetDrawTempCnt()++;
		}
	}
};

void FiniteElement2D::DrawElement()
{
	mbs->SetColor(col);

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

	int res = GetMBS()->GetDrawResolution() - 1;
	if (colormode == 0 && NS() <= 4)
		res = 0;		// linear elements, not colored
	if (res < 0) res = 0;

	double tilex = pow(2.,res);

	double tiley = tilex;
	int tileyn = (int)tiley;

	double v = 0;

	TArray<Vector3D> points((int)(tilex+1)*(int)(tiley+1));
	TArray<double> vals((int)(tilex+1)*(int)(tiley+1));
	points.SetLen(0);
	vals.SetLen(0);

	if (GetMBS()->GetIOption(118) && NNodes() == 4)
	{
		// plotting linear element interpolated - can be done quickly
		// add rectangular grid of points and values
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
		// full version
		Vector2D p0, vx, vy;

		double shrink = GetMBS()->GetDOption(106);//0.995;
		p0 = Vector2D(-1*shrink,-1*shrink);
		vx = Vector2D(2./tilex*shrink,0);
		vy = Vector2D(0,2./tileyn*shrink);

		for (double iy = 0; iy <= tileyn+1e-8; iy++)
		{
			for (double ix = 0; ix <= tilex+1e-8; ix++)
			{
				Vector2D ploc = p0 + ix * vx + iy * vy;
				Vector3D pg = ToP3D(GetPos2DD(ploc, 1));
				points.Add(pg);
				if (colormode)
					v = GetFieldVariableValue(*GetMBS()->GetActualPostProcessingFieldVariable(), ploc, true);
				vals.Add(v);
			}
		}
		mbs->DrawColorQuads(points,vals,(int)tilex+1,(int)tileyn+1,colormode,linemode);
	}
}

void FiniteElement2D::GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables)
{
	// add all the field variables of the parent class
	Body2D::GetAvailableFieldVariables(variables);
	//FVT_position,FVT_displacement, FVT_velocity: implemented in Element
	// possibility to remove entries does not exist yet, just do not use the Function of the parent class
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_stress,
		FieldVariableDescriptor::FVCI_y, true);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_stress_mises);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_total_strain,
		FieldVariableDescriptor::FVCI_y, true);
	if(GetMaterial().IsInelasticMaterial())
	{
		FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_inelastic_strain,
			FieldVariableDescriptor::FVCI_y, true);
		if(GetMaterial().GetInelasticityType() == 1) 
		{
			FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_hardening_parameter);
			FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_yield_function);
		}
	}
}

void FiniteElement2D::UpdateFieldVariableComputationPoint(const Vector2D & local_position_original, Vector2D & local_position_actual, Matrix2D & initial_strain, Vector & inelastic_variables)
{
	assert(0); // !!!!!!!!!!!!!!!
	/*
	IntPointsStiffnessIterator ip(this);
	//plastic_strain.SetAll(0);			// no plastic strains yet
	if(ip.Data().initialStrainComponents == NULL && !IsInelasticMaterial())
	{
		// no initial strains - the point remains the same
		local_position_actual = local_position_original;
		initial_strain.SetAll(0.);
	}
	else
	{
		// initial strain is not known at the actual point;
		// therefore we simply find the nearest integration point and return the initial strain in it;
		// same holds for the vector of inelastic variables;
		double dist = 10;
		const Vector * isc;
		bool has_initial_strain = false;
		while(!ip.IsEnd())
		{
			double current_dist = (local_position_original - ip.Point2D()).Norm();
			if(current_dist < dist)
			{
				dist = current_dist;
				if(ip.Data().initialStrainComponents != 0)
				{
					isc = ip.Data().initialStrainComponents;
					has_initial_strain = true;
				}
				if(IsInelasticMaterial())
					DataToInelasticVariables(inelastic_variables, ip.GetIndex());
				local_position_actual = ip.Point2D();
			}
			++ip;
		}
		if(has_initial_strain)
		{
			initial_strain = Matrix2D(	(*isc)(1), (*isc)(3),
																 	(*isc)(3), (*isc)(2)   );
		}
	}
	*/
}

double FiniteElement2D::GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector2D & local_position, bool flagD)
{
	switch(fvd.VariableType())
	{
	case FieldVariableDescriptor::FVT_displacement:
		{
			if(flagD)
				return fvd.GetComponent(GetDisplacement2DD(local_position));
			else
				return fvd.GetComponent(GetDisplacement2D(local_position));
		}
	case FieldVariableDescriptor::FVT_position:
		{
			if(flagD)
				return fvd.GetComponent(GetPos2DD(local_position));
			else
				return fvd.GetComponent(GetPos2D(local_position));
		}
	case FieldVariableDescriptor::FVT_velocity:
		{
			if(flagD)
				return fvd.GetComponent(GetVel2DD(local_position));
			else
				return fvd.GetComponent(GetVel2D(local_position));
		}
	case FieldVariableDescriptor::FVT_stress:
	case FieldVariableDescriptor::FVT_stress_mises:
		{ //EK 2012-03-02 updated according to FiniteElement3D (needs to be tested)
						
			//Matrix2D initial_strain;
			Vector inelastic_variables;
			//EK -> talked with PG -> old version
			//UpdateFieldVariableComputationPoint(local_position, local_position_updated, initial_strain, inelastic_variables);
			//according to YV 2d has to be adapted as well

			if(IsInelasticMaterial())
			{
				IntPointsStiffnessIterator ip(this);
				Vector3D local_position3D(local_position.X(), local_position.Y(),0);
				ip.GoClosestTo(local_position3D);
				if(flagD)
				{
					DataToInelasticVariablesD(inelastic_variables, ip.GetIndex());
				}
				else
				{
					DataToInelasticVariables(inelastic_variables, ip.GetIndex());
				}
			}

			
			Matrix2D gradu;
			if(flagD) 
				GraduD(local_position, gradu);
			else
			{
				ConstVector<FE2DmaxDOF> xgCache;
				SetComputeCoordinates(xgCache);
				Gradu(local_position, xgCache, gradu);
			}
			
			Matrix2D strain;
			GetMaterial().ComputeStrain2D(strain, gradu, GetGeometricNonlinearityStatus() == GNS_NonlinearLargeStrain);
			
			//Matrix3D stress;
			Matrix2D stress;
			if(HasInitialStrains())
			{
				Matrix2D initial_strain;
				GetInitialStrain(initial_strain, local_position, flagD);
				GetMaterial().ComputeStressFromStrain2D(strain, stress, inelastic_variables, initial_strain);
			}
			else
			{
				GetMaterial().ComputeStressFromStrain2D(strain, stress, inelastic_variables);
			}

			if(fvd.VariableType() == FieldVariableDescriptor::FVT_stress_mises)
			{
				// mises is not defined for matrix 2D
				double tr = (1./3.)*stress.Trace();
				stress(1,1) -= tr;
				stress(2,2) -= tr;
				return sqrt(3./2.*(stress.InnerProduct(stress)));
			}
			else
				return fvd.GetComponent(stress);
		}
	case FieldVariableDescriptor::FVT_total_strain:
		{
			Matrix2D gradu;
			if(flagD) 
				GraduD(local_position, gradu);
			else
			{
				ConstVector<FE2DmaxDOF> xgCache;
				SetComputeCoordinates(xgCache);
				Gradu(local_position, xgCache, gradu);
			}
			Matrix2D strain;
			GetMaterial().ComputeStrain2D(strain, gradu, GetGeometricNonlinearityStatus() == GNS_NonlinearLargeStrain);
			return fvd.GetComponent(strain);
		}
	case FieldVariableDescriptor::FVT_initial_strain:
		{
			//EK initial strain not implemented so far
			mbs->UO() << "ERROR : Initial Strain not implemented so far\n";
			assert(0);
			return 0;
		}
	case FieldVariableDescriptor::FVT_inelastic_strain:
	case FieldVariableDescriptor::FVT_hardening_parameter:
	case FieldVariableDescriptor::FVT_yield_function:
		{
			Vector inelastic_variables;
			IntPointsStiffnessIterator ip(this);
			if(IsInelasticMaterial())
			{
				Vector3D local_position3D(local_position.X(), local_position.Y(),0);
				ip.GoClosestTo(local_position3D);
				if(flagD)
				{
					DataToInelasticVariablesD(inelastic_variables, ip.GetIndex());
				}
				else
				{
					DataToInelasticVariables(inelastic_variables, ip.GetIndex());
				}
			}

			switch(fvd.VariableType())
			{
			case FieldVariableDescriptor::FVT_inelastic_strain:
				{
					Matrix3D inelastic_strain = 0;
					if (IsInelasticMaterial())
					{
						if (GetMaterial().IsPlanarMaterial())
						{
							inelastic_strain(1,1) = inelastic_variables(1);
							inelastic_strain(2,2) = inelastic_variables(2);
							inelastic_strain(1,2) = inelastic_strain(1,2) = .5*inelastic_variables(3);
						}
						else
						{
							StrainVectorToMatrix3D(inelastic_strain, inelastic_variables);
						}
					}
					return fvd.GetComponent(inelastic_strain);
				}
			case FieldVariableDescriptor::FVT_hardening_parameter:
				{
					return inelastic_variables(4);
				}
			case FieldVariableDescriptor::FVT_yield_function:
				{
					//$ EK 2012-09-04 yield function 
					if (GetMaterial().GetInelasticityType() == Material::IT_ElastoPlastic)
					{
						Vector2D new_local_position = ip.Point2D();
						Matrix2D gradu;
						if(flagD) 
							GraduD(new_local_position, gradu);
						else
						{
							ConstVector<FE2DmaxDOF> xgCache;
							SetComputeCoordinates(xgCache);
							Gradu(new_local_position, xgCache, gradu);
						}
						Matrix2D strain;
						GetMaterial().ComputeStrain2D(strain, gradu, GetGeometricNonlinearityStatus() == GNS_NonlinearLargeStrain);
						if (HasInitialStrains())
						{
							Matrix2D initial_strain;
							GetInitialStrain(initial_strain, new_local_position, flagD);
							strain -= initial_strain;
						}
						return GetMaterial().ComputeYieldFunction2D(strain, inelastic_variables);
					}
					break;
				}
			}
		}
	}
	return FIELD_VARIABLE_NO_VALUE;
}