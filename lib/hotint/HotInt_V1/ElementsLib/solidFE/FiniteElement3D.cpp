//#**************************************************************
//#
//# filename:             FiniteElement3D.cpp
//#
//# project:              
//#
//# author:               Gerstmayr Johannes, Aigner Larissa
//#
//# generated:						May 11, 2009
//# description:          general 3D FiniteElement
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


// this helper class will simplify browsing through all integration points for stiffness
// with the direct access to the associated data units
class IntPointsStiffnessIterator : public IntegrationPointsIterator
{
	FiniteElement3D * fe;
public:
	IntPointsStiffnessIterator(FiniteElement3D * fe_) :
			IntegrationPointsIterator(fe_->integrationRuleStiffness),
			fe(fe_)
			{
			}
	FiniteElement3D::IntegrationPointStiffnessMatrixLocalData & Data()
	{
		return *fe->integrationPointStiffnessMatrixLocalData(GetIndex());
	}
};

void FiniteElement3D::SetFiniteElement3D(int bodyindi, const TArray<int>& nodelist, int material_num, const Vector3D& coli, int CMSelementi)
{
	FiniteElementGeneric<Body3D>::SetFiniteElementGeneric(bodyindi,nodelist,material_num,coli,CMSelementi);
	assert(!GetMaterial().IsPlanarMaterial());
	ResetOuterFaceFlags();
};

void FiniteElement3D::CopyFrom(const Element& e)
{
	FiniteElementGeneric<Body3D>::CopyFrom(e);
	const FiniteElement3D& ce = (const FiniteElement3D&)e;
	outer_face = ce.outer_face;
}

void FiniteElement3D::Initialize() 
{
	Body3D::Initialize();
	ComputeMass();
}

//compute element transformation matrix (jacobian)
void FiniteElement3D::GetJacobi(Matrix3D& jac, const Vector3D& ploc) const
{
	jac.SetSize(3,3);
	int ns = NS();
	int dim = Dim();
	jac.SetAll(0.);
	ConstMatrix<FEmaxDOF*FEmaxDOF> dsMatrix(FlexDOF(),FlexDOF());
	GetDSMatrix0(ploc,dsMatrix);
	// double aa[2][2] = {	{1,1}, {2,2} };		// *YV - commented out
	for (int j = 1; j <= dim; j++)
	{
		for (int k=1; k <= ns; k++)
		{ 
			Vector3D p0 = GetNodeRefPos(k);
			double sum = 0;
			for (int i = 1; i <= dim; i++)
			{ 
				jac(i,j) += dsMatrix(j,k)*p0(i);
			}
		}
	}
}

//compute DSMAtrix and apply element jacobian inverse; also compte determinant of jacobian
double FiniteElement3D::GetJacInvDS(const Vector3D& ploc, Matrix& jacinvDS) const
{
	//jacinvDS should have already Dim() x NS()
	jacinvDS.SetSize(Dim(), NS());

	Matrix3D jac0, jacinv;
	jacinv.SetSize(3,3);

	GetJacobi(jac0,ploc);

	jac0.GetInverse(jacinv);
	jacinv.TpYs();

	GetDSMatrix0(ploc,jacinvDS);
	for (int i=1; i <= NS(); i++)
	{
		Vector3D v(jacinvDS(1,i),jacinvDS(2,i),jacinvDS(3,i));
		v = jacinv*v;
		jacinvDS(1,i) = v(1);
		jacinvDS(2,i) = v(2);
		jacinvDS(3,i) = v(3);
	}
	//Mult(jacinv, DS, grad(i));

	return jac0.Det();
}

void FiniteElement3D::Gradu(const Vector& u, const Matrix& jacinvDS, Matrix3D& gradu) const
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

void FiniteElement3D::Gradu(const Vector3D& ploc, const Vector& u, Matrix3D& gradu) const
{
	ConstMatrix<FEmaxDOF> jacinvDS(Dim(),NS());

	GetJacInvDS(ploc, jacinvDS);
	Gradu(u, jacinvDS, gradu);
}

//ploc -1 ... +1
void FiniteElement3D::GraduD(const Vector3D& ploc, Matrix3D& gradu) const
{
	ConstVector<FEmaxDOF> xgd(SOS(),1);		//*YV,JG
	
	for (int i=1; i<=SOS(); i++)
	{
		xgd(i) = XGD(i);
	}
	Gradu(ploc, xgd, gradu);
}

//set all outer faces active, which have no opposite face from an other element
void FiniteElement3D::SetOuterFaces()
{
	TArray<int> neighbours(20);
	GetNeighbourElements(neighbours);

	ResetOuterFaceFlags();

	for (int i=1; i<=NFaces(); i++)
	{
		FEFace face = GetFace(i);
		face.Invert();

		int hasface = 0;
		for (int j=1; j<=neighbours.Length(); j++)
		{
			int no_surfaceupdate = mbs->GetIOption(119);

			if (GetElement(neighbours(j)).IsFiniteElement() && GetElement(neighbours(j)).Dim()==3 && 
				(GetElement(neighbours(j)).DrawElementFlag()!=0 || no_surfaceupdate!=0) ) //Hex/Tet, drawbody flag is set if it matters
			{
				FiniteElement3D* fe = (FiniteElement3D*)GetElementPtr(neighbours(j));
				//const FiniteElement3D& h=(const FiniteElement3D&)(GetElement(neighbours(j)));
				if (fe->HasFace(face)) 
				{
					hasface = 1;
					break;
				}
			}
		}
		if (!hasface) SetOuterFaceFlag(i);
		//if (GetOuterFaceFlag(i) ) UO() << "face" << i << "=" << GetOuterFaceFlag(i) << "\n";
	}
}

//set all outer faces active, which have no opposite face from an other element (this version also checks visibility due to cutting plane)
void FiniteElement3D::SetOuterFacesCuttingPlane()
{
	TArray<int> neighbours(20);
	GetNeighbourElements(neighbours);

// check if any neighbor has drawflag set but is not drawn due to cutting plane
	int needrecalc = 0;
	for (int i=1; i<=neighbours.Length(); i++)
	{
		if(GetElement(neighbours(i)).DrawElementFlag() != 0)
		{
			if( !mbs->CuttingPlanesAllow(GetElement(neighbours(i)).GetRefPosD()))
			{
				// element has a neighbor that should be drawn due to flag, but is not drawn due to cutting plane
				needrecalc = 1;
				break;
			}
		}
	}
	if (!needrecalc)
		return;

	ResetOuterFaceFlags();

	for (int i=1; i<=NFaces(); i++)
	{
		FEFace face = GetFace(i);
		face.Invert();

		int hasface = 0;
		for (int j=1; j<=neighbours.Length(); j++)
		{
			if (GetElement(neighbours(j)).IsFiniteElement() && GetElement(neighbours(j)).Dim()==3 && 
				GetElement(neighbours(j)).DrawElementFlag()!=0 && mbs->CuttingPlanesAllow(GetElement(neighbours(j)).GetRefPosD())) //Hex/Tet, drawbody flag is set, neighbour is drawn
			{
				FiniteElement3D* fe = (FiniteElement3D*)GetElementPtr(neighbours(j));
				//const FiniteElement3D& h=(const FiniteElement3D&)(GetElement(neighbours(j)));
				if (fe->HasFace(face)) 
				{
					hasface = 1;
					break;
				}
			}
		}
		if (!hasface) SetOuterFaceFlag(i);
		//if (GetOuterFaceFlag(i) ) UO() << "face" << i << "=" << GetOuterFaceFlag(i) << "\n";
	}
}

//-1..+1 based!!!
//FFRF: deformed position relative to floating frame; else: absolute deformed position
Vector3D FiniteElement3D::GetPosRel(const Vector3D& p_loc) const
{
	Vector3D p; //initialized with zeros

	for (int j = 1; j <= NS(); j++)
	{
		double s0 = GetS0(p_loc,j);
		Vector3D p0 = GetNodeRefPos(j);
		for (int i = 1; i <= Dim(); i++)
		{
			p(i) += s0*(XG((j-1)*Dim()+i) + p0(i));
		}
	}
	return p;
};

//-1..+1 based!!!
Vector3D FiniteElement3D::GetDisplacementRel(const Vector3D& p_loc) const
{
	Vector3D p; //initialized with zeros

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
Vector3D FiniteElement3D::GetVelRel(const Vector3D& p_loc) const
{
	Vector3D p; //initialized with zeros

	for (int j = 1; j <= NS(); j++)
	{
		double s0 = GetS0(p_loc,j);
		for (int i = 1; i <= Dim(); i++)
		{
			double xgp = XGP((j-1)*Dim()+i);
			p(i) += s0*XGP((j-1)*Dim()+i);
		}
	}
	return p;
};

//-1..+1 based!!!
Vector3D FiniteElement3D::GetPosRelD(const Vector3D& p_loc, int use_magnification) const
{
	Vector3D p; //initialized with zeros
	double factor = GetMBS()->GetDOption(105);
	if (!use_magnification) factor = 1;

	for (int j = 1; j <= NS(); j++)
	{
		double s0 = GetS0(p_loc,j);
		Vector3D p0 = GetNodeRefPos(j);
		for (int i = 1; i <= Dim(); i++)
		{
			p(i) += s0*(XGD((j-1)*Dim()+i)*factor + p0(i));
		}
	}
	return p;
};

//-1..+1 based!!!
Vector3D FiniteElement3D::GetDisplacementRelD(const Vector3D& p_loc) const
{
	Vector3D p; //initialized with zeros

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
Vector3D FiniteElement3D::GetVelRelD(const Vector3D& p_loc) const
{
	Vector3D p;

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
Vector3D FiniteElement3D::GetRefConfPos(const Vector3D& p_loc) const
{
	Vector3D p;

	for (int j = 1; j <= NS(); j++)
	{
		Vector3D p0 = GetNodeRefPos(j);
		for (int i = 1; i <= Dim(); i++)
		{
			p(i) += GetS0(p_loc,j)*p0(i);
		}
	}
	return p;
}

Vector3D FiniteElement3D::GetNodeRefPos(int i) const //return nodal position of node i in reference configuration
{
	return GetNode(i).RefConfPos();
}

Vector3D FiniteElement3D::GetAngularVel(const Vector3D &ploc) const
{
	Vector3D omega;
	Vector3D ploc2 = ploc + Vector3D(0.5, 0.5, 0.5);
	Vector3D pglob = GetPos(ploc);
	Vector3D pglob2 = GetPos(ploc2);
	Vector3D r(0.), axis(0.);
	Vector3D v1, v2, v;
	double normr, invnormr;

	for (int i=1; i<=3; i++)
	{
		r = pglob2 - pglob;
		normr = r.Norm();
		invnormr = 1./normr;
		// axis is ith unit vector
		axis.Set(0.,0.,0.);
		axis(i) = 1.;
		int j=0;
		// find ploc2 such that pglob2+pglob is not parallel to axis
		while ( fabs(r*axis) > (0.55)*normr && j < 3 ) 
		{
			//global_uo << "not using " << axis << " and " << r << "\n";
			j++;
			Vector3D add(0.);
			add(j) = 1.;
			ploc2 = ploc + add;
			pglob = GetPos(ploc);
			pglob2 = GetPos(ploc2);
			r = pglob2 - pglob;
			normr = r.Norm();
			invnormr = 1./normr;
		}
		if (j==3)
		{
			UO() << "Error in FiniteElement3D::GetAngularVel: could not find direction not parallel to axis " << axis << "\n";
		}
		//double rax = r*axis;
		//r -= rax * axis;
		r(i) = 0;
		v1 = GetVel(ploc);
		v2 = GetVel(ploc2);
		v = v2-v1;
		//double vax = v*axis;
		//v -= vax*axis;
		v(i) = 0;
		double vr = -1./r.Norm2()*(v*r);
		v += vr * r;
		omega(i) = invnormr*v.Norm();
		//global_uo << "Pglob = " << pglob << ", axis = " << axis << " r = " << r << ", vel = " << v << "\n";
	}
	return omega;
}

void FiniteElement3D::GetNodedPosdqT(int node, Matrix& dpdqi)
{
	dpdqi.SetSize(SOS(),Dim());
	dpdqi.FillWithZeros();

	dpdqi((node-1)*3+1,1) = 1; //dnx/dnx
	dpdqi((node-1)*3+2,2) = 1; //dny/dny
	dpdqi((node-1)*3+3,3) = 1; //dny/dny
}

Vector3D FiniteElement3D::GetNodePos(int i) const
{
	return Vector3D(XG((i-1)*3+1), XG((i-1)*3+2), XG((i-1)*3+3)) + GetNodeRefPos(i);
}

Vector3D FiniteElement3D::GetNodePosD(int i) const 
{
	return Vector3D(XGD((i-1)*3+1), XGD((i-1)*3+2), XGD((i-1)*3+3)) + GetNodeRefPos(i);
}

Vector3D FiniteElement3D::GetNodeVel(int i) const
{
	return Vector3D(XGP((i-1)*3+1), XGP((i-1)*3+2), XGP((i-1)*3+3));
}

void FiniteElement3D::ComputeMass()
{
	mass = 0;
	Matrix3D jac;
	for (IntegrationPointsIterator ip(integrationRuleMass); !ip.IsEnd(); ++ip)
	{
		GetJacobi(jac, ip.Point());
		double jacdet = jac.Det();
		mass += fabs(jacdet) * Rho() * ip.Weight();
	}
}

// fill only FlexibleDof entries of massmatrix, other entries remain unchanged
void FiniteElement3D::EvalMff(Matrix& m, double t)
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

		ConstVector<FEmaxDOF> SF(ns);
		ConstMatrix<FEmaxDOF*FEmaxDOF> SST(flexdof,flexdof);
		Matrix3D jac;
		for (IntegrationPointsIterator ip(integrationRuleMass); !ip.IsEnd(); ++ip)
		{
			GetJacobi(jac, ip.Point());
			double jacdet = jac.Det();
			// compute SF = vector containing shape functions
			for (int i=1; i<=ns; i++)
			{
				SF(i)=GetS0(ip.Point(), i);
			}

			// compute S * S^T  where S is the shape function matrix 
			//     [ SF(1) 0 0 SF(2) 0 0 ...
			// S = [ 0 SF(1) 0 0 SF(2) 0 ..
			//     [ 0 0 SF(1) 0 0 SF(2) ..
			SST.SetAll(0);
			for (int i=0; i<ns; i++)
				for (int j=0; j<ns; j++)
					for (int di=1; di<=dim; di++)
							SST(i*dim+di,j*dim+di) += SF(i+1)*SF(j+1);

			m.AddSubmatrix(SST,1,1,1,1,flexdof,flexdof,
				fabs(jacdet) * GetMaterial().Density() * ip.Weight());
			//massmatrix += fabs (jacdet) * GetMaterial().Density() * w(i1) * (HL*HL.GetTp());
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

void FiniteElement3D::GetIntDuDqFCentrifugal(Matrix& H, const Vector3D& omega, const Vector3D& r0)
{	
	int dim = Dim();
	int ns = NS();
	int sos = SOS();

	H.SetSize(sos,dim);
	H.SetAll(0);
	Matrix3D jac;
	for (IntegrationPointsIterator ip(integrationRuleLoad); !ip.IsEnd(); ++ip)
	{
		GetJacobi(jac, ip.Point());
		double jacdet = jac.Det();
		double fact = fabs (jacdet) * ip.Weight();
		Vector3D r_rel = r0 - GetPos(ip.Point());
		Vector3D force = omega.Cross(omega.Cross(r_rel));

		for (int i = 0; i < ns; i++)
		{
			double combinedfactor = fact*GetS0(ip.Point(), i+1)*Rho();
			H(i*dim+1,1) += combinedfactor * force.X();
			H(i*dim+2,1) += combinedfactor * force.Y();
			H(i*dim+3,1) += combinedfactor * force.Z();
		}
	}
}

void FiniteElement3D::GetH(Matrix& H) 
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
		ConstMatrix<(FEmaxDOF+FFRFsize)*3> dudq;
		Matrix3D jac;
		// flexible part Hf
		for (IntegrationPointsIterator ip(integrationRuleLoad); !ip.IsEnd(); ++ip)
		{
			GetJacobi(jac, ip.Point());
			double jacdet = jac.Det();
			double fact = fabs(jacdet) * ip.Weight();
			for (int i=0; i<ns; i++)
			{
				for (int j=1; j<=dim; j++)
				{
					H(i*dim+j,j) += fact * GetS0(ip.Point(), i+1);
				}
			}
		}
		Hmatrix = H;
	}
}


void FiniteElement3D::GetdPosdqT(const Vector3D& ploc, Matrix& d)
{
	d.SetSize(SOS(),Dim());
	d.FillWithZeros();

	for (int i = 1; i <= NS(); i++)
	{
		double sv =	GetS0(ploc, i);
		d((i-1)*Dim()+1,1) = sv;
		d((i-1)*Dim()+2,2) = sv;
		d((i-1)*Dim()+3,3) = sv;
	}
}


//fill in sos x sos components, m might be larger
//compute stiffness matrix
void FiniteElement3D::StiffnessMatrix(Matrix& m) 
{
	int small_strain_assumption = (GetGeometricNonlinearityStatus() != GNS_NonlinearLargeStrain);

	//linear stiffness matrix: for linear problems or component mode synthesis:
	int dim = Dim(); 
	int ns = NS();

	if (
		mbs->GetSolSet().store_FE_matrices &&
		small_strain_assumption &&
		!GetMaterial().IsInelasticMaterial() &&      //$ PG 2013-8-14: assemble every stiffness matrix in case of inelastic material
		stiffnessmatrix.Getcols() == FlexDOF()
		)
	{
		m.SetSubmatrix(stiffnessmatrix,1,1);
		return;
	}

	m.SetAll(0);

	// B = deps/dq in the linear case, in the nonlinear case B = dE/dq
	ConstMatrix<6*FEmaxDOF> B((dim-1)*3,dim*ns), CB((dim-1)*3,dim*ns); //filled with zeros
	// intermediate quantities
	ConstMatrix<FEmaxDOF*FEmaxDOF> BCB(dim*ns,dim*ns); 

	// Elasticity Matrix
	ConstMatrix<36> C;
	GetMaterial().ComputeElasticityMatrix(C);

	// needed only in case of geometric nonlinear:
	Matrix3D F, sigma, strain;  // deformation gradient F, second piola kirchhoff tensor, nonlinear strain tensor 1/2(F^T F - I)
	// dB_nl/dq^T sigma (B_nl is the nonlinear part of B in case of geometric nonlinear problem, see Skript p. 38
	ConstMatrix<FEmaxDOF*FEmaxDOF> dB_nldq_sigma; // Contains dB_nldq_sigma(3*(alpha-1)+j, 3*(beta-1)+k) = dB_nl(i, 3*(alpha-1)+j)/dq(3*(beta-1)+k)*stress(i)

	ConstVector<FEmaxDOF> xgCache;
	if (GetGeometricNonlinearityStatus() != GNS_Linear || GetMaterial().IsInelasticMaterial())
	{
		SetComputeCoordinates(xgCache);
	}

	// Integration routine
	for (IntPointsStiffnessIterator ip(this); !ip.IsEnd(); ++ip)
	{
		//grad contains derivatives in compressed form:
		//grad(1,1..ns) = da/dx     where a = du/dqi or dv/dqi
		//grad(2,1..ns) = da/dy
		Matrix dShapei_dxj;
		GetGrad(ip, dShapei_dxj);

		if(small_strain_assumption) // Small strain
		{
			B.SetAll(0);

			// Geometrically linearized B matrix
			for (int alpha = 1; alpha <= ns; alpha++)
			{
				B(1,3*(alpha-1)+1) = dShapei_dxj(1,alpha); //du/dx
				B(2,3*(alpha-1)+2) = dShapei_dxj(2,alpha); //dv/dy
				B(3,3*(alpha-1)+3) = dShapei_dxj(3,alpha); //dw/dz

				B(4,3*(alpha-1)+2) = dShapei_dxj(3,alpha); //dv/dz, eps_yz
				B(4,3*(alpha-1)+3) = dShapei_dxj(2,alpha); //dw/dy
				B(5,3*(alpha-1)+3) = dShapei_dxj(1,alpha); //dw/dx, eps_xz
				B(5,3*(alpha-1)+1) = dShapei_dxj(3,alpha); //du/dz
				B(6,3*(alpha-1)+1) = dShapei_dxj(2,alpha); //du/dy, eps_xy
				B(6,3*(alpha-1)+2) = dShapei_dxj(1,alpha); //dv/dx
			}
			// compute consistent tangent stiffness matrix instead of pure elasticity matrix C
			if (GetMaterial().IsInelasticMaterial() && GetMaterial().GetInelasticitySolutionMethod() == Material::ISM_ConsistentTangentStiffness)
			{
				// compute linearized strain tensor
				ConstVector<6> strain_vector = (ConstVector<6>&) (B*xgCache);   // (eps_xx, eps_yy, eps_zz, 2 eps_yz, 2 eps_zx, 2 eps_xy)

				// get inelastic variables
				ConstVector<MAXInelasticVariablesCount> inelastic_variables;
				DataLastStepToInelasticVariables(inelastic_variables, ip.GetIndex());   // (eps^p_xx, eps^p_yy, eps^p_zz, 2 eps^p_yz, 2 eps^p_zx, 2 eps^p_xy, hardening_parameter)
				GetMaterial().ComputeConsistentTangentStiffnessMatrix(C, strain_vector, inelastic_variables);  // C as input: elasticity tensor, as output: consistent tangent matrix C = C (I - 2 mu d actual_plastic_strain/d strain)
			}
		}
		else // Large Strain
		{
			B.SetAll(0.);
			dB_nldq_sigma.SetSize(dim*ns,dim*ns);  
			dB_nldq_sigma.SetAll(0);

			// compute actual displacement gradient
			Gradu(xgCache, dShapei_dxj, F);
			// add I matrix to obtain deformation gradient
			F(1,1) += 1; F(2,2) += 1; F(3,3) += 1;

			// Green-Lagrange strain tensor
			//strain = 0.5 * (F.GetTp() * F);
			F.GetATA2(strain);
			strain(1,1) -= 0.5; strain(2,2) -= 0.5; strain(3,3) -= 0.5;

			if (IsInelasticMaterial())
			{
				mbs->UO(UO_LVL_warn) << "WARNING: Large strain inelasticity not implemented!\n";
			}
			
			GetMaterial().ComputeStressFromStrain(strain, sigma);

			for (int alpha = 1; alpha <= ns; alpha++)
			{
				// nonlinear B matrix
				// Variation of nonlinear Green Lagrange strain tensor E
				for (int k=1; k<=dim; k++)
				{
					int alphak = 3*(alpha-1)+k;
					B(1,alphak) = F(k,1)*dShapei_dxj(1,alpha);   // duk/dx1 * dShape_alpha/dx1
					B(2,alphak) = F(k,2)*dShapei_dxj(2,alpha);   // duk/dx2 * dShape_alpha/dx2
					B(3,alphak) = F(k,3)*dShapei_dxj(3,alpha);   // duk/dx3 * dShape_alpha/dx3
					B(4,alphak) = F(k,2)*dShapei_dxj(3,alpha)+F(k,3)*dShapei_dxj(2,alpha);  // duk/dx2 * dShape_alpha/dx3 + duk/dx3 * dShape_alpha/dx2
					B(5,alphak) = F(k,1)*dShapei_dxj(3,alpha)+F(k,3)*dShapei_dxj(1,alpha);  
					B(6,alphak) = F(k,2)*dShapei_dxj(1,alpha)+F(k,1)*dShapei_dxj(2,alpha);  
				}

				//// entries of dB_nl/dq^T sigma
				// dB_nldq^T sigma (alpha j, beta j) = Sum_k dBnl(k, alphaj) / dq(betaj) * sigma(k) = add_dBdq_sigma is independent of j
				for (int beta=1; beta<=ns; beta++)
				{
					double add_dBdq_sigma = 0;
					add_dBdq_sigma += dShapei_dxj(1,beta)*dShapei_dxj(1,alpha) * sigma(1,1);
					add_dBdq_sigma += dShapei_dxj(2,beta)*dShapei_dxj(2,alpha) * sigma(2,2);
					add_dBdq_sigma += dShapei_dxj(3,beta)*dShapei_dxj(3,alpha) * sigma(3,3);
					add_dBdq_sigma += (dShapei_dxj(2,beta)*dShapei_dxj(3,alpha)+dShapei_dxj(3,beta)*dShapei_dxj(2,alpha)) * sigma(3,2);
					add_dBdq_sigma += (dShapei_dxj(3,beta)*dShapei_dxj(1,alpha)+dShapei_dxj(1,beta)*dShapei_dxj(3,alpha)) * sigma(1,3);
					add_dBdq_sigma += (dShapei_dxj(1,beta)*dShapei_dxj(2,alpha)+dShapei_dxj(2,beta)*dShapei_dxj(1,alpha)) * sigma(1,2);

					for (int j=1; j<=dim; j++)
					{
						int alphaj = 3*(alpha-1)+j;
						int betaj = 3*(beta-1)+j;
						dB_nldq_sigma(alphaj, betaj) += add_dBdq_sigma; 
					}
				}
			}
			if (GetMaterial().IsInelasticMaterial() && GetMaterial().GetInelasticitySolutionMethod() == Material::ISM_ConsistentTangentStiffness)
			{
				GetMBS()->UO(UO_LVL_warn) << "WARNING: FiniteElement3D: Computation of Consistent Tangent Stiffness Matrix not implemented for large strain case!\n";
			}
		}

		// compute B^T C B
		Mult(C, B, CB);
		MultTp(B, CB, BCB);

		// Stiffness matrix linear:    K = B^T C B,
		//                  nonlinear: K = B_hat ^T C B_hat + dB_hat/dq sigma

		// integration weight factor
		double fact = -fabs(ip.Data().jacdet) * ip.Weight();

		// Add part  B^T C B 
		// AddSubmatrix(mat, offset rows submat, offset cols submat, offset rows mat, offset cols mat, # rows, # cols, factor)
		m.AddSubmatrix(BCB, 1, 1, 1, 1, ns*dim, ns*dim, fact);

		// geometrically nonlinear case:
		// Add part  dB/dq^T  sigma
		if (!small_strain_assumption)
		{
			m.AddSubmatrix(dB_nldq_sigma, 1, 1, 1, 1, ns*dim, ns*dim, fact);
		}

	}

	if (mbs->GetSolSet().store_FE_matrices && small_strain_assumption && !(GetMaterial().IsInelasticMaterial() && GetMaterial().GetInelasticitySolutionMethod() == Material::ISM_ConsistentTangentStiffness)) // if stiffness matrix is stored
	{
		stiffnessmatrix.SetSize(FlexDOF(), FlexDOF());
		for (int i = 1; i <= FlexDOF(); i++)
			for (int j = 1; j <= FlexDOF(); j++)
				stiffnessmatrix(i,j) = m(i,j);
	}
}



//$EK 2013-02-12: compute volume of Finite element via integration of 1
double FiniteElement3D::ComputeVolume()
{
	double volume = 0.;
	// Integration routine
	for (IntPointsStiffnessIterator ip(this); !ip.IsEnd(); ++ip)
	{
		double fact = fabs(ip.Data().jacdet) * ip.Weight();
		volume += fact;
	}
	return fabs(volume);
}

// $EK 2013-03-21 compute center of mass via integration
void FiniteElement3D::ComputeCenterOfMass(Vector3D& center)
{
	center.Set(0,0,0);
	// Integration routine
	for (IntPointsStiffnessIterator ip(this); !ip.IsEnd(); ++ip)
	{
		double fact = fabs(ip.Data().jacdet) * ip.Weight();
		center += fact* GetPosRel(ip.Point());
		
	}
	center *= 1./this->GetVolume();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//COMPUTE RESIDUAL       COMPUTE RESIDUAL       COMPUTE RESIDUAL       COMPUTE RESIDUAL       COMPUTE RESIDUAL       COMPUTE RESIDUAL
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void FiniteElement3D::EvalF2(Vector& f, double t) 
{
	//add loads, constraint forces and other things in parent function:
	Body3D::EvalF2(f,t);
	TMStartTimer(22);

	int sos = SOS();
	int ns = NS();
	int dim = Dim();

	ConstVector<FEmaxDOF> fadd;
	fadd.SetLen(sos);
	fadd.SetAll(0);

	if ( GetGeometricNonlinearityStatus() == GNS_Linear || IsFFRF())
	{
		EvalF2GeomLin(fadd, t);
	}
	else
	{
		EvalF2GeomNonlin(fadd, t);
	}

	if (IsFFRF())
	{
		AddQuadraticVelocityVector(fadd, t);
	}	

	f -= fadd;
	TMStopTimer(22);

	// +++++++ damping: +++++++
	if (GetMaterial().DampingM() != 0)
	{
		// set velocities into a temporary vector, modified by YV:
		ConstVector<FEmaxDOF> velocities;
		velocities.SetLen(sos);
		for (int i = 1; i <= sos; i++)
			velocities(i) = XGP(i);

		ConstVector<FEmaxDOF> temp;
		temp.SetLen(sos);
		temp.SetAll(0.);

		if (mbs->GetSolSet().store_FE_matrices) // mass matrix is stored
		{
			if (massmatrix.Getcols() != sos) {UO() << "ERROR: finite element: mass matrix not built!\n"; return;}
			Mult(massmatrix,velocities,temp);
		}
		else // mass matrix is not stored
		{
			ConstMatrix<FEmaxDOF*FEmaxDOF> m(FlexDOF(),FlexDOF());
			EvalMff(m,t);
			Mult(m,velocities,temp);
		}
		temp *= GetMaterial().DampingM();
		f -= temp;
	}
}; 


void FiniteElement3D::EvalF2GeomNonlin(Vector& fadd, double t) 
{
	//add loads, constraint forces and other things in parent function:
	int sos = SOS();
	int flexdof = FlexDOF();

	int ns = NS();
	int dim = Dim();

	ConstVector<FEmaxDOF> xgCache;
	SetComputeCoordinates(xgCache);

	Matrix3D strain, piola1, F;

	ConstVector<FEmaxDOF> temp;
	temp.SetLen(sos);
	temp.SetAll(0.);

	//compute deformation gradient, compute strain, compute stress and compute elastic forces:
	for (IntPointsStiffnessIterator ip(this); !ip.IsEnd(); ++ip)
	{
		//int_V (P : d F / d q) dV
		//P = F * S
		//S = D : E
		//E = 1/2*(F^T F-I)

		Matrix agrad;
		GetGrad(ip, agrad);

		Gradu(xgCache, agrad, F);
		F(1,1) += 1; F(2,2) += 1; F(3,3) += 1;

		// Green-Lagrange strain tensor
		//strain = 0.5 * (F.GetTp() * F);
		F.GetATA2(strain);
		strain(1,1) -= 0.5; strain(2,2) -= 0.5; strain(3,3) -= 0.5;

		// get last valid inelastic variables
		// either from last time step (consistent tangent stiffness matrix) or from last nonlin step (return mapping)
		Vector inelastic_variables(0);
		if (IsInelasticMaterial())
		{
			Material::InelasticitySolutionMethod solmeth = GetMaterial().GetInelasticitySolutionMethod();
			if (solmeth == Material::ISM_ConsistentTangentStiffness || solmeth == Material::ISM_ReturnMapping)
			{
				DataLastStepToInelasticVariables(inelastic_variables, ip.GetIndex());
			}
			else if (solmeth == Material::ISM_Default || solmeth == Material::ISM_FixedPoint)
			{
				DataToInelasticVariables(inelastic_variables, ip.GetIndex());
			}
			else
			{
				mbs->UO(UO_LVL_warn) << "WARNING: FiniteElement3D: Inelasticity-solution method " << solmeth << " not specified for this Element!\n";
			}
		}

		//{ initial strain - YV
		if(HasInitialStrains())
		{
			Matrix3D initial_strain;
			GetInitialStrain(initial_strain, ip.Point(), false);
			GetMaterial().ComputeStressFromStrain(strain, piola1, inelastic_variables, initial_strain);
		}
		//}
		else
		{
			GetMaterial().ComputeStressFromStrain(strain, piola1, inelastic_variables);
		}
		piola1 = F * piola1;  //piola1 = F * piola2

		for (int j=1; j <= Dim(); j++)
		{
			for (int i = 1; i <= ns; i++)
			{
				temp(3*(i-1)+j) = agrad(1, i)*piola1(j,1)
					+ agrad(2, i)*piola1(j,2)
					+ agrad(3, i)*piola1(j,3);
			}
		}
		fadd.MultAdd(fabs(ip.Data().jacdet) * ip.Weight(), temp);
	}
}; 


void FiniteElement3D::EvalF2GeomLin(Vector& fadd, double t) 
{
	//add loads, constraint forces and other things in parent function:
	int sos = SOS();
	int flexdof = FlexDOF();

	int ns = NS();
	int dim = Dim();

	ConstVector<FEmaxDOF> xgCache;
	SetComputeCoordinates(xgCache);

	//double u;

	Matrix3D strain, stress, gradu, F;

	ConstVector<FEmaxDOF> temp;
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
		//EK 2012-06-27 is it necessary in the linear case to compute F ??? -> NO!
		//F = gradu;
		//F(1,1) += 1.; F(2,2) += 1.; F(3,3) += 1.;

		// Green-Lagrange strain tensor
		strain = 0.5 * (gradu.GetTp() + gradu);

		Vector inelastic_variables(0);
		if (IsInelasticMaterial())
		{
			Material::InelasticitySolutionMethod solmeth = GetMaterial().GetInelasticitySolutionMethod();
			if (solmeth == Material::ISM_ConsistentTangentStiffness || solmeth == Material::ISM_ReturnMapping)
			{
				DataLastStepToInelasticVariables(inelastic_variables, ip.GetIndex());    //$ PG 2012-10-24: replaced DataToInelasticVariables --> DataLastStepToInelasticVariables
			}
			else if (solmeth == Material::ISM_Default || solmeth == Material::ISM_FixedPoint)
			{
				DataToInelasticVariables(inelastic_variables, ip.GetIndex());
			}
			else
			{
				mbs->UO(UO_LVL_warn) << "WARNING: FiniteElement3D: Inelasticity-solution method " << solmeth << " not specified for this Element!\n";
			}
		}
		
		//{ initial strain - YV
		if(HasInitialStrains())
		{
			Matrix3D initial_strain;
			GetInitialStrain(initial_strain, ip.Point(), false);
			GetMaterial().ComputeStressFromStrain(strain, stress, inelastic_variables, initial_strain);
		}
		//}
		else
		{
			GetMaterial().ComputeStressFromStrain(strain, stress, inelastic_variables);
		}

		for (int j=1; j <= Dim(); j++)
		{
			for (int i = 1; i <= ns; i++)
			{
				int k = Dim()*(i-1)+j;
				temp(k) = agrad(1, i)*stress(j,1)
					+ agrad(2, i)*stress(j,2)
					+ agrad(3, i)*stress(j,3);
				fadd(k) += fabs(ip.Data().jacdet) * ip.Weight() * temp(k);
			}
		}
	}
}; 

double FiniteElement3D::PostNewtonStep(double t)
{

	if (!IsInelasticMaterial()) return 0;
	Material::InelasticitySolutionMethod solmeth = GetMaterial().GetInelasticitySolutionMethod();


	double err = 0;
	double fe_volume = 0;

	int sos = FlexDOF();
	int ns = NS();
	int dim = Dim();

	ConstVector<FEmaxDOF> xgCache;
	SetComputeCoordinates(xgCache);

	Matrix3D F,jac;
	
	ConstVector<MAXInelasticVariablesCount> inelastic_variables(GetMaterial().GetInelasticVariablesCount());
	ConstVector<MAXInelasticVariablesCount> inelastic_variables_last_step(GetMaterial().GetInelasticVariablesCount());
	ConstVector<MAXInelasticVariablesCount> inelastic_variables_update(GetMaterial().GetInelasticVariablesCount());
	//----------------------
	// loop over integration points
	for (IntPointsStiffnessIterator ip(this); !ip.IsEnd(); ++ip)
	{
		// precompute average value used for plastic strain update (yield criterion)
		//if (solmeth == Material::ISM_FixedPoint || solmeth == Material::ISM_ReturnMapping)
		//{
		//	GetMaterial().SetAveragedYieldParamter(PostNewtonAveragingStep(t));
		//}

		// Jacobi-matrix and Jacobi-determinant for integration point p
		GetJacobi(jac,ip.Point());
		double jacdet = jac.Det();

		// compute F = grad u + I
		Matrix agrad;
		GetGrad(ip, agrad);
		Gradu(xgCache, agrad, F);
		F(1,1) += 1; F(2,2) += 1; F(3,3) += 1;

		// compute strain from F
		Matrix3D strain;
		if ( GetGeometricNonlinearityStatus() == GNS_Linear || IsFFRF())
		{
			strain = 0.5*(F.GetTp() + F) - Matrix3D(1.);
		}
		else
		{
			strain = 0.5*(F.GetTp()*F - Matrix3D(1.));
		}

		// potentially subtract initial strain
		Matrix3D initial_strain;
		if(HasInitialStrains())
		{
			GetInitialStrain(initial_strain, ip.Point(), false);
		}

		// get old inelastic strain (from last time step)
		DataLastStepToInelasticVariables(inelastic_variables_last_step, ip.GetIndex());  // read old inelastic strain from last computed time/load step
		
		// and last computed update of inelastic strain (either from last nonlin step (return mapping) or also from last time step (consistent tangent stiffness matrix))
		switch (solmeth)
		{
		case Material::ISM_ConsistentTangentStiffness:
		case Material::ISM_ReturnMapping:
			DataLastStepToInelasticVariables(inelastic_variables, ip.GetIndex());  // read old inelastic strain from last computed time/load step
			break;
		case Material::ISM_FixedPoint:
			DataToInelasticVariables(inelastic_variables, ip.GetIndex());  // read old inelastic strain from last nonlinear (discontinuous) step
			break;
		default:
			mbs->UO(UO_LVL_warn) << "WARNING: in FiniteElement3D::PostNewtonStep(double t) not implemented for inelasticity-solution-method " << solmeth << "\n";
		}

		// compute actual inelastic strain as function from old inelastic strain and elastic strain (minus initial strain)
		err += fabs(jacdet) * ip.Weight() * GetMaterial().ComputeInelasticVariablesUpdateFromStrain(inelastic_variables_update, inelastic_variables, strain - initial_strain);
		fe_volume += fabs(jacdet) * ip.Weight();

		inelastic_variables = inelastic_variables_last_step + inelastic_variables_update;
		InelasticVariablesToData(inelastic_variables, ip.GetIndex());
	}

	return err/fe_volume;
}

void FiniteElement3D::PostprocessingStep()
{
}


void FiniteElement3D::DrawElementPreProc() 
{
	if (!GetMBS()->GetIOption(118)) return;

	if (HasOuterFaces() || !GetMBS()->GetIOption(146)||GetMBS()->GetIOption(149)) //do not average stress at nodes if element does not have visible faces //$ DR 2012-09-20: bugfix, added GetMBS()->GetIOption(149)
	{
		FieldVariableDescriptor * fvd = GetMBS()->GetActualPostProcessingFieldVariable();
		if(fvd != NULL)
		{
			for (int i=1; i <= NNodes(); i++)
			{
				GetNode(i).GetDrawTemp() += GetFieldVariableValue(*fvd, (1./sqrt(3.))*GetNodeLocPos(i), true);
				GetNode(i).GetDrawTempCnt()++;
			}
		}
	}
};

void FiniteElement3D::DrawElement() 
{
	//$ YV 2012-12-11: the treatment of cutting planes and outer faces was placed here from mbs
	char backUpOfOuterFaces = Outer_Face();

	int no_surfaceupdate = mbs->GetIOption(119);
	if (mbs->UseCuttingPlanes() && no_surfaceupdate == 0)

	SetOuterFacesCuttingPlane();

	//$ DR 2011-12-02:[ do not draw element if no surface at all
	bool draw_this_element = 0;
	if(!(GetMBS()->GetIOption(146))) 
	{
		draw_this_element = 1;		// not only surface elements are drawn
	}
	else												// only surface elements are drawn
	{
		for (int side=1; (side <= NFaces())&&(!draw_this_element) ; side++)
		{
			if (GetOuterFaceFlag(side)) 
			{
				draw_this_element = 1;
			}
		}
	}
	//$ DR 2011-12-02:] do not draw element if no surface at all
	if(draw_this_element)
	{

		mbs->SetColor(col);

		//SetDrawCoordinates(); //for computation of GraduD; 12.4.9 ***not necessary any more!

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

		//if(col==Vector3D(0.6,0.6,0.6)) colormode = 0;		//$!DR 2011-05-26: Hack: if body is grey: draw body color, for all the other bodies: draw mises or what ever else is defined

		int res = GetMBS()->GetDrawResolution()-1;
		if (res == -1) res = 0;

		double tilex = pow(2.,res);
		static TArray<Vector3D> points;
		points.SetLen((int)(tilex+1)*(int)(tilex+1));
		static TArray<double> vals;
		vals.SetLen((int)(tilex+1)*(int)(tilex+1));

		double v=0;
		Vector3D p0, vx, vy;

		double shrink = GetMBS()->GetDOption(106);//0.995;

		for (int side=1; side <= NFaces(); side++)
		{
			if (GetOuterFaceFlag(side) || !GetMBS()->GetIOption(146)) //draw surface elements only?
			{
				points.SetLen(0); vals.SetLen(0);
				Vector3D p0, vx, vy, vxy;

				GetLocalFaceCoordinates(side, p0, vx, vy, vxy);
				vx /= tilex;
				vy /= tilex;
				vxy /= tilex;

				if (shrink != 1) 
				{
					if (p0.Norm() != 0)
						p0 *= shrink; //hexahedral
					else
						p0 = Vector3D(0.5*(1.-shrink)); //tetrahedral

					vx *= shrink;
					vy *= shrink;
				}

				for (double iy = 0; iy <= tilex; iy++)
				{
					for (double ix = 0; ix <= tilex; ix++)
					{
						Vector3D ploc;
						if (GetElementType() == TFE_Tetrahedral)
						{
							ploc = (p0+ix*vx*(1.-iy/tilex)+iy*vy);
							if (side == 4) ploc = (p0+(1.-ix/tilex)*tilex*vx*(1.-iy/tilex)+ix*vxy*(1.-iy/tilex)+iy*vy);
						}
						// $EK 2013-03-05 special rule for prism
						else if (GetElementType() == TFE_Prism)
						{
							if (side <= 3) // quadrilateral sides
								ploc = (p0+ix*vx+iy*vy);
							else // triangluar sides
								//ploc = (p0+ix*vx*(1.-iy/tilex)+iy*vy);
								ploc = (p0+ix*vx*(1.-iy/tilex)+iy*vy);
						}
						// $EK 2013-03-05 special rule for pyramid
						else if (GetElementType() == TFE_Pyramid)
						{
							if (side <= 1) // quadrilateral sides
								ploc = (p0+ix*vx+iy*vy);
							else // triangluar sides
								//ploc = (p0+ix*vx*(1.-iy/tilex)+iy*vy);
								ploc = (p0+ix*vx*(1.-iy/tilex)+iy*vy);
						}
						else
						{
							ploc = (p0+ix*vx+iy*vy);
						}
						points.Add(GetPosD(ploc));


						if (colormode)
						{
							if (!GetMBS()->GetIOption(118)) //compute stress at every local position
							{
								if (colormode)
									v = GetFieldVariableValue(*GetMBS()->GetActualPostProcessingFieldVariable(), ploc, true);
							}
							else //stress interpolated from nodes!
							{
								v = 0;
								for (int i=1; i <= NNodes(); i++)
								{
									v += GetS0(ploc, i)*GetNode(i).GetDrawTemp();
								}
							}
						}
						vals.Add(v);
					}
				}

				FitPointsInDrawingRegion(points, (int)tilex+1, (int)tilex+1);		// in case the derived class wishes to control the drawing region

				mbs->DrawColorQuads(points,vals,(int)tilex+1,(int)tilex+1,colormode,linemode);
			}
		}

	} // end: if (draw_this_element)

	Outer_Face() = backUpOfOuterFaces;
};

void FiniteElement3D::SetOuterFaceFlag(int face, int set /* = 1 */)
{
	if (set == 1)
	{
		outer_face = outer_face|GetCharBit(face);
	}
	else
	{
		if (GetOuterFaceFlag(face))
		{
			outer_face = (char)((int)outer_face - GetIntBit(face));
		}
	}
}

void FiniteElement3D::GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables)
{
	// add all the field variables of the parent class
	Body3D::GetAvailableFieldVariables(variables);
	//FVT_position,FVT_displacement, FVT_velocity: implemented in Element
	// possibility to remove entries does not exist yet, just do not use the Function of the parent class

	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_stress,
		FieldVariableDescriptor::FVCI_z, true);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_stress_mises);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_total_strain,
		FieldVariableDescriptor::FVCI_z, true);
	if(GetMaterial().IsInelasticMaterial())
	{
		FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_inelastic_strain,
			FieldVariableDescriptor::FVCI_z, true);
		if(GetMaterial().GetInelasticityType() == 1) 
		{
			FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_hardening_parameter);
			FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_yield_function);
		}
	}
}

double FiniteElement3D::GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & local_position, bool flagD)
{
	switch(fvd.VariableType())
	{
	case FieldVariableDescriptor::FVT_displacement:
		{
			if(flagD)
				return fvd.GetComponent(GetDisplacementD(local_position));
			else
				return fvd.GetComponent(GetDisplacement(local_position));
		}
	case FieldVariableDescriptor::FVT_position:
		{
			if(flagD)
				return fvd.GetComponent(GetPosD(local_position));
			else
				return fvd.GetComponent(GetPos(local_position));
		}
	case FieldVariableDescriptor::FVT_velocity:
		{
			if(flagD)
				return fvd.GetComponent(GetVelD(local_position));
			else
				return fvd.GetComponent(GetVel(local_position));
		}
	case FieldVariableDescriptor::FVT_stress:
	case FieldVariableDescriptor::FVT_stress_mises:
		{
			Vector inelastic_variables;
			if(IsInelasticMaterial())
			{
				IntPointsStiffnessIterator ip(this);
				ip.GoClosestTo(local_position);
				if(flagD)
				{
					DataToInelasticVariablesD(inelastic_variables, ip.GetIndex());
				}
				else
				{
					DataToInelasticVariables(inelastic_variables, ip.GetIndex());
				}
			}

			Matrix3D gradu;
			if(flagD) 
				GraduD(local_position, gradu);
			else
			{
				ConstVector<FEmaxDOF> xgCache;
				SetComputeCoordinates(xgCache);
				Gradu(local_position, xgCache, gradu);
			}
			
			Matrix3D strain;
			GetMaterial().ComputeStrain(strain, gradu, GetGeometricNonlinearityStatus() == GNS_NonlinearLargeStrain);
			
			Matrix3D stress;
			if(HasInitialStrains())
			{
				Matrix3D initial_strain;
				GetInitialStrain(initial_strain, local_position, flagD);
				GetMaterial().ComputeStressFromStrain(strain, stress, inelastic_variables, initial_strain);
			}
			else
			{
				GetMaterial().ComputeStressFromStrain(strain, stress, inelastic_variables);
			}
			
			if(fvd.VariableType() == FieldVariableDescriptor::FVT_stress_mises)
				return stress.Mises();
			else
				return fvd.GetComponent(stress);
		}
	case FieldVariableDescriptor::FVT_total_strain:
		{
			Matrix3D gradu;
			if(flagD) 
				GraduD(local_position, gradu);
			else
			{
				ConstVector<FEmaxDOF> xgCache;
				SetComputeCoordinates(xgCache);
				Gradu(local_position, xgCache, gradu);
			}
			Matrix3D strain;
			GetMaterial().ComputeStrain(strain, gradu, GetGeometricNonlinearityStatus() == GNS_NonlinearLargeStrain);
			return fvd.GetComponent(strain);
		}
	case FieldVariableDescriptor::FVT_initial_strain:
		{
			Matrix3D initial_strain;
			GetInitialStrain(initial_strain, local_position, flagD);
			return fvd.GetComponent(initial_strain);
		}
	case FieldVariableDescriptor::FVT_inelastic_strain:
	case FieldVariableDescriptor::FVT_hardening_parameter:
	case FieldVariableDescriptor::FVT_yield_function:
		{
			if (!IsInelasticMaterial())
			{
				GetMBS()->UO(UO_LVL_err) << "ERROR: Inelastic field variables requested. Material is not declared to be inelastic, though!\n";
				return 0.;
			}

			Vector inelastic_variables;
			IntPointsStiffnessIterator ip(this);
			ip.GoClosestTo(local_position);
			if (flagD)
			{
				DataToInelasticVariablesD(inelastic_variables, ip.GetIndex());
			}
			else
			{
				DataToInelasticVariables(inelastic_variables, ip.GetIndex());
			}

			switch(fvd.VariableType())
			{
			case FieldVariableDescriptor::FVT_inelastic_strain:
				{
					Matrix3D inelastic_strain;
					if (IsInelasticMaterial())
						StrainVectorToMatrix3D(inelastic_strain, inelastic_variables);
					return fvd.GetComponent(inelastic_strain);
				}
			case FieldVariableDescriptor::FVT_hardening_parameter:
				{
					return inelastic_variables(7);
				}
			case FieldVariableDescriptor::FVT_yield_function:
				{
					Vector3D new_local_position = ip.Point();
					Matrix3D gradu;
					if(flagD) 
						GraduD(new_local_position, gradu);
					else
					{
						ConstVector<FEmaxDOF> xgCache;
						SetComputeCoordinates(xgCache);
						Gradu(new_local_position, xgCache, gradu);
					}
					Matrix3D strain;
					GetMaterial().ComputeStrain(strain, gradu, GetGeometricNonlinearityStatus() == GNS_NonlinearLargeStrain);
					Matrix3D initial_strain;
					GetInitialStrain(initial_strain, new_local_position, flagD);
					return GetMaterial().ComputeYieldFunction(strain-initial_strain, inelastic_variables);
				}
			}
		}
	}
	
	return FIELD_VARIABLE_NO_VALUE;
}

void FiniteElement3D::PreAssemble()
{
	FiniteElementGeneric<Body3D>::PreAssemble();
	SetOuterFaces();
}