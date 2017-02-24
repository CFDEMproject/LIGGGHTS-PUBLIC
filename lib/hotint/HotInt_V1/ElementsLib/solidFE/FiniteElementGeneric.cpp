//#**************************************************************
//#
//# filename:             FiniteElementGeneric.cpp
//#
//# author:               Gerstmayr Johannes, Aigner Larissa, YV
//#
//# generated:						October 2010
//# description:          common finite element base class - for 2D & 3D finite elements
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
#include "body2d.h"
#include "material.h"
#include "node.h"
#include "finiteelementgeneric.h"
#include "elementdataaccess.h"
#include "solversettings_auto.h"

// explicit instantatiation of the template
template class FiniteElementGeneric<Body3D>;
template class FiniteElementGeneric<Body2D>;

template<class BodyXD>
void FiniteElementGeneric<BodyXD>::CopyFrom(const Element& e)
{
	BodyXD::CopyFrom(e);
	const FiniteElementGeneric<BodyXD>& ce = (const FiniteElementGeneric<BodyXD>&)e;
	massmatrix = ce.massmatrix;
	Hmatrix = ce.Hmatrix;
	stiffnessmatrix = ce.stiffnessmatrix;

	// integration rules
	integrationRuleStiffness = ce.integrationRuleStiffness;
	integrationRuleMass = ce.integrationRuleMass;
	integrationRuleLoad = ce.integrationRuleLoad;

	SetGeometricNonlinearityStatus(ce.GetGeometricNonlinearityStatus());

	nodes.SetLen(ce.NNodes());
	for (int i = 1; i <= ce.NNodes(); i++)
	{
		NodeNum(i) = ce.NodeNum(i);
		//concentratedmass[i-1] = ce.concentratedmass[i-1];
	}

	bodyind = ce.bodyind;
}

template<class BodyXD>
void FiniteElementGeneric<BodyXD>::GetElementData(ElementDataContainer& edc)
{
	BodyXD::GetElementData(edc);
	ElementData ed;
	IVector v;
	for (int i=1; i <= NNodes(); i++)
	{
		v.Add(NodeNum(i));
	}
	SetElemDataIVector(edc, v, "Node_list");
	ed.SetInt(bodyind, "Body_index");
	edc.Add(ed);
}

template<class BodyXD>
int FiniteElementGeneric<BodyXD>::SetElementData(ElementDataContainer& edc)
{
	int rv = BodyXD::SetElementData(edc);
	ElementData ed;
	IVector nodelist(NNodes());
	GetElemDataIVector(GetMBS(), edc, "Node_list", nodelist, 1);
	GetElemDataInt(GetMBS(), edc, "Body_index", bodyind);
	SetFiniteElementGeneric(bodyind, nodelist, GetMaterialNum(), GetCol());
	return rv;
}

template<class BodyXD>
void FiniteElementGeneric<BodyXD>::LinkToElements()
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
		// remaining, if exist
		for (int i=1; i <= Nnonode; i++)
		{
			AddLTG(storeltg(i+Nnonode));
		}
	}
}

template<class BodyXD>
void FiniteElementGeneric<BodyXD>::SetMaterialNum(int mnum)
{
	materialnum = mnum;
	//$ YV 2011-08.02 moved the following line to PreAssemble():
	//if (GetMBS()->GetMaterial(materialnum).IsInelasticMaterial())
}

template<class BodyXD>
const IntegrationRule * FiniteElementGeneric<BodyXD>::GetIntegrationRule(IntegrationRule::IntegratedValueType ivt) const
{
	switch(ivt)
	{
	case IntegrationRule::IVT_Stiffness: return integrationRuleStiffness;
	case IntegrationRule::IVT_Mass: return integrationRuleMass;
	case IntegrationRule::IVT_Load: return integrationRuleLoad;
	}
	assert(0 && "unknown entity to integrate");
	return NULL;
}

template<class BodyXD>
void FiniteElementGeneric<BodyXD>::BuildDSMatrices()
{
	// precompute DShape Matrices and multiply with element transformation
	// in every integration point of the stiffness matrix
	for (IntegrationPointsIterator ip(integrationRuleStiffness); !ip.IsEnd(); ++ip)
	{
		IntegrationPointStiffnessMatrixLocalData * newData = new IntegrationPointStiffnessMatrixLocalData;
		integrationPointStiffnessMatrixLocalData(ip.GetIndex()) = newData;
		if(mbs->GetSolSet().store_FE_matrices)
		{
			newData->grad = new Matrix;
			newData->jacdet = GetJacInvDS(ip.Point(), *(newData->grad));
		}
		else
		{
			static Matrix gradTemp;
			newData->jacdet = GetJacInvDS(ip.Point(), gradTemp);
		}
	}
}

template<class BodyXD>
void FiniteElementGeneric<BodyXD>::SetComputeCoordinates(Vector & xgCache)
	{
		xgCache.SetLen(XGLength());
		for (int i = 1; i <= xgCache.Length(); i++)
			xgCache(i) = XG(i);
		
		//int nnodes = NNodes();
		//int dof_per_node = FlexDOF()/nnodes;
		//Vector3D mean;
		//for (int i=0; i<nnodes; i++)
		//{
		//	for (int j=1; j<=Dim(); j++)
		//	{
		//		mean(j) += xgCache(i*dof_per_node + j);
		//	}
		//}
		//mean *= 1./nnodes;
		//for (int i=0; i<nnodes; i++)
		//{
		//	for (int j=1; j<=Dim(); j++)
		//	{
		//		xgCache(i*dof_per_node + j) -= mean(j);
		//	}
		//}
	}

template<class BodyXD>
void FiniteElementGeneric<BodyXD>::GetNeighbourElements(TArray<int>& neighbours)
{
	neighbours.SetLen(0);

	for (int i=1; i<=NNodes(); i++)
	{
		Node& n = GetNode(i);
		for (int j=1; j<=n.NElements(); j++)
		{
			int el = n.ElemNum(j);
			if (!neighbours.Find(el) && el != GetOwnNum()) neighbours.Add(el);
		}
	}
	//UO() << "elem " << GetOwnNum() << " has " << neighbours.Length() << "neighbours\n";
}

template<class BodyXD>
void FiniteElementGeneric<BodyXD>::SetFiniteElementGeneric(int bodyindi, const TArray<int>& nodelist, int material_num, const Vector3D& coli, int CMSelementi)
{
	SetGeometricNonlinearityStatus(GNS_NonlinearSmallStrain);		// YV

	bodyind = bodyindi;
	nodes = nodelist;
	elementname = GetElementSpec();

	// materials
	SetMaterialNum(material_num);
	col = coli;

	x_init = Vector(2*FlexDOF());

	if (CMSelementi)
	{
		SetFFRFElement(CMSelementi);
	}
	//else //JG 2011-03-21: initialization not needed in case of CMS element
	//{
		// (AD)
		// get initial displacements & velocities from nodes (nodelist) if x_init of (all) nodes are initialized (Length != 0) 
		for (int i=1; i<=NNodes(); i++)
		{
			// AP: Changed mbs->GetNode(nodes(i))  to  FiniteElement3D::GetNode(i),
			//         otherwise, does not work for CMS elements
			// Node node(mbs->GetNode(nodes(i)));
			Node& node(GetNode(i));


			// {x1,y1,z1,x2,y2,z2,... vx1,...} -> twice the size
			if (node.X_Init().GetLen())
			{
				x_init.Copy(node.X_Init(), 1,								(i-1)*DOFPerNode()+1,          DOFPerNode());
				x_init.Copy(node.X_Init(), 1+DOFPerNode(),	(i-1+NNodes())*DOFPerNode()+1, DOFPerNode());
			}
		}
	//}
// (AD) end

	elementname = GetElementSpec();
}

template<class BodyXD>
void FiniteElementGeneric<BodyXD>::PreAssemble()
{
	// we set up the integration rules
	IntegrationRule::IntegrationRuleSettings settings(
		GetElementType(), GetActualInterpolationOrder(), IntegrationRule::IVT_Stiffness, GetGeometricNonlinearityStatus(), 0
		);
	integrationRuleStiffness = GetIntegrationRulesLibrary()->GetIntegrationRule(settings, this);
	settings.integratedValueType = IntegrationRule::IVT_Mass;
	integrationRuleMass = GetIntegrationRulesLibrary()->GetIntegrationRule(settings, this);
	settings.integratedValueType = IntegrationRule::IVT_Load;
	integrationRuleLoad = GetIntegrationRulesLibrary()->GetIntegrationRule(settings, this);
	assert(
		integrationRuleStiffness != NULL &&
		integrationRuleMass != NULL &&
		integrationRuleLoad != NULL
		);

	// Build matrices
	BuildDSMatrices();

	if (DataS() != 0)
	{
		Vector datainit(DataS()); //automatically initialized with zeros
		SetDataInit(datainit); //initialize data variables with zero = initial inelastic strain == zero
	}
}

template<class BodyXD>
void FiniteElementGeneric<BodyXD>::GetGrad(IntegrationPointsIterator & it, Matrix & grad) const
{
	if(
		integrationPointStiffnessMatrixLocalData.Length() < it.GetIndex() ||
		integrationPointStiffnessMatrixLocalData(it.GetIndex())->grad == NULL
		)
		GetJacInvDS(it.Point(), grad);		// if nothing is pre-computed, let's do it now
	else
	{
		Matrix * preComputedGrad = integrationPointStiffnessMatrixLocalData(it.GetIndex())->grad;
		grad.LinkWith(preComputedGrad->Getrows(), preComputedGrad->Getcols(), preComputedGrad->GetMatPtr());
	}
}

template<class BodyXD>
void FiniteElementGeneric<BodyXD>::AddMSparse(SparseMatrix& ms, double t) //add sparse mass matrix into full system matrix
{
	if (mbs->GetSolSet().store_FE_matrices || massmatrix.Getcols() != FlexDOF())
	{
		ConstMatrix<FEmaxDOF*FEmaxDOF> m(FlexDOF(),FlexDOF());
		EvalMff(m, t);
		ms.AddMatrix(ltg,ltg,FlexDOF(),FlexDOF(),m); //ltg may be longer than FlexDOF()
	}
	else
	{
		//add dense matrix
		ms.AddMatrix(ltg,ltg,FlexDOF(),FlexDOF(),massmatrix); //ltg may be longer than FlexDOF()
	}
}

template<class BodyXD>
void FiniteElementGeneric<BodyXD>::EraseInelasticVariables()
{
	// here we assume that inelastic variables are stored in each integration point for the stiffness matrix
	int n = GetMaterial().GetInelasticVariablesCount();
	for(int nip = 1; nip <= integrationRuleStiffness->GetNumberOfIntegrationPoints(); nip++)
		for (int i = 1; i <= n; i++)
			XData((nip-1)*n+i) = 0;
}

template<>
Vector3D FiniteElementGeneric<Body3D>::GetAngularMomentum(const Vector3D& p_ref) const
{
	Vector3D D(0.,0.,0.);
	double rho = GetMaterial().Density();

	/*Vector3D v = GetVel(Vector3D(1.,1.,1.));
	Vector3D p = GetPos(Vector3D(1.,1.,1.));

	Vector xgp(NNodes()*Dim()), xg(NNodes()*Dim());*/

	/*for (int i = 1; i <= NNodes()*Dim(); i++)
	{
		xg(i) = XG(i);
		xgp(i) = XGP(i);
	}*/

	/*mbs->UO() << xg << "\n";
	mbs->UO() << xgp << "\n";
	mbs->UO() << v << "\n";
	mbs->UO() << p << "\n\n";*/

	for (IntegrationPointsIterator ip(integrationRuleMass); !ip.IsEnd(); ++ip)
	{
		Vector3D p = ip.Point();
		double det = this->GetPrecomputedJacDet(ip.GetIndex());

		Vector3D r = GetPos(p) - p_ref;
		Vector3D v = GetVel(p);
		D += r.Cross(v) * rho * ip.Weight() * det;
	}
	return D;
}

template<>
Vector3D FiniteElementGeneric<Body2D>::GetAngularMomentum(const Vector3D& p_ref) const
{
	Vector3D D(0.,0.,0.);
	// AH: to be implemented
	assert(0);
	return D;
}