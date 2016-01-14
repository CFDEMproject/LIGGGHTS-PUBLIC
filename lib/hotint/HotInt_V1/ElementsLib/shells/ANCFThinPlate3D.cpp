//#**************************************************************
//#
//# filename:             ANCFThinPlate3D.cpp
//#
//# author:               Gerstmayr Johannes, Vetyukov Yury
//#
//# generated:						first: 24. October 2006, rewritten: January-February 2011
//# description:          Kirchhoff Plate model, large deformation, 36 DOF
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
#include "material.h"
#include "ANCFThinPlate3D.h"
#include "Node.h"
//#include "graphicsconstants.h"
//#include "elementdata.h"
//#include "elementdataaccess.h"
#include "solversettings_auto.h"

PlasticVariables & PlasticVariables::operator*=(double k)
{
	strain *= k;
	hardening_parameter *= k;
	yield_function *= k;
	return *this;
}

PlasticVariables & PlasticVariables::operator+=(const PlasticVariables & T)
{
	strain += T.strain;
	hardening_parameter += T.hardening_parameter;
	yield_function += T.yield_function;
	return *this;
}

PlasticVariables & PlasticVariables::operator-=(const PlasticVariables & T)
{
	strain -= T.strain;
	hardening_parameter -= T.hardening_parameter;
	yield_function -= T.yield_function;
	return *this;
}


// this helper class will simplify browsing through all integration points for stiffness
// with the direct access to the associated data units
class IntPointsStiffnessIterator : public IntegrationPointsIterator
{
	ANCFThinPlate3D * fe;
public:
	IntPointsStiffnessIterator(ANCFThinPlate3D * fe_) :
			IntegrationPointsIterator(fe_->integrationRuleStiffness),
			fe(fe_)
			{
			}
	ANCFThinPlate3D::IntegrationPointStiffnessMatrixLocalData & Data()
	{
		return *fe->integrationPointStiffnessMatrixLocalData(GetIndex());
	}
};

#pragma region initialization

void ANCFThinPlate3D::ANCFThinPlate3DInitialize()
{
	size1 = 1;
	size2 = 1;
	thickness_layers = 1;
	nodes.SetLen(4);
	xgReferenceState.SetLen(FlexDOF());
	xgInit.SetXGProvider(this, &xgReferenceState);
	xgCompute.SetXGProvider(this, &xgReferenceState);
	xgDraw.SetXGProvider(this, &xgReferenceState);
}

void ANCFThinPlate3D::CopyFrom(const Element& e)
{
	FiniteElementGeneric<Body3D>::CopyFrom(e);
	const ANCFThinPlate3D & ce = (const ANCFThinPlate3D&)e;
	size1 = ce.size1;
	size2 = ce.size2;
	nodes = ce.nodes;
	thickness_layers = ce.thickness_layers;
	plateCrossSection = ce.plateCrossSection;
	xgReferenceState.SetLen(FlexDOF());
	xgReferenceState.Copy(ce.xgReferenceState, 1, 1, FlexDOF());
}

void ANCFThinPlate3D::Initialize() 
{
	Body3D::Initialize();
	ComputeMass();
}

// sets the characteristic sizes of the finite element in the reference configuration
// in the directions of both coordinate lines;
// the sizes are set to 2 by default
void ANCFThinPlate3D::SetSizes(double size1, double size2)
{
	this->size1 = size1;
	this->size2 = size2;
	//$ YV 2011-08: BuildDSMatrices() is now called in PreAssemble() anyway
	//BuildDSMatrices();
}

// sets the number of thickness layers ( needs to be >= 1, by default = 4 )
void ANCFThinPlate3D::SetThicknessLayers(int layers)
{
	assert(layers >= 1);
	thickness_layers = layers;
}

// this version of the set function requires that the nodes are already added to mbs;
// the nodes must have 12 degrees of freedom and be initialized according to the reference configuration
void ANCFThinPlate3D::SetANCFThinPlate3D(
	int bodyindi,									// identifier of the body
	const TArray<int>& node_list,	// list of the four indices of the nodes of the element in the global array
	int material_num,							// index of the material of the element in the global array
	double thickness,							// thickness of the plate
	const Vector3D& coli					// color of the graphical representation
	)
{
	plateCrossSection.thickness = thickness;
	if(node_list.Length() != 4)
	{
		GetMBS()->UO() << "Warning: ANCFThinPlate3D, number of nodes is invalid.\n";
		assert(0);
		return;
	}

	// we need to initialize xgReferenceState before calling to SetFiniteElementGeneric(..),
	// becase building the matrices for the reference state requires the functionality of xgInit
	nodes = node_list;
	xgReferenceState.SetLen(FlexDOF());
	for (int i = 1; i <= NNodes(); i++)
	{
		Node& node = GetNode(i);
		xgReferenceState.Copy(node.X_Init(), 1, (i-1)*DOFPerNode()+1, DOFPerNode());
	}

	FiniteElementGeneric<Body3D>::SetFiniteElementGeneric(bodyindi,node_list,material_num,coli);

	// we redefine the default behavior: as the reference state is stored in nodes,
	// and xg are displacements, we set the first half of x_init of the element (positions) to 0,
	// and keep the velocities as they are;
	// the original reference values are already stored in xgReferenceState
	for (int i = 1; i <= FlexDOF(); i++)
		x_init(i) = 0;
}

// this version of the set function creates and adds nodes to mbs;
// each node is defined through its 9 degrees of freedom in the reference configuration
void ANCFThinPlate3D::SetANCFThinPlate3D(
	int bodyindi,									// identifier of the body
	const TArray<Vector> nodes,		// degrees of freedom of the four nodes
	int material_num,							// index of the material of the element in the global array
	double thickness,							// thickness of the plate
	const Vector3D& coli					// color of the graphical representation
	)
{
	TArray<int> node_list;
	for (int i=1; i <= nodes.Length(); i++)
	{
		Vector3D pos(nodes(i)(1), nodes(i)(2), nodes(i)(3));
		Node n(9, bodyindi, pos);
		node_list.Add(AddBodyNode(n)); // adds the node either to the reference frame or to the bode and checks if node exists
	}
	SetANCFThinPlate3D(bodyindi, node_list, material_num, thickness, coli);
}

#pragma endregion

#pragma region shape functions and their derivatives

// these shape functions for the elements with 36 degrees of freedom do not provide
// full C1 continuity over the element boundaries (only at nodes)

// degrees of freedom are spatial components of the nodal values of the position vector,
// its derivative with respect to q1 and then its derivative with respect to q2 (i.e. 9 nodal dof);
// the nodes are numbered as follows:
//          q2
//          ^
//   +------+------+
//   |4     |     3|
//   |      |      |
//   +------+------+-> q1
//   |      |      |
//   |1     |     2|
//   +------+------+

// shape functions as they are, accessible per index nShapeFunction = 1 .. 12
// local coordinates q1, q2 vary from -1 to +1
double ANCFThinPlate3D::GetS0(const Vector2D& ploc, int nShapeFunction) const
{
	double q1 = ploc.X();
	double q2 = ploc.Y();

	switch(nShapeFunction)
	{
	case 1: return q1*q2/2.0-q1*q1*q1*q2/8.0-q1*q2*q2*q2/8.0+1.0/4.0+q1*q1*q1/8.0+q2*q2
		*q2/8.0-3.0/8.0*q2-3.0/8.0*q1;
	case 2: return size1*q1*q1*q2/16.0+size1/16.0+size1*q1*q2/16.0-size1*q1/16.0-size1*q1*q1/16.0-size1*
		q1*q1*q1*q2/16.0+size1*q1*q1*q1/16.0-size1*q2/16.0;
	case 3: return size2/16.0-size2*q1*q2*q2*q2/16.0-size2*q2/16.0-size2*q2*q2/16.0+size2*q2*q2*q2/
		16.0-size2*q1/16.0+size2*q1*q2/16.0+size2*q1*q2*q2/16.0;
	case 4: return -q1*q2/2.0+q1*q1*q1*q2/8.0+q1*q2*q2*q2/8.0+1.0/4.0-q1*q1*q1/8.0+q2*
		q2*q2/8.0-3.0/8.0*q2+3.0/8.0*q1;
	case 5: return -size1/16.0-size1*q1/16.0+size1*q1*q1/16.0+size1*q1*q1*q1/16.0+size1*q2/16.0+size1*q1*
		q2/16.0-size1*q1*q1*q2/16.0-size1*q1*q1*q1*q2/16.0;
	case 6: return size2/16.0-size2*q2/16.0+size2*q1/16.0-size2*q1*q2/16.0-size2*q2*q2/16.0-size2*q1*q2*
		q2/16.0+size2*q2*q2*q2/16.0+size2*q1*q2*q2*q2/16.0;
	case 7: return q1*q2/2.0-q1*q1*q1*q2/8.0-q1*q2*q2*q2/8.0+1.0/4.0-q1*q1*q1/8.0-q2*q2
		*q2/8.0+3.0/8.0*q2+3.0/8.0*q1;
	case 8: return -size1/16.0-size1*q2/16.0-size1*q1/16.0-size1*q1*q2/16.0+size1*q1*q1/16.0+size1*q1*q1*
		q2/16.0+size1*q1*q1*q1/16.0+size1*q1*q1*q1*q2/16.0;
	case 9: return -size2/16.0-size2*q2/16.0+size2*q2*q2/16.0-size2*q1/16.0-size2*q1*q2/16.0+size2*q1*q2*
		q2/16.0+size2*q2*q2*q2/16.0+size2*q1*q2*q2*q2/16.0;
	case 10: return -q1*q2/2.0+q1*q1*q1*q2/8.0+q1*q2*q2*q2/8.0+1.0/4.0+q1*q1*q1/8.0-q2*
		q2*q2/8.0+3.0/8.0*q2-3.0/8.0*q1;
	case 11: return size1/16.0+size1*q2/16.0-size1*q1/16.0-size1*q1*q2/16.0-size1*q1*q1/16.0-size1*q1*q1*
		q2/16.0+size1*q1*q1*q1/16.0+size1*q1*q1*q1*q2/16.0;
	case 12: return -size2/16.0-size2*q2/16.0+size2*q2*q2/16.0+size2*q2*q2*q2/16.0+size2*q1/16.0+size2*q1*
			q2/16.0-size2*q1*q2*q2/16.0-size2*q1*q2*q2*q2/16.0;
	}
	assert(0);
	return 0;
}

// derivatives of the shape functions with respect to the local coordinates; alpha = 1,2
double ANCFThinPlate3D::GetDS0(const Vector2D& ploc, int nShapeFunction, int alpha) const
{
	double q1 = ploc.X();
	double q2 = ploc.Y();

	if(alpha == 1)
	{
		// derivatives with respect to q1
		switch(nShapeFunction)
		{
		case 1: return q2/2.0-3.0/8.0*q1*q1*q2-q2*q2*q2/8.0+3.0/8.0*q1*q1-3.0/8.0;
		case 2: return size1*q1*q2/8.0+size1*q2/16.0-size1/16.0-size1*q1/8.0-3.0/16.0*size1*q1*q1*q2+3.0/16.0*size1*q1*q1;
		case 3: return -size2*q2*q2*q2/16.0-size2/16.0+size2*q2/16.0+size2*q2*q2/16.0;
		case 4: return -q2/2.0+3.0/8.0*q1*q1*q2+q2*q2*q2/8.0-3.0/8.0*q1*q1+3.0/8.0;
		case 5: return -size1/16.0+size1*q1/8.0+3.0/16.0*size1*q1*q1+size1*q2/16.0-size1*q1*q2/8.0-3.0/16.0*size1*q1*q1*q2;
		case 6: return size2/16.0-size2*q2/16.0-size2*q2*q2/16.0+size2*q2*q2*q2/16.0;
		case 7: return q2/2.0-3.0/8.0*q1*q1*q2-q2*q2*q2/8.0-3.0/8.0*q1*q1+3.0/8.0;
		case 8: return -size1/16.0-size1*q2/16.0+size1*q1/8.0+size1*q1*q2/8.0+3.0/16.0*size1*q1*q1+3.0/16.0*size1*q1*q1*q2;
		case 9: return -size2/16.0-size2*q2/16.0+size2*q2*q2/16.0+size2*q2*q2*q2/16.0;
		case 10: return -q2/2.0+3.0/8.0*q1*q1*q2+q2*q2*q2/8.0+3.0/8.0*q1*q1-3.0/8.0;
		case 11: return -size1/16.0-size1*q2/16.0-size1*q1/8.0-size1*q1*q2/8.0+3.0/16.0*size1*q1*q1+3.0/16.0*size1*q1*q1*q2;
		case 12: return size2/16.0+size2*q2/16.0-size2*q2*q2/16.0-size2*q2*q2*q2/16.0;
		}
	}
	else if(alpha == 2)
	{
		// derivatives with respect to q2
		switch(nShapeFunction)
		{
		case 1: return q1/2.0-q1*q1*q1/8.0-3.0/8.0*q1*q2*q2+3.0/8.0*q2*q2-3.0/8.0;
		case 2: return size1*q1*q1/16.0+size1*q1/16.0-size1*q1*q1*q1/16.0-size1/16.0;
		case 3: return -3.0/16.0*size2*q1*q2*q2-size2/16.0-size2*q2/8.0+3.0/16.0*size2*q2*q2+size2*q1/16.0+size2*q1*q2/8.0;
		case 4: return -q1/2.0+q1*q1*q1/8.0+3.0/8.0*q1*q2*q2+3.0/8.0*q2*q2-3.0/8.0;
		case 5: return size1/16.0+size1*q1/16.0-size1*q1*q1/16.0-size1*q1*q1*q1/16.0;
		case 6: return -size2/16.0-size2*q1/16.0-size2*q2/8.0-size2*q1*q2/8.0+3.0/16.0*size2*q2*q2+3.0/16.0*size2*q1*q2*q2;
		case 7: return q1/2.0-q1*q1*q1/8.0-3.0/8.0*q1*q2*q2-3.0/8.0*q2*q2+3.0/8.0;
		case 8: return -size1/16.0-size1*q1/16.0+size1*q1*q1/16.0+size1*q1*q1*q1/16.0;
		case 9: return -size2/16.0+size2*q2/8.0-size2*q1/16.0+size2*q1*q2/8.0+3.0/16.0*size2*q2*q2+3.0/16.0*size2*q1*q2*q2;
		case 10: return -q1/2.0+q1*q1*q1/8.0+3.0/8.0*q1*q2*q2-3.0/8.0*q2*q2+3.0/8.0;
		case 11: return size1/16.0-size1*q1/16.0-size1*q1*q1/16.0+size1*q1*q1*q1/16.0;
		case 12: return -size2/16.0+size2*q2/8.0+3.0/16.0*size2*q2*q2+size2*q1/16.0-size2*q1*q2/8.0-3.0/16.0*size2*q1*q2*q2;
		}
	}
	assert(0);
	return 0;
}

// second order derivatives of the shape functions with respect to the local coordinates
// alphaBeta = 1 (11), 2 (12), 3 (22)
double ANCFThinPlate3D::GetDDS0(const Vector2D& ploc, int nShapeFunction, int alphaBeta) const
{
	double q1 = ploc.X();
	double q2 = ploc.Y();

	if(alphaBeta == 1)
	{
	// derivatives d^2 S / d q1^2
		switch(nShapeFunction)
		{
		case 1: return -3.0/4.0*q1*q2+3.0/4.0*q1;
		case 2: return -3.0/8.0*size1*q1*q2+size1*q2/8.0+3.0/8.0*size1*q1-size1/8.0;
		case 3: return 0.0;
		case 4: return 3.0/4.0*q1*q2-3.0/4.0*q1;
		case 5: return size1/8.0+3.0/8.0*size1*q1-size1*q2/8.0-3.0/8.0*size1*q1*q2;
		case 6: return 0.0;
		case 7: return -3.0/4.0*q1*q2-3.0/4.0*q1;
		case 8: return size1/8.0+size1*q2/8.0+3.0/8.0*size1*q1+3.0/8.0*size1*q1*q2;
		case 9: return 0.0;
		case 10: return 3.0/4.0*q1*q2+3.0/4.0*q1;
		case 11: return -size1/8.0-size1*q2/8.0+3.0/8.0*size1*q1+3.0/8.0*size1*q1*q2;
		case 12: return 0.0;
		}
	}
	else if(alphaBeta == 2)
	{
	// derivatives d^2 S / (d q1 d q2)
		switch(nShapeFunction)
		{
		case 1: return 1.0/2.0-3.0/8.0*q1*q1-3.0/8.0*q2*q2;
		case 2: return -3.0/16.0*size1*q1*q1+size1*q1/8.0+size1/16.0;
		case 3: return size2/16.0-3.0/16.0*size2*q2*q2+size2*q2/8.0;
		case 4: return -1.0/2.0+3.0/8.0*q1*q1+3.0/8.0*q2*q2;
		case 5: return size1/16.0-size1*q1/8.0-3.0/16.0*size1*q1*q1;
		case 6: return -size2/16.0-size2*q2/8.0+3.0/16.0*size2*q2*q2;
		case 7: return 1.0/2.0-3.0/8.0*q1*q1-3.0/8.0*q2*q2;
		case 8: return -size1/16.0+size1*q1/8.0+3.0/16.0*size1*q1*q1;
		case 9: return -size2/16.0+size2*q2/8.0+3.0/16.0*size2*q2*q2;
		case 10: return -1.0/2.0+3.0/8.0*q1*q1+3.0/8.0*q2*q2;
		case 11: return -size1/16.0-size1*q1/8.0+3.0/16.0*size1*q1*q1;
		case 12: return size2/16.0-size2*q2/8.0-3.0/16.0*size2*q2*q2;
		}
	}
	else if(alphaBeta == 3)
	{
	// derivatives d^2 S / d q2^2
		switch(nShapeFunction)
		{
		case 1: return -3.0/4.0*q1*q2+3.0/4.0*q2;
		case 2: return 0.0;
		case 3: return -3.0/8.0*size2*q1*q2+size2*q1/8.0-size2/8.0+3.0/8.0*size2*q2;
		case 4: return 3.0/4.0*q1*q2+3.0/4.0*q2;
		case 5: return 0.0;
		case 6: return -size2/8.0-size2*q1/8.0+3.0/8.0*size2*q2+3.0/8.0*size2*q1*q2;
		case 7: return -3.0/4.0*q1*q2-3.0/4.0*q2;
		case 8: return 0.0;
		case 9: return size2/8.0+size2*q1/8.0+3.0/8.0*size2*q2+3.0/8.0*size2*q1*q2;
		case 10: return 3.0/4.0*q1*q2-3.0/4.0*q2;
		case 11: return 0.0;
		case 12: return size2/8.0+3.0/8.0*size2*q2-size2*q1/8.0-3.0/8.0*size2*q1*q2;
		}
	}
	assert(0);
	return 0;
}

#pragma endregion

#pragma region integration rules

void ANCFThinPlate3D::DefineIntegrationRule(IntegrationRule & integrationRule)
{
	assert(integrationRule.settings.elementType == TFE_ThinPlate);

	int ruleOrder = 0;
	// here the particular rule order will be chosen depending on the settings
	if(integrationRule.settings.integratedValueType == IntegrationRule::IVT_Stiffness)
	{
		if(integrationRule.settings.geometricNonlinearityStatus == GNS_NonlinearLargeStrain)
			ruleOrder = 8;
		else
			ruleOrder = 6;
	}
	else
		ruleOrder = 4;
	assert(ruleOrder != 0);

	if(IsInelasticMaterial())
	{
		//IntegrationRule::DefineIntegrationRuleLobattoSquare(integrationRule, ruleOrder);
		IntegrationRule::DefineIntegrationRuleSquare(integrationRule, ruleOrder);
	}
	else
	{
		//IntegrationRule::DefineIntegrationRuleLobattoSquare(integrationRule, ruleOrder);
		IntegrationRule::DefineIntegrationRuleSquare(integrationRule, ruleOrder);
	}
}

#pragma endregion

#pragma region geometry of the element

// covariant basis
Vector3D ANCFThinPlate3D::GetRAlpha(const Vector2D& ploc, int alpha, const XGProvider & xg) const
{
	Vector3D RAlpha;
	for (int i = 1; i <= Dim(); i++)
		for (int j = 1; j <= NS(); j++)
			RAlpha(i) += GetDS0(ploc,j,alpha) * xg.XGcoord((j-1)*Dim()+i);
	return RAlpha;
}

Vector3D ANCFThinPlate3D::GetRAlphaBeta(const Vector2D& ploc, int alphaBeta, const XGProvider & xg) const
{
	Vector3D RAlphaBeta;
	for (int i = 1; i <= Dim(); i++)
		for (int j = 1; j <= NS(); j++)
			RAlphaBeta(i) += GetDDS0(ploc,j,alphaBeta) * xg.XGcoord((j-1)*Dim()+i);
	return RAlphaBeta;
}

Vector3D ANCFThinPlate3D::GetNormal(const Vector2D& ploc, const XGProvider & xg) const
{
	Vector3D n = GetRAlpha(ploc, 1, xg).Cross(GetRAlpha(ploc, 2, xg));
	n.Normalize();
	return n;
}

Vector3D ANCFThinPlate3D::GetDerivativeRAlphaDOF(const Vector2D& ploc, int nDOF, int alpha) const
{
	Vector3D v;
	v((nDOF - 1) % 3 + 1) = GetDS0(ploc, (nDOF - 1) / 3 + 1, alpha);
	return v;
}

Vector3D ANCFThinPlate3D::GetDerivativeRAlphaBetaDOF(const Vector2D& ploc, int nDOF, int alphaBeta) const
{
	Vector3D v;
	v((nDOF - 1) % 3 + 1) = GetDDS0(ploc, (nDOF - 1) / 3 + 1, alphaBeta);
	return v;
}

Vector3D ANCFThinPlate3D::GetDerivativeNormalDOF(const Vector2D& ploc, int nDOF, const XGProvider & xg) const
{
	Vector3D r1 = GetRAlpha(ploc, 1, xg);
	Vector3D r2 = GetRAlpha(ploc, 2, xg);
	Vector3D n = r1.Cross(r2);
	double norm = n.Norm();
	n.Normalize();

	// v = dr1/dei x r2 + r1 x dr2/dei
	Vector3D v = GetDerivativeRAlphaDOF(ploc, nDOF, 1).Cross(r2) + r1.Cross(GetDerivativeRAlphaDOF(ploc, nDOF, 2));
	// projection of v onto the tangent plane
	v -= n * (v * n);

	return v * (1 / norm);
}

Vector3D ANCFThinPlate3D::GetPosD(const Vector3D& ploc, int use_magnification) const
{
	Vector3D p0 = GetRefConfPos(ploc);
	double factor = GetMBS()->GetDOption(105);
	if (!use_magnification)
		factor = 1;
	Vector3D u = GetDisplacement(ploc, xgDraw);
	return p0 + factor * u;
}

Vector3D ANCFThinPlate3D::GetPos(const Vector3D& ploc, const XGProvider & xg) const
{
	Vector3D p;
	for (int i = 1; i <= Dim(); i++)
		for (int j = 1; j <= NS(); j++)
			p(i) += GetS0((const Vector2D &)ploc, j) * xg.XGcoord((j-1)*Dim()+i);

	if(ploc.Z() == 0)
		return p;		// time saving

	double thickness = GetThicknessAtPoint((const Vector2D &)ploc, xg.IsDrawConfiguration());
	return p + GetNormal((const Vector2D&)ploc, xg) * (ploc.Z() * thickness / 2);
}

Vector3D ANCFThinPlate3D::GetDisplacement(const Vector3D& ploc, const XGProvider & xg) const
{
	if(GetGeometricNonlinearityStatus() == GNS_Linear)
	{
		// linearized formulas
		Vector3D u0;	// displacement on the middle surface
		for (int i = 1; i <= Dim(); i++)
			for (int j = 1; j <= NS(); j++)
				u0(i) += GetS0((const Vector2D &)ploc, j) * xg.XGdispl((j-1)*Dim()+i);
		if(ploc.Z() == 0)
			return u0;		// time saving
		// change of the normal vector
		Vector3D dn;
		for(int i = 1; i <= DOFPerNode() * NNodes(); i++)
			dn += GetDerivativeNormalDOF((const Vector2D&)ploc, i, xgInit) * xg.XGdispl(i);
		double thickness = GetThicknessAtPoint((const Vector2D &)ploc, xg.IsDrawConfiguration());
		return u0 + dn * (ploc.Z() * thickness / 2);
	}
	else
		return GetPos(ploc,xg) - GetRefConfPos(ploc);
}

Vector3D ANCFThinPlate3D::GetVel(const Vector3D& ploc, const XGProvider & xg) const
{
	Vector3D v0;	// velocity on the middle surface
	for (int i = 1; i <= Dim(); i++)
		for (int j = 1; j <= NS(); j++)
			v0(i) += GetS0((const Vector2D &)ploc, j) * xg.XGP((j-1)*Dim()+i);
	if(ploc.Z() == 0)
		return v0;		// time saving
	// time derivative of the normal vector
	Vector3D dndt;
	for(int i = 1; i <= DOFPerNode() * NNodes(); i++)
		dndt += GetDerivativeNormalDOF((const Vector2D&)ploc, i, xg) * xg.XGP(i);
	double thickness = GetThicknessAtPoint((const Vector2D &)ploc, xg.IsDrawConfiguration());
	return v0 + dndt * (ploc.Z() * thickness / 2);
}

Vector3D ANCFThinPlate3D::GetNodeLocPos(int local_node_number) const
{
	switch(local_node_number)
	{
	case 1: return Vector3D(-1,-1,0);
	case 2: return Vector3D(1,-1,0);
	case 3: return Vector3D(1,1,0);
	case 4: return Vector3D(-1,1,0);
	}
	assert(0);
	return 0;
}

Vector3D ANCFThinPlate3D::GetDOFPosD(int idof) const
{
	int nodeStartDOF = ((idof-1)/9) * 9 + 1;
	return Vector3D(
		xgDraw.XGcoord(nodeStartDOF + 0),
		xgDraw.XGcoord(nodeStartDOF + 1),
		xgDraw.XGcoord(nodeStartDOF + 2)
		);
}

Vector3D ANCFThinPlate3D::GetDOFDirD(int idof) const
{
	int dir = (idof-1)%3;
	if (dir == 0) return Vector3D(1.,0.,0.);
	else if (dir == 1) return Vector3D(0.,1.,0.);
	else return Vector3D(0.,0.,1.);
}

#pragma endregion

#pragma region metric tensors and strain tensors

// compute the components of the first metric tensor
PlaneSymmetricTensorComponents ANCFThinPlate3D::GetMetricTensorAComponents(const Vector2D& ploc, const XGProvider & xg) const
{
	Vector3D r1 = GetRAlpha(ploc, 1, xg);
	Vector3D r2 = GetRAlpha(ploc, 2, xg);
	return PlaneSymmetricTensorComponents( r1 * r1, r1 * r2, r2 * r2 );
}

// compute the components of the second metric tensor
PlaneSymmetricTensorComponents ANCFThinPlate3D::GetMetricTensorBComponents(const Vector2D& ploc, const XGProvider & xg) const
{
	Vector3D n = GetNormal(ploc, xg);
	PlaneSymmetricTensorComponents B;
	for(int alphaBeta = 1; alphaBeta <= 3; alphaBeta++)
		B.SetComponent(alphaBeta, GetRAlphaBeta(ploc, alphaBeta, xg) * n);
	return B;
}

// compute the components of the first strain tensor (in-plane deformations)
PlaneSymmetricTensorComponents ANCFThinPlate3D::GetFirstStrainTensorCComponents(
												const Vector2D & ploc, const Matrix & grad, const XGProvider & xg) const
{
	// CAlphaBeta = (AAlphaBeta - aAlphaBeta) / 2
	PlaneSymmetricTensorComponents C = GetMetricTensorAComponents(ploc, xg);
	PlaneSymmetricTensorComponents a;
	// we read the components of the first metric tensor in the reference configuration from the saved data
	for(int alphaBeta = 1; alphaBeta <= 3; alphaBeta++)
		a.SetComponent(alphaBeta, grad(2, alphaBeta));
	C -= a;
	C *= 0.5;
	return C;
}

// compute the components of the second strain tensor (bending - change of curvature)
PlaneSymmetricTensorComponents ANCFThinPlate3D::GetSecondStrainTensorKComponents(
												const Vector2D & ploc, const Matrix & grad, const XGProvider & xg) const
{
	// KAlphaBeta = BAlphaBeta - bAlphaBeta
	PlaneSymmetricTensorComponents K = GetMetricTensorBComponents(ploc, xg);
	PlaneSymmetricTensorComponents b;
	// we read the components of the second metric tensor in the reference configuration from the saved data
	for(int alphaBeta = 1; alphaBeta <= 3; alphaBeta++)
		b.SetComponent(alphaBeta, grad(3, alphaBeta));
	K -= b;
	return K;
}

#pragma endregion

#pragma region computation

// this is the central point of the class: the differential geometry of the surface in the Gauss point,
// which is used later on for computing the parameters of the actually deformed surface.
// what the function returns is the jacobian (integration coefficient), in our case
// it is the determinant of the matrix of components of the first metric tensor;
// the data, which need to be pre-computed in each Gauss point, is stored in jacinvDS.
// the argument is Vector3D, as the function is called from the Generic class
double ANCFThinPlate3D::GetJacInvDS(const Vector3D& ploc, Matrix& jacinvDS) const
{
	// in each integration point we pre-compute the following data:
	// 1. inverse of the components of the first metric tensor a^{\alpha\beta} (3 values - diagonal and non-diaginal) (contravariant)
	// 2. components of the first metric tensor a_{\alpha\beta} (covariant)
	// 3. components of the second metric tensor b_{\alpha\beta} (covariant)
	PlaneSymmetricTensorComponents a = GetMetricTensorAComponents((const Vector2D&)ploc, xgInit);
	PlaneSymmetricTensorComponents b = GetMetricTensorBComponents((const Vector2D&)ploc, xgInit);
	PlaneSymmetricTensorComponents aContravariant = a.Inverse();
	jacinvDS = Matrix(		aContravariant.t11,	aContravariant.t12,	aContravariant.t22,
												a.t11,							a.t12,							a.t22,
												b.t11,							b.t12,							b.t22				);
	return sqrt(a.Det());
}

double ANCFThinPlate3D::GetJacobianDeterminant(const Vector2D& ploc) const
{
	return sqrt(GetMetricTensorAComponents(ploc, xgInit).Det());
}

void ANCFThinPlate3D::ComputeMass()
{
	mass = 0;
	for (IntegrationPointsIterator ip(integrationRuleMass); !ip.IsEnd(); ++ip)
	{
		double thickness = GetThicknessAtPoint((const Vector2D &)ip.Point2D(), false);
		mass += GetJacobianDeterminant(ip.Point2D()) * ip.Weight() * thickness * Rho();
	}
}

void ANCFThinPlate3D::GetH(Matrix& H) 
{
	int ns = NS();
	int sos = SOS();
	int dim = Dim();
	if (Hmatrix.Getrows() == ns*dim)
	{
		H = Hmatrix;
		return;
	}
	else
	{
		H.SetSize(ns * dim, dim);
		H.FillWithZeros();
		ConstMatrix<FEmaxDOF*3> dudq;
		for (IntegrationPointsIterator ip(integrationRuleLoad); !ip.IsEnd(); ++ip)
		{
			double thickness = GetThicknessAtPoint((const Vector2D &)ip.Point2D(), false);
			double fact = GetJacobianDeterminant(ip.Point2D()) * ip.Weight() * thickness;
			for (int i = 0; i < ns; i++)
				for (int j = 1; j <= dim; j++)
					H(i*dim+j,j) += fact * GetS0(ip.Point2D(), i+1);
		}
		Hmatrix = H;
	}
}

// here we ignore the "moment" effect of a force, shifted away from the mid-surface
void ANCFThinPlate3D::GetdPosdqT(const Vector3D& ploc, Matrix& d)
{
	d.SetSize(FlexDOF(),Dim());
	d.FillWithZeros();

	for (int i = 1; i <= NS(); i++)
	{
		double sv =	GetS0((const Vector2D&)ploc, i);
		d((i-1)*Dim()+1,1) = sv;
		d((i-1)*Dim()+2,2) = sv;
		d((i-1)*Dim()+3,3) = sv;
	}
}

// fill only FlexibleDof entries of massmatrix, other entries remain unchanged
void ANCFThinPlate3D::EvalMff(Matrix& m, double t)
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
		for (IntegrationPointsIterator ip(integrationRuleMass); !ip.IsEnd(); ++ip)
		{
			// compute SF = vector containing shape functions
			for (int i = 1; i <= ns; i++)
				SF(i) = GetS0(ip.Point2D(), i);

			// compute S * S^T  where S is the shape function matrix 
			//     [ SF(1) 0 0 SF(2) 0 0 ...
			// S = [ 0 SF(1) 0 0 SF(2) 0 ..
			//     [ 0 0 SF(1) 0 0 SF(2) ..
			SST.SetAll(0);
			for (int i = 0; i < ns; i++)
				for (int j = 0; j < ns; j++)
					for (int di = 1; di <= dim; di++)
							SST(i*dim+di,j*dim+di) += SF(i+1)*SF(j+1);

			double thickness = GetThicknessAtPoint((const Vector2D &)ip.Point2D(), false);
			m.AddSubmatrix(SST,1,1,1,1,flexdof,flexdof,
				GetJacobianDeterminant(ip.Point2D()) * ip.Weight() * thickness * Rho());
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

void ANCFThinPlate3D::EvalF2(Vector& f, double t) 
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

	if (GetGeometricNonlinearityStatus() == GNS_Linear)
	{
		EvalF2GeomLin(fadd, t);
	}
	else
	{
		EvalF2GeomNonlin(fadd, t);
	}

	f -= fadd;
}; 


void ANCFThinPlate3D::EvalF2GeomNonlin(Vector& fadd, double t) 
{
	//add loads, constraint forces and other things in parent function:
	int sos = SOS();
	int flexdof = FlexDOF();
	int ns = NS();
	int dim = Dim();

	XGProviderCoordCached<FEmaxDOF> xg(xgCompute, XGLength());

	// compute elastic forces
	for (IntPointsStiffnessIterator ip(this); !ip.IsEnd(); ++ip)
	{
		// metric of the reference configuration in this point
		Matrix grad;
		GetGrad(ip, grad);
		PlaneSymmetricTensorComponents aInv( grad(1,1), grad(1,2), grad(1,3) );
		const Vector2D & ploc = ip.Point2D();

		// actual strains in this point
		PlaneSymmetricTensorComponents C = GetFirstStrainTensorCComponents(ploc, grad, xg);
		PlaneSymmetricTensorComponents K = GetSecondStrainTensorKComponents(ploc, grad, xg);

		// initial strains are simply subtracted from the first strain tensor
		// setting the initial strains can yet only be done in the derived classes
		if(HasInitialStrains())
		{
			Matrix2D initial_strain_cartesian;
			GetInitialStrainComponentsAtPoint(initial_strain_cartesian, ip.Point2D(), false);
			PlaneSymmetricTensorComponents initial_strain = CartesianComponentsToLocal(ip, initial_strain_cartesian);
			C -= initial_strain;
		}

		Vector3D r1 = GetRAlpha(ploc, 1, xg);
		Vector3D r2 = GetRAlpha(ploc, 2, xg);
		Vector3D rAlphaBeta[4];
		for(int alphaBeta = 1; alphaBeta <= 3; alphaBeta++)
		{
			rAlphaBeta[alphaBeta] = GetRAlphaBeta(ploc, alphaBeta, xg);
		}
		Vector3D n = GetNormal(ploc, xg);

		plateCrossSection.thickness = GetThicknessAtPoint(ploc, false);
		plateCrossSection.ComputeStiffnesses(Em(), Nu());

#if 0
		// TODO
		if (GetMaterial().IsInelasticMaterial())
		{
			PlaneSymmetricTensorComponents a = GetMetricTensorAComponents(ip.Point2D(), xg);
			PlaneSymmetricTensorComponents b = GetMetricTensorBComponents(ip.Point2D(), xg);
			PlaneSymmetricTensorComponents Ep(0.,0.,0.);

			if (thickness_layers == 1)
			{
				Ep += GetPlasticVariables(ip, 1, xg, aInv, a, b).strain;
			}
			else
			{
				Vector w(thickness_layers, 0, 1.);   //weights
				w(1) = 0.5;
				w(thickness_layers) = 0.5;

				for (int layer = 1; layer <= thickness_layers; layer++)
				{
					// integration by trapeciod rule (uniform distribution of integration points along thickness)
					PlaneSymmetricTensorComponents Ep_at_layer = GetPlasticVariables(ip, layer, xg, aInv, a, b).strain;
					Ep_at_layer *= w(layer);
					Ep += Ep_at_layer;
				}
				Ep *= (thickness_layers - 1) / plateCrossSection.thickness;
			}

			C -= Ep;
			// todo:
			// K -= z* Ep;
		}
#endif

		// a loop over the local degrees of freedom e
		for(int i = 1; i <= FlexDOF(); i++)
		{
			// derivatives of the strain measure components with respect to the current degree of freedom
			// CAlphaBeta = (AAlphaBeta - aAlphaBeta) / 2
			Vector3D dr1de = GetDerivativeRAlphaDOF(ploc, i, 1);
			Vector3D dr2de = GetDerivativeRAlphaDOF(ploc, i, 2);
			PlaneSymmetricTensorComponents dCde(r1 * dr1de, (r1 * dr2de + dr1de * r2) / 2, r2 * dr2de);
			// KAlphaBeta = BAlphaBeta - bAlphaBeta
			Vector3D dnde = GetDerivativeNormalDOF(ploc, i, xg);
			PlaneSymmetricTensorComponents dKde;
			for(int alphaBeta = 1; alphaBeta <= 3; alphaBeta++)
			{
				dKde.SetComponent(
									alphaBeta,
									rAlphaBeta[alphaBeta] * dnde + GetDerivativeRAlphaBetaDOF(ploc, i, alphaBeta) * n
									);
			}

			// derivative of the strain energy per unit area in the reference configuration
			double dUde = C.Convolute(dCde, plateCrossSection.E1, plateCrossSection.E2, aInv) +
										K.Convolute(dKde, plateCrossSection.D1, plateCrossSection.D2, aInv);
			// contribution to the force vector
			fadd(i) += dUde * ip.Data().jacdet * ip.Weight();
		}
	}
}

double ANCFThinPlate3D::PostNewtonStep(double t)
{

	if (!IsInelasticMaterial()) return 0;

	//----------------------
	// Variable declaration and Parameters
	double err = 0; // error for nonlinear iteration
	//XGProviderCoordCached<FEmaxDOF> xg(xgInit, XGLength());
	XGProviderCoordCached<FEmaxDOF> xg(xgCompute, XGLength());

	//----------------------
	// loop over integration points
#if 0
	// TODO
	for (IntPointsStiffnessIterator ip(this); !ip.IsEnd(); ++ip)
	{
		Matrix grad;
		GetGrad(ip, grad);
		PlaneSymmetricTensorComponents aInv( grad(1,1), grad(1,2), grad(1,3) );
		PlaneSymmetricTensorComponents a( grad(2,1), grad(2,2), grad(2,3) );        //A = GetMetricTensorAComponents(ip.Point2D(), xg);
		PlaneSymmetricTensorComponents b( grad(3,1), grad(3,2), grad(3,3) );        //B = GetMetricTensorBComponents(ip.Point2D(), xg);

		PlasticVariables plastic_variables;
		for (int layer = 1; layer <= thickness_layers; layer++)
		{
			//DataToPlasticVariables(plastic_variables, ip.GetIndex(), layer);
			plastic_variables = GetPlasticVariables(ip, layer, xg, aInv, a, b);
			PlasticVariablesToData(plastic_variables, ip.GetIndex(), layer);
		}
	}
#endif
	return err;
}

// this method computes the threedimensional position regarding the specific thickness layer
Vector3D ANCFThinPlate3D::GetLocalPositionAtThicknessLayer(const Vector2D & ploc, int thickness_layer, bool flagDraw) const
{
	double zpos;
	
	if (thickness_layers == 1)
	{
		zpos = 0.;
	}
	else
	{
		zpos = GetThicknessAtPoint(ploc, flagDraw) * ((thickness_layer - 1) / (thickness_layers - 1) - .5);
	}

	return Vector3D(ploc.X(), ploc.Y(), zpos);
}

#if 0

// save plastic variables to xdata (back and forth + drawmode)
void ANCFThinPlate3D::DataToPlasticVariables(PlasticVariables & plastic_variables, int ip_idx, int layer_idx)
{	
	int n = 5;
	int block = ((ip_idx - 1) * thickness_layers + layer_idx - 1) * n;

	plastic_variables.strain.SetComponent(1, XData(block + 1));
	plastic_variables.strain.SetComponent(2, XData(block + 2));
	plastic_variables.strain.SetComponent(3, XData(block + 3));
	plastic_variables.hardening_parameter = XData(block + 4);
	plastic_variables.yield_function = XData(block + 5);
}
void ANCFThinPlate3D::DataToPlasticVariablesD(PlasticVariables & plastic_variables, int ip_idx, int layer_idx)
{	
	int n = 5;
	int block = ((ip_idx - 1) * thickness_layers + layer_idx - 1) * n;

	plastic_variables.strain.SetComponent(1, XDataD(block + 1));
	plastic_variables.strain.SetComponent(2, XDataD(block + 2));
	plastic_variables.strain.SetComponent(3, XDataD(block + 3));
	plastic_variables.hardening_parameter = XDataD(block + 4);
	plastic_variables.yield_function = XDataD(block + 5);
}
void ANCFThinPlate3D::PlasticVariablesToData(const PlasticVariables & plastic_variables, int ip_idx, int layer_idx)
{	
	int n = 5;
	int block = ((ip_idx - 1) * thickness_layers + layer_idx - 1) * n;

	XData(block + 1) = plastic_variables.strain.t11;
	XData(block + 2) = plastic_variables.strain.t12;
	XData(block + 3) = plastic_variables.strain.t22;
	XData(block + 4) = plastic_variables.hardening_parameter;
	XData(block + 5) = plastic_variables.yield_function;
}
void ANCFThinPlate3D::PlasticVariablesToDataD(const PlasticVariables & plastic_variables, int ip_idx, int layer_idx)
{	
	int n = 5;
	int block = ((ip_idx - 1) * thickness_layers + layer_idx - 1) * n;

	XDataD(block + 1) = plastic_variables.strain.t11;
	XDataD(block + 2) = plastic_variables.strain.t12;
	XDataD(block + 3) = plastic_variables.strain.t22;
	XDataD(block + 4) = plastic_variables.hardening_parameter;
	XDataD(block + 5) = plastic_variables.yield_function;
}

// computes local components of the 3D plastic strain tensor
PlasticVariables ANCFThinPlate3D::GetPlasticVariables(const IntegrationPointsIterator & ip, int thickness_layer,
																											const XGProvider & xg,
																											const PlaneSymmetricTensorComponents & tainv,
																											const PlaneSymmetricTensorComponents & ta,
																											const PlaneSymmetricTensorComponents & tb)
{
	PlasticVariables plastic_variables;

PlaneSymmetricTensorComponents a(1,0,1);
PlaneSymmetricTensorComponents ainv(1,0,1);
PlaneSymmetricTensorComponents b(0,0,0);
	if (1) // new version
	{
		DataToPlasticVariables(plastic_variables, ip.GetIndex(), thickness_layer);

		PlaneSymmetricTensorComponents trial_stress =	GetStressComponentsD(GetLocalPositionAtThicknessLayer(ip.Point2D(), thickness_layer, xg.IsDrawConfiguration()), xg, ainv, a, b);
		double trace = trial_stress.Trace(ainv);
		trial_stress.BuildDeviator(a, ainv);
		double norm_trial_stress = sqrt(trial_stress.SelfConvolute(0, 1, ainv) + trace*trace/9.);
		double scaled_yield_stress = sqrt(2./3.) * GetMaterial().YieldStress();
		double tangent_module = GetMaterial().TangentModule();

		double beta = norm_trial_stress - scaled_yield_stress*(1. + tangent_module * plastic_variables.hardening_parameter);
		if (beta > 0)  // new plasticity occurring
		{
			double mu = ElasticMu();
			double xi = 2. * mu + pow(scaled_yield_stress * tangent_module, 2);
			PlaneSymmetricTensorComponents update = trial_stress;
			update *= beta / xi / norm_trial_stress;
			plastic_variables.strain += update;

			double norm_update = sqrt(update.SelfConvolute(0, 1, ainv) + pow(update.Trace(ainv),2));
			plastic_variables.hardening_parameter += scaled_yield_stress * tangent_module * norm_update;

			// compute yield function due to new plastic strain
			// therefore, compute deviator of new stress, due to new elastic strain (= total strain - new plastic strain)
			PlaneSymmetricTensorComponents eps = GetStrainComponentsD(GetLocalPositionAtThicknessLayer(ip.Point2D(), thickness_layer, xg.IsDrawConfiguration()), xg, a, b, ST_total_strain);
			eps -= plastic_variables.strain;
			PlaneSymmetricTensorComponents sigma1 = eps;
			sigma1 *= 2 * mu;
			PlaneSymmetricTensorComponents sigma2 = a;
			sigma2 *= eps.Trace(ainv) * 2 * mu * ElasticCoefficient();
			sigma1 += sigma2;

			trace = sigma1.Trace(ainv);
			sigma1.BuildDeviator(a,ainv);
			norm_trial_stress = sqrt(sigma1.SelfConvolute(0, 1, ainv) + trace*trace/9.);
		}
		plastic_variables.yield_function = norm_trial_stress - scaled_yield_stress * (1 + tangent_module * plastic_variables.hardening_parameter);
		plastic_variables.yield_function *= sqrt(3./2.);
	}
	else // old version
	{
		DataToPlasticVariables(plastic_variables, ip.GetIndex(), thickness_layer);

		PlaneSymmetricTensorComponents trial_stress =	plastic_variables.strain;
		double mu = ElasticMu();
		trial_stress *= -2. * mu;
		trial_stress += GetStressComponentsD(GetLocalPositionAtThicknessLayer(ip.Point2D(), thickness_layer, xg.IsDrawConfiguration()), xg, ainv, a, b);
		trial_stress.BuildDeviator(a, ainv);
		double norm_trial_stress = sqrt(trial_stress.SelfConvolute(1, 1, ainv));
		double scaled_yield_stress = sqrt(2./3.) * GetMaterial().YieldStress();
		double tangent_module = GetMaterial().TangentModule();

		double beta = norm_trial_stress - scaled_yield_stress*(1. + tangent_module * plastic_variables.hardening_parameter);
		if (beta > 0)  // new plasticity occurring
		{
			double xi = 2. * mu + pow(scaled_yield_stress * tangent_module, 2);
			PlaneSymmetricTensorComponents update = trial_stress;
			update *= beta / xi / norm_trial_stress;
			plastic_variables.strain += update;

			double norm_update = sqrt(update.SelfConvolute(1, 1, ainv));
			plastic_variables.hardening_parameter += scaled_yield_stress * tangent_module * norm_update;

			update *= 2. * mu;
			trial_stress -= update;
			norm_trial_stress = sqrt(trial_stress.SelfConvolute(1, 1, ainv));
		}
		plastic_variables.yield_function = norm_trial_stress - scaled_yield_stress * (1 + tangent_module * plastic_variables.hardening_parameter);
		plastic_variables.yield_function *= sqrt(3./2.);
	}

	return plastic_variables;
}

PlasticVariables ANCFThinPlate3D::GetOldPlasticVariables(const IntegrationPointsIterator & ip, int thickness_layer)
{
	PlasticVariables plastic_variables;
	DataToPlasticVariables(plastic_variables, ip.GetIndex(), thickness_layer);
	return plastic_variables;
}

// return plastic variables at x,y of integration point and z (in [-1,1]) linearly interpolated between values in thickness layers
PlasticVariables ANCFThinPlate3D::GetOldPlasticVariables(const IntegrationPointsIterator& ip, double z, bool flagDraw)
{
	PlasticVariables inelastic_variables;
	double thickness = GetThicknessAtPoint(ip.Point2D(), flagDraw);
	double scaled_z_position = .5*(z + 1);   // in [0,1]
	double scaled_mesh_size_in_height = 1./(thickness_layers - 1);
	double scaled_z_position_of_layer = 0.;
	double epsilon = (1e-16)*(thickness_layers - 1);
	int layer = 1;

	while ( (layer < thickness_layers) && (scaled_z_position_of_layer <= scaled_z_position + epsilon) )
	{
		scaled_z_position_of_layer += scaled_mesh_size_in_height;
		layer++;
	}

	if (layer == 1)
	{
		// return plastic variables at lowest thickness layer
		inelastic_variables = GetOldPlasticVariables(ip, layer);
	}
	else   // layer > 1
	{
		//return a linear combination of plastic variables of upper and lower thickness layer
		double lambda = (scaled_z_position_of_layer - scaled_z_position)/(thickness_layers - 1);
		PlasticVariables lower_inelastic_variables(GetOldPlasticVariables(ip, layer-1));
		lower_inelastic_variables *= 1 - lambda;
		inelastic_variables = GetOldPlasticVariables(ip, layer);
		inelastic_variables *= lambda;
		inelastic_variables += lower_inelastic_variables;
	}

	return inelastic_variables;
}

#endif

void ANCFThinPlate3D::EvalF2GeomLin(Vector& fadd, double t) 
{
	assert(0 && "Geometrically linear shells are not yet implemented");
}

//fill in sos x sos components, m might be larger
//compute stiffness matrix
void ANCFThinPlate3D::StiffnessMatrix(Matrix& m) 
{
	if(
		mbs->GetSolSet().store_FE_matrices &&
		GetGeometricNonlinearityStatus() == GNS_Linear &&
		stiffnessmatrix.Getcols() == FlexDOF()
		)
	{
		m.SetSubmatrix(stiffnessmatrix,1,1);
		return;
	}

	m.SetAll(0);

	XGProviderCoordCached<FEmaxDOF> xg(xgCompute, XGLength());
	ConstMatrix<FEmaxDOF*FEmaxDOF> d2Ude2(FlexDOF(), FlexDOF()); 

	// Integration routine
	for (IntPointsStiffnessIterator ip(this); !ip.IsEnd(); ++ip)
	{
		// metric of the reference configuration in this point
		Matrix grad;
		GetGrad(ip, grad);
		PlaneSymmetricTensorComponents aInv( grad(1,1), grad(1,2), grad(1,3) );
		const Vector2D & ploc = ip.Point2D();

		// actual strain measures in this point
		PlaneSymmetricTensorComponents C = GetFirstStrainTensorCComponents(ploc, grad, xg);
		PlaneSymmetricTensorComponents K = GetSecondStrainTensorKComponents(ploc, grad, xg);

		// initial strains are simply subtracted from the first strain tensor
		// YV: for in-plane problem this subtraction leads somehow to worse convergence,
		// but for problems with bending it is nesessary
		if(HasInitialStrains())
		{
			Matrix2D initial_strain_cartesian;
			GetInitialStrainComponentsAtPoint(initial_strain_cartesian, ip.Point2D(), false);
			PlaneSymmetricTensorComponents initial_strain = CartesianComponentsToLocal(ip, initial_strain_cartesian);
			C -= initial_strain;
		}

		// actual kinematic characteristics in this point
		Vector3D r1 = GetRAlpha(ploc, 1, xg);
		Vector3D r2 = GetRAlpha(ploc, 2, xg);
		Vector3D rAlphaBeta[4];
		for(int alphaBeta = 1; alphaBeta <= 3; alphaBeta++)
			rAlphaBeta[alphaBeta] = GetRAlphaBeta(ploc, alphaBeta, xg);
		Vector3D v = r1.Cross(r2);		// simple vector product - unnormed normal
		double vNormInv = 1 / v.Norm();
		Vector3D n = v * vNormInv;

		plateCrossSection.thickness = GetThicknessAtPoint(ip.Point2D(), false);
		plateCrossSection.ComputeStiffnesses(Em(), Nu());

		// first we pre-compute the first order derivatives of the normal N
		// and of the strain measures with respect to the degrees of freedom of the element;
		// we don't pre-compute the derivarives of RAlpha and RAlphaBeta as they can be fast computed;
		// we reserve memory on the stack
		Vector3D dvde[FEmaxDOF + 1];
		Vector3D dnde[FEmaxDOF + 1];
		PlaneSymmetricTensorComponents dCde[FEmaxDOF + 1];
		PlaneSymmetricTensorComponents dKde[FEmaxDOF + 1];
		for(int i = 1; i <= FlexDOF(); i++)
		{
			Vector3D dr1dei = GetDerivativeRAlphaDOF(ploc, i, 1);
			Vector3D dr2dei = GetDerivativeRAlphaDOF(ploc, i, 2);
			Vector3D drAlphaBetadei[4];
			for(int alphaBeta = 1; alphaBeta <= 3; alphaBeta++)
				drAlphaBetadei[alphaBeta] = GetDerivativeRAlphaBetaDOF(ploc, i, alphaBeta);
			dCde[i].SetComponents(r1 * dr1dei, (r1 * dr2dei + dr1dei * r2) / 2, r2 * dr2dei);
			Vector3D & dvdei = dvde[i] = dr1dei.Cross(r2) + dr2dei.Cross(r1);
			Vector3D & dndei = dnde[i] = (dvdei - n * (n * dvdei)) * vNormInv;
			for(int alphaBeta = 1; alphaBeta <= 3; alphaBeta++)
				dKde[i].SetComponent(
									alphaBeta,
									rAlphaBeta[alphaBeta] * dndei + drAlphaBetadei[alphaBeta] * n
									);
			// in a nested loop we compute the second order derivatives
			// and, finally, the contributions in the components of the stiffness matrix;
			// the entries of the arrays for i and j are already computed
			for(int j = 1; j <= i; j++)
			{
				Vector3D dr1dej = GetDerivativeRAlphaDOF(ploc, j, 1);
				Vector3D dr2dej = GetDerivativeRAlphaDOF(ploc, j, 2);
				Vector3D drAlphaBetadej[4];
				for(int alphaBeta = 1; alphaBeta <= 3; alphaBeta++)
					drAlphaBetadej[alphaBeta] = GetDerivativeRAlphaBetaDOF(ploc, j, alphaBeta);
				Vector3D d2vdeidej = dr1dei.Cross(dr2dej) + dr2dei.Cross(dr1dej);
        Vector3D d2ndeidej = d2vdeidej - n * (n * d2vdeidej) -
										dndei * (n * dvde[j]) - dnde[j] * (n * dvdei) - n * (dndei * dvde[j]);
				PlaneSymmetricTensorComponents d2Cdeidej(
												dr1dei * dr1dej,
												(dr1dei * dr2dej + dr1dej * dr2dei)/2,
												dr2dei * dr2dej
												);
				PlaneSymmetricTensorComponents d2Kdeidej;
				for(int alphaBeta = 1; alphaBeta <= 3; alphaBeta++)
					d2Kdeidej.SetComponent(
												alphaBeta,
												drAlphaBetadei[alphaBeta] * dnde[j] + drAlphaBetadej[alphaBeta] * dndei + rAlphaBeta[alphaBeta] * d2ndeidej
												);
				double contribution =	C.Convolute(d2Cdeidej, plateCrossSection.E1, plateCrossSection.E2, aInv) +
															dCde[i].Convolute(dCde[j], plateCrossSection.E1, plateCrossSection.E2, aInv) +
															K.Convolute(d2Kdeidej, plateCrossSection.D1, plateCrossSection.D2, aInv) +
															dKde[i].Convolute(dKde[j], plateCrossSection.D1, plateCrossSection.D2, aInv);
				contribution *= ip.Data().jacdet * ip.Weight();
				d2Ude2(i, j) += contribution;
				if (i != j)
					d2Ude2(j, i) += contribution;
			}
		}
	}

	m.AddSubmatrix(d2Ude2, 1, 1, 1, 1, FlexDOF(), FlexDOF(), -1);

	if (mbs->GetSolSet().store_FE_matrices && GetGeometricNonlinearityStatus() == GNS_Linear) // if stiffness matrix is stored
	{
		stiffnessmatrix.SetSize(FlexDOF(), FlexDOF());
		for (int i = 1; i <= FlexDOF(); i++)
			for (int j = 1; j <= FlexDOF(); j++)
				stiffnessmatrix(i,j) = d2Ude2(i,j);
	}

	/////////////////////////////DEBUG
#if 0
	std::ofstream os("d:\\matrix_hotint.txt");
	os << "{\n";
	for (int i=1; i<=d2Ude2.Getrows(); i++)
	{
		os << "{ ";
		for (int j=1; j<=d2Ude2.Getcols(); j++)
		{
			char str[32];
			sprintf(str,"%1.16g",d2Ude2(i,j));
			os << str;
			if (j!=d2Ude2.Getcols())
				os << ", ";
		}
		os << " }";
		if (i!=d2Ude2.Getrows())
				os << ",";
		os << "\n";
	}
	os << "}";
	os.flush();
#endif
	//////////////////////////////////
}

void ANCFThinPlate3D::PlateCrossSection::ComputeStiffnesses(double YoungModulus, double PoissonRatio)
{
	E1 = YoungModulus * thickness * PoissonRatio / (1 -PoissonRatio * PoissonRatio);
	E2 = YoungModulus * thickness / (1 + PoissonRatio);
	D1 = E1 * thickness * thickness / 12;
	D2 = E2 * thickness * thickness / 12;
}

#pragma endregion

#pragma region field variables computation and drawing

void ANCFThinPlate3D::GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables)
{
	// add all the field variables of the parent class
	Body3D::GetAvailableFieldVariables(variables);
	//FVT_position,FVT_diplacement, FVT_velocity: implemented in Element

	// for testing the models it is practical to be able to visualize the local coordinates
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_acceleration,
		FieldVariableDescriptor::FVCI_z, false);

	// all entities are defined in the local coordinate system on the surface,
	// such that by now we let them be post-processed in some Cartesian frame in the tangent surface;
	// a recomputation into a fixed spatial three-dimensional frame may be implemented later
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_stress,
		FieldVariableDescriptor::FVCI_y, true);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_stress_mises);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_total_strain,
		FieldVariableDescriptor::FVCI_y, true);
	
	if(GetMaterial().IsInelasticMaterial())
	{
		FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_inelastic_strain,
			FieldVariableDescriptor::FVCI_y, true);
		if(GetMaterial().GetInelasticityType() == 1) {
			FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_hardening_parameter);
			FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_yield_function);
		}
	}

	IntPointsStiffnessIterator ip(this);
	PlaneSymmetricTensorComponents initial_strain;
	if(HasInitialStrains())
		FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_initial_strain,
							FieldVariableDescriptor::FVCI_y, true);
}

// plane part of the total Lagrangian strain tensor,
// covariant components in the local basis of the reference configuration;
// a point is identified by the integration point in the plane and by the local thickness coordinate, which varies from -1 to 1
PlaneSymmetricTensorComponents ANCFThinPlate3D::ComputeTotalStrainLocalComponents(IntegrationPointsIterator & ip, double localZ, const XGProvider & xg)
{
	Matrix grad;
	GetGrad(ip, grad);
	PlaneSymmetricTensorComponents aInv( grad(1,1), grad(1,2), grad(1,3) );
	const Vector2D & ploc = ip.Point2D();

	// actual strains in terms of the shell theory in this point
	PlaneSymmetricTensorComponents C = GetFirstStrainTensorCComponents(ploc, grad, xg);
	PlaneSymmetricTensorComponents K = GetSecondStrainTensorKComponents(ploc, grad, xg);

	// the plane part of the three-dimensional strain tensor is linearly distributed over the thickness
	double z = localZ * GetThicknessAtPoint(ploc, xg.IsDrawConfiguration()) / 2;
	K *= z;
	K += C;
	return K;
}

// converts covariant components of an in-plane tensor into a matrix
// of components in a particular Cartesian basis in the tangent plane
Matrix2D ANCFThinPlate3D::LocalComponentsToCartesian(IntegrationPointsIterator & ip, PlaneSymmetricTensorComponents & localComponents)
{
	PlaneSymmetricTensorComponents a;
	// we read the components of the first metric tensor in the reference configuration from the saved data
	Matrix grad;
	GetGrad(ip, grad);
	for(int alphaBeta = 1; alphaBeta <= 3; alphaBeta++)
		a.SetComponent(alphaBeta, grad(2, alphaBeta));
	// and now apply the derived formulas, which were presented in "Cartesian components for shell elements.docx"
	double adet = a.Det();
	double epsx = localComponents.t11 / a.t11;
	double epsy = (localComponents.t11 * a.t12 * a.t12 / a.t11 + localComponents.t22 * a.t11 - 2 * localComponents.t12 * a.t12) / adet;
	double epsxy = (localComponents.t12 - localComponents.t11 * a.t12 / a.t11) / sqrt(adet);
	return Matrix2D(epsx, epsxy, epsxy, epsy);
}

// converts Cartesian components of an in-plane tensor into local basis (produces covariant components)
PlaneSymmetricTensorComponents ANCFThinPlate3D::CartesianComponentsToLocal(IntegrationPointsIterator & ip, Matrix2D & cartesianComponents)
{
	PlaneSymmetricTensorComponents a;
	// we read the components of the first metric tensor in the reference configuration from the saved data
	Matrix grad;
	GetGrad(ip, grad);
	for(int alphaBeta = 1; alphaBeta <= 3; alphaBeta++)
		a.SetComponent(alphaBeta, grad(2, alphaBeta));
	// and now apply the derived formulas, which were presented in "Cartesian components for shell elements.docx"
	double adet = a.Det();
	return PlaneSymmetricTensorComponents(
		a.t11 * cartesianComponents(1,1),
		a.t12 * cartesianComponents(1,1) + sqrt(adet) * cartesianComponents(1,2),
		(a.t12 * a.t12 * cartesianComponents(1,1) + 2 * sqrt(adet) * a.t12 * cartesianComponents(1,2) + adet * cartesianComponents(2,2)) / a.t11
		);
}

#if 0


// computes local covariant components of the plane part of the 3D strain tensor
// (similar to Green-Lagrange strain components)
// for postprocessing purposes only (searches for nearest integration point around ploc)
PlaneSymmetricTensorComponents ANCFThinPlate3D::GetStrainComponentsD(const Vector3D & ploc, const XGProvider & xg,
												const PlaneSymmetricTensorComponents & a, const PlaneSymmetricTensorComponents & b, StrainType strainType)
{
	// geometrically nonlinear version
	PlaneSymmetricTensorComponents initial_strain, inelastic_strain;
	if(strainType != ST_total_strain)
	{
		// first we seek whether initial strains or plasticity take place
		if (HasInitialStrains() || IsInelasticMaterial())
		{
			if(HasInitialStrains())
			{
				GetInitialStrainComponentsAtPoint(initial_strain, (Vector2D&) ploc, xg.IsDrawConfiguration());
				if(strainType == ST_initial_strain)
				{
					return initial_strain;
				}
			}
			if(IsInelasticMaterial())
			{
				// inelastic strain is not known at the actual point;
				// therefore we simply find the nearest integration point and use the inealstic strain in it;
				// the same holds for the vector of inelastic variables;
				IntPointsStiffnessIterator ip(this);
				ip.GoClosestTo(ploc);

				if (thickness_layers == 1)
				{
					inelastic_strain = GetOldPlasticVariables(ip, 1).strain;
				}
				else   //thickness_layers > 1
				{
					inelastic_strain = GetOldPlasticVariables(ip, ploc.Z(), xg.IsDrawConfiguration()).strain;
				}

				if(strainType == ST_inelastic_strain)
				{
					return inelastic_strain;
				}
			}
		}
	}

	// here either the total strain or the elastic part of strain is computed
	PlaneSymmetricTensorComponents C = GetMetricTensorAComponents((const Vector2D&)ploc, xg);
	C -= a;
	C *= 0.5;
	C -= initial_strain;
	C -= inelastic_strain;
	PlaneSymmetricTensorComponents K = GetMetricTensorBComponents((const Vector2D&)ploc, xg);
	K -= b;
	double thickness = GetThicknessAtPoint((const Vector2D &)ploc, xg.IsDrawConfiguration());
	//$ YV 2013-11-14: instead of dividing by thickness, we shall multiply by thickness
	K *= ploc.Z() * thickness * 2;
	//$ PG 2011-9-21: yet to be done: curvatoric relaxation due to inplane inelastic_strain at ploc.Z() (level-wise)
	C += K;
	return C;
}



double ANCFThinPlate3D::ElasticCoefficient()
{
	return Nu() / (1 - Nu());
}

double ANCFThinPlate3D::ElasticMu()
{
	return Em() / 2 / (1 + Nu());
}


// computes local components of the 3D stress tensor (similar to second Piola-Kirchhoff stress tensor)
// for drawing purposes only
PlaneSymmetricTensorComponents ANCFThinPlate3D::GetStressComponentsD(const Vector3D & ploc, const XGProvider & xg,
												 const PlaneSymmetricTensorComponents & aInv,
												 const PlaneSymmetricTensorComponents & a, const PlaneSymmetricTensorComponents & b)
{
	PlaneSymmetricTensorComponents eps = GetStrainComponentsD(ploc, xg, a, b, ST_elastic_strain);
	PlaneSymmetricTensorComponents sigma1 = eps;
	double mu = ElasticMu();
	sigma1 *= 2 * mu;
	PlaneSymmetricTensorComponents sigma2 = a;
	sigma2 *= eps.Trace(aInv) * 2 * mu * ElasticCoefficient();
	sigma1 += sigma2;
	return sigma1;
}

// recomputes the components in the local basis into a corresponding physical component of the Cartesian basis;
// the tensor may have an out-of-plane part Tnn (no "shear" components are needed for a homogeneous material)
double ANCFThinPlate3D::GetCartesianComponent(
							const Vector3D & ploc,
							const PlaneSymmetricTensorComponents & T, double Tnn,
							FieldVariableDescriptor::FieldVariableComponentIndex component1,
							FieldVariableDescriptor::FieldVariableComponentIndex component2,
							const PlaneSymmetricTensorComponents & aInv
							)
{
	if(component1 == FieldVariableDescriptor::FVCI_magnitude)	// sqrt(T..T)
		return sqrt(T.SelfConvolute(0, 1, aInv) + Tnn * Tnn);

	Vector3D r1 = GetRAlpha((const Vector2D&)ploc, 1, xgInit);
	Vector3D r2 = GetRAlpha((const Vector2D&)ploc, 2, xgInit);

	// components of the local basis vectors of the reference configuration in the Cartesian basis
	double r1i = r1(component1);
	double r2i = r2(component1);
	double r1j = r1(component2);
	double r2j = r2(component2);

	// components of the contravariant local basis vectors in the Cartesian basis
	double rc1i = aInv.t11 * r1i + aInv.t12 * r2i;
	double rc2i = aInv.t12 * r1i + aInv.t22 * r2i;
	double rc1j = aInv.t11 * r1j + aInv.t12 * r2j;
	double rc2j = aInv.t12 * r1j + aInv.t22 * r2j;

	double result = T.t11 * rc1i * rc1j + T.t12 * rc1i * rc2j +
									T.t12 * rc2i * rc1j + T.t22 * rc2i * rc2j;

	if(Tnn != 0)
	{
		Vector3D n = r1.Cross(r2);
		n.Normalize();
		result += Tnn * n(component1) * n(component2);
	}

	return result;
}

#endif

double ANCFThinPlate3D::GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & local_position, bool flagD)
{
	if(flagD)
		return GetFieldVariableValue(fvd, local_position, xgDraw);
	return GetFieldVariableValue(fvd, local_position, xgCompute);
}

double ANCFThinPlate3D::GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & drawing_local_position, const XGProvider & xg)
{
	Vector3D local_position(drawing_local_position);

	////$ PG 2012-10-16: [ test!!!
	//IntPointsStiffnessIterator ip(this);
	//ip.GoClosestTo(drawing_local_position);
	//local_position = ip.Point();
	////$ PG 2012-10-16: ]
	
	switch(fvd.VariableType())
	{
	case FieldVariableDescriptor::FVT_displacement: return fvd.GetComponent(GetDisplacement(local_position, xg));
	case FieldVariableDescriptor::FVT_position: return fvd.GetComponent(GetPos(local_position, xg));
	case FieldVariableDescriptor::FVT_velocity: return fvd.GetComponent(GetVel(local_position, xg));
	case FieldVariableDescriptor::FVT_acceleration: return fvd.GetComponent(local_position);
	case FieldVariableDescriptor::FVT_stress:
	case FieldVariableDescriptor::FVT_stress_mises:
	case FieldVariableDescriptor::FVT_total_strain:
	case FieldVariableDescriptor::FVT_inelastic_strain:
	case FieldVariableDescriptor::FVT_initial_strain:
		{
			PlaneSymmetricTensorComponents a = GetMetricTensorAComponents((const Vector2D&)local_position, xgInit);
			PlaneSymmetricTensorComponents b = GetMetricTensorBComponents((const Vector2D&)local_position, xgInit);
			PlaneSymmetricTensorComponents aInv = a.Inverse();

			/////////////////////////////////////////////////// TEST - testing the cartesian components
#if 0
			if(mbs->GetSimulationStatus().GetStatusFlag(TSimulationEndedRegularly))
			{
				static bool flagPrintedOut = false;
				if(!flagPrintedOut)
				{
					flagPrintedOut = true;
					mbs->UO(UO_LVL_0) << "testing the conversion between the cartesian and local components\n";
					IntPointsStiffnessIterator ip(this);
					ip.GoClosestTo(local_position);
					PlaneSymmetricTensorComponents strainLocal = ComputeTotalStrainLocalComponents(ip, xg);
					Matrix2D strainGlobal = LocalComponentsToCartesian(ip, strainLocal);
					PlaneSymmetricTensorComponents strainLocalBack = CartesianComponentsToLocal(ip, strainGlobal);
					mbs->UO(UO_LVL_0) << "local 1: " << strainLocal.t11 << " " << strainLocal.t12 << " " << strainLocal.t22 << "\n";
					mbs->UO(UO_LVL_0) << "local 2: " << strainLocalBack.t11 << " " << strainLocalBack.t12 << " " << strainLocalBack.t22 << "\n";
					mbs->UO(UO_LVL_0) << "global: " << strainGlobal(1,1) << " " << strainGlobal(1,2) << " " << strainGlobal(2,2) << "\n";
				}
			}
#endif

			if(fvd.VariableType() == FieldVariableDescriptor::FVT_stress || fvd.VariableType() == FieldVariableDescriptor::FVT_stress_mises)
			{
				// TODO: by now we ignore the inelastic strains, which needs to be implemented
				IntPointsStiffnessIterator ip(this);
				ip.GoClosestTo(local_position);
				Matrix2D strain = LocalComponentsToCartesian(ip, ComputeTotalStrainLocalComponents(ip, local_position.Z(), xg));
				Matrix2D stress;
				GetMaterial().ComputeStressFromStrain2D(strain, stress);
				if(fvd.VariableType() == FieldVariableDescriptor::FVT_stress_mises)
				{
					double trace = stress.Trace();
					stress(1,1) -= trace / 3;
					stress(2,2) -= trace / 3;
					// in the inner product of the deviator of the stress tensor with itself we are taking the out-of-plane part into account
					return sqrt(1.5 * (stress.InnerProduct(stress) + trace * trace / 9));
				}
				return fvd.GetComponent(stress);
			}
			if(fvd.VariableType() == FieldVariableDescriptor::FVT_total_strain)
			{
				IntPointsStiffnessIterator ip(this);
				ip.GoClosestTo(local_position);
				Matrix2D strain = LocalComponentsToCartesian(ip, ComputeTotalStrainLocalComponents(ip, local_position.Z(), xg));
				return fvd.GetComponent(strain);
			}
			if(fvd.VariableType() == FieldVariableDescriptor::FVT_initial_strain)
			{
				IntPointsStiffnessIterator ip(this);
				ip.GoClosestTo(local_position);
				Matrix2D initial_strain_cartesian;
				GetInitialStrainComponentsAtPoint(initial_strain_cartesian, ip.Point2D(), true);
				return fvd.GetComponent(initial_strain_cartesian);
			}
			if(fvd.VariableType() == FieldVariableDescriptor::FVT_inelastic_strain)
			{
				// TODO!!!!!!!
				return 0;
			}
		}
	case FieldVariableDescriptor::FVT_hardening_parameter:
		{
#if 0		// TODO
			IntPointsStiffnessIterator ip(this);
			ip.GoClosestTo(local_position);
			if (thickness_layers == 1)
			{
				return GetOldPlasticVariables(ip, 1).hardening_parameter;
			}
			return GetOldPlasticVariables(ip, local_position.Z(), xg.IsDrawConfiguration()).hardening_parameter;
#endif
			return 0;
		}
	case FieldVariableDescriptor::FVT_yield_function:
		{
#if 0		// TODO
			IntPointsStiffnessIterator ip(this);
			ip.GoClosestTo(local_position);
			if (thickness_layers == 1)
			{
				return GetOldPlasticVariables(ip, 1).yield_function;
			}
			return GetOldPlasticVariables(ip, local_position.Z(), xg.IsDrawConfiguration()).yield_function;
#endif
			return 0;
		}
	}
	return FIELD_VARIABLE_NO_VALUE;
}

void ANCFThinPlate3D::DrawElement() 
{
	XGProviderCached<FEmaxDOF> xgCache(xgDraw, XGLength());

	mbs->SetColor(col);
	double lx1 = GetMBS()->GetDOption(106); 
	double ly1 = GetMBS()->GetDOption(106); 
	double lz1 = GetMBS()->GetDOption(106)*GetMBS()->GetMagnifyYZ();
	int chladni = 0;//GetMBS()->GetIOption(116); //for chladni figures in eigenmodes
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
		linemode = 3;
	}

	int colormode = 0;
	if (GetMBS()->GetActualPostProcessingFieldVariable() != NULL)
		colormode = 1;

	double tilex = GetMBS()->GetIOption(137);
	if (!colormode)
		tilex = GetMBS()->GetIOption(136);

	double tiley = tilex;

	static TArray<Vector3D> points;
	points.SetLen((int)(tilex+1)*(int)(tiley+1));
	static TArray<Vector3D> normals;
	normals.SetLen((int)(tilex+1)*(int)(tiley+1));
	static TArray<double> vals;
	vals.SetLen((int)(tilex+1)*(int)(tiley+1));

	int starts = 1;
	int ends = 6;
	if (chladni || drawflat)
	{
		starts = 3;
		ends = 3;
	}

	for (int side = starts; side <= ends; side++)
	{
		points.SetLen(0);
		vals.SetLen(0);
		normals.SetLen(0);
		Vector3D p0, vx, vy;
		int tileyn = (int)tiley;
		if (side <= 2 || side >= 5)
			tileyn = GetMBS()->GetIOption(138);

		switch(side)
		{
		case 1:
			{ //bottom
				p0 = Vector3D(-lx1,-ly1,-lz1);
				vx = Vector3D(2.*lx1/tilex,0,0);
				vy = Vector3D(0,0,2.*lz1/tileyn);
				break;
			}
		case 2:
			{ //top
				p0 = Vector3D(-lx1, ly1, lz1);
				vx = Vector3D(2.*lx1/tilex,0,0);
				vy = Vector3D(0,0,-2.*lz1/tileyn);
				break;
			}
		case 3:
			{ //front
				p0 = Vector3D(-lx1,-ly1, lz1);
				vx = Vector3D(2.*lx1/tilex,0,0);
				vy = Vector3D(0,2.*ly1/tileyn,0);
				if (chladni || drawflat)
					p0 = Vector3D(-lx1,-ly1, 0);
				break;
			}
		case 4:
			{ //back
				p0 = Vector3D(-lx1, ly1,-lz1);
				vx = Vector3D(2.*lx1/tilex,0,0);
				vy = Vector3D(0,-2.*ly1/tileyn,0);
				break;
			}
		case 5:
			{ //left
				p0 = Vector3D(-lx1, ly1,-lz1);
				vx = Vector3D(0,-2.*ly1/tilex,0);
				vy = Vector3D(0,0,2.*lz1/tileyn);
				break;
			}
		case 6:
			{ //right
				p0 = Vector3D( lx1,-ly1,-lz1);
				vx = Vector3D(0,2.*ly1/tilex,0);
				vy = Vector3D(0,0,2.*lz1/tileyn);
				break;
			}
		}

		for (double iy = 0; iy <= tileyn+1e-10; iy++)
		{
			for (double ix = 0; ix <= tilex+1e-10; ix++)
			{
				if (!chladni || !colormode)
				{
					Vector3D ploc = (p0+ix*vx+iy*vy);
					points.Add(GetPosD(ploc));

					if (side == 3)
						normals.Add(GetNormalD(ploc));
					if (side == 4)
						normals.Add(-1.*GetNormalD(ploc));

					if (colormode)
					{
						Vector3D plocActual = ploc;
						//plocActual.Z() = 0;		// if this line is uncommented, then the post-processing values are referred to the mid-surface
						vals.Add(GetFieldVariableValue(*GetMBS()->GetActualPostProcessingFieldVariable(), ploc, xgCache));
					}
					else
						vals.Add(0);
				}
				else
				{
					//draw chladni figures:
					Vector3D ploc = (p0+ix*vx+iy*vy);
					points.Add(GetPos(ploc, xgInit));
					normals.Add(GetNormal((const Vector2D&)ploc, xgInit));
					vals.Add(GetDisplacementD(ploc).Norm());
				}
			}
		}
		for (int i=1; i <= normals.Length(); i++)
			points.Add(normals(i));

		FitPointsInDrawingRegion(points, (int)tilex+1, (int)tileyn+1);		// in case the derived class wishes to control the drawing region

		mbs->DrawColorQuads(points,vals,(int)tilex+1,(int)tileyn+1,colormode,linemode);
	}
}

#pragma endregion