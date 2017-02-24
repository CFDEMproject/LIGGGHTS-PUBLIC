//#**************************************************************
//#
//# filename:             IntegrationRule.cpp
//#
//# author:               YV
//#
//# generated:						october 2010
//# description:          Management of the integration points for different finite elements
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
 
#include "mbs_interface.h"
#include "IntegrationRule.h"
#include "femathhelperfunctions.h"

//$ YV 2013-01-12: the library of integration rules is common for all providers of integration rules
IntegrationRulesLibrary IntegrationRule::IntegrationRuleProvider::integrationRulesLibrary;

bool IntegrationRule::IntegrationRuleSettings::operator==(const IntegrationRuleSettings & settings) const
{
	return
		elementType == settings.elementType &&
		geometricNonlinearityStatus == settings.geometricNonlinearityStatus &&
		integratedValueType == settings.integratedValueType &&
		interpolationOrder == settings.interpolationOrder &&
		flags == settings.flags;
}

IntegrationRule::IntegrationRuleSettings::IntegrationRuleSettings(
			TFiniteElementType elementType,
			int interpolationOrder,
			IntegratedValueType integratedValueType,
			GeometricNonlinearityStatus geometricNonlinearityStatus,
			short flags
			)
{
	this->elementType = elementType;
	this->interpolationOrder = interpolationOrder;
	this->integratedValueType = integratedValueType;
	this->geometricNonlinearityStatus = geometricNonlinearityStatus;
	this->flags = flags;
}

IntegrationRulesLibrary::~IntegrationRulesLibrary()
{
	for(int i = 1; i <= integrationRules.Length(); i++)
		delete integrationRules(i);
}

IntegrationRule * IntegrationRulesLibrary::GetIntegrationRule(
		const IntegrationRule::IntegrationRuleSettings & settings,
		IntegrationRule::IntegrationRuleProvider * pIntegrationRuleProvider /*= NULL*/
		)
{
	// if this search once gets slow because of too many integration rules,
	// it can be accelerated by introducing an integer hash value for IntegrationRuleSettings
	// and making a dictionary (map) of the kind hash => integr. rule
	for(int i = 1; i <= integrationRules.Length(); i++)
		if(settings == integrationRules(i)->settings)
			return integrationRules(i);
	if(pIntegrationRuleProvider == NULL)
		return NULL;
	IntegrationRule * newIntegrationRule = new IntegrationRule();
	newIntegrationRule->settings = settings;
	pIntegrationRuleProvider->DefineIntegrationRule(*newIntegrationRule);
	integrationRules.Add(newIntegrationRule);		// AH: add new integration rule to set of existing rules
	return newIntegrationRule;
}


///////////////////////////////////////////////////////////////////////////////////
// Iterator

/// This version is fast, but not safe;
/// it is based on the particular structure of Vector3D and started working
/// as the virtual destructor was removed from Vector3D.
/// The structure of Vector3D is checked at the run time.
const Vector3D & IntegrationPointsIterator::Point() const
{
	return (Vector3D &)(pIntegrationRule->integrationPoints(index));
}

const Vector2D & IntegrationPointsIterator::Point2D() const
{
	return (Vector2D &)(pIntegrationRule->integrationPoints(index));
}

double IntegrationPointsIterator::Weight()
{
	return pIntegrationRule->integrationPoints(index).weight;
}


//////////////////////////////////////
// helper functions

void IntegrationPointsIterator::GoClosestTo(const Vector3D& ploc)
{
	// bring (*this) to integration point which is closest to ploc
	double dist = 10;
	index = 1;
	int index_closest = index;
	
	while(!IsEnd())
	{
		double current_dist = (ploc - Point()).Norm();
		if(current_dist < dist)
		{
			dist = current_dist;
			index_closest = index;
		}
		++index;
	}
	index = index_closest;
}

void IntegrationRule::DefineIntegrationRuleSquare(IntegrationRule & integrationRule, int ruleOrder)
{
	// we use a ready function, which provides a one-dimensional integration rule
	Vector x;
	Vector w;
	GetIntegrationRule(x, w, ruleOrder);

	// and now we create a 2D integration rule
	integrationRule.integrationPoints.SetLen(x.Length() * x.Length());
	int cnt = 1;
	for (int i1 = 1; i1 <= x.GetLen(); i1++)
	{
		for (int i2 = 1; i2 <= x.GetLen(); i2++)
		{
			integrationRule.integrationPoints(cnt).x = x(i1);
			integrationRule.integrationPoints(cnt).y = x(i2);
			integrationRule.integrationPoints(cnt).z = 0;
			integrationRule.integrationPoints(cnt).weight = w(i1) * w(i2);
			cnt++;
		}
	}
}

void IntegrationRule::DefineIntegrationRuleLobattoSquare(IntegrationRule & integrationRule, int ruleOrder)
{
	// we use a ready function, which provides a one-dimensional integration rule
	Vector x;
	Vector w;
	GetIntegrationRuleLobatto(x, w, ruleOrder);

	// and now we create a 2D integration rule
	integrationRule.integrationPoints.SetLen(x.Length() * x.Length());
	int cnt = 1;
	for (int i1 = 1; i1 <= x.GetLen(); i1++)
	{
		for (int i2 = 1; i2 <= x.GetLen(); i2++)
		{
			integrationRule.integrationPoints(cnt).x = x(i1);
			integrationRule.integrationPoints(cnt).y = x(i2);
			integrationRule.integrationPoints(cnt).z = 0;
			integrationRule.integrationPoints(cnt).weight = w(i1) * w(i2);
			cnt++;
		}
	}
}

void IntegrationRule::DefineIntegrationRuleSquareAnisotropic(IntegrationRule & integrationRule, int ruleOrder_x, int ruleOrder_y)
{
	// we use a ready function, which provides a one-dimensional integration rule
	// we use different integration rules for x and y direction
	Vector xx, xy;
	Vector wx, wy;
	GetIntegrationRule(xx, wx, ruleOrder_x);
	GetIntegrationRule(xy, wy, ruleOrder_y);

	// and now we create a 2D integration rule
	integrationRule.integrationPoints.SetLen(xx.Length() * xy.Length());
	int cnt = 1;
	for (int i1 = 1; i1 <= xx.GetLen(); i1++)
	{
		for (int i2 = 1; i2 <= xy.GetLen(); i2++)
		{
			integrationRule.integrationPoints(cnt).x = xx(i1);
			integrationRule.integrationPoints(cnt).y = xy(i2);
			integrationRule.integrationPoints(cnt).z = 0;
			integrationRule.integrationPoints(cnt).weight = wx(i1) * wy(i2);
			cnt++;
		}
	}
}

void IntegrationRule::DefineIntegrationRuleTriangle(IntegrationRule & integrationRule, int ruleOrder)
{
	// we use a ready function, which provides a triangular integration rule
	Vector x1, x2, w;
	GetIntegrationRuleTrig(x1, x2, w, ruleOrder);

	// and now we create a 3D integration rule
	integrationRule.integrationPoints.SetLen(x1.Length());
	assert(x1.Length() == x2.Length());
	assert(x1.Length() == w.Length());
	for (int i1 = 1; i1 <= x1.GetLen(); i1++)
	{
		integrationRule.integrationPoints(i1).x = x1(i1);
		integrationRule.integrationPoints(i1).y = x2(i1);
		integrationRule.integrationPoints(i1).z = 0;
		integrationRule.integrationPoints(i1).weight = w(i1);
	}
}

void IntegrationRule::DefineIntegrationRuleCube(IntegrationRule & integrationRule, int ruleOrder)
{
	// we use a ready function, which provides a one-dimensional integration rule
	Vector x;
	Vector w;
	GetIntegrationRule(x ,w, ruleOrder);

	// and now we create a 3D integration rule
	integrationRule.integrationPoints.SetLen(x.Length() * x.Length() * x.Length());
	int cnt = 1;
	for (int i1 = 1; i1 <= x.GetLen(); i1++)
	{
		for (int i2 = 1; i2 <= x.GetLen(); i2++)
		{
			for (int i3 = 1; i3 <= x.GetLen(); i3++)
			{
				integrationRule.integrationPoints(cnt).x = x(i1);
				integrationRule.integrationPoints(cnt).y = x(i2);
				integrationRule.integrationPoints(cnt).z = x(i3);
				integrationRule.integrationPoints(cnt).weight = w(i1) * w(i2) * w(i3);
				cnt++;
			}
		}
	}
}

void IntegrationRule::DefineIntegrationRuleTetrahedron(IntegrationRule & integrationRule, int ruleOrder)
{
	// we use a ready function, which provides a tetrahedral integration rule
	Vector x1, x2, x3, w;
	GetIntegrationRuleTet(x1, x2, x3, w, ruleOrder);

	// and now we create a 3D integration rule
	integrationRule.integrationPoints.SetLen(x1.Length());
	assert(x1.Length() == x2.Length());
	assert(x1.Length() == x3.Length());
	assert(x1.Length() == w.Length());
	for (int i1 = 1; i1 <= x1.GetLen(); i1++)
	{
		integrationRule.integrationPoints(i1).x = x1(i1);
		integrationRule.integrationPoints(i1).y = x2(i1);
		integrationRule.integrationPoints(i1).z = x3(i1);
		integrationRule.integrationPoints(i1).weight = w(i1);
	}
}


//$EK 2013-03-05 added
void IntegrationRule::DefineIntegrationRulePrism(IntegrationRule & integrationRule, int ruleOrder)
{
	// we use a ready function, which provides a tetrahedral integration rule
	Vector x1, x2, xline, w, wline;
	GetIntegrationRuleTrig(x1, x2, w, ruleOrder);
	GetIntegrationRule(xline, wline, ruleOrder);
	for (int i = 1; i <= xline.Length(); ++i)
		xline(i) = (xline(i) + 1.)/2.;
	
	// and now we create a 3D integration rule
	integrationRule.integrationPoints.SetLen(x1.Length()*xline.Length());
	assert(x1.Length() == x2.Length());
	assert(x1.Length() == w.Length());
	assert(xline.Length() == wline.Length());
	int cnt = 1;
	for (int i1 = 1; i1 <= x1.GetLen(); i1++)
	{
		for (int i2 = 1; i2 <= xline.GetLen(); ++i2)
		{
			integrationRule.integrationPoints(cnt).x = x1(i1);
			integrationRule.integrationPoints(cnt).y = x2(i1);
			integrationRule.integrationPoints(cnt).z = xline(i2);
			//the /2. is due to the transformation of [-1,1] to [0,1]
			integrationRule.integrationPoints(cnt).weight = w(i1)*wline(i2)/2.;
			cnt++;
	  }
	}
}

//$EK 2013-03-05 added
void IntegrationRule::DefineIntegrationRulePyramid(IntegrationRule & integrationRule, int ruleOrder)
{
	// we use a ready function, which provides a tetrahedral integration rule
	Vector x_z, xline, w_z, wline;
	//$EK I don't know why in z-direction one needs in increased order???
	//-------------------------------------------------------------------------
	//FROM NGSolve ... prismatic integration rule (in intrule.cpp)
	//		const IntegrationRule & quadrule = SelectIntegrationRule (ET_QUAD, order);
	//	const IntegrationRule & segrule = SelectIntegrationRule (ET_SEGM, order+2);
	//const IntegrationPoint & ipquad = quadrule[i]; // .GetIP(i);
	//const IntegrationPoint & ipseg = segrule[j]; // .GetIP(j);
	//point[0] = (1-ipseg(0)) * ipquad(0);
	//point[1] = (1-ipseg(0)) * ipquad(1);
	//point[2] = ipseg(0);
	//weight = ipseg.Weight() * sqr (1-ipseg(0)) * ipquad.Weight();
	//--------------------------------------------------------------------------
	GetIntegrationRule(x_z, w_z, ruleOrder+2);
	GetIntegrationRule(xline, wline, ruleOrder);
	//adapt integration points to interval [0,1] instead of [-1,1]
	for (int i = 1; i <= xline.Length(); ++i)
		xline(i) = (xline(i) + 1.)/2.;
	for (int i = 1; i <= x_z.Length(); ++i)
		x_z(i) = (x_z(i) + 1.)/2.;
	
	// and now we create a 3D integration rule
	integrationRule.integrationPoints.SetLen(x_z.Length()*sqr(xline.Length()));
	assert(x_z.Length() == w_z.Length());
	assert(xline.Length() == wline.Length());
	int cnt = 1;
	for (int i = 1; i <= xline.GetLen(); i++)
		for (int j = 1; j <= xline.GetLen(); j++)
			for (int k = 1; k <= x_z.GetLen(); ++k)
			{
				double z = x_z(k);
				integrationRule.integrationPoints(cnt).x = xline(i)*(1-z);
				integrationRule.integrationPoints(cnt).y = xline(j)*(1-z);
				integrationRule.integrationPoints(cnt).z = z;
				//the /8. is due to the transformation of [-1,1] to [0,1] 
				integrationRule.integrationPoints(cnt).weight = w_z(k) * sqr(1- z) * wline(i) * wline(j)/8.;
				cnt++;
			}
}


void IntegrationRule::DefineIntegrationRuleLine(IntegrationRule & integrationRule, int ruleOrder)
{
	// we use a ready function, which provides a one-dimensional integration rule
	Vector x, w;
	GetIntegrationRule(x, w, ruleOrder);

	integrationRule.integrationPoints.SetLen(x.Length());
	for (int i = 1; i <= x.GetLen(); i++)
	{
		integrationRule.integrationPoints(i).x = x(i);
		integrationRule.integrationPoints(i).y = 0.;
		integrationRule.integrationPoints(i).z = 0.;
		integrationRule.integrationPoints(i).weight = w(i);
	}
}

#if 0
///////////////////////////////////////////////////////////////////////////////////
// Populating the tables with ready-made data.
// Possibly, such a strategy is not optimal as we create much redundant information

class IntegrationRulesCreator
{
	IntegrationRule::IntegrationRuleInfo * currentIntegrationRuleInfo;
	void AddCurrentRuleInfo();

	/// Adds all hexahedral rules.
	void AddHexahedralRules();
	/// Adds a single hexahedral rule of the specified order (ruleOrder) for the given settings.
	void AddHexahedralRule(int interpolationOrder,
		IntegrationRule::IntegratedValueType ivt, GeometricNonlinearityStatus gns, int ruleOrder);
	/// Actually creates a particular rule for the given order.
	void CreateHexahedralRule(int order);

	/// Adds all tertahedral rules.
	void AddTetrahedralRules();
	/// Adds a single tertahedral rule of the specified order (ruleOrder) for the given settings.
	void AddTetrahedralRule(int interpolationOrder,
		IntegrationRule::IntegratedValueType ivt, GeometricNonlinearityStatus gns, int ruleOrder);
	/// Actually creates a particular rule for the given order.
	void CreateTetrahedralRule(int order);

	/// Adds all quadrilateral rules.
	void AddQuadrilateralRules();
	/// Adds a single quadrilateral rule of the specified order (ruleOrder) for the given settings.
	void AddQuadrilateralRule(int interpolationOrder,
		IntegrationRule::IntegratedValueType ivt, GeometricNonlinearityStatus gns, int ruleOrder);
	/// Actually creates a particular rule for the given order.
	void CreateQuadrilateralRule(int order);

	/// Adds all thin plate rules.
	void AddThinPlateRules();
	/// Adds a single thin plate rule of the specified order (ruleOrder) for the given settings.
	void AddThinPlateRule(int interpolationOrder,
		IntegrationRule::IntegratedValueType ivt, GeometricNonlinearityStatus gns, int ruleOrder);

public:
	IntegrationRulesCreator();
};

IntegrationRulesCreator dummy;	// creating this variable, the system populates the array with the available data

IntegrationRulesCreator::IntegrationRulesCreator()
{
	// these tests are needed to ensure the functionality of the fast access
	// of the integration point coordinates in the class IntegrationPointsIterator.
	assert(sizeof(Vector3D) == 3 * sizeof(double));
	assert(sizeof(Vector2D) == 2 * sizeof(double));

	AddHexahedralRules();
	AddTetrahedralRules();
	AddQuadrilateralRules();
	AddThinPlateRules();
}

void IntegrationRulesCreator::AddCurrentRuleInfo()
{
	IntegrationRule::integrationRules.Add(currentIntegrationRuleInfo);
}

#pragma region Hexahedral rules

void IntegrationRulesCreator::CreateHexahedralRule(int order)
{
	// we use a ready function, which provides a one-dimensional integration rule
	Vector x;
	Vector w;
	GetIntegrationRule(x,w,order);

	// and now we create a 3D integration rule
	currentIntegrationRuleInfo = new IntegrationRule::IntegrationRuleInfo();
	currentIntegrationRuleInfo->integrationRule.SetLen(x.Length() * x.Length() * x.Length());
	int cnt = 1;
	for (int i1 = 1; i1 <= x.GetLen(); i1++)
	{
		for (int i2 = 1; i2 <= x.GetLen(); i2++)
		{
			for (int i3 = 1; i3 <= x.GetLen(); i3++)
			{
				currentIntegrationRuleInfo->integrationRule(cnt).x = x(i1);
				currentIntegrationRuleInfo->integrationRule(cnt).y = x(i2);
				currentIntegrationRuleInfo->integrationRule(cnt).z = x(i3);
				currentIntegrationRuleInfo->integrationRule(cnt).weight = w(i1) * w(i2) * w(i3);
				cnt++;
			}
		}
	}
}

void IntegrationRulesCreator::AddHexahedralRule(int interpolationOrder,
		IntegrationRule::IntegratedValueType ivt, GeometricNonlinearityStatus gns, int ruleOrder)
{
	CreateHexahedralRule(ruleOrder);
	currentIntegrationRuleInfo->settings = IntegrationRule::IntegrationRuleSettings(
		TFE_Hexahedral,interpolationOrder,ivt,gns,false
		);
	AddCurrentRuleInfo();
}

void IntegrationRulesCreator::AddHexahedralRules()
{
	// we start with cubic interpolation - 8 nodes hexahedral
	// geometric linear case
	AddHexahedralRule(3, IntegrationRule::IVT_Stiffness, GNS_Linear, 3);
	AddHexahedralRule(3, IntegrationRule::IVT_Mass, GNS_Linear, 2);
	AddHexahedralRule(3, IntegrationRule::IVT_Load, GNS_Linear, 2);
	// geometric nonlinear, small strain - same rule orders
	AddHexahedralRule(3, IntegrationRule::IVT_Stiffness, GNS_NonlinearSmallStrain, 3);
	AddHexahedralRule(3, IntegrationRule::IVT_Mass, GNS_NonlinearSmallStrain, 2);
	AddHexahedralRule(3, IntegrationRule::IVT_Load, GNS_NonlinearSmallStrain, 2);
	// geometric nonlinear, large strain - same rule orders
	AddHexahedralRule(3, IntegrationRule::IVT_Stiffness, GNS_NonlinearLargeStrain, 3);
	AddHexahedralRule(3, IntegrationRule::IVT_Mass, GNS_NonlinearLargeStrain, 2);
	AddHexahedralRule(3, IntegrationRule::IVT_Load, GNS_NonlinearLargeStrain, 2);

	// fourth order interpolation - 20 nodes hehaxedral
	// geometric linear case
	AddHexahedralRule(4, IntegrationRule::IVT_Stiffness, GNS_Linear, 4);
	AddHexahedralRule(4, IntegrationRule::IVT_Mass, GNS_Linear, 3);
	AddHexahedralRule(4, IntegrationRule::IVT_Load, GNS_Linear, 2);
	// geometric nonlinear, small strain - same rule orders
	AddHexahedralRule(4, IntegrationRule::IVT_Stiffness, GNS_NonlinearSmallStrain, 4);
	AddHexahedralRule(4, IntegrationRule::IVT_Mass, GNS_NonlinearSmallStrain, 3);
	AddHexahedralRule(4, IntegrationRule::IVT_Load, GNS_NonlinearSmallStrain, 2);
	// geometric nonlinear, large strain
	AddHexahedralRule(4, IntegrationRule::IVT_Stiffness, GNS_NonlinearLargeStrain, 7);
	AddHexahedralRule(4, IntegrationRule::IVT_Mass, GNS_NonlinearLargeStrain, 4);
	AddHexahedralRule(4, IntegrationRule::IVT_Load, GNS_NonlinearLargeStrain, 2);
}

#pragma endregion

#pragma region Tetrahedral rules

void IntegrationRulesCreator::CreateTetrahedralRule(int order)
{
	// we use a ready function, which provides a one-dimensional integration rule
	Vector x1, x2, x3, w;
	GetIntegrationRuleTet(x1, x2, x3, w, order);

	// and now we create a 3D integration rule
	currentIntegrationRuleInfo = new IntegrationRule::IntegrationRuleInfo();
	currentIntegrationRuleInfo->integrationRule.SetLen(x1.Length());
	assert(x1.Length() == x2.Length());
	assert(x1.Length() == x3.Length());
	assert(x1.Length() == w.Length());
	for (int i1 = 1; i1 <= x1.GetLen(); i1++)
	{
		currentIntegrationRuleInfo->integrationRule(i1).x = x1(i1);
		currentIntegrationRuleInfo->integrationRule(i1).y = x2(i1);
		currentIntegrationRuleInfo->integrationRule(i1).z = x3(i1);
		currentIntegrationRuleInfo->integrationRule(i1).weight = w(i1);
	}
}

void IntegrationRulesCreator::AddTetrahedralRule(int interpolationOrder,
		IntegrationRule::IntegratedValueType ivt, GeometricNonlinearityStatus gns, int ruleOrder)
{
	CreateTetrahedralRule(ruleOrder);
	currentIntegrationRuleInfo->settings = IntegrationRule::IntegrationRuleSettings(
		TFE_Tetrahedral,interpolationOrder,ivt,gns,false
		);
	AddCurrentRuleInfo();
}

void IntegrationRulesCreator::AddTetrahedralRules()
{
	// we start with cubic interpolation - 4 nodes tetrahedral
	// geometric linear case
	AddTetrahedralRule(1, IntegrationRule::IVT_Stiffness, GNS_Linear, 3);
	AddTetrahedralRule(1, IntegrationRule::IVT_Mass, GNS_Linear, 2);
	AddTetrahedralRule(1, IntegrationRule::IVT_Load, GNS_Linear, 2);
	// geometric nonlinear, small strain - same rule orders
	AddTetrahedralRule(1, IntegrationRule::IVT_Stiffness, GNS_NonlinearSmallStrain, 3);
	AddTetrahedralRule(1, IntegrationRule::IVT_Mass, GNS_NonlinearSmallStrain, 2);
	AddTetrahedralRule(1, IntegrationRule::IVT_Load, GNS_NonlinearSmallStrain, 2);
	// geometric nonlinear, large strain - same rule orders
	AddTetrahedralRule(1, IntegrationRule::IVT_Stiffness, GNS_NonlinearLargeStrain, 3);
	AddTetrahedralRule(1, IntegrationRule::IVT_Mass, GNS_NonlinearLargeStrain, 2);
	AddTetrahedralRule(1, IntegrationRule::IVT_Load, GNS_NonlinearLargeStrain, 2);

	// fourth order interpolation - 8 nodes tetraxedral
	// geometric linear case
	AddTetrahedralRule(2, IntegrationRule::IVT_Stiffness, GNS_Linear, 4);
	AddTetrahedralRule(2, IntegrationRule::IVT_Mass, GNS_Linear, 3);
	AddTetrahedralRule(2, IntegrationRule::IVT_Load, GNS_Linear, 2);
	// geometric nonlinear, small strain - same rule orders
	AddTetrahedralRule(2, IntegrationRule::IVT_Stiffness, GNS_NonlinearSmallStrain, 4);
	AddTetrahedralRule(2, IntegrationRule::IVT_Mass, GNS_NonlinearSmallStrain, 3);
	AddTetrahedralRule(2, IntegrationRule::IVT_Load, GNS_NonlinearSmallStrain, 2);
	// geometric nonlinear, large strain
	AddTetrahedralRule(2, IntegrationRule::IVT_Stiffness, GNS_NonlinearLargeStrain, 5);
	AddTetrahedralRule(2, IntegrationRule::IVT_Mass, GNS_NonlinearLargeStrain, 4);
	AddTetrahedralRule(2, IntegrationRule::IVT_Load, GNS_NonlinearLargeStrain, 2);
}

#pragma endregion

#pragma region Quadrilateral rules

void IntegrationRulesCreator::AddQuadrilateralRules()
{
	// we start with quadratic interpolation - 4 nodes quadrilateral
	// geometric linear case
	AddQuadrilateralRule(2, IntegrationRule::IVT_Stiffness, GNS_Linear, 2);
	AddQuadrilateralRule(2, IntegrationRule::IVT_Mass, GNS_Linear, 2);
	AddQuadrilateralRule(2, IntegrationRule::IVT_Load, GNS_Linear, 2);
	// geometric nonlinear, small strain - same rule orders
	AddQuadrilateralRule(2, IntegrationRule::IVT_Stiffness, GNS_NonlinearSmallStrain, 2);
	AddQuadrilateralRule(2, IntegrationRule::IVT_Mass, GNS_NonlinearSmallStrain, 2);
	AddQuadrilateralRule(2, IntegrationRule::IVT_Load, GNS_NonlinearSmallStrain, 2);
	// geometric nonlinear, large strain - same rule orders
	AddQuadrilateralRule(2, IntegrationRule::IVT_Stiffness, GNS_NonlinearLargeStrain, 2);
	AddQuadrilateralRule(2, IntegrationRule::IVT_Mass, GNS_NonlinearLargeStrain, 2);
	AddQuadrilateralRule(2, IntegrationRule::IVT_Load, GNS_NonlinearLargeStrain, 2);

	// fourth order interpolation - 9 nodes quadrilateral
	// geometric linear case
	AddQuadrilateralRule(4, IntegrationRule::IVT_Stiffness, GNS_Linear, 4);
	AddQuadrilateralRule(4, IntegrationRule::IVT_Mass, GNS_Linear, 3);
	AddQuadrilateralRule(4, IntegrationRule::IVT_Load, GNS_Linear, 2);
	// geometric nonlinear, small strain - same rule orders
	AddQuadrilateralRule(4, IntegrationRule::IVT_Stiffness, GNS_NonlinearSmallStrain, 4);
	AddQuadrilateralRule(4, IntegrationRule::IVT_Mass, GNS_NonlinearSmallStrain, 3);
	AddQuadrilateralRule(4, IntegrationRule::IVT_Load, GNS_NonlinearSmallStrain, 2);
	// geometric nonlinear, large strain
	AddQuadrilateralRule(4, IntegrationRule::IVT_Stiffness, GNS_NonlinearLargeStrain, 5);
	AddQuadrilateralRule(4, IntegrationRule::IVT_Mass, GNS_NonlinearLargeStrain, 4);
	AddQuadrilateralRule(4, IntegrationRule::IVT_Load, GNS_NonlinearLargeStrain, 2);
}

void IntegrationRulesCreator::AddQuadrilateralRule(int interpolationOrder,
	IntegrationRule::IntegratedValueType ivt, GeometricNonlinearityStatus gns, int ruleOrder)
{
	CreateQuadrilateralRule(ruleOrder);
	currentIntegrationRuleInfo->settings = IntegrationRule::IntegrationRuleSettings(
		TFE_Quadrilateral,interpolationOrder,ivt,gns,false
		);
	AddCurrentRuleInfo();
}

void IntegrationRulesCreator::CreateQuadrilateralRule(int order)
{
	// we use a ready function, which provides a one-dimensional integration rule
	Vector x;
	Vector w;
	GetIntegrationRule(x,w,order);

	// and now we create a 2D integration rule
	currentIntegrationRuleInfo = new IntegrationRule::IntegrationRuleInfo();
	currentIntegrationRuleInfo->integrationRule.SetLen(x.Length() * x.Length());
	int cnt = 1;
	for (int i1 = 1; i1 <= x.GetLen(); i1++)
	{
		for (int i2 = 1; i2 <= x.GetLen(); i2++)
		{
				currentIntegrationRuleInfo->integrationRule(cnt).x = x(i1);
				currentIntegrationRuleInfo->integrationRule(cnt).y = x(i2);
				currentIntegrationRuleInfo->integrationRule(cnt).z = 0;
				currentIntegrationRuleInfo->integrationRule(cnt).weight = w(i1) * w(i2);
				cnt++;
		}
	}
}

#pragma endregion

#pragma region Thin plate

void IntegrationRulesCreator::AddThinPlateRules()
{
	// fourth order interpolation - 4 nodes thin plate
	// geometric linear case
	AddThinPlateRule(4, IntegrationRule::IVT_Stiffness, GNS_Linear, 6);
	AddThinPlateRule(4, IntegrationRule::IVT_Mass, GNS_Linear, 4);
	AddThinPlateRule(4, IntegrationRule::IVT_Load, GNS_Linear, 4);
	// geometric nonlinear, small strain - higher rule order for the stiffness matrix
	AddThinPlateRule(4, IntegrationRule::IVT_Stiffness, GNS_NonlinearSmallStrain, 6);
	AddThinPlateRule(4, IntegrationRule::IVT_Mass, GNS_NonlinearSmallStrain, 4);
	AddThinPlateRule(4, IntegrationRule::IVT_Load, GNS_NonlinearSmallStrain, 4);
	// geometric nonlinear, large strain - higher rule order for the stiffness matrix
	AddThinPlateRule(4, IntegrationRule::IVT_Stiffness, GNS_NonlinearLargeStrain, 8);
	AddThinPlateRule(4, IntegrationRule::IVT_Mass, GNS_NonlinearLargeStrain, 4);
	AddThinPlateRule(4, IntegrationRule::IVT_Load, GNS_NonlinearLargeStrain, 4);
}

void IntegrationRulesCreator::AddThinPlateRule(int interpolationOrder,
	IntegrationRule::IntegratedValueType ivt, GeometricNonlinearityStatus gns, int ruleOrder)
{
	CreateQuadrilateralRule(ruleOrder);
	currentIntegrationRuleInfo->settings = IntegrationRule::IntegrationRuleSettings(
		TFE_ThinPlate,interpolationOrder,ivt,gns,false
		);
	AddCurrentRuleInfo();
}

#pragma endregion

#endif