//#**************************************************************
//#
//# filename:             FE2DTriQuad.cpp
//#
//# author:               Gerstmayr Johannes, YV
//#
//# generated:						December 2010
//# description:          2D particular finite elements
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
#include "Material.h"
#include "solversettings_auto.h"
#include "femathhelperfunctions.h"
#include "fe2dtriquad.h"

void Quadrilateral::SetQuadrilateral(int bodyindi, const TArray<int>& nodelist, int material_num, double thickness, const Vector3D& coli)
{
	assert(nodelist.Length() == 4 || nodelist.Length() == 9);		// Quadrilateral, number of nodes in the set function is invalid
	FiniteElement2D::SetFiniteElement2D(bodyindi, nodelist, material_num, thickness, coli);
	SetGeometricNonlinearityStatus(GNS_Linear);
}

double Quadrilateral::GetS0(const Vector2D& ploc, int shape) const
{
	// nodes are always sorted counterclock-wise in x-y (r-s) plane
	double r = ploc.X();
	double s = ploc.Y();

	if(NNodes() == 4)
	{
		switch(shape)
		{
		case 1:	return 0.25*(1+r)*(1+s);
		case 2:	return 0.25*(1-r)*(1+s);
		case 3:	return 0.25*(1-r)*(1-s);
		case 4:	return 0.25*(1+r)*(1-s);
		}
	}
	else
	{
		double r2 = r*r;
		double s2 = s*s;

		switch(shape)
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
		}
	}
	assert(0);
	return 0;
}

void Quadrilateral::GetDSMatrix0(const Vector2D& ploc, Matrix& sf) const
{
	double r = ploc.X();
	double s = ploc.Y();

	if(NNodes() == 4)
	{
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
	else
	{
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
}

Vector2D Quadrilateral::GetNodeLocPos2D(int i) const
{
	switch(i)
	{ //1..4 for linear element, 1..9 for quadratic element
	case 1: return Vector2D( 1., 1.); break;
	case 2: return Vector2D(-1., 1.); break;
	case 3: return Vector2D(-1.,-1.); break;
	case 4: return Vector2D( 1.,-1.); break;
	case 5: return Vector2D( 0., 1.); break;
	case 6: return Vector2D(-1., 0.); break;
	case 7: return Vector2D( 0.,-1.); break;
	case 8: return Vector2D( 1., 0.); break;
	case 9: return Vector2D( 0., 0.); break;
	default: return Vector2D(0,0);
	}
}

int Quadrilateral::GetActualInterpolationOrder() const
{
	switch(NNodes())
	{
	case 4: return 2;
	case 9: return 4;
	}
	assert(0);
	return 0;
}

void Quadrilateral::DefineIntegrationRule(IntegrationRule & integrationRule)
{
	assert(integrationRule.settings.elementType == TFE_Quadrilateral);

	int ruleOrder = 0;
	// here the particular rule order will be chosen depending on the settings
	if(
		integrationRule.settings.interpolationOrder == 2 ||
		integrationRule.settings.integratedValueType == IntegrationRule::IVT_Load
		)
		ruleOrder = 2;
	else
	{
		if(integrationRule.settings.integratedValueType == IntegrationRule::IVT_Stiffness)
			ruleOrder = 4;
		else
			ruleOrder = 3;
		if(integrationRule.settings.geometricNonlinearityStatus == GNS_NonlinearLargeStrain)
			ruleOrder++;
	}
	assert(ruleOrder != 0);

	IntegrationRule::DefineIntegrationRuleSquare(integrationRule, ruleOrder);
}