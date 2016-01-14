//#**************************************************************
//#
//# filename:             IntegrationRule.h
//#
//# author:               YV
//#
//# generated:						october 2010; refactored in august 2011
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
 
#pragma once

#include "finite_element_definitions.h"

/// This class encapsulates the concept of an integration rule for integration over 2D/3D finite elements.
/// Depending on
/// - exact element type (TFiniteElementType)
/// - interpolation order (quadrilinear, ...)
/// - entity type to be computed (mass matrix, stiffness matrix, H matrix, etc.)
/// - status of geometric nonlinearity
/// - possible user-defined flags (reduced integration, locking compensation, etc.)
/// a corresponding set of integration points and weights must be supplied by the library of integration rules.
/// Integration rules are stored separately from finite elements to optimize the memory usage.
/// "Browsing" or "iterating" the integration points in a loop should be done
/// with a helper class IntegrationPointsIterator.
/// Particular definitions of the integraion rules must be implemented in the client finite element classes.
class IntegrationRule
{
public:
#pragma region general definitions

	// data of a particular integraion point
	struct IntegrationPoint
	{
		double x;
		double y;
		double z;
		double weight;
	};
	/// Type of the entity to integrate.
	enum IntegratedValueType
	{
		IVT_Stiffness,		///< K matrix
		IVT_Mass,					///< M matrix
		IVT_Load					///< H matrix
	};
	/// "Input data", depending on which a particular integration rule will be provided.
	struct IntegrationRuleSettings
	{
		TFiniteElementType elementType;		// defined in finiteelementgeneric.h
		int interpolationOrder;						// 2 - quadrilinear, ... - order of the interpolating polinomials
		IntegratedValueType integratedValueType;
		GeometricNonlinearityStatus geometricNonlinearityStatus;	// defined in femesh.h
		short flags;

		/// default constructor, required by TArray
		IntegrationRuleSettings() {}
		/// full constructor
		IntegrationRuleSettings(
			TFiniteElementType elementType,
			int interpolationOrder,
			IntegratedValueType integratedValueType,
			GeometricNonlinearityStatus geometricNonlinearityStatus,
			short flags
			);
		// comparison of two settings
		bool operator==(const IntegrationRuleSettings & settings) const;
	};
	// each finite element must be able to create own integration rules, i.e. must implement this interface
	// finite elements may access the library of integration rules with given settings;
	// if there is no corresponding entry in the library, then the element is asked to create it by calling DefineIntegrationRule()
	struct IntegrationRuleProvider
	{
		// the element obtains a "half-finished product":
		// the settings of the desired integration rule are set,
		// and the integration points need to be computed
		virtual void DefineIntegrationRule(IntegrationRule & integrationRule) = 0;

		//$ YV 2013-01-12: the library of integration rules is common for all providers of integration rules
	private:
		static IntegrationRulesLibrary integrationRulesLibrary;
	protected:
		static IntegrationRulesLibrary * GetIntegrationRulesLibrary() { return &integrationRulesLibrary; }
	};

#pragma endregion

	// data of the class: an array of integration points and the particular settings
	TArray<IntegrationPoint> integrationPoints;
	IntegrationRuleSettings settings;

	int GetNumberOfIntegrationPoints() const { return integrationPoints.Length(); }

	// helper functions for the creation of "typical" integration rules
	static void DefineIntegrationRuleSquare(IntegrationRule & integrationRule, int ruleOrder);
	static void DefineIntegrationRuleLobattoSquare(IntegrationRule & integrationRule, int ruleOrder);
	static void DefineIntegrationRuleSquareAnisotropic(IntegrationRule & integrationRule, int ruleOrder_x, int ruleOrder_y);
	static void DefineIntegrationRuleTriangle(IntegrationRule & integrationRule, int ruleOrder);
	static void DefineIntegrationRuleCube(IntegrationRule & integrationRule, int ruleOrder);
	static void DefineIntegrationRuleTetrahedron(IntegrationRule & integrationRule, int ruleOrder);
	//$EK 2013-03-05 added int rule for prisms
	static void DefineIntegrationRulePrism(IntegrationRule & integrationRule, int ruleOrder);
	//$EK 2013-03-05 added int rule for pyramid
	static void DefineIntegrationRulePyramid(IntegrationRule & integrationRule, int ruleOrder);
	static void DefineIntegrationRuleLine(IntegrationRule & integrationRule, int ruleOrder); // AH: integration on a line
};

class IntegrationRulesLibrary
{
	TArray<IntegrationRule*> integrationRules;
public:
	~IntegrationRulesLibrary();
	// finds a suitable integration rule in the library and returns a pointer to it;
	// if the rule is not found, then the supplied IntegrationRuleProvider is invoked
	// to produce a new entry in the library, which is then returned;
	// if the rule is not found and no IntegrationRuleProvider is supplied, then NULL is returned.
	IntegrationRule * GetIntegrationRule(
		const IntegrationRule::IntegrationRuleSettings & settings,
		IntegrationRule::IntegrationRuleProvider * pIntegrationRuleProvider = NULL
		);
};

/// With this helper class the client code
/// can iterate over all integration points of an integration rule
class IntegrationPointsIterator
{
public:
	/// Creation of an iterator and setting it to the first integration point.
	IntegrationPointsIterator(const IntegrationRule * pIntegrationRule_) : pIntegrationRule(pIntegrationRule_), index(1) {}
	/// Moving the iterator to the next integration point.
	void operator++() { ++index; }
	/// Testing if we have iterated through all integration point.
	bool IsEnd() { return index > pIntegrationRule->integrationPoints.Length(); }
	/// Getting the actual integration point in 3D case.
	const Vector3D & Point() const;
	/// Getting the actual integration point in 2D case.
	const Vector2D & Point2D() const;
	/// Find the closest Integration Point to ploc and set this-pointer to it
	void GoClosestTo(const Vector3D& ploc);
	/// Getting the weight of the actual integration point.
	double Weight();
	/// Getting the 1-based index of the actual integration point -
	/// is often needed as certain data may be  stored in external arrays (inelastic properties, plastic strains, etc.)
	int GetIndex() const { return index; }
	/// In certain cituations it is needed to access elements directly by index,
	/// so that we allow explicit 1-based index setting
	void SetIndex(int index) { this->index = index; }


protected:
	const IntegrationRule * pIntegrationRule;
	int index;
};