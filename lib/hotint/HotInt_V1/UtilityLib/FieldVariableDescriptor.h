//#**************************************************************
//#
//# filename:             FieldVariableDescriptor.h
//#
//# author:               Yury Vetyukov
//#
//# generated:						November 2010
//# description:          FieldVariableDescriptor
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
//#**************************************************************

#pragma once

#include "mystring.h"
#include "tarray.h"

class Vector3D;
class Vector2D;
class Matrix3D;
class Matrix2D;

// functions, which compute values of field variables,
// should return this constant to indicate that the value could not actuallu be computed;
// in the future this may affect processing of the values
const double FIELD_VARIABLE_NO_VALUE = 0;

// - this structure identifies a single post-processing field value for contour plotting and sensors.
// - elements (mainly finite elements) provide lists of relevant FieldVariableDescriptors
//   depending on their settings (i.e. inelastic_strain will appear only
//   when a corresponding material is selected). this affects user interface
//   (only the relevant options will appear in the selection).
// - while drawing, the elements have access to the actually selected FieldVariableDescriptor,
//   and plot the corresponding variable.
// - each value is defined by its type and, possibly, one or two identifiers of the components.
// - there exists a possibility to define a non-standard variable for post-processing.
struct FieldVariableDescriptor
{
public:
	// these are the most common (standard) variable types for post-processing;
	// adding a new standard type should be completed by adding its textual identifier in the constructor
	enum FieldVariableType
	{
		// kinematic value types
		FVT_displacement = 0,						// vector, 3 components in the reference basis
		FVT_position,								// vector, 3 components in the reference basis, not for plotting
		FVT_velocity,								// vector, 3 components in the reference basis
		FVT_acceleration,						// vector, 3 components in the reference basis
		// kinematic values in the local basis
		FVT_displacement_local_basis,						// vector, 3 components in the local basis
		FVT_velocity_local_basis,								// vector, 3 components in the local basis
		FVT_acceleration_local_basis,						// vector, 3 components in the local basis

		FVT_bryant_angle, 							//$ DR 2012-12-12 added according to JG (rotation x-y-z) //$ DR 2012-12-14 renamed
		// r = A*r_bar with A = A3*A2*A1 and Ai is the rotation matrix around axis i

		FVT_angular_velocity,							//$ DR 2012-12-12 added according to JG
		FVT_angular_velocity_local_basis, //$ DR 2012-12-12 added according to JG

		FVT_angular_acceleration,		//$ SW 2013-10-21: added

		// density
		FVT_density,								// scalar       //PG newly added
		FVT_pressure,								// scalar       //PG newly added

		// deformations
		FVT_total_strain,						// tensor, 6 symmetric components in the reference basis

		FVT_inelastic_strain,				// tensor, 6 symmetric components in the reference basis
		FVT_hardening_parameter,		// scalar value, (linear isotropic) hardening law parameter
		FVT_yield_function,					// scalar value, always <= 0, defines plastic flow law
		FVT_initial_strain,					// tensor, 6 symmetric components in the reference basis

		FVT_beam_curvature_Y,					// Y component in the local (rotating) basis
		FVT_beam_curvature_Z,					// Z component in the local (rotating) basis
		FVT_beam_torsion,						// scalar value

		FVT_beam_shear_Y,							// Y component in the local (rotating) basis
		FVT_beam_shear_Z,							// Z component in the local (rotating) basis
		FVT_beam_axial_extension,		// scalar value

		//+++++++++++++++++++++++++++++++++++++++++++
		//depricated:
		FVT_beam_curvature,					// component in the local (rotating) basis
		FVT_beam_shear,							// component in the local (rotating) basis
		FVT_beam_moment_bending,	  // component in the local (rotating) basis
		FVT_beam_force_transversal,	// component in the local (rotating) basis

		//+++++++++++++++++++++++++++++++++++++++++++

		// stresses in 3D
		FVT_stress,									// tensor, 6 symmetric components in the reference basis
		FVT_stress_mises,						// scalar value - synonim for stress magnitude
		// stress resultants for a 1D continuum (beam, rod)
		FVT_beam_force,							// vector, 3 components in the reference basis
		FVT_beam_force_axial,				// scalar value
		FVT_beam_force_transversal_Y,	// Y component in the local (rotating) basis
		FVT_beam_force_transversal_Z,	// Z component in the local (rotating) basis
		FVT_beam_moment,						// vector, 3 components in the reference basis
		FVT_beam_moment_torsional,	// scalar value
		FVT_beam_moment_bending_Y,	   // Y component in the local (rotating) basis
		FVT_beam_moment_bending_Z,	   // Z component in the local (rotating) basis
		// stress resultants for a 2D continuum (plate, shell)
		FVT_shell_force_in_plane,		// in-plane tensor with components xx, xy, yy in the local 
		FVT_shell_force_transversal,	// in-plane vector with 2 components in the local basis
		FVT_shell_moment,						// in-plane tensor with components xx, xy, yy in the local basis
		// "user-defined" value type for particular problems (contact stress, etc.)
		FVT_problem_specific
	};
	// identifier of the components of a variable
	enum FieldVariableComponentIndex
	{
		FVCI_none = 0,	// this component is not relevant
		FVCI_x = 1,
		FVCI_y = 2,
		FVCI_z = 3,
		FVCI_magnitude = 4		// norm of a vector or a tensor (root of the inner product with itself), should not be used for the second component
	};
	// possible dimensions of the variables - will be used in the plot legend
	enum FieldVariableDimensionality
	{
		FVD_none,				// non-dimensional, for strains
		FVD_length,			// for displacements etc., length
		FVD_density,
		FVD_pressure,
		FVD_velocity,
		FVD_acceleration,
		FVD_force,
		FVD_force_per_length,		// force / length
		FVD_force_per_length_square,		// for stresses etc., force / length / length
		FVD_force_length,			// moment, force*length
		FVD_1_per_length		// curvature etc., 1 / length
	};

private:
	// actual data
	FieldVariableType variable_type;
	FieldVariableComponentIndex component_index_1;
	FieldVariableComponentIndex component_index_2;
	bool is_not_for_plotting;	// this variable should not be available for the contour plotting
	//$EK 2012-11-07 dangerous use due to wrong use of copy constructer (only pointer is set -> also in operator=)
	const char * problem_specific_textual_identifier_without_components;

public:
	// creation of a standard value info
	FieldVariableDescriptor(
		FieldVariableType variable_type,
		FieldVariableComponentIndex component_index_1 = FVCI_none,
		FieldVariableComponentIndex component_index_2 = FVCI_none
		);
	// creation of a problem specific value info;
	// the textual identifier will be saved; the string should not contain the components;
	// the string pointer will be saved, so that it must be constant
	// and valid during the life cycle of this FieldVariableDescriptor
	FieldVariableDescriptor(
		const char * textual_identifier,
		FieldVariableComponentIndex component_index_1 = FVCI_none,
		FieldVariableComponentIndex component_index_2 = FVCI_none
		);
	// copy constructor
	FieldVariableDescriptor(const FieldVariableDescriptor & fvd);

	// when this function is called with true, the variable is marked to be unavailable for contout plotting
	void SetNotForPlotting(bool is_not_for_plotting_ = true) { is_not_for_plotting = is_not_for_plotting_; }

	// actual variable type
	FieldVariableType VariableType() const { return variable_type; }
	// actual components
	FieldVariableComponentIndex ComponentIndex1() const { return component_index_1; }
	FieldVariableComponentIndex ComponentIndex2() const { return component_index_2; }
	// whether this variable should not be available for the contour plotting
	bool IsNotForPlotting() const { return is_not_for_plotting; }

	// these functions provide a short informative text description
	// for the user interface - with and without the notion of components
	mystr GetTextualIdentifier() const;		// full name - with components
	const char * GetTextualIdentifierWithoutComponents() const;
	const char * GetTextualIdentifierComponentsOnly() const;
	void InitializeFromTextualIdentifier(mystr & textualIdentifier);		// restores the state from the textual identifier

	// dimensionality of the variable - evaluated depending on the FieldVariableType;
	// problem specific types are not supported in the present implementation
	virtual FieldVariableDimensionality GetDimensionality() const;

private:
	// comparison for sorting to the desired order,
	// in which variables should appear in the user interface.
	// the functions returns -1 if "this" variable should go first, 0 if the variables are the same and +1 otherwise.
	int Compare(const FieldVariableDescriptor & fvd_another) const;

public:
	bool operator==(const FieldVariableDescriptor & fvd) const
	{
		int rv = Compare(fvd);
		if(rv==0)
		{
			return true;
		}
		return false;
		//return
		//	variable_type == FVD.variable_type &&
		//	component_index_1 == FVD.component_index_1 && 
		//	component_index_2 == FVD.component_index_2; 
	}

public:
	// helper functions, which should be used to simplify the usage of this class:

	// the contents of the array source will be merged to the contents of the array dest such,
	// that dest will contain only the unique elements in the order, determined by the Compare() function.
	static void MergeArrays(TArray<FieldVariableDescriptor> & dest, const TArray<FieldVariableDescriptor> & source);
	// adds to the array family of variable descriptors with the same main type;
	// when components are relevant, then the magnitude of the variable will be added.
	static void AddTypeIntoArray(
		TArray<FieldVariableDescriptor> & dest,														// array
		FieldVariableType variable_type,																	// type of the added variables
		FieldVariableComponentIndex component_index_range = FVCI_none,		// range of the components (if relevant)
		bool two_components_symbol = false,																// true - tensor (2-indexed-symbol), false - vector (1-indexed-symbol)
		bool symmetric = true																							// for tensors (2-indexed-symbols): if true, only the symmetric components will be added
		);
	// the same as before for problem specific types
	static void AddTypeIntoArray(
		TArray<FieldVariableDescriptor> & dest,														// array
		//$EK changed to const 
		const char * textual_identifier,																				// string of text, which identifies the problem specific type
		FieldVariableComponentIndex component_index_range = FVCI_none,		// range of the components (if relevant)
		bool two_components_symbol = false,																// true - tensor (2-indexed-symbol), false - vector (1-indexed-symbol)
		bool symmetric = true																							// for tensors (2-indexed-symbols): if true, only the symmetric components will be added
		);
	// looks, if a type of variable descriptor is found in the array;
	// returns position of first occurence in the array, and returns 0 if not found.
	static int FindTypeInArray(
		TArray<FieldVariableDescriptor> & fvd_array,											// array
		FieldVariableType variable_type,																	// type of the added variables
		FieldVariableComponentIndex component_index_1,										// component 1
		FieldVariableComponentIndex component_index_2 = FVCI_none					// component 2 (if relevant)
		);

	// the following functions compute corresponding components or magnitudes of a vector or a matrix
	double GetComponent(const Vector3D & v) const;
	double GetComponent(const Vector2D & v) const;
	double GetComponent(const Matrix3D & m) const;
	double GetComponent(const Matrix2D & m) const;

private:
	// initializer
	void Initialize(FieldVariableType variable_type, FieldVariableComponentIndex component_index_1, FieldVariableComponentIndex component_index_2);

public:
	void operator=(const FieldVariableDescriptor & fvd);
	FieldVariableDescriptor() {}

	static char GetComponentsDelimiter();
};