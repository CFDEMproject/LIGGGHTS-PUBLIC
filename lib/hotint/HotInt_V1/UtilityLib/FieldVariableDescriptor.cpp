//#**************************************************************
//#
//# filename:             FieldVariableDescriptor.cpp
//#
//# author:               YV
//#
//# generated:						November 2010
//# description:          
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
//#**************************************************************


#include "mbs_interface.h"
#include "FieldVariableDescriptor.h"

FieldVariableDescriptor::FieldVariableDescriptor(
		FieldVariableType variable_type,
		FieldVariableComponentIndex component_index_1,
		FieldVariableComponentIndex component_index_2
		)
{
	Initialize(variable_type,component_index_1,component_index_2);
}

FieldVariableDescriptor::FieldVariableDescriptor(
		const char * textual_identifier,
		FieldVariableComponentIndex component_index_1,
		FieldVariableComponentIndex component_index_2
		)
{
	Initialize(FVT_problem_specific,component_index_1,component_index_2);
	problem_specific_textual_identifier_without_components = textual_identifier;
}

void FieldVariableDescriptor::Initialize(FieldVariableType variable_type, FieldVariableComponentIndex component_index_1, FieldVariableComponentIndex component_index_2)
{
	this->component_index_1 = component_index_1;
	this->component_index_2 = component_index_2;
	this->variable_type = variable_type;

	switch(variable_type)
	{
	case FVT_position: SetNotForPlotting(true); break;
	default: SetNotForPlotting(false);
	}
	problem_specific_textual_identifier_without_components = GetTextualIdentifierWithoutComponents();
}

FieldVariableDescriptor::FieldVariableDescriptor(const FieldVariableDescriptor & fvd)
{
	//$EK 2012-11-07 dangerous, since this means a set of a pointer with a pointer
	*this = fvd;
}

// this is a helper function, which is actually doing the work for standard types
const char * GetTextualIdentifierWithoutComponents(FieldVariableDescriptor::FieldVariableType vt)
{
	// no spaces are allowed
	switch(vt)
	{
	case FieldVariableDescriptor::FVT_displacement:						return "displacement";
	case FieldVariableDescriptor::FVT_position:								return "position";
	case FieldVariableDescriptor::FVT_velocity:								return "velocity";
	case FieldVariableDescriptor::FVT_acceleration:						return "acceleration";
	case FieldVariableDescriptor::FVT_displacement_local_basis:	return "displacement_local_basis";
	case FieldVariableDescriptor::FVT_velocity_local_basis:		return "velocity_local_basis";
	case FieldVariableDescriptor::FVT_acceleration_local_basis:	return "acceleration_local_basis";
	case FieldVariableDescriptor::FVT_bryant_angle:					return "bryant_angle";	//$ DR 2012-12-12 added according to JG
	case FieldVariableDescriptor::FVT_angular_velocity:				return "angular_velocity"; //$ DR 2012-12-12 added according to JG
	case FieldVariableDescriptor::FVT_angular_velocity_local_basis:	return "angular_velocity_local_basis"; //$ DR 2012-12-12 added according to JG
	case FieldVariableDescriptor::FVT_angular_acceleration:		return "angular_acceleration";//$ SW 2013-10-21: added
	case FieldVariableDescriptor::FVT_density:									return "density"; //PG newly added
	case FieldVariableDescriptor::FVT_pressure:								return "pressure"; //PG newly added
	case FieldVariableDescriptor::FVT_total_strain:						return "total_strain";
	case FieldVariableDescriptor::FVT_inelastic_strain:				return "inelastic_strain";
	case FieldVariableDescriptor::FVT_hardening_parameter:			return "hardening_parameter";
	case FieldVariableDescriptor::FVT_yield_function:					return "yield_function";
//	case FieldVariableDescriptor::FVT_plastic_strain:					return "plastic strain";
	case FieldVariableDescriptor::FVT_initial_strain:					return "initial_strain";

	case FieldVariableDescriptor::FVT_beam_torsion:						return "beam_torsion";
	case FieldVariableDescriptor::FVT_beam_shear:							return "beam_shear";


	//+++++++++++++++++++++++++++++++++++++++++++
	//depricated, use single values Y/Z; JG2013-09
	case FieldVariableDescriptor::FVT_beam_curvature:					return "beam_curvature";
	case FieldVariableDescriptor::FVT_beam_axial_extension:		return "beam_axial_extension";
	case FieldVariableDescriptor::FVT_beam_moment_bending:			return "beam_moment_bending";
	case FieldVariableDescriptor::FVT_beam_force_transversal:	return "beam_force_transversal";
	//+++++++++++++++++++++++++++++++++++++++++++

	//new ones
	case FieldVariableDescriptor::FVT_beam_curvature_Y:					return "beam_curvature_Y";
	case FieldVariableDescriptor::FVT_beam_shear_Y:							return "beam_shear_Y";
	case FieldVariableDescriptor::FVT_beam_moment_bending_Y:			return "beam_moment_bending_Y";
	case FieldVariableDescriptor::FVT_beam_force_transversal_Y:	return "beam_force_transversal_Y";

	case FieldVariableDescriptor::FVT_beam_curvature_Z:					return "beam_curvature_Z";
	case FieldVariableDescriptor::FVT_beam_shear_Z:							return "beam_shear_Z";
	case FieldVariableDescriptor::FVT_beam_moment_bending_Z:			return "beam_moment_bending_Z";
	case FieldVariableDescriptor::FVT_beam_force_transversal_Z:	return "beam_force_transversal_Z";

	case FieldVariableDescriptor::FVT_stress:									return "stress";
	case FieldVariableDescriptor::FVT_stress_mises:						return "stress_mises";
	case FieldVariableDescriptor::FVT_beam_force:							return "beam_force";
	case FieldVariableDescriptor::FVT_beam_force_axial:				return "beam_force_axial";
	case FieldVariableDescriptor::FVT_beam_moment:							return "beam_moment";
	case FieldVariableDescriptor::FVT_beam_moment_torsional:		return "beam_moment_torsional";
	case FieldVariableDescriptor::FVT_shell_force_in_plane:		return "shell_force_in_plane";
	case FieldVariableDescriptor::FVT_shell_force_transversal:	return "shell_force_transversal";
	case FieldVariableDescriptor::FVT_shell_moment:						return "shell_moment";
	}
	return NULL;
}

const char * FieldVariableDescriptor::GetTextualIdentifierWithoutComponents() const
{
	const char * textualIdentifier = ::GetTextualIdentifierWithoutComponents(variable_type);
	if(textualIdentifier != NULL)
		return textualIdentifier;
	if(variable_type == FVT_problem_specific)
		return problem_specific_textual_identifier_without_components;
	assert(0);
	return "unknown_variable_type";
}

static const char * component_names1[] = {"x", "y", "z", "magnitude"};
static const char * component_names21[] = {"xx", "xy", "xz"};
static const char * component_names22[] = {"yx", "yy", "yz"};
static const char * component_names23[] = {"zx", "zy", "zz"};
static const char ** component_names2[] = {component_names21, component_names22, component_names23};
static const char componentsDelimiter = '.';

char FieldVariableDescriptor::GetComponentsDelimiter(){return componentsDelimiter;}

const char * FieldVariableDescriptor::GetTextualIdentifierComponentsOnly() const
{
	if(component_index_1 == FVCI_none)
		return "";
	if(component_index_2 == FVCI_none || component_index_1 == FVCI_magnitude)
		return component_names1[component_index_1 - 1];
	return component_names2[component_index_1 - 1][component_index_2 - 1];
}

mystr FieldVariableDescriptor::GetTextualIdentifier() const
{
	return mystr(GetTextualIdentifierWithoutComponents()) + mystr(componentsDelimiter) + GetTextualIdentifierComponentsOnly();
}

void FieldVariableDescriptor::InitializeFromTextualIdentifier(mystr & textualIdentifier)
{
	int delimiterPos = textualIdentifier.Find(componentsDelimiter);
	mystr name = delimiterPos == -1 ? textualIdentifier : textualIdentifier.SubString(0, delimiterPos - 1);
	variable_type = FVT_problem_specific;
	for(int vt = 0; vt < FVT_problem_specific; vt++)
		if(name.Compare(::GetTextualIdentifierWithoutComponents((FieldVariableType)vt)))
		{
			variable_type = (FieldVariableType)vt;
			break;
		}
	if(variable_type == FVT_problem_specific)
		problem_specific_textual_identifier_without_components = "unknown_variable_type";

	mystr components = delimiterPos == -1 ? "" : textualIdentifier.SubString(delimiterPos + 1, textualIdentifier.Length() - 1);
	component_index_1 = FVCI_none;
	component_index_2 = FVCI_none;
	if(components.Compare("magnitude"))
		component_index_1 = FVCI_magnitude;
	else if(components.Compare("x"))
		component_index_1 = FVCI_x;
	else if(components.Compare("y"))
		component_index_1 = FVCI_y;
	else if(components.Compare("z"))
		component_index_1 = FVCI_z;
	else if(components.Compare("xx"))
	{
		component_index_1 = FVCI_x;
		component_index_2 = FVCI_x;
	}
	else if(components.Compare("xy"))
	{
		component_index_1 = FVCI_x;
		component_index_2 = FVCI_y;
	}
	else if(components.Compare("xz"))
	{
		component_index_1 = FVCI_x;
		component_index_2 = FVCI_z;
	}
	else if(components.Compare("yx"))
	{
		component_index_1 = FVCI_y;
		component_index_2 = FVCI_x;
	}
	else if(components.Compare("yy"))
	{
		component_index_1 = FVCI_y;
		component_index_2 = FVCI_y;
	}
	else if(components.Compare("yz"))
	{
		component_index_1 = FVCI_y;
		component_index_2 = FVCI_z;
	}
	else if(components.Compare("zx"))
	{
		component_index_1 = FVCI_z;
		component_index_2 = FVCI_x;
	}
	else if(components.Compare("zy"))
	{
		component_index_1 = FVCI_z;
		component_index_2 = FVCI_y;
	}
	else if(components.Compare("zz"))
	{
		component_index_1 = FVCI_z;
		component_index_2 = FVCI_z;
	}
}

int FieldVariableDescriptor::Compare(const FieldVariableDescriptor & fvd_another) const
{
	// first we compare by the variable type
	if(variable_type != fvd_another.variable_type)
		return variable_type < fvd_another.variable_type ? -1 : 1;
	// the variable types are the same, let's see if this is a problem specific variable
	if(variable_type == FVT_problem_specific)
	{
		int cmp = strcmp(problem_specific_textual_identifier_without_components, fvd_another.problem_specific_textual_identifier_without_components);
		if(cmp != 0)
			return cmp;
	}
	// now we compare the components
	if(component_index_1 != fvd_another.component_index_1)
		return component_index_1 < fvd_another.component_index_1 ? -1 : 1;
	if(component_index_2 != fvd_another.component_index_2)
		return component_index_2 < fvd_another.component_index_2 ? -1 : 1;
	return 0;			// the variable definitions are identical
}

void FieldVariableDescriptor::MergeArrays(TArray<FieldVariableDescriptor> & dest, const TArray<FieldVariableDescriptor> & source)
{
	// in case of a performance loss this function should be rewritten with dictionaries (maps) and hash codes
	for(int j = 1; j <= source.Length(); j++)
	{
		const FieldVariableDescriptor & fvd = source(j);
		int i;
		int cmp = 1;
		for(i = 1; i <= dest.Length(); i++)
		{
			cmp = fvd.Compare(dest(i));
			if(cmp != 1)
				break;	// found either the same entry or the point where this entry is to be inserted
		}
		if(cmp != 0)
			dest.Insert(i, fvd);
	}
}

void FieldVariableDescriptor::AddTypeIntoArray(
		TArray<FieldVariableDescriptor> & dest,
		FieldVariableType variable_type,
		FieldVariableComponentIndex component_index_range,
		bool two_components_symbol, bool symmetric
		)
{
	if(component_index_range == FVCI_none)
	{
		dest.Add(FieldVariableDescriptor(variable_type));
	}
	else
	{
		for(FieldVariableComponentIndex component1 = FVCI_x; component1 <= component_index_range; component1 = (FieldVariableComponentIndex)(component1 + 1))
		{
			if(two_components_symbol)
			{
				for(FieldVariableComponentIndex component2 = FVCI_x; component2 <= component_index_range; component2 = (FieldVariableComponentIndex)(component2 + 1))
				{
					if(!symmetric || component1 <= component2)
						dest.Add(FieldVariableDescriptor(variable_type, component1, component2));
				}
			}
			else
			{
				dest.Add(FieldVariableDescriptor(variable_type, component1));
			}
		}
		dest.Add(FieldVariableDescriptor(variable_type, FVCI_magnitude));
	}
}

void FieldVariableDescriptor::AddTypeIntoArray(
		TArray<FieldVariableDescriptor> & dest,
		const char * textual_identifier,
		FieldVariableComponentIndex component_index_range,
		bool two_components_symbol, bool symmetric
		)
{
	if(component_index_range == FVCI_none)
	{
		dest.Add(FieldVariableDescriptor(textual_identifier));
	}
	else
	{
		for(FieldVariableComponentIndex component1 = FVCI_x; component1 <= component_index_range; component1 = (FieldVariableComponentIndex)(component1 + 1))
		{
			if(two_components_symbol)
			{
				for(FieldVariableComponentIndex component2 = FVCI_x; component2 <= component_index_range; component2 = (FieldVariableComponentIndex)(component2 + 1))
				{
					if(!symmetric || component1 <= component2)
						dest.Add(FieldVariableDescriptor(textual_identifier, component1, component2));
				}
			}
			else
			{
				dest.Add(FieldVariableDescriptor(textual_identifier, component1));
			}
		}
		dest.Add(FieldVariableDescriptor(textual_identifier, FVCI_magnitude));
	}
}

int FieldVariableDescriptor::FindTypeInArray(
		TArray<FieldVariableDescriptor> & fvd_array,											// array
		FieldVariableType variable_type,																	// type of the added variables
		FieldVariableComponentIndex component_index_1,										// component 1
		FieldVariableComponentIndex component_index_2											// component 2 (if relevant)
		)
{
	int rv = 0;
	FieldVariableDescriptor fvd(variable_type, component_index_1, component_index_2);

	for(int i=1; i<=fvd_array.Length(); i++)
	{
		if (!fvd.Compare(fvd_array(i)))
		{
			return i;
		}
	}
	
	return rv;
}

void FieldVariableDescriptor::operator=(const FieldVariableDescriptor & fvd)
{
	Initialize(fvd.variable_type, fvd.component_index_1, fvd.component_index_2);
	is_not_for_plotting = fvd.is_not_for_plotting;
	//$EK 2012-11-07 bug fix: textual identifier has to be set (!! dangerous -> pointer = pointer)
	problem_specific_textual_identifier_without_components = fvd.problem_specific_textual_identifier_without_components;
}

double FieldVariableDescriptor::GetComponent(const Vector3D & v) const
{
	assert(component_index_1 != FVCI_none && component_index_2 == FVCI_none);
	if(component_index_1 == FVCI_magnitude)
		return v.Norm();
	return v(component_index_1);
}

double FieldVariableDescriptor::GetComponent(const Vector2D & v) const
{
	assert(component_index_1 != FVCI_none && component_index_2 == FVCI_none);
	if(component_index_1 == FVCI_magnitude)
		return v.Norm();
	if(component_index_1 == FVCI_z)
		return FIELD_VARIABLE_NO_VALUE;		// this may happen if there are 3D and 2D elements mixed in the system
	return v(component_index_1);
}

double FieldVariableDescriptor::GetComponent(const Matrix3D & m) const
{
	assert(component_index_1 != FVCI_none && component_index_2 != FVCI_magnitude);
	assert(component_index_2 != FVCI_none || component_index_1 == FVCI_magnitude);
	if(component_index_1 == FVCI_magnitude)
		return sqrt(m.InnerProduct(m));
	return m(component_index_1, component_index_2);
}

double FieldVariableDescriptor::GetComponent(const Matrix2D & m) const
{
	//EK 2012-03-02 assert bug fix - if component_index_1 == FVCI_magnitude, component_index_2 does not matter
	assert(component_index_1 != FVCI_none && component_index_2 != FVCI_magnitude);
	assert(component_index_2 != FVCI_none || component_index_1 == FVCI_magnitude);
	if(component_index_1 == FVCI_magnitude)
		return sqrt(m.InnerProduct(m));
	if(component_index_1 == FVCI_z || component_index_2 == FVCI_z)
		return FIELD_VARIABLE_NO_VALUE;		// this may happen if there are 3D and 2D elements mixed in the system
	return m(component_index_1, component_index_2);
}

FieldVariableDescriptor::FieldVariableDimensionality FieldVariableDescriptor::GetDimensionality() const
{
	switch(variable_type)
	{
	case FVT_displacement:
	case FVT_position:										return FVD_length;
	case FVT_velocity:										return FVD_velocity;
	case FVT_acceleration:								return FVD_acceleration;
	case FVT_density:											return FVD_density;
	case FVT_pressure:										return FVD_pressure;
	case FVT_total_strain:
	case FVT_inelastic_strain:
	case FVT_hardening_parameter:
	case FVT_initial_strain:

	case FVT_beam_axial_extension:
	case FVT_beam_shear:									return FVD_none;
	case FVT_beam_curvature:
	case FVT_beam_torsion:                return FVD_1_per_length;

	case FVT_beam_shear_Y:								return FVD_none;
	case FVT_beam_curvature_Y:							return FVD_1_per_length;
	case FVT_beam_force_transversal_Y:			return FVD_force;
	case FVT_beam_moment_bending_Y:					return FVD_force_length;

	case FVT_beam_shear_Z:								return FVD_none;
	case FVT_beam_curvature_Z:							return FVD_1_per_length;
	case FVT_beam_force_transversal_Z:			return FVD_force;
	case FVT_beam_moment_bending_Z:					return FVD_force_length;

	case FVT_stress:
	case FVT_stress_mises:
	case FVT_yield_function:							return FVD_force_per_length_square;
	case FVT_beam_force:
	case FVT_beam_force_axial:
	case FVT_beam_force_transversal:			return FVD_force;
	case FVT_beam_moment:
	case FVT_beam_moment_torsional:
	case FVT_beam_moment_bending:					return FVD_force_length;
	case FVT_shell_force_in_plane:
	case FVT_shell_force_transversal:			return FVD_force_per_length;
	case FVT_shell_moment:								return FVD_force;
		// still missing:
		//case FVT_bryant_angle:							
		//case FVT_angular_velocity:						
		//case FVT_angular_velocity_local_basis:
	}
	return FVD_none;
}