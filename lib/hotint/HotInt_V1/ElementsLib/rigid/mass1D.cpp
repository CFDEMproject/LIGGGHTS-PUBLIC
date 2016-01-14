//#**************************************************************
//#
//# filename:             Mass1D.cpp           
//#
//# author:               Daniel Reischl
//#
//# generated:						June 2013
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
//#***************************************************************************************
 
#include "mass1d.h"
#include "material.h"
#include "Node.h"
#include "rigid3D.h"
#include "femathhelperfunctions.h"
//#include "graphicsconstants.h"
#include "elementdataaccess.h"
#include "solversettings_auto.h"
#include "options_class_auto.h"
#include "rendercontext.h"


int Mass1D::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = Element::CheckConsistency(errorstr);
	if (rv) return rv;

	if(!mass)
	{
		rv = 1;
		errorstr = mystr("No mass is defined for the Mass1D!\n");
	}
	return rv;
}

//this function assigns default values to the element variables
void Mass1D::ElementDefaultConstructorInitialization()
{
	type = TBody;
	drawres = 8;
	x_init = Vector(2*SOS()); //zero initialized
	elementname = GetElementSpec();
	radius = 0.1;						// DR just some default value
	pref3D = Vector3D(0.);
	rotref3D  = Matrix3D(1.);
}

void Mass1D::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Element::GetElementData(edc);
}

int Mass1D::SetElementData(ElementDataContainer& edc) //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
{
	int rv = Element::SetElementData(edc);
	return rv;
}


//---------------------------------------------------------------
//for visualization of displacements and velocities:
//---------------------------------------------------------------
double Mass1D::GetFieldVariableValue(const FieldVariableDescriptor & fvd, const double & local_position, bool flagD)
{		
	double p0=local_position;

	if (flagD)
	{
		switch (fvd.VariableType())
		{
		case FieldVariableDescriptor::FVT_position:
			return GetPos1DD(p0);
		case FieldVariableDescriptor::FVT_velocity:
			return GetVel1DD(p0);
		case FieldVariableDescriptor::FVT_displacement:
			return GetPos1DD(p0) - x_init(1);
		}
	}
	else
	{
		switch (fvd.VariableType())
		{
		case FieldVariableDescriptor::FVT_position:
			return GetPos1D(p0);
		case FieldVariableDescriptor::FVT_velocity:
			return GetVel1D(p0);
		case FieldVariableDescriptor::FVT_displacement:
			return GetPos1D(p0) - x_init(1);
		case FieldVariableDescriptor::FVT_acceleration:	//$ DR 2013-07-01 added
			return XGPP(1); // $ MSax 2013-07-16 : changed from GetMBS()->GetAcceleration(LTG(1 + SOS())) to XGPP(1)
		}
	}
	return FIELD_VARIABLE_NO_VALUE;
}


void Mass1D::GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables)
{
	// add all the field variables of the parent class
	Element::GetAvailableFieldVariables(variables);
	//FVT_position,FVT_displacement, FVT_velocity: implemented in Element
	// possibility to remove entries does not exist yet, just do not use the Function of the parent class

	FieldVariableDescriptor::FieldVariableComponentIndex max_component = FieldVariableDescriptor::FieldVariableComponentIndex::FVCI_x;
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_acceleration, max_component); //$ DR 2013-07-01 added

}

void Mass1D::DrawElement()
{
	Element::DrawElement();

	mbs->SetColor(col);
	
	double fac = mbs->GetOptions()->PostProcOptions()->BodiesDeformationScaleFactor();
	if (!GetMBS()->GetIOption(151)) fac = 1;

	double p = GetRefPos1DD(); 

	if (GetMBS()->GetActualPostProcessingFieldVariable() != NULL)
	{
		double v = GetFieldVariableValue(*GetMBS()->GetActualPostProcessingFieldVariable(), 0., true);
		mbs->DrawColorSphere(v, ToP3D(p), radius, drawres);
	}
	else
	{
		mbs->DrawSphere(ToP3D(p), radius, drawres);
	}
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Rotor1D Rotor1D	Rotor1D	Rotor1D	Rotor1D	Rotor1D	Rotor1D	Rotor1D	Rotor1D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int Rotor1D::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = Element::CheckConsistency(errorstr);
	if (rv) return rv;

	if(!mass)
	{
		rv = 1;
		errorstr = mystr("No moment of inertia is defined for the Rotor1D!\n");
	}
	return rv;
}

//this function assigns default values to the element variables
void Rotor1D::ElementDefaultConstructorInitialization()
{
	type = TBody;
	drawres = 8;
	x_init = Vector(2*SOS()); //zero initialized
	elementname = GetElementSpec();
	radius = 0.1;						// DR just some default value
	length = 0.2;						// DR just some default value
	pref3D = Vector3D(0.);
	rotref3D  = Matrix3D(1.);
}


//---------------------------------------------------------------
//for visualization of displacements and velocities:
//---------------------------------------------------------------
double Rotor1D::GetFieldVariableValue(const FieldVariableDescriptor & fvd, const double & local_position, bool flagD)
{		
	double p0=local_position;

	if (flagD)
	{
		switch (fvd.VariableType())
		{
		case FieldVariableDescriptor::FVT_bryant_angle:
			return GetPos1DD(p0);
		case FieldVariableDescriptor::FVT_angular_velocity:
			return GetVel1DD(p0);
		}
	}
	else
	{
		switch (fvd.VariableType())
		{
		case FieldVariableDescriptor::FVT_bryant_angle:
			return GetPos1D(p0);
		case FieldVariableDescriptor::FVT_angular_velocity:
			return GetVel1D(p0);
		case FieldVariableDescriptor::FVT_angular_acceleration:
			return XGPP(1);
		}
	}
	return FIELD_VARIABLE_NO_VALUE;
}


void Rotor1D::GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables)
{
	FieldVariableDescriptor::FieldVariableComponentIndex max_component = FieldVariableDescriptor::FieldVariableComponentIndex::FVCI_x;

	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_bryant_angle, max_component);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_angular_velocity, max_component);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_angular_acceleration, max_component);
}


void Rotor1D::DrawElement()
{
	Element::DrawElement();
	mbs->SetColor(col);

	if(length<=0)
	{
		mbs->DrawArrow(GetRefPos(),GetPosD(Vector3D(0,radius,0)));
	}
	else
	{
		Vector3D pos1 = GetRefPos();
		Vector3D pos2 = GetRefPos() + rotref3D * Vector3D(length,0,0);
		mbs->DrawZyl(pos1, pos2, radius,64);

		Vector3D v1,v2,v3;
		v1 = pos2-pos1;																//axis
		v2 = GetPosD(Vector3D(0,radius,0)) - pos1;		// radius
		v3 = v1.Cross(v2);
		double l = 0.1*radius/(v3.Norm());
		v3 = l*v3;
		mbs->SetColor(Vector3D(0.,0.,0.));
		pos1 = pos1 - 0.025*v1;
		mbs->DrawCube(pos1,1.05*v1,v2,v3); 
	}
}