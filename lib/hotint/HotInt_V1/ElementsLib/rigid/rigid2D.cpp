//#**************************************************************
//#
//# filename:             rigid2D.cpp
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						17.October 2004
//# description:          Driver and model for timeintegration
//#                       Model of a rigid arm and a hydraulic zylinder
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
#include "material.h"
#include "Node.h"
#include "rigid3D.h"
#include "femathhelperfunctions.h"
#include "graphicsconstants.h"
#include "elementdataaccess.h"
#include "solversettings_auto.h"
#include "options_class_auto.h"
#include "rendercontext.h"


void Body2D::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Element::GetElementData(edc);

	ElementData ed;

	//SetElemDataVector3D(edc, pref3D, "ReferencePoint3D"); edc.Last().SetToolTipText("Reference point for transformation of planar objects to 3D; p = [X, Y, Z]");

	//Vector v(9);
	//v(1) = rotref3D(1,1); v(2) = rotref3D(1,2); v(3) = rotref3D(1,3);
	//v(4) = rotref3D(2,1); v(5) = rotref3D(2,2); v(6) = rotref3D(2,3);
	//v(7) = rotref3D(3,1); v(8) = rotref3D(3,2); v(9) = rotref3D(3,3);
	//ed.SetVector(v.GetVecPtr(), v.Length(), "ReferenceRot3D");    ed.SetToolTipText("Rotation matrix for transformation of planar objects to 3D; [R11, R12,..., R33]"); edc.Add(ed);

	//Vector3D size must be set in derived class!!!	
}

int Body2D::SetElementData(ElementDataContainer& edc) //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
{
	int rv = Element::SetElementData(edc); 

	//GetElemDataVector3D(GetMBS(), edc, "ReferencePoint3D", pref3D);

	//Vector v(9);
	//GetElemDataVector(GetMBS(), edc, "ReferenceRot3D", v); 
	//if (v.Length() == 9) 
	//{
	//	rotref3D(1,1) = v(1); rotref3D(1,2) = v(2); rotref3D(1,3) = v(3);
	//	rotref3D(2,1) = v(4); rotref3D(2,2) = v(5); rotref3D(2,3) = v(6);
	//	rotref3D(3,1) = v(7); rotref3D(3,2) = v(8); rotref3D(3,3) = v(9);
	//}
	//else
	//{
	//	GetMBS()->LoadError("ReferenceRot3D must have 9 components!");
	//}

	return rv;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//$ DR 2012-07: CheckConsistency added
int Mass2D::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = Element::CheckConsistency(errorstr);
	if (rv) return rv;

	if(!mass)
	{
		rv = 1;
		errorstr = mystr("No mass is defined for the Mass2D!\n");
	}

	return rv;
}

//this function assigns default values to the element variables
void Mass2D::ElementDefaultConstructorInitialization()
{
	drawres = 8;
	x_init = Vector(2*SOS()); //zero initialized
	elementname = GetElementSpec();
	size.X() = 0.1;						// DR just some default value
}

void Mass2D::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Body2D::GetElementData(edc);

	ElementData ed;

	ElementData* edf = edc.TreeFind("Physics.density");
	if (edf) edf->SetLocked(1);

	//ed.SetDouble(mass, "mass"); edc.TreeAdd("Physics",ed);

	edf = edc.TreeFind("Geometry.size");
	if (edf) edf->SetLocked(1);

	//ed.SetDouble(size.X(), "radius"); edc.TreeAdd("Geometry",ed);

	//ed.SetVector2D(x_init(1), x_init(2), "initial_position"); edc.TreeAdd("Initialization",ed);
	//ed.SetVector2D(x_init(3), x_init(4), "initial_velocity"); edc.TreeAdd("Initialization",ed);

	//ed.SetInt(drawres, "Draw_resolution", 1, 1000); ed.SetToolTipText("Number of segments used to approximate circumference of circle"); edc.Add(ed);
}

int Mass2D::SetElementData(ElementDataContainer& edc) //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
{
	int rv = Body2D::SetElementData(edc);

	//GetElemDataDouble(GetMBS(), edc, "Physics.mass", mass);
	//GetElemDataDouble(GetMBS(), edc, "Geometry.radius", size.X());

	//Vector2D v;
	//GetElemDataVector2D(GetMBS(), edc, "Initialization.initial_position", v); x_init(1) = v.X(); x_init(2) = v.Y(); 
	//GetElemDataVector2D(GetMBS(), edc, "Initialization.initial_velocity", v); x_init(3) = v.X(); x_init(4) = v.Y(); 

	//GetElemDataInt(GetMBS(), edc, "Draw_resolution", drawres);

	return rv;
}


//---------------------------------------------------------------
//for visualization of displacements and velocities:
//---------------------------------------------------------------
double Mass2D::GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector2D & local_position, bool flagD)
{		
	Vector2D p0(local_position);

	if (flagD)
	{
		switch (fvd.VariableType())
		{
		case FieldVariableDescriptor::FVT_position:
			return fvd.GetComponent(GetPos2DD(p0));
		case FieldVariableDescriptor::FVT_velocity:
			return fvd.GetComponent(GetVel2DD(p0));
		case FieldVariableDescriptor::FVT_displacement:
			return fvd.GetComponent(GetPos2DD(p0) - Vector2D(x_init(1), x_init(2)));
		case FieldVariableDescriptor::FVT_acceleration:
			return fvd.GetComponent(GetAcceleration(p0));
		}
	}
	else
	{
		switch (fvd.VariableType())
		{
		case FieldVariableDescriptor::FVT_position:
			return fvd.GetComponent(GetPos2D(p0));
		case FieldVariableDescriptor::FVT_velocity:
			return fvd.GetComponent(GetVel2D(p0));
		case FieldVariableDescriptor::FVT_displacement:
			return fvd.GetComponent(GetPos2D(p0) - Vector2D(x_init(1), x_init(2)));
		case FieldVariableDescriptor::FVT_acceleration:
			return fvd.GetComponent(GetAcceleration(p0));
		}
	}
	return FIELD_VARIABLE_NO_VALUE;
}


void Mass2D::GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables)
{
	// add all the field variables of the parent class
	Body2D::GetAvailableFieldVariables(variables);
	//FVT_position,FVT_displacement, FVT_velocity: implemented in Element
	// possibility to remove entries does not exist yet, just do not use the Function of the parent class
	//$ SW 2013-10-18: add acceleration
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_acceleration, FieldVariableDescriptor::FieldVariableComponentIndex::FVCI_y);
}

void Mass2D::DrawElement()
{
	Body2D::DrawElement();

	mbs->SetColor(col);
	
	double drawrad = size.X();
	double fac = mbs->GetOptions()->PostProcOptions()->BodiesDeformationScaleFactor();
	if (!GetMBS()->GetIOption(151)) fac = 1;

	if (type & TParticle)
	{
		drawrad *= mbs->GetOptions()->PostProcOptions()->BodiesParticlesDrawSizeFactor();
		fac = mbs->GetOptions()->PostProcOptions()->BodiesParticlesDisplacementScaleFactor();
	}

	Vector2D p = GetRefPos2DD(); // fac*GetRefPos2DD() + (1-fac)*GetRefPosInit2D();    // p_scaled = (p-p_init)*fac + p_0 

	if (GetMBS()->GetActualPostProcessingFieldVariable() != NULL)
	{
		double v = GetFieldVariableValue(*GetMBS()->GetActualPostProcessingFieldVariable(), Vector2D(), true);
		mbs->DrawColorSphere(v, ToP3D(p), drawrad, drawres);
	}
	else
	{
		mbs->DrawSphere(ToP3D(p), drawrad, drawres);
	}
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Rigid2D::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Body2D::GetElementData(edc);

	ElementData ed;

	//int pos = edc.Find("Density");
	//if (pos) edc.Get(pos).SetLocked(1);
	ElementData* edf = edc.TreeFind("Physics.density");
	if (edf) edf->SetLocked(1);

	double frad = 1;
	if (GetMBS()->GetIOption(120) && !GetMBS()->IsLoadSaveMode()) {frad = 180./MY_PI;}

	//$RE 2013-08-20:[ moved to GetElementDataAuto
	//ed.SetVector3D(size.X(), size.Y(), size.Z(), "Body_dimensions"); ed.SetToolTipText("Dimensions of a regular cube [L_x, L_y, L_z]"); edc.Add(ed);
	//ed.SetDouble(mass, "mass"); ed.SetToolTipText("Used if 'Auto_compute_inertia' is deactivated only"); edc.TreeAdd("Physics",ed);
	
	//ed.SetDouble(Iphi, "principal_inertia");
	//ed.SetToolTipText("With respect to center of mass; used if 'Auto_compute_inertia' is deactivated"); edc.TreeAdd("Physics",ed);

	//ed.SetBool(0, "auto_compute_inertia"); ed.SetToolTipText("Principal inertia is auto-computed from mass and body_dimensions"); edc.TreeAdd("Physics",ed);

	////initial position / angle
	//ed.SetVector2D(x_init(1), x_init(2), "initial_position"); edc.TreeAdd("Initialization",ed);
	//ed.SetDouble(x_init(3)*frad, "initial_angle"); 
	//if (GetMBS()->GetIOption(120)) ed.SetToolTipText("° (Degrees)"); //alternatively: GetRotUnitStr(GetMBS()->GetIOption(120)).c_str()
	//else ed.SetToolTipText("(Radiant)"); 
	//edc.TreeAdd("Initialization",ed);

	////initial velocity / angular vel
	//ed.SetVector2D(x_init(4), x_init(5), "initial_velocity"); edc.TreeAdd("Initialization",ed);
	//ed.SetDouble(x_init(6)*frad, "initial_angular_velocity"); 
	//if (GetMBS()->GetIOption(120)) ed.SetToolTipText("(Degrees/s)"); 
	//else ed.SetToolTipText("(Radiant/s)");
	//edc.TreeAdd("Initialization",ed);
	//$RE 2013-08-20:]
}

int Rigid2D::SetElementData(ElementDataContainer& edc) //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
{
	int rv = Body2D::SetElementData(edc);

	double frad = 1;
	if (GetMBS()->GetIOption(120) && !GetMBS()->IsLoadSaveMode()) {frad = 180./MY_PI;}

	//GetElemDataVector3D(GetMBS(), edc, "Body_dimensions", size, 1);
	GetElemDataDouble(GetMBS(), edc, "Physics.mass", mass, 1);
	
	double V = size.X()*size.Y()*size.Z();
	
	rho = 0; //$RE 2013-08-20: 'rho' is not in the equations, possibly used in access-functions 'Rho()' => possibly redundant!
	if (V!=0) rho = mass/V;
	//$RE 2013-08-20:[ moved to SetElementDataAuto
	//int flag=0;
	//GetElemDataBool(GetMBS(), edc, "Physics.auto_compute_inertia", flag);
	//if (flag)
	//{
	//	Iphi = 1./12.*mass*(Sqr(size.X())+Sqr(size.Y()));
	//}
	//else
	//{
	//	GetElemDataDouble(GetMBS(), edc, "Physics.principal_inertia", Iphi, 1);
	//}
	//
	//Vector2D v;
	//GetElemDataVector2D(GetMBS(), edc, "Initialization.initial_position", v); x_init(1) = v.X(); x_init(2) = v.Y(); 
	//GetElemDataDouble(GetMBS(), edc, "Initialization.initial_angle", x_init(3), 1); x_init(3) /= frad;
	//GetElemDataVector2D(GetMBS(), edc, "Initialization.initial_velocity", v); x_init(4) = v.X(); x_init(5) = v.Y(); 
	//GetElemDataDouble(GetMBS(), edc, "Initialization.initial_angular_velocity", x_init(6), 1); x_init(6) /= frad;
	//$RE 2013-08-20:] moved to SetElementDataAuto

	return rv;
}

//$RE 2013-08-20:[ 
//this function assigns default values to the element variables
void Rigid2D::ElementDefaultConstructorInitialization() 
{	
	SetRigid2D(Vector(0.,0.,0., 0.,0.,0.), 1., Vector3D(0.1,0.1,0.01), colblue);	
}
//$RE 2013-08-20:]

void Rigid2D::GetdPosdqT(const Vector2D& ploc, Matrix& dpdqi)
{
	dpdqi.SetSize(3,2);
	dpdqi(1,1) = 1;
	dpdqi(1,2) = 0;
	dpdqi(2,1) = 0;
	dpdqi(2,2) = 1;
	// dR/dphi*ploc:		
	double phi = XG(3);
	double sphi = sin(phi);
	double cphi = cos(phi);
	dpdqi(3,1) = -ploc.X()*sphi-ploc.Y()*cphi;
	dpdqi(3,2) =  ploc.X()*cphi-ploc.Y()*sphi;
}

void Rigid2D::GetdRotvdqT(const Vector2D& vloc, const Vector2D& ploc, Matrix& d) 
{
	d.SetSize(3,2);
	d(1,1) = 0;
	d(1,2) = 0;
	d(2,1) = 0;
	d(2,2) = 0;
	// d(A*vloc)/dphi*ploc:		
	double phi = XG(3);
	double sphi = sin(phi);
	double cphi = cos(phi);
	d(3,1) = -vloc.X()*sphi - vloc.Y()*cphi;
	d(3,2) =  vloc.X()*cphi - vloc.Y()*sphi;
}

void Rigid2D::GetdAngle2DdqT(const Vector2D& ploc, Matrix& d)
{
	d.SetSize(3,1);
	d(1,1) = 0;
	d(2,1) = 0;
	d(3,1) = 1;
}


//---------------------------------------------------------------
//for visualization of displacements and velocities:
//---------------------------------------------------------------
double Rigid2D::GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector2D & local_position, bool flagD)
{	
	//$RE 2013-08-21:[ added according to "//$ DR 2012-12-12 added according to JG"
	Vector2D p0(local_position);

	if (flagD)
	{
		switch (fvd.VariableType())
		{
		case FieldVariableDescriptor::FVT_position:
			return fvd.GetComponent(GetPos2DD(p0));
		case FieldVariableDescriptor::FVT_velocity:
			return fvd.GetComponent(GetVel2DD(p0));
		case FieldVariableDescriptor::FVT_displacement:
			return fvd.GetComponent(GetPos2DD(p0) - Vector2D(x_init(1), x_init(2)));
		case FieldVariableDescriptor::FVT_bryant_angle:
			return XGD(3); //$ DR 2013-08-26
		case FieldVariableDescriptor::FVT_angular_velocity:	
			return XGPD(3); //$ DR 2013-08-26
		case FieldVariableDescriptor::FVT_acceleration: //$ SW 2013-10-18: added
			return fvd.GetComponent(GetAcceleration(p0));
		//case FieldVariableDescriptor::FVT_angular_velocity_local_basis://$ DR 2013-08-26 removed, does not make sense 
			//return fvd.GetComponent(GetAngularVelLocal());
		}
	}
	else
	{
		switch (fvd.VariableType())
		{
		case FieldVariableDescriptor::FVT_position:
			return fvd.GetComponent(GetPos2D(p0));
		case FieldVariableDescriptor::FVT_velocity:
			return fvd.GetComponent(GetVel2D(p0));
		case FieldVariableDescriptor::FVT_displacement:
			return fvd.GetComponent(GetPos2D(p0) - Vector2D(x_init(1), x_init(2)));	
		case FieldVariableDescriptor::FVT_bryant_angle:
			//return fvd.GetComponent(ToP3D(local_position)); 
			return XG(3);//$ DR 2013-08-26
		case FieldVariableDescriptor::FVT_angular_velocity:	
			return XGP(3);//$ DR 2013-08-26
		case FieldVariableDescriptor::FVT_acceleration://$ SW 2013-10-18: added
			return fvd.GetComponent(GetAcceleration(p0));
			//return fvd.GetComponent(GetAngularVel(ToP3D(local_position)));
		//case FieldVariableDescriptor::FVT_angular_velocity_local_basis: //$ DR 2013-08-26 removed, does not make sense 
			//return fvd.GetComponent(GetAngularVelLocal());
		}
	}
	return FIELD_VARIABLE_NO_VALUE;
	//$RE 2013-08-21:]
}


void Rigid2D::GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables)
{
	//$RE 2013-08-21:[
	Body2D::GetAvailableFieldVariables(variables);
	// add some specific ones of this class
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_bryant_angle, FieldVariableDescriptor::FVCI_x);		//$ DR 2013-08-26 changed to FVCI_x
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_angular_velocity, FieldVariableDescriptor::FVCI_x);  //$ DR 2013-08-26 changed to FVCI_x
	//FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_angular_velocity_local_basis, FieldVariableDescriptor::FVCI_z);  //$ DR 2013-08-26 removed, does not make sense
	//$RE 2013-08-21:]
	//$ SW 2013-10-18:[
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_acceleration, FieldVariableDescriptor::FVCI_y);
	//$ SW 2013-10-18:]
}

void Rigid2D::DrawElement() 
{
	Body2D::DrawElement();

	mbs->SetColor(col);

	//draw local frame at origin:
	if (GetMBS()->GetIOption(125))
	{
		double s = GetMBS()->GetDOption(104);

		GetMBS()->ChooseColor(0.3f,0.3f,0.3f);

		Vector3D v1 = ToP3D(GetPos2DD(Vector2D(-0.2*s, 0)));
		Vector3D v2 = ToP3D(GetPos2DD(Vector2D( s, 0)));
		Vector3D v3 = ToP3D(GetPos2DD(Vector2D( 0,-0.2*s)));
		Vector3D v4 = ToP3D(GetPos2DD(Vector2D( 0, s)));
		double d = GetMBS()->GetDOption(114);
		GetMBS()->MyDrawLine(v1,v2,d);
		GetMBS()->MyDrawLine(v3,v4,d);

		char str[20];
		sprintf(str, "X%d", GetOwnNum());
		GetMBS()->GetRC()->PrintText3D((float)v2.X(), (float)v2.Y(), (float)v2.Z(), str);
		sprintf(str, "Y%d", GetOwnNum());
		GetMBS()->GetRC()->PrintText3D((float)v4.X(), (float)v4.Y(), (float)v4.Z(), str);
	}




	Vector3D p8(GetPosD(Vector3D(-0.5*size.X(),-0.5*size.Y(),-0.5*size.Z())));
	Vector3D p7(GetPosD(Vector3D(-0.5*size.X(),-0.5*size.Y(), 0.5*size.Z())));
	Vector3D p6(GetPosD(Vector3D( 0.5*size.X(),-0.5*size.Y(),-0.5*size.Z())));
	Vector3D p5(GetPosD(Vector3D( 0.5*size.X(),-0.5*size.Y(), 0.5*size.Z())));
	Vector3D p4(GetPosD(Vector3D(-0.5*size.X(), 0.5*size.Y(),-0.5*size.Z())));
	Vector3D p3(GetPosD(Vector3D(-0.5*size.X(), 0.5*size.Y(), 0.5*size.Z())));
	Vector3D p2(GetPosD(Vector3D( 0.5*size.X(), 0.5*size.Y(),-0.5*size.Z())));
	Vector3D p1(GetPosD(Vector3D( 0.5*size.X(), 0.5*size.Y(), 0.5*size.Z())));

	mbs->DrawHex(p1,p2,p3,p4,p5,p6,p7,p8);

};



