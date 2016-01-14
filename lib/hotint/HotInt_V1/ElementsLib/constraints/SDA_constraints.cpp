//#**************************************************************
//#
//# filename:             SDA_constraints.cpp
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
 
#include "element.h"
#include "body3d.h"
#include "body2d.h"
//#include "femathhelperfunctions.h"
#include "material.h"
#include "node.h"
#include "rigid3d.h"
//#include "rigid3dkardan.h"
#include "constraint.h"
#include "sensors.h"
#include "control.h"
#include "sda_constraints.h"
//#include "myfile.h"
#include "elementdataaccess.h"
#include "geomelements.h"
#include "graphicsconstants.h"


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Spring Damper Actuator Element 3D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SDActor::SetSDActor(int en1, int en2, const Vector3D& lc1, const Vector3D& lc2,	double spring_stiffness,
												 double spring_init_len, double damping_coeff, double actuator_force, Vector3D cdim, const Vector3D& coli)
{	
	x_init = Vector(0);
	GetCol() = coli;
	draw_dim = cdim;
	AddElementCoord(en1, lc1);
	AddElementCoord(en2, lc2);
	k = spring_stiffness;
	l0 = spring_init_len;
	d = damping_coeff;
	fa = actuator_force;
	forcemode = 0;
	groundjoint = 0;
	InitForce();
	elementname = GetElementSpec();
	forcedirected = 0;
};

void SDActor::SetSDActor(int en1, const Vector3D& lc1, const Vector3D& gc2, double spring_stiffness, double spring_init_len, 
												 double damping_coeff, double actuator_force,Vector3D cdim, const Vector3D& coli)
{	
	x_init = Vector(0);
	GetCol() = coli;
	draw_dim = cdim;
	AddElementCoord(en1, lc1);
	k = spring_stiffness;
	l0 = spring_init_len;
	d = damping_coeff;
	fa = actuator_force;
	p_global = gc2;
	forcemode = 0;
	groundjoint = 1;
	InitForce();
	elementname = GetElementSpec();
	forcedirected = 0;
};

void SDActor::SetSDActor(int en1, const Vector3D& lc1, const Vector3D& gc2, double spring_stiffness, double spring_init_len, Vector3D& direction, 
												 double damping_coeff, double actuator_force,Vector3D cdim, const Vector3D& coli)
{	
	SetSDActor(en1, lc1,  gc2, spring_stiffness, spring_init_len, damping_coeff, actuator_force, cdim, coli);
	SetForceDirected(direction);
};

void SDActor::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Element::GetElementData(edc); //do not use Constraint::, because of element number
	InputOutputElement::GetElementData(edc); //data from InputOutputelement --> menu

	ElementData ed;

	//from Constraint::	
	int NEeff = NE();
	if (NEeff == 3-IsGroundJoint()) NEeff = 2-IsGroundJoint();

	Vector elnum(NEeff);
	for (int i=1; i <= NEeff; i++) {elnum(i) = GetElnum(i);}
	ed.SetVector(elnum.GetVecPtr(), elnum.Length(), "Constraint_element_numbers"); ed.SetToolTipText("Only valid element numbers permitted!"); edc.Add(ed);
	//++++++++++++++++++++

	//drawdim: X=spring drawsize, Y=windings, Z=damper size

	ed.SetDouble(draw_dim.X(), "Draw_size"); ed.SetToolTipText("Diameter of spring"); edc.Add(ed);
	ed.SetDouble(draw_dim.Y(), "Draw_windings"); ed.SetToolTipText("Number of windings in spring"); edc.Add(ed);
	ed.SetDouble(draw_dim.Z(), "Draw_damper_size"); ed.SetToolTipText("Diameter of damper, for drawing"); edc.Add(ed);

	ed.SetDouble(k, "Spring_stiffness"); edc.Add(ed);
	ed.SetDouble(l0, "Spring_initial_length"); edc.Add(ed);
	ed.SetDouble(d, "Damping"); edc.Add(ed);
	ed.SetDouble(fa, "Actuator_force"); edc.Add(ed);

	int fm = forcemode; if (fm == 0) fm = 2;

	ed.SetBoolGroup(fm==2, 3, "Linear-nonlinear"); edc.Add(ed);
	ed.SetBoolGroup(fm==1, 3, "Tension_only"); edc.Add(ed);
	ed.SetBoolGroup(fm==3, 3, "Clearance"); edc.Add(ed);
	ed.SetBoolGroup(fm==4, 3, "Controlled"); edc.Add(ed);

	if (groundjoint)
	{
		SetElemDataVector3D(edc, loccoords(1), "Local_joint_pos_body"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector3D(edc, p_global, "Global_joint_pos"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");

	} else
	{
		SetElemDataVector3D(edc, loccoords(1), "Local_joint_pos_body1"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector3D(edc, loccoords(2), "Local_joint_pos_body2"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
	}

	SetElemDataVector2D(edc, Vector2D(k2, k3), "Nonlinear_stiffness"); edc.Last().SetToolTipText("Quadratic and cubic stiffness parameters, f = k*x+k2*Sgn(x)*x^2+k3*x^3, [k2 k3]");
	SetElemDataVector2D(edc, Vector2D(x_range1, x_range2), "Clearance_range"); edc.Last().SetToolTipText("Range for x=l-l0 where no spring stiffness is active, [xmin xmax]");

	int celnum = 0;
	if (elements.Length() == 3-IsGroundJoint())
	{
		celnum = elements(3-IsGroundJoint());
	}
	ed.SetInt(celnum, "Controller_elem_num", 0, GetMBS()->NE()); 
	ed.SetToolTipText("Element number of controller element"); edc.Add(ed);

	int mode;
	Matrix data;
	mathfunc_k.GetData(mode, data);
	SetElemDataMatrix(edc, data, "Stiffness_curve"); 
	//mathfunc_k_flag ... not tunable in options yet 
	if(mathfunc_k_flag == 0) edc.Last().SetToolTipText("Set piecewise linear force values f versus displacements x=l-l0");
	else if(mathfunc_k_flag == 1)edc.Last().SetToolTipText("Set piecewise linear stiffness values c12(l); force results to c12(l)*(l-l0)");
	edc.Last().SetVariableLength();

	mathfunc_d.GetData(mode, data);
	SetElemDataMatrix(edc, data, "Damping_curve"); edc.Last().SetToolTipText("Set piecewise linear force values f versus velocity dx/dt");
	edc.Last().SetVariableLength();
}

int SDActor::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Element::SetElementData(edc); //do not use Constraint:: !!!!
    rv += InputOutputElement::SetElementData(edc);
	InitForce(); //initialize data, if some items are not read

	GetElemDataDouble(GetMBS(), edc, "Draw_size", draw_dim.X());
	GetElemDataDouble(GetMBS(), edc, "Draw_windings", draw_dim.Y());
	GetElemDataDouble(GetMBS(), edc, "Draw_damper_size", draw_dim.Z());

	GetElemDataDouble(GetMBS(), edc, "Spring_stiffness", k);
	GetElemDataDouble(GetMBS(), edc, "Spring_initial_length", l0);
	GetElemDataDouble(GetMBS(), edc, "Damping", d);
	GetElemDataDouble(GetMBS(), edc, "Actuator_force", fa);


	if (IsGroundJoint())
	{
		GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body", loccoords(1));
		GetElemDataVector3D(GetMBS(), edc, "Global_joint_pos", p_global);
	} else
	{
		GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body1", loccoords(1));
		GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body2", loccoords(2));
	}

	int flag = 0;
	forcemode = 0;

	GetElemDataBool(GetMBS(), edc, "Linear-nonlinear", flag,0);
	if (flag) forcemode = 2;
	flag = 0;
	GetElemDataBool(GetMBS(), edc, "Tension_only", flag,0);
	if (flag) forcemode = 1;
	flag = 0;
	GetElemDataBool(GetMBS(), edc, "Clearance", flag,0);
	if (flag) forcemode = 3;
	flag = 0;
	GetElemDataBool(GetMBS(), edc, "Controlled", flag,0);
	if (flag) forcemode = 4;

	int en = 0;
	GetElemDataInt(GetMBS(), edc, "Controller_elem_num", en, 0);
	if (en && forcemode == 4)
	{
		//if (GetMBS()->GetElement(en).IsType(TController))
		{
			elements.SetLen(3-IsGroundJoint());
			SetController(en);
		}
		/*else
		{
			GetMBS()->EDCError("Invalid controller element number in SDActor");
			forcemode = 0;
		}*/
	}
	else
	{
		elements.SetLen(2-IsGroundJoint());
		if (forcemode == 4) forcemode = 0;
	}

	//+++++++++++++++++++++++++++++++++++++++++++++++
	//from Constraint::
	Vector ennew;
	GetElemDataVector(GetMBS(), edc, "Constraint_element_numbers", ennew, 1);
	
	for (int i=1; i <= ennew.Length(); i++) 
	{
		int en = ennew(i);
		if (en <= 0 || en > GetMBS()->GetNElements()) 
		{
			GetMBS()->EDCError(mystr("Constraint: Element number ")+mystr(en)+mystr(" is out of range"));
			en = 1;
		}
		SetElnum(i, en);
	}
	//+++++++++++++++++++++++++++++++++++++++++++++++


	GetElemDataVector2D(GetMBS(), edc, "Nonlinear_stiffness", k2, k3, 0);
	GetElemDataVector2D(GetMBS(), edc, "Clearance_range", x_range1, x_range2, 0);


	//mathfunc_k_flag ... not tunable in options yet
	Matrix data;
	int pos = GetElemDataMatrix(GetMBS(), edc, "Stiffness_curve", data, 0);
	if (pos) mathfunc_k.SetData(TMFpiecewiselinear, data);

	pos = GetElemDataMatrix(GetMBS(), edc, "Damping_curve", data, 0);
	if (pos) mathfunc_d.SetData(TMFpiecewiselinear, data);

	return rv;
}

Vector3D SDActor::ComputeForceDirection() const   //return global force direction, normalized
{
	double l = 0;
	double xp = 0;
	Vector3D v12;

	if (forcedirected)
	{	
		v12 = GetBody3D(1).GetRotMatrix(loccoords(1))*locaxis1; // see SetForceDirected
	}
	else
	{
		if (IsGroundJoint())
		{
			v12 = GetBody3D(1).GetPos(loccoords(1)) - p_global;
		}
		else
		{
			v12 = GetBody3D(1).GetPos(loccoords(1)) - GetBody3D(2).GetPos(loccoords(2));
		}
	}
	v12.Normalize();
	
	return v12;
}

Vector3D SDActor::ComputeForceDirectionD() const   //return global force direction, normalized
{
	double l = 0;
	double xp = 0;
	Vector3D v12;

	if (forcedirected)
	{	
		v12 = GetBody3D(1).GetRotMatrixD(loccoords(1))*locaxis1; // see SetForceDirected
	}
	else
	{
		if (IsGroundJoint())
		{
			v12 = GetBody3D(1).GetPosD(loccoords(1)) - p_global;
		}
		else
		{
			v12 = GetBody3D(1).GetPosD(loccoords(1)) - GetBody3D(2).GetPosD(loccoords(2));
		}
	}
	v12.Normalize();
	return v12;
}
                                  
double SDActor::GetOutput(double t, int i) const 
{
	if (i==1 || i==2)  // First output = Spring Deflection
										 // Second output = Total length
	{
		double l = 0; 	

		Vector3D dir = ComputeForceDirection(); 

		double l0act = l0; 	

		if( GetNInputs() == 1 && ioElement_actionIndex(1) == (TConstraintIOType)	TSRefPos)
		{
			InputOutputElement* ioe = (InputOutputElement*)GetMBS()->GetElementPtr(inputs(1));
			l0act = ioe->GetOutput(t, 1); 
		}
		if (IsGroundJoint())
		{
			Vector3D v12 = GetBody3D(1).GetPos(loccoords(1)) - p_global;
			l = v12*dir;
		}
		else
		{
			Vector3D v12 = GetBody3D(1).GetPos(loccoords(1)) - GetBody3D(2).GetPos(loccoords(2));
			l = v12*dir;
		}

		double x = (l-l0act); 	

		if (i==1) return x; else return l;
	}
	else if(i==3)
	{
		Vector3D dir = ComputeForceDirection(); 
		return ComputeForce(t)*dir;
	}	
	else
	{
		assert(0 && "output not defined.");
		return 0.;
	}
}

Vector3D SDActor::ComputeForce(double t)  const
{
	double l = 0;
	double xp = 0;
	Vector3D dir = ComputeForceDirection();
	double fact = fa;
	double l0act = l0;
	if( GetNInputs() == 1 && ioElement_actionIndex(1) == (TConstraintIOType)TSRefForce)
	{
		InputOutputElement* ioe = (InputOutputElement*)GetMBS()->GetElementPtr(inputs(1));
		fact = ioe->GetOutput(t, 1); 
	}

	if( GetNInputs() == 1 && ioElement_actionIndex(1) == (TConstraintIOType)	TSRefPos)
	{
		InputOutputElement* ioe = (InputOutputElement*)GetMBS()->GetElementPtr(inputs(1));
		l0act = ioe->GetOutput(t, 1); 
	}

	if (IsGroundJoint())
	{
		Vector3D v12 = GetBody3D(1).GetPos(loccoords(1)) - p_global;
		l = v12*dir;

		if (d != 0)
			xp = (GetBody3D(1).GetVel(loccoords(1)))*dir;
	}
	else
	{
		Vector3D v12 = GetBody3D(1).GetPos(loccoords(1)) - GetBody3D(2).GetPos(loccoords(2));
		l = v12*dir;

		if (d != 0) 
			xp = (GetBody3D(1).GetVel(loccoords(1)) - GetBody3D(2).GetVel(loccoords(2)))*dir;
	}
	
	double x = (l-l0act); 


	if (forcemode <= 2) //linear/only tension/nonlinear/piecewise
	{
		if (forcemode == 1 && x < 0) {x = 0; xp = 0;}

		double f = x*k + Sqr(x)*Sgn(x)*k2 + Cub(x)*k3;
		double fd = xp*d;
	
		if (forcemode != 1)
		{
			if (mathfunc_k.GetFuncMode() == TMFpiecewiselinear && mathfunc_k_flag == 1)f = mathfunc_k.Evaluate(l)*x; //mathfunc_k ... stiffness depending on 'l'
			else if(mathfunc_k.GetFuncMode() == TMFpiecewiselinear)f = mathfunc_k.Evaluate(x); //mathfunc_k ... force
			if (mathfunc_d.GetFuncMode() == TMFpiecewiselinear) fd = mathfunc_d.Evaluate(xp);
		}

		return (f + fd + fact)*dir;
	}
	else if (forcemode == 3) //"clearance"
	{
		double fd = xp*d;
		if (x >= x_range1 && x <= x_range2) return (fact + fd)*dir;

		if (x <= x_range1) x -= x_range1;
		else if (x >= x_range2) x -= x_range2;

		double f = x*k + Sqr(x)*Sgn(x)*k2 + Cub(x)*k3;
		return (f + fd + fact)*dir;
	}
	else if (forcemode == 4) //controller
	{
		GetMBS()->UO() << "Error: MBSController not available any more!!!\n";
		return 0;
		//if (NE() != 3-groundjoint) GetMBS()->UO() << "Error: SDActor, no MBSController defined!!!\n";
		//else
		//{
		//	if (GetElem(3-groundjoint).IsType(TController))
		//	{
		//		MBSController& c = (MBSController&)GetElem(3-groundjoint);
		//		fact = c.ComputeControlValue();
		//		//UO() << "ma=" << ma << "\n";
		//	}
		//	else GetMBS()->UO() << "Error: SDActor, Element is no MBSController!!!\n";
		//}
		//return fact*dir;
	}
	else return Vector3D(0); //should not happen!
}

double SDActor::GetActorForce(double computation_time, int dir) const 
{
	if (dir >=1 && dir <= 3) return (ComputeForce(computation_time))(dir);
	else return ComputeForce(computation_time)*ComputeForceDirection();
} 
int SDActor::GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables)
{
	// call base class routines
	InputOutputElement::GetAvailableSpecialValues(available_variables);

	// Automatic entries for this class 
	SDActor::GetAvailableSpecialValuesAuto(available_variables);

	// Manual entries for this class
	available_variables.Add(ReadWriteElementDataVariableType("Internal.force_axial",0,0,0.,mystr("force in axial direction, applied to the kinematic pairs due to the constraint"))) ;
	available_variables.Add(ReadWriteElementDataVariableType("Internal.force",3,0,0.,mystr("force applied to the kinematic pairs due to the constraint. range: 1-3 corresponds to force in x-y-z direction"))) ;

	return 0;
}

int SDActor::ReadSingleElementData(ReadWriteElementDataVariableType& RWdata) 		
{
	// call base class routine 	
	int rv = Constraint::ReadSingleElementData(RWdata);
	if (rv == 1) return 1;

	// manual things to read 
	double t = GetMBS()->GetTime();
	Vector3D f = ComputeForce(t);


 if( RWdata.variable_name == mystr("Internal.force_axial") )
	{
		RWdata.value = f*ComputeForceDirection(); 
		return 1; 
	}
	else	if( RWdata.variable_name == mystr("Internal.force") )
	{
		if (RWdata.comp1 > 0 && RWdata.comp1 <= 3) //range check
		{
			RWdata.value = ComputeForce(t)(RWdata.comp1); 
			return 1; 
		}
		else return -2; 
	}
	return ReadSingleElementDataAuto(RWdata);
}


void SDActor::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	//-C_q^T \lambda = [dC/dx; dC/dy; dC/dphi]^T [\lambda1;\lambda2] 
	//C = p_ref+R*p_loc
	if (locelemind >= 3-IsGroundJoint()) return; //do not add lagrange multipliers to controller

	double sign = 1;
	if (locelemind == 2) sign = -1;

	if (locelemind > elements.Length()) {UO() << "ERROR: inconsistency in SDActor::AddElementCqTLambda\n";}

	GetBody3D(locelemind).GetdPosdqT(loccoords(locelemind),dpdq);
	Vector3D F(ComputeForce(t));
	for (int i=1; i <= f.Length(); i++)
	{
		f(i) -= sign*(dpdq(i,1)*F(1)+dpdq(i,2)*F(2)+dpdq(i,3)*F(3));
	}
};

Vector3D SDActor::GetRefPosD()	const 
{
	return GetBody3D(1).GetPosD(loccoords(1));
}

void SDActor::DrawElement() 
{
	Constraint::DrawElement();

	mbs->SetColor(GetCol());

	Vector3D pa = GetBody3D(1).GetPosD(loccoords(1));
	Vector3D pb;
	if (IsGroundJoint())
		pb = p_global;
	else
		pb = GetBody3D(2).GetPosD(loccoords(2));

	if (IsGroundJoint()) Swap(pa,pb);

	if (draw_dim.Y() == 0)
	{
		double l = (pa - pb).Norm();
		if (l < l0 && forcemode == 1) mbs->SetColor(Vector3D(0.95,0.95,0.95));

		mbs->DrawZyl(pa, pb, 0.5*draw_dim.X(),12);
	}
	else
	{
		//draw spring:
		double phi,phi2,xx;
		Vector3D vf;
		double steps = 200*0.5; //0.2 ->> very coarse
		double rots = draw_dim.Y();
		double r = draw_dim.X()/2.;
		double rz = draw_dim.X()/10.;
		int ztile = 6;

		vf = pa-pb;
		double lenvf = vf.Norm();
		Vector3D dir = vf;
		vf = ComputeForceDirectionD();
		vf.Normalize();
		vf *= lenvf;
		if (vf*dir < 0) vf *= -1;

		Vector3D vfn = vf;
		Vector3D vz,vy;
		vfn.SetNormalBasis(vz,vy);
		Vector3D pz1, pz2;

		double off = steps*0.1;
		for (xx = off; xx <= steps-off; xx++)
		{
			phi = xx*2*MY_PI/steps*rots;
			phi2 = (xx+1.2)*2*MY_PI/steps*rots;

			Vector3D pr1=cos(phi)*r*vz+sin(phi)*r*vy;
			Vector3D pr2=cos(phi2)*r*vz+sin(phi2)*r*vy;

			pz1=pb+xx/steps*vf+pr1;
			pz2=pb+(xx+1)/steps*vf+pr2;
			mbs->DrawZyl(pz1,pz2,rz,ztile);
			if (xx==off) mbs->DrawZyl(pb,pz1,rz,ztile);
		}
		mbs->DrawZyl(pa,pz2,rz,ztile);

		if (draw_dim.Z() != 0 && d != 0)
		{
			mbs->DrawZyl(pa, pa-0.5*l0*vfn, 0.5*draw_dim.Z(),12);
			mbs->SetColor(colgrey2);
			mbs->DrawZyl(pa-0.5*l0*vfn, pb, 0.5*draw_dim.Z()*0.6,12);

			mbs->DrawSphere(pa,0.65*draw_dim.Z(),12);
			mbs->DrawSphere(pb,0.65*draw_dim.Z(),12);
		}
	}

};

// ++++++++++++++++++++++++
//b: for InputOutputelement
Vector2D SDActor::GetInputPosD(int i) const //dummy, must be programmed for 2D-Actors
{
	assert(0 && "SDActor is not 2D!!");
	return Vector2D(0.);
}


Vector3D SDActor::GetInputPos3DD(int i) const //return absolute position of input #i
{
	return GetPosD(Vector3D(0.));	
}

// similar to InputOutputElement
Vector2D SDActor::GetOutputPosD(int i) const //return absolute position of input #i
{
	Vector3D vec=GetPosD(Vector3D(0.));
	return Vector2D(vec(1),vec(2));
}

Vector3D SDActor::GetOutputPos3DD(int i) const //return absolute position of input #i
{
	return GetPosD(Vector3D(0.));
}
//e: for InputOutputelement
// ++++++++++++++++++++++++



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Rotational Spring Damper Actuator Element
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void SDRotActor::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Element::GetElementData(edc); //do not use Constraint::, because of element number
	ElementData ed;

	double frad = 1; //factor for radiant to degree
	if (GetMBS()->GetIOption(120) && !GetMBS()->IsLoadSaveMode()) {frad = 180./MY_PI;}

	//from Constraint::	
	int NEeff = NE();
	if (NEeff == 3-IsGroundJoint()) NEeff = 2-IsGroundJoint();

	Vector elnum(NEeff);
	for (int i=1; i <= NEeff; i++) {elnum(i) = GetElnum(i);}
	ed.SetVector(elnum.GetVecPtr(), elnum.Length(), "Constraint_element_numbers"); ed.SetToolTipText("Only valid element numbers permitted!"); edc.Add(ed);
	//++++++++++++++++++++

	//drawdim: X=spring drawsize, Y=windings, Z=axis radius (cylinder)
	ed.SetDouble(draw_dim.X(), "Draw_size"); ed.SetToolTipText("Radius of torsional spring"); edc.Add(ed);
	ed.SetDouble(draw_dim.Y(), "Draw_windings"); ed.SetToolTipText("Number of windings of torsional spring"); edc.Add(ed);
	ed.SetDouble(draw_dim.Z(), "Draw_axis_radius"); ed.SetToolTipText("Radius of torsional spring axis"); edc.Add(ed);

	ed.SetDouble(k, "Spring_stiffness"); ed.SetToolTipText("Torsional stiffness K of the spring, the torque is K*(phi-phi0)"); edc.Add(ed);

	ed.SetDouble(phi0*frad, "Spring_angle_offset");
	ed.SetToolTipText((mystr("angle offset phi0 ")+GetRotUnitStr(GetMBS()->GetIOption(120))).c_str()); edc.Add(ed);

	ed.SetDouble(d, "Damping"); ed.SetToolTipText("Torsional damping D of the spring, the torque is D times angular velocity"); edc.Add(ed);
	ed.SetDouble(ma, "Actuator_torque"); ed.SetToolTipText("Constant torque of an actuator"); edc.Add(ed);

	int fm = forcemode; if (fm == 0) fm = 2;

	ed.SetBoolGroup(fm==2, 3, "Linear-nonlinear"); edc.Add(ed);
	ed.SetBoolGroup(fm==3, 3, "Clearance"); edc.Add(ed);
	ed.SetBoolGroup(fm==4, 3, "Controlled"); edc.Add(ed);

	if (groundjoint)
	{
		SetElemDataVector3D(edc, loccoords(1), "Local_joint_pos_body"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector3D(edc, loccoords(2), "Global_rotation_axis"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
	} 
	else
	{
		SetElemDataVector3D(edc, loccoords(1), "Local_joint_pos_body1"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector3D(edc, loccoords(2), "Local_joint_pos_body2"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector3D(edc, loccoords(3), "Global_rotation_axis"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
	}

	SetElemDataVector2D(edc, Vector2D(k2, k3), "Nonlinear_stiffness"); edc.Last().SetToolTipText("Quadratic and cubic stiffness parameters, f = k*x+k2*Sgn(x)*x^2+k3*x^3, [k2 k3]");
	SetElemDataVector2D(edc, Vector2D(phi_range1, phi_range2), "Clearance_range"); edc.Last().SetToolTipText("Range for phi [rad] where no spring stiffness is active, [phi_min phi_max]");

	int celnum = 0;
	if (elements.Length() == 3-IsGroundJoint())
	{
		celnum = elements(3-IsGroundJoint());
	}
	ed.SetInt(celnum, "Controller_elem_num", 0, GetMBS()->NE()); 
	ed.SetToolTipText("Element number of controller element"); edc.Add(ed);

	int mode;
	Matrix data;
	mathfunc_k.GetData(mode, data);
	SetElemDataMatrix(edc, data, "Stiffness_curve"); edc.Last().SetToolTipText("Set piecewise linear torque values m versus relative rotation");
	edc.Last().SetVariableLength();

	mathfunc_d.GetData(mode, data);
	SetElemDataMatrix(edc, data, "Damping_curve"); edc.Last().SetToolTipText("Set piecewise linear torque values m versus angular velocity dphi/dt");
	edc.Last().SetVariableLength();
}

int SDRotActor::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Element::SetElementData(edc); //do not use Constraint:: !!!!

	double frad = 1;
	if (GetMBS()->GetIOption(120) && !GetMBS()->IsLoadSaveMode()) {frad = 180./MY_PI;}

	InitForce(); //initialize data, if some items are not read

	GetElemDataDouble(GetMBS(), edc, "Draw_size", draw_dim.X());
	GetElemDataDouble(GetMBS(), edc, "Draw_windings", draw_dim.Y());
	GetElemDataDouble(GetMBS(), edc, "Draw_axis_radius", draw_dim.Z());

	GetElemDataDouble(GetMBS(), edc, "Spring_stiffness", k);
	GetElemDataDouble(GetMBS(), edc, "Spring_angle_offset", phi0); phi0 /= frad;
	GetElemDataDouble(GetMBS(), edc, "Damping", d);
	GetElemDataDouble(GetMBS(), edc, "Actuator_torque", ma);


	if (IsGroundJoint())
	{
		GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body", loccoords(1));
		GetElemDataVector3D(GetMBS(), edc, "Global_rotation_axis", loccoords(2));
	} else
	{
		GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body1", loccoords(1));
		GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body2", loccoords(2));
		GetElemDataVector3D(GetMBS(), edc, "Global_rotation_axis", loccoords(3));
	}


	int flag = 0;
	forcemode = 0;

	GetElemDataBool(GetMBS(), edc, "Linear-nonlinear", flag,0);
	if (flag) forcemode = 2;
	flag = 0;
	GetElemDataBool(GetMBS(), edc, "Clearance", flag,0);
	if (flag) forcemode = 3;
	flag = 0;
	GetElemDataBool(GetMBS(), edc, "Controlled", flag,0);
	if (flag) forcemode = 4;

	int en = 0;
	GetElemDataInt(GetMBS(), edc, "Controller_elem_num", en, 0);
	if (en && forcemode == 4)
	{
		//if (GetMBS()->GetElement(en).IsType(TController))
		{
			elements.SetLen(3-IsGroundJoint());
			SetController(en);
		}
		/*else
		{
			GetMBS()->EDCError("Invalid controller element number in SDActor");
			forcemode = 0;
		}*/
	}
	else
	{
		elements.SetLen(2-IsGroundJoint());
		if (forcemode == 4) forcemode = 0;
	}

	//+++++++++++++++++++++++++++++++++++++++++++++++
	//from Constraint::
	Vector ennew;
	GetElemDataVector(GetMBS(), edc, "Constraint_element_numbers", ennew, 1);
	
	for (int i=1; i <= ennew.Length(); i++) 
	{
		int en = ennew(i);
		if (en <= 0 || en > GetMBS()->GetNElements()) 
		{
			GetMBS()->EDCError(mystr("Constraint: Element number ")+mystr(en)+mystr(" is out of range"));
			en = 1;
		}
		SetElnum(i, en);
	}
	//+++++++++++++++++++++++++++++++++++++++++++++++


	GetElemDataVector2D(GetMBS(), edc, "Nonlinear_stiffness", k2, k3, 0);
	GetElemDataVector2D(GetMBS(), edc, "Clearance_range", phi_range1, phi_range2, 0);


	Matrix data;
	int pos = GetElemDataMatrix(GetMBS(), edc, "Stiffness_curve", data, 0);
	if (pos) mathfunc_k.SetData(TMFpiecewiselinear, data);

	pos = GetElemDataMatrix(GetMBS(), edc, "Damping_curve", data, 0);
	if (pos) mathfunc_d.SetData(TMFpiecewiselinear, data);


	return rv;
}




void SDRotActor::Initialize()
{
	if (IsGroundJoint())
	{
		Vector3D lpos = loccoords(1);
		Vector3D grot = loccoords(2);
		grot.Normalize();

		Vector3D gn1, gt1;
		grot.SetNormalBasis(gn1,gt1);

		Matrix3D RT=GetBody3D(1).GetRotMatrix(lpos).GetTp();

		loccoords(3)=RT*grot; //local (body fixed) rotation axis
		loccoords(4)=RT*gn1;	//local (body fixed) normal to rotation axis
		loccoords(5)=gn1;			//global normal to rotation axis
	} 
	else
	{
		Vector3D lp1 = loccoords(1);
		Vector3D lp2 = loccoords(2);
		Vector3D grot = loccoords(3);
		grot.Normalize();

		Vector3D gn1, gt1;
		grot.SetNormalBasis(gn1,gt1);

		Matrix3D RT1=GetBody3D(1).GetRotMatrix(lp1).GetTp();
		Matrix3D RT2=GetBody3D(2).GetRotMatrix(lp2).GetTp();

		loccoords(4)=RT1*grot; //local (body fixed) rotation axis
		loccoords(5)=RT2*grot; //local (body fixed) rotation axis
		loccoords(6)=RT1*gn1;	 //local (body fixed) normal to rotation axis
		loccoords(7)=RT2*gn1;	 //local (body fixed) normal to rotation axis
	}
}

double SDRotActor::GetOutput(double t, int i) const 
{

	if (i==1 || i==2)  // First output = Spring Deflection
										 // Second output = Total length
	{

		double phi0act = phi0;

		if( GetNInputs() == 1 && ioElement_actionIndex(1) == (TConstraintIOType)	TSRefAngle)
		{
			InputOutputElement* ioe = (InputOutputElement*)GetMBS()->GetElementPtr(inputs(1));
			phi0act = ioe->GetOutput(t, 1); 
		}

		double phi, phip;
		if (IsGroundJoint())
		{
			const Vector3D& lp1 = loccoords(1);
			Matrix3D R1=GetBody3D(1).GetRotMatrix(lp1);

			const Vector3D& rot = loccoords(2);		//==R1 * loccoords(3)
			Vector3D n1 = R1*loccoords(4); //rotation axis normal body1
			const Vector3D& n2 = loccoords(5);		//rotation axis normal global

			phi = NormalizedVectorAngle(n2,n1); //phi should always be positive
			if (phi < 0) phi += 2.*MY_PI; //should not happen!

			if (rot*(n2.Cross(n1)) < 0) phi = -phi;
		}
		else
		{
			const Vector3D& lp1 = loccoords(1);
			const Vector3D& lp2 = loccoords(2);
			Matrix3D R1=GetBody3D(1).GetRotMatrix(lp1);
			Matrix3D R2=GetBody3D(2).GetRotMatrix(lp2);

			Vector3D rot= R1*loccoords(4);		//actual global rotation axis body 1 == body 2
			Vector3D rot2= R2*loccoords(5);		//actual global rotation axis body 1 == body 2
			Vector3D n1 = R1*loccoords(6);		//rotation axis normal body1
			Vector3D n2 = R2*loccoords(7);		//rotation axis normal body2

			n1.Normalize(); //not necessary because Rigid3D Rotationmatrix is now consistent
			n2.Normalize();

			phi = NormalizedVectorAngle(n1,n2); //phi should always be positive
			if (phi < 0) phi += 2.*MY_PI; //should not happen!

			if ((0.5*(rot+rot2))*(n1.Cross(n2)) < 0) phi = -phi;

		}

		if (i==1) return phi-phi0; else return phi;
	}
	else
	{
		assert(0 && "output not defined.");
		return 0.;
	}
}

double SDRotActor::ComputeTorque(double t) const
{

	double mact = ma;
	double phi0act = phi0;
	if( GetNInputs() == 1 && ioElement_actionIndex(1) == (TConstraintIOType)TSRefMom)
	{
		InputOutputElement* ioe = (InputOutputElement*)GetMBS()->GetElementPtr(inputs(1));
		mact = ioe->GetOutput(t, 1); 
	}

	if( GetNInputs() == 1 && ioElement_actionIndex(1) == (TConstraintIOType)	TSRefAngle)
	{
		InputOutputElement* ioe = (InputOutputElement*)GetMBS()->GetElementPtr(inputs(1));
		phi0act = ioe->GetOutput(t, 1); 
	}
	//compute reacton moment:

	double phi, phip;
	if (IsGroundJoint())
	{
		const Vector3D& lp1 = loccoords(1);
		Matrix3D R1=GetBody3D(1).GetRotMatrix(lp1);

		const Vector3D& rot = loccoords(2);		//==R1 * loccoords(3)
		Vector3D n1 = R1*loccoords(4); //rotation axis normal body1
		const Vector3D& n2 = loccoords(5);		//rotation axis normal global
		//Vector3D t1 = rot.Cross(n1);		//rotation axis second normal body 1

		Vector3D angvel = GetBody3D(1).GetAngularVel(lp1);
		phip = -(rot*angvel);

		phi = NormalizedVectorAngle(n2,n1); //phi should always be positive
		if (phi < 0) phi += 2.*MY_PI; //should not happen!

		if (rot*(n2.Cross(n1)) < 0) phi = -phi;
		//if (!GetMBS()->IsJacobianComputation())
		//{
		//	 GetMBS()->UO() << "t=" << GetMBS()->GetTime() << ", phi=" << phi << "\n";
		//}
	}
	else
	{
		const Vector3D& lp1 = loccoords(1);
		const Vector3D& lp2 = loccoords(2);
		Matrix3D R1=GetBody3D(1).GetRotMatrix(lp1);
		Matrix3D R2=GetBody3D(2).GetRotMatrix(lp2);

		Vector3D rot= R1*loccoords(4);		//actual global rotation axis body 1 == body 2
		Vector3D rot2= R2*loccoords(5);		//actual global rotation axis body 1 == body 2
		Vector3D n1 = R1*loccoords(6);		//rotation axis normal body1
		Vector3D n2 = R2*loccoords(7);		//rotation axis normal body2
		//Vector3D t1 = rot.Cross(n1);			//rotation axis second normal body 1

		Vector3D angvel1 = GetBody3D(1).GetAngularVel(lp1);
		Vector3D angvel2 = GetBody3D(2).GetAngularVel(lp2);
		
		//phip = (rot*angvel1 - (R2*loccoords(5))*angvel2);
		phip = (rot*(angvel1 - angvel2));

		n1.Normalize(); //not necessary because Rigid3D Rotationmatrix is now consistent
		n2.Normalize();

		phi = NormalizedVectorAngle(n1,n2); //phi should always be positive
		if (phi < 0) phi += 2.*MY_PI; //should not happen!

		if ((0.5*(rot+rot2))*(n1.Cross(n2)) < 0) phi = -phi;

	}

	if (forcemode == 0 || forcemode == 2) //linear/nonlinear/piecewise
	{
		double x = phi-phi0act;
		double f = x*k + Sqr(x)*Sgn(x)*k2 + Cub(x)*k3;
		double fd = phip*d;
		if (mathfunc_k.GetFuncMode() == TMFpiecewiselinear)
		{
			f = mathfunc_k.Evaluate(x);
		}
		if (mathfunc_d.GetFuncMode() == TMFpiecewiselinear) 
		{
			fd = mathfunc_d.Evaluate(phip);
		}
		return mact + f - fd;
	}
	else if (forcemode == 3) //"clearance"
	{
		double x = phi-phi0act;
		if (x >= phi_range1 && x <= phi_range2) 
		{
			//GetMBS()->UO() << "clearance=" << x << "\n";
			return mact - phip * d;
		}

		if (x <= phi_range1) x -= phi_range1;
		else if (x >= phi_range2) x -= phi_range2;

		return mact + x*k + Sqr(x)*Sgn(x)*k2 + Cub(x)*k3 - phip * d;
	}
	else if (forcemode == 4) //controller
	{
		GetMBS()->UO() << "Error: MBSController not available any more!!!\n";
		return 0;
		//if (NE() != 3-groundjoint) GetMBS()->UO() << "Error: SDRotActor, no MBSController defined!!!\n";
		//else
		//{
		//	if (GetElem(3-groundjoint).IsType(TController))
		//	{
		//		MBSController& c = (MBSController&)GetElem(3-groundjoint);
		//		mact = c.ComputeControlValue();
		//		//UO() << "mact=" << mact << "\n";
		//	}
		//	else GetMBS()->UO() << "Error: SDRotActor, Element is no MBSController!!!\n";
		//}
		//return mact;
	}
	else return 0; //should not happen!
}

void SDRotActor::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{

	double sign = 1;
	if (locelemind == 2) sign *= -1;
	if (locelemind > elements.Length()) {UO() << "ERROR: inconsistency in SDRotActor::AddElementCqTLambda\n";}

	if (locelemind == 3-groundjoint) return;

	Vector3D rot;
	if (IsGroundJoint())
	{
		sign = -1;
		rot = loccoords(2);		//global rotation axis
	}
	else
	{
		rot = GetBody3D(2).GetRotMatrix(loccoords(2))*loccoords(5);	//global rotation axis body2
		//Vector3D rot1 = GetBody3D(1).GetRotMatrix(loccoords(1))*loccoords(4);	//global rotation axis body1
		//rot = 0.5*(rot+rot1);
	}

	Vector3D torque((sign* ComputeTorque(t)) * rot);

	//Apply torque to body 'locelemind', see MBSload
	if (GetBody3D(locelemind).IsRigid())
	{
		static Vector fadd;
		fadd.SetLen(f.Length()); fadd.SetAll(0);

		GetBody3D(locelemind).ApplyDrotrefdq(fadd, torque); 

		f += fadd;
	}
	else
	{
		//see also MBSLoad: Tpointmoment!!!
		//UO() << "ERROR: SDRotActor not implemented for flexible elements!\n";

		static Matrix H;
		H.SetSize(f.Length(), GetBody3D(locelemind).Dim()); //done also in GetdRotdqT
		static Vector fadd;
		fadd.SetLen(f.Length());

		GetBody3D(locelemind).GetdRotdqT(loccoords(locelemind), H);
		Mult(H, torque, fadd);

		f += fadd;

		//GetBody3D(locelemind).GetdRotvdqT(rot,loccoords(locelemind),dpdq); //not implemented for ANCFbeam!
		
	}

};

Vector3D SDRotActor::GetRefPosD()	const 
{
	return GetBody3D(1).GetPosD(loccoords(1));
}

void SDRotActor::DrawElement() 
{
	Constraint::DrawElement();

	double r = draw_dim.X(); //radius of torsional spring
	double rev = draw_dim.Y(); //number of revolutions for torsional spring
	double t = draw_dim.Z(); //radius of axis cylinder

	mbs->SetColor(GetCol());

	Vector3D p1 = GetBody3D(1).GetPosD(loccoords(1));
	Vector3D p2, rot, n1, t1, t2;
	if (IsGroundJoint())
	{
		Matrix3D R1=GetBody3D(1).GetRotMatrixD(loccoords(1));
		rot = loccoords(2);
		n1 = R1*loccoords(4);
		t1 = rot.Cross(n1);
		t2 = rot.Cross(loccoords(5));
		p2 = p_global;
	}
	else
	{
		Matrix3D R1=GetBody3D(1).GetRotMatrixD(loccoords(1));
		Matrix3D R2=GetBody3D(2).GetRotMatrixD(loccoords(2));

		rot = R1*loccoords(4);
		n1 = R1*loccoords(6);
		t1 = rot.Cross(n1);
		t2 = rot.Cross(R2*loccoords(7));

		p2 = GetBody3D(2).GetPosD(loccoords(2));
	}

	double dphi = NormalizedVectorAngle(t2,n1);

	double phioff = 0;
	dphi=-MY_PI*0.5+dphi; //only 90°+90°
	if (IsGroundJoint()) dphi += MY_PI;
	else 
	{
		phioff = dphi;
		dphi = -dphi + MY_PI;
	}

	if (r != 0) 
	{
		int mode = 2; //mode=1: zig-zag, mode=2: revolutions
		double ntile = 16*mode; //tiling
		if (mode == 2) dphi += rev*MY_PI; //6 revolutions
		for (double i = 1; i <= ntile; i++)
		{
			double phia = (i-1.)/ntile*dphi+phioff;
			double phib = (i)/ntile*dphi+phioff;
			double ra = r;
			double rb = r;

			if (mode == 1)
			{
				if ((i-1) > 1 && (i-1) < ntile-1 && ((int)i-1)% 2 == 0) ra = r*0.85;
				if (i > 1 && i < ntile-1 && (int)i % 2 == 0) rb = r*0.85;
			}
			else
			{
				ra *= (1.-0.33*phia/dphi);
				rb *= (1.-0.33*phib/dphi);
			}
			Vector3D va = p1 + ra*cos(phia)*n1 + ra*sin(phia)*t1;
			Vector3D vb = p1 + rb*cos(phib)*n1 + rb*sin(phib)*t1;

			mbs->MyDrawLine(va,vb,2);
			//mbs->DrawZyl(va,vb,t,6);
		}
	}

	mbs->DrawZyl(p1-(0.5*r)*rot,p1+(0.5*r)*rot,t,8);

	//double rs = 0.65*draw_dim.Z();//original
	double rs = 0.02*draw_dim.X();
	mbs->DrawSphere(p1,rs,12);
	mbs->SetColor(colred);
	mbs->DrawSphere(p2,rs,12);
	/*
	//draw n1, t1, t2:
	mbs->SetColor(colblue2);
	mbs->DrawZyl(p1-(0.5*r)*n1,p1+(0.5*r)*n1,t,8);
	mbs->SetColor(colred);
	mbs->DrawZyl(p1-(0.5*r)*t1,p1+(0.5*r)*t1,t,8);

	mbs->SetColor(colbrown);
	mbs->DrawZyl(p1-(0.5*r)*t2,p1+(0.5*r)*t2,t,8);
	*/

};



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Spring Damper Actuator Element 2D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Vector2D SDActor2D::ComputeForce(double t) 
{
	double l = 0;
	double lp = 0;
	Vector2D v12;
	if (elements.Length() == 1)
	{
		v12 = GetBody2D(1).GetPos2D(loccoords(1)) - p_global;
		l = v12.Norm();
		if (l!=0) v12.Normalize();
		if (d != 0) 
			lp = (GetBody2D(1).GetVel2D(loccoords(1)))*v12;
	}
	else
	{
		v12 = GetBody2D(1).GetPos2D(loccoords(1)) - GetBody2D(2).GetPos2D(loccoords(2));
		l = v12.Norm();
		if (l!=0) v12.Normalize();
		if (d != 0) 
			lp = (GetBody2D(1).GetVel2D(loccoords(1)) - GetBody2D(2).GetVel2D(loccoords(2)))*v12;
	}
	double deltad = (l-l0); 
	if (forcemode == 1 && deltad < 0) {deltad = 0; lp = 0;}
  double kstiff = k;
  if (k_delay > -1)
  {
    double tt = mbs->GetTime();
    double funct_fact = kparams(2);    
    double layer = kparams(3);
    double k_stiff_max = kparams(4);
    kstiff = k*exp(-funct_fact*(tt-(layer-1)*k_delay));
    if (kstiff > k_stiff_max) kstiff = k_stiff_max;
    //if (mbs->GetTime() > 0.02) mbs->UO() << "ElNr = " << GetBody2D(1).GetOwnNum() << "; LayerNr = " << layer << "; Spring stiffness = " << kstiff << "\n";
  }
	return (kstiff*deltad+d*lp+EvaluateActorForce())*v12;
}

void SDActor2D::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	//-C_q^T \lambda = [dC/dx; dC/dy; dC/dphi]^T [\lambda1;\lambda2] 
	//C = p_ref+R*p_loc
	double sign = 1;
	if (locelemind == 2) sign = -1;
	if (locelemind > elements.Length()) {UO() << "ERROR: inconsistency in SDActor2D::AddElementCqTLambda\n";}

	GetBody2D(locelemind).GetdPosdqT(loccoords(locelemind),dpdq);
	Vector2D F(ComputeForce(t));

	//UO() << "F=" << F << "\n";
	for (int i=1; i <= f.Length(); i++)
	{
		f(i) -= sign*(dpdq(i,1)*F(1)+dpdq(i,2)*F(2));
	}
};

Vector3D SDActor2D::GetRefPosD()	const 
{
	return GetBody2D(1).ToP3D(GetBody2D(1).GetPos2DD(loccoords(1)));
}


double SDActor2D::EvaluateActorForce()
{
	if( GetNInputs() == 1)
	{
		InputOutputElement* ioe = (InputOutputElement*)GetMBS()->GetElementPtr(inputs(1));
		return ioe->GetOutput(1); 
	}
	else
	{
		return fa;
	}
}


void SDActor2D::DrawElement() 
{
	Constraint::DrawElement();

	mbs->SetColor(GetCol());

	Vector3D pa = GetBody2D(1).ToP3D(GetBody2D(1).GetPos2DD(loccoords(1)));
	Vector3D pb;
	if (elements.Length() == 1)
		pb = GetBody2D(1).ToP3D(p_global);
	else
		pb = GetBody2D(2).ToP3D(GetBody2D(2).GetPos2DD(loccoords(2)));

	if (elements.Length() == 1) Swap(pa,pb);

	if (draw_dim.Y() == 0)
	{
		double l = (pa-pb).Norm();
		if (l < l0 && forcemode == 1) mbs->SetColor(Vector3D(0.95,0.95,0.95)); //forcemode = 1 -> no compression force

		mbs->DrawZyl(pa, pb, 0.5*draw_dim.X(),12);
	}
	else
	{
		//draw spring:
		double phi,phi2,xx;
		Vector3D vf;
		double steps = 150;
		double rots = draw_dim.Y();
		double r = draw_dim.X()/2;
		double rz = draw_dim.X()/10;
		int ztile = 6;

		vf = pa-pb;
		Vector3D vfn = vf;
		Vector3D vz,vy;
		vfn.SetNormalBasis(vz,vy);
		Vector3D pz1, pz2;

		double off = steps*0.1;
		for (xx = off; xx <= steps-off; xx++)
		{
			phi = xx*2*MY_PI/steps*rots;
			phi2 = (xx+1.2)*2*MY_PI/steps*rots;

			Vector3D pr1=cos(phi)*r*vz+sin(phi)*r*vy;
			Vector3D pr2=cos(phi2)*r*vz+sin(phi2)*r*vy;

			pz1=pb+xx/steps*vf+pr1;
			pz2=pb+(xx+1)/steps*vf+pr2;
			mbs->DrawZyl(pz1,pz2,rz,ztile);
			if (xx==off) mbs->DrawZyl(pb,pz1,rz,ztile);
		}
		mbs->DrawZyl(pa,pz2,rz,ztile);

		if (draw_dim.Z() != 0 && d != 0)
		{
			mbs->DrawZyl(pa, pa-0.5*l0*vfn, 0.5*draw_dim.Z(),12);
			mbs->SetColor(colgrey2);
			mbs->DrawZyl(pa-0.5*l0*vfn, pb, 0.5*draw_dim.Z()*0.6,12);

			mbs->DrawSphere(pa,0.65*draw_dim.Z(),12);
			mbs->DrawSphere(pb,0.65*draw_dim.Z(),12);
		}
	}

};




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Hydraulic Actuator Element - single acting cylinder
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Vector3D HydraulicActor::ComputeForce(double t) 
{
	//XG(1) is pressure in bar
	double l = 0;
	Vector3D v12 = GetBody3D(1).GetPos(loccoords(1)) - GetBody3D(2).GetPos(loccoords(2));
	l = v12.Norm();
	if (l!=0) v12*=1./l;;

	double damp = 0.; //damping does not improve a lot!!!
	double lp = (GetBody3D(1).GetVel(loccoords(1)) - GetBody3D(2).GetVel(loccoords(2)))*v12;

	return (-XG(1)*1e5*hyddata.A_zyl)*v12-damp*lp*1e5;
}

void HydraulicActor::EvalF(Vector& f, double t)
{
	double l = 0;
	Vector3D v12 = GetBody3D(1).GetPos(loccoords(1)) - GetBody3D(2).GetPos(loccoords(2));
	l = v12.Norm();
	if (l!=0) v12.Normalize();
	double lp = (GetBody3D(1).GetVel(loccoords(1)) - GetBody3D(2).GetVel(loccoords(2)))*v12;

	double s = l-hyddata.l_Zyl0; //position of piston, s=0 is closed

	double A_des = hyddata.R_P*(hyddata.s_end-s)+hyddata.R_D*(-lp);

	double A_v = Sgn(A_des)*(1.0-exp(-fabs(A_des)));

	double taccel = 0.;
	if (t < t_on)
	{
		A_v = 0;
	}
	else if (t < t_on+taccel) 
	{
		if (taccel == 0) taccel = 1e-10;
		double p = (t-t_on)/taccel; 
		A_v *= 0.5*(1-cos(p*MY_PI));
	}


	double s0 = s;
	if (s0 < 0) s0 = 0; //no negative volume ...

	double Vtot = hyddata.V_0+hyddata.A_zyl*s;
	if (Vtot < hyddata.V_0) Vtot = hyddata.V_0;

	double term;
	if (A_v > 0)
		term = Sgn(hyddata.p_s-XG(1))*hyddata.Q_n/sqrt(hyddata.p_n)*A_v*sqrt(fabs(hyddata.p_s-XG(1)));
	else
		term = Sgn(XG(1)-hyddata.p_T)*hyddata.Q_n/sqrt(hyddata.p_n)*A_v*sqrt(fabs(XG(1)-hyddata.p_T));

	double pp_new = hyddata.E_oil/Vtot*(-hyddata.A_zyl*lp+term);

	//if (this->GetMBS()->GetJacCol() == 0)	{		UO() << t << ", s=" << s << ", Av=" << A_v << ", sp=" << lp << ", p=" << XG(1) << ", pp=" << pp_new << "\n";}
	f(1) = pp_new;

};

void HydraulicActor::DrawElement() 
{
	Constraint::DrawElement();

	mbs->SetColor(GetCol());

	const Vector3D& pa = GetBody3D(1).GetPosD(loccoords(1));
	const Vector3D& pb = GetBody3D(2).GetPosD(loccoords(2));

	Vector3D vfn = pa-pb;
	vfn.Normalize();

	double l0 = hyddata.l_Zyl0-hyddata.l_0off;
	mbs->DrawZyl(pa, pa-l0*vfn, 0.5*draw_dim.X(),16);
	mbs->SetColor(colgrey2);
	mbs->DrawZyl(pa-l0*vfn, pb, 0.5*draw_dim.X()*0.6,16);

	mbs->DrawSphere(pa,0.5*draw_dim.Y(),12);
	mbs->DrawSphere(pb,0.5*draw_dim.Y(),12);
	/*
	char str[100];
	sprintf(str," s=%8.6g, p=%8.6g   ",
	(pb-pa).Norm()-hyddata.l_Zyl0, XGD(1));
	this->GetMBS()->GetRC()->PrintTextStruct(5,-1,str);
	*/

};




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Hydraulic Actuator Element - DOUBLE acting cylinder
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void HydraulicActorDA::DrawElement() 
{
	InputOutputElement::DrawElement();

	mbs->SetColor(GetCol());

	const Vector3D& pa = GetBody3D(1).GetPosD(loccoords(1));
	const Vector3D& pb = GetBody3D(2).GetPosD(loccoords(2));

	Vector3D vfn = pa-pb;
	vfn.Normalize();

	double l0d = hyddata.l_Zyl0-hyddata.l_0off;
	mbs->DrawZyl(pa, pa-l0d*vfn, 0.5*draw_dim.X(),16);
	mbs->SetColor(colgrey2);
	mbs->DrawZyl(pa-l0d*vfn, pb, 0.5*draw_dim.X()*0.6,16);

	mbs->DrawSphere(pa,0.5*draw_dim.Y(),12);
	mbs->DrawSphere(pb,0.5*draw_dim.Y(),12);

};


Vector3D HydraulicActorDA::ComputeForce(double t) const
{
	//XG(1) is pressure in bar
	double l = 0;
	Vector3D v12 = GetBody3D(1).GetPos(loccoords(1)) - GetBody3D(2).GetPos(loccoords(2));
	l = v12.Norm();
	if (l!=0) v12*=1./l;;

	double damp = 0.; //damping does not improve a lot!!!
	double lp = (GetBody3D(1).GetVel(loccoords(1)) - GetBody3D(2).GetVel(loccoords(2)))*v12;

	return (-XG(1)*1e5*hyddata.A_zyl+XG(2)*1e5*hyddata.A_zyl2)*v12-damp*lp*1e5;
}

void HydraulicActorDA::EvalF(Vector& f, double t)
{
	double l = 0;
	Vector3D v12 = GetBody3D(1).GetPos(loccoords(1)) - GetBody3D(2).GetPos(loccoords(2));
	l = v12.Norm();
	if (l!=0) v12.Normalize();
	double lp = (GetBody3D(1).GetVel(loccoords(1)) - GetBody3D(2).GetVel(loccoords(2)))*v12;

	double s = l-hyddata.l_Zyl0; //position of piston, s=0 is closed

	double A_des = hyddata.R_P*(hyddata.s_end-s)+hyddata.R_D*(-lp);

	double A_v = Sgn(A_des)*(1.0-exp(-fabs(A_des)));
	if( GetNInputs() == 1 )
	{
		InputOutputElement* ioe = (InputOutputElement*)GetMBS()->GetElementPtr(inputs(1));
		A_v = ioe->GetOutput(t, 1); 
		A_ext = A_v;
	}

	double s0 = s;
	if (s0 < 0) 
	{
		UO() << "Hydraulic-Warning: s negative ->leads to negative Volume! \n";
		s0 = 0; //no negative volume ...
	}

	double Vtot = hyddata.V_0+hyddata.A_zyl*s;
	double Vtot2 = hyddata.V_02-hyddata.A_zyl2*s;
	if (Vtot2 < 1e-6)
	{
		Vtot2 = 1e-6;
		UO() << "Warning: Vtot2 in hydraulic zylinder smaller than 0.001 liter!!!\n";
	}

	double term, term2;
	//cylinder 1:
	if (A_v > 0)
		term = Sgn(hyddata.p_s-XG(1))*hyddata.Q_n/sqrt(hyddata.p_n)*A_v*sqrt(fabs(hyddata.p_s-XG(1)));
	else
		term = Sgn(XG(1)-hyddata.p_T)*hyddata.Q_n/sqrt(hyddata.p_n)*A_v*sqrt(fabs(XG(1)-hyddata.p_T));

	double pp_new1 = hyddata.E_oil/Vtot*(-hyddata.A_zyl*lp+term);

	//cylinder 2:

	A_v = - A_v;
	if (A_v < 0)
		term2 = Sgn(XG(2)-hyddata.p_T)*hyddata.Q_n/sqrt(hyddata.p_n)*A_v*sqrt(fabs(XG(2)-hyddata.p_T));
	else
		term2 = Sgn(hyddata.p_s-XG(2))*hyddata.Q_n/sqrt(hyddata.p_n)*A_v*sqrt(fabs(hyddata.p_s-XG(2)));

	double pp_new2 = hyddata.E_oil/Vtot2*(hyddata.A_zyl2*lp+term2);

	f(1) = pp_new1;
	f(2) = pp_new2;

};






//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Hydraulic Actuator Element - DOUBLE acting cylinder
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


Vector2D HydraulicActorDA2D::ComputeForce(double t)
{
	double barfact = 1.e-6;
	//XG(1) is pressure in bar
	double damp = 0.; //damping does not improve a lot!!!

	double l = 0;
	Vector2D v12 = GetBody2D(1).GetPos2D(loccoords(1)) - GetBody2D(2).GetPos2D(loccoords(2));
	l = v12.Norm();
	if (l!=0) v12*=1./l;;

	double lp = (GetBody2D(1).GetVel2D(loccoords(1)) - GetBody2D(2).GetVel2D(loccoords(2)))*v12;

	return (-XG(1)*1e5*barfact*hyddata.A_zyl+XG(2)*1e5*barfact*hyddata.A_zyl2)*v12-damp*lp*1e5;
}

void HydraulicActorDA2D::EvalF(Vector& f, double t)
{
	double barfact = 1.e-6;
	//if (GetMBS()->GetTime() > 6) out = 1;

	double l = 0;
	Vector2D v12 = GetBody2D(1).GetPos2D(loccoords(1)) - GetBody2D(2).GetPos2D(loccoords(2));
	l = v12.Norm();
	if (l!=0) v12*=1./l;;

	//if (out) UO() << "l=" << l << "\n";

	double lp = (GetBody2D(1).GetVel2D(loccoords(1)) - GetBody2D(2).GetVel2D(loccoords(2)))*v12;

	//if (out) UO() << "lp=" << lp << "\n";

	double s = l-hyddata.l_Zyl0; //position of piston, s=0 is closed

	//if (out) UO() << "s=" << s << "\n";

	//TODO: same as 3D
	if( GetNInputs() == 1 )
	{
		InputOutputElement* ioe = (InputOutputElement*)GetMBS()->GetElementPtr(inputs(1));
		A_ext = ioe->GetOutput(t, 1); 
	}
	double A_des = A_ext; // value from connected InputOutputElement //old: hyddata.R_P*(hyddata.s_end-s)+hyddata.R_D*(-lp);
	double A_v = Sgn(A_des)*(1.0-exp(-fabs(A_des)));

	double taccel = 0.;
	if (t < t_on) A_v = 0;
	else if (taccel != 0 && t < t_on+taccel) 
	{
		double p = (t-t_on)/taccel; 
		A_v *= 0.5*(1-cos(p*MY_PI));
	}
	//A_v = 1;
	if (A_v < 0) A_v = 0;

	//if (out) UO() << "A_v=" << A_v << "\n";

	//double s0 = s;
	//if (s0 < 0) s0 = 0; //no negative volume ...

	double Vtot = hyddata.V_0+hyddata.A_zyl*s;
	if (Vtot < hyddata.V_0) Vtot = hyddata.V_0;

	double Vtot2 = hyddata.V_02-hyddata.A_zyl2*s;
	if (Vtot2 < 1e-6)
	{
		Vtot2 = 1e-6;
		UO() << "Warning: Vtot2 in hydraulic zylinder smaller than 0.001 liter!!!\n";
	}

	double term, term2;
	//cylinder 1:

	if (A_v > 0)
		term = Sgn(hyddata.p_s-XG(1)*barfact)*hyddata.Q_n/sqrt(hyddata.p_n)*A_v*sqrt(fabs(hyddata.p_s-XG(1)*barfact));
	else
		term = Sgn(XG(1)*barfact-hyddata.p_T)*hyddata.Q_n/sqrt(hyddata.p_n)*A_v*sqrt(fabs(XG(1)*barfact-hyddata.p_T));

	double pp_new1 = hyddata.E_oil/Vtot*(-hyddata.A_zyl*lp+term);

	//cylinder 2:

	A_v = - A_v;
	if (A_v < 0)
		term2 = Sgn(XG(2)*barfact-hyddata.p_T)*hyddata.Q_n/sqrt(hyddata.p_n)*A_v*sqrt(fabs(XG(2)*barfact-hyddata.p_T));
	else
		term2 = Sgn(hyddata.p_s-XG(2)*barfact)*hyddata.Q_n/sqrt(hyddata.p_n)*A_v*sqrt(fabs(hyddata.p_s-XG(2)*barfact));

	double pp_new2 = hyddata.E_oil/Vtot2*(hyddata.A_zyl2*lp+term2);

/*
	if (A_v > 0)  
		term = Sgn(hyddata.p_s-XG(1)*barfact)*hyddata.Q_n/sqrt(hyddata.p_n)*A_v*sqrt(fabs(hyddata.p_s-XG(1)*barfact));
	else
		term = Sgn(XG(1)*barfact-hyddata.p_T)*hyddata.Q_n/sqrt(hyddata.p_n)*A_v*sqrt(fabs(XG(1)*barfact-hyddata.p_T));

	double pp_new1 = hyddata.E_oil/Vtot*(-hyddata.A_zyl*lp+term);

	//UO() << "Eoil = " << hyddata.E_oil << "\n";

	//cylinder 2:
	A_v = - A_v;
	if (A_v < 0)  
		term2 = Sgn(XG(2)*barfact-hyddata.p_T)*hyddata.Q_n/sqrt(hyddata.p_n)*A_v*sqrt(fabs(XG(2)*barfact-hyddata.p_T));
	else
		term2 = Sgn(hyddata.p_s-XG(2)*barfact)*hyddata.Q_n/sqrt(hyddata.p_n)*A_v*sqrt(fabs(hyddata.p_s-XG(2)*barfact));

	double pp_new2 = hyddata.E_oil/Vtot2*(hyddata.A_zyl2*lp+term2);
*/

	//if (out) UO() << "pp1=" << pp_new1 << "\n";
	//if (out) UO() << "pp2=" << pp_new2 << "\n";

	//if (out) UO() << "p1=" << XG(1) << "\n";
	//if (out) UO() << "p2=" << XG(2) << "\n";

	f(1) = pp_new1;
	f(2) = pp_new2;

};

void HydraulicActorDA2D::DrawElement() 
{
	Constraint::DrawElement();

	mbs->SetColor(GetCol());

	const Vector3D& pa = GetBody2D(1).GetPosD(GetBody2D(1).ToP3D(loccoords(1)));
	const Vector3D& pb = GetBody2D(2).GetPosD(GetBody2D(2).ToP3D(loccoords(2)));

	Vector3D vfn = pa-pb;
	vfn.Normalize();

	double l0 = hyddata.l_Zyl0-hyddata.l_0off;
	mbs->DrawZyl(pa, pa-l0*vfn, 0.5*draw_dim.X(),16);
	mbs->SetColor(colgrey2);
	mbs->DrawZyl(pa-l0*vfn, pb, 0.5*draw_dim.X()*0.6,16);

	mbs->DrawSphere(pa,0.5*draw_dim.Y(),12);
	mbs->DrawSphere(pb,0.5*draw_dim.Y(),12);

};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Hydraulic Actuator Element - DOUBLE acting cylinder 2D NEW
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Vector2D HydraulicActorDA2Dn::ComputeForce(double t) 
{
	//XG(1), XG(2) are pressures in bar
	// 1 bar = 10^5 N/m²
	
	double l = 0;
	Vector2D v12 = GetBody2D(1).GetPos2D(loccoords(1)) - GetBody2D(2).GetPos2D(loccoords(2));
	l = v12.Norm();
	if (l!=0) v12*=1./l;;

	double lp = (GetBody2D(1).GetVel2D(loccoords(1)) - GetBody2D(2).GetVel2D(loccoords(2)))*v12;

	return (-XG(1)*1e5*hyddata.A_zyl+XG(2)*1e5*hyddata.A_zyl2)*v12;
}

void HydraulicActorDA2Dn::EvalF(Vector& f, double t)
{
	int out = 0;
	if (GetMBS()->GetTime() > 100) out = 1;

	double l = 0;
	Vector2D v12 = GetBody2D(1).GetPos2D(loccoords(1)) - GetBody2D(2).GetPos2D(loccoords(2));
	l = v12.Norm();
	if (l!=0) v12*=1./l;;

	if (out) UO() << "l=" << l << "\n";

	double lp = (GetBody2D(1).GetVel2D(loccoords(1)) - GetBody2D(2).GetVel2D(loccoords(2)))*v12;

	if (out) UO() << "lp=" << lp << "\n";

	double s = l-hyddata.l_Zyl0; //position of piston, s=0 is closed

	if (out) UO() << "s=" << s << "\n";

	if( GetNInputs() == 1 )
	{
		InputOutputElement* ioe = (InputOutputElement*)GetMBS()->GetElementPtr(inputs(1));
		A_ext = ioe->GetOutput(t, 1); 
	}

	double A_v = A_ext; // value from connected InputOutputElement //old: hyddata.R_P*(hyddata.s_end-s)+hyddata.R_D*(-lp);

	//if (out) UO() << "A_v=" << A_v << "\n";
	
	if (A_v < 0) A_v = 0;

	

	//double s0 = s;
	//if (s0 < 0) s0 = 0; //no negative volume ...

	double Vtot = hyddata.V_0+hyddata.A_zyl*s;
	if (Vtot < hyddata.V_0) Vtot = hyddata.V_0;

	double Vtot2 = hyddata.V_02-hyddata.A_zyl2*s;
	if (Vtot2 < 1e-6)
	{
		Vtot2 = 1e-6;
		UO() << "Warning: Vtot2 in hydraulic zylinder smaller than 0.001 liter!!!\n";
	}

	double term, term2;
	//cylinder 1:

	if (A_v > 0)
		term = Sgn(hyddata.p_s-XG(1))*hyddata.Q_n/sqrt(hyddata.p_n)*A_v*sqrt(fabs(hyddata.p_s-XG(1)));
	else
		term = Sgn(XG(1)-hyddata.p_T)*hyddata.Q_n/sqrt(hyddata.p_n)*A_v*sqrt(fabs(XG(1)-hyddata.p_T));

	double pp_new1 = hyddata.E_oil/Vtot*(-hyddata.A_zyl*lp+term);

	//cylinder 2:

	A_v = - A_v;
	if (A_v < 0)
		term2 = Sgn(XG(2)-hyddata.p_T)*hyddata.Q_n/sqrt(hyddata.p_n)*A_v*sqrt(fabs(XG(2)-hyddata.p_T));
	else
		term2 = Sgn(hyddata.p_s-XG(2))*hyddata.Q_n/sqrt(hyddata.p_n)*A_v*sqrt(fabs(hyddata.p_s-XG(2)));

	double pp_new2 = hyddata.E_oil/Vtot2*(hyddata.A_zyl2*lp+term2);


	if (out) UO() << "pp1=" << pp_new1 << "\n";
	if (out) UO() << "pp2=" << pp_new2 << "\n";

	//if (out) UO() << "p1=" << XG(1) << "\n";
	//if (out) UO() << "p2=" << XG(2) << "\n";

	f(1) = pp_new1;
	f(2) = pp_new2;

};

void HydraulicActorDA2Dn::DrawElement() 
{
	Constraint::DrawElement();

	mbs->SetColor(GetCol());

	const Vector3D& pa = GetBody2D(1).GetPosD(GetBody2D(1).ToP3D(loccoords(1)));
	const Vector3D& pb = GetBody2D(2).GetPosD(GetBody2D(2).ToP3D(loccoords(2)));

	Vector3D vfn = pa-pb;
	vfn.Normalize();

	double l0 = hyddata.l_Zyl0-hyddata.l_0off;
	mbs->DrawZyl(pa, pa-l0*vfn, 0.5*draw_dim.X(),16);
	mbs->SetColor(colgrey2);
	mbs->DrawZyl(pa-l0*vfn, pb, 0.5*draw_dim.X()*0.6,16);

	mbs->DrawSphere(pa,0.5*draw_dim.Y(),12);
	mbs->DrawSphere(pb,0.5*draw_dim.Y(),12);

	/*
	InputOutputElement* ioe = (InputOutputElement*)GetMBS()->GetElementPtr(IOGetElementNum(1));	
	Vector2D po2D = ioe->GetOutputPosD(IOGetOutputNum(1));
	Vector3D po = ioe->ToP3D(ioe->GetOutputPosD(IOGetOutputNum(1)));
	GetMBS()->MyDrawLine(po, pa, 1., colblack);
*/
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// InertialLinearSpringDamper
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void InertialLinearSpringDamper::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	InputOutputElement::GetElementData(edc);

	ElementData ed;
	ed.SetDouble(draw_dim.X(), "Draw_size"); ed.SetToolTipText("General drawing size of constraint"); edc.Add(ed);
	ed.SetInt((int)draw_dim.Z(), "Draw_resolution", 0, 1000); ed.SetToolTipText("Number of quads to approximate round geometries, suggested = 16"); edc.Add(ed);
	

	Matrix kmat(3,3),dmat(3,3);
	for(int i=1; i<=kmat.Getrows(); i++)
	{
		for(int j=1; j<=kmat.Getcols(); j++)
		{
			kmat(i,j) = k(i,j);
			dmat(i,j) = d(i,j);
		}
	}
	ed.SetMatrix(kmat.GetMatPtr(), kmat.Getrows(), kmat.Getcols(), "spring_stiffness"); ed.SetToolTipText("Spring stiffness w.r.t. inertial coordinate system"); edc.Add(ed);
	ed.SetMatrix(dmat.GetMatPtr(), dmat.Getrows(), dmat.Getcols(), "damping_coefficient"); ed.SetToolTipText("Damping coefficient w.r.t. inertial coordinate system"); edc.Add(ed);
	

	if (elements.Length()==1)
	{
		SetElemDataVector3D(edc, loccoords(1), "Local_joint_pos_body"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector3D(edc, p_global, "Global_joint_pos"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
	} else
	{
		SetElemDataVector3D(edc, loccoords(1), "Local_joint_pos_body1"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector3D(edc, loccoords(2), "Local_joint_pos_body2"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
	}
	
	Vector3D v = init_disp;
	ed.SetVector3D(v(1), v(2), v(3), "Initial_spring_displacement"); edc.Add(ed);
	int pos = edc.Find("Initial_spring_displacement");
	if (pos) edc.Get(pos).SetLocked(1);

	Vector3D vp;
	if (elements.Length() == 2) vp = GetBody3D(1).GetVel(loccoords(1)) - GetBody3D(2).GetVel(loccoords(2));
	else vp = GetBody3D(1).GetVel(loccoords(1));
	ed.SetVector3D(vp(1), vp(2), vp(3), "Initial_spring_velocity"); edc.Add(ed);
	pos = edc.Find("Initial_spring_velocity");
	if (pos) edc.Get(pos).SetLocked(1);
}

int InertialLinearSpringDamper::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = InputOutputElement::SetElementData(edc);

	GetElemDataDouble(mbs, edc, "Draw_size", draw_dim.X(), 0);
	int dd;
	if (GetElemDataInt(mbs, edc, "Draw_resolution", dd, 0)) draw_dim.Z() = dd;
	

	Matrix kmat(3,3),dmat(3,3);
	GetElemDataMatrix(GetMBS(), edc, "spring_stiffness", kmat, 1);
	GetElemDataMatrix(GetMBS(), edc, "damping_coefficient", kmat, 1);
	for(int i=1; i<=3; i++)
	{
		for(int j=1; j<=3; j++)
		{
			k(i,j) = kmat(i,j);
			d(i,j) = dmat(i,j);
		}
	}

	if (elements.Length()==1)
	{
		GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body", loccoords(1));
		GetElemDataVector3D(GetMBS(), edc, "Global_joint_pos", p_global);
	} else
	{
		GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body1", loccoords(1));
		GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body2", loccoords(2));
	}

	return rv;
}


Vector3D InertialLinearSpringDamper::ComputeForce(double t)  const
{
	Vector3D vp(0.,0.,0.);
	Vector3D v(0.,0.,0.);

	if (IsGroundJoint())
	{
		v = GetBody3D(1).GetPos(loccoords(1)) - p_global;
		
		if (d.AbsNorm() != 0) 
			vp = (GetBody3D(1).GetVel(loccoords(1)));
	}
	else
	{
		v = GetBody3D(1).GetPos(loccoords(1)) - GetBody3D(2).GetPos(loccoords(2));

		if (d.AbsNorm() != 0) 
		{
			vp = GetBody3D(1).GetVel(loccoords(1)) - GetBody3D(2).GetVel(loccoords(2));
		}
	}

	// k, d ... stiffness and damping matrix
	Vector3D delta = (v-l0);
	return k*delta + d*vp;
}

double InertialLinearSpringDamper::GetActorForce(double computation_time, int dir) const
{
	if (dir >=1 && dir <= 3) return (ComputeForce(computation_time))(dir);
	else return (ComputeForce(computation_time)).Norm();
}

// function no longer needed, calculation now in EvalF2
void InertialLinearSpringDamper::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	return; // (AD) function no longer needed, calculation now in EvalF2
// old code
	////-C_q^T \lambda = [dC/dx; dC/dy; dC/dphi]^T [\lambda1;\lambda2] 
	////e.g.: C = p_ref+R*p_loc
	//double sign = 1;
	//if (locelemind == 2) sign = -1;

	//Vector3D F3D(ComputeForce(t));


	//GetBody3D(locelemind).GetdPosdqT(loccoords(locelemind),dpdq);
	////UO() << "dpdq=" << dpdq;

	//for (int i=1; i <= f.Length(); i++)
	//{
	//	f(i) -= sign*(dpdq(i,1)*F3D.X()+dpdq(i,2)*F3D.Y()+dpdq(i,3)*F3D.Z());
	//}
}

int InertialLinearSpringDamper::SOS() const
{
	int nsos = 0;
	for (int i=1; i <= NE(); i++)
	{
		nsos += GetElem(i).SOS();
	}
	return nsos;
}

void InertialLinearSpringDamper::EvalF2(Vector& f, double t)
{
	InputOutputElement::EvalF2(f,t);

	for (int i=1; i <= NE(); i++)
	{
		double sign = 1.0;
		int offset = 0;
		if (i==2) 
		{
			sign = -1.0;
		  offset = GetBody3D(1).SOS();
		}

		Vector3D F3D(ComputeForce(t));
		GetBody3D(i).GetdPosdqT(loccoords(i),dpdq);

		for (int j=1; j<=GetBody3D(i).SOS(); j++)
		{
			f(j+offset) -= sign*(dpdq(j,1)*F3D.X() + dpdq(j,2)*F3D.Y() + dpdq(j,3)*F3D.Z());
		}
	}
}

void InertialLinearSpringDamper::LinkToElements()
{
		LTGreset();

	// add all SOS dofs from the elements
	//Position(first SOS) 
	for (int k=1; k <= NE(); k++)
	{
		for (int i=1; i <= GetBody3D(k).SOS(); i++)
		{
			AddLTG(GetBody3D(k).LTG(i));
		}
	}
	//and Velocity (second SOS):
	for (int k=1; k <= NE(); k++)
	{
		for (int i=1; i <= GetBody3D(k).SOS(); i++)
		{
			AddLTG(GetElem(k).LTG(i+GetBody3D(k).SOS()));
		}
	}
}


Vector3D InertialLinearSpringDamper::GetRefPosD()	const 
{
	return GetBody3D(1).GetPosD(loccoords(1));
}

void InertialLinearSpringDamper::DrawElement() 
{
	InputOutputElement::DrawElement();

	int res = (int)draw_dim.Z();
	if (res < 2) res = 3;

	if (draw_dim.X() != 0)
	{
		mbs->SetColor(GetCol());
		Vector3D p1 = GetBody3D(1).GetPosD(loccoords(1));
		mbs->DrawSphere(p1,0.5*draw_dim.X(),res);

		mbs->SetColor(Vector3D(0.8,0.1,0.1));
		Vector3D p2;
		if (elements.Length()==1)
		{
			p2 = p_global;
		}
		else
		{
			p2 = GetBody3D(2).GetPosD(loccoords(2));
		}
		mbs->DrawSphere(p2,0.5*draw_dim.X(),res-1);
		if (draw_dim.Y() > 0) mbs->DrawZyl(p1,p2, draw_dim.Y(), res);
	}
};


Vector3D InertialSpringDamperNL::ComputeForce(double t) const
{
	Vector3D vp(0.,0.,0.);
	Vector3D v(0.,0.,0.);

	if (IsGroundJoint())
	{
		v = GetBody3D(1).GetPos(loccoords(1)) - p_global;
		
		if (d.AbsNorm() != 0) 
			vp = (GetBody3D(1).GetVel(loccoords(1)));
	}
	else
	{
		v = GetBody3D(1).GetPos(loccoords(1)) - GetBody3D(2).GetPos(loccoords(2));

		if (d.AbsNorm() != 0) 
		{
			vp = GetBody3D(1).GetVel(loccoords(1)) - GetBody3D(2).GetVel(loccoords(2));
		}
	}

	// k, d ... stiffness and damping matrix
	Vector3D delta = (v-l0);
	if (usestiffmfunc)
	{
	
 
		//int push = 0;
		//if (delta.Z() > 0)
		//{
		//	push = 0.;
		//}
		//else
		//{
		//	//GetMBS()->UO() << "Time = " << GetMBS()->GetTime() << ", deltaZ = " << delta.Z() << "\n";
		//	push = 1;
		//}

		double actstiffX = XData(7);
		double actstiffY = XData(8);
		double actstiffZ = XData(9);
		
		Vector3D force = Vector3D(actstiffX*delta.X(),actstiffY*delta.Y(),actstiffZ*delta.Z())+ d*vp;
		//Vector3D force = Vector3D(push*actstiffX*delta.X(),push*actstiffY*delta.Y(),actstiffZ*delta.Z())+ push*d*vp;
		//GetMBS()->UO() << "Delta = " << delta << " Kraft = " << force << "\n";
		return force;
	}
	else
	{
		GetMBS()->UO() << "Error: No stiffness Math-function definied for NL-IntertialSpringDamper! \n";
		return k*delta + d*vp;
	}
}


double InertialSpringDamperNL::PostNewtonStep(double t)
{
	if (usestiffmfunc)
	{
		if (TMFpiecewiselinear || TMFpiecewiseconst)
		{
			Vector3D vp(0.,0.,0.);
			Vector3D v(0.,0.,0.);

			if (IsGroundJoint())
			{
				v = GetBody3D(1).GetPos(loccoords(1)) - p_global;

				if (d.AbsNorm() != 0) 
					vp = (GetBody3D(1).GetVel(loccoords(1)));
			}
			else
			{
				v = GetBody3D(1).GetPos(loccoords(1)) - GetBody3D(2).GetPos(loccoords(2));

				if (d.AbsNorm() != 0) 
				{
					vp = GetBody3D(1).GetVel(loccoords(1)) - GetBody3D(2).GetVel(loccoords(2));
				}
			}

			Vector3D delta = (v-l0);

			XData(7) = stiffmfuncX.Evaluate(delta.X());
			XData(8) = stiffmfuncY.Evaluate(delta.Y());
			XData(9) = stiffmfuncZ.Evaluate(delta.Z());

			double errorx = 0; double errory = 0; double errorz = 0;

			int index1 = XData(1);
			int index2 = stiffmfuncX.FindIndexPiecewise(delta.X());
			double delta1 = XData(2);
			double delta2 = delta.X();
			XData(1) = index2;
			XData(2) = delta2;
			
			if ((index1!=index2) && (index1!=0) && (index2!=0))
			{
				Matrix xdata;
				int mode = stiffmfuncX.GetFuncMode();
				stiffmfuncX.GetData(mode,xdata);

				double x1 = xdata.GetColVec(1).Get(index1);
				double x2 = xdata.GetColVec(1).Get(index2);
				double x = 0;
				if ((Minimum(delta1,delta2) <= x1) && (x1 <= Maximum(delta1,delta2)))
				{
					x = x1;
				}
				else
				{
					x = x2;
				}

				double b = delta2 - x;
				double a = x - delta1;
				double dt = GetMBS()->GetStepSize();
				double vel = (delta2-delta1)/dt;
				double newstep_t = a/vel;

				errorx = b;
				
				if (newstep_t < GetMBS()->GetStepRecommendation() && newstep_t < GetMBS()->GetStepSize()) 
				{
					//GetMBS()->UO() << "Discont.Step Change! : error_x = " << errorx << ", time = " << mbs->GetTime() << ", oldstep = " << mbs->GetStepSize() << ", newstep = " << newstep_t << "\n";
					GetMBS()->SetStepRecommendation(newstep_t);	
				}

			}

			index1 = XData(3);
			index2 = stiffmfuncY.FindIndexPiecewise(delta.Y());
			delta1 = XData(4);
			delta2 = delta.Y();
			XData(3) = index2;
			XData(4) = delta2;
			
			if ((index1!=index2) && (index1!=0) && (index2!=0))
			{
				Matrix ydata;
				int mode = stiffmfuncY.GetFuncMode();
				stiffmfuncY.GetData(mode,ydata);

				double y1 = ydata.GetColVec(1).Get(index1);
				double y2 = ydata.GetColVec(1).Get(index2);
				double y = 0;
				if ((Minimum(delta1,delta2) <= y1) && (y1 <= Maximum(delta1,delta2)))
				{
					y = y1;
				}
				else
				{
					y = y2;
				}

				double b = delta2 - y;
				double a = y - delta1;
				double dt = GetMBS()->GetStepSize();
				double vel = (delta2-delta1)/dt;
				double newstep_t = a/vel;

				errory = b;
				
				if (newstep_t < GetMBS()->GetStepRecommendation() && newstep_t < GetMBS()->GetStepSize()) 
				{
					//GetMBS()->UO() << "Discont.Step Change! : error_y = " << errory << ", time = " << mbs->GetTime() << ", oldstep = " << mbs->GetStepSize() << ", newstep = " << newstep_t << "\n";
					GetMBS()->SetStepRecommendation(newstep_t);	
				}

			}

			index1 = XData(5);
			index2 = stiffmfuncZ.FindIndexPiecewise(delta.Z());
			delta1 = XData(6);
			delta2 = delta.Z();
			XData(5) = index2;
			XData(6) = delta2;
			
			if ((index1!=index2) && (index1!=0) && (index2!=0))
			{
				Matrix Zdata;
				int mode = stiffmfuncZ.GetFuncMode();
				stiffmfuncZ.GetData(mode,Zdata);

				double z1 = Zdata.GetColVec(1).Get(index1);
				double z2 = Zdata.GetColVec(1).Get(index2);
				double z = 0;
				if ((Minimum(delta1,delta2) <= z1) && (z1 <= Maximum(delta1,delta2)))
				{
					z = z1;
				}
				else
				{
					z = z2;
				}

				double b = delta2 - z;
				double a = z - delta1;
				double dt = GetMBS()->GetStepSize();
				double vel = (delta2-delta1)/dt;
				double newstep_t = a/vel;

				errorz = b;
				
				if (newstep_t < GetMBS()->GetStepRecommendation() && newstep_t < GetMBS()->GetStepSize()) 
				{
					//GetMBS()->UO() << "Discont.Step Change! : error_z = " << errorz << ", time = " << mbs->GetTime() << ", oldstep = " << mbs->GetStepSize() << ", newstep = " << newstep_t << "\n";
					GetMBS()->SetStepRecommendation(newstep_t);	
				}

			}
			return errorx+errory+errorz;
		}
		else
		{
			return 0;
		}
	}
	else
	{
		return 0;
	}
}







//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// InertialLinearSpringDamper2D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void InertialLinearSpringDamper2D::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Constraint::GetElementData(edc);

	ElementData ed;
	ed.SetDouble(draw_dim.X(), "Draw_size"); ed.SetToolTipText("General drawing size of constraint"); edc.Add(ed);
	ed.SetInt((int)draw_dim.Z(), "Draw_resolution", 0, 1000); ed.SetToolTipText("Number of quads to approximate round geometries, suggested = 16"); edc.Add(ed);
	

	Matrix kmat(2,2),dmat(2,2);
	for(int i=1; i<=kmat.Getrows(); i++)
	{
		for(int j=1; j<=kmat.Getcols(); j++)
		{
			kmat(i,j) = k(i,j);
			dmat(i,j) = d(i,j);
		}
	}
	ed.SetMatrix(kmat.GetMatPtr(), kmat.Getrows(), kmat.Getcols(), "spring_stiffness"); ed.SetToolTipText("Spring stiffness w.r.t. inertial coordinate system"); edc.Add(ed);
	ed.SetMatrix(dmat.GetMatPtr(), dmat.Getrows(), dmat.Getcols(), "damping_coefficient"); ed.SetToolTipText("Damping coefficient w.r.t. inertial coordinate system"); edc.Add(ed);
	
	if (elements.Length()==1)
	{
		SetElemDataVector2D(edc, loccoords(1), "Local_joint_pos_body"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector2D(edc, p_global, "Global_joint_pos"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
	} 
	else
	{
		SetElemDataVector2D(edc, loccoords(1), "Local_joint_pos_body1"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector2D(edc, loccoords(2), "Local_joint_pos_body2"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
	}

	Vector2D v = init_disp;
	ed.SetVector2D(v(1), v(2), "Initial_spring_displacement"); edc.Add(ed);
	int pos = edc.Find("Initial_spring_displacement");
	if (pos) edc.Get(pos).SetLocked(1);

	Vector2D vp;
	if (elements.Length() == 2) vp = GetBody2D(1).GetVel2D(loccoords(1)) - GetBody2D(2).GetVel2D(loccoords(2));
	else vp = GetBody2D(1).GetVel2D(loccoords(1));
	ed.SetVector2D(vp(1), vp(2), "Initial_spring_velocity"); edc.Add(ed);
	pos = edc.Find("Initial_spring_velocity");
	if (pos) edc.Get(pos).SetLocked(1);
}

int InertialLinearSpringDamper2D::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Constraint::SetElementData(edc);

	GetElemDataDouble(mbs, edc, "Draw_size", draw_dim.X(), 0);
	int dd;
	if (GetElemDataInt(mbs, edc, "Draw_resolution", dd, 0)) draw_dim.Z() = dd;
	

	Matrix kmat(2,2),dmat(2,2);
	GetElemDataMatrix(GetMBS(), edc, "spring_stiffness", kmat, 1);
	GetElemDataMatrix(GetMBS(), edc, "damping_coefficient", dmat, 1);
	for(int i=1; i<=2; i++)
	{
		for(int j=1; j<=2; j++)
		{
			k(i,j) = kmat(i,j);
			d(i,j) = dmat(i,j);
		}
	}

	if (elements.Length()==1)
	{
		GetElemDataVector2D(GetMBS(), edc, "Local_joint_pos_body", loccoords(1));
		GetElemDataVector2D(GetMBS(), edc, "Global_joint_pos", p_global);
	} 
	else
	{
		GetElemDataVector2D(GetMBS(), edc, "Local_joint_pos_body1", loccoords(1));
		GetElemDataVector2D(GetMBS(), edc, "Local_joint_pos_body2", loccoords(2));
	}
	
	return rv;
}

void InertialLinearSpringDamper2D::Initialize() 
{
	if (!autocompute2ndcoord) return;

	if (elements.Length() == 1) // ground joint: 
	{
		Matrix3D A = GetBody2D(1).GetRotMatrix2D();
		Vector2D cg = GetBody2D(2).GetRefPos2D();
		
		if(autocompute2ndcoord == 1) // calculate global from local1
		{
//			p_global = A * loccoords(1) + cg ;
			p_global = A * loccoords(1) + cg - init_disp; //(AD) - hack initial displacement from rest position
		}
		if(autocompute2ndcoord == 2) // calculate local from global
		{
//			loccoords(1) = A.GetTp() * (p_global - cg);
			loccoords(1) = A.GetTp() * (p_global - cg + init_disp); //(AD) - hack initial displacement from rest position
		}
	}
	else // regular 2-element joint
	{
		Matrix3D A1 = GetBody2D(1).GetRotMatrix2D();
		Matrix3D A2 = GetBody2D(2).GetRotMatrix2D();
		Vector2D cg1 = GetBody2D(1).GetRefPos2D();
		Vector2D cg2 = GetBody2D(2).GetRefPos2D();
		
		if(autocompute2ndcoord == 1) // calculate global & local2 from local1
		{
			p_global = A1 * loccoords(1) + cg1;
//			loccoords(2) = A2.GetTp() * (p_global - cg2);
			loccoords(2) = A2.GetTp() * (p_global - cg2 + init_disp); //(AD) - hack initial displacement from rest position
		}
		if(autocompute2ndcoord == 2) // calculate local1 & local2 from global
		{
			loccoords(1) = A1.GetTp() * (p_global - cg1);
//			loccoords(2) = A2.GetTp() * (p_global - cg2);
			loccoords(2) = A2.GetTp() * (p_global - cg2 + init_disp); //(AD) - hack initial displacement from rest position
		}
	}
//  autocompute2ndcoord = 0; // compute only once
};

Vector2D InertialLinearSpringDamper2D::ComputeForce(double t)  const
{
	Vector2D vp(0.,0.);
	Vector2D v(0.,0.);

	if (IsGroundJoint())
	{
		v = GetBody2D(1).GetPos2D(loccoords(1)) - p_global;
		
		if (d.AbsNorm() != 0) 
			vp = (GetBody2D(1).GetVel2D(loccoords(1)));
	}
	else
	{
		Vector2D p1 = GetBody2D(1).GetPos2D(loccoords(1));
		Vector2D p2 = GetBody2D(2).GetPos2D(loccoords(2));

		v = GetBody2D(1).GetPos2D(loccoords(1)) - GetBody2D(2).GetPos2D(loccoords(2));

		if (d.AbsNorm() != 0) 
			vp = GetBody2D(1).GetVel2D(loccoords(1)) - GetBody2D(2).GetVel2D(loccoords(2));
	}

	// k, d ... stiffness and damping matrix
	return k*v + d*vp;
}

double InertialLinearSpringDamper2D::GetActorForce(double computation_time, int dir) const
{
	if (dir >=1 && dir <= 2) return (ComputeForce(computation_time))(dir);
	else return (ComputeForce(computation_time)).Norm();
}

// function no longer needed, calculation now in EvalF2
void InertialLinearSpringDamper2D::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	return; // (AD) function no longer needed, calculation now in EvalF2
// old code
	////-C_q^T \lambda = [dC/dx; dC/dy; dC/dphi]^T [\lambda1;\lambda2] 
	////e.g.: C = p_ref+R*p_loc
	//double sign = 1;
	//if (locelemind == 2) sign = -1;

	//Vector2D F2D(ComputeForce(t));

	//GetBody2D(locelemind).GetdPosdqT(loccoords(locelemind),dpdq);
	////UO() << "dpdq=" << dpdq;

	//for (int i=1; i <= f.Length(); i++)
	//{
	//	f(i) -= sign*(dpdq(i,1)*F2D.X()+dpdq(i,2)*F2D.Y());
	//}
};


int InertialLinearSpringDamper2D::SOS() const
{
// return;
	int nsos = 0;
	for (int i=1; i <= NE(); i++)
	{
		nsos += GetElem(i).SOS();
	}
	return nsos;
}

void InertialLinearSpringDamper2D::LinkToElements()
{
	LTGreset();

	// add all SOS dofs from the elements
	//Position(first SOS) 
	for (int k=1; k <= NE(); k++)
	{
		for (int i=1; i <= GetBody2D(k).SOS(); i++)
		{
			AddLTG(GetBody2D(k).LTG(i));
		}
	}
	//and Velocity (second SOS):
	for (int k=1; k <= NE(); k++)
	{
		for (int i=1; i <= GetBody2D(k).SOS(); i++)
		{
			AddLTG(GetElem(k).LTG(i+GetBody2D(k).SOS()));
		}
	}
}

void InertialLinearSpringDamper2D::EvalF2(Vector &f, double t)
{
	Constraint::EvalF2(f,t);

	for (int i=1; i <= NE(); i++)
	{
		double sign = 1.0;
		int offset = 0;
		if (i==2) 
		{
			sign = -1.0;
		  offset = GetBody2D(1).SOS();
		}

		Vector2D F2D(ComputeForce(t));
		GetBody2D(i).GetdPosdqT(loccoords(i),dpdq);

		for (int j=1; j<=GetBody2D(i).SOS(); j++)
		{
			f(j+offset) -= sign*(dpdq(j,1)*F2D.X() + dpdq(j,2)*F2D.Y());
		}
	}
}


Vector2D InertialLinearSpringDamper2D::GetRefPos2DD()	const 
{
	return GetBody2D(1).GetPos2DD(loccoords(1));
}

void InertialLinearSpringDamper2D::DrawElement() 
{
	Constraint::DrawElement();

	int res = (int)draw_dim.Z();
	if (res < 2) res = 3;

	if (draw_dim.X() != 0)
	{
		mbs->SetColor(GetCol());
		Vector3D p1 = GetBody2D(1).ToP3D(GetBody2D(1).GetPos2DD(loccoords(1)));

		mbs->DrawSphere(p1,0.5*draw_dim.X(),res);

		mbs->SetColor(Vector3D(0.8,0.1,0.1));
		Vector3D p2;
		if (elements.Length()==1)
		{
			p2 = GetBody2D(1).ToP3D(p_global);
		}
		else
		{
			p2 = GetBody2D(2).ToP3D(GetBody3D(2).GetPos2DD(loccoords(2)));
		}

		mbs->DrawSphere(p2,0.5*draw_dim.X(),res-1);

		if (draw_dim.Y() > 0) mbs->DrawZyl(p1,p2, draw_dim.Y(), res);
	}
};








//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// LinearRotationalSpringDamper
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void LinearRotationalSpringDamper::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Constraint::GetElementData(edc);

	ElementData ed;
	ed.SetDouble(draw_dim.X(), "Draw_size"); ed.SetToolTipText("General drawing size of constraint"); edc.Add(ed);
	ed.SetInt((int)draw_dim.Z(), "Draw_resolution", 0, 1000); ed.SetToolTipText("Number of quads to approximate round geometries, suggested = 16"); edc.Add(ed);
	

	Matrix kmat(3,3),dmat(3,3);
	for(int i=1; i<=kmat.Getrows(); i++)
	{
		for(int j=1; j<=kmat.Getcols(); j++)
		{
			kmat(i,j) = k(i,j);
			dmat(i,j) = d(i,j);
		}
	}
	ed.SetMatrix(kmat.GetMatPtr(), kmat.Getrows(), kmat.Getcols(), "spring_stiffness"); ed.SetToolTipText("Spring stiffness w.r.t. inertial coordinate system"); edc.Add(ed);
	ed.SetMatrix(dmat.GetMatPtr(), dmat.Getrows(), dmat.Getcols(), "damping_coefficient"); ed.SetToolTipText("Damping coefficient w.r.t. inertial coordinate system"); edc.Add(ed);
	


	if (elements.Length()==1)
	{
		SetElemDataVector3D(edc, loccoords(1), "Local_joint_pos_body"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector3D(edc, p_global, "Global_joint_pos"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
	} else
	{
		SetElemDataVector3D(edc, loccoords(1), "Local_joint_pos_body1"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector3D(edc, loccoords(2), "Local_joint_pos_body2"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
	}
}

int LinearRotationalSpringDamper::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Constraint::SetElementData(edc);

	GetElemDataDouble(mbs, edc, "Draw_size", draw_dim.X(), 0);
	int dd;
	if (GetElemDataInt(mbs, edc, "Draw_resolution", dd, 0)) draw_dim.Z() = dd;
	

	Matrix kmat(3,3),dmat(3,3);
	GetElemDataMatrix(GetMBS(), edc, "spring_stiffness", kmat, 1);
	GetElemDataMatrix(GetMBS(), edc, "damping_coefficient", kmat, 1);
	for(int i=1; i<=3; i++)
	{
		for(int j=1; j<=3; j++)
		{
			k(i,j) = kmat(i,j);
			d(i,j) = dmat(i,j);
		}
	}



	if (elements.Length()==1)
	{
		GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body", loccoords(1));
		GetElemDataVector3D(GetMBS(), edc, "Global_joint_pos", p_global);
	} else
	{
		GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body1", loccoords(1));
		GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body2", loccoords(2));
	}

	return rv;
}
                     
//this results in a torque
Vector3D LinearRotationalSpringDamper::ComputeForce(double t)  const
{
	Matrix3D Ap(0.);
	Matrix3D A(0.);
	Matrix3D A1(1.);
	Matrix3D A2(1.);

	//assume: small rotations!

	if (IsGroundJoint())
	{
		A = GetBody3D(1).GetRotMatrix(loccoords(1));

		if (d.AbsNorm() != 0) 
			Ap = A.GetTp()*GetBody3D(1).GetRotMatrixP(loccoords(1));
	}
	else
	{
		A1 = GetBody3D(1).GetRotMatrix(loccoords(1));
		A2 = GetBody3D(2).GetRotMatrix(loccoords(2));
		A = A1 * A2.GetTp();

		if (d.AbsNorm() != 0) 
		{
			Matrix3D Ap1 = GetBody3D(1).GetRotMatrixP(loccoords(1));
			Matrix3D Ap2 = GetBody3D(2).GetRotMatrixP(loccoords(2));
			Ap = Ap1*A1.GetTp() - Ap2*A2.GetTp();
		}
	}

	//compute local linearized angle and angular velocity
	Vector3D phi  (-A(2,3),A(1,3),-A(1,2));
	Vector3D omega(-Ap(2,3),Ap(1,3),-Ap(1,2));

	//compute torque in local coordinates and transform back into global coordinate system
	// k, d ... stiffness and damping matrix
	//mbs->UO() << "k:" << k << "\n";
	return A2*k*phi + A2*d*A2.GetTp()*omega;
}

double LinearRotationalSpringDamper::GetActorForce(double computation_time, int dir) const
{
	if (dir >=1 && dir <= 3) return (ComputeForce(computation_time))(dir);
	else return (ComputeForce(computation_time)).Norm();
}


void LinearRotationalSpringDamper::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	//-C_q^T \lambda = [dC/dx; dC/dy; dC/dphi]^T [\lambda1;\lambda2] 
	//e.g.: C = p_ref+R*p_loc
	double sign = 1;
	if (locelemind == 2) sign = -1;

	Vector3D M3D(ComputeForce(t)); //compute torque in global coordinates


	if (GetBody3D(locelemind).IsRigid())
	{
		if (f.Length() > max_rigid3D_coordinates)
		{
			assert(0 && "LinearRotationalSpringDamper::AddElementCqTLambda");
		}
		ConstVector<max_rigid3D_coordinates>  fadd(f.Length());
		fadd.SetAll(0);

		GetBody3D(locelemind).ApplyDrotrefdq(fadd, M3D); 
		fadd *= sign;

		f -= fadd;
	}
	else
	{
		GetBody3D(locelemind).GetdRotdqT(loccoords(locelemind),dpdq);
		//UO() << "dpdq=" << dpdq;

		for (int i=1; i <= f.Length(); i++)
		{
			f(i) -= sign*(dpdq(i,1)*M3D.X()+dpdq(i,2)*M3D.Y()+dpdq(i,3)*M3D.Z());
		}
	}
};

Vector3D LinearRotationalSpringDamper::GetRefPosD()	const 
{
	return GetBody3D(1).GetPosD(loccoords(1));
}

void LinearRotationalSpringDamper::DrawElement() 
{
	Constraint::DrawElement();

	int res = (int)draw_dim.Z();
	if (res < 2) res = 3;

	if (draw_dim.X() != 0)
	{
		mbs->SetColor(GetCol());
		Vector3D p1 = GetBody3D(1).GetPosD(loccoords(1));
		mbs->DrawSphere(p1,0.5*draw_dim.X(),res);

		mbs->SetColor(Vector3D(0.8,0.1,0.1));
		Vector3D p2;
		if (elements.Length()==1)
		{
			p2 = p_global;
		}
		else
		{
			p2 = GetBody3D(2).GetPosD(loccoords(2));
		}
		mbs->DrawSphere(p2,0.5*draw_dim.X(),res-1);
		if (draw_dim.Y() > 0) mbs->DrawZyl(p1,p2, draw_dim.Y(), res);
	}
};

Vector3D SpringDamperBearing::ComputeForce(double t)  const
{
	Vector3D p1;
	Vector3D p2;
	Vector3D v1;
	Vector3D v2(0.);

	Vector3D vp(0.,0.,0.);
	Vector3D v(0.,0.,0.);

	if (IsGroundJoint())
	{
		p1 = GetBody3D(1).GetPos(loccoords(1));
		p2 = p_global;
		v = p1 - p2;

		v1 = GetBody3D(1).GetVel(loccoords(1));
		
		if (d.AbsNorm() != 0) 
			vp = v1;
	}
	else
	{
		p1 = GetBody3D(1).GetPos(loccoords(1));
		p2 = GetBody3D(2).GetPos(loccoords(2));
		v = p1 - p2;

		if (d.AbsNorm() != 0) 
		{
			v1 = GetBody3D(1).GetVel(loccoords(1));
			v2 = GetBody3D(2).GetVel(loccoords(2));
			vp = v1 - v2;
		}
	}

	Vector3D f(0.);
	if (k_nonlinear)
	{
		//this model only uses diagonal components!
		double kv11 = 0;
		double kv22 = 0;
		double kv33 = 0;

		if (rotcomp == 3)
		{
			kv11 = mf_kxx.Evaluate(v.X());
			kv22 = mf_kyy.Evaluate(v.Y());
			kv33 = k(3,3) * v.Z();
		}
		else
		{
			GetMBS()->UO() << "Error: SpringDamperBearing:rotcomp must be 3 for nonlinear bearing!\n";
		}
		f += Vector3D(kv11,kv22,kv33) + d*vp;
	}
	else
	{
		f += k*v + d*vp;
	}
	if (use_bearing_fault)
	{
		int thisrotcomp = 3;
		if (rotcomp != 3) 
			GetMBS()->UO() << "Error: SpringDamperBearing: rotcomp has been set to 3!!!\n";

		double phi_inner = 0;
		double phi_outer = 0;
		double phi_cage;

		if (!GetBody3D(1).IsRigid() || GetBody3D(1).SOS() != 6) 
		{
			GetMBS()->UO() << "Error: SpringDamperBearing:no rigid/kardan body1!\n";
		}
		else 
		{
			phi_inner = GetBody3D(1).XG(3+thisrotcomp);
		}
		if (!IsGroundJoint())
		{
			//first body must be on inner ring
			//second body must be on outer ring
			if (!GetBody3D(2).IsRigid() || GetBody3D(2).SOS() != 6) 
			{
				GetMBS()->UO() << "Error: SpringDamperBearing:no rigid/kardan body2!\n";
			}
			else 
			{
				phi_outer = GetBody3D(2).XG(3+thisrotcomp);
			}
		}

		Vector3D dir;
		double f_fault;
		ComputeBearingGeometry(p1, p2, v1, v2, phi_inner, phi_outer, phi_cage, dir, f_fault);

		//Vector3D dir(cos(bearing_fault_dir*phidir),sin(bearing_fault_dir*phidir),0.);
		//double f_fault = bearing_fault_amp * pow((sin(delta_phi*bearing_fault_freqfact)+1.)*0.5,bearing_fault_power);

		f += f_fault * dir;
	}
	if (random_noise_factor != 0.)
	{
		f += random_noise_factor * random_value;
	}

	// k, d ... stiffness and damping matrix
	return f;
}

void SpringDamperBearing::StartTimeStep()
{
	random_value = (double)rand()/(double)RAND_MAX;
}

void SpringDamperBearing::ComputeBearingGeometry(const Vector3D& p1, const Vector3D& p2, const Vector3D& v1, const Vector3D& v2, 
		const double& phi_inner, const double& phi_outer, double& phi_cage, Vector3D& bearing_fault_vec, double& f_bearing_fault) const
{
	int thisrotcomp = rotcomp;
	if (rotcomp != 3)
	{
		thisrotcomp = 3;
		GetMBS()->UO() << "Error: SpringDamperBearing: rotcomp has been set to 3!!!\n";
	}

	double phidir = phi_outer;

	if (bearing_fault_innerring)
	{
		phidir = phi_inner;
	}

	//distange cage has traveled (proportional to v_cage)
	double x_cage = 0.5*(inner_roll_radius * phi_inner + outer_roll_radius * phi_outer);
	phi_cage = x_cage / (0.5*(inner_roll_radius + outer_roll_radius ));
	//wrong: phi_cage = 0.5*(phi_inner + phi_outer);

	bearing_fault_vec = Vector3D(1.,0.,0.);
	f_bearing_fault = 0;

	if (use_bearing_fault)
	{

		double phidir = phi_outer;
		double draw_rad = outer_roll_radius;


		double phi_fault;

		if (bearing_fault_innerring)
		{
			phidir = phi_inner;
			draw_rad = inner_roll_radius;
			phi_fault = phi_inner - phi_cage;
		}
		else
		{
			phi_fault = phi_outer - phi_cage;
		}

		bearing_fault_vec = Vector3D(cos(phidir),sin(phidir),0.);

		f_bearing_fault = bearing_fault_amp * pow((cos(phi_fault*n_bearing_balls)+1.)*0.5,bearing_fault_power);

	}

}



void SpringDamperBearing::DrawElement() 
{
	Constraint::DrawElement();

	if (!n_bearing_balls || !draw_bearing)
	{
		InertialLinearSpringDamper::DrawElement();
		return;
	}

	int res = (int)draw_dim.Z();
	if (res < 2) res = 3;

	mbs->SetColor(GetCol());
	Vector3D p1 = GetBody3D(1).GetPosD(loccoords(1));
	Vector3D p2;
	//mbs->DrawSphere(p1,0.5*draw_dim.X(),res);

	Vector3D v1(0.); //not needed for drawing
	Vector3D v2(0.); //not needed for drawing

	Matrix3D A1 = GetBody3D(1).GetRotMatrixD(Vector3D(0.,0.,0.));
	Matrix3D A2(1.);

	if (elements.Length()==1)
	{
		p2 = p_global;
	}
	else
	{
		p2 = GetBody3D(2).GetPosD(loccoords(2));
		A2 = GetBody3D(2).GetRotMatrixD(Vector3D(0.,0.,0.));
	}

	//draw inner ring of bearing:
	{
		TArray<Vector3D> cols;
		cols.Add(col);	cols.Add(col);	cols.Add(col);	cols.Add(col);	cols.Add(col);
		TArray<Vector2D> segs;
		segs.Add(Vector2D(-t_ring_draw*0.5 ,r_ring_inner_draw));	
		segs.Add(Vector2D( t_ring_draw*0.5 ,r_ring_inner_draw));	
		segs.Add(Vector2D( t_ring_draw*0.5 ,inner_roll_radius));	
		segs.Add(Vector2D(-t_ring_draw*0.5 ,inner_roll_radius));	
		segs.Add(Vector2D(-t_ring_draw*0.5 ,r_ring_inner_draw));	

		GeomRotObject3D zyl1(GetMBS(), p1, segs, cols, A1, 32);
		zyl1.DrawYourself();

		mbs->SetColor(Vector3D(0.8,0.2,0.2));
		mbs->DrawSphere(p1+A1*Vector3D(r_ring_inner_draw,0.,0.),t_ring_draw*0.2,res-1);
	}
	//draw outer ring of bearing
	{
		TArray<Vector3D> cols;
		cols.Add(col);	cols.Add(col);	cols.Add(col);	cols.Add(col);	cols.Add(col);
		TArray<Vector2D> segs;
		segs.Add(Vector2D(-t_ring_draw*0.5 ,outer_roll_radius));	
		segs.Add(Vector2D( t_ring_draw*0.5 ,outer_roll_radius));	
		segs.Add(Vector2D( t_ring_draw*0.5 ,r_ring_outer_draw));	
		segs.Add(Vector2D(-t_ring_draw*0.5 ,r_ring_outer_draw));	
		segs.Add(Vector2D(-t_ring_draw*0.5 ,outer_roll_radius));	

		GeomRotObject3D zyl1(GetMBS(), p2, segs, cols, A2, 32);
		zyl1.DrawYourself();

		mbs->SetColor(Vector3D(0.8,0.2,0.2));
		mbs->DrawSphere(p2+A2*Vector3D(r_ring_outer_draw,0.,0.),t_ring_draw*0.2,res-1);
	}

	if (rotcomp == 3) //only draw in this case
	{
		double phi_inner = 0;
		double phi_outer = 0;
		double phi_cage;

		if (!GetBody3D(1).IsRigid() || GetBody3D(1).SOS() != 6) 
		{
			GetMBS()->UO() << "Error: SpringDamperBearing:no rigid/kardan body1!\n";
		}
		else 
		{
			phi_inner = GetBody3D(1).XGD(3+rotcomp);
		}
		if (!IsGroundJoint())
		{
			//first body must be on inner ring
			//second body must be on outer ring
			if (!GetBody3D(2).IsRigid() || GetBody3D(2).SOS() != 6) 
			{
				GetMBS()->UO() << "Error: SpringDamperBearing:no rigid/kardan body2!\n";
			}
			else 
			{
				phi_outer = GetBody3D(2).XGD(3+rotcomp);
			}
		}

		Vector3D dir;
		double f_fault;
		//compute phi_cage, dir(ection) of fault and, f_fault:
		ComputeBearingGeometry(p1, p2, v1, v2, phi_inner, phi_outer, phi_cage, dir, f_fault);
		//if (random_noise_factor != 0.)
		//{
		//	f_fault += random_noise_factor * random_value;
		//}


		//draw balls, rotating with cage frequency:
		double r_ball = 0.5*(outer_roll_radius - inner_roll_radius);
		double r_ball_run = r_ball + inner_roll_radius;
		double nb = n_bearing_balls;

		double x_outer = outer_roll_radius * phi_outer; //velocity on outer ring
		double x_ball = (x_outer-phi_cage*outer_roll_radius); //velocity on ball surface, relative to ball
		double phi_ball = 2.*x_ball/r_ball + phi_cage;

		for (int i=1; i <= n_bearing_balls; i++)
		{
			double phi = (double)i/nb*2*MY_PI + phi_cage;
			Vector3D dir(cos(phi), sin(phi), 0.);
			Vector3D p = r_ball_run * dir + p1;

			mbs->SetColor(Vector3D(0.7,0.7,0.7));
			mbs->DrawSphere(p,r_ball,res-1);

			if (i == 1)
			{
				Vector3D r(cos(phi_ball),sin(phi_ball),0.);
				r *= r_ball;
				mbs->SetColor(Vector3D(0.8,0.2,0.2));
				mbs->DrawSphere(p + r,r_ball*0.15,res-1);
			}
		}

		if (use_bearing_fault)
		{
			double sign = 1;
			double draw_rad = outer_roll_radius;

			if (bearing_fault_innerring)
			{
				sign = -1;
				draw_rad = inner_roll_radius;
			}

			double fsize = f_fault * GetMBS()->GetDOption(120);
			double linethickness = fsize * 0.04;
			double headsize = fsize * 0.2;

			Vector3D pf = p1 + draw_rad * dir;
			mbs->MyDrawArrow(pf +  (-1.* sign * fsize) * dir, pf,
			Vector3D(mbs->GetDOption(121),mbs->GetDOption(122),mbs->GetDOption(123)), linethickness, headsize, 6);

		}
	}
};




