//#**************************************************************
//#
//# filename:             SpecialConstraints.cpp
//#
//# author:               Gerstmayr Johannes, Reischl Daniel
//#
//# generated:						20.April 2011
//# description:          In this file all implemented constraints, that are not lower kinematic pairs or rigid joints are collected.
//#												Lower kinematic pairs and rigid joints are in KinematicPairs.h
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
 

#include "element.h"
#include "constraint.h"
#include "kinematicpairs.h"
#include "specialconstraints.h"
#include "geomelements.h"
#include "elementdataaccess.h"
#include "graphicsconstants.h"
#include "rendercontext.h"
#include "material.h"
#include "finiteelement3d.h"
#include "..\UtilityLib\femathHelperFunctions.h"	//$ DR 2013-03-25 needed for Quicksort in Rope3D

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// CoordConstraint
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//$ MSax 2013-01-15: added
//void CoordConstraint::Initialize() 
//{
//	//if (MaxIndex()==3 || UsePenaltyFormulation()) // position constraint
//	//{
//	//	if (elements.Length()==2)
//	//	{
//	//		absolute_coord_offset = coord_gain_factor*(GetElem(1).GetXInit()(loccoords(1)))-GetElem(2).GetXInit()(loccoords(2))+coord_offset;
//	//	}
//	//	else 
//	//	{
//	//		absolute_coord_offset = coord_gain_factor*(GetElem(1).GetXInit()(loccoords(1)))+coord_offset;
//	//	}
//	//}
//	//else // velocity constraint
//	//{
//	//	if (elements.Length()==2)
//	//	{
//	//		absolute_coord_offset = coord_gain_factor*(GetElem(1).GetXInit()(loccoords(1)+GetElem(1).SOS()))-GetElem(2).GetXInit()(loccoords(2)+GetElem(2).SOS());
//	//	}
//	//	else 
//	//	{
//	//		absolute_coord_offset = coord_gain_factor*(GetElem(1).GetXInit()(loccoords(1)+GetElem(1).SOS()));
//	//	}
//	//}
//	// position constraint
//	if (elements.Length()==2)
//	{
//		absolute_coord_offset = coord_gain_factor*(GetElem(1).GetXInit()(loccoords(1)))-GetElem(2).GetXInit()(loccoords(2))+coord_offset;
//	}
//	else 
//	{
//		absolute_coord_offset = coord_gain_factor*(GetElem(1).GetXInit()(loccoords(1)))+coord_offset;
//	}
//
//	// velocity constraint
//	if (elements.Length()==2)
//	{
//		absolute_coord_offset_vel = coord_gain_factor*(GetElem(1).GetXInit()(loccoords(1)+GetElem(1).SOS()))-GetElem(2).GetXInit()(loccoords(2)+GetElem(2).SOS());
//	}
//	else 
//	{
//		absolute_coord_offset_vel = coord_gain_factor*(GetElem(1).GetXInit()(loccoords(1)+GetElem(1).SOS()));
//	}
//};

int CoordConstraint::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = Constraint::CheckConsistency(errorstr);
	if (rv) return rv;

	if(!UsePenaltyFormulation())
	{
		if (coord_offset!=0 && relative_to_inital_values) // position constraint
		{
			mbs->UO() << "ERROR: CoordConstraint::Invalid coord_offset value. For a lagrange constraint the coord_offset value must be set to 0 if flag relative_to_initial_values is set.\n";
			rv=1;
		}
	}

	// AP 26-04-2013: bugfix, check XInit only for elements which exist, not automatically for element 1 and 2
	for (int i=1; i<=elements.Length(); i++)
	{
		//if (loccoords(i) <= GetElem(i).SOS()) //do not check this, if constraint applied to velocity constraints!
		//{
		//	if (GetElem(i).GetXInit()(loccoords(i)+GetElem(i).SOS()) != 0)
		//	{
		//		mbs->UO() << "WARNING: CoordConstraint::Constraint coord. have initial velocity unequal 0, which might lead to erroneous results. Use VelocityCoordinateConstraint instead.\n";		
		//		rv=1;
		//	}
		//}

		//if (GetElem(i).GetXInit()(loccoords(i)+GetElem(i).SOS()) != 0)
		//{
		//	mbs->UO() << "WARNING: CoordConstraint::Constraint coord. have initial velocity unequal 0, which might lead to erroneous results. Use VelocityCoordinateConstraint instead.\n";		
		//	rv=1;
		//}
	}


	return rv;
}




void CoordConstraint::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Constraint::GetElementData(edc);

	ElementData ed;	
	//ed.SetDouble(coord_offset, "coord_offset"); ed.SetToolTipText("coordinate offset, see documentation section equation"); edc.Add(ed);
	//ed.SetDouble(coord_gain_factor, "coord_gain_factor"); ed.SetToolTipText("coordinate gain factor, see documentation section equation"); edc.Add(ed);

	ed.SetDouble(GetDrawSizeScalar(), "drawing_size"); ed.SetToolTipText("General drawing size of constraint"); edc.TreeAdd("Graphics",ed);
	//ed.SetInt((int)draw_dim.Z(), "Draw_resolution", 0, 1000); ed.SetToolTipText("Number of quads to approximate round geometries, suggested = 16"); edc.Add(ed);

	int en1 = 0;
	int en2 = 0;
	int lc1 = 0;
	int lc2 = 0;

	if (NE() >= 1) en1 = GetElnum(1); 
	if (NE() >= 2) en2 = GetElnum(2); 
	if (loccoords.Length() >= 1) lc1 = loccoords(1); 
	if (loccoords.Length() >= 2) lc2 = loccoords(2); 

	ed.SetInt(en1, "element_number"); ed.SetToolTipText("element number for coordinate 1"); edc.TreeAdd("Coordinate1",ed);
	ed.SetInt(lc1, "local_coordinate"); ed.SetToolTipText("Local coordinate of element 1 to be constrained"); edc.TreeAdd("Coordinate1",ed);

	ed.SetInt(en2, "element_number"); ed.SetToolTipText("element number for coordinate 2; for ground joint, set element number to zero"); edc.TreeAdd("Coordinate2",ed);
	ed.SetInt(lc2, "local_coordinate"); ed.SetToolTipText("Local coordinate of element 2 to be constrained"); edc.TreeAdd("Coordinate2",ed);

	//Vector elnum(NE());
	//for (int i=1; i <= NE(); i++) {elnum(i) = GetElnum(i);}
	//if (NE())
	//{
	//	ed.SetVector(elnum.GetVecPtr(), elnum.Length(), "Constraint_element_numbers"); ed.SetToolTipText("Only valid element numbers permitted!"); edc.Add(ed);
	//}

	//ed.SetInt(loccoords(1), "Element_Coordinate1"); ed.SetToolTipText("Local coordinate of element 1 to be constrained"); edc.Add(ed);
	//if (loccoords.Length() == 2)
	//{
	//	ed.SetInt(loccoords(2), "Element_Coordinate2"); ed.SetToolTipText("Local coordinate of element 2 to be constrained"); edc.Add(ed);
	//}

	//if(UsePenaltyFormulation())
	//{
	//	ed.SetDouble(GetPenaltyStiffness(), "Penalty_stiffness"); ed.SetToolTipText("Penalty stiffness of coordinate constraint"); edc.Add(ed);
	//}
}

int CoordConstraint::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Constraint::SetElementData(edc);
	double tmp;
	GetElemDataDouble(mbs, edc, "Graphics.drawing_size", tmp, 1);
	SetDrawSizeScalar(tmp);
	//int dd;
	//if (GetElemDataInt(mbs, edc, "Draw_resolution", dd, 0)) draw_dim.Z() = dd;


	//SetPenaltyFormulation(1); //done in constraint
	//SetPenaltyStiffness(spring_stiffnessi); //done in constraint

	int en1 = 0;
	int en2 = 0;
	int lc1 = 0;
	int lc2 = 0;

	loccoords.SetLen(0);
	elements.SetLen(0);

	GetElemDataInt(mbs, edc, "Coordinate1.local_coordinate", lc1, 1);
	GetElemDataInt(mbs, edc, "Coordinate1.element_number", en1, 1);

	GetElemDataInt(mbs, edc, "Coordinate2.local_coordinate", lc2, 1);
	GetElemDataInt(mbs, edc, "Coordinate2.element_number", en2, 1);

	AddElementCoord(en1, lc1);
	if (en2!=0 && lc2 != 0) {AddElementCoord(en2, lc2);}

	x_init = Vector(SS()); //this must be done after penalty formulation has been set

	//GetElemDataDouble(mbs, edc, "coord_offset", coord_offset, 1);
	//GetElemDataDouble(mbs, edc, "coord_gain_factor", coord_gain_factor, 1);

	//GetElemDataInt(mbs, edc, "Element_Coordinate1", loccoords(1), 1);
	//if (loccoords(1) < 1 || loccoords(1) > GetElem(1).SS())
	//{
	//	loccoords(1) = 1;
	//	GetMBS()->EDCError("Illegal element coordinate in CoordConstraint");
	//	rv = 0;
	//}

	//if (loccoords.Length() == 2)
	//{
	//	GetElemDataInt(mbs, edc, "Element_Coordinate2", loccoords(2), 1);

	//	if (loccoords(2) < 1 || loccoords(2) > GetElem(2).SS())
	//	{
	//		loccoords(2) = 1;
	//		GetMBS()->EDCError("Illegal element coordinate in CoordConstraint");
	//		rv = 0;
	//	}

	//}

	//if(UsePenaltyFormulation())
	//{
	//	GetElemDataDouble(mbs, edc, "Penalty_stiffness", tmp);
	//	SetPenaltyStiffness(tmp);
	//}

	return rv;
}

void CoordConstraint::EvalG(Vector& f, double t) 
{
	if(UsePenaltyFormulation()) return; //(RL)

	if (elements.Length()<1 || elements.Length()>2)
	{
		mbs->UO() << "ERROR: CoordConstraint::EvalG, number of elements != 1 or 2\n";
	}
	double init_val = 0;

	if (MaxIndex()==3)
	{
		if (elements.Length()>=1) 
		{
			if (relative_to_inital_values)
			{
				init_val = GetElem(1).GetXInit()(loccoords(1));
			}

			//f(1) = GetElem(1).XG(loccoords(1))-GetElem(1).GetXInit()(loccoords(1))-coord_offset;
			f(1) = coord_gain_factor*(GetElem(1).XG(loccoords(1))-init_val)-coord_offset;  // $ MSax 2013-04-22 ==> added gain factor
		}
		if (elements.Length()==2) 
		{
			if (relative_to_inital_values)
			{
				init_val = GetElem(2).GetXInit()(loccoords(2));
			}
			//f(1)-= GetElem(2).XG(loccoords(2))-GetElem(2).GetXInit()(loccoords(2))/*-coord_offset*/; //changed by JG 2013-01-11 ==> this cannot work, coord_offset would be useless
			f(1)-= GetElem(2).XG(loccoords(2))-init_val; //changed by JG 2013-01-11 ==> this cannot work, coord_offset would be useless
		}
	}
	else
	{ //velocity constraints:
		if (elements.Length()>=1) 
		{
			//f(1) = GetElem(1).XGP(loccoords(1))-GetElem(1).GetXInit()(loccoords(1)+GetElem(1).SOS());
			f(1) = coord_gain_factor*(GetElem(1).XGP(loccoords(1)));  // $ MSax 2013-04-22 ==> added gain factor, removed init velocity
		}
		if (elements.Length()==2) 
		{
			//f(1)-= GetElem(2).XGP(loccoords(2))-GetElem(2).GetXInit()(loccoords(2)+GetElem(2).SOS());
			f(1)-= GetElem(2).XGP(loccoords(2));
		}
	}
};

void CoordConstraint::EvalF2(Vector& f, double t)
{
	//f = [f1 f2], where f1 is the residual vector of constraint element1 and f2 of constraint element 2

	if(!UsePenaltyFormulation()) return; //(RL)

	double force = ComputeForce(t);
	for (int i=1; i <= NE(); i++)
	{
		double sign = coord_gain_factor; //sign = du/dq; //$ MSax 2013-01-14: ==> added gain factor
		int offset = 0;
		if (i==2) 
		{
			sign = -1.; //sign = du/dq;
			offset = GetElem(1).SOS();
		}
		int lc = loccoords(i);
		if(lc > GetElem(i).SOS() && lc <=GetElem(i).SOS()*2)	//$ JG+DR 2012-12-13 is constraint of a velocity coordinate
		{
			lc -= GetElem(i).SOS();
		}
		f(lc + offset) -= sign*force;
	}
}

void CoordConstraint::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	if (UsePenaltyFormulation()) return; //(RL)

	//-C_q^T \lambda = -1*\lambda at q_loccoord (coordinate which is constrained)
	//double sign = 1;
	//if (locelemind == 2) sign = -1;
	//f(loccoords(locelemind)) -= sign*ComputeForce(t);

	//$ MSax 2013-01-14: ==> added gain factor
	double signFactor;
	if (locelemind == 1) 
	{
		signFactor = coord_gain_factor;
	}
	else // locelemind = 2
	{
		signFactor = -1;
	}
	f(loccoords(locelemind)) -= signFactor*ComputeForce(t);

	//if (!UsePenaltyFormulation())
	//{
	//	f(loccoords(locelemind)) -= sign*XG(1);
	//}
	//else
	//{
	//	if (elements.Length()==1) 
	//	{
	//		double u1 = GetElem(1).XG(loccoords(1))-GetElem(1).GetXInit()(loccoords(1));
	//		f(loccoords(locelemind)) -= sign * spring_stiffness * u1;
	//	}
	//	else
	//	{
	//		double u1 = GetElem(1).XG(loccoords(1))-GetElem(1).GetXInit()(loccoords(1));
	//		double u2 = GetElem(2).XG(loccoords(2))-GetElem(2).GetXInit()(loccoords(2));
	//		f(loccoords(locelemind)) -= sign * spring_stiffness * (u1-u2);
	//	}
	//}
	//GetMBS()->UO() << "fc(" << loccoords(locelemind) << ")=" << sign*XG(1) << "\n";
};

double CoordConstraint::ComputeForce(double t) const
{
	if (UsePenaltyFormulation())
	{
		double f = 0;
		double u1 = GetElem(1).XG(loccoords(1));
		if (relative_to_inital_values)
		{
			u1 = u1-GetElem(1).GetXInit()(loccoords(1));
		}
		double u2 = 0;
		double u1p = GetElem(1).XGP(loccoords(1));
		double u2p = 0;
		if (elements.Length()==2) 
		{
			u2 = GetElem(2).XG(loccoords(2)); //***JG //double u2
			if (relative_to_inital_values)
			{
				u2 = u2-GetElem(2).GetXInit()(loccoords(2));
			}
			u2p = GetElem(2).XGP(loccoords(2));
		}
		f = GetPenaltyStiffness() * (coord_gain_factor*u1-u2-coord_offset);
		f += GetPenaltyDamping()*(coord_gain_factor*u1p-u2p); // $ MSax 2013-04-22
		return f;
	}
	else
	{
		return XG(1); //Lagrange parameter
	}
}


void CoordConstraint::DrawElement() 
{
	Constraint::DrawElement();
	//

	if (GetDrawSizeScalar() == 0) return;

	mbs->SetColor(GetCol());

	Vector3D dir = GetElem(1).GetDOFDirD(loccoords(1));
	if (RelApproxi(dir.Norm(),1.))
	{
		Vector3D p = GetElem(1).GetDOFPosD(loccoords(1));
		mbs->DrawCone(p + (-GetDrawSizeScalar()*0.75)*dir,p,(GetDrawSizeScalar()*0.6),6,1);
	}
	else if (RelApproxi(dir.Norm(),2.))
	{
		dir *= 0.5;
		Vector3D p = GetElem(1).GetDOFPosD(loccoords(1));
		mbs->DrawCone(p + (-GetDrawSizeScalar()*1.)*dir,p,(GetDrawSizeScalar()*0.5),6,1);
		mbs->DrawCone(p + (-GetDrawSizeScalar()*1.9)*dir,p + (-GetDrawSizeScalar()*0.9)*dir,(GetDrawSizeScalar()*0.5),8,1);
	}
	else
	{
		mbs->DrawSphere(GetElem(1).GetDOFPosD(loccoords(1)),GetDrawSizeScalar()*0.5);
	}
	if (loccoords.Length() == 2)
	{
		Vector3D dir = GetElem(2).GetDOFDirD(loccoords(2));
		if (RelApproxi(dir.Norm(),1.))
		{
			Vector3D p = GetElem(2).GetDOFPosD(loccoords(2));
			mbs->DrawCone(p + (-GetDrawSizeScalar()*0.75)*dir,p,(GetDrawSizeScalar()*0.6),6,1);
		}
		else if (RelApproxi(dir.Norm(),2.))
		{
			dir *= 0.5;
			Vector3D p = GetElem(2).GetDOFPosD(loccoords(2));
			mbs->DrawCone(p + (-GetDrawSizeScalar()*1.)*dir,p,(GetDrawSizeScalar()*0.5),6,1);
			mbs->DrawCone(p + (-GetDrawSizeScalar()*1.9)*dir,p + (-GetDrawSizeScalar()*0.9)*dir,(GetDrawSizeScalar()*0.5),8,1);
		}
		else
		{
			mbs->DrawSphere(GetElem(2).GetDOFPosD(loccoords(2)),GetDrawSizeScalar()*0.5);
		}
	}

};

void CoordConstraint::GetGlobalConstrainedDOFs(TArray<int>& dofs) const
{
	int GlobalNr;
	const Element& elem = GetMBS()->GetElement(elements(1));
	GlobalNr = elem.LTG(loccoords(1));
	dofs.Add(GlobalNr);
}


void VelCoordConstraint::EvalG(Vector& f, double t) 
{
	if(UsePenaltyFormulation()) return;

	if (elements.Length()<1 || elements.Length()>2)
	{
		mbs->UO() << "ERROR: VelCoordConstraint::EvalG, number of elements != 1 or 2\n";
	}
	double init_val = 0;

	if (elements.Length()>=1) 
	{
		if (relative_to_inital_values)
		{
			init_val = GetElem(1).GetXInit()(loccoords(1)+GetElem(1).SOS());
		}
		f(1) = coord_gain_factor*(GetElem(1).XGP(loccoords(1))-init_val)-coord_offset;
	}
	if (elements.Length()==2) 
	{
		if (relative_to_inital_values)
		{
			init_val = GetElem(2).GetXInit()(loccoords(2)+GetElem(2).SOS());
		}
		f(1)-= GetElem(2).XGP(loccoords(2))-init_val;
	}
};

//void VelCoordConstraint::EvalF2(Vector& f, double t)
//{
//
//	//if(!UsePenaltyFormulation()) return;
//
//	//mbs->UO() << "ERROR: VelCoordConstraint::EvalF2, no penalty formulation implemented.\n";
//
//	if(!UsePenaltyFormulation()) return; //(RL)
//
//	double force = ComputeForce(t);
//	for (int i=1; i <= NE(); i++)
//	{
//		double sign = coord_gain_factor; //sign = du/dq; //$ MSax 2013-01-14: ==> added gain factor
//		int offset = 0;
//		if (i==2) 
//		{
//			sign = -1.; //sign = du/dq;
//			offset = GetElem(1).SOS();
//		}
//		int lc = loccoords(i);
//		//if(lc > GetElem(i).SOS() && lc <=GetElem(i).SOS()*2)	//$ JG+DR 2012-12-13 is constraint of a velocity coordinate
//		//{
//		//	lc -= GetElem(i).SOS();
//		//}
//		f(lc + offset) -= sign*force;
//	}
//}

double VelCoordConstraint::ComputeForce(double t) const // $ MSax 2013-07-01 : added penalty formulation
{
	if (UsePenaltyFormulation())
	{
		double f = 0;
		double v1 = GetElem(1).XGP(loccoords(1));
		if (relative_to_inital_values)
		{	
			v1 = v1-GetElem(1).GetXInit()(loccoords(1)+GetElem(1).SOS());
		}
		double v2 = 0;
		//double u1p = GetElem(1).XGP(loccoords(1));
		//double u2p = 0;
		if (elements.Length()==2) 
		{
			v2 = GetElem(2).XGP(loccoords(2)); //***JG //double u2
			if (relative_to_inital_values)
			{	
				v2 = v2-GetElem(2).GetXInit()(loccoords(2)+GetElem(2).SOS());
			}
			//u2p = GetElem(2).XGP(loccoords(2));
		}
		f = GetPenaltyStiffness() * (coord_gain_factor*v1-v2-coord_offset);
		//if(1)
		//{
		//	f += GetPenaltyDamping()*(u1p-u2p); // $ MSax 2013-04-22
		//}
		return f;
	}
	else
	{
		return XG(1); //Lagrange parameter
	}
}

void VelCoordConstraint::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	if (UsePenaltyFormulation()) return;

	double signFactor;
	if (locelemind == 1) 
	{
		signFactor = coord_gain_factor;
	}
	else // locelemind = 2
	{
		signFactor = -1;
	}
	f(loccoords(locelemind)) -= signFactor*ComputeForce(t);
};

int VelCoordConstraint::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = Constraint::CheckConsistency(errorstr);
	if (rv) return rv;

	if(!UsePenaltyFormulation())
	{
		if (coord_offset!=0 && relative_to_inital_values) // wrong initial configuration
		{
			mbs->UO() << "ERROR: VelCoordConstraint::Invalid coord_offset value. For a lagrange constraint the coord_offset value must be set to 0 if flag relative_to_initial_values is set.\n";
			rv=1;
		}
	}
	//else
	//{
	//		mbs->UO() << "ERROR: VelCoordConstraint::Penalty formulation not implemented.\n";		
	//		rv=1;
	//}

	return rv;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ArcLengthConstraint
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//void ArcLengthConstraint::SetArcLengthConstraint(double load, int elem, int loccoord)
//{
//	SetCoordConstraint(elem, loccoord, -1);
//	this->targetload = load;
//}
void ArcLengthConstraint::EvalG(Vector& f, double t)
{
	double L, psi;
	L = 1000.; 
	psi = 0.5;

	int dof = mbs->GetSecondOrderSize();

	double g = 0.;

	double du = this->GetElem(1).XG(this->loccoords(1)) - mbs->GetLastSolVector().Get(this->GetElem(1).LTG(this->loccoords(1)));
	double dlambda = XG(1) - mbs->GetLastSolVector().Get(LTG(1));

	g = du*du + psi*psi * dlambda*dlambda*targetload*targetload - L*L;	// equ. 5 in Neto, Feng: On the determinaiton of the path direction for arc-length methods

	f(1) += g;
}

void ArcLengthConstraint::AddElementCqTLambda(int locelemind, Vector& f)
{
	f(loccoords(locelemind)) -= 2*ComputeForce(0);
	//GetMBS()->UO() << "Lagrange_Multiplier_at_loccoord_" << loccoords(locelemind) << " = " << XG(1) << "\n";
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// PrescribedCoordConstraint
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void PrescribedCoordConstraint::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	CoordConstraint::GetElementData(edc);

	ElementData ed;

	//SetElemDataMathFunc(edc, mathfunc, "Function");
	ElementDataContainer edc_mf;
	mathfunc.GetElementData(edc_mf);
	ed.SetEDC(&edc_mf, "Prescribed_function"); ed.SetToolTipText("mathfunction which describes the prescribed coordinate"); edc.TreeAdd("Kinematics",ed);

	//position/velocity:
	
	//ed.SetBoolGroup(1-constraintmode, 2, "Position_level"); ed.SetToolTipText("Constraint is applied directly at position level"); edc.Add(ed);
	//ed.SetBoolGroup(constraintmode, 2, "Velocity_level"); ed.SetToolTipText("Constraint is applied at velocity level (reduced index)"); edc.Add(ed);
}

int PrescribedCoordConstraint::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = CoordConstraint::SetElementData(edc);

	//rv = rv*GetElemDataMathFunc(mbs, edc, "Function", mathfunc, 1);
	//const ElementData* ed = edc.TreeFind("Kinematics.Prescribed_function");
	ElementData* ed = edc.TreeFind("Kinematics.Prescribed_function");
	if (ed != 0 && ed->IsEDC())
	{
		//mathfunc.SetElementData(GetMBS(), *ed->GetEDC());
		mathfunc.SetElementData(GetMBS(), *ed->GetEDC()); //$ DR 2012-12-12
	}
	else
	{
		GetMBS()->EDCError("PrescribedCoordConstraint::SetElementData: Kinematics.Prescribed_function not found, or is not an EDC!\n");
		rv = 0;
	}
	/*
	//+++++++++++++++++++++++++++
	//mathfunc
	int mode;
	Matrix data;

	GetElemDataInt(GetMBS(), edc, "Function_type", mode, 1); 
	GetElemDataMatrix(GetMBS(), edc, "Function_data", data, 1);

	rv = rv*mathfunc.SetData(mode, data); //set data only accepts valid data!
	//+++++++++++++++++++++++++++
	*/

	//constraintmode = 0; //position level
	//int flag = 0;
	//GetElemDataBool(mbs, edc, "Velocity_level", flag, 0); 
	//if (flag) constraintmode = 1; //velocity level

	return rv;
}

void PrescribedCoordConstraint::EvalG(Vector& f, double t) 
{
	if(IsActive())
	{
		if (constraintmode == 0 && MaxIndex()==3)
		{ 
			//position constraint:
			if (elements.Length()==1)
			{
				f(1) = mathfunc.Evaluate(t) - GetElem(1).XG(loccoords(1));
			}
			else
			{
				f(1) = mathfunc.Evaluate(t) - (GetElem(1).XG(loccoords(1))-GetElem(2).XG(loccoords(2))); //x1-x2 = prescribed value
			}

		}
		// AP: not necessary to have MaxIndex <= 2 for velocity constraint
		else if (constraintmode == 1 /*&& MaxIndex()<=2*/)
		{
			//velocity constraint:
			if (elements.Length()==1)
			{
				f(1) = mathfunc.Evaluate(t) - GetElem(1).XGP(loccoords(1)); //first derivative w.r.t. time
			}
			else
			{
				f(1) = mathfunc.Evaluate(t) - (GetElem(1).XGP(loccoords(1))-GetElem(2).XGP(loccoords(2))); //x1-x2 = prescribed value
			}
		}
		else if (constraintmode == 0 && MaxIndex()<=2)
		{
			//position constraint differentiated:
			if (elements.Length()==1)
			{
				f(1) = mathfunc.Evaluate(t, 1) - GetElem(1).XGP(loccoords(1)); //first derivative w.r.t. time
			}
			else
			{
				f(1) = mathfunc.Evaluate(t, 1) - (GetElem(1).XGP(loccoords(1))-GetElem(2).XGP(loccoords(2))); //x1-x2 = prescribed value
			}
		}
	}
	else
	{
		f(1) = XG(1);
	}
};

//To be replaced in derived class
void PrescribedCoordConstraint::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	//-C_q^T \lambda = -1*\lambda at q_loccoord (coordinate which is constrained)
	double sign = 1;
	if (locelemind == 2) sign = -1;

	f(loccoords(locelemind)) -= sign*XG(1);
};


////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//// ArcLengthConstraint
////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//void ArcLengthConstraint::SetArcLengthConstraint(double load, double L, int elem, int loccoord)
//{
//	SetPenaltyFormulation(0); //JG
//	loccoords.SetLen(0);
//	elements.SetLen(0);
//	x_init = Vector(SS());
//	AddElementCoord(elem, loccoord);
//	this->load = load;
//	this->L = L,
//	elementname = GetElementSpec();
//
//	Vector v(1);
//	v(1) = 1;
//	
//	x_init = v;
//}
//
//void ArcLengthConstraint::EvalG(Vector& f, double t)
//{
//	double psi = 0.5;
//	double u = GetElem(1).XG(loccoords(1));
//	double u_old = mbs->GetLastSolVector().Get(GetElem(1).LTG(loccoords(1)));
//	double du = u - u_old;
//	double dlambda = XG(1) - mbs->GetLastSolVector().Get(LTG(1));
//	//f(1) += du*du + psi*psi * dlambda*dlambda * load*load - L*L;	// equ. 5 in Neto, Feng: On the determinaiton of the path direction for arc-length methods
//	f(1) += /*du*du +*/ dlambda - 1;
//}
//
//void ArcLengthConstraint::AddElementCqTLambda(double t, int locelemind, Vector& f)
//{
//	double u = GetElem(1).XG(loccoords(1));
//	double u_old = mbs->GetLastSolVector().Get(GetElem(1).LTG(loccoords(1)));
//	double du = u - u_old;
//
//	f(loccoords(locelemind)) += /*2*du*/XG(1);
//	//f(loccoords(locelemind)) += 2*(GetElem(1).XG(loccoords(1)) - mbs->GetLastSolVector().Get(GetElem(1).LTG(loccoords(1)))) * XG(1);
//	//GetMBS()->UO() << "Lagrange_Multiplier_at_loccoord_" << loccoords(locelemind) << " = " << XG(1) << "\n";
//}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Pos2DConstraint
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Pos2DConstraint::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Constraint::GetElementData(edc);

	ElementData ed;
	ed.SetDouble(draw_dim.X(), "Draw_size"); ed.SetToolTipText("General drawing size of constraint"); edc.Add(ed);

	ed.SetVector2D(loccoords(1).X(),loccoords(1).Y(), "Element1_Local_Pos"); ed.SetToolTipText("Local position of element 1 to be constrained"); edc.Add(ed);
	if (loccoords.Length() == 2)
	{
		ed.SetVector2D(loccoords(2).X(),loccoords(2).Y(), "Element2_Local_Pos"); ed.SetToolTipText("Local position of element 2 to be constrained"); edc.Add(ed);
	}
	else
	{
		ed.SetVector2D(p_global.X(), p_global.Y(), "Global_Pos"); ed.SetToolTipText("Global position to be constrained"); edc.Add(ed);
	}
}

int Pos2DConstraint::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Constraint::SetElementData(edc);

	GetElemDataDouble(mbs, edc, "Draw_size", draw_dim.X(), 0);

	GetElemDataVector2D(mbs, edc, "Element1_Local_Pos", loccoords(1), 1);

	if (loccoords.Length() == 2)
	{
		GetElemDataVector2D(mbs, edc, "Element2_Local_Pos", loccoords(2), 1);
	}
	else
	{
		GetElemDataVector2D(mbs, edc, "Global_Pos", p_global, 1);
	}

	return rv;
}



void Pos2DConstraint::EvalG(Vector& f, double t) 
{
	if (elements.Length()<1 || elements.Length()>2)
	{
		mbs->UO() << "ERROR: Pos2DConstraint::EvalG, number of elements != 1 or 2\n"; return;
	}

	if (MaxIndex()==3)
	{
		if (elements.Length()==1)
		{
			if (!IsPrescribed())
			{
				Vector2D v;
				if (!DOFmode)
				{
					v = GetBody2D(1).GetPos2D(loccoords(1));
				}
				else
				{
					v = GetBody2D(1).GetNodePos2D((int)loccoords(1).X());
				}

				v -= p_global;
				f(1) = v(1);
				f(2) = v(2);
			} 
			else if (IsPrescribed() == 1)
			{
				//prescribed curve
				double phi = p4;
				if (t < p3) phi += 0.5*Sqr(t)/p3*p2;
				else phi += 0.5*p3*p2+(t-p3)*p2;

				Vector2D p = p_global+Vector2D(p1*cos(phi),p1*sin(phi));

				Vector2D v;
				if (!DOFmode)
				{
					v = GetBody2D(1).GetPos2D(loccoords(1));
				}
				else
				{
					v = GetBody2D(1).GetNodePos2D((int)loccoords(1).X());
				}
				v -= p;
				f(1) = v(1);
				f(2) = v(2);
				//UO() << "pos2dc v=" << v << "\n";
			}
			else if (IsPrescribed() == 2)
			{
				int i = 1;

				while (i <= prescribedtime.Length() && t > prescribedtime(i)) {i++;}

				if (i > prescribedtime.Length()) {i = prescribedtime.Length();}

				if (i != 0)
				{
					Vector2D v;
					double t1 = 0;
					Vector2D vel1(0,0);
					if (i >= 2) 
					{
						t1 = prescribedtime(i-1);
						vel1 = prescribedvel(i-1);
					}
					double t2 = prescribedtime(i);
					double dt = t2-t1;
					if (dt <= 0) dt = 1;

					Vector2D vel = (1.-(t-t1)/dt)*vel1 + (1.-(t-t2)/dt)*prescribedvel(i);

					if (!DOFmode)
					{
						v = GetBody2D(1).GetVel2D(loccoords(1)) - vel;
					}
					else
					{
						v = GetBody2D(1).GetNodeVel2D((int)loccoords(1).X()) - vel;
					}
					f(1) = v(1);
					f(2) = v(2);
				}
			}
		}
		else
			if (elements.Length()==2) 
			{
				Vector2D v;
				if (DOFmode == 0)
				{
					v = GetBody2D(1).GetPos2D(loccoords(1))-GetBody2D(2).GetPos2D(loccoords(2));
				}
				else if (DOFmode == 1)
				{
					v = GetBody2D(1).GetNodePos2D((int)loccoords(1).X()) - 
						GetBody2D(2).GetPos2D(loccoords(2));
				}
				else if (DOFmode == 2)
				{
					v = GetBody2D(1).GetNodePos2D((int)loccoords(1).X()) - 
						GetBody2D(2).GetNodePos2D((int)loccoords(2).X());
				}

				f(1) = v(1);
				f(2) = v(2);
			}
	}
	else if (MaxIndex()<=2)
	{
		if (elements.Length()==1)
		{ 
			if (!IsPrescribed())
			{
				Vector2D v;
				if (!DOFmode)
				{
					v = GetBody2D(1).GetVel2D(loccoords(1));
				}
				else
				{
					v = GetBody2D(1).GetNodeVel2D((int)loccoords(1).X());
				}
				f(1) = v(1);
				f(2) = v(2);
			}
			else if (IsPrescribed() == 1)
			{
				//prescribed curve
				double phi = p4;
				if (t < p3) phi += 0.5*Sqr(t)/p3*p2;
				else phi += 0.5*p3*p2+(t-p3)*p2;

				double phip;
				if (t < p3) phip = t/p3*p2;
				else phip = p2;

				Vector2D vpre = Vector2D(-p1*phip*sin(phi),p1*phip*cos(phi));

				Vector2D v;
				if (!DOFmode)
				{
					v = GetBody2D(1).GetVel2D(loccoords(1))-vpre;
				}
				else
				{
					v = GetBody2D(1).GetNodeVel2D((int)loccoords(1).X())-vpre;
				}
				f(1) = v(1);
				f(2) = v(2);
			}
			else if (IsPrescribed() == 2)
			{
				//prescribed velocity list
				int i = 1;
				while (i <= prescribedtime.Length() && t > prescribedtime(i)) {i++;}

				if (i > prescribedtime.Length()) {i = prescribedtime.Length();}

				if (i != 0)
				{
					Vector2D v;
					double t1 = 0;
					Vector2D vel1(0,0);
					if (i >= 2) 
					{
						t1 = prescribedtime(i-1);
						vel1 = prescribedvel(i-1);
					}
					double t2 = prescribedtime(i);
					double dt = t2-t1;
					if (dt <= 0) dt = 1;

					//linear interpolation between times and velocities
					Vector2D vel = (1.-(t-t1)/dt)*vel1 + (1.-(t2-t)/dt)*prescribedvel(i);
					//if (t>=0.25) UO() << "prevel" << i << "=" << vel << ", v3=" << prescribedvel(3) << ", t3=" << prescribedtime(3) << ", len=" << prescribedtime.Length() << "\n";

					if (!DOFmode)
					{
						v = GetBody2D(1).GetVel2D(loccoords(1)) - vel;
					}
					else
					{
						v = GetBody2D(1).GetNodeVel2D((int)loccoords(1).X()) - vel;
					}
					f(1) = v(1);
					f(2) = v(2);
				}
			}
		}
		else
			if (elements.Length()==2) 
			{
				Vector2D v;
				if (DOFmode == 0)
				{
					v = GetBody2D(1).GetVel2D(loccoords(1))-GetBody2D(2).GetVel2D(loccoords(2));
				}
				else if (DOFmode == 1) //only first node has node number
				{
					v = GetBody2D(1).GetNodeVel2D((int)loccoords(1).X()) - 
						GetBody2D(2).GetVel2D(loccoords(2));
				}
				else if (DOFmode == 2)
				{
					v = GetBody2D(1).GetNodeVel2D((int)loccoords(1).X()) - 
						GetBody2D(2).GetNodeVel2D((int)loccoords(2).X());
				}

				f(1) = v(1);
				f(2) = v(2);
			}
	}

};

void Pos2DConstraint::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	//-C_q^T \lambda = [dC/dx; dC/dy; dC/dphi]^T [\lambda1;\lambda2] 
	//C = p_ref+R*p_loc
	double sign = 1;
	if (locelemind == 2) sign = -1;

	if (DOFmode == 0 || (DOFmode == 1 && locelemind == 2))
	{
		GetBody2D(locelemind).GetdPosdqT(loccoords(locelemind),dpdq);
		//UO() << "pos2dc dpdq1=" << dpdq << "\n";
		for (int i=1; i <= f.Length(); i++)
		{
			f(i) -= sign*(dpdq(i,1)*XG(1)+dpdq(i,2)*XG(2));
		}
	}
	else
	{
		Vector2D lambda(-sign*XG(1), -sign*XG(2)); //negative sign: needs to be subtracted from f!!!!
		GetBody2D(locelemind).AddNodedPosdqTLambda((int)loccoords(locelemind).X(),lambda,f);

		/*
		//large matrices dpdq:
		GetBody2D(locelemind).GetNodedPosdqT((int)loccoords(locelemind).X(),dpdq);
		for (int i=1; i <= f.Length(); i++)
		{
		f(i) -= sign*(dpdq(i,1)*XG(1)+dpdq(i,2)*XG(2));
		}
		*/	
	}

};

Vector3D Pos2DConstraint::GetRefPosD() const 
{
	Vector2D v;
	if (DOFmode == 0)
	{
		v = GetBody2D(1).GetPos2DD(loccoords(1));
	}
	else
	{
		v = GetBody2D(1).GetNodePos2DD((int)loccoords(1).X());
	}
	return Vector3D(v.X(),v.Y(),0);
}

void Pos2DConstraint::DrawElement() 
{
	Constraint::DrawElement();
	//

	int drawmode = 2; //1=conventional with spheres, 2=advanced like revolute joint

	if (drawmode == 1)
	{
		Vector2D v;
		if (DOFmode == 0)
		{
			v = GetBody2D(1).GetPos2DD(loccoords(1));
		}
		else
		{
			v = GetBody2D(1).GetNodePos2DD((int)loccoords(1).X());
			//GetMBS()->UO() << "Node" << (int)loccoords(1).X() << "=" << v << "->" << GetBody2D(1).ToP3D(v)<< "\n";
		}

		mbs->SetColor(GetCol());
		mbs->DrawSphere(GetBody2D(1).ToP3D(v),0.5*draw_dim.X(),6);
		if (IsPrescribed() && 0)
		{
			mbs->SetColor(Vector3D(0.3,0.3,0.1));
			mbs->DrawZyl(GetBody2D(1).ToP3D(v), GetBody2D(1).ToP3D(p_global), 0.5*draw_dim.X(), 8);
		}

		if (elements.Length()==2)
		{
			Vector2D v2;
			if (DOFmode != 2)
			{
				v2 = ((Body2D&)GetBody2D(2)).GetPos2DD(loccoords(2));
			}
			else
			{
				v2 = GetBody2D(2).GetNodePos2DD((int)loccoords(2).X());
			}
			mbs->SetColor(Vector3D(0.8,0.2,0.2));
			mbs->DrawSphere(GetBody2D(2).ToP3D(v2),0.5*draw_dim.X(),8);
		}
		else
		{
			mbs->SetColor(Vector3D(0.8,0.2,0.2));
			mbs->DrawSphere(GetBody2D(1).ToP3D(p_global),0.5*draw_dim.X(),8);
		}
	}
	else if (drawmode == 2) //draw as 3D joint
	{
		
		mbs->SetColor(GetCol());
		int res = 16;

		if (elements.Length()==1) 
		{

			Vector3D p1 = GetBody2D(1).GetPosD(Vector3D(loccoords(1).X(),loccoords(1).Y(),0.));
			Vector3D p2 = GetBody2D(1).ToP3D(p_global);
			Vector3D rot1(0.,0.,draw_dim.X()*0.5);
			Vector3D rot1r = GetBody2D(1).ToP3D(rot1);

			Vector3D pg = GetBody2D(1).ToP3D(p_global);

			mbs->DrawZyl(p1+rot1,p1+rot1*0.1,0.5*draw_dim.X(),res);
			mbs->DrawZyl(p1-rot1,p1-rot1*0.1,0.5*draw_dim.X(),res);
			mbs->SetColor(colgrey2);
			mbs->DrawZyl(pg+1.2*rot1,pg-1.2*rot1,0.4*draw_dim.X(),res);			
		} 
		else
		{
			Vector3D p1 = GetBody2D(1).GetPosD(Vector3D(loccoords(1).X(),loccoords(1).Y(),0.));
			Vector3D p2 = GetBody2D(2).GetPosD(Vector3D(loccoords(2).X(),loccoords(2).Y(),0.));

			Vector3D rot1(0.,0.,draw_dim.X()*0.5);
			Vector3D rot2(0.,0.,draw_dim.X()*0.5);

			Vector3D rot1r = GetBody2D(1).ToP3D(rot1);
			Vector3D rot2r = GetBody2D(2).ToP3D(rot2);

			mbs->DrawZyl(p1+rot1,p1+rot1*0.1,0.5*draw_dim.X(),res);
			mbs->DrawZyl(p1-rot1,p1-rot1*0.1,0.5*draw_dim.X(),res);
			mbs->SetColor(colgrey2);
			mbs->DrawZyl(p2+1.2*rot2,p2-1.2*rot2,0.4*draw_dim.X(),res);
		}
	}
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Angle2DConstraint
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Angle2DConstraint::EvalG(Vector& f, double t) 
{
	if (MaxIndex()==3)
	{
		if (elements.Length()==1)
		{
			f(1) = GetBody2D(1).GetAngle2D(loccoords(1)) - angle_global;
		}
		else if (elements.Length()==2)
		{
			f(1) = GetBody2D(1).GetAngle2D(loccoords(1)) - GetBody2D(2).GetAngle2D(loccoords(2));
		}
	}
	else if (MaxIndex()<=2)
	{
		if (elements.Length()==1)
		{
			f(1) = GetBody2D(1).GetAngle2DP(loccoords(1));
		}
		else if (elements.Length()==2) 
		{
			f(1) = GetBody2D(1).GetAngle2DP(loccoords(1)) - GetBody2D(2).GetAngle2DP(loccoords(2));
		}
	}
};

void Angle2DConstraint::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{

	//-C_q^T \lambda = [dC/dx; dC/dy; dC/dphi]^T [\lambda1;\lambda2] 
	//C = p_ref+R*p_loc
	double sign = 1;
	if (locelemind == 2) sign = -1;

	GetBody2D(locelemind).GetdAngle2DdqT(loccoords(locelemind),dpdq);

	for (int i=1; i <= f.Length(); i++)
	{
		f(i) -= sign*dpdq(i,1)*XG(1);
	}
};

Vector3D Angle2DConstraint::GetRefPosD() const 
{
	return GetBody2D(1).ToP3D(GetBody2D(1).GetPos2DD(loccoords(1)));
}

void Angle2DConstraint::DrawElement() 
{
	Constraint::DrawElement();
	//

	Vector2D v;
	v = GetBody2D(1).GetPos2DD(loccoords(1));

	mbs->SetColor(GetCol());
	mbs->DrawSphere(GetBody2D(1).ToP3D(v),0.5*draw_dim.X(),4);

	if (elements.Length()==2)
	{
		Vector2D v2;
		v2 = ((Body2D&)GetBody2D(2)).GetPos2DD(loccoords(2));

		mbs->SetColor(Vector3D(0.8,0.2,0.2));
		mbs->DrawSphere(GetBody2D(2).ToP3D(v2),0.5*draw_dim.X(),4);
	}
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Rotational Spring-Damper-Actor Element 2D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double SDRotActor2D::ComputeForce(double t) 
{
	double splitangle;
	double angspeed;
	if (elements.Length() == 1)
	{
		splitangle = GetBody2D(1).GetAngle2D(loccoords(1));
		angspeed = GetBody2D(1).GetAngle2DP(loccoords(1));
	}
	else
	{
		splitangle = GetBody2D(1).GetAngle2D(loccoords(1)) - GetBody2D(2).GetAngle2D(loccoords(2));
		angspeed = GetBody2D(1).GetAngle2DP(loccoords(1)) - GetBody2D(2).GetAngle2DP(loccoords(2));
	}
	return GetPenaltyStiffness()*splitangle + d*angspeed + ma;
}


void SDRotActor2D::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{

	//-C_q^T \lambda = [dC/dx; dC/dy; dC/dphi]^T [\lambda1;\lambda2] 
	//C = p_ref+R*p_loc
	double sign = 1;
	if (locelemind == 2) sign = -1;

	GetBody2D(locelemind).GetdAngle2DdqT(loccoords(locelemind),dpdq);
	double M = ComputeForce(t);
	
	for (int i=1; i <= f.Length(); i++)
	{
		f(i) -= sign*dpdq(i,1)*M;
	}
};

Vector3D SDRotActor2D::GetRefPosD() const 
{
	return GetBody2D(1).ToP3D(GetBody2D(1).GetPos2DD(loccoords(1)));
}

void SDRotActor2D::DrawElement() 
{
	Constraint::DrawElement();
	//

	Vector2D v;
	v = GetBody2D(1).GetPos2DD(loccoords(1));

	mbs->SetColor(GetCol());
	mbs->DrawSphere(GetBody2D(1).ToP3D(v),0.5*GetDrawSizeScalar(),4);

	if (elements.Length()==2)
	{
		Vector2D v2;
		v2 = ((Body2D&)GetBody2D(2)).GetPos2DD(loccoords(2));

		mbs->SetColor(Vector3D(0.8,0.2,0.2));
		mbs->DrawSphere(GetBody2D(2).ToP3D(v2),0.5*GetDrawSizeScalar(),4);
	}
};









//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// UniversalJointOLD
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void UniversalJointOLD::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Constraint::GetElementData(edc);

	ElementData ed;

	ed.SetDouble(draw_dim.X(), "Draw_size"); ed.SetToolTipText("General drawing size of constraint"); edc.Add(ed);
	ed.SetInt((int)draw_dim.Z(), "Draw_resolution", 0, 1000); ed.SetToolTipText("Number of quads to approximate round geometries, suggested = 16"); edc.Add(ed);
	ed.SetDouble(draw_dim.Y(), "Draw_sphere_radius"); ed.SetToolTipText("Draw dimension of sphere in joint"); edc.Add(ed);
	//ed.SetDouble(draw_dim.Y(), "Draw_axis_thickness"); ed.SetToolTipText("Draw dimension of axis"); edc.Add(ed);

	if (elements.Length()==1)
	{
		SetElemDataVector3D(edc, loccoords(1), "Local_joint_pos_body"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector3D(edc, p_global, "Global_joint_pos"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector3D(edc, loccoords(2), "Global_joint_axis1"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector3D(edc, loccoords(3), "Global_joint_axis2"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
	} else
	{
		SetElemDataVector3D(edc, loccoords(1), "Local_joint_pos_body1"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector3D(edc, loccoords(2), "Local_joint_pos_body2"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector3D(edc, loccoords(3), "Global_joint_axis1"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector3D(edc, loccoords(4), "Global_joint_axis2"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
	}
}

int UniversalJointOLD::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Constraint::SetElementData(edc);

	GetElemDataDouble(mbs, edc, "Draw_size", draw_dim.X(), 0);
	int dd;
	if (GetElemDataInt(mbs, edc, "Draw_resolution", dd, 0)) draw_dim.Z() = dd;

	GetElemDataDouble(GetMBS(), edc, "Draw_sphere_radius", draw_dim.Y());
	//GetElemDataDouble(GetMBS(), edc, "Draw_axis_thickness", draw_dim.Y());

	if (elements.Length()==1)
	{
		GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body", loccoords(1));
		GetElemDataVector3D(GetMBS(), edc, "Global_joint_pos", p_global);
		GetElemDataVector3D(GetMBS(), edc, "Global_joint_axis1", loccoords(2));
		GetElemDataVector3D(GetMBS(), edc, "Global_joint_axis2", loccoords(3));
	} else
	{
		GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body1", loccoords(1));
		GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body2", loccoords(2));
		GetElemDataVector3D(GetMBS(), edc, "Global_joint_axis1", loccoords(3));
		GetElemDataVector3D(GetMBS(), edc, "Global_joint_axis2", loccoords(4));
	}

	return rv;
}


void UniversalJointOLD::Initialize() 
{
	if (elements.Length()==1)
	{
		Vector3D lpos = loccoords(1);
		Vector3D vrot1 = loccoords(2);
		vrot1.Normalize();

		//if vrot1 and vrot2 not orthogonal --> orthogonalize vrot2!!!
		vrot1.GramSchmidt(loccoords(3));

		Matrix3D RT=GetBody3D(1).GetRotMatrix(lpos).GetTp();

		loccoords(4) = RT*vrot1; //local rot1 body1
		loccoords(5) = RT*loccoords(3); //local rot2 in body1, only for drawing
	} else
	{
		const Vector3D& lp1 = loccoords(1);
		const Vector3D& lp2 = loccoords(2);
		Vector3D vrot1 = loccoords(3);
		Vector3D vrot2 = loccoords(4);
		Matrix3D RT=GetBody3D(1).GetRotMatrix(lp1).GetTp();
		vrot1.Normalize();
		loccoords(5) = RT*vrot1;

		Matrix3D RT2=GetBody3D(2).GetRotMatrix(lp2).GetTp();
		vrot2.Normalize();
		loccoords(6) =  RT2*vrot2;

		//for drawing:
		loccoords(7) =  RT2*vrot1;
		loccoords(8) =  RT*vrot2;

	}
};

void UniversalJointOLD::EvalG(Vector& f, double t) 
{
	if (elements.Length()<1 || elements.Length()>2)
	{
		mbs->UO() << "ERROR: UniversalJointOLD::EvalG, number of elements != 1 or 2\n"; return;
	}

	if (MaxIndex()==3)
	{
		if (elements.Length()==1)
		{
			const Vector3D& lpos  = loccoords(1);
			const Vector3D& vrot2 = loccoords(3);
			const Vector3D& lrot1 = loccoords(4);

			Vector3D v = GetBody3D(1).GetPos(lpos)-p_global;
			f(1) = v(1); f(2) = v(2); f(3) = v(3);

			//vrot2^T*A(1)*lrot1=0
			Matrix3D A=GetBody3D(1).GetRotMatrix(lpos); 
			f(4) = vrot2*(A*lrot1);
		}
		else if (elements.Length()==2) 
		{
			const Vector3D& lp1   = loccoords(1);
			const Vector3D& lp2   = loccoords(2);
			const Vector3D& lrot1 = loccoords(5); //body1
			const Vector3D& lrot2 = loccoords(6); //body2

			Vector3D v = GetBody3D(1).GetPos(lp1)-GetBody3D(2).GetPos(lp2);
			f(1) = v(1); f(2) = v(2); f(3) = v(3);

			//lrot1^T*A(1)^T*A(2)*lrot2=0
			Matrix3D A1=GetBody3D(1).GetRotMatrix(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Matrix3D A2=GetBody3D(2).GetRotMatrix(lp2); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			v = A1*lrot1;
			f(4) = v*(A2*lrot2);
		}
	}
	else if (MaxIndex()<=2)
	{
		if (elements.Length()==1)
		{
			const Vector3D& lpos  = loccoords(1);
			const Vector3D& vrot2 = loccoords(3);
			const Vector3D& lrot1 = loccoords(4);

			Vector3D v = GetBody3D(1).GetVel(lpos);
			f(1) = v(1); f(2) = v(2); f(3) = v(3);

			//vrot2^T*Ap(1)*lrot1=0
			Matrix3D Ap=GetBody3D(1).GetRotMatrixP(lpos); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			f(4) = vrot2*(Ap*lrot1);
		}
		else if (elements.Length()==2) 
		{
			const Vector3D& lp1   = loccoords(1);
			const Vector3D& lp2   = loccoords(2);
			const Vector3D& lrot1 = loccoords(5); //body1
			const Vector3D& lrot2 = loccoords(6); //body2

			Vector3D v = GetBody3D(1).GetVel(lp1)-GetBody3D(2).GetVel(lp2);
			f(1) = v(1); f(2) = v(2); f(3) = v(3);

			//lrot1^T*Ap(1)^T*A(2)*lrot2 + lrot1^T*A(1)^T*Ap(2)*lrot2 = 0
			Matrix3D A1=GetBody3D(1).GetRotMatrix(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Matrix3D A2=GetBody3D(2).GetRotMatrix(lp2); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Matrix3D A1p=GetBody3D(1).GetRotMatrixP(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Matrix3D A2p=GetBody3D(2).GetRotMatrixP(lp2); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			v = A1*lrot1; 
			Vector3D vp = A1p*lrot1; 
			f(4) = v*(A2p*lrot2)+vp*(A2*lrot2);
		}
	}
};

void UniversalJointOLD::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{

	hmat.SetSize(f.Length(),IS());

	if (elements.Length() == 1) //--> ignore locelemind, only first element!
	{
		const Vector3D& lpos  = loccoords(1);
		const Vector3D& vrot2 = loccoords(3);
		const Vector3D& lrot1 = loccoords(4);

		GetBody3D(1).GetdPosdqT(lpos,dpdq);
		hmat.SetSubmatrix(dpdq,1,1,-1);

		GetBody3D(1).GetdRotvdqT(lrot1,lpos,dpdq);
		Mult(dpdq,vrot2,hvec); hvec*=-1;
		hmat.SetColVec(hvec,4);
	}
	else
	{
		const Vector3D& lp1   = loccoords(1);
		const Vector3D& lp2   = loccoords(2);
		const Vector3D& lrot1 = loccoords(5); //body1
		const Vector3D& lrot2 = loccoords(6); //body2

		double sign = 1;
		if (locelemind==2) sign = -1;

		GetBody3D(locelemind).GetdPosdqT(loccoords(locelemind),dpdq);
		hmat.SetSubmatrix(dpdq,1,1, sign);


		if (locelemind==1)
		{
			Vector3D vrot2 = GetBody3D(2).GetRotMatrix(lp2)*lrot2;
			GetBody3D(1).GetdRotvdqT(lrot1,lp1,dpdq);
			Mult(dpdq,vrot2,hvec);
			hmat.SetColVec(hvec,4);
		}
		else
		{
			Vector3D vrot1 = GetBody3D(1).GetRotMatrix(lp1)*lrot1;
			GetBody3D(2).GetdRotvdqT(lrot2,lp2,dpdq);
			Mult(dpdq,vrot1,hvec);
			hmat.SetColVec(hvec,4);
		}
	}

	hvec.SetLen(IS());
	hvec(1) = -XG(1);
	hvec(2) = -XG(2);
	hvec(3) = -XG(3);
	hvec(4) = -XG(4);

	for (int i=1; i <= f.Length(); i++)
	{
		f(i) -= hmat(i,1)*hvec(1)+hmat(i,2)*hvec(2)+hmat(i,3)*hvec(3)+hmat(i,4)*hvec(4);
	}
};

Vector3D UniversalJointOLD::GetRefPosD()	const 
{
	return GetBody3D(1).GetPosD(loccoords(1));
}

void UniversalJointOLD::DrawElement() 
{
	Constraint::DrawElement();
	//

	mbs->SetColor(GetCol());
	int res = (int)draw_dim.Z();
	if (res == 0) res = 16;

	double r = 0.5*draw_dim.X();
	double d = 0.25*draw_dim.Y();
	double t = 0.5*draw_dim.Y();
	double rsphere = draw_dim.Y();

	if (elements.Length()==1) 
	{
		Vector3D vrot1b = loccoords(2);
		Vector3D vrot2 = loccoords(3);
		Vector3D lrot1 = loccoords(4);
		Vector3D lrot2 = loccoords(5);
		Vector3D vrot1 = GetBody3D(1).GetRotMatrixD(loccoords(1))*lrot1;
		Vector3D vrot2b = GetBody3D(1).GetRotMatrixD(loccoords(1))*lrot2;

		Vector3D p = GetBody3D(1).GetPosD(loccoords(1));
		mbs->DrawSphere(p,1.*draw_dim.Y(),res);

		mbs->SetColor(colgrey3);
		mbs->DrawZyl(p+1.001*(r-d)*vrot1,p-1.001*(r-d)*vrot1,0.5*t,res);
		mbs->DrawZyl(p+1.001*r*vrot2,p-1.001*r*vrot2,0.5*t,res);

		//rings:
		GeomPipe3D pipe(GetMBS(),  p-0.5*t*vrot1b, t*vrot1b, vrot2, r-d, r, res, colgrey4, MY_PI, 2.*MY_PI); 
		pipe.DrawYourself();
		GeomPipe3D pipe2(GetMBS(), p-0.5*t*vrot2b, t*vrot2b, vrot1, r-2*d, r-d, res, colgrey4, MY_PI, 2.*MY_PI); 
		pipe2.DrawYourself();
	}
	else
	{
		Vector3D p1 = GetBody3D(1).GetPosD(loccoords(1));
		Vector3D p2 = GetBody3D(2).GetPosD(loccoords(2));
		const Vector3D& lrot1 = loccoords(5); //body1
		const Vector3D& lrot2 = loccoords(6); //body2
		const Vector3D& lrot1b = loccoords(7); //body2
		const Vector3D& lrot2b = loccoords(8); //body1

		Matrix3D A1 = GetBody3D(1).GetRotMatrixD(loccoords(1));
		Matrix3D A2 = GetBody3D(2).GetRotMatrixD(loccoords(2));

		Vector3D vrot1 = A1*lrot1;
		Vector3D vrot2 = A2*lrot2;
		Vector3D vrot1b= A2*lrot1;
		Vector3D vrot2b= A1*lrot2;

		mbs->DrawSphere(p1,1.*draw_dim.Y()*0.99,res*0.3);
		mbs->DrawSphere(p2,1.*draw_dim.Y(),res);

		mbs->SetColor(colgrey3);
		mbs->DrawZyl(p1+1.001*(r-d)*vrot1,p1-1.001*(r-d)*vrot1,0.5*t,res);
		mbs->DrawZyl(p2+1.001*r*vrot2,p2-1.001*r*vrot2,0.5*t,res);

		//rings:
		GeomPipe3D pipe(GetMBS(),  p1-0.5*t*vrot1b, t*vrot1b, vrot2, r-d, r, res, colgrey4, 0*MY_PI, 1.*MY_PI); 
		pipe.DrawYourself();
		GeomPipe3D pipe2(GetMBS(), p2-0.5*t*vrot2b, t*vrot2b, vrot1, r-2*d, r-d, res, colgrey4, 0*MY_PI, 1.*MY_PI); 
		pipe2.DrawYourself();
	}

};







//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// PrismaticJoint2D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void PrismaticJoint2D::Initialize() 
{
	//find orthogonal vectors:
	const Vector2D& lp1 = loccoords(1);
	const Vector2D& lp2 = loccoords(2);
	Vector2D t1, n2;
	Vector2D gp1 = GetBody2D(1).GetPos2D(lp1);
	Vector2D gp2 = lp2; //for ground joints, local = global
	if (elements.Length() == 2) 
		gp2 = GetBody2D(2).GetPos2D(lp2);

	t1 = gp2 - gp1;
	n2 = Vector2D(-t1.Y(), t1.X());


	Matrix3D RT1=GetBody2D(1).GetRotMatrix2D().GetTp();
	Matrix3D RT2;
	RT2.Set22(1,0,0,1);
	if (elements.Length() == 2) 
		RT2=GetBody2D(2).GetRotMatrix2D().GetTp();

	t1 = RT1*t1;
	n2 = RT2*n2;

	loccoords(3) = t1; 
	loccoords(4) = n2;

};

void PrismaticJoint2D::EvalG(Vector& f, double t) 
{
	if (elements.Length()<1 || elements.Length()>2)
	{
		mbs->UO() << "ERROR: PrismaticJoint2D::EvalG, number of elements != 1 or 2\n"; return;
	}

	if (MaxIndex()==3)
	{
		if (elements.Length()==1)
		{
			const Vector2D& lp1 = loccoords(1);
			const Vector2D& gp2 = loccoords(2);
			const Vector2D& lt1 = loccoords(3);
			const Vector2D& ln2 = loccoords(4);

			Vector2D rpij = GetBody2D(1).GetPos2D(lp1)-gp2;

			Matrix3D A = GetBody2D(1).GetRotMatrix2D(); 
			Vector2D t1 = A*lt1;

			f(1) = rpij * ln2;
			f(2) = t1 * ln2;
		}
		else
			if (elements.Length()==2) 
			{
				const Vector2D& lp1 = loccoords(1);
				const Vector2D& lp2 = loccoords(2);
				const Vector2D& lt1 = loccoords(3);
				const Vector2D& ln2 = loccoords(4);

				Vector2D rpij = GetBody2D(1).GetPos2D(lp1) - GetBody2D(2).GetPos2D(lp2);

				Matrix3D A1 = GetBody2D(1).GetRotMatrix2D(); 
				Matrix3D A2 = GetBody2D(2).GetRotMatrix2D(); 
				Vector2D t1 = A1*lt1;
				Vector2D n2 = A2*ln2;

				f(1) = rpij * n2;
				f(2) = t1 * n2;

			}
	}
	else if (MaxIndex()<=2)
	{
		if (elements.Length()==1)
		{
			const Vector2D& lp1 = loccoords(1);
			const Vector2D& gp2 = loccoords(2);
			const Vector2D& lt1 = loccoords(3);
			const Vector2D& ln2 = loccoords(4);

			//Vector2D rpij = GetBody2D(1).GetPos2D(lp1)-gp2;
			Vector2D vpij = GetBody2D(1).GetVel2D(lp1);

			//Matrix3D A = GetBody2D(1).GetRotMatrix2D(); 
			Matrix3D Ap = GetBody2D(1).GetRotMatrix2DP();

			Vector2D tp1 = Ap*lt1;

			f(1) = vpij * ln2;
			f(2) = tp1 * ln2;
		}
		else if (elements.Length()==2) 
		{
			const Vector2D& lp1 = loccoords(1);
			const Vector2D& lp2 = loccoords(2);
			const Vector2D& lt1 = loccoords(3);
			const Vector2D& ln2 = loccoords(4);

			Vector2D rpij = GetBody2D(1).GetPos2D(lp1) - GetBody2D(2).GetPos2D(lp2);
			Vector2D vpij = GetBody2D(1).GetVel2D(lp1) - GetBody2D(2).GetVel2D(lp2);

			Matrix3D A1 = GetBody2D(1).GetRotMatrix2D(); 
			Matrix3D A1p = GetBody2D(1).GetRotMatrix2DP(); 
			Matrix3D A2 = GetBody2D(2).GetRotMatrix2D(); 
			Matrix3D A2p = GetBody2D(2).GetRotMatrix2DP(); 
			Vector2D t1 = A1*lt1;
			Vector2D t1p = A1p*lt1;
			Vector2D n2 = A2*ln2;
			Vector2D n2p = A2p*ln2;

			f(1) = vpij * n2 + rpij * n2p;
			f(2) = t1p * n2 + t1 * n2p;
			//UO() << "f=" << f << "\n";
		}
	}
};

void PrismaticJoint2D::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	hmat.SetSize(f.Length(),2);

	if (elements.Length() == 1) //--> ignore locelemind, only first element!
	{
		const Vector2D& lp1 = loccoords(1);
		const Vector2D& gp2 = loccoords(2);
		const Vector2D& lt1 = loccoords(3);
		const Vector2D& ln2 = loccoords(4);

		Vector2D rpij = GetBody2D(1).GetPos2D(lp1) - gp2;

		Matrix3D A = GetBody2D(1).GetRotMatrix2D(); 

		Vector2D t1 = A*lt1;

		GetBody2D(1).GetdPosdqT(lp1,dpdq);
		Mult(dpdq,ln2,hvec);
		hmat.AddColVec(1,hvec); //for first constraint

		GetBody2D(1).GetdRotvdqT(lt1,lp1,dpdq);
		Mult(dpdq,ln2,hvec);
		hmat.SetColVec(hvec,2); //for second constraint
	}
	else
	{
		const Vector2D& lp1 = loccoords(1);
		const Vector2D& lp2 = loccoords(2);
		const Vector2D& lt1 = loccoords(3);
		const Vector2D& ln2 = loccoords(4);

		if (locelemind==1)
		{
			Vector2D rpij = GetBody2D(1).GetPos2D(lp1) - GetBody2D(2).GetPos2D(lp2);

			Matrix3D A1 = GetBody2D(1).GetRotMatrix2D(); 
			Matrix3D A2 = GetBody2D(2).GetRotMatrix2D(); 

			Vector2D t1 = A1*lt1;
			Vector2D n2 = A2*ln2;

			GetBody2D(1).GetdPosdqT(lp1,dpdq);
			Mult(dpdq,n2,hvec);
			hmat.SetColVec(hvec,1); //for first constraint

			GetBody2D(1).GetdRotvdqT(lt1,lp1,dpdq);
			Mult(dpdq,n2,hvec);
			hmat.SetColVec(hvec,2); //for second constraint

			//f(1) = rpij * n2;
			//f(2) = t1 * n2;
		}
		else
		{
			Vector2D rpij = GetBody2D(1).GetPos2D(lp1) - GetBody2D(2).GetPos2D(lp2);

			Matrix3D A1 = GetBody2D(1).GetRotMatrix2D(); 
			Matrix3D A2 = GetBody2D(2).GetRotMatrix2D(); 

			Vector2D t1 = A1*lt1;
			Vector2D n2 = A2*ln2;

			GetBody2D(2).GetdPosdqT(lp2,dpdq);
			Mult(dpdq,n2,hvec); hvec *= -1.;
			hmat.SetColVec(hvec,1); //for first constraint

			GetBody2D(2).GetdRotvdqT(ln2,lp2,dpdq);
			Mult(dpdq,rpij,hvec);
			hmat.AddColVec(1,hvec); //for first constraint
			Mult(dpdq,t1,hvec);
			hmat.SetColVec(hvec,2); //for second constraint
		}
	}

	for (int i=1; i <= f.Length(); i++)
	{
		f(i) -= hmat(i,1)*XG(1)+hmat(i,2)*XG(2);
	}
};

Vector3D PrismaticJoint2D::GetRefPosD()	const 
{
	return GetBody2D(1).ToP3D(GetBody2D(1).GetPos2DD(loccoords(1)));
}

void PrismaticJoint2D::DrawElement() 
{
	Constraint::DrawElement();
	//


	mbs->SetColor(GetCol());
	//...
	if (draw_dim.X() == 0) return;
	{

		Vector3D p1 = GetBody2D(1).ToP3D(GetBody2D(1).GetPos2DD(loccoords(1)));

		Vector3D rot1, rot2;

		rot1 = GetBody2D(1).ToP3D(GetBody2D(1).GetRotMatrix2DD()*loccoords(3));
		rot1.Normalize();
		rot1 *= draw_dim.Y();
		mbs->DrawZyl(p1,p1+rot1,0.5*draw_dim.X(),12);

		if (elements.Length()==2)
		{
			Vector3D p2 = GetBody2D(2).ToP3D(GetBody2D(2).GetPos2DD(loccoords(2)));
			rot2 = GetBody2D(2).ToP3D(GetBody2D(2).GetRotMatrix2DD()*loccoords(4));
			rot2.Normalize();
			rot2 *= draw_dim.Y();

			mbs->DrawZyl(p2,p2+rot2,0.5*draw_dim.X(),12);
		}
	}
};













//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// RollingJoint3D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void RollingJoint3D::ElementDefaultConstructorInitialization()
{
	Constraint::ElementDefaultConstructorInitialization();

	dpdq.SetSize(0,0);
	x_init = Vector(SS());

	Vector xdatainit(0.,0.); //sticking+is_contact
	SetDataInit(xdatainit);

	//standard setting uses idealized rolling!
	rolling_plane_point = Vector3D(0.,0.,0.);
	rolling_plane_normal = Vector3D(0.,0.,1.); //default: rolling in xy-plane

	wheel_local_center_point = Vector3D(0.,0.,0.);
	wheel_local_axis = Vector3D(1.,0.,0.);
	wheel_radius = 0; //radius of the idealized wheel

	use_penalty_inplane_dir = 0;
	use_penalty_contact = 0;

	penalty_stiffness_inplane = 1e3; //penalty stiffness in direction of rolling
	penalty_stiffness_contact = 1e3; //penalty stiffness in contact direction (contact wheel-ground)

	rolling_friction = 0;

	use_friction = 0;
	friction_coeff_inplane = 0; //friction coeffiction in rolling direction

	contact_damping = 0;
	use_contact_condition = 0;

	elements.Set1(1); // $ MSax 2013-03-08: added

	draw_dim.X() = 0.05; // $ MSax 2013-05-06: draw size sphere
	draw_dim.Y() = 0.00001; // $ MSax 2013-05-06: force scaling factor; 0.01mm/N
}

void RollingJoint3D::Initialize() 
{
	wheel_local_axis.Normalize(); // $ MSax 2013-03-08: added
	rolling_plane_normal.Normalize(); // $ MSax 2013-03-08: added

	if (!use_penalty_contact) // $ MSax 2013-06-17: lagrange method for contact: always in contact, always sticking
	{
		IsContact() = 1;
		IsSticking() = 1;
		use_contact_condition = 0;
		use_friction = 0;
	}
}

int RollingJoint3D::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute // $ MSax 2013-05-07: added
{
	int rv = Constraint::CheckConsistency(errorstr);
	if (rv) return rv;

	//if (TODO:check if initial states are correct)
	//{
	//		errorstr = mystr("ERROR: RollingJoint3D: wrong initial states for 'sticking' or 'contact'!\n");
	//		rv = 1;
	//}
	return rv;
}

int RollingJoint3D::GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) //$ MSax 2013-05-06: added
{
	// Automatic entries for this class 
	RollingJoint3D::GetAvailableSpecialValuesAuto(available_variables);

	// Manual READ entries for this class
	available_variables.Add(ReadWriteElementDataVariableType("Connector.force_forward",0,0,0.,mystr("force component in forward direction"), TRWElementDataRead));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.force_lateral",0,0,0.,mystr("force component in lateral direction"), TRWElementDataRead));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.force_normal",0,0,0.,mystr("contact force normal to plane"), TRWElementDataRead));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.sticking",0,0,0.,mystr("returns 1 if sticking, else 0"), TRWElementDataRead));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.contact",0,0,0.,mystr("returns 1 if in contact, else 0"), TRWElementDataRead));

	return 0;
}

int RollingJoint3D::ReadSingleElementData(ReadWriteElementDataVariableType& RWdata) //$ MSax 2013-05-06: added	
{
	// call base class routine 	
	int rv = Constraint::ReadSingleElementData(RWdata);
	if (rv == 1) return 1;

	//manual things to read  
	if(RWdata.variable_name.CStrCompare("Connector.force_forward"))
	{
		RWdata.value = XG(1);
		return 1;
	}
	if(RWdata.variable_name.CStrCompare("Connector.force_lateral"))
	{
		RWdata.value = XG(2);
		return 1;
	}
	if(RWdata.variable_name.CStrCompare("Connector.force_normal"))
	{
		RWdata.value = abs(GetContactForce()); // $ MSax 2013-06-18: changed to absolute value of contact force
		return 1;
	}
	if(RWdata.variable_name.CStrCompare("Connector.contact"))
	{
		RWdata.value = IsContact()*5000;
		return 1;
	}
	if(RWdata.variable_name.CStrCompare("Connector.sticking"))
	{
		RWdata.value = IsSticking()*5000;
		return 1;
	}
	return ReadSingleElementDataAuto(RWdata);
}

void RollingJoint3D::EvalG(Vector& f, double t) 
{

	//2 Bedingungen in Vorwärtsrichtung?
	//if (MaxIndex()<=3)
	//{
	//	Matrix3D AT = GetBody3D(1).GetRotMatrix();
	//	AT.TpYs();
	//	Vector3D locpos = AT*loccoords(1);

	//	Vector3D v = GetBody3D(1).GetVel(locpos);
	//	f = v;
	//}

	//int oo = 0; //show output

	//contact must be distiguished between index 2 and 3!
	Vector3D lcp;
	Vector3D forward_dir;
	Vector3D lateral_dir;
	
	double act_radius = WheelLocalContactPosition(lcp, forward_dir, lateral_dir);
	//UO() << "lcp=" << lcp << "\n";
	//if (oo) UO() << "forward_dir=" << forward_dir << "\n";
	//if (oo) UO() << "lateral_dir=" << lateral_dir << "\n";
	Vector3D up_dir = forward_dir.Cross(lateral_dir);

	Vector3D wheelcenterp = GetBody3D(1).GetPos(wheel_local_center_point);
	Vector3D wheelv = lcp - wheelcenterp; //vector which points from wheel center point to contact point
	wheelv.Normalize();
	Vector3D wheelp = wheelcenterp + wheel_radius * wheelv;	//point at wheel at original radius (penetrated)

	Vector3D pp = wheelp;
	ProjectInPlane(rolling_plane_point, rolling_plane_normal, pp);

	double contact_penetration = (pp-wheelp) * rolling_plane_normal;

	//if (oo) UO() << "wheelcenterp=" << wheelcenterp << "\n";
	//if (oo) UO() << "wheelv=" << wheelv << "\n";
	//if (oo) UO() << "wheelp=" << wheelp << "\n";
	//if (oo) UO() << "pp=" << pp << "\n";
	//if (oo) UO() << "contact_penetration=" << contact_penetration << "\n";

	//compute local position of contact point:
	Vector3D v = lcp-GetBody3D(1).GetRefPos();
	Vector3D vloc = GetBody3D(1).GetRotMatrix().GetTp()*v;
	//if (oo) UO() << "v=" << v << "\n";
	//if (oo) UO() << "vloc=" << vloc << "\n";

	//compute velocity of contact point ==> should be zero because of friction condition!
	Vector3D contact_point_vel = GetBody3D(1).GetVel(vloc);
	//Vector3D contact_point_acc = GetBody3D(1).GetAcceleration(vloc);

	double contact_velocity = contact_point_vel * rolling_plane_normal;

	//UO() << "contact_point_vel=" << contact_point_vel << "\n";

	//Vector3D wheel_ang_vel = GetBody3D(1).GetAngularVel(wheel_local_center_point);
	//Vector3D forward_vel = wheel_radius * wheel_ang_vel.Cross(rolling_plane_normal)
	//	- GetBody3D(1).GetVel(Vector3D(0.,0.,0.));

	//UO() << "contact_point_vel=" << contact_point_vel << "\n";
	//UO() << "forward_vel_CP=   " << forward_vel << "\n";

	if (IsContact() || !use_contact_condition)
	{
		if (IsSticking() || !use_friction)
		{
			if (use_penalty_inplane_dir)
			{
				f(1) = XG(1) - (contact_point_vel * forward_dir) * penalty_stiffness_inplane;
			}
			else
			{
				f(1) = contact_point_vel * forward_dir;
				//f(1) = forward_vel * forward_dir;
			}
			if (use_penalty_inplane_dir)
			{
				f(2) = XG(2) - (contact_point_vel * lateral_dir) * penalty_stiffness_inplane;
			}
			else
			{
				f(2) = contact_point_vel * lateral_dir;
			}
		}
		else //not sticking ==> friction force inplane against sliding velocity
		{
			
			if (contact_point_vel.Norm() != 0)
			{
				Vector3D friction_dir = contact_point_vel;
				friction_dir.Z() = 0; // $ MSax 2013-06-17 : added
				friction_dir.Normalize();
				double FN = abs(GetContactForce());  // $ MSax 2013-06-18 : changed to absolut value of contact force
				
				Vector3D friction_force = (-FN*friction_coeff_inplane)*friction_dir;

				f(1) = XG(1) - friction_force*forward_dir;
				f(2) = XG(2) - friction_force*lateral_dir;
			}
			else //no sliding velocity ==> zero force
			{
				f(1) = XG(1);
				f(2) = XG(2);
			}

		}

		if (use_penalty_contact)
		{
			f(3) = XG(3) + (contact_penetration * penalty_stiffness_contact - contact_velocity * contact_damping);
		}
		else
		{
			if (MaxIndex()==3)
			{
				f(3) = -contact_penetration;
			}
			else if (MaxIndex()<=2)
			{
				f(3) = contact_point_vel * up_dir;//contact_penetration;
			}
			else
			{
				UO() << "ERROR: RollingJoint3D: illegal index! \n";
			}
		}
	}
	else //no contact ==> contact forces are zero
	{
		f(1) = XG(1);
		f(2) = XG(2);
		f(3) = XG(3);
	}

	//f(1) = XG(1);
	//f(2) = XG(2);
	//f(3) = XG(3);
	//UO() << "f=" << f << "\n";
};

void RollingJoint3D::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	Vector3D force_local = ComputeForce(); //in rolling coordinates!

	Vector3D lcp;
	Vector3D forward_dir;
	Vector3D lateral_dir;
	
	double act_radius = WheelLocalContactPosition(lcp, forward_dir, lateral_dir);

	Vector3D contact_dir = rolling_plane_normal;

	Vector3D force = force_local.X()*forward_dir + force_local.Y()*lateral_dir + force_local.Z()*contact_dir;

	//UO() << "force=" << force << "\n";
	//UO() << "force_local=" << force_local << "\n";

	//compute local position of contact point:
	Vector3D v = lcp-GetBody3D(1).GetRefPos();
	Vector3D vloc = GetBody3D(1).GetRotMatrix().GetTp()*v;
	//UO() << "v=" << v << ", ";
	//UO() << "vloc=" << vloc << "\n";

	//vloc = Vector3D(0.,0.,0.); //test
	GetBody3D(1).GetdPosdqT(vloc,dpdq);
	//GetBody3D(1).GetdPosdqT(Vector3D(0.,0.,0.),dpdq);

	for (int i=1; i <= f.Length(); i++)
	{
		f(i) -= (dpdq(i,1)*force.X()+dpdq(i,2)*force.Y()+dpdq(i,3)*force.Z());
	}
	//UO() << "f=" << f << "\n";
	//UO() << "dpdq=" << dpdq << "\n";

	//rolling friction:
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//approximated: force against rolling direction, applied to wheel centerpoint

	//Vector3D wheel_vel = GetBody3D(1).GetVel(wheel_local_center_point);
	//double friction_force = rolling_friction * (forward_dir * wheel_vel);

	//force = forward_dir*Vector3D(friction_force, 0., 0.);
	////vloc = Vector3D(0.,0.,-0.2);
	//GetBody3D(1).GetdPosdqT(wheel_local_center_point,dpdq);

	//for (int i=1; i <= f.Length(); i++)
	//{
	//	f(i) -= (dpdq(i,1)*force.X()+dpdq(i,2)*force.Y()+dpdq(i,3)*force.Z());
	//}

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//BETTER: apply moment against angular velocity at wheelpoint
	// $ MSax 2013-06-17 : removed rolling friction
	//GetBody3D(1).GetdRotdqT(wheel_local_center_point, dpdq);
	//Vector3D wheel_ang_vel = GetBody3D(1).GetAngularVel(wheel_local_center_point);

	//Vector3D torque = rolling_friction * wheel_ang_vel;

	//for (int i=1; i <= f.Length(); i++)
	//{
	//	f(i) -= (dpdq(i,1)*torque.X()+dpdq(i,2)*torque.Y()+dpdq(i,3)*torque.Z());
	//}


};

double RollingJoint3D::PostNewtonStep(double t)
{
	double error = 0;
	
	if (use_penalty_contact)
	{
		Vector3D force_local = ComputeForce();

		Vector3D lcp;
		Vector3D forward_dir;
		Vector3D lateral_dir;

		double act_radius = WheelLocalContactPosition(lcp, forward_dir, lateral_dir);
		Vector3D contact_dir = rolling_plane_normal;
		Vector3D force = force_local.X()*forward_dir + force_local.Y()*lateral_dir + force_local.Z()*contact_dir;

		//compute local position of contact point:
		Vector3D v = lcp-GetBody3D(1).GetRefPos();
		Vector3D vloc = GetBody3D(1).GetRotMatrix().GetTp()*v;

		Vector3D contact_point_vel = GetBody3D(1).GetVel(vloc);

		//Vector3D up_dir = forward_dir.Cross(lateral_dir);

		Vector3D wheelcenterp = GetBody3D(1).GetPos(wheel_local_center_point);
		Vector3D wheelv = lcp - wheelcenterp; //vector which points from wheel center point to contact point
		wheelv.Normalize();
		Vector3D wheelp = wheelcenterp + wheel_radius * wheelv;	//point at wheel at original radius (penetrated)

		Vector3D pp = wheelp;
		ProjectInPlane(rolling_plane_point, rolling_plane_normal, pp);

		double contact_penetration = (pp-wheelp) * rolling_plane_normal;

		// change the state of contact if necessary
		if (contact_penetration >= 0) 
		{
			if (!(int)IsContact()) error = fabs(contact_penetration); //if first step with contact
			IsContact() = 1;
		}
		else
			//better: GetContactForce() < 0
		{
			if ((int)IsContact()) error = fabs(contact_penetration);
			IsContact() = 0;
		}

		// change the state of sticking if necessary
		Vector3D surf_vel = contact_point_vel;
		surf_vel.Z() = 0;

		if(surf_vel.Norm2() < 1e-3 && IsContact())
		{
			//if(!(int)IsSticking()) error = error + fabs(surf_vel.Norm2());
			IsSticking() = 1;
		}

		if(!(int)IsContact())
		{
			IsSticking() = 0;
		}
		else //contact
		{
			double FN = abs(GetContactForce());
			double max_friction_force = (FN*friction_coeff_inplane);
			Vector3D friction_force_local(XG(1),XG(2),0.);

			if (friction_force_local.Norm() > max_friction_force)
			{
				//if((int)IsSticking()) error = error + fabs(friction_force_local.Norm() - max_friction_force);
				IsSticking() = 0; // $ MSax 2013-05-06: changed from 1 to 0, otherwise it makes no sense
			}
		}

		//	Vector3D surf_vel = contact_point_vel;
		//	surf_vel.Z() = 0;
		//	Vector3D friction_force_local(XG(1),XG(2),0.);

		//	double val_cross = (surf_vel.Cross(friction_force_local)).Norm2();
		//	double val = friction_force_local.Norm2();

		//	if((val_cross > 1e-10 || val < 1e-10) && IsContact())
		//	{
		//		IsSticking() = 1;
		//	}
	}

	return error;
}

void RollingJoint3D::PostprocessingStep()
{
}

Vector3D RollingJoint3D::GetRefPosD()	const 
{
	return GetBody3D(1).GetRefPosD();
}

void RollingJoint3D::DrawElement() 
{
	Constraint::DrawElement();
	//

	mbs->SetColor(GetCol());

	Vector3D lcp;
	Vector3D forward_dir;
	Vector3D lateral_dir;
	
	WheelLocalContactPositionD(lcp, forward_dir, lateral_dir);
	Vector3D up_dir = forward_dir.Cross(lateral_dir);

	Vector3D v = lcp-GetBody3D(1).GetRefPosD();
	Vector3D vloc = GetBody3D(1).GetRotMatrixD().GetTp()*v;
	mbs->DrawSphere(lcp,0.5*draw_dim.X(),6);

	//double loadsizefactor = mbs->GetDOption(120);

	Vector3D F_forward = draw_dim(2)*XGD(1)*forward_dir;
	Vector3D F_lateral = draw_dim(2)*XGD(2)*lateral_dir;
	Vector3D F_contact = draw_dim(2)*GetContactForceD()*up_dir;

	mbs->SetColor(colblue);
	mbs->MyDrawArrow(lcp, lcp+F_forward, GetCol());
	mbs->SetColor(colgreen);
	mbs->MyDrawArrow(lcp, lcp+F_lateral, GetCol());
	mbs->SetColor(colred);
	mbs->MyDrawArrow(lcp, lcp+F_contact, GetCol());

	//mbs->SetColor(colgrey3);
	//mbs->MyDrawArrow(lcp, lcp+contact_point_vel, GetCol());


	mbs->SetColor(colgrey2);
	mbs->DrawSphere(GetBody3D(1).GetPosD(vloc),0.5*draw_dim.X(),5);

	//mbs->DrawSphere(GetBody3D(1).GetRefPosD() + loccoords(1),0.5*draw_dim.X(),16);
};

Vector3D RollingJoint3D::ComputeForce() const
{
	Vector3D force(-1*XG(1),-1*XG(2),XG(3)); //x=rolling direction, y=lateral direction, z=contact direction

	return force;
}

double RollingJoint3D::WheelLocalContactPosition(Vector3D& contact_point, Vector3D& forward_dir, Vector3D& lateral_dir) const
{
	int oo = 0;
  Matrix3D A = GetBody3D(1).GetRotMatrix();
	if (oo) UO() << "A=" << A << "\n";

	//transform wheel coordinates to global coordinates:
	Vector3D wheelp = GetBody3D(1).GetPos(wheel_local_center_point);
	Vector3D wheelaxis = A*wheel_local_axis;

	if (oo) UO() << "rolling_plane_normal=" << rolling_plane_normal << "\n";
	if (oo) UO() << "rolling_plane_point=" << rolling_plane_point << "\n";
	if (oo) UO() << "wheelp=" << wheelp << "\n";
	if (oo) UO() << "wheelaxis=" << wheelaxis << "\n";

	Vector3D linepoint, linevector;

	int rv = PlaneIntersection(wheelaxis, wheelp, rolling_plane_normal, rolling_plane_point, linepoint, linevector);
	if (oo) UO() << "linepoint=" << linepoint << "\n";
	if (oo) UO() << "linevector=" << linevector << "\n";
	if (oo) UO() << "rv=" << rv << "\n";

	if (rv == 1)
	{
		double act_radius = DistToLine(linepoint, linevector, wheelp, contact_point);  

		forward_dir = rolling_plane_normal.Cross(wheelaxis);
		forward_dir.Normalize();
		//??forward_dir should point in direction of movement (this is advantageous for friction force determination)
		//not so good: if (forward_dir * GetBody3D(1).GetVel(contact_point_local) < 0.) {forward_dir *= -1;}

		lateral_dir = rolling_plane_normal.Cross(forward_dir);
		lateral_dir.Normalize();

		return act_radius;
	}
	else
	{
		UO() << "RollingJoint: wheel is lying\n";
	}

	//wheel axis and plane normal are colinear; wheel is lying; contact point=wheelp
	contact_point = Vector3D(wheelp);
	wheelaxis.SetNormalBasis(forward_dir, lateral_dir);

	return 0.; //wheel is lying
}

double RollingJoint3D::WheelLocalContactPositionD(Vector3D& contact_point, Vector3D& forward_dir, Vector3D& lateral_dir) const
{
	Matrix3D A = GetBody3D(1).GetRotMatrixD();

	//transform wheel coordinates to global coordinates:
	Vector3D wheelp = GetBody3D(1).GetPosD(wheel_local_center_point);
	Vector3D wheelaxis = A*wheel_local_axis;

	Vector3D linepoint, linevector;

	int rv = PlaneIntersection(wheelaxis, wheelp, rolling_plane_normal, rolling_plane_point, linepoint, linevector);

	if (rv == 1)
	{
		double act_radius = DistToLine(linepoint, linevector, wheelp, contact_point);  

		forward_dir = rolling_plane_normal.Cross(wheelaxis);
		//forward_dir should point in direction of movement (this is advantageous for friction force determination)
		//if (forward_dir * GetBody3D(1).GetVelD(contact_point) < 0.) {forward_dir *= -1;}
		forward_dir.Normalize();

		lateral_dir = rolling_plane_normal.Cross(forward_dir);
		lateral_dir.Normalize();

		return act_radius;
	}
	contact_point = wheelp;
	wheelaxis.SetNormalBasis(forward_dir, lateral_dir);

	return 0.; //wheel is lying
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// RollingSegment3D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void RollingSegment3D::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	Vector3D proll = GetBody3D(1).GetPos(loccoords(1));
	Vector3D pCOM = GetBody3D(1).GetRefPos(); //center of mass
	Vector3D r = Vector3D(0.,0.,-rollrad);

	Vector3D v = pCOM-proll;

	double phi = acos((v*r)/(rollrad*v.Norm()));

	if (phi > MY_PI) phi -= 2.*MY_PI;
	if (phi < -MY_PI) phi += 2.*MY_PI;

	if (fabs(phi) <= 0.5*rollphi)
	{
		Vector3D t = r-v;


		Matrix3D AT = GetBody3D(1).GetRotMatrix();
		AT.TpYs();
		Vector3D tloc = AT*t;

		GetBody3D(1).GetdPosdqT(tloc,dpdq);
		Vector3D lambda = (rollstiff) * GetBody3D(1).GetVel(tloc);


		for (int i=1; i <= f.Length(); i++)
		{
			f(i) -= (dpdq(i,1)*lambda.X()+dpdq(i,2)*lambda.Y());
		}
	}

};

Vector3D RollingSegment3D::GetRefPosD()	const 
{
	return GetBody3D(1).GetRefPosD();
}

void RollingSegment3D::DrawElement() 
{
	Constraint::DrawElement();
	//


	mbs->SetColor(GetCol());

	Vector3D v = GetBody3D(1).GetPosD(loccoords(1));
	Vector3D v2 = Vector3D(0.,0.,-rollrad);

	mbs->DrawSphere(v+v2,0.5*draw_dim.X(),16);
};





//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// FlexibleRollingCable3D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void FlexibleRollingCable3D::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{

	Vector3D l1 = loccoords(1);
	Vector3D l2 = loccoords(2);

	for (int i=1; i <= nseg; i++)
	{
		double fact = 1./(double)nseg;
		double t = ((double)i-0.5)/(double)nseg;
		Vector3D lp = (1.-t)*l1 + t*l2;
		Vector3D v = GetBody3D(1).GetPos(lp);

		if (v.Z() <= rollz) 
		{
			Body3D& ce = GetBody3D(1);


			ce.GetdPosdqT(lp,dpdq);
			Vector3D lambda = fact * rollstiff * ce.GetVel(lp);


			for (int i=1; i <= f.Length(); i++)
			{
				f(i) -= (dpdq(i,1)*lambda.X()+dpdq(i,2)*lambda.Y());
			}

		}
	}

};

Vector3D FlexibleRollingCable3D::GetRefPosD()	const 
{
	return GetBody3D(1).GetRefPosD();
}

void FlexibleRollingCable3D::DrawElement() 
{
	Constraint::DrawElement();
	//


	mbs->SetColor(GetCol());

	Vector3D l1 = loccoords(1);
	Vector3D l2 = loccoords(2);

	for (int i=1; i <= nseg; i++)
	{
		double t = ((double)i-0.5)/(double)nseg;
		Vector3D v = GetBody3D(1).GetPosD((1.-t)*l1 + t*l2);

		if (v.Z() <= rollz) 
		{
			mbs->DrawSphere(v,0.5*draw_dim.X(),16);
		}
	}
};





//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// RollingJoint2D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void RollingJoint2D::Initialize() 
{
};

void RollingJoint2D::EvalG(Vector& f, double t) 
{
	if (elements.Length()<1 || elements.Length()>2)
	{
		mbs->UO() << "ERROR: RollingJoint2D::EvalG, number of elements != 1 or 2\n"; return;
	}

	if (MaxIndex()<=3)
	{
		const Vector2D& lp1 = loccoords(1);
		const Vector2D& lt1 = loccoords(2);

		Matrix3D A = GetBody2D(1).GetRotMatrix2D(); 
		Vector2D v = GetBody2D(1).GetVel2D(lp1);
		Vector2D t = A*lt1;
		Vector2D n(-t.Y(), t.X());

		f(1) = v * n;
	}
};

void RollingJoint2D::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	const Vector2D& lp1 = loccoords(1);
	const Vector2D& lt1 = loccoords(2);

	Matrix3D A = GetBody2D(1).GetRotMatrix2D(); 
	Vector2D tv = A*lt1;
	Vector2D n(-tv.Y(), tv.X());

	GetBody2D(1).GetdPosdqT(lp1,dpdq);

	for (int i=1; i <= f.Length(); i++)
	{
		f(i) -= (dpdq(i,1)*n.X() + dpdq(i,2)*n.Y())*XG(1);
	}
};

Vector3D RollingJoint2D::GetRefPosD()	const 
{
	return GetBody2D(1).ToP3D(GetBody2D(1).GetPos2DD(loccoords(1)));
}

void RollingJoint2D::DrawElement() 
{
	Constraint::DrawElement();
	//


	mbs->SetColor(GetCol());
	//...
	if (draw_dim.X() == 0) return;
	{

		Vector3D p1 = GetBody2D(1).ToP3D(GetBody2D(1).GetPos2DD(loccoords(1)));

		Vector3D rot1, rot2;

		rot1 = GetBody2D(1).ToP3D(GetBody2D(1).GetRotMatrix2DD()*loccoords(2));
		rot1.Normalize();
		rot1 *= draw_dim.Y()*0.5;

		mbs->DrawZyl(p1-rot1,p1+rot1,0.5*draw_dim.X(),12);
	}
};









//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// SlidingJoint
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SlidingJoint::SlidingJoint(MBS* mbsi, int en1, int en2, 
													 const Vector3D& lc1, const Vector3D& lc2, double cdim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq() 
{	
	ElementDefaultConstructorInitialization();

	x_init = Vector(SS()); // s, lambda1, lambda2, lambda3
	x_init(1) = lc2.X();
	GetCol() = coli;
	draw_dim.X() = cdim;
	AddElementCoord(en1, lc1);
	AddElementCoord(en2, lc2);
	elemind = 1; //sliding element, 1 means elements(2)
	xoff = 0;
	ind_init = elemind;
};

SlidingJoint::SlidingJoint(MBS* mbsi, int en1, const TArray<int>& en2, int in2, 
													 const Vector3D& lc1, const Vector3D& lc2, double cdim, 
													 const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq() 
{	
	ElementDefaultConstructorInitialization();

	x_init = Vector(SS()); // s, lambda1, lambda2, lambda3
	x_init(1) = lc2.X();
	GetCol() = coli;
	draw_dim.X() = cdim;
	AddElementCoord(en1, lc1);
	AddElementCoord(en2(1), lc2);
	for (int i=2; i <= en2.Length(); i++)
	{
		AddElementCoord(en2(i), lc2);
	}
	elemind = in2;
	ind_init = elemind;
	xoff = 0;
	double x = lc2.X();
	int n = 1;
	while (n < in2) 
	{
		xoff += 0.5*GetBody3D(n+1).GetSize().X(); 
		xoff += 0.5*GetBody3D(n+2).GetSize().X(); 
		n++;
	}

};

void SlidingJoint::ElementDefaultConstructorInitialization()
{
	noforcey = 0;
	nlstepcnt = 0;
	iscontact = 1;
	iscontact2 = 1;
	driven = 0;
	driveendtime = 1e10;
	drivenextstart = 0;
	driveflag = 0;
	drivevel = 0;
	drivelen = 0;
	driveacceltime = 0;
	elementname = GetElementSpec();
}

void SlidingJoint::CopyFrom(const Element& e)
{
	Constraint::CopyFrom(e);
	const SlidingJoint& ce = (const SlidingJoint&)e;
	loccoords = ce.loccoords;
	dpdq = ce.dpdq;
	elemind = ce.elemind;
	xoff = ce.xoff;
	nlstepcnt = ce.nlstepcnt;
	ind_init = ce.ind_init;
	iscontact = ce.iscontact;
	iscontact2 = ce.iscontact2;
	noforcey = ce.noforcey;

	driven = ce.driven;
	drivestarttime = ce.drivestarttime;
	driveendtime = ce.driveendtime;
	drivevel = ce.drivevel;
	drivelen = ce.drivelen;
	driveflag = ce.driveflag;
	drivenextstart = ce.drivenextstart;
	driveacceltime = ce.driveacceltime;
}

void SlidingJoint::EvalF(Vector& f, double t) 
{
	Element::EvalF(f,t);
	f(ns) = XG(nsp);
}; 

void SlidingJoint::PostprocessingStep() 
{
	if (!driven) return;

	if (drivenextstart == 0)
	{
		if (XG(ns) >= drivelen && !driveflag)
		{
			driveflag = 1;
		}
	}
	else 
	{
		if (XG(ns) >= drivelen && !driveflag)
		{
			driveflag = 1;
		}

		//very special for fluid problem!!!!!! WARNING ??????
		//UO() << "drivestarttime=" << drivestarttime << "\n";
		if (GetMBS()->GetTime()+1*GetMBS()->GetStepSize() >= drivestarttime && driveacceltime == 0 
			&& fabs(GetBody3D(1).XG(4)) < 0.1*drivevel)
		{
			Vector3D dpdx, gv2;
			Vector3D p2 = loccoords(2);
			p2.X() = XG(ns)-xoff;
			GetBody3D(elemind+1).GetdPosdx(p2, dpdx);

			double dpdxn = dpdx.Norm();
			if (dpdxn == 0) dpdxn = 1;
			dpdxn = 1./dpdxn;

			gv2 = dpdx*(dpdxn*XG(nsp));

			GetBody3D(1).XG(4) = gv2.X(); 
			GetBody3D(1).XG(5) = gv2.Y();
			GetBody3D(1).XG(6) = gv2.Z();

			//GetBody3D(1).XG(4) = drivevel; 
			XG(nsp) = drivevel;
		}

		if (/*driven &&*/ (/*XG(ns) >= drivelen*/ GetMBS()->GetTime() >= drivestarttime + drivenextstart) /*&& !driveflag*/)
		{

			XG(l1) = 0;
			XG(l2) = 0;
			XG(l3) = 0;
			XG(ns) = loccoords(2).X();
			XG(nsp) = drivevel;
			//driven = 2;

			drivestarttime += drivenextstart;
			//drivestarttime = 0;

			elemind = 1;
			xoff = 0;

			Vector3D dpdx, gv2;
			//Vector3D gp2;
			//ComputeDpDx(dpdx, gp2, gv2);


			Vector3D p2 = loccoords(2);
			p2.X() = XG(ns)-xoff;

			GetBody3D(elemind+1).GetdPosdx(p2, dpdx);

			double dpdxn = dpdx.Norm();
			if (dpdxn == 0) dpdxn = 1;
			dpdxn = 1./dpdxn;

			gv2 = dpdx*(dpdxn*XG(nsp));

			GetBody3D(1).XG(1) = 0;
			GetBody3D(1).XG(2) = 0;
			GetBody3D(1).XG(3) = 0;

			GetBody3D(1).XG(4) = gv2.X(); 
			GetBody3D(1).XG(5) = gv2.Y();
			GetBody3D(1).XG(6) = gv2.Z();

			driveflag = 0;
		}
	}

	//UO() << "x=" <<	GetBody3D(1).XG(1) << "\n";


	nlstepcnt = 0;
}
int contactz = 0;

double SlidingJoint::PostNewtonStep(double t) 
{
	//return 0;
	nlstepcnt++;
	if (driveflag || nlstepcnt>=2) return 0;
	double chg = 0;
	double epstol = 1e-10;

	if (XG(ns)-xoff > 0.5*(1.+epstol)*GetBody3D(elemind+1).GetSize().X())
	{
		chg = 1;
		//chg = fabs(XG(ns)-xoff - 0.5*GetBody3D(elemind+1).GetSize().X());
		if (elemind+1 < elements.Length()) 
		{
			xoff += 0.5*(GetBody3D(elemind+1).GetSize()).X(); 
			xoff += 0.5*(GetBody3D(elemind+2).GetSize()).X(); 
			elemind++;
		}
		else chg = 0;
		//UO() << "->xoff=" << xoff << ", elemind=" << elemind << "\n";
	} 
	else if (XG(ns)-xoff < -0.5*(1.+epstol)*GetBody3D(elemind+1).GetSize().X())
	{
		//UO() << "-:XG=" << XG(ns) << "\n";
		//UO() << "  elemind=" << elemind << "\n";
		chg = 1;
		//chg=fabs(XG(ns)-xoff + 0.5*GetBody3D(elemind+1).GetSize().X());
		if (elemind+1 > 2) 
		{
			xoff -= 0.5*GetBody3D(elemind+1).GetSize().X(); 
			xoff -= 0.5*GetBody3D(elemind).GetSize().X(); 
			elemind--;
		}
		else chg = 0;
		//UO() << "->xoff=" << xoff << "\n";
	}

	if (contactz)
	{
		Vector3D p2 = loccoords(2);
		p2.X() = XG(ns)-xoff;
		Vector3D d = GetBody3D(1).GetPos(loccoords(1))-GetBody3D(elemind+1).GetPos(p2);

		int iscontactold = iscontact2;
		iscontact2 = iscontact;

		if (d.Z() > -1e-6 && !iscontact) 
		{
			iscontact = 1;
			chg += fabs(d.Z());
		}
		if (XG(l3) > 1e-12 && iscontact)
		{
			iscontact = 0;
			chg += fabs(XG(l3));
		}
	}

	return chg;
}

double lastwrote = 0;
void SlidingJoint::ComputeDpDx(Vector3D& dpdx, Vector3D& gp2, Vector3D& gv2)
{
	Vector3D p2 = loccoords(2);
	p2.X() = XG(ns)-xoff;

	GetBody3D(elemind+1).GetdPosdx(p2, dpdx);

	double dpdxn = dpdx.Norm();
	if (dpdxn == 0) dpdxn = 1;
	dpdxn = 1./dpdxn;

	gp2 = GetBody3D(elemind+1).GetPos(p2);
	gv2 = GetBody3D(elemind+1).GetVel(p2)+dpdx*(dpdxn*XG(nsp));


	int dosmoothing = 0;

	if (dosmoothing)
	{
		//const double smooth = 0.995;
		const double smoothdist = 0.1;

		if (XG(ns)-xoff > 0.5*GetBody3D(elemind+1).GetSize().X()-smoothdist && elemind+1 < elements.Length())
		{
			double fact = 1+0.5*(XG(ns)-xoff - 0.5*GetBody3D(elemind+1).GetSize().X()-smoothdist)/smoothdist;

			Vector3D dpdx2;
			double xoff2;
			Vector3D p22 = loccoords(2);

			xoff2 = xoff;
			xoff2 += 0.5*(GetBody3D(elemind+1).GetSize()).X(); 
			xoff2 += 0.5*(GetBody3D(elemind+2).GetSize()).X(); 

			p22.X() = XG(ns)-xoff2;

			GetBody3D(elemind+1+1).GetdPosdx(p22, dpdx2);
			double dpdxn2 = dpdx2.Norm();
			if (dpdxn2 == 0) dpdxn2 = 1;
			dpdxn2 = 1./dpdxn2;

			dpdx = (1-fact)*dpdx + fact*dpdx2;
			gp2 = (1-fact)*gp2 + fact*GetBody3D(elemind+1+1).GetPos(p22);
			gv2 = (1-fact)*gv2 + fact*(GetBody3D(elemind+1+1).GetVel(p22)+dpdx2*dpdxn2*XG(nsp));
		} 
		else if (XG(ns)-xoff < -0.5*GetBody3D(elemind+1).GetSize().X()+smoothdist && elemind+1 > 2)
		{
			double fact = 0.5*(1-(-(XG(ns)-xoff) - 0.5*GetBody3D(elemind+1).GetSize().X()-smoothdist)/smoothdist)-1;

			Vector3D dpdx2;
			double xoff2;
			Vector3D p22 = loccoords(2);

			xoff2 = xoff;
			xoff2 -= 0.5*GetBody3D(elemind+1).GetSize().X(); 
			xoff2 -= 0.5*GetBody3D(elemind).GetSize().X(); 

			p22.X() = XG(ns)-xoff2;

			GetBody3D(elemind+1-1).GetdPosdx(p22, dpdx2);
			double dpdxn2 = dpdx2.Norm();
			if (dpdxn2 == 0) dpdxn2 = 1;
			dpdxn2 = 1./dpdxn2;

			dpdx = (0.5+fact)*dpdx + (0.5-fact)*dpdx2;
			gp2 = (0.5+fact)*gp2 + (0.5-fact)*GetBody3D(elemind+1-1).GetPos(p22);
			gv2 = (0.5+fact)*gv2 + (0.5-fact)*(GetBody3D(elemind+1-1).GetVel(p22)+dpdx2*dpdxn2*XG(nsp));
		}
	}
}

void SlidingJoint::EvalG(Vector& f, double t) 
{
	// s, lambda1, lambda2, lambda3
	if (MaxIndex()==3)
	{
		/*
		Vector3D p2 = loccoords(2);
		p2.X() = XG(ns)-xoff;
		Vector3D v = GetBody3D(1).GetPos(loccoords(1))-GetBody3D(elemind+1).GetPos(p2);
		Vector3D dpdx;
		GetBody3D(elemind+1).GetdPosdx(p2, dpdx);
		*/

		Vector3D dpdx, gp2, gv2;
		ComputeDpDx(dpdx, gp2, gv2);
		Vector3D v = GetBody3D(1).GetPos(loccoords(1)) - gp2;

		if (!driveflag)
		{
			f(1) = v(1);
			if (!noforcey)
			{
				f(2) = v(2);			
			}
			else
			{
				f(2) = XG(l2);
			}

			if (!contactz)
			{
				f(3) = v(3);
			}
			else
			{
				if (iscontact)
				{
					f(3) = v(3);
				}
				else
				{
					f(3) = XG(l3);
				}
			}
			if (driven)
			{
				double v = drivevel;
				if (driveacceltime != 0) UO() << "ERROR in Slinding Joint: driveacceltime must be zero for index 3 formulation\n";
				if (t < drivestarttime)
				{
					v = 0;
				}
				f(4) = XG(ns) - (t-drivestarttime)*v;
			}
			else
			{
				f(4) = dpdx(1)*XG(l1)+dpdx(2)*XG(l2)+dpdx(3)*XG(l3);
			}
		}
		else
		{
			f(1) = XG(l1);
			f(2) = XG(l2);
			f(3) = XG(l3);
			f(4) = XG(ns)-(drivelen*(1.+1e-14));
		}
		//UO() << "evalfg=" << f << "\n";
	}
	else if (MaxIndex()<=2)
	{

		/*
		Vector3D p2 = loccoords(2);
		p2.X() = XG(ns)-xoff;
		Vector3D dpdx;
		GetBody3D(elemind+1).GetdPosdx(p2, dpdx);
		*/

		if (!driveflag)
		{
			Vector3D dpdx, gp2, gv2;
			ComputeDpDx(dpdx, gp2, gv2);
			Vector3D v = GetBody3D(1).GetVel(loccoords(1))-gv2;

			f(1) = v(1);
			if (!noforcey)
			{
				f(2) = v(2);			
			}
			else
			{
				f(2) = XG(l2);
			}

			if (!contactz)
			{
				f(3) = v(3);
			}
			else
			{
				if (iscontact)
				{
					f(3) = v(3);
				}
				else
				{
					f(3) = XG(l3);
				}
			}

			if (0 && driven && t <= driveendtime)
			{
				double v = drivevel;
				if (t < drivestarttime)
				{
					v = 0;
				} else if (t <= drivestarttime+driveacceltime)
				{
					double p = (t-drivestarttime)/driveacceltime; //from 0 to 1
					v = v*(0.5-0.5*cos(p*MY_PI)); //smooth otherwise lambdas start to oscillate heavily!
					//v = v*p; //linear
				}

				f(4) = XG(nsp) - v;
				//if (t > drivestarttime+driveacceltime) f(4) = dpdx(1)*XG(l1)+dpdx(2)*XG(l2)+dpdx(3)*XG(l3); //not driven any more after acceleration!!!
			}
			/*			else if (driven == 2) //driven mass, fluid, ???
			{
			f(4) = XG(nsp) - drivevel;
			}*/
			else
			{
				f(4) = dpdx(1)*XG(l1)+dpdx(2)*XG(l2)+dpdx(3)*XG(l3);
			}

		}
		else if (driveflag == 1)
		{
			//UO() << "out!\n";
			f(1) = XG(l1);
			f(2) = XG(l2);
			f(3) = XG(l3);
			f(4) = XG(nsp);
		}
	}

};

void SlidingJoint::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	//-C_q^T \lambda = [dC/dq;]^T [\lambda_i] 

	if (driveflag) return;

	if (locelemind == 1)
	{
		GetBody3D(1).GetdPosdqT(loccoords(1),dpdq);
		//UO() << "dpdq1=" << dpdq << "\n";

		for (int i=1; i <= f.Length(); i++)
		{
			f(i) -= (dpdq(i,1)*XG(l1)+dpdq(i,2)*XG(l2)+dpdq(i,3)*XG(l3));
		}
	}
	else if (locelemind == elemind+1)
	{
		Vector3D p = loccoords(elemind+1);
		p.X() = XG(ns)-xoff;
		GetBody3D(elemind+1).GetdPosdqT(p,dpdq);
		//UO() << "dpdq2=" << dpdq << "\n";

		for (int i=1; i <= f.Length(); i++)
		{
			f(i) -= -(dpdq(i,1)*XG(l1)+dpdq(i,2)*XG(l2)+dpdq(i,3)*XG(l3));
		}
	}
};

Vector3D SlidingJoint::GetRefPosD()	const 
{
	return GetBody3D(1).GetPosD(loccoords(1));
}

Vector3D SlidingJoint::GetRefPos()	const 
{
	return GetBody3D(1).GetPos(loccoords(1));
}

Vector3D SlidingJoint::GetRefVelD()	const 
{
	return GetBody3D(1).GetVelD(loccoords(1));
}

Vector3D SlidingJoint::GetRefVel()	const 
{
	return GetBody3D(1).GetVel(loccoords(1));
}

double SlidingJoint::GetDrift() const
{
	Vector3D p2 = loccoords(2);
	p2.X() = XG(ns)-xoff;
	Vector3D v = GetBody3D(1).GetPos(loccoords(1))-GetBody3D(elemind+1).GetPos(p2);
	return v.Norm();
}

double SlidingJoint::GetDriftP()
{
	Vector3D p2 = loccoords(2);
	p2.X() = XG(ns)-xoff;
	Vector3D dpdx;
	GetBody3D(elemind+1).GetdPosdx(p2, dpdx);

	double dpdxn = dpdx.Norm();
	if (dpdxn == 0) dpdxn = 1;
	dpdxn = 1./dpdxn;

	Vector3D v = GetBody3D(1).GetVel(loccoords(1))-(GetBody3D(elemind+1).GetVel(p2)+dpdx*dpdxn*XG(nsp));
	return v.Norm();
}

Vector3D SlidingJoint::GetSlidingPos()	const 
{
	Vector3D p2 = loccoords(2);
	p2.X() = XG(ns)-xoff;
	return GetBody3D(elemind+1).GetPos(p2);
}

void SlidingJoint::DrawElement() 
{
	Constraint::DrawElement();
	//

	Vector3D lc1=loccoords(1);
	Vector3D lc2=loccoords(2);

	mbs->SetColor(GetCol());
	Vector3D p1 = GetBody3D(1).GetPosD(lc1);
	mbs->DrawSphere(p1,0.5*draw_dim.X(),16);
	double x = XGD(ns);
	int n = 1;
	while (x > GetBody3D(n+1).GetSize().X()*0.5 && n+1 < elements.Length()) 
	{
		x -= 0.5*GetBody3D(n+1).GetSize().X(); 
		x -= 0.5*GetBody3D(n+2).GetSize().X(); 
		n++;
	}
	Vector3D p2 = GetBody3D(n+1).GetPosD(Vector3D(x,lc2.Y(),lc2.Z()));

	if (1)
	{
		mbs->SetColor(colred);
		mbs->DrawSphere(p2,0.5*draw_dim.X(),12);
	}
	//draw lambda-vector:
	//Vector3D v(XGD(l1),XGD(l2),XGD(l3));
	//v.Normalize();
	//mbs->DrawZyl(p2,0.2*v+p2,0.01,8);

	if(0)
	{
		char str[100];
		sprintf(str,"    drift = %18.15g   ", (p1-p2).Norm());
		GetMBS()->GetRC()->PrintTextStruct(7,-1,str);
	}
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// SlidingPointJoint
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SlidingPointJoint::SetSlidingPointJoint(int en1, int en2, 
													 const Vector3D& lc1, const Vector3D& lc2, double cdim, const Vector3D& coli)
{	
	GetCol() = coli;
	draw_dim.X() = cdim;
	elements(1) = en1;
	elements(2) = en2;
	loccoords(1) = lc1;
	loccoords(2) = lc2;
	elemind = 1; //sliding element, 1 means elements(2)
	xoff = 0;
};

void SlidingPointJoint::SetSlidingPointJoint(int en1, const TArray<int>& en2, int in2, 
													 const Vector3D& lc1, const Vector3D& lc2, double cdim, 
													 const Vector3D& coli)
{	
	GetCol() = coli;
	draw_dim.X() = cdim;
	elements.SetLen(1+en2.Length());
	elements(1) = en1;
	for (int i=1; i <= en2.Length(); i++)
	{
		elements(i+1) = en2(i);
	}
	loccoords(1) = lc1;
	loccoords(2) = lc2;
	elemind = in2;
	xoff = 0;
};

void SlidingPointJoint::SetGlobalInitConditions(Vector& x_glob)
{
	// set s (x_init(1)) as the x component of lc2 and the other entries are already zero
	x_init(1) = loccoords(2).X();
	for (int i=1; i<=SS(); i++)
	{
		x_glob(ltg.Get(i)) = x_init(i);
	}
}

void SlidingPointJoint::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Element::GetElementData(edc);
}

int SlidingPointJoint::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	return Element::SetElementData(edc);
}

int SlidingPointJoint::GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables)
{
	return 0;
}

int SlidingPointJoint::ReadSingleElementData(ReadWriteElementDataVariableType& RWdata) 		
{
	// call base class routine 	
	int rv = Constraint::ReadSingleElementData(RWdata);
	if (rv == 1) return 1;

	return ReadSingleElementDataAuto(RWdata);
}

void SlidingPointJoint::ElementDefaultConstructorInitialization()
{	
	elemind = 1;
	loccoords.SetLen(2);
	elements.SetLen(2);
	elements(1) = 1;
	elements(2) = 2;
	elementname = GetElementSpec();
	draw_dim.X() = 0.01;
	xoff = 0;
	x_init = Vector(SS()); // s, lambda1, lambda2, lambda3
}

void SlidingPointJoint::CopyFrom(const Element& e)
{
	Constraint::CopyFrom(e);
	const SlidingPointJoint& ce = (const SlidingPointJoint&)e;
	loccoords = ce.loccoords;
	dpdq = ce.dpdq;
	elemind = ce.elemind;
	xoff = ce.xoff;
}

void SlidingPointJoint::Initialize() 
{
	int n = 1;
	xoff = 0;
	while (n < elemind) 
	{
		xoff += 0.5*GetBody3D(n+1).GetSize().X(); 
		xoff += 0.5*GetBody3D(n+2).GetSize().X(); 
		n++;
	}
};

void SlidingPointJoint::GetNecessaryKinematicAccessFunctions(TArray<int> &KinAccFunc, int numberOfKinematicPair)
{
	KinAccFunc.SetLen(0);
	int kaf = (int)(TKinematicsAccessFunctions(TKAF_none));

	kaf = (int)(TKinematicsAccessFunctions(TKAF_position+TKAF_velocity+TKAF_D_pos_D_x+TKAF_D_pos_D_q));	
	KinAccFunc.Add(kaf);
}
void SlidingPointJoint::EvalF(Vector& f, double t) 
{
	Element::EvalF(f,t);
	f(ns) = XG(nsp); // d(s)/dt = sp
}; 

double SlidingPointJoint::PostNewtonStep(double t) 
{
	// chg is returned: 1 if index and xoff changed, else 0
	// xoff and elemind are updated if necessary

	double chg = 0;
	double epstol = 1e-10;

	if (XG(ns)-xoff > 0.5*(1.+epstol)*GetBody3D(elemind+1).GetSize().X()) // if sliding point changes the index of en2 elements ==> correct xoff and elemind
	{
		chg = 1;
		if (elemind+1 < elements.Length()) 
		{
			xoff += 0.5*(GetBody3D(elemind+1).GetSize()).X(); 
			xoff += 0.5*(GetBody3D(elemind+2).GetSize()).X(); 
			elemind++;
		}
		else chg = 0;
	} 
	else if (XG(ns)-xoff < -0.5*(1.+epstol)*GetBody3D(elemind+1).GetSize().X())
	{
		chg = 1;
		if (elemind+1 > 2) 
		{
			xoff -= 0.5*GetBody3D(elemind+1).GetSize().X(); 
			xoff -= 0.5*GetBody3D(elemind).GetSize().X(); 
			elemind--;
		}
		else chg = 0;
	}

	return chg;
}


void SlidingPointJoint::ComputeDpDx(Vector3D& dpdx, Vector3D& gp2, Vector3D& gv2) 
{
	Vector3D p2 = loccoords(2);
	p2.X() = XG(ns)-xoff; // XG(ns) ==> sliding parameter s
	// p2 is now the vector from the cog from the actual element en2_x to the sliding point

	GetBody3D(elemind+1).GetdPosdx(p2, dpdx); // +1 because body 1 comes first

	double dpdxn = dpdx.Norm();
	if (dpdxn == 0) dpdxn = 1;
	dpdxn = 1./dpdxn;

	gp2 = GetBody3D(elemind+1).GetPos(p2);
	gv2 = GetBody3D(elemind+1).GetVel(p2)+dpdx*(dpdxn*XG(nsp));
	// formulas see SlidingJoints.pdf document, eq. (5)
}

void SlidingPointJoint::EvalG(Vector& f, double t) 
{
	// s, lambda1, lambda2, lambda3
	if (MaxIndex()==3)
	{

		Vector3D dpdx, gp2, gv2;
		ComputeDpDx(dpdx, gp2, gv2);
		Vector3D p = GetBody3D(1).GetPos(loccoords(1)) - gp2;

		f(1) = p(1);
		f(2) = p(2);
		f(3) = p(3);
		f(4) = dpdx(1)*XG(l1)+dpdx(2)*XG(l2)+dpdx(3)*XG(l3); // force acts normal to the sliding direction, see eq. (5) in SlidingJoints.pdf document
	}
	else if (MaxIndex()<=2)
	{
		Vector3D dpdx, gp2, gv2;
		ComputeDpDx(dpdx, gp2, gv2);
		Vector3D v = GetBody3D(1).GetVel(loccoords(1))-gv2;

		f(1) = v(1);
		f(2) = v(2);
		f(3) = v(3);
		f(4) = dpdx(1)*XG(l1)+dpdx(2)*XG(l2)+dpdx(3)*XG(l3);
	}
};

void SlidingPointJoint::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	//-C_q^T \lambda = [dC/dq;]^T [\lambda_i] 

	if (locelemind == 1)
	{
		GetBody3D(1).GetdPosdqT(loccoords(1),dpdq);

		for (int i=1; i <= f.Length(); i++)
		{
			f(i) -= (dpdq(i,1)*XG(l1)+dpdq(i,2)*XG(l2)+dpdq(i,3)*XG(l3));
		}
	}
	else if (locelemind == elemind+1)
	{
		Vector3D p = loccoords(2);
		p.X() = XG(ns)-xoff;
		GetBody3D(elemind+1).GetdPosdqT(p,dpdq);

		for (int i=1; i <= f.Length(); i++)
		{
			f(i) -= -(dpdq(i,1)*XG(l1)+dpdq(i,2)*XG(l2)+dpdq(i,3)*XG(l3));
		}
	}
};

Vector3D SlidingPointJoint::GetRefPosD()	const 
{
	return GetBody3D(1).GetPosD(loccoords(1));
}

Vector3D SlidingPointJoint::GetRefPos()	const 
{
	return GetBody3D(1).GetPos(loccoords(1));
}

Vector3D SlidingPointJoint::GetRefVelD()	const 
{
	return GetBody3D(1).GetVelD(loccoords(1));
}

Vector3D SlidingPointJoint::GetRefVel()	const 
{
	return GetBody3D(1).GetVel(loccoords(1));
}

double SlidingPointJoint::GetDrift() const
{
	// returns the vector length (vector is the difference in sliding point positions)
	Vector3D p2 = loccoords(2);
	p2.X() = XG(ns)-xoff;
	Vector3D v = GetBody3D(1).GetPos(loccoords(1))-GetBody3D(elemind+1).GetPos(p2);
	return v.Norm();
}

double SlidingPointJoint::GetDriftP()
{
	// returns the drift velocity between the sliding points
	Vector3D p2 = loccoords(2);
	p2.X() = XG(ns)-xoff;
	Vector3D dpdx;
	GetBody3D(elemind+1).GetdPosdx(p2, dpdx);

	double dpdxn = dpdx.Norm();
	if (dpdxn == 0) dpdxn = 1;
	dpdxn = 1./dpdxn;

	Vector3D v = GetBody3D(1).GetVel(loccoords(1))-(GetBody3D(elemind+1).GetVel(p2)+dpdx*dpdxn*XG(nsp));
	return v.Norm();
}

Vector3D SlidingPointJoint::GetSlidingPos()	const 
{
	// returns the global vector to the sliding point on element 2
	Vector3D p2 = loccoords(2);
	p2.X() = XG(ns)-xoff;
	return GetBody3D(elemind+1).GetPos(p2);
}

void SlidingPointJoint::DrawElement() 
{
	Constraint::DrawElement();

	Vector3D lc1=loccoords(1);
	Vector3D lc2=loccoords(2);

	mbs->SetColor(GetCol());
	Vector3D p1 = GetBody3D(1).GetPosD(lc1);
	mbs->DrawSphere(p1,0.5*draw_dim.X(),16);
	double x = XGD(ns);
	int n = 1;
	while (x > GetBody3D(n+1).GetSize().X()*0.5 && n+1 < elements.Length()) 
	{
		x -= 0.5*GetBody3D(n+1).GetSize().X(); 
		x -= 0.5*GetBody3D(n+2).GetSize().X(); 
		n++;
	}
	Vector3D p2 = GetBody3D(n+1).GetPosD(Vector3D(x,lc2.Y(),lc2.Z()));

	if (1)
	{
		mbs->SetColor(colred);
		mbs->DrawSphere(p2,0.5*draw_dim.X(),12);
	}
	if(0)
	{
		char str[100];
		sprintf(str,"    drift = %18.15g   ", (p1-p2).Norm());
		GetMBS()->GetRC()->PrintTextStruct(7,-1,str);
	}
};


void SlidingPrismaticJoint::SetSlidingPrismaticJoint(int en1, int en2, 
													 const Vector3D& lc1, const Vector3D& lc2, double cdim, const Vector3D& coli, const Vector3D& stiffness, const Vector3D& damping)
{	
	GetCol() = coli;
	draw_dim.X() = cdim;
	elements(1) = en1;
	elements(2) = en2;
	loccoords(1) = lc1;
	loccoords(2) = lc2;
	elemind = 1; //sliding element, 1 means elements(2)
	xoff = 0;
	k1 = stiffness(1);
	k2 = stiffness(2);
	k3 = stiffness(3);
	d1 = damping(1);
	d2 = damping(2);
	d3 = damping(3);
};

void SlidingPrismaticJoint::SetSlidingPrismaticJoint(int en1, const TArray<int>& en2, int in2, 
													 const Vector3D& lc1, const Vector3D& lc2, double cdim, 
													 const Vector3D& coli, const Vector3D& stiffness, const Vector3D& damping)
{	
	GetCol() = coli;
	draw_dim.X() = cdim;
	elements.SetLen(1+en2.Length());
	elements(1) = en1;
	for (int i=1; i <= en2.Length(); i++)
	{
		elements(i+1) = en2(i);
	}
	loccoords(1) = lc1;
	loccoords(2) = lc2;
	elemind = in2;
	xoff = 0;
	k1 = stiffness(1);
	k2 = stiffness(2);
	k3 = stiffness(3);
	d1 = damping(1);
	d2 = damping(2);
	d3 = damping(3);
};

void SlidingPrismaticJoint::LinkToElements() 
{
	for (int i = 1; i <= NE(); i++) 
	{
		GetElem(i).AddConstraint(this,i);
	}
}

void SlidingPrismaticJoint::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Element::GetElementData(edc);
}

int SlidingPrismaticJoint::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	return Element::SetElementData(edc);
}

int SlidingPrismaticJoint::GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables)
{
	Constraint::GetAvailableSpecialValues(available_variables);

	// Automatic entries for this class 
	Constraint::GetAvailableSpecialValuesAuto(available_variables);

	// Manual READ entries for this class
	available_variables.Add(ReadWriteElementDataVariableType("Connector.sliding_parameter",0,0,0.,mystr("internal sliding parameter s"), TRWElementDataRead));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.sliding_parameter_p",0,0,0.,mystr("internal time derivative of sliding parameter s"), TRWElementDataRead));

	return 0;
}

int SlidingPrismaticJoint::ReadSingleElementData(ReadWriteElementDataVariableType& RWdata) 		
{
	// call base class routine 	
	int rv = Constraint::ReadSingleElementData(RWdata);
	if (rv == 1) return 1;

	//manual things to read  
	if(RWdata.variable_name.CStrCompare("Connector.sliding_parameter"))
	{
		RWdata.value = XG(ns);
		return 1;
	}
	if(RWdata.variable_name.CStrCompare("Connector.sliding_parameter_p"))
	{
		RWdata.value = XG(nsp);
		return 1;
	}
	return ReadSingleElementDataAuto(RWdata);
}

void SlidingPrismaticJoint::ElementDefaultConstructorInitialization()
{	
	elemind = 1;
	loccoords.SetLen(8);
	elements.SetLen(2);
	elements(1) = 1;
	elements(2) = 2;
	use_penalty_formulation = 1;
	k1 = k2 = k3 = 1e5;
	d1 = d2 = d3 = 100;
	elementname = GetElementSpec();
	draw_dim.X() = 0.01;
	xoff = 0;
	x_init = Vector(SS()); // s, lambda1, lambda2, lambda3
}

void SlidingPrismaticJoint::CopyFrom(const Element& e)
{
	Constraint::CopyFrom(e);
	const SlidingPrismaticJoint& ce = (const SlidingPrismaticJoint&)e;
	loccoords = ce.loccoords;
	dpdq = ce.dpdq;
	drotdq = ce.drotdq;
	elemind = ce.elemind;
	xoff = ce.xoff;
	ind = ce.ind;
	k1 = ce.k1;
	k2 = ce.k2;
	k3 = ce.k3;
	d1 = ce.d1;
	d2 = ce.d2;
	d3 = ce.d3;
}

void SlidingPrismaticJoint::Initialize() 
{
	ind = 3;
	if (k1 == 0 && k2 == 0 && k3 == 0) 
	{
		SetPenaltyFormulation(0);
	} else {
		SetPenaltyFormulation(1);
	}

	int n = 1;
	xoff = 0;
	while (n < elemind) 
	{
		xoff += 0.5*GetBody3D(n+1).GetSize().X(); 
		xoff += 0.5*GetBody3D(n+2).GetSize().X(); 
		n++;
	}
	//find orthogonal vectors:
	const Vector3D& lp1 = loccoords(1);
	const Vector3D& lp2 = loccoords(2);
	Vector3D vi1,vi2,vi3;
	Vector3D vj1,vj2,vj3;
	
	vi1 = Vector3D(1.,0.,0.);
	vj1 = vi1;
	vi2 = Vector3D(0.,1.,0.);
	vj2 = vi2;
	vi3 = Vector3D(0.,0.,1.);
	vj3 = vi3;

	Matrix3D RTi = GetBody3D(1).GetRotMatrix(lp1).GetTp();
	Matrix3D RTj = GetBody3D(elemind+1).GetRotMatrix(lp2).GetTp();
	vi1 = RTi*vi1;
	vi2 = RTi*vi2;
	vi3 = RTi*vi3;
	vj1 = RTj*vj1;
	vj2 = RTj*vj2;
	vj3 = RTj*vj3;

	loccoords(ind) = vi1; //body i
	loccoords(ind+1) = vi2;	//body i
	loccoords(ind+2) = vi3;	//body i
	loccoords(ind+3) = vj1; //body j
	loccoords(ind+4) = vj2; //body j
	loccoords(ind+5) = vj3; //body j
};

void SlidingPrismaticJoint::GetNecessaryKinematicAccessFunctions(TArray<int> &KinAccFunc, int numberOfKinematicPair)
{
	KinAccFunc.SetLen(0);
	int kaf = (int)(TKinematicsAccessFunctions(TKAF_none));

	kaf = (int)(TKinematicsAccessFunctions(TKAF_position+TKAF_velocity+TKAF_D_pos_D_x+TKAF_D_pos_D_q+TKAF_rotation_matrix+TKAF_rotation_matrix_P+TKAF_D_rot_D_q));	
	KinAccFunc.Add(kaf);
}

void SlidingPrismaticJoint::EvalF(Vector& f, double t) 
{
	Element::EvalF(f,t);
	f(ns) = XG(nsp);
}; 

double SlidingPrismaticJoint::PostNewtonStep(double t) 
{
	// chg is returned: 1 if index and xoff changed, else 0
	// xoff and elemind are updated if necessary

	double chg = 0;
	double epstol = 1e-10;

	if (XG(ns)-xoff > 0.5*(1.+epstol)*GetBody3D(elemind+1).GetSize().X())
	{
		chg = 1;
		if (elemind+1 < elements.Length()) 
		{
			xoff += 0.5*(GetBody3D(elemind+1).GetSize()).X(); 
			xoff += 0.5*(GetBody3D(elemind+2).GetSize()).X(); 
			elemind++;
		}
		else chg = 0;
	} 
	else if (XG(ns)-xoff < -0.5*(1.+epstol)*GetBody3D(elemind+1).GetSize().X())
	{
		chg = 1;
		if (elemind+1 > 2) 
		{
			xoff -= 0.5*GetBody3D(elemind+1).GetSize().X(); 
			xoff -= 0.5*GetBody3D(elemind).GetSize().X(); 
			elemind--;
		}
		else chg = 0;
	}
	return chg;
}


void SlidingPrismaticJoint::ComputeDpDx(Vector3D& dpdx, Vector3D& gp2, Vector3D& gv2)
{
	Vector3D p2 = loccoords(2);
	p2.X() = XG(ns)-xoff;

	GetBody3D(elemind+1).GetdPosdx(p2, dpdx);

	double dpdxn = dpdx.Norm();
	if (dpdxn == 0) dpdxn = 1;
	dpdxn = 1./dpdxn;

	gp2 = GetBody3D(elemind+1).GetPos(p2);
	gv2 = GetBody3D(elemind+1).GetVel(p2)+dpdx*(dpdxn*XG(nsp));
}

void SlidingPrismaticJoint::EvalG(Vector& f, double t) 
{
	// s, lambda1, lambda2, lambda3
	if (MaxIndex()==3)
	{
		Vector3D dpdx, gp2, gv2;
		ComputeDpDx(dpdx, gp2, gv2);
		Vector3D v = GetBody3D(1).GetPos(loccoords(1)) - gp2;

		f(1) = v(1);
		f(2) = v(2);
		f(3) = v(3);
		f(4) = dpdx(1)*XG(l1)+dpdx(2)*XG(l2)+dpdx(3)*XG(l3);

		Vector3D p1 = loccoords(1);
		Vector3D p2 = loccoords(2);
		p2.X() = XG(ns)-xoff;
		// p2 is now the vector from the cog from the actual element en2_x to the sliding point

		Matrix3D Roti = GetBody3D(1).GetRotMatrix(p1); // rot matrix at local position
		Matrix3D Rotj = GetBody3D(elemind+1).GetRotMatrix(p2); // rot matrix at local position
		Matrix3D Rotip= GetBody3D(1).GetRotMatrixP(p1);
		Matrix3D Rotjp= GetBody3D(elemind+1).GetRotMatrixP(p2);


		Vector3D vi1g,vi2g,vi3g,vi1pg,vi2pg,vi3pg;
		Vector3D vj1g,vj2g,vj3g,vj1pg,vj2pg,vj3pg;

		const Vector3D& vi1 = loccoords(ind);
		const Vector3D& vi2 = loccoords(ind+1);
		const Vector3D& vi3 = loccoords(ind+2);
		const Vector3D& vj1 = loccoords(ind+3);
		const Vector3D& vj2 = loccoords(ind+4);
		const Vector3D& vj3 = loccoords(ind+5);

		vi1g = Roti*vi1;
		vi2g = Roti*vi2;
		vi3g = Roti*vi3;
		vj1g = Rotj*vj1;
		vj2g = Rotj*vj2;
		vj3g = Rotj*vj3;
		vi1pg = Rotip*vi1;
		vi2pg = Rotip*vi2;
		vi3pg = Rotip*vi3;
		vj1pg = Rotjp*vj1;
		vj2pg = Rotjp*vj2;
		vj3pg = Rotjp*vj3;


		// formulas see SlidingJoints.pdf document
		if (!UsePenaltyFormulation())
		{
			f(5) = vj2g*vi3g;
			f(6) = vj3g*vi1g;
			f(7) = vj2g*vi1g;
		} else // these are global conditions / rotations
		{
			double vj2gvi3g = vj2g*vi3g;
			double vj3gvi1g = vj3g*vi1g;
			double vj2gvi1g = vj2g*vi1g;

			if(vj2gvi3g > 1) {
				vj2gvi3g = 1;
			} else if(vj2gvi3g < -1) {
				vj2gvi3g = -1;
			}
			if(vj3gvi1g > 1) {
				vj3gvi1g = 1;
			} else if(vj3gvi1g < -1) {
				vj3gvi1g = -1;
			}
			if(vj2gvi1g > 1) {
				vj2gvi1g = 1;
			} else if(vj2gvi1g < -1) {
				vj2gvi1g = -1;
			}

			f(5) = asin(vj2gvi3g)*k1+1/sqrt(1-Sqr(fabs(vj2gvi3g)))*(vj2g*vi3pg + vj2pg*vi3g)*d1+XG(l4);
			f(6) = asin(vj3gvi1g)*k2+1/sqrt(1-Sqr(fabs(vj3gvi1g)))*(vj3g*vi1pg + vj3pg*vi1g)*d2+XG(l5);
			f(7) = -asin(vj2gvi1g)*k3-1/sqrt(1-Sqr(fabs(vj2gvi1g)))*(vj2g*vi1pg + vj2pg*vi1g)*d3+XG(l6);
		}

	}
	else if (MaxIndex()<=2)
	{
		Vector3D dpdx, gp2, gv2;
		ComputeDpDx(dpdx, gp2, gv2);
		Vector3D v = GetBody3D(1).GetVel(loccoords(1))-gv2;

		f(1) = v(1);
		f(2) = v(2);
		f(3) = v(3);
		f(4) = dpdx(1)*XG(l1)+dpdx(2)*XG(l2)+dpdx(3)*XG(l3);

		Vector3D p1 = loccoords(1);
		Vector3D p2 = loccoords(2);
		p2.X() = XG(ns)-xoff;

		Matrix3D Roti = GetBody3D(1).GetRotMatrix(p1); // rot matrix at local position
		Matrix3D Rotj = GetBody3D(elemind+1).GetRotMatrix(p2); // rot matrix at local position
		Matrix3D Rotip= GetBody3D(1).GetRotMatrixP(p1);
		Matrix3D Rotjp= GetBody3D(elemind+1).GetRotMatrixP(p2);


		Vector3D vi1g,vi2g,vi3g,vi1pg,vi2pg,vi3pg;
		Vector3D vj1g,vj2g,vj3g,vj1pg,vj2pg,vj3pg;

		const Vector3D& vi1 = loccoords(ind);
		const Vector3D& vi2 = loccoords(ind+1);
		const Vector3D& vi3 = loccoords(ind+2);
		const Vector3D& vj1 = loccoords(ind+3);
		const Vector3D& vj2 = loccoords(ind+4);
		const Vector3D& vj3 = loccoords(ind+5);

		vi1g = Roti*vi1;
		vi2g = Roti*vi2;
		vi3g = Roti*vi3;
		vj1g = Rotj*vj1;
		vj2g = Rotj*vj2;
		vj3g = Rotj*vj3;
		vi1pg = Rotip*vi1;
		vi2pg = Rotip*vi2;
		vi3pg = Rotip*vi3;
		vj1pg = Rotjp*vj1;
		vj2pg = Rotjp*vj2;
		vj3pg = Rotjp*vj3;

		if (!UsePenaltyFormulation())
		{
			f(5) = vj2g*vi3pg + vj2pg*vi3g;
			f(6) = vj3g*vi1pg + vj3pg*vi1g;
			f(7) = vj2g*vi1pg + vj2pg*vi1g;
		} else
		{
			double vj2gvi3g = vj2g*vi3g;
			double vj3gvi1g = vj3g*vi1g;
			double vj2gvi1g = vj2g*vi1g;

			if(vj2gvi3g > 1) {
				vj2gvi3g = 1;
			} else if(vj2gvi3g < -1) {
				vj2gvi3g = -1;
			}
			if(vj3gvi1g > 1) {
				vj3gvi1g = 1;
			} else if(vj3gvi1g < -1) {
				vj3gvi1g = -1;
			}
			if(vj2gvi1g > 1) {
				vj2gvi1g = 1;
			} else if(vj2gvi1g < -1) {
				vj2gvi1g = -1;
			}

			f(5) = asin(vj2gvi3g)*k1+1/sqrt(1-Sqr(fabs(vj2gvi3g)))*(vj2g*vi3pg + vj2pg*vi3g)*d1+XG(l4);
			f(6) = asin(vj3gvi1g)*k2+1/sqrt(1-Sqr(fabs(vj3gvi1g)))*(vj3g*vi1pg + vj3pg*vi1g)*d2+XG(l5);
			f(7) = -asin(vj2gvi1g)*k3-1/sqrt(1-Sqr(fabs(vj2gvi1g)))*(vj2g*vi1pg + vj2pg*vi1g)*d3+XG(l6);
		}
	}
};

void SlidingPrismaticJoint::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	//-C_q^T \lambda = [dC/dq;]^T [\lambda_i] 

	if (locelemind == 1)
	{
		GetBody3D(1).GetdPosdqT(loccoords(1),dpdq);
		GetBody3D(1).GetdRotdqT(loccoords(1),drotdq);
		//UO() << "dpdq1=" << dpdq << "\n";

		for (int i=1; i <= f.Length(); i++)
		{
			f(i) -= (dpdq(i,1)*XG(l1)+dpdq(i,2)*XG(l2)+dpdq(i,3)*XG(l3)+drotdq(i,1)*XG(l4)+drotdq(i,2)*XG(l5)+drotdq(i,3)*XG(l6));
		}

	}
	else if (locelemind == elemind+1)
	{
		Vector3D p = loccoords(2);
		p.X() = XG(ns)-xoff;
		GetBody3D(elemind+1).GetdPosdqT(p,dpdq);
		GetBody3D(elemind+1).GetdRotdqT(p,drotdq);
		//UO() << "dpdq2=" << dpdq << "\n";

		for (int i=1; i <= f.Length(); i++)
		{
			f(i) -= -(dpdq(i,1)*XG(l1)+dpdq(i,2)*XG(l2)+dpdq(i,3)*XG(l3)+drotdq(i,1)*XG(l4)+drotdq(i,2)*XG(l5)+drotdq(i,3)*XG(l6));
		}
	}
	//mbs->UO() << "Fx:" << XG(l1) << "\n";
	//mbs->UO() << "Fy:" << XG(l2) << "\n";
	//mbs->UO() << "Fz:" << XG(l3) << "\n";
	//mbs->UO() << "Mx:" << XG(l4) << "\n";
	//mbs->UO() << "My:" << XG(l5) << "\n";
	//mbs->UO() << "Mz:" << XG(l6) << "\n";
};

Vector3D SlidingPrismaticJoint::GetRefPosD()	const 
{
	return GetBody3D(1).GetPosD(loccoords(1));
}

Vector3D SlidingPrismaticJoint::GetRefPos()	const 
{
	return GetBody3D(1).GetPos(loccoords(1));
}

Vector3D SlidingPrismaticJoint::GetRefVelD()	const 
{
	return GetBody3D(1).GetVelD(loccoords(1));
}

Vector3D SlidingPrismaticJoint::GetRefVel()	const 
{
	return GetBody3D(1).GetVel(loccoords(1));
}

double SlidingPrismaticJoint::GetDrift() const
{
	Vector3D p2 = loccoords(2);
	p2.X() = XG(ns)-xoff;
	Vector3D v = GetBody3D(1).GetPos(loccoords(1))-GetBody3D(elemind+1).GetPos(p2);
	return v.Norm();
}

double SlidingPrismaticJoint::GetDriftP()
{
	Vector3D p2 = loccoords(2);
	p2.X() = XG(ns)-xoff;
	Vector3D dpdx;
	GetBody3D(elemind+1).GetdPosdx(p2, dpdx);

	double dpdxn = dpdx.Norm();
	if (dpdxn == 0) dpdxn = 1;
	dpdxn = 1./dpdxn;

	Vector3D v = GetBody3D(1).GetVel(loccoords(1))-(GetBody3D(elemind+1).GetVel(p2)+dpdx*dpdxn*XG(nsp));
	return v.Norm();
}

Vector3D SlidingPrismaticJoint::GetSlidingPos()	const 
{
	Vector3D p2 = loccoords(2);
	p2.X() = XG(ns)-xoff;
	return GetBody3D(elemind+1).GetPos(p2);
}

void SlidingPrismaticJoint::DrawElement() 
{
	Constraint::DrawElement();

	Vector3D lc1=loccoords(1);
	Vector3D lc2=loccoords(2);

	mbs->SetColor(GetCol());
	Vector3D p1 = GetBody3D(1).GetPosD(lc1);
	mbs->DrawSphere(p1,0.5*draw_dim.X(),16);
	double x = XGD(ns);
	int n = 1;
	while (x > GetBody3D(n+1).GetSize().X()*0.5 && n+1 < elements.Length()) 
	{
		x -= 0.5*GetBody3D(n+1).GetSize().X(); 
		x -= 0.5*GetBody3D(n+2).GetSize().X(); 
		n++;
	}

	Vector3D p2 = GetBody3D(n+1).GetPosD(Vector3D(x,lc2.Y(),lc2.Z()));

	if (1)
	{
		mbs->SetColor(colred);
		mbs->DrawSphere(p2,0.5*draw_dim.X(),12);
	}

	if(0)
	{
		char str[100];
		sprintf(str,"    drift = %18.15g   ", (p1-p2).Norm());
		GetMBS()->GetRC()->PrintTextStruct(7,-1,str);
	}
};





//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// SlidingJoint2D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SlidingJoint2D::SlidingJoint2D(MBS* mbsi, int en1, const TArray<int>& en2, int in2, 
															 const Vector2D& lc1, const Vector2D& lc2, double cdim, 
															 const Vector3D& coli, int constrainrotationI):Constraint(mbsi), loccoords(), dpdq() 
{	
	constrainrotation = constrainrotationI;

	x_init = Vector(SS()); // s, lambda1, lambda2
	x_init(1) = lc2.X(); //initial value for sliding parameter
	
	GetCol() = coli;
	draw_dim.X() = cdim;
	AddElementCoord(en1, lc1);
	AddElementCoord(en2(1), lc2);
	for (int i=2; i <= en2.Length(); i++)
	{
		AddElementCoord(en2(i), lc2);
	}
	elemind = in2;
	ind_init = elemind;
	xoff = 0;
	double x = lc2.X();
	int n = 1;
	while (n < in2) 
	{
		xoff += 0.5*GetBody2D(n+1).GetSize().X(); 
		xoff += 0.5*GetBody2D(n+2).GetSize().X(); 
		n++;
	}
	x_init(1) = xoff+lc2.X(); //initial value for sliding parameter

	//UO() << "elemind=" << elemind << ", xoff=" << xoff << "\n";
	nlstepcnt = 0;
	driven = 0;
	driveendtime = 1e10;
	drivenextstart = 0;
	driveflag = 0;
	drivelen = 1e30;

	elementname = GetElementSpec();
};

void SlidingJoint2D::EvalF(Vector& f, double t) 
{
	Element::EvalF(f,t);
	f(ns) = XG(nsp);
}; 

void SlidingJoint2D::PostprocessingStep() 
{
	if (!driven) return;

	if (drivenextstart == 0)
	{
		if (XG(ns) >= drivelen && !driveflag)
		{
			driveflag = 1;
		}
	}
	else 
	{
		if (XG(ns) >= drivelen && !driveflag)
		{
			driveflag = 1;
		}

		//very special for fluid problem!!!!!! WARNING ??????
		//UO() << "drivestarttime=" << drivestarttime << "\n";
		if (0 && GetMBS()->GetTime()+1*GetMBS()->GetStepSize() >= drivestarttime && driveacceltime == 0 
			&& fabs(GetBody2D(1).XGP(1)) < 0.01*drivevel)
		{
			Vector2D dpdx, gv2;
			Vector2D p2 = loccoords(2);
			p2.X() = XG(ns)-xoff;
			GetBody2D(elemind+1).GetdPosdx(p2, dpdx);

			double dpdxn = dpdx.Norm();
			if (dpdxn == 0) dpdxn = 1;
			dpdxn = 1./dpdxn;

			gv2 = dpdx*(dpdxn*XG(nsp));

			//UO() << "v0=" << gv2.X() << "," << gv2.Y() << "\n";
			//UO() << "xgp=" << GetBody2D(1).XGP(1) << ", drivevel=" << drivevel << ", dst=" << drivestarttime << "\n";

			GetBody2D(1).XGP(1) = gv2.X(); 
			GetBody2D(1).XGP(2) = gv2.Y();
			XG(nsp) = drivevel;
		}

		if ((GetMBS()->GetTime() >= drivestarttime + drivenextstart) && driven)
		{

			XG(l1) = 0;
			XG(l2) = 0;
			XG(ns) = loccoords(2).X(); //0.5*drivevel*GetMBS()->GetStepSize(); //correction for trapezoidal rule

			drivestarttime += drivenextstart;

			elemind = 1;
			xoff = 0;

			Vector2D dpdx, gv2;

			Vector2D p2 = loccoords(2);
			p2.X() = XG(ns)-xoff; // - 0.5*drivevel*GetMBS()->GetStepSize();//correction for trapezoidal rule

			GetBody2D(elemind+1).GetdPosdx(p2, dpdx);

			double dpdxn = dpdx.Norm();
			if (dpdxn == 0) dpdxn = 1;
			dpdxn = 1./dpdxn;

			gv2 = dpdx*(dpdxn*drivevel);

			//reset position:
			if (GetMBS()->GetTime() >= drivenextstart)
			{
				GetBody2D(1).XG(1) = 0;
				GetBody2D(1).XG(2) = 0;
			}

			//if (driveacceltime == 0)
			{
				GetBody2D(1).XGP(1) = gv2.X(); 
				GetBody2D(1).XGP(2) = gv2.Y();
				XG(nsp) = drivevel;
			}

			//UO() << "v0=" << gv2.X() << "," << gv2.Y() << "\n";
			//UO() << "time=" << GetMBS()->GetTime() << gv2 << "\n";

			driveflag = 0;
		}
	}

	nlstepcnt = 0;
}

double SlidingJoint2D::PostNewtonStep(double t) 
{
	//return 0;
	nlstepcnt++;
	if (driveflag || nlstepcnt>=2) return 0;
	double chg = 0;
	double epstol = 1e-10;

	if (XG(ns)-xoff > 0.5*(1.+epstol)*GetBody2D(elemind+1).GetSize().X())
	{
		chg = 1;

		if (elemind+1 < elements.Length()) 
		{
			xoff += 0.5*(GetBody2D(elemind+1).GetSize()).X(); 
			xoff += 0.5*(GetBody2D(elemind+2).GetSize()).X(); 
			elemind++;
		}
		else chg = 0;
	} 
	else if (XG(ns)-xoff < -0.5*(1.+epstol)*GetBody2D(elemind+1).GetSize().X())
	{
		chg = 1;
		if (elemind+1 > 2) 
		{
			xoff -= 0.5*GetBody2D(elemind+1).GetSize().X(); 
			xoff -= 0.5*GetBody2D(elemind).GetSize().X(); 
			elemind--;
		}
		else chg = 0;
		//UO() << "->xoff=" << xoff << "\n";
	}

	//UO() << "elemind=" << elemind << " ";
	//UO() << "->xoff=" << xoff << "\n";
	return chg;
}

void SlidingJoint2D::ComputeDpDx(Vector2D& dpdx, Vector2D& gp2, Vector2D& gv2)
{
	Vector2D p2 = loccoords(2);
	p2.X() = XG(ns)-xoff;

	GetBody2D(elemind+1).GetdPosdx(p2, dpdx);

	//double dpdxn = dpdx.Norm();
	//if (dpdxn == 0) dpdxn = 1;
	//dpdxn = 1./dpdxn; //leads to wrong results with curved elements!

	gp2 = GetBody2D(elemind+1).GetPos2D(p2);
	gv2 = GetBody2D(elemind+1).GetVel2D(p2) + XG(nsp)*dpdx; //1./dpdx.Norm()*dpdx*XG(nsp) leads to wrong results in curved elements!
	//UO() << "dpdx=" << dpdx << "\n";


	int dosmoothing = 0;

	if (dosmoothing)
	{
		//const double smooth = 0.995;
		const double smoothdist = 0.1;

		if (XG(ns)-xoff > 0.5*GetBody2D(elemind+1).GetSize().X()-smoothdist && elemind+1 < elements.Length())
		{
			double fact = 1+0.5*(XG(ns)-xoff - 0.5*GetBody2D(elemind+1).GetSize().X()-smoothdist)/smoothdist;

			Vector2D dpdx2;
			double xoff2;
			Vector2D p22 = loccoords(2);

			xoff2 = xoff;
			xoff2 += 0.5*(GetBody2D(elemind+1).GetSize()).X(); 
			xoff2 += 0.5*(GetBody2D(elemind+2).GetSize()).X(); 

			p22.X() = XG(ns)-xoff2;

			GetBody2D(elemind+1+1).GetdPosdx(p22, dpdx2);
			double dpdxn2 = dpdx2.Norm();
			if (dpdxn2 == 0) dpdxn2 = 1;
			dpdxn2 = 1./dpdxn2;

			dpdx = (1-fact)*dpdx + fact*dpdx2;
			gp2 = (1-fact)*gp2 + fact*GetBody2D(elemind+1+1).GetPos2D(p22);
			gv2 = (1-fact)*gv2 + fact*(GetBody2D(elemind+1+1).GetVel2D(p22)+dpdx2*dpdxn2*XG(nsp));
		} 
		else if (XG(ns)-xoff < -0.5*GetBody2D(elemind+1).GetSize().X()+smoothdist && elemind+1 > 2)
		{
			double fact = 0.5*(1-(-(XG(ns)-xoff) - 0.5*GetBody2D(elemind+1).GetSize().X()-smoothdist)/smoothdist)-1;

			Vector2D dpdx2;
			double xoff2;
			Vector2D p22 = loccoords(2);

			xoff2 = xoff;
			xoff2 -= 0.5*GetBody2D(elemind+1).GetSize().X(); 
			xoff2 -= 0.5*GetBody2D(elemind).GetSize().X(); 

			p22.X() = XG(ns)-xoff2;

			GetBody2D(elemind+1-1).GetdPosdx(p22, dpdx2);
			double dpdxn2 = dpdx2.Norm();
			if (dpdxn2 == 0) dpdxn2 = 1;
			dpdxn2 = 1./dpdxn2;

			dpdx = (0.5+fact)*dpdx + (0.5-fact)*dpdx2;
			gp2 = (0.5+fact)*gp2 + (0.5-fact)*GetBody2D(elemind+1-1).GetPos2D(p22);
			gv2 = (0.5+fact)*gv2 + (0.5-fact)*(GetBody2D(elemind+1-1).GetVel2D(p22)+dpdx2*dpdxn2*XG(nsp));
		}
	}
}

void SlidingJoint2D::EvalG(Vector& f, double t) 
{
	// s, lambda1, lambda2, lambda3
	if (MaxIndex()==3)
	{

		Vector2D dpdx, gp2, gv2;
		ComputeDpDx(dpdx, gp2, gv2);
		Vector2D v = GetBody2D(1).GetPos2D(loccoords(1)) - gp2;

		if (!driveflag)
		{
			f(1) = v(1);
			f(2) = v(2);
			if (constrainrotation)
			{
				Vector2D p2 = loccoords(2);
				p2.X() = XG(ns)-xoff;
				f(4) = GetBody2D(1).GetAngle2D(loccoords(1)) - 	GetBody2D(elemind+1).GetAngle2D(p2);
			}

			if (driven)
			{
				double v = drivevel;
				if (driveacceltime != 0) UO() << "ERROR in Sliding Joint: driveacceltime must be zero for index 3 formulation\n";
				if (t < drivestarttime)
				{
					v = 0;
				}
				f(3) = XG(ns) - (t-drivestarttime)*v;
			}
			else
			{
				f(3) = dpdx(1)*XG(l1)+dpdx(2)*XG(l2);
			}
		}
		else
		{
			f(1) = XG(l1);
			f(2) = XG(l2);
			f(3) = XG(ns)-(drivelen*(1.+1e-15));
			if (constrainrotation)
			{
				f(4) = XG(lphi);
			}
		}
		//UO() << "evalfg=" << f << "\n";
	}
	else if (MaxIndex()<=2)
	{

		if (!driveflag)
		{
			Vector2D dpdx, gp2, gv2;
			ComputeDpDx(dpdx, gp2, gv2);
			Vector2D v = GetBody2D(1).GetVel2D(loccoords(1))-gv2;

			f(1) = v(1);
			f(2) = v(2);

			if (constrainrotation)
			{
				Vector2D p2 = loccoords(2);
				p2.X() = XG(ns)-xoff;
				double dphidx;
				GetBody2D(elemind+1).GetdAngle2Ddx(p2,dphidx);

				f(4) = GetBody2D(1).GetAngle2DP(loccoords(1)) - 
					(GetBody2D(elemind+1).GetAngle2DP(p2) + dphidx * XG(nsp));
			}

			if (1 && driven /*&& XG(nsp) > 0.01*drivevel && t <= drivestarttime+drivenextstart*/)
			{
				double v = drivevel;

				if (t < drivestarttime)
				{
					v = 0;
				} 
				else if (t <= drivestarttime+driveacceltime && driveacceltime!=0)
				{
					double p = (t-drivestarttime)/driveacceltime; //from 0 to 1
					//v = v*(0.5-0.5*cos(p*MY_PI)); //smooth otherwise lambdas start to oscillate heavily!
					v = v*p; //linear
				}

				f(3) = XG(nsp) - v;
				//f(3) = GetBody2D(1).GetVel2D(loccoords(1)).Norm() - fabs(v); //alternatively?
			}
			else
			{
				f(3) = dpdx(1)*XG(l1)+dpdx(2)*XG(l2);
			}

		}
		else if (driveflag == 1)
		{
			//UO() << "out!\n";
			f(1) = XG(l1);
			f(2) = XG(l2);
			f(3) = XG(nsp);
			if (constrainrotation)
			{
				f(4) = XG(lphi);
			}
		}
	}

};

void SlidingJoint2D::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	//-C_q^T \lambda = [dC/dq;]^T [\lambda_i] 

	if (driveflag) return;

	if (locelemind == 1)
	{
		GetBody2D(1).GetdPosdqT(loccoords(1),dpdq);

		for (int i=1; i <= f.Length(); i++)
		{
			f(i) -= (dpdq(i,1)*XG(l1)+dpdq(i,2)*XG(l2));
		}

		if (constrainrotation)
		{
			GetBody2D(1).GetdAngle2DdqT(loccoords(1),dpdq);
			for (int i=1; i <= f.Length(); i++)
			{
				f(i) -= (dpdq(i,1)*XG(lphi));
			}
		}
	}
	else if (locelemind == elemind+1)
	{
		Vector2D p = loccoords(elemind+1);
		p.X() = XG(ns)-xoff;
		GetBody2D(elemind+1).GetdPosdqT(p,dpdq);

		for (int i=1; i <= f.Length(); i++)
		{
			f(i) -= -(dpdq(i,1)*XG(l1)+dpdq(i,2)*XG(l2));
		}

		if (constrainrotation)
		{
			GetBody2D(elemind+1).GetdAngle2DdqT(p,dpdq);
			for (int i=1; i <= f.Length(); i++)
			{
				f(i) -= -(dpdq(i,1)*XG(lphi));
			}
		}

	}
};

Vector2D SlidingJoint2D::GetRefPos2DD()	const 
{
	return GetBody2D(1).GetPos2DD(loccoords(1));
}

Vector2D SlidingJoint2D::GetRefPos2D()	const 
{
	return GetBody2D(1).GetPos2D(loccoords(1));
}

double SlidingJoint2D::GetDrift() const
{
	Vector2D p2 = loccoords(2);
	p2.X() = XG(ns)-xoff;
	Vector2D v = GetBody2D(1).GetPos2D(loccoords(1))-GetBody2D(elemind+1).GetPos2D(p2);
	return v.Norm();
}

Vector2D SlidingJoint2D::GetSlidingPos()	const 
{
	Vector2D p2 = loccoords(2);
	p2.X() = XG(ns)-xoff;
	return GetBody2D(elemind+1).GetPos2D(p2);
}

void SlidingJoint2D::DrawElement() 
{
	Constraint::DrawElement();
	//


	mbs->SetColor(GetCol());
	Vector2D p1 = GetBody2D(1).GetPos2DD(loccoords(1));

	mbs->SetColor(colred);
	double x = XGD(ns);
	int n = 1;
	while (x > GetBody2D(n+1).GetSize().X()*0.5 && n+1 < elements.Length()) 
	{
		x -= 0.5*GetBody2D(n+1).GetSize().X(); 
		x -= 0.5*GetBody2D(n+2).GetSize().X(); 
		n++;
	}
	Vector2D p2 = GetBody2D(n+1).GetPos2DD(Vector2D(x,0));
	mbs->DrawSphere(GetBody2D(n+1).ToP3D(p2),0.5*draw_dim.X(),6);

	mbs->SetColor(Vector3D(0.1,0.8,0.1));
	mbs->DrawSphere(GetBody2D(1).ToP3D(p1),0.5*draw_dim.X(),5);

	if(0) 
	{
		char str[100];
		sprintf(str,"    drift = %18.15g, %18.15g   ", (p1-p2).X(), (p1-p2).Y());
		GetMBS()->GetRC()->PrintTextStruct(7+GetOwnNum(),-1,str);
	}
};










//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// CylindricalContact
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CylindricalContact::EvalF(Vector& f, double t) 
{
	Element::EvalF(f,t);
	//f(1) = XG(2);
}; 

void CylindricalContact::PostprocessingStep() 
{
	nlstepcnt = 0;
}

double CylindricalContact::PostNewtonStep(double t) 
{
	if (nlstepcnt++ > 2) return 0;
	double nlerror = 0;

	Vector2D x0 = MinDist();
	Vector3D p1 = loccoords(1);
	p1.X() = x0(1)-xoff;
	Vector3D p2 = loccoords(2);
	p2.X() = x0(2)-xoff;

	int searchnewcontact = 0;
	double tol = 1+1e-8;
	if (x0(2) > 0.5*tol*GetBody3D(elemind+1).GetSize().X() ||
		x0(2) <-0.5*tol*GetBody3D(elemind+1).GetSize().X()) 
	{
		searchnewcontact = 1;
	}

	//Contact-Vector:
	Vector3D v = GetBody3D(1).GetPos(p1)-GetBody3D(elemind+1).GetPos(p2);
	//Penetration:
	double gap = v.Norm()-cdist;

	if ((gap > 0 && switchvar == 0) || searchnewcontact)
	{
		//search for contact on other elements
		int oldind = elemind;

		int foundcontact = 0;
		double mingap = 0;
		int found_ind = 1;
		for (elemind=1; elemind < elements.Length(); elemind++)
		{
			Vector2D x0 = MinDist();
			Vector3D p1 = loccoords(1);
			p1.X() = x0(1)-xoff;
			Vector3D p2 = loccoords(2);
			p2.X() = x0(2)-xoff;
			//Contact-Vector:
			Vector3D v = GetBody3D(1).GetPos(p1)-GetBody3D(elemind+1).GetPos(p2);
			//Penetration:
			double gap = v.Norm()-cdist;
			if (gap <= mingap) 
			{
				mingap = gap;
				foundcontact = 1;
				nlerror = fabs(gap);
				found_ind = elemind;
			}
		}
		elemind = found_ind;
		if (!foundcontact) elemind = oldind;
	}
	double tolgap = 0;
	if (gap <= tolgap && switchvar == 0) 
	{
		switchvar = 1; //1
		nlerror = fabs(gap);
	}
	else
		if (gap > tolgap && switchvar == 1) 
		{
			switchvar = 0;
			nlerror = fabs(gap);
		}
		if (XG(1) > 0 && switchvar == 1) 
		{
			switchvar = 0; //0 !!!!!!!
			nlerror = fabs(XG(1));
		}
		return nlerror;
}

void CylindricalContact::ResMinDist(const Vector2D& x, Vector2D& res)
{
	Vector3D p1 = loccoords(1);
	p1.X() = x(1)-xoff;
	Vector3D p2 = loccoords(2);
	p2.X() = x(2)-xoff;
	Vector3D v = GetBody3D(1).GetPos(p1)-GetBody3D(elemind+1).GetPos(p2);

	Vector3D dpdx1, dpdx2;
	GetBody3D(1).GetdPosdx(p1, dpdx1);
	GetBody3D(elemind+1).GetdPosdx(p2, dpdx2);
	res(1) = v*dpdx1;
	res(2) = v*dpdx2;
}

void CylindricalContact::JacMinDist(const Vector2D& x, Matrix3D& jac)
{
	Vector2D x0, res0, res;
	x0 = x;
	ResMinDist(x0,res0);
	double eps = 1e-8;
	x0(1)+=eps;
	ResMinDist(x0,res);
	x0(1)-=eps;
	jac.SetSize(2,2);
	jac(1,1) = 1./eps*(res(1)-res0(1));
	jac(2,1) = 1./eps*(res(2)-res0(2));

	x0(2)+=eps;
	ResMinDist(x0,res);
	x0(2)-=eps;
	jac(1,2) = 1./eps*(res(1)-res0(1));
	jac(2,2) = 1./eps*(res(2)-res0(2));
}

//compute points which are nearest
Vector2D CylindricalContact::MinDist()
{
	Vector2D x0 = nlfinit;
	Vector2D res, xold;
	Matrix3D jac;
	//ResMinDist(x0,res);
	//double initerr = res.Norm();
	//if (initerr < 1e-6) {initerr = 1;}
	double tol = 1e-12;
	xold = x0+Vector2D(1000*tol,1000*tol);
	int maxit = 6; int it = 0;
	while(it < maxit && (xold-x0).Norm() > tol)
	{
		xold = x0;
		ResMinDist(x0,res);
		it++;
		JacMinDist(x0,jac);
		jac.Invert();
		//x0 -= jac*res;
		x0(1) -= jac(1,1)*res(1)+jac(1,2)*res(2);
		x0(2) -= jac(2,1)*res(1)+jac(2,2)*res(2);

		if (x0(1) > 0.5*GetBody3D(1).GetSize().X()) {x0(1) = 0.5*GetBody3D(1).GetSize().X();}
		if (x0(1) <-0.5*GetBody3D(1).GetSize().X()) {x0(1) =-0.5*GetBody3D(1).GetSize().X();}
		if (x0(2) > 0.51*GetBody3D(elemind+1).GetSize().X()) {x0(2) = 0.51*GetBody3D(elemind+1).GetSize().X();}
		if (x0(2) <-0.51*GetBody3D(elemind+1).GetSize().X()) {x0(2) =-0.51*GetBody3D(elemind+1).GetSize().X();}
		//UO() << "Mindist err = " << res.Norm() << ", it=" << it << "\n";
	}
	if (it < maxit) {nlfinit = x0;}


	return x0;
}

void CylindricalContact::EvalG(Vector& f, double t) 
{
	if (MaxIndex()>=1) 
	{
		Vector2D x0 = MinDist();
		Vector3D p1 = loccoords(1);
		p1.X() = x0(1)-xoff;
		Vector3D p2 = loccoords(2);
		p2.X() = x0(2)-xoff;

		//Contact-Vector:
		Vector3D v = GetBody3D(1).GetPos(p1)-GetBody3D(elemind+1).GetPos(p2);
		//Penetration:
		double gap = v.Norm()-cdist;
		Vector3D vv = GetBody3D(1).GetVel(p1)-GetBody3D(elemind+1).GetVel(p2);
		v.Normalize();
		//UO() << "gap=" << gap << "\n";
		if (switchvar)
		{
			double cdamp = cstiff*0.01;
			f(1) = XG(1) - gap*cstiff-vv*v*cdamp; //XG(1) = gap * cstiff  -> pressure is negative, XG(1) always negative!
		}
		else
		{
			f(1) = XG(1);
		}
	}
	else
	{
		UO() << "Index0 not yet implemented\n";
	}

};

void CylindricalContact::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	//-C_q^T \lambda = [dC/dq;]^T [\lambda_i] 

	if (!switchvar) return;

	if (locelemind == elemind+1)
	{
		Vector2D x0 = MinDist();
		Vector3D p1 = loccoords(1);
		p1.X() = x0(1)-xoff;
		Vector3D p2 = loccoords(2);
		p2.X() = x0(2)-xoff;

		//Contact-Vector:
		Vector3D v = GetBody3D(1).GetPos(p1)-GetBody3D(elemind+1).GetPos(p2);
		v.Normalize();
		GetBody3D(elemind+1).GetdPosdqT(p2,dpdq);

		double fc = XG(1);

		for (int i=1; i <= f.Length(); i++)
		{
			f(i) -= -(dpdq(i,1)*fc*v(1)+dpdq(i,2)*fc*v(2)+dpdq(i,3)*fc*v(3));
		}
	}
	else if (locelemind == 1)
	{
		Vector2D x0 = MinDist();
		Vector3D p1 = loccoords(1);
		p1.X() = x0(1)-xoff;
		Vector3D p2 = loccoords(2);
		p2.X() = x0(2)-xoff;

		//Contact-Vector:
		Vector3D v = GetBody3D(1).GetPos(p1)-GetBody3D(elemind+1).GetPos(p2);
		v.Normalize();

		GetBody3D(1).GetdPosdqT(p1,dpdq);
		double fc = XG(1);

		for (int i=1; i <= f.Length(); i++)
		{
			f(i) -= (dpdq(i,1)*fc*v(1)+dpdq(i,2)*fc*v(2)+dpdq(i,3)*fc*v(3));
		}
	}
};

Vector3D CylindricalContact::GetRefPosD()	const 
{
	return GetBody3D(1).GetPosD(loccoords(1));
}

Vector3D CylindricalContact::GetRefPos()	const 
{
	return GetBody3D(1).GetPos(loccoords(1));
}

void CylindricalContact::DrawElement() 
{
	Constraint::DrawElement();


	//How or Why draw Contact?
	//mbs->SetColor(col);
	//mbs->DrawSphere(GetBody3D(1).GetPosD(loccoords(1)),0.5*draw_dim.X(),16);
};



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// GeneralizedAngleConstraint
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


void GeneralizedAngleConstraint::EvalG(Vector& f, double t) 
{
	if (penalty) return;

	Body3D& b1 = GetBody3D(1);
	Body3D& b2 = GetBody3D(2);
	Vector3D r1(b1.XG((c1-1)*3+1),b1.XG((c1-1)*3+2),b1.XG((c1-1)*3+3));
	Vector3D r2(b2.XG((c2-1)*3+1),b2.XG((c2-1)*3+2),b2.XG((c2-1)*3+3));

	if (MaxIndex()>=3) 
	{
		f(1) = (r1*r2-calpha);
		//UO() << "f1=" << f(1) << "\n";
	}
	else if (MaxIndex()>=2) 
	{
		Vector3D r1p(b1.XGP((c1-1)*3+1),b1.XGP((c1-1)*3+2),b1.XGP((c1-1)*3+3));
		Vector3D r2p(b2.XGP((c2-1)*3+1),b2.XGP((c2-1)*3+2),b2.XGP((c2-1)*3+3));
		f(1) = r1p*r2+r1*r2p;
	}
	else
	{
		UO() << "Index0 not yet implemented\n";
	}

};

void GeneralizedAngleConstraint::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	//-C_q^T \lambda = [dC/dq;]^T [\lambda_i] 

	Body3D& b1 = GetBody3D(1);
	Body3D& b2 = GetBody3D(2);
	Vector3D r1(b1.XG((c1-1)*3+1),b1.XG((c1-1)*3+2),b1.XG((c1-1)*3+3));
	Vector3D r2(b2.XG((c2-1)*3+1),b2.XG((c2-1)*3+2),b2.XG((c2-1)*3+3));
	dpdq.SetSize(f.GetLen(),1);
	dpdq.SetAll(0);

	if (locelemind == 1)
	{
		dpdq((c1-1)*3+1,1) = r2.X();
		dpdq((c1-1)*3+2,1) = r2.Y();
		dpdq((c1-1)*3+3,1) = r2.Z();
	}
	else
	{
		dpdq((c2-1)*3+1,1) = r1.X();
		dpdq((c2-1)*3+2,1) = r1.Y();
		dpdq((c2-1)*3+3,1) = r1.Z();
	}

	double lambda;
	if (!penalty) lambda = XG(1);
	else lambda = (r1*r2-calpha)*2*penaltyfact;

	for (int i=1; i <= f.Length(); i++)
	{
		f(i) -= dpdq(i,1)*lambda;
	}
};

Vector3D GeneralizedAngleConstraint::GetRefPosD()	const 
{
	return GetBody3D(1).GetPosD(0);
}

void GeneralizedAngleConstraint::DrawElement() 
{
	Constraint::DrawElement();
	//


	//position hard to determine ... not stored ... c1 < or > 4
	mbs->SetColor(GetCol());
	double side = -1;
	if (c1 > 4) {side = 1;}
	mbs->DrawSphere(GetBody3D(1).GetPosD(Vector3D(side*0.5*GetBody3D(1).GetSize().X(),0,0)),0.5*draw_dim.X(),16);
};






//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// PrescribedAngleConstraint
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void PrescribedAngleConstraint::EvalG(Vector& f, double t)
{
	if (mode == 1)
	{
		Body3D& b1 = GetBody3D(1);
		Vector3D lp1 = loccoords(1);
		Vector3D p1 = GetBody3D(1).GetPos(lp1);
		if (MaxIndex()>=3)
		{
			Vector3D vrot = GetRotVector(Angle(GetMBS()->GetTime()));
			//Matrix3D A=GetBody3D(1).GetRotMatrix(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Vector3D Axv = GetBody3D(1).GetRotMatv(lp1,GetNormalRotVector());
			f(1) = vrot*(Axv);
			//UO() << "Evalg=" << f(1) << "\n";
		}
		else
			if (MaxIndex()>=2)
			{
				Vector3D vrot = GetRotVector(Angle(GetMBS()->GetTime()));
				Vector3D vrotp = GetRotVectorp(Angle(GetMBS()->GetTime()),Anglep(GetMBS()->GetTime()));
				//Matrix3D A=GetBody3D(1).GetRotMatrix(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
				//Matrix3D Ap=GetBody3D(1).GetRotMatrixP(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
				Vector3D Axv = GetBody3D(1).GetRotMatv(lp1,GetNormalRotVector());
				Vector3D Apxv = GetBody3D(1).GetRotMatPv(lp1,GetNormalRotVector());
				f(1) = vrot*(Apxv)+vrotp*(Axv);
				//UO() << "Evalg=" << f(1) << "\n";
			}
			else
			{
				UO() << "This Index not yet implemented\n";
			}
	}
	else
	{
		Body2D& b1 = GetBody2D(1);
		Vector2D lp1(loccoords(1).X(), loccoords(1).Y());
		if (MaxIndex()>=3)
		{
			f(1) = b1.GetAngle2D(lp1) - Angle(GetMBS()->GetTime());
		}
		else if (MaxIndex()>=2)
		{
			f(1) = b1.GetAngle2DP(lp1) - Anglep(GetMBS()->GetTime());
		}
	}
};
void PrescribedAngleConstraint::AddElementCqTLambda(double t, int locelemind, Vector& f)
{
	//-C_q^T \lambda = [dC/dq;]^T [\lambda_i]
	if (mode == 1) //3D
	{
		dpdq.SetSize(f.GetLen(),3);
		dpdq.SetAll(0);
		Body3D& b1 = GetBody3D(1);
		Vector3D lp1 = loccoords(1);
		Vector3D p1 = b1.GetPos(lp1);
		Vector3D vrot = GetRotVector(Angle(GetMBS()->GetTime()));
		b1.GetdRotvdqT(GetNormalRotVector(),lp1,dpdq);
		for (int i=1; i <= f.Length(); i++)
		{
			f(i) -= (dpdq(i,1)*vrot(1)+dpdq(i,2)*vrot(2)+dpdq(i,3)*vrot(3))*XG(1);
		}
	}
	else //planar
	{
		dpdq.SetSize(f.GetLen(),1);
		dpdq.SetAll(0);
		Body2D& b1 = GetBody2D(1);
		Vector2D lp1(loccoords(1).X(), loccoords(1).Y());
		b1.GetdAngle2DdqT(lp1,dpdq);
		for (int i=1; i <= f.Length(); i++)
		{
			f(i) -= dpdq(i,1) * XG(1);
		}
	}
};
Vector3D PrescribedAngleConstraint::GetRefPosD() const
{
	if (mode == 1)
		return GetBody3D(1).GetPosD(Vector3D(0.,0.,0.));
	else
		return GetBody2D(1).GetPosD(Vector3D(0.,0.,0.));
}
void PrescribedAngleConstraint::DrawElement()
{
	Constraint::DrawElement();
	//


	mbs->SetColor(GetCol());
	Vector3D p1;
	if (mode == 1) //3D
	{
		p1 = (GetBody3D(1).GetPosD(loccoords(1)));
		Vector3D v1 = (2.*draw_dim.X())*GetRotVector(Angle(GetMBS()->GetDrawTime())-MY_PI*0.5);
		mbs->DrawZyl(p1,p1+v1,0.25*draw_dim.X(),8);
	}
	else //mode==2, planar
	{
		p1 = GetBody2D(1).ToP3D(GetBody2D(1).GetPos2DD(Vector2D(loccoords(1).X(),loccoords(1).Y())));
	}
	mbs->DrawSphere(p1,0.5*draw_dim.X(),16);
};




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// PrescribedAngularVel
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void PrescribedAngularVel::GetPrescribedOmega(Vector3D& prescribed_omega, double t) const
{

	Vector3D f(1.,1.,1.);

	switch(ptype)
	{
	case 1:
		{
			if (t*data1(1) < MY_PI) f.X() = 0.5;
			if (t*data2(1) < MY_PI) f.Y() = 0.5;
			if (t*data3(1) < MY_PI) f.Z() = 0.5;

			prescribed_omega.X() = f.X()*data1(3)*sin(data1(1)*t + data1(2));
			prescribed_omega.Y() = f.Y()*data2(3)*sin(data2(1)*t + data2(2));
			prescribed_omega.Z() = f.Z()*data3(3)*sin(data3(1)*t + data3(2)); 
			//if (t*data2(1) > MY_PI*0.5) prescribed_omega.Y() = data2(3); //prescribed, special
			break;
		}
	default: prescribed_omega = Vector3D(0.,0.,0.);
	}
}

void PrescribedAngularVel::EvalG(Vector& f, double t) 
{
	Vector3D lp1 = loccoords(1);

	if (MaxIndex()>=2) 
	{
		Vector3D omega = GetBody3D(1).GetAngularVel(lp1);

		Vector3D prescribed_omega;
		GetPrescribedOmega(prescribed_omega, t);

		f(1) = omega.X() - prescribed_omega.X();
		f(2) = omega.Y() - prescribed_omega.Y();
		f(3) = omega.Z() - prescribed_omega.Z();
		/*
		f(1) = XG(1);
		f(2) = XG(2);
		f(3) = XG(3);*/
	}
	else
	{
		UO() << "This Index not yet implemented\n";
	}
};

void PrescribedAngularVel::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	//-C_q^T \lambda = [dC/dq;]^T [\lambda_i] 

	dpdq.SetSize(f.GetLen(),3);
	dpdq.SetAll(0);

	Vector3D lp1 = loccoords(1);

	GetBody3D(1).GetdAngVeldqpT(lp1, dpdq);

	for (int i=1; i <= f.Length(); i++)
	{
		f(i) -= (dpdq(i,1)*XG(1) + dpdq(i,2)*XG(2) + dpdq(i,3)*XG(3));
	}
	/*
	Vector ff1, ff2, ff3;
	ff1.SetLen(f.Length());
	ff2.SetLen(f.Length());
	ff3.SetLen(f.Length());
	GetBody3D(1).ApplyDrotrefdq(ff1, Vector3D(1.,0.,0.));
	GetBody3D(1).ApplyDrotrefdq(ff2, Vector3D(0.,1.,0.));
	GetBody3D(1).ApplyDrotrefdq(ff3, Vector3D(0.,0.,1.));
	*/
	/*
	UO() << "G=" << dpdq << "\n";
	UO() << "l1=" << XG(1) << ", l2=" << XG(2) << ", l3=" << XG(3) << "\n";
	UO() << "angv=" << GetBody3D(1).GetAngularVel(Vector3D(0.,0.,0.)) << "\n";
	*/
};

Vector3D PrescribedAngularVel::GetRefPosD()	const 
{
	return GetBody3D(1).GetPosD(Vector3D(0.,0.,0.));
}

void PrescribedAngularVel::DrawElement() 
{
	Constraint::DrawElement();
	//


	//position hard to determine ... not stored ... c1 < or > 4????
	mbs->SetColor(GetCol());
	Vector3D p1(GetBody3D(1).GetPosD(loccoords(1)));
	mbs->DrawSphere(p1,0.5*draw_dim.X(),16);
};




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// GravityConstraint
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void GravityConstraint::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	//-C_q^T \lambda = [dC/dq;]^T [\lambda_i] 

	Body3D& b1 = GetBody3D(1);
	Body3D& b2 = GetBody3D(2);
	Vector3D r1(b1.XG(1),b1.XG(2),b1.XG(3));
	Vector3D r2(b2.XG(1),b2.XG(2),b2.XG(3));
	dpdq.SetSize(f.GetLen(),1);

	double sign = 1;
	if (locelemind == 2)
	{
		sign = -1;
	}

	Vector3D r12 = r2-r1;
	double l12 = r12.Norm();
	double crit = 1e9;
	if (l12 < crit) l12 = (l12+1./crit*Sqr(crit-l12));
	Vector3D F = (sign*gconstant*b1.GetMass()*b2.GetMass()/Cub(l12))*r12;

	for (int i=1; i <= f.Length(); i++)
	{
		f(i) += F(i); //external force -> plus sign
	}
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// GroundContact
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//compute points which are nearest, compute Normal vector to ground
double GroundContact::MinDist(Vector3D& n, int& nquad)
{

	//assume Body to be of type sphere
	double r = cdist;
	double mindist = 1e100;

	if (!ground->NQuads()) return 0;

	if (nquad == 0)
	{
		int foundquad = 0;
		Box3D b(GetBody3D(1).GetElementBox());
		//UO() << "box=" << b << ", r=" << r << "\n";
		b.Increase(b.Radius()*0.5);

		//UO() << "box2=" << b << "\n";
		mystatic TArray<int> quadnums;

		ground->GetQuadsInBox(b,quadnums);
		//UO() << "quads=" << quadnums << "\n";

		//for (int i = 1; i <= ground->NQuads(); i++)
		for (int ii = 1; ii <= quadnums.Length(); ii++)
		{
			int i = quadnums(ii);
			if (b.Intersect(ground->Boxes(i)))
			{
				int4 quad = ground->GetQuad(i);
				const Vector3D& p1 = ground->PT(quad.Get(1));
				const Vector3D& p2 = ground->PT(quad.Get(2));
				const Vector3D& p3 = ground->PT(quad.Get(3));
				const Vector3D& p4 = ground->PT(quad.Get(4));

				Vector3D p = GetBody3D(1).GetPos(loccoords(1));

				double dist = DistToQuad(p1,p2,p3,p4,p)-r;
				if (dist < mindist)
				{
					mindist = dist;
					foundquad = i;
				}
			}
		}
		//UO() << "found=" << foundquad << "\n";

		//UO() << "midist=" << mindist << "\n";
		//UO() << "foundquad=" << foundquad << "\n";

		if (!foundquad) 
		{	
			foundquad = 1;
			mindist = r;
		}
		int4 quad = ground->GetQuad(foundquad);
		const Vector3D& p1 = ground->PT(quad.Get(1));
		const Vector3D& p2 = ground->PT(quad.Get(2));
		const Vector3D& p4 = ground->PT(quad.Get(4));
		Normal3D(p1,p2,p4,n);

		nquad = foundquad;
	}
	else
	{
		int4 quad = ground->GetQuad(nquad);
		const Vector3D& p1 = ground->PT(quad.Get(1));
		const Vector3D& p2 = ground->PT(quad.Get(2));
		const Vector3D& p3 = ground->PT(quad.Get(3));
		const Vector3D& p4 = ground->PT(quad.Get(4));

		Vector3D p = GetBody3D(1).GetPos(loccoords(1));

		mindist = DistToQuad(p1,p2,p3,p4,p)-r;
		Normal3D(p1,p2,p4,n);
	}

	return mindist;
}

//ofstream testf("..\\..\\output\\test.txt");
void GroundContact::PostprocessingStep() 
{
	nlstepcnt = 0;


	Vector3D n;
	int q=0;
	double gap = MinDist(n, usedquad);
	//testf << GetMBS()->GetTime() << " " << gap << " " << q << endl;

}

double GroundContact::PostNewtonStep(double t) 
{
	double lambdatol = 1e-16*cstiff;
	if (nlstepcnt++ > 2) 
	{
		return 0;
		if (switchvar == 1 && XG(1) < lambdatol) return 0;
	}
	double nlerror = 0;

	Vector3D n;

	int oldquad = usedquad;
	usedquad = 0;
	double gap = MinDist(n, usedquad);

	//if (t>12.) UO() << "step=" << nlstepcnt << ", gap=" << gap << ", usedquad=" << usedquad << ", oldquad=" << oldquad << "switchvar=" << switchvar << "lambda=" << XG(1) << "\n";

	double tolgap = 1e-6;
	if (oldquad != usedquad) 
	{
		nlerror = 1e10;

		if (gap <= 0) 
			switchvar = 1;
		else
			switchvar = 0;
	}
	else
	{
		if (gap <= 0 && switchvar == 0) 
		{
			switchvar = 1;
			nlerror = fabs(gap);
		}
		else if ((XG(1) > lambdatol /*|| gap > tolgap*/) && switchvar == 1) 
		{
			switchvar = 0;
			nlerror = Minimum(fabs(XG(1))/cstiff, fabs(gap));
		}
	}
	return nlerror;
}

void GroundContact::EvalG(Vector& f, double t) 
{
	if (MaxIndex()>=1) 
	{
		Vector3D n;
		double gap = MinDist(n, usedquad);
		int hardcontact = 1;

		if (switchvar)
		{
			Vector3D lp1 = loccoords(1);
			Vector3D vv = GetBody3D(1).GetVel(lp1);

			if (hardcontact)
			{
				f(1) = vv*n;
			}
			else
			{
				double fdamp = vv*n*cstiff*0.01;
				//double sgn = Sgn(fdamp);
				//if (fabs(fdamp)  > fabs(gap*cstiff)) fdamp = sgn*fabs(gap*cstiff);
				f(1) = XG(1) - (gap*cstiff-fdamp); //XG(1) = gap * cstiff  -> pressure is negative, XG(1) always negative!
			}
		}
		else
		{
			f(1) = XG(1);
		}
	}
	else
	{
		UO() << "Index0 not yet implemented\n";
	}

};

void GroundContact::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	//-C_q^T \lambda = [dC/dq;]^T [\lambda_i] 

	if (!switchvar) return;

	if (locelemind == 1)
	{
		Vector3D n;
		double gap = MinDist(n, usedquad);
		Vector3D lp1 = loccoords(1);

		double friccoeff = 0.1*0;
		Vector3D ffric(0,0,0);

		if (friccoeff)
		{
			Vector3D t1,t2;
			n.SetNormalBasis(t1,t2);
			Vector3D vv = GetBody3D(1).GetVel(lp1);
			vv = (-(vv*t1)*t1-(vv*t2)*t2);

			double limitv = 0.5;
			if (vv.Norm() < limitv) friccoeff *= To4(vv.Norm()/(limitv));
			//if (vv.Norm() < 1) friccoeff = 0;

			vv.Normalize();
			ffric = friccoeff*vv*XG(1);
		}


		GetBody3D(1).GetdPosdqT(lp1,dpdq);
		double fc = -XG(1);

		for (int i=1; i <= f.Length(); i++)
		{
			f(i) -= (dpdq(i,1)*(fc*n(1)+ffric(1))+dpdq(i,2)*(fc*n(2)+ffric(2))+dpdq(i,3)*(fc*n(3)+ffric(3)));
		}
	}
};

Vector3D GroundContact::GetRefPosD()	const 
{
	return GetBody3D(1).GetPosD(loccoords(1));
}

Vector3D GroundContact::GetRefPos()	const 
{
	return GetBody3D(1).GetPos(loccoords(1));
}

void GroundContact::DrawElement() 
{
	Constraint::DrawElement();
	//


	//UO() << "ground-box=" << ground->GetBoundingBox() << "\n";
	//ground->DrawYourself();
	mbs->SetColor(GetCol());
	mbs->DrawSphere(GetBody3D(1).GetPosD(loccoords(1)),0.5*draw_dim.X(),8);

};

Box3D GroundContact::GetElementBoxD() const
{
	return ground->GetBoundingBox();
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// CircleContact2D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CircleContact2D::PostprocessingStep() 
{
	if (!iscontact) gapp_init = gapp_init_min;

	nlstepcnt = 0;
}

double CircleContact2D::PostNewtonStep(double t) 
{
	if (contactmode == 1 || contactmode == 2) //always contact
	{
		if (!iscontact)
		{
			iscontact = 1;
			return 1;
		}
		else return 0;
	}

	if (nlstepcnt++ > 5) return 0;
	double nlerror = 0;

	Vector2D v;
	if (!DOFmode)
	{
		v = GetBody2D(1).GetPos2D(loccoords(1))-p_global;
	}
	else
	{
		v = GetBody2D(1).GetNodePos2D((int)loccoords(1).X()) - p_global;
	}

	double gap = v.Norm() - r;
	v.Normalize();

	double gapp;
	if (!DOFmode)
	{
		gapp = -(GetBody2D(1).GetVel2D(loccoords(1))*v);
	}
	else
	{
		gapp = -(GetBody2D(1).GetNodeVel2D((int)loccoords(1).X())*v);
	}

	if (gap < 0 && !iscontact)
	{
		//GetMBS()->UO() << "contact!\n";
		iscontact = 1;

		nlerror += fabs(gap);
	} 
	else 	if (gap > 0 && iscontact)
	{
		iscontact = 0;
		nlerror += fabs(XG(1))/GetPenaltyStiffness();
	}
	if (gapp_init < fabs(gapp) && iscontact)
	{
		nlerror = fabs(gapp_init - fabs(gapp))/Maximum(gapp_init,fabs(gapp));
		gapp_init = fabs(gapp);
	}


	return nlerror;
}


void CircleContact2D::EvalG(Vector& f, double t) 
{
	if (elements.Length()!=1)
	{
		mbs->UO() << "ERROR: CircleContact2D::EvalG, number of elements != 1\n"; return;
	}

	if (MaxIndex()<=3)
	{
		if (elements.Length()==1)
		{

			if (!iscontact)
			{
				f(1) = XG(1);
			}
			else 
			{
				Vector2D v;
				if (!DOFmode)
				{
					v = GetBody2D(1).GetPos2D(loccoords(1))-p_global;
				}
				else
				{
					v = GetBody2D(1).GetNodePos2D((int)loccoords(1).X()) - p_global;
				}

				double r1 = v.Norm();
				v.Normalize();
				double gap = r1-r;
				double gapp;
				if (!DOFmode)
				{
					gapp = -(GetBody2D(1).GetVel2D(loccoords(1))*v);
				}
				else
				{
					gapp = -(GetBody2D(1).GetNodeVel2D((int)loccoords(1).X())*v);
				}


				double e = 0.95; //coeff of restitution; 0.95 ususally
				double n = 1; //coeff for Hertzian contact; 1 works good
				double gi = gapp_init;
				if (fabs(gapp) > gapp_init) gi = fabs(gapp);


				//f(1) = -gap*c_contact-XG(1); //XG(1) is compression force, linear contact model
				switch (contactmode)
				{
				case 0: //Hertian contact
					f(1) = -Sgn(gap)*GetPenaltyStiffness()*(1.+0.75*(1.-Sqr(e))*gapp/gi)*pow(fabs(gap),n)-XG(1); //XG(1) is compression force, Hertzian contact model
					break;
				case 1: //tied, linear elastic
					f(1) = -GetPenaltyStiffness()*gap - XG(1);
					break;
				case 2: //tied, rigid, Lagrange multiplier
					f(1) = gapp;
					break;
				case 3: //linear contact model
					f(1) = -GetPenaltyStiffness()*gap - XG(1);
					break;
				default: UO() << "ERROR: contactmode\n"; f(1) = 0;
				}
			}
		}
	}

};

void CircleContact2D::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	//-C_q^T \lambda = [dC/dx; dC/dy; dC/dphi]^T [\lambda1;\lambda2] 
	//C = p_ref+R*p_loc
	double sign = 1;
	if (locelemind == 2) sign = -1;

	if (iscontact)
	{
		Vector2D v;
		if (!DOFmode)
		{
			v = GetBody2D(1).GetPos2D(loccoords(1))-p_global;
			GetBody2D(locelemind).GetdPosdqT(loccoords(locelemind),dpdq);
			v.Normalize();

			for (int i=1; i <= f.Length(); i++)
			{
				f(i) -= -sign*(dpdq(i,1)*XG(1)*v.X()+dpdq(i,2)*XG(1)*v.Y());
			}
		}
		else
		{
			v = GetBody2D(1).GetNodePos2D((int)loccoords(1).X()) - p_global;
			v.Normalize();

			Vector2D lambda(-sign*XG(1)*v.X(), -sign*XG(1)*v.Y()); //negative sign: needs to be subtracted from f!!!!
			GetBody2D(locelemind).AddNodedPosdqTLambda((int)loccoords(locelemind).X(),lambda,f);

			/*
			GetBody2D(locelemind).GetNodedPosdqT((int)loccoords(locelemind).X(),dpdq);
			for (int i=1; i <= f.Length(); i++)
			{
			f(i) -= -sign*(dpdq(i,1)*XG(1)*v.X()+dpdq(i,2)*XG(1)*v.Y());
			}*/

			//old:
			/*
			v = GetBody2D(1).GetNodePos2D((int)loccoords(1).X()) - p_global;
			v.Normalize();

			if (0)
			{
			Vector2D lambda(sign*XG(1)*v.X(), sign*XG(1)*v.Y());
			GetBody2D(locelemind).AddNodedPosdqTLambda((int)loccoords(locelemind).X(),lambda,f);
			}
			else
			{
			GetBody2D(locelemind).GetNodedPosdqT((int)loccoords(locelemind).X(),dpdq);
			for (int i=1; i <= f.Length(); i++)
			{
			f(i) -= -sign*(dpdq(i,1)*XG(1)*v.X()+dpdq(i,2)*XG(1)*v.Y());
			}
			}*/


		}
	}
};

Vector3D CircleContact2D::GetRefPosD() const 
{
	Vector2D v;
	if (!DOFmode)
	{
		v = GetBody2D(1).GetPos2DD(loccoords(1));
	}
	else
	{
		v = GetBody2D(1).GetNodePos2DD((int)loccoords(1).X());
	}
	return Vector3D(v.X(),v.Y(),0);
}

void CircleContact2D::DrawElement() 
{
	Constraint::DrawElement();
	//


	Vector2D v;
	if (!DOFmode)
	{
		v = GetBody2D(1).GetPos2DD(loccoords(1));
	}
	else
	{
		//int ni = elements(1);
		v = GetBody2D(1).GetNodePos2DD((int)loccoords(1).X());
	}

	//draw contact node:
	if (GetMBS()->GetIOption(112))
	{
		mbs->SetColor(GetCol());
		mbs->DrawSphere(GetBody2D(1).ToP3D(v),0.5*GetDrawSizeScalar(),2);
	}

	/*
	if (elements.Length()==2)
	{
	mbs->SetColor(Vector3D(0.8,0.2,0.2));
	mbs->DrawSphere(GetBody2D(2).ToP3D(GetBody2D(2).GetPos2DD(loccoords(2))),0.5*draw_dim.X(),4);
	}*/

	Vector3D pz1 = GetBody2D(1).ToP3D(p_global);
	//mbs->DrawZyl(pz1,pz1+Vector3D(0,0,0.001), r, draw_dim.Y());

	if (GetMBS()->GetIOption(113))
	{

		//show fixed circle (rigid contact body)
		pz1+=Vector3D(0,0,-0.0004);
		for (int i = 1; i <= GetDrawSizeResolution(); i++)
		{
			double phi1 = (double)(i-1)/GetDrawSizeResolution()*2.*MY_PI;
			double phi2 = (double)(i)/GetDrawSizeResolution()*2.*MY_PI;
			mbs->DrawQuad(pz1,pz1,pz1+GetBody2D(1).ToP3D(Vector2D(r*sin(phi1),r*cos(phi1))),
				pz1+GetBody2D(1).ToP3D(Vector2D(r*sin(phi2),r*cos(phi2))));
		}

		//show fixed circle (rigid contact body)
		mbs->SetColor(Vector3D(0.7,0.7,0.1));
		pz1+=Vector3D(0,0,0.0002);
		for (int i = 1; i <= GetDrawSizeResolution(); i++)
		{
			double phi1 = (double)(i-1)/GetDrawSizeResolution()*2.*MY_PI;
			double phi2 = (double)(i)/GetDrawSizeResolution()*2.*MY_PI;
			mbs->DrawQuad(pz1,pz1,pz1+GetBody2D(1).ToP3D(Vector2D(r*0.99*sin(phi1),r*0.99*cos(phi1))),
				pz1+GetBody2D(1).ToP3D(Vector2D(r*0.99*sin(phi2),r*0.99*cos(phi2))));
		}
	}
};








//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// IntegralPosConstraint2D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void IntegralPosConstraint2D::EvalG(Vector& f, double t) 
{
	if (MaxIndex()==3) //position constraint
	{
		Vector2D v(0,0);

		for (int i=1; i <= loccoords.Length(); i++)
		{
			if (DOFmode == 0)
			{
				v += weights(i)*GetBody2D(i).GetPos2D(loccoords(i));
			}
			else
			{
				UO() << "not implemented integral constraint 2D, index 3, #DOFmode > 0\n";
				v += weights(i)*GetBody2D(i).GetNodePos2D((int)loccoords(i).X());
			}
		}
		v /= (double)loccoords.Length(); //could be eliminated
		v -= p_global;

		f(1) = v(1);
		f(2) = v(2);
	}
	else if (MaxIndex()<=2) //velocity constraint
	{
		Vector2D v(0,0);

		if (DOFmode <= 3)
		{
			for (int i=1; i <= loccoords.Length(); i++)
			{
				if (DOFmode == 0)
				{
					v += weights(i) * GetBody2D(i).GetVel2D(loccoords(i));
				}
				else
				{
					if (i == loccoords.Length() && DOFmode == 3) 
					{
						v -= weights(i) * GetBody2D(i).GetVel2D(loccoords(i));
					}
					else
					{
						v += weights(i) * GetBody2D(i).GetNodeVel2D((int)loccoords(i).X());
					}
				}
			}
			//v /= (double)loccoords.Length(); //could be eliminated
			f(1) = v(1);
			f(2) = v(2);
		}
		else if (DOFmode == 4)
		{
			if (maxmoment == 0)
			{
				for (int i=1; i <= npoints1; i++)
				{
					v += weights(i) * GetBody2D(i).GetNodeVel2D((int)loccoords(i).X());
				}
				for (int i=npoints1+1; i <= npoints2+npoints1; i++)
				{
					v -= weights(i) * GetBody2D(i).GetVel2D(loccoords(i));
				}
				f(1) = v(1);
				f(2) = v(2);
			}
			else
			{
				double fact;
				for (int p = 0; p <= maxmoment; p++)
				{
					v = Vector2D(0.,0.);
					for (int i=1; i <= npoints1; i++)
					{
						if (p == 0) fact = 1;
						else fact = pow(s_body1(i), p)/length1;
						v += fact * weights(i) * GetBody2D(i).GetNodeVel2D((int)loccoords(i).X());
					}
					for (int i=npoints1+1; i <= npoints2+npoints1; i++)
					{
						if (p == 0) fact = 1;
						else fact = pow(s_body2(i-npoints1), p)/length2;
						v -= fact * weights(i) * GetBody2D(i).GetVel2D(loccoords(i));
					}
					f(1+p*2) = v(1);
					f(2+p*2) = v(2);
				}
			}
		}
	}
};

void IntegralPosConstraint2D::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	//double sign = 1;
	//if (locelemind == 2) sign = -1;

	if (DOFmode == 0) //all elements different, local coordinates
	{
		GetBody2D(locelemind).GetdPosdqT(loccoords(locelemind),dpdq);
		Vector2D lam(XG(1)/(double)loccoords.Length(), XG(2)/(double)loccoords.Length());
		for (int i=1; i <= f.Length(); i++)
		{
			f(i) -= weights(i)*(dpdq(i,1)*lam(1)+dpdq(i,2)*lam(2));
		}
	}
	else if (DOFmode == 1) //all elements different, node numbers
	{
		Vector2D lam(-weights(locelemind)*XG(1)/(double)loccoords.Length(), 
			-weights(locelemind)*XG(2)/(double)loccoords.Length()); //negative sign: needs to be subtracted from f!!!!
		GetBody2D(locelemind).AddNodedPosdqTLambda((int)loccoords(locelemind).X(),lam,f);
	}
	else if (DOFmode == 2) //all elements are the same, different nodes
	{
		Vector2D lam(-XG(1)/(double)loccoords.Length(), -XG(2)/(double)loccoords.Length()); //negative sign: needs to be subtracted from f!!!!
		Vector2D lam2;
		for (int i=1; i <= loccoords.Length(); i++)
		{
			lam2 = weights(i)*lam; 
			GetBody2D(i).AddNodedPosdqTLambda((int)loccoords(i).X(),lam2,f);
		}
	}
	else if (DOFmode == 3) //all elements are the same except last element, different nodes
	{
		if (locelemind == 1)
		{
			Vector2D lam(XG(1), XG(2)); //negative sign: needs to be subtracted from f!!!!
			Vector2D lam2;
			for (int i=1; i <= loccoords.Length()-1; i++)
			{
				lam2 = (-weights(i))*lam;
				GetBody2D(i).AddNodedPosdqTLambda((int)loccoords(i).X(),lam2,f);
			}
		}
		else if (locelemind == loccoords.Length())
		{
			GetBody2D(loccoords.Length()).GetdPosdqT(loccoords(loccoords.Length()),dpdq);
			//Vector2D lam(XG(1), XG(2));
			//GetMBS()->UO() << "lam=" << lam << "\n";

			for (int i=1; i <= f.Length(); i++)
			{
				f(i) -= -(dpdq(i,1)*XG(1) + dpdq(i,2)*XG(2)); //negative sign: because subtracted from weighted integral
			}

		}
	}
	else if (DOFmode == 4)
	{
		double fact;
		if (locelemind <= npoints1)
		{
			for (int p = 0; p <= maxmoment; p++)
			{
				//nodal points, first body:
				Vector2D lam(-XG(1+2*p), -XG(2+2*p)); //negative sign: needs to be subtracted from f!!!!
				if (p == 0) fact = 1;
				else fact = pow(s_body1(locelemind), p)/length1;
				lam *= fact * weights(locelemind);
				GetBody2D(locelemind).AddNodedPosdqTLambda((int)loccoords(locelemind).X(),lam,f);
			}

		}
		else
		{
			//locelemind == npoints1+1
			//UO() << "s_body1=" << s_body1 << "\n";
			//UO() << "s_body2=" << s_body2 << "\n";

			for (int p = 0; p <= maxmoment; p++)
			{
				//local coordinates, second body:
				Vector2D lam(XG(1+2*p), XG(2+2*p));
				Vector2D lam2;
				for (int k = 1; k <= npoints2; k++)
				{
					if (p == 0) fact = 1;
					else fact = pow(s_body2(k), p)/length2;

					lam2 = fact/length2 * weights(npoints1 + k)*lam;

					GetBody2D(locelemind).GetdPosdqT(loccoords(npoints1 + k),dpdq);

					for (int i=1; i <= f.Length(); i++)
					{
						f(i) -= (dpdq(i,1)*lam2.X()+dpdq(i,2)*lam2.Y());
					}
				}
			}
		}
	}

};

Vector3D IntegralPosConstraint2D::GetRefPosD() const 
{
	Vector2D v;
	if (!DOFmode)
	{
		v = GetBody2D(1).GetPos2DD(loccoords(1));
	}
	else
	{
		v = GetBody2D(1).GetNodePos2DD((int)loccoords(1).X());
	}
	return Vector3D(v.X(),v.Y(),0);
}

void IntegralPosConstraint2D::DrawElement() 
{
	Constraint::DrawElement();
	//


	Vector2D v;
	if (GetMBS()->GetIOption(112))
	{
		mbs->SetColor(GetCol());
		if (!DOFmode)
		{
			for (int i=1; i <= loccoords.Length(); i++)
			{
				v = GetBody2D(i).GetPos2DD(loccoords(i));

				//draw joint node:
				mbs->DrawSphere(GetBody2D(i).ToP3D(v),0.5*draw_dim.X(),4);
			}

		}
		else
		{
			Vector2D vmid(0.,0.);
			Vector2D vmid2(0.,0.);
			if (DOFmode <= 3)
			{
				int l = loccoords.Length();
				if (DOFmode == 3) 
				{
					v = GetBody2D(l).GetPos2DD(loccoords(l));
					//UO() << "vl=" << v << ",bn=" << GetElnum(l) << "\n";
					mbs->DrawSphere(GetBody2D(l).ToP3D(v),0.5*draw_dim.X(),4);
					l -= 1;
				}

				for (int i=1; i <= l; i++)
				{
					//int ni = elements(1);
					v = GetBody2D(i).GetNodePos2DD((int)loccoords(i).X());
					vmid += v;
					vmid2 += weights(i)*v;

					//draw contact node:
					mbs->DrawSphere(GetBody2D(i).ToP3D(v),0.5*draw_dim.X(),4);
				}

				if (l != 0)
				{
					//vmid2 *= 1./l;
					mbs->SetColor(Vector3D(0.1,0.1,0.8));
					//mbs->DrawSphere(GetBody2D(1).ToP3D(vmid2),0.5*draw_dim.X(),6);

					/*vmid *= 1./l;
					mbs->SetColor(Vector3D(0.8,0.1,0.1));
					mbs->DrawSphere(GetBody2D(1).ToP3D(vmid),0.5*draw_dim.X(),6);*/
				}
			}
			else if (DOFmode == 4)
			{
				for (int i=1; i <= npoints1; i++)
				{
					v = GetBody2D(i).GetNodePos2DD((int)loccoords(i).X());
					vmid += weights(i)*v;
					mbs->DrawSphere(GetBody2D(i).ToP3D(v),0.5*draw_dim.X(),5);
				}
				for (int i=1+npoints1; i <= npoints1+npoints2; i++)
				{
					mbs->SetColor(Vector3D(0.4,0.4,0.8));
					v = GetBody2D(i).GetPos2DD(loccoords(i));
					vmid2 += weights(i)*v;
					mbs->DrawSphere(GetBody2D(i).ToP3D(v),0.55*draw_dim.X(),3);
				}
				mbs->SetColor(Vector3D(0.8,0.1,0.1));
				mbs->DrawSphere(GetBody2D(1).ToP3D(vmid),0.6*draw_dim.X(),6);
				mbs->SetColor(Vector3D(0.1,0.1,0.8));
				mbs->DrawSphere(GetBody2D(npoints1+1).ToP3D(vmid2),0.7*draw_dim.X(),3);

			}
		}
	}
};





//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Rotational Spring Damper Actuator Element
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void KardanSDRotActor::GetElementData(ElementDataContainer& edc) 		//fill in all element data
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
	//ed.SetVector3D(spring_stiffness.X(), spring_stiffness.Y(), spring_stiffness.Z(), "Spring_stiffness"); ed.SetToolTipText("Nm/rad");edc.Add(ed);	
	ed.SetVector3D(GetPenaltyStiffness3(1), GetPenaltyStiffness3(2), GetPenaltyStiffness3(3), "Spring_stiffness"); ed.SetToolTipText("Nm/rad");edc.Add(ed);	
}

int KardanSDRotActor::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Element::SetElementData(edc); //do not use Constraint:: !!!!

	double frad = 1;
	if (GetMBS()->GetIOption(120) && !GetMBS()->IsLoadSaveMode()) {frad = 180./MY_PI;}
	Vector3D tmp;
	GetElemDataVector3D(mbs, edc, "Spring_stiffness", tmp, 1);	
	SetPenaltyStiffness3(tmp);
	return rv;
}


Vector3D KardanSDRotActor::ComputeTorque(double t) const   
{
	//compute reacton moment:
	if (IsGroundJoint())
	{
		const Vector3D& lp1 = loccoords(1);
		Matrix3D R1=GetBody3D(1).GetRotMatrix(lp1);

		double b0,b1,b2,b3;
		RotMatToQuaternions(R1,b0,b1,b2,b3);
		Vector3D phi;
		QuaternionsToKardanAngles(b0,b1,b2,b3,phi);

		//phi contains kardan angles (only correct in linearized case!!!)
		if (phi.X() < -MY_PI) phi.X() += 2.*MY_PI;
		if (phi.X() >  MY_PI) phi.X() -= 2.*MY_PI;
		if (phi.Y() < -MY_PI) phi.Y() += 2.*MY_PI;
		if (phi.Y() >  MY_PI) phi.Y() -= 2.*MY_PI;
		if (phi.Z() < -MY_PI) phi.Z() += 2.*MY_PI;
		if (phi.Z() >  MY_PI) phi.Z() -= 2.*MY_PI;

		//
		//?test M^t.G.theta = M^t.A.G_bar.theta = M^t.G.theta, wenn M globales Moment ist
		//return Vector3D(spring_stiffness.X() * phi.X(), spring_stiffness.Y() * phi.Y(), spring_stiffness.Z() * phi.Z());
		return Vector3D(GetPenaltyStiffness3(1) * phi.X(), GetPenaltyStiffness3(2) * phi.Y(), GetPenaltyStiffness3(3) * phi.Z());
	}
	assert(0);
	return Vector3D(0.);
}

void KardanSDRotActor::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	double sign = 1;
	if (locelemind == 2) sign *= -1;
	if (locelemind > elements.Length()) {UO() << "ERROR: inconsistency in KardanSDRotActor::AddElementCqTLambda\n";}

	if (locelemind == 3-groundjoint) return;

	Vector3D rot;
	if (IsGroundJoint())
	{
		sign = -1;
	}

	Vector3D torque(sign*ComputeTorque(t));

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

Vector3D KardanSDRotActor::GetRefPosD()	const 
{
	return GetBody3D(1).GetPosD(loccoords(1));
}

void KardanSDRotActor::DrawElement() 
{
	Constraint::DrawElement();
	//

	double l = GetDrawSizeCylinderLength(); //length of cylinders
	double r = GetDrawSizeScalar(); //radius of cylinders

	mbs->SetColor(GetCol());

	Vector3D p0 = GetBody3D(1).GetPosD(loccoords(1));

	if (IsGroundJoint())
	{
		Vector3D p1x = p0-Vector3D(-0.5*l,0.,0.);
		Vector3D p2x = p0+Vector3D(-0.5*l,0.,0.);
		Vector3D p1y = p0-Vector3D(0.,-0.5*l,0.);
		Vector3D p2y = p0+Vector3D(0.,-0.5*l,0.);
		Vector3D p1z = p0-Vector3D(0.,0.,-0.5*l);
		Vector3D p2z = p0+Vector3D(0.,0.,-0.5*l);
		mbs->DrawZyl(p1x,p2x,r,12);
		mbs->DrawZyl(p1y,p2y,r,12);
		mbs->DrawZyl(p1z,p2z,r,12);
	}
};



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// RigidLinkConstraint (RL)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void RigidLinkConstraint::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Constraint::GetElementData(edc);

	ElementData ed;
	ed.SetDouble(draw_dim.X(), "Draw_size_sphere"); ed.SetToolTipText("General drawing size of constraint"); edc.Add(ed);
	ed.SetDouble(draw_dim.Y(), "Draw_size_cylinder"); ed.SetToolTipText("General drawing size of constraint"); edc.Add(ed);
	ed.SetInt((int)draw_dim.Z(), "Draw_resolution", 0, 1000); ed.SetToolTipText("Number of quads to approximate round geometries, suggested = 16"); edc.Add(ed);
	
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

int RigidLinkConstraint::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Constraint::SetElementData(edc);

	GetElemDataDouble(mbs, edc, "Draw_size_sphere", draw_dim.X(), 0);
	GetElemDataDouble(mbs, edc, "Draw_size_cylinder", draw_dim.Y(), 0);
	int dd;
	if (GetElemDataInt(mbs, edc, "Draw_resolution", dd, 0)) draw_dim.Z() = dd;
	
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
                                  
Vector3D RigidLinkConstraint::ComputeForce(double t)  const//  const war nicht möglich wegen SDActor::fa (forcemode == 4), jetzt wird lokale variable verwendet
{
	Vector3D vp(0.,0.,0.);
	Vector3D v(0.,0.,0.);

	if (IsGroundJoint())
	{
		v = GetBody3D(1).GetPos(loccoords(1)) - p_global;
		v.Normalize();
	}
	else
	{
		v = GetBody3D(1).GetPos(loccoords(1)) - GetBody3D(2).GetPos(loccoords(2));
		v.Normalize();	
	}
	return XG(1) * v;
}

void RigidLinkConstraint::EvalG(Vector& f, double t) 
{
	if (elements.Length()<1 || elements.Length()>2)
	{
		mbs->UO() << "ERROR: RigidLinkConstraint::EvalG, number of elements != 1 or 2\n";
	}
	if (MaxIndex()==3)
	{
		//position constraints: sqrt((p1-p2)^2)-dist) = 0
		if (elements.Length()==1) 
		{
			f(1) = (GetElem(1).GetPos(loccoords(1))-p_global).Norm() - dist;
		}
		else if (elements.Length()==2) 
		{
			f(1) = (GetElem(1).GetPos(loccoords(1))-GetElem(2).GetPos(loccoords(2))).Norm() - dist;
		}
	}
	else
	{ 
		//velocity constraints: d/dt(sqrt((p1-p2)^2)-dist) = 0
		if (elements.Length()==1) 
		{
			f(1) = (GetElem(1).GetPos(loccoords(1))-p_global)*(GetElem(1).GetVel(loccoords(1))); 
		}
		if (elements.Length()==2) 
		{
			f(1) = (GetElem(1).GetPos(loccoords(1))-GetElem(2).GetPos(loccoords(2)))*(GetElem(1).GetVel(loccoords(1))-GetElem(2).GetVel(loccoords(2))); 
		}
	}
};

double RigidLinkConstraint::GetActorForce(double computation_time, int dir) const
{
	if (dir >=1 && dir <= 3) return (ComputeForce(computation_time))(dir);
	else return (ComputeForce(computation_time)).Norm();
}


void RigidLinkConstraint::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	//-C_q^T \lambda = [dC/dx; dC/dy; dC/dphi]^T [\lambda1;\lambda2] 
	//e.g.: C = p_ref+R*p_loc
	double sign = 1;
	if (locelemind == 2) sign = -1;

	Vector3D F3D(ComputeForce(t));

	GetBody3D(locelemind).GetdPosdqT(loccoords(locelemind),dpdq);
	//UO() << "dpdq=" << dpdq;

	for (int i=1; i <= f.Length(); i++)
	{
		f(i) -= sign*(dpdq(i,1)*F3D.X()+dpdq(i,2)*F3D.Y()+dpdq(i,3)*F3D.Z());
	}
};

Vector3D RigidLinkConstraint::GetRefPosD()	const 
{
	return GetBody3D(1).GetPosD(loccoords(1));
}

void RigidLinkConstraint::DrawElement() 
{
	Constraint::DrawElement();
	//


	int res = (int)draw_dim.Z();
	if (res < 2) res = 3;

	if (draw_dim.X() != 0)
	{
		Vector3D p1 = GetBody3D(1).GetPosD(loccoords(1));
		
		mbs->SetColor(GetCol());
		Vector3D p2;
		if (elements.Length()==1)
		{
			p2 = p_global;
		}
		else
		{
			p2 = GetBody3D(2).GetPosD(loccoords(2));
		}

		// cylinder
		if (draw_dim.Y() > 0) mbs->DrawZyl(p1,p2, draw_dim.Y(), res);

		// spheres
		mbs->SetColor(colgrey4);
		mbs->DrawSphere(p1,draw_dim.X(),res);
		mbs->DrawSphere(p2,draw_dim.X(),res);
	}
};







/*
void RevoluteJoint::EvalG(Vector& f, double t) 
{
	//given local vectors:
	//lp1 ... local position vector in body 1
	//lp2 ... local position vector in body 2
	//lr1 ... local axis of rotation in body 1
	//ln2 ... local vector n normal to rotation axis in body 2
	//lt2 ... local second vector t normal to rotation axis in body 2

	if (MaxIndex()==3)
	{
		//compute relative vector between joint position in body 1 and 2:
		Vector3D v = GetBody3D(1).GetPos(lp1) - GetBody3D(2).GetPos(lp2);

		//compute local to global transformation for both bodies
		Matrix3D A1=GetBody3D(1).GetRotMatrix(lp1);
		Matrix3D A2=GetBody3D(2).GetRotMatrix(lp2); 
		Vector rot = A1*lr1; //compute global rotation axis

		//return residuals of constraint equations:
		f(1) = v(1); 
		f(2) = v(2); 
		f(3) = v(3);
		f(4) = rot*(A2*ln2);
		f(5) = rot*(A2*lt2);
	}
	else if (MaxIndex()<=2)
	{
		//compute relative velocity between joint velocity in body 1 and 2:
		Vector3D v = GetBody3D(1).GetVel(lp1) - GetBody3D(2).GetVel(lp2);

		//compute local to global transformation for both bodies
		Matrix3D A1=GetBody3D(1).GetRotMatrix(lp1);   
		Matrix3D A2=GetBody3D(2).GetRotMatrix(lp2);   
		//compute time derivative of transformation
		Matrix3D A1p=GetBody3D(1).GetRotMatrixP(lp1); 
		Matrix3D A2p=GetBody3D(2).GetRotMatrixP(lp2);

		Vector3D rot = A1*lr1; //compute global rotation axis
		Vector3D rotp = A1p*lr1; //compute time derivative of rotation axis
		f(1) = v(1); 
		f(2) = v(2); 
		f(3) = v(3);
		f(4) = rot*(A2p*ln2)+rotp*(A2*ln2);
		f(5) = rot*(A2p*lt2)+rotp*(A2*lt2);
	}
};
*/


/*
//add the generalized constraint forces to vector f
//locelemind is the local index of the constraint
//locelemind==1 ==> first body, locelemind==2 ==> second body
void RevoluteJoint::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	//f has the size of the second order equations of the constrained body
	hmat.SetSize(f.Length(),5);

	double sign = 1;
	if (locelemind==2) sign = -1;

	GetBody3D(locelemind).GetdPosdqT(loccoords(locelemind),dpdq);
	hmat.SetSubmatrix(dpdq,1,1, sign);

	if (locelemind==1)
	{
		//compute global vectors n and t, 
		Vector3D vn2_glob = GetBody3D(2).GetRotMatrix(lp2)*ln2;
		Vector3D vt2_glob = GetBody3D(2).GetRotMatrix(lp2)*lt2;
		
		//compute derivative of lr1 with respect to element DOF at position lp1
		GetBody3D(1).GetdRotvdqT(lr1,lp1,dpdq);

		Mult(dpdq,vn2_glob,hvec);
		hmat.SetColVec(hvec,4);

		Mult(dpdq,vt2_glob,hvec);
		hmat.SetColVec(hvec,5);
	}
	else
	{
		//compute global rotation axis
		Vector3D vrglob = GetBody3D(1).GetRotMatrix(lp1)*lr1;

		//compute derivative of ln2 with respect to element DOF at position lp2
		GetBody3D(2).GetdRotvdqT(ln2,lp2,dpdq);
		Mult(dpdq,vrglob,hvec);
		hmat.SetColVec(hvec,4);

		//compute derivative of lt2 with respect to element DOF at position lp2
		GetBody3D(2).GetdRotvdqT(lt2,lp2,dpdq);
		Mult(dpdq,vrglob,hvec);
		hmat.SetColVec(hvec,5);
	}

	//Vector lambda(5) contain the lagrange multipliers for the revolute joint
	for (int i=1; i <= f.Length(); i++)
		for (int j=1; j <= 5; j++)
			f(i) -= hmat(i,j)*lambda(j);
};
*/









//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Surface
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Surface::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Constraint::GetElementData(edc);

	ElementData ed;
	//ed.SetDouble(draw_dim.X(), "Draw_size"); ed.SetToolTipText("General drawing size of constraint"); edc.Add(ed);
	//ed.SetInt((int)draw_dim.Z(), "Draw_resolution", 0, 1000); ed.SetToolTipText("Number of quads to approximate round geometries, suggested = 16"); edc.Add(ed);

	//ed.SetInt(loccoords(1), "Element_Coordinate1"); ed.SetToolTipText("Local coordinate of element 1 to be constrained"); edc.Add(ed);
	//if (loccoords.Length() == 2)
	//{
	//	ed.SetInt(loccoords(2), "Element_Coordinate2"); ed.SetToolTipText("Local coordinate of element 2 to be constrained"); edc.Add(ed);
	//}

	//if(GetPenalty())
	//{
	//	ed.SetDouble(spring_stiffness, "Penalty_stiffness"); ed.SetToolTipText("Penalty stiffness of coordinate constraint"); edc.Add(ed);
	//}
}

int Surface::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Constraint::SetElementData(edc);

	//GetElemDataDouble(mbs, edc, "Draw_size", draw_dim.X(), 0);
	//int dd;
	//if (GetElemDataInt(mbs, edc, "Draw_resolution", dd, 0)) draw_dim.Z() = dd;

	//GetElemDataInt(mbs, edc, "Element_Coordinate1", loccoords(1), 1);
	//if (loccoords(1) < 1 || loccoords(1) > GetElem(1).SS())
	//{
	//	loccoords(1) = 1;
	//	GetMBS()->EDCError("Illegal element coordinate in CoordConstraint");
	//	rv = 0;
	//}

	//if (loccoords.Length() == 2)
	//{
	//	GetElemDataInt(mbs, edc, "Element_Coordinate2", loccoords(2), 1);

	//	if (loccoords(2) < 1 || loccoords(2) > GetElem(2).SS())
	//	{
	//		loccoords(2) = 1;
	//		GetMBS()->EDCError("Illegal element coordinate in CoordConstraint");
	//		rv = 0;
	//	}

	//}

	//if(GetPenalty())
	//{
	//	GetElemDataDouble(mbs, edc, "Penalty_stiffness", spring_stiffness); 		
	//}

	return rv;
}
//
//void CoordConstraint::EvalG(Vector& f, double t) 
//{
//	if(penalty) return; //(RL)
//
//	if (elements.Length()<1 || elements.Length()>2)
//	{
//		mbs->uout << "ERROR: CoordConstraint::EvalG, number of elements != 1 or 2\n";
//	}
//	if (MaxIndex()==3)
//	{
//		if (elements.Length()>=1) 
//		{
//			f(1) = GetElem(1).XG(loccoords(1))-GetElem(1).GetXInit()(loccoords(1));
//		}
//		if (elements.Length()==2) 
//		{
//			f(1)-= GetElem(2).XG(loccoords(2))-GetElem(2).GetXInit()(loccoords(2));
//		}
//	}
//	else
//	{ //velocity constraints:
//		if (elements.Length()>=1) 
//		{
//			f(1) = GetElem(1).XGP(loccoords(1))-GetElem(1).GetXInit()(loccoords(1)+GetElem(1).SOS());
//		}
//		if (elements.Length()==2) 
//		{
//			f(1)-= GetElem(2).XGP(loccoords(2))-GetElem(2).GetXInit()(loccoords(2)+GetElem(2).SOS());
//		}
//	}
//};
//
////To be replaced in derived class
//void CoordConstraint::AddElementCqTLambda(double t, int locelemind, Vector& f)


void Surface::UpdateSurfNodeNrList(int en, int fn)
{
	//get global node numbers
	FiniteElement3D& fel = dynamic_cast<FiniteElement3D&> (GetMBS()->GetElement(en));

	FEFace felf = fel.GetFace(fn);
	int locnodenr=0, globnodenr=0;

	for(int i=1; i < felf.NFaceNodes()+1; i++)
	{
		locnodenr=felf.Node(i);
		globnodenr=fel.NodeNum(locnodenr);
		//check if already in list
		if(surf_node_nrs_nodouble.Find(globnodenr)==0)
			surf_node_nrs_nodouble.Add(globnodenr); //add global node numbers to list
	}
};

void Surface::DrawElement() 
{
	if (draw_dim.X() == 0) return;

	mbs->SetColor(GetCol());

};



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Spherical: MultiNodalSphericalJoint
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: MultiNodalSphericalJoint (Constraint)
// short description: 
// available formulations: Penalty Formulation, Lagrange Formulation
// type of constraint equation: nonlinear constraint equations
// index formulation available: 
// ground constraint: yes
// element-element constraint: yes
// development status: in progress
// long description:	
//								
// class variables:
//  
//**end(ued)**
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


int MultiNodalSphericalJoint::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = Constraint::CheckConsistency(errorstr);
	if (rv) return rv;
	if(!UsePenaltyFormulation())
	{
			errorstr = mystr("ERROR: MultiNodalSphericalJoint: only penalty formulation is possible\n");
			rv = 1;
	}
	
	return rv;
}

void MultiNodalSphericalJoint::ElementDefaultConstructorInitialization()
{
	SetPenaltyFormulation(1);
	SetUseLocalCoordinateSystem(0);
	loccoords.Set2(0,0);
	nodes.SetLen(0);
	nodes1.Set1(1);
	nodes2.SetLen(0);
	velocity_constraint = 0;
	SetDampingCoeff(0.0);
	auto_comp_ground = 0;
	stiffness_in_joint_local_frame = 0;
	elements.Set2(0,0);
	elementname = GetElementSpec();
	draw_dim.Y() = GetDrawSizeScalar();
	draw_dim.Z() = 1;
}

void MultiNodalSphericalJoint::CopyFrom(const Element& e)
{
	BasePointJoint::CopyFrom(e);
	const MultiNodalSphericalJoint& ce = (const MultiNodalSphericalJoint&)e;

	nodes1 = ce.nodes1;
	nodes2 = ce.nodes2;
}

// average position of all nodes of the specified kinematic pair
Vector3D MultiNodalSphericalJoint::GetAveragePosition(int kinpair) const
{
	const TArray<int> *nodes_tmp;
	if(kinpair==1)	{ nodes_tmp = &nodes1;}
	else			
	{ 
		if(nodes2.Length())
		{
			nodes_tmp = &nodes2;
		}
		else	// ground
		{
			return loccoords(2);
		}
	}

	Vector3D pos(0);
	for(int n=1; n<=nodes_tmp->Length(); n++)
	{
		pos +=(GetMBS()->GetNode(nodes_tmp->Get(n))).GetPos();
	}
	double fact = 1./(nodes_tmp->Length());
	pos *= fact;
	return pos;
}

Vector3D MultiNodalSphericalJoint::GetAverageDrawPosition(int kinpair)
{
	const TArray<int> *nodes_tmp;
	if(kinpair==1)	{ nodes_tmp = &nodes1;}
	else			
	{ 
		if(nodes2.Length())
		{
			nodes_tmp = &nodes2;
		}
		else	// ground
		{
			return loccoords(2);
		}
	}

	Vector3D pos(0);
	for(int n=1; n<=nodes_tmp->Length(); n++)
	{
		pos +=(GetMBS()->GetNode(nodes_tmp->Get(n))).GetPosD();
	}
	double fact = 1./(nodes_tmp->Length());
	pos *= fact;
	return pos;
}

// average velocity of all nodes of the specified kinematic pair
Vector3D MultiNodalSphericalJoint::GetAverageVelocity(int kinpair) const
{
	const TArray<int> *nodes_tmp;
	if(kinpair==1)	{ nodes_tmp = &nodes1;}
	else			
	{
		if(nodes2.Length())
		{
			nodes_tmp = &nodes2;
		}
		else	//ground
		{	
			return Vector3D(0);
		}
	}

	Vector3D vel(0);
	for(int n=1; n<=nodes_tmp->Length(); n++)	
	{
		vel +=(GetMBS()->GetNode(nodes_tmp->Get(n))).GetVel();
	}
	double fact = 1./(nodes_tmp->Length());
	vel *= fact;
	return vel;
}

//position of the node i of the specified kinematic pair
Vector3D MultiNodalSphericalJoint::GetDiscretePosition(int kinpair, int i) const
{
	Vector3D pos;
	//Body3D body(GetMBS());

	if (kinpair <= NKinPairs())					// not ground
	{
		if(kinpair == 1)	{	pos = (GetMBS()->GetNode(nodes1(i))).GetPos();}
		else							{	pos = (GetMBS()->GetNode(nodes2(i))).GetPos();}
	}
	else if (kinpair==1)								// pos1 is ground	--> not possible
	{
		GetMBS()->UO(UO_LVL_err) << "ERROR: MultiNodalSphericalJoint: no element number defined for body 1";
		return Vector3D(0.0);
	}
	else if (kinpair==2)								// pos2 is ground	
	{
		pos = loccoords(2);					// if SetPos2ToGlobalCoord is used, then loccords(2) contains the global(!) coords of ground
		if(displ.Length() != 0)
		{
			double t_glob = Constraint::GetGlobalTime();
			for(int j=1; j<=displ.Length(); j++)
			{
				if(displ(j) != NULL)
				{
					pos(j) = pos(j)+ (displ(j)->Evaluate(t_glob));
				}
			}
		}			
	}
	return pos;
}

//velocity of the node i of the specified kinematic pair
Vector3D MultiNodalSphericalJoint::GetDiscreteVelocity(int kinpair, int i) const
{
	Vector3D vel;
	if (kinpair <= NKinPairs())					// not ground
	{
		if(kinpair == 1)	{ vel = (GetMBS()->GetNode(nodes1(i))).GetVel();}
		else							{ vel = (GetMBS()->GetNode(nodes2(i))).GetVel();}
	}
	else if (kinpair==1)								// pos1 is ground	--> not possible
	{
    GetMBS()->UO(UO_LVL_err) << "ERROR: MultiNodalSphericalJoint: no element number defined for body 1";
		return Vector3D(0.0);
	}
	else if (kinpair==2)								// pos2 is ground	
	{
		return Vector3D(0.0);
	}
	return vel;
}
void MultiNodalSphericalJoint::SetPosToGlobalNodes(int kinpair,  TArray<int> glob_node_nrs)
{
	// for global nodes the element number '0' is inserted
	// if i > elements.length resize is done automatically by TArray
	if(kinpair == 1)	
	{
		nodes1 = glob_node_nrs;	
		elements(1) = 0;
	}
	else
	{
		nodes2 = glob_node_nrs;	
		elements(2) = 0;
	}
}


void MultiNodalSphericalJoint::SetPos2ToGlobalCoord(Vector3D ground)
{
	loccoords(2) = ground;
	elements(2) = 0;
	nodes2.SetLen(0);
}

//void MultiNodalSphericalJoint::GetElementData(ElementDataContainer& edc) 		//fill in all element data
//{
//	Element::GetElementData(edc);
//}
//
//int MultiNodalSphericalJoint::SetElementData(const ElementDataContainer& edc) //set element data according to ElementDataContainer
//{
//	return Element::SetElementData(edc);
//}


Vector3D MultiNodalSphericalJoint::ComputeForce(double t) const
{
	Vector3D u; // displacement
	Vector3D v; // velocity
	Vector3D f; // resulting force
	Vector3D k = this->GetPenaltyStiffness3(t); // spring stiffness
	Matrix3D rot; // reference rotation
	Matrix3D k_mat;

	if (UsePenaltyFormulation())
	{
		u = GetPos1() - GetPos2();
		if(UseDamping()) {	v = GetVel1() - GetVel2();}

		k_mat = Matrix3D(k.X(),k.Y(),k.Z());

		// F = A*K*A'*u
		if(UseLocalCoordinateSystem()||stiffness_in_joint_local_frame)
		{
			rot = GetRotMati();
			rot.Transpose();
			k_mat = k_mat*rot;
			rot = GetRotMati();
			k_mat = rot*k_mat;
		}

		f = k_mat*u;	

		if(UseDamping()) { f += GetDampingCoeff() * v; }

		return f;
	}
	else
	{
		GetMBS()->UO(UO_LVL_err) << "ERROR: MultiNodalSphericalJoint: ComputeForce just implemented for PenaltyFormulation";
		return Vector3D(0.0);
	}

};


void MultiNodalSphericalJoint::EvalF2(Vector& f, double t)
{
	//f = [f1_1 f1_2 .. f1_n1 f2_1 f2_2 .. f2_n2], where f1_i is the residual vector of the node i of constraint element1 
	//$ YV 2012.10.10 - modified the computation of the factor; all nodes are assumed to be 3D, so that there are 3 force components per node
	if(!UsePenaltyFormulation())
		return;	// no penalty method --> Lagrange multiplier --> no EvalF2
	Vector3D force = ComputeForce(t);
	int offset;
	int length;
	for (int i=1; i <= NKinPairs(); i++)
	{
		if (i==2) 
		{
			length = nodes2.Length();
			force *= -1. / length;
			offset = 3 * nodes1.Length();
		}
		else
		{
			length = nodes1.Length();
			force *= 1. / length;
			offset = 0;
		}
		for(int n = 0; n <= length - 1; n++)
			for(int k = 1; k <= 3; k++)
				f(3 * n + k + offset) -= force(k);
	}	
};

void MultiNodalSphericalJoint::LinkToElementsPenalty()
{
	// add all SOS dofs from the elements
	//Position(first SOS) 
	for (int k=1; k <= NKinPairs(); k++)
	{
		TArray<int> *nodes_tmp;
		if(k==1)	{ nodes_tmp = &nodes1;}
		else			{ nodes_tmp = &nodes2;}
		for(int n=1; n<=nodes_tmp->Length(); n++)
		{
			const Node& node = GetMBS()->GetNode(nodes_tmp->Get(n));
			for (int i=1; i <= node.Dim(); i++)
			{
				AddLTG(node.Get(i));
			}
		}
	}
	//and Velocity (second SOS):
	for (int k=1; k <= NKinPairs(); k++)
	{
		TArray<int> *nodes_tmp;
		if(k==1)	{ nodes_tmp = &nodes1;}
		else			{ nodes_tmp = &nodes2;}
		for(int n=1; n<=nodes_tmp->Length(); n++)
		{
			const Node& node = GetMBS()->GetNode(nodes_tmp->Get(n));
			for (int i=1; i <= node.Dim(); i++)
			{
				AddLTG(node.Get(i+node.SOS()));
			}
		}
	}
}


int MultiNodalSphericalJoint::SOS() const 
{
	int nsos = 3*nodes1.Length()+3*nodes2.Length();
	return nsos;

};  // explicit size, number of constrained dofs


void MultiNodalSphericalJoint::DrawElement() 
// 3 options available depending on the entries in draw_dim
// if only draw_dim.X() is set --> a sphere in the center of the nodes is drawn
// X and Y are set --> additionally a sphere for every node in nodelist is drawn
// X, Y and Z are set --> additionally lines from the center sphere to the nodes is drawn
{
	Constraint::DrawElement();

	if (GetDrawSizeScalar() == 0) return;		// do not draw MultiNodalSphericalJoint at all

	Vector3D dir(0.,0.,0.);
	Vector3D p(0.,0.,0.);
	Vector3D center(0.,0.,0.);
	int res = 8;

	for (int j=1; j <= NKinPairs(); j++)
	{
		//mbs->SetColor(col);
		if (j == 1) {	mbs->SetColor(Vector3D(0.,0.,0.8));}
		if (j == 2) {	mbs->SetColor(Vector3D(0.8,0.1,0.1)); res = 7;}

		// draw the sphere in the center
		center = GetAverageDrawPosition(j);
		mbs->DrawSphere(center,GetDrawSizeScalar(),res);

		if(GetDrawSizeCircleSpheres()!=0)		// draw a sphere for every node in nodelist
		{
			const TArray<int> *nodes_tmp;
			if(j==1)	{ nodes_tmp = &nodes1;}
			else			{ nodes_tmp = &nodes2;}

			for (int i=1; i<=nodes_tmp->Length(); i++)
			{	
				p = (GetMBS()->GetNode(nodes_tmp->Get(i))).GetPosD();
				mbs->DrawSphere(p,GetDrawSizeCircleSpheres(),res);
				if(draw_dim.Z()!=0)		// draw lines from the center sphere to the nodes
				{
					mbs->MyDrawLine(center,p,GetDrawSizeLineThicknessStar());
				}
			}
		}
	}
};







//$ DR 2011-04:[ MultiNodalConstraint added
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// MultiNodalConstraint
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//initialize the offset between the geometric center and the center computed using the nodelist
	//constrain 2 elements: center1 = center2 + center_offset
	//constrain 1 element to ground: center1 = center_offset
void MultiNodalConstraint::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Constraint::GetElementData(edc);

	ElementData ed;
	ed.SetDouble(GetDrawSizeScalar(), "Draw_size"); ed.SetToolTipText("size of the one sphere in the center"); edc.Add(ed);
	ed.SetDouble(GetDrawSizeCircleSpheres(), "DrawSizeCircleSpheres"); ed.SetToolTipText("size of the spheres"); edc.Add(ed);
	ed.SetInt(GetDrawSizeLineThicknessStar(), "LineThicknessStar"); ed.SetToolTipText("line thickness of the lines to the center sphere"); edc.Add(ed);

	int i;
	for (i=1; i <= loccoords.Length(); i++)
	{
		ed.SetInt(loccoords(i), mystr("Node_Coordinate")+mystr(i)); ed.SetToolTipText("Local coordinate of node to be constrained"); edc.Add(ed);
	}

	if (UsePenaltyFormulation()) //$ DR 2011-04-20: old code: if(GetPenalty())
	{
		ed.SetDouble(GetPenaltyStiffness(), "Penalty_stiffness"); ed.SetToolTipText("Penalty stiffness of coordinate constraint"); edc.Add(ed);
		ed.SetDouble(damping_coeff, "Penalty_damping"); ed.SetToolTipText("Penalty damping of coordinate constraint"); edc.Add(ed);
	}
}

int MultiNodalConstraint::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Constraint::SetElementData(edc);
	double tmp;
	GetElemDataDouble(mbs, edc, "Draw_size", tmp, 0);
	SetDrawSizeScalar(tmp);
	GetElemDataDouble(mbs, edc, "DrawSizeCircleSpheres", tmp, 0);
	SetDrawSizeCircleSpheres(tmp);
	int th;
	GetElemDataInt(mbs, edc, "LineThicknessStar", th, 0);
	SetDrawSizeLineThicknessStar(th);

	int i;

	for (i=1; i <= loccoords.Length(); i++)
	{
		GetElemDataInt(mbs, edc, mystr("Node_Coordinate")+mystr(i), loccoords(i), 1);
		if (loccoords(i) < 0 )
		{
			loccoords(i) = 0;
			GetMBS()->EDCError("Illegal element coordinate in MultiNodalConstraint");
			rv = 0;
		}
	}

	if (UsePenaltyFormulation()) //$ DR 2011-04-20: old code: if(GetPenalty())
	{
		GetElemDataDouble(mbs, edc, "Penalty_stiffness", tmp);
		SetPenaltyStiffness(tmp);
	  GetElemDataDouble(mbs, edc, "Penalty_damping", damping_coeff);
	}

	return rv;
}


void MultiNodalConstraint::Initialize()
	{
		center_offset = ComputeCenterOfNodes(1,nodes_per_element[0]);
		mbs->UO(UO_LVL_dbg1) << "MultiNodalConstraint::Initialize: center of element 1 =" << center_offset << "\n";					
		if(NodeNum(2,1)!=0) // constrain 2 elements 
		{
			mbs->UO(UO_LVL_dbg1) << "MultiNodalConstraint::Initialize: center of element 2 =" << ComputeCenterOfNodes(2,nodes_per_element[1]) << "\n";					
			center_offset -= ComputeCenterOfNodes(2,nodes_per_element[1]);			
		}
		mbs->UO(UO_LVL_all) << "MultiNodalConstraint::Initialize: center_offset =" << center_offset << "\n";					


		if(GetUseConstantNodeDpdq()){				// compute constant average matrix dpdq, e.g. when GCMS is used
			for (int j=1; j <= NE(); j++) 
			{
				Matrix& mat = constant_node_dpdqT[j-1];
				for (int i=1; i <= NodesPerElementLength(j); i++)
				{
					GetBody3D(j).GetNodedPosdqT(NodeNum(j,i), dpdq);
					if (i == 1) {	mat = GetWeight(j,i)*dpdq;	}
					else {mat += GetWeight(j,i)*dpdq;}
				}
				//mbs->UO(UO_LVL_all) << "constant_node_dpdqT[" << j-1	<<	"]= " <<constant_node_dpdqT[j-1] << "\n";	
				Matrix dpdq = constant_node_dpdqT[j-1];							// get the dpdqT matrix
				dpdq.TpYs();																				// Transpose, result: dpdq matrix
				average_node_dpdq[j-1]=dpdq;												// save dpdq matrix
			}
		}
	}

// return position of node (nodenumbers of constrain are stored in nodes_per_element)
Vector3D MultiNodalConstraint::GetElementNodePos(int elem_nr, int node_nr, int flagD) const
{
	if (!flagD)
	{
		return GetBody3D(elem_nr).GetNodePos(NodeNum(elem_nr, node_nr));
	}
	else
	{
		return GetBody3D(elem_nr).GetNodePosD(NodeNum(elem_nr, node_nr));
	}
}

// return velocity of node (nodenumbers of constrain are stored in nodes_per_element)
Vector3D MultiNodalConstraint::GetElementNodeVel(int elem_nr, int node_nr, int flagD) const
{
	if (!flagD)
	{
		return GetBody3D(elem_nr).GetNodeVel(NodeNum(elem_nr, node_nr));
	}
	else
	{
		return GetBody3D(elem_nr).GetNodeVelD(NodeNum(elem_nr, node_nr));
	}
}

// get the distance of the 2 centers of nodelists
// without displacement (initial positions): center1 = center2 + center_offset
Vector3D MultiNodalConstraint::GetDistance() const 
{
	Vector3D dist;
	if(GetUseConstantNodeDpdq()){															//e.g. when GCMS is used
		Vector center, dist_vec;
		for (int i=1; i <= NE(); i++) 
		{
			const Matrix& const_dpdq = average_node_dpdq[i-1];		
			Vector q;																							// get vector q
			q.SetLen(GetBody3D(i).SOS());
			for(int j=1; j <=GetBody3D(i).SOS(); j++)
			{
				q(j) = GetBody3D(i).XG(j);
			}
			Mult(const_dpdq,q,center);															// u = N*q
			if(i==1) dist_vec = center;															// dist = r1 - (r2 + offset) = u1 - u2
			else dist_vec -= center;
		}
		dist = Vector3D(dist_vec(1),dist_vec(2),dist_vec(3));
	}
	else
	{
		Vector3D center[2];
		for (int i=1; i <= NE(); i++) 
		{
			center[i-1] = ComputeCenterOfNodes(i,nodes_per_element[i-1]);
		}
		dist = center[0] - (center[1] + center_offset);						// dist = r1 - (r2 + offset)
	}
	return dist;
}

// get the relatice velocity of the 2 centers of nodelists
Vector3D MultiNodalConstraint::GetRelVel() const 
{
	Vector3D relvel;
	if(GetUseConstantNodeDpdq()){																//e.g. when GCMS is used
		Vector vel, vel_vec;
		for (int i=1; i <= NE(); i++) 
		{
			const Matrix& const_dpdq = average_node_dpdq[i-1];		
			Vector qp;																							
			qp.SetLen(GetBody3D(i).SOS());													// get vector qp
			for(int j=1; j <=GetBody3D(i).SOS(); j++)
			{
				qp(j) = GetBody3D(i).XGP(j);
			}
			Mult(const_dpdq,qp,vel);																// v = N*qp													
			if(i==1) vel_vec = vel;															
			else vel_vec -= vel;
		}
		relvel = Vector3D(vel_vec(1),vel_vec(2),vel_vec(3));
	}
	else
	{
		Vector3D vel[2];
		for (int i=1; i <= NE(); i++) 
		{
			vel[i-1] = ComputeCenterOfNodesVel(i,nodes_per_element[i-1]);
		}
		relvel = vel[0] - vel[1];
	}
	return relvel;
} 

void MultiNodalConstraint::EvalF2(Vector& f, double t)
{
	// always assume that it is local nodal constraint
	// just implemented for penalty
	// without displacement (initial positions): center1 = center2 + center_offset

	//TMStartTimer(28);
	Vector3D dist = GetDistance();		
	Vector3D relvel = GetRelVel();		
	//TMStopTimer(28);
	double tempS, tempD;
	int lc;
	
	if(GetUseConstantNodeDpdq()){															//e.g. when GCMS is used
		for (int j=1; j <= NE(); j++) 
		{
			int sosoff = 0;
			double flag = 1;
			if (j==2) 
			{
				flag = -1;
				sosoff = GetBody3D(1).SOS();
			}
			Matrix& const_dpdqT = constant_node_dpdqT[j-1];		

			tempS = flag *  GetPenaltyStiffness();
			tempD = flag *  damping_coeff;
			for (int k=1; k<=const_dpdqT.Getrows(); k++)		
			{
				lc=loccoords(j);
				f(k + sosoff) -= tempS*dist(lc)*const_dpdqT(k,lc);
				f(k + sosoff) -= tempD*relvel(lc)*const_dpdqT(k,lc);

				if(loccoords.Length()>3)
				{
					lc=loccoords(j+2);
					f(k + sosoff) -= tempS*dist(lc)*const_dpdqT(k,lc);
					f(k + sosoff) -= tempD*relvel(lc)*const_dpdqT(k,lc);

					if(loccoords.Length()>5)
					{
						lc=loccoords(j+4);
						f(k + sosoff) -= tempS*dist(lc)*const_dpdqT(k,lc);
						f(k + sosoff) -= tempD*relvel(lc)*const_dpdqT(k,lc);
					}
				} 
			}
		}
	}
	else
	{
		for (int j=1; j <= NE(); j++) 
		{
			for (int i=1; i <= NodesPerElementLength(j); i++)
			{

				int sosoff = 0;
				GetBody3D(j).GetNodedPosdqT(NodeNum(j,i), dpdq);		

				double flag = 1;
				if (j==2) 
				{
					flag = -1;
					sosoff = GetBody3D(1).SOS();
				}

				tempS = flag * GetWeight(j,i) * GetPenaltyStiffness();
				tempD = flag * GetWeight(j,i) * damping_coeff;
				for (int k=1; k<=dpdq.Getrows(); k++)		
				{
					lc=loccoords(j);
					f(k + sosoff) -= tempS*dist(lc)*dpdq(k,lc);
					f(k + sosoff) -= tempD*relvel(lc)*dpdq(k,lc);

					if(loccoords.Length()>3)
					{
						lc=loccoords(j+2);
						f(k + sosoff) -= tempS*dist(lc)*dpdq(k,lc);
						f(k + sosoff) -= tempD*relvel(lc)*dpdq(k,lc);

						if(loccoords.Length()>5)
						{
							lc=loccoords(j+4);
							f(k + sosoff) -= tempS*dist(lc)*dpdq(k,lc);
							f(k + sosoff) -= tempD*relvel(lc)*dpdq(k,lc);
						}
					} 
				}
			}
		}
	}
}

void MultiNodalConstraint::DrawElement() 
// draws MultiNodalConstraint
// 3 options available depending on the entries in draw_dim
// if only draw_dim.X() is set --> a sphere in the center of the nodes is drawn
// X and Y are set --> additionally a sphere for every node in nodelist is drawn
// X, Y and Z are set --> additionally lines from the center sphere to the nodes is drawn
{
	Constraint::DrawElement();
	//

	if (GetDrawSizeScalar() == 0) return;		// do not draw MultiNodalConstraint at all

	Vector3D dir(0.,0.,0.);
	Vector3D p(0.,0.,0.);
	Vector3D center(0.,0.,0.);
	int res = 8;

	for (int j=1; j <= NE(); j++)
	{
		//mbs->SetColor(col);
		if (j == 1) {	mbs->SetColor(Vector3D(0.,0.,0.8));}
		if (j == 2) {	mbs->SetColor(Vector3D(0.8,0.1,0.1)); res = 7;}

		// draw the sphere in the center
		if(nodes_per_element[j-1](1)!=0) center = ComputeCenterOfNodes(j,nodes_per_element[j-1], 1);
		else center = center_offset;	// if it is a constrait to ground use the initial center_offset as position
		mbs->DrawSphere(center,GetDrawSizeScalar(),res);

		if(GetDrawSizeCircleSpheres()!=0)		// draw a sphere for every node in nodelist
		{
			for (int i=1; i<=nodes_per_element[j-1].Length(); i++)
			{	
				p = GetBody3D(j).GetNodePosD(NodeNum(j,i));
				mbs->DrawSphere(p,GetDrawSizeCircleSpheres(),res);
				if(draw_dim.Z()!=0)		// draw lines from the center sphere to the nodes
				{
					mbs->MyDrawLine(center,p,GetDrawSizeLineThicknessStar());
				}
			}
		}
	}
};
//$ DR 2011-04:] MultiNodalConstraint added


// ##############################################################################################
//$ DR 2012: FrictionConstraintUniDir added

void FrictionConstraintUniDir::SetDefaultValues()
	{
		SetPenaltyFormulation(1);
		SetUseLocalCoordinateSystem(0);
		SetUseConstantNormalForce(1);
		SetVelocityTolerance(1e-5);

		Vector datainit(DataS());
		datainit.SetAll(0.);
		datainit(1)=1.;
		SetDataInit(datainit);
		Fn0 = 0.;
		keep_sliding = 0;
	}

void FrictionConstraintUniDir::SetFrictionConstraintUniDir_LocalNode_to_LocalNode(int elem1, int node1, int elem2, int node2, Vector3D normal_vec, double stiffness, double damping, double friction_coeff_st,  double friction_coeff_kin, Vector3D axial_dir)
	{
		SetPos1ToLocalNode(elem1, node1);
		SetPos2ToLocalNode(elem2, node2);
		SetAxialDir(axial_dir);
		SetPenaltyStiffness(stiffness);
		SetFrictionCoeff_st(friction_coeff_st);
		SetFrictionCoeff_kin(friction_coeff_kin);
		SetDampingCoeff(damping);

		SetDefaultValues();
	}
void FrictionConstraintUniDir::SetFrictionConstraintUniDir_LocalPos_to_LocalPos(int elem1, Vector3D lc1, int elem2, Vector3D lc2, Vector3D normal_vec, double stiffness, double damping, double friction_coeff_st,  double friction_coeff_kin, Vector3D axial_dir)
	{
		SetPos1ToLocalCoord(elem1, lc1);
		SetPos2ToLocalCoord(elem2, lc2);
		SetAxialDir(axial_dir);
		SetPenaltyStiffness(stiffness);
		SetFrictionCoeff_st(friction_coeff_st);
		SetFrictionCoeff_kin(friction_coeff_kin);
		SetDampingCoeff(damping);

		SetDefaultValues();
	}

// Get force that acts during stick-phase. Returns force in axial direction.
Vector3D FrictionConstraintUniDir::GetAxialForce(Vector3D u, double stiffness) const
	{
		Vector3D vec_ax = GetAxialDir();
		Vector3D f = (u*vec_ax)*stiffness*vec_ax; // |vec_ax| = 1
		Vector3D v = 	SlidingVelAx()*vec_ax;
		f += v*GetDampingCoeff();
		return f;
	}

// NormalForce Fn is positve for pressure!
double FrictionConstraintUniDir::GetNormalForce(Vector3D u, double stiffness)const
	{
		return Fn0;
	}

	// Friction Force when body is sliding: F = µk*Fn
Vector3D FrictionConstraintUniDir::GetSlidingFrictionForce(double coef_kin)const
	{
		Vector3D dir = 	SlidingVelAx()*GetAxialDir();
		dir.Normalize();
		return dir*coef_kin*(Fn0); 	
	}

	// maximal Force at which sliding starts: Fmax = µs*Fn
double FrictionConstraintUniDir::GetTraction(double t, Vector3D u) const
{
		return 	GetFrictionCoeff_st()*(Fn0);	
}

void FrictionConstraintUniDir::PostprocessingStep()
	{
		double lastXData = GetMBS()->GetLastDataVector()(LTGdata(1));

		//Vector3D v = GetVel1() - GetVel2();

		if (IsSticking() == 1 && lastXData != 1) //switch from sliding to stick ==> remember actual sticking position
		{
			Vector3D u = GetPos1() - GetPos2();
			Vector3D vec_ax = GetAxialDir();
			StickingPosAx(u*vec_ax);				// origin of spring is last sticking point
			SlidingVelAxOLD(SlidingVelAx());							// last sliding velocity
			if(mbs->UO(UO_LVL_dbg1).PrintMsg())
			{	
				GetMBS()->UO() << "change slide to STICK: elementNum = " << GetOwnNum() <<", t = " << GetMBS()->GetTime() <<", u_ax = " << StickingPosAx() <<", u_tan = " << StickingPosTan() << "\n"; 
			}
		}
		else if (IsSticking() == 0) //during sliding ==> remember last sliding velocity
		{
			//XData(2) = 0;					// origin of spring is reset
			//XData(3) = 0;
			SlidingVelAxOLD(SlidingVelAx());		// last sliding velocity
			if(lastXData && (mbs->UO(UO_LVL_dbg1).PrintMsg()) )
			{	
				GetMBS()->UO() << "change stick to SLIDE: elementNum = " << GetOwnNum() <<", t = " << GetMBS()->GetTime() <<", v_ax= " << SlidingVelAx() <<", v_tan= " << SlidingVelTan() << "\n";
			}
		}
	}

double FrictionConstraintUniDir::PostNewtonStep(double t)
	{
		//return 0; //$!DR Hack
		double error = 0.;

		Vector3D u; // displacement
		Vector3D v; // velocity
		Vector3D vec_ax = GetAxialDir();

		double stiffness = GetPenaltyStiffness();

		u = GetPos1() - GetPos2();
		v = GetVel1() - GetVel2(); 

		//double vax = v*vec_ax;
		SlidingVelAx(v*vec_ax);

		u = u - StickingPosAx()*vec_ax;	// origin of spring is last sticking point

		if(IsSticking())	// stick
		{
			double Fa = GetAxialForce(u,stiffness).Norm();
			double threshhold = GetTraction(t,u);

			if(fabs(Fa) <= threshhold) // keep sticking
			{
				error = 0.;
			}
			else									// start sliding
			{
				error = fabs(Fa-threshhold);
				if (threshhold!=0) {error /= threshhold;}
				//GetMBS()->UO() << "change stick to SLIDE: t = " << t <<", Fa = " << Fa <<" > threshhold = " << threshhold <<", error = " << error << "\n";
				IsSticking(0);
			}
		}
		else					// sliding
		{
			if(keep_sliding)				// no sticking again
			{
				error = 0;
			}
			else
			{
				double tol = GetVelocityTolerance(); // if v < tol than it is assumed that v = 0	

				if( SlidingVelAx()*SlidingVelAx() < tol*tol) // start sticking when new velocity is very small (assumed to be zero)
				{
					IsSticking(1);
					error = fabs((SlidingVelAx()*SlidingVelAx())/tol);
					//GetMBS()->UO() << "v_new = " << SlidingVelAx() << " < " << tol<< "\n";
				}
				else
				{
					if(SlidingVelAxOLD()*SlidingVelAxOLD() < tol*tol) {error = 0;}// keep sliding when old velocity is very small (assumed to be zero) and new velocity is not
					else
					{
						if(SlidingVelAx()*SlidingVelAxOLD()	> 0) {error = 0;}// keep sliding when old and new velocity have similar directions
						else			// the velocity changed the direction 
						{
							IsSticking(1);
							error = fabs(SlidingVelAx()-SlidingVelAxOLD());		// this could be changed to error w.r.t. vmax
							//GetMBS()->UO() << "v changed sign\n";
						}
					}
				}
			}
		}
		return error;
	}
	


Vector3D FrictionConstraintUniDir::ComputeForce(double t) const
	{
		Vector3D u; // displacement
		Vector3D v; // velocity
		Vector3D f; // resulting force

		double stiffness = GetPenaltyStiffness();

		Vector3D vec_ax = GetAxialDir();
		
		if (UsePenaltyFormulation())
		{
			u = GetPos1() - GetPos2();
			//v = GetVel1() - GetVel2();

			u = u - StickingPosAx()*vec_ax ;	// origin of spring is last sticking point
			double Fn = GetNormalForce(u,stiffness);

			if(IsSticking())	// stick
			{
				f = GetAxialForce(u,stiffness);
			}

			else					// slip
			{
				f = GetSlidingFrictionForce(GetFrictionCoeff_kin());	// constant normal force for computation of friction force
			}
			//if(f.Norm()) GetMBS()->UO() << "f = " << f << "\n";
			return f;
		}
		else
		{
			GetMBS()->UO() << "ERROR: FrictionConstraintUniDir: ComputeForce just implemented for PenaltyFormulation";
			return Vector3D(0.0);
		}

	};

void FrictionConstraintUniDir::DrawElement() 
	{
		//Constraint::DrawElement();
		//double t_draw = mbs->GetDrawTime();

		if (GetDrawSizeScalar() == 0) return;

		Vector3D color_green(0.0,1.0,0.0);
		Vector3D color_pink(1.0,0.5,0.5);
		double scale = 1;

		if(XDataD(1)) 
		{
			mbs->SetColor(color_pink);
		}
		else
		{
			mbs->SetColor(color_green);
			scale = 0.5;
		}

		TArray<Vector3D> constr_dir;
		Vector3D p;

		constr_dir.Add(GetAxialDir()); 

		for (int i=1; i <= NKinPairs(); i++)
		{
			p = GetDrawPosition(i);

			if ( i== 2) 
			{ 
				constr_dir(1) *= -1;
			}

			for (int j=1; j<=constr_dir.Length(); j++)
			{
				mbs->DrawCone(p + (-GetDrawSizeScalar()*0.75*scale)*constr_dir(j),p,(GetDrawSizeScalar()*0.6*scale),6,1);
			}
		}
	};

int FrictionConstraintUniDir::GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables)
{
	// call base class routines
	BasePointJoint::GetAvailableSpecialValues(available_variables);

	// Automatic entries for this class 
	FrictionConstraintUniDir::GetAvailableSpecialValuesAuto(available_variables);

	// Manual entries for this class
	available_variables.Add(ReadWriteElementDataVariableType("Internal.stick_slip",0,0,0.,mystr("1 if sticking, 0 if sliding"))) ;
	available_variables.Add(ReadWriteElementDataVariableType("Internal.force_axial",0,0,0.,mystr("force in axial direction, applied to the kinematic pairs due to the constraint"))) ;
	available_variables.Add(ReadWriteElementDataVariableType("Internal.force",3,0,0.,mystr("force applied to the kinematic pairs due to the constraint. range: 1-3 corresponds to force in x-y-z direction"))) ;
	//JG 2013-01-11: only range of components is necessary!!!!

	return 0;
}

int FrictionConstraintUniDir::ReadSingleElementData(ReadWriteElementDataVariableType& RWdata) 		
{
	// call base class routine 	
	int rv = Constraint::ReadSingleElementData(RWdata);
	if (rv == 1) return 1;

	// manual things to read 
	double t = GetMBS()->GetTime();
	Vector3D f = ComputeForce(t);


	if( RWdata.variable_name == mystr("Internal.stick_slip") )
	{
		RWdata.value = IsSticking();
		return 1; 	
	}
	else if( RWdata.variable_name == mystr("Internal.force_axial") )
	{
		RWdata.value = f*GetAxialDir(); 
		return 1; 
	}
	else	if( RWdata.variable_name == mystr("Internal.force") )
	{
		if (RWdata.comp1 > 0 && RWdata.comp1 <= 3) //range check
		{
			RWdata.value = ComputeForce(t)(RWdata.comp1); ; 
			return 1; 
		}
		else return -2; 
	}
	return ReadSingleElementDataAuto(RWdata);
}

	// get force projected in special direction
	double FrictionConstraintUniDir::GetSpecialSensorValue(int nr, double time) const
	{
		double t = time; //GetMBS()->GetTime();
		Vector3D f = ComputeForce(t);
		double value = 0;	
		switch(nr)
		{
		case 1:		// force in axial direction
			value = f*GetAxialDir();
			break;
		case 2:		// force in tangential direction
			value = 0;
			break;
		case 3:		// force in normal direction
			value = Fn0;
			break;
		default:	// error handling
			return 0;	break;
		}
		return value;
	}


/*
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Controller
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int MBSController::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = Constraint::CheckConsistency(errorstr);
	if (rv) return rv;
	
	//check sensor numbers
	for (int i=1; i <= sensors.Length(); i++)
	{
		if (sensors(i) <= 0 || sensors(i) > GetMBS()->NSensors())
		{
			errorstr += mystr("    MBSController has invalid sensor number: ") + mystr(sensors(i)) + mystr("!\n");
			rv = 1;
		}
	}

	return rv;
}



void MBSController::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Constraint::GetElementData(edc);

	ElementData ed;

	Vector offsetsV(offsets.Length());
	for (int i=1; i <= offsets.Length(); i++) offsetsV(i) = offsets(i);
	Vector factorsV(factors.Length());
	for (int i=1; i <= factors.Length(); i++) factorsV(i) = factors(i);

	SetElemDataIVector(edc, sensors, "Sensor_numbers"); edc.Last().SetToolTipText("[Sensor number1, Sensor number2, ...]");
	edc.Last().SetVariableLength();
	SetElemDataVector(edc, offsetsV, "Sensor_offsets"); edc.Last().SetToolTipText("[Sensor offset1, Sensor offset2, ...]");
	edc.Last().SetVariableLength();
	SetElemDataVector(edc, factorsV, "Sensor_factors"); edc.Last().SetToolTipText("[Sensor factor1, Sensor factor2, ...]");
	edc.Last().SetVariableLength();
	SetElemDataIVector(edc, options, "Sensor_options"); edc.Last().SetToolTipText("Option=1 means integrate; [Sensor option1, Sensor option2, ...]");
	edc.Last().SetVariableLength();

	SetElemDataMathFunc(edc, mathfunc, "Ctrl_limit_func"); edc.Last().SetToolTipText("Function that modifies the control output (limitation, activation, etc.)");
}

int MBSController::SetElementData(const ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Constraint::SetElementData(edc);

	//store data in temporary array first:
	TArray<int> sensorsI;
	Vector factorsI;
	Vector offsetsI;
	TArray<int> optionsI;

	GetElemDataIVector(mbs, edc, "Sensor_numbers", sensorsI, 1);
	GetElemDataVector(mbs, edc, "Sensor_offsets", offsetsI, 1);
	GetElemDataVector(mbs, edc, "Sensor_factors", factorsI, 1);
	GetElemDataIVector(mbs, edc, "Sensor_options", optionsI, 1);

	int len = sensorsI.Length();
	if (offsetsI.Length() != len) {offsetsI.SetLen(len); GetMBS()->EDCError("Wrong number of offsets in Controller"); rv = 0;}
	if (factorsI.Length() != len) {factorsI.SetLen(len); GetMBS()->EDCError("Wrong number of factors in Controller"); rv = 0;}
	if (optionsI.Length() != len) {optionsI.SetLen(len); GetMBS()->EDCError("Wrong number of options in Controller"); rv = 0;}

	FlushSensors();
	for (int i=1; i <= sensorsI.Length(); i++)
	{
		AddSensor(sensorsI(i), offsetsI(i), factorsI(i), optionsI(i));
		//Error checking not possible because of file loading
		//if (sensorsI(i) >= 1 && sensorsI(i) <= GetMBS()->NSensors())
		//{
		//	AddSensor(sensorsI(i), offsetsI(i), factorsI(i), optionsI(i));
		//}
		//else
		//{
		//	GetMBS()->EDCError("Controller: invalid sensor number");
		//	rv = 0;
		//}
	}

	rv = rv*GetElemDataMathFunc(mbs, edc, "Ctrl_limit_func", mathfunc, 1);

	return rv;
}

void MBSController::EvalF(Vector& f, double t)
{
	int j = 1;
	for (int i=1; i <= sensors.Length(); i++)
	{
		if (options(i)&1) //integrate
		{
			f(j) = GetSensor(i).GetValues()(1);
			j++;
		}
	}
}

int MBSController::AddSensor(int sensor_num, double offset, double factor, int option) 
{
	//UO() << "Sensor: " << sensors.Length()+1 << ", num=" << sensor_num << ", off=" << offset << ", fac=" << factor << ", option=" << option << "\n";
	offsets.Add(offset);
	factors.Add(factor);
	options.Add(option);

	if (option&1)  //if sensor integrated --> add variable for explicit size
	{
		nvariables++;
		x_init = Vector(SS());
	}

	return sensors.Add(sensor_num);
}

double MBSController::ComputeControlValue()
{
	double t = GetMBS()->GetTime();

	double rv = 0;
	double val;
	int j=1;
	for (int i=1; i <= sensors.Length(); i++)
	{
		if (options(i)&1) //integrate
		{
			val = XG(j); 
			j++;
		}
		else
		{
			val = GetSensor(i).GetValues()(1);
		}

		rv += (val + offsets(i))*factors(i);
	}
	rv *= mathfunc.Evaluate(t);
	//UO() << "val=" << rv << "\n";

	return rv;
}
*/


// ##############################################################################################
//$ DR 2013: Rope3D added

int Rope3D::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = Constraint::CheckConsistency(errorstr);
	if (rv) return rv;
	if(!UsePenaltyFormulation())
	{
			errorstr = mystr("ERROR: Rope3D: only penalty formulation is possible\n");
			rv = 1;
	}
	
	return rv;
}

// for constraints only, these are the necessary (!) access functions!	//$ DR 2013-02-13
void Rope3D::GetNecessaryKinematicAccessFunctions(TArray<int> &KinAccFunc, int numberOfKinematicPair)
{
	KinAccFunc.SetLen(0);
	int kaf = (int)(TKinematicsAccessFunctions(TKAF_none));
	kaf = (int)(TKinematicsAccessFunctions(TKAF_position+TKAF_D_pos_D_q));	
	KinAccFunc.Add(kaf);
}

void Rope3D::ElementDefaultConstructorInitialization()
{
	SetPenaltyFormulation(1);
	SetUseLocalCoordinateSystem(0);
	loccoords.SetLen(2);
	loccoords(1)=Vector3D(0,0,0);
	loccoords(1)=Vector3D(1,0,0);
	initial_length = 1.;
	coiled_length = 0;
	velocity_constraint = 0;
	SetDampingCoeff(0.0);
	elements.Set2(0,0);
	elementname = GetElementSpec();
	nodes.SetLen(0); // not used

	Vector datainit(DataS());
	datainit.SetAll(0.);
	SetDataInit(datainit);
}

void Rope3D::Initialize() 
{
	//do not use BasePointJoint::Initialize because then the wrong SetPos2ToGlobalCoord is called
	Constraint::Initialize();

	initial_length = GetLengthOfRope();
	GetMBS()->UO(UO_LVL_ext) << "Info: Rope3D: initial length of rope = " << initial_length << "\n";
};

void Rope3D::CopyFrom(const Element& e)
{
	BasePointJoint::CopyFrom(e);
	const Rope3D& ce = (const Rope3D&)e;
	initial_length = ce.initial_length;
	coiled_length = ce.coiled_length;
}

void Rope3D::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	BasePointJoint::GetElementData(edc);
	ElementData ed;
	Matrix m;
	m.SetSize(loccoords.Length(),3);
	for(int i=1; i<=loccoords.Length();i++)
	{
		Vector v;
		v.Set3D(loccoords(i).X(),loccoords(i).Y(),loccoords(i).Z());
		m.SetRowVec(v,i);
	}
	int cols = 3;
	int rows = loccoords.Length();
	ed.SetMatrix(m.GetMatPtr(),rows,cols,"positions");
	ed.SetToolTipText("(local) positions of the suspension points");
	edc.TreeAdd("Geometry",ed);

	//ed.SetDouble(GetPenaltyStiffness()/initial_length,"rope_stiffness",0,0,1,0); ed.SetToolTipText("[N] stiffness parameter c of the rope, F = c * (l-l0)/l0");
	ed.SetDouble(GetPenaltyStiffness(),"rope_stiffness",0,0,1,0); ed.SetToolTipText("[N] stiffness parameter c of the rope, F = c * (l-l0)/l0");
	edc.TreeAdd("Physics.Penalty",ed);

	ed.SetDouble(GetPenaltyStiffness()*(initial_length+coiled_length),"spring_stiffness",0,0,1,0); ed.SetToolTipText("total stiffness c1 of the rope F = c1 * (l-l0)"); ed.SetLocked(1);
	edc.TreeAdd("Physics.Penalty",ed);

}

int Rope3D::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = BasePointJoint::SetElementData(edc);

	int rows,cols;
	double *mp;

	if(edc.TreeGetMatrix("Geometry.positions",&mp,rows,cols))
	{
		Matrix m(rows,cols,mp);
		loccoords.SetLen(rows);
		for(int i=1; i<=rows; i++)
		{
			loccoords(i) = Vector3D(m(i,1),m(i,2),m(i,3));
		}
	}

	double c = edc.TreeGetDouble("Physics.Penalty.rope_stiffness");		
	SetPenaltyStiffness(c);														

	return rv;
}

int Rope3D::GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables)
{
	//BasePointJoint::GetAvailableSpecialValues(available_variables);

	// Automatic entries for this class 
	Rope3D::GetAvailableSpecialValuesAuto(available_variables);

	// Manual READ entries for this class
	available_variables.Add(ReadWriteElementDataVariableType("Internal.actor_force",1,0,0.,mystr("force in the rope"),TRWElementDataRead)) ;
	available_variables.Add(ReadWriteElementDataVariableType("Internal.rope_length",1,0,0.,mystr("length of the rope"),TRWElementDataRead)) ;

	// Manual WRITE entries for this class
	available_variables.Add(ReadWriteElementDataVariableType("Internal.coiled_length",1,0,0.,mystr("(additional) length of the rope that is provided by a coil. length = rope_length + coiled_length"),TRWElementDataReadWrite)) ;

	return 0;
}

int Rope3D::ReadSingleElementData(ReadWriteElementDataVariableType& RWdata) 		
{
	// call base class routine 	
	//int rv = BasePointJoint::ReadSingleElementData(RWdata);
	//if (rv == 1) return 1;

	// manual things to read  
	if(RWdata.variable_name.CStrCompare("Internal.actor_force") )
	{
		double t = mbs->GetTime();
		RWdata.value = ComputeRopeForce(mbs->GetTime()); 
		return 1; 
	}
	if(RWdata.variable_name.CStrCompare("Internal.rope_length") )
	{
		RWdata.value = GetLengthOfRope();
		return 1; 
	}
	if(RWdata.variable_name.CStrCompare("Internal.coiled_length") )
	{
		RWdata.value = coiled_length;
		return 1; 
	}

	return ReadSingleElementDataAuto(RWdata);
}

int Rope3D::WriteSingleElementData(const ReadWriteElementDataVariableType& RWdata) 		
{
	// call base class routine ( not required in Element )	
	int rv = BasePointJoint::WriteSingleElementData(RWdata);
	if (rv == 1) return 1;
	//manual things to write
	if(RWdata.variable_name.CStrCompare("Internal.coiled_length"))
	{
		coiled_length = RWdata.value;
		return 1; 
	}

	return WriteSingleElementDataAuto(RWdata);
}

void Rope3D::SetPositions(TArray<int> element_numbers,  TArray<Vector3D> coordinates)
{
	// for global positions (=ground) the element number '0' is inserted
	// if i > elements.length resize is done automatically by TArray

	if(element_numbers.Length()==0) {assert(0);}
	if(element_numbers.Length()!= coordinates.Length()) {assert(0);}
	elements = element_numbers;
	loccoords = coordinates;
}

Vector3D Rope3D::GetPosition(int i) const
{
	Vector3D pos;
	if(elements.Length()!=loccoords.Length())
	{
		GetMBS()->UO(UO_LVL_err) << "ERROR: Rope3D: The number of elements must be equal to the number of positions!";
		return Vector3D(0.0);
	}
	else
	{
		if(elements(i)==0)				// ground position
		{
			pos = loccoords(i);
		}
		else											
		{
			pos = GetBody3D(i).GetPos(loccoords(i));
		}
	}

	return pos;
}

Vector3D Rope3D::GetDrawPosition(int i) const
{
	Vector3D pos;
	if(elements.Length()!=loccoords.Length())
	{
		GetMBS()->UO(UO_LVL_err) << "ERROR: Rope3D: The number of elements must be equal to the number of positions!";
		return Vector3D(0.0);
	}
	else
	{
		if(elements(i)==0)				// ground position
		{
			pos = loccoords(i);
		}
		else											
		{
			pos = GetBody3D(i).GetPosD(loccoords(i));
		}
	}

	return pos;
}

void Rope3D::PostprocessingStep()
{
	double t = mbs->GetTime();

	double dl = GetLengthOfRope()-OldLength();
	double dt = t-OldTime();
	
	if(dt)
		SetVelocityOfRope(dl/dt);
	else
		SetVelocityOfRope(0);

	OldLength(GetLengthOfRope());
	OldTime(t);
}

double Rope3D::GetLengthOfRope() const
{
	double l=0;
	Vector3D v;
	for(int i=1; i<=elements.Length()-1; i++)
	{
		v = GetPosition(i+1)-GetPosition(i);
		l += v.Norm();
	}
	return l/*+coiled_length*/;
}

//double Rope3D::GetVelocityOfRope() const	//$ DR 2013-08-29 moved to postprocessing step
//{
//	double dl = GetLengthOfRope()-OldLength();
//	double dt = (mbs->GetTime())-OldTime();
//	if(dt)
//		return dl/dt;
//	else
//		return 0;
//}

double Rope3D::ComputeRopeForce(double t) const
{
	double f; // resulting force
	if (UsePenaltyFormulation())
	{
		double l = GetLengthOfRope();
		double l0 = initial_length+coiled_length;
		double u = l-l0;
		if( (l0==0) || (l0<0))
		{
			GetMBS()->UO(UO_LVL_warn) << "WARNING: initial length of rope is equal to or smaller than 0, this is not possible and therefore the force is set to 0.\n";
			return 0;
		}
		f = GetPenaltyStiffness()*u/l0;	
		if(GetDampingCoeff())
		{
			double v = GetVelocityOfRope();
			f += GetDampingCoeff()*v;
		}
		if(f < 0)
		{
			GetMBS()->UO(UO_LVL_warn) << "WARNING: force in rope is negative, this is not possible and therefore set to 0.\n";
			return 0;
		}
		return f;
	}
	else
	{
		GetMBS()->UO(UO_LVL_err) << "ERROR: Rope3D: ComputeForce just implemented for PenaltyFormulation";
		return 0;
	}
};

Vector3D Rope3D::ComputeForceDirection(int i) const
{
	Vector3D v1 = Vector3D(0);
	Vector3D v2 = v1;

	if(i != NKinPairs())
	{
		v1 = GetPosition(i+1)-GetPosition(i);
		v1.Normalize();
	}
	if(i != 1)
	{
		v2 = GetPosition(i-1)-GetPosition(i);
		v2.Normalize();
	}
	Vector3D dir = v1 + v2;
	dir.Normalize();
	return dir;
};


void Rope3D::EvalF2(Vector& f, double t)
{
	//f = [f1 f2 ... fn], where fi is the residual vector of element i
	if(!UsePenaltyFormulation())
		return;	// no penalty method --> Lagrange multiplier --> no EvalF2

	double rope_force = ComputeRopeForce(t);
	Vector3D force;

	int offset = 0;
	for(int i=1; i<=NKinPairs(); i++)
	{
		int el = elements(i);
		if(el)	// no force applied to ground
		{
			// compute force in rope
			Vector3D dir = ComputeForceDirection(i);
			force = rope_force * dir;
			if((i != 1) && (i != NKinPairs()))		// not the first and not the last position
			{
				force = 2*force;
			}

			// enter force vector in f-vector
			GetBody3D(i).GetdPosdqT(loccoords(i),dpdq);
			for (int j=1; j<=GetElem(i).SOS(); j++)
			{
				f(j+offset) += dpdq(j,1)*force.X() + dpdq(j,2)*force.Y() + dpdq(j,3)*force.Z();
			}

			offset += GetElem(i).SOS();	// update offset in f-vector for next element
		}
	}
};

void Rope3D::LinkToElementsPenalty()
{
	LTGreset();

	// add all SOS dofs from the elements
	//Position(first SOS) 
	for (int k=1; k <= NKinPairs(); k++)
	{
		if(elements(k)!=0)						
		{
			for (int i=1; i <= GetElem(k).SOS(); i++)
			{
				AddLTG(GetElem(k).LTG(i));
			}
		}
	}
	//and Velocity (second SOS):
	for (int k=1; k <= NKinPairs(); k++)
	{
		if(elements(k)!=0)						
		{
			for (int i=1; i <= GetElem(k).SOS(); i++)
			{
				AddLTG(GetElem(k).LTG(i+GetElem(k).SOS()));
			}
		}
	}
};


int Rope3D::SOS() const 
{
	int nsos = 0;
	for(int i=1; i<=elements.Length(); i++)
	{
		if(elements(i))
		{
			nsos += GetElem(i).SOS();
		}
	}
	return nsos;
};  // explicit size, number of constrained dofs


void Rope3D::DrawElement() 
{
	mbs->SetColor(GetCol());

	Vector3D p1,p2;
	for (int i=1; i <= NKinPairs(); i++)
	{
		p1 = GetDrawPosition(i);
		mbs->DrawSphere(p1, GetDrawSizeScalar());
		if(i < NKinPairs())
		{
			p2 = GetDrawPosition(i+1);
			mbs->MyDrawLine(p1,p2,1.);
		}
	}
};


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void FrictionConstraint::ElementDefaultConstructorInitialization()
{
	SetPenaltyFormulation(1); 
	SetUseLocalCoordinateSystem(0);

	loccoords.SetLen(2);
	loccoords.SetAll(1);
	elements.SetLen(2);
	elements(1)=0;
	elements(2)=0;
	elementname = GetElementSpec();

	SetVelocityTolerance(1e-5);
	SetFrictionCoeff_st(0.15);
	SetFrictionCoeff_kin(0.1);

	Vector datainit(DataS());
	datainit.SetAll(0.);
	datainit(1)=1.;	// start with sticking
	SetDataInit(datainit);
	Fn0 = 0.;
	keep_sliding = 0;

	damping_coeff = 0;

	//x_init = Vector(SS());
}

int FrictionConstraint::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = Constraint::CheckConsistency(errorstr);
	if (rv) return rv;

	if(GetVelocityTolerance() == 0)
	{
			mbs->UO() << "ERROR: FrictionConstraint:: velocity_tolerance MUST be greater than zero.\n";
			rv=1;
	}

	return rv;
}


void FrictionConstraint::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Constraint::GetElementData(edc);
	ElementData ed;
	ed.SetInt(GetElnum(1), "element_number"); ed.SetToolTipText("element number for coordinate 1"); edc.TreeAdd("Coordinate1",ed);
	ed.SetInt(loccoords(1), "local_coordinate"); ed.SetToolTipText("Local coordinate of element 1 to be constrained"); edc.TreeAdd("Coordinate1",ed);

	ed.SetInt(GetElnum(2), "element_number"); ed.SetToolTipText("element number for coordinate 2; for ground joint, set element number to zero"); edc.TreeAdd("Coordinate2",ed);
	ed.SetInt(loccoords(2), "local_coordinate"); ed.SetToolTipText("Local coordinate of element 2 to be constrained"); edc.TreeAdd("Coordinate2",ed);

}

int FrictionConstraint::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Constraint::SetElementData(edc);

	int tmp;
	GetElemDataInt(mbs, edc, "Coordinate1.local_coordinate", loccoords(1), 1);
	GetElemDataInt(mbs, edc, "Coordinate1.element_number", tmp, 1);
	elements(1) = tmp;

	GetElemDataInt(mbs, edc, "Coordinate2.local_coordinate", loccoords(2), 1);
	GetElemDataInt(mbs, edc, "Coordinate2.element_number", tmp, 1);
	elements(2) = tmp;

	return rv;
}


int FrictionConstraint::WriteSingleElementData(const ReadWriteElementDataVariableType& RWdata) 
{
	// call base class routine ( not required in Element )	
	int rv = Constraint::WriteSingleElementData(RWdata);
	if (rv == 1) return 1;
	//manual things to write

	if(RWdata.variable_name.CStrCompare("Internal.normal_force"))
	{
		Fn0 = RWdata.value;
	}

	return WriteSingleElementDataAuto(RWdata);
}

int FrictionConstraint::ReadSingleElementData(ReadWriteElementDataVariableType& RWdata) 		
{
	// call base class routine 	
	int rv = Constraint::ReadSingleElementData(RWdata);
	if (rv == 1) return 1;

	// manual things to read 
	double t = GetMBS()->GetTime();
	double f = ComputeForce(t);

	if( RWdata.variable_name == mystr("Internal.stick_slip") )
	{
		RWdata.value = IsSticking();
		return 1; 	
	}
	else if( RWdata.variable_name == mystr("Internal.force") )
	{
		RWdata.value = f; 
		return 1; 
	}
	else if( RWdata.variable_name == mystr("Internal.force_abs") )
	{
		RWdata.value = abs(f); 
		return 1; 
	}
	else if( RWdata.variable_name == mystr("Internal.normal_force") )
	{
		RWdata.value = Fn0; 
		return 1; 
	}
	return ReadSingleElementDataAuto(RWdata);
}

int FrictionConstraint::GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables)
{
	// no parent class, flush the array here ( and only here )
	available_variables.Flush();

	// Automatic entries for this class 
	Element::GetAvailableSpecialValuesAuto(available_variables);

	// Manual READ entries for this class
	available_variables.Add(ReadWriteElementDataVariableType("Internal.stick_slip",0,0,0.,mystr("1 if sticking, 0 if sliding"))) ;
	available_variables.Add(ReadWriteElementDataVariableType("Internal.force",0,0,0.,mystr("force, applied to the kinematic pairs due to the constraint"))) ;
	available_variables.Add(ReadWriteElementDataVariableType("Internal.force_abs",0,0,0.,mystr("absolute value of the force, applied to the kinematic pairs due to the constraint"))) ;

	// Manual READ and WRITE entries for this class
	available_variables.Add(ReadWriteElementDataVariableType("Internal.normal_force",0,0,0.,mystr("force, applied to the kinematic pairs due to the constraint"),TRWElementDataReadWrite)) ;

	return 0;
}

// Get force that acts during stick-phase. Returns force in axial direction.
double FrictionConstraint::GetAxialForce(double u, double stiffness) const
{
	return u*stiffness + SlidingVel()*GetDampingCoeff();
}

// NormalForce Fn is positve for pressure!
double FrictionConstraint::GetNormalForce(double u, double stiffness)const
{
	return Fn0;
}

// Friction Force when body is sliding: F = µk*Fn
double FrictionConstraint::GetSlidingFrictionForce(double coef_kin)const
{
	int sign = 1;
	if(SlidingVel()<0) { sign = -1;}
	double F = sign*coef_kin*(Fn0);
	
	if(keep_sliding)
	{
		if(abs(SlidingVel()) < GetVelocityTolerance())
		{
			F = F * abs(SlidingVel())/GetVelocityTolerance();
		}
	}

	return F;
}

// maximal Force at which sliding starts: Fmax = µs*Fn
double FrictionConstraint::GetTraction(double t, double u) const
{
	return 	GetFrictionCoeff_st()*(Fn0);	
}

double FrictionConstraint::ComputeForce(double t) const
{
	double u; // displacement
	double v; // velocity
	double f; // resulting force

	double stiffness = GetPenaltyStiffness();

	if (UsePenaltyFormulation())
	{
		u = GetPos1() - GetPos2();
		//v = GetVel1() - GetVel2();

		u = u - StickingPos();	// origin of spring is last sticking point
		double Fn = GetNormalForce(u,stiffness);

		if(IsSticking())	// stick
		{
			f = GetAxialForce(u,stiffness);
		}

		else					// slip
		{
			f = GetSlidingFrictionForce(GetFrictionCoeff_kin());	// constant normal force for computation of friction force
		}
		return f;
	}
	else
	{
		GetMBS()->UO() << "ERROR: FrictionConstraint: ComputeForce just implemented for PenaltyFormulation";
		return 0;
	}

};

void FrictionConstraint::EvalF2(Vector& f, double t)
{
	//f = [f1 f2], where f1 is the residual vector of constraint element1 and f2 of constraint element 2

	if(!UsePenaltyFormulation()) return; 

	double force = ComputeForce(t);
	for (int i=1; i <= NE(); i++)
	{
		double sign = 1.;
		int offset = 0;
		if (i==2) 
		{
			sign = -1.; //sign = du/dq;
			offset = GetElem(1).SOS();
		}
		f(loccoords(i)+offset) -= sign*force;
	}
}


void FrictionConstraint::PostprocessingStep()
{
	double lastXData = GetMBS()->GetLastDataVector()(LTGdata(1));

	if (IsSticking() == 1 && lastXData != 1) //switch from sliding to stick ==> remember actual sticking position
	{
		double u = GetPos1() - GetPos2();

		// adjust u, such that jumps of force at slide-stick-transition are reduced
		double du = abs(GetSlidingFrictionForce(GetFrictionCoeff_kin()))/GetPenaltyStiffness();
		if(SlidingVelOLD() < 0) u+=du;
		else u-=du;

		StickingPos(u);														// origin of spring is last sticking point
		SlidingVelOLD(SlidingVel());							// last sliding velocity
		if(mbs->UO(UO_LVL_dbg1).PrintMsg())
		{	
			GetMBS()->UO() << "change slide to STICK: elementNum = " << GetOwnNum() <<", t = " << GetMBS()->GetTime() <<", u = " << StickingPos() << "\n"; 
		}
	}
	else if (IsSticking() == 0) //during sliding ==> remember last sliding velocity
	{
		SlidingVelOLD(SlidingVel());		// last sliding velocity
		if(lastXData && (mbs->UO(UO_LVL_dbg1).PrintMsg()) )
		{	
			GetMBS()->UO() << "change stick to SLIDE: elementNum = " << GetOwnNum() <<", t = " << GetMBS()->GetTime() <<", v = " << SlidingVel() << "\n";
		}
	}
}

double FrictionConstraint::PostNewtonStep(double t)
{
	double error = 0.;

	double u; // displacement
	double v; // velocity

	double stiffness = GetPenaltyStiffness();

	u = GetPos1() - GetPos2();
	v = GetVel1() - GetVel2(); 

	SlidingVel(v);

	u = u - StickingPos();	// origin of spring is last sticking point

	if(IsSticking())	// stick
	{
		double Fa = GetAxialForce(u,stiffness);
		double threshhold = GetTraction(t,u);

		if(fabs(Fa) <= threshhold) // keep sticking
		{
			error = 0.;
		}
		else									// start sliding
		{
			error = fabs(Fa-threshhold);
			if (threshhold!=0) {error /= threshhold;}
			IsSticking(0);
		}
	}
	else					// sliding
	{
		if(keep_sliding)				// no sticking again
		{
			error = 0;
		}
		else
		{
			double tol = GetVelocityTolerance(); // if v < tol than it is assumed that v = 0	

			if( SlidingVel()*SlidingVel() < tol*tol) // start sticking when new velocity is very small (assumed to be zero)
			{
				IsSticking(1);
				error = fabs((SlidingVel()*SlidingVel())/tol);
			}
			else
			{
				if(SlidingVelOLD()*SlidingVelOLD() < tol*tol) {error = 0;}// keep sliding when old velocity is very small (assumed to be zero) and new velocity is not
				else
				{
					if(SlidingVel()*SlidingVelOLD()	> 0) {error = 0;}// keep sliding when old and new velocity have similar directions
					else			// the velocity changed the direction 
					{
						IsSticking(1);
						error = fabs(SlidingVel()-SlidingVelOLD());		// this could be changed to error w.r.t. vmax
						//GetMBS()->UO() << "v changed sign\n";
					}
				}
			}
		}
	}
	return error;
}


void FrictionConstraint::DrawElement() 
{
	Constraint::DrawElement();

	if (GetDrawSizeScalar() == 0) return;

	Vector3D color_green(0.0,1.0,0.0);
	Vector3D color_pink(1.0,0.5,0.5);
	double scale = 1;

	if(XDataD(1)) 
	{
		mbs->SetColor(color_pink);
	}
	else
	{
		mbs->SetColor(color_green);
		scale = 0.5;
	}

	Vector3D p;
	Vector3D v0(0);
	for (int i=1; i <= NE(); i++)
	{
		p = GetElem(i).GetPosD(v0);
		mbs->DrawSphere(p, scale*GetDrawSizeScalar());
	}

};


//----------------------------------------------------------------------------------------------------------------------------------------
//$ DR 2013-09-30 added Contact1D
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// set functions - for element-based formulations
void Contact1D::SetContact1DElement(int en1, int lc1) 
{
	ElementDefaultConstructorInitialization();
	elements(1) = en1;
	lcoord1 = lc1;
}

void Contact1D::SetContact1DElement(int en1, int en2, int lc1, int lc2)
{	
	ElementDefaultConstructorInitialization();

	elements(1) = en1;
	lcoord1 = lc1;
	elements(2) = en2;
	lcoord2 = lc2;
}

// set functions - for node-based formulations
void Contact1D::SetContact1DNode(int nn1, int lc1) 
{
	ElementDefaultConstructorInitialization();
	elements(1) = 0;	//$ DR 2013-12-09 default value 1. element=0 for nodal formulation
	nodes(1) = nn1;
	lcoord1 = lc1;
}

void Contact1D::SetContact1DNode(int nn1, int nn2, int lc1, int lc2)
{	
	ElementDefaultConstructorInitialization();
	elements(1) = 0;	//$ DR 2013-12-09 default value 1. element=0 for nodal formulation
	elements(2) = 0;	//$ DR 2013-12-09 default value 1. element=0 for nodal formulation

	nodes(1) = nn1;
	lcoord1 = lc1;
	nodes(2) = nn2;
	lcoord2 = lc2;
}

void Contact1D::ElementDefaultConstructorInitialization()
{
	SetPenaltyFormulation(1); 
	SetUseLocalCoordinateSystem(0);

	lcoord1 = 1;
	lcoord2 = 1;
	lpos1 =		0;
	lpos2 =		0;

	contact_direction = 1;

	elements.SetLen(2);
	elements(1)=1;		//$ DR 2013-12-09 default value 1 necessary for autogenerated docu
	elements(2)=0;
	nodes.SetLen(2);
	nodes(1)=0;
	nodes(2)=0;
	elementname = GetElementSpec();
	
	damping_coeff = 0;

	mode = 1;
}

int Contact1D::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = Constraint::CheckConsistency(errorstr);
	if (rv) return rv;

	return rv;
}


void Contact1D::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Constraint::GetElementData(edc);
	ElementData ed;
	ed.SetInt(GetElnum(1), "element_number"); ed.SetToolTipText("element number for coordinate 1; set to zero if you use nodal coordinates!"); edc.TreeAdd("Coordinate1",ed);
	ed.SetInt(GetElnum(2), "element_number"); ed.SetToolTipText("element number for coordinate 2; for ground joint or nodal coordinates, set element number to zero"); edc.TreeAdd("Coordinate2",ed);
	ed.SetInt(nodes(1), "node_number"); ed.SetToolTipText("(just used if element number > 0) node number for coordinate 1"); edc.TreeAdd("Coordinate1",ed);
	ed.SetInt(nodes(2), "node_number"); ed.SetToolTipText("(just used if element number > 0) node number for coordinate 2; for ground joint, set node number to zero"); edc.TreeAdd("Coordinate2",ed);
}

int Contact1D::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Constraint::SetElementData(edc);
	elements(1)=edc.TreeGetInt("Coordinate1.element_number");
	elements(2)=edc.TreeGetInt("Coordinate2.element_number");
	nodes(1)=edc.TreeGetInt("Coordinate1.node_number");
	nodes(2)=edc.TreeGetInt("Coordinate2.node_number");

	if(mode==1)
	{
		SetPenaltyFormulation(1); 
	}
	if(mode==2)
	{
		SetPenaltyFormulation(0); 
	}
	return rv;
}

void Contact1D::CopyFrom(const Element& e)
{
	Constraint::CopyFrom(e);
	const Contact1D& ce = (const Contact1D&)e;
	lcoord1 = ce.lcoord1;
	lcoord2 = ce.lcoord2;
	lpos1 = ce.lpos1;
	lpos2 = ce.lpos2;
	damping_coeff=ce.damping_coeff;			
	mode=ce.mode;
	nodes = ce.nodes;
	contact_direction = ce.contact_direction;
}

double Contact1D::GetPos1() const
{
	if(elements(1))
		return GetElem(1).XG(lcoord1) + lpos1;
	return GetNode(1).XG(lcoord1) + lpos1;
}

double Contact1D::GetPos2() const 
{
	if(elements(2))
		return GetElem(2).XG(lcoord2) + lpos2;
	if(nodes(2))
		return GetNode(2).XG(lcoord2) + lpos2;
	return lpos2;		// ground
}

double Contact1D::GetVel1() const
{
	if(elements(1))
		return GetElem(1).XGP(lcoord1);
	return GetNode(1).XGP(lcoord1);
}

double Contact1D::GetVel2() const
{
	if(elements(2))
		return GetElem(2).XGP(lcoord2);
	if(nodes(2))
		return GetNode(2).XGP(lcoord2);
	return 0;		// ground
}

int Contact1D::NE() const 
{
	if(elements(2))
		return 2;
	else if(elements(1))
		return 1;
	else
		return 0; // when adding the empty constraint or if the constraint is node-based
}

int Contact1D::WriteSingleElementData(const ReadWriteElementDataVariableType& RWdata) 
{
	// call base class routine ( not required in Element )	
	int rv = Constraint::WriteSingleElementData(RWdata);
	if (rv == 1) return 1;
	//manual things to write

	//if(RWdata.variable_name.CStrCompare("Internal.normal_force"))
	//{
	//	Fn0 = RWdata.value;
	//}

	return WriteSingleElementDataAuto(RWdata);
}

int Contact1D::ReadSingleElementData(ReadWriteElementDataVariableType& RWdata) 		
{
	// call base class routine 	
	int rv = Constraint::ReadSingleElementData(RWdata);
	if (rv == 1) return 1;

	//// manual things to read 
	//double t = GetMBS()->GetTime();
	//double f = ComputeForce(t);

	//if( RWdata.variable_name == mystr("Internal.stick_slip") )
	//{
	//	RWdata.value = IsSticking();
	//	return 1; 	
	//}
	//else if( RWdata.variable_name == mystr("Internal.force") )
	//{
	//	RWdata.value = f; 
	//	return 1; 
	//}
	//else if( RWdata.variable_name == mystr("Internal.force_abs") )
	//{
	//	RWdata.value = abs(f); 
	//	return 1; 
	//}
	//else if( RWdata.variable_name == mystr("Internal.normal_force") )
	//{
	//	RWdata.value = Fn0; 
	//	return 1; 
	//}
	return ReadSingleElementDataAuto(RWdata);
}

int Contact1D::GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables)
{
	// no parent class, flush the array here ( and only here )
	available_variables.Flush();

	// Automatic entries for this class 
	Element::GetAvailableSpecialValuesAuto(available_variables);

	// Manual READ entries for this class
	//available_variables.Add(ReadWriteElementDataVariableType("Internal.stick_slip",0,0,0.,mystr("1 if sticking, 0 if sliding"))) ;
	//available_variables.Add(ReadWriteElementDataVariableType("Internal.force",0,0,0.,mystr("force, applied to the kinematic pairs due to the constraint"))) ;
	//available_variables.Add(ReadWriteElementDataVariableType("Internal.force_abs",0,0,0.,mystr("absolute value of the force, applied to the kinematic pairs due to the constraint"))) ;

	// Manual READ and WRITE entries for this class
	//available_variables.Add(ReadWriteElementDataVariableType("Internal.normal_force",0,0,0.,mystr("force, applied to the kinematic pairs due to the constraint"),TRWElementDataReadWrite)) ;

	return 0;
}

// the base-class implementation needs to be redefined in case of node-based contact
void Contact1D::LinkToElements()
{
	// the implementation below is valid only for the penalty formulation
	// add all relevant dofs from the elements or from the nodes
	// Position (first SOS)
	for (int k=1; k <= 2; k++)
	{
		int lc = k == 1 ? lcoord1 : lcoord2;
		if(elements(k))
			AddLTG(GetElem(k).LTG(lc));
		else
		{
			if(nodes(k))
				AddLTG(GetNode(k).LTG(lc));
		}
	}
	// and Velocity (second SOS)
	for (int k=1; k <= 2; k++)
	{
		int lc = k == 1 ? lcoord1 : lcoord2;
		if(elements(k))
			AddLTG(GetElem(k).LTG(lc + GetElem(k).SOS()));
		else
		{
			if(nodes(k))
				AddLTG(GetNode(k).LTG(lc + GetNode(k).SOS()));
		}
	}
}

int Contact1D::SOS() const
{
	if (!UsePenaltyFormulation()) 
	{
		return 0;
	}
	else
	{
		int nsos = 0;
		for (int k=1; k <= 2; k++)
		{
			if(elements(k) || nodes(k))
				nsos++;
		}
		return nsos;
	}
}


double Contact1D::ComputeForce(double t) const
{
	double u; // displacement
	double v; // velocity
	double f=0; // resulting force

	double stiffness = GetPenaltyStiffness();

	u = contact_direction * (GetPos1() - GetPos2());
	v = contact_direction * (GetVel1() - GetVel2());

	if(u > 0.)
	{
		// do nothing, no contact
		return 0.;
	}
	else
	{
		switch(mode)
		{
		case 1: // linear elastic
			{
				f = contact_direction * (stiffness*u+damping_coeff*v);
				break;
			}
		case 2: // hard contact
			{
				f = XG(1);
				break;
			}
			//case 3: 
			//	break;
		}
	}
	return f;
}

void Contact1D::EvalF2(Vector& f, double t)
{
	//f = [f1 f2], where fi are the corresponding residual forces of the contacting degrees of freedom of elements or nodes

	if(!UsePenaltyFormulation()) return; 

	double force = ComputeForce(t);

	f(1) -= force;
	if(elements(2) || nodes(2)) 
	{
		f(2) += force;
	}
}

void Contact1D::EvalG(Vector& f, double t)
{
	if(UsePenaltyFormulation()) return; 

	double u = GetPos1() - GetPos2();
	if(u > 0.)
	{
		// do nothing, no contact
		f(1) -= 0;
	}
	else
	{
		f(1) += XG(1);
	}
}

void Contact1D::AddElementCqTLambda(double t, int locelemind, Vector& f)
{	//-C_q^T \lambda = [dC/dq;]^T [\lambda_i] 
	
	//int dim = mbs->GetElement(elements(1)).Dim();
	//Matrix dpdq;

	//switch(dim)
	//	{
	//	case 1: 
	//		{
	//			dpdq.SetSize(1,1);
	//			dpdq(1,1)=1;
	//			break;
	//		}
	//	case 2: 
	//		{
	//			Body2D* b2 =(Body2D*)mbs->GetElementPtr(elements(1));
	//			Vector2D lc2(0.);
	//			lc2(lcoord1) = lpos1;
	//			b2->GetdPosdqT(lc2,dpdq);
	//			break;
	//		}
	//	case 3: 
	//		{
	//
	//			break;
	//		}
	//	}
	//
	////GetBody3D(1).GetdPosdqT(lp1,dpdq);
	//double fc = -XG(1);

	//Vector ftmp(dim);
	//ftmp(lcoord1)=fc;

	//f -= dpdq*ftmp;
	
	if(UsePenaltyFormulation()) return; 

	double u = contact_direction * (GetPos1() - GetPos2());
	if(u > 0.)
	{
		return;
	}

	double test = XG(1);
	if(locelemind == 1)
	{
		//f(lcoord1) -= ComputeForce(t);
		f(lcoord1) += contact_direction * XG(1);
	}
	else
	{
		//f(lcoord2) += ComputeForce(t);
		f(lcoord2) -= contact_direction * XG(1);
	}
}
void Contact1D::PostprocessingStep()
{
	//double lastXData = GetMBS()->GetLastDataVector()(LTGdata(1));

	//if (IsSticking() == 1 && lastXData != 1) //switch from sliding to stick ==> remember actual sticking position
	//{
	//	double u = GetPos1() - GetPos2();

	//	// adjust u, such that jumps of force at slide-stick-transition are reduced
	//	double du = abs(GetSlidingFrictionForce(GetFrictionCoeff_kin()))/GetPenaltyStiffness();
	//	if(SlidingVelOLD() < 0) u+=du;
	//	else u-=du;

	//	StickingPos(u);														// origin of spring is last sticking point
	//	SlidingVelOLD(SlidingVel());							// last sliding velocity
	//	if(mbs->UO(UO_LVL_dbg1).PrintMsg())
	//	{	
	//		GetMBS()->UO() << "change slide to STICK: elementNum = " << GetOwnNum() <<", t = " << GetMBS()->GetTime() <<", u = " << StickingPos() << "\n"; 
	//	}
	//}
	//else if (IsSticking() == 0) //during sliding ==> remember last sliding velocity
	//{
	//	SlidingVelOLD(SlidingVel());		// last sliding velocity
	//	if(lastXData && (mbs->UO(UO_LVL_dbg1).PrintMsg()) )
	//	{	
	//		GetMBS()->UO() << "change stick to SLIDE: elementNum = " << GetOwnNum() <<", t = " << GetMBS()->GetTime() <<", v = " << SlidingVel() << "\n";
	//	}
	//}
}

double Contact1D::PostNewtonStep(double t)
{
	if(mode==1) return 0;		// linear elastic

	double u; // displacement
	u = contact_direction * (GetPos1() - GetPos2());

	double error=0;
	if(u > 0.)
	{
		error = 0.;
	}
	else
	{
		error = fabs(u);
	}
	return error;
}


void Contact1D::DrawElement() 
{
	Constraint::DrawElement();

	if (GetDrawSizeScalar() == 0) return;

	mbs->SetColor(colgreen);
	Vector3D p;
	Vector3D v0(0.);

	for (int i=1; i <= 2; i++)
	{
		if(elements(i))
			mbs->DrawSphere(GetElem(i).GetPosD(v0), GetDrawSizeScalar());
		else
		{
			if(nodes(i))
				mbs->DrawSphere(GetNode(i).GetPosD(), GetDrawSizeScalar());
		}
	}
};