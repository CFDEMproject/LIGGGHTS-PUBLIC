//#**************************************************************
//#
//# filename:             constraint.cpp
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						17.October 2004
//# description:          class constraint
//#                       
//# remarks:						  
//#
///# Copyright (c) 2003-2013 Johannes Gerstmayr, Linz Center of Mechatronics GmbH, Austrian
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
#include "body2d.h"
#include "body3d.h"
#include "constraint.h"

Constraint::Constraint(MBS* mbsi):Element(mbsi), steps(0)
{
	ElementDefaultConstructorInitialization();
};

//To be overwritten in derived class:
Element* Constraint::GetCopy()
{
	Element* ec = new Constraint(*this);
	return ec;
}
//To be overwritten in derived class:
void Constraint::CopyFrom(const Element& e)
{
	Element::CopyFrom(e);
	const Constraint& ce = (const Constraint&)e;
	draw_dim = ce.draw_dim;
	use_local_coordinate_system = ce.use_local_coordinate_system;
	spring_stiffness = ce.spring_stiffness;
	use_penalty_formulation = ce.use_penalty_formulation;
	steps = ce.steps;
	col_ext = ce.col_ext;
}


void Constraint::BuildDependencies() 
{
	Element::BuildDependencies();
}

int Constraint::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = Element::CheckConsistency(errorstr);
	if (rv) return rv;
	
	return rv;
}

// test if the provided KinematicAccessFunctions of an element can be used with this constraint
// if there are differences whether the element is the first or second element of the constraint, use numberOfKinematicPair
bool Constraint::IsSuitableElement(TKinematicsAccessFunctions KAF_of_element, int numberOfKinematicPair)
{
#ifdef __RELEASE_VERSION__
	//$ DR 2013-02-13
	TArray<int> KAF_of_joint;
	GetNecessaryKinematicAccessFunctions(KAF_of_joint, numberOfKinematicPair);

	int max_bit = (int)(TKinematicsAccessFunctions(TKAF_maximum)/2);	//the highest available KinematixAccessFunction

	for(int i=1; i<=KAF_of_joint.Length();i++)
	{
		bool suitable=1;
		int KAF_i_of_joint = KAF_of_joint(i);
		for(int bit = 1; bit <= max_bit; bit=bit*2)
		{
			if(KAF_i_of_joint &bit)	// this AccessFunction is needed by the joint
			{
				if(!(KAF_of_element &bit)) // the AccessFunction is not provided by the element
				{
					suitable = 0;
					break;				//break the inner for-loop
				}
			}
		}
		if(suitable) return 1;	//found one suitable combination of AccessFunctions
	}

	return false;	// the element can be used with the constraint
#else
	return true;	// if you are not in release version, than this check is not performed
#endif
}

void Constraint::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Element::GetElementData(edc);

	ElementData ed;

	//SDActor and SDRotActor do not inherit the following commands:

	//element numbers should be managed by constraints, because elements can appear in different folders
	//Vector elnum(NE());
	//for (int i=1; i <= NE(); i++) {elnum(i) = GetElnum(i);}
	//if (NE())
	//{
	//	ed.SetVector(elnum.GetVecPtr(), elnum.Length(), "Constraint_element_numbers"); ed.SetToolTipText("Only valid element numbers permitted!"); edc.Add(ed);
	//}
}

int Constraint::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Element::SetElementData(edc);

	//element numbers should be managed by constraints, because elements can appear in different folders
	//int pos;

	//pos = edc.Find("Constraint_element_numbers");
	//if (pos)
	//{
	//	const ElementData& ed = edc.Get(pos);
	//	if (ed.IsVector())
	//	{
	//		for (int i=1; i <= ed.GetVectorLen(); i++)
	//		{
	//			int en = (int)ed.GetVectorVal(i);
	//			if (en <= 0 || en > GetMBS()->GetNElements()) 
	//			{
	//				GetMBS()->EDCError(mystr("Constraint: Element number ")+mystr(en)+mystr(" is out of range"));
	//				en = 1;
	//			}
	//			SetElnum(i, en);
	//		}
	//	}
	//} //not all controllers have element numbers!
	return rv;
}


void Constraint::LinkToElements() 
{
	if (UsePenaltyFormulation()) //in penalty formulation, the SOS-DOF are hired from constrained elements
	{
		LinkToElementsPenalty();
	}
	else
	{
		for (int i = 1; i <= NE(); i++) //$ DR 2012-11-02: changed elements.Length() to NE(), because some constraints add element(2) = 0 for ground joints, NE() returns the correct value for these constraints 
		{
			GetElem(i).AddConstraint(this,i);
		}
	}
}

void Constraint::LinkToElementsPenalty()
{
	if (IS()!=0 && SOS()!=0)
	{
		UO(UO_LVL_err) << "Error: Constraint::LinkToElementsPenalty() is not possible for mixed penalty/Lagrange elements\n";
	}
	LTGreset();

	// add all SOS dofs from the elements
	//Position(first SOS) 
	for (int k=1; k <= NE(); k++)
	{
		for (int i=1; i <= GetElem(k).SOS(); i++)
		{
			AddLTG(GetElem(k).LTG(i));
		}
	}
	//and Velocity (second SOS):
	for (int k=1; k <= NE(); k++)
	{
		for (int i=1; i <= GetElem(k).SOS(); i++)
		{
			AddLTG(GetElem(k).LTG(i+GetElem(k).SOS()));
		}
	}
}

//$ LA 2011-04-21: Get/Set draw dimensions
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//				get- and set-functions for drawing
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double Constraint::GetDrawSizeScalar()
{
	if(draw_dim.X()==-1)
		return GetMBS()->GetDOption(171);			
	else
		return draw_dim.X();
}
void Constraint::SetDrawSizeScalar(double drawdim_scalar)
{
	draw_dim.X()=drawdim_scalar;
}

void Constraint::StandardJointDrawing(Matrix3D& A, Vector3D& p, Vector3D& d_dir, Vector3D& d_rot, int flag, double draw_factor) 
{
	//double draw_factor = GetDrawSizeScalar();
	Vector3D draw_dir;

	for (int i=1; i<=3; i++)
	{
		draw_dir = draw_factor*Vector3D(A.Get0(0,i-1),A.Get0(1,i-1),A.Get0(2,i-1));
		if(d_dir(i)!=0)
		{
			mbs->SetColor(GetCol());
			mbs->DrawCone(p + flag*draw_dir,p,draw_factor*1.1,6,1);
		}
		if(d_rot(i)!=0)
		{
			//GetMBS()->ChooseColor(0.0f,0.6f,0.0f);
			mbs->SetColor(GetColExt());
			mbs->DrawCone(p + 2*flag*draw_dir,p,draw_factor,6,1);
			mbs->DrawCone(p + 4*flag*draw_dir,p + 2*flag*draw_dir,draw_factor,6,1);
			//mbs->SetColor(GetCol());
		}
	}
}


//$ AD 2011-09-07: Constraint Steps
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//              functions for Constraint - Steps
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int const Constraint::NSteps() const 
{
	return steps.Length();
}

// this function is an exact copy of MBSLoad::GetConstraintStepFact
double Constraint::GetCStepsFact(double t) const
{
	if(t <= 0) 
		return 0.;
	
  double steptime = GetCStepTime(t);    // step number (computation step number from TimeInt)
	int stepnumber = GetCStepNumber(t);   // step time   (computation step time from TimeInt)

// return value when time is after last computation step endtime...
	if(stepnumber < 0 && t > 0.) 
		return GetFinalValue(steps.Length());

// return value when too few load steps are defined
	if(stepnumber > NSteps()) 
		return GetFinalValue(steps.Length());


	double finalvalue = GetFinalValue(stepnumber); // value at the end of this interval
// return value for STEP ( finalvalue is valid for entire interval )
	if(GetRampMode(stepnumber) == TCSRStep)
		return finalvalue; 

	double startvalue = 1.;  // value at beginning of this interval, =1. means that all constraints start ON !!!
	if (stepnumber != 1) 
    startvalue = GetFinalValue(stepnumber-1);
// return value for LINEAR ( linear interpolation between beginning and end of the interval )
  if(GetRampMode(stepnumber) == TCSRLinear)
	{
		double stepfactor = startvalue + (finalvalue-startvalue)*steptime;
		return stepfactor;
	}

	const double eps_t = 1e-12;
// return value for EXPONENTIAL ( eponential interpolation between beginning and end of the interval )
	if(steptime < eps_t) 
		return startvalue;
	else if(steptime == 1.)
		return finalvalue;
	else
	{
		double factor_growth_decay = 1.;
		if(finalvalue < startvalue) 
			factor_growth_decay=-1;

// exponential interpolation:  y(x) = y(0) * (y(1) / y(0))^x ... with y(0) and y(1) may not be 0.
// but this always works
// with y(x) = a*e^+x + d , x[0..1] --> y(0) = a+d, y(1) = ae+d    --> a = (y(1)-y(0)) / (e-1),    d = y(0) - a  (growth)
// with y(x) = a*e^-x + d , x[0..1] --> y(0) = a+d, y(1) = ae^-1+d --> a = (y(1)-y(0)) / (e^-1-1), d = y(0) - a  (decay)
		double a = (finalvalue-startvalue) / (exp(factor_growth_decay)-1.);
		double d = startvalue - a;
		double stepfactor = a * exp(steptime*factor_growth_decay) + d;

////// exponential interpolation Mk.II  factor 10 in 1/10 th step (  10^0 for 1, 10^-1 for 0.9, 10^-4 for 0.5, cutoff 10^-20,...
////// results:  y(x)  = (y(1)-y(0)) * max(10^-(10*(1-x)), cutoff) + y(0)    (growth)
//////	         y(x)  = (y(1)-y(0)) * max(10^-(10*(x)), cutoff) + y(0)
////	  const double eps_factor = 1e-20; // <-- minimum for 
////    double exponent = (1.-steptime) * 10.;
////		if (finalvalue < startvalue)
////			exponent = (-steptime) * 10;
////		double stepfactor = (finalvalue - startvalue) * Maximum(eps_factor, pow(10.,exponent)) + startvalue;


		return stepfactor;
	}
  return 0.;
}
 
int Constraint::GetCStepNumber(double t) const
{
	return ((MBS*) GetMBS()) ->GetCSNumber(t); // do indeed calculate
//	return GetMBS()->ComputationStepNumber(); // access-function: gives wrong values for t=0 during assemble
}

double Constraint::GetCStepTime(double t) const
{ 
	return ((MBS*) GetMBS()) ->GetCSTime(t); // do indeed calculate
//	return GetMBS()->ComputationStepTime(); // access-function: gives wrong values for t=0 during assemble
}

int Constraint::GetRampMode(int i) const
{ 
	return steps(i).GetRampMode(); 
}

double Constraint::GetFinalValue(int i) const
{ 
	return steps(i).GetLoadFactor(); 
}

void Constraint::SetConstraintSteps(const TArray<StepSettings>& settings)
{
	steps = settings;
}
