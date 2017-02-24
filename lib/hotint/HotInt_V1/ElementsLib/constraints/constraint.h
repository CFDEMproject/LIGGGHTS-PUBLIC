//#**************************************************************
//#
//# filename:             constraint.h
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						17.October 2004
//# description:          class constraint
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
 
#ifndef CONSTRAINT__H
#define CONSTRAINT__H

#include "body2d.h"
#include "body3d.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//   Constraint   Constraint   Constraint   Constraint   Constraint
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

typedef enum {TSvalve = 0, TSRefAngle = 1, TSRefPos = 2, TSRefForce = 3, TSRefMom = 4} TConstraintIOType; // unique enumeration for defining type of output of an InputOutputElement used in a specific constraint
const Vector3D inherited_draw_dim  = Vector3D(-1,-1,-1);
const Vector3D inherited_color = Vector3D(0.3,0.8,0.3);
const Vector3D inherited_color1 = Vector3D(0.7,0.8,0.3);

class Constraint: public Element //$EDC$[beginclass,classname=Constraint,parentclassname=Element]
{
public:
	//Constraint() {mbs = NULL;};
	Constraint(MBS* mbsi);
	

	//To be overwritten in derived class:
	virtual Element* GetCopy();
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e);

		//this function assigns default values to the element variables
	virtual void ElementDefaultConstructorInitialization()
	{
		draw_dim = inherited_draw_dim;
		col = inherited_color;
		type = TConstraint;
		use_local_coordinate_system=0;
		SetPenaltyFormulation(0);
		elementname = GetElementSpec();
		spring_stiffness = 0;
		steps.SetLen(0);
		col_ext = inherited_color1;

	}

	virtual void Initialize() 
	{
	};
	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	
	// for constraints, no accessfunctions are provided usually //$ DR 2013-02-13
	virtual TKinematicsAccessFunctions GetKinematicsAccessFunctions(int mode) const
	{
		return TKinematicsAccessFunctions(TKAF_none);
	}

	// for constraints only, these are the necessary (!) access functions!	//$ DR 2013-02-13
	virtual void GetNecessaryKinematicAccessFunctions(TArray<int> &KAF, int numberOfKinematicPair)
	{
		KAF.SetLen(0);
		KAF.Add((int)(TKinematicsAccessFunctions(TKAF_none)));
	}
	
	// test if the provided KinematicAccessFunctions of an element can be used with this constraint
	// if there are differences whether the element is the first or second element of the constraint, use numberOfKinematicPair
	bool IsSuitableElement(TKinematicsAccessFunctions KAF_of_element, int numberOfKinematicPair=0);

	virtual const char* GetElementSpec() const {return "Constraint";}
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer	
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter
	virtual void EvalF(Vector& f, double t) {};  //first order equations: \dot q=f(q,t), return vector:  len(q)
	virtual void EvalG(Vector& f, double t) {UO() << "ERROR1_constraint\n";};  //evaluate constraints: len(z)
	virtual void EvalM(Matrix& m, double t) {};  //evaluate mass matrix
	virtual void EvalF2(Vector& f, double t) {}; //second order equations: M \ddot u = F2(u,\dot u,t), len(u)

	virtual int IS() const {return 1;};  //implicit (algebraic) size
	//virtual int SOS() const {return 0;}; //size of stiffness/mass matrix
	virtual int SOS() const 
	{
		if (!UsePenaltyFormulation()) 
		{
			return 0;
		}
		else
		{
			int nsos = 0;
			for (int i=1; i <= NE(); i++)
			{
				nsos += GetElem(i).SOS();
			}
			return nsos;
		}
	}; //size of stiffness/mass matrix
	virtual int SOSowned() const 
	{
		if (UsePenaltyFormulation()) return 0;
		else return SOS();
	}; //number of second order unknowns added by constraint
	virtual int ES() const {return 0;};  //size of first order explicit equations
	virtual int SS() const {return 2*SOS()+ES()+IS();};  //system size

	virtual int IsGroundJoint() const {return (elements.Length() == 1);} //needs to be overwritten if more/less elements

	//# Timeint specific derived functions: for discontinuities
	virtual double PostNewtonStep(double t) {return 0;};
	virtual void PostprocessingStep() {};
	virtual double GetError() const 
	{
		double err = 0;
		for (int i=1; i<=SS(); i++)
		{
			err += fabs(XG(i));
		}
		return err;
	};

	virtual void SetGlobalInitConditions(Vector& x_glob)
	{
		if (!UsePenaltyFormulation())
		{
			Element::SetGlobalInitConditions(x_glob);
		}
	}

  virtual void ComputationFinished() {}; //function is called when computation is finished (e.g. in order to free memory, write results, close connections, etc.)

	virtual int Dim() const {return 2;} //default value
	//
	virtual void SetInitConditions(const Vector& x0) {x_init = x0;}

	virtual void LinkToElements();

  virtual void LinkToElementsPenalty();

	virtual void BuildDependencies();

	//$ PG 2013-6-14:[ ElementsShareDOFs(): 
	// for multi-point constraints, which constrain a set of elements
	// to ground or to another set of elements (e.g. AverageConstraint),
	// elements of the same set might share the same global dofs. 
	// if this is the case, this function returns 1, else zero. 
	// this information is used in the assembling of the global system 
	// see also void Element::JacobianG(double t, Matrix& m, IVector& colref)
	// called by void MultiBodySystem::LocalJacobianG(SparseMatrix& m, Vector& x).
	virtual int ElementsShareDOFs() { return 0; }  //$ PG 2013-6-14: to be overwritten in derived classes of constraint, which act on not just one but a set of elements per kinematic pair

	//elements not doubled: only for certain elements (contact)!!!
	virtual int NE_nodouble() const {return elements.Length();}
	virtual const Element& GetElem_nodouble(int i) const
	{
		return mbs->GetElement(elements(i));
	}
	virtual Element& GetElem_nodouble(int i)
	{
		return mbs->GetElement(elements(i));
	}


	//do not try to get bodies which are not Body2D or Body3D!!!!
	virtual const Body2D& GetBody2D(int i) const {return (const Body2D&)mbs->GetElement(elements(i));}
	virtual Body2D& GetBody2D(int i) {return (Body2D&)mbs->GetElement(elements(i));}
	virtual const Body3D& GetBody3D(int i) const {return (const Body3D&)mbs->GetElement(elements(i));}
	virtual Body3D& GetBody3D(int i) {return (Body3D&)mbs->GetElement(elements(i));}

	//Function, which adds the constraint force and direction (action on the element DOF) to the element residual vector
	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f) {UO() << "ERROR2_constraint\n";};

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const 
	{
		return GetElem(1).GetRefPosD();
	}
	virtual Vector3D GetRefPos() const 
	{
		return GetElem(1).GetRefPos();
	}
	virtual Vector3D GetRefVelD() const 
	{
		return GetElem(1).GetRefVelD();
	}
	virtual Vector3D GetRefVel() const 
	{
		return GetElem(1).GetRefVel();
	}

	virtual void DrawElement() 
	{
		if (!GetMBS()->GetIOption(121)) return;
		Element::DrawElement();

		//mbs->SetColor(col);
		//mbs->DrawSphere(GetRefPosD(),draw_dim.X());
	};

	virtual void StandardJointDrawing(Matrix3D& A, Vector3D& p, Vector3D& d_dir, Vector3D& d_rot, int flag, double draw_factor);

	virtual const Vector3D& GetColExt() const {return col_ext;} 
	virtual Vector3D& GetColExt() {return col_ext;} 

	virtual void SetDrawDim(const Vector3D& drawdim) {draw_dim = drawdim;}
  
	//virtual double GetActorForce(double computation_time, int dir = 0) const {return 0;} //get actor force, with optional direction (x,y,z, etc.)
  virtual double GetActorForce(double computation_time, int dir = 0) const {return 0;} //get actor force, with optional direction (x,y,z, etc.) and computation time (for time dependent forces)

	//$ LA 2011-4-21: Get/Set draw dimensions
	virtual double GetDrawSizeScalar();
	virtual void SetDrawSizeScalar(double drawdim_scalar);

	virtual double GetDrawSizeAxisLength(){UO() << "Warning: GetDrawSizeAxisLength is not implemented in your Constraint\n"; return 0;};
	virtual void SetDrawSizeAxisLength(double drawdim_axislength){UO() << "Warning: SetDrawSizeAxisLength is not implemented in your Constraint\n";};
	virtual int GetDrawSizeResolution(){return (int)(GetMBS()->GetDOption(173));};
	//virtual void SetDrawSizeResolution(int drawdim_resolution){UO() << "Warning: SetDrawSizeResolution is not implemented in your Constraint\n"; };


	// functions to catch the input - implemented in derived class
	virtual int RespondToKey(int key) {return false;}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//access function, which determines, if constraint is in penalty mode;  //$ JG 2011-02
	virtual int UsePenaltyFormulation() const {return use_penalty_formulation;}
	//change penalty mode of constraint
	virtual void SetPenaltyFormulation(int flag) {use_penalty_formulation = flag;};
// 1D spring stifffness
	//return stiffness parameter for penalty formulation
	virtual double GetPenaltyStiffness() const {return spring_stiffness;}
	//set stiffness parameter for penalty formulation
	virtual void SetPenaltyStiffness(double stiffness) {spring_stiffness = stiffness;}
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//////////////// to be overwritten in derived class - return 3 stiffness parameters for penalty formulation
	//////////////virtual Vector3D GetPenaltyStiffness3() const {return 0;}
	//////////////virtual double GetPenaltyStiffness3(int i) const {return spring_stiffness3[i-1];}
	//////////////// to be overwritten in derived class - set 3 stiffness parameters for penalty formulation
	//////////////virtual void SetPenaltyStiffness3(Vector3D stiffness3i) {}
	//////////////// to be overwritten in derived class - return 3 damping parameters for penalty formulation
	//////////////virtual Vector3D GetPenaltyDamping3() const {return 0;}
	//////////////// to be overwritten in derived class - set 3 damping parameters for penalty formulation
	//////////////virtual void SetPenaltyDamping3(Vector3D damping3i) {}

	//access function, which determines, if constraint is applied in local coordinate system 
	virtual int UseLocalCoordinateSystem() const {return use_local_coordinate_system;}
	//change used coordinate system for constraint
	virtual void SetUseLocalCoordinateSystem(int flag) {use_local_coordinate_system = flag;}

//$ AD 2011-09-08
	// Constraint Steps
	virtual int const NSteps() const;
	virtual double GetCStepsFact(double t) const;
	virtual int GetCStepNumber(double t) const;
	virtual double GetCStepTime(double t) const;
	virtual int GetRampMode(int i) const;
	virtual double GetFinalValue(int i) const;
	virtual void SetConstraintSteps(const TArray<StepSettings>& settings);
	virtual double GetGlobalTime() const { return ((MBS*) GetMBS())->GetTime(); }

protected:
	Vector3D draw_dim; //EDC access should be done in derived element, in order to specify the necessary drawing dimensions!!! ==> e.g.
										 //EDC Vector3D draw_dim; EDC[varaccess,EDCvarname="drawing_dimensions",EDCfolder="Graphics",tooltiptext="general drawing dimensions of constraint"]
	
	//$ DR 2013-05-17: not readonly anymore but still int_bool. 
	int use_local_coordinate_system;	//$EDC$[varaccess,EDCvarname="use_local_coordinate_system",EDCfolder="Geometry",int_bool,tooltiptext="0=use global coordinates, 1=use local coordinate system of Body 1"]
																		// flag to dedicate in which coordinate system the stiffness is defined (just used in some Constraints)
																		// 0 .. use global coordinate system: k_glob = [k_x, k_y, k_z]
																		// 1 .. use local coordinate system of body 1: k_loc = [k_x1, k_y1, k_z1]
																		// n .. use local coordinate system of body n: k_loc = [k_xn, k_yn, k_zn]	//$ DR 2013-05-17
	int use_penalty_formulation; //$EDC$[varaccess,EDCvarname="use_penalty_formulation",EDCfolder="Physics",int_bool,tooltiptext="0 = use lagrange multipliers (index 3 DAE, exact), 1 = use penalty formulation (no additional equation added, approximate constraint)"]
	double spring_stiffness; //$EDC$[varaccess,EDCvarname="spring_stiffness",EDCfolder="Physics.Penalty",tooltiptext="general or penalty stiffness parameter"]

	Vector3D col_ext;

//$ AD 2011-09-08
	// Constraint Steps
	TArray<StepSettings> steps; // array containing values of loadfunction for the intervals defined in TimeInt::CSEndtimes

	//Vector3D inherited_draw_dim(-1,-1,-1);
	//Vector3D inherited_color(0.3,0.8,0.3);

	////DEPRECIATED:
	//TArray<int> ioElement_number;       // elememt number of InputOutputElement
	//TArray<int> ioElement_actionIndex;  // unique type number (e.g.: valve position for hydraulic, ...)
	//TArray<int> ioElement_outputNumber; // elememt number of InputOutputElement
	
	//EDC	int draw_element; //$EDC$[varaccess,remove,EDCvarname="show_element",EDCfolder="Graphics"]
	//EDC	int draw_element; //$EDC$[varaccess,EDCvarname="show_connector",EDCfolder="Graphics",int_bool,tooltiptext="Flag to draw connector"]

};//$EDC$[endclass,Constraint]


#endif
