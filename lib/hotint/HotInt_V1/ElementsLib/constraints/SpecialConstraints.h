//#**************************************************************
//#
//# filename:             SpecialConstraints.h
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
 
#ifndef SPECIALCONSTRAINT__H
#define SPECIALCONSTRAINT__H

#include "KinematicPairs.h"

class GeomGround;


//EXAMPLE OF UNIFIED CONSTRAINT DESCRIPTION (UED):
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: CoordConstraint (Constraint)
// short description: Direct coordinate constraint based on generalized coordiantes
// available formulations: Lagrange multiplier, penalty
// type of constraint equation: linear
// index formulation available: 2, 3
// ground constraint: yes
// element-element constraint: yes
// development status: incomplete (EvalF2 for penalty missing)
// long description:	Constraint for constraining a specific generalized (local) coordinate of one or two bodies; 
//										in the case of one body ==> ground joint; in the case of two bodies ==> direct coordinate constraint between two bodies
// class variables:
//   - spring_stiffness: specifies the penalty factor / spring in case of penalty formulation only!
//**end(ued)**
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class CoordConstraint: public Constraint //$EDC$[beginclass,classname=CoordConstraint,parentclassname=Constraint,addelementtype=TAEconstraint+TAEspecial_connector,addelementtypename=CoordinateConstraint,texdescription="
//The CoordinateConstraint constrains two elements by constraining a single coordinate of each element, e.g. the x-displacement of two different elements. If the second element number is zero, a ground
//joint can be realized. The CoordinateConstraint uses the lagrange multiplier formulation by default, which means that there is no constraint violation at all. For static problems, the lagrange multiplier 
//constraint formulation is applied directly, by adding the kinematical conditions to the nonlinear system equations. In dynamic (time dependent) simulations, the constraint is solved on the position (displacement) level
//with index 3 solvers and on the velocity level with index 2 solvers. 
//Alternatively, the penalty formulation can be used, which means that a certain (very high) spring stiffness is used instead of lagrange multipliers. Thus, no additional equation is added, however, the system
//equations may become unsolvable stiff (ill conditioned) in case of static problems; for dynamical problems, the very high stiffness might lead to high-frequency oscillations, inaccurate solutions or no convergence.",
//texdescriptionEquations="\textbf{Lagrange formulation:}\\ \\
//position constraint (index 3 solver) \\ 
//2 elements (coordinate to coordinate): 
//$C=k\left(q^{el1}_i-q^{el1}_{i,0}\right)-\left(q^{el2}_j-q^{el2}_{j,0}\right)-d=0$ \\
//1 element (coordinate to ground): 
//$C=k\left(q^{el1}_i-q^{el1}_{i,0}\right)-d=0$ \\ \\
//velocity constraint - index reduction (index 2 solver) \\
//2 elements (coordinate to coordinate): 
//$C=k\,\dot{q}^{el1}_i-\dot{q}^{el2}_j=0$, \\
//1 element (coordinate to ground): 
//$C=\dot{q}^{el1}_i=0$ \\ \\
//Langrange multiplier \\
//$\frac{\partial{C}}{\partial{\mathbf{q}^{el1}}}^T= \left[0 \, ......... \, 0 \, , k \, , 0 \, ... \, 0\right]$ ... with $k$ at index i \\
//$\frac{\partial{C}}{\partial{\mathbf{q}^{el2}}}^T= \left[0 \, ... \, 0 \, , -1 \, , 0 \, ......... \, 0\right]$ ... with $-1$ at index j \\ \\
//\textbf{Penalty formulation:}\\ \\
//2 elements (coordinate to coordinate): 
//$f=S_P\left(k\left(q^{el1}_i-q^{el1}_{i,0}\right)-\left(q^{el2}_j-q^{el2}_{j,0}\right)-d\right)+D_P\left(k \,\dot{q}^{el1}_i-\dot{q}^{el2}_j\right)$ \\
//1 element (coordinate to ground): 
//$f=S_P\left(k\left(q^{el1}_i-q^{el1}_{i,0}\right))-d\right)+D_P \, k \, \dot{q}^{el1}_i$ \\ \\
//\textbf{Description:} \\
//$k$ ... coordinate gain factor \\
//$d$ ... coordinate offset (for index 2 solvers not used) \\
//$q^{el1}_i$ ... $i^{th}$ coordinate of element 1 \\
//$q^{el1}_{i,0}=q^{el1}_i\left(t=0\right)$ ... $i^{th}$ coordinate of element 1 at initialization \\
//$q^{el2}_j$ ... $j^{th}$ coordinate of element 2 \\
//$q^{el2}_{j,0}=q^{el2}_j\left(t=0\right)$ ... $j^{th}$ coordinate of element 2 at initialization \\
//$S_P$ ... spring stiffness \\
//$D_P$ ... damping \\
//$C$ ... Lagrange equation \\
//$f$ ... force vector (penalty formulation)",
//modus="{coordinate to ground}{Coordinate2.element\_number AND Coordinate2.local\_coordinate have to be equal to 0}",
//modus="{coordinate to coordinate}{Coordinate2.element\_number and/or Coordinate2.local\_coordinate must not be equal to 0}",
//modus="{Lagrange}{For Physics.use\_penalty\_formulation = 0 no stiffness parameter is used.}",
//modus="{relative or absolute to initial values}{Only important for max index 3 solvers. \newline If relative\_to\_inital\_values is set to 1: Equation above is used. \newline If set to 0: Simplified equation is used ($q^{el1}_{i,0} = q^{el2}_{j,0} = 0$).}",example="CoordinateConstraint.txt"]

{
public:
	//CoordConstraint():loccoords() {mbs = NULL;};
	CoordConstraint(MBS* mbsi):Constraint(mbsi), loccoords() 
	{	ElementDefaultConstructorInitialization();};

	CoordConstraint(const CoordConstraint& ct): Constraint(ct.mbs), loccoords(0)
	{
		CopyFrom(ct);
	};

	CoordConstraint(MBS* mbsi, int en1, int lc1, double r=-1):Constraint(mbsi), loccoords()
	{	
		SetCoordConstraint(en1, lc1, r);
	};

	CoordConstraint(MBS* mbsi, int en1, int en2, int lc1, int lc2, double r=-1):Constraint(mbsi), loccoords() 
	{	
		SetCoordConstraint(en1, en2, lc1, lc2, r);
	}

	// penalty formulation instead of lagrange parameter (=> ground joint with spring)
	CoordConstraint(MBS* mbsi, int en1, int lc1, double r, double spring_stiffnessi):Constraint(mbsi), loccoords() 
	{	
		SetCoordConstraintSpring(en1, lc1, r, spring_stiffnessi);
	}

	// set functions
	virtual void SetCoordConstraint(int en1, int lc1, double r=-1) // $ MSax 2013-04-22 removed coord_gain_factor from set function
	{
		ElementDefaultConstructorInitialization();

		AddElementCoord(en1, lc1);
		if(r!=-1)
			SetDrawSizeScalar(r); // draw_dim.X() = r;

	};
	virtual void SetCoordConstraint(int en1, int en2, int lc1, int lc2, double r=-1) // $ MSax 2013-04-22 removed coord_gain_factor from set function
	{	
		ElementDefaultConstructorInitialization();

		AddElementCoord(en1, lc1);
		AddElementCoord(en2, lc2);
		if(r!=-1)
			SetDrawSizeScalar(r); // draw_dim.X() = r;

	}
	virtual void SetCoordConstraintSpring(int en1, int lc1, double r, double spring_stiffnessI) // $ MSax 2013-04-22 removed coord_gain_factor and coord_offset from set function
	{
		ElementDefaultConstructorInitialization();

		SetPenaltyFormulation(1); //JG
		SetPenaltyStiffness(spring_stiffnessI);

		AddElementCoord(en1, lc1);
		SetDrawSizeScalar(r); //draw_dim.X() = r;
		x_init = Vector(SS()); //$JG2012-02 this is necessary, because in penalty constraint, it changes from default

	}
	virtual void SetCoordConstraintSpring(int en1, int en2, int lc1,  int lc2, double r, double spring_stiffnessI) // $ MSax 2013-04-22 removed coord_gain_factor from set function
	{
		SetPenaltyFormulation(1); //JG
		SetPenaltyStiffness(spring_stiffnessI);

		//x_init = Vector(SS());
		AddElementCoord(en1, lc1);
		AddElementCoord(en2, lc2);
		SetDrawSizeScalar(r); //draw_dim.X() = r;

		x_init = Vector(SS()); //$JG2012-02 this is necessary, because in penalty constraint, it changes from default

	}

	virtual void SetCoordOffset(double coord_offsetI)
	{
		coord_offset = coord_offsetI;
	}

	virtual void SetCoordGainFactor(double coord_gain_factorI)
	{
		coord_gain_factor = coord_gain_factorI;
	}

	//virtual void Initialize();

	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new CoordConstraint(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const CoordConstraint& ce = (const CoordConstraint&)e;
		loccoords = ce.loccoords;
		coord_offset = ce.coord_offset;
		coord_gain_factor = ce.coord_gain_factor;
		relative_to_inital_values = ce.relative_to_inital_values;
		damping_coeff = ce.damping_coeff;
	}
	
	virtual void ElementDefaultConstructorInitialization()
	{
		SetPenaltyFormulation(0); //JG
		loccoords.SetLen(0);
		elements.SetLen(0);
		elementname = GetElementSpec();
		x_init = Vector(SS());
		coord_offset = 0;
		coord_gain_factor = 1;
		relative_to_inital_values = 1;
		damping_coeff = 0;
	}

	virtual void AddElementCoord(int en, int loccoord)
	{
		AddElement(en);
		loccoords.Add(loccoord);
	}

	virtual const char* GetElementSpec() const {return "CoordinateConstraint";}
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data according to ElementDataContainer
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	//evaluate constraints: len(z)
	virtual void EvalG(Vector& f, double t);
	virtual void EvalF2(Vector& f, double t); //second order equations for constraints: M \ddot u = F2(u,\dot u,t), len(u)

	//To be replaced in derived class
	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	virtual int IS() const {return 1-UsePenaltyFormulation();};  //implicit (algebraic) size (RL)
	virtual int Dim() const 
	{
		if (NE()>=1 && GetElnum(1) != 0) 
		{
			return GetElem(1).Dim();
		} else 
		{
			return 3; //default dimension
		}
	}  //has 3D position???

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//DEPRECIATED FUNCTIONS: //$ JG 2011-02
	//virtual void SetPenalty(int flag) {SetPenaltyFormulation(flag);} // use penalty coordinate constraint (RL)
	//virtual int GetPenalty() {return UsePenaltyFormulation();} // use penalty coordinate constraint (RL)

	//virtual void SetSpringStiffness(double val){SetPenaltyStiffness(val);} // use penalty coordinate constraint (RL)
	//virtual double GetSpringStiffness(){return GetPenaltyStiffness();} // use penalty coordinate constraint (RL)
	//END DEPRECIATED FUNCTIONS
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//do not add initial conditions for 2nd order size DOF, because they are set in constrained elements
	virtual void SetGlobalInitConditions(Vector& x_glob) {}

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const 
	{
		return GetElem(1).GetRefPosD();
	}
	
	virtual void GetGlobalConstrainedDOFs(TArray<int>& dofs) const;

	virtual void DrawElement();

	//compute constraint force, for penalty or Lagrange approach
	virtual double ComputeForce(double t) const;

	virtual int WriteSingleElementData(const ReadWriteElementDataVariableType& RWdata) 
	{
		// call base class routine ( not required in Element )	
		int rv = Constraint::WriteSingleElementData(RWdata);
		if (rv == 1) return 1;
		//manual things to write

		if(RWdata.variable_name.CStrCompare("Connector.coordinate_offset"))
		{
			coord_offset = RWdata.value;
		}

		if(RWdata.variable_name.CStrCompare("Connector.gain_factor"))
		{
			coord_gain_factor = RWdata.value;
		}

		return WriteSingleElementDataAuto(RWdata);
	}

	virtual int ReadSingleElementData(ReadWriteElementDataVariableType& RWdata) 		
	{
		// call base class routine 	
		int rv = Constraint::ReadSingleElementData(RWdata);
		if (rv == 1) return 1;

		// manual things to read  
		if(RWdata.variable_name.CStrCompare("Connector.constraint_force"))
		{
			//if (RWdata.comp1 >= 0 && RWdata.comp1 <= 1) //range check
			//{
			RWdata.value = ComputeForce(mbs->GetTime());//ConstraintForce(mbs->GetTime(),RWdata.comp1); 
			return 1; 
			//}
		}

		else if(RWdata.variable_name.CStrCompare("Connector.coordinate_difference"))
		{
			//if (RWdata.comp1 == 1) //range check
			//{
			RWdata.value = /*XG(loccoords(1));*/ GetElem(1).XG(loccoords(1));
			if (elements.Length()==2) 
			{
				RWdata.value -= /*XG(loccoords(2));*/ GetElem(2).XG(loccoords(2));
			}
			return 1; 
			//}
		}
		else if(RWdata.variable_name.CStrCompare("Connector.coordinate_offset"))
		{
			RWdata.value = coord_offset; 
			return 1; 
		}

		else if(RWdata.variable_name.CStrCompare("Connector.gain_factor"))
		{
			RWdata.value = coord_gain_factor; 
			return 1; 
		}

		else return -2; 
	}

	//virtual double ConstraintForce(double computation_time, int dir=0) const //get actor force, with optional direction (x,y,z, etc.)
	//{
	//	if (dir == 0) return abs(ComputeForce(computation_time));
	//	if (dir > 1) assert(0 && "CoordConstraint::ConstraintForce");

	//	return ComputeForce(computation_time); // dir == 1
	//}

	virtual int GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables)
	{
		// no parent class, flush the array here ( and only here )
		available_variables.Flush();

		// Automatic entries for this class 
		Element::GetAvailableSpecialValuesAuto(available_variables);

		// Manual READ entries for this class
		available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_force",0,0,0.,mystr("internal force of connector"), TRWElementDataRead)); //test only!!!!

		available_variables.Add(ReadWriteElementDataVariableType("Connector.coordinate_difference",0,0,0.,mystr("difference between the coordinates"), TRWElementDataRead)); //test only!!!!

		// Manual WRITE/READWRITE entries for this class
		available_variables.Add(ReadWriteElementDataVariableType("Connector.coordinate_offset",0,0,0.,mystr("coordiante offset for CoordinateConstraint (w.r.t. ground or between two element coordinates); offset is ignored for Index 2 (setting of time integration) velocity level constraint"), TRWElementDataReadWrite)); //test only!!!!

		available_variables.Add(ReadWriteElementDataVariableType("Connector.gain_factor",0,0,0.,mystr("coordiante gain factor for CoordinateConstraint"), TRWElementDataReadWrite)); //test only!!!!

		return 0;
	}


	//get a component of the constraint force ==> for sensors
	//deprecated function for old sensors
	virtual double GetActorForce(double computation_time, int dir=0) const {return ComputeForce(computation_time);} //get actor force, with optional direction (x,y,z, etc.)

	virtual void SetPenaltyDamping(double damping_factor) {damping_coeff=damping_factor;};
	virtual double GetPenaltyDamping() const {return damping_coeff;};

	//see users Documentation on the action of this flag; JG
	virtual void SetRelativeToInitialValues(int flag) {relative_to_inital_values = flag;}
	virtual int GetRelativeToInitialValues() {return relative_to_inital_values;}

protected:
	TArray<int> loccoords;
	double coord_offset;						//$EDC$[varaccess,EDCvarname="coord_offset",tooltiptext="coordinate offset d, see documentation section equation"]
	double coord_gain_factor;				//$EDC$[varaccess,EDCvarname="coord_gain_factor",tooltiptext="coordinate gain factor k, see documentation section equation"]
	int relative_to_inital_values;  //$EDC$[varaccess,EDCvarname="relative_to_inital_values",int_bool,tooltiptext="flag == 1: full equation is used, see documentation; flag == 0: the init state values qi0 and qj0 are neglected."]
	double damping_coeff;						//$EDC$[varaccess,EDCvarname="damping",minval=0,EDCfolder="Physics.Penalty",tooltiptext="damping coefficient Dp for viscous damping"]

	//EDC	int use_local_coordinate_system;	//$EDC$[varaccess,remove,EDCvarname="use_local_coordinate_system",EDCfolder="Geometry"]

	//EDC double spring_stiffness; //$EDC$[varaccess,remove,EDCvarname="spring_stiffness",EDCfolder="Physics.Penalty"]
	//EDC double spring_stiffness; //$EDC$[varaccess,EDCvarname="spring_stiffness",EDCfolder="Physics.Penalty",tooltiptext="general or penalty stiffness parameter Sp"]

};//$EDC$[endclass,CoordConstraint]

class VelCoordConstraint: public CoordConstraint //$EDC$[beginclass,classname=VelCoordConstraint,parentclassname=CoordConstraint,addelementtype=TAEconstraint+TAEspecial_connector,addelementtypename=VelocityCoordinateConstraint,texdescription="Similar to CoordinateConstraint. Lagrangian constraint implemented for index 3 and index 2 solvers. A penalty formulation is also implemented.",
//texdescriptionEquations="\textbf{Lagrange formulation:}\\ \\
//velocity constraint (index 2 and 3 solvers) \\ 
//2 elements (coordinate to coordinate): 
//$C=k\left(\dot{q}^{el1}_i-\dot{q}^{el1}_{i,0}\right)-\left(\dot{q}^{el2}_j-\dot{q}^{el2}_{j,0}\right)-d=0$, \\
//1 element (coordinate to ground): 
//$C=k\left(\dot{q}^{el1}_i-\dot{q}^{el1}_{i,0}\right)-d=0$ \\ \\
//Langrange multiplier \\
//$\frac{\partial{C}}{\partial{\mathbf{\dot{q}}^{el1}}}^T= \left[0 \, ......... \, 0 \, , k \, , 0 \, ... \, 0\right]$ ... with $k$ at index i \\
//$\frac{\partial{C}}{\partial{\mathbf{\dot{q}}^{el2}}}^T= \left[0 \, ... \, 0 \, , -1 \, , 0 \, ......... \, 0\right]$ ... with $-1$ at index j \\ \\
//\textbf{Penalty formulation:}\\ \\
//2 elements (coordinate to coordinate): 
//$f=S_P\left(k\left(\dot{q}^{el1}_i-\dot{q}^{el1}_{i,0}\right)-\left(\dot{q}^{el2}_j-\dot{q}^{el2}_{j,0}\right)-d\right)$ \\
//1 element (coordinate to ground): 
//$f=S_P\left(k\left(\dot{q}^{el1}_i-\dot{q}^{el1}_{i,0}\right))-d\right)$ \\ \\
//\textbf{Description:} \\
//$k$ ... coordinate velocity gain factor \\
//$d$ ... coordinate velocity offset \\
//$\dot{q}^{el1}_i$ ... $i^{th}$ coordinate velocity of element 1 \\
//$\dot{q}^{el1}_{i,0}=\dot{q}^{el1}_i\left(t=0\right)$ ... $i^{th}$ coordinate velocity of element 1 at initialization \\
//$\dot{q}^{el2}_j$ ... $j^{th}$ coordinate velocity of element 2 \\
//$\dot{q}^{el2}_{j,0}=\dot{q}^{el2}_j\left(t=0\right)$ ... $j^{th}$ coordinate velocity of element 2 at initialization \\
//$S_P$ ... spring stiffness \\
//$C$ ... Lagrange equation \\
//$f$ ... force vector (penalty formulation)",
//modus="{coordinate to ground}{Coordinate2.element\_number AND Coordinate2.local\_coordinate have to be equal to 0}",
//modus="{coordinate to coordinate}{Coordinate2.element\_number and/or Coordinate2.local\_coordinate must not be equal to 0}",
//modus="{relative or absolute to initial values}{If relative\_to\_inital\_values is set to 1: Equation above is used. \newline If set to 0: Simplified equation is used ($\dot{q}^{el1}_{i,0} = \dot{q}^{el2}_{j,0} = 0$).}",example="VelocityCoordinateConstraint.txt"]

{
public:
	VelCoordConstraint(MBS* mbsi):CoordConstraint(mbsi)
	{	ElementDefaultConstructorInitialization();};

	VelCoordConstraint(const VelCoordConstraint& ct): CoordConstraint(ct.mbs)
	{
		CopyFrom(ct);
	};

	// penalty formulation instead of lagrange parameter (=> ground joint with spring)
	VelCoordConstraint(MBS* mbsi, int en1, int lc1, double r, double spring_stiffnessi):CoordConstraint(mbsi)
	{	
		mbs->UO() << "WARNING: VelCoordConstraint::Wrong constructor used. Penalty formulation not implemented.\n";		
		assert(0);
	}

private:
	virtual void SetCoordConstraintSpring(int en1, int lc1, double r, double spring_stiffnessI) {assert(0);} // $ MSax 2013-04-23 not used in derived class (no penalty formulation)
	virtual void SetCoordConstraintSpring(int en1, int en2, int lc1,  int lc2, double r, double spring_stiffnessI){assert(0);} // $ MSax 2013-04-23 not used in derived class (no penalty formulation)

public:
	virtual void EvalG(Vector& f, double t);
	//virtual void EvalF2(Vector& f, double t); //second order equations for constraints: M \ddot u = F2(u,\dot u,t), len(u)
	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	virtual const char* GetElementSpec() const {return "VelocityCoordinateConstraint";}

	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute

	virtual Element* GetCopy()
	{
		Element* ec = new VelCoordConstraint(*this);
		return ec;
	}

	virtual void CopyFrom(const Element& e)
	{
		CoordConstraint::CopyFrom(e);
		const VelCoordConstraint& ce = (const VelCoordConstraint&)e;
	}
	
	virtual void ElementDefaultConstructorInitialization()
	{
		CoordConstraint::ElementDefaultConstructorInitialization();
		elementname = GetElementSpec();
	}

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data according to ElementDataContainer
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual double ComputeForce(double t) const;
protected:
	//EDC int use_penalty_formulation; //!EDC$[varaccess,remove,EDCvarname="use_penalty_formulation",EDCfolder="Physics"]
	//EDC double spring_stiffness; //!EDC$[varaccess,remove,EDCvarname="spring_stiffness",EDCfolder="Physics.Penalty"]
	//EDC double damping_coeff;						//$EDC$[varaccess,remove,EDCvarname="damping",EDCfolder="Physics.Penalty"]
	//EDC int relative_to_inital_values;  //$EDC$[varaccess,remove,EDCvarname="relative_to_inital_values"]
	//EDC int relative_to_inital_values;  //$EDC$[varaccess,EDCvarname="relative_to_inital_values",int_bool,tooltiptext="flag == 1: full equation is used, see documentation; flag == 0: the init state derivatives d(qi0)/dt and d(qj0)/dt are neglected."]

};//$EDC$[endclass,VelCoordConstraint]

class ArcLengthConstraint : public CoordConstraint
{
	public:

	ArcLengthConstraint(MBS* mbsi) : CoordConstraint(mbsi) 
	{	
	};

	//ArcLengthConstraint(MBS* mbsi, double loadi, int elem, int loccoordi) : CoordConstraint(mbsi) 
	//{	
	//	load = loadi;
	//	loccoord = loccoordi;
	//	AddElement(elem);
	//};

	ArcLengthConstraint(const ArcLengthConstraint& ct): CoordConstraint(ct.mbs)
	{
		CopyFrom(ct);
	};

	~ArcLengthConstraint()
	{
	};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new ArcLengthConstraint(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		CoordConstraint::CopyFrom(e);
		const ArcLengthConstraint& ce = (const ArcLengthConstraint&)e;
		/*this->*/
		targetload = ce.targetload;
		//loccoord= ce.loccoord;
		//elem= ce.elem;
	}

	virtual const char* GetElementSpec() const {return "ArcLengthConstraint";}

	virtual void SetArcLengthConstraint(double load, int elem, int loccoord)
	{
		SetPenaltyFormulation(0); //JG
		loccoords.SetLen(0);
		elements.SetLen(0);
		x_init = Vector(SS());
		AddElementCoord(elem, loccoord);
		targetload = load;

		elementname = GetElementSpec();
	}

	// virtual void SetArcLengthConstraint();

	virtual int IS() const { return 1; } //implicit size -> number of constraint equations

	virtual void EvalG(Vector& f, double t);
	virtual void AddElementCqTLambda(int locelemind, Vector& f);

protected:
	mutable double targetload;
	//int			loccoord;
	//int			elem;
	//double	load;
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Prescribed coordinate constraint
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class PrescribedCoordConstraint: public CoordConstraint
{
public:
	//PrescribedCoordConstraint():loccoords() {mbs = NULL;};
	PrescribedCoordConstraint(MBS* mbsi):CoordConstraint(mbsi) {};

	PrescribedCoordConstraint(MBS* mbsi, int en1, int lc1, double r=-1):CoordConstraint(mbsi)
	{	
		int ss = SS();
		x_init = Vector(SS());
		AddElementCoord(en1, lc1);
		if(r!=-1)
			SetDrawSizeScalar(r); //draw_dim.X() = r;

		mathfunc.SetConstant(0);
		constraintmode = 0;

		elementname = GetElementSpec();

		endtime = -1.;
	};

	PrescribedCoordConstraint(MBS* mbsi, int en1, int en2, int lc1, int lc2, double r=-1):CoordConstraint(mbsi)
	{	
		x_init = Vector(SS());
		AddElementCoord(en1, lc1);
		AddElementCoord(en2, lc2);
		if (r!=-1)
			SetDrawSizeScalar(r); //draw_dim.X() = r;

		mathfunc.SetConstant(0);
		constraintmode = 0;

		elementname = GetElementSpec();
		
		endtime = -1.;
	};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new PrescribedCoordConstraint(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		CoordConstraint::CopyFrom(e);
		const PrescribedCoordConstraint& ce = (const PrescribedCoordConstraint&)e;

		mathfunc = ce.mathfunc;
		constraintmode = ce.constraintmode;
		endtime = ce.endtime;
	}
	virtual const char* GetElementSpec() const {return "PrescribedCoordConstraint";}
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	//evaluate constraints: len(z)
	virtual void EvalG(Vector& f, double t);
	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	virtual void SetPrescribedHarmonic(double omega, double phase, double amplitude)
	{
		mathfunc.SetHarmonic(Vector(omega), Vector(phase), Vector(amplitude));
	}
	virtual MathFunction& GetMathFunction() {return mathfunc;}
	virtual const MathFunction& GetMathFunction() const {return mathfunc;}

	virtual void SetConstraintMode(int cm) {constraintmode = cm;}

	//$MaSch 2012-11-19
	virtual void SetEndTime(double time) {endtime=time;}
	virtual int IsActive() {return (endtime>0 && mbs->GetTime()>endtime) ? 0 : 1;}

protected:

	MathFunction mathfunc;
	int constraintmode; //$EDC[varaccess,EDCvarname="prescribe_velocity",EDCfolder="Kinematics",int_bool,tooltiptext="1 = prescribe function at velocity level, 0 = prescribe velocity at position/displacement level"] //0=position level, 1=velocity level

	double endtime; //$MaSch 2012-11-19: contraint is only active during 0<=t<=end_time; if end_time < 0 (default is -1), constraint is always active;
									//see also member functions SetEndTime and IsActive, and EvalG

	//int ptype;		//1=harmonic functions for velocity
	//Vector data;	//ptype=1: omega, phaseshift and amplitude of velocity
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Position constraint 2D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//Constraint for a position measured relative within one or two bodies
//leads to nonlinear constraint equations
class Pos2DConstraint: public Constraint
{
public:
	//Pos2DConstraint():loccoords(), dpdq() {mbs = NULL; };
	Pos2DConstraint(MBS* mbsi):Constraint(mbsi), loccoords(), dpdq(), prescribedtime(), prescribedvel()
	{
	};
	Pos2DConstraint(MBS* mbsi, int en1, const Vector2D& lc1, 
		const Vector2D& pglob, double cdim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq(), prescribedtime(), prescribedvel()
	{	
		DOFmode = 0;
		prescribed = 0;
		x_init = Vector(SS());
		GetCol() = coli;
		draw_dim.X() = cdim;
		p_global = pglob;
		AddElementCoord(en1, lc1);
		elementname = GetElementSpec();
	};
	Pos2DConstraint(MBS* mbsi, int en1, int en2, 
		const Vector2D& lc1, const Vector2D& lc2, double cdim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq(), prescribedtime(), prescribedvel() 
	{	
		DOFmode = 0;
		prescribed = 0;
		x_init = Vector(SS());
		GetCol() = coli;
		draw_dim.X() = cdim;
		AddElementCoord(en1, lc1);
		AddElementCoord(en2, lc2);
		elementname = GetElementSpec();
	};

	Pos2DConstraint(MBS* mbsi, int en1, int node1, 
		const Vector2D& pglob, double cdim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq(), prescribedtime(), prescribedvel() 
	{	
		Vector2D lc1((double)node1,0.);
		DOFmode = 1;
		prescribed = 0;
		x_init = Vector(SS());
		GetCol() = coli;
		draw_dim.X() = cdim;
		p_global = pglob;
		AddElementCoord(en1, lc1);
		elementname = GetElementSpec();
	};
	Pos2DConstraint(MBS* mbsi, int en1, int node1, int en2,
		const Vector2D& lc2, double cdim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq(), prescribedtime(), prescribedvel() 
	{	
		Vector2D lc1((double)node1,0.);
		DOFmode = 1; //only first node is node number
		prescribed = 0;
		x_init = Vector(SS());
		GetCol() = coli;
		draw_dim.X() = cdim;
		AddElementCoord(en1, lc1);
		AddElementCoord(en2, lc2);
		elementname = GetElementSpec();
	};
	Pos2DConstraint(MBS* mbsi, int en1, int en2, int node1, int node2, 
		double cdim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq(), prescribedtime(), prescribedvel() 
	{	
		Vector2D lc1((double)node1,0.);
		Vector2D lc2((double)node2,0.);
		DOFmode = 2; //both nodes are node numbers!!!
		prescribed = 0;
		x_init = Vector(SS());
		GetCol() = coli;
		draw_dim.X() = cdim;
		AddElementCoord(en1, lc1);
		AddElementCoord(en2, lc2);
		elementname = GetElementSpec();
	};



	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new Pos2DConstraint(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const Pos2DConstraint& ce = (const Pos2DConstraint&)e;
		loccoords = ce.loccoords;
		p_global = ce.p_global;
		dpdq = ce.dpdq;
		prescribed = ce.prescribed;
		p1 = ce.p1;
		p2 = ce.p2;
		p3 = ce.p3;
		p4 = ce.p4;
		DOFmode = ce.DOFmode;
		prescribedtime = ce.prescribedtime;
		prescribedvel = ce.prescribedvel;
	}

	virtual const char* GetElementSpec() const {return "Pos2DConstraint";}
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	virtual void AddElementCoord(int en, const Vector2D& loccoord)
	{
		AddElement(en);
		loccoords.Add(loccoord);
	}

	virtual void EvalG(Vector& f, double t);

	//To be replaced in derived class
	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size
	virtual int IS() const 
	{
		return 2;
	};

	virtual int Dim() const {return GetElem(1).Dim();}  //has 3D position???

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;

	virtual int IsPrescribed() const {return prescribed;}
	virtual void SetPrescribedCurve(double radius, double ang_vel, double t_accel, double init_angle)
	{
		prescribed = 1;
		p1 = radius;
		p2 = ang_vel;
		p3 = t_accel;
		p4 = init_angle;
	}
	virtual void SetPrescribedVel(TArray<double>& time, TArray<Vector2D> &vel)
	{
		prescribed = 2;
		prescribedtime = time;
		prescribedvel = vel;

		prescribedtime.Add(1e20); //continue last interval forever ...
		prescribedvel.Add(prescribedvel.Last());
	}

	virtual void DrawElement();

protected:
	TArray<Vector2D> loccoords;
	Vector2D p_global;

	Matrix dpdq; //temporary element, no need to copy???
	int prescribed;  // prescribed curves
	double p1;			 // parameter 1; 
	double p2;			 // parameter 2; 
	double p3;			 // parameter 3; 
	double p4;			 // parameter 4; 

	TArray<double> prescribedtime;
	TArray<Vector2D> prescribedvel;


	int DOFmode; //if DOFmode = 1, then loccoords(1).X() is the local node number in the element
};

/*class ArcLengthConstraint : public CoordConstraint
{
	public:
	ArcLengthConstraint(MBS* mbsi) : CoordConstraint(mbsi) {};

	ArcLengthConstraint(const ArcLengthConstraint& ct): CoordConstraint(ct.mbs)
	{
		CopyFrom(ct);
	};

	~ArcLengthConstraint() {};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new ArcLengthConstraint(*this);
		return ec;
	}

	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		CoordConstraint::CopyFrom(e);
		const ArcLengthConstraint& ce = (const ArcLengthConstraint&)e;
		load = ce.load;
		L = ce.L;
	}

	virtual const char* GetElementSpec() const {return "ArcLengthConstraint";}

	virtual void SetArcLengthConstraint(double load, double L, int elem, int loccoord);

	virtual int IS() const { return 1; } //implicit size -> number of constraint equations

	virtual void EvalG(Vector& f, double t);

	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

protected:
	double load;
	double L;
};*/

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Angle constraint 2D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//Constraint for a position measured relative within one or two bodies
//leads to nonlinear constraint equations
class Angle2DConstraint: public Constraint
{
public:
	//Angle2DConstraint():loccoords(), dpdq() {mbs = NULL; };
	Angle2DConstraint(MBS* mbsi):Constraint(mbsi), loccoords(), dpdq()
	{
	};
	Angle2DConstraint(MBS* mbsi, int en1, const Vector2D& lc1, 
		const double angle_glob, double cdim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq()
	{	
		x_init = Vector(SS());
		GetCol() = coli;
		draw_dim.X() = cdim;
		angle_global = angle_glob;
		AddElementCoord(en1, lc1);
		elementname = GetElementSpec();
	};
	Angle2DConstraint(MBS* mbsi, int en1, int en2, 
		const Vector2D& lc1, const Vector2D& lc2, double cdim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq()
	{	
		x_init = Vector(SS());
		GetCol() = coli;
		draw_dim.X() = cdim;
		AddElementCoord(en1, lc1);
		AddElementCoord(en2, lc2);
		elementname = GetElementSpec();
	};



	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new Angle2DConstraint(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const Angle2DConstraint& ce = (const Angle2DConstraint&)e;
		loccoords = ce.loccoords;
		angle_global = ce.angle_global;
		dpdq = ce.dpdq;
	}
	virtual const char* GetElementSpec() const {return "Angle2DConstraint";}

	virtual void AddElementCoord(int en, const Vector2D& loccoord)
	{
		AddElement(en);
		loccoords.Add(loccoord);
	}

	virtual void EvalG(Vector& f, double t);

	//To be replaced in derived class
	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size
	virtual int IS() const 
	{
		return 1;
	};

	virtual int Dim() const {return GetElem(1).Dim();}  //has 3D position???

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;

	virtual void DrawElement();

protected:
	TArray<Vector2D> loccoords;
	double angle_global;

	Matrix dpdq; //temporary element, no need to copy???
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Rotational Spring-Damper-Actor Element 2D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//Rotational Spring-Damper-Actuator connects two elements, with spring (kphi...stiffness, phi0...initial rotation),
// damper (d ...  damping coefficient) and a constant moment (ma ... constant actuator force)
class SDRotActor2D: public Constraint
{
public:
	//SDRotActor2D():loccoords(), dpdq() {mbs = NULL; };
	SDRotActor2D(MBS* mbsi):Constraint(mbsi), loccoords(), dpdq()
	{
	};
	SDRotActor2D(MBS* mbsi, int en1, const Vector2D& lc1, const Vector2D& gc2,
		double spring_stiffness, double spring_init_angle, double damping_coeff, double actuator_moment, double cdim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq()
	{	
		x_init = Vector(SS());
		GetCol() = coli;
		//draw_dim.X() = cdim;
		SetDrawSizeScalar(cdim);
		//k = spring_stiffness;
		SetPenaltyStiffness(spring_stiffness);
		phi0 = spring_init_angle;
		d = damping_coeff;
		ma = actuator_moment;
		AddElementCoord(en1, lc1);
		elementname = GetElementSpec();
	};
	SDRotActor2D(MBS* mbsi, int en1, int en2, 
		const Vector2D& lc1, const Vector2D& lc2,double spring_stiffness, double spring_init_angle, double damping_coeff, double actuator_moment, double cdim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq()
	{	
		x_init = Vector(SS());
		GetCol() = coli;
		//draw_dim.X() = cdim;
		SetDrawSizeScalar(cdim);
		//k = spring_stiffness;
		SetPenaltyStiffness(spring_stiffness);
		phi0 = spring_init_angle;
		d = damping_coeff;
		ma = actuator_moment;
		AddElementCoord(en1, lc1);
		AddElementCoord(en2, lc2);
		elementname = GetElementSpec();
	};

	virtual void EvalG(Vector& f, double t) {};
	virtual void EvalF(Vector& f, double t) {};
	virtual void EvalF2(Vector& f, double t) {};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new SDRotActor2D(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const SDRotActor2D& ce = (const SDRotActor2D&)e;
		loccoords = ce.loccoords;
		dpdq = ce.dpdq;
		//k = ce.k;
		d = ce.d;
		phi0 = ce.phi0;
		ma = ce.ma;
	}

	virtual void AddElementCoord(int en, const Vector2D& loccoord)
	{
		AddElement(en);
		loccoords.Add(loccoord);
	}

	//To be replaced in derived class
	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	virtual double ComputeForce(double t);
	
	//implicit (algebraic) size
	virtual int IS() const 
	{
		return 0;
	};

	virtual int Dim() const {return GetElem(1).Dim();}  //has 3D position???

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;

	virtual void DrawElement();

protected:
	TArray<Vector2D> loccoords;
	double phi0,d,ma;			//$ DR 2011-04-22: removed k
	Matrix dpdq; //temporary element, no need to copy???
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// UniversalJointOLD 3D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//Constraint for a position measured relative within one or two bodies
//leads to nonlinear constraint equations
//draw_dim.X()==diameter
//draw_dim.Y()==axis length
//draw_dim.Z()==draw resolution
class UniversalJointOLD: public Constraint
{
public:
	UniversalJointOLD(MBS* mbsi):Constraint(mbsi), loccoords(), dpdq(), hmat(), hvec()
	{
	};
	UniversalJointOLD(MBS* mbsi, int en1, const Vector3D& locp, const Vector3D& pglob, const Vector3D& globrot1,  const Vector3D& globrot2,
		const Vector3D& ddim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq(), hmat(), hvec() 
	{	
		x_init = Vector(SS());
		GetCol() = coli;
		draw_dim = ddim;
		p_global = pglob;

		Vector3D vrot1;
		vrot1 = globrot1;
		vrot1.Normalize();

		Vector3D vrot2;
		vrot2 = globrot2;
		vrot2.Normalize();

		AddElement(en1);
		loccoords.Add(locp);
		loccoords.Add(vrot1);
		loccoords.Add(vrot2); //global rotation axis
		loccoords.Add(Vector3D(0)); //lrot1body1
		loccoords.Add(Vector3D(0)); //lrot1body2, only for drawing

		elementname = GetElementSpec();
	};
	UniversalJointOLD(MBS* mbsi, int en1, int en2, 
		const Vector3D& lp1, const Vector3D& lp2, const Vector3D& globrot1,  const Vector3D& globrot2,
		const Vector3D& ddim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq(), hmat(), hvec() 
	{
		x_init = Vector(SS());
		GetCol() = coli;
		draw_dim = ddim;

		Vector3D vrot1;
		vrot1 = globrot1;
		vrot1.Normalize();

		Vector3D vrot2;
		vrot2 = globrot2;
		vrot2.Normalize();

		AddElement(en1);
		AddElement(en2);
		loccoords.Add(lp1);
		loccoords.Add(lp2);
		loccoords.Add(vrot1);
		loccoords.Add(vrot2);
		loccoords.Add(Vector3D(0)); //lrot1body1
		loccoords.Add(Vector3D(0));	//lrot2body2
		loccoords.Add(Vector3D(0)); //lrot1body2, only for drawing
		loccoords.Add(Vector3D(0)); //lrot2body1, only for drawing

		elementname = GetElementSpec();
	};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new UniversalJointOLD(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const UniversalJointOLD& ce = (const UniversalJointOLD&)e;
		loccoords = ce.loccoords;
		p_global = ce.p_global;
		dpdq = ce.dpdq;
		hmat = ce.hmat;
		hvec = ce.hvec;
	}

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	virtual const char* GetElementSpec() const {return "UniversalJointOLD";} //this name names MUST coincide with the type names defined in MBS::MBS() !!!!
	virtual void Initialize();

	virtual void EvalG(Vector& f, double t);

	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size
	virtual int IS() const {return 4;};

	virtual int Dim() const {return 3;}

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;

	virtual void DrawElement();

protected:
	TArray<Vector3D> loccoords;
	Vector3D p_global;

	Matrix dpdq; //temporary element, no need to copy???
	Matrix hmat;
	Vector hvec;
};



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// PrismaticJoint2D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//Constraint for a position measured relative within one or two bodies
//leads to nonlinear constraint equations
class PrismaticJoint2D: public Constraint
{
public:
	PrismaticJoint2D(MBS* mbsi):Constraint(mbsi), loccoords(), dpdq() 
	{
	};
	//relative motion of body and ground along axis lp1-gp2, lp1 and gp2 must never be equal!!
	PrismaticJoint2D(MBS* mbsi, int en1, 
		const Vector2D& lp1, const Vector2D& gp2,
		const Vector3D& ddim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq(), hmat(), hvec() 
	{	
		x_init = Vector(SS());
		GetCol() = coli;
		draw_dim = ddim;
		AddElement(en1);
		loccoords.Add(lp1);
		loccoords.Add(gp2);
		loccoords.Add(Vector2D(0)); //t1
		loccoords.Add(Vector2D(0)); //n2
		elementname = GetElementSpec();
	};
	//relative motion of two bodies along axis lp1-lp2, lp1 and lp2 must never be equal!!
	PrismaticJoint2D(MBS* mbsi, int en1, int en2, 
		const Vector2D& lp1, const Vector2D& lp2,
		const Vector3D& ddim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq(), hmat(), hvec() 
	{	
		x_init = Vector(SS());
		GetCol() = coli;
		draw_dim = ddim;
		AddElement(en1);
		AddElement(en2);
		loccoords.Add(lp1);
		loccoords.Add(lp2);
		loccoords.Add(Vector2D(0)); //t1
		loccoords.Add(Vector2D(0)); //n2
		elementname = GetElementSpec();
	};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new PrismaticJoint2D(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const PrismaticJoint2D& ce = (const PrismaticJoint2D&)e;
		loccoords = ce.loccoords;
		dpdq = ce.dpdq;
		hmat = ce.hmat;
		hvec = ce.hvec;
	}

	virtual void Initialize();

	virtual void EvalG(Vector& f, double t);

	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size
	virtual int IS() const {return 2;};

	virtual int Dim() const {return 2;}

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;

	virtual void DrawElement();

protected:
	TArray<Vector2D> loccoords;

	Matrix dpdq; //temporary element, no need to copy???
	Matrix hmat;
	Vector hvec;
};




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// RollingJoint3D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//A point at fixed relative position (in global coordinates) with respect to the center of the body must have zero velocity (or nearly zero)
//The attached body must be a rigid body!

// $ MSax 2013-03-07: [added to scriptlanguage
class RollingJoint3D: public Constraint //$EDC$[beginclass,classname=RollingJoint3D,parentclassname=Constraint,addelementtype=TAEconstraint+TAENotInRelease,addelementtypename=RollingJoint3D,texdescription="
//A point at fixed relative position (in global coordinates) with respect to the center of the body must have zero velocity (or nearly zero). The attached body must be a rigid body!",
//texdescriptionEquations="
//	contact point / normal directions: \\
//	$F_z = k_c\,z + d_c\,v$ \quad (a) \\
//	$k_c$ ... penalty\_stiffness\_contact (user defined parameter) \\
//	$d_c$ ... contact\_damping (user defined parameter) \\
//	$z$ ... contact\_penetration \\
//	$v$ ... contact\_velocity \\ \\
//	contact point - sticking / inplane directions: \\
//	$F_x = k_i\,x$ \quad (b.1) \\
//	$F_y = k_i\,y$ \quad (b.2) \\
//	$k_i$ ... penalty\_stiffness\_inplane (user defined parameter) \\
//	$x,y$ ... contact\_point\_veloctiy (relative velocity between ground and wheel contact point which is zero for exact sticking \\ \\
//	contact point - sliding / inplane directions: \\
//	$\mathbf{F_R} = -F_N\,f_c\,\mathbf{d}$ \\
//	$F_x = \mathbf{F_R}\,\mathbf{e_x}$ \quad (c.1) \\
//	$F_y = \mathbf{F_R}\,\mathbf{e_y}$ \quad (c.2) \\
//	$F_N$ ... contact force \\
//	$f_c$ ... friction\_coeff\_inplane (user defined parameter) \\
//	$\mathbf{d}$ ... direction vector of sliding (contact point) velocity \\
//	$\mathbf{e_x}, \left(\mathbf{e_y}\right)$ ... unit vector in forward (lateral) direction w.r.t gobal axes \\ \\
//	wheel center point / rolling friction: \\
//	$\mathbf{M_r} = f_{c,r}\,\mathbf{\omega}$ \quad (d) \\
//	$f_{c,r}$ ... rolling\_friction (user defined parameter) \\
//	$\mathbf{\omega}$ ... angular velocity vector \\
//	\begin{figure}[H]
//		\begin{center}
//			\includegraphics[width=0.95\textwidth]{D:/cpp/HotInt_V1/documentation/EDCauto_documentation/figures/RollingJoint3DStructure.jpg}
//			\caption{RollingJoint3D}
//		\end{center}
//	\end{figure}"]
{

public:

	RollingJoint3D(MBS* mbsi):Constraint(mbsi) // $ MSax 2013-03-07: added
	{	
		ElementDefaultConstructorInitialization();
	};

	RollingJoint3D(const RollingJoint3D& ct): Constraint(ct.mbs) // $ MSax 2013-03-07: added
	{
		ElementDefaultConstructorInitialization();
		CopyFrom(ct);
	};

	~RollingJoint3D() // $ MSax 2013-03-07: added
	{
	};

	void SetRollingJoint3D(int element_number, double wheel_radiusI, const Vector3D& wheel_local_center_pointI, const Vector3D& wheel_local_axisI, 
		const Vector3D& rolling_plane_pointI, const Vector3D& rolling_plane_normalI,
		double cdim, const Vector3D& coli)
	{	
		dpdq.SetSize(0,0);
		x_init = Vector(SS());

		Vector xdatainit(1.,1.); //sticking+is_contact
		SetDataInit(xdatainit);

		GetCol() = coli;
		draw_dim.X() = cdim;
		elements.Set1(element_number);
		//AddElement(element_number); //already has 1 element due to default constructor initialization!
		elementname = GetElementSpec();

		wheel_radius = wheel_radiusI;

		rolling_plane_point = rolling_plane_pointI;
		rolling_plane_normal = rolling_plane_normalI;

		wheel_local_center_point = wheel_local_center_pointI;
		wheel_local_axis = wheel_local_axisI;

		wheel_local_axis.Normalize();
		rolling_plane_normal.Normalize();
	};

	virtual void GetNecessaryKinematicAccessFunctions(TArray<int> &KinAccFunc, int numberOfKinematicPair) // MSax 2013-08-05: added
	{
		KinAccFunc.SetLen(0);
		//int kaf = (int)(TKinematicsAccessFunctions(TKAF_none));

		int kaf = (int)(TKinematicsAccessFunctions(TKAF_position+TKAF_ref_position+TKAF_rotation_matrix+TKAF_velocity+TKAF_D_pos_D_q));	
		KinAccFunc.Add(kaf);
	}

	void SetForceScalingFactor(double force_scaling_factor)
	{
		draw_dim.Y() = force_scaling_factor;
	}

	void SetContactFriction(int use_frictionI, double friction_coeff_inplaneI, int use_contact_conditionI)
	{
		use_friction = 1;
		friction_coeff_inplane = friction_coeff_inplaneI; //friction coeffiction in rolling direction
		use_contact_condition = use_contact_conditionI;

		//UO() << "ERROR: RollingJoint: friction not yet implemented!\n";
	}
	void SetRollingFrictionCoefficient(double rolling_frictionI) {rolling_friction = rolling_frictionI;}

	void SetPenaltyStiffness(int use_penalty_inplane_dirI, double penalty_stiffness_inplaneI,
		int use_penalty_contactI, double penalty_stiffness_contactI, double contact_dampingI)
	{
		use_penalty_inplane_dir = use_penalty_inplane_dirI;
		use_penalty_contact     = use_penalty_contactI;

		penalty_stiffness_inplane = penalty_stiffness_inplaneI; 
		penalty_stiffness_contact = penalty_stiffness_contactI; 
		contact_damping = contact_dampingI;

	}

	virtual const char* GetElementSpec() const {return "RollingJoint3D";} // $ MSax 2013-03-08: added

	virtual void ElementDefaultConstructorInitialization(); // $ MSax 2013-03-07: added

	virtual void Initialize(); // $ MSax 2013-03-07: added

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new RollingJoint3D(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const RollingJoint3D& ce = (const RollingJoint3D&)e;
		//loccoords = ce.loccoords;
		dpdq = ce.dpdq;

		rolling_plane_point = ce.rolling_plane_point;
		rolling_plane_normal = ce.rolling_plane_normal;

		wheel_local_center_point = ce.wheel_local_center_point;
		wheel_local_axis = ce.wheel_local_axis;

		wheel_radius = ce.wheel_radius; 

		use_penalty_inplane_dir = ce.use_penalty_inplane_dir;
		use_penalty_contact = ce.use_penalty_contact;

		penalty_stiffness_inplane = ce.penalty_stiffness_inplane; 
		penalty_stiffness_contact = ce.penalty_stiffness_contact; 

		use_friction = ce.use_friction;
		friction_coeff_inplane = ce.friction_coeff_inplane;

		contact_damping = ce.contact_damping;
		use_contact_condition = ce.use_contact_condition;

		rolling_friction = ce.rolling_friction;
	}

	//virtual void AddElementCoord(int en, const Vector3D& loccoord)
	//{
	//	AddElement(en);
	//	loccoords.Add(loccoord);
	//}

	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute

  virtual int GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //for special value sensor
	virtual int ReadSingleElementData(ReadWriteElementDataVariableType& RWdata); 		//for special value sensor

	virtual void EvalG(Vector& f, double t);

	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	virtual double PostNewtonStep(double t);

	virtual void PostprocessingStep();

	//implicit (algebraic) size
	virtual int IS() const {return 3 /*-use_penalty_contact-2*use_penalty_inplane_dir*/;};
	virtual int DataS() const {return 2;}	  //1=(is_sticking=1, is_sliding=0), 2=(1=is_contact,0=is_free)

	virtual const double& IsSticking() const {return XData(1);}
	virtual const double& IsContact() const {return XData(2);}
	virtual double& IsSticking() {return XData(1);}
	virtual double& IsContact() {return XData(2);}

	virtual int Dim() const {return 3;}

	virtual Vector3D ComputeForce() const;

	virtual double GetContactForce() const {return -XG(3);} //contact pressure is positive
	virtual double GetContactForceD() const {return -XGD(3);} //contact pressure is positive

	//compute distance to plane, contact point on plane
	virtual double WheelLocalContactPosition(Vector3D& contact_point, Vector3D& forward_dir, Vector3D& lateral_dir) const;
	virtual double WheelLocalContactPositionD(Vector3D& contact_point, Vector3D& forward_dir, Vector3D& lateral_dir) const;

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		
	virtual int SetElementDataAuto(ElementDataContainer& edc);
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;

	virtual void DrawElement();

protected:
	//TArray<Vector3D> loccoords; //not needed?//loccoords(1) is the wheel (local) center point, loccoords(2) is the wheel (local) axis of rotation
	
	Vector3D rolling_plane_point;  //$EDC$[varaccess,EDCvarname="rolling_plane_point",EDCfolder="Geometry",tooltiptext="any point at plane on that wheel is rolling"]//any point at plane on that wheel is rolling
	Vector3D rolling_plane_normal; //$EDC$[varaccess,EDCvarname="rolling_plane_normal",EDCfolder="Geometry",tooltiptext="normal vector of plane on that wheel is rolling"]//normal vector of plane on that wheel is rolling

	Vector3D wheel_local_center_point;  //$EDC$[varaccess,EDCvarname="wheel_local_center_point",EDCfolder="Geometry",tooltiptext="wheel center point in local coordinates of wheel"]//wheel center point in local coordinates of wheel
	Vector3D wheel_local_axis; //$EDC$[varaccess,EDCvarname="wheel_local_axis",EDCfolder="Geometry",tooltiptext="wheel axis in local coordinates of wheel"]//wheel axis in local coordinates of wheel

	double wheel_radius;  //$EDC$[varaccess,EDCvarname="wheel_radius",EDCfolder="Geometry",tooltiptext="radius of the idealized wheel"]//radius of the idealized wheel

	int use_penalty_inplane_dir;//$EDC$[varaccess,EDCvarname="use_penalty_inplane",int_bool,EDCfolder="Physics.Penalty",tooltiptext="1 = enabled, 0 = disabled"]
	int use_penalty_contact; //$EDC$[varaccess,EDCvarname="use_penalty_contact",int_bool,EDCfolder="Physics.Penalty",tooltiptext="1 = enabled, 0 = disabled"]

	double penalty_stiffness_inplane; //$EDC$[varaccess,EDCvarname="penalty_stiffness_inplane",EDCfolder="Physics.Penalty",minval=0,tooltiptext="penalty stiffness in rolling plane"]
	double penalty_stiffness_contact; //$EDC$[varaccess,EDCvarname="penalty_stiffness_contact",EDCfolder="Physics.Penalty",minval=0,tooltiptext="penalty stiffness in contact direction (contact wheel-ground)"] //penalty stiffness in contact direction (contact wheel-ground)
	double contact_damping; //$EDC$[varaccess,EDCvarname="contact_damping",EDCfolder="Physics.Penalty",minval=0,tooltiptext="linear contact damping coefficient"] //penalty stiffness in contact direction (contact wheel-ground)//damping in contact direction (contact wheel-ground)

	int use_contact_condition; //$EDC$[varaccess,EDCvarname="use_contact_condition",EDCfolder="Physics.Penalty",int_bool,tooltiptext="consider contact condition of wheel"] //consider contact condition of wheel
	int use_friction; //$EDC$[varaccess,EDCvarname="use_friction",EDCfolder="Physics.Penalty",int_bool,tooltiptext="friction force inplane against sliding velocity (in case if not sticking)"]
	double friction_coeff_inplane; //$EDC$[varaccess,EDCvarname="friction_coeff_inplane",EDCfolder="Physics.Penalty",tooltiptext="friction coefficient in rolling forward/lateral direction (must be combined)"] //friction coeffiction in rolling forward/lateral direction (must be combined)
	
	double rolling_friction;  //$EDC$[varaccess,EDCvarname="rolling_friction",EDCfolder="Physics.Penalty",tooltiptext="velocity proportional friction"]//velocity proportional friction

	//EDC double draw_dim(1); //$EDC$[varaccess,EDCvarname="draw_dim",EDCfolder="Graphics",tooltiptext="size of sphere (not used for computation)"]
	//EDC double draw_dim(2); //$EDC$[varaccess,EDCvarname="force_dim",EDCfolder="Graphics",tooltiptext="force scaling factor in [m/N] (not used for computation)"]

	//EDC int elements(1); //$EDC$[varaccess,EDCvarname="element_number",minval=1,EDCfolder="Element",tooltiptext="number of constrained element"]

	//EDC double data_init(1); //$EDC$[varaccess,EDCvarname="sticking",EDCfolder="Initialization",int_bool,tooltiptext="sticking condition at initial time step (1 if sticking, 0 if not sticking)"]
	//EDC double data_init(2); //$EDC$[varaccess,EDCvarname="contact",EDCfolder="Initialization",int_bool,tooltiptext="contact condition at initial time step (1 if in contact, 0 if not in contact)"]

	//remove unused data
	//==========================
	//EDC	int use_penalty_formulation; //$EDC$[varaccess,remove,EDCvarname="use_penalty_formulation",EDCfolder="Physics"]
	//EDC int use_local_coordinate_system;	//$EDC$[varaccess,remove,EDCvarname="use_local_coordinate_system",EDCfolder="Geometry"]
	//EDC	double spring_stiffness; //$EDC$[varaccess,remove,EDCvarname="spring_stiffness",EDCfolder="Physics.Penalty"]

	Matrix dpdq; //temporary element
};//$EDC$[endclass,RollingJoint3D]
// MSax 2013-03-07:]

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// RollingSegment3D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//A segment with defined angle +/- phi/2 can roll, about a relative point locpos1
//The segment rolls if a point of the segment has the rolling z-position rollz
//works only with a penalized approach in the case of rigid bodies

class RollingSegment3D: public Constraint
{
public:
	RollingSegment3D(MBS* mbsi):Constraint(mbsi), loccoords(), dpdq() {	};

	RollingSegment3D(MBS* mbsi, int en1, const Vector3D& relpos, double rollradI, double rollphiI, 
		double rollzI, double rollstiffI, double rolldampI, double cdim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq() 
	{	
		x_init = Vector(SS());
		GetCol() = coli;
		draw_dim.X() = cdim;

		rollz = rollzI;
		rollphi = rollphiI;
		rollrad = rollradI;
		rollstiff = rollstiffI;
		rolldamp = rolldampI;


		AddElementCoord(en1, relpos);
		elementname = GetElementSpec();
	};


	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new RollingSegment3D(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const RollingSegment3D& ce = (const RollingSegment3D&)e;
		loccoords = ce.loccoords;
		dpdq = ce.dpdq;

		rollphi = ce.rollphi;
		rollz = ce.rollz;
		rollrad = ce.rollrad;
		rollstiff = ce.rollstiff;
		rolldamp = ce.rolldamp;
	}

	virtual void AddElementCoord(int en, const Vector3D& loccoord)
	{
		AddElement(en);
		loccoords.Add(loccoord);
	}

	virtual void EvalG(Vector& f, double t) {};

	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size
	virtual int IS() const {return 0;};

	virtual int Dim() const {return 3;}

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;

	virtual void DrawElement();

protected:
	TArray<Vector3D> loccoords;

	Matrix dpdq; //temporary element, no need to copy???

	double rollphi;
	double rollz;
	double rollrad;
	double rollstiff;
	double rolldamp;
};




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// FlexibleRollingCable3D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//A segment with between two points p1 and p2 can roll on a surface z=rollz;
//The segment rolls if a point of the segment has the rolling z-position rollz
//works only with a penalized approach in the case of rigid bodies
class FlexibleRollingCable3D: public Constraint
{
public:
	FlexibleRollingCable3D(MBS* mbsi):Constraint(mbsi), loccoords(), dpdq() {	};

	FlexibleRollingCable3D(MBS* mbsi, int en1, const Vector3D& p1, const Vector3D& p2, 
		double rollzI, double rollstiffI, int nsegI, double cdim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq() 
	{	
		x_init = Vector(SS());
		GetCol() = coli;
		draw_dim.X() = cdim;

		rollz = rollzI;
		rollstiff = rollstiffI;
		nseg = nsegI;

		AddElementCoord(en1, p1);
		AddElementCoord(en1, p2);
		elementname = GetElementSpec();
	};


	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new FlexibleRollingCable3D(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const FlexibleRollingCable3D& ce = (const FlexibleRollingCable3D&)e;
		loccoords = ce.loccoords;
		dpdq = ce.dpdq;

		rollz = ce.rollz;
		rollstiff = ce.rollstiff;
		nseg = ce.nseg;
	}

	virtual void AddElementCoord(int en, const Vector3D& loccoord)
	{
		AddElement(en);
		loccoords.Add(loccoord);
	}

	virtual void EvalG(Vector& f, double t) {};

	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size
	virtual int IS() const {return 0;};

	virtual int Dim() const {return 3;}

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;

	virtual void DrawElement();

protected:
	TArray<Vector3D> loccoords;

	Matrix dpdq; //temporary element, no need to copy???

	double rollstiff;
	double rollz;
	int nseg;
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// RollingJoint2D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//Constraint for a position measured relative within one or two bodies
//leads to nonlinear constraint equations
class RollingJoint2D: public Constraint
{
public:
	RollingJoint2D(MBS* mbsi):Constraint(mbsi), loccoords(), dpdq() 
	{
	};
	//relative motion of body and ground along axis lp1-gp2, lp1 and gp2 must never be equal!!
	RollingJoint2D(MBS* mbsi, int en1, 
		const Vector2D& lp1, const Vector2D& lt1,
		const Vector3D& ddim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq() 
	{	
		x_init = Vector(SS());
		GetCol() = coli;
		draw_dim = ddim;
		AddElement(en1);
		loccoords.Add(lp1);
		loccoords.Add(lt1);
		elementname = GetElementSpec();
	};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new RollingJoint2D(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const RollingJoint2D& ce = (const RollingJoint2D&)e;
		loccoords = ce.loccoords;
		dpdq = ce.dpdq;
	}

	virtual void Initialize();

	virtual void EvalG(Vector& f, double t);

	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size
	virtual int IS() const {return 1;};

	virtual int Dim() const {return 2;}

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;

	virtual void DrawElement();

protected:
	TArray<Vector2D> loccoords;

	Matrix dpdq; //temporary element, no need to copy???
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: SlidingJoint (Constraint)
// short description: 
// available formulations: Lagrange multiplier
// type of constraint equation: nonlinear constraint equations
// index formulation available: 
// ground constraint: no
// element-element constraint: yes
// development status: 
// long description:	a point of body 2 may slide along x-axis of body 1
//										
// class variables:
//  
//**end(ued)**
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class SlidingJoint: public Constraint
{
public:
	SlidingJoint(MBS* mbsi):Constraint(mbsi), loccoords(), dpdq() 
	{
		ElementDefaultConstructorInitialization();
	};
	SlidingJoint(const SlidingJoint& e):Constraint(e.mbs) {CopyFrom(e);};

	//lc2 is a proper initial condition!!!
	SlidingJoint(MBS* mbsi, int en1, int en2, 
		const Vector3D& lc1, const Vector3D& lc2, double cdim, const Vector3D& coli);

	SlidingJoint(MBS* mbsi, int en1, const TArray<int>& en2, int in2, 
		const Vector3D& lc1, const Vector3D& lc2, double cdim, const Vector3D& coli);

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new SlidingJoint(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e);

	virtual void ElementDefaultConstructorInitialization();
	virtual const char* GetElementSpec() const {return "SlidingJoint";}

	virtual void AddElementCoord(int en, const Vector3D& loccoord)
	{
		AddElement(en);
		loccoords.Add(loccoord);
	}

	virtual void EvalG(Vector& f, double t);
	virtual void EvalF(Vector& f, double t);

	virtual double PostNewtonStep(double t);
	virtual void PostprocessingStep();

	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size
	virtual int IS() const {return 4;};
	virtual int ES() const {return 1;};

	virtual int Dim() const {return 3;}

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;
	virtual Vector3D GetRefPos() const;
	virtual Vector3D GetRefVelD() const;
	virtual Vector3D GetRefVel() const;

	virtual Vector3D GetSlidingPos() const;

	virtual void DrawElement();
	virtual double GetDrift() const;
	virtual double GetDriftP();
	virtual void SetNoForceY(int val) {noforcey = val;}

	//compute dpdx, global position of point at sliding axis and global velocity at sliding axis: 
	virtual void ComputeDpDx(Vector3D& dpdx, Vector3D& gp2, Vector3D& gv2);

	virtual void SetDriveParameters(int p_driven, double p_drivestarttime, double p_drivevel, double p_drivelen, double p_driveacceltime = 0.1)
	{
		driven = p_driven;
		drivestarttime = p_drivestarttime;
		drivevel = p_drivevel;
		drivelen = p_drivelen;
		driveacceltime = p_driveacceltime;
		driveflag = 0;
	}

	virtual void SetDriveEndTime(double te) {driveendtime = te;}
	virtual void SetDriveNextStart(double te) {drivenextstart = te;}

protected:
	TArray<Vector3D> loccoords;
	int elemind;
	double xoff;
	int nlstepcnt;
	int ind_init;

	//especially for 3D pantograph:
	int noforcey;

	Matrix dpdq; //temporary element, no need to copy???

	static const int ns = 1;
	static const int nsp = 2;
	static const int l1 = 3;
	static const int l2 = 4;
	static const int l3 = 5;

	int iscontact, iscontact2; //state of contact, old contact state

	//for driven mass along body:
	int driven;            //1=yes, driven
	double drivestarttime; //time when driving is started
	double driveendtime;   //time when driving is finished
	double drivevel;       //driving velocity
	double drivelen;			 //curve length, after which sliding is stopped (body falls down)
	double driveacceltime;			 //curve length, after which sliding is stopped (body falls down)
	int driveflag;         //if set to 1, driving is changed into completly free motion
	double drivenextstart; //time after which body is started again
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//# author:               Saxinger Martin
//#
//# generated:						19. Juli 2012
//# description:          This file contains the SlidingPointJoint and the PrismaticJoint
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: SlidingPointJoint (Constraint)
// short description: 
// available formulations: Lagrange multiplier
// type of constraint equation: nonlinear constraint equations
// index formulation available: 
// ground constraint: no
// element-element constraint: yes
// development status: 
// long description:	a point of body 2 may slide along x-axis of body 1
//										
// class variables:
//  
//**end(ued)**
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//MSax: For more information see document SlidingJoints.pdf
class SlidingPointJoint: public Constraint //$EDC$[beginclass,classname=SlidingPointJoint,parentclassname=Constraint,addelementtype=TAEconstraint,addelementtypename=SlidingPointJoint,figure="SlidingPointJoint",texdescription="
//This joint enables sliding of a fixed point of a body i along the x - axis of another body j. Both body i and body j can be flexible or rigid. Body j can contain more than one elements. No rotations are constrained at all. Only a Lagrangian formulation is implemented, the penalty formulation is not implemented yet. A MaxIndex 2 and 3 formulation exists.",
//modus="{sliding along a single body}{The vector Geomety.element\_numbers is equal to $[en1,en2]$. Index Geomety.elemind must be $1$.}",
//modus="{sliding along more than one body}{Geomety.element\_numbers has to be set to $[en1,en2_1,en2_2,...,en2_n]$. Geomety.elemind is the body j index of the element in inital configuration, e.g. for $en2_2$ the elemind is $2$.}",
//texdescriptionDOF="
//The vector of the DOF contains the sliding parameter $s$, its time derivative $\dot{s}$ and the vector of the Lagrangian parameters $\mathbf{\lambda}=\left[\lambda_1 \lambda_2 \lambda_3\right]^\mathbf{T}$. The Lagrange parameters $\lambda_1$ to $\lambda_3$ are representing the sliding forces in the global coordinate system.
//\begin{eqnarray}
//	\mathbf{q} &=& \left[
//	\begin{array}{ccccc}
//		s & \dot{s} & \lambda_1 & \lambda_2 & \lambda_3 
//	\end{array}
//	\right]^\mathbf{T} \quad
//\end{eqnarray}",
//texdescriptionEquations="
//positions:
//\begin{eqnarray}
//	\mathbf{x}^i &=& \left[
//	\begin{array}{ccc}
//		x_{1}^i & x_{2}^i & x_{3}^i 
//	\end{array}
//	\right]^\mathbf{T} 
//\end{eqnarray} \quad
//\begin{eqnarray}
//	\mathbf{x}^j &=& \left[
//	\begin{array}{ccc}
//		x_{1}^j = s & x_{2}^j & x_{3}^j 
//	\end{array}
//	\right]^\mathbf{T}
//\end{eqnarray}
//constraint equation - position level
//\begin{eqnarray}
//	\mathbf{C} &=& \left[
//	\begin{array}{c}
//		\mathbf{r}^i\left(\mathbf{x}^i\right)-\mathbf{r}^j\left(\mathbf{x}^j\right)  \medskip \\
//		\frac{\partial \mathbf{r}^j\left(\mathbf{x}^j\right)}{\partial x_1^j} \mathbf{\lambda}
//	\end{array}
//	\right] = \mathbf{0} \quad
//\end{eqnarray}
//The first three constraints restrict the motion of the sliding point on body $i$ and $j$. The fourth constraint equation ensures, that there is no force in the sliding direction. \\ \\
//constraint equation - velocity level:
//\begin{eqnarray}
//	\mathbf{C} &=& \left[
//	\begin{array}{c}
//		\frac{\partial \mathbf{r}^i\left(\mathbf{x}^i\right)}{\partial t}-\frac{\mathbf{r}^j\left(\mathbf{x}^j\right)}{\partial t}-\frac{\mathbf{r}^j\left(\mathbf{x}^j\right)}{\partial x_1^j} \dot{s}  \medskip \\
//		\frac{\partial \mathbf{r}^j\left(\mathbf{x}^j\right)}{\partial x_1^j} \mathbf{\lambda}
//	\end{array}
//	\right] = \mathbf{0} \quad
//\end{eqnarray}
//To obtain the constraints for velocity level, the first three equations are differentiated with respect to time. The sliding parameter $s$ is also a function of time. The fourth constraint equation is equal to the position level equation.",
//example="SlidingPointJoint.txt"]
{
public:
	SlidingPointJoint(MBS* mbsi):Constraint(mbsi)
	{
		ElementDefaultConstructorInitialization();
	};
	SlidingPointJoint(const SlidingPointJoint& e):Constraint(e.mbs) {CopyFrom(e);};

	//lc2 is a proper initial condition!!!
	void SetSlidingPointJoint(int en1, int en2, 
		const Vector3D& lc1, const Vector3D& lc2, double cdim, const Vector3D& coli);

	void SetSlidingPointJoint(int en1, const TArray<int>& en2, int in2, 
		const Vector3D& lc1, const Vector3D& lc2, double cdim, const Vector3D& coli);

	virtual void ElementDefaultConstructorInitialization();

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new SlidingPointJoint(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e);

	virtual void Initialize();

	virtual const char* GetElementSpec() const {return "SlidingPointJoint";}

	// for constraints only, these are the necessary (!) access functions!	//$ DR 2013-02-18
	virtual void GetNecessaryKinematicAccessFunctions(TArray<int> &KAF, int numberOfKinematicPair);

	virtual void EvalG(Vector& f, double t);
	virtual void EvalF(Vector& f, double t);

	virtual double PostNewtonStep(double t);

	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size
	virtual int IS() const {return 4;};
	virtual int ES() const {return 1;};

	virtual void SetGlobalInitConditions(Vector& x_glob);

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer
  virtual int GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //for special value sensor
	virtual int ReadSingleElementData(ReadWriteElementDataVariableType& RWdata); 		//for special value sensor

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		
	virtual int SetElementDataAuto(ElementDataContainer& edc);
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual int Dim() const {return 3;}

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;
	virtual Vector3D GetRefPos() const;
	virtual Vector3D GetRefVelD() const;
	virtual Vector3D GetRefVel() const;

	virtual Vector3D GetSlidingPos() const;

	virtual void DrawElement();
	virtual double GetDrift() const;
	virtual double GetDriftP();

	//compute dpdx, global position of point at sliding axis and global velocity at sliding axis: 
	virtual void ComputeDpDx(Vector3D& dpdx, Vector3D& gp2, Vector3D& gv2); // $MSax: gp...global position, gv...global velocity

protected:
	TArray<Vector3D> loccoords;
	int elemind; //$EDC$[varaccess,minval=1,EDCvarname="elemind",EDCfolder="Geometry",tooltiptext="Index of the initial sliding body."] 

	//EDC Vector3D loccoords(1);  //$EDC$[varaccess,EDCvarname="position_1",EDCfolder="Geometry",tooltiptext="Vector from the center of body number 1 (en1) to the sliding point in the local body 1 coordinate system."]
	//EDC Vector3D loccoords(2);	//$EDC$[varaccess,EDCvarname="position_2",EDCfolder="Geometry",tooltiptext="Vector from the center of the first body of en2 array to the sliding point in the local body 2 coordinate system."]
	//EDC TArray<int> elements;  //$EDC$[varaccess,EDCvarname="element_numbers",variable_length_vector,minval=1,EDCfolder="Geometry",tooltiptext="Element numbers: [en1,en2_1,en2_2,...,en2_n]."]

	double xoff; // position offset in x direction (el 2) between cog el 2_1 to cog el 2_elemind

	Matrix dpdq;

	static const int ns = 1; // index of sliding parameter s in XG vector
	static const int nsp = 2; // index of time derivative of sliding parameter s in XG vector
	static const int l1 = 3;  // index of lagrange paramter 1 in XG vector
	static const int l2 = 4;
	static const int l3 = 5;

	//EDC int use_local_coordinate_system;	//$EDC$[varaccess,remove,EDCvarname="use_local_coordinate_system",EDCfolder="Geometry"]
	//EDC int use_penalty_formulation;	//$EDC$[varaccess,remove,EDCvarname="use_penalty_formulation",EDCfolder="Physics"]
	//EDC double spring_stiffness; //$EDC$[varaccess,remove,EDCvarname="spring_stiffness",EDCfolder="Physics.Penalty"]

	//EDC double draw_dim(1);	 //$EDC$[varaccess,EDCvarname="draw_size",EDCfolder="Graphics",tooltiptext="Drawing dimensions of constraint. If set to -1, than global_draw_scalar_size is used."]

};  //$EDC$[endclass,SlidingPointJoint]

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: SlidingPrismaticJoint (Constraint)
// short description: 
// available formulations: Lagrange multiplier
// type of constraint equation: nonlinear constraint equations
// index formulation available: 
// ground constraint: no
// element-element constraint: yes
// development status: 
// long description:	a point of body 2 may slide along x-axis of body 1
//										
// class variables:
//  
//**end(ued)**
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//MSax: For more information see document SlidingJoints.pdf
class SlidingPrismaticJoint: public Constraint //$EDC$[beginclass,classname=SlidingPrismaticJoint,parentclassname=Constraint,addelementtype=TAEconstraint,addelementtypename=SlidingPrismaticJoint,texdescription="
//This joint enables sliding of a fixed point of a body i along the x - axis of another body j. Both body i and body j can be flexible or rigid. Body j can contain more than one element. The difference to the SlidingPointJoint is that the relative rotation between the bodies is also constrained. A Lagrangian formulation is used for both stiff and springy constrained rotation. For the position constraint only a stiff formulation exists. A penalty formulation is not implemented yet. There is a MaxIndex 2 and 3 formulation implemented.",
//modus="{sliding along a single body}{The vector Geomety.element\_numbers is equal to $[en1,en2]$. Index Geomety.elemind must be $1$.}",
//modus="{sliding along more than one body}{Geomety.element\_numbers has to be set to $[en1,en2_1,en2_2,...,en2_n]$. Geomety.elemind is the body j index of the element in inital configuration, e.g. for $en2_2$ the elemind is $2$.}",
//modus="{stiff constrained rotation}{Physics.use\_penalty\_formulation is set to $0$.}",
//modus="{springy constrained rotation}{Physics.use\_penalty\_formulation is set to $1$. The values for stiffness and damping must be set in Physics.Penalty folder.}",
//texdescriptionDOF="
//The vector of the DOF contains the sliding parameter $s$, its time derivative $\dot{s}$ and the vector of the Lagrangian parameters $\mathbf{\lambda}=\left[\lambda_1 \lambda_2 \lambda_3\right]^\mathbf{T}$. The Lagrange parameters $\lambda_1$ to $\lambda_3$ are representing the sliding forces in the global coordinate system. The three Lagrangian parameters $\lambda_4$ to $\lambda_6$ are the sliding moments about the global coordinate system axes.
//\begin{eqnarray}
//	\mathbf{q} &=& \left[
//	\begin{array}{cccccccc}
//		s & \dot{s} & \lambda_1 & \lambda_2 & \lambda_3 & \lambda_4 & \lambda_5 & \lambda_6
//	\end{array}
//	\right]^\mathbf{T} \quad
//\end{eqnarray}",
//texdescriptionEquations="
//At initialization the unit vectors of the global coordinate system are transformed to the local coordinate system of each body and the vectors $\mathbf{v}_1^i$, $\mathbf{v}_2^i$ and $\mathbf{v}_3^i$ for body $i$ and $\mathbf{v}_1^j$, $\mathbf{v}_2^j$ and $\mathbf{v}_3^j$ for body $j$ are obtained. The vectors are fixed in the body coordinate system. The position vectors are the same as for the SlidingPointJoint. \\ \\
//constraint equation - position level (stiff connection)
//\begin{eqnarray}
//	\mathbf{C} &=& \left[
//	\begin{array}{c}
//		\mathbf{r}^i\left(\mathbf{x}^i\right)-\mathbf{r}^j\left(\mathbf{x}^j\right)  \medskip \\
//		\frac{\partial \mathbf{r}^j\left(\mathbf{x}^j\right)}{\partial x_1^j} \mathbf{\lambda}  \medskip \\
//		\mathbf{v}_2^j\mathbf{v}_3^i  \medskip \\
//		\mathbf{v}_3^j\mathbf{v}_1^i  \medskip \\
//		\mathbf{v}_2^j\mathbf{v}_1^i
//	\end{array}
//	\right] = \mathbf{0} \quad
//\end{eqnarray}
//constraint equation - position level (springy connection)
//\begin{eqnarray}
//	\mathbf{C} &=& \left[
//	\begin{array}{c}
//		\mathbf{r}^i\left(\mathbf{x}^i\right)-\mathbf{r}^j\left(\mathbf{x}^j\right)  \medskip \\
//		\frac{\partial \mathbf{r}^j\left(\mathbf{x}^j\right)}{\partial x_1^j} \mathbf{\lambda}  \medskip \\
//		\mathbf{v}_2^j\mathbf{v}_3^ik_1+\left(\mathbf{\dot{v}}_2^j\mathbf{v}_3^i+\mathbf{v}_2^j\mathbf{\dot{v}}_3^i\right)d_1+\lambda_4  \medskip \\
//		\mathbf{v}_3^j\mathbf{v}_1^ik_1+\left(\mathbf{\dot{v}}_3^j\mathbf{v}_1^i+\mathbf{v}_3^j\mathbf{\dot{v}}_1^i\right)d_2+\lambda_5  \medskip \\
//		-\mathbf{v}_2^j\mathbf{v}_1^ik_1-\left(\mathbf{\dot{v}}_2^j\mathbf{v}_1^i+\mathbf{v}_2^j\mathbf{\dot{v}}_1^i\right)d_3+\lambda_6  \medskip \\
//	\end{array}
//	\right] = \mathbf{0} \quad,
//\end{eqnarray}
//constraint equation - velocity level (stiff connection)
//\begin{eqnarray}
//	\mathbf{C} &=& \left[
//	\begin{array}{c}
//		\frac{\partial \mathbf{r}^i\left(\mathbf{x}^i\right)}{\partial t}-\frac{\mathbf{r}^j\left(\mathbf{x}^j\right)}{\partial t}-\frac{\mathbf{r}^j\left(\mathbf{x}^j\right)}{\partial x_1^j} \dot{s}  \medskip \\
//		\frac{\partial \mathbf{r}^j\left(\mathbf{x}^j\right)}{\partial x_1^j} \mathbf{\lambda}  \medskip \\
//		\mathbf{\dot{v}}_2^j\mathbf{v}_3^i+\mathbf{v}_2^j\mathbf{\dot{v}}_3^i  \medskip \\
//		\mathbf{\dot{v}}_3^j\mathbf{v}_1^i+\mathbf{v}_3^j\mathbf{\dot{v}}_1^i  \medskip \\
//		\mathbf{\dot{v}}_2^j\mathbf{v}_1^i+\mathbf{v}_2^j\mathbf{\dot{v}}_1^i
//	\end{array}
//	\right] = \mathbf{0} \quad
//\end{eqnarray}
//constraint equation - velocity level (springy connection)
//\begin{eqnarray}
//	\mathbf{C} &=& \left[
//	\begin{array}{c}
//		\frac{\partial \mathbf{r}^i\left(\mathbf{x}^i\right)}{\partial t}-\frac{\mathbf{r}^j\left(\mathbf{x}^j\right)}{\partial t}-\frac{\mathbf{r}^j\left(\mathbf{x}^j\right)}{\partial x_1^j} \dot{s}  \medskip \\
//		\frac{\partial \mathbf{r}^j\left(\mathbf{x}^j\right)}{\partial x_1^j} \mathbf{\lambda}  \medskip \\		\mathbf{v}_2^j\mathbf{v}_3^ik_1+\left(\mathbf{\dot{v}}_2^j\mathbf{v}_3^i+\mathbf{v}_2^j\mathbf{\dot{v}}_3^i\right)d_1+\lambda_4  \medskip \\
//		\mathbf{v}_3^j\mathbf{v}_1^ik_1+\left(\mathbf{\dot{v}}_3^j\mathbf{v}_1^i+\mathbf{v}_3^j\mathbf{\dot{v}}_1^i\right)d_2+\lambda_5  \medskip \\
//		-\mathbf{v}_2^j\mathbf{v}_1^ik_1-\left(\mathbf{\dot{v}}_2^j\mathbf{v}_1^i+\mathbf{v}_2^j\mathbf{\dot{v}}_1^i\right)d_3+\lambda_6  \medskip \\
//	\end{array}
//	\right] = \mathbf{0} \quad
//\end{eqnarray}",example="SlidingPrismaticJoint.txt"]

{
public:
	SlidingPrismaticJoint(MBS* mbsi):Constraint(mbsi)
	{ 
		ElementDefaultConstructorInitialization();
	};
	SlidingPrismaticJoint(const SlidingPrismaticJoint& e):Constraint(e.mbs) {CopyFrom(e);};

	//lc2 is a proper initial condition!!!
	void SetSlidingPrismaticJoint(int en1, int en2, 
		const Vector3D& lc1, const Vector3D& lc2, double cdim, const Vector3D& coli, const Vector3D& stiffness=Vector3D(0), const Vector3D& damping=Vector3D(0));
	
	void SetSlidingPrismaticJoint(int en1, const TArray<int>& en2, int in2,
		const Vector3D& lc1, const Vector3D& lc2, double cdim, const Vector3D& coli, const Vector3D& stiffness=Vector3D(0), const Vector3D& damping=Vector3D(0));

	virtual void LinkToElements(); 

	virtual void ElementDefaultConstructorInitialization();

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new SlidingPrismaticJoint(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e);

	virtual const char* GetElementSpec() const {return "SlidingPrismaticJoint";}

	virtual void Initialize();

	// for constraints only, these are the necessary (!) access functions!	//$ DR 2013-02-18
	virtual void GetNecessaryKinematicAccessFunctions(TArray<int> &KAF, int numberOfKinematicPair);

	virtual void EvalG(Vector& f, double t);
	virtual void EvalF(Vector& f, double t);

	virtual double PostNewtonStep(double t);

	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size
	virtual int IS() const {return 7;};
	virtual int ES() const {return 1;};
	virtual int SOS() const {return 0;}; //size of stiffness/mass matrix
	virtual int SOSowned() const {return SOS();}; //number of second order unknowns added by constraint

	virtual void SetGlobalInitConditions(Vector& x_glob)
	{
		// set s (x_init(1)) as the x component of lc2 and the other entries are already zero
		x_init(1) = loccoords(2).X();
		for (int i=1; i<=SS(); i++)
		{
			x_glob(ltg.Get(i)) = x_init(i);
		}
	}

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer
  virtual int GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //for special value sensor
	virtual int ReadSingleElementData(ReadWriteElementDataVariableType& RWdata); 		//for special value sensor

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		
	virtual int SetElementDataAuto(ElementDataContainer& edc);
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual int Dim() const {return 3;}

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;
	virtual Vector3D GetRefPos() const;
	virtual Vector3D GetRefVelD() const;
	virtual Vector3D GetRefVel() const;

	virtual Vector3D GetSlidingPos() const;

	virtual void DrawElement();
	virtual double GetDrift() const;
	virtual double GetDriftP();

	//compute dpdx, global position of point at sliding axis and global velocity at sliding axis: 
	virtual void ComputeDpDx(Vector3D& dpdx, Vector3D& gp2, Vector3D& gv2);

protected:
	TArray<Vector3D> loccoords;
	//EDC Vector3D loccoords(1);  //$EDC$[varaccess,EDCvarname="position_1",EDCfolder="Geometry",tooltiptext="Vector from the center of body number 1 (en1) to the sliding point in the local body 1 coordinate system."]
	//EDC Vector3D loccoords(2);	//$EDC$[varaccess,EDCvarname="position_2",EDCfolder="Geometry",tooltiptext="Vector from the center of the first body of en2 array to the sliding point in the local body 2 coordinate system."]
	//EDC TArray<int> elements;  //$EDC$[varaccess,EDCvarname="element_numbers",variable_length_vector,minval=1,EDCfolder="Geometry",tooltiptext="Element numbers: [en1,en2_1,en2_2,...,en2_n]."]

	int elemind; //$EDC$[varaccess,minval=1,EDCvarname="elemind",EDCfolder="Geometry",tooltiptext="Index of the initial sliding body."] 
	double xoff;
	int ind;

	double k1; //$EDC$[varaccess,minval=0,EDCvarname="k1",EDCfolder="Physics.Penalty",tooltiptext="Stiffness for rotation about global x - axis."] 
	double k2; //$EDC$[varaccess,minval=0,EDCvarname="k2",EDCfolder="Physics.Penalty",tooltiptext="Stiffness for rotation about global y - axis."] 
	double k3; //$EDC$[varaccess,minval=0,EDCvarname="k3",EDCfolder="Physics.Penalty",tooltiptext="Stiffness for rotation about global z - axis."] 
	double d1; //$EDC$[varaccess,minval=0,EDCvarname="d1",EDCfolder="Physics.Penalty",tooltiptext="Damping of rotation about global x - axis."] 
	double d2; //$EDC$[varaccess,minval=0,EDCvarname="d2",EDCfolder="Physics.Penalty",tooltiptext="Damping of rotation about global x - axis."] 
	double d3; //$EDC$[varaccess,minval=0,EDCvarname="d3",EDCfolder="Physics.Penalty",tooltiptext="Damping of rotation about global x - axis."] 

	Matrix dpdq;
	Matrix drotdq;

	static const int ns = 1;
	static const int nsp = 2;
	static const int l1 = 3; // Fx
	static const int l2 = 4; // Fy
	static const int l3 = 5; // Fz
	static const int l4 = 6; // Mx
	static const int l5 = 7; // My
	static const int l6 = 8; // Mz

	//EDC int use_local_coordinate_system;	//$EDC$[varaccess,remove,EDCvarname="use_local_coordinate_system",EDCfolder="Geometry"]
	//EDC double spring_stiffness; //$EDC$[varaccess,remove,EDCvarname="spring_stiffness",EDCfolder="Physics.Penalty"]

	//EDC double draw_dim(1);	 //$EDC$[varaccess,EDCvarname="draw_size",EDCfolder="Graphics",tooltiptext="Drawing dimensions of constraint. If set to -1, than global_draw_scalar_size is used."]

};  //$EDC$[endclass,SlidingPrismaticJoint]


//Constraint: a point of body 2 may slide along x-axis of body 1
//leads to nonlinear constraint equations
class SlidingJoint2D: public Constraint
{
public:
	SlidingJoint2D(MBS* mbsi):Constraint(mbsi), loccoords(), dpdq() 
	{
		constrainrotation = 0; //for SS()
	};
	//lc2 is a proper initial condition!!!
	SlidingJoint2D(MBS* mbsi, int en1, int en2, 
		const Vector2D& lc1, const Vector2D& lc2, double cdim, const Vector3D& coli, int constrainrotationI = 0):Constraint(mbsi), loccoords(), dpdq() 
	{	
		constrainrotation = constrainrotationI; //must be first, for x_init!!!!

		x_init = Vector(SS()); // s, lambda1, lambda2, lambda3
		x_init(1) = lc2.X();

		GetCol() = coli;
		draw_dim.X() = cdim;
		AddElementCoord(en1, lc1);
		AddElementCoord(en2, lc2);
		elemind = 1; //sliding element, 1 means elements(2)
		xoff = 0;
		nlstepcnt = 0;
		ind_init = elemind;
		driven = 0;
		driveendtime = 1e10;
		drivenextstart = 0;
		driveflag = 0;
		elementname = GetElementSpec();
	};
	SlidingJoint2D(MBS* mbsi, int en1, const TArray<int>& en2, int in2, 
		const Vector2D& lc1, const Vector2D& lc2, double cdim, const Vector3D& coli, int constrainrotationI = 0);

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new SlidingJoint2D(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const SlidingJoint2D& ce = (const SlidingJoint2D&)e;
		loccoords = ce.loccoords;
		dpdq = ce.dpdq;
		elemind = ce.elemind;
		xoff = ce.xoff;
		nlstepcnt = ce.nlstepcnt;
		ind_init = ce.ind_init;

		constrainrotation = ce.constrainrotation;

		driven = ce.driven;
		drivestarttime = ce.drivestarttime;
		driveendtime = ce.driveendtime;
		drivevel = ce.drivevel;
		drivelen = ce.drivelen;
		driveflag = ce.driveflag;
		drivenextstart = ce.drivenextstart;
	}

	virtual void AddElementCoord(int en, const Vector2D& loccoord)
	{
		AddElement(en);
		loccoords.Add(loccoord);
	}

	virtual void EvalG(Vector& f, double t);
	virtual void EvalF(Vector& f, double t);

	virtual double PostNewtonStep(double t);
	virtual void PostprocessingStep();

	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size
	virtual int IS() const {return 3+constrainrotation;};
	virtual int ES() const {return 1;};

	virtual int Dim() const {return 2;}

	//get a reference position of the body in 3d
	virtual Vector2D GetRefPos2D() const;
	virtual Vector2D GetRefPos2DD() const;

	virtual Vector2D GetSlidingPos() const;

	virtual void DrawElement();
	virtual double GetDrift() const;

	//compute dpdx, global position of point at sliding axis and global velocity at sliding axis: 
	virtual void ComputeDpDx(Vector2D& dpdx, Vector2D& gp2, Vector2D& gv2);

	virtual void SetDriveParameters(int p_driven, double p_drivestarttime, double p_drivevel, double p_drivelen, double p_driveacceltime = 0.1)
	{
		driven = p_driven;
		drivestarttime = p_drivestarttime;
		drivevel = p_drivevel;
		drivelen = p_drivelen;
		driveacceltime = p_driveacceltime;
		driveflag = 0;
	}

	virtual void SetDriveEndTime(double te) {driveendtime = te;}
	virtual void SetDriveNextStart(double te) {drivenextstart = te;}

protected:
	TArray<Vector2D> loccoords;
	int elemind;
	double xoff;
	int nlstepcnt;
	int ind_init;

	Matrix dpdq; //temporary element, no need to copy???

	static const int ns = 1;
	static const int nsp = 2;
	static const int l1 = 3;
	static const int l2 = 4;
	static const int lphi = 5;

	//for driven mass along body:
	int constrainrotation; //additionally constrain the rotation of the attached body (sliding cylindrical joint)
	int driven;            //1=yes, driven
	double drivestarttime; //time when driving is started
	double driveendtime;   //time when driving is finished
	double drivevel;       //driving velocity
	double drivelen;			 //curve length, after which sliding is stopped (body falls down)
	double driveacceltime;			 //curve length, after which sliding is stopped (body falls down)
	int driveflag;         //if set to 1, driving is changed into completly free motion
	double drivenextstart; //time after which body is started again
};





//Constraint: two Cylindrical bodies (rigid/flexible, axis = X) may be in contact at 1 point
//leads to nonlinear constraint equations
class CylindricalContact: public Constraint
{
public:
	CylindricalContact(MBS* mbsi):Constraint(mbsi), loccoords(), dpdq() 
	{
	};
	//body 1 (element number en1), bodies 2-n (element numbers in en2, initial index in2)
	//lc1.X() and lc2.X() are proper initial X-coords!!
	//elemind from 1 to en2.Length()
	//contactdist is the contact distance where the linear contact force starts, usually the sum of both radii
	//Has one contact parameter, lc1.X() and lc2.X() are internally resolved due to nonlinear equations
	CylindricalContact(MBS* mbsi, int en1, const TArray<int>& en2, int in2, 
		const Vector3D& lc1, const Vector3D& lc2, double contactdist, double contactstiff,
		double cdim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq() 
	{	
		x_init = Vector(SS()); //lambda-contact
		GetCol() = coli;
		draw_dim.X() = cdim;
		AddElementCoord(en1, lc1);
		AddElementCoord(en2(1), lc2);
		for (int i=2; i <= en2.Length(); i++)
		{
			AddElementCoord(en2(i), lc2);
		}
		elemind = in2;
		xoff = 0;
		nlstepcnt = 0;
		nlfinit(1) = lc1.X();
		nlfinit(2) = lc2.X();
		cdist = contactdist;
		cstiff = contactstiff;
		switchvar = 0; //1 is contact, 0 is free
		elementname = GetElementSpec();
	};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new CylindricalContact(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const CylindricalContact& ce = (const CylindricalContact&)e;
		loccoords = ce.loccoords;
		dpdq = ce.dpdq;
		elemind = ce.elemind;
		xoff = ce.xoff;
		nlstepcnt = ce.nlstepcnt;
		nlfinit = ce.nlfinit;
		cdist = ce.cdist;
		cstiff = ce.cstiff;
		switchvar = ce.switchvar;
	}

	virtual void AddElementCoord(int en, const Vector3D& loccoord)
	{
		AddElement(en);
		loccoords.Add(loccoord);
	}

	virtual void EvalG(Vector& f, double t);
	virtual void EvalF(Vector& f, double t);

	virtual double PostNewtonStep(double t);
	virtual void PostprocessingStep();

	virtual Vector2D MinDist();
	//compute the residual of the mindist function
	virtual void ResMinDist(const Vector2D& x, Vector2D& res);
	virtual void JacMinDist(const Vector2D& x, Matrix3D& jac);

	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size
	virtual int IS() const {return 1;};
	virtual int ES() const {return 0;};

	virtual int Dim() const {return 3;}

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;
	virtual Vector3D GetRefPos() const;

	virtual void DrawElement();

protected:
	TArray<Vector3D> loccoords;
	int elemind;
	double xoff;
	int nlstepcnt;
	Vector2D nlfinit;
	double cdist, cstiff;
	int switchvar;

	Matrix dpdq; //temporary element, no need to copy???

};


//Constraint: constrain 3 coordinates x1 to other 3 coordinates x2 to be x1*x2=y
//leads to nonlinear constraint equations
class GeneralizedAngleConstraint: public Constraint
{
public:
	GeneralizedAngleConstraint(MBS* mbsi):Constraint(mbsi), loccoords(), dpdq() 
	{
	};
	//element1: en1, coordinates at (c1-1)*3+1
	//element2: en2, coordinates at (c2-1)*3+1
	//calpha means the cosine of the included angle
	//calpha=0 -> both vectors are perpendicular
	//calpha=1 and en1=en2 -> norm of vector is 1 (length is 1) 
	GeneralizedAngleConstraint(MBS* mbsi, int en1, int en2, int c1i, int c2i, double calphai,
		double cdim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq() 
	{	
		penalty = 0;
		penaltyfact = 1e-4;

		x_init = Vector(SS());
		GetCol() = coli;
		draw_dim.X() = cdim;
		AddElement(en1);
		AddElement(en2);
		loccoords.Add(Vector3D(0,0,0));
		loccoords.Add(Vector3D(0,0,0));

		c1 = c1i;
		c2 = c2i;
		calpha = calphai;

		elementname = GetElementSpec();
	};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new GeneralizedAngleConstraint(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const GeneralizedAngleConstraint& ce = (const GeneralizedAngleConstraint&)e;
		loccoords = ce.loccoords;
		dpdq = ce.dpdq;
		c1 = ce.c1;
		c2 = ce.c2;
		calpha = ce.calpha;
		penalty = ce.penalty;
		penaltyfact = ce.penaltyfact;
	}

	virtual void EvalG(Vector& f, double t);

	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size
	virtual int IS() const {if (!penalty) return 1; else return 0;};
	virtual int Dim() const {return 3;}

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;

	virtual void DrawElement();

protected:
	TArray<Vector3D> loccoords;
	int c1, c2;
	double calpha;

	int penalty;
	double penaltyfact;

	Matrix dpdq; //temporary element, no need to copy???

};




//Prescribed Angle constraint for example of spin-up

class PrescribedAngleConstraint: public Constraint
{
public:
	PrescribedAngleConstraint(MBS* mbsi):Constraint(mbsi), loccoords(), dpdq()
	{
		mode = 1;
	};
	PrescribedAngleConstraint(MBS* mbsi, int en1, const Vector3D& lp1,
		double cdim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq()
	{
		mode = 1;
		x_init = Vector(SS());
		GetCol() = coli;
		draw_dim.X() = cdim;
		AddElement(en1);
		loccoords.Add(lp1);
		elementname = GetElementSpec();
	};
	PrescribedAngleConstraint(MBS* mbsi, int en1, const Vector2D& lp1,
		double cdim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq()
	{
		mode = 2;
		x_init = Vector(SS());
		GetCol() = coli;
		draw_dim.X() = cdim;
		AddElement(en1);
		Vector3D l3d(lp1.X(), lp1.Y(), 0.);
		loccoords.Add(l3d);
		elementname = GetElementSpec();
	};
	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new PrescribedAngleConstraint(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const PrescribedAngleConstraint& ce = (const PrescribedAngleConstraint&)e;
		loccoords = ce.loccoords;
		dpdq = ce.dpdq;
		mode = ce.mode;
	}
	virtual void EvalG(Vector& f, double t);
	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);
	//implicit (algebraic) size
	virtual int IS() const {return 1;};
	virtual int Dim() const {if (mode == 1) return 3; else return 2;}
	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;
	virtual void DrawElement();
	virtual double Angle(double t)
	{
		const double Ts = 15;
		const double os = 4;
		//double theta;
		if (t < Ts)
		{
			return os/Ts*(Sqr(t)/2.+Sqr(Ts/(2.*MY_PI))*(cos(2.*MY_PI*t/Ts)-1.));
		}
		else
		{
			return os*(t-Ts/2.);
		}
	}
	virtual double Anglep(double t)
	{
		const double Ts = 15;
		const double os = 4;
		//double theta;
		if (t < Ts)
		{
			return os/Ts*(t+Sqr(Ts/(2.*MY_PI))*(-2.*MY_PI/Ts*sin(2.*MY_PI*t/Ts)));
		}
		else
		{
			return os;
		}
	}
	//rotation around y-axis: initially (phi=0) normal to x
	virtual Vector3D GetRotVector(double phi)
	{
		//return Vector3D(cos(phi),sin(phi),0); //start vector along x-axis, does not work for cable
		return Vector3D(-sin(phi),cos(phi),0); //old, start along y-axis
	}
	virtual Vector3D GetRotVectorp(double phi, double phip)
	{
		//return Vector3D(-sin(phi)*phip,cos(phi)*phip,0); //does not work for cable
		return Vector3D(-cos(phi)*phip,-sin(phi)*phip,0); //old, y-axis
	}
	virtual Vector3D GetNormalRotVector() const {return Vector3D(1,0,0);}
	//virtual Vector3D GetNormalRotVector() const {return Vector3D(0,1,0);} //does not work for cable
protected:
	TArray<Vector3D> loccoords;
	Matrix dpdq; //temporary element, no need to copy???
	int mode; //1=3D, 2=2D
};



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Angular velocity
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//Prescribed Angular Velocity for 3 axes parallel to global coordinate system axes
class PrescribedAngularVel: public Constraint
{
public:
	PrescribedAngularVel(MBS* mbsi):Constraint(mbsi), loccoords(), dpdq() 
	{
	};
	PrescribedAngularVel(MBS* mbsi, int en1, const Vector3D& lp1,
		double cdim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq() 
	{	
		ptype = 0;
		x_init = Vector(SS());
		GetCol() = coli;
		draw_dim.X() = cdim;
		AddElement(en1);
		loccoords.Add(lp1);
		elementname = GetElementSpec();
	};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new PrescribedAngularVel(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const PrescribedAngularVel& ce = (const PrescribedAngularVel&)e;
		loccoords = ce.loccoords;
		dpdq = ce.dpdq;

		ptype = ce.ptype;
		data1 = ce.data1;
		data2 = ce.data2;
		data3 = ce.data3;
	}

	virtual void EvalG(Vector& f, double t);

	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size
	virtual int IS() const {return 3;};
	virtual int Dim() const {return 3;}

	virtual void SetPrescribedHarmonic(const Vector3D& omega, const Vector3D& phase,const Vector3D& amplitude)
	{
		ptype = 1;
		data1.SetLen(3);
		data2.SetLen(3);
		data3.SetLen(3);
		data1(1) = omega.X(); data1(2) = phase.X(); data1(3) = amplitude.X();
		data2(1) = omega.Y(); data2(2) = phase.Y(); data2(3) = amplitude.Y();
		data3(1) = omega.Z(); data3(2) = phase.Z(); data3(3) = amplitude.Z();
	}

	virtual void GetPrescribedOmega(Vector3D& prescribed_omega, double t) const;

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;

	virtual void DrawElement();

protected:
	TArray<Vector3D> loccoords;
	Matrix dpdq; //temporary element, no need to copy???

	int ptype;				//1=harmonic functions
	Vector data1;	//omega, phaseshift and amplitude, ang_vel1
	Vector data2;	//omega, phaseshift and amplitude, ang_vel2
	Vector data3;	//omega, phaseshift and amplitude, ang_vel3
};





//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Gravity constraint ...
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class GravityConstraint: public Constraint
{
public:
	GravityConstraint(MBS* mbsi):Constraint(mbsi), loccoords(), dpdq() 
	{
	};
	//element1: en1, element2: en2
	GravityConstraint(MBS* mbsi, int en1, int en2, double Gi):Constraint(mbsi), loccoords(), dpdq() 
	{	
		x_init = Vector(SS());

		AddElement(en1);
		AddElement(en2);
		loccoords.Add(Vector3D(0,0,0));
		loccoords.Add(Vector3D(0,0,0));

		gconstant = Gi;

		elementname = GetElementSpec();
	};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new GravityConstraint(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const GravityConstraint& ce = (const GravityConstraint&)e;
		loccoords = ce.loccoords;
		dpdq = ce.dpdq;
		gconstant = ce.gconstant;
	}

	virtual void EvalG(Vector& f, double t) {};

	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size
	virtual int IS() const {return 0;};
	virtual int Dim() const {return 3;}

	virtual void DrawElement() {};

protected:
	TArray<Vector3D> loccoords;
	double gconstant;

	Matrix dpdq; //temporary element, no need to copy???

};


//Constraint: Contact with ground
//may lead to nonlinear constraint equations
class GroundContact: public Constraint
{
public:
	GroundContact(MBS* mbsi):Constraint(mbsi), loccoords(), dpdq() 
	{
	};
	GroundContact(MBS* mbsi, int en1, GeomGround* ground1, double contactdist, double contactstiff,
		double cdim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq() 
	{	
		x_init = Vector(SS()); //lambda-contact
		GetCol() = coli;
		draw_dim.X() = cdim;
		AddElementCoord(en1, Vector3D(0,0,0));

		ground = ground1;

		cdist = contactdist;
		cstiff = contactstiff;
		switchvar = 0; //1 is contact, 0 is free
		usedquad = 0;
		elementname = GetElementSpec();
	};

	GroundContact(MBS* mbsi, int en1, Vector3D locpos, GeomGround* ground1, double contactdist, double contactstiff,
		double cdim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq() 
	{	
		x_init = Vector(SS()); //lambda-contact
		GetCol() = coli;
		draw_dim.X() = cdim;
		AddElementCoord(en1, locpos);

		ground = ground1;

		cdist = contactdist;
		cstiff = contactstiff;
		switchvar = 0; //1 is contact, 0 is free
		usedquad = 0;
		elementname = GetElementSpec();
	};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new GroundContact(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const GroundContact& ce = (const GroundContact&)e;
		loccoords = ce.loccoords;
		dpdq = ce.dpdq;

		ground = ce.ground;

		cdist = ce.cdist;
		cstiff = ce.cstiff;
		switchvar = ce.switchvar;
		nlstepcnt = ce.nlstepcnt;
		usedquad = ce.usedquad;
	}

	virtual void AddElementCoord(int en, const Vector3D& loccoord)
	{
		AddElement(en);
		loccoords.Add(loccoord);
	}

	virtual double PostNewtonStep(double t);
	virtual void PostprocessingStep();

	virtual void EvalG(Vector& f, double t);

	//compute minimal distance to ground, negative if penetration ...
	virtual double MinDist(Vector3D& n, int& nquad);

	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size
	virtual int IS() const {return 1;};
	virtual int ES() const {return 0;};

	virtual int Dim() const {return 3;}

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;
	virtual Vector3D GetRefPos() const;

	virtual void DrawElement();

	virtual Box3D GetElementBoxD() const;

protected:
	TArray<Vector3D> loccoords;

	int nlstepcnt;
	double cdist, cstiff;
	int switchvar;
	int usedquad;

	GeomGround* ground;

	Matrix dpdq; //temporary element, no need to copy???

};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// CircleContact 2D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//a point can not move inside a given circle in space
class CircleContact2D: public Constraint
{
public:
	//Pos2DConstraint():loccoords(), dpdq() {mbs = NULL; };
	CircleContact2D(MBS* mbsi):Constraint(mbsi), loccoords(), dpdq() 
	{
	};
	CircleContact2D(MBS* mbsi, int en1, const Vector2D& lc1, 
		const Vector2D& pglob, double radius, double contact_stiffness, const Vector3D& cdim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq() 
	{	
		DOFmode = 0;
		r = radius;
		//c_contact = contact_stiffness;
		SetPenaltyStiffness(contact_stiffness);
		iscontact = 0;
		nlstepcnt = 0;
		contactmode = 0;
		gapp_init_min = 1e-6;
		gapp_init = gapp_init_min;

		x_init = Vector(SS());
		GetCol() = coli;
		//draw_dim = cdim;
		SetDrawSizeScalar(cdim.X());
		SetDrawSizeResolution((int)(cdim.Y()));
		p_global = pglob;
		AddElementCoord(en1, lc1);
		elementname = GetElementSpec();
	};

	CircleContact2D(MBS* mbsi, int en1, int nodenum, 
		const Vector2D& pglob, double radius, double contact_stiffness, const Vector3D& cdim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq() 
	{	
		DOFmode = 1;
		Vector2D lc1(nodenum,0.);

		r = radius;
		//c_contact = contact_stiffness;
		SetPenaltyStiffness(contact_stiffness);
		iscontact = 0;
		contactmode = 0;
		nlstepcnt = 0;
		gapp_init_min = 1e-6;
		gapp_init = gapp_init_min;

		x_init = Vector(SS());
		GetCol() = coli;
		//draw_dim = cdim;
		SetDrawSizeScalar(cdim.X());
		SetDrawSizeResolution((int)(cdim.Y()));
		p_global = pglob;
		AddElementCoord(en1, lc1);
		elementname = GetElementSpec();
	};


	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new CircleContact2D(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const CircleContact2D& ce = (const CircleContact2D&)e;
		loccoords = ce.loccoords;
		p_global = ce.p_global;
		dpdq = ce.dpdq;
		r = ce.r;
		//c_contact = ce.c_contact;
		iscontact = ce.iscontact;
		nlstepcnt = ce.nlstepcnt;
		gapp_init = ce.gapp_init;
		gapp_init_min = ce.gapp_init_min;

		DOFmode = ce.DOFmode;
		contactmode = ce.contactmode;
	}

	virtual void AddElementCoord(int en, const Vector2D& loccoord)
	{
		AddElement(en);
		loccoords.Add(loccoord);
	}

	virtual double PostNewtonStep(double t);
	virtual void PostprocessingStep();

	virtual void EvalG(Vector& f, double t);

	//To be replaced in derived class
	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size
	virtual int IS() const 
	{
		return 1;
	};
	virtual void SetContactMode(int i) {contactmode = i;}

	virtual int Dim() const {return GetElem(1).Dim();}  //has 3D position???

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;

	virtual void DrawElement();
	int CircleContact2D::GetDrawSizeResolution()
	{
		if(draw_dim.Z()==-1)
			return (int)(GetMBS()->GetDOption(173));
		else
			return (int)draw_dim.Z();
	}

	void CircleContact2D::SetDrawSizeResolution(int drawdim_resolution)
	{
		draw_dim.Z()= (double)drawdim_resolution;
	}

protected:
	TArray<Vector2D> loccoords;
	int DOFmode; //if DOFmode = 1, then loccoords(1).X() is the local node number in the element

	Vector2D p_global; //position of fixed circle
	double r; //radius of fixed circle
	//double c_contact; //contact stiffness
	double gapp_init; //initial contact velocity
	double gapp_init_min;
	int contactmode; //0= contact, 1=no loss of contact (always contact stiffness), 2=no loss of contact, lagrange mult., no stiffness

	int iscontact;
	int nlstepcnt;

	Matrix dpdq; //temporary element, no need to copy???
};



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// IntegralPosConstraint2D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//a point can not move inside a given circle in space
class IntegralPosConstraint2D: public Constraint
{
public:
	//Pos2DConstraint():loccoords(), dpdq() {mbs = NULL; };
	IntegralPosConstraint2D(MBS* mbsi):Constraint(mbsi), loccoords(), dpdq(), weights(), elements_nodouble(), s_body1(), s_body2()
	{
		maxmoment = 0;
	};

	//contraint int. local points with ground
	//must be different elements!!!
	IntegralPosConstraint2D(MBS* mbsi, const TArray<int>& elnums1, const TArray<Vector2D>& loccoords1, const TArray<double>& weights1, 
		const Vector2D& pglob, const Vector3D& cdim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq(), weights(), elements_nodouble(), s_body1(), s_body2()
	{	
		DOFmode = 0; //with local coordinates
		maxmoment = 0;
		x_init = Vector(SS());
		GetCol() = coli;
		draw_dim = cdim;
		p_global = pglob;

		for (int i=1; i <= elnums1.Length(); i++)
		{
			AddElementCoord(elnums1(i), loccoords1(i));
			weights.Add(weights1(i));
		}

		//for element jacobian, no double elements:
		elements_nodouble.SetLen(0);
		for (int i = 1; i <= elnums1.Length(); i++)
		{	
			if (!Find(elnums1(i), elements_nodouble))
			{
				elements_nodouble.Add(elnums1(i));
			}
		}
		elementname = GetElementSpec();
	};

	//integral constraint, node numbers with ground
	//must be different elements!!!!
	IntegralPosConstraint2D(MBS* mbsi, const TArray<int>& elnums1, const TArray<int>& nodenums1, const TArray<double>& weights1,
		const Vector2D& pglob, const Vector3D& cdim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq(), weights(), elements_nodouble(), s_body1(), s_body2() 
	{	
		DOFmode = 1; //with node numbers
		maxmoment = 0;
		x_init = Vector(SS());
		GetCol() = coli;
		draw_dim = cdim;
		p_global = pglob;

		for (int i=1; i <= elnums1.Length(); i++)
		{
			AddElementCoord(elnums1(i), Vector2D(nodenums1(i),0));
			weights.Add(weights1(i));
		}

		//for element jacobian, no double elements:
		elements_nodouble.SetLen(0);
		for (int i = 1; i <= elnums1.Length(); i++)
		{	
			if (!Find(elnums1(i), elements_nodouble))
			{
				elements_nodouble.Add(elnums1(i));
			}
		}
		elementname = GetElementSpec();
	};

	//int. constraints, node numbers with ground
	//same elements
	IntegralPosConstraint2D(MBS* mbsi, int elnum, const TArray<int>& nodenums1, const TArray<double>& weights1,
		const Vector2D& pglob, const Vector3D& cdim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq(), weights(), elements_nodouble(), s_body1(), s_body2() 
	{
		DOFmode = 2; //with node numbers, only one element
		maxmoment = 0;
		x_init = Vector(SS());
		GetCol() = coli;
		draw_dim = cdim;
		p_global = pglob;

		for (int i=1; i <= nodenums1.Length(); i++)
		{
			AddElementCoord(elnum, Vector2D(nodenums1(i),0));
			weights.Add(weights1(i));
		}

		//for element jacobian, no double elements:
		elements_nodouble.Add(elnum);
	};

	//integral constraint, first element with node numbers, second element with local position
	//same elements
	IntegralPosConstraint2D(MBS* mbsi, int elnum, int elnum2, const TArray<int>& nodenums1, const TArray<double>& weights1,
		const Vector2D& ploc2, const Vector3D& cdim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq(), weights(), elements_nodouble(), s_body1(), s_body2() 
	{	
		DOFmode = 3; //with node numbers, 2 elements, second element with local position
		maxmoment = 0;
		x_init = Vector(SS());
		GetCol() = coli;
		draw_dim = cdim;

		for (int i=1; i <= nodenums1.Length(); i++)
		{
			AddElementCoord(elnum, Vector2D(nodenums1(i),0));
			weights.Add(weights1(i));
		}
		AddElementCoord(elnum2, ploc2);
		weights.Add(1);

		//for element jacobian, no double elements:
		elements_nodouble.Add(elnum);
		elements_nodouble.Add(elnum2);
	};


	//integral constraint, first elements with node numbers, second element with local positions
	//nodes are sorted, weigths are computed from distance of nodes
	IntegralPosConstraint2D(MBS* mbsi, const TArray<int>& elnums1, int elnum2, const TArray<int>& nodenums1,
		const TArray<Vector2D>& nodepos1, const TArray<Vector2D>& loccoords2, int maxmomentI,
		const Vector3D& cdim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq(), weights(), elements_nodouble(), s_body1(), s_body2() 
	{
		DOFmode = 4; //2 element sets, first element set with node numbers, second element set with local position
		x_init = Vector(SS());
		GetCol() = coli;
		draw_dim = cdim;

		npoints1 = nodenums1.Length();
		npoints2 = loccoords2.Length();
		maxmoment = maxmomentI;

		assert(elnums1.Length() == nodenums1.Length());
		weights.SetLen(npoints1 + npoints2); //add second element later


		length1 = 0;
		s_body1.SetLen(npoints1);
		s_body1(1) = 0;
		for (int i=1; i < nodenums1.Length(); i++)
		{
			length1 += Dist(nodepos1(i), nodepos1(i+1));
			s_body1(i+1) = length1;
			weights(i) = 0;
		}
		weights(nodenums1.Length()) = 0;

		if (length1 == 0) length1 = 1;

		for (int i=1; i <= nodenums1.Length(); i++)
		{
			s_body1(i) -= 0.5*length1;
			AddElementCoord(elnums1(i), Vector2D(nodenums1(i),0));
			if (i < nodenums1.Length())
			{
				double w = 0.5*Dist(nodepos1(i), nodepos1(i+1))/length1;
				weights(i) += w;
				weights(i+1) += w;
			}
		}

		//second body:
		length2 = 0;
		int noff = nodenums1.Length();
		s_body2.SetLen(npoints2);
		s_body2(1) = 0;
		for (int i=1; i < loccoords2.Length(); i++)
		{
			length2 += Dist(loccoords2(i), loccoords2(i+1));
			s_body2(i+1) = length2;
			weights(i+noff) = 0;
		}

		weights(loccoords2.Length()+noff) = 0;
		if (length2 == 0)
		{
			weights(loccoords2.Length()+noff) = 1; //only one pos!!!
			length2 = 1;
		}

		for (int i=1; i <= loccoords2.Length(); i++)
		{
			s_body2(i) -= 0.5*length2;
			AddElementCoord(elnum2, loccoords2(i));
			if (i < loccoords2.Length())
			{
				double w = 0.5*Dist(loccoords2(i), loccoords2(i+1))/length2;
				weights(i+noff) += w;
				weights(i+1+noff) += w;
			}
		}

		//for element jacobian, no double elements:
		elements_nodouble.SetLen(0);
		for (int i = 1; i <= elnums1.Length(); i++)
		{	
			if (!Find(elnums1(i), elements_nodouble))
			{
				elements_nodouble.Add(elnums1(i));
			}
		}
		elements_nodouble.Add(elnum2);

	};


	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new IntegralPosConstraint2D(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const IntegralPosConstraint2D& ce = (const IntegralPosConstraint2D&)e;

		elements_nodouble = ce.elements_nodouble;

		loccoords = ce.loccoords;
		weights = ce.weights;
		p_global = ce.p_global;
		dpdq = ce.dpdq;
		DOFmode = ce.DOFmode;

		npoints1 = ce.npoints1;
		npoints2 = ce.npoints2;
		maxmoment = ce.maxmoment;
		length1 = ce.length1;
		length2 = ce.length2;

		s_body1 = ce.s_body1;
		s_body2 = ce.s_body2;

	}

	virtual void AddElementCoord(int en, const Vector2D& loccoord)
	{
		AddElement(en);
		loccoords.Add(loccoord);
	}

	virtual void LinkToElements() 
	{
		//for integral constraint???

		if (DOFmode <= 3) 
		{
			GetElem(1).AddConstraint(this,1);
			if (elements.Length() >= 2 && DOFmode == 3)
			{
				GetElem(elements.Length()).AddConstraint(this, elements.Length());
			}
		}
		else if (DOFmode == 4) 
		{
			for (int i = 1; i <= npoints1; i++)
			{
				GetElem(i).AddConstraint(this,i);
			}
			GetElem(npoints1+1).AddConstraint(this,npoints1+1);
		}

	}

	//elements not doubled: only for certain elements!!!
	virtual int NE_nodouble() const {return elements_nodouble.Length();}
	virtual const Element& GetElem_nodouble(int i) const {return mbs->GetElement(elements_nodouble(i));}
	virtual Element& GetElem_nodouble(int i) {return mbs->GetElement(elements_nodouble(i));}

	virtual void EvalG(Vector& f, double t);

	//To be replaced in derived class
	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size
	virtual int IS() const {return 2+2*maxmoment;};

	virtual int Dim() const {return GetElem(1).Dim();}  //has 3D position???

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;

	virtual void DrawElement();

protected:
	TArray<int> elements_nodouble; //every element is added only once, for automatic jacobian!!!

	TArray<Vector2D> loccoords;
	TArray<double> weights; //weights for integral, if points are not equally distributed
	int DOFmode; //if DOFmode = 1, then loccoords(1).X() is the local node number in the element

	int npoints1;  //for DOFmode == 4, number of points body 1
	int npoints2;  //for DOFmode == 4, number of points body 2
	int maxmoment; //how many moments to take?, DOFmode == 4
	double length1, length2; //length of integral for body 1 and 2
	TArray<double> s_body1; //body 1: curve-distances from center of moment, for moments > 0!
	TArray<double> s_body2; //body 2: curve-distances from center of moment, for moments > 0!

	Vector2D p_global; //position of fixed circle

	Matrix dpdq; //temporary element, no need to copy???
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Kardan Rotational Spring-Damper-Actor Element, acts on Kardan angles!!! 3D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//description
class KardanSDRotActor: public Constraint
{
public:
	KardanSDRotActor(MBS* mbsi):Constraint(mbsi), loccoords(), dpdq()
	{
		groundjoint = 0;
	};


	void SetKardanSDRotActor(int en1, const Vector3D& lc1, const Vector3D& spring_stiffnessi, 
				Vector3D cdim, const Vector3D& coli)
	{	
		groundjoint = 1;
		x_init = Vector(0);
		GetCol() = coli;
		//draw_dim = cdim;
		SetDrawSizeScalar(cdim.Y());
		SetDrawSizeCylinderLength(cdim.X());
		AddElementCoord(en1, lc1);
		//spring_stiffness = spring_stiffnessi;
		SetPenaltyStiffness3(spring_stiffnessi);
		elementname = GetElementSpec();
	};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new KardanSDRotActor(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const KardanSDRotActor& ce = (const KardanSDRotActor&)e;
		loccoords = ce.loccoords;
		//p_global = ce.p_global;
		dpdq = ce.dpdq;

		spring_stiffness3 = ce.spring_stiffness3;
		groundjoint = ce.groundjoint;
	}

	virtual const char* GetElementSpec() const {return "KardanRotActor";}
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	virtual void AddElementCoord(int en, const Vector3D& loccoord)
	{
		AddElement(en);
		loccoords.Add(loccoord);
	}

	virtual void EvalG(Vector& f, double t) {};
	virtual void EvalF(Vector& f, double t) {};
	virtual void EvalF2(Vector& f, double t) {};

	virtual int IsGroundJoint() const {return groundjoint;}

	virtual Vector3D ComputeTorque(double t) const;

	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size
	virtual int IS() const {return 0;};

	virtual int Dim() const {return 3;}

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;

	virtual void DrawElement();

	double KardanSDRotActor::GetDrawSizeCylinderLength()
	{
		if(draw_dim.Y()==-1)
			return GetMBS()->GetDOption(172);
		else
			return draw_dim.Y();
	}
	void KardanSDRotActor::SetDrawSizeCylinderLength(double drawdim_cyllength)
	{
		draw_dim.Y()=drawdim_cyllength;
	}

	virtual double GetActorForce(double computation_time, int dir=0) const 
	{
		return (ComputeTorque(computation_time))(dir);
	} 
	// 3D spring stiffness
	//return 3 stiffness parameters for penalty formulation
	virtual Vector3D GetPenaltyStiffness3() const {return spring_stiffness3;}
	virtual double GetPenaltyStiffness3(int i) const {return spring_stiffness3[i-1];}
	//set 3 stiffness parameters for penalty formulation
	virtual void SetPenaltyStiffness3(Vector3D stiffness3i) {spring_stiffness3 = stiffness3i;}

protected:
	TArray<Vector3D> loccoords;

	Matrix dpdq; //temporary element, no need to copy???

	int groundjoint; //==1 if groundjoint
	Vector3D spring_stiffness3;
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Surface
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//no equations
class Surface: public Constraint
{
public:
	Surface(MBS* mbsi):Constraint(mbsi), face_nrs(), surf_node_nrs_nodouble()
	{	};

	Surface(const Surface& sur): Constraint(sur.mbs), face_nrs(0), surf_node_nrs_nodouble(0)
	{
		CopyFrom(sur);
	};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new Surface(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const Surface& ce = (const Surface&)e;
		face_nrs = ce.face_nrs;
		surf_node_nrs_nodouble = ce.surf_node_nrs_nodouble;
	}

	// set functions
	///*Übergabeparameter anpassen*/
	//void SetSurface(int en1, int fn1, double r) //calculates TArray<int> surf_node_nrs
	
	//{
	//	surf_node_nrs_nodouble.SetLen(0);
	//  

	//	x_init = Vector(SS());
	//	AddElementCoord(en1, lc1);
	//	draw_dim.X() = r;

	//	{
	//		UpdateSurfNodeNrList(en, fn);
	//	}
	//};

	virtual const char* GetElementSpec() const {return "Surface";}
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	virtual int IS() const {return 0;};  //implicit (algebraic) size (RL) anzahl lagrangeparameter=0 und sos=0
	virtual int Dim() const {return GetElem(1).Dim();}  //has 3D position???

	virtual void AddFaceToSurface(int en, int fn)
	{
		//element_nrs.Add(en);
		elements.Add(en);
		face_nrs.Add(fn);
		UpdateSurfNodeNrList(en, fn);
	}

	virtual void ComputeSurfNodeNrList()
	{
		int en=0, fn=0;
		for (int i=1; i<=face_nrs.Length(); i++)
		{
			en=elements(i);/*element_nrs(i);*/
			fn=face_nrs(i);
			UpdateSurfNodeNrList(en, fn);
		}
	}

	//adds global nodenumbers if not already contained
	virtual void UpdateSurfNodeNrList(int en, int fn); 

	virtual TArray<int> GetSurfaceNodeList(){return surf_node_nrs_nodouble;}

	//virtual void SetPenalty(int flag){penalty = flag;} // use penalty coordinate constraint (RL)
	//virtual int GetPenalty(){return penalty;} // use penalty coordinate constraint (RL)

	//virtual void SetSpringStiffness(double val){spring_stiffness = val;} // use penalty coordinate constraint (RL)
	//virtual double GetSpringStiffness(){return spring_stiffness;} // use penalty coordinate constraint (RL)
	//
	////get a reference position of the body in 3d
	//virtual Vector3D GetRefPosD() const 
	//{
	//	return GetElem(1).GetRefPosD();
	//}
	
	////wandelt lokale in globale Knotennummer um
	//virtual int GetGlobalConstrainedDOF() const;

	virtual void DrawElement();

protected:
	//TArray<int> element_nrs;
	TArray<int> face_nrs;
	TArray<int> surf_node_nrs_nodouble; //list of global node numbers of involved nodes

};



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// RigidLinkConstraint 3D (RL)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// constraint for constant distance between two points (global/local body point or two points of different bodies)
class RigidLinkConstraint: public Constraint
{
public:
	RigidLinkConstraint(MBS* mbsi):Constraint(mbsi), loccoords(), dpdq() 
	{
	};
	void SetRigidLinkConstraint(int en1, const Vector3D& lc1, 
		const Vector3D& pglob, const Vector3D& drawdimi, const Vector3D& coli)
	{	
		x_init = Vector(SS());
		GetCol() = coli;
		draw_dim = drawdimi;
		p_global = pglob;
		AddElementCoord(en1, lc1);
		elementname = GetElementSpec();
	};

	void SetRigidLinkConstraint(int en1, int en2, 
		const Vector3D& lc1, const Vector3D& lc2, const Vector3D& drawdimi, const Vector3D& coli) 
	{	
		x_init = Vector(SS());
		GetCol() = coli;
		draw_dim = drawdimi;
		AddElementCoord(en1, lc1);
		AddElementCoord(en2, lc2);
		elementname = GetElementSpec();
	};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new RigidLinkConstraint(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const RigidLinkConstraint& ce = (const RigidLinkConstraint&)e;
		loccoords = ce.loccoords;
		p_global = ce.p_global;
		dpdq = ce.dpdq;
		dist = ce.dist;
	}

	virtual void Initialize() 
	{
		if (elements.Length()==1) 
		{
			dist = (GetElem(1).GetPos(loccoords(1))-p_global).Norm(); 
		}
		if (elements.Length()==2) 
		{
			dist = (GetElem(1).GetPos(loccoords(1))-GetElem(2).GetPos(loccoords(2))).Norm(); 
		}	
	};

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	virtual const char* GetElementSpec() const {return "RigidLinkConstraint";} //this name names MUST coincide with the type names defined in MBS::MBS() !!!!

	virtual void AddElementCoord(int en, const Vector3D& loccoord)
	{
		AddElement(en);
		loccoords.Add(loccoord);
	}

	virtual void EvalG(Vector& f, double t);

	virtual Vector3D ComputeForce(double t) const;    
	virtual Vector3D ComputeForceDirection() const{assert(0 && "RigidLinkConstraint: ComputeForceDirection called (niy)!");return Vector3D(0.);}; 
	virtual Vector3D ComputeForceDirectionD() const{assert(0 && "RigidLinkConstraint: ComputeForceDirectionD called (niy)!");return Vector3D(0.);}; 
	virtual double GetActorForce(double computation_time, int dir=0) const;
	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size
	virtual int IS() const {return 1;};

	virtual int Dim() const {return 3;}

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;

	virtual void DrawElement();

protected:
	TArray<Vector3D> loccoords;
	Vector3D p_global;
	Matrix dpdq; // temporary element, no need to copy???
	double dist; // distance between bodies
};



//$ DR 2011-04:[ MultiNodalConstraint added
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// MultiNodalConstraint
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// name: MultiNodalConstraint (Constraint)
// short description: two sets of nodes are constrained to each other
// available formulations: penalty
// type of constraint equation: 
// index formulation available: 
// ground constraint: yes
// element-element constraint: yes
// development status: complete, tested for CMS + GCMS
// long description: two sets of nodes are constrained to each other, a center of each set of nodes is computed and the displacement of this two center nodes is used to compute the force.
// class variables:
//	- nodes_per_element					these are the nodes of elements, the TArray nodes of NodalConstraint is not used in this class
//	- weight_per_node[2]				these are the weights of the nodes to compute the center and to distribute the force
//	- center_offset							constraint is set relative to offset (not corotating); for ground-joint, this is the ground position
//	- constant_node_dpdqT[2]		temporary transposed node_dpdq matrix
//	- average_node_dpdq[2]			temporary node_dpdq matrix
//	- use_constant_node_dpdq		flag which indicates that dpdq is constant (e.g. for GCMS or solid FE) do not use for CMS !!!
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


class MultiNodalConstraint: public NodalConstraint
{
public:
	MultiNodalConstraint(MBS* mbsi):NodalConstraint(mbsi) {}

	MultiNodalConstraint(const MultiNodalConstraint& ct): NodalConstraint(ct.mbs)
	{
		CopyFrom(ct);
	};

	// a set of nodes (nodelist1) is constrained to ground, cdirlist1 defines the constrained (global) degrees of freedom
	// the displacement of the center of the nodes in nodelist1 is computed; due to the displacement a force is computed; this force is distributed to the nodes
	// aground is not in use up to now --> position of ground is the initial center of nodelist1
	void SetLocalMultiNodalConstraintSpringDamper(int elem1, TArray<int> &nodelist1, TArray<int> &cdirlist1, Vector3D drawsize, int penaltyi = 0, Vector3D aground=Vector3D(0.0), double spring_stiffnessi=0, double damping_coeffi=0)
	{
		int elem2 = elem1;
		TArray<int> nodelist2; nodelist2.Add(0);
		TArray<int> cdirlist2; cdirlist2.Add(0);

		if(cdirlist1.Length()>1) cdirlist2.Add(0);
		if(cdirlist1.Length()>2) cdirlist2.Add(0);

		Vector3D offset=aground;

		SetLocalMultiNodalConstraintSpringDamper(elem1, nodelist1, cdirlist1, elem2, nodelist2, cdirlist2, drawsize, penaltyi, offset, spring_stiffnessi, damping_coeffi);
	}

	// two sets of nodes (nodelist1+nodelist2) are constrained to each other, cdirlist1 defines the constrained (global) degrees of freedom
	// the displacement of the centers of the nodes is computed; due to the displacement a force is computed; this force is distributed to the nodes
	// offset is not in use up to now --> offset is the initial difference of the centers of nodelist1 and nodelist2
	void SetLocalMultiNodalConstraintSpringDamper(int elem1, TArray<int> &nodelist1, TArray<int> &cdirlist1, int elem2, TArray<int> &nodelist2, TArray<int> &cdirlist2, Vector3D drawsize, int penaltyi = 0, Vector3D offset=0, double spring_stiffnessi=0, double damping_coeffi=0)
	{
		InitNC();

		if (penaltyi == 0) GetMBS()->UO() << "***** WARNING: penalty == 0 does not work in 'SetLocalMultiNodalConstraintSpringDamper'! ******************\n";

		SetPenaltyFormulation(1);

		nodes_per_element[0].Merge(nodelist1);
		nodes_per_element[1].Merge(nodelist2);

		loccoords.Add(cdirlist1(1));
		loccoords.Add(cdirlist2(1));

		elements.Add(elem1);

		if(cdirlist1.Length()>1) loccoords.Add(cdirlist1(2));
		if(cdirlist2.Length()>1) loccoords.Add(cdirlist2(2));
		if(cdirlist1.Length()>2) loccoords.Add(cdirlist1(3));
		if(cdirlist2.Length()>2) loccoords.Add(cdirlist2(3));

		if(nodes_per_element[1](1)!=0) elements.Add(elem2);

		//draw_dim = drawsize;
		SetDrawSizeScalar(drawsize.X());
		SetDrawSizeCircleSpheres(drawsize.Y());
		SetDrawSizeLineThicknessStar((int)(drawsize.Z()));

		SetPenaltyStiffness(spring_stiffnessi); //$ DR 2011-04-21: old code: spring_stiffness = spring_stiffnessi;
		damping_coeff = damping_coeffi;

		elementname = GetElementSpec();
		x_init = Vector(SS());

		velocity_constraint = 0; //default

		center_offset = offset;
	}

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new MultiNodalConstraint(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		NodalConstraint::CopyFrom(e);
		const MultiNodalConstraint& ce = (const MultiNodalConstraint&)e;

		nodes_per_element[0] = ce.nodes_per_element[0];
		nodes_per_element[1] = ce.nodes_per_element[1];

		weight_per_node[0] = ce.weight_per_node[0];
		weight_per_node[1] = ce.weight_per_node[1];

		center_offset = ce.center_offset;

		constant_node_dpdqT[0] = ce.constant_node_dpdqT[0];			
		constant_node_dpdqT[1] = ce.constant_node_dpdqT[1];	

		average_node_dpdq[0] = ce.average_node_dpdq[0];			
		average_node_dpdq[1] = ce.average_node_dpdq[1];	

		use_constant_node_dpdq = ce.use_constant_node_dpdq;	
	}

	virtual void Initialize();
	virtual const char* GetElementSpec() const {return "MultiNodalConstraint";}
	virtual void EvalF2(Vector& f, double t);
	virtual void DrawElement(); 

	double MultiNodalConstraint::GetDrawSizeCircleSpheres()
	{
		if(draw_dim.Y()==-1)
			return GetMBS()->GetDOption(172);
		else
			return draw_dim.Y();
	}
	void MultiNodalConstraint::SetDrawSizeCircleSpheres(double drawdim_circlespheres)
	{
		draw_dim.Y()=drawdim_circlespheres;
	}

	int MultiNodalConstraint::GetDrawSizeLineThicknessStar()
	{
		if(draw_dim.Z()==-1)
			return (int)(GetMBS()->GetDOption(173));
		else
			return (int)draw_dim.Z();
	}

	void MultiNodalConstraint::SetDrawSizeLineThicknessStar(int drawdim_line)
	{
		draw_dim.Z()= (double)drawdim_line;
	}

	// when computing the center of nodes, the nodes are weighted with this values
	// up to now just one option is implemented: equal weight for every node in nodelist
	virtual double GetWeight(int element_nr, int node_nr) const
	{
		if(!NodesPerElementLength(element_nr)) {assert(0);}
		return 1.0/(double)NodesPerElementLength(element_nr);
	}

	// return nodenumber stored in nodes_per_element
	virtual int NodeNum(int elem_nr, int node_nr) const 
	{
		return nodes_per_element[elem_nr-1](node_nr);
	}

	virtual int NodesPerElementLength(int elem_nr) const 
	{
		return nodes_per_element[elem_nr-1].Length();
	}

	// return node (nodenumbers of constrain are stored in nodes_per_element)
	virtual const Node& GetNode(int elem_nr, int node_nr)  const
	{
		return GetElem(elem_nr).GetNode(NodeNum(elem_nr,node_nr)); 
	}

	virtual Vector3D GetRefPosD() const 
	{
		return GetNode(1,1).Pos();	//GetPosD
	}

	// return position of node (nodenumbers of constrain are stored in nodes_per_element)
	virtual Vector3D GetElementNodePos(int elem_nr, int node_nr, int flagD=0) const;
	// return velocity of node (nodenumbers of constrain are stored in nodes_per_element)
	virtual Vector3D GetElementNodeVel(int elem_nr, int node_nr, int flagD=0) const;

	// compute center of nodes: all nodes are from element with defined elem_nr, nodes are stored in nodelist
	virtual Vector3D ComputeCenterOfNodes(int elem_nr, const TArray<int>& nodelist, int flagD=0) const
	{
		if (nodelist.Length()<1) { assert(0);}
		Vector3D center(0.,0.,0.);
		for(int i=1; i<= nodelist.Length(); i++)
		{
			center += GetWeight(elem_nr,i)*GetElementNodePos(elem_nr, i, flagD);		
		}	
		return center;
	}

	// compute average velocity of nodes: all nodes are from element with defined elem_nr, nodes are stored in nodelist
	virtual Vector3D ComputeCenterOfNodesVel(int elem_nr, const TArray<int>& nodelist, int flagD=0) const
	{
		if (nodelist.Length()<1) { assert(0);}
		Vector3D center(0.,0.,0.);
		for(int i=1; i<= nodelist.Length(); i++)
		{
			center += GetWeight(elem_nr,i)*GetElementNodeVel(elem_nr, i, flagD);		
		}
		return center;
	}

	// get the distance of the 2 centers of nodelists
	// without displacement (initial positions): center1 = center2 + center_offset
	virtual Vector3D GetDistance() const ;


	// get the relatice velocity of the 2 centers of nodelists
	virtual Vector3D GetRelVel() const ;


	virtual int GetUseConstantNodeDpdq() const { return use_constant_node_dpdq;}
	virtual void SetUseConstantNodeDpdq(int flag) { use_constant_node_dpdq=flag;}
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

protected:
	TArray<int> nodes_per_element[2];			//these are the nodes of elements, the TArray nodes of NodalConstraint is not used in this class
	TArray<int> weight_per_node[2];				//these are the weights of the nodes to compute the center and to distribute the force
	Vector3D center_offset;								//constraint is set relative to offset (not corotating); for ground-joint, this is the ground position
	Matrix constant_node_dpdqT[2];				//temporary transposed node_dpdq matrix
	Matrix average_node_dpdq[2];					//temporary node_dpdq matrix
	int use_constant_node_dpdq;						//flag which indicates that dpdq is constant (e.g. for GCMS or solid FE) do not use for CMS !!!
};

////////$ DR 2011-04:] MultiNodalConstraint added
//////


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Spherical: AveragingJoint
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: AveragingJoint (Constraint)
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

// DR and PG are responsible for this constraint
class MultiNodalSphericalJoint: public BasePointJoint  //$EDC$[beginclass,classname=MultiNodalSphericalJoint,parentclassname=BasePointJoint,texdescription="
// MultiNodalSphericalJoint is used like SphericalJoint, but each of the kinematic pairs is a set of body points (element number and local position) or global nodes.
// The global positions and velocities of these points are averaged.
// MultiNodalSphericalJoint constrains the weighted mid point M = \sum_j r_j w_j of a set of weights w_j and global positions r_j
// Such weighted mid point can be either constrained to a prescribed goal m (where m=0 would mean a ground constraint),
// or to a second weighted mid point M2 (prescribed in the same manner, as the original weighted mid point M).
//"]
{
public: 

	MultiNodalSphericalJoint(MBS* mbsi):BasePointJoint(mbsi)
	{	
		ElementDefaultConstructorInitialization();
	};

	MultiNodalSphericalJoint(const MultiNodalSphericalJoint& ct): BasePointJoint(ct.mbs)
	{
		ElementDefaultConstructorInitialization();
		CopyFrom(ct);
	};

	~MultiNodalSphericalJoint()
	{
	};

	void SetMultiNodalSphericalJoint_GlobalNodes_to_GlobalNodes(TArray<int> nodes1, TArray<int> nodes2, Vector3D stiffness, double damping = 0.0)
	{
		SetPos1ToGlobalNodes(nodes1);
		SetPos2ToGlobalNodes(nodes2);
		SetPenaltyStiffness3(stiffness);
		SetDampingCoeff(damping);
	}

	void SetMultiNodalSphericalJoint_GlobalNodes_to_GlobalPos(TArray<int> nodes1, Vector3D stiffness, double damping = 0.0, TArray<MathFunction*>& displacement = TArray<MathFunction*>(0) )
	{
		SetPos1ToGlobalNodes(nodes1);
		AutoComputeGround();

		if(displacement.Length() != 0) SetDisplacement(displacement);
		SetPenaltyStiffness3(stiffness);
		SetDampingCoeff(damping);
	}

	virtual void ElementDefaultConstructorInitialization();

	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute

	virtual Element* GetCopy()
	{
		Element* ec = new MultiNodalSphericalJoint(*this);
		return ec;
	}

	virtual void CopyFrom(const Element& e);

	virtual const char* GetElementSpec() const {return "MultiNodalSphericalJoint";}

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		
	//virtual int SetElementDataAuto(const ElementDataContainer& edc);
	virtual int SetElementDataAuto(ElementDataContainer& edc);
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual void DrawElement();

	virtual void EvalF2(Vector& f, double t);
	virtual Vector3D ComputeForce(double t) const;	

	virtual void LinkToElementsPenalty();

	//do not use BasePointJoint::Initialize because then the wrong SetPos2ToGlobalCoord is called
	void Initialize() 
	{
		if(auto_comp_ground)
		{
			SetPos2ToGlobalCoord(GetPos1());
		}

		// if general (scalar) penalty stiffness is used, set Stiffness-Vector
		if(GetPenaltyStiffness())
		{
			SetPenaltyStiffness3(Vector3D(GetPenaltyStiffness()));
		}
	}

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const 
	{
		return GetPos1();		// <--------- change it
	}

	Vector3D GetPos1() const {return GetAveragePosition(1);}
	Vector3D GetPos2() const {return GetAveragePosition(2);}
	Vector3D GetVel1() const {return GetAverageVelocity(1);}
	Vector3D GetVel2() const {return GetAverageVelocity(2);}
	void SetPos1ToGlobalNodes(TArray<int> &glob_node_nrs){SetPosToGlobalNodes(1, glob_node_nrs);}
	void SetPos2ToGlobalNodes(TArray<int> &glob_node_nrs){SetPosToGlobalNodes(2, glob_node_nrs);}

	// do not use this function!
	void SetPos1ToLocalNode(int elem_nr, int loc_node_nr)		{GetMBS()->UO(UO_LVL_err) << "ERROR: MultiNodalSphericalJoint: SetPos1ToLocalNode not implemented.";}
	// do not use this function!
	void SetPos2ToLocalNode(int elem_nr, int loc_node_nr)		{GetMBS()->UO(UO_LVL_err) << "ERROR: MultiNodalSphericalJoint: SetPos2ToLocalNode not implemented.";}
	// do not use this function!
	void SetPos1ToLocalCoord(int elem_nr, Vector3D loccoord){GetMBS()->UO(UO_LVL_err) << "ERROR: MultiNodalSphericalJoint: SetPos1ToLocalCoord not implemented.";}
	// do not use this function!
	void SetPos2ToLocalCoord(int elem_nr, Vector3D loccoord){GetMBS()->UO(UO_LVL_err) << "ERROR: MultiNodalSphericalJoint: SetPos2ToLocalCoord not implemented.";}
	// define ground position
	void SetPos2ToGlobalCoord(Vector3D ground);

	virtual int GetNumberOfConstrainedCoords() const 
	{
			return 3;				// <- ??? is it called at all?
	}	
	virtual int SOS() const;
	virtual int NE_nodouble() const {  return 0; }		
	virtual int NKinPairs() const 
	{
		if(nodes2.Length())
		{
			return 2;
		}
		else
		{
			return 1; //ground joint
		}
	}
	Vector3D GetAverageDrawPosition(int kinpair);
	Vector3D GetDiscreteDrawPosition(int kinpair, int i);

	double GetDrawSizeCircleSpheres()		{ return draw_dim.Y();}
	int GetDrawSizeLineThicknessStar() 	{	return (int)draw_dim.Z();}

protected:
	//EDC double draw_dim(1);						//$EDC$[varaccess,EDCvarname="draw_size_center_sphere",EDCfolder="Graphics",tooltiptext="drawing dimensions of the sphere indicating the average position of the constraint. If set to -1, than global_draw_scalar_size is used."]
	//EDC double draw_dim(2);						//$EDC$[varaccess,EDCvarname="draw_size_circle_spheres",EDCfolder="Graphics",tooltiptext="drawing dimensions of the spheres indicating the positions of the multiple nodes."]
	//EDC double draw_dim(3);						//$EDC$[varaccess,EDCvarname="draw_size_line_thickness",EDCfolder="Graphics",tooltiptext="thickness of the lines from the outer spheres to the center sphere"]
	//EDC Vector3D loccoords(2);				//$EDC$[varaccess,EDCvarname="ground",EDCfolder="Position2",tooltiptext="global (average) position, just used if node_numbers.Length()==0"]
	TArray<int> nodes1;									//$EDC$[varaccess,EDCvarname="node_numbers",EDCfolder="Position1",tooltiptext="global node numbers of kinematic pair 1"]
	TArray<int>	nodes2;									//$EDC$[varaccess,EDCvarname="node_numbers",EDCfolder="Position2",tooltiptext="global node numbers of kinematic pair 2"]

	//// remove entries from EDC //$ DR 2012-08-22
	//EDC double draw_dim(1);						//$EDC$[varaccess,remove,EDCvarname="draw_size",EDCfolder="Graphics"]
	//EDC int nodes1(1);								//$EDC$[varaccess,remove,EDCvarname="node_number",EDCfolder="Position1"]
	//EDC int nodes1(2);								//$EDC$[varaccess,remove,EDCvarname="node_number",EDCfolder="Position2"]
	//EDC int elements(1);							//$EDC$[varaccess,remove,EDCvarname="element_number",EDCfolder="Position1"]
	//EDC int elements(2);							//$EDC$[varaccess,remove,EDCvarname="element_number",EDCfolder="Position2"]
	//EDC Vector3D loccoords(1);				//$EDC$[varaccess,remove,EDCvarname="position",EDCfolder="Position1"]
	//EDC Vector3D loccoords(2);				//$EDC$[varaccess,remove,EDCvarname="position",EDCfolder="Position2"]
	//EDC int use_penalty_formulation;	//$EDC$[varaccess,remove,EDCvarname="use_penalty_formulation",EDCfolder="Physics"]
	//EDC TArray<int> dir;							//$EDC$[varaccess,remove,EDCvarname="constrained_directions",EDCfolder="Physics.Lagrange"]

public:
	Vector3D GetAveragePosition(int kinpair) const;						// position of kinematic pair kinpair (average position)
	Vector3D GetAverageVelocity(int kinpair) const;						// velocity of kinematic pair kinpair (average velocity)
	Vector3D GetDiscretePosition(int kinpair, int i) const;		// position of node i of kinematic pair kinpair
	Vector3D GetDiscreteVelocity(int kinpair, int i) const;		// velocity of node i of kinematic pair kinpair
	//void SetPosToLocalNode(int i, int elem_nr, int loc_node_nr);
	//void SetPosToLocalCoord(int i, int elem_nr, Vector3D loccoord);
	void SetPosToGlobalNodes(int i,  TArray<int> glob_node_nrs);
};//$EDC$[endclass,MultiNodalSphericalJoint]





// ##############################################################################################
//$ DR 2012: FrictionConstraintUniDir added
class FrictionConstraintUniDir: public BasePointJoint//$EDC$[beginclass,classname=FrictionConstraintUniDir,parentclassname=BasePointJoint,texdescription="
//The FrictionConstraint constrains two elements at a local position or node each. Different friction coefficients for static and dynamic case are considered as well as stick/slip transition.
//"]
{
// DR: this will be the entries in far future: addelementtype=TAEconstraint,addelementtypename=FrictionConstraintUniDir
public:

	FrictionConstraintUniDir(MBS* mbsi):BasePointJoint(mbsi)
	{
	};

	virtual Element* GetCopy()
	{
		Element* ec = new FrictionConstraintUniDir(*this);
		return ec;
	}

	virtual void CopyFrom(const Element& e)
	{
		BasePointJoint::CopyFrom(e);
		const FrictionConstraintUniDir& ce = (const FrictionConstraintUniDir&)e;
		Fn0 = ce.Fn0;
		ax_vec = ce.ax_vec;
		constantNormalF = ce.constantNormalF;

		//$ 2012-03-21 Test
		keep_sliding = ce.keep_sliding;
	}

	virtual int GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //for special value sensor
	virtual int ReadSingleElementData(ReadWriteElementDataVariableType& RWdata); 		//for special value sensor

	virtual void GetElementDataAuto(ElementDataContainer& edc);   
	//virtual int SetElementDataAuto(const ElementDataContainer& edc);
	virtual int SetElementDataAuto(ElementDataContainer& edc);
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual const char* GetElementSpec() const {return "FrictionConstraintUniDir";}
	virtual void SetDefaultValues();

	virtual void SetFrictionConstraintUniDir_LocalNode_to_LocalNode(int elem1, int node1, int elem2, int node2, Vector3D normal_vec, double stiffness, double damping, double friction_coeff_st,  double friction_coeff_kin, Vector3D axial_dir = Vector3D(0.,0.,1.));
	virtual void SetFrictionConstraintUniDir_LocalPos_to_LocalPos(int elem1, Vector3D lc1, int elem2, Vector3D lc2, Vector3D normal_vec, double stiffness, double damping, double friction_coeff_st,  double friction_coeff_kin, Vector3D axial_dir = Vector3D(0.,0.,1.));

	//1=stick, 2/3=pos_ax/pos_tan, 4/5=vel_ax/vel_tan 6/7= OLD vel_ax/vel_tan
	virtual int DataS() const {return 7;}
	
	virtual void PostprocessingStep();
	virtual double PostNewtonStep(double t);
	virtual Vector3D ComputeForce(double t) const;

	// Get force that acts during stick-phase
	virtual Vector3D GetAxialForce(Vector3D u, double stiffness) const;
	virtual double GetNormalForce(Vector3D u,double stiffness)const;
	virtual Vector3D GetAxialDir() const {	return ax_vec;}
	virtual void SetAxialDir(Vector3D v)
	{
		v.Normalize();
		ax_vec = v;
	}

	// Friction Force when body is sliding: F = µk*Fn
	virtual Vector3D GetSlidingFrictionForce(double coef_kin)const;
	// maximal Force at which sliding starts: Fmax = µs*Fn
	virtual double GetTraction(double t, Vector3D u) const;
	
	virtual void SetNormalForce(double Fnormal) const {	Fn0 = Fnormal;}
	virtual void SetUseConstantNormalForce(bool flag_constant){	constantNormalF = flag_constant;}

	virtual void DrawElement() ;
	
	// get actor force, with optional direction TSX, TSY, TSZ
	// if neither TSX, TSY nor TSZ is used, than the flag stick/slip is returned
	virtual double GetActorForce(double computation_time, int dir) const
	{
		if (dir == 0) return IsSticking();
		if (dir > 3) assert(0 && "FrictionConstraintUniDir::GetActorForce");

		return ComputeForce(computation_time)(dir);
	}

	// get force projected in special direction
	virtual double GetSpecialSensorValue(int nr, double time) const;

	virtual double GetFrictionCoeff_st() const {	return fr_coeff_st;	}
	virtual double GetFrictionCoeff_kin() const{	return fr_coeff_kin;}
	virtual void SetFrictionCoeff_st(double coeff) { fr_coeff_st = coeff; }
	virtual void SetFrictionCoeff_kin(double coeff){ fr_coeff_kin = coeff; }
	virtual void SetVelocityTolerance(double tol) {	velocity_tolerance = tol;	}
	virtual double GetVelocityTolerance() const 	{	return velocity_tolerance; }

	//XData(): 1=stick, 2/3=pos_ax/pos_tan, 4/5=vel_ax/vel_tan
	virtual bool IsSticking() const		{	return XData(1) != 0;}
	virtual void IsSticking(bool flag){ XData(1) = flag;}
	virtual double StickingPosAx() const		{	return XData(2);}
	virtual void StickingPosAx(double pos){ XData(2) = pos;}
	virtual double StickingPosTan() const		{	return XData(3);}
	virtual void StickingPosTan(double pos){ XData(3) = pos;}
	virtual double SlidingVelAx() const		{	return XData(4);}
	virtual void SlidingVelAx(double vel){ XData(4) = vel;}
	virtual double SlidingVelTan() const		{	return XData(5);}
	virtual void SlidingVelTan(double vel){ XData(5) = vel;}

	virtual double SlidingVelAxOLD() const		{	return XData(6);}
	virtual void SlidingVelAxOLD(double vel){ XData(6) = vel;}
	virtual double SlidingVelTanOLD() const		{	return XData(7);}
	virtual void SlidingVelTanOLD(double vel){ XData(7) = vel;}

	virtual void SetKeepSliding(int flag) {keep_sliding=flag;}

protected:
	mutable double Fn0;					//$EDC$[varaccess,EDCvarname="Fn",EDCfolder="Physics",tooltiptext="constant normal force Fn"]
	Vector3D ax_vec;						// vector that defines axis of rotor or vector in sliding surface

	int constantNormalF;				//$EDC$[varaccess,EDCvarname="use_constant_Fn",EDCfolder="Physics",int_bool,tooltiptext="if this flag is set, normal force is assumed to be constant"]
	double velocity_tolerance;	//$EDC$[varaccess,EDCvarname="velocity_tolerance",EDCfolder="Physics",tooltiptext="if velocity is below this value, than velocity is assumed to be zero"]
	double fr_coeff_st;					//$EDC$[varaccess,EDCvarname="fr_coeff_st",EDCfolder="Physics",tooltiptext="friction coefficient µs, used when elements are sticking (F = µs*Fn)"]
	double fr_coeff_kin;				//$EDC$[varaccess,EDCvarname="fr_coeff_kin",EDCfolder="Physics",tooltiptext="friction coefficient µk, used when elements are sliding (F = µk*Fn)"]
	int keep_sliding;
	// all other informations are stored in XData --> see DataS() for documentation
};//$EDC$[endclass,FrictionConstraintUniDir]


/*
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// MBSController
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// The controller is treated as an element in order to be able to include differential equations
// However, for simple purposes the controller might just compute a force
// The controller input is based on sensor values!
class MBSController: public Constraint
{
public:
	//MBSController(MBS* mbsi):Constraint(mbsi), offsets(), factors(), options() 
	//{	};

	MBSController(MBS* mbsi):Constraint(mbsi), offsets(), factors(), options()
	{	
		mbs = mbsi; //is set in ::Contraint::Element(mbsi)
		nvariables = 0;
		x_init = Vector(SS());
		mathfunc.SetConstant(1);
		type = TMBSElement(TConstraint + TController);
		elementname = GetElementSpec();
		elements.SetLen(0);
	};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new MBSController(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const MBSController& ce = (const MBSController&)e;

		offsets = ce.offsets;
		factors = ce.factors;
		options = ce.options;

		nvariables = ce.nvariables;
		mathfunc = ce.mathfunc;
	}

	virtual const char* GetElementSpec() const {return "Controller";}
	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(const ElementDataContainer& edc); //set element data according to ElementDataContainer

	//evaluate constraints: len(z)
	virtual void EvalF(Vector& f, double t);
	virtual void EvalG(Vector& f, double t) {};

	//To be replaced in derived class
	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f) {};

	virtual int IS() const {return 0;};
	virtual int ES() const {return nvariables;};
	virtual int Dim() const {return 3;}  //is always 3D

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const 
	{
		return Vector3D(0.,0.,0.);
	}

	virtual void DrawElement() 
	{
	};

	virtual int AddSensor(int sensor_num, double offset, double factor, int option);

	virtual void FlushSensors()
	{
		sensors.SetLen(0);
		offsets.SetLen(0);
		factors.SetLen(0);
		options.SetLen(0);
		nvariables = 0;
	}
	virtual double GetSensorValue(int i) {return GetSensor(i).GetValues()(1);}

	virtual double ComputeControlValue(); //evaluate controller

	virtual void SetInitialValues(const Vector& xi) {x_init = xi;}

	virtual void SetMathFunction(const MathFunction& mf) {mathfunc = mf;}

protected:
	//TArray<int> sensors;			//sensor number i ==> now in element
	TArray<double> offsets;		//offset for sensor number i
	TArray<double> factors;		//factor for sensor number i
	TArray<int> options;			//option == 1: integrate sensor number i

	int nvariables;						//number of integrators or other operations that require variables to integrate

	MathFunction mathfunc;		//generally manipulating the output value, e.g. by linear interpolated table ...
};
*/

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Rope3D (DR)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// elastic rope that is always under tension and can be fixed to multiple bodies and ground
class Rope3D: public BasePointJoint//$EDC$[beginclass,addelementtypename=Rope3D,classname=Rope3D,parentclassname=BasePointJoint,addelementtype=TAEconstraint,
//texdescription="Elastic rope that is always under tension and can be fixed to multiple bodies and ground. There are 2 different kinds of suspensions points. Suspension points fixed on the ground are defined with the element number 0 and the global position. Suspension points on bodies are defined with the element number and the corresponding local position.",
//texdescriptionLimitations="The rope is assumed to be straight between 2 suspension points. No negative forces can be transmitted by a rope. The computation of the time derivative of the length of the rope is just an approximation. Therefore the damping of the rope may be represented slightly incorrect.",
//example="Rope3D.txt",
//figure="Rope3D,Point mass with rope"]
{
public:
	Rope3D(MBS* mbsi):BasePointJoint(mbsi)
	{	
		ElementDefaultConstructorInitialization();
	};

	Rope3D(const Rope3D& ct): BasePointJoint(ct.mbs)
	{
		ElementDefaultConstructorInitialization();
		CopyFrom(ct);
	};

	~Rope3D()
	{
	};

	virtual void ElementDefaultConstructorInitialization();
	virtual void Initialize();

	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	virtual void GetNecessaryKinematicAccessFunctions(TArray<int> &KAF, int numberOfKinematicPair);

	virtual Element* GetCopy()
	{
		Element* ec = new Rope3D(*this);
		return ec;
	}

	virtual void CopyFrom(const Element& e);

	virtual const char* GetElementSpec() const {return "Rope3D";}

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer
  virtual int GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //for special value sensor
	virtual int ReadSingleElementData(ReadWriteElementDataVariableType& RWdata); 		//for special value sensor
	virtual int WriteSingleElementData(const ReadWriteElementDataVariableType& RWdata); // for IOElementDataModifier

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		
	virtual int SetElementDataAuto(ElementDataContainer& edc);
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual void DrawElement();

	virtual double GetLengthOfRope() const;
	//virtual double GetVelocityOfRope() const;		//$ DR 2013-08-29 bugfix for sensor of actor force

	virtual Vector3D GetPosition(int i) const;
	virtual Vector3D GetDrawPosition(int i) const;
	virtual double ComputeRopeForce(double t) const;
	virtual Vector3D ComputeForceDirection(int i) const;
	virtual void EvalF2(Vector& f, double t);

	// compute the velocity of the rope out of the change of the length of the rope
	virtual int DataS() const {return 3;}
	virtual double OldLength() const		{	return XData(1);}
	virtual void OldLength(double l)		{ XData(1) = l;}
	virtual double OldTime() const			{	return XData(2);}
	virtual void OldTime(double t)			{ XData(2) = t;}
	virtual double GetVelocityOfRope() const			{	return XData(3);} //$ DR 2013-08-29 bugfix for sensor of actor force
	virtual void SetVelocityOfRope(double v)			{ XData(3) = v;}		//$ DR 2013-08-29 bugfix for sensor of actor force
	virtual void PostprocessingStep();

	virtual void LinkToElementsPenalty();

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const { return GetPosition(1);}

	// do not use these functions!
	void SetPos1ToLocalNode(int elem_nr, int loc_node_nr)		{GetMBS()->UO(UO_LVL_err) << "ERROR: Rope3D: SetPos1ToLocalNode not implemented.";}
	void SetPos2ToLocalNode(int elem_nr, int loc_node_nr)		{GetMBS()->UO(UO_LVL_err) << "ERROR: Rope3D: SetPos2ToLocalNode not implemented.";}
	void SetPos1ToLocalCoord(int elem_nr, Vector3D loccoord){GetMBS()->UO(UO_LVL_err) << "ERROR: Rope3D: SetPos1ToLocalCoord not implemented.";}
	void SetPos2ToLocalCoord(int elem_nr, Vector3D loccoord){GetMBS()->UO(UO_LVL_err) << "ERROR: Rope3D: SetPos2ToLocalCoord not implemented.";}
	void SetPos2ToGlobalCoord(Vector3D ground) {GetMBS()->UO(UO_LVL_err) << "ERROR: Rope3D: SetPos2ToGlobalCoord not implemented.";}
	void SetPositions(TArray<int> element_numbers,  TArray<Vector3D> coordinates);

	virtual int SOS() const;
	virtual int NE_nodouble() const {  return NKinPairs(); }		
	virtual int NKinPairs() const {return elements.Length();}

protected:
	double initial_length;	//$EDC$[varaccess, readonly, EDCvarname="rope_length", EDCfolder="Geometry", tooltiptext="initial length l0 of rope (computed automatically)"]
	double coiled_length;

	// remove entries
	//EDC double draw_local_frame_size;					//$EDC$[varaccess,remove,EDCvarname="draw_size_joint_local_frame",EDCfolder="Graphics"]
	//EDC Vector3D spring_stiffness3;						//$EDC$[varaccess,remove,EDCvarname="spring_stiffness_vector",EDCfolder="Physics.Penalty"]	
	//EDC	int use_local_coordinate_system;			//$EDC$[varaccess,remove,EDCvarname="use_local_coordinate_system",EDCfolder="Geometry"]
	//EDC	int stiffness_in_joint_local_frame;		//$EDC$[varaccess,remove,EDCvarname="use_joint_local_frame",EDCfolder="Geometry"]
	//EDC	Matrix3D JA0i;												//$EDC$[varaccess,remove,EDCvarname="joint_local_frame",EDCfolder="Geometry"]
  //EDC int use_penalty_formulation;					//$EDC$[varaccess,remove,EDCvarname="use_penalty_formulation",EDCfolder="Physics"]
	//EDC	TArray<int> dirRot;										//$EDC$[varaccess,remove,EDCvarname="constrained_directions",EDCfolder="Physics.Lagrange"]
	//EDC int elements(1);											//$EDC$[varaccess,remove,EDCvarname="element_number",EDCfolder="Position1"]
	//EDC int elements(2);											//$EDC$[varaccess,remove,EDCvarname="element_number",EDCfolder="Position2"]
	//EDC Vector3D loccoords(1);								//$EDC$[varaccess,remove,EDCvarname="position",EDCfolder="Position1"]
	//EDC Vector3D loccoords(2);								//$EDC$[varaccess,remove,EDCvarname="position",EDCfolder="Position2"]
	//EDC int nodes(1);													//$EDC$[varaccess,remove,EDCvarname="node_number",EDCfolder="Position1"]
	//EDC int nodes(2);													//$EDC$[varaccess,remove,EDCvarname="node_number",EDCfolder="Position2"]
	//EDC	double spring_stiffness;							//$EDC$[varaccess,remove,EDCvarname="spring_stiffness",EDCfolder="Physics.Penalty"]

	// new access
	//EDC TArray<int> elements;									//$EDC$[varaccess,EDCvarname="element_numbers",EDCfolder="Geometry",variable_length_vector,tooltiptext="element numbers of the suspension points"]

};//$EDC$[endclass,Rope3D]


class FrictionConstraint: public Constraint //$EDC$[beginclass,classname=FrictionConstraint,parentclassname=Constraint,addelementtype=TAEconstraint+TAEspecial_connector,addelementtypename=FrictionConstraint,texdescription="
//The FrictionConstraint is acting on an arbitraty coordinate, including rotations. 
//It can be used to connect two elements to each other or one element to ground. 
//Up to a specified threshold of the force, the constraint is sticking, which is realized by a spring-damper formulation. 
//Above this threshold, a constant friction force is applied during the sliding phase. 
//Alternatively sticking can be switched off and a coulomb friction force, with a transition region for very small velocities, can be applied.",
//texdescriptionEquations="
//\begin{equation}
//	F_{st} = c x + d v
//\end{equation} \begin{equation}
//	F_{sl} = \mu_{kin} F_n
//\end{equation}",
//texdescriptionComments="The switching from sticking phase to sliding phase is done automatically, as soon as $F_{st} > \mu_{st} F_n$. The switching to sticking phase is performed when the absolute value of the velocity $v$ is smaller than the specified velocity\_tolerance.\\
//If the solver does not converge close to the switching points, set the solver option SolverOptions.Discontinuous.ignore\_max\_iterations = 1.\\
//If you are using the FrictionConstraint in order to constrain rotations, problems may occur when the change of the angle is discontinous, e.g. if it exceeds pi/2.",
//modus="{sticking}{During sticking phase, the constraint is implemented as spring-damper, with the force $F_{st}$, the spring stiffness $c$ and the damping coefficient $d$.}",
//modus="{sliding}{During sliding phase, a constant friction force $F_{sl}$ is applied. $F_{sl}$ depends on the normal force $F_n$.
//If the flag keep\_sliding is active, than a transition region for small velocities is used.}",
//example="FrictionConstraint.txt",
//figure="FrictionConstraint,FrictionConstraint with friction forces $F_{st}$ and $F_{sl}$, with sticking (left figure, keep\_sliding = 0) and without sticking (right figure, keep\_sliding = 1)."
//]
{
public:
	FrictionConstraint(MBS* mbsi):Constraint(mbsi), loccoords() 
	{	ElementDefaultConstructorInitialization();};

	FrictionConstraint(const FrictionConstraint& ct): Constraint(ct.mbs), loccoords(0)
	{
		CopyFrom(ct);
	};

	// set functions
	virtual void SetFrictionConstraint(int en1, int lc1) 
	{
		ElementDefaultConstructorInitialization();
		elements(1) = en1;
		loccoords(1) = lc1;
	};
	virtual void SetFrictionConstraint(int en1, int en2, int lc1, int lc2)
	{	
		ElementDefaultConstructorInitialization();

		elements(1) = en1;
		loccoords(1) = lc1;
		elements(2) = en2;
		loccoords(2) = lc2;
	}

	//virtual void Initialize();

	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new FrictionConstraint(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const FrictionConstraint& ce = (const FrictionConstraint&)e;
		loccoords = ce.loccoords;
		Fn0=ce.Fn0;				
		velocity_tolerance=ce.velocity_tolerance;
		fr_coeff_st=ce.fr_coeff_st;				
		fr_coeff_kin=ce.fr_coeff_kin;			
		damping_coeff=ce.damping_coeff;			
		keep_sliding=ce.keep_sliding;
	}
	
	virtual void ElementDefaultConstructorInitialization();

	virtual void Initialize() 
	{
		if(keep_sliding)
		{
			IsSticking(0);	// start with sliding
		}
	};


	virtual double GetPos1() const {return GetElem(1).XG(loccoords(1));}
	virtual double GetPos2() const 
	{
		if (elements(2)) return GetElem(2).XG(loccoords(2));
		else return GetElem(1).GetXInit()(loccoords(1));		// initial value = ground
	}
	virtual double GetVel1() const {return GetElem(1).XGP(loccoords(1));}
	virtual double GetVel2() const
	{
		if (elements(2)) return GetElem(2).XGP(loccoords(2));
		else return 0; // initial value = ground
	}
	virtual int NE() const 
	{
		if(elements(2)) return 2;
		else if(elements(1)) return 1;
		else return 0; // when adding the empty constraint
	}		
	
	virtual const char* GetElementSpec() const {return "FrictionConstraint";}
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data according to ElementDataContainer
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter


	//1=stick, 2=pos, 3=vel ,4= OLD vel
	virtual int DataS() const {return 4;}
	virtual void PostprocessingStep();
	virtual double PostNewtonStep(double t);
	virtual void EvalF2(Vector& f, double t); //second order equations for constraints: M \ddot u = F2(u,\dot u,t), len(u)

	virtual double GetFrictionCoeff_st() const {	return fr_coeff_st;	}
	virtual double GetFrictionCoeff_kin() const{	return fr_coeff_kin;}
	virtual void SetFrictionCoeff_st(double coeff) { fr_coeff_st = coeff; }
	virtual void SetFrictionCoeff_kin(double coeff){ fr_coeff_kin = coeff; }
	virtual void SetVelocityTolerance(double tol) {	velocity_tolerance = tol;	}
	virtual double GetVelocityTolerance() const 	{	return velocity_tolerance; }

	//XData(): 1=stick, 2=pos, 3=vel ,4= OLD vel
	virtual bool IsSticking() const		{	return XData(1);}
	virtual void IsSticking(bool flag){ XData(1) = flag;}
	virtual double StickingPos() const		{	return XData(2);}
	virtual void StickingPos(double pos){ XData(2) = pos;}
	virtual double SlidingVel() const		{	return XData(3);}
	virtual void SlidingVel(double vel){ XData(3) = vel;}
	virtual double SlidingVelOLD() const		{	return XData(4);}
	virtual void SlidingVelOLD(double vel){ XData(4) = vel;}

	virtual void SetKeepSliding(int flag) {keep_sliding=flag;}

	virtual double GetAxialForce(double u, double stiffness) const;
	virtual double GetNormalForce(double u,double stiffness)const;
	virtual double GetSlidingFrictionForce(double coef_kin)const;
	virtual double GetTraction(double t, double u) const;

	virtual int IS() const {return 0;}; 
	virtual int Dim() const 
	{
		if(elements(1)) {return GetElem(1).Dim();}
		else return 3;
	}

	//do not add initial conditions for 2nd order size DOF, because they are set in constrained elements
	virtual void SetGlobalInitConditions(Vector& x_glob) {}

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const 
	{
		return GetElem(1).GetRefPosD();
	}
	
	virtual void DrawElement();

	//compute constraint force, for penalty or Lagrange approach
	virtual double ComputeForce(double t) const;

	virtual int WriteSingleElementData(const ReadWriteElementDataVariableType& RWdata);
	virtual int ReadSingleElementData(ReadWriteElementDataVariableType& RWdata);
	virtual int GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables);


	//get a component of the constraint force ==> for sensors
	//deprecated function for old sensors
	virtual double GetActorForce(double computation_time, int dir=0) const {return ComputeForce(computation_time);} //get actor force, with optional direction (x,y,z, etc.)

	virtual double GetDampingCoeff() const {return damping_coeff;}
	virtual void SetDampingCoeff(double damp) {	damping_coeff = damp;	}

protected:
	TArray<int> loccoords;
	mutable double Fn0;					//$EDC$[varaccess,EDCvarname="normal_force",EDCfolder="Physics",tooltiptext="constant normal force Fn"]

	double velocity_tolerance;	//$EDC$[varaccess,EDCvarname="velocity_tolerance",EDCfolder="Physics",minval=0,tooltiptext="If velocity is below this value, sticking starts, or if 'keep sliding' is active, the transition region is used."]
	double fr_coeff_st;					//$EDC$[varaccess,EDCvarname="fr_coeff_st",EDCfolder="Physics",tooltiptext="static friction coefficient, used to determine the threshold when sliding starts."]
	double fr_coeff_kin;				//$EDC$[varaccess,EDCvarname="fr_coeff_kin",EDCfolder="Physics",tooltiptext="kinematic friction coefficient, used to calculate the constant force during sliding phase."]
	int keep_sliding;						//$EDC$[varaccess,EDCvarname="keep_sliding",int_bool,EDCfolder="Physics",tooltiptext="The constraint will never go to modus 'stick'."]

	// all other informations are stored in XData --> see DataS() for documentation

	//EDC	int use_local_coordinate_system;			//$EDC$[varaccess,remove,EDCvarname="use_local_coordinate_system",EDCfolder="Geometry"]
  //EDC int use_penalty_formulation;					//$EDC$[varaccess,remove,EDCvarname="use_penalty_formulation",EDCfolder="Physics"]
	//EDC   Vector3D col;                       //$EDC$[varaccess,remove,EDCvarname="RGB_color",EDCfolder="Graphics"]
	//EDC double draw_dim(1);										//$EDC$[varaccess,EDCvarname="draw_size",EDCfolder="Graphics",tooltiptext="Drawing dimensions of constraint. If set to -1, than global_draw_scalar_size is used."]

	//EDC	double spring_stiffness;							//$EDC$[varaccess,remove,EDCvarname="spring_stiffness",EDCfolder="Physics.Penalty"]
	//EDC	double spring_stiffness;							//$EDC$[varaccess,EDCvarname="spring_stiffness",EDCfolder="Physics.Penalty",tooltiptext="spring stiffness c, only used during sticking phase!"]
	double damping_coeff;				//$EDC$[varaccess,EDCvarname="damping",minval=0,EDCfolder="Physics.Penalty",tooltiptext="damping coefficient d for viscous damping, only used during sticking phase!"]

};//$EDC$[endclass,FrictionConstraint]



//----------------------------------------------------------------------------------------------------------------------------------------
//$ DR 2013-09-30 added Contact1D
class Contact1D: public Constraint //$EDC$[beginclass,classname=Contact1D,parentclassname=Constraint,addelementtype=TAEconstraint+TAEspecial_connector,addelementtypename=Contact1D,texdescription="
//Contact1D realizes a contact formulation between two elements or one element and ground. Only one coordinate (direction) is considered per element.",
//figure="Contact1DGround, Description of the geometry options in the case of a ground constraint.",
//figure="Contact1D, Description of the geometry options in the case of 2 elements.",
//texdescriptionGeometry="Figure \ref{Contact1Dfigure1} shows the meaning of the values local coordinate and position in the case of a ground constraint. The only direction which is considered is that defined by Coordinate1.local coordinate. Figure \ref{Contact1Dfigure1} shows the case for 2 elements.\\
//ATTENTION: Be carefull when using coordinates which do not represent a position!",
//texdescriptionEquations="Some general definitions:\begin{equation}
//pos = coordinate + local position \end{equation} \begin{equation}
//u = pos_1 - pos_2 \end{equation} \begin{equation}
//v = vel_1 - vel_2 \end{equation}
//\textbf{Mode 1:}\\
//if $u \geq 0$: \begin{equation}
//  F = 0 \end{equation}
//else: \begin{equation}
//  F = c u + d v \end{equation}",
//example="Contact1D.txt",
//modus="{Mode 1}{Penalty Formulation with spring and damper. The bodies will penetrate slightly according to the spring stiffness. Results may depend on chosen step size!}",
//modus="{Mode 2}{Lagrange Formulation (not implemented yet)}"]
{
public:
	Contact1D(MBS* mbsi):Constraint(mbsi)  
	{	ElementDefaultConstructorInitialization();};

	Contact1D(const Contact1D& ct): Constraint(ct.mbs)
	{
		CopyFrom(ct);
	};

	// set functions - for element-based formulations
	virtual void SetContact1DElement(int en1, int lc1);
	virtual void SetContact1DElement(int en1, int en2, int lc1, int lc2);
	// set functions - for node-based formulations
	virtual void SetContact1DNode(int nn1, int lc1);
	virtual void SetContact1DNode(int nn1, int nn2, int lc1, int lc2);

	void SetLPos1(double lpos1) { this->lpos1 = lpos1; }
	void SetDirection(double direction) { this->contact_direction = direction; }

	virtual void Initialize() 
	{

	};

	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new Contact1D(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e);
	
	virtual void ElementDefaultConstructorInitialization();

	// the base-class implementation needs to be redefined in case of node-based contact
	virtual void LinkToElements();
	virtual int SOS() const;

	virtual double GetPos1() const;
	virtual double GetPos2() const;
	virtual double GetVel1() const;
	virtual double GetVel2() const;
	virtual int NE() const;
	virtual const Node& GetNode(int i) const
	{
		return GetMBS()->GetNode(nodes(i));
	}
	
	virtual const char* GetElementSpec() const {return "Contact1D";}
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data according to ElementDataContainer
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual int DataS() const {return 0;}
	virtual void PostprocessingStep();
	virtual double PostNewtonStep(double t);
	virtual void EvalF2(Vector& f, double t); //second order equations for constraints: M \ddot u = F2(u,\dot u,t), len(u)
	virtual void EvalG(Vector& f, double t);
	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	virtual int IS() const 
	{
		if(mode==2) {return 1;}
		else				{return 0;}
	}; 
	virtual int Dim() const 
	{
		if(elements(1)) {return GetElem(1).Dim();}
		return GetNode(1).Dim();
	}

	//do not add initial conditions for 2nd order size DOF, because they are set in constrained elements
	virtual void SetGlobalInitConditions(Vector& x_glob) {}

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const 
	{
		if(elements(1))
			return GetElem(1).GetRefPosD();
		return GetNode(1).GetPosD();
	}
	
	virtual void DrawElement();

	//compute constraint force, for penalty or Lagrange approach
	virtual double ComputeForce(double t) const;

	virtual int WriteSingleElementData(const ReadWriteElementDataVariableType& RWdata);
	virtual int ReadSingleElementData(ReadWriteElementDataVariableType& RWdata);
	virtual int GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables);

	virtual double GetDampingCoeff() const {return damping_coeff;}
	virtual void SetDampingCoeff(double damp) {	damping_coeff = damp;	}

protected:
	int lcoord1;																//$EDC$[varaccess,EDCvarname="local_coordinate",EDCfolder="Coordinate1",tooltiptext="Local coordinate of element 1 to be constrained"]
	int lcoord2;																//$EDC$[varaccess,EDCvarname="local_coordinate",EDCfolder="Coordinate2",tooltiptext="Local coordinate of element 2 to be constrained (not used if ground constraint)"]

	double lpos1;																//$EDC$[varaccess,EDCvarname="position",EDCfolder="Coordinate1",tooltiptext="Local position at which contact occurs"]
	double lpos2;																//$EDC$[varaccess,EDCvarname="position",EDCfolder="Coordinate2",tooltiptext="Local (or global if ground) position at which contact occurs"]
	double contact_direction;													//$EDC$[varaccess,EDCvarname="direction",EDCfolder="Physics",tooltiptext="Direction of the contact: +1 if the first body is on top, or else -1"]

	TArray<int> nodes;

	int mode;																		//$EDC$[varaccess,EDCvarname="mode",EDCfolder="Physics",tooltiptext="mode of computation"]

	//EDC	int use_local_coordinate_system;			//$EDC$[varaccess,remove,EDCvarname="use_local_coordinate_system",EDCfolder="Geometry"]
  //EDC int use_penalty_formulation;					//$EDC$[varaccess,remove,EDCvarname="use_penalty_formulation",EDCfolder="Physics"]
	//EDC   Vector3D col;                       //$EDC$[varaccess,remove,EDCvarname="RGB_color",EDCfolder="Graphics"]
	//EDC double draw_dim(1);										//$EDC$[varaccess,EDCvarname="draw_size",EDCfolder="Graphics",tooltiptext="Drawing dimensions of constraint. If set to -1, than global_draw_scalar_size is used."]

	//EDC	double spring_stiffness;							//$EDC$[varaccess,remove,EDCvarname="spring_stiffness",EDCfolder="Physics.Penalty"]
	//EDC	double spring_stiffness;							//$EDC$[varaccess,EDCvarname="spring_stiffness",EDCfolder="Physics.Mode1",tooltiptext="spring stiffness c"]
	double damping_coeff;				//$EDC$[varaccess,EDCvarname="damping",minval=0,EDCfolder="Physics.Mode1",tooltiptext="damping coefficient d for viscous damping"]

};//$EDC$[endclass,Contact1D]




#endif
