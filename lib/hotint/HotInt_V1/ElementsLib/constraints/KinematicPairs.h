//#**************************************************************
//#
//# filename:             KinematicPairs.h
//#
//# author:               Gerstmayr Johannes, Reischl Daniel
//#
//# generated:						20.April 2011
//# description:          Kinematic pairs can be divided in lower (surface contact) and higher (point or line contact) kinematic pairs. 
//#												In this file the lower kinematic pairs and the rigid joints (all d.o.f. are constrained) are provided.
//#												Other implemented constraints are in SpecialConstraints.h and some old kinematic pairs can be found in
//#												KinematicPairsDeprecated.h.
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
 
//$ SW 2013-08-29: Moved CylindricalJoint, PrismaticJoint, RevoluteJoint and RigidJoint to KinematicPairsDeprecated.h
// and renamed them into CylindricalJointOLD, PrismaticJointOLD, RevoluteJointOLD and RigidJointOLD respectively.
// The classes CylindricalJoint_, PrismaticJoint_, RevoluteJoint_ and RigidJoint_ have been renamed to CylindricalJoint,
// PrismaticJoint, RevoluteJoint and RigidJoint (these classes can be found in RigidBodyJoints.h)

#ifndef KINEMATICPAIRS__H
#define KINEMATICPAIRS__H

#include "constraint.h"
#include "node.h"

//# Rigid:				all d.o.f. are constrained, e.g. welding

//# lower kinematic pairs:

//# Spherical:		all translatory d.o.f. are constrained
//# Rotary:				all translatory and 2 rotatory d.o.f. are constrained, only rotation about 1 axis possible
//# Cylindrical:	like Rotary, but additionally translation along rotational-axis possible
//# Translatory:	all rotary and 2 translatory d.o.f are constrained, e.g. prismatic joint
//# Plane:				not implemented
//# Skrew-type:		not implemented




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Spherical	Spherical	Spherical	Spherical	Spherical	Spherical	Spherical	
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Spherical: NodalConstraint	
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: NodalConstraint (Constraint)
// short description: constrains global nodes (global node number, no element attached)
// available formulations: Lagrange multiplier (only ground constraint), Penalty Formulation
// type of constraint equation: 
// index formulation available: 
// ground constraint: yes
// element-element constraint: yes
// development status: 
// long description:	all translatory d.o.f. are constrained
//								
// class variables:
//  
//**end(ued)**
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class NodalConstraint: public Constraint
{
public:
	//$ DR 2011-04-20:[ old code
	//NodalConstraint(MBS* mbsi):Constraint(mbsi), penalty(0), loccoords(), ground(), velocity_constraint(0)
	//{	
	//	spring_stiffness = 0.0;
	//	damping_coeff = 0.0;
	//};

	//NodalConstraint(const NodalConstraint& ct): Constraint(ct.mbs), loccoords(), ground(), velocity_constraint(0)
	//{
	//	CopyFrom(ct);
	//};
	//$ DR 2011-04-20:] old code

	NodalConstraint(MBS* mbsi):Constraint(mbsi)
	{	
		InitNC();
	};

	NodalConstraint(const NodalConstraint& ct): Constraint(ct.mbs)
	{
		InitNC();
		CopyFrom(ct);
	};

	~NodalConstraint()
	{
		for(int i=1; i<=ground.Length(); i++)
		{
			delete ground(i);
		}
	};

	virtual void InitNC()
	{
		//penalty = 0;
		SetPenaltyFormulation(0);
		elements.SetLen(0);
		loccoords.SetLen(0);
		ground.SetLen(0);
		nodes.SetLen(0);
	}

	// set functions
	void SetNodalConstraint(int nn1, int lc1, double ground1, double r=-1, int velconstraint=0)
	{
		SetNodalConstraintSpringDamper(nn1, lc1, ground1, r, 0, 0, velconstraint);

		//penalty = 0;							//$ DR 2011-04-20: old code, now in SetNodalConstraintSpringDamper if stiffness = 0

		//InitNC();

		//node = nn1;
		//MathFunction* mathfun= new MathFunction();
		//mathfun->SetConstant(ground1);
		//AddCoordGround(lc1, mathfun);
		//draw_dim.X() = r;

		//x_init = Vector(SS());
		//elementname = GetElementSpec();
		//velocity_constraint = velconstraint;
		//if (!velconstraint)
		//	x_init(1) = ground1;
		//else
		//	x_init(SOS()+1) = ground1;
	};

	void SetNodalConstraint(int nn1, TArray<int>& lc, TArray<double>& aground, double r=-1, int velconstraint=0)
	{
		SetNodalConstraintSpringDamper(nn1, lc, aground, r, 0, 0, velconstraint);

		//penalty = 0;							//$ DR 2011-04-20: old code, now in SetNodalConstraintSpringDamper if stiffness = 0

		//InitNC();

		////ground = aground;
		//ground.SetLen(aground.Length());
		//for(int i=1; i<=ground.Length(); i++)
		//{
		//	MathFunction* mathfun= new MathFunction;
		//	mathfun->SetConstant(aground(i));
		//	ground(i)= mathfun;
		//}
		//node = nn1;
		//draw_dim.X() = r;

		//x_init = Vector(SS());
		//elementname = GetElementSpec();
		//velocity_constraint = velconstraint;
		//if (!velconstraint)
		//	for(int i=1; i<=ground.Length(); i++)
		//		x_init(i) = ground(i)->Evaluate(0);
		//else
		//	for(int i=1; i<=ground.Length(); i++)
		//		x_init(SOS()+i) = ground(i)->Evaluate(0);
	};

	void SetNodalConstraint(int nn1, TArray<int>& lc, TArray<MathFunction*>& aground, double r=-1, int velconstraint=0)
	{
		InitNC();

		loccoords = lc;
		//ground = aground;
		ground.SetLen(aground.Length());
		for(int i=1; i<=ground.Length(); i++)
		{
			MathFunction* mathfun= new MathFunction;
			mathfun=(aground(i));
			ground(i)= mathfun;
			//delete aground(i);
		}
		nodes.Add(nn1);
		SetDrawSizeScalar(r);		//draw_dim.X() = r;

		x_init = Vector(SS());
		elementname = GetElementSpec();
		velocity_constraint = velconstraint;
		if (!velconstraint)
			for(int i=1; i<=ground.Length(); i++)
				x_init(i) = ground(i)->Evaluate(0);
		else
			for(int i=1; i<=ground.Length(); i++)
				x_init(SOS()+i) = ground(i)->Evaluate(0);
	};

	void SetNodalConstraintSpringDamper(int nn1, int lc1, double ground1, double r, double spring_stiffnessi, double damping_coeffi = 0, int velconstraint=0)
	{
		TArray<int> lcs; lcs.Add(lc1);
		TArray<double> grounds; grounds.Add(ground1);

		SetNodalConstraintSpringDamper(nn1, lcs, grounds, r, spring_stiffnessi, damping_coeffi, velconstraint);
	}

	void SetNodalConstraintSpringDamper(int nn1, TArray<int> &lc, TArray<double> &aground, double r, double spring_stiffnessi, double damping_coeffi = 0, int velconstraint=0)
	{
		InitNC();

		//penalty = 1; 		//$ DR 2011-04-20: old code
		if(spring_stiffnessi) { SetPenaltyFormulation(1); } //$ DR 2011-04-20

		loccoords = lc;
		//ground = aground;
		ground.SetLen(aground.Length());
		for(int i=1; i<=ground.Length(); i++)
		{
			MathFunction* mathfun= new MathFunction();
			mathfun->SetConstant(aground(i));
			ground(i)= mathfun;
		}
		nodes.Add(nn1);
		SetDrawSizeScalar(r);	//draw_dim.X() = r;
		SetPenaltyStiffness(spring_stiffnessi); //$ DR 2011-04-21: old code: spring_stiffness = spring_stiffnessi;
		damping_coeff = damping_coeffi;

		elementname = GetElementSpec();
		x_init = Vector(SS());
		velocity_constraint = velconstraint;
		if (!velconstraint)
			for(int i=1; i<=ground.Length(); i++)
				x_init(i) = ground(i)->Evaluate(0);
		else
			for(int i=1; i<=ground.Length(); i++)
				x_init(SOS()+i) = ground(i)->Evaluate(0); 
	}


	//constraint local nodes (nn1,nn2) of 2 elements (elem1, elem2) and nodal direction (lc1,lc2)
	void SetLocalNodalConstraintSpringDamper(int elem1, int nn1, int lc1, int elem2, int nn2, int lc2,double r=-1, int penaltyi = 0, double spring_stiffnessi=0, double damping_coeffi=0, int velconstraint=0)
	{
		TArray<int> lcs1; lcs1.Add(lc1);
		TArray<int> lcs2; lcs2.Add(lc2);

		SetLocalNodalConstraintSpringDamper(elem1, nn1, lcs1, elem2, nn2, lcs2, r, penaltyi, spring_stiffnessi, damping_coeffi, velconstraint);

		//$ DR 2011-03-09: [ old code

		//InitNC();

		//if (penaltyi == 0) GetMBS()->UO() << "***** WARNING: penalty == 0 does not work in 'SetLocalNodalConstraintSpringDamper'! ******************\n";

		//penalty = penaltyi;

		//nodes.Add(nn1);
		//loccoords.Add(lc1);
		//elements.Add(elem1);

		//nodes.Add(nn2);
		//loccoords.Add(lc2);
		//elements.Add(elem2);

		//draw_dim.X() = r;
		//spring_stiffness = spring_stiffnessi;
		//damping_coeff = damping_coeffi;

		//elementname = GetElementSpec();
		//x_init = Vector(SS());
		//velocity_constraint = velconstraint;

		//$ DR 2011-03-09: ] old code

	}

		void SetLocalNodalConstraintSpringDamper(int elem1, int nn1, TArray<int> &lc1, int elem2, int nn2, TArray<int> &lc2,double r, int penaltyi = 0, double spring_stiffnessi=0, double damping_coeffi=0, int velconstraint=0)
	{
		InitNC();

		if (penaltyi == 0) GetMBS()->UO() << "***** WARNING: penalty == 0 does not work in 'SetLocalNodalConstraintSpringDamper'! ******************\n";

		//penalty = penaltyi; 		//$ DR 2011-04-20: old code
		if(penaltyi) { SetPenaltyFormulation(1); } //$ DR 2011-05-11

		nodes.Add(nn1);
		//loccoords.Add(lc1);
		loccoords.Add(lc1(1));
		elements.Add(elem1);

		nodes.Add(nn2);
		//loccoords.Add(lc2);
		loccoords.Add(lc2(1));
		elements.Add(elem2);

		if(lc1.Length()>1) loccoords.Add(lc1(2));
		if(lc2.Length()>1) loccoords.Add(lc2(2));
		if(lc1.Length()>2) loccoords.Add(lc1(3));
		if(lc2.Length()>2) loccoords.Add(lc2(3));

		//draw_dim.X() = r;
		SetDrawSizeScalar(r);
		SetPenaltyStiffness(spring_stiffnessi);			//$ DR 2011-04-21: old code: spring_stiffness = spring_stiffnessi;
		damping_coeff = damping_coeffi;

		elementname = GetElementSpec();
		x_init = Vector(SS());
		velocity_constraint = velconstraint;

	}



	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new NodalConstraint(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const NodalConstraint& ce = (const NodalConstraint&)e;
		nodes = ce.nodes;
		loccoords = ce.loccoords;
		//spring_stiffness = ce.spring_stiffness; //$ DR 2011-04-21: old code
		damping_coeff = ce.damping_coeff;
		//penalty = ce.penalty; //$ DR 2011-04-20: old code
		ground.SetLen(ce.ground.Length());

		dpdq = ce.dpdq;

		for(int i=1; i<=ground.Length(); i++)
		{
			MathFunction* mathfun= new MathFunction(*ce.ground(i));
			ground(i)= mathfun;
		}
		velocity_constraint = ce.velocity_constraint;

	}

	virtual void AddCoordGround(int loccoord, MathFunction* aground)
	{
		loccoords.Add(loccoord);
		ground.Add(aground);
	}

	virtual const char* GetElementSpec() const {return "NodalConstraint";}
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	//evaluate constraints: len(z)
	virtual void EvalG(Vector& f, double t);
	virtual void EvalF2(Vector& f, double t);

	virtual int SOS() const 
	{
		if (IsLocalNodeConstraint())
		{
			int nsos = 0;
			for (int i=1; i <= NE(); i++)
			{
				nsos += GetElem(i).SOS();
			}
			return nsos;
		}
		else return loccoords.Length();
	};  // explicit size, number of constrained dofs
	virtual int SOSowned() const {return 0;}  // number of explicit equations added by element
	//virtual int IS() const {return (1-penalty)*loccoords.Length();};  //implicit (algebraic) size (RL)
	virtual int IS() const 
	{ 
		if(UsePenaltyFormulation()) //$ DR 2011-04-20: old code: if (penalty)
		{
			return 0; 
		}
		else 
		{
			if (IsLocalNodeConstraint())
			{
				return 1; //only one coordinate is constrained!
			}
			else
			{
				return loccoords.Length();
			}
		}
	};  //more penalties added (valid values of "penalty" now 0..3) (AD)

	virtual int Dim() const 
	{
		if (elements.Length() != 0) 
		{
			return GetElem(1).Dim();
		}
		else if (nodes.Length() != 0) 
		{
			return GetMBS()->GetNode(nodes(1)).Dim();
		}
		else return 3;

	}  //has 3D position???

	////$ DR 2011-04:[ old code
	virtual void SetPenalty(int flag){SetPenaltyFormulation(flag); UO() << "Warning: SetPenalty will be replaced by SetPenaltyFormulation!!! \n";} // use penalty coordinate constraint (RL)
	//virtual int GetPenalty(){return penalty;} // use penalty coordinate constraint (RL)	
	virtual void SetSpringStiffness(double val){SetPenaltyStiffness(val); UO() << "Warning: SetSpringStiffness will be replaced by SetPenaltyStiffness!!! \n";} // use penalty coordinate constraint (RL)
	//virtual double GetSpringStiffness(){return spring_stiffness;} // use penalty coordinate constraint (RL)
	////$ DR 2011-04:] old code

	virtual void SetDampingCoeff(double val){damping_coeff = val;} // use penalty coordinate constraint (AD)
	virtual double GetDampingCoeff(){return damping_coeff;} // use penalty coordinate constraint (AD)

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const 
	{
		return GetNode(1).Pos();
	}

	virtual Node& GetNode(int i) 
	{
		if (IsLocalNodeConstraint()) 
		{ return GetElem(i).GetNode(NodeNum(i));} 
		else 
		{ return GetMBS()->GetNode(NodeNum(i)); }
	}
	virtual const Node& GetNode(int i) const
	{
		if (IsLocalNodeConstraint()) 
		{ return GetElem(i).GetNode(NodeNum(i));} 
		else 
		{ return GetMBS()->GetNode(NodeNum(i)); }
	}
	

	virtual void GetGlobalConstrainedDOFs(TArray<int>& dofs) const
	{
		//$ YV 2013-01-03: old contents of the list remains in there - we simply add new constrained dofs
		//dofs.Flush();
		int GlobalNr;
		if (IsLocalNodeConstraint()) // element & "local" node number
		{
			const Element& elem = GetMBS()->GetElement(elements(1));
			for(int i=1; i<= loccoords.Length(); i++)
			{
				GlobalNr = elem.LTG(loccoords(i));
				dofs.Add(GlobalNr);
			}
		}
		else // global node number
		{
			int nodenr = nodes(1);
			const Node& node = GetMBS()->GetNode(nodenr);

			int offset = 0;
			if(IsVelocityConstraint()) 
				offset = node.SOS(); // run test with a velocity constraint...

			for(int i=1; i<= loccoords.Length(); i++)
			{
				GlobalNr = node.LTG(loccoords(i+offset));
				dofs.Add(GlobalNr);
			}
		}
	}

	virtual void DrawElement();


	//do not use NNodes, this is assumed to be part of a finite element
	virtual int NConstrainedNodes() const {return nodes.Length(); };

	virtual const int& NodeNum(int i) const { return nodes(i); }
	virtual int& NodeNum(int i) { return nodes(i); }

	virtual void LinkToElements() ;

	virtual void SetGlobalInitConditions(Vector& x_glob)
	{
	}

	virtual int IsVelocityConstraint() const {return velocity_constraint;}

	virtual int IsLocalNodeConstraint() const {return NE() != 0;}

	//$ YV 2011-05-03[
	// it is sometimes necessary to access/modify the ground during computation
	MathFunction * & Ground(int i) { return ground(i); }
	//$ YV 2011-05-03]

	//$ AD 2011-09-08
	double GetPenaltyStiffness(double t=-1.);  // time dependent stiffness factor (computation steps)

protected:
	// node degrees of freedom
	TArray<int> loccoords;
	TArray<int> nodes;
	// values to which node position is set
	TArray<MathFunction*> ground;
	// number of node in CMS element
	//int node;

	Matrix dpdq; //temporary

	//$ DR 2011-04-20: old code
	//// use penalty method
	//int penalty; //use penalty formulation instead of lagrange parameter (RL)

	// stiffness in case of penalty method
	//double spring_stiffness;
	// damping of spring
	double damping_coeff;
	// apply constraint to velocities, not positions
	int velocity_constraint;
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Spherical: BasePointJoint
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: BasePointJoint (Constraint)
// short description: 
// available formulations: Lagrange multiplier, Penalty Formulation
// type of constraint equation: nonlinear constraint equations
// index formulation available: 
// ground constraint: yes
// element-element constraint: yes
// development status: in progress
// long description:	should replace spherical joint and nodal constraint in future
//								
// class variables:
//  
//**end(ued)**
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// DR is responsible for this constraint
class BasePointJoint: public Constraint  //$EDC$[beginclass,classname=BasePointJoint,parentclassname=Constraint,addelementtype=TAEconstraint,addelementtypename=PointJoint,
//texdescription="The PointJoint constrains two elements at a local position or node each. If only one element is specified (second element 0), a ground PointJoint is realized.",
//modus="{element to ground}{Position2.element\_number AND Position2.node\_number have to be equal to 0}",
//modus="{element to element}{Position2.element\_number and/or Position2.node\_number must not be equal to 0}",
//modus="{Lagrange}{If Physics.use\_penalty\_formulation = 0, than no stiffness and no damping parameters are used.}",example="PointJointShort.txt"]
{
public:
	//$ PG 2013-4-17: [ remove these functions again!!!  just for testing
	virtual TArray<int>& GetElementsArray() { return elements; }
	virtual TArray<Vector3D>& GetLoccoordsArray() { return loccoords; }
	//$ PG 2013-4-17: ]

	BasePointJoint(MBS* mbsi):Constraint(mbsi)
	{	
		ElementDefaultConstructorInitialization();
	};

	BasePointJoint(const BasePointJoint& ct): Constraint(ct.mbs)
	{
		ElementDefaultConstructorInitialization();
		CopyFrom(ct);
	};

	~BasePointJoint()
	{
	};

	//virtual void InitBasePointJoint();	//$ DR 2012-07-26: renamed to: ElementDefaultConstructorInitialization
	virtual void ElementDefaultConstructorInitialization();

	//$ DR 2012-07: CheckConsistency added
	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute

	// for constraints only, these are the necessary (!) access functions!	//$ DR 2013-02-13
	virtual void GetNecessaryKinematicAccessFunctions(TArray<int> &KAF, int numberOfKinematicPair);

	virtual Element* GetCopy()
	{
		Element* ec = new BasePointJoint(*this);
		return ec;
	}

	virtual void CopyFrom(const Element& e);

	virtual const char* GetElementSpec() const {return "PointJoint";}

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer
  virtual int GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //for special value sensor
	virtual int ReadSingleElementData(ReadWriteElementDataVariableType& RWdata); 		//for special value sensor


	virtual void GetElementDataAuto(ElementDataContainer& edc); 		
	//virtual int SetElementDataAuto(const ElementDataContainer& edc);
	virtual int SetElementDataAuto(ElementDataContainer& edc);
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	void BasePointJoint::DrawElement();

	virtual void EvalG(Vector& f, double t);
	virtual void EvalF2(Vector& f, double t);
	virtual Vector3D ComputeForce(double t) const;
	virtual double GetActorForce(double computation_time, int dir=0) const; //get actor force, with optional direction (x,y,z, etc.)
	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);


	// spherical verwendet das von constraint (wenn lagrange dann Element::SetGlobalInitConditions(x_glob);)
	// nodal hat das leere wie hier
	virtual void SetGlobalInitConditions(Vector& x_glob)
	{
	}

	//this function is called after assembly and AFTER filysis
	void Initialize() 
	{
		if(auto_comp_ground)
		{
			SetPos2ToGlobalCoord(GetPos1());
		}

		// if general (scalar) penalty stiffness is used, set Stiffness-Vector
		if(GetPenaltyStiffness() && UsePenaltyFormulation())
		{
			SetPenaltyStiffness3(Vector3D(GetPenaltyStiffness()));
		}
	}

	virtual	void GetdRotTvdqT(const Vector3D& vloc, const Vector3D& ploc, Matrix& d, int bodyindex); // MSax: GetdRotTvdqT is derived from GetdRotvdqT

	virtual int SOS() const;
	virtual int SOSowned() const {return 0;}  // number of explicit equations added by element

	virtual int IS() const 
	{ 
		if(UsePenaltyFormulation())
		{
			return 0; 
		}
		else 
		{
			return GetNumberOfConstrainedCoords();
		}
	}; 

	virtual int Dim() const {return 3;}		// what is it for? constraint = 2, nodal: 3 or Dim of Element, spherical: 3

	////// wenn nur local Nodes/Coords --> wie in constraint --> kann man entfernen

	virtual void LinkToElements() 
	{
		if (UsePenaltyFormulation()) 
		{
			LinkToElementsPenalty();
		}
		else
		{
			LinkToElementsLagrange();		// wie in constraint --> kann man entfernen
		}
		//UO(UO_LVL_dbg1) << "LTG: \n" << ltg << "\n";
	}

	virtual void LinkToElementsPenalty();
	virtual void LinkToElementsLagrange();


		//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const 
	{
		if(elements(1)) 
		{
			if((mbs->GetElementPtr(elements(1)))->GetType() >= TCMS)	// CMS-Element
			{
				return GetPos1();
			}
			else
			{
				return (GetMBS()->GetElement(elements(1))).GetRefPos();
			}
		}
		else		// global node
		{
			return (GetMBS()->GetNode(nodes(1))).Pos();
		}
	}

	Vector3D GetPos1() const {return GetPosition(1);}
	Vector3D GetPos2() const {return GetPosition(2);}
	Vector3D GetVel1() const {return GetVelocity(1);}
	Vector3D GetVel2() const {return GetVelocity(2);}
	void SetPos1ToGlobalNode(int glob_node_nr){SetPosToGlobalNode(1, glob_node_nr);}
	void SetPos2ToGlobalNode(int glob_node_nr){SetPosToGlobalNode(2, glob_node_nr);}
	void SetPos1ToLocalNode(int elem_nr, int loc_node_nr){SetPosToLocalNode(1, elem_nr, loc_node_nr);}
	void SetPos2ToLocalNode(int elem_nr, int loc_node_nr){SetPosToLocalNode(2, elem_nr, loc_node_nr);}
	void SetPos1ToLocalCoord(int elem_nr, Vector3D loccoord){SetPosToLocalCoord(1, elem_nr, loccoord);}
	void SetPos2ToLocalCoord(int elem_nr, Vector3D loccoord){SetPosToLocalCoord(2, elem_nr, loccoord);}
	void SetPos2ToGlobalCoord(Vector3D ground);
	void AutoComputeGround();
	void SetDisplacement(TArray<MathFunction*>& disp);
	

	// velocity constraint, even when Index 3 solver is used
	virtual int IsVelocityConstraint() const {return velocity_constraint;}	
	virtual void SetConstrainedDirections(IVector directions) { dir = directions;}	//$ MSax 2013-01: new SetFunction
	virtual void SetConstrainedDirections(Vector3D directions)	//$ DR 2013-02-19 bugfix, keep the old set function alive
	{ 
		dir(1) = (int)directions.X();
		dir(2) = (int)directions.Y();
		dir(3) = (int)directions.Z();
	}
	virtual IVector GetConstrainedDirections() {return dir;}

	virtual int GetNumberOfConstrainedCoords() const 
	{
		if(UsePenaltyFormulation())		// not used at all ?
		{
			return 3;					
		}
		else
		{
			int counter = 0;
			for(int i=1; i<=3; i++)
			{
				if(dir(i)!=0) counter++;
			}
			return counter;
		}
	}	

	virtual double GetDampingCoeff() const {return damping_coeff;}
	virtual void SetDampingCoeff(double damp) {	damping_coeff = damp;	}
	virtual int UseDamping() const {
		if(damping_coeff!=0.) {return 1;}
		else {return 0;};
	}

	//return 3 stiffness parameters for penalty formulation
	virtual Vector3D GetPenaltyStiffness3(double t) const;

	//set 3 stiffness parameters for penalty formulation
	virtual void SetPenaltyStiffness3(Vector3D stiffness3i) 
	{
		spring_stiffness3 = stiffness3i; 
		SetPenaltyFormulation(1);
	}

	// if flag stiffness_in_joint_local_frame is set:	use rotated (local or global) coordinate system to define stiffness
	virtual void SetJointLocalFrame(Matrix3D A)
	//virtual void SetStiffnessInRotatedGlobCoordSys(Matrix3D A)
	{
		JA0i = A;
		stiffness_in_joint_local_frame = 1;
	}

	// drawing dimensions of joint local frame. If set to -1, than global_draw_scalar_size is used. If set to 0, than no joint local frame is drawn.
	virtual void SetDrawSize_JointLocalFrame (double size)
	{
		draw_local_frame_size = size;
	}

	virtual Matrix3D GetRotMati() const;
	virtual Matrix3D GetRotMatiP() const;
	virtual Matrix3D GetRotMatiD() const;

	virtual int IsLocalNodeConstraint() const 
	{
		////return NE() != 0;
		//int flag = 1;
		//for(int i=1; i<=elements.Length(); i++)
		//{
		//	if(elements(i)==0) flag = 0;				// global nodes are added with element-nr = 0
		//}
		//return flag;
		if(elements(1))				// local node or local coordinate
		{
			return 1;
		}
		else									// global node
		{
			return 0;
		}
	}

	virtual int NE_nodouble() const {
		int counter = 0;
		for(int i=1; i<=elements.Length(); i++)
		{
			if(elements(i)!=0) counter++;				// global nodes are added with element-nr = 0
		}
		return counter;
	}

	virtual int NE() const {return NE_nodouble();}		// returns 0 if global nodes to ground
	virtual int NKinPairs() const 
	{/*return elements.Length();*/
		int counter = 0;
		for(int i=1; i<=elements.Length(); i++)
		{
			if(elements(i)!=0) counter++;		// local node or local coordinate	
			else														// global nodes are added with element-nr = 0
			{
				if((nodes(i)!=0)) counter++;
			}
		}
		return counter;
	}

	virtual void SetDrawingOptions(double cdim = -1, const Vector3D& coli = inherited_color)
	{
		SetDrawSizeScalar(cdim);
		GetCol() = coli;
	}
	Vector3D GetDrawPosition(int i);

protected:
	//	=== flags ===
	//int use_local_coordinate_system;		// flag, which is used to define the stiffness w.r.t. local coordinate system	--> already in constraint
	//int use_penalty_formulation;				// flag, which is used to switch constraints into penalty mode (if available)--> already in constraint
	int velocity_constraint;						// flag, which is used to apply constraint to velocities, not positions
	int auto_comp_ground;								// flag to define if position of ground is computed automatically
	int stiffness_in_joint_local_frame;	//$EDC$[varaccess,EDCvarname="use_joint_local_frame",EDCfolder="Geometry",int_bool,tooltiptext="Use a special joint local frame"]

	//  === data ===
	//EDC int elements(1);							//$EDC$[varaccess,EDCvarname="element_number",minval=1,EDCfolder="Position1",tooltiptext="Number of constrained element"]
	//EDC int elements(2);							//$EDC$[varaccess,EDCvarname="element_number",minval=0,EDCfolder="Position2",tooltiptext="Number of constrained element"]

	TArray<Vector3D> loccoords;					// local coordinates of point in element, or (if ground joint) global coordinate of ground
	//EDC Vector3D loccoords(1);				//$EDC$[varaccess,EDCvarname="position",EDCfolder="Position1",tooltiptext="local position. Only used if node_number == 0!"]
	//EDC Vector3D loccoords(2);				//$EDC$[varaccess,EDCvarname="position",EDCfolder="Position2",tooltiptext="local or global (if element_number == 0) position. Only used if node_number == 0!"]
	
	TArray<int> nodes;									// node numbers
	//EDC int nodes(1);									//$EDC$[varaccess,EDCvarname="node_number",minval=1,EDCfolder="Position1",tooltiptext="local or global (if element_number == 0) node number."]
	//EDC int nodes(2);									//$EDC$[varaccess,EDCvarname="node_number",minval=1,EDCfolder="Position2",tooltiptext="local or global (if element_number == 0) node number."]

	//Vector3D dir;												//!EDC$[varaccess,EDCvarname="constrained_directions",minval=0,maxval=1,EDCfolder="Physics.Lagrange",tooltiptext="[x,y,z]...(1 = constrained, 0 = free), can be defined as local or global directions (see Geometry)"]
	TArray<int> dir;											//$EDC$[varaccess,EDCvarname="constrained_directions",minval=0,maxval=1,EDCfolder="Physics.Lagrange",tooltiptext="[x,y,z]...(1 = constrained, 0 = free), can be defined as local or global directions (see Geometry)"]


	TArray<MathFunction*> displ;				// global displacement of ground if ground is not constant; (1 mathfunction per direction)		
	Vector3D spring_stiffness3;					//$EDC$[varaccess,EDCvarname="spring_stiffness_vector",EDCfolder="Physics.Penalty",tooltiptext="penalty stiffness parameter [kx,ky,kz]. Just used if scalar spring_stiffness == 0, otherwise kx=ky=kz=spring_stiffness"]
	Matrix3D JA0i;											//$EDC$[varaccess,EDCvarname="joint_local_frame",EDCfolder="Geometry",tooltiptext="Prerotate stiffness vector w.r.t. global coordinate system or local coordinate system of body 1. Just used if use_joint_local_frame == 1"]

	//EDC int use_local_coordinate_system;	//$EDC$[varaccess,remove,EDCvarname="use_local_coordinate_system",EDCfolder="Geometry"]
	//EDC int use_local_coordinate_system;	//$EDC$[varaccess,EDCvarname="use_local_coordinate_system",EDCfolder="Geometry",int_bool,tooltiptext="0=use global coordinates, 1=use local coordinate system of Body 1"]

	double damping_coeff;								//$EDC$[varaccess,EDCvarname="damping",minval=0,EDCfolder="Physics.Penalty",tooltiptext="damping coefficient for viscous damping (F = d*v), applied in all constrained directions"]
	Matrix dpdq;
	Matrix drotTvdq;

	//EDC double draw_dim(1);						//$EDC$[varaccess,EDCvarname="draw_size",EDCfolder="Graphics",tooltiptext="drawing dimensions of constraint. If set to -1, than global_draw_scalar_size is used."]
	double	draw_local_frame_size;			//$EDC$[varaccess,EDCvarname="draw_size_joint_local_frame",EDCfolder="Graphics",tooltiptext="drawing dimensions of joint local frame. If set to -1, than global_draw_scalar_size is used. If set to 0, than no joint local frame is drawn."]

public:
	//Vector3D GetDrawPosition(int i);		//$ DR 2011-11-2: not private anymore
	Vector3D GetPosition(int i) const;
	Vector3D GetVelocity(int i) const;
	void SetPosToLocalNode(int i, int elem_nr, int loc_node_nr);
	void SetPosToLocalCoord(int i, int elem_nr, Vector3D loccoord);
	void SetPosToGlobalNode(int i,  int glob_node_nr);

};//$EDC$[endclass,BasePointJoint]

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Spherical: TSphericalJoint
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: TSphericalJoint (Constraint)
// short description: 
// available formulations: Lagrange multiplier, Penalty Formulation
// type of constraint equation: nonlinear constraint equations
// index formulation available: 
// ground constraint: yes
// element-element constraint: yes
// development status: finished
// long description:	constructors and set functions for the old spherical joint
//										computation is done by basepointjoint
//								
// class variables: none
//  
//**end(ued)**
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <class T>
class TSphericalJoint: public T		
{
	public:

		TSphericalJoint(MBS *mbs):T(mbs)
		{
		}

		// ###################################
		// constructors of old SphericalJoint:

		// ---- element to element : Lagrange -----------
		TSphericalJoint(MBS* mbs, int en1, int en2, const Vector3D& lc1, const Vector3D& lc2, double cdim = -1 , const Vector3D& coli = inherited_color):T(mbs)
		{
			SetSphericalJoint_LocalPos_to_LocalPos(en1, lc1, en2, lc2);
			SetDrawingOptions(cdim, coli);
		}

		// ---- element to element : Penalty  -----------
		TSphericalJoint(MBS *mbs, int en1, int en2, const Vector3D& lc1, const Vector3D& lc2, Vector3D spring_stiffness, double cdim = -1, const Vector3D& coli = inherited_color):T(mbs)
		{
			SetSphericalJoint_LocalPos_to_LocalPos(en1, lc1, en2, lc2, Vector3D(spring_stiffness));
			SetDrawingOptions(cdim, coli);
		}

		// ---- element to ground  : Lagrange -----------
		TSphericalJoint(MBS* mbs, int en1, const Vector3D& lc1, const Vector3D& pglob, double cdim = -1, const Vector3D& coli = inherited_color):T(mbs)
		{
			SetSphericalJoint_LocalPos_to_GlobalPos(en1, lc1, pglob); 
			SetDrawingOptions(cdim, coli);
		}

		// ---- element to ground  : Penalty  -----------
		TSphericalJoint(MBS* mbs, int en1, const Vector3D& lc1, const Vector3D& pglob, Vector3D spring_stiffness, double cdim = -1, const Vector3D& coli = inherited_color):T(mbs)
		{
			SetSphericalJoint_LocalPos_to_GlobalPos(en1, lc1, pglob,spring_stiffness); 
			SetDrawingOptions(cdim, coli);
		}

		// ###################################
		// SetFunctions

		// Penalty + Lagrange: Element to Element 
		void SetSphericalJoint_LocalNode_to_LocalNode(int elem1, int node1, int elem2, int node2, Vector3D stiffness = NULL, int islocalStiffness=0, double damping = 0.0)
		{
			SetPos1ToLocalNode(elem1, node1);
			SetPos2ToLocalNode(elem2, node2);
			if(!(stiffness==NULL)) SetPenaltyStiffness3(stiffness);
			else 
			{
				IVector dir;
				dir.Set3(1,1,1);
				SetConstrainedDirections(dir);
			}
			SetUseLocalCoordinateSystem(islocalStiffness);
			SetDampingCoeff(damping);
		}

		// Penalty + Lagrange: Element to Ground
		void SetSphericalJoint_LocalNode_to_GlobalPos(int elem1, int node1, Vector3D ground = NULL, Vector3D stiffness = NULL, int islocalStiffness=0, double damping = 0.0, TArray<MathFunction*>& displacement = TArray<MathFunction*>(0))
		{
			SetPos1ToLocalNode(elem1, node1);
			if(!(ground == NULL))	SetPos2ToGlobalCoord(ground);
			else	AutoComputeGround();

			if(!(stiffness==NULL)) SetPenaltyStiffness3(stiffness);
			else 
			{
				IVector dir;
				dir.Set3(1,1,1);
				SetConstrainedDirections(dir);
			}
			if(displacement.Length() != 0) SetDisplacement(displacement);
			SetUseLocalCoordinateSystem(islocalStiffness);
			SetDampingCoeff(damping);
		}

		// Penalty + Lagrange: Element to Ground
		void SetSphericalJoint_LocalPos_to_GlobalPos(int elem1, Vector3D lc1, Vector3D ground = NULL, Vector3D stiffness = NULL, int islocalStiffness=0, double damping = 0.0,TArray<MathFunction*>& displacement = TArray<MathFunction*>(0))
		{
			SetPos1ToLocalCoord(elem1, lc1);
			if(!(ground == NULL))	SetPos2ToGlobalCoord(ground);
			else	AutoComputeGround();

			if(!(stiffness==NULL)) SetPenaltyStiffness3(stiffness);
			else
			{
				IVector dir;
				dir.Set3(1,1,1);
				SetConstrainedDirections(dir);
			}
			if(displacement.Length() != 0) SetDisplacement(displacement);
			SetUseLocalCoordinateSystem(islocalStiffness);
			SetDampingCoeff(damping);
		}

		// Penalty + Lagrange: Element to Element
		void SetSphericalJoint_LocalPos_to_LocalPos(int elem1, Vector3D lc1, int elem2, Vector3D lc2, Vector3D stiffness = NULL, int islocalStiffness=0, double damping = 0.0)
		{
			SetPos1ToLocalCoord(elem1, lc1);
			SetPos2ToLocalCoord(elem2, lc2);
			if(!(stiffness==NULL)) SetPenaltyStiffness3(stiffness);
			else
			{
				IVector dir;
				dir.Set3(1,1,1);
				SetConstrainedDirections(dir);
			}
			SetUseLocalCoordinateSystem(islocalStiffness);
			SetDampingCoeff(damping);
		}
		
		// Penalty + Lagrange: Element to Element
		void SetSphericalJoint_LocalPos_to_LocalNode(int elem1, Vector3D lc1, int elem2, int node2, Vector3D stiffness = NULL, int islocalStiffness=0, double damping = 0.0)
		{
			SetPos1ToLocalCoord(elem1, lc1);
			SetPos2ToLocalNode(elem2, node2);
			if(!(stiffness==NULL)) SetPenaltyStiffness3(stiffness);
			else 
			{
				IVector dir;
				dir.Set3(1,1,1);
				SetConstrainedDirections(dir);
			}
			SetUseLocalCoordinateSystem(islocalStiffness);
			SetDampingCoeff(damping);
		}

		// Penalty + Lagrange: Element to Element
		void SetSphericalJoint_LocalNode_to_LocalPos(int elem1, int node1, int elem2, Vector3D lc2, Vector3D stiffness = NULL, int islocalStiffness=0, double damping = 0.0)
		{
			SetPos1ToLocalNode(elem1, node1);
			SetPos2ToLocalCoord(elem2, lc2);
			if(!(stiffness==NULL)) SetPenaltyStiffness3(stiffness);
			else 
			{
				IVector dir;
				dir.Set3(1,1,1);
				SetConstrainedDirections(dir);
			}
			SetUseLocalCoordinateSystem(islocalStiffness);
			SetDampingCoeff(damping);
		}

		// Penalty + Lagrange: Element to Element
		void SetSphericalJoint_GlobalNode_to_GlobalNode(int node1, int node2, Vector3D stiffness = NULL , double damping = 0.0)
		// not implemented yet for Lagrange
		{
			SetPos1ToGlobalNode(node1);
			SetPos2ToGlobalNode(node2);
			if(!(stiffness==NULL)) SetPenaltyStiffness3(stiffness);
			else 
			{
				IVector dir;
				dir.Set3(1,1,1);
				SetConstrainedDirections(dir);
			}
			//SetUseLocalCoordinateSystem(islocalStiffness);	// not possible for global nodes, because no element is specified
			SetDampingCoeff(damping);
		}

		// Penalty + Lagrange: Element to Ground
		void SetSphericalJoint_GlobalNode_to_GlobalPos(int node1, Vector3D ground = NULL, Vector3D stiffness = NULL, double damping = 0.0, TArray<MathFunction*>& displacement = TArray<MathFunction*>(0) )
		{
			//// DR moved to CheckConsistency
			//if(stiffness == NULL)
			//{
			//	GetMBS()->UO(UO_LVL_err) << "ERROR: SphericalJoint: global node to ground not implemented yet for Lagrange multiplier\n";
			//	assert(0);
			//}

			SetPos1ToGlobalNode(node1);
			if(!(ground == NULL))	SetPos2ToGlobalCoord(ground);
			else	AutoComputeGround();

			if(displacement.Length() != 0) SetDisplacement(displacement);
			if(!(stiffness==NULL)) SetPenaltyStiffness3(stiffness);
			else
			{
				IVector dir;
				dir.Set3(1,1,1);
				SetConstrainedDirections(dir);
			}
			//SetUseLocalCoordinateSystem(islocalStiffness);	// not possible for global nodes, because no element is specified
			SetDampingCoeff(damping);
		}


};




// ========================================================================================================================


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Spherical: Spherical Joint 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: SphericalJointDEPRECATED (Constraint)
// short description: 
// available formulations: Lagrange multiplier, Penalty Formulation
// type of constraint equation: nonlinear constraint equations
// index formulation available: 
// ground constraint: yes
// element-element constraint: yes
// development status: 
// long description:	all translatory d.o.f. are constrained
//								
// class variables:
//  
//**end(ued)**
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//#define use_sphericaljoint_new	//$!DR 2011-07-06: removed old code
//#ifdef use_sphericaljoint_old		//$!DR 2011-07-06: removed old code

class SphericalJointDEPRECATED: public Constraint
{
public:

	// Spherical Joint: Ground Joint & Lagrange
	SphericalJointDEPRECATED(MBS* mbsi, int en1, const Vector3D& lc1, const Vector3D& pglob, 
		double cdim = -1, const Vector3D& coli = inherited_color):Constraint(mbsi), loccoords(), dpdq()
	{	
		p_global = pglob;
		AddElementCoord(en1, lc1);
		SetSphericalJointParameters(0.,cdim,coli);
	};

	// Spherical Joint: Ground Joint & penalty_3directons
	SphericalJointDEPRECATED(MBS* mbsi, int en1, const Vector3D& lc1, const Vector3D& pglob, Vector3D spring_stiffness, 
		double cdim = -1, const Vector3D& coli = inherited_color):Constraint(mbsi), loccoords(), dpdq()
	{	
		p_global = pglob;
		AddElementCoord(en1, lc1);
		SetSphericalJointParameters(spring_stiffness,cdim,coli);
	};

	// Spherical Joint: 2 Elements & Lagrange
	SphericalJointDEPRECATED(MBS* mbsi, int en1, int en2, const Vector3D& lc1, const Vector3D& lc2, 
		double cdim = -1 , const Vector3D& coli = inherited_color):Constraint(mbsi), loccoords(), dpdq()
	{	
		AddElementCoord(en1, lc1);
		AddElementCoord(en2, lc2);
		SetSphericalJointParameters(0.,cdim,coli);
	};

	// Spherical Joint: Ground Joint & penalty_3directons
	SphericalJointDEPRECATED(MBS* mbsi, int en1, int en2, const Vector3D& lc1, const Vector3D& lc2, Vector3D spring_stiffness, 
		double cdim = -1, const Vector3D& coli = inherited_color):Constraint(mbsi), loccoords(), dpdq()
	{	
		AddElementCoord(en1, lc1);
		AddElementCoord(en2, lc2);
		SetSphericalJointParameters(spring_stiffness,cdim,coli);
	};

	//$!DR 2011-06-30:[ === old constructors: do not use them any more!

	//SphericalJoint(MBS* mbsi, int en1, const Vector3D& lc1, const Vector3D& pglob):Constraint(mbsi), loccoords(), dpdq()
	//{
	//	SphericalJoint(mbsi,  en1, lc1, pglob, GetDrawSizeScalar(), GetCol());
	//};

	//// Spherical Joint: Ground Joint & penalty_uniform
	//	SphericalJoint(MBS* mbsi, int en1, const Vector3D& lc1, 
	//		const Vector3D& pglob, double spring_stiffness, double cdim = -1, const Vector3D& coli = inherited_color):Constraint(mbsi), loccoords(), dpdq()
	//	{	
	//		p_global = pglob;
	//		AddElementCoord(en1, lc1);
	//		SetSphericalJointParameters(Vector3D(spring_stiffness,spring_stiffness,spring_stiffness),cdim,coli);
	//	};

	//SphericalJoint(MBS* mbsi, int en1, const Vector3D& lc1, const Vector3D& pglob, double spring_stiffness):Constraint(mbsi), loccoords(), dpdq()
	//{
	//	SphericalJoint(mbsi,  en1, lc1, pglob, spring_stiffness, GetDrawSizeScalar(), GetCol());
	//};

	//SphericalJoint(MBS* mbsi, int en1, const Vector3D& lc1, const Vector3D& pglob, Vector3D spring_stiffness):Constraint(mbsi), loccoords(), dpdq()
	//{
	//	SphericalJoint(mbsi, en1, lc1, pglob, spring_stiffness, GetDrawSizeScalar(), GetCol());
	//};

	//SphericalJoint(MBS* mbsi, int en1, int en2, const Vector3D& lc1, const Vector3D& lc2):Constraint(mbsi), loccoords(), dpdq()
	//{
	//	//SphericalJoint(mbsi, en1, en2, lc1, lc2, GetDrawSizeScalar(), GetCol());		//$ DR 2011-06-08: bugfix, replaced this (old) line by the following 3 (new) ones
	//	AddElementCoord(en1, lc1);
	//	AddElementCoord(en2, lc2);
	//	SetSphericalJointParameters(0.,GetDrawSizeScalar(),GetCol());
	//};


	//// Spherical Joint: Ground Joint & penalty_uniform
	//	SphericalJoint(MBS* mbsi, int en1, int en2, 
	//		const Vector3D& lc1, const Vector3D& lc2, double spring_stiffness, double cdim = -1, const Vector3D& coli = inherited_color):Constraint(mbsi), loccoords(), dpdq()
	//	{	
	//		AddElementCoord(en1, lc1);
	//		AddElementCoord(en2, lc2);
	//		SetSphericalJointParameters(Vector3D(spring_stiffness,spring_stiffness,spring_stiffness),cdim,coli);
	//	};

	//SphericalJoint(MBS* mbsi, int en1, int en2, const Vector3D& lc1, const Vector3D& lc2, double spring_stiffness):Constraint(mbsi), loccoords(), dpdq()
	//{
	//	SphericalJoint(mbsi,  en1,  en2, lc1, lc2, spring_stiffness, GetDrawSizeScalar(), GetCol());
	//};

	//SphericalJoint(MBS* mbsi, int en1, int en2, const Vector3D& lc1, const Vector3D& lc2, Vector3D spring_stiffness):Constraint(mbsi), loccoords(), dpdq()
	//{
	//	SphericalJoint(mbsi, en1, en2, lc1, lc2, spring_stiffness, GetDrawSizeScalar(), GetCol());
	//};

//$!DR 2011-06-30:] === old constructors: do not use them any more!

	void SetSphericalJointParameters(Vector3D spring_stiffness, double cdim=-1., const Vector3D& coli=NULL)
	{
		mbs->UO() << "######################################################################## \n";
		mbs->UO() << "########## WARNING: you are using SphericalJointDEPRECATED ############# \n";
		mbs->UO() << "###### change typedef in KinematicPairs.h to use new version ########### \n";
		mbs->UO() << "######################################################################## \n";

		if (cdim != -1.)
			SetDrawSizeScalar(cdim);		//draw_dim.X() = cdim;
		if (&coli != NULL)
			GetCol() = coli;

		elementname = GetElementSpec();
		if (spring_stiffness.X())
		{
			SetPenaltyFormulation(1);
			SetPenaltyStiffness3(spring_stiffness);
		}
		else
			SetPenaltyFormulation(0);

		x_init = Vector(SS());  // correct System size
	};

	//void SetSphericalJointParameters(Vector3D spring_stiffness)
	//{
	//	SetSphericalJointParameters(spring_stiffness, GetDrawSizeScalar(), GetCol());
	//};

// 3D spring stiffness
	//return 3 stiffness parameters for penalty formulation
	virtual Vector3D GetPenaltyStiffness3(double t=-1.) const //$ AD 2011-09-08 computation steps for spherical joint {return spring_stiffness3;}
	{
		if(t < 0.)
			return spring_stiffness3;                // case no t passed for evaluation 

		if(NSteps()==0)
			return spring_stiffness3;
		else
			return spring_stiffness3 * Constraint::GetCStepsFact(t);
	}

	//set 3 stiffness parameters for penalty formulation
	virtual void SetPenaltyStiffness3(Vector3D stiffness3i) {spring_stiffness3 = stiffness3i;}
// 3D damping
	//return 3 damping parameters for penalty formulation
	virtual Vector3D GetPenaltyDamping3() const {return damping3;}
	//set 3 damping parameters for penalty formulation
	virtual void SetPenaltyDamping3(Vector3D damping3i) { damping3 = damping3i;}

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new SphericalJointDEPRECATED(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const SphericalJointDEPRECATED& ce = (const SphericalJointDEPRECATED&)e;
		loccoords = ce.loccoords;
		p_global = ce.p_global;
		spring_stiffness3 = ce.spring_stiffness3;
		damping3 = ce.damping3;
		dpdq = ce.dpdq;
	}

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	virtual const char* GetElementSpec() const {return "SphericalJointDEPRECATED";} //this name names MUST coincide with the type names defined in MBS::MBS() !!!!

	virtual void AddElementCoord(int en, const Vector3D& loccoord)
	{
		AddElement(en);
		loccoords.Add(loccoord);
	}

	virtual void EvalG(Vector& f, double t);
	virtual Vector3D ComputeForce(double t=-1.) const;
	virtual void EvalF2(Vector& f, double t); //second order equations for constraints: M \ddot u = F2(u,\dot u,t), len(u)

	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size
	//virtual int IS() const {return 3;};
	virtual int IS() const {return 3*(1-this->UsePenaltyFormulation());}
	virtual int Dim() const {return 3;}

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;

	virtual void DrawElement();

	int SphericalJointDEPRECATED::GetDrawSizeResolution()
	{
		if(draw_dim.Z()==-1)
			return (int)(GetMBS()->GetDOption(173));
		else
			return (int)draw_dim.Z();
	}

	void SphericalJointDEPRECATED::SetDrawSizeResolution(int drawdim_resolution)
	{
		draw_dim.Z()= (double)drawdim_resolution;
	}

protected:
	TArray<Vector3D> loccoords;
	Vector3D p_global;
	Vector3D spring_stiffness3; //penalty or general stiffness parameters for each axis (if available)
	Vector3D damping3; //penalty or general damping parameters for each axis (if available)

	Matrix dpdq; //temporary element, no need to copy???

};

//=========== use SphericalJointDEPRECATED ======================
// do not use SphericalJointDEPRECATED anymore, it will be removed soon! 
// report bugs of TSphericalJoint or BasePointJoint to DR

typedef TSphericalJoint<BasePointJoint> SphericalJoint;
//typedef SphericalJointDEPRECATED SphericalJoint;
//==========================================================


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Cylindrical	Cylindrical	Cylindrical	Cylindrical	Cylindrical	Cylindrical			
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Cylindrical: CylindricalPointJoint
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: CylindricalPointJoint (Constraint)
// short description: 
// available formulations: Lagrange multiplier
// type of constraint equation: leads to nonlinear constraint equations
// index formulation available: 
// ground constraint: no
// element-element constraint: yes
// development status: get+setelementdata not implemented yet
// long description:	first body must have rotation matrix, must be rigid
//										second body only needs local to global transformation, can be deformable
// class variables:
//  
//**end(ued)**
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//Constraint for a position measured relative within one or two bodies

class CylindricalPointJoint: public Constraint
{
public:
	CylindricalPointJoint(MBS* mbsi):Constraint(mbsi), loccoords(), dpdq() 
	{
	};

	//relative motion of two bodies along axis lp1-lp2, lp1 and lp2 must never be equal!!
	CylindricalPointJoint(MBS* mbsi, int en1, int en2, 
		const Vector3D& lp1, const Vector3D& lp2,
		const Vector3D& ddim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq(), hmat(), hvec() 
	{	
		x_init = Vector(SS());
		GetCol() = coli;
		//draw_dim = ddim;
		SetDrawSizeScalar(ddim.X());
		SetDrawSizeAxisLength(ddim.Y());
		AddElement(en1);
		AddElement(en2);
		loccoords.Add(lp1);
		loccoords.Add(lp2);
		loccoords.Add(Vector3D(0));
		loccoords.Add(Vector3D(0));
		elementname = GetElementSpec();
	};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new CylindricalPointJoint(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const CylindricalPointJoint& ce = (const CylindricalPointJoint&)e;
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

	virtual int Dim() const {return 3;}

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;

	virtual void DrawElement();

	double CylindricalPointJoint::GetDrawSizeAxisLength()
	{
		if(draw_dim.Y()==-1)
			return GetMBS()->GetDOption(172);
		else
			return draw_dim.Y();
	}
	void CylindricalPointJoint::SetDrawSizeAxisLength(double drawdim_axislength)
	{
		draw_dim.Y()=drawdim_axislength;
	}


protected:
	TArray<Vector3D> loccoords;

	Matrix dpdq; //temporary element, no need to copy???
	Matrix hmat;
	Vector hvec;
};


#endif
