//#**************************************************************
//#
//# filename:             KinematicPairsDeprecated.h
//#
//# author:               Gerstmayr Johannes, Reischl Daniel
//#
//# generated:						29. August 2013

//# description:          This class contains only DEPRECATED kinematic pairs. This file exists only for backward compatibilty. DO NOT USE THE CLASSES FOR NEW MODELS.
//#												Kinematic pairs can be divided in lower (surface contact) and higher (point or line contact) kinematic pairs. 
//#												In this file the lower kinematic pairs and the rigid joints (all d.o.f. are constrained) are provided.
//#												Other implemented constraints are in SpecialConstraints.h
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
 

#ifndef KINEMATICPAIRSDEPRECATED__H
#define KINEMATICPAIRSDEPRECATED__H

#include "constraint.h"


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Cylindrical	Cylindrical	Cylindrical	Cylindrical	Cylindrical	Cylindrical			
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Cylindrical: CylindricalJointOLD 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: CylindricalJointOLD (Constraint)
// short description: 
// available formulations: Lagrange multiplier
// type of constraint equation: 4 nonlinear constraint equations
// index formulation available: 
// ground constraint: yes
// element-element constraint: yes
// development status: 
// long description:	relative rotation with respect to rotation axis and movement along rotation axis is possible
//								
// class variables:
//  
//**end(ued)**
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// drawing options
// draw_dim.X()==diameter
// draw_dim.Y()==axis length
// draw_dim.Z()==draw resolution
class CylindricalJointOLD: public Constraint
{
public:
	CylindricalJointOLD(MBS* mbsi):Constraint(mbsi), loccoords(), dpdq(), hmat(), hvec()
	{
	};
	CylindricalJointOLD(MBS* mbsi, int en1, const Vector3D& locp, const Vector3D& pglob,  const Vector3D& globrot,
		const Vector3D& ddim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq(), hmat(), hvec() 
	{	
		x_init = Vector(SS());
		GetCol() = coli;
		//draw_dim = ddim;
		SetDrawSizeScalar(0.5*ddim.X());
		SetDrawSizeAxisLength(ddim.Y());
		SetDrawSizeResolution((int)(ddim.Z()));

		p_global = pglob;

		Vector3D vrot;
		vrot = globrot;
		vrot.Normalize();

		AddElement(en1);
		loccoords.Add(locp);
		loccoords.Add(vrot);
		loccoords.Add(Vector3D(0)); //dummy for later on
		loccoords.Add(Vector3D(0));

		elementname = GetElementSpec();
	};
	CylindricalJointOLD(MBS* mbsi, int en1, int en2, 
		const Vector3D& lp1, const Vector3D& lp2, const Vector3D& globrot,
		const Vector3D& ddim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq(), hmat(), hvec() 
	{
		x_init = Vector(SS());
		GetCol() = coli;
		//draw_dim = ddim;
		SetDrawSizeScalar(0.5*ddim.X());
		SetDrawSizeAxisLength(ddim.Y());
		SetDrawSizeResolution((int)(ddim.Z()));

		Vector3D vrot;
		vrot = globrot;
		vrot.Normalize();

		AddElement(en1);
		AddElement(en2);
		loccoords.Add(lp1);
		loccoords.Add(lp2);
		loccoords.Add(vrot);
		loccoords.Add(Vector3D(0)); //ln2
		loccoords.Add(Vector3D(0));	//lt2
		loccoords.Add(Vector3D(0)); //lrot2

		elementname = GetElementSpec();
	};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new CylindricalJointOLD(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const CylindricalJointOLD& ce = (const CylindricalJointOLD&)e;
		loccoords = ce.loccoords;
		p_global = ce.p_global;
		dpdq = ce.dpdq;
		hmat = ce.hmat;
		hvec = ce.hvec;
	}

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	virtual const char* GetElementSpec() const {return "CylindricalJointOLD";} //this name names MUST coincide with the type names defined in MBS::MBS() !!!!
	virtual void Initialize();

	virtual void EvalG(Vector& f, double t);

	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size
	virtual int IS() const {return 4;};

	virtual int Dim() const {return 3;}

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;

	virtual void DrawElement();

	double CylindricalJointOLD::GetDrawSizeAxisLength()
	{
		if(draw_dim.Y()==-1)
			return GetMBS()->GetDOption(172);
		else
			return draw_dim.Y();
	}
	void CylindricalJointOLD::SetDrawSizeAxisLength(double drawdim_axislength)
	{
		draw_dim.Y()=drawdim_axislength;
	}

	int CylindricalJointOLD::GetDrawSizeResolution()
	{
		if(draw_dim.Z()==-1)
			return (int)(GetMBS()->GetDOption(173));
		else
			return (int)draw_dim.Z();
	}

	void CylindricalJointOLD::SetDrawSizeResolution(int drawdim_resolution)
	{
		draw_dim.Z()= (double)drawdim_resolution;
	}

protected:
	TArray<Vector3D> loccoords;
	Vector3D p_global;

	Matrix dpdq; //temporary element, no need to copy???
	Matrix hmat;
	Vector hvec;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Translatory	Translatory	Translatory	Translatory	Translatory	Translatory				
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Translatory: PrismaticJointOLD
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: PrismaticJointOLD (Constraint)
// short description: 
// available formulations: Lagrange multiplier
// type of constraint equation: nonlinear constraint equations
// index formulation available: 
// ground constraint: yed
// element-element constraint: yes
// development status: 
// long description:	
//										
// class variables:
//  
//**end(ued)**
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//Constraint for a position measured relative within one or two bodies

class PrismaticJointOLD: public Constraint
{
public:
	PrismaticJointOLD(MBS* mbsi):Constraint(mbsi), loccoords(), dpdq() 
	{
	};
	//relative motion of body and ground along axis lp1-gp2, lp1 and gp2 must never be equal!!
	//relative rotation around axis lp1-lp2 is constrained
	PrismaticJointOLD(MBS* mbsi, int en1, 
		const Vector3D& lp1, const Vector3D& gp2,
		const Vector3D& ddim, const Vector3D& coli);

	//relative motion of body and ground along local body1 axis laxis1, lp1 and gp2 may be equal!!
	//relative rotation around axis lp1-lp2 is constrained
	PrismaticJointOLD(MBS* mbsi, int en1, 
		const Vector3D& lp1, const Vector3D& gp2, const Vector3D& laxis1,
		const Vector3D& ddim, const Vector3D& coli);

	//relative motion of two bodies along axis lp1-lp2, lp1 and lp2 must never be equal!!
	//relative rotation around axis lp1-lp2 is constrained
	PrismaticJointOLD(MBS* mbsi, int en1, int en2, 
		const Vector3D& lp1, const Vector3D& lp2,
		const Vector3D& ddim, const Vector3D& coli);

	//relative motion of two bodies along local (body1) axis laxis1, lp1 and lp2 may be equal!!
	//relative rotation around axis lp1-lp2 is constrained
	PrismaticJointOLD(MBS* mbsi, int en1, int en2, 
		const Vector3D& lp1, const Vector3D& lp2, const Vector3D& laxis1, 
		const Vector3D& ddim, const Vector3D& coli);

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new PrismaticJointOLD(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const PrismaticJointOLD& ce = (const PrismaticJointOLD&)e;
		loccoords = ce.loccoords;
		dpdq = ce.dpdq;
		hmat = ce.hmat;
		hvec = ce.hvec;
	}

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer
	virtual const char* GetElementSpec() const {return "PrismaticJointOLD";} //this name names MUST coincide with the type names defined in MBS::MBS() !!!!

	virtual void Initialize();

	virtual void EvalG(Vector& f, double t);

	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size
	virtual int IS() const {return 5;};

	virtual int Dim() const {return 3;}

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;

	virtual void DrawElement();

protected:
	TArray<Vector3D> loccoords;

	Matrix dpdq; //temporary element, no need to copy???
	Matrix hmat;
	Vector hvec;
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Rotary		Rotary		Rotary		Rotary		Rotary		Rotary		Rotary		
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Rotary: RevoluteJointOLD
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: RevoluteJointOLD (Constraint)
// short description: 
// available formulations: Lagrange multiplier
// type of constraint equation: nonlinear constraint equations
// index formulation available: 
// ground constraint: yes
// element-element constraint: yes
// development status: 
// long description:	all translatory and 2 rotatory d.o.f. are constrained, only rotation about 1 axis possible
//								
// class variables:
//  
//**end(ued)**
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// drawing options
// draw_dim.X()==diameter
// draw_dim.Y()==axis length
// draw_dim.Z()==draw resolution
class RevoluteJointOLD: public Constraint
{
public:
	RevoluteJointOLD(MBS* mbsi):Constraint(mbsi), loccoords(), dpdq(), hmat(), hvec()
	{
	};
	RevoluteJointOLD(MBS* mbsi, int en1, const Vector3D& locp, const Vector3D& pglob,  const Vector3D& globrot,
		const Vector3D& ddim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq(), hmat(), hvec() 
	{	
		x_init = Vector(SS());
		GetCol() = coli;
		//draw_dim = ddim;
		SetDrawSizeScalar(ddim.X());
		SetDrawSizeAxisLength(ddim.Y());
		SetDrawSizeResolution((int)ddim.Z());

		p_global = pglob;

		Vector3D vrot;
		vrot = globrot;
		vrot.Normalize();

		AddElement(en1);
		loccoords.Add(locp);
		loccoords.Add(vrot);
		loccoords.Add(Vector3D(0)); //dummy
		loccoords.Add(Vector3D(0));

		elementname = GetElementSpec();
	};
	RevoluteJointOLD(MBS* mbsi, int en1, int en2, 
		const Vector3D& lp1, const Vector3D& lp2, const Vector3D& globrot,
		const Vector3D& ddim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq(), hmat(), hvec() 
	{
		x_init = Vector(SS());
		GetCol() = coli;
		//draw_dim = ddim;
		SetDrawSizeScalar(ddim.X());
		SetDrawSizeAxisLength(ddim.Y());
		SetDrawSizeResolution((int)ddim.Z());

		Vector3D vrot;
		vrot = globrot;
		vrot.Normalize();

		AddElement(en1);
		AddElement(en2);
		loccoords.Add(lp1);
		loccoords.Add(lp2);
		loccoords.Add(vrot);
		loccoords.Add(Vector3D(0)); //ln2
		loccoords.Add(Vector3D(0));	//lt2
		loccoords.Add(Vector3D(0)); //lrot2

		elementname = GetElementSpec();
	};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new RevoluteJointOLD(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const RevoluteJointOLD& ce = (const RevoluteJointOLD&)e;
		loccoords = ce.loccoords;
		p_global = ce.p_global;
		dpdq = ce.dpdq;
		hmat = ce.hmat;
		hvec = ce.hvec;
	}

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	virtual const char* GetElementSpec() const {return "RevoluteJointOLD";} //this name names MUST coincide with the type names defined in MBS::MBS() !!!!
	virtual void Initialize();

	virtual void EvalG(Vector& f, double t);

	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size
	virtual int IS() const {return 5;};

	virtual int Dim() const {return 3;}

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;

	virtual void DrawElement();

	double RevoluteJointOLD::GetDrawSizeAxisLength()
	{
		if(draw_dim.Y()==-1)
			return GetMBS()->GetDOption(172);
		else
			return draw_dim.Y();
	}
	void RevoluteJointOLD::SetDrawSizeAxisLength(double drawdim_axislength)
	{
		draw_dim.Y()=drawdim_axislength;
	}

	int RevoluteJointOLD::GetDrawSizeResolution()
	{
		if(draw_dim.Z()==-1)
			return (int)(GetMBS()->GetDOption(173));
		else
			return (int)draw_dim.Z();
	}

	void RevoluteJointOLD::SetDrawSizeResolution(int drawdim_resolution)
	{
		draw_dim.Z()= (double)drawdim_resolution;
	}

protected:
	TArray<Vector3D> loccoords;
	Vector3D p_global;

	Matrix dpdq; //temporary element, no need to copy???
	Matrix hmat;
	Vector hvec;
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Rigid	Rigid	Rigid	Rigid	Rigid	Rigid	Rigid	Rigid	Rigid	Rigid	Rigid	Rigid	
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Rigid: RigidJointOLD
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: RigidJointOLD (Constraint)
// short description: 
// available formulations: Lagrange multiplier
// type of constraint equation: nonlinear constraint equations
// index formulation available: 
// ground constraint: yes
// element-element constraint: yes
// development status: 
// long description:	Constraint for constraining all d.o.f. 
//								
// class variables:
//  
//**end(ued)**
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class RigidJointOLD: public Constraint //!EDC$[beginclass,classname=RigidJointOLD,parentclassname=Constraint,addelementtype=TAEconstraint,addelementtypename=RigidJointOLD,texdescription="
//The RigidJointOLD constrains the position and relative angles of an element at a specified position. If only one element is specified, a ground joint is realized. 
//"]
{
public:
	RigidJointOLD(MBS* mbsi):Constraint(mbsi), loccoords(), dpdq() 
	{
		ElementDefaultConstructorInitialization();
	};
	//Position and relative angle of body is constrained at local point lp1 and gp2
	RigidJointOLD(MBS* mbsi, int en1, 
		const Vector3D& lp1, const Vector3D& gp2,
		const Vector3D& ddim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq(), hmat(), hvec() 
	{	
		ElementDefaultConstructorInitialization();
		SetRigidJointOLDElementToGround(en1, lp1, gp2);
		//x_init = Vector(SS());	//$ DR 2012-10: moved to ElementDefaultConstructorInitialization
		GetCol() = coli;
		SetDrawSizeScalar(ddim.X());
		//AddElement(en1); //$ DR 2012-10: moved to ElementDefaultConstructorInitialization and SetFunction
		//loccoords.Add(lp1); //$ DR 2012-10: moved to ElementDefaultConstructorInitialization and SetFunction
		//loccoords.Add(gp2); //$ DR 2012-10: moved to ElementDefaultConstructorInitialization and SetFunction
		//loccoords.Add(Vector3D(0)); //$ DR 2012-10: moved to ElementDefaultConstructorInitialization
		//loccoords.Add(Vector3D(0)); //$ DR 2012-10: moved to ElementDefaultConstructorInitialization
		//loccoords.Add(Vector3D(0)); //$ DR 2012-10: moved to ElementDefaultConstructorInitialization
		//loccoords.Add(Vector3D(0)); //$ DR 2012-10: moved to ElementDefaultConstructorInitialization
		//loccoords.Add(Vector3D(0)); //$ DR 2012-10: moved to ElementDefaultConstructorInitialization
		//elementname = GetElementSpec(); //$ DR 2012-10: moved to ElementDefaultConstructorInitialization
	};
	//Position and relative angle of both bodies is constrained in local points lp1 and lp2
	RigidJointOLD(MBS* mbsi, int en1, int en2, 
		const Vector3D& lp1, const Vector3D& lp2,
		const Vector3D& ddim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq(), hmat(), hvec() 
	{	
		ElementDefaultConstructorInitialization();
		SetRigidJointOLDElementToElement(en1, en2, lp1, lp2);
		//x_init = Vector(SS()); //$ DR 2012-10: moved to ElementDefaultConstructorInitialization
		GetCol() = coli;
		SetDrawSizeScalar(ddim.X());
		//AddElement(en1); //$ DR 2012-10: moved to ElementDefaultConstructorInitialization and SetFunction
		//AddElement(en2); //$ DR 2012-10: moved to ElementDefaultConstructorInitialization and SetFunction
		//loccoords.Add(lp1); //$ DR 2012-10: moved to ElementDefaultConstructorInitialization and SetFunction
		//loccoords.Add(lp2); //$ DR 2012-10: moved to ElementDefaultConstructorInitialization and SetFunction
		//loccoords.Add(Vector3D(0)); //$ DR 2012-10: moved to ElementDefaultConstructorInitialization
		//loccoords.Add(Vector3D(0)); //$ DR 2012-10: moved to ElementDefaultConstructorInitialization
		//loccoords.Add(Vector3D(0)); //$ DR 2012-10: moved to ElementDefaultConstructorInitialization
		//loccoords.Add(Vector3D(0)); //$ DR 2012-10: moved to ElementDefaultConstructorInitialization
		//loccoords.Add(Vector3D(0)); //$ DR 2012-10: moved to ElementDefaultConstructorInitialization
		//elementname = GetElementSpec(); //$ DR 2012-10: moved to ElementDefaultConstructorInitialization
		
	};

	virtual void ElementDefaultConstructorInitialization();

	virtual void SetRigidJointOLDElementToElement(int en1, int en2, const Vector3D& lp1, const Vector3D& lp2)
	{
		elements.SetLen(2);
		elements(1)=en1;
		elements(2)=en2;
		loccoords(1)=lp1;
		loccoords(2)=lp2;
	}

	virtual void SetRigidJointOLDElementToGround(int en1, const Vector3D& lp1, const Vector3D& gp2)
	{
		elements.SetLen(2);
		elements(1)=en1;
		elements(2)=0;
		loccoords(1)=lp1;
		loccoords(2)=gp2;
	}

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new RigidJointOLD(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const RigidJointOLD& ce = (const RigidJointOLD&)e;
		loccoords = ce.loccoords;
		dpdq = ce.dpdq;
		hmat = ce.hmat;
		hvec = ce.hvec;
	}
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer
	
	//virtual void GetElementDataAuto(ElementDataContainer& edc); 		//automatically generated function from EDC_converter
	//virtual int SetElementDataAuto(ElementDataContainer& edc); //automatically generated function from EDC_converter
	//virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	//virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
 // virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); // returns all available directly (ReadWrite-) accessable variables

	virtual const char* GetElementSpec() const {return "RigidJointOLD";} //this name names MUST coincide with the type names defined in MBS::MBS() !!!!

	virtual int NE_nodouble() const //$ DR 2012-11-02: returns 1 if ground joint and 2 otherwise
	{
		int cnt=0;
		for(int i=1; i<=elements.Length(); i++)
		{
			if(elements(i)) {cnt++;}
		}
		return cnt;
	}
	virtual int NE() const {return NE_nodouble();}		

	virtual void Initialize();

	virtual void EvalG(Vector& f, double t);

	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size
	virtual int IS() const {return 6;};

	virtual int Dim() const {return 3;}

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;

	virtual void DrawElement();

protected:
	TArray<Vector3D> loccoords;

	Matrix dpdq; //temporary element, no need to copy???
	Matrix hmat;
	Vector hvec;
	//EDC double draw_dim(1);						//!EDC$[varaccess,EDCvarname="draw_size",EDCfolder="Graphics",tooltiptext="drawing dimensions of constraint. If set to -1, than global_draw_scalar_size is used."]
	//EDC int elements(1);							//!EDC$[varaccess,EDCvarname="element_number",EDCfolder="Position1",tooltiptext="Number of constrained element"]
	//EDC int elements(2);							//!EDC$[varaccess,EDCvarname="element_number",EDCfolder="Position2",tooltiptext="Number of constrained element"]
	//EDC Vector3D loccoords(1);				//!EDC$[varaccess,EDCvarname="position",EDCfolder="Position1",tooltiptext="local position"]
	//EDC Vector3D loccoords(2);				//!EDC$[varaccess,EDCvarname="position",EDCfolder="Position2",tooltiptext="local or global (if element_number == 0) position."]

	//// remove entries from EDC //$ DR 2012-11-02
	//EDC	int use_penalty_formulation;			//!EDC$[varaccess,remove,EDCvarname="use_penalty_formulation",EDCfolder="Physics"]
	//EDC	double spring_stiffness;   				//!EDC$[varaccess,remove,EDCvarname="spring_stiffness",EDCfolder="Physics.Penalty"]
	//EDC	int use_local_coordinate_system;	//!EDC$[varaccess,remove,EDCvarname="use_local_coordinate_system",EDCfolder="Geometry"]

};//!EDC$[endclass,RigidJointOLD]

#endif