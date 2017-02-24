//#**************************************************************
//#
//# filename:             Mass1D.h
//#
//# project:              
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
 
#ifndef Mass1D__H
#define Mass1D__H

#include "element.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Mass1D	Mass1D	Mass1D	Mass1D	Mass1D	Mass1D	Mass1D	Mass1D	Mass1D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Mass1D: public Element //$EDC$[beginclass,classname=Mass1D,parentclassname=Element,addelementtype=TAEBody+TAE1D,addelementtypename=Mass1D,texdescription="A point mass in one dimensions with 1 position coordinate. The computation of the dynamics of the point mass is extremely simple. The Mass1D can be used for a lot of applications which can be represented by the same type of equations. If you interpret the 'mass' to be 'moment of inertia' and the 'position' to be 'angle', then you can realize a 1D rotatory element as well.",
//texdescriptionDOF="1 degree of freedom: the position in x-direction",
//texdescriptionLimitations="The mass has no rotations, thus external moments can not be applied. The transformation of local to global coordinates is based on a translation, e.g. the global mass position is added to the local coordinates.",
//texdescriptionEquations="
//\begin{equation}
//m \ddot{x} = F	
//\end{equation}
// with the mass $m$ and the force $F$.",
//texdescriptionGeometry="The global position $p_{glob}$ of a local point p is computed as
//\begin{equation}
// p_{glob} = p_0 + A \left( \left( \begin{array}{c}x \\ 0 \\ 0 \end{array} \right) + p \right)
//\end{equation}
//with the reference\_position $p_0$ and the rotation\_matrix $A$.",
//example="Mass1D.txt",
//figure="Mass1D, Mass1D"]
{
public:
	Mass1D(MBS* mbsi):Element(mbsi) {ElementDefaultConstructorInitialization();};
	Mass1D(const Mass1D& e):Element(e.mbs) {CopyFrom(e);};

	//this function assigns default values to the element variables
	virtual void ElementDefaultConstructorInitialization();		

	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute

	virtual Element* GetCopy()
	{
		Element* ec = new Mass1D(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Element::CopyFrom(e);
		const Mass1D& ce = (const Mass1D&)e;
		drawres = ce.drawres;
		pref3D = ce.pref3D;  
		rotref3D = ce.rotref3D;
		radius = ce.radius;
	}

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//for load/save/edit:
	virtual const char* GetElementSpec() const {return "Mass1D";}
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	virtual void Initialize() 
	{
		Element::Initialize();
	};
	
	virtual void EvalM(Matrix& m, double t) 
	{
		m(1,1) = mass;
	}; 
	virtual void EvalF2(Vector& f, double t) 
	{
		Element::EvalF2(f,t);

		if (GetMassDamping() != 0)
		{
			f(1) -= GetMassDamping()*mass*XGP(1);
		}
	}; 

	virtual int SOS() const {return 1;}; //size of second order equations, len(u)
	virtual int Dim() const {return 1;} //default value
	virtual int IsRigid() const {return 1;} //default value

	virtual const double& Rho() const {assert(0); double dummy = -1.; return dummy; } 
	virtual double& Rho() {assert(0); double dummy = -1.; return dummy;} 

	// ............ 1D-access functions ............
	virtual double GetRefPos1D() const {return XG(1);};
	virtual double GetRefPos1DD() const 
	{
		if (GetMBS()->GetIOption(151) && IsRigid() && GetMBS()->GetDOption(105) != 1.)
		{
			double p = XGD(1);
			double p0 = x_init(1);
			return (p - p0)*GetMBS()->GetDOption(105) + p0;
		}
		else
		{
			return XGD(1);
		}
	};
	virtual double GetRefVel1D() const {return XGP(1);};
	virtual double GetRefVel1DD() const {return XGPD(1);};

	virtual double GetPos1D(const double& p_loc) const
  {
		return p_loc+GetRefPos1D();
	};
	virtual double GetPos1DD(const double& p_loc) const
  {
		return p_loc+GetRefPos1DD();
	};
	virtual double GetVel1D(const double& p_loc) const
  {
		return p_loc+GetRefVel1D();
	};
	virtual double GetVel1DD(const double& p_loc) const
  {
		return p_loc+GetRefVel1DD();
	};

	// ............ 3D-access functions ............
	virtual Vector3D GetRefPos() const
	{
		double p = GetRefPos1D();
		return ToP3D(p);
	}
	virtual Vector3D GetRefPosD() const
	{
		double p = GetRefPos1DD();
		return ToP3D(p);
	}
	virtual Vector3D GetPos(const Vector3D& p_loc) const
	{
		double p = GetPos1D(p_loc.X());
		return ToP3D(p);
	}
	virtual Vector3D GetPosD(const Vector3D& p_loc) const
	{
		double p = GetPos1DD(p_loc.X());
		return ToP3D(p);
	}
	virtual Vector3D GetVel(const Vector3D& p_loc) const
	{
		double p = GetVel1D(p_loc.X());
		return ToV3D(p);
	}

	virtual Matrix3D GetRotMatrix() const
	{
		Matrix3D rot(1.);
		return rot;
	}
	virtual Matrix3D GetRotMatrixD() const {return GetRotMatrix();};

	virtual Matrix3D GetRotMatrixP() const
	{
		Matrix3D rot(0.);
		return rot;
	}
	virtual Matrix3D GetRotMatrixPD() const {return GetRotMatrixP();};

	// transform points and velocities
	//set offset for 3D drawing and computation
	virtual void SetPRef3D(const Vector3D& p) {pref3D = p;}
	//set rotation for 3D transformation
	virtual void SetRotRef3D(const Matrix3D& rot) {rotref3D = rot;}
	virtual Vector3D ToP3D(const double& p) const 
	{
		return Vector3D(p,0,0)*rotref3D+pref3D;
	}
	virtual Vector3D ToV3D(const double& v) const 
	{
		return Vector3D(v,0,0)*rotref3D;
	}

	// ............ for loads:
	virtual void ApplyDprefdq(Vector& f, const double& x) 
	{
		f(1)=x;
	};
	virtual void GetDuxDq(Vector& dudq) 
	{
		dudq.SetAll(0);
		dudq(1) = 1;
	};
	virtual void GetDuyDq(Vector& dudq) {dudq.SetAll(0);};
	virtual void GetDuzDq(Vector& dudq) {dudq.SetAll(0);};
	// the rest is not done yet
	// ............ end for loads:

	virtual void DrawElement();
	virtual Box3D GetElementBox() const
	{
		//return Box3D(ToP3D(GetPos1D(-0.5*radius)),ToP3D(GetPos1D(0.5*radius)));
		return Box3D(Vector3D(GetPos1D(-radius))+pref3D,Vector3D(GetPos1D(radius))+pref3D); 		//$ DR 2013-09-23
	}

	virtual Box3D GetElementBoxD() const
	{
		//return Box3D(ToP3D(GetPos1DD(-0.5*radius)),ToP3D(GetPos1DD(0.5*radius)));
		return Box3D(Vector3D(GetPos1DD(-radius))+pref3D,Vector3D(GetPos1DD(radius))+pref3D); 		//$ DR 2013-09-23
	}

	
	//for visualization and sensoring of displacements and velocities:
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const double & local_position, bool flagD);
	virtual void GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables);

protected:
	int drawres; //$EDC$[varaccess,EDCvarname="drawing_tiling",EDCfolder="Graphics",tooltiptext="tiling of circle/sphere to represent Mass1D;
	//the drawing_tiling should be set small in order to improve efficiency,
	//but large for nice graphical represenations"]
	
  double radius;			 //$EDC$[varaccess,EDCvarname="radius",EDCfolder="Graphics",tooltiptext="drawing radius of mass"]
	//EDC Vector x_init; //$EDC$[varaccess,EDCvarname="initial_position",EDCfolder="Initialization",vecstart=1,vecend=1,tooltiptext="initial values for position [x]"]
	//EDC Vector x_init; //$EDC$[varaccess,EDCvarname="initial_velocity",EDCfolder="Initialization",vecstart=2,vecend=2,tooltiptext="initial values for velocity [v]"]

	//EDC double mass; //$EDC$[varaccess,EDCvarname="mass",EDCfolder="Physics",tooltiptext="total mass of point mass"] 

	Vector3D pref3D;   //$EDC$[varaccess,EDCvarname="reference_position",EDCfolder="Graphics",tooltiptext="Reference point for transformation of 1D objects to 3D; p = [X, Y, Z]"]
	Matrix3D rotref3D; //$EDC$[varaccess,EDCvarname="rotation_matrix",EDCfolder="Graphics",tooltiptext="Rotation matrix for transformation of 1D objects to 3D"]


}; //$EDC$[endclass,Mass1D]




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Rotor1D Rotor1D	Rotor1D	Rotor1D	Rotor1D	Rotor1D	Rotor1D	Rotor1D	Rotor1D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Rotor1D: public Mass1D //$EDC$[beginclass,classname=Rotor1D,parentclassname=Mass1D,addelementtype=TAEBody+TAE1D,addelementtypename=Rotor1D,texdescription="A rotor with 1 degree of freedom (the rotation). Mathematically implemented like Mass1D but different geometric representation.",
//texdescriptionDOF="1 degree of freedom: the rotation",
//texdescriptionEquations="
//\begin{equation}
//I \ddot{\varphi} = M	
//\end{equation}
//with the moment of inertia $I$ and the torque $M$.",
//texdescriptionGeometry="The global position $p_{glob}$ of a local point p is computed as
//\begin{equation}
// p_{glob} = p_0 + A_0 A p 
//\end{equation}
//with the reference\_position $p_0$, the constant rotation\_matrix $A_0$ and the non-constant rotation matrix 
//\begin{equation}
// A = \left( \begin{array}{ccc} 1 & 0 & 0 \\ 0  & \cos \varphi & -\sin \varphi \\ 0 &  \sin \varphi & \cos \varphi \end{array}\right)
//\end{equation}",
//example="Rotor1D.txt",
//figure="Rotor1D, Rotor1D is represented as rotating disc."]
{
public:
	Rotor1D(MBS* mbsi):Mass1D(mbsi) {ElementDefaultConstructorInitialization();};
	Rotor1D(const Rotor1D& e):Mass1D(e.mbs) {CopyFrom(e);};

	//this function assigns default values to the element variables
	virtual void ElementDefaultConstructorInitialization();		

	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute

	virtual Element* GetCopy()
	{
		Element* ec = new Rotor1D(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Mass1D::CopyFrom(e);
		const Rotor1D& ce = (const Rotor1D&)e;
		length = ce.length;
	}

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//for load/save/edit:
	virtual const char* GetElementSpec() const {return "Rotor1D";}
	//virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	//virtual int SetElementData(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// ............ 1D-access functions ............
	// all implemented in Mass1D

	// ............ 3D-access functions ............
	virtual Vector3D GetRefPos() const	{ return pref3D;	}
	virtual Vector3D GetRefPosD() const	{	return pref3D;	}
	virtual Vector3D GetPos(const Vector3D& p_loc) const
	{
		return pref3D + GetRotMatrix()*p_loc;
	}
	virtual Vector3D GetPosD(const Vector3D& p_loc) const
	{
		return pref3D + GetRotMatrixD()*p_loc;
	}
	virtual Vector3D GetVel(const Vector3D& p_loc) const
	{
		return GetRotMatrixP()*p_loc;
	}
	virtual Matrix3D GetRotMatrix() const
	{
		Matrix3D rot = RotMatrix1(GetRefPos1D());	//rotation around x-axis
		return rotref3D*rot;
	}
	virtual Matrix3D GetRotMatrixD() const 
	{
		Matrix3D rot = RotMatrix1(GetRefPos1DD());	//rotation around x-axis
		return rotref3D*rot;
	};

	virtual Matrix3D GetRotMatrixP() const
	{
		double p = GetRefPos1D();
		Matrix3D Ap(0.,0.,0.,
		0.,-sin(p),-cos(p),
		0.,cos(p), -sin(p));
		return GetRefVel1D()*Ap;
	}
	virtual Matrix3D GetRotMatrixPD() const 
	{
		double p = GetRefPos1DD();
		Matrix3D Ap(0.,0.,0.,
		0.,-sin(p),-cos(p),
		0.,cos(p), -sin(p));
		return GetRefVel1DD()*Ap;
	};

	virtual void DrawElement();
	virtual Box3D GetElementBox() const
	{
		return GetElementBoxD();
	}

	virtual Box3D GetElementBoxD() const
	{
		//return Box3D(ToP3D(GetPos1DD(-0.5*radius)),ToP3D(GetPos1DD(0.5*radius)));
		return Box3D(pref3D-Vector3D(radius+length)*1.1,pref3D+Vector3D(radius+length)*1.1);
	}
	
	//for visualization and sensoring of displacements and velocities:
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const double & local_position, bool flagD);
	virtual void GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables);

protected:
	//EDC Vector x_init; //$EDC$[varaccess,EDCvarname="initial_position",EDCfolder="Initialization",remove]
	//EDC Vector x_init; //$EDC$[varaccess,EDCvarname="initial_velocity",EDCfolder="Initialization",remove]

	//EDC Vector x_init; //$EDC$[varaccess,EDCvarname="initial_rotation",EDCfolder="Initialization",vecstart=1,vecend=1,tooltiptext="initial value for rotation"]
	//EDC Vector x_init; //$EDC$[varaccess,EDCvarname="initial_angular_velocity",EDCfolder="Initialization",vecstart=2,vecend=2,tooltiptext="initial value for angular velocity"]

	//EDC double mass; //$EDC$[varaccess,EDCvarname="mass",EDCfolder="Physics",remove] 
	//EDC double mass; //$EDC$[varaccess,EDCvarname="moment_of_inertia",EDCfolder="Physics",tooltiptext="mass moment of inertia in kg*m*m"] 

	//EDC int drawres;				 //$EDC$[varaccess,EDCvarname="drawing_tiling",EDCfolder="Graphics",remove]
	//EDC double radius;			 //$EDC$[varaccess,EDCvarname="radius",EDCfolder="Graphics",remove]
	//EDC double radius;			 //$EDC$[varaccess,EDCvarname="radius",EDCfolder="Graphics",tooltiptext="radius of rotor"]
	double length;						 //$EDC$[varaccess,EDCvarname="length",EDCfolder="Graphics",tooltiptext="length of rotor"]

}; //$EDC$[endclass,Rotor1D]

#endif