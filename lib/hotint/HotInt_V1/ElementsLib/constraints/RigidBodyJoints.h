//#**************************************************************
//#
//# filename:             RigidBodyJoints.h
//#
//# author:               Saxinger Martin
//#
//# generated:						19. December 2012
//# description:          
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
 
#ifndef RIGIDBODYJOINTS__H
#define RIGIDBODYJOINTS__H

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Rigid: BaseBodyJoint
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: BaseBodyJoint (Constraint)
// short description: 
// available formulations: Lagrange multiplier, Penalty Formulation
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

#include "element.h"
#include "constraint.h"
#include "KinematicPairs.h"

// SM is responsible for this constraint
class BaseBodyJoint: public BasePointJoint  //$EDC$[beginclass,classname=BaseBodyJoint,parentclassname=BasePointJoint,addelementtype=TAEconstraint,addelementtypename=GenericBodyJoint,figure="genericBodyJoint",texdescription="
//The GenericBodyJoint constrains two elements at a local position each. If only one element is specified (second element 0), a ground GenericBodyJoint is realized. A penalty and lagrange formulation is available.", 
//texdescriptionLimitations="It is strongly recommended to prefer the Lagrangian method for free rotation instead of penalty formulation to avoid simulation problems.",
//modus="{element to ground}{Position2.element\_number AND Position2.node\_number have to be equal to 0}",
//modus="{element to element}{Position2.element\_number and/or Position2.node\_number must not be equal to 0}",
//modus="{Lagrange}{Physics.use\_penalty\_formulation must be set to 0. \newline 
//Set the vector of constrained directions in Physics.Lagrange.constrained\_directions ($\left[x,y,z\right]$, 1 = constrained, 0 = free). The directions are w.r.t the local body 1 joint coordinate system. \newline 
//Set the vector of constrained rotations in Physics.Lagrange.constrained\_rotations ($\left[\phi_x,\phi_y,\phi_z\right]$, 1 = constrained, 0 = free). The rotations are about the axes of local body 1 joint coordinate system.}",
//modus="{Penalty}{Physics.use\_penalty\_formulation must be set to 1. \newline 
//In Physics.Penalty.stiffness\_matrix and Physics.Penalty.damping\_matrix all parameters for translational stiffness and damping w.r.t. local body 1 coordinate system can be set. \newline 
//In Physics.Penalty.stiffness\_matrix\_rotation and Physics.Penalty.damping\_matrix\_rotation all parameters for rotational stiffness and damping w.r.t. local body 1 coordinate system can be set.}",example="GenericBodyJointShort.txt"]
{

public:

	BaseBodyJoint(MBS* mbsi):BasePointJoint(mbsi)
	{	
		ElementDefaultConstructorInitialization();
	};

	BaseBodyJoint(const BaseBodyJoint& ct): BasePointJoint(ct.mbs)
	{
		ElementDefaultConstructorInitialization();
		CopyFrom(ct);
	};

	~BaseBodyJoint()
	{
	};

	void SetBaseBodyJoint_LocalPos_to_LocalPos(int elem1, Vector3D lc1, int elem2, Vector3D lc2, IVector constrDir = IVector(3), IVector constrDirRot = IVector(3), Vector3D bryantAnglesI = Vector3D(0.), Matrix3D stiffness = Matrix3D(0), Matrix3D stiffnessRot = Matrix3D(0), Matrix3D damping = Matrix3D(0), Matrix3D dampingRot = Matrix3D(0));

	virtual void ElementDefaultConstructorInitialization();

	int IsNullMatrix(Matrix3D mat) const; // TODO $DR+MSax: verschieben

	virtual void Initialize();

	int PenaltyRotModus() const; // only for penalty formulation: returns the rotation axis: 0 = rigid constraint, 1 = rot about x-axis, 2 = rot about y-axis, 3 = rot about z-axis, 4 = free rotations

	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute

	virtual Element* GetCopy()
	{
		Element* ec = new BaseBodyJoint(*this);
		return ec;
	}

	virtual void CopyFrom(const Element& e);

	virtual const char* GetElementSpec() const {return "GenericBodyJoint";}

	// for constraints only, these are the necessary (!) access functions!	//$ MSax 2013-02-18
	virtual void GetNecessaryKinematicAccessFunctions(TArray<int> &KAF, int numberOfKinematicPair);

  virtual int GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //for special value sensor
	virtual int ReadSingleElementData(ReadWriteElementDataVariableType& RWdata); 		//for special value sensor
	virtual int WriteSingleElementData(const ReadWriteElementDataVariableType& RWdata);

	virtual int GetSignOfLagrangeParam(int dir);  // only for lagrange formulation

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		
	virtual int SetElementDataAuto(ElementDataContainer& edc);
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual int GetNumberOfConstrainedCoords() const;

	virtual void EvalG(Vector& f, double t);
	virtual void EvalF2(Vector& f, double t);
	virtual Vector3D ComputeForce(double t) const;
	virtual Vector3D ComputeMoment(double t) const;
	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	Matrix3D GetStiffnessMatrix() const {return penaltyStiffness;};
	void SetStiffnessMatrix(Matrix3D stiffness3i) 
	{
		penaltyStiffness = stiffness3i; 
		SetDampingMatrix(Matrix3D(0.));
	}
	Matrix3D GetDampingMatrix() const {return penaltyDamping;};
	void SetDampingMatrix(Matrix3D damping3i) {penaltyDamping = damping3i;};
	
	Matrix3D GetStiffnessMatrixRot() const {return penaltyStiffnessRot;};
	void SetStiffnessMatrixRot(Matrix3D stiffness3i) 
	{
		penaltyStiffnessRot = stiffness3i; 
		SetDampingMatrixRot(Matrix3D(0.));
	}
	Matrix3D GetDampingMatrixRot() const {return penaltyDampingRot;};
	void SetDampingMatrixRot(Matrix3D damping3i) {penaltyDampingRot = damping3i;};

	void SetConstrainedDirectionsRot(IVector directions) {dirRot = directions;}

	virtual Matrix3D GetRotMati() const; // rot matrix from joint coordinate system i to global coordinate system
	virtual Matrix3D GetRotMatj() const;
	virtual Matrix3D GetRotMatiD() const;
	virtual Matrix3D GetRotMatjD() const;

	Matrix3D GetRotMatBodyi() const {return GetRotMatrixBody(1);}; // rot matrix from body coordinate system i to global coordinate system
	Matrix3D GetRotMatBodyj() const {return GetRotMatrixBody(2);};

	Matrix3D GetRotMatrixBody(int i) const;
	Matrix3D GetRotMatrixBodyD(int i) const;

	virtual Matrix3D GetRotMatiP() const; // time derivative of rot matrix from joint coordinate system i to global coordinate system
	virtual Matrix3D GetRotMatjP() const;

	Matrix3D GetRotMatBodyiP() const {return GetRotMatrixBodyP(1);};
	Matrix3D GetRotMatBodyjP() const {return GetRotMatrixBodyP(2);};

	Matrix3D GetRotMatrixBodyP(int i) const;

	virtual void DrawElement();
	virtual double GetDrawSizeCone();

protected:

	// IVector dir; // already defined in BasePointJoint
	TArray<int> dirRot;	 //$EDC$[varaccess,EDCvarname="constrained_rotations",minval=0,maxval=1,EDCfolder="Physics.Lagrange",tooltiptext="[angle about x axis,angle about y axis,angle about z axis]...(1 = constrained, 0 = free), can be defined as local or global directions (see Geometry)"]

	Matrix3D JA0j;	// tranformation matrix from joint coordinate system to body coordinate system
	Matrix3D penaltyStiffness; //$EDC$[varaccess,EDCvarname="stiffness_matrix",EDCfolder="Physics.Penalty",tooltiptext="3x3 matrix with stiffness parameters"]
	Matrix3D penaltyDamping; //$EDC$[varaccess,EDCvarname="damping_matrix",EDCfolder="Physics.Penalty",tooltiptext="3x3 matrix with damping parameters"]
	Matrix3D penaltyStiffnessRot; //$EDC$[varaccess,EDCvarname="stiffness_matrix_rotation",EDCfolder="Physics.Penalty",tooltiptext="3x3 matrix with stiffness parameters for rotation"]
	Matrix3D penaltyDampingRot; //$EDC$[varaccess,EDCvarname="damping_matrix_rotation",EDCfolder="Physics.Penalty",tooltiptext="3x3 matrix with damping parameters for rotation"]

	//EDC	Matrix3D JA0i;		//$EDC$[varaccess,remove,EDCvarname="joint_local_frame",EDCfolder="Geometry"]
	//EDC	Matrix3D JA0i;		//$EDC$[varaccess,readonly,EDCvarname="joint_local_frame",EDCfolder="Geometry"]
	Vector3D bryantAngles;	//$EDC$[varaccess,EDCvarname="joint_local_frame_in_bryant_angles",EDCfolder="Geometry",tooltiptext="Prerotate joint coordinate system w.r.t. local coordinate system of body 1 [phi x, phi y, phi z]. Rot. sequence: JA0i=A(phi z)A(phi y)A(phi x)"]

	int freeRotPenalty; // rotation modus used for penalty formulation: 0 = rigid constraint, 1 = rot about x-axis, 2 = rot about y-axis, 3 = rot about z-axis, 4 = free rotations

	// remove unused data from class
	//EDC	int use_local_coordinate_system;	//$EDC$[varaccess,remove,EDCvarname="use_local_coordinate_system",EDCfolder="Geometry"]
	//EDC	int stiffness_in_joint_local_frame;	//$EDC$[varaccess,remove,EDCvarname="use_joint_local_frame",EDCfolder="Geometry"]
	//EDC	Vector3D spring_stiffness3;	//$EDC$[varaccess,remove,EDCvarname="spring_stiffness_vector",EDCfolder="Physics.Penalty"]
	//EDC	double damping_coeff;	//$EDC$[varaccess,remove,EDCvarname="damping",EDCfolder="Physics.Penalty"]
	//EDC	double spring_stiffness; //$EDC$[varaccess,remove,EDCvarname="spring_stiffness",EDCfolder="Physics.Penalty"]
	//EDC int nodes(1); //$EDC$[varaccess,remove,EDCvarname="node_number",EDCfolder="Position1"]
	//EDC int nodes(2);	//$EDC$[varaccess,remove,EDCvarname="node_number",EDCfolder="Position2"]
	//EDC double draw_dim(1);	//$EDC$[varaccess,remove,EDCvarname="draw_size",EDCfolder="Graphics"]

	int standard_joint_drawing; // flag for drawing mode; 1 ==> show constrained directions and rotations; 2 ==> draw constraint element
	double draw_cone_size; //$EDC$[varaccess,EDCvarname="cone_size",EDCfolder="Graphics",tooltiptext="cone size for standard joint drawing"]

	//EDC	Vector3D col; //$EDC$[varaccess,remove,EDCvarname="RGB_color",EDCfolder="Graphics"] //color for body
	//EDC	Vector3D col; //$EDC$[varaccess,EDCvarname="color_body1",EDCfolder="Graphics",tooltiptext="[red, green, blue] first color of constraint, range = 0..1, use default color:[-1,-1,-1]"] //color for constraint
	//EDC	Vector3D col_ext; //$EDC$[varaccess,EDCvarname="color_body2",EDCfolder="Graphics",tooltiptext="[red, green, blue] second color of constraint, range = 0..1, use default color:[-1,-1,-1]"] //color for constraint
private:
	Matrix tmpmat1;
	Matrix tmpmat2;
	Vector tmpvec1;
};//$EDC$[endclass,BaseBodyJoint]



// SM is responsible for this constraint
//$ SW 2013-08-29: changed class name from RevoluteJoint_ to RevoluteJoint
class RevoluteJoint: public BaseBodyJoint  //$EDC$[beginclass,classname=RevoluteJoint,parentclassname=BaseBodyJoint,addelementtype=TAEconstraint,addelementtypename=RevoluteJoint,figure="revoluteJoint",texdescription="
//The RevoluteJoint constrains all relative degrees of freedom between two bodies except the rotation about a local rotation axis. A penalty formulation exists, which replaces the exact lagrange constraint by a approximation with joint stiffness and damping.",
//example="RevoluteJointShort.txt"]
{

public:

	RevoluteJoint(MBS* mbsi):BaseBodyJoint(mbsi)
	{	
		ElementDefaultConstructorInitialization();
	};

	RevoluteJoint(const RevoluteJoint& ct): BaseBodyJoint(ct.mbs)
	{
		ElementDefaultConstructorInitialization();
		CopyFrom(ct);
	};

	~RevoluteJoint()
	{
	};

	// MSax: Set function not testet yet
	void SetRevoluteJointLocalPos_to_LocalPos(int elem1, Vector3D lc1, int elem2, Vector3D lc2, Vector3D rot_axisI, double stiffnessI = 0., double stiffnessRotI = 0., double dampingI = 0., double dampingRotI = 0.)
	{
		rotAxis = rot_axisI;
		SetPos1ToLocalCoord(elem1, lc1);
		SetPos2ToLocalCoord(elem2, lc2);
	}

	virtual void ElementDefaultConstructorInitialization();

	virtual void Initialize();

	virtual int CheckConsistency(mystr& errorstr);

	virtual Element* GetCopy()
	{
		Element* ec = new RevoluteJoint(*this);
		return ec;
	}

	virtual void CopyFrom(const Element& e);

	virtual const char* GetElementSpec() const {return "RevoluteJoint";}

  virtual int GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //for special value sensor

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		
	virtual int SetElementDataAuto(ElementDataContainer& edc);
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual void DrawElement();

	virtual double GetDrawSizeAxisLength();
	virtual double GetDrawSizeScalar();

	virtual void SetDrawSize(double axis_length, double diameter) {draw_axis_length = axis_length; draw_diameter = diameter;}

protected:
	//EDC	TArray<int> dir; //$EDC$[varaccess,remove,EDCvarname="constrained_directions",EDCfolder="Physics.Lagrange"]
	//EDC	TArray<int> dirRot; //$EDC$[varaccess,remove,EDCvarname="constrained_rotations",EDCfolder="Physics.Lagrange"]
	//EDC	TArray<int> dir; //$EDC$[varaccess,readonly,EDCvarname="constrained_directions",EDCfolder="Physics.Lagrange",tooltiptext="constrained directions cannot be changed"]
	//EDC	TArray<int> dirRot; //$EDC$[varaccess,readonly,EDCvarname="constrained_rotations",EDCfolder="Physics.Lagrange",tooltiptext="constrained rotations cannot be changed"]

	Vector3D rotAxis;		 //$EDC$[varaccess,EDCvarname="rotation_axis",EDCfolder="Physics",tooltiptext="local rotation axis w.r.t body 1 coordinate system"]

	//EDC	int standard_joint_drawing; //$EDC$[varaccess,EDCvarname="standard_joint_drawing",EDCfolder="Graphics",tooltiptext="flag for drawing mode; 1 == draw constraint element; 0 == show constrained directions and rotations;",int_bool]
	double draw_diameter; //$EDC$[varaccess,EDCvarname="diameter",EDCfolder="Graphics",tooltiptext="diameter of the revolute joint (for drawing)"]
	double draw_axis_length; //$EDC$[varaccess,EDCvarname="axis_length",EDCfolder="Graphics",tooltiptext="axis length of the revolute joint (for drawing)"]

	//EDC	double damping_coeff;	//$EDC$[varaccess,minval=0,EDCvarname="damping",EDCfolder="Physics.Penalty",tooltiptext="damping parameter used for translation and rotation"]
	//EDC	double spring_stiffness; //$EDC$[varaccess,minval=0,EDCvarname="stiffness",EDCfolder="Physics.Penalty",tooltiptext="stiffness parameter used for translation and rotation"]

	//remove unused data
	//EDC	Matrix3D JA0i;	//$EDC$[varaccess,remove,EDCvarname="joint_local_frame",EDCfolder="Geometry"]
	//EDC	Vector3D bryantAngles;	//$EDC$[varaccess,remove,EDCvarname="joint_local_frame_in_bryant_angles",EDCfolder="Geometry"]
	//EDC	Matrix3D penaltyStiffness; //$EDC$[varaccess,remove,EDCvarname="stiffness_matrix",EDCfolder="Physics.Penalty"]
	//EDC	Matrix3D penaltyDamping; //$EDC$[varaccess,remove,EDCvarname="damping_matrix",EDCfolder="Physics.Penalty"]
	//EDC	Matrix3D penaltyStiffnessRot; //$EDC$[varaccess,remove,EDCvarname="stiffness_matrix_rotation",EDCfolder="Physics.Penalty"]
	//EDC	Matrix3D penaltyDampingRot; //$EDC$[varaccess,remove,EDCvarname="damping_matrix_rotation",EDCfolder="Physics.Penalty"]

};//$EDC$[endclass,RevoluteJoint]


// SM is responsible for this constraint
//$ SW 2013-08-29: changed class name from PrismaticJoint_ to PrismaticJoint
class PrismaticJoint: public BaseBodyJoint  //$EDC$[beginclass,classname=PrismaticJoint,parentclassname=BaseBodyJoint,addelementtype=TAEconstraint,addelementtypename=PrismaticJoint,figure="prismaticJoint",texdescription="
//The PrismaticJoint constrains all relative degrees of freedom between two bodies except the translation along a local sliding axis. A penalty formulation exists, which replaces the exact lagrange constraint by a approximation with joint stiffness and damping.",
//example="PrismaticJointShort.txt"]
{
public:

	PrismaticJoint(MBS* mbsi):BaseBodyJoint(mbsi)
	{	
		ElementDefaultConstructorInitialization();
	};

	PrismaticJoint(const PrismaticJoint& ct): BaseBodyJoint(ct.mbs)
	{
		ElementDefaultConstructorInitialization();
		CopyFrom(ct);
	};

	~PrismaticJoint()
	{
	};

	void SetPrismaticJointLocalPos_to_LocalPos(int elem1, Vector3D lc1, int elem2, Vector3D lc2, Vector3D sliding_dirI, double stiffnessI = 0., double stiffnessRotI = 0., double dampingI = 0., double dampingRotI = 0.)
	{
		slidingDirection = sliding_dirI;
		slidingDirection.Normalize();
		SetPos1ToLocalCoord(elem1, lc1);
		SetPos2ToLocalCoord(elem2, lc2);
	}

	virtual void ElementDefaultConstructorInitialization();

	virtual void Initialize();

	virtual int CheckConsistency(mystr& errorstr);

	virtual Element* GetCopy()
	{
		Element* ec = new PrismaticJoint(*this);
		return ec;
	}

	virtual void CopyFrom(const Element& e);

	virtual const char* GetElementSpec() const {return "PrismaticJoint";}

  virtual int GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //for special value sensor

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		
	virtual int SetElementDataAuto(ElementDataContainer& edc);
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual void DrawElement();

	virtual double GetDrawSizeAxisLength();
	virtual Vector3D GetDrawSize();

protected:
	//EDC	TArray<int> dir; //$EDC$[varaccess,remove,EDCvarname="constrained_directions",EDCfolder="Physics.Lagrange"]
	//EDC	TArray<int> dirRot; //$EDC$[varaccess,remove,EDCvarname="constrained_rotations",EDCfolder="Physics.Lagrange"]
	//EDC	TArray<int> dir; //$EDC$[varaccess,readonly,EDCvarname="constrained_directions",EDCfolder="Physics.Lagrange",tooltiptext="constrained directions cannot be changed"]
	//EDC	TArray<int> dirRot; //$EDC$[varaccess,readonly,EDCvarname="constrained_rotations",EDCfolder="Physics.Lagrange",tooltiptext="constrained rotations cannot be changed"]

	Vector3D slidingDirection;		 //$EDC$[varaccess,EDCvarname="sliding_direction",EDCfolder="Physics",tooltiptext="local sliding direction w.r.t body 1 coordinate system"]

	//EDC	int standard_joint_drawing; //$EDC$[varaccess,EDCvarname="standard_joint_drawing",EDCfolder="Graphics",tooltiptext="flag for drawing mode; 1 == draw constraint nicely; 0 == show constrained directions and rotations;",int_bool]
	double draw_length; //$EDC$[varaccess,EDCvarname="rail_length",EDCfolder="Graphics",tooltiptext="length of the prismatic joint (for drawing)"]
	Vector3D draw_cube_size; //$EDC$[varaccess,EDCvarname="joint_cube_size",EDCfolder="Graphics",tooltiptext="cube dimension of prismatic joint (for drawing); [lx (in sl. dir.),ly (normal to sl. dir.),lz (normal to sl. dir.)]"]

	//EDC	double damping_coeff;	//$EDC$[varaccess,minval=0,EDCvarname="damping",EDCfolder="Physics.Penalty",tooltiptext="damping parameter used for translation and rotation"]
	//EDC	double spring_stiffness; //$EDC$[varaccess,minval=0,EDCvarname="stiffness",EDCfolder="Physics.Penalty",tooltiptext="stiffness parameter used for translation and rotation"]

	//remove unused data
	//EDC	Matrix3D JA0i;	//$EDC$[varaccess,remove,EDCvarname="joint_local_frame",EDCfolder="Geometry"]
	//EDC	Vector3D bryantAngles;	//$EDC$[varaccess,remove,EDCvarname="joint_local_frame_in_bryant_angles",EDCfolder="Geometry"]
	//EDC	Matrix3D penaltyStiffness; //$EDC$[varaccess,remove,EDCvarname="stiffness_matrix",EDCfolder="Physics.Penalty"]
	//EDC	Matrix3D penaltyDamping; //$EDC$[varaccess,remove,EDCvarname="damping_matrix",EDCfolder="Physics.Penalty"]
	//EDC	Matrix3D penaltyStiffnessRot; //$EDC$[varaccess,remove,EDCvarname="stiffness_matrix_rotation",EDCfolder="Physics.Penalty"]
	//EDC	Matrix3D penaltyDampingRot; //$EDC$[varaccess,remove,EDCvarname="damping_matrix_rotation",EDCfolder="Physics.Penalty"]

};//$EDC$[endclass,PrismaticJoint]


//$ SW 2013-9-18: added UniversalJoint
class UniversalJoint: public BaseBodyJoint  //$EDC$[beginclass,classname=UniversalJoint,
//parentclassname=BaseBodyJoint,addelementtype=TAEconstraint,addelementtypename=UniversalJoint,
//texdescription="The UniversalJoint constains the local position of two elements and keeps two axes, one on each body, perpendicular to each other.",
//figure="UniversalJoint_with_beams",
//texdescriptionDOF="The vector of the DOF contains the Lagrangian parameters $\lambda = 
//	\left[\lambda_1\ \lambda_2\ \lambda_3\ \lambda_4 \right]^{\mathbf{T}}$, where $\lambda_1, \lambda_2, \lambda_3$ are measures for the violation of the displacement condition and $\lambda_4$ is a measure for the violation of the orthogonality condition of the two axes.",
//texdescriptionGeometry="For this constraint one needs to specify the axes of the cross and the directions in which the hinges are drawn. The direction of the hinge and the axis connected to this hinge have to be given in the local coordinate system of the respective body. See figure \ref{UniversalJointfigure2}",
//figure="UniversalJoint_detail",
//texdescriptionEquations="The positions and axes are given in local coordinates of body 1 respectively body 2. However the calculations are done internally in global coordinates.\\Let
	//\[ \mathbf{x}^i = \left[
	//	\begin{array}{ccc}
	//		x_1^i & x_2^i & x_3^i
	//	\end{array}
	//\right]^{\mathbf{T}}\]
	//be the position (in global coordinates) where the joint is connected to the first body and let
	//\[ \mathbf{x}^j = \left[
	//	\begin{array}{ccc}
	//		x_1^j & x_2^j & x_3^j
	//	\end{array}
	//\right]^{\mathbf{T}}\]
	//be the position (in global coordinates) where the joint is connected to the second body.\\
	//Let
	//\[ \mathbf{a}^i = 
	//	\left[\begin{array}{ccc}
	//		a_1^i & a_2^i & a_3^i
	//	\end{array}\right]^{\mathbf{T}}\]
	//be the axis (in global coordinates) connected to the first body and let
	//\[ \mathbf{a}^j = 
	//	\left[\begin{array}{ccc}
	//		a_1^j & a_2^j & a_3^j
	//	\end{array}\right]^{\mathbf{T}}\]
	//be the axis (in global coordinates) connected to the second body.
	//Then the constraint equations at position level are
	//\[\mathbf{C} =
	//\left[\begin{array}{c}
	//	\mathbf{x}^i-\mathbf{x}^j \\
	//	{\mathbf{a}^i}^{\mathbf{T}}\cdot\mathbf{a}^j
	//\end{array}\right]=\mathbf 0.\]\\
	//The first three constraints restrict the position of the connection points of body 1 and 2. The fourth equation ensures that the two axes of the cross are perpendicular to each other.\\
	//The constraint equations at velocity level are
	//\[ \mathbf C =
	//\left[ \begin{array}{c}
	//	\frac{\partial \mathbf x^i}{\partial t}-\frac{\partial \mathbf x^j}{\partial t} \\
	//	{\frac{\mathbf{a}^i}{\partial t}}^{\mathbf{T}}\cdot\frac{\mathbf{a}^j}{\partial t}
	//\end{array} \right] = 0. \]",
//texdescriptionLimitations="No penalty formulation is available.",
//example="UniversalJoint.txt"
//]
{
public:


	UniversalJoint(MBS* mbsi):BaseBodyJoint(mbsi)
	{	
		ElementDefaultConstructorInitialization();
	};

	UniversalJoint(const UniversalJoint& ct): BaseBodyJoint(ct.mbs)
	{
		ElementDefaultConstructorInitialization();
		CopyFrom(ct);
	};

	~UniversalJoint()
	{
	};

	virtual int IS() const;

	virtual void ElementDefaultConstructorInitialization();

	virtual void SetUniversalJoint(int elem1, Vector3D lc1, Vector3D axis_1, int elem2, Vector3D lc2, Vector3D axis_2); 

	virtual void Initialize();

	virtual int CheckConsistency(mystr& errorstr);

	virtual void GetNecessaryKinematicAccessFunctions(TArray<int> &KAF, int numberOfKinematicPair);
	
	virtual Element* GetCopy()
	{
		Element* ec = new UniversalJoint(*this);
		return ec;
	}

	virtual void CopyFrom(const Element& e);

	virtual const char* GetElementSpec() const {return "UniversalJoint";}

  virtual int GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //for special value sensor

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		
	virtual int SetElementDataAuto(ElementDataContainer& edc);
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType&  RWdata);
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata);
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter


	virtual void DrawElement();

	virtual double GetDrawSizeAxisLength();
	virtual Vector3D GetDrawSize();

	
	//for Lagrange formulation only
	virtual void EvalG(Vector& f, double t);

	//for Penalty formulation only
	virtual void EvalF2(Vector& f, double t);

	//for Lagrange formulation only
	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);
	

protected:
	//EDC	TArray<int> dir; //$EDC$[varaccess,remove,EDCvarname="constrained_directions",EDCfolder="Physics.Lagrange"]
	//EDC	TArray<int> dirRot; //$EDC$[varaccess,remove,EDCvarname="constrained_rotations",EDCfolder="Physics.Lagrange"]

	//EDC	Vector3D col; //$EDC$[varaccess, remove, EDCvarname="color_body1",EDCfolder="Graphics"] //color for constraint
	//EDC	Vector3D col_ext; //$EDC$[varaccess, remove, EDCvarname="color_body2",EDCfolder="Graphics"] //color for constraint
	//EDC	Vector3D col; //$EDC$[varaccess, EDCvarname="color_body1", EDCfolder="Graphics",tooltiptext="[red, green, blue] color of the hinge connected to the first body, range = 0..1"] //color for constraint
	//EDC	Vector3D col_ext; //$EDC$[varaccess, EDCvarname="color_body2",EDCfolder="Graphics",tooltiptext="[red, green, blue] color of the hinge connected to the first body, range = 0..1"] //color for constraint
	Vector3D col_t; //$EDC$[varaccess,EDCvarname="color_cross",EDCfolder="Graphics",tooltiptext="[red, green, blue] color of the cross shaft"]

	//EDC Vector3D loccoords(1); //$EDC$[varaccess, remove, EDCvarname="position",EDCfolder="Position1"]
	//EDC Vector3D loccoords(2); //$EDC$[varaccess, remove, EDCvarname="position",EDCfolder="Position2"]
	//EDC Vector3D loccoords(1); //$EDC$[varaccess, EDCvarname="position",EDCfolder="Position1", tooltiptext="local position"]
	//EDC Vector3D loccoords(2); //$EDC$[varaccess, EDCvarname="position",EDCfolder="Position2", tooltiptext="local or global (if element_number == 0) position"]

	Vector3D axis_1;		 //$EDC$[varaccess,EDCvarname="axis",EDCfolder="Position1",tooltiptext="the axis of the cross connected to body 1 in local coordinates"]
	Vector3D axis_2;		 //$EDC$[varaccess,EDCvarname="axis",EDCfolder="Position2",tooltiptext="the axis of the cross connected to body 2 in local coordinates"]

	double draw_length; //$EDC$[varaccess,EDCvarname="draw_length",EDCfolder="Graphics",tooltiptext="length of the universal joint (for drawing)"]
	double draw_width;  //$EDC$[varaccess,EDCvarname="draw_width",EDCfolder="Graphics",tooltiptext="width of the universal joint (for drawing)"]
	Vector3D draw_direction_1; //$EDC$[varaccess,EDCvarname="draw_direction_1",EDCfolder="Graphics",tooltiptext="direction from body 1 to joint (for drawing)"]
	Vector3D draw_direction_2; //$EDC$[varaccess,EDCvarname="draw_direction_2",EDCfolder="Graphics",tooltiptext="direction from body 2 to joint (for drawing)"]
	
	//remove unused data

	//EDC int use_penalty_formulation;		//$EDC$[varaccess,remove,EDCvarname="use_penalty_formulation",EDCfolder="Physics"]
	static const int use_penalty_formulation = false;
	//EDC double	draw_local_frame_size; //$EDC$[varaccess, remove, EDCvarname="draw_size_joint_local_frame",EDCfolder="Graphics"]

	//EDC	double damping_coeff;	//$EDC$[varaccess,remove,EDCvarname="damping",EDCfolder="Physics.Penalty"]
	//EDC	double spring_stiffness; //$EDC$[varaccess,remove,EDCvarname="stiffness",EDCfolder="Physics.Penalty"]

	//EDC double draw_cone_size; //$EDC$[varaccess,remove, EDCvarname="cone_size",EDCfolder="Graphics"]
	//EDC	Matrix3D JA0i;	//$EDC$[varaccess,remove,EDCvarname="joint_local_frame",EDCfolder="Geometry"]
	//EDC	Vector3D bryantAngles;	//$EDC$[varaccess,remove,EDCvarname="joint_local_frame_in_bryant_angles",EDCfolder="Geometry"]
	//EDC	Matrix3D penaltyStiffness; //$EDC$[varaccess,remove,EDCvarname="stiffness_matrix",EDCfolder="Physics.Penalty"]
	//EDC	Matrix3D penaltyDamping; //$EDC$[varaccess,remove,EDCvarname="damping_matrix",EDCfolder="Physics.Penalty"]
	//EDC	Matrix3D penaltyStiffnessRot; //$EDC$[varaccess,remove,EDCvarname="stiffness_matrix_rotation",EDCfolder="Physics.Penalty"]
	//EDC	Matrix3D penaltyDampingRot; //$EDC$[varaccess,remove,EDCvarname="damping_matrix_rotation",EDCfolder="Physics.Penalty"]

	//member variables used as temporary matrices
	Vector tmp_Vec1, tmp_Vec2;
	Matrix tmp_dRotvdqT;
	Vector tmp_v1, tmp_vp1, tmp_v2, tmp_vp2;
};//$EDC$[endclass,UniversalJoint]

// SM is responsible for this constraint
//$ SW 2013-08-29: changed class name from RigidJoint_ to RigidJoint
class RigidJoint: public BaseBodyJoint  //$EDC$[beginclass,classname=RigidJoint,parentclassname=BaseBodyJoint,addelementtype=TAEconstraint,addelementtypename=RigidJoint,figure="rigidJoint",texdescription="
//The RigidJoint constrains the position and relative angles of an element at a specified local position. If only one element is specified, a ground joint is realized. A penalty formulation exists, which replaces the exact lagrange constraint by a approximation with joint stiffness and damping.",
//example="RigidJointShort.txt"]
{
public:

	RigidJoint(MBS* mbsi):BaseBodyJoint(mbsi)
	{	
		ElementDefaultConstructorInitialization();
	};

	RigidJoint(const RigidJoint& ct): BaseBodyJoint(ct.mbs)
	{
		ElementDefaultConstructorInitialization();
		CopyFrom(ct);
	};

	~RigidJoint()
	{
	};

	void SetRigidJoint_LocalPos_to_LocalPos(int elem1, Vector3D lc1, int elem2, Vector3D lc2)
	{
		SetPos1ToLocalCoord(elem1, lc1);
		SetPos2ToLocalCoord(elem2, lc2);
	}

	virtual void ElementDefaultConstructorInitialization();

	virtual void Initialize();

	virtual int CheckConsistency(mystr& errorstr);

	virtual Element* GetCopy()
	{
		Element* ec = new RigidJoint(*this);
		return ec;
	}

	virtual void CopyFrom(const Element& e);

	virtual const char* GetElementSpec() const {return "RigidJoint";}

  virtual int GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //for special value sensor

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		
	virtual int SetElementDataAuto(ElementDataContainer& edc);
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual void DrawElement();

	virtual double GetDrawSizeScalar();

protected:
	//EDC	TArray<int> dir; //$EDC$[varaccess,remove,EDCvarname="constrained_directions",EDCfolder="Physics.Lagrange"]
	//EDC	TArray<int> dirRot; //$EDC$[varaccess,remove,EDCvarname="constrained_rotations",EDCfolder="Physics.Lagrange"]
	//EDC	TArray<int> dir; //$EDC$[varaccess,readonly,EDCvarname="constrained_directions",EDCfolder="Physics.Lagrange",tooltiptext="constrained directions cannot be changed"]
	//EDC	TArray<int> dirRot; //$EDC$[varaccess,readonly,EDCvarname="constrained_rotations",EDCfolder="Physics.Lagrange",tooltiptext="constrained rotations cannot be changed"]

	//EDC	int standard_joint_drawing; //$EDC$[varaccess,EDCvarname="standard_joint_drawing",EDCfolder="Graphics",tooltiptext="flag for drawing mode; 1 == draw constraint element; 0 == show constrained directions and rotations;",int_bool]
	double draw_dimension; //$EDC$[varaccess,EDCvarname="cube_length",EDCfolder="Graphics",tooltiptext="rigid joint dimension (for drawing)"]

	//EDC	double damping_coeff;	//$EDC$[varaccess,minval=0,EDCvarname="damping",EDCfolder="Physics.Penalty",tooltiptext="damping parameter used for translation and rotation"]
	//EDC	double spring_stiffness; //$EDC$[varaccess,minval=0,EDCvarname="stiffness",EDCfolder="Physics.Penalty","tooltiptext="stiffness parameter used for translation and rotation"]

	//remove unused data
	//EDC	Matrix3D JA0i;	//$EDC$[varaccess,remove,EDCvarname="joint_local_frame",EDCfolder="Geometry"]
	//EDC	Vector3D bryantAngles;	//$EDC$[varaccess,remove,EDCvarname="joint_local_frame_in_bryant_angles",EDCfolder="Geometry"]
	//EDC	Matrix3D penaltyStiffnessRot; //$EDC$[varaccess,remove,EDCvarname="stiffness_matrix_rotation",EDCfolder="Physics.Penalty"]
	//EDC	Matrix3D penaltyDampingRot; //$EDC$[varaccess,remove,EDCvarname="damping_matrix_rotation",EDCfolder="Physics.Penalty"]
	//EDC	Matrix3D penaltyStiffness; //$EDC$[varaccess,remove,EDCvarname="stiffness_matrix",EDCfolder="Physics.Penalty"]
	//EDC	Matrix3D penaltyDamping; //$EDC$[varaccess,remove,EDCvarname="damping_matrix",EDCfolder="Physics.Penalty"]

};//$EDC$[endclass,RigidJoint]


// SM is responsible for this constraint
//$ SW 2013-08-29: changed class name from CylindricalJoint_ to CylindricalJoint
class CylindricalJoint: public BaseBodyJoint  //$EDC$[beginclass,classname=CylindricalJoint,parentclassname=BaseBodyJoint,addelementtype=TAEconstraint,addelementtypename=CylindricalJoint,figure="cylindricalJoint",texdescription="
//The CylindricalJoint constrains like the RevoluteJoint, but allows additionally translation along the rotational axis. A penalty formulation exists, which replaces the exact lagrange constraint by a approximation with joint stiffness and damping.",
//example="CylindricalJointShort.txt"]
{
public:

	CylindricalJoint(MBS* mbsi):BaseBodyJoint(mbsi)
	{	
		ElementDefaultConstructorInitialization();
	};

	CylindricalJoint(const CylindricalJoint& ct): BaseBodyJoint(ct.mbs)
	{
		ElementDefaultConstructorInitialization();
		CopyFrom(ct);
	};

	~CylindricalJoint()
	{
	};

	virtual void ElementDefaultConstructorInitialization();

	virtual void Initialize();

	virtual int CheckConsistency(mystr& errorstr);

	virtual Element* GetCopy()
	{
		Element* ec = new CylindricalJoint(*this);
		return ec;
	}

	virtual void CopyFrom(const Element& e);

	virtual const char* GetElementSpec() const {return "CylindricalJoint";}

  virtual int GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //for special value sensor

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		
	virtual int SetElementDataAuto(ElementDataContainer& edc);
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual void DrawElement();

	virtual double GetDrawSizeAxisLength();
	virtual Vector2D GetDrawSize();

protected:
	//EDC	TArray<int> dir; //$EDC$[varaccess,remove,EDCvarname="constrained_directions",EDCfolder="Physics.Lagrange"]
	//EDC	TArray<int> dirRot; //$EDC$[varaccess,remove,EDCvarname="constrained_rotations",EDCfolder="Physics.Lagrange"]
	//EDC	TArray<int> dir; //$EDC$[varaccess,readonly,EDCvarname="constrained_directions",EDCfolder="Physics.Lagrange",tooltiptext="constrained directions cannot be changed"]
	//EDC	TArray<int> dirRot; //$EDC$[varaccess,readonly,EDCvarname="constrained_rotations",EDCfolder="Physics.Lagrange",tooltiptext="constrained rotations cannot be changed"]

	Vector3D rotSlideAxis;		 //$EDC$[varaccess,EDCvarname="rotation_sliding_axis",EDCfolder="Physics",tooltiptext="local rotation/sliding axis w.r.t body 1 coordinate system"]

	//EDC	int standard_joint_drawing; //$EDC$[varaccess,EDCvarname="standard_joint_drawing",EDCfolder="Graphics",tooltiptext="flag for drawing mode; 1 == draw constraint element; 0 == show constrained directions and rotations;",int_bool]
	Vector2D draw_cylinder_size;  //$EDC$[varaccess,EDCvarname="joint_cylinder_size",EDCfolder="Graphics",tooltiptext="cylinder dimension of cylindrical joint (for drawing); [lx (cyl. length, in sl. dir.),d (cylinder diameter)]"]
	double draw_axis_length; //$EDC$[varaccess,EDCvarname="axis_length",EDCfolder="Graphics",tooltiptext="axis length of the revolute joint (for drawing)"]

	//EDC	double damping_coeff;	//$EDC$[varaccess,minval=0,EDCvarname="damping",EDCfolder="Physics.Penalty",tooltiptext="damping parameter used for translation and rotation"]
	//EDC	double spring_stiffness; //$EDC$[varaccess,minval=0,EDCvarname="stiffness",EDCfolder="Physics.Penalty",tooltiptext="stiffness parameter used for translation and rotation"]

	//remove unused data
	//EDC	Matrix3D JA0i;	//$EDC$[varaccess,remove,EDCvarname="joint_local_frame",EDCfolder="Geometry"]
	//EDC	Vector3D bryantAngles;	//$EDC$[varaccess,remove,EDCvarname="joint_local_frame_in_bryant_angles",EDCfolder="Geometry"]
	//EDC	Matrix3D penaltyStiffnessRot; //$EDC$[varaccess,remove,EDCvarname="stiffness_matrix_rotation",EDCfolder="Physics.Penalty"]
	//EDC	Matrix3D penaltyDampingRot; //$EDC$[varaccess,remove,EDCvarname="damping_matrix_rotation",EDCfolder="Physics.Penalty"]
	//EDC	Matrix3D penaltyStiffness; //$EDC$[varaccess,remove,EDCvarname="stiffness_matrix",EDCfolder="Physics.Penalty"]
	//EDC	Matrix3D penaltyDamping; //$EDC$[varaccess,remove,EDCvarname="damping_matrix",EDCfolder="Physics.Penalty"]

};//$EDC$[endclass,CylindricalJoint]

//***************************************************************************************************************************************************
// SpringDamperActuator
//***************************************************************************************************************************************************

//Spring-Damper-Actuator connects two elements, with spring (k...stiffness, l0...initial length),
// damper (d ...  damping coefficient) and a constant force (fa ... constant actuator force)
class SpringDamperActuator: public BasePointJoint  //$EDC$[beginclass,classname=SpringDamperActuator,parentclassname=BasePointJoint,addelementtype=TAEconstraint,addelementtypename=SpringDamperActuator,figure="springDamperActuator",texdescription="
//The Spring-Damper-Actuator connects two points with a spring, a damper and a actor element, in which actuator force fa remains constant. The resultant force is applied in the connection line of these points. 
//There are different modes available, how the spring and damper force is calculated. It is also possible to change the neutral spring length. This joint is realized in Penalty formulation only.",
//modus="{element to ground}{Position2.element\_number AND Position2.node\_number have to be equal to 0}",
//modus="{element to element}{Position2.element\_number and/or Position2.node\_number must not be equal to 0}",
//modus="{forcemode}{
//\textbf{Physics.forcemode = 0:} \newline Force is computed as (a) with constant stiffness and damping factors $k$ and $d$. The factors can be defined in the two fields in Physics.Linear. \newline 
//\textbf{Physics.forcemode = 1:} \newline 2 MathFunctions are used to describe piecewise linear stiffness $k\left(\Delta x\right)$ and damping $d\left(v\right)$, see formula (b) and Physics.MathFunction. \newline 
//\textbf{Physics.forcemode = 2:} \newline 2 IOElementDataModifiers describe the force (c) due to stiffness and damping. You should use this mode if full nonlinear behavior is required, e.g. $f_k = f_k\left(t,l,v,...\right)$ and $f_d=d\left(t,l,v,... \right)$.  \newline
//\textbf{Physics.forcemode = 3:} \newline 2 MathFunctions are used to describe piecewise linear spring force $f_k\left(\Delta x\right)$ and damping force $f_d\left(v\right)$, see formula (d) and Physics.MathFunction. \newline \newline 
//modifier value names for forcemode == 2: \newline $f_k$: Connector.spring\_force' \newline $f_d$: Connector.damper\_force'}",
//modus="{spring length offset}{It is possible to change the spring length $l_0$ (neutral length of the spring) during the simulation, e.g. for the usage of the SpringDamperActuator as a linear actuator. In 
//standard mode the value in the field Physics.spring\_length remains constant. This value can be modified by a IOElementDataModifier via 'Connector.spring\_length\_offset'.}",
//modus="{additional actor force}{In Physics.actor\_force a constant offset force $f_a$ can be added.}",
//texdescriptionLimitations="If the 2 end points of the spring are the same point in the initial configuration, this may lead to problems! The direction of the spring can not be determined in that case!",
//texdescriptionEquations="
//point positions: $\mathbf{p}^{(1)} = \left[p_x^{(1)} p_y^{(1)} p_z^{(1)}\right]^T$; \qquad $\mathbf{p}^{(2)} = \left[p_x^{(2)} p_y^{(2)} p_z^{(2)}\right]^T$. \\ \\
//point velocities: $\mathbf{\dot{p}}^{(1)} = \left[\dot{p}_x^{(1)} \dot{p}_y^{(1)} \dot{p}_z^{(1)}\right]^T$; \qquad $\mathbf{\dot{p}}^{(2)} = \left[\dot{p}_x^{(2)} \dot{p}_y^{(2)} \dot{p}_z^{(2)}\right]^T$. \\ \\
//spring length: $l_0$ \\ \\
//direction vector: $\mathbf{dir} = \frac{\mathbf{p}^{(1)}-\mathbf{p}^{(2)}}{\sqrt{\left(p_x^{(1)}-p_x^{(2)}\right)^2+\left(p_y^{(1)}-p_y^{(2)}\right)^2+\left(p_z^{(1)}-p_z^{(2)}\right)^2}}$ \\
//spring deflection: $\Delta x = l - l_0 = \left(\mathbf{p}^{(1)}-\mathbf{p}^{(2)}\right)^T\,\mathbf{dir} - l_0$ \\
//spring velocity: $v = \left(\mathbf{\dot{p}}^{(1)}-\mathbf{\dot{p}}^{(2)}\right)^T\,\mathbf{dir}$ \\ \\
//resultant force (see section forcemode): \\
//forcemode 0: $f = k\,\Delta x+d\,v+f_a$ (a) \\ forcemode 1: $f = k\left(\Delta x\right)\Delta x+d\left(v\right)v+f_a$ (b) \\ forcemode 2: $f = f_k+f_d+f_a$ (c) \\ forcemode 3: $f=f_k\left(\Delta x\right)+f_k\left(v\right)+f_a$ (d)",example="SpringDamperActuator.txt"]
{
public:

	SpringDamperActuator(MBS* mbsi):BasePointJoint(mbsi)
	{	
		ElementDefaultConstructorInitialization();
	};

	SpringDamperActuator(const SpringDamperActuator& ct): BasePointJoint(ct.mbs)
	{
		ElementDefaultConstructorInitialization();
		CopyFrom(ct);
	};

	~SpringDamperActuator()
	{
	};

	virtual void Initialize();

	virtual Element* GetCopy()
	{
		Element* ec = new SpringDamperActuator(*this);
		return ec;
	}

	virtual void CopyFrom(const Element& e);

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		
	virtual int SetElementDataAuto(ElementDataContainer& edc);
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

  virtual int GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //for special value sensor
	virtual int ReadSingleElementData(ReadWriteElementDataVariableType& RWdata); 		//for special value sensor
	virtual int WriteSingleElementData(const ReadWriteElementDataVariableType& RWdata);

	virtual void ElementDefaultConstructorInitialization();
	virtual void SetStiffnessDamping();
	virtual void SetPiecewiseLinearSpring(const Vector& xpos, const Vector& force);
	virtual void SetPiecewiseLinearDamper(const Vector& vel, const Vector& force);
	//sets forcemode to 3 and sets the piecewise points and the corresponsing force values.
	virtual void SetPiecewiseLinearSpringForce(const Vector& xpos, const Vector& force);
	virtual void SetPiecewiseLinearDamperForce(const Vector& vel, const Vector& force);

	virtual const char* GetElementSpec() const {return "SpringDamperActuator";}

	// for constraints only, these are the necessary (!) access functions!	//$ MSax 2013-02-18
	virtual void GetNecessaryKinematicAccessFunctions(TArray<int> &KAF, int numberOfKinematicPair);

	virtual void DrawElement();

	virtual Vector3D ComputeForce(double t) const;

	virtual Vector3D ComputeForceDirection() const; 
	virtual Vector3D ComputeForceDirectionD(); 

	// Penalty + Lagrange: Element to Element
	void SetSpringDamperActuator_LocalPos_to_LocalPos(int elem1, Vector3D lc1, int elem2, Vector3D lc2, double stiffness, double damping = 0.0, double spring_length=0.0)
	{
		SetPos1ToLocalCoord(elem1, lc1);
		SetPos2ToLocalCoord(elem2, lc2);
		SetDampingCoeff(damping);
		SetPenaltyStiffness(stiffness);
		SetPenaltyFormulation(1);
		l0 = spring_length;
	}

protected:
	double l0; //$EDC$[varaccess,EDCvarname="spring_length",EDCfolder="Physics",tooltiptext="length of the spring in the initial configuration"]
	double fa; //$EDC$[varaccess,EDCvarname="actor_force",EDCfolder="Physics",tooltiptext="constant force acting on the spring"]
	int forcemode; //$EDC$[varaccess,EDCvarname="forcemode",EDCfolder="Physics",tooltiptext="defines how the spring and damper force is computed: 0..constant coefficient, 1..MathFunction (stiffness and damping), 2..IOElementDataModifier, 3..MathFunction (spring force and damping force)"]
	
	MathFunction mathfunc_k; //evaluate(x) gives the stiffness as function of x
	MathFunction mathfunc_d; //evaluate(x) gives damping as function of xp
	MathFunction mathfunc_fk; //evaluate(x) gives the spring force as function of x
	MathFunction mathfunc_fd; //evaluate(x) gives the damping force as function of x

	double force_k;
	double force_d;

	//EDC	Vector3D col; //$EDC$[varaccess,remove,EDCvarname="RGB_color",EDCfolder="Graphics"] //color for body
	//EDC	Vector3D col; //$EDC$[varaccess,EDCvarname="color_body1",EDCfolder="Graphics",tooltiptext="[red, green, blue] first color of constraint (spring), range = 0..1, use default color:[-1,-1,-1]"] //color for constraint
	//EDC	Vector3D col_ext; //$EDC$[varaccess,EDCvarname="color_body2",EDCfolder="Graphics",tooltiptext="[red, green, blue] second color of constraint (damper), range = 0..1, use default color:[-1,-1,-1]"] //color for constraint

	//EDC double spring_stiffness; //$EDC$[varaccess,remove,EDCvarname="spring_stiffness",EDCfolder="Physics.Penalty"]
	//EDC double spring_stiffness; //$EDC$[varaccess,EDCvarname="spring_stiffness",minval=0,EDCfolder="Physics.Linear",tooltiptext="stiffness coefficient of the linear spring. Only used if forcemode is 0."]
	//EDC double damping_coeff; //$EDC$[varaccess,remove,EDCvarname="damping",EDCfolder="Physics.Penalty"]
	//EDC double damping_coeff; //$EDC$[varaccess,EDCvarname="damping",minval=0,EDCfolder="Physics.Linear",tooltiptext="damping coefficient for viscous damping. Only used if forcemode is 0."]

	// remove unused components from class
  //EDC int use_penalty_formulation;				//$EDC$[varaccess,remove,EDCvarname="use_penalty_formulation",EDCfolder="Physics"]
	//EDC int stiffness_in_joint_local_frame;	//$EDC$[varaccess,remove,EDCvarname="use_joint_local_frame",EDCfolder="Geometry"]
	//EDC TArray<int> dir;												//$EDC$[varaccess,remove,EDCvarname="constrained_directions",minval=0,maxval=1,EDCfolder="Physics.Lagrange"]
	//EDC Vector3D spring_stiffness3;					//$EDC$[varaccess,remove,EDCvarname="spring_stiffness_vector",EDCfolder="Physics.Penalty"]
	//EDC Matrix3D JA0i;											//$EDC$[varaccess,remove,EDCvarname="joint_local_frame",EDCfolder="Geometry"]
	//EDC int use_local_coordinate_system;	//$EDC$[varaccess,remove,EDCvarname="use_local_coordinate_system",EDCfolder="Geometry"]
	//EDC double draw_dim(1);						//$EDC$[varaccess,remove,EDCvarname="draw_size",EDCfolder="Graphics"]
	//EDC double draw_local_frame_size;			//$EDC$[varaccess,remove,EDCvarname="draw_size_joint_local_frame",EDCfolder="Graphics"]

	//EDC double draw_dim(1);						//$EDC$[varaccess,EDCvarname="spring_diameter",EDCfolder="Graphics",minval=0,tooltiptext="spring diameter used for drawing only."]
	//EDC double draw_dim(2);						//$EDC$[varaccess,EDCvarname="spring_coils",EDCfolder="Graphics",minval=0,tooltiptext="spring coils used for drawing. If set to 0, then a cylinder with the value 'spring_diameter' as diameter is shown instead of the coils."]
	double spring_res;									//$EDC$[varaccess,EDCvarname="spring_resolution",EDCfolder="Graphics",minval=1,maxval=10,tooltiptext="spring resolution used for drawing (very coarse = 1, very smooth = 10)."]
	//EDC double draw_dim(3);						//$EDC$[varaccess,EDCvarname="damper_diameter",EDCfolder="Graphics",minval=0,tooltiptext="damper diameter used for drawing only. If set to 0, then the damper is not shown. It's recommended to choose the value smaller then the spring diameter."]

};//$EDC$[endclass,SpringDamperActuator]

class RigidLink: public BasePointJoint  //$EDC$[beginclass,classname=RigidLink,parentclassname=BasePointJoint,addelementtype=TAEconstraint,addelementtypename=RigidLink,figure="rigidLink",texdescription="
//A rigid link is a rigid constraint element that provides a stiff connection between nodes or positions in the model. In standard mode the distance between the connected points remains constant. In extended mode it is possible to change the distance as a function of time or input. 
// There is only a Lagrange formulation implemented.",
//modus="{element to ground}{Position2.element\_number AND Position2.node\_number have to be equal to 0}",
//modus="{element to element}{Position2.element\_number and/or Position2.node\_number must not be equal to 0}",
//modus="{distancemode}{
//\textbf{Physics.distancemode = 0:} \newline The distance remains constant. The value can be defined in the field Physics.Constant.link\_length. \newline 
//\textbf{Physics.distancemode = 1:} \newline A MathFunction is used to describe piecewise linear distance or velocity development over time t, e.g. for a rigid link actuator. See Physics.MathFunction. \newline 
//\textbf{Physics.distancemode = 2:} \newline A IOElementDataModifier describes the developing distance or velocity over time t, e.g. for a rigid link actuator. See section limitations.}",
//texdescriptionLimitations="For a position constraint (index 3 solver) with variable distance it is necessary to define the link length $l_0$ as a function of time. In this case the velocity input $v$ (the derivative of the distance with respect to time) is not considered, see formula (a).
//Reverse conditions apply to the velocity constraint with formula (b).",
//texdescriptionEquations="
//point positions: $\mathbf{p}^{(1)} = \left[p_x^{(1)} p_y^{(1)} p_z^{(1)}\right]^T$; \qquad $\mathbf{p}^{(2)} = \left[p_x^{(2)} p_y^{(2)} p_z^{(2)}\right]^T$. \\ 
//point velocities: $\mathbf{\dot{p}}^{(1)} = \left[\dot{p}_x^{(1)} \dot{p}_y^{(1)} \dot{p}_z^{(1)}\right]^T$; \qquad $\mathbf{\dot{p}}^{(2)} = \left[\dot{p}_x^{(2)} \dot{p}_y^{(2)} \dot{p}_z^{(2)}\right]^T$. \\ 
//link length: $l_0$ \\ 
//time derivative of link length: $v$ (equates $\dot{l}_0$)\\ 
//direction vector: $\mathbf{dir} = \frac{\mathbf{p}^{(1)}-\mathbf{p}^{(2)}}{\sqrt{\left(p_x^{(1)}-p_x^{(2)}\right)^2+\left(p_y^{(1)}-p_y^{(2)}\right)^2+\left(p_z^{(1)}-p_z^{(2)}\right)^2}}$ \\
//position constraint: $\mathbf{C}= \left(\mathbf{p}^{(1)}-\mathbf{p}^{(2)}\right)^T\,\mathbf{dir}-l_0 = 0$ (a)\\
//velocity constraint: $\mathbf{\dot{C}}= \left(\mathbf{\dot{p}}^{(1)}-\mathbf{\dot{p}}^{(2)}\right)^T\,\mathbf{dir}-v = 0$ (b)\\
//$\mathbf{\frac{\partial{C}}{\partial{q}}}^T= \left(\mathbf{\frac{\partial{p^{(1)}}}{\partial{q}}}-\mathbf{\frac{\partial{p^{(2)}}}{\partial{q}}}\right)^T\mathbf{dir}$",
//example="RigidLink.txt"]
{
public:

	RigidLink(MBS* mbsi):BasePointJoint(mbsi)
	{	
		ElementDefaultConstructorInitialization();
	};

	RigidLink(const RigidLink& ct): BasePointJoint(ct.mbs)
	{
		ElementDefaultConstructorInitialization();
		CopyFrom(ct);
	};

	~RigidLink()
	{
	};

	virtual void Initialize();

	virtual Element* GetCopy()
	{
		Element* ec = new RigidLink(*this);
		return ec;
	}

	virtual void CopyFrom(const Element& e);

	virtual const char* GetElementSpec() const {return "RigidLink";}

	// for constraints only, these are the necessary (!) access functions!	//$ MSax 2013-02-18
	virtual void GetNecessaryKinematicAccessFunctions(TArray<int> &KAF, int numberOfKinematicPair);

	virtual void ElementDefaultConstructorInitialization();

	virtual int SOS() const;

	virtual int IS() const {return 1;}

	virtual void EvalG(Vector& f, double t);
	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	virtual void DrawElement();

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		
	virtual int SetElementDataAuto(ElementDataContainer& edc);
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

  virtual int GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //for special value sensor
	virtual int WriteSingleElementData(const ReadWriteElementDataVariableType& RWdata);

	virtual Vector3D ComputeForceDirection() const; 
	virtual Vector3D ComputeForceDirectionD(); 

protected:
	int distancemode; //$EDC$[varaccess,EDCvarname="distancemode",EDCfolder="Physics",tooltiptext="defines the distance: 0..constant distance, 1..MathFunction, 2..IOElementDataModifier"]
	
	double distance; //$EDC$[varaccess,EDCvarname="link_length",EDCfolder="Physics.Constant",tooltiptext="constant distance is used, when distancemode = 0"]

	MathFunction mathfunc_l; //evaluate(t) gives the distance at time t
	MathFunction mathfunc_v; //evaluate(t) gives the velocity at time t

	double IO_dist;
	double IO_vel;

	//EDC	Vector3D col; //$EDC$[varaccess,remove,EDCvarname="RGB_color",EDCfolder="Graphics"] //color for body
	//EDC	Vector3D col; //$EDC$[varaccess,EDCvarname="color_body1",EDCfolder="Graphics",tooltiptext="[red, green, blue] first color of constraint, range = 0..1, use default color:[-1,-1,-1]"] //color for constraint
	//EDC	Vector3D col_ext; //$EDC$[varaccess,EDCvarname="color_body2",EDCfolder="Graphics",tooltiptext="[red, green, blue] second color of constraint, range = 0..1, use default color:[-1,-1,-1]"] //color for constraint

	// remove unused components from class
  //EDC int use_penalty_formulation;				//$EDC$[varaccess,remove,EDCvarname="use_penalty_formulation",EDCfolder="Physics"]
	//EDC int stiffness_in_joint_local_frame;	//$EDC$[varaccess,remove,EDCvarname="use_joint_local_frame",EDCfolder="Geometry"]
	//EDC TArray<int> dir;												//$EDC$[varaccess,remove,EDCvarname="constrained_directions",EDCfolder="Physics.Lagrange"]
	//EDC Vector3D spring_stiffness3;					//$EDC$[varaccess,remove,EDCvarname="spring_stiffness_vector",EDCfolder="Physics.Penalty"]
	//EDC Matrix3D JA0i;											//$EDC$[varaccess,remove,EDCvarname="joint_local_frame",EDCfolder="Geometry"]
	//EDC int use_local_coordinate_system;	//$EDC$[varaccess,remove,EDCvarname="use_local_coordinate_system",EDCfolder="Geometry"]
	//EDC double draw_dim(1);						//$EDC$[varaccess,remove,EDCvarname="draw_size",EDCfolder="Graphics"]
	//EDC double draw_local_frame_size;			//$EDC$[varaccess,remove,EDCvarname="draw_size_joint_local_frame",EDCfolder="Graphics"]
	//EDC double spring_stiffness; //$EDC$[varaccess,remove,EDCvarname="spring_stiffness",EDCfolder="Physics.Penalty"]
	//EDC double damping_coeff; //$EDC$[varaccess,remove,EDCvarname="damping",EDCfolder="Physics.Penalty"]

	//EDC double draw_dim(1);						//$EDC$[varaccess,EDCvarname="cylinder1_diameter",EDCfolder="Graphics",minval=0,tooltiptext="cylinder one diameter (drawing only)."]
	//EDC double draw_dim(2);						//$EDC$[varaccess,EDCvarname="cylinder2_diameter",EDCfolder="Graphics",minval=0,tooltiptext="cylinder two diameter (drawing only). Only used if distance not constant = distancemode 1 or 2."]
	//EDC double draw_dim(3);						//$EDC$[varaccess,EDCvarname="cylinder1_length",EDCfolder="Graphics",minval=0,tooltiptext="cylinder one length (drawing only). Only used if distance not constant = distancemode 1 or 2."]

};//$EDC$[endclass,RigidLink]


class RotatorySpringDamperActuator: public BaseBodyJoint  //$EDC$[beginclass,classname=RotatorySpringDamperActuator,parentclassname=BaseBodyJoint,addelementtype=TAEconstraint,addelementtypename=RotatorySpringDamperActuator,figure="rotatorySpringDamperActuator",texdescription="
//The RotatorySpringDamperActuator connects two elements with rotatory spring, damper and a constant actuator moment ma. Positive rotation around rotation axis according to right hand rule. 
//There are different modes available, how the spring and damper moment is calculated. It is also possible to change the neutral spring angle. This joint is realized in Penalty formulation only.",
//texdescriptionLimitations="The RotatorySpringDamperActuator should be used together with a RevoluteJoint to avoid useless simulation results. It is important to ensure that the relative angle of rotation between the two bodies must never be greater than $\pm\pi$. This has to be taken into accout when using an offset angle $\phi_0$.",
//modus="{element to ground}{Position2.element\_number AND Position2.node\_number have to be equal to 0}",
//modus="{element to element}{Position2.element\_number and/or Position2.node\_number must not be equal to 0}",
//modus="{forcemode}{
//\textbf{Physics.forcemode = 0:} \newline Moment is computed as (a) with constant stiffness and damping factors $k$ and $d$. The factors can be defined in the two fields in Physics.Linear. \newline 
//\textbf{Physics.forcemode = 1:} \newline A MathFunction is used to describe piecewise linear stiffness $k\left(\Delta \phi\right)$ and damping $d\left(\omega\right)$, see formula (b) and Physics.MathFunction. \newline 
//\textbf{Physics.forcemode = 2:} \newline 2 IOElementDataModifiers describe the moment (c) due to stiffness and damping. You should use this mode if full nonlinear behavior is required, e.g. $m_k = m_k\left(t,\phi,\omega,...\right)$ and $m_d=d\left(t,\phi,\omega,... \right)$.}",
//modus="{spring angle offset}{It is possible to change the spring angle $\phi_0$ (neutral angle of the spring) during the simulation, e.g. for the usage of the RotatorySpringDamperActuator as a rotational actuator. In 
//standard mode the offset remains constant. The value can be defined in the field Physics.spring\_angle\_offset. This offset can be modified by a IOElementDataModifier via 'Connector.angle\_offset'.}",
//modus="{additional actuator moment}{In Physics.actuator\_torque a constant offset moment $m_a$ can be added.}",
//texdescriptionEquations="spring angular deflection $\Delta \phi = \phi-\phi_0$ \\ spring angular velocity $\omega$ \\ \\ resultant moment (see section forcemode): \\ forcemode 0: $m = k\,\Delta \phi+d\,\omega+m_a$ (a) \\ forcemode 1: $m = k\left(\Delta \phi\right)\Delta \phi+d\left(\omega\right)\omega+m_a$ (b) \\ forcemode 2: $m = m_k+m_d+m_a$ (c)",example="RotationalSpringDamperActuator.txt"]
{
public:

	RotatorySpringDamperActuator(MBS* mbsi):BaseBodyJoint(mbsi)
	{	
		ElementDefaultConstructorInitialization();
	};

	RotatorySpringDamperActuator(const RotatorySpringDamperActuator& ct): BaseBodyJoint(ct.mbs)
	{
		ElementDefaultConstructorInitialization();
		CopyFrom(ct);
	};

	~RotatorySpringDamperActuator()
	{
	};

	virtual void ElementDefaultConstructorInitialization();

	virtual void Initialize();

	virtual Element* GetCopy()
	{
		Element* ec = new RotatorySpringDamperActuator(*this);
		return ec;
	}

	virtual void CopyFrom(const Element& e);

	virtual const char* GetElementSpec() const {return "RotatorySpringDamperActuator";}

	// for constraints only, these are the necessary (!) access functions!	//$ MSax 2013-02-18
	virtual void GetNecessaryKinematicAccessFunctions(TArray<int> &KAF, int numberOfKinematicPair);

	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		
	virtual int SetElementDataAuto(ElementDataContainer& edc);
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

  virtual int GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //for special value sensor
	virtual int ReadSingleElementData(ReadWriteElementDataVariableType& RWdata); 		//for special value sensor
	virtual int WriteSingleElementData(const ReadWriteElementDataVariableType& RWdata);

	virtual Vector3D ComputeForce(double t) const;
	virtual Vector3D ComputeMoment(double t) const;

	virtual void DrawElement();

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

protected:
	double phi0; //$EDC$[varaccess,EDCvarname="spring_angle_offset",EDCfolder="Physics",tooltiptext="spring angle offset is used if constant_spring_angle_offset is enabled. A positive offset equates a positve angle about the rotation axis."]
	double ma; //$EDC$[varaccess,EDCvarname="actuator_torque",EDCfolder="Physics",tooltiptext="constant torque of an actuator. A positive torque is acting about the rotation axis in a positive sense."]
	Vector3D glob_rot_axis; //$EDC$[varaccess,EDCvarname="rotation_axis",EDCfolder="Physics",tooltiptext="local axis of rotation w.r.t. body 1 coordinate system in inital configuration"]
	int forcemode; //$EDC$[varaccess,EDCvarname="forcemode",EDCfolder="Physics",tooltiptext="defines how the spring and damper moment is computed: 0..constant coefficient, 1..MathFunction, 2..IOElementDataModifier"]
	
	MathFunction mathfunc_k; //evaluate(x) gives the stiffness as function of x
	MathFunction mathfunc_d; //evaluate(x) gives damping as function of xp

	double moment_k;
	double moment_d;

	//EDC double spring_stiffness; //$EDC$[varaccess,EDCvarname="spring_stiffness",minval=0,EDCfolder="Physics.Linear",tooltiptext="stiffness parameter of the rotatory spring. Only used if forcemode is 0."]
	//EDC double damping_coeff; //$EDC$[varaccess,EDCvarname="damping",minval=0,EDCfolder="Physics.Linear",tooltiptext="damping coefficient for viscous damping. Only used if forcemode is 0."]

	// remove unused components from class
  //EDC int use_penalty_formulation;				//$EDC$[varaccess,remove,EDCvarname="use_penalty_formulation",EDCfolder="Physics"]
	//EDC int stiffness_in_joint_local_frame;	//$EDC$[varaccess,remove,EDCvarname="use_joint_local_frame",EDCfolder="Geometry"]
	//EDC TArray<int> dir;												//$EDC$[varaccess,remove,EDCvarname="constrained_directions",EDCfolder="Physics.Lagrange"]
	//EDC Vector3D spring_stiffness3;					//$EDC$[varaccess,remove,EDCvarname="spring_stiffness_vector",EDCfolder="Physics.Penalty"]
	//EDC Matrix3D JA0i;											//$EDC$[varaccess,remove,EDCvarname="joint_local_frame",EDCfolder="Geometry"]
	//EDC int use_local_coordinate_system;	//$EDC$[varaccess,remove,EDCvarname="use_local_coordinate_system",EDCfolder="Geometry"]
	//EDC double draw_dim(1);						//$EDC$[varaccess,remove,EDCvarname="draw_size",EDCfolder="Graphics"]
	//EDC double draw_local_frame_size;			//$EDC$[varaccess,remove,EDCvarname="draw_size_joint_local_frame",EDCfolder="Graphics"]
	//EDC TArray<int> dirRot; //$EDC$[varaccess,remove,EDCvarname="constrained_rotations",EDCfolder="Physics.Lagrange"]
	//EDC Matrix3D penaltyStiffness; //$EDC$[varaccess,remove,EDCvarname="stiffness_matrix",EDCfolder="Physics.Penalty"]
	//EDC Matrix3D penaltyDamping; //$EDC$[varaccess,remove,EDCvarname="damping_matrix",EDCfolder="Physics.Penalty"]
	//EDC Matrix3D penaltyStiffnessRot; //$EDC$[varaccess,remove,EDCvarname="stiffness_matrix_rotation",EDCfolder="Physics.Penalty"]
	//EDC Matrix3D penaltyDampingRot; //$EDC$[varaccess,remove,EDCvarname="damping_matrix_rotation",EDCfolder="Physics.Penalty"]
	//EDC Vector3D bryantAngles;	//$EDC$[varaccess,remove,EDCvarname="joint_local_frame_in_bryant_angles",EDCfolder="Geometry"]
	//EDC double draw_cone_size; //$EDC$[varaccess,remove,EDCvarname="cone_size",EDCfolder="Graphics"]

	//drawdim: X=spring drawsize, Y=windings, Z=axis radius (cylinder)
	//EDC double draw_dim(1);						//$EDC$[varaccess,EDCvarname="spring_size",EDCfolder="Graphics",minval=0,tooltiptext="radius of torsional spring. This parameter is used for drawing only."]
	//EDC double draw_dim(2);						//$EDC$[varaccess,EDCvarname="windings",EDCfolder="Graphics",minval=0,tooltiptext="number of windings of torsional spring. This parameter is used for drawing only."]
	double spring_res;									//$EDC$[varaccess,EDCvarname="spring_resolution",EDCfolder="Graphics",minval=1,maxval=10,tooltiptext="spring resolution used for drawing (very coarse = 1, very smooth = 10)."]
	//EDC double draw_dim(3);						//$EDC$[varaccess,EDCvarname="axis_radius",EDCfolder="Graphics",minval=0,tooltiptext="radius of torsional spring axis (cylinder). This parameter is used for drawing only."]

};//$EDC$[endclass,RotatorySpringDamperActuator]




class SpringDamperActuator2D: public Constraint  //$EDC$[beginclass,classname=SpringDamperActuator2D,parentclassname=Constraint,addelementtype=TAEconstraint,addelementtypename=SpringDamperActuator2D,
//texdescription="The SpringDamperActuator2D is a simplified version of the SpringDamperActuator for 2D elements. Nodes are not supported in the 2D version. Apart from that the constraint has the same functionality as the 3D version. See  the SpringDamperActuator documentation for more information.",
//modus="{element to ground}{Position2.element\_number has to be equal to 0}",
//modus="{element to element}{Position2.element\_number must not be equal to 0}",
//modus="{Lagrange}{If Physics.use\_penalty\_formulation = 0, than no stiffness and no damping parameters are used.}",
//example="SpringDamperActuator2D.txt"]
{
public:
	SpringDamperActuator2D(MBS* mbsi):Constraint(mbsi)
	{	
		ElementDefaultConstructorInitialization();
	};

	SpringDamperActuator2D(const SpringDamperActuator2D& ct): Constraint(ct.mbs)
	{
		ElementDefaultConstructorInitialization();
		CopyFrom(ct);
	};

	~SpringDamperActuator2D()
	{
	};

	virtual void Initialize() 
	{
	};

	virtual void ElementDefaultConstructorInitialization();

	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute

	virtual void GetNecessaryKinematicAccessFunctions(TArray<int> &KAF, int numberOfKinematicPair);

	virtual Element* GetCopy()
	{
		Element* ec = new SpringDamperActuator2D(*this);
		return ec;
	}

	virtual void CopyFrom(const Element& e);

	virtual const char* GetElementSpec() const {return "SpringDamperActuator2D";}

	virtual int NE() const {
		int counter = 0;
		for(int i=1; i<=elements.Length(); i++)
		{
			if(elements(i)!=0) counter++;
		}
		return counter;
	}

	virtual int IS() const {return 0;};  //implicit (algebraic) size

	virtual int SOS() const;
	virtual int SOSowned() const {return 0;}  // number of explicit equations added by element

	virtual void EvalF2(Vector& f, double t); //second order equations: M \ddot u = F2(u,\dot u,t), len(u)

	//virtual void EvalF2(Vector& f, double t) {}; //second order equations: M \ddot u = F2(u,\dot u,t), len(u)

	virtual Vector2D ComputeForce(double t) const;

	virtual int Dim() const {return 2;}

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		
	virtual int SetElementDataAuto(ElementDataContainer& edc);
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

  virtual int GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //for special value sensor
	virtual int ReadSingleElementData(ReadWriteElementDataVariableType& RWdata); 		//for special value sensor
	virtual int WriteSingleElementData(const ReadWriteElementDataVariableType& RWdata);

	virtual void SetStiffnessDamping();
	virtual void SetPiecewiseLinearSpring(const Vector& xpos, const Vector& force);
	virtual void SetPiecewiseLinearDamper(const Vector& vel, const Vector& force);

	virtual void DrawElement();

	virtual Vector2D ComputeForceDirection() const; 
	virtual Vector3D ComputeForceDirectionD(); 

	Vector2D GetPosition(int i) const;
	Vector2D GetVelocity(int i) const;
	Vector3D GetDrawPosition(int i) const;

	virtual double GetDampingCoeff() const {return damping_coeff;}
	virtual void SetDampingCoeff(double damp) {	damping_coeff = damp;	}
	virtual int UseDamping() const {
		if(damping_coeff!=0.) {return 1;}
		else {return 0;};
	}
protected:
	//  === data ===
	//EDC int elements(1);							//$EDC$[varaccess,EDCvarname="element_number",minval=1,EDCfolder="Position1",tooltiptext="Number of constrained element"]
	//EDC int elements(2);							//$EDC$[varaccess,EDCvarname="element_number",minval=0,EDCfolder="Position2",tooltiptext="Number of constrained element (0 if ground joint)"]

	TArray<Vector2D> loccoords;					// local coordinates of point in element, or (if ground joint) global coordinate of ground
	//EDC Vector2D loccoords(1);				//$EDC$[varaccess,EDCvarname="position",EDCfolder="Position1",tooltiptext="local position 1"]
	//EDC Vector2D loccoords(2);				//$EDC$[varaccess,EDCvarname="position",EDCfolder="Position2",tooltiptext="local or global position 2"]

	Matrix dpdq;

	//EDC	Vector3D col; //$EDC$[varaccess,remove,EDCvarname="RGB_color",EDCfolder="Graphics"] //color for body
	//EDC	Vector3D col; //$EDC$[varaccess,EDCvarname="color_body1",EDCfolder="Graphics",tooltiptext="[red, green, blue] first color of constraint, range = 0..1, use default color:[-1,-1,-1]"] //color for constraint
	//EDC	Vector3D col_ext; //$EDC$[varaccess,EDCvarname="color_body2",EDCfolder="Graphics",tooltiptext="[red, green, blue] second color of constraint, range = 0..1, use default color:[-1,-1,-1]"] //color for constraint

	// remove unused components from class
  //EDC int use_penalty_formulation;				//$EDC$[varaccess,remove,EDCvarname="use_penalty_formulation",EDCfolder="Physics"]
	//EDC int use_local_coordinate_system;	//$EDC$[varaccess,remove,EDCvarname="use_local_coordinate_system",EDCfolder="Geometry"]

	double l0; //$EDC$[varaccess,EDCvarname="spring_length",EDCfolder="Physics",tooltiptext="length of the spring in the initial configuration"]
	double fa; //$EDC$[varaccess,EDCvarname="actor_force",EDCfolder="Physics",tooltiptext="constant force acting on the spring"]
	int forcemode; //$EDC$[varaccess,EDCvarname="forcemode",EDCfolder="Physics",tooltiptext="defines how the spring and damper force is computed: 0..constant coefficient, 1..MathFunction, 2..IOElementDataModifier"]
	
	MathFunction mathfunc_k; //evaluate(x) gives the stiffness as function of x
	MathFunction mathfunc_d; //evaluate(x) gives damping as function of xp

	double force_k;
	double force_d;

	//EDC double spring_stiffness; //$EDC$[varaccess,remove,EDCvarname="spring_stiffness",EDCfolder="Physics.Penalty"]
	//EDC double spring_stiffness; //$EDC$[varaccess,EDCvarname="spring_stiffness",minval=0,EDCfolder="Physics.Linear",tooltiptext="stiffness coefficient of the linear spring. Only used if forcemode is 0."]

	double damping_coeff;		//$EDC$[varaccess,EDCvarname="damping",minval=0,EDCfolder="Physics.Linear",tooltiptext="damping coefficient for viscous damping. Only used if forcemode is 0."]

	//EDC double draw_dim(1);						//$EDC$[varaccess,EDCvarname="spring_diameter",EDCfolder="Graphics",minval=0,tooltiptext="spring diameter used for drawing only."]
	//EDC double draw_dim(2);						//$EDC$[varaccess,EDCvarname="spring_coils",EDCfolder="Graphics",minval=0,tooltiptext="spring coils used for drawing. If set to 0, then a cylinder with the value 'spring_diameter' as diameter is shown instead of the coils."]
	double spring_res;									//$EDC$[varaccess,EDCvarname="spring_resolution",EDCfolder="Graphics",minval=1,maxval=10,tooltiptext="spring resolution used for drawing (very coarse = 1, very smooth = 10)."]
	//EDC double draw_dim(3);						//$EDC$[varaccess,EDCvarname="damper_diameter",EDCfolder="Graphics",minval=0,tooltiptext="damper diameter used for drawing only. If set to 0, then the damper is not shown. It's recommended to choose the value smaller then the spring diameter."]

};//$EDC$[endclass,SpringDamperActuator2D]

class PointJoint2D: public Constraint  //$EDC$[beginclass,classname=PointJoint2D,parentclassname=Constraint,addelementtype=TAEconstraint,addelementtypename=PointJoint2D,
//texdescription="The PointJoint2D is a simplified version of the PointJoint for 2D elements. It constrains two elements at at local element positions. If only one element is specified (second element 0), a ground PointJoint is realized. It provides both Lagrangian and penalty formulation.",
//modus="{element to ground}{Position2.element\_number has to be equal to 0}",
//modus="{element to element}{Position2.element\_number must not be equal to 0}",
//modus="{Lagrange}{If Physics.use\_penalty\_formulation = 0, than no stiffness and no damping parameters are used.}"
//,example="PointJoint2D.txt"]
{
public:
	PointJoint2D(MBS* mbsi):Constraint(mbsi)
	{	
		ElementDefaultConstructorInitialization();
	};

	PointJoint2D(const PointJoint2D& ct): Constraint(ct.mbs)
	{
		ElementDefaultConstructorInitialization();
		CopyFrom(ct);
	};

	~PointJoint2D()
	{
	};

	virtual void Initialize() 
	{
		// if general (scalar) penalty stiffness is used, set Stiffness-Vector
		if(GetPenaltyStiffness() && UsePenaltyFormulation())
		{
			SetPenaltyStiffness2(Vector2D(GetPenaltyStiffness()));
		}
	};

	virtual void ElementDefaultConstructorInitialization();

	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute

	virtual void GetNecessaryKinematicAccessFunctions(TArray<int> &KAF, int numberOfKinematicPair);

	virtual Element* GetCopy()
	{
		Element* ec = new PointJoint2D(*this);
		return ec;
	}

	virtual void CopyFrom(const Element& e);

	virtual const char* GetElementSpec() const {return "PointJoint2D";}

	virtual int NE() const {
		int counter = 0;
		for(int i=1; i<=elements.Length(); i++)
		{
			if(elements(i)!=0) counter++;
		}
		return counter;
	}

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

	virtual int GetNumberOfConstrainedCoords() const 
	{
		int counter = 0;
		for(int i=1; i<=2; i++)
		{
			if(dir(i)!=0) counter++;
		}
		return counter;
	}	

	virtual int SOS() const;
	virtual int SOSowned() const {return 0;}  // number of explicit equations added by element

	virtual double GetDampingCoeff() const {return damping_coeff;}
	virtual void SetDampingCoeff(double damp) {	damping_coeff = damp;	}
	virtual int UseDamping() const {
		if(damping_coeff!=0.) {return 1;}
		else {return 0;};
	}

	//return 2 stiffness parameters for penalty formulation
	virtual Vector2D GetPenaltyStiffness2(double t) const;

	//set 2 stiffness parameters for penalty formulation
	virtual void SetPenaltyStiffness2(Vector2D stiffness2i) 
	{
		spring_stiffness2 = stiffness2i; 
		SetPenaltyFormulation(1);
	}

	virtual void EvalF2(Vector& f, double t); //second order equations: M \ddot u = F2(u,\dot u,t), len(u)
	virtual void EvalG(Vector& f, double t);
	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	virtual Vector2D ComputeForce(double t) const;

	virtual int Dim() const {return 2;}

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		
	virtual int SetElementDataAuto(ElementDataContainer& edc);
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

  virtual int GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //for special value sensor
	virtual int ReadSingleElementData(ReadWriteElementDataVariableType& RWdata); 		//for special value sensor

	virtual void DrawElement();

	virtual	void GetdRotTvdqT(const Vector2D& vloc, const Vector2D& ploc, Matrix& d, int bodyindex);

	virtual Matrix3D GetRotMati() const;
	virtual Matrix3D GetRotMatiP() const;
	virtual Matrix3D GetRotMatiD() const;

	Vector2D GetPosition(int i) const;
	Vector2D GetVelocity(int i) const;
	Vector3D GetDrawPosition(int i) const;

protected:
	//  === data ===
	//EDC int elements(1);							//$EDC$[varaccess,EDCvarname="element_number",minval=1,EDCfolder="Position1",tooltiptext="Number of constrained element"]
	//EDC int elements(2);							//$EDC$[varaccess,EDCvarname="element_number",minval=0,EDCfolder="Position2",tooltiptext="Number of constrained element (0 if ground joint)"]

	TArray<Vector2D> loccoords;					// local coordinates of point in element, or (if ground joint) global coordinate of ground
	//EDC Vector2D loccoords(1);				//$EDC$[varaccess,EDCvarname="position",EDCfolder="Position1",tooltiptext="local position 1"]
	//EDC Vector2D loccoords(2);				//$EDC$[varaccess,EDCvarname="position",EDCfolder="Position2",tooltiptext="local or global position 2"]

	Vector2D spring_stiffness2;					//$EDC$[varaccess,EDCvarname="spring_stiffness_vector",EDCfolder="Physics.Penalty",tooltiptext="penalty stiffness parameter [kx,ky]. Just used if scalar spring_stiffness == 0, otherwise kx=ky=spring_stiffness"]

	TArray<int> dir;											//$EDC$[varaccess,EDCvarname="constrained_directions",minval=0,maxval=1,EDCfolder="Physics.Lagrange",tooltiptext="[x,y]...(1 = constrained, 0 = free), can be defined as local or global directions (see Geometry)"]
	
	double damping_coeff;								//$EDC$[varaccess,EDCvarname="damping",minval=0,EDCfolder="Physics.Penalty",tooltiptext="damping coefficient for viscous damping (F = d*v), applied in all constrained directions"]

	Matrix dpdq;
	Matrix drotTvdq;
	Matrix tmpmat;
	

	int stiffness_in_joint_local_frame;	//$EDC$[varaccess,EDCvarname="use_joint_local_frame",EDCfolder="Geometry",int_bool,tooltiptext="Use a special joint local frame"]

	double phi_z;												//$EDC$[varaccess,EDCvarname="joint_local_frame",EDCfolder="Geometry",tooltiptext="Prerotate stiffness vector w.r.t. global coordinate system or local coordinate system of body 1 with angle phi_z about the z axis. Just used if use_joint_local_frame == 1"]

	double	draw_local_frame_size;			//$EDC$[varaccess,EDCvarname="draw_size_joint_local_frame",EDCfolder="Graphics",tooltiptext="drawing dimensions of joint local frame. If set to -1, than global_draw_scalar_size is used. If set to 0, than no joint local frame is drawn."]

	//EDC double draw_dim(1);						//$EDC$[varaccess,EDCvarname="draw_size",EDCfolder="Graphics",tooltiptext="drawing dimensions of constraint. If set to -1, than global_draw_scalar_size is used."]

	//EDC int use_local_coordinate_system;	//$EDC$[varaccess,remove,EDCvarname="use_local_coordinate_system",EDCfolder="Geometry"]
	//EDC int use_local_coordinate_system;	//$EDC$[varaccess,EDCvarname="use_local_coordinate_system",EDCfolder="Geometry",int_bool,tooltiptext="0=use global coordinates, 1=use local coordinate system of Body 1"]

};//$EDC$[endclass,PointJoint2D]



#endif
