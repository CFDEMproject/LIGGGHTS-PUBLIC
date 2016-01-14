//#***************************************************************************************
//#
//# filename:             element.h
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						July 2004
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

#ifndef ELEMENT__H
#define ELEMENT__H

#include "mbs_interface.h"

#include "FieldVariableDescriptor.h"
#include "mbsload.h"		// elements need to deal with loads

class ReadWriteElementDataVariableType;

class Constraint;


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                                  ELEMENT
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

typedef enum {TElement=1, TConstraint=2, TBody=4,/* TFFRF=8, TBody2D=16, TBody3D=32,*/
TCMS = 64, TCMSflag = 128, TController = 256, TFiniteElement = 512, TGCMS = 1024, TParticle = 2048, TIOElementDataModifier = 4096} TMBSElement;
// TCMSflag ... set if the Element is a BaseCMSElement
// TCMS ....... set if the Element is a finite element belonging to a BaseCMSElement (i.e. CMS *OR* GCMS)
// TGCMS ...... set if the Element is a finite element belonging to a GCMSElement
// TParticle... set if the Element is a particle (e.g., of a fluid, see the classes SPHParticle2D, SPHParticle3D, SPHMass2D, SPHMass3D)

//definition of element equations structure:
//this information is for element description and for decision of solvers which method to use
typedef enum {TET_algebraic=1, TET_algebraic_linear=2,
TET_first_order_ODE=4, TET_first_order_ODE_linear=8,
TET_second_order_ODE=16, TET_second_order_ODE_linear=32,
TET_discontinuous=64, TET_fixedpoint_iteration=128,
TET_delay_equation=256, TET_stochastic_process=512,
TET_constant_mass_matrix=1024, TET_constant_stiffness_matrix=2048,
TET_Lagrange_multipliers=4096	} TEquationType;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//new position, displacement, velocity ... access functions
typedef enum{TKAF_none=0,															//dummy for initialization
TKAF_position=1,																			//GetPos	
TKAF_displacement=2,																	//GetDisplacement
TKAF_ref_conf_position=131072,												//GetRefConfPos
TKAF_ref_position=4194304,														//GetRefPos
TKAF_node_position = 4,																//GetNodePos
TKAF_node_ref_conf_position = 262144,									//GetNodeRefPos
TKAF_global_to_local_position=8,
TKAF_velocity=16,																			//GetVel
TKAF_node_velocity=32,																//GetNodeVel
TKAF_acceleration=64,		
TKAF_angular_velocity=128,														//GetAngularVel
TKAF_rotation_matrix=256,															//GetRotMatrix
TKAF_rotation_matrix_P=512,														//GetRotMatrixP
TKAF_D_pos_D_q=1024,																	//GetdPosdqT
TKAF_D_pos_D_x=2048,																	//GetdPosdx
TKAF_D_node_pos_D_q = 4096,														//GetNodedPosdqT
TKAF_D_ang_vel_D_q =8192,															//GetdAngVeldqpT
TKAF_D_rot_D_q=16384,																	//GetdRotdqT
TKAF_D_rot_v_D_q=32768,																//GetdRotvdqT
TKAF_int_D_u_D_q =65536,															//GetIntDuDq
TKAF_position_2D =524288,															//GetPos2D
TKAF_velocity_2D =1048576,														//GetVel2D
TKAF_D_pos_D_q_2D =2097152,														//GetdPosdqT2D


//GetIntDuDqFCentrifugal
TKAF_maximum=8388608} TKinematicsAccessFunctions;		//always update maximum, such that it is still the maximum, used in Constraint:IsSuitableElement

//flag T_compute, T_draw, T_reference_configuration and T_initial_values decides, where to retrieve the system coordinate
	//flag T_draw_magnified decides, whether the displacment is magnified or not in the drawing function (e.g. GetPos_dc)
	//flag T_reference_configuration returns zero under the assumption that in reference configuration the coordinates are zero (especially for rigid bodies and classical displacement based finite elements)
	//flag  are not implemented yet
typedef enum {TCD_compute=1, TCD_draw=2, TCD_reference_configuration=4, TCD_initial_values=8, TCD_draw_magnified=16, TCD_cached=32} TComputeDrawInitFlag;

class Element //$EDC$[beginclass,classname=Element]
{

public:
	Element():ltg(), ltgdata(), constraintindices(), loads(), constraints(), constraints_nodouble(), drawelements(), sensors(), elements(), x_init(), data_init()
	{
		Element::ElementDefaultConstructorInitialization();
		//damping_m=0;
		//altshape = 0;
		//type = TElement;
		//elementname = GetElementSpec();
		//rho = 0;
		//materialnum = 0;
		//draw_element = 1;
	};
	Element(MBS* mbsi):ltg(), ltgdata(), constraintindices(), loads(), constraints(), constraints_nodouble(), drawelements(), sensors(), elements(), x_init(), data_init()
	{
		Element::ElementDefaultConstructorInitialization();
		mbs=mbsi;
		//damping_m=0;
		//altshape = 0;
		//type = TElement;
		//elementname = GetElementSpec();
		//rho = 0;
		//materialnum = 0;
		//draw_element = 1;
	};
	Element(const Element& e):ltg(), ltgdata(), constraintindices(), loads(), constraints(), constraints_nodouble(), drawelements(), sensors(), elements(), x_init(), data_init()
	{
		CopyFrom(e);
	};
	Element& operator=(const Element& e) 
	{
		if (this == &e) {return *this;}
		CopyFrom(e);
		return *this;
	}
	//To be overwritten in derived class:
	virtual Element* GetCopy();

	virtual ~Element();

	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e);


	//this function assigns default values to the element variables
	virtual void ElementDefaultConstructorInitialization()
	{
		damping_m=0;
		altshape = 0;
		type = TElement;
		elementname = GetElementSpec();
		//rho = 0;
		materialnum = 0;
		draw_element = 1;
		mass = 0;
		elnum = 0;
		col = Vector3D(0.1,0.1,0.8); //default body color
		//col = Vector3D(-1,-1,-1); //default body color

		x_init = Vector(2*SOS()+ES()+IS());   //initial conditions
		data_init = Vector(DataS());   //initial conditions for data variables
		
		loads.SetLen(0);	//$ DR 2012-10

		//not initilized: TArray<int> ltg; //local to global dof reference vector;
		//not initilized: TArray<int> ltgdata; //local to global data reference vector;
		//not initilized: TArray<int> constraintindices; //index of this element in the constraint (mostly 1 or 2)
		//not initilized: TArray<char> dependencies; //contains dependencies to coordinates
		//not initilized: TArray<Constraint*> constraints; //this constraint list might be double!
		//not initilized: TArray<Constraint*> constraints_nodouble; //one constraint is only added once to an element
		//not initilized: TArray<int> sensors; //EDC[varaccess,EDCvarname="sensors",EDCfolder="Info",tooltiptext="attached sensors",readonly]//if the element equations of motion are dependent on sensor values (e.g. as an input), then add the sensor numbers here
		//not initilized: TArray<int> elements; //dependent elements (mainly for contraints, but also for follower loads!)
		//not initilized: TArray<int> drawelements; 
	}

	virtual void Initialize() //initialize after initial conditions are set!!!
	{
	};

	virtual const MBS* GetMBS() const {return mbs;}
	virtual MBS* GetMBS() {return mbs;}
	//virtual UserOutput& UO() {return mbs->uout;};
	//virtual UserOutput& UO(int message_level = UO_LVL_all) { mbs->uout.SetLocalMessageLevel(message_level); return mbs->uout;}; //(AD)
	//$ YV 2012-11-28
	//virtual UserOutput& UO(int message_level = UO_LVL_all, int output_prec = -1) { mbs->uout.SetLocalMessageLevel(message_level); mbs->uout.SetOutputPrec(output_prec); return mbs->uout;}; //$ AD 2011-02 output_prec
	virtual UserOutputInterface& UO(int message_level = UO_LVL_all, int output_prec = -1) const { return mbs->UO(message_level,output_prec); }
	virtual int MaxIndex() const {return mbs->MaxIndex();}
	
	//$ JG 2011-02: compute index 3 constraint drift, only for postprocessing! (evaluate sensors)
	virtual double GetConstraintDrift(double t) const;

	virtual int GetOwnNum() const {return elnum;}
	virtual void SetOwnNum(int i) {elnum = i;}

	virtual const mystr& GetElementName() const {return elementname;}
	virtual mystr& GetElementName() {return elementname;}
	virtual void SetElementName(const char* name) {elementname = name;}
	virtual const char* GetElementSpec() const {return "Element";}	//$ DR this is implemented already in many elements and used at multiple places
	virtual mystr GetElementSpecification() //$EDC$[funcaccess,EDCvarname="element_type",tooltiptext="specification of element type. Once the element is added to the mbs, you MUST NOT change this type anymore!"]
	{
		mystr type = GetElementSpec();	//$ DR this is necessary for the skript language
		return type;
	}	

	virtual int IsType(TMBSElement te) const 
	{
		return (te&type) != 0;
	}
	virtual void SetType(TMBSElement te) {type = te;}
	
	virtual TMBSElement GetType() const {return type;} //return element type (e.g. for knowing if element is a constraint, etc.)
	//the following string returns a list of strings which contains main specification of the element
	//virtual MyStrList GetType_String() const
	//{
	//	MyStrList strlist;
	//	if (IsType(TElement)) {strlist.Add("connector element");}
	//	if (IsType(TConstraint)) {strlist.Add("connector element");}
	//	if (IsType(TBody) && !IsType(TFiniteElement)) {strlist.Add("body");}
	//	if (IsType(TFiniteElement)) {strlist.Add("finite element");}
	//	if (IsType(TCMS)) {strlist.Add("finite element belonging to a CMS reference frame");}
	//	if (IsType(TCMSflag)) {strlist.Add("reference frame for component mode synthesis (CMS)\n");}
	//	if (IsType(TGCMS)) {strlist.Add("GCMS modally reduced element\n");}
	//	if (IsType(TController)) {strlist.Add("control element (with input and output)\n");}
	//	if (IsType(TParticle)) {strlist.Add("particle for sph simulation\n");}

	//	return str;
	//}

	virtual void AddType(TMBSElement te) 
	{
		type = (TMBSElement)(type|te);
	}
	virtual int IsRigid() const {return 0;} //$EDC$[funcaccess,readonly,EDCvarname="is_rigid",EDCfolder="Info",tooltiptext="flag if body is a rigid body"]//default value
	virtual int IsFiniteElement() const {return 0;} //$EDC$[funcaccess,readonly,EDCvarname="is_finite_element",EDCfolder="Info",tooltiptext="flag if body is a finite element"] //if is of type "FiniteElement" 2D or 3D

	////return the information about the equations type of the element
	//virtual TElementEquationType GetEquationsType() const
	//{
	//	TElementEquationType eet = 0;
	//	if (IS() != 0) eet += TET_algebraic;
	//	if (ES() != 0) eet += TET_first_order_ODE;
	//	if (SOS() != 0) eet += TET_second_order_ODE;
	//	if (IS() != 0) eet += TET_Lagrange_multipliers;

	//	return eet;
	//}
	//return stringlist containing the possible types of the equation that apply
	//virtual MyStrList GetEquationsType_String() const
	//{
	//	MyStrList strlist;
	//	if ((GetEquationsType()&TET_algebraic != 0) && (GetEquationsType()&TET_algebraic_linear != 0)) {strlist.Add("linear algebraic equations");}
	//	else if ((GetEquationsType()&TET_algebraic != 0)) {strlist.Add("nonlinear algebraic equations");}

	//	if ((GetEquationsType()&TET_first_order_ODE != 0) && (GetEquationsType()&TET_first_order_ODE_linear != 0)) {strlist.Add("linear first order differential equations");}
	//	else if ((GetEquationsType()&TET_first_order_ODE != 0)) {strlist.Add("nonlinear first order differential equations");}

	//	if ((GetEquationsType()&TET_second_order_ODE != 0) && (GetEquationsType()&TET_second_order_ODE_linear != 0)) {strlist.Add("linear second order differential equations");}
	//	else if ((GetEquationsType()&TET_second_order_ODE != 0)) {strlist.Add("nonlinear second order differential equations");}

	//	if (GetEquationsType()&TET_discontinuous != 0) {strlist.Add("discontinuous equations (switches, jumps, etc.)");}
	//	if (GetEquationsType()&TET_fixedpoint_iteration != 0) {strlist.Add("fixed point iteration for discontinuous equations");}

	//	if (GetEquationsType()&TET_delay_equation != 0) {strlist.Add("delay equations (storage and retrieving of data)");}
	//	if (GetEquationsType()&TET_stochastic_process != 0) {strlist.Add("stochastic process (solution might be non-deterministic)");}
	//	if (GetEquationsType()&TET_constant_mass_matrix != 0) {strlist.Add("second order equations contain a constant mass matrix");}
	//	if (GetEquationsType()&TET_constant_stiffness_matrix != 0) {strlist.Add("second order equations contain a constant stiffness matrix");}
	//	if (GetEquationsType()&TET_Lagrange_multipliers != 0) {strlist.Add("Lagrange multipliers are used and forces added to other elements");}
	//}

#pragma endregion

	virtual int PerformNodeCheck() const {return 1;} //node number are not checked for CMS-elements ...
	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	
	//the following functions are for reading/writing the whole element data including the necessary initialization after setting the element data
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!


	//the following functions are for read/write access to single variables or functions (mapped to EDC-style-variables)
	//these functions should be implemented in all derived classes
	virtual int ReadSingleElementData(ReadWriteElementDataVariableType& RWdata); 		//retrieve value of element variable/function named RWdata.variable_name with specified components (vector/matrix); return 1, if variable found and read; 0 if not found, -1 if no readaccess, -2 if range-fault
	virtual int WriteSingleElementData(const ReadWriteElementDataVariableType& RWdata); //write value of element variable/function named RWdata.variable_name with specified components (vector/matrix); return 1, if variable found and written; 0 if not found, -1 if no writeaccess, -2 if range-fault
  virtual int GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //$ DR 2012-10

	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	//virtual int IsConstraint() const {return 0;}

	virtual void PrecomputeEvalFunctions() {}; //PG: precomputation for EvalF(), EvalF2(), EvalM(), and also for EvalG() and AddElementCqTLambda() in case of constraints
	virtual void EvalF(Vector& f, double t) {};  //first order equations: \dot q=f(q,t), return vector:  len(q)
	virtual void EvalG(Vector& f, double t) {};  //evaluate constraints: len(z)
	//second order equations: M \ddot u = F2(u,\dot u,t), len(u)
	virtual void EvalM(Matrix& m, double t) {};  //evaluate mass matrix
	//Element-elastic forces+external forces+Cq^T*lambda
	virtual void EvalF2(Vector& f, double t);
	virtual void EvalMinvF2(Vector& f, double t);
	//Jacobian for F2, compute element Jacobian m, LTGVector ref
	virtual void JacobianF2(double t, Matrix& m, IVector& colref);
	virtual void JacobianG(double t, Matrix& m, IVector& colref);
	virtual double GetKineticEnergy() {return 0;};   //compute kinetic energy of system
	virtual double GetPotentialEnergy(); //compute potential (strain) energy of system

	virtual int FastStiffnessMatrix() const {return 0;}
	virtual void StiffnessMatrix(Matrix& m) {assert(0 && "ERROR: StiffnessMatrix() not defined"); }; //fill in sos x sos components, m might be larger

	virtual void GyroscopicMatrix(SparseMatrix& gy)const  // $ MSax 2013-07-25 : added
	{
		// gyroscopic matrix for rotordynamic elements. Used e.g. for computation of eigenmodes (campbell)
		gy = SparseMatrix();
		gy.SetSize(SOS(),SOS());
		gy.FillWithZeros();
	};

	//to be activated, MBS->SetTransformJacApply must be set:
	virtual int TransformJacApply() const {return 0;}; //transform jacobian in solver (x = x_0 + A*J^(-1)*B * f)
	virtual void ApplyTransform(const Vector& v, Vector& Av, int mode) {}; //compute Av=A^T*v in mode==0 and Av=A*v in mode==1

	//sparse jacobian operations:
	virtual void AddMSparse(SparseMatrix& m, double t) {};  //add sparse mass matrix into full system matrix
	virtual void AddKSparse(SparseMatrix& m, double t) {};  //add sparse stiffness matrix into full system matrix
	virtual void AddDSparse(SparseMatrix& m, double t) {};  //add sparse damping matrix into full system matrix (at sosmbs)
	virtual int UseSparseM() const {return 0;};
	virtual int UseSparseK() const {return 0;};
	virtual void SetUseSparseM(int i=1) {};
	virtual void SetUseSparseK(int i=1) {};

	//floating frame of reference formulation: for FFRF elements
	virtual void GetI1(Vector& I1) {assert(0);}; //Shabana p. 209-211
	virtual double GetIkl(int k, int l) {assert(0);return 0;};
	virtual void GetIbarkl(int k, int l, Vector& I1) {assert(0);};
	virtual void GetSbar(Matrix& Sbar) {assert(0);};
	virtual void GetSbarkl(int k, int l, Matrix& Sbar) {assert(0);};
	virtual void EvalMff(Matrix& m, double t) {assert(0);}; //matrix without rigid body motion, only flexible part
	virtual void GetH(Matrix& H)  {assert(0);}; 
	virtual int NCMSNodes() const {mbs->UO() << "Element::NCMSNodes() called, should be called for FFRFElement only!!\n"; return 0;}

	//dir: normal == 0=X, 1=-X, 2=Y, 3=-Y, 4=Z, 5=-Z
	virtual void AddSurfacePressure(Vector& f, double pressure, int dir) {assert(0);}
	virtual Vector3D GetSurfaceNormalD(int dir) { assert(0); return Vector3D(); }		// AH: get normal to a specfied surface in order to plot surface pressure

	virtual int IS() const {return 0;};  //$EDC$[funcaccess,readonly,EDCvarname="n_algebraic_equations",EDCfolder="Info",tooltiptext="number of algebraic equations"] //implicit (algebraic) size
	virtual int SOS() const {return 0;}; //$EDC$[funcaccess,readonly,EDCvarname="n_second_order_ODEs",EDCfolder="Info",tooltiptext="number of second order ODEs in which the element contributes with some terms to the mass matrix and RHS"] //size of second order components (size of M and EvalF2=^=K)
	virtual int SOSowned() const {return SOS();}; //$EDC$[funcaccess,readonly,EDCvarname="n_second_order_ODE_unknowns",EDCfolder="Info", tooltiptext="number of second order ODE unknowns, which are owned (caused) by one element (not borrowed e.g. from another element or node)] //size of second order unknowns, len(u) or len(v)
	virtual int ES() const {return 0;};  //$EDC$[funcaccess,readonly,EDCvarname="n_first_order_ODEs",EDCfolder="Info",tooltiptext="number of first order ODEs"] //size of first order explicit equations
	virtual int SS() const {return 2*SOS()+ES()+IS();};  //$EDC$[funcaccess,readonly,EDCvarname="element_equation_size",EDCfolder="Info",tooltiptext="element system size=2*SOS()+ES()+IS()"]//system size
	virtual int DataS() const {return 0;} //$EDC$[funcaccess,readonly,EDCvarname="data_size",EDCfolder="Info",tooltiptext="number of data variables (e.g. plastic strains), for which there are no separate equations"] //Data size for non-state variables (contact condition, plastic strains, discrete states, etc.)

  virtual int IS_RS() const {return IS();};  //implicit (algebraic) size
  virtual int SOSowned_RS() const {return SOSowned();}; //for resort, the implicit element variables belong to SOS2
  virtual int FlexDOF() const {return SOS();}

   //# Timeint specific derived functions: for discontinuities
  virtual void StartTimeStep() {}; //function is called when computation of time step is started
  virtual void EndTimeStep() {}; //function is called when computation of time step is started
  virtual void ComputationFinished() {}; //function is called when computation is finished (e.g. in order to free memory, write results, close connections, etc.)

	//$!LA 2011-02-16: NonlinStep was replaced by PostNewtonStep
	// virtual double NonlinStep(double t) {return 0;};
  virtual double PostNewtonStep(double t) {return 0;};
	//$!LA 2011-02-16: FixNonlinStep was replaced by PostprocessingStep
	// virtual void FixNonlinStep() {};
	virtual void PostprocessingStep() {};
	virtual double GetError() const 
	{
		double err = 0;
		for (int i=1; i<=SS()-IS(); i++)
		{
			//err += fabs(XG(i));
			err += Sqr(XG(i));
		}
		return err;
	};

	virtual void SetGlobalInitConditions(Vector& x_glob)
	{
		for (int i=1; i<=SS(); i++)
		{
// (AD) changed () to .Get()
			x_glob(ltg.Get(i)) = x_init(i);
//			x_glob(ltg(i)) = x_init(i);
		}
	}

	virtual void SetGlobalInitData(Vector& data_glob)
	{
		for (int i=1; i<=DataS(); i++)
		{
// (AD) changed () to .Get()
			data_glob(ltgdata.Get(i)) = data_init(i);
//			data_glob(ltgdata(i)) = data_init(i);
		}
	}

	virtual int Dim() const {return 2;} //$EDC$[funcaccess,readonly,EDCvarname="element_dimension",EDCfolder="Info",tooltiptext="dimensionality of element (2D/3D)"]//default value
	virtual const double& GetXact(int i) const {return mbs->GetXact(i);}
	virtual double& GetXact(int i) {return mbs->GetXact(i);}
	//virtual const Vector& GetXact() const {return mbs->GetXact();}
	//virtual Vector& GetXact() {return mbs->GetXact();}
	//
	virtual void SetInitConditions(const Vector& x0) {x_init = x0;}
	virtual const Vector& GetXInit() const {return x_init;}

	virtual const double& GetXInit(int i) const {return x_init.Get(i);} // $ MSax 2013-07-12 : added
	
	virtual void SetDataInit(const Vector& data_initI) {data_init = data_initI;}
	virtual const Vector& GetDataInit() const {return data_init;}

	//compute things which should be actualized after every change of xact
	virtual void Update() {};
	//local to global transformation, for positions, velocities, first order ODE and algebraic variables

	virtual double GetOutput(double t, int i=1) const {return 0;}

	virtual int ElementBandwidth() const {return SS();}


	//local to global transformation variables
// (AD) changed () to .Get() & Elem()
	virtual int LTG(int iloc) const {return ltg.Get(iloc);}
	virtual int& LTG(int iloc) {return ltg.Elem(iloc);}
//	virtual int LTG(int iloc) const {return ltg(iloc);}
//	virtual int& LTG(int iloc) {return ltg(iloc);}
	virtual const TArray<int>& GetLTGArray() const {return ltg;} 
	virtual void AddLTG(int gi)	{ltg.Add(gi);}
	virtual int LTGlength() const {return ltg.Length();}
	virtual void LTGreset() {ltg.Flush();}

	//local to global transformation variables for data elements:
// (AD) changed () to .Get() & Elem()
	virtual int LTGdata(int iloc) const {return ltgdata.Get(iloc);}
	virtual int& LTGdata(int iloc) {return ltgdata.Elem(iloc);}
//	virtual int LTGdata(int iloc) const {return ltgdata(iloc);}
//	virtual int& LTGdata(int iloc) {return ltgdata(iloc);}
	virtual const TArray<int>& GetLTGdataArray() const {return ltgdata;} 
	virtual void AddLTGdata(int gi)	{ltgdata.Add(gi);}
	virtual void LTGdataReset() {ltgdata.Flush();}

	//tell loads the element reference pointer
	virtual void LinkLoads()
	{
		for (int i=1; i <= loads.Length(); i++)
		{
			//loads(i)->SetElement(this); //$ DR 2012-10: loads moved from element to mbs, old code
			//mbs->GetLoad(loads(i)).SetElement(this); //$ DR 2012-10: loads moved from element to mbs, old code
			mbs->GetLoad(loads(i)).SetMBS(mbs); //$ DR 2012-10: loads moved from element to mbs
			mbs->GetLoad(loads(i)).SetDim(this->Dim()); //$ DR 2012-10: loads moved from element to mbs
		}
	}
	//Link element to other elements if necessary (e.g. constraints)
	virtual void LinkToElements() {};

	virtual void BuildDependencies(); 

	//$ YV 2011-08-02:
	// the following function can be implemented by elements,
	// which wish to perform some actions just before the system is assembled
	virtual void PreAssemble() {}

	virtual int IsDependent(int i) const 
	{
		if (i == 0 || !GetMBS()->UseDependencies()) return 1;
#ifdef _DEBUG
		if (i > dependencies.Length()) 
		{
			mbs->UO() << "Error in IsDependent!!!!\n";
			mbs->UO() << "deplen=" << dependencies.Length() << ", i=" << i << "\n";
			return 1;
		}
#endif
		return dependencies(i);
	};

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//get system coordinate of actual state (must be set in EvalF, EvalG, EvalF2, EvalM)
	virtual const double XG_dc(int iloc, TComputeDrawInitFlag flag) const 
	{
		if(flag & TCD_compute)
		{
			return GetXact(ltg.Get(iloc));
		}
		else if(flag & TCD_draw)
		{
			return GetDrawValue(ltg.Get(iloc));
		}
		else if(flag & TCD_reference_configuration)
		{
			return 0;
		}
		else if(flag & TCD_initial_values)
		{
			return GetXInit(iloc);
		}
		else 
		{
			assert("XG_dc called with illegal flag" && 0);
			return 0;
		}
	}


	//get system coordinate of actual state (must be set in EvalF, EvalG, EvalF2, EvalM)
	virtual const double& XG(int iloc) const {return GetXact(ltg.Get(iloc));}
	//EDC Vector XG(); //$EDC[funcaccess,readonly,EDCvarname="SOS_ODE_GC_poslevel",EDCfolder="Debug",condition=SOS()>0,vecstart=1,vecend=SOS(),tooltiptext="second order size state variables at position level (XG(i))"]//default value
	virtual double& XG(int iloc) {return GetXact(ltg(iloc));}
	virtual const double& XGP(int iloc) const {return GetXact(ltg.Get(iloc+SOS()));}
	//EDC Vector XGP(); //$EDC[funcaccess,readonly,EDCvarname="SOS_ODE_GC_vellevel",EDCfolder="Debug",condition=SOS()>0,vecstart=1,vecend=SOS(),tooltiptext="second order size state variables at velocity level (XGP(i))"]//default value
	virtual double& XGP(int iloc) {return GetXact(ltg(iloc+SOS()));}
	//EDC Vector XG(); //$EDC[funcaccess,readonly,EDCvarname="FOS_ODE_GC",EDCfolder="Debug",condition=ES()>0,vecstart=2*SOS()+1,vecend=2*SOS()+ES(),tooltiptext="first order size state variables (XG(2*SOS()+i))"]//default value
	//EDC Vector XG(); //$EDC[funcaccess,readonly,EDCvarname="AE_lagrange_multipliers",EDCfolder="Debug",condition=IS()>0,vecstart=2*SOS()+ES()+1,vecend=SS(),tooltiptext="first order size state variables (XG(2*SOS()+ES()+i))"]//default value
	//virtual const double& XGPP(int iloc) const {return mbs->GetVelocityAndAcceleration(ltg.Get(iloc+SOS()));} //this retrieves the acceleration values; only with implicit time integration! // $ MSax 2013-07-16 : removed
	virtual double XGPP(int iloc) const  // $ MSax 2013-07-16 : added, returns the acceleration in dynamic case, returns zero in static case
	{
		if(mbs->DoStaticComputation()) return 0;
		double xgpp = mbs->GetVelocityAndAcceleration(ltg.Get(iloc+SOS())); //this retrieves the acceleration values; only with implicit time integration!
		return xgpp;
	}


	//global functions
	virtual const double& XGG(int iglob) const {return GetXact(iglob);}
	virtual double& XGG(int iglob) {return GetXact(iglob);}

	//functions with values for drawing:
	//get coordinate of actual state (must be set in EvalF, EvalG, EvalF2, EvalM)
// (AD) changed () to .Get()
	virtual const double& XGD(int iloc) const {return GetDrawValue(ltg.Get(iloc));}
	virtual double& XGD(int iloc) {return GetDrawValue(ltg.Get(iloc));}
//	virtual const double& XGD(int iloc) const {return GetDrawValue(ltg(iloc));}
	//virtual double& XGD(int iloc) {return GetDrawValue(ltg(iloc));}
// (AD) changed () to .Get()
	virtual const double& XGPD(int iloc) const {return GetDrawValue(ltg.Get(iloc+SOS()));}
///	virtual const double& XGPD(int iloc) const {return GetDrawValue(ltg(iloc+SOS()));}
	//virtual double& GetDrawValue(int iloc) {return mbs->GetDrawValue(iloc);}
	virtual const double& GetDrawValue(int iloc) const {return mbs->GetDrawValue(iloc);}
	virtual double& GetDrawValue(int iloc) {return mbs->GetDrawValue(iloc);}

// (AD) changed () to .Get()
	virtual const double& XData(int iloc) const {return GetMBS()->GetDataAct(ltgdata.Get(iloc));}
//	virtual const double& XData(int iloc) const {return GetMBS()->GetDataAct(ltgdata(iloc));}
	virtual double& XData(int iloc) {return GetMBS()->GetDataAct(ltgdata(iloc));}

// (AD) changed () to .Get()
	virtual const double& XDataD(int iloc) const {return GetMBS()->GetDataDraw(ltgdata.Get(iloc));}
//	virtual const double& XDataD(int iloc) const {return GetMBS()->GetDataDraw(ltgdata(iloc));}
	virtual double& XDataD(int iloc) {return GetMBS()->GetDataDraw(ltgdata(iloc));}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//Get rot matrix in 2D and 3D
	//rot matrix A describes the transformation of the local(body) to the
	//global coordinate system, such that r_g=R_g+A_bg*u_b
	virtual Matrix3D GetRotMatrix() const {return Matrix3D(1);};
	virtual Matrix3D GetRotMatrixP() const {return Matrix3D();};
	virtual Matrix3D GetRotMatrix(const Vector3D& ploc) const {return Matrix3D(1);};
	//virtual const Matrix& GetRotMatrix2D() {return rot;};

	//Get roation matrix in 2D or 3D
	virtual Matrix3D GetRotMatrixD() const {return Matrix3D(1);};

	virtual Matrix3D GetRotMatrixInit() const {return Matrix3D(1);};
	//virtual const Matrix& GetRotMatrix2DD() {return rot;};
	//get a reference position of the body in 3d

	// in specialized classes: keep this function up to date for consistency check with constraints, see Constraint::IsSuitableElement(..)
	// add all kinematic access functions available to those of the base class
	virtual TKinematicsAccessFunctions GetKinematicsAccessFunctions(int mode = 1) const
	{
		return TKinematicsAccessFunctions(TKAF_position+TKAF_displacement+TKAF_velocity);
	}


	//virtual Vector3D GetPos_dc(const Vector3D& ploc, TComputeDrawInitFlag flag) const 
	//{
	//	return GetRefPos_dc(flag)+ploc;
	//}

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	virtual Vector3D GetRefPos() const 
	{
		if (Dim()==3)
		{
			return Vector3D(XG(1),XG(2),XG(3));
		}
		else
		{
			return Vector3D(XG(1),XG(2),0.);
		}
	}
	virtual Vector3D GetRefPosInit() const 
	{
		if (Dim()==3)
		{
			return Vector3D(x_init(1),x_init(2),x_init(3));
		}
		else
		{
			return Vector3D(x_init(1),x_init(2),0.);
		}
	}
	virtual Vector3D GetRefVel() const 
	{
		if (Dim()==3)
		{
			return Vector3D(XGP(1),XGP(2),XGP(3));
		}
		else
		{
			return Vector3D(XGP(1),XGP(2),0.);
		}
	}
	//$ PG 2013-3-28: added GetRefConfPos to Element (for calling Constraint::GetBody3D(..).GetRefConfPos(...))
	virtual Vector3D GetRefConfPos(const Vector3D& p_loc) const
	{
		assert(0);
		return Vector3D(0);  //maby no assert and return GetRefPosInit()+ploc; ???
	}
	virtual Vector3D GetPos(const Vector3D& ploc) const
	{
		return GetRefPos()+ploc;
	}
	virtual Vector3D GetVel(const Vector3D& ploc) const
	{
		return GetRefVel();
	}
	//$ PG 2013-6-25: added for availablility in Element::GetFieldVariableValue(..), since displacement is added via Element::GetAvailableFieldVariables(..)
	virtual Vector3D GetDisplacement(const Vector3D& ploc) const 
	{
		return ploc;
	}
	//$ YV: sensors should use field variables, and the function below should be eliminated
	virtual Vector3D GetAcceleration(const Vector3D& ploc) const // global acceleration vector used in MBSSensor
	{
		assert(0);
		return Vector3D(0.);
	}

	//draw-functions:
	virtual Vector3D GetPosD(const Vector3D& ploc) const 
	{
		return GetRefPosD()+ploc;
	}
	virtual Vector3D GetVelD(const Vector3D& ploc) const
	{
		return GetRefVelD();
	}
	//$ PG 2013-6-25: added for availablility in Element::GetFieldVariableValue(..), since displacement is added via Element::GetAvailableFieldVariables(..)
	virtual Vector3D GetDisplacementD(const Vector3D& ploc) const 
	{
		return ploc;
	}

	virtual Vector3D GetRefPosD() const 
	{
		if (Dim()==3)
		{
			return Vector3D(XGD(1),XGD(2),XGD(3));
		}
		else
		{
			return Vector3D(XGD(1),XGD(2),0.);
		}
	}
	virtual Vector3D GetRefVelD() const 
	{
		if (Dim()==3)
		{
			return Vector3D(XGPD(1),XGPD(2),XGPD(3));
		}
		else
		{
			return Vector3D(XGPD(1),XGPD(2),0.);
		}
	}
	virtual Vector3D GetDOFPosD(int idof) const //returns postion of i-th DOF
	{
		return GetRefPosD();
	}
	virtual Vector3D GetDOFDirD(int idof) const //returns direction of action of i-th DOF
	{
		return Vector3D(0.,0.,0.);
	}


	virtual int NNodes() const {return 0; }; //$EDC$[funcaccess,readonly,EDCvarname="number_of_nodes",EDCfolder="Info",tooltiptext="number of nodes in finite elements"]
	virtual const int& NodeNum(int i) const {assert(0); return altshape;}
	virtual int& NodeNum(int i) {assert(0); return altshape;}

	virtual const Node& GetNode(int i) const {return GetMBS()->GetNode(NodeNum(i));} //get local node number i
	virtual Node& GetNode(int i) {
		int j = NodeNum(i);
		return GetMBS()->GetNode(j);
	} //get local node number i
	virtual Vector3D GetNodeLocPos(int i) const {assert(0); return Vector3D(0.,0.,0.);} //local position of node

	virtual Vector3D GetNodePos(int i) const {return GetPos(GetNodeLocPos(i));}

	virtual Vector3D GetNodePosD(int i) const {return GetPosD(GetNodeLocPos(i));}

	virtual Vector3D GetNodeVel(int i) const {return GetVel(GetNodeLocPos(i));}

	virtual Vector3D GetNodeVelD(int i) const {return GetVelD(GetNodeLocPos(i));}

	virtual Vector2D GetNodeLocPos2D(int i) const {assert(0); return Vector2D(0.,0.);} //local position of node

	virtual Vector2D GetNodePos2D(int i) const {return GetPos2D(GetNodeLocPos2D(i));}

	virtual Vector2D GetNodePos2DD(int i) const {return GetPos2DD(GetNodeLocPos2D(i));}

	virtual Vector2D GetNodeVel2D(int i) const {return GetVel2D(GetNodeLocPos2D(i));}

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//2D:
	virtual Vector2D GetRefPosInit2D() const 
	{
		return Vector2D(x_init(1),x_init(2));
	}
	virtual Vector2D GetRefPos2D() const 
	{
		return Vector2D(XG(1),XG(2));
	}
	virtual Vector2D GetRefVel2D() const 
	{
		return Vector2D(XGP(1),XGP(2));
	}

	virtual Vector2D GetPos2D(const Vector2D& ploc) const 
	{
		return GetRefPos2D()+ploc;
	}
	virtual Vector2D GetVel2D(const Vector2D& ploc) const
	{
		return GetRefVel2D();
	}
	//$ PG 2013-6-25: added for availablility in Element::GetFieldVariableValue(..), since displacement is added via Element::GetAvailableFieldVariables(..)
	virtual Vector2D GetDisplacement2D(const Vector2D& ploc) const 
	{
		return ploc;
	}

	//draw-functions 2D:
	virtual Vector2D GetRefPos2DD() const 
	{
		return Vector2D(XGD(1),XGD(2));
	}
	virtual Vector2D GetRefVel2DD() const 
	{
		return Vector2D(XGPD(1),XGPD(2));
	}

	virtual Vector2D GetPos2DD(const Vector2D& ploc) const 
	{
		return GetRefPos2DD()+ploc;
	}
	virtual Vector2D GetVel2DD(const Vector2D& ploc) const
	{
		return GetRefVel2DD();
	}
	//$ PG 2013-6-25: added for availablility in Element::GetFieldVariableValue(..), since displacement is added via Element::GetAvailableFieldVariables(..)
	virtual Vector2D GetDisplacement2DD(const Vector2D& ploc) const 
	{
		return ploc;
	}


	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	//Add load to Element
	//$ DR 2012-10:[ loads moved from element to mbs
	// old code:
	//virtual void AddLoad(const MBSLoad& li);
	//virtual int NLoads() const {return loads.Length();} //EDC$[funcaccess,readonly,EDCvarname="number_of_loads",EDCfolder="Info"]
	//virtual const MBSLoad& GetLoad(int i) const {return *loads(i);}
	//virtual MBSLoad& GetLoad(int i) {return *loads(i);}
	//virtual void DeleteLoad(int i)
	//{
	//	delete loads(i);
	//	for (int j=i; j < NLoads(); j++)
	//	{
	//		loads(j) = loads(j+1);
	//	}
	//	loads.SetLen(loads.Length()-1);
	//}

	// new code:
	virtual void AddLoad(const MBSLoad& li);
	virtual int NLoads() const {return loads.Length();} //$EDC$[funcaccess,readonly,EDCvarname="number_of_loads",EDCfolder="Info"]
	virtual const MBSLoad& GetLoad(int i) const {return mbs->GetLoad(loads(i));}
	virtual MBSLoad& GetLoad(int i) {return mbs->GetLoad(loads(i));}
	virtual void DeleteLoad(int i);
	virtual int GetLoadNr(int i) {return loads(i);}
	virtual TArray<int> GetLoadNrs() {return loads;} 
	virtual void SetLoadNrs(TArray<int> load_nrs) {loads = load_nrs;} 

	//$ DR 2012-10:] loads moved from element to mbs

	

	//Add constraint, do not copy!!!
	virtual void AddConstraint(Constraint* c, int lind)
	{
		constraints.Add(c);
		constraintindices.Add(lind);

		int found = 0;
		for (int i=1; i <= constraints_nodouble.Length(); i++)
		{
			if (constraints_nodouble(i) == c) {found = 1; break;}
		}
		if (!found) constraints_nodouble.Add(c);
	}
	virtual void RemoveConstraints()
	{
		constraints.SetLen(0);
		constraints_nodouble.SetLen(0);
		constraintindices.SetLen(0);
	}
	virtual int AddSensor(int sensornum) {return sensors.Add(sensornum);}
	virtual int NSensors() const {return sensors.Length();}

	virtual int GetSensorNum(int i) const {return sensors(i);}
	virtual void SetSensorNum(int i, int sensnum) {sensors(i) = sensnum;}
	virtual const Sensor& GetSensor(int i) const {return mbs->GetSensor(sensors(i));}
	virtual Sensor & GetSensor(int i) {return mbs->GetSensor(sensors(i));}
	virtual double GetSpecialSensorValue(int nr, double time) const {return 0;}	// DR 2012-01-12: has to be overwritten in derived class

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//access functions for constraints
	virtual void AddElement(int en);
	virtual int NE() const {return elements.Length();}

	virtual int GetElnum(int i) const {return elements(i);}
	virtual void SetElnum(int i, int elnum) {elements(i) = elnum;}

	virtual const Element& GetElem(int i) const {return mbs->GetElement(elements(i));}
	virtual Element& GetElem(int i) {return mbs->GetElement(elements(i));}

	virtual void GetDirectFeedThroughElements(TArray<int>& elnums) const {} //add all elements which are depending on this element by direct feed through to list

	// functions to catch the input - implemented in derived class
	virtual int RespondToKey(int key) {return false;}

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	virtual int NC() const {return constraints.Length();}
	virtual Constraint* GetConstraint(int i) const {return constraints(i);}
	virtual void ExternalSetRefCoord(double val) {} //assigns an external input value to an element data --> used in rotactor and in LinearTransformation element

	virtual const Vector3D& GetSize() const {assert(0); return col;}

	virtual const int& GetMaterialNum() const {return materialnum;}
	virtual int& GetMaterialNum() {return materialnum;}
	virtual void SetMaterialNum(int mnum);
	virtual void AddMaterial(const Material& m);
	virtual const Material& GetMaterial() const {return GetMBS()->GetMaterial(materialnum);}
	virtual Material& GetMaterial() {return GetMBS()->GetMaterial(materialnum);}

	//virtual double GetRho() const {return rho;} //old version, should vanish in new or revised elements!
	virtual const double& Rho() const; //link to material data
	virtual const double& Em() const;  //link to material data
	virtual const double& Nu() const;  //link to material data
	virtual double& Rho(); //link to material data
	virtual double& Em();  //link to material data
	virtual double& Nu();  //link to material data

	virtual double GetMass() const {return mass;}
	virtual double GetVolume() const {return mass/Rho();}
	virtual double GetMassDamping() const {return damping_m;}
	virtual void SetMassDamping(double d) {damping_m = d;}

	//for loads, rigid bodies:
	virtual void GetDuxDq(Vector& dudq) {assert(0);};
	virtual void GetDuyDq(Vector& dudq) {assert(0);};
	virtual void GetDuzDq(Vector& dudq) {assert(0);};
	virtual void GetDrotxDq(Vector& dudq) {assert(0);};
	virtual void GetDrotyDq(Vector& dudq) {assert(0);};
	virtual void GetDrotzDq(Vector& dudq) {assert(0);};
	//volume loads for flexible bodies:
	virtual void GetIntDuDq(Matrix& dudq) {assert(0);}; //in fact it is DuDq Transposed
	// gravity load for flexible bodies
	virtual void GetIntRhoDuDq(Matrix& rhodudq)  //in fact it is Int_V Rho DuDq Transposed
	{
		GetIntDuDq(rhodudq);
		rhodudq *= Rho();
	}
	virtual void GetIntDkappaDq2D(Vector& dudq) {assert(0);};
	// centrifugal load for flexible bodies
	virtual void GetIntDuDqFCentrifugal(Matrix& dudq, const Vector3D& omega, const Vector3D& r0) {assert(0);};

	virtual void SetAltShape(int truefalse) {altshape = truefalse;}
	virtual int GetAltShape() const {return altshape;}

	virtual int NGeomElements() {return drawelements.Length();}
	virtual int AddGeomElement(int i) {return drawelements.Add(i);}
	virtual const int& GetGeomElementNum(int i) const {return drawelements(i);}
	virtual int& GetGeomElementNum(int i) {return drawelements(i);}
	virtual void DeleteGeomElement(int i) {drawelements.Erase(i);}
	virtual GeomElement* GetGeomElement(int i) {return GetMBS()->GetDrawElement(drawelements(i));}
	virtual int Add(const GeomElement& de) //adds GeomElement to MBS and adds index to drawelments list
	{
		int i = GetMBS()->Add(de, GetOwnNum());
		return drawelements.Add(i);

		//GeomElement* den = de.GetCopy(); //old version!
		//drawelements.Add(den);
	}

	virtual Box3D GetBoundingBox() const;
	virtual Box3D GetBoundingBoxD() const;
	virtual Box2D GetBoundingBox2D() const;
	virtual Box2D GetBoundingBox2DD() const;

	virtual Box3D GetElementBox() const
	{
		return Box3D(GetRefPos(),GetRefPos());
	}

	virtual Box3D GetElementBoxD() const
	{
		return Box3D(GetRefPosD(),GetRefPosD());
	}

	virtual Box2D GetElementBox2D() const
	{
		return Box2D((const Vector2D&) GetRefPos(), (const Vector2D&) GetRefPos());
	}

	virtual Box2D GetElementBox2DD() const
	{
		return Box2D((const Vector2D&) GetRefPosD(), (const Vector2D&) GetRefPosD());
	}

	virtual void DrawElementAdd();

	// the derived classes may fill in this list to inform the system
	// which field variables they can produce for contour plotting and for sensors
	virtual void GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables);
	// computes a value of a given field variable at a given point;
	// flagD is true for drawing (contour plotting) and false for sensing (computation),
	// and the actual computation should be based either on XG() or XGD()
	// three versions are available with different definitions of the point, in which the field variable is computed:
	// 3D local coordinates, 2D local coordinates, and local node number

	//$ PG 2013-6-25:[ return position by GetPos(), displacement by GetDisplacemenet(), and velocity by GetVel(), those functions have to be implemented in the inherited element classes, though.
	//virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & local_position, bool flagD) { assert(0); return 0; }
	//virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector2D & local_position, bool flagD) { assert(0); return 0; }
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & local_position, bool flagD);
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector2D & local_position, bool flagD);
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const double  & local_position, bool flagD) { assert(0); return 0; }; //$ DR 2013-06: added
	//$ PG 2013-6-25:]
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, int local_node_number, bool flagD) { assert(0); return 0; }			//$ YV 2012-06: added
	

	virtual void DrawElementPreProc() {};

	virtual void DrawElement();
	virtual void DrawElement2D();
	virtual void DrawElementLocalFrame();
	virtual void DrawElementVelocityVector();
	virtual const int& DrawElementFlag() const {return draw_element;}
	virtual int& DrawElementFlag() {return draw_element;}

	virtual const Vector3D& GetCol() const {return col;} 
	virtual Vector3D& GetCol() {return col;} 

	virtual Vector3D GetAngularMomentum() const { return GetAngularMomentum(Vector3D(0.,0.,0.)); }
	virtual Vector3D GetAngularMomentum(const Vector3D& p_ref) const { assert(0); return Vector3D(0.,0.,0.); }

	//$ YV 2013-01-03: for constraint elements this function adds the global numbers of constrained
	// dofs to the list (i.e. the old contents of the list remains in there!);
	// the function is used by the eigenmodes solver
	virtual void GetGlobalConstrainedDOFs(TArray<int>& dofs) const {}

//$ AD 2013-01-20: added functions to write name of ouput (load) to list available in InputOutputElement
	virtual int GetNOutputs() const {return 0;}
	virtual void SetOutputName(int output_nr, mystr& str) {}

	//$ AD 2013-07-08: added functions to maniputlate the connection lines between the IOElements 
	virtual void InsertConNode(int list_idx, int input_nr, Vector2D& node) {}
	virtual void DeleteConNode(int list_idx) {}
//$AD 2013-07-10: Change Position of a single element (MBS element, conNode, ...) 
	virtual void MoveConNode2D(int list_idx, double delta_x, double delta_y) {}
	virtual void MoveElement(double delta_x, double delta_y, double delta_z) {}

//$ AD 2013-02-22: added functions to manipulate an internal variable of InputOutputElement - increase or decrease the current value of IOGUIResponse
	virtual void Increase() {};
	virtual void Decrease() {};

protected:
	mystr elementname; //$EDC$[varaccess,EDCvarname="name",EDCfolder="",tooltiptext="name of the element"]
	mutable MBS* mbs; //$!JG 2011-02: switched to mutable in order to modify e.g. in GetConstraintDrift
	Vector x_init;   //initial conditions
	Vector data_init;   //initial conditions for data variables
	//Matrix rot;			//temporal matrix, 2D

	TArray<int> ltg; //local to global dof reference vector;
	TArray<int> ltgdata; //local to global data reference vector;

	TArray<int> constraintindices; //index of this element in the constraint (mostly 1 or 2)
	TArray<char> dependencies; //contains dependencies to coordinates
	int elnum; //$EDC$[varaccess,EDCvarname="element_number",EDCfolder="",readonly,tooltiptext="number of the element in the mbs"] //number in MBS-system
	TMBSElement type;


	//TArray<MBSLoad*> loads; //$ DR 2012-10: loads moved from element to mbs
	TArray<int> loads;  //$EDC$[varaccess,EDCvarname="loads",variable_length_vector,condition=!IsType(TConstraint),tooltiptext="Set loads attached to this element: 'nr_load1, nr_load2, ...' or empty"]
	TArray<Constraint*> constraints; //this constraint list might be double!
	TArray<Constraint*> constraints_nodouble; //one constraint is only added once to an element

	TArray<int> sensors; //$EDC$[varaccess,EDCvarname="sensors",EDCfolder="Info",variable_length_vector,tooltiptext="attached sensors",readonly]//if the element equations of motion are dependent on sensor values (e.g. as an input), then add the sensor numbers here
											 //this improves the Jacobian; otherwise, convergence might fail
	TArray<int> elements; //dependent elements (mainly for contraints, but also for follower loads!)

	//density for gravitational force ==> erase!!!:
	//double rho;   	 //JG+DR: should not be accessible! EDC$[varaccess,EDCvarname="density",EDCfolder="Physics",condition=GetMaterialNum()==0,tooltiptext="density of simple bodies, which do not obtain density from material"]

	int materialnum; //JG+DR: should not be accessible! EDC$[varaccess,EDCvarname="material_number",EDCfolder="Physics",tooltiptext="material number which contains the main material properties of the body or element"] //material number in MBS which contains material information

	double mass;			//JG+DR: should not be accessible! EDC$[varaccess,EDCvarname="mass",EDCfolder="Physics",tooltiptext="for rigid bodies and FFRF bodies, not available for constraints"]	//mass, for rigid bodies and FFRF bodies	// DR 2012-07
	double damping_m; //JG+DR: should not be accessible! EDC[varaccess,EDCvarname="mass_prop_damping",EDCfolder="Physics",tooltiptext="only for bodies, not available for constraints"] //damping constant for mass proportional damping

	Vector3D col; //$EDC$[varaccess,EDCvarname="RGB_color",EDCfolder="Graphics",tooltiptext="[red, green, blue] color of element, range = 0..1, use default color:[-1,-1,-1]"] //color for body
	
	//geomelements are managed in MBS main class
	TArray<int> drawelements; //$EDC$[varaccess,EDCvarname="geom_elements",EDCfolder="Graphics",variable_length_vector,condition=!IsType(TConstraint),tooltiptext="Set Geometric elements to represent body 'geomelem1, geomelem2, ...' or empty"]
	int altshape; //$EDC$[varaccess,EDCvarname="use_alternative_shape",EDCfolder="Graphics",int_bool,condition=!IsType(TConstraint),tooltiptext="Graphical representation of element with geom-objects that are attached to the element"]

	int draw_element; //$EDC$[varaccess,EDCvarname="show_element",EDCfolder="Graphics",int_bool,tooltiptext="Flag to draw element"]
}; //$EDC$[endclass,Element]


typedef enum {TRWElementDataNoAccess=0, TRWElementDataRead=1, TRWElementDataWrite=2, TRWElementDataReadWrite=3} TReadWriteElementDataVariableAccess;

class ReadWriteElementDataVariableType
{
public:
	ReadWriteElementDataVariableType():variable_name(),comp1(0),comp2(0),value(0.),tooltiptext()/*, RWaccess(1)*/ {RWaccess = TRWElementDataRead;}; // default constructor
	ReadWriteElementDataVariableType(const mystr& full_string) { GetVariableNameAndComponents(full_string); };
	ReadWriteElementDataVariableType(const mystr& name, int c1, int c2, double val, const mystr& tooltip, TReadWriteElementDataVariableAccess RWaccessI = TRWElementDataRead) // used in function GetAvailableSpecialValuesAuto
	{
		RWaccess = RWaccessI;
	  variable_name = name;
		comp1 = c1;
		comp2 = c2;
		value = 0.;
		tooltiptext = tooltip; 
	}
	//possibly needed in future: copy-constructor, destructor
  // convert an full string that CAN contain [] into desired parts ( variable_name / component1 / component2 )
	// example  "MassMatrix[2,2]" is converted to: { variable_name = "MassMatrix", comp1 = 2, comp2 = 2 }
public:
	int GetVariableNameAndComponents(const mystr& full_string) //initialization from string (for access by sensor)
	{
		//initialize with defaults:
		value = 0.;
		RWaccess = TRWElementDataRead;
		//tooltiptext = ""; //not needed!

// --> MYSTR function
		variable_name = full_string;
		int n_comp = variable_name.Trim_And_Get_Indices(comp1,comp2);

		if(n_comp>=0)	{return 1;}		// parsed successfully
		else {return 0;}						// error
	}

	// YV 2012-11-5: comparison, e.g. for checking consistensy of a sensor
	// JG 2013-01-11: changed to function which checks range
	//int operator==(const ReadWriteElementDataVariableType & RWdata) const
	//{
	//	return variable_name == RWdata.variable_name && comp1 == RWdata.comp1 && comp2 == RWdata.comp2;
	//}
	//this function checks, if the RWdata provided fits into the name and range (comp1 for maximum range of comp1 and comp2 for maximum range of comp2)
	//and if RWaccess fits
	int IsAvailable(const ReadWriteElementDataVariableType & RWdata) const
	{
		if (variable_name == RWdata.variable_name)
		{
			if (comp1 != 0)
			{
				if (RWdata.comp1 </*=*/ 0 || RWdata.comp1 > comp1) return 0;
			}
			if (comp2 != 0)
			{
				if (RWdata.comp2 </*=*/ 0 || RWdata.comp2 > comp2) return 0;
			}
			
			if (((int)RWdata.RWaccess & (int)RWaccess) == 0) return 0; //does not work for case that (int)RWdata.RWaccess == 3, but should not occur

			return 1;
		}

		return 0;
	}

public:
	mystr variable_name; //name of variable for read/write access
	int comp1; //component for vector variables or matrix row
	int comp2; //second component for matrix variable (column)
	TReadWriteElementDataVariableAccess RWaccess; //binary flags: 0==no access, 1==read access, 2==write access, 3==read and write access
	double value; //read/write value
	mystr tooltiptext; //short description matching the variable name - to be filled  ONLY by function GetAvailableSpecialValues
	double time; // time for evaluation
};

#endif