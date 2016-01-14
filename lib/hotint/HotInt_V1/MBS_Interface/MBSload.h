//#**************************************************************
//#
//# filename:             MBSload.h
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
//#**************************************************************

#ifndef MBSLOAD__H
#define MBSLOAD__H

#include "mathfunc.h"
#include "stepsettings.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  Load   Load   Load   Load   Load   Load   Load   Load   Load   Load   
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//TGCload ... load applied to generalized coordinate
//Tgravity ... load must be multiplied with mass (val=gravity)
//Tspecload ... specific load, e.g. applied at certain point at beam, surface load at FE
//Tspecloadtype .... eventually needed ... specify surface load in FE, ...
//Trotatingforce ... force applied to a point in direction normal to axis
//                   two such forces can be used to simulate moment (Kraeftepaar)
//                   2D: only one axis point needed, axis is in direction Z, node_ax2 can be dummy
//the following list does not coincide with the type names defined in MBS::MBS() !!!!
typedef enum {TGCload=1, Tbodyload=2, Tsurfaceload=3,
				      Tgravity=4, Tforcevector2D=5, Tforcevector3D=6, TcontactZ=7, Tsurfacepressure=8,
							Tdistributedmoment=9, Tcentrifugal=10, Trotatingforce=11,Tmomentvector2D=12, Tmomentvector3D=13,
							Tparameterint= 101, Tparameterdouble=102, Tparametervector3d=103} TMBSLoad;
//$ DR 2012-10: to become independent from element: "Tpointload=5, Tpointmoment=6," changed to "Tforcevector2D=5, Tforcevector3D=6" and "Tmomentvector2D=12, Tmomentvector3D=13,"

//$ AD 2013-01: added enum for load function
typedef enum {TLoadValue = 0, TLoadMathFunc = 1, TLoadIOElement = 2 } TMBSLoadFunc;



class MBSLoad 
{
public:
	MBSLoad();

	MBSLoad(const MBSLoad& l): /*data(),*/ steps()
	{
		CopyFrom(l);
	};
	virtual void CopyFrom(const MBSLoad& l)
	{
		val = l.val;
		gc = l.gc;
		loadtype = l.loadtype;
		specpos = l.specpos;
		forcepos = l.forcepos;
		vecval = l.vecval;
		dudq = l.dudq;
		//el = l.el;
		loadfunc = l.loadfunc;
		steps = l.steps;
		loadname = l.loadname;
		loadnum = l.loadnum;

		//$ DR 2012-10: changed to MathFunction
		//for (int i=1; i <= l.data.Length(); i++)
		//	data.Add(l.data(i));
		mathfunc = l.mathfunc;

		IOelement = l.IOelement;
		IOelement_outputnum = l.IOelement_outputnum;

		node_f = l.node_f;
		node_1 = l.node_1;
		node_2 = l.node_2;

		p_int = l.p_int;
		p_double = l.p_double;
		p_vector3d = l.p_vector3d;

		mbs = l.mbs; //$ DR 2012-10
		dim = l.dim; //$ DR 2012-10

		localBodyAxis = l.localBodyAxis; // SM 2013-01-11
	};

	virtual MBSLoad* GetCopy()
	{
		MBSLoad* l = new MBSLoad(*this);
		return l;
	}

	virtual	MBSLoad* GetCopy() const
	{
		MBSLoad* l = new MBSLoad(*this);
		return l;
	}

	virtual TMBSLoad LoadType() const {return loadtype;};
	virtual int Dim() const{return dim;};
	virtual void SetDim(int element_dim){dim = element_dim;}; //$ DR 2012-10

	virtual mystr LoadName() const;
	virtual void SetLoadName(const char* name) {loadname = name;}

	virtual int GetOwnNum() const {return loadnum;}
	virtual void SetOwnNum(int i) {loadnum = i;}

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer
	//virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	//virtual int SetElementDataAuto(const ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int CheckConsistency(Element& el, mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute

	//get force
	virtual void GetGCLoad(double& vali, int& gci);
	virtual void GetBodyLoad(double& vali, int& dir);
	virtual void GetGravity(double& vali, int& dir);
	virtual void GetForceVector3D(Vector3D & force, Vector3D& pos);
	virtual void GetMomentVector3D(Vector3D & moment, Vector3D& pos);
	virtual void GetForceVector2D(Vector2D & force, Vector2D& pos);
	virtual void GetMomentVector2D(double& moment, Vector2D& pos);
	virtual void GetCentrifugal(Vector3D & omega, Vector3D& pos);
	virtual void GetLoadFunc(TMBSLoadFunc& i) {loadfunc = (TMBSLoadFunc) i;}
	virtual void GetContactZ(double& zpos);
	virtual void GetSurfacePressure(double& pressure, int& surfacenum);	//surfacenum: normal: 0=X, 1=-X, 2=Y, 3=-Y, 4=Z, 5=-Z
	virtual void GetRotatingForce(double& vali, int& node_force, int& node_ax1, int& node_ax2);
	virtual void GetParameterInt(int* & p_inti);
	virtual void GetParameterDouble(double* & p_doublei);
	virtual void GetParameterVector3D(double* & p_vector3di);
	//virtual void GetHarmonic(double& omega_i, double& phase_i) ; 		//$ DR 2012-10: changed to MathFunction
	//virtual void GetTimeSpan(double& ts, double& te); 		//$ DR 2012-10: changed to MathFunction
	//virtual void GetTimeRamp(double& mag, double& t1on, double& t2on, double& t1off, double& t2off); 		//$ DR 2012-10: changed to MathFunction
	virtual void GetMathFunction(MathFunction& userfunction);
	
	virtual int const NLoadSteps() const;
	virtual double GetCStepsFact(double t) const;
	virtual int GetCStepNumber(double t) const;
	virtual double GetCStepTime(double t) const;
	virtual int GetRampMode(int i) const;
	virtual double GetFinalValue(int i) const;
	virtual int GetLoadType() const {return (int)loadtype;};
	virtual mystr GetLoadTypeString() const;
	virtual mystr GetLoadTypeTexDescription() const;


	//set force type and value:
	virtual void SetGCLoad(double vali, int gci);
	virtual void SetBodyLoad(double vali, int dir);
	virtual void SetGravity(double vali, int dir);
	virtual void SetForceVector3D(const Vector3D & force, const Vector3D& pos, const int useLocalBodyAxis = 0);
	virtual void SetForceVector3D(const Vector3D & force, const int& nodenr, const int useLocalBodyAxis = 0); // DR 2011-05-17
	virtual void SetMomentVector3D(const Vector3D & moment, const Vector3D& pos, const int useLocalBodyAxis = 0);
	virtual void SetForceVector2D(const Vector2D & force, const Vector2D& pos, const int useLocalBodyAxis = 0);
	virtual void SetMomentVector2D(double moment, const Vector2D& pos);
	virtual void SetSurfacePressure(double pressure, int surfacenum);	//surfacenum: normal: 0=X, 1=-X, 2=Y, 3=-Y, 4=Z, 5=-Z
	virtual void SetDistributedMoment2D(double moment); //for planar beam
  virtual void SetCentrifugalLoad(const Vector3D& omega, const Vector3D& pos);
	virtual void SetRotatingForce(double vali, int node_force, int node_ax1, int node_ax2);
	virtual void SetParameterInt(int vali, int* p_doublei);
	virtual void SetParameterDouble(double vali, double* p_doublei);
	virtual void SetParameterVector3D(Vector3D& vecvali, Vector3D* p_vertor3di);

	virtual void SetLoadFunc(int i) {loadfunc = (TMBSLoadFunc) i;}
	virtual void SetContactZ(double zpos);
	virtual void SetHarmonic(double omega_i, double phase_i);
	virtual void SetTimeSpan(double ts, double te);
	virtual void SetTimeRamp(double mag, double t1on, double t2on, double t1off=-1, double t2off=-1);
	virtual void SetMathFunction(MathFunction userfunction); // sets any mathfunction
	virtual void SetIOElement(int IOelementI, int IOelement_outputnumI);
	virtual int HasIOElement() const {return HasElementDependence(); }
	virtual void GetIOElementNr(int& IOelement_, int& IOelement_outputnum_)
	{
		 IOelement_ = IOelement;
		 IOelement_outputnum_ = IOelement_outputnum;
	}

	virtual void SetLoadSteps(const TArray<StepSettings>& settings);

	virtual double Evaluate(double t) const;
	virtual int GC() const {return gc;};
	//virtual void AddElementLoad(Vector& f, double t); //$ DR 2012-10: loads moved from element to mbs, old code
	virtual void AddElementLoad(Vector& f, double t, Element* el);

	//virtual void SetElement(Element* eli) {el = eli;}
	virtual void SetMBS(MBS* mbsi) {mbs = mbsi;}	//$ DR 2012-10
	virtual MBS* GetMBS() const {	return mbs;}	//$ DR 2012-10
	
	//virtual void DrawLoad(MBS* mbs); //$ DR 2012-10: loads moved from element to mbs, old code
	virtual void DrawLoad(MBS* mbs, Element* el);

	virtual int HasElementDependence() const //return element number to which the load is dependent, otherwise return 0
	{
		if (loadfunc == TLoadIOElement) return IOelement; 		//$ DR 2012-10: changed the concept of loadfunc (IOElement changed from 5 to 2)
		else return 0;
	}

private:
	mystr loadname;			//$ DR 2012-10
	int loadnum;				//$ DR 2012-12-05
	double val;
	int gc;
	TMBSLoad loadtype;  //gc-force, gravity, body load, surface load, point load, ...
	int specpos;        //specifies the direction (body load), surface load, ...
	Vector3D forcepos;  //specifies the position (point load)
	Vector3D vecval;		//Force vector for point load

	int localBodyAxis;	// flag which describes, if local or global coordinate system is used //$ MSax 2013-01

	Vector dudq;			  // Int_V \partial u_i/partial q dV
	Matrix H;
	//Element* el;

	TMBSLoadFunc loadfunc; //0=val, >0 : specific function //$ DR 2012-10: new concept 0 = val, 1 = MathFunction, 2 = IOElement

	TArray<StepSettings> steps; // array containing values of loadfunction for the intervals defined in TimeInt::CSEndtimes

	//TArray<double> data; //for time-dependent forces
	MathFunction mathfunc; // link a math function to load

	int IOelement;			//loadfunc = TLoadIOElement; Load Size is computed from the output of an IOElement
	int IOelement_outputnum; //local output number of IOelement

	int node_f;					//specifies the position (e.g. rotating force)
	int node_1;					
	int node_2;					

	int* p_int;
  double* p_double;
	Vector3D* p_vector3d;

	int dim; //$ DR 2012-10
	MBS* mbs;		//$ DR 2012-10
};

#endif
