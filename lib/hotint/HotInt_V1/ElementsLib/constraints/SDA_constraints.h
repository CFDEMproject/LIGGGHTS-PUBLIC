//#**************************************************************
//#
//# filename:             SDA_constraint.h
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
//#***************************************************************************************
 
#ifndef SDA_CONSTRAINTS__H
#define SDA_CONSTRAINTS__H

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//   SDA_Constraints   SDA_Constraints   SDA_Constraints   SDA_Constraints
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//SDActor
//SDRotActor
//SDActor2D
//SDRotActor2D
//HydraulicActor
//HydraulicActorDA
//InertialLinearSpringDamper
//InertialLinearSpringDamper2D
//LinearRotationalSpringDamper
//SpringDamperBearing

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Spring-Damper-Actor Element 3D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//Spring-Damper-Actuator connects two elements, with spring (k...stiffness, l0...initial length),
// damper (d ...  damping coefficient) and a constant force (fa ... constant actuator force)
class SDActor: public InputOutputElement
{
public:
	SDActor(MBS* mbsi):InputOutputElement(mbsi), loccoords(), dpdq(), ioElement_actionIndex()
	{
		InitForce();
		groundjoint = 0;
		forcedirected = 0;
		//--------------------------------
        //b: same as in InputOutputElement
		InitConstructor();
		//SetNStates(0);
		SetNOutputs(0);
        //e: same as in InputOutputElement
		//--------------------------------
	};

	SDActor(MBS* mbsi, int en1, int en2, const Vector3D& lc1, const Vector3D& lc2,double spring_stiffness, 
		double spring_init_len, double damping_coeff, double actuator_force, Vector3D cdim, const Vector3D& coli):InputOutputElement(mbsi), loccoords(), dpdq(), ioElement_actionIndex()
	{	
		SetSDActor(en1, en2, lc1, lc2, spring_stiffness, spring_init_len, damping_coeff, actuator_force, cdim, coli);
	};

	SDActor(MBS* mbsi, int en1, const Vector3D& lc1, const Vector3D& gc2, double spring_stiffness, double spring_init_len, 
		double damping_coeff, double actuator_force, Vector3D cdim, const Vector3D& coli):InputOutputElement(mbsi), loccoords(), dpdq(), ioElement_actionIndex()
	{	
		SetSDActor(en1, lc1, gc2, spring_stiffness, spring_init_len, damping_coeff, actuator_force, cdim, coli);
	};
  // set-functions
	virtual void SetSDActor(int en1, int en2, const Vector3D& lc1, const Vector3D& lc2,	double spring_stiffness,
		double spring_init_len, double damping_coeff, double actuator_force, Vector3D cdim, const Vector3D& coli);

	virtual void SetSDActor(int en1, const Vector3D& lc1, const Vector3D& gc2, double spring_stiffness, double spring_init_len, 
		double damping_coeff, double actuator_force,Vector3D cdim, const Vector3D& coli);

	virtual void SetActionIndex(int action_type){assert(action_type == TSRefPos || action_type == TSRefForce);ioElement_actionIndex.Add(action_type);} //type of input (valve position/force/...)
	// set functions without constructor (use constructor with mbs for initialization)
	virtual void SetSDActor(int en1, const Vector3D& lc1, const Vector3D& gc2, double spring_stiffness, double spring_init_len, Vector3D& direction,
		double damping_coeff, double actuator_force, Vector3D cdim, const Vector3D& coli);

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new SDActor(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		InputOutputElement::CopyFrom(e);
		const SDActor& ce = (const SDActor&)e;
		loccoords = ce.loccoords;
		//p_global = ce.p_global;
		dpdq = ce.dpdq;
		k = ce.k;
		d = ce.d;
		l0 = ce.l0;
		fa = ce.fa;
		p_global = ce.p_global;
		forcemode = ce.forcemode;
		groundjoint = ce.groundjoint;

		k2 = ce.k2;
		k3 = ce.k3;
		x_range1 = ce.x_range1;
		x_range2 = ce.x_range2;

		mathfunc_k_flag = ce.mathfunc_k_flag;
		mathfunc_k = ce.mathfunc_k;
		mathfunc_d = ce.mathfunc_d;
		elementname = GetElementSpec();

		forcedirected = ce.forcedirected;
		locaxis1 = ce.locaxis1;
		ioElement_actionIndex = ce.ioElement_actionIndex;
	}

	//-----------------------------
    //b: same as in InputOutputElement
	virtual void InitConstructor()
	{
		InputOutputElement::InitConstructor();
		elementname = SymbolText();
	}
	//virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	//virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	//virtual int SetElementData(const ElementDataContainer& edc); //set element data according to ElementDataContainer

	virtual double GetOutput(double t, int i=1) const;
	
	virtual int IsDirectFeedThrough() const 
	{
		return 0;
	}


	//e: same as in InputOutputElement
	//-----------------------------


	virtual void InitForce()
	{
		forcemode = 0;

		k2 = 0;
		k3 = 0;
		x_range1 = 0;
		x_range2 = 0;
		mathfunc_d.SetConstant(0);
		mathfunc_k_flag = 0;
		mathfunc_k.SetConstant(0);
	}

	virtual const char* GetElementSpec() const {return "SpringDamperActor";}
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

	virtual int GetForceMode() const {return forcemode;}
	virtual void SetForceMode(int m) {forcemode = m;}

	virtual int IsGroundJoint() const {return groundjoint;}

	virtual void SetForceDirected(const Vector3D& locaxis1I)
	{
		locaxis1 = locaxis1I;
		locaxis1.Normalize();
		forcedirected = 1;
	}

	virtual void SetPolynomial(double k2I, double k3I) //for sgn(x)*Sqr(x)*k2, etc.
	{
		if (forcemode == 0) forcemode = 2;

		k2 = k2I;
		k3 = k3I;
	}
	virtual void SetRange(double range1, double range2) 
	{
		forcemode = 3;

		x_range1 = range1;
		x_range2 = range2;
	}
	virtual void SetPiecewiseLinearSpring(const Vector& xpos, const Vector& force)
	{
		forcemode = 2;
		mathfunc_k.SetPiecewise(xpos, force, 1);
		mathfunc_k_flag = 0;
	}
	virtual void SetPiecewiseLinearSpringStiffness(const Vector& xpos, const Vector& stiffness)
	{
		forcemode = 2;
		mathfunc_k.SetPiecewise(xpos, stiffness, 1);
		mathfunc_k_flag = 1;
	}
	virtual void SetPiecewiseLinearDamper(const Vector& vel, const Vector& force) 
	{
		forcemode = 2;
		mathfunc_d.SetPiecewise(vel, force, 1);
	}
	virtual void SetController(int elem_ind)
	{
		if (elements.Length() == 3-IsGroundJoint())
			elements(3-IsGroundJoint()) = elem_ind; 
		else
			elements.Add(elem_ind); 
		forcemode = 4;
	}

	virtual Vector3D ComputeForce(double t) const;    

	virtual Vector3D ComputeForceDirection() const; 
	virtual Vector3D ComputeForceDirectionD() const; 
	
	virtual double GetActorForce(double computation_time, int dir=0) const;
	virtual int GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //for special value sensor
	virtual int ReadSingleElementData(ReadWriteElementDataVariableType& RWdata); 		//for special value sensor

  //old: virtual double GetActorForce() {return ComputeForce()*ComputeForceDirection();} // scalar product
	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size
	virtual int IS() const {return 0;};

	virtual int Dim() const {return 3;}

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;

	virtual void DrawElement();

	virtual Vector2D GetInputPosD(int i) const; //return absolute position of input #i
	virtual Vector3D GetInputPos3DD(int i) const; //return absolute position of input #i (3D)

	virtual Vector2D GetOutputPosD(int i) const; //return absolute position of input #i	
	virtual Vector3D GetOutputPos3DD(int i) const;

protected:
	TArray<Vector3D> loccoords; //locp1, locp2
	Vector3D p_global;
	Vector3D locaxis1;

	Matrix dpdq; //temporary element, no need to copy???

	double k,l0,d,fa;
	int groundjoint;	//1==ground
	int forcedirected; //0==direction given by attached points; 1==direction of spring given by local vector of body 1

	int forcemode; //0=polynomial, 1=tension force only+polynomial, 2=polynomial, 3=with zero zone+polynomial,
	//							 4=controller for Ma with control elements!!!, 5=piecewise linear by mathfunc

	MathFunction mathfunc_k; //evaluate(x) gives a force in case mathfunc_k_flag == 0
	MathFunction mathfunc_d; //evaluate(x) gives a force
	int mathfunc_k_flag; // 0...mathfunc_k is force F(l-l0); 1...mathfunc_k is stiffness depending on position 'l': F(l,l0)=c12(l)*(l-l0)

	double x_range1, x_range2;
	double k2, k3; //coefficients for square and cubic dependence on angle

	TArray<int> ioElement_actionIndex; // number for definition type of input signal, if InputOutputElement is used in ComputeForce  (see TConstraintIOType)
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Rotational Spring-Damper-Actor Element 3D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//Rotational Spring-Damper-Actuator connects two elements, with spring (k...stiffness, phi0...initial rotation),
// damper (d ...  damping coefficient) and a constant moment (ma ... constant actuator force)
//positive rotation around rotation axis according to right hand rule!!
//rot x n1 = t1, rotation is positive if rot*(n1 x n2) >= 0
class SDRotActor: public InputOutputElement
{
public:
	SDRotActor(MBS* mbsi):InputOutputElement(mbsi), loccoords(), dpdq(), ioElement_actionIndex()
	{
		groundjoint = 0;
		InitForce();
	};

	SDRotActor(MBS* mbsi, int en1, int en2, 
		const Vector3D& lc1, const Vector3D& lc2, const Vector3D& globrot, 
		double spring_stiffness, double spring_init_angle, 
		double damping_coeff, double actuator_moment,
		Vector3D cdim, const Vector3D& coli):InputOutputElement(mbsi), loccoords(), dpdq(), ioElement_actionIndex() 
	{	
		groundjoint = 0;
		x_init = Vector(0);
		GetCol() = coli;
		draw_dim = cdim;
		AddElementCoord(en1, lc1);
		AddElementCoord(en2, lc2);
		loccoords.Add(globrot);
		loccoords.Add(0); //local rotation axis body1
		loccoords.Add(0); //local rotation axis body2
		loccoords.Add(0); //local rotation normal body1
		loccoords.Add(0); //local rotation normal body2
		k = spring_stiffness;
		phi0 = spring_init_angle;
		d = damping_coeff;
		ma = actuator_moment;

		InitForce();
		elementname = GetElementSpec();
	};

	SDRotActor(MBS* mbsi, int en1, 
		const Vector3D& lc1, const Vector3D& gc2, const Vector3D& globrot, 
		double spring_stiffness, double spring_init_angle, 
		double damping_coeff, double actuator_moment,
		Vector3D cdim, const Vector3D& coli):InputOutputElement(mbsi), loccoords(), dpdq(), ioElement_actionIndex()
	{	
		groundjoint = 1;
		x_init = Vector(0);
		GetCol() = coli;
		draw_dim = cdim;
		p_global = gc2;
		AddElementCoord(en1, lc1);
		loccoords.Add(globrot);
		loccoords.Add(0); //local rotation axis
		loccoords.Add(0); //local rotation normal
		loccoords.Add(0); //global rotation normal
		k = spring_stiffness;
		phi0 = spring_init_angle;
		d = damping_coeff;
		ma = actuator_moment;

		InitForce();
		elementname = GetElementSpec();
	};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new SDRotActor(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const SDRotActor& ce = (const SDRotActor&)e;
		loccoords = ce.loccoords;
		//p_global = ce.p_global;
		dpdq = ce.dpdq;
		k = ce.k;
		d = ce.d;
		phi0 = ce.phi0;
		ma = ce.ma;
		p_global = ce.p_global;
		forcemode = ce.forcemode;
		groundjoint = ce.groundjoint;

		k2 = ce.k2;
		k3 = ce.k3;
		phi_range1 = ce.phi_range1;
		phi_range2 = ce.phi_range2;

		mathfunc_k = ce.mathfunc_k;
		mathfunc_d = ce.mathfunc_d;
		
		ioElement_actionIndex = ce.ioElement_actionIndex;  // unique type number (e.g.: valve position for hydraulic, ...)
	}

	virtual void SetActionIndex(TConstraintIOType action_type)
	{
		assert(action_type == TSRefAngle || action_type == TSRefMom);
		ioElement_actionIndex.Add(action_type);
	}

	virtual const char* GetElementSpec() const {return "RotSpringDamperActor";}
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	virtual void AddElementCoord(int en, const Vector3D& loccoord)
	{
		AddElement(en);
		loccoords.Add(loccoord);
	}

	virtual void Initialize();
	virtual void InitForce()
	{
		forcemode = 0;

		k2 = 0;
		k3 = 0;
		phi_range1 = 0;
		phi_range2 = 0;
		mathfunc_d.SetConstant(0);
		mathfunc_k.SetConstant(0);
	}



	virtual void EvalG(Vector& f, double t) {};
	virtual void EvalF(Vector& f, double t) {};
	virtual void EvalF2(Vector& f, double t) {};

	virtual int IsGroundJoint() const {return groundjoint;}
	virtual int GetForceMode() const {return forcemode;}
	virtual void SetForceMode(int m) {forcemode = m;}

	virtual void SetPolynomial(double k2I, double k3I) //for sgn(phi)*Sqr(phi)*k2, etc.
	{
		if (forcemode == 0) forcemode = 2;

		k2 = k2I;
		k3 = k3I;
	}
	virtual void SetRange(double range1, double range2) 
	{
		forcemode = 3;

		phi_range1 = range1;
		phi_range2 = range2;
	}
	virtual void SetPiecewiseLinearSpring(const Vector& xpos, const Vector& force) 
	{
		forcemode = 2;
		mathfunc_k.SetPiecewise(xpos, force, 1);
	}
	virtual void SetPiecewiseLinearDamper(const Vector& vel, const Vector& force) 
	{
		forcemode = 2;
		mathfunc_d.SetPiecewise(vel, force, 1);
	}

	virtual double ComputeTorque(double t) const;  //?J const is possible now,  ma local defined 
	//virtual double ComputeTorque(); 

	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size
	virtual int IS() const {return 0;};

	virtual int Dim() const {return 3;}

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;

	virtual void DrawElement();

	virtual void SetController(int elem_ind)
	{
		if (elements.Length() == 3-IsGroundJoint())
			elements(3-IsGroundJoint()) = elem_ind; 
		else
			elements.Add(elem_ind); 
		forcemode = 4;
	}
	
	// overwrite initial spring lenght with external value
	//virtual void ExternalSetRefCoord(double phi0I) {phi0 = phi0I;}

	//// add input output element e.g.: for definition of spring length
	virtual double GetActorForce(double computation_time, int dir=0) const 
	{
		if (dir != 3) return ComputeTorque(computation_time);
		else return phi0;
	} 

	virtual double GetOutput(double t, int i=1) const;

protected:
	TArray<Vector3D> loccoords;
	Vector3D p_global;

	Matrix dpdq; //temporary element, no need to copy???

	int groundjoint; //==1 if groundjoint
	double k,phi0,d,ma;
	int forcemode; //0=linear,, 3=with zero zone, then polynomial,
	//							 4=controller for Ma with control elements!!!, 2=given by mathfunc
	

	MathFunction mathfunc_k, mathfunc_d; //Evaluate(x) gives already a force!!!!
	double phi_range1, phi_range2;
	double k2, k3; //coefficients for square and cubic dependence on angle
	TArray<int> ioElement_actionIndex;  // unique type number (e.g.: valve position for hydraulic, ...)
};



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Spring-Damper-Actor Element 2D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//Spring-Damper-Actuator connects two elements, with spring (k...stiffness, l0...initial length),
// damper (d ...  damping coefficient) and a constant force (fa ... constant actuator force)
class SDActor2D: public InputOutputElement
{
public:
	SDActor2D(MBS* mbsi):InputOutputElement(mbsi), loccoords(), dpdq()
	{
		forcemode = 0;
	};

	SDActor2D(MBS* mbsi, int en1, int en2, 
		const Vector2D& lc1, const Vector2D& lc2, 
		double spring_stiffness, double spring_init_len, 
		double damping_coeff, double actuator_force,
		Vector3D cdim, const Vector3D& coli):InputOutputElement(mbsi), loccoords(), dpdq() 
	{	
		x_init = Vector(0);
		GetCol() = coli;
		draw_dim = cdim;
		AddElementCoord(en1, lc1);
		AddElementCoord(en2, lc2);
		k = spring_stiffness;
    k_delay = -1;
		l0 = spring_init_len;
		d = damping_coeff;
		fa = actuator_force;
		forcemode = 0;
		elementname = GetElementSpec();
	};

	SDActor2D(MBS* mbsi, int en1, 
		const Vector2D& lc1, const Vector2D& gc2, 
		double spring_stiffness, double spring_init_len, 
		double damping_coeff, double actuator_force,
		Vector3D cdim, const Vector3D& coli):InputOutputElement(mbsi), loccoords(), dpdq() 
	{	
		x_init = Vector(0);
		GetCol() = coli;
		draw_dim = cdim;
		AddElementCoord(en1, lc1);
		k = spring_stiffness;
    k_delay = -1;
		l0 = spring_init_len;
		d = damping_coeff;
		fa = actuator_force;
		p_global = gc2;
		forcemode = 0;
		elementname = GetElementSpec();
	};

  SDActor2D(MBS* mbsi, int en1, 
		const Vector2D& lc1, const Vector2D& gc2, 
    double spring_stiffness, Vector stiffness_params, double spring_init_len, 
		double damping_coeff, double actuator_force,
		Vector3D cdim, const Vector3D& coli):InputOutputElement(mbsi), loccoords(), dpdq() 
	{	
		x_init = Vector(0);
		GetCol() = coli;
		draw_dim = cdim;
		AddElementCoord(en1, lc1);

		k = spring_stiffness;
    //int kparams_len = stiffness_params.GetLen();
    kparams = stiffness_params; //.SetLen(kparams_len);
    k_delay = kparams(1);
    
		l0 = spring_init_len;
		d = damping_coeff;
		fa = actuator_force;
		p_global = gc2;
		forcemode = 0;
		elementname = GetElementSpec();
	};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new SDActor2D(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const SDActor2D& ce = (const SDActor2D&)e;
		loccoords = ce.loccoords;
		//p_global = ce.p_global;
		dpdq = ce.dpdq;
		k = ce.k;
    k_delay = ce.k_delay;
		d = ce.d;
		l0 = ce.l0;
		fa = ce.fa;
		p_global = ce.p_global;
		forcemode = ce.forcemode;

    if (k_delay > -1) kparams = ce.kparams;
	}

	virtual void AddElementCoord(int en, const Vector2D& loccoord)
	{
		AddElement(en);
		loccoords.Add(loccoord);
	}

	virtual void EvalG(Vector& f, double t) {};
	virtual void EvalF(Vector& f, double t) {};
	virtual void EvalF2(Vector& f, double t) {};

	virtual int GetForceMode() const {return forcemode;}
	virtual void SetForceMode(int m) {forcemode = m;}

	virtual Vector2D ComputeForce(double t);

	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size
	virtual int IS() const {return 0;};

	virtual int Dim() const {return 2;}

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;

	virtual void DrawElement();

	virtual double EvaluateActorForce();

protected:
	TArray<Vector2D> loccoords;
	Vector2D p_global;

	Matrix dpdq; //temporary element, no need to copy???

	double k,l0,d,fa;
  double k_delay;
  Vector kparams;
	int forcemode; //0=linear, 1=tension force only
};

class HydraulicData
{
	//everything is public
public:
	HydraulicData()
	{
		E_oil=1.4e9;//# oil stiffness
		V_0=0.001;   //# oil remaining volume
		A_zyl=0.; //# cross-section of cylinder
		p_s=200;    //# system pressure (bar)
		Q_n=2e-3;   //# nominal flow
		p_n=35.;     //# nominal pressure (bar)
		p_T=1.;      //# outside pressure (bar)
		l_Zyl0=0.;
		l_0off=0.;

		R_P=5.;    //# control parameter P
		R_D=0.;    //# control parameter D
		s_end=0.;  //# control endposition (final extension of cylinder)

		//double acting:
		A_zyl2=0.; //# cross-section of cylinder
		initp1 = 0;
		initp2 = 0;
		V_02 = 0.0135;   //# oil remaining volume
	}

	double E_oil;//# oil stiffness
	double V_0;   //# oil remaining volume
	double A_zyl; //# cross-section of cylinder
	double A_zyl2; //# cross-section of cylinder
	double p_s;    //# system pressure (bar)
	double Q_n;   //# nominal flow
	double p_n;     //# nominal pressure (bar)
	double p_T;      //# outside pressure (bar)
	double R_P;      //# control parameter P
	double R_D;     //# control parameter D
	double s_end;  //# final extension of cylinder
	double l_Zyl0;	//# Length where cylinder is closed (V=V_0)
	double l_0off;	//# length of piston outside cylinder when closed (for drawing)

	double initp1;	//# initial pressure p1
	double initp2;	//# initial pressure p2
	double V_02;   //# oil remaining volume when closed, second chamber
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Hydraulic-Actor
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Hydraulik Actor (k...stiffness, l0...initial length),
// damper (d ...  damping coefficient) and a constant force (fa ... constant actuator force)
class HydraulicActor: public SDActor
{
public:
	HydraulicActor(MBS* mbsi):SDActor(mbsi)
	{
	};

	HydraulicActor(MBS* mbsi, int en1, int en2, 
		const Vector3D& lc1, const Vector3D& lc2, 
		const HydraulicData& hydraulic_data,
		Vector3D cdim, const Vector3D& coli):SDActor(mbsi) 
	{	
		x_init = Vector(ES());
		x_init(1) = 0;
		GetCol() = coli;
		draw_dim = cdim;
		AddElementCoord(en1, lc1);
		AddElementCoord(en2, lc2);
		hyddata = hydraulic_data;
		t_on = 0;
		elementname = GetElementSpec();
	};


	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new HydraulicActor(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const HydraulicActor& ce = (const HydraulicActor&)e;
		hyddata = ce.hyddata;
		t_on = ce.t_on;
	}

	virtual void EvalF(Vector& f, double t);

	virtual Vector3D ComputeForce(double t);

	//implicit (algebraic) size
	virtual int IS() const {return 0;};
	virtual int ES() const {return 1;};

	virtual int Dim() const {return 3;}
	virtual void SetTon(int t) {t_on = t;}

	virtual void DrawElement();

protected:
	HydraulicData hyddata;
	double t_on;
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Hydraulic-Actor - double acting
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Hydraulik Actor (k...stiffness, l0...initial length),
// damper (d ...  damping coefficient) and a constant force (fa ... constant actuator force)
class HydraulicActorDA: public SDActor
{
public:
	HydraulicActorDA(MBS* mbsi):SDActor(mbsi)
	{
		elementname = GetElementSpec();
	};

	void SetHydraulicActorDA(int en1, int en2, 
		const Vector3D& lc1, const Vector3D& lc2, 
		const HydraulicData& hydraulic_data,
		Vector3D cdim, const Vector3D& coli)
	{	
		hyddata = hydraulic_data;

		x_init = Vector(ES());
		x_init(1) = hyddata.initp1; //initial pressure chamber 1
		x_init(2) = hyddata.initp2; //initial pressure chamber 2
		GetCol() = coli;
		draw_dim = cdim;
		AddElementCoord(en1, lc1);
		AddElementCoord(en2, lc2);
		t_on = 0;
		A_ext = 0;
	};


	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new HydraulicActorDA(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const HydraulicActorDA& ce = (const HydraulicActorDA&)e;
		hyddata = ce.hyddata;
		t_on = ce.t_on;
	}

	virtual void EvalF(Vector& f, double t);

	//implicit (algebraic) size
	virtual int IS() const {return 0;};
	virtual int ES() const {return 2;};
	virtual Vector3D ComputeForce(double t) const;

	virtual int Dim() const {return 3;}
	virtual void SetTon(int t) {t_on = t;}

	virtual void DrawElement();

protected:
	HydraulicData hyddata;
	double t_on;
	double A_ext;     // signal from InputOutputElement 
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Hydraulic-Actor2D - double acting
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Hydraulik Actor (k...stiffness, l0...initial length),
// damper (d ...  damping coefficient) and a constant force (fa ... constant actuator force)
// reference trajectory (A_des/A_ext) is defined in an InputOutputElement
class HydraulicActorDA2D: public SDActor2D
{
public:
	HydraulicActorDA2D(MBS* mbsi):SDActor2D(mbsi)
	{
		elementname = GetElementSpec();
	};

	void SetHydraulicActorDA2D(int en1, int en2, 
		const Vector2D& lc1, const Vector2D& lc2, 
		const HydraulicData& hydraulic_data,
		Vector3D cdim, const Vector3D& coli)
	{	
		hyddata = hydraulic_data;

		x_init = Vector(ES());
		x_init(1) = hyddata.initp1; //initial pressure chamber 1
		x_init(2) = hyddata.initp2; //initial pressure chamber 2
		GetCol() = coli;
		draw_dim = cdim;
		AddElementCoord(en1, lc1);
		AddElementCoord(en2, lc2);
		t_on = 0;
		A_ext = 0;

	};


	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new HydraulicActorDA2D(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const HydraulicActorDA2D& ce = (const HydraulicActorDA2D&)e;
		hyddata = ce.hyddata;
		t_on = ce.t_on;
		A_ext = ce.A_ext;
	}

	virtual void EvalF(Vector& f, double t);

	virtual Vector2D ComputeForce(double t);

	//implicit (algebraic) size
	virtual int IS() const {return 0;};
	virtual int ES() const {return 2;}; 

	virtual void SetTon(int t) {t_on = t;}

	virtual void DrawElement();

	//virtual void ExternalSetRefCoord(double A_des) {A_ext = A_des;}

	////            add input output element with io_type "TSvalve"
	//virtual void AddInputOutputElement(int io_elnum, int io_type, int io_outnum = 1){Constraint::AddInputOutputElement(io_elnum, io_type, io_outnum);};
	
	//function is called when computation of time step is started ... read output of InputOutputElement and overwrite phi0 (one spring end follows reference trajectory)
	//virtual void StartTimeStep(){if( IOGetNElements() == 1 && IOGetOutputType(1) == (TConstraintIOType)TSvalve)ExternalSetRefCoord(IOGetOutputValue(1, IOGetOutputNum(1)));};

protected:
	HydraulicData hyddata;
	double t_on;

	double A_ext;     // signal from InputOutputElement 
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Hydraulic-Actor2D - double acting NEW
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Hydraulik Actor (k...stiffness, l0...initial length),
// damper (d ...  damping coefficient) and a constant force (fa ... constant actuator force)
// reference trajectory (A_des/A_ext) should be defined in an InputOutputElement
class HydraulicActorDA2Dn: public SDActor2D
{
public:
	HydraulicActorDA2Dn(MBS* mbsi):SDActor2D(mbsi)
	{
		elementname = GetElementSpec();
	};

	void SetHydraulicActorDA2Dn(int en1, int en2, 
		const Vector2D& lc1, const Vector2D& lc2, 
		const HydraulicData& hydraulic_data,
		Vector3D cdim, const Vector3D& coli)
	{	
		hyddata = hydraulic_data;

		x_init = Vector(ES());
		x_init(1) = hyddata.initp1; //initial pressure chamber 1
		x_init(2) = hyddata.initp2; //initial pressure chamber 2
		GetCol() = coli;
		draw_dim = cdim;
		AddElementCoord(en1, lc1);
		AddElementCoord(en2, lc2);
		t_on = 0;
		A_ext = 0;

	};


	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new HydraulicActorDA2Dn(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const HydraulicActorDA2Dn& ce = (const HydraulicActorDA2Dn&)e;
		hyddata = ce.hyddata;
		t_on = ce.t_on;
		A_ext = ce.A_ext;
	}

	virtual void EvalF(Vector& f, double t);

	virtual Vector2D ComputeForce(double t);

	//implicit (algebraic) size
	virtual int IS() const {return 0;};
	virtual int ES() const {return 2;}; 

	virtual void SetTon(int t) {t_on = t;}

	virtual void DrawElement();

	//virtual void ExternalSetRefCoord(double A_des) {A_ext = A_des;}

	////            add input output element with io_type "TSvalve"
	//virtual void AddInputOutputElement(int io_elnum, int io_type, int io_outnum = 1){Constraint::AddInputOutputElement(io_elnum, io_type, io_outnum);};
	
	//function is called when computation of time step is started ... read output of InputOutputElement and overwrite phi0 (one spring end follows reference trajectory)
	//virtual void StartTimeStep(){if( IOGetNElements() == 1 && IOGetOutputType(1) == (TConstraintIOType)TSvalve)ExternalSetRefCoord(IOGetOutputValue(1, IOGetOutputNum(1)));};

protected:
	HydraulicData hyddata;
	double t_on;

	double A_ext;     // signal from InputOutputElement 
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// InertialLinearSpringDamper 3D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//add spring/damper force due to position and velocity in direction of global (inertial) system
class InertialLinearSpringDamper: public InputOutputElement
{
public:
	InertialLinearSpringDamper(MBS* mbsi):InputOutputElement(mbsi), loccoords(), dpdq() 
	{

		//--------------------------------
        //b: same as in InputOutputElement
		InitConstructor();
		//SetNStates(0);
		SetNOutputs(0);
        //e: same as in InputOutputElement
		//--------------------------------
	};


	void SetInertialLinearSpringDamperParameters(Matrix3D& spring_stiffness, Matrix3D& damping_coeff, const Vector3D& drawdimi, const Vector3D& coli)
	{
    SetXInit();

		GetCol() = coli;
		draw_dim = drawdimi;
		k = spring_stiffness;
		d = damping_coeff;

		elementname = GetElementSpec();
	}
		
	void SetXInit()
	{
		//not necessary, just ignore SetGlobalInitConditions
		int soselem1 = GetElem(1).SOS();

		Vector elem1xinitpos = GetElem(1).GetXInit().SubVector(1,soselem1);
		Vector elem1xinitvel = GetElem(1).GetXInit().SubVector(soselem1+1,2*soselem1);
		Vector elem2xinitpos;
		Vector elem2xinitvel;
		
		if (NE() == 2)
		{
			int soselem2 = GetElem(2).SOS();

			elem2xinitpos = GetElem(2).GetXInit().SubVector(1,soselem2);
			elem2xinitvel = GetElem(2).GetXInit().SubVector(soselem2+1,2*soselem2);
		}
		x_init = elem1xinitpos;
		x_init = x_init.Append(elem2xinitpos);
		x_init = x_init.Append(elem1xinitvel);
		x_init = x_init.Append(elem2xinitvel);
	}

	void SetInertialLinearSpringDamper(int en1, const Vector3D& lc1, 
		const Vector3D& pglob, Matrix3D& spring_stiffness, Matrix3D& damping_coeff, const Vector3D& drawdimi, const Vector3D& coli, const Vector3D& spring_init_len = Vector3D(0.))
	{	
		p_global = pglob;
		AddElementCoord(en1, lc1);

		SetInertialLinearSpringDamperParameters(spring_stiffness,damping_coeff,drawdimi,coli);
		l0 = spring_init_len;
	};

	void SetInertialLinearSpringDamper(int en1, int en2, 
		const Vector3D& lc1, const Vector3D& lc2, Matrix3D& spring_stiffness, Matrix3D& damping_coeff, const Vector3D& drawdimi, const Vector3D& coli, const Vector3D& spring_init_len = Vector3D(0.)) 
	{	
		AddElementCoord(en1, lc1);
		AddElementCoord(en2, lc2);
		SetInertialLinearSpringDamperParameters(spring_stiffness,damping_coeff,drawdimi,coli);
		l0 = spring_init_len;
	};

	
	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new InertialLinearSpringDamper(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		InputOutputElement::CopyFrom(e);
		const InertialLinearSpringDamper& ce = (const InertialLinearSpringDamper&)e;
		loccoords = ce.loccoords;
		p_global = ce.p_global;
		dpdq = ce.dpdq;
		k = ce.k;
		d = ce.d;
		l0 = ce.l0;
	}

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	virtual const char* GetElementSpec() const {return "InertialLinearSpringDamper";} //this name names MUST coincide with the type names defined in MBS::MBS() !!!!

	virtual void AddElementCoord(int en, const Vector3D& loccoord)
	{
		AddElement(en);
		loccoords.Add(loccoord);
	}

	virtual void EvalG(Vector& f, double t){assert(0 && "InertialLinearSpringDamper: EvalG called - should not happen!");};

	virtual Vector3D ComputeForce(double t) const;    
	virtual Vector3D ComputeForceDirection() const{assert(0 && "InertialLinearSpringDamper: ComputeForceDirection called (niy)!");return Vector3D(0.);}; 
	virtual Vector3D ComputeForceDirectionD() const{assert(0 && "InertialLinearSpringDamper: ComputeForceDirectionD called (niy)!");return Vector3D(0.);}; 
	virtual double GetActorForce(double computation_time, int dir=0) const;
	
	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f); // function no longer needed, calculation now in EvalF2
// switch back to evaluation with AddElementCqTLambda:
// *remove functions: EvalF2, LinkToElements, SOS, SOS2
// *in SetInertialLinearSpringDamperParameters: use  x_init = Vector(SS()); for initialization

	virtual void EvalF2(Vector& f, double t);
	virtual void LinkToElements();

	//implicit (algebraic) size
	virtual int IS() const {return 0;};
	virtual int Dim() const {return 3;}
	virtual int SOS() const;
	virtual int SOSowned() const {return 0;} // not equal SOS because of penalty

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;

	virtual void DrawElement();

protected:
	TArray<Vector3D> loccoords;
	Vector3D p_global;
	Vector3D init_disp; // (AD) added for spring with initial displacement, use with autocompute
  Matrix3D k, d; // (intertial) spring stiffness, damping
	Matrix dpdq;   // temporary element, no need to copy???
	Vector3D l0;     // l0 ... initial spring length
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// InertialSpringDamperNL 3D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//add non-linear spring/damper force due to position and velocity in direction of global (inertial) system
class InertialSpringDamperNL: public InertialLinearSpringDamper
{
public:
	InertialSpringDamperNL(MBS* mbsi):InertialLinearSpringDamper(mbsi) 
	{
		Vector datainit(DataS());
		datainit.SetAll(0.);
		SetDataInit(datainit);
		usestiffmfunc = 0; //Use MatFunc instead of stiffness matrix
	};


	void SetInertialSpringDamperNLParameters(Matrix3D& spring_stiffness, Matrix3D& damping_coeff, const Vector3D& drawdimi, const Vector3D& coli)
	{
		InertialLinearSpringDamper::SetInertialLinearSpringDamperParameters(spring_stiffness, damping_coeff, drawdimi,coli);
	}
		
	void SetXInit()
	{
		InertialLinearSpringDamper::SetXInit();
	}

	void SetInertialSpringDamperNL(int en1, const Vector3D& lc1, 
		const Vector3D& pglob, Matrix3D& spring_stiffness, Matrix3D& damping_coeff, const Vector3D& drawdimi, const Vector3D& coli, const Vector3D& spring_init_len = Vector3D(0.))
	{	
		InertialLinearSpringDamper::SetInertialLinearSpringDamper(en1,lc1, pglob, spring_stiffness, damping_coeff, drawdimi, coli, spring_init_len);
	};

	void SetInertialSpringDamperNL(int en1, int en2, 
		const Vector3D& lc1, const Vector3D& lc2, Matrix3D& spring_stiffness, Matrix3D& damping_coeff, const Vector3D& drawdimi, const Vector3D& coli, const Vector3D& spring_init_len = Vector3D(0.)) 
	{	
		InertialLinearSpringDamper::SetInertialLinearSpringDamper(en1,en2, lc1, lc2, spring_stiffness,damping_coeff, drawdimi, coli, spring_init_len); 
	};

	// MS, RL
	// Use MatFunc instead of stiffness matrix
	virtual void SetStiffnessMathFunction(const MathFunction& stiffmfuncXI,const MathFunction& stiffmfuncYI,const MathFunction& stiffmfuncZI)
	{
		usestiffmfunc = 1;
		stiffmfuncX = stiffmfuncXI;
		stiffmfuncY = stiffmfuncYI;
		stiffmfuncZ = stiffmfuncZI;
	}

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new InertialSpringDamperNL(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		InertialLinearSpringDamper::CopyFrom(e);
		const InertialSpringDamperNL& ce = (const InertialSpringDamperNL&)e;
		usestiffmfunc = ce.usestiffmfunc;
		stiffmfuncX = ce.stiffmfuncX;
		stiffmfuncY = ce.stiffmfuncY;
		stiffmfuncZ = ce.stiffmfuncZ;
	}

	virtual Vector3D ComputeForce(double t) const;
	
	virtual const char* GetElementSpec() const {return "InertialSpringDamperNL";} //this name names MUST coincide with the type names defined in MBS::MBS() !!!!

	virtual void EvalG(Vector& f, double t){assert(0 && "InertialSpringDamperNL: EvalG called - should not happen!");};
  
	virtual Vector3D ComputeForceDirection() const{assert(0 && "InertialSpringDamperNL: ComputeForceDirection called (niy)!");return Vector3D(0.);}; 
	virtual Vector3D ComputeForceDirectionD() const{assert(0 && "InertialSpringDamperNL: ComputeForceDirectionD called (niy)!");return Vector3D(0.);}; 

	virtual double PostNewtonStep(double t);
	virtual int DataS() const {return 9;} //Data size for non-state variables (length of XData-vector)

protected:

	MathFunction stiffmfuncX, stiffmfuncY, stiffmfuncZ; //MatFunc for stiffness
	int usestiffmfunc; //Use MatFunc instead of stiffness matrix

};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// InertialLinearSpringDamper 2D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//add spring/damper force due to position and velocity in direction of global (inertial) system

//use constraint similar to finite element ==> for local differentiation (modified newton)!
#define useEvalF2

class InertialLinearSpringDamper2D: public Constraint
{
public:
	InertialLinearSpringDamper2D(MBS* mbsi):Constraint(mbsi), loccoords(), dpdq() 
	{
		init_disp = Vector2D(0.,0.);
		autocompute2ndcoord = 0;
	};

	void SetInertialLinearSpringDamper2DParameters(Matrix3D& spring_stiffness, Matrix3D& damping_coeff, const Vector3D& drawdimi, const Vector3D& coli)
	{
		SetXInit();

		GetCol() = coli;
		draw_dim = drawdimi;
		k = spring_stiffness;
		d = damping_coeff;

		elementname = GetElementSpec();
	}

	void SetXInit()
	{
		int soselem1 = GetElem(1).SOS();

		Vector elem1xinitpos = GetElem(1).GetXInit().SubVector(1,soselem1);
		Vector elem1xinitvel = GetElem(1).GetXInit().SubVector(soselem1+1,2*soselem1);
		Vector elem2xinitpos;
		Vector elem2xinitvel;
		
		if (NE() == 2)
		{
			int soselem2 = GetElem(2).SOS();

			elem2xinitpos = GetElem(2).GetXInit().SubVector(1,soselem2);
			elem2xinitvel = GetElem(2).GetXInit().SubVector(soselem2+1,2*soselem2);
		}
		x_init = elem1xinitpos;
		x_init = x_init.Append(elem2xinitpos);
		x_init = x_init.Append(elem1xinitvel);
		x_init = x_init.Append(elem2xinitvel);
	}

	void SetInertialLinearSpringDamper2D(int en1, const Vector2D& lc1, 
		const Vector2D& pglob, Matrix3D& spring_stiffness, Matrix3D& damping_coeff, const Vector3D& drawdimi, const Vector3D& coli)
	{	
		p_global = pglob;
		AddElementCoord(en1, lc1);
		SetInertialLinearSpringDamper2DParameters(spring_stiffness,damping_coeff,drawdimi,coli);
	
	};

	//Ground Joint with local coordinates
	void SetInertialLinearSpringDamper2D1Loc(int en1, const Vector2D& lc1, 
		Matrix3D& spring_stiffness, Matrix3D& damping_coeff, const Vector3D& drawdimi, const Vector3D& coli)
	{	
		p_global = Vector2D(0.,0.);
		AddElementCoord(en1, lc1);
		SetInertialLinearSpringDamper2DParameters(spring_stiffness,damping_coeff,drawdimi,coli);

		autocompute2ndcoord = 1; // calculate from lc1
	};

//Ground Joint with global coordinates
	void SetInertialLinearSpringDamper2D1Glob(int en1, const Vector2D& pglob, 
		Matrix3D& spring_stiffness, Matrix3D& damping_coeff, const Vector3D& drawdimi, const Vector3D& coli)
	{
		p_global = pglob;
		AddElementCoord(en1, Vector2D(0.0,0.0));
		autocompute2ndcoord = 2; // calculate from pglob
		SetInertialLinearSpringDamper2DParameters(spring_stiffness,damping_coeff,drawdimi,coli);
	}

//2-Body Joint with 2 local coordinates
	void SetInertialLinearSpringDamper2D(int en1, int en2, 
		const Vector2D& lc1, const Vector2D& lc2, Matrix3D& spring_stiffness, Matrix3D& damping_coeff, const Vector3D& drawdimi, const Vector3D& coli) 
	{	
		AddElementCoord(en1, lc1);
		AddElementCoord(en2, lc2);
		SetInertialLinearSpringDamper2DParameters(spring_stiffness,damping_coeff,drawdimi,coli);
	};

//2-Body Joint with 1 local coordinate
	void SetInertialLinearSpringDamper2D1Loc(int en1, int en2, 
		const Vector2D& lc1, Matrix3D& spring_stiffness, Matrix3D& damping_coeff, const Vector3D& drawdimi, const Vector3D& coli) 
	{	
		p_global = Vector2D(0.0,0.0);
		AddElementCoord(en1, lc1);
		AddElementCoord(en2, Vector2D(0.0));	
		autocompute2ndcoord = 1; // calculate from lc1
		SetInertialLinearSpringDamper2DParameters(spring_stiffness,damping_coeff,drawdimi,coli);
	};

//2-Body Joint with global coordinates
	void SetInertialLinearSpringDamper2D1Glob(int en1, int en2, 
		const Vector2D& pglob, Matrix3D& spring_stiffness, Matrix3D& damping_coeff, const Vector3D& drawdimi, const Vector3D& coli) 
	{	
		p_global = pglob;
		AddElementCoord(en1, Vector2D(0.0,0.0));
		AddElementCoord(en2, Vector2D(0.0,0.0));
		autocompute2ndcoord = 2; // calculate from pglob
		SetInertialLinearSpringDamper2DParameters(spring_stiffness,damping_coeff,drawdimi,coli);
	}

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new InertialLinearSpringDamper2D(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const InertialLinearSpringDamper2D& ce = (const InertialLinearSpringDamper2D&)e;
		loccoords = ce.loccoords;
		p_global = ce.p_global;
		init_disp = ce.init_disp;
		autocompute2ndcoord = ce.autocompute2ndcoord;
		dpdq = ce.dpdq;
		k = ce.k;
		d = ce.d;
	}

	virtual void Initialize();

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	virtual const char* GetElementSpec() const {return "InertialLinearSpringDamper2D";} //this name names MUST coincide with the type names defined in MBS::MBS() !!!!

	virtual void AddElementCoord(int en, const Vector2D& loccoord)
	{
		AddElement(en);
		loccoords.Add(loccoord);
	}

	virtual void EvalG(Vector& f, double t){assert(0 && "InertialLinearSpringDamper2D: EvalG called - should not happen!");};

	virtual Vector2D ComputeForce(double t) const;    
	virtual Vector2D ComputeForceDirection() const{assert(0 && "InertialLinearSpringDamper2D: ComputeForceDirection called (niy)!");return Vector2D(0.);}; 
	virtual Vector2D ComputeForceDirectionD() const{assert(0 && "InertialLinearSpringDamper2D: ComputeForceDirectionD called (niy)!");return Vector2D(0.);}; 
	virtual double GetActorForce(double computation_time, int dir=0) const;
	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f); // function no longer needed, calculation now in EvalF2
// switch back to evaluation with AddElementCqTLambda:
// *remove functions: EvalF2, LinkToElements, SOS, SOS2
// *in SetInertialLinearSpringDamper2DParameters: use  x_init = Vector(SS()); for initialization

	virtual void EvalF2(Vector& f, double t);
	virtual void LinkToElements();

	//implicit (algebraic) size
	virtual int IS() const {return 0;};
	virtual int Dim() const {return 2;}
	virtual int SOS() const;
	virtual int SOSowned() const {return 0;}

	//get a reference position of the body in 3d
	virtual Vector2D GetRefPos2DD() const;

	virtual void SetInitDisp(const Vector2D& x0) {init_disp = x0;}

	virtual void DrawElement();

	virtual void SetAutoComputeFromGlobal() { autocompute2ndcoord = 2; }
	virtual void SetAutoComputeFromLocal() { autocompute2ndcoord = 1; }
	virtual Vector2D& GetPGlobal() { return p_global; }

protected:
	TArray<Vector2D> loccoords;
	Vector2D p_global;
	Vector2D init_disp; // (AD) added for spring with initial displacement, use with autocompute
	int autocompute2ndcoord; // automatically compute all other coordinates from (=1): local coord 1, (=2): global coord
  Matrix3D k, d; // (intertial) spring stiffness, damping
	Matrix dpdq; // temporary element, no need to copy???
};



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// LinearRotationalSpringDamper 3D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//add spring/damper force due to relative rotation of two bodies
class LinearRotationalSpringDamper: public Constraint
{
public:
	LinearRotationalSpringDamper(MBS* mbsi):Constraint(mbsi), loccoords(), dpdq() 
	{
	};
	//void SetLinearRotationalSpringDamper(int en1, const Vector3D& lc1, 
	//	const Vector3D& pglob, Matrix3D& spring_stiffness, Matrix3D& damping_coeff, const Vector3D& drawdimi, const Vector3D& coli)
	//{	
	//	x_init = Vector(SS());
	//	col = coli;
	//	draw_dim = drawdimi;
	//	p_global = pglob;
	//	k = spring_stiffness;
	//	d = damping_coeff;
	//	AddElementCoord(en1, lc1);

	//	elementname = GetElementSpec();
	//};

	void SetLinearRotationalSpringDamper(int en1, int en2, 
		const Vector3D& lc1, const Vector3D& lc2, Matrix3D& spring_stiffness, Matrix3D& damping_coeff, const Vector3D& drawdimi, const Vector3D& coli) 
	{	
		x_init = Vector(SS());
		GetCol() = coli;
		draw_dim = drawdimi;
		AddElementCoord(en1, lc1);
		AddElementCoord(en2, lc2);
		k = spring_stiffness;
		d = damping_coeff;
		elementname = GetElementSpec();
	};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new LinearRotationalSpringDamper(*this);
		//ec.CopyFrom(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const LinearRotationalSpringDamper& ce = (const LinearRotationalSpringDamper&)e;
		loccoords = ce.loccoords;
		p_global = ce.p_global;
		dpdq = ce.dpdq;
		k = ce.k;
		d = ce.d;
	}

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	virtual const char* GetElementSpec() const {return "LinearRotationalSpringDamper";} //this name names MUST coincide with the type names defined in MBS::MBS() !!!!

	virtual void AddElementCoord(int en, const Vector3D& loccoord)
	{
		AddElement(en);
		loccoords.Add(loccoord);
	}

	virtual void EvalG(Vector& f, double t){assert(0 && "LinearRotationalSpringDamper: EvalG called - should not happen!");};

	virtual Vector3D ComputeForce(double t) const;    
	virtual Vector3D ComputeForceDirection() const{assert(0 && "LinearRotationalSpringDamper: ComputeForceDirection called (niy)!");return Vector3D(0.);}; 
	virtual Vector3D ComputeForceDirectionD() const{assert(0 && "LinearRotationalSpringDamper: ComputeForceDirectionD called (niy)!");return Vector3D(0.);}; 
	virtual double GetActorForce(double computation_time, int dir=0) const;
	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size
	virtual int IS() const {return 0;};

	virtual int Dim() const {return 3;}

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;

	virtual void DrawElement();

protected:
	TArray<Vector3D> loccoords;
	Vector3D p_global;
  Matrix3D k, d; // (intertial) spring stiffness, damping
	Matrix dpdq; // temporary element, no need to copy???
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// SpringDamperBearing 3D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//add spring/damper force due to position and velocity in direction of global (inertial) system
class SpringDamperBearing: public InertialLinearSpringDamper
{
public:
	SpringDamperBearing(MBS* mbsi):InertialLinearSpringDamper(mbsi)
	{
	};
	void SetSpringDamperBearing(int en1, const Vector3D& lc1, 
		const Vector3D& pglob, Matrix3D& spring_stiffness, Matrix3D& damping_coeff, const Vector3D& drawdimi, const Vector3D& coli)
	{	
		InertialLinearSpringDamper::SetInertialLinearSpringDamper(en1, lc1, pglob, spring_stiffness, damping_coeff, drawdimi, coli);
		elementname = GetElementSpec();
		k_nonlinear = 0;
		use_bearing_fault = 0;
		n_bearing_balls = 0;
		draw_bearing = 0;
		random_value = 0;
		random_noise_factor = 0;
	};

	void SetSpringDamperBearing(int en1, int en2, 
		const Vector3D& lc1, const Vector3D& lc2, Matrix3D& spring_stiffness, Matrix3D& damping_coeff, const Vector3D& drawdimi, const Vector3D& coli) 
	{	
		InertialLinearSpringDamper::SetInertialLinearSpringDamper(en1, en2, lc1, lc2, spring_stiffness, damping_coeff, drawdimi, coli);
		elementname = GetElementSpec();
		k_nonlinear = 0;
		use_bearing_fault = 0;
		n_bearing_balls = 0;
		draw_bearing = 0;
		random_value = 0;
		random_noise_factor = 0;
	};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new SpringDamperBearing(*this);

		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		InertialLinearSpringDamper::CopyFrom(e);
		const SpringDamperBearing& ce = (const SpringDamperBearing&)e;

		rotcomp = ce.rotcomp;
		mf_kxx = ce.mf_kxx;
		mf_kyy = ce.mf_kyy;
		k_nonlinear = ce.k_nonlinear;

		bearing_fault_freqfact = ce.bearing_fault_freqfact;
		bearing_fault_amp = ce.bearing_fault_amp;		
		bearing_fault_power = ce.bearing_fault_power;	
		use_bearing_fault = ce.use_bearing_fault;
		//bearing_fault_dir = ce.bearing_fault_dir;

		bearing_fault_innerring = ce.bearing_fault_innerring;

		draw_bearing = ce.draw_bearing;

		n_bearing_balls = ce.n_bearing_balls;  
		inner_roll_radius = ce.inner_roll_radius;
		outer_roll_radius = ce.outer_roll_radius;
		r_ring_outer_draw = ce.r_ring_outer_draw;
		r_ring_inner_draw = ce.r_ring_inner_draw;
		t_ring_draw = ce.t_ring_draw;      

		random_value = ce.random_value;
		random_noise_factor = ce.random_noise_factor;
	}

	void SetNonlinearK(MathFunction mf_kxxi, MathFunction mf_kyyi)
	{
		mf_kxx = mf_kxxi;
		mf_kyy = mf_kyyi;
		k_nonlinear = 1;
	}
	//1,2,3 .. rotation about global x,y,z-axis
	void SetRotComp(int rotcompi) {rotcomp = rotcompi;}
	void SetDrawBearing(int draw) {draw_bearing = draw;}

	void SetRandomNoiseFactor(double random_noise_factor_i) {random_noise_factor = random_noise_factor_i;}

	//needs rotcomp and rotelements!!!!
	//sets a harmonic fault, higher frequency than relative rotation speed!
	void SetBearingFault(double bearing_fault_freqfact_i, double bearing_fault_amp_i, double bearing_fault_power_i, int innerring = 0)
	{
		bearing_fault_freqfact = bearing_fault_freqfact_i;
		bearing_fault_amp = bearing_fault_amp_i;
		bearing_fault_power = bearing_fault_power_i;
		use_bearing_fault = 1;
		bearing_fault_innerring = innerring;
	}

	void SetBearingDimensions(int n_bearing_balls_i, double inner_roll_radius_i, double outer_roll_radius_i, double r_ring_outer_draw_i, double r_ring_inner_draw_i, double t_ring_draw_i)
	{
		n_bearing_balls =   n_bearing_balls_i;  
		inner_roll_radius = inner_roll_radius_i;
		outer_roll_radius = outer_roll_radius_i;
		r_ring_outer_draw = r_ring_outer_draw_i;
		r_ring_inner_draw = r_ring_inner_draw_i;
		t_ring_draw = t_ring_draw_i;  
	}

	virtual const char* GetElementSpec() const {return "SpringDamperBearing";} //this name names MUST coincide with the type names defined in MBS::MBS() !!!!

	virtual void ComputeBearingGeometry(const Vector3D& p1, const Vector3D& p2, const Vector3D& v1, const Vector3D& v2, 
		const double& phi_inner, const double& phi_outer, double& phi_cage, Vector3D& bearing_fault_vec, double& f_bearing_fault) const;

	virtual Vector3D ComputeForce(double t) const;    
	virtual void DrawElement();

  virtual void StartTimeStep(); //function is called when computation of time step is started; do random computations

protected:
	int rotcomp; //component of Rigid3DKardan, which is rotation respectively for bearing, must be same for all rotors, (e.g. 3==rotation around z-axis)
	
	int use_bearing_fault; //activate
	double bearing_fault_freqfact; //frequency factor for harmonic bearing fault
	double bearing_fault_amp; //amplitude for harmonic bearing fault
	double bearing_fault_power; //power for bearing fault of type f= amp * ((sin(omega*t)+1)/2)^power
	//not needed: double bearing_fault_dir; //+1 or -1 for direction of rotation of bearing
	int bearing_fault_innerring; //0..bearing fault on outer ring, 1..on inner ring

	MathFunction mf_kxx; //mathfunction: nonlinearity kxx
	MathFunction mf_kyy; //mathfunction: nonlinearity kyy
	int k_nonlinear; //use nonlinear kxx,kyy

	int n_bearing_balls; //number of balls of bearing: if 0 ==> bearing is not drawn!
	double inner_roll_radius; //inner rolling radius of bearing (omega*r_i = v_ball_i at inner side)
	double outer_roll_radius; //outer rolling radius of bearing (omega*r_o = v_ball_o at outer side), d_balls = outer_roll_radius - inner_roll_radius
	double r_ring_outer_draw; //outer radius for drawing of outer ring
	double r_ring_inner_draw; //inner radius for drawing of inner ring
	double t_ring_draw; //thickness of ring for drawing

	int draw_bearing; //flag, if bearing shall be drawn
	double random_value; //computed in StartTimeStep
	double random_noise_factor; //size of random noise
};



#endif
