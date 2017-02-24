//#**************************************************************
//#
//# filename:             GCMSElement.h
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						
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
 
#ifndef GeneralizedCMSElement3D__H
#define GeneralizedCMSElement3D__H



//this element contains the whole structure of the body (finite element numbers, nodes) in this class.
//all nodes have local indices in the CMSElement
template<class RIGID>
class GCMSElement: public BaseCMSElement<RIGID>
{
public:

	// Constructor, TCMSflag is set
	GCMSElement(MBS* mbsi):BaseCMSElement<RIGID>(mbsi), refnode1(0), refnode2(0), refnode3(0),	
		dofs_per_local_mode(0), flag_rotation_axis(0), rotation_axis(0.), flexible_dofs(0),
		use_rb_constraints(0), use_mode_constraints(0)
	{
		type = (TMBSElement)(type|TCMSflag);
	};

	// Copy constructor
	GCMSElement(const GCMSElement& e):BaseCMSElement<RIGID>(e.mbs)
		{CopyFrom(e);};

	// Constructor,
	// p ..... initial position
	// v ..... initial velocity
	// phi ... initial angle
	// phip .. initial angular velocity
	// nimodesi .. number of internal modes
	// sizei ..... size of frame (for drawing)
	// coli ...... color (for drawing, currently not used)
	GCMSElement(MBS* mbsi, const Vector3D& p, const Vector3D& v, Vector3D phi, Vector3D phip, const int nimodesi,
		const Vector3D& sizei, const Vector3D& coli) : BaseCMSElement<RIGID>(mbsi)
	{
		type = (TMBSElement)(type|TCMSflag);
		solverparameters.DefaultInitialize();
		EVfile = mystr("");
		SetGCMSElement(p, v, phi, phip, nimodesi, sizei, coli);
		use_rb_constraints = 0;
		use_mode_constraints = 0;
	};

	// Set-Function,
	// p ..... initial position
	// v ..... initial velocity
	// phi ... initial angle
	// phip .. initial angular velocity
	// nimodesi .. number of internal modes
	// sizei ..... size of frame (for drawing)
	// coli ...... color (for drawing, currently not used)
	void SetGCMSElement(const Vector3D& p, const Vector3D& v, Vector3D phi, Vector3D phip, const int nimodesi,
		const Vector3D& sizei, const Vector3D& coli);

	// Set-Function using a reference point P to define v and phip
	// p									initial position of origin
	// v_refP							initial velocity of reference point P
	// phi								initial angle
	// phip_refP					initial angular velocity of reference point P
	// ref_node1TOref_pos	vector from RefNode1 to reference point P
	// nimodesi						number of internal modes
	// sizei							size of frame (for drawing)
	// coli								color (for drawing, currently not used)
	void SetGCMSElement(const Vector3D& p, const Vector3D& v_refP, Vector3D phi, Vector3D phip_refP, const Vector3D& ref_node1TOref_pos, const int nimodesi,
		const Vector3D& sizei, const Vector3D& coli);


	// destructor
	~GCMSElement()
	{
	}


	virtual Element* GetCopy()
	{
		Element* ec = new GCMSElement(*this);
		return ec;
	}
	
	virtual void CopyFrom(const Element& e)
	{
		BaseCMSElement<RIGID>::CopyFrom(e);
		const GCMSElement& ce = (const GCMSElement&)e;
		refnode1 = ce.refnode1;
		refnode2 = ce.refnode2;
		refnode3 = ce.refnode3;
		dofs_per_local_mode = ce.dofs_per_local_mode;
		flag_rotation_axis = ce.flag_rotation_axis;
		rotation_axis = ce.rotation_axis;
		flexible_dofs = ce.flexible_dofs;
		use_rb_constraints =  ce.use_rb_constraints;
		use_mode_constraints = ce.use_mode_constraints;
	}

	// initialization of  FFRF matrices is done here
	virtual void Initialize();
	virtual int PerformNodeCheck() const {return 0;}
	virtual const char* GetElementSpec() const {return "GCMSElement";}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//must be called before Assemble!!!!!
	//link elements, ltg-elements, compute nbmodes, compute M and K, modal analysis, transformation matrix, modal matrices
	virtual void DoModalAnalysis(const TArray<int2>& fixednodes);
	//final matrices of FFRF-formulation:
	// compute mass matrix of CMS element (not constant!)
	virtual void EvalM(Matrix& m, double t); 
	// stiffness+quadratic velocity vector:
	virtual void EvalF2(Vector& f, double t);
	// quadratic velocity vector
	//virtual void AddQuadraticVelocityVector(Vector& f, double t, Vector& temp);
	virtual void EvalG(Vector& g, double t);
		

	// can be set to 2 if the stiffness matrix routine works fine (quadratic velocity vector missing currently)
	virtual int FastStiffnessMatrix() const;
	// returns Kr
	virtual void StiffnessMatrix(Matrix& m); //fill in sos x sos components only of stiffness matrix, m might be larger


	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// number of (local) mode shapes
	virtual int NModes() const {return NBModes()+NIModes();}
	// number of dofs induced by each local mode (9 for general rotation, 3 for rotation around axis
	virtual int& GetNDofPerLocalMode() {return dofs_per_local_mode; }
	virtual const int& GetNDofPerLocalMode() const {return dofs_per_local_mode; }

	// number of reduced dofs
	// number of rigid body motion dofs --> 12 or 5
	virtual int SOSRigid() const {if (flag_rotation_axis) {return 5;} else {return 12;} };
	// total number of reduced dofs
	//size of second order equations, len(u), number of modes*number of dofs per local mode + number of rigid body dofs
	virtual int SOS() const {return FlexDOF()+SOSRigid();};
	//size of second order dofs, len(u), which are not taken from external mbs nodes (all dofs in this case)
	virtual int SOSowned() const {return SOS();}; 
	virtual int SOSowned_RS() const {return SOS()-SOSRigid();}; //size of second order equations for resorting???
	virtual int IS_RS() const {return SOSRigid();}; // implicit size  for resorting???
	virtual int IS() const //implicit (algebraic) size
	{
		if (use_rb_constraints && !use_mode_constraints) 
		{return 6;} 
		else if (!use_rb_constraints && !use_mode_constraints) 
		{return 0;} 
		else if (!use_rb_constraints && use_mode_constraints) 
		{return NModes()*(dofs_per_local_mode-1); }
		else if (use_rb_constraints && use_mode_constraints) 
		{return 6+NModes()*(dofs_per_local_mode-1); }
		assert(0); return 0; //JG 2011-03-22

	};  
	virtual int ES() const {return 0;};  //size of first order explicit equations


	// flag if 6 constraints for rigid body motion shall be used
	virtual const int& UseRigidBodyConstraints() const {return use_rb_constraints;}
	virtual void SetUseRigidBodyConstraints(int use_rbc);
	// flag if 8 constraints for each 9 modes stemming from one local mode shall be used
	virtual const int& UseFlexibleModesConstraints() const {return use_mode_constraints;}
	virtual void SetUseFlexibleModeConstraints(int use_fmc);


	// number of flexible degrees of freedom
	virtual int FlexDOF() const {return flexible_dofs;}

	virtual int IsCMS() const {return 1;}
	virtual int IsGCMS() const {return 1;}

	virtual int& RefNode1() {return refnode1;}
	virtual int& RefNode2() {return refnode2;}
	virtual int& RefNode3() {return refnode3;}
	virtual const int& RefNode1() const {return refnode1;}
	virtual const int& RefNode2() const {return refnode2;}
	virtual const int& RefNode3() const {return refnode3;}

	virtual void ApplyRotation(const Matrix3D& rot, Vector& u) const;
	virtual void ApplyRotationFromLeft(const Matrix3D& rot, Matrix& K) const;
	virtual void ApplyRotationFromRight(const Matrix3D& rot, Matrix& K) const;
	virtual void ComputeURigid(const Matrix3D& A, const Vector3D& uref1, Vector& urigid) const;
// compute derivative of rotation matrix entries with respect to all reduced degrees of freedom
	// analytically
	virtual void ComputeDADq(ConstVector<CMSmaxDOF> dAijdq[3][3]) const;
// compute derivative of rotation matrix entries with respect to reduced degree of freedom q_k
// using numerical differentiation
	virtual void ComputeDADqk_NumericDiff(int k, Matrix3D& dAdq);
// compute rotation matrix for actual set of reduced degrees of freedom
	virtual Matrix3D GetRotMatrix() const;
// compute rotation matrix for given set of reduced degrees of freedom xgc
	virtual Matrix3D GetRotMatrix(const Vector& xgc) const;
	// compute rotation matrix for drawing reduced degrees of freedom
	virtual Matrix3D GetRotMatrixD() const;
	// compute translation vector for given set of reduced degrees of freedom xgc
	virtual Matrix3D GetRotMatrixP() const;
	// compute angular velocity vector
	virtual Vector3D GetAngularVel(const Vector3D& p_loc) const ;
	virtual Vector3D GetAngularVelLocal(const Vector3D& p_loc) const ;
	virtual Vector3D GetAngularVelLocal() const { return GetAngularVelLocal(Vector3D(0.,0.,0.)); }	//AH
	virtual Vector3D GetTranslation(const Vector& xgc) const;
	virtual Vector3D GetRefPosD() const;
	virtual Vector3D GetRefPos() const;
	virtual Vector3D GetRefPosInit() const; 

	virtual const double& XGLastStep(int iloc) const {return mbs->GetLastSolVector()(ltg.Get(iloc));}
	//virtual double& XGLastStep(int iloc) {return mbs->GetLastSolVector()(ltg.Get(iloc));}

	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// position / velocity of node in the CMS element
// compute position vector of node for actual/drawing set of reduced degrees of freedom
	virtual Vector3D GetNodePos(int i) const ;
	virtual Vector3D GetNodePosD(int i) const;
// compute position vector of node for given set of reduced degrees of freedom xgc
	virtual Vector3D GetNodePos(int i, const Vector& xgc) const ;
	// compute position vector of node due to 12 RB-Modes, other modes are neglected
	virtual Vector3D GetNodePosRBModes(int i, const Vector& xgc) const ;
	virtual Vector3D GetNodeLocPos(int i) const ;
	virtual Vector3D GetNodeVel(int i, const Vector& xgp) const ;
	virtual Vector3D GetNodeVel(int i) const ;
	virtual void DrawElement();

	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// position of point in the GCMS elemnt according to rigid body motion only
	virtual Vector3D GetPos(const Vector3D& p_loc) const;
	virtual Vector3D GetVel(const Vector3D& p_loc) const;
	virtual Vector3D GetPosD(const Vector3D& p_loc) const;
	virtual Vector3D GetVelD(const Vector3D& p_loc) const;
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// d = d pos / dq(ploc)
	// GetdPosdqt: only for rigid body position of rigid frame!!
	virtual void GetdPosdqT(const Vector3D& ploc, Matrix& d); 
	// GetdPosdqt: only for rigid body rotation of rigid frame!!
	virtual void GetdRotvdqT(const Vector3D& vloc, const Vector3D& ploc, Matrix& d);
	virtual void GetdRotdqT(const Vector3D& ploc, Matrix& d);

	// differentiation of global position of Node node with respect to REDUCED degrees of freedom!
	virtual void GetNodedPosdqT(int node, Matrix& dpdqi);
	// differentiation of position of Node node due to the 12 RB-Modes
	virtual void GetNodedPosdqTRBModes(int node, Matrix& dpdqi) const;
	virtual void AddNodedPosdqTLambda(int node, const Vector3D& lambda, Vector& f);   // f += dpdq*lambda

	// int du/dq dV for the full cms element, nonlinear, depends on xg
	virtual void GetIntDuDq(Matrix& dudq); //in fact it is DuDq Transposed
	// int rho du/dq dV for the full cms element, nonlinear, depends on xg
	virtual void GetIntRhoDuDq(Matrix& rhodudq); //in fact it is DuDq Transposed
	// Phi_CB^T H, where H is the assembled H-matrix of the FFRFelements
	virtual void GetH(Matrix& H)
	{
		H = Hr;
	}

	virtual Vector3D GetAngularMomentum(const Vector3D& p_ref) const;
	virtual double GetKineticEnergy();
	virtual double GetPotentialEnergy();

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	virtual void EigenModesToDataManager();

	protected:
	// dofs per local mode in theory, in reality less dofs may be used due to linear dependence of modes
	int dofs_per_local_mode;

	int refnode1, refnode2, refnode3;

	// use 6 constraints for rigid body motion
	int use_rb_constraints;
	// use 8 constraints for each 9 modes generated from one local mode
	int use_mode_constraints;

	// number of flexible modes
	int flexible_dofs;

	int flag_rotation_axis;
	Vector3D rotation_axis;

	Vector xg;
};



#endif

