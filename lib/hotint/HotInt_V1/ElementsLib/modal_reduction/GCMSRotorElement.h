//#**************************************************************
//#
//# filename:             GCMSRotorElement.h
//#
//# project:              
//#
//# author:               Pascal Ziegler, Alexander Humer
//#
//# generated:						2012
//# description:          Reference frame element for FFRF formulation
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
 
#pragma once

#include "GCMSElement.h"

template<class RIGID>
class GCMSRotorElement: public GCMSElement<RIGID>
{
public:
	// Constructor, TCMSflag is set
	GCMSRotorElement(MBS* mbsi): GCMSElement<RIGID>(mbsi)
	{
		dofs_per_local_mode = 3;
	};

	// Copy constructor
	GCMSRotorElement(const GCMSRotorElement& e): GCMSElement<RIGID>(e.mbs) {	CopyFrom(e);	};

	// Constructor,
	// p ..... initial position
	// v ..... initial velocity
	// phi ... initial angle
	// phip .. initial angular velocity
	// nimodesi .. number of internal modes
	// sizei ..... size of frame (for drawing)
	// coli ...... color (for drawing, currently not used)
	/*GCMSRotorElement(MBS* mbsi, const Vector3D& p, const Vector3D& v, Vector3D phi, Vector3D phip, const int nimodesi,
		const Vector3D& sizei, const Vector3D& coli) : GCMSElement<RIGID>(mbsi)
	{
		type = (TMBSElement)(type|TCMSflag);
		solverparameters.DefaultInitialize();
		EVfile = mystr("");
		SetGCMSElement(p, v, phi, phip, nimodesi, sizei, coli);
		use_rb_constraints = 0;
		use_mode_constraints = 0;
	};*/

	// Set-Function,
	// p ..... initial position
	// v ..... initial velocity
	// phi ... initial angle
	// phip .. initial angular velocity
	// nimodesi .. number of internal modes
	// sizei ..... size of frame (for drawing)
	// coli ...... color (for drawing, currently not used)
	void SetGCMSRotorElement(const Vector3D& p, const Vector3D& v, Vector3D phi, Vector3D rotation_axis, double omega, const int nimodesi,
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
	/*void SetGCMSRotorElement(const Vector3D& p, const Vector3D& v_refP, Vector3D phi, Vector3D phip_refP, const Vector3D& ref_node1TOref_pos, const int nimodesi,
		const Vector3D& sizei, const Vector3D& coli);*/


	// destructor
	~GCMSRotorElement()
	{
	}


	virtual Element* GetCopy()
	{
		Element* ec = new GCMSRotorElement(*this);
		return ec;
	}
	
	virtual void CopyFrom(const Element& e)
	{
		GCMSElement<RIGID>::CopyFrom(e);
		const GCMSRotorElement& ce = (const GCMSRotorElement&)e;
		this->rotation_axis = ce.rotation_axis;

		/*BaseCMSElement<RIGID>::CopyFrom(e);
		const GCMSElement& ce = (const GCMSElement&)e;
		refnode1 = ce.refnode1;
		refnode2 = ce.refnode2;
		refnode3 = ce.refnode3;
		dofs_per_local_mode = ce.dofs_per_local_mode;
		flag_rotation_axis = ce.flag_rotation_axis;
		rotation_axis = ce.rotation_axis;
		flexible_dofs = ce.flexible_dofs;
		use_rb_constraints =  ce.use_rb_constraints;
		use_mode_constraints = ce.use_mode_constraints;*/
	}

	virtual const char* GetElementSpec() const {return "GCMSRotorElement";}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//must be called before Assemble!!!!!
	//link elements, ltg-elements, compute nbmodes, compute M and K, modal analysis, transformation matrix, modal matrices
	virtual void DoModalAnalysis(const TArray<int2>& fixednodes);
	//compute residuum
	virtual void EvalF2(Vector& f, double t);

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// number of (local) mode shapes
	virtual int NModes() const {return NBModes()+NIModes();}
	// number of dofs induced by each local mode (9 for general rotation, 3 for rotation around axis
	virtual int& GetNDofPerLocalMode() {return dofs_per_local_mode; }
	virtual const int& GetNDofPerLocalMode() const {return dofs_per_local_mode; }

	// number of reduced dofs
	// number of rigid body motion dofs --> 12 or 5
	virtual int SOSRigid() const { return 5; };
	// total number of reduced dofs
	//size of second order equations, len(u), number of modes*number of dofs per local mode + number of rigid body dofs
	virtual int SOS() const {return FlexDOF()+SOSRigid();};
	//size of second order dofs, len(u), which are not taken from external mbs nodes (all dofs in this case)
	virtual int SOSowned() const {return SOS();}; 
	virtual int SOSowned_RS() const {return SOS()-SOSRigid();}; //size of second order equations for resorting???
	virtual int IS_RS() const {return SOSRigid();}; // implicit size  for resorting???
	virtual int IS() const //implicit (algebraic) size
	{
		return 0;
	};  
	virtual int ES() const {return 0;};  //size of first order explicit equations

	// rotation matrix A = A1 + A2 cos(phi) + A3 sin(phi)
	virtual Matrix3D GetA1();
	virtual Matrix3D GetA2();
	virtual Matrix3D GetA3();

	virtual Matrix3D GetRotMatrix() const;
	virtual Matrix3D GetRotMatrix(const Vector& xgc) const;
	virtual Matrix3D GetAbar(const Matrix3D& A) const;

	
	virtual int FastStiffnessMatrix() const;

	virtual void ApplyRotation(const Matrix3D& A, Vector& u) const;
	virtual void ApplyRotationFromLeft(const Matrix3D& A, Matrix& K) const;
	virtual void ApplyRotationFromRight(const Matrix3D& A, Matrix& K) const;
	virtual void ApplyDADAkl(int k, int l, Vector& u) const;

	virtual void GetSinCosPhi(const Matrix3D& A, double* sinphi, double* cosphi) const;
	virtual void ComputeURigid(const Matrix3D& A, const Vector3D& uref1, Vector& urigid) const;

	virtual void ComputeDADq(ConstVector<CMSmaxDOF> dAijdq[3][3]) const;
	virtual void ComputeDADqk_NumericDiff(int k, Matrix3D& dAdq) const;

	virtual void GetNodedPosdqT(int node, Matrix& dpdqi);
	virtual Vector3D GetNodePos(int i, const Vector& xgc) const; 
	virtual Vector3D GetNodeVel(int i, const Vector& xgp) const;

protected:
	Vector3D rotation_axis;
	Matrix PhiCBallmodes;
	Matrix project2indep;	// projects dependent local modes to linear independent set
};