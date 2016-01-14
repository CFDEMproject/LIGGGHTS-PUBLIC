//#**************************************************************
//#
//# filename:             FiniteElement2D.h
//#
//# author:               Gerstmayr Johannes, Aigner Larissa, YV, PG
//#
//# generated:						November 2010
//# description:          2D - general finite element
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

#include "FiniteElementGeneric.h"


/// Base class for finite elements in 2D plane.
/// Local coordinates in the element vary in the range -1...+1
class FiniteElement2D : public FiniteElementGeneric<Body2D>
{
public:
	///@name Creation, initialization and copying of an element
	//@{
	/// Default constructor
	FiniteElement2D(MBS * mbs) : FiniteElementGeneric<Body2D>(mbs) {}
	/// Initialization
	void SetFiniteElement2D(
		int bodyindi,					///< [in] identifier of the body
		const TArray<int>& nodelist,	///< [in] list of the indices of the nodes of the element in the global array
		int material_num,				///< [in] index of the material of the element in the global array
		double thickness,					// thickness of the plate
		const Vector3D& coli			///< [in] color of the graphical representation
		);
	/// Filling in the data from an existing element
	virtual void CopyFrom(const Element& e);
	//@}

	///@name Late initialization
	//@{
	virtual void Initialize(); ///< After assembly and after initial conditions have been set or loaded
	//@}

	///@name FE system matrices and gradients for computation and drawing
	//@{
	/// Get the geometric finite element jacobian \b J for the point with given local coordinates
	virtual void GetJacobi(Matrix2D& jac, const Vector2D& ploc) const;
	/// Compute DS matrix and apply element jacobian inverse for the point with given local coordinates;\n
	/// also compute determinate of jacobian;\n
	/// jacinvDS should have already the dimensions Dim() x NS()
	///@return Determinant of the jacobian
	virtual double GetJacInvDS(const Vector2D& ploc, Matrix& jacinvDS) const;
	virtual double GetJacInvDS(const Vector3D& ploc, Matrix& jacinvDS) const { return GetJacInvDS((const Vector2D &) ploc, jacinvDS); }
	/// Compute displacement gradient gradu for given result of GetJacInvDS()
	virtual void Gradu(
		const Vector& u,		///< Nodal coordinates
		const Matrix& jacinvDS,	///< Pre-computed result of GetJacInvDS()
		Matrix2D& gradu			///< [out] result of computation
		) const;
	/// Compute displacement gradient at a point with given local coordinates
	virtual void Gradu(
		const Vector2D& ploc,   ///< Local coordinates
		const Vector& u,		///< Nodal coordinates
		Matrix2D& gradu) const;	///< [out] result of computation
	/// Compute displacement gradient at a point with given local coordinates
	/// for the actual "drawing" values of degrees of freedom of the element, see Element::XGD().
	virtual void GraduD(const Vector2D& ploc, Matrix2D& gradu) const;
	//@}

	//use sparse methods to write mass/stiffness matrix --> may be faster for sparse systems!
	
	//virtual int UseSparseK() const {return 0;}; //stiffness matrix is not calculated yet
	
	virtual int Dim() const { return 2; }

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// SHAPE FUNCTIONS
	virtual double GetS0(const Vector2D& ploc, int shape) const = 0;
	//get derivative of shape functions w.r.t. all local coordinates:
	virtual void GetDSMatrix0(const Vector2D& ploc, Matrix& sm) const = 0;
	//get "dxj"-th derivative of "shape"-th shape function
	//virtual double GetDS0(const Vector3D& ploc, int shape, int dxj) const = 0;
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	virtual int FastStiffnessMatrix() const { return 2; }
	virtual void StiffnessMatrix(Matrix& m); //fill in sos x 2*sos components, m might be larger
	
	virtual void ComputeMass();

	// in the functions, which compute positions/displacements/velocities,
	// we preserve the structure of the 3D case - YV
	// positions (determined by the own deformation of the element),
	//-1..+1 based!!!
	virtual Vector2D GetDisplacement2DRel(const Vector2D& p_loc) const;
	virtual Vector2D GetDisplacement2DRelD(const Vector2D& p_loc) const;
	virtual Vector2D GetPos2DRel(const Vector2D& p_loc) const;
	virtual Vector2D GetVel2DRel(const Vector2D& p_loc) const;
	virtual Vector2D GetPos2DRelD(const Vector2D& p_loc, int use_magnification) const;
	virtual Vector2D GetVel2DRelD(const Vector2D& p_loc) const;

	//position in the reference configuration:
	//-1..+1 based!!!
	virtual Vector2D GetRefConfPos2D(const Vector2D& p_loc) const;
	virtual Vector2D GetPos2D(const Vector2D& p_loc) const { return GetPos2DRel(p_loc); }
	virtual Vector2D GetDisplacement2D(const Vector2D& p_loc) const { return GetDisplacement2DRel(p_loc); }
	virtual Vector2D GetVel2D(const Vector2D& p_loc) const { return GetVel2DRel(p_loc); }
	virtual Vector2D GetPos2DD(const Vector2D& p_loc, int use_magnification) const { return GetPos2DRelD(p_loc, use_magnification); }
	virtual Vector2D GetPos2DD(const Vector2D& p_loc) const { return GetPos2DD(p_loc, 1); }
	//virtual Vector3D GetPosD(const Vector3D& p_loc, double factor) const;
	virtual Vector2D GetDisplacement2DD(const Vector2D& p_loc) const { return GetDisplacement2DRelD(p_loc); }
	virtual Vector2D GetVel2DD(const Vector2D& p_loc) const { return GetVel2DRelD(p_loc); }
	//reference position (for estimates), -1..+1 based!!!
	virtual Vector2D GetRefPos2D() const {return GetPos2D(Vector2D(0.));};
	//reference position (for graphical issues), -1..+1 based!!!
	virtual Vector2D GetRefPos2DD() const {return GetPos2DD(Vector2D(0.));};

	//assign a position to each DOF (for drawing of DOF-constraints)
	virtual Vector3D GetDOFPosD(int idof) const //returns postion of i-th DOF
	{
		int node = (idof-1)/2+1;
		return ToP3D(GetPos2DD(GetNodeLocPos2D(node)));
	}
	//assign a direction of action to each DOF (for drawing of DOF-constraints)
	virtual Vector3D GetDOFDirD(int idof) const //returns direction of action of i-th DOF
	{
		int dir = (idof-1)%2;
		if (dir == 0)
			return Vector3D(1.,0.,0.);
		else
			return Vector3D(0.,1.,0.);
	}

	//return box in which the element fits
	virtual Box3D GetElementBox() const
	{
		Box3D b;
		for (int i=1;i <= NNodes(); i++)
		{
			b.Add(ToP3D(GetNodePos2D(i)));
		}
		return b;
	}

	virtual Box3D GetElementBoxD() const
	{
		Box3D b;
		for (int i=1;i <= NNodes(); i++)
		{
			b.Add(ToP3D(GetNodePos2DD(i)));
		}
		return b;
	}

	/// Integral over the element of the displacement with respect to the element coordinates.
	/// Integral over the element of the displacement with respect to the element coordinates \f$\bf{q}\f$
	/// and Jacobian \f$\bf{J}\f$.\n
	/// For finite elements wish shape matrix \f$\bf{S}\f$ and displacements \f$\bf{u}=\bf{Sq}\f$
	/// this reads \f$\bf{H}=\int_{V_{elem}}\bf{S}\,dV_{elem}.\f$
	virtual void GetH(Matrix& H);

	/// Evaluate Mass matrix
	virtual void EvalM(Matrix& m, double t) { EvalMff(m, t); }			// will be overwritten in FFRF case
	// Evaluate Mass matrix for flexible dofs, Attention: Size of mass matrix is not changed, is not initialized with zeros!
	// Size of mass matrix is SOS(), also for Floating Frame formulation!
	virtual void EvalMff(Matrix& m, double t);

	virtual void EvalF2(Vector& f, double t);

	// for volumeloads (gravity ...)
	//in fact it is DuDq Transposed
	virtual void GetIntDuDq(Matrix& dudq) { GetH(dudq); }

	// for centrifugal force as MBSLoad - parameters omega, r0 is point on rotation axis (AD 2012-07-09)
	virtual void GetIntDuDqFCentrifugal(Matrix& dudq, const Vector3D& omega, const Vector3D& r0); 

	virtual void GetdPosdqT(const Vector2D& ploc, Matrix& d);

	virtual void GetNodedPosdqT(int node, Matrix& dpdqi);

	virtual void AddNodedPosdqTLambda(int node, const Vector2D& lambda, Vector& f) 
	{
		f((node-1)*2+1) += lambda(1);
		f((node-1)*2+2) += lambda(2);
	};   // f += dpdq*lambda

	//local node position is generated in the actual element classes (Tri, Quad, ...)
	virtual Vector2D GetNodeLocPos2D(int i) const = 0;

	virtual Vector2D GetNodeRefPos2D(int i) const //return nodal position of node i in reference configuration
	{
		return GetNode(i).RefConfPos2D();
	}

	virtual Vector2D GetNodePos2D(int i) const;
	virtual Vector2D GetNodePos2DD(int i) const;
	virtual Vector2D GetNodeVel2D(int i) const;

	// variables, available for post-processing and sensing
	virtual void GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables);
	// computation of the variables
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector2D & local_position, bool flagD);
	// these are the helper functions for the computation of the field variables
	// certain entities may be defined only in the integration points of the stiffness matrix -
	// - if necessary, this function shifts the original point to the nearest integration point
	// and computes initial strains and plastic strains there
	void UpdateFieldVariableComputationPoint(const Vector2D & local_position_original, Vector2D & local_position_actual,
		Matrix2D & initial_strain, Vector & inelastic_variables);

	//for interpolated stress plots!
	virtual void DrawElementPreProc();

	virtual void DrawElement();


	// NONLINEAR MATERIAL LAWS
	virtual double PostNewtonStep(double t);
	virtual void PostprocessingStep();

	// this function tests whether there are real initial strains on the element - may be redefined in specialized derived classes
	virtual bool HasInitialStrains() { return false; }
	// the initial strains may be interpolated over the nodes in a derived class
	virtual void GetInitialStrain(Matrix2D & initial_strain, const Vector2D & ploc, bool flagDraw) 
	{
		mbs->UO() << "***Warning: NO inital strain computation implemented so far!\n";
	}
	
	//EK 2012-07-10 integrate the function given by the values of the shape functions in f
	double L2InnerProduct(Vector& f1, Vector &f2);

	// $EK 20130408 added in order to be consistent to constraints
	virtual TKinematicsAccessFunctions GetKinematicsAccessFunctions(int mode = 1) const
	{
		// $EK 20130416 not correct to use inheritage here -> addition twice results in the next higher bit to be set!!!!
		TKinematicsAccessFunctions tkf = FiniteElementGeneric<Body2D>::GetKinematicsAccessFunctions(mode);
		//TKAF_displacement .... only GetDisplacement2D available
		//TKAF_node_position ... only GetNodePos2D available
		//tkf = TKinematicsAccessFunctions(tkf);
		//TKAF_node_velocity ... only GetNodeVel2D available
		//tkf = TKinematicsAccessFunctions(TKAF_node_velocity + tkf);
		tkf = TKinematicsAccessFunctions(TKAF_node_position + TKAF_ref_conf_position + TKAF_D_pos_D_q + TKAF_D_node_pos_D_q + TKAF_int_D_u_D_q + tkf);
		//TKAF_position+TKAF_displacement+TKAEF_velocity bereits definiert in Element (Konflikt ????)
		return tkf;
	}
protected:
	double thickness;
}; 
