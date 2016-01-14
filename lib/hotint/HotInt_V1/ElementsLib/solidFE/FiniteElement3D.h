//#**************************************************************
//#
//# filename:             FiniteElement3D.h
//#
//# author:               Gerstmayr Johannes, Aigner Larissa, YV
//#
//# generated:						March 2007
//# description:          3D - general finite element
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

#include "body3d.h"
#include "FiniteElementGeneric.h"


/// Face of a spatial finite element.
/// A face of a finite element is defined by a set of node indices
/// and may have either 4 or 3 nodes.
class FEFace
{
public:
	/** @name Constructors */
	//@{
	FEFace() {};
	/// constructor for a face with four nodes
	FEFace(int v1, int v2, int v3, int v4)
	{
		f[0] = v1; f[1] = v2; f[2] = v3; f[3] = v4; 
		nfacenodes = 4;
	}
	/// constructor for a face with three nodes
	FEFace(int v1, int v2, int v3)
	{
		f[0] = v1; f[1] = v2; f[2] = v3;
		nfacenodes = 3;
	}
	/// copy constructor
	FEFace(const FEFace& face)
	{
		f[0] = face.f[0];
		f[1] = face.f[1];
		f[2] = face.f[2];
		f[3] = face.f[3];
		nfacenodes = face.nfacenodes;
	}
	//@

	/** @name Extracting information concerning the nodes */
	//@{
	virtual const int& Node(int i) const {return f[i-1];}
	virtual int& Node(int i) {return f[i-1];}
	virtual const int& NodeMod(int i) const {return f[(i-1)%NFaceNodes()];}
	virtual int& NodeMod(int i) {return f[(i-1)%NFaceNodes()];}
	//@}

	/** @name Obtaining general information concerning the face */
	//@{
	/// number of nodes
	virtual int NFaceNodes() const {return nfacenodes;}
	virtual int IsTrig() const {return nfacenodes == 3;}
	virtual int IsQuad() const {return nfacenodes == 4;}
	//@}

	/** @name Other general purpose functions */
	//@{

	/// inverts the order of nodes
	virtual void Invert()
	{
		if (NFaceNodes() == 4)
		{
			int ii = f[1];
			f[1] = f[3];
			f[3] = ii;
		}
		else if (NFaceNodes() == 3)
		{
			Swap(f[1],f[2]);
		}
	}

	/// returns 1 if both faces are exactly the same
	virtual int IsEqual(const FEFace& other) const
	{
		//if faces have different number of nodes, they can not be equal
		if (NFaceNodes() != other.NFaceNodes()) return 0;

		for (int i=0; i<NFaceNodes(); i++)
		{
			if (f[i] != other.f[i]) return 0;
		}
		return 1;
	}

	/// returns 1 if both faces are the same or if any cyclic permutation is the same
	virtual int IsCyclicEqual(const FEFace& other) const
	{
		//if faces have different number of nodes, they can not be equal
		if (NFaceNodes() != other.NFaceNodes()) return 0;

		for (int j=0; j<NFaceNodes(); j++)
		{
			int equal = 1;
			for (int i=1; i<=NFaceNodes(); i++)
			{
				if (NodeMod(i+j) != other.NodeMod(i)) 
				{
					equal = 0;
					break;
				}
			}
			if (equal) return 1;
		}
		return 0;
	}

	//@}

private:
	int nfacenodes;
	static const int FEMaxFaceNodes = 4;
	int f[FEMaxFaceNodes];
};

/// Base class for finite elements in 3D space.
/// Provides basic functionality of ANCF elements
/// and serves as a base class for FFRF and CMS elements.
/// Local coordinates in the element vary in the range -1...+1
class FiniteElement3D : public FiniteElementGeneric<Body3D>
{
public:
	///@name Creation, initialization and copying of an element
	//@{
	/// Default constructor
	FiniteElement3D(MBS * mbs) : FiniteElementGeneric<Body3D>(mbs) {}
	/// Initialization
	void SetFiniteElement3D(
		int bodyindi,					///< [in] identifier of the body
		const TArray<int>& nodelist,	///< [in] list of the indices of the nodes of the element in the global array
		int material_num,				///< [in] index of the material of the element in the global array
		const Vector3D& coli,			///< [in] color of the graphical representation
		int CMSelnumi=0					///< [in] [optional] number of the corresponding CMS element
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
	virtual void GetJacobi(Matrix3D& jac, const Vector3D& ploc) const;
	/// Compute DS matrix and apply element jacobian inverse for the point with given local coordinates;\n
	/// also compute determinate of jacobian;\n
	/// jacinvDS should have already the dimensions Dim() x NS()
	///@return Determinant of the jacobian
	virtual double GetJacInvDS(const Vector3D& ploc, Matrix& jacinvDS) const;
	/// Compute displacement gradient gradu for given result of GetJacInvDS()
	virtual void Gradu(
		const Vector& u,		///< Nodal coordinates
		const Matrix& jacinvDS,	///< Pre-computed result of GetJacInvDS()
		Matrix3D& gradu			///< [out] result of computation
		) const;
	/// Compute displacement gradient at a point with given local coordinates
	virtual void Gradu(
		const Vector3D& ploc,   ///< Local coordinates
		const Vector& u,		///< Nodal coordinates
		Matrix3D& gradu) const;	///< [out] result of computation
	/// Compute displacement gradient at a point with given local coordinates
	/// for the actual "drawing" values of degrees of freedom of the element, see Element::XGD().
	virtual void GraduD(const Vector3D& ploc, Matrix3D& gradu) const;
	//@}

	//use sparse methods to write mass/stiffness matrix --> may be faster for sparse systems!
	
	//virtual int UseSparseK() const {return 0;}; //stiffness matrix is not calculated yet
	
	virtual int Dim() const { return 3; }

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// SHAPE FUNCTIONS
	virtual double GetS0(const Vector3D& ploc, int shape) const = 0;
	//get derivative of shape functions w.r.t. all local coordinates:
	virtual void GetDSMatrix0(const Vector3D& ploc, Matrix& sm) const = 0;
	//get "dxj"-th derivative of "shape"-th shape function
	//virtual double GetDS0(const Vector3D& ploc, int shape, int dxj) const = 0;
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	virtual int FastStiffnessMatrix() const { return 2; }
	
	virtual void StiffnessMatrix(Matrix& m); //fill in sos x 2*sos components, m might be larger
	
	virtual void ComputeMass();
	// quadratic velocity vector - geometric stiffening effect is not present unless it is FFRF
	virtual void AddQuadraticVelocityVector(Vector& fadd, double t) { assert(0); }

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//floating frame of reference formulation: for FFRF elements
	//Position relative to frame, undeformed
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	// positions (determined by the own deformation of the element),
	// which are possibly relative to the reference frame:
	//-1..+1 based!!!
	virtual Vector3D GetDisplacementRel(const Vector3D& p_loc) const;
	virtual Vector3D GetDisplacementRelD(const Vector3D& p_loc) const;
	virtual Vector3D GetPosRel(const Vector3D& p_loc) const;
	virtual Vector3D GetVelRel(const Vector3D& p_loc) const;
	virtual Vector3D GetPosRelD(const Vector3D& p_loc, int use_magnification) const;
	virtual Vector3D GetVelRelD(const Vector3D& p_loc) const;

	//position in the reference configuration:
	//-1..+1 based!!!
	virtual Vector3D GetRefConfPos(const Vector3D& p_loc) const;
	virtual Vector3D GetPos(const Vector3D& p_loc) const { return GetPosRel(p_loc); }
	virtual Vector3D GetDisplacement(const Vector3D& p_loc) const { return GetDisplacementRel(p_loc); }
	virtual Vector3D GetVel(const Vector3D& p_loc) const { return GetVelRel(p_loc); }
	virtual Vector3D GetPosD(const Vector3D& p_loc, int use_magnification) const { return GetPosRelD(p_loc, use_magnification); }
	virtual Vector3D GetPosD(const Vector3D& p_loc) const { return GetPosD(p_loc, 1); }
	//virtual Vector3D GetPosD(const Vector3D& p_loc, double factor) const;
	virtual Vector3D GetDisplacementD(const Vector3D& p_loc) const { return GetDisplacementRelD(p_loc); }
	virtual Vector3D GetVelD(const Vector3D& p_loc) const { return GetVelRelD(p_loc); }
	//reference position (for estimates), -1..+1 based!!!
	virtual Vector3D GetRefPos() const {return GetPos(Vector3D(0.));};
	//reference position (for graphical issues), -1..+1 based!!!
	virtual Vector3D GetRefPosD() const {return GetPosD(Vector3D(0.));};

	// angular velocity of the finite element
	virtual Vector3D GetAngularVel(const Vector3D& ploc) const;
	//virtual double GetAngularVel(const Vector3D &ploc, const Vector3D &ploc2, const Vector3D& glob_axis) const;

	//assign a position to each DOF (for drawing of DOF-constraints)
	virtual Vector3D GetDOFPosD(int idof) const //returns postion of i-th DOF
	{
		int node = (idof-1)/3+1;
		return GetPosD(GetNodeLocPos(node));
	}
	//assign a direction of action to each DOF (for drawing of DOF-constraints)
	virtual Vector3D GetDOFDirD(int idof) const //returns direction of action of i-th DOF
	{
		int dir = (idof-1)%3;
		if (dir == 0) return Vector3D(1.,0.,0.);
		else if (dir == 1) return Vector3D(0.,1.,0.);
		else return Vector3D(0.,0.,1.);
	}

	//$ PG 2011-3-2: utilize Body3D::GetElementBox() and Body3D::GetElementBoxD()
	//return box in which the element fits
	virtual Box3D GetElementBox() const
	{
		Box3D b = Body3D::GetElementBox();
		b.Increase(2);
		/*for (int i=1;i <= NNodes(); i++)
		{
			b.Add(GetNodePos(i));
		}*/
		return b;
	}

	virtual Box3D GetElementBoxD() const
	{
		Box3D b = Body3D::GetElementBoxD();
		b.Increase(2);
		/*for (int i=1;i <= NNodes(); i++)
		{
			b.Add(GetNodePosD(i));
		}*/
		return b;
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//FACES   FACES   FACES   FACES   FACES   FACES   FACES   FACES   FACES   FACES   FACES   FACES   FACES   FACES
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	virtual int NFaces() const = 0;
	virtual const FEFace & GetLocalFace(int i) const = 0;
	virtual FEFace GetFace(int i) const 
	{
		FEFace face = GetLocalFace(i);
		for (int j=1; j<=face.NFaceNodes(); j++)
		{
			face.Node(j) = NodeNum(face.Node(j));
		}
		return face;
	}
	//return if this element has the face 'face' in any cyclic manner
	virtual int HasFace(const FEFace& face) const
	{
		for (int i=1; i <= NFaces(); i++)
		{
			if (face.IsCyclicEqual(GetFace(i))) return 1;
		}
		return 0;
	}
	//compute local face coordinates for drawing or for accessing face surface
	//the vectors v1 and v2 are local vectors which point from refpos to the other nodes
	virtual void GetLocalFaceCoordinates(int face, Vector3D& refpos, Vector3D& v1, Vector3D& v2, Vector3D& v3) const = 0;
	virtual void SetOuterFaces();
	virtual void SetOuterFacesCuttingPlane();
	virtual void ResetOuterFaceFlags() {outer_face = (char)0;}
	char& Outer_Face() {return outer_face;}

	virtual int GetOuterFaceFlag(int face) const
	{
		return (outer_face&GetCharBit(face)) != 0;
	}
	virtual int HasOuterFaces() const {return (int)outer_face != 0;}

	virtual void SetOuterFaceFlag(int face, int set = 1);
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
	virtual void EvalF2GeomLin(Vector& fadd, double t); 
	virtual void EvalF2GeomNonlin(Vector& fadd, double t);

	// for volumeloads (gravity ...)
	//in fact it is DuDq Transposed
	virtual void GetIntDuDq(Matrix& dudq) { GetH(dudq); }		// will be overwritten in FFRF case

	// for centrifugal force as MBSLoad - parameters omega, r0 is point on rotation axis (AD)
	virtual void GetIntDuDqFCentrifugal(Matrix& dudq, const Vector3D& omega, const Vector3D& r0); 

	virtual void GetdPosdqT(const Vector3D& ploc, Matrix& d);

	virtual void GetNodedPosdqT(int node, Matrix& dpdqi);

	virtual void AddNodedPosdqTLambda(int node, const Vector3D& lambda, Vector& f) 
	{
		f((node-1)*3+1) += lambda(1);
		f((node-1)*3+2) += lambda(2);
		f((node-1)*3+3) += lambda(3);
	};   // f += dpdq*lambda

	virtual Vector3D GetNodeRefPos(int i) const; //return nodal position of node i in reference configuration
	virtual Vector3D GetNodePos(int i) const;
	virtual Vector3D GetNodePosD(int i) const;
	virtual Vector3D GetNodeVel(int i) const;

	// variables, available for post-processing and sensing
	virtual void GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables);
	// computation of the variables
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & local_position, bool flagD);

	//for interpolated stress plots!
	virtual void DrawElementPreProc();

	virtual void DrawElement();


	// NONLINEAR MATERIAL LAWS
	virtual double PostNewtonStep(double t);
	virtual void PostprocessingStep();

	// this function tests whether there are real initial strains on the element - may be redefined in specialized derived classes
	virtual bool HasInitialStrains() { return false; }
	// the initial strains may be interpolated over the nodes in a derived class
	virtual void GetInitialStrain(Matrix3D & initial_strain, const Vector3D & ploc, bool flagDraw) {}

	virtual void PreAssemble();

	// $EK 2013-02-12 compute volume via numerical integration of 1
	virtual double ComputeVolume();
	// $EK 2013-03-21 compute center of mass via integration
	virtual void ComputeCenterOfMass(Vector3D& center);
	
	// $EK 20130408 added in order to be consistent to constraints
	virtual TKinematicsAccessFunctions GetKinematicsAccessFunctions(int mode = 1) const
	{
		// $EK 20130416 not correct to use inheritage here -> addition twice results in the next higher bit to be set!!!!
		//TKinematicsAccessFunctions tkf = FiniteElementGeneric<Body3D>::GetKinematicsAccessFunctions(mode);
		TKinematicsAccessFunctions tkf = TKinematicsAccessFunctions(TKAF_position +TKAF_ref_conf_position+ TKAF_displacement + TKAF_node_position+TKAF_node_ref_conf_position);
		tkf = TKinematicsAccessFunctions(TKAF_velocity + TKAF_node_velocity + TKAF_angular_velocity + tkf);
		tkf = TKinematicsAccessFunctions(TKAF_D_pos_D_q + TKAF_D_node_pos_D_q + TKAF_int_D_u_D_q + tkf);
		//TKAF_position+TKAF_displacement+TKAEF_velocity bereits definiert in Element (Konflikt ????)
		return tkf;
	}

protected:
	//drawing:
	char outer_face; //bitwise flags if side 1 - side 6 is an outer (visible) face	
}; 
