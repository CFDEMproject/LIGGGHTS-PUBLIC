 
//#***************************************************************************************
//# filename:     ANCFThinPlate3D.h
//# authors:      Gerstmayr Johannes, Vetyukov Yury
//# generated:    first: 24. October 2006, 
//# rewritten:    January-February 2011
//# description:  Kirchhoff Plate model, large deformation, 36 DOF
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

#include "FiniteElement3D.h"
#include "XGProvider.h"
#include "PlaneSymmetricTensorComponents.h"
	
struct PlasticVariables
{
	PlaneSymmetricTensorComponents strain;
	double hardening_parameter, yield_function;

	PlasticVariables & operator*=(double k);
	PlasticVariables & operator+=(const PlasticVariables & T);
	PlasticVariables & operator-=(const PlasticVariables & T);
	void SetAllZero() { strain.t11=0; strain.t22=0; strain.t12=0; hardening_parameter=0; yield_function=0; }
};

const int ANCFThinPlateMaxIP=25; //5*5 for plate ??

// thin plate/shell element with initial curvarure,
// modeled on the basis of the theory of classical Kirchhoff-Love material surface.
// each element has 4 nodes; the topology of the mesh must be regular,
// as the coordinate lines on the element must be continuous over the element interfaces.
// each node has 9 mechanical degrees of freedom:
// - 3 degrees of freedom - Cartesian components of the position vector
// - 3 components of the derivative of the position vector with respect to the first local coordinate
// - 3 components --//-- with respect to the second local coordinate
// two sizes of the elements can be given, which affect scaling of the "slope" degrees
// of freedom (derivatives of the position vector) and which should improve the convergence
// in the case of variable element sizes;
// the sizes should equal to the actual lengths of the coordinate lines in the reference configuration;
// if both sizes are set to 1 (by default), then the local derivatives in nodes are considered
// to be relative to the local coordinates on the element, which vary from -1 to +1.
// the nodes, which constitute the element, are treated in a non-standard way:
// pos is ignored, and x_init is understood not as displacements, but as absolute coordinates
// in the reference configuration.
// what concerns the elements, XG() are the displacements (0 in the reference state) and absolute velocities,
// which is achieved through x_init; this is better for the geometrically linear analysis.
class ANCFThinPlate3D: public FiniteElementGeneric<Body3D>
{
public:
	ANCFThinPlate3D(MBS* mbsi) : FiniteElementGeneric<Body3D>(mbsi) { ANCFThinPlate3DInitialize(); }
	ANCFThinPlate3D(const ANCFThinPlate3D& e) : FiniteElementGeneric<Body3D>(e.mbs)
	{
		ANCFThinPlate3DInitialize();
		CopyFrom(e);
	}

	// this version of the set function requires that the nodes are already added to mbs;
	// the nodes must have 9 degrees of freedom and be initialized according to the reference configuration
	// (x_init of the node must be set to the position vector and local derivatives)
	void SetANCFThinPlate3D(
		int bodyindi,									// identifier of the body
		const TArray<int>& node_list,	// list of the four indices of the nodes of the element in the global array
		int material_num,							// index of the material of the element in the global array
		double thickness,							// thickness of the plate
		const Vector3D& coli					// color of the graphical representation
		);

	// this version of the set function creates and adds nodes to mbs;
	// each node is defined through its 9 degrees of freedom in the reference configuration
	void SetANCFThinPlate3D(
		int bodyindi,									// identifier of the body
		const TArray<Vector> nodes,		// degrees of freedom of the four nodes
		int material_num,							// index of the material of the element in the global array
		double thickness,							// thickness of the plate
		const Vector3D& coli					// color of the graphical representation
		);

	// sets the characteristic sizes of the finite element in the reference configuration
	// in the directions of both coordinate lines;
	// the sizes are set to 2 by default
	void SetSizes(double size1, double size2);

	// sets the number of thickness layers ( needs to be >= 1, by default = 4 )
	void SetThicknessLayers(int layers);

	// standard copying procedures
	virtual Element* GetCopy()
	{
		Element * ec = new ANCFThinPlate3D(*this);
		return ec;
	}
	virtual void CopyFrom(const Element& e);

	// late initialization
	virtual void Initialize();

	// identification of the element type
	virtual const char* GetElementSpec() const { return "ANCFThinPlate3D"; }
	virtual TFiniteElementType GetElementType() const { return TFE_ThinPlate; }

	// spatial dimension of the element
	virtual int Dim() const { return 3; }
	// number of degrees of freedom in a node
	virtual int DOFPerNode() const { return 9; }
	// number of shape functions
	virtual int NS() const {return 12;}

	// this is the central point of the class: the differential geometry of the surface in the Gauss point,
	// which is used later on for computing the parameters of the actually deformed surface.
	// what the function returns is the jacobian (integration coefficient), in our case
	// it is the determinant of the matrix of components of the first metric tensor;
	// the data, which need to be pre-computed in each Gauss point, is stored in jacinvDS.
	// the argument is Vector3D, as the function is called from the Generic class
	virtual double GetJacInvDS(const Vector3D& ploc, Matrix& jacinvDS) const;
	// this function produces the determinant of the jacobian matrix
	// (J in the relation J*dq1*dq2 = dOmega, Omega being the surface area)
	// for the reference configuration - the same as GetJacInvDS, but without overhead
	double GetJacobianDeterminant(const Vector2D& ploc) const;

	// accessing the element geometric state from outside of the class
	// the local coordinates are -1..+1 based;
	// the first two local coordinates are on the element surface,
	// the third one over the thickness (in the direction of the normal)
	virtual Vector3D GetRefConfPos(const Vector3D& ploc) const { return GetPos(ploc, xgInit); }
	virtual Vector3D GetPos(const Vector3D& ploc) const { return GetPos(ploc, xgCompute); }

	virtual Vector3D GetPosD(const Vector3D& ploc) const { return GetPos(ploc, xgDraw); }
	virtual Vector3D GetPosD(const Vector3D& ploc, int use_magnification) const;
	// reference positions of the element for painting
	virtual Vector3D GetRefPos() const { return GetPos(Vector3D(0.)); }
	virtual Vector3D GetRefPosD() const { return GetPosD(Vector3D(0.)); }
	// the displacements are computed based on the difference from the reference configuration
	virtual Vector3D GetDisplacement(const Vector3D& ploc) const { return GetDisplacement(ploc, xgCompute); }
	virtual Vector3D GetDisplacementD(const Vector3D& ploc) const { return GetDisplacement(ploc, xgDraw); }
	// the velocities are taken just on the middle surface, no rotation with the normal vector
	virtual Vector3D GetVel(const Vector3D& ploc) const { return GetVel(ploc, xgCompute); }
	virtual Vector3D GetVelD(const Vector3D& ploc) const { return GetVel(ploc, xgDraw); }
	// this method computes the threedimensional position regarding the specific thickness layer
	Vector3D GetLocalPositionAtThicknessLayer(const Vector2D & ploc, int thickness_layer, bool flagDraw) const;
	// unit normal vector
	Vector3D GetNormal(const Vector3D& ploc) const { return GetNormal((const Vector2D&)ploc, xgCompute); }
	Vector3D GetNormalD(const Vector3D& ploc) const { return GetNormal((const Vector2D&)ploc, xgDraw); }
	virtual Vector3D GetNodeLocPos(int local_node_number) const;		//$ YV 2012-06

	// positions and effective directions of degrees of freedom (for drawing of DOF-constraints)
	virtual Vector3D GetDOFPosD(int idof) const;
	virtual Vector3D GetDOFDirD(int idof) const;

	virtual TKinematicsAccessFunctions GetKinematicsAccessFunctions(int mode) const
	{
		TKinematicsAccessFunctions kaf = Body3D::GetKinematicsAccessFunctions(mode);
		return TKinematicsAccessFunctions(kaf+TKAF_D_pos_D_q+TKAF_int_D_u_D_q+TKAF_ref_conf_position);
	}

	virtual int GetActualInterpolationOrder() const { return 4; }

	// computation functions

	void ComputeMass();
	/// Integral over the element of the displacement with respect to the element coordinates.
	virtual void GetH(Matrix& H);
	/// Evaluate Mass matrix
	virtual void EvalM(Matrix& m, double t) { EvalMff(m, t); }
	// Evaluate Mass matrix for flexible dofs, Attention: Size of mass matrix is not changed, is not initialized with zeros!
	virtual void EvalMff(Matrix& m, double t);
	virtual void EvalF2(Vector& f, double t); 
	virtual void EvalF2GeomLin(Vector& fadd, double t); 
	virtual void EvalF2GeomNonlin(Vector& fadd, double t);
	virtual void GetdPosdqT(const Vector3D& ploc, Matrix& d);
	virtual int FastStiffnessMatrix() const { return 2; }
	virtual void StiffnessMatrix(Matrix& m); //fill in sos x 2*sos components, m might be larger

	// for volumeloads (gravity ...)
	//in fact it is DuDq Transposed
	virtual void GetIntDuDq(Matrix& dudq) { GetH(dudq); }		// will be overwritten in FFRF case

	// specific plastic computation
	virtual double PostNewtonStep(double t);

	PlasticVariables GetPlasticVariables(const IntegrationPointsIterator & ip, int thickness_layer, const XGProvider & xg,
							const PlaneSymmetricTensorComponents & ainv, const PlaneSymmetricTensorComponents & a, const PlaneSymmetricTensorComponents & b);
	PlasticVariables GetOldPlasticVariables(const IntegrationPointsIterator & ip, int thickness_layer);
	// return plastic variables at x,y of integration point and z linearly interpolated between values in thickness layers
	PlasticVariables GetOldPlasticVariables(const IntegrationPointsIterator& ip, double z, bool flagDraw);
	
	void DataToPlasticVariables(PlasticVariables & plastic_variables, int ip_idx, int layer_idx);
	void DataToPlasticVariablesD(PlasticVariables & plastic_variables, int ip_idx, int layer_idx);
	void PlasticVariablesToData(const PlasticVariables & plastic_variables, int ip_idx, int layer_idx);
	void PlasticVariablesToDataD(const PlasticVariables & plastic_variables, int ip_idx, int layer_idx);
	virtual int DataS() const { return thickness_layers * 5 * integrationRuleStiffness->GetNumberOfIntegrationPoints(); }

	// drawing
	void DrawElement();

protected:
	void ANCFThinPlate3DInitialize();
	double size1, size2;
	int thickness_layers;

	void GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables);
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & local_position, bool flagD);		// for using from the outside
	// for using from the class itself
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & local_position, const XGProvider & xg);

	// "state" providers
	Vector xgReferenceState;
	XGProviderInit xgInit;					// initial positions and velocities
	XGProviderCompute xgCompute;		// "compute" positions and velocities
	XGProviderDraw xgDraw;					// "drawing" positions and velocities

	Vector3D GetNormal(const Vector2D& ploc, const XGProvider & xg) const;
	Vector3D GetPos(const Vector3D& ploc, const XGProvider & xg) const;
	Vector3D GetDisplacement(const Vector3D& ploc, const XGProvider & xg) const;
	Vector3D GetVel(const Vector3D& ploc, const XGProvider & xg) const;

	// shape functions as they are, accessible per index nShapeFunction = 1 .. 12
	// local coordinates q1, q2 vary from -1 to +1
	double GetS0(const Vector2D& ploc, int nShapeFunction) const;
	// derivatives of the shape functions with respect to the local coordinates; alpha = 1,2
	double GetDS0(const Vector2D& ploc, int nShapeFunction, int alpha) const;
	// second order derivatives of the shape functions with respect to the local coordinates
	// alphaBeta = 1 (11), 2 (12), 3 (22)
	double GetDDS0(const Vector2D& ploc, int nShapeFunction, int alphaBeta) const;

	// compute the local covariant basis vector Ralpha (dr/dq^alpha)
	Vector3D GetRAlpha(const Vector2D& ploc, int alpha, const XGProvider & xg) const;
	// compute the second order derivative vector Ralpha (d^2r / (dq^alpha dq^\beta)), alphaBeta as in GetDDS0(...)
	Vector3D GetRAlphaBeta(const Vector2D& ploc, int alphaBeta, const XGProvider & xg) const;
	// compute the derivarive of the local basis vector RAlpha with respect to a particular degree of freedom
	// nDOF = 1..36, alpha = 1,2, the result is independent from the considered configuration
	Vector3D GetDerivativeRAlphaDOF(const Vector2D& ploc, int nDOF, int alpha) const;
	// compute the derivarive of the vector RAlphaBeta with respect to a particular degree of freedom
	// nDOF = 1..36, alphaBeta as in GetDDS0(...), the result is independent from the considered configuration
	Vector3D GetDerivativeRAlphaBetaDOF(const Vector2D& ploc, int nDOF, int alphaBeta) const;
	// compute the derivative of the unit normal vector with respect to a particular degree of freedom
	Vector3D GetDerivativeNormalDOF(const Vector2D& ploc, int nDOF, const XGProvider & xg) const;

	// compute the covariant components of the first metric tensor
	PlaneSymmetricTensorComponents GetMetricTensorAComponents(const Vector2D& ploc, const XGProvider & xg) const;
	// compute the components of the second metric tensor
	PlaneSymmetricTensorComponents GetMetricTensorBComponents(const Vector2D& ploc, const XGProvider & xg) const;
	// compute the components of the first strain tensor (in-plane deformations)
	PlaneSymmetricTensorComponents GetFirstStrainTensorCComponents(
								const Vector2D & ploc, const Matrix & grad, const XGProvider & xg) const;
	// compute the components of the second strain tensor (bending - change of curvature)
	PlaneSymmetricTensorComponents GetSecondStrainTensorKComponents(
								const Vector2D & ploc, const Matrix & grad, const XGProvider & xg) const;	


	// for additional information concerning the recomputation between local and Cartesian
	// components see the document "Cartesian components for shell elements.docx"
	// plane part of the total Lagrangian strain tensor,
	// covariant components in the local basis of the reference configuration;
	// a point is identified by the integration point in the plane and by the local thickness coordinate, which varies from -1 to 1
	PlaneSymmetricTensorComponents ComputeTotalStrainLocalComponents(IntegrationPointsIterator & ip, double localZ, const XGProvider & xg);
	// converts covariant components of an in-plane tensor into a matrix
	// of components in a particular Cartesian basis in the tangent plane;
	// the integration point is required as the recomputation involves the metric of the surface
	Matrix2D LocalComponentsToCartesian(IntegrationPointsIterator & ip, PlaneSymmetricTensorComponents & localComponents);
	// converts Cartesian components of an in-plane tensor into local basis (produces covariant components)
	PlaneSymmetricTensorComponents CartesianComponentsToLocal(IntegrationPointsIterator & ip, Matrix2D & cartesianComponents);

	// thickness and stiffnesses of the cross-section for bending and in plane deformations
	struct PlateCrossSection
	{
		double thickness, E1, E2, D1, D2;
		void ComputeStiffnesses(double YoungModulus, double PoissonRatio);
	};
	PlateCrossSection plateCrossSection;

	// actual access to the plate thickness should be performed via this function,
	// which would allow to derive a class with variable thickness;
	// thickness should be specified for the "compute" or "draw" states
	virtual double GetThicknessAtPoint(const Vector2D & ploc, bool flagDraw) const { return plateCrossSection.thickness; }

	// this function tests whether there are real initial strains on the element - may be redefined in specialized derived classes
	virtual bool HasInitialStrains() { return false; }
	// the initial strains may be interpolated over the nodes in a derived class
	virtual void GetInitialStrainComponentsAtPoint(Matrix2D & isc, const Vector2D & ploc, bool flagDraw) {}

	// Elastic material parameters
	double ElasticCoefficient();		// lambda/(lambda+2mu) == nu/(1-nu)
	double ElasticMu();

	// implementation of IntegrationRuleProvider
	virtual void DefineIntegrationRule(IntegrationRule & integrationRule);
};