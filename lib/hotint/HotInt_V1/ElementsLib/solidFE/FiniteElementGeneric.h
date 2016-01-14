//#**************************************************************
//#
//# filename:             FiniteElementGeneric.h
//#
//# author:               Gerstmayr Johannes, Aigner Larissa, YV
//#
//# generated:						October 2010
//# description:          common finite element base class - for 2D & 3D finite elements
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

#include "IntegrationRule.h"

// these numbers are needed as a constant -
// - it will affect the creation of arrays with fixed size on the stack (ConstVector, ConstMatrix)
const int FEmaxDOF = 81;	//27 nodes hexahedral
const int FFRFsize = 7;		// additional variables for FFRF formulation
const int FE2DmaxDOF = 18; // 3*3 nodes for second order quadriangle, 2 DOFs for each

template<class BodyXD>
class FiniteElementGeneric :
	public BodyXD,		// it is either a 3D or a 2D body
	public IntegrationRule::IntegrationRuleProvider		// each finite element must be able to provide its own integration rules
{
public:
	FiniteElementGeneric(MBS* mbsi) : BodyXD(mbsi), massmatrix(), Hmatrix(), stiffnessmatrix()
	{
		InitConstructor();
	}
	/// Copy constructor
	FiniteElementGeneric(const FiniteElementGeneric& e) : BodyXD(e.mbs), massmatrix(), Hmatrix(), stiffnessmatrix()
	{
		InitConstructor();
		CopyFrom(e);
	};

	virtual int UseSparseM() const { return 1; }

	/// Filling in the data from an existing element
	virtual void CopyFrom(const Element& e);

	/// For a given mesh topology, computes the mapping between the local and the global degrees of freedom
	virtual void LinkToElements();

	// integration rules and ds matrices are initialized here
	virtual void PreAssemble();

	///@name Element data container
	//@{
	/// fill in all element data
	virtual void GetElementData(ElementDataContainer& edc);
	
	/// set element data according to ElementDataContainer
	virtual int SetElementData(ElementDataContainer& edc);
	//@}

	// assigns the material
	virtual void SetMaterialNum(int mnum);

	///@name Identification of the element type
	//@{
	/// Text form id
	virtual const char* GetElementSpec() const = 0;
	/// Type of this particular finite element
	virtual TFiniteElementType GetElementType() const = 0;
	//@}

	/// Get/set the actual status of geometric nonlinearityFlag to set if computation of element is geometrically nonlinear
	GeometricNonlinearityStatus GetGeometricNonlinearityStatus() const { return geometricNonlinearityStatus; }
	void SetGeometricNonlinearityStatus(GeometricNonlinearityStatus gns) { geometricNonlinearityStatus = gns; }

	///@name Basic properties of the finite element
	//@{
	// set to (1-cms)*(FlexDof+ffrfdim)
	virtual int SOS() const { return FlexDOF(); } ///<Size of K and M (stiffness matrix and mass matrix)
	virtual int SOSowned() const { return 0; } ///< Number of internal degrees of freedom
	virtual int ES() const   {return 0; }  ///< Size of first order explicit equations
	virtual int IS() const  { return 0; }  ///< Implicit (algebraic) size

	virtual int IsRigid() const { return 0; }	///< The element is not rigid
	virtual int IsFiniteElement() const { return 1; } ///< It is a finite element
	//virtual void SetConcentratedMass(double cm, int i) {concentratedmass[i-1] = cm;}
	//@}

	/// Returns an integration rule for the given integrated value or NULL if this integrated value is not supported
	const IntegrationRule * GetIntegrationRule(IntegrationRule::IntegratedValueType ivt) const;

	//change for other shape functions and dimension:
	// spatial dimension of the element
	virtual int Dim() const = 0;
	// number of degrees of freedom in a node
	virtual int DOFPerNode() const { return Dim(); }
	//number of shape functions:
	virtual int NS() const { return nodes.Length(); }
	//number of flexible degrees of freedom:
	virtual int FlexDOF() const { return DOFPerNode() * NNodes(); }		//*YV: used to be Dim()*NS()
	//number of nodes:
	virtual int NNodes() const { return nodes.Length(); }
	//size of K and M for not reduced system, number of total degrees of freedom
	virtual int XGLength() const { return FlexDOF(); }

	// ACCESS TO GENERALIZED COORDINATES of the element itself and of the reference frame
	virtual const double& XGP(int iloc) const {return GetXact(ltg.Get(iloc+XGLength()));}
	virtual const double& XGPD(int iloc) const {return GetDrawValue(ltg.Get(iloc+XGLength()));}
	// ACCESS TO GENERALIZED COORDINATES of the element itself and (possibly) of the reference frame
	virtual const double& GetXact(int i) const { return mbs->GetXact(i); }
	virtual const double& GetDrawValue(int iloc) const { return mbs->GetDrawValue(iloc); }

	/// Estimate bandwidth for sparse matrices of assembled matrices.
	/// Improved for sparse matrix computation, *2 for overlapping.
	virtual int ElementBandwidth() const { return SOS()*2; }


	////return box in which the element fits
	//virtual Box3D GetElementBox() const
	//{
	//	Box3D b;
	//	for (int i=1;i <= NNodes(); i++)
	//	{
	//		b.Add(GetNodePos(i));
	//	}
	//	return b;
	//}

	//virtual Box3D GetElementBoxD() const
	//{
	//	Box3D b;
	//	for (int i=1;i <= NNodes(); i++)
	//	{
	//		b.Add(GetNodePosD(i));
	//	}
	//	return b;
	//}




	// this element is a part of a floating frame of reference formulation model
	virtual int IsFFRF() const { return 0; }
	// this element is a part of a component mode synthesis model
	virtual int IsCMS() const { return 0; }

	virtual const Element& GetElement(int globind) const 
	{
		return GetMBS()->GetElement(globind);
	}
	virtual Element* GetElementPtr(int globind)
	{
		return GetMBS()->GetElementPtr(globind);
	}

	/// Compute DS matrix and apply element jacobian inverse for the point with given local coordinates;\n
	/// also compute determinate of jacobian;\n
	/// for simple spatial elements jacinvDS has the dimensions Dim() x NS()
	///@return Determinant of the jacobian
	// this function explicitly operates with Vector3D as we do not parametrize with respect to VectorXD;
	// in 2D case it will be enough to convert this argument to const Vector2D &.
	virtual double GetJacInvDS(const Vector3D& ploc, Matrix& jacinvDS) const = 0;

	/// Precompute DS Matrices and apply element transformation matrix
	/// in all integration points of the stiffness matrix.
	virtual void BuildDSMatrices();
	/// Precomputed jacobian at a Gauss point with the given number
	virtual double GetPrecomputedJacDet(int ip) const
	{
		return integrationPointStiffnessMatrixLocalData(ip)->jacdet;
	}

	/// Copies global element coordinates into some local "compute" cache
	void SetComputeCoordinates(Vector & xgCache);

	/// Provides a (possibly) pre-computed matrix of gradients of shape functions
	/// at i-th integration point for the stiffness matrix
	void GetGrad(int ip, Matrix & grad) const
	{
		IntegrationPointsIterator it(integrationRuleStiffness);
		it.SetIndex(ip);
		GetGrad(it, grad);
	}

	// access to the node indices
	virtual const int& NodeNum(int i) const { return nodes(i); }
	virtual int& NodeNum(int i) { return nodes(i); }
	virtual const Node& GetNode(int i) const 
	{
		return GetMBS()->GetNode(NodeNum(i));
	}
	virtual Node& GetNode(int i) 
	{
		return GetMBS()->GetNode(NodeNum(i));
	}

	virtual void GetNeighbourElements(TArray<int>& neighbours);		///< get a list of elements which are neighbours of this element
	virtual void SetFFRFElement(int CMSelementI) { assert(0 && "Non-FFRF element was constructed with CMSelement index."); }
	virtual int AddBodyNode(Node & n) { return mbs->AddBodyNode(&n); }

	virtual int DataS() const
	{
		return integrationRuleStiffness->GetNumberOfIntegrationPoints() * GetMaterial().GetInelasticVariablesCount();
	}

	// inelasticity
	// {
	virtual int IsInelasticMaterial() const {return GetMaterial().IsInelasticMaterial();}

	// Vector notation
	//( eps(1,1), eps(2,2), eps(3,3), 2*eps(3,2), 2*eps(3,1), 2*eps(1,2), param_1, param_2, ... )
	virtual void InelasticVariablesToData(const Vector& inelastic_variables_IP, int nip)
	{
		int n = GetMaterial().GetInelasticVariablesCount();
		assert(inelastic_variables_IP.Length() == n);
		for (int i=1; i<=n; i++)
		{
			XData((nip-1)*n+i) = inelastic_variables_IP(i);
		}
	}
	virtual void InelasticVariablesToDataD(const Vector& inelastic_variables_IP, int nip)
	{
		int n = GetMaterial().GetInelasticVariablesCount();
		assert(inelastic_variables_IP.Length() == n);
		for (int i=1; i<=n; i++)
		{
			XDataD((nip-1)*n+i) = inelastic_variables_IP(i);
		}
	}
	virtual void DataToInelasticVariables(Vector& inelastic_variables_IP, int nip)
	{
		int n = GetMaterial().GetInelasticVariablesCount();
		inelastic_variables_IP.SetLen(n);
		for (int i=1; i<=n; i++)
		{
			inelastic_variables_IP(i) = XData((nip-1)*n+i);
		}
	}

	virtual void DataToInelasticVariablesD(Vector& inelastic_variables_IP, int nip)
	{
		int n = GetMaterial().GetInelasticVariablesCount();
		inelastic_variables_IP.SetLen(n);
		for (int i=1; i<=n; i++)
		{
			inelastic_variables_IP(i) = XDataD((nip-1)*n+i);
		}
	}

	virtual void DataLastStepToInelasticVariables(Vector& inelastic_variables_IP, int nip) //read inelastic strain from last computed step
	{
		int n=GetMaterial().GetInelasticVariablesCount();
		inelastic_variables_IP.SetLen(n);
		for (int i=1; i<=n; i++)
		{
			inelastic_variables_IP(i) = GetMBS()->GetLastDataVector()(LTGdata((nip-1)*n + i));
		}
	}

	virtual void DataLastNLitToInelasticVariables(Vector& inelastic_variables_IP, int nip) //read inelastic strain from last computed nonlinear iteration
	{
		int n=GetMaterial().GetInelasticVariablesCount();
		inelastic_variables_IP.SetLen(n);
		for (int i=1; i<=n; i++)
		{
			inelastic_variables_IP(i) = GetMBS()->GetLastNLItDataVector()(LTGdata((nip-1)*n + i));
		}
	}

	// resets the accumulated inelastic state
	virtual void EraseInelasticVariables();
	// }

	virtual void AddMSparse(SparseMatrix& m, double t);  //add sparse mass matrix into full system matrix

	virtual Vector3D GetAngularMomentum(const Vector3D& p_ref) const;

	//$ YV 2012-06: an element needs to know local coordinates of its nodes
	// to be able to compute field variables there; the computation is to be implemented in derived finite element classes
	virtual Vector3D GetNodeLocPos(int local_node_number) const { assert(0); return 0; }
	
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & local_position, bool flagD) {assert(0); 	return 0; }	//$ DR 2012-09-12 bugfix: is to be implemented in derived finite element classes
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, int local_node_number, bool flagD)
	{
		assert(!IsCMS());		// for CMS elements the function needs to be redefined
		return GetFieldVariableValue(fvd, GetNodeLocPos(local_node_number), flagD);	//$ DR 2012-09-12 bugfix: "return GetField..." instead of "return Element::GetField..."
	}

	// $EK 2013-02-12 computes volume of the finite element (for FiniteElement3D and FiniteElement2D)
	virtual double ComputeVolume() { return 0.; }
	// $EK 2013-03-21 compute center of mass via integration (for FiniteElement3D)
	virtual void ComputeCenterOfMass(Vector3D& center) { center.Set(0,0,0); }
protected:
	/// General initialization of internal variables
	void InitConstructor()
	{
		massmatrix.Init();
		stiffnessmatrix.Init();
		Hmatrix.Init();
		SetGeometricNonlinearityStatus(GNS_NonlinearSmallStrain);
	}

	// common part of the set-functions
	void SetFiniteElementGeneric(int bodyindi, const TArray<int>& nodelist, int material_num, const Vector3D& coli, int CMSelementi = 0);

	TArray<int> nodes;

	int bodyind;

	//temporary storage for mass matrix, damping, etc:
	Matrix massmatrix;
	Matrix Hmatrix;
	Matrix stiffnessmatrix;
	// information for each integration point of the stiffness matrix
	struct IntegrationPointStiffnessMatrixLocalData
	{
		// when the flag "store_FE_matrices" is switched on,
		// then we pre-compute and store here the matrix of components of the
		// gradients of shape functions with respect to the spatial coordinates
		// in the reference configuration
		Matrix * grad;
		// determinant of the jacobian of the transformation matrix between the
		// local (element) and spatial coordinates
		double jacdet;
		//$ YV 2011-09-22
		// the initial strains were removed from the integration points
		IntegrationPointStiffnessMatrixLocalData() : grad(NULL) {}
		~IntegrationPointStiffnessMatrixLocalData()
		{
			if(grad != NULL)
				delete grad;
		}
		// copy operator is needed because of the structure, in which these objects are stored
		IntegrationPointStiffnessMatrixLocalData & operator=(const IntegrationPointStiffnessMatrixLocalData & data)
		{
			jacdet = data.jacdet;
			if(data.grad == NULL)
				grad = NULL;
			else
				grad = new Matrix(*data.grad);
			return *this;
		}
	};
	// either provides a pre-computed gradient, or computes it "on the fly"
	void GetGrad(IntegrationPointsIterator & integrationPointStiffnessIterator, Matrix & grad) const;

	/// integration rules for stiffness K, mass M and H-matrices
	/// (H-matrix - simple integral over shape functions, needed to apply distributed forces)
	IntegrationRule * integrationRuleStiffness;		///< K
	TArray<IntegrationPointStiffnessMatrixLocalData*> integrationPointStiffnessMatrixLocalData;
	friend class IntPointsStiffnessIterator;
	
	//$ EK 2012-09-07, renamed AD 2013-09-30: include extended/custom elements
	friend class IntPointsStiffnessIterator_Quadrilateral_custom;
	// $EK 2012-12-13 now Hexaherdal_extended is inherited from Hexahedral and FiniteElement3D_extended
	friend class IntPointsStiffnessIterator_FiniteElement3D_custom;

	IntegrationRule * integrationRuleMass;				///< M
	IntegrationRule * integrationRuleLoad;				///< H

	/// maximal power of the shape functions of the present element;
	/// determines the choice of the integration rules
	virtual int GetActualInterpolationOrder() const = 0;
	/// current geometric nonlinearity status, which affects different aspects of the behavior of the element
	GeometricNonlinearityStatus geometricNonlinearityStatus;

	//$ YV 2012-04-10: a derived class may wish to control the region, in which the finit element is plotted, then it may implement this function to modify the array of points
	virtual void FitPointsInDrawingRegion(TArray<Vector3D> & points, int n1, int n2) {}

		// $EK 20130408 added in order to be consistent to constraints
	virtual TKinematicsAccessFunctions GetKinematicsAccessFunctions(int mode = 1) const
	{
		return BodyXD::GetKinematicsAccessFunctions(mode);
	}
public:
	~FiniteElementGeneric()
	{
		for(int i = 1; i <= integrationPointStiffnessMatrixLocalData.Length(); i++)
			delete integrationPointStiffnessMatrixLocalData(i);
	}
};
