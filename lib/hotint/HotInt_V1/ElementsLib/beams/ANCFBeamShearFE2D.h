//#**************************************************************
//# filename:             ANCFBeamShearFE2D.h
//#
//# author:               Gerstmayr, Vetyukov
//#
//# generated:						
//# description:          
//# comments:
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
#include "FiniteElementGenericBeam2D.h"

const int ANCFBeamShearFE2DGenericmaxDOF = 12;
const int ANCFBeamShearFE2DGenericmaxIP = 10;

class ANCFBeamShearFE2DGeneric : public FiniteElementGenericBeam2D //$EDC$[beginclass,classname=ANCFBeamShearFE2DGeneric,parentclassname=FiniteElementGenericBeam2D]
{
public:
	ANCFBeamShearFE2DGeneric(MBS* mbs) : FiniteElementGenericBeam2D(mbs), use_reduced_integration(0), use_reduced_integration_poisson(0),
		use_contmech_formulation(0), shear_correction_factor(1)
	{
	};
	ANCFBeamShearFE2DGeneric(const ANCFBeamShearFE2DGeneric& e) : FiniteElementGenericBeam2D(e.mbs) 
	{ 

		CopyFrom(e); 
	}
	~ANCFBeamShearFE2DGeneric(void) { ; }

	// sets material number and checks if for inelastic material the Poisson ratio is zero (necessary for theory so far)
	virtual void SetMaterialNum(int mnum);

	virtual void CopyFrom(const Element &e)
	{
		FiniteElementGenericBeam2D::CopyFrom(e);
		const ANCFBeamShearFE2DGeneric& beam = dynamic_cast<const ANCFBeamShearFE2DGeneric&> (e);
		this->use_reduced_integration = beam.use_reduced_integration;
		this->use_reduced_integration_poisson = beam.use_reduced_integration_poisson;
		this->use_contmech_formulation = beam.use_contmech_formulation;
		this->shear_correction_factor = beam.shear_correction_factor;
	};

	virtual const char* GetElementSpec() const { return "ANCFBeamShearFE2DGeneric"; }
	virtual TFiniteElementType GetElementType() const { return TFE_Beam2D; }
	virtual int GetActualInterpolationOrder() const { return NNodes()*2; }

	int UseReducedIntegration() const {return use_reduced_integration; }
	void SetReducedIntegration(int use_red_int) {use_reduced_integration = use_red_int; }
	int UseReducedIntegrationPoisson() const { return use_reduced_integration_poisson; }
	void SetReducedIntegrationPoisson(int use_red_int_p) {use_reduced_integration_poisson = use_red_int_p; }
	int UseContinuumMechanicsFormulation() const {return use_contmech_formulation; }
	void SetUseContinuumMechanicsFormulation(int use_contm_f) { use_contmech_formulation = use_contm_f; }
	double ShearCorrectionFactor() const { return shear_correction_factor; }
	void SetShearCorrectionFactor(double ks) { shear_correction_factor = ks; }
	virtual void DefineIntegrationRule(IntegrationRule& integrationRule, int use_reduced_integration_y);
	virtual void DefineIntegrationRule(IntegrationRule& integrationRule) { DefineIntegrationRule(integrationRule, 0); }
	virtual void DefineIntegrationRulesStructuralMech(IntegrationRule& ruleGamma, IntegrationRule& ruleTheta, IntegrationRule& ruleThickness) = 0;

	virtual void GetElasticityMatrix(Matrix3D& Dm, int use_red_integration=0, int poisson_case=0);
	virtual void EvalF2(Vector& f, double t);
	virtual void EvalF2StructuralMech(Vector& f, double t);
	virtual void EvalF2ContMech(Vector& f, double t);
	virtual void EvalM(Matrix& m, double t);

	virtual void GetH(Matrix& H);
	virtual void GetIntDuDq(Matrix& dudq) {	GetH(dudq);	}

	virtual Vector2D GetRefPosy2D(const Vector2D& p_loc) const;
	virtual Vector2D GetPosy2D(const Vector2D& p_loc, const XGProvider& xg) const;
	virtual Vector2D GetPosx2D(const Vector2D& p_loc, const XGProvider& xg) const;
	virtual Vector2D GetVely2D(const Vector2D& p_loc, const XGProvider& xg) const;

	virtual Vector2D GetRefPos2D(const Vector2D& p_loc) const;
	virtual Vector2D GetDisplacement2D(const Vector2D& p_loc) const;
	virtual Vector2D GetDisplacement2DD(const Vector2D& p_loc) const;
	virtual Vector2D GetVel2D(const Vector2D& p_loc) const;
	virtual Vector2D GetVel2DD(const Vector2D& p_loc) const;
	virtual Vector2D GetPos2D(const Vector2D& p_loc)	const;
	virtual Vector2D GetPos2DD(const Vector2D& p_loc)	const;
	virtual Vector2D GetPos2DD(const Vector2D& p_loc, double def_scale) const;


	virtual Vector2D GetDisplacement2D(const Vector2D& p_loc, const XGProvider& xg) const;
	virtual Vector2D GetVel2D(const Vector2D& p_loc, const XGProvider& xg) const;
	virtual Vector2D GetPos2D(const Vector2D& p_loc, const XGProvider& xg)	const;
	virtual Vector2D GetPos2D(const Vector2D& p_loc, double def_scale, const XGProvider& xg) const;



	virtual void GetdPosdqT(const Vector2D& p_loc, Matrix& dpdqi);
// get jacobian [du/dx, du/dy] 
// continuum mechanics formulation: for u corresponding to dof-vector xg0, x, y in reference coordinates, works for pre-curved elements
// reissner formulation: simple scaling ((lx/2, 0), (0, ly/2)), in this formulation pre-curvature is not possible up to now anyway
	virtual void GetJacobi(Matrix3D& jac, const Vector2D& ploc, const Vector& xg0) const;

	// number of shape functions - 4 for linear element, 6 for quadratic element
	virtual int NS() const { return 2*NNodes(); };
	// always four dofs for each node - posx, posy, slopex, slopey
	virtual int DOFPerNode() const { return 4; }
	virtual Vector3D GetDOFDirD(int idof) const;
	virtual Vector3D GetDOFPosD(int idof) const;

	virtual Vector2D GetInplaneUnitVector2D(const double& p_loc) const;
	virtual Vector2D GetInplaneUnitVector2DD(const double& p_loc) const;
	virtual Vector2D GetRefInplaneUnitVector2D(const double& p_loc) const;
	virtual Vector2D GetInplaneUnitVectorP2D(const double& p_loc) const;

#pragma region shapefunction routines
	virtual double GetS0(const Vector2D& ploc, int shape) const = 0;
	virtual double GetS0x(const Vector2D& ploc, int shape) const = 0;
	virtual double GetS0xx(const Vector2D& ploc, int shape) const = 0;
	virtual double GetS0y(const Vector2D& ploc, int shape) const = 0;
	virtual double GetS0xy(const Vector2D& ploc, int shape) const = 0;
#pragma endregion

	virtual double GetEpsAxial(const double& p_loc) const { return GetGamma12D(p_loc, xgCompute); }
	virtual double GetEpsAxialD(const double& p_loc) const { return GetGamma12D(p_loc, xgDraw); }
	virtual double GetKappa(const double& p_loc) const{ return GetThetax2D(p_loc, xgCompute); }
	virtual double GetKappaD(const double& p_loc) const { return GetThetax2D(p_loc, xgDraw); }

	// compute 2D strain tensor in Matrix2D format
	virtual void ComputeStrainMatrix2D(Vector2D ploc, Matrix2D& strain, const XGProvider& xg) const;
	
	virtual void GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables);
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector2D& p_loc, bool flagD);

	double GetThetax2D(double x0, const XGProvider& xg) const;
	void GetDeltaThetax2D(double x0, double& kappa, Vector& DeltaThetax, const XGProvider& xg) const;
	//Gamma1, Gamma2 (and corresponding delta)
	double GetGamma12D(double x0, const XGProvider& xg) const;
	double GetGamma22D(double x0, const XGProvider& xg) const;
	void GetDeltaGamma1Gamma22D(double x0, double& Gamma1, double& Gamma2, Vector& DeltaGamma1, Vector& DeltaGamma2, const XGProvider& xg) const;
	void GetDeltaEyy(double x0, double& Eyy, Vector& DeltaEyy, const XGProvider& xg) const;
	double GetEyy(double x0, const XGProvider& xg) const;
	// Stretchy ... sqrt(r_y^T r_y)
	void GetDeltaStretchy(double x0, double& stretchy, Vector& DeltaStretchy, const XGProvider& xg) const;
	double GetStretchy(double x0, const XGProvider& xg) const;

	// plastic material law is essentially two-dimensional
	virtual int DimensionOfPlasticMaterialLaw() const { return 2;};
	virtual double PostNewtonStep(double t);

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		
	virtual int SetElementDataAuto(ElementDataContainer& edc);
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata);
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); 

protected:
	// EDC access for q0 Vector
	//EDC Vector q0;		//$EDC$[varaccess,readonly,EDCvarname="node1_initial_position",EDCfolder="Initialization",vecstart=1,vecend=2,tooltiptext="initial values for position of node 1 (left). r1 = [x1,y1]"]
	//EDC Vector q0;		//$EDC$[varaccess,readonly,EDCvarname="node1_initial_slope_eta",EDCfolder="Initialization",vecstart=3,vecend=4,tooltiptext="initial values for slope vector of node 1 (left). [d(r1)/d(eta)]"]
	//EDC Vector q0;		//$EDC$[varaccess,readonly,EDCvarname="node2_initial_position",EDCfolder="Initialization",vecstart=5,vecend=6,tooltiptext="initial values for position of node 2 (right). r2 = [x2,y2]"]
	//EDC Vector q0;		//$EDC$[varaccess,readonly,EDCvarname="node2_initial_slope_eta",EDCfolder="Initialization",vecstart=7,vecend=8,tooltiptext="initial values for slope vector of node 2 (right). [d(r2)/d(eta)]"]

	// EDC access for size Vector - lx, ly, lz
	//EDC double size(1); //$EDC$[varaccess,EDCvarname="lx",EDCfolder="Geometry",tooltiptext="length of undeformed beam element l_x"]
	//EDC double size(2); //$EDC$[varaccess,EDCvarname="ly",EDCfolder="Geometry",tooltiptext="length of undeformed beam element l_y"]
	//EDC double size(3); //$EDC$[varaccess,EDCvarname="lz",EDCfolder="Geometry",tooltiptext="length of undeformed beam element l_z"]

	// use reduced integration along axis
	int use_reduced_integration; //$EDC$[varaccess,EDCvarname="use_reduced_integration",EDCfolder="Initialization",tooltiptext="use reduced integration along the beam axis"]
	// use reduced integration in thickness direction for Poisson locking in continuum mechanics formulation
	int use_reduced_integration_poisson;  //$EDC$[varaccess,EDCvarname="use_reduced_integration_poisson",EDCfolder="Initialization",tooltiptext="use reduced integration across beam thickness to avoid Poisson locking"]
	// use continuum mechanics formulation based on stress-strain relationship
	int use_contmech_formulation;  //$EDC$[varaccess,EDCvarname="use_continuum_mechanics",EDCfolder="Initialization",tooltiptext="use continuum mechanics formulation [1] or structural mechanics formulation [0]"]

	double shear_correction_factor;  //$EDC$[varaccess,EDCvarname="shear_correction_factor",EDCfolder="Initialization",tooltiptext="shear correction factor"]
};//$EDC$[endclass,ANCFBeamShearFE2DGeneric]


class ANCFBeamShearFE2DLinear : public ANCFBeamShearFE2DGeneric //$EDC$[beginclass,classname=ANCFBeamShearFE2DLinear,parentclassname=ANCFBeamShearFE2DGeneric,addelementtype=TAEBody+TAENotInRelease,addelementtypename=ANCFBeamShearFE2DLinear]
{
public:
	ANCFBeamShearFE2DLinear(MBS* mbs) : ANCFBeamShearFE2DGeneric(mbs)
	{
		nodes.SetLen(NNodes());
		nodes.SetAll(0);
		q0.SetLen(2*SOS());
	};

	ANCFBeamShearFE2DLinear(const ANCFBeamShearFE2DLinear& e) : ANCFBeamShearFE2DGeneric(e.mbs) 
	{ 
		CopyFrom(e); 
	}
	~ANCFBeamShearFE2DLinear(void) { ; }

	virtual Element* GetCopy()
	{
		Element* e = new ANCFBeamShearFE2DLinear(*this);
		return e;
	};
	virtual void CopyFrom(const Element &e)
	{
		ANCFBeamShearFE2DGeneric::CopyFrom(e);
	};

	void SetANCFBeamShearFE2DLinear(int bodyind,
		const Vector& xc1, 
		const Vector& xc2, 
		int n1, 
		int n2, 
		int material_num,
		const Vector3D& size,
		const Vector3D& color);

	virtual void SetANCFBeamShearFE2DLinear(int n1, int n2, 
																					int material_num,	const Vector3D& asize,	const Vector3D& color);

	virtual const char* GetElementSpec() const { return "ANCFBeamShearFE2DLinear"; }
	//number of nodes:
	virtual int NNodes() const { return 2; }

#pragma region shapefunction routines
	virtual double GetS0(const Vector2D& ploc, int shape) const;
	virtual double GetS0x(const Vector2D& ploc, int shape) const;
	virtual double GetS0xx(const Vector2D& ploc, int shape) const;
	virtual double GetS0y(const Vector2D& ploc, int shape) const;
	virtual double GetS0xy(const Vector2D& ploc, int shape) const;
#pragma endregion

	virtual void DefineIntegrationRulesStructuralMech(IntegrationRule& ruleGamma, IntegrationRule& ruleTheta, IntegrationRule& ruleThickness);

	virtual int SetElementData(ElementDataContainer&edc);
	virtual void GetElementData(ElementDataContainer &edc);
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		
	virtual int SetElementDataAuto(ElementDataContainer& edc);
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata);
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); 

};//$EDC$[endclass,ANCFBeamShearFE2DLinear]

class ANCFBeamShearFE2DQuadratic : public ANCFBeamShearFE2DGeneric //$EDC$[beginclass,classname=ANCFBeamShearFE2DQuadratic,parentclassname=ANCFBeamShearFE2DGeneric,addelementtype=TAEBody+TAENotInRelease,addelementtypename=ANCFBeamShearFE2DQuadratic]
{
public:
	ANCFBeamShearFE2DQuadratic(MBS* mbs) : ANCFBeamShearFE2DGeneric(mbs)
	{
		nodes.SetLen(NNodes());
		nodes.SetAll(0);
		q0.SetLen(2*SOS());

	};

	ANCFBeamShearFE2DQuadratic(const ANCFBeamShearFE2DQuadratic& e) : ANCFBeamShearFE2DGeneric(e.mbs) 
	{ 
		CopyFrom(e); 
	}
	~ANCFBeamShearFE2DQuadratic(void) { ; }

	virtual Element* GetCopy()
	{
		Element* e = new ANCFBeamShearFE2DQuadratic(*this);
		return e;
	};
	virtual void CopyFrom(const Element &e)
	{
		ANCFBeamShearFE2DGeneric::CopyFrom(e);
	};

	void SetANCFBeamShearFE2DQuadratic(int bodyind,
		const Vector& xc1, 
		const Vector& xc2, 
		const Vector& xc3,
		int n1, 
		int n2, 
		int n3,
		int material_num,
		const Vector3D& size,
		const Vector3D& color);

	virtual void SetANCFBeamShearFE2DQuadratic(int n1, int n2, int n3,
																						 int material_num,	const Vector3D& asize,	const Vector3D& color);

	virtual const char* GetElementSpec() const { return "ANCFBeamShearFE2DQuadratic"; }
	//number of nodes:
	virtual int NNodes() const { return 3; }

	#pragma region shapefunction routines
	virtual double GetS0(const Vector2D& ploc, int shape) const;
	virtual double GetS0x(const Vector2D& ploc, int shape) const;
	virtual double GetS0xx(const Vector2D& ploc, int shape) const;
	virtual double GetS0y(const Vector2D& ploc, int shape) const;
	virtual double GetS0xy(const Vector2D& ploc, int shape) const;
#pragma endregion

	virtual void DefineIntegrationRulesStructuralMech(IntegrationRule& ruleGamma, IntegrationRule& ruleTheta, IntegrationRule& ruleThickness);

	virtual int SetElementData(ElementDataContainer&edc);
	virtual void GetElementData(ElementDataContainer &edc);
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		
	virtual int SetElementDataAuto(ElementDataContainer& edc);
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata);
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); 

	// EDC access for values corresponding to third node
	//EDC Vector q0;		//$EDC$[varaccess,readonly,EDCvarname="node3_initial_position",EDCfolder="Initialization",vecstart=9,vecend=10,tooltiptext="initial values for position of node 3 (middle). r3 = [x3,y3]"]
	//EDC Vector q0;		//$EDC$[varaccess,readonly,EDCvarname="node3_initial_slope_eta",EDCfolder="Initialization",vecstart=11,vecend=12,tooltiptext="initial values for slope vector of node 3 (right). [d(r3)/d(eta)]"]

};//$EDC$[endclass,ANCFBeamShearFE2DQuadratic]
