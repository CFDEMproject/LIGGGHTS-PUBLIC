//#**************************************************************
//# filename:             FiniteElementGenericBeam2D.h
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
#include "FiniteElementGeneric.h"
#include "XGProvider.h"


class FiniteElementGenericBeam2D : public FiniteElementGeneric<Body2D> //$EDC$[beginclass,classname=FiniteElementGenericBeam2D]
{
public:
enum {MAX_PLAST_GRIDP_Y = 50}; // maximal number of plastic grid points in y direction
enum {MAX_PLASTIC_QUANTITIES = 8}; // maximal number of plastic quantities
enum {MAX_IP = 10}; // maximum number of integration points in one direction
enum {MAX_NS = 6}; // maximum number of shape functions
enum {DIM = 2};

	FiniteElementGenericBeam2D(MBS* mbs) : FiniteElementGeneric<Body2D>(mbs), plasticip_x(0), plasticip_y(0) 
	{
		this->geometricNonlinearityStatus = GNS_NonlinearLargeStrain;

		xgReferenceState.SetLen(FlexDOF());
		xgInit.SetXGProvider(this, &xgReferenceState);
		xgCompute.SetXGProvider(this, &xgReferenceState);
		xgDraw.SetXGProvider(this, &xgReferenceState);
	}	// default constructor
	FiniteElementGenericBeam2D(const FiniteElementGenericBeam2D& e) : FiniteElementGeneric<Body2D>(e.mbs) 
	{	
		xgReferenceState.SetLen(FlexDOF());
		xgInit.SetXGProvider(this, &xgReferenceState);
		xgCompute.SetXGProvider(this, &xgReferenceState);
		xgDraw.SetXGProvider(this, &xgReferenceState);

		CopyFrom(e);	
	}

	void SetFiniteElementGenericBeam2D(
		int bodyind,
		const Vector& xc1, 
		const Vector& xc2, 
		int n1, 
		int n2, 
		int material_num,
		const Vector3D& size,
		const Vector3D& color);

	virtual void CopyFrom(const Element& e);
	virtual void PreAssemble();

	virtual mystr GetElementSpecification() //$EDC$[funcaccess,EDCvarname="element_type",tooltiptext="specification of element type. Once the element is added to the mbs, you MUST NOT change this type anymore!"]
	{
		mystr type = GetElementSpec();	//$ DR this is necessary for the skript language
		return type;
	}	

	virtual int Dim() const { return 2; }
	virtual Vector3D GetDOFDirD(int idof) const = 0;

	virtual int UseSparseM() const { return 0; }

	virtual double GetJacInvDS(const Vector3D& p_loc, Matrix& jacinvDS) const { assert(0); return 0; }

	virtual double GetBeamEIy() const {return GetMaterial().BeamEIy(); }
	virtual double GetBeamEA() const {return GetMaterial().BeamEA(); }
	virtual double GetBeamGAky() const {return GetMaterial().BeamGAky(); }
	virtual double GetBeamRhoA() const {return GetMaterial().BeamRhoA(); }
	virtual double GetBeamRhoIz() const {return GetMaterial().BeamRhoIx(); }

	virtual double& GetLx() { return size.X(); }
	virtual const double& GetLx() const { return size.X(); }
	virtual double& GetLy() { return size.Y(); }
	virtual const double& GetLy() const { return size.Y(); }
	virtual double& GetLz() { return size.Z(); }
	virtual const double& GetLz() const { return size.Z(); }

	// AH: keep it simple, stupid
	// current pos. vectors + deriv.
	// p_loc in (-1, 1) always
	virtual Vector2D GetPos2D(const double& p_loc)		const	{	return GetRefPos2D(p_loc)	+	GetDisplacement2D(p_loc);	}
	virtual Vector2D GetPos2DD(const double& p_loc)		const	{	return GetRefPos2D(p_loc)	+	GetDisplacement2DD(p_loc);	}
	virtual Vector2D GetPosx2D(const double& p_loc)		const	{	return GetRefPosx2D(p_loc)	+	GetDisplacementx2D(p_loc);	}
	virtual Vector2D GetPosx2DD(const double& p_loc)	const	{	return GetRefPosx2D(p_loc)	+	GetDisplacementx2DD(p_loc);	}
	virtual Vector2D GetPosxx2D(const double& p_loc)	const	{	return GetRefPosxx2D(p_loc)	+	GetDisplacementxx2D(p_loc);	}
	virtual Vector2D GetPosxx2DD(const double& p_loc)	const	{	return GetRefPosxx2D(p_loc)	+	GetDisplacementxx2DD(p_loc);	}

	virtual Vector2D GetPos2D(const Vector2D& p_loc)	const	{	return GetPos2D(p_loc.X())	+	0.5*GetLy()*p_loc.Y() * GetInplaneUnitVector2D(p_loc.X());	}
	virtual Vector2D GetPos2DD(const Vector2D& p_loc)	const	{	return GetPos2DD(p_loc.X())	+	0.5*GetLy()*p_loc.Y() * GetInplaneUnitVector2DD(p_loc.X());	}
	virtual Vector2D GetPos2DD(const Vector2D& p_loc, double def_scale) const {	return GetRefPos2D(p_loc) + def_scale * GetDisplacement2DD(p_loc); }

	// AH: this method is called by mbsload drawing routine, implementation in rigid2d returns zero in case of displacement dofs
	virtual Vector3D GetRefPosD() const 
	{
		Vector2D p = GetPos2DD(0.);
		return Vector3D(p.X(), p.Y(), 0.);
	}

	// pos. vectors in reference configuration + deriv.
	virtual Vector2D GetRefPos2D(const double& p_loc) const;
	virtual Vector2D GetRefPosx2D(const double& p_loc) const;
	virtual Vector2D GetRefPosxx2D(const double& p_loc) const;

	virtual Vector2D GetRefPos2D(const Vector2D& p_loc) const;

	// displacement vectors + deriv.
	virtual Vector2D GetDisplacement2D(const double& p_loc) const;
	virtual Vector2D GetDisplacement2DD(const double& p_loc) const;
	virtual Vector2D GetDisplacementx2D(const double& p_loc) const;
	virtual Vector2D GetDisplacementx2DD(const double& p_loc) const;
	virtual Vector2D GetDisplacementxx2D(const double& p_loc) const;
	virtual Vector2D GetDisplacementxx2DD(const double& p_loc) const;

	virtual Vector2D GetDisplacement2D(const Vector2D& p_loc) const		{	return GetPos2D(p_loc)	- GetRefPos2D(p_loc);	}
	virtual Vector2D GetDisplacement2DD(const Vector2D& p_loc) const	{	return GetPos2DD(p_loc) - GetRefPos2D(p_loc);	}

	// unit vector in direction of cross-section
	virtual Vector2D GetInplaneUnitVector2D(const double& p_loc) const = 0;
	virtual Vector2D GetInplaneUnitVector2DD(const double& p_loc) const = 0;
	virtual Vector2D GetRefInplaneUnitVector2D(const double& p_loc) const = 0;

	virtual Vector2D GetInplaneUnitVectorP2D(const double& p_loc) const = 0;

	// veloctiy vectors
	virtual Vector2D GetVel2D(const double& p_loc) const;
	virtual Vector2D GetVel2D(const Vector2D& p_loc) const;
	virtual Vector2D GetVelx2D(const double& p_loc) const;

	virtual double GetEpsAxial(const double& p_loc) const = 0;
	virtual double GetEpsAxialD(const double& p_loc) const = 0;
	virtual double GetKappa(const double& p_loc) const = 0;
	virtual double GetKappaD(const double& p_loc) const = 0;

	virtual void ComputeStrainMatrix2D(Vector2D ploc, Matrix2D& strain, const XGProvider& xg) const;
			
#pragma region shape-functions
	virtual double GetS0(const Vector2D& p_loc, int shape) const = 0;
	virtual double GetS0x(const Vector2D& p_loc, int shape) const = 0;
	virtual double GetS0xx(const Vector2D& p_loc, int shape) const = 0;
#pragma endregion

	virtual void GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables);
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector2D & local_position, bool flagD);

	virtual void DrawElement();

	~FiniteElementGenericBeam2D(void);
	
#pragma region plasticity
	virtual void SetPlasticStrainGridsize(int ip_x, int ip_y) {
		plasticip_x = ip_x; plasticip_y = ip_y; 
		Vector datainit(DataS()); //automatically initialized with zeros
		SetDataInit(datainit); //initialize data variables with zero = initial inelastic strain == zero
	}
	// get plastic variables at gridpoint
	virtual void DataToPlasticVariables(Vector & plastic_variables, int ip_x, int ip_y) const;
	virtual void DataToPlasticVariablesLastTimeStep(Vector& plastic_variables, int ip_x, int ip_y) const;
	virtual void DataToPlasticVariablesD(Vector& plastic_variables, int ip_x, int ip_y) const;
	// interpolate plastic variables to ploc
	virtual void DataToPlasticVariables(Vector& plastic_variables, Vector2D& ploc);
	virtual void DataToPlasticVariablesD(Vector& plastic_variables, const Vector2D ploc);
	// interpolate plastic strain only, not hardening parameter and yield function
	virtual void DataToPlasticStrain(Vector& plasticstrain, Vector2D& ploc);
	// set plastic variables to internal XData
	virtual void PlasticVariablesToData(const Vector & plastic_variables, int ip_x, int ip_y);
	virtual void PlasticVariablesToDataD(const Vector & plastic_variables, int ip_x, int ip_y);
	// stored are maximal 5 quantities per gridpoint (like ANCFThinPlate3D):
	// in case of 1D-Plasticity: 
	// plasticstrain_xx, hardeningparameter, actual value of yield function
	// plasticstrain_yy and plasticstrain_xy are always zero
	// in case of 2D-Plasticity:
	// plasticstrain_xx, plasticstrain_yy, plasticstrain_xy,
	// hardeningparameter, actual value of yield function
	virtual int DataS() const { return NPlasticParameters() * NPlasticGridpoints(); }
	virtual int NPlasticParameters() const 
	{ 
		if (DimensionOfPlasticMaterialLaw()==1) { return 3; } 
		else if (DimensionOfPlasticMaterialLaw()==2) { return 5; } 
		else return 0;
	}
	virtual int NPlasticGridpoints() const { return plasticip_x * plasticip_y; }
	// dimension of plastic material law:
	// return 1 for 1D-Plasticity as in cable element
	// return 2 for 2D-Plasticity as in shear-deformable beam
	virtual int DimensionOfPlasticMaterialLaw() const { return 1;};

	// link XData corresponding to the plastic quantity quantitynr with Matrix mat
	// 1 .. epsp_xx
	// 2 .. epsp_yy
	// 3 .. epsp_xy
	// 4 .. hardeningparameter
	// 5 .. yield function
	virtual void GetPlasticQuantityMatrix(int quantitynr, const Matrix& mat, int flagD) const;
	virtual void GetPlasticQuantityMatrix(int quantitynr, Matrix& mat, int flagD);

	// routines for retrieving integral values of plastic strain
	// plastic axial strain = 1/h int_{-h/2}^h/2 epsp_xx dy
	virtual double GetPlasticAxialStrain(double x0, int flagD=0);
	// plastic material curvature = 12/h^3 int_{-h/2}^h/2 y epsp_xx dy
	virtual double GetPlasticKappa(double x0, int flagD=0);
	// plastic shear strain = 1/h int_{-h/2}^h/2 epsp_xy dy
	virtual double GetPlasticShearStrain(double x0, int flagD=0);
	// plastic thickness strain = 1/h int_{-h/2}^h/2 epsp_yy dy
	virtual double GetPlasticThicknessStrain(double x0, int flagD=0);

	// help routine to integrate some plastic strain over height
	// return 1/h int_{-h/2}^h/2 epsp dy
	double IntegralMean(double x0, Matrix& plasticstrain_mat);
	// help routine to evaluate some plastic strain on center line
	// return plastic strain evaluated at y=0
	double CenterLineValue(double x0, Matrix& plasticstrain_mat);

	virtual double PostNewtonStep(double t);

#pragma endregion

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		
	virtual int SetElementDataAuto(ElementDataContainer& edc);
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata);
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); 

protected:
	//double lx, ly, lz; -> use GetLx(), GetLy(), GetLz() instead
	Vector q0;	// a vector for the initial positions and velocities is needed as x_init is the initial diplacement
				// one could maybe use x_init of the node instead (see yuri's ancfthinplate3d)
	
	// number of plasticstrain-gridpoints
	int plasticip_x; //$EDC$[varaccess,EDCvarname="number_plastic_intpoints_x",EDCfolder="Geometry",tooltiptext="Number of integration points for plastic variable in x direction"]
	int plasticip_y; //$EDC$[varaccess,EDCvarname="number_plastic_intpoints_y",EDCfolder="Geometry",tooltiptext="Number of integration points for plastic variable in y direction"]

	// EDC access for elnum
	//EDC	int elnum; //$EDC$[varaccess,EDCvarname="element_number",EDCfolder="",readonly,tooltiptext="number of the element in the mbs"] //number in MBS-system

	// EDC access for name of the element
	//EDC	mystr elementname; //$EDC$[varaccess,EDCvarname="name",EDCfolder="",tooltiptext="name of the element"]

	// EDC access for nodes array
	//EDC TArray<int> nodes; //$EDC$[varaccess,EDCvarname="node_numbers",EDCfolder="Geometry",tooltiptext="Global node numbers"]

	// EDC access for mat_num
	//EDC int materialnum; //$EDC$[varaccess,EDCvarname="material_number",EDCfolder="Physics",tooltiptext="material number which contains the main material properties of the body or element"] //material number in MBS which contains material information

	//EDC Vector3D col; //$EDC$[varaccess,EDCvarname="RGB_color",EDCfolder="Graphics",tooltiptext="[red, green, blue] color of element, range = 0..1, use default color:[-1,-1,-1]"] //color for body

	// EDC access for loads
	//EDC	TArray<int> loads;  //$EDC$[varaccess,EDCvarname="loads",variable_length_vector,condition=!IsType(TConstraint),tooltiptext="Set loads attached to this element: 'nr_load1, nr_load2, ...' or empty"]

	// EDC access for sensors
	//EDC	TArray<int> sensors; //$EDC$[varaccess,EDCvarname="sensors",EDCfolder="Info",variable_length_vector,tooltiptext="attached sensors",readonly]//if the element equations of motion are dependent on sensor values (e.g. as an input), then add the sensor numbers here

	Vector xgReferenceState;
	XGProviderInit xgInit;					// initial positions and velocities
	XGProviderCompute xgCompute;		// "compute" positions and velocities
	XGProviderDraw xgDraw;					// "drawing" positions and velocities
}; //$EDC$[endclass,FiniteElementGenericBeam2D]
