//#**************************************************************
//#
//# filename:             ANCFBeamShear3D.h
//#
//# author:               Astrid und Karin
//#
//# generated:						Nov. 2010
//# description:          3D ANCF beam element with shear deformation
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
 
// Info for KN:
// lx, ly, lz are no variables of ANCFBeamShear3DGeneric anymore!
// this values are stored in size of Element

#pragma once


#pragma region ANCFBeamShear3DGeneric
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ANCFBeamShear3DGeneric ANCFBeamShear3DGeneric ANCFBeamShear3DGeneric ANCFBeamShear3DGeneric
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// base class for class ANCFBeamShear3DLinear and class ANCFBeamShear3DQuadratic

class ANCFBeamShear3DGeneric: public Body3D    //$EDC$[beginclass,classname=ANCFBeamShear3DGeneric,parentclassname=Body3D]
{
public:
	ANCFBeamShear3DGeneric(MBS* mbsi):Body3D(mbsi), massmatrix(), Hmatrix(), 
		xg(), xgd(), q0() 
	{
		perform_reduced_integration = 0;    // flag for reduced integration (formerly "RI"), used in the 3 variations of cont.mech.formulation (CMF), see BeamFormulation
		BeamFormulation = 4;                // default beam formulations
		                                    //   1: continuum mechanics approach (CMF)
		                                    //   2: CMF, elasticity matrix is splitted into a sum of a diagonal part and a shear part, the diagonal part is integrated over the cross section, the shear part is integrated over the beam's axis
		                                    //   3: CMF, like beam_formulation=2, but shear part is neglected here (see paper Duffa)
		                                    //   4: structural mechanics approach (SMF), similar to Reissner's beam theory + shear deformation + crosssection-deformation (which is thickness deformation and "cross-sectional shear")
		calclinear = 0;                     // flag, if linearized strain tensor is used in cont.mech.formulation (CMF)
		order_axial = 4;                    // integration order in axial direction
		order_crosssectional = 2;           // integration order on cross section
	}

	ANCFBeamShear3DGeneric(const ANCFBeamShear3DGeneric& e):Body3D(e.mbs),massmatrix(), Hmatrix(), 
		xg(), xgd(), q0() 
	{
		CopyFrom(e);
	};

	virtual Element* GetCopy() = 0;   // to be defined in derived classes

	virtual void CopyFrom(const Element& e)  //new variable -> copy it here
	{
		Body3D::CopyFrom(e);
		const ANCFBeamShear3DGeneric& ce = (const ANCFBeamShear3DGeneric&)e;
		xg = ce.xg;           // temporary vector for faster computation
		xgd = ce.xgd;         // temporary vector for faster computation
		massmatrix = ce.massmatrix; // storage for mass matrix
		Hmatrix = ce.Hmatrix;       // storage for matrix H

		//integration order
		order_axial = ce.order_axial;
		order_crosssectional = ce.order_crosssectional;

		q0 = ce.q0;            //coordinates in reference configuration
		penalty=ce.penalty;    //penalty is used for an artificial weighting of the summation terms of the energy potential in structural mechanics formulation (SMF)

		nnodes=ce.nnodes;                   //number of nodes of this element (2 for ANCFBeamShear3DLinear, 3 for ANCFBeamShear3DQuadratic)
		n1 = ce.n1; n2 = ce.n2; n3=ce.n3;   //ANCFBeamShear3DLinear owns only node n1 and n2, whereas ANCFBeamShear3DQuadratic has all three nodes
		materialnum = ce.materialnum;       //assign material to element by materialnumber (stored in the mbs)

		perform_reduced_integration=ce.perform_reduced_integration;      //flag for reduced integration (formerly "RI"), used in the 3 variations of cont.mech.formulation (CMF), see BeamFormulation
		BeamFormulation=ce.BeamFormulation; //different beam formulations:
		                                    //1: continuum mechanics approach (CMF)
		                                    //2: CMF, elasticity matrix is splitted into a sum of a diagonal part and a shear part, the diagonal part is integrated over the cross section, the shear part is integrated over the beam's axis
		                                    //3: CMF, like beam_formulation=2, but shear part is neglected here (see paper Duffa)
		                                    //4: structural mechanics approach (SMF), similar to Reissner's beam theory + shear deformation + crosssection-deformation (which is thickness deformation and "cross-sectional shear")
		calclinear=ce.calclinear;    // flag, if linearized strain tensor is used in cont.mech.formulation (CMF)
		is_straight_beam_in_reference_configuration = ce.is_straight_beam_in_reference_configuration;
	}

	virtual void Initialize()
	{
		// TODO for speedup: all precomputations at integration points (see XGProviderCached<data_size>)
		Body3D::Initialize();
	}

	virtual TKinematicsAccessFunctions GetKinematicsAccessFunctions(int mode) const
	{
		TKinematicsAccessFunctions kaf = Body3D::GetKinematicsAccessFunctions(mode);
		return TKinematicsAccessFunctions(kaf+TKAF_node_position+TKAF_angular_velocity+TKAF_rotation_matrix+TKAF_rotation_matrix_P+TKAF_D_pos_D_q+TKAF_D_pos_D_x+TKAF_D_rot_D_q+TKAF_D_rot_v_D_q+TKAF_int_D_u_D_q);
	}

	virtual void LinkToElements() = 0;
	virtual void BuildDSMatrices();
	virtual const char* GetElementSpec() const {return "ANCFBeamShear3DGeneric";}

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual int NS() const {return 3*nnodes;};     // number of shape functions
	virtual int SOS() const {return Dim()*NS();};  // size of K and M  
	virtual int SOSowned() const {return 0;};          // len(u)
	virtual int ES() const  {return 0;};           // size of first order explicit equations
	virtual int IS() const  {return 0;};           // implicit (algebraic) size

	virtual int DataS() const {return 0;}          // number of internal variables
	
	virtual int NNodes() const 
	{
		if (SOSowned() == 0)
		{
			return nnodes;
		}
		else
		{
			return 0;
		}
	};

	virtual void SetNNodes(int nnodes_i)
	{
		nnodes = nnodes_i;
	}

	virtual const int& NodeNum(int i) const = 0;
	virtual int& NodeNum(int i) = 0;

	virtual Node& GetNode(int i) {return GetMBS()->GetNode(NodeNum(i));}

	virtual void SetMaterialNum(int matnr)
	{
		materialnum = matnr;
	}

	virtual const Vector& GetQ0() {return q0;}

	//void SetBeamParameters(double beamEIi, double beamEAi, double beamRhoAi)
	//{
	//	GetMaterial().BeamEIy() = beamEIi;
	//	GetMaterial().BeamEA() = beamEAi;
	//	GetMaterial().BeamRhoA() = beamRhoAi;
	//}
	//int IsBeamParameters() const {if (GetBeamEIy() < 0) return 0; else return 1;}

	virtual int IsRigid() const {return 0;}

	// (AD) changed () to .Get()
	virtual const double& XGP(int iloc) const {return GetXact(ltg.Get(iloc+SOS()));}
	virtual const double& XGPD(int iloc) const {return mbs->GetDrawValue(ltg.Get(iloc+SOS()));}

	//Shapefunctions
	virtual void GetShapes(Vector& sf, const Vector3D& ploc) const = 0;
	virtual void GetShapesx(Vector& sfy, const Vector3D& ploc) const = 0;
	virtual void GetShapesy(Vector& sfy, const Vector3D& ploc) const = 0;
	virtual void GetShapesz(Vector& sfy, const Vector3D& ploc) const = 0;
	virtual void GetShapesxy(Vector& sfy, const Vector3D& ploc) const = 0;
	virtual void GetShapesxz(Vector& sfy, const Vector3D& ploc) const = 0;
	virtual double GetSF(int i, const Vector3D& ploc) const = 0;
	virtual double GetSFx(int i, const Vector3D& ploc) const = 0;
	virtual double GetSFy(int i, const Vector3D& ploc) const = 0;
	virtual double GetSFxy(int i, const Vector3D& ploc) const = 0;
	virtual double GetSFz(int i, const Vector3D& ploc) const = 0;
	virtual double GetSFxz(int i, const Vector3D& ploc) const = 0;



	//Dimensions
	virtual double GetLx() const {return size.X();}
	virtual double GetLy() const {return size.Y();}
	virtual double GetLz() const {return size.Z();}

	//for reduced-integration-parameter:
	virtual int GetReducedIntegration() const {return perform_reduced_integration;}
	virtual void SetReducedIntegration(int val)  //Get-Set-Fct
	{
		perform_reduced_integration = val;
	}
	virtual void SetPenalty(int i, double val)
	{
		penalty(i) = val;
	}
	virtual void SetBeamFormulation(int nr)
	{
		BeamFormulation = nr;
	}
	virtual bool InContinuumMechanicsFormulation()
	{
		return (BeamFormulation < 4);
	}
	virtual void SetCalcLinear(int Calc_linear)
	{
		calclinear = Calc_linear;
	}
	virtual void SetOrderAxial(int int_order_dx)
	{
		order_axial = int_order_dx;
	}
	virtual void SetOrderCrossSectional(int int_order_dydz)
	{
		order_crosssectional = int_order_dydz;
	}

	//Velocity, Displacement
	
	//typedef enum {TCDRI_compute = 0, TCDRI_draw = 1, TCDRI_ref_conf = 2, TCDRI_init = 2} TCDRI_flag;
	//
	//userelement.h:
	//virtual Vector3D GetPos_CDRI(const Vector3D& p_loc, TCDRI_flag flag = TCDRI_compute) const;
	//element.h:
	//virtual Vector3D GetPos(const Vector3D& p_loc) const;
	//virtual Vector3D GetPosD(const Vector3D& p_loc) const;
	//virtual Vector3D GetInitPos(const Vector3D& p_loc) const;
	//virtual Vector3D GetRefConfPos(const Vector3D& p_loc) const;


	virtual Vector3D GetVel3D(const Vector3D& p_loc, int flagD=0) const;
	virtual Vector3D GetVel(const Vector3D& p_loc) const {return GetVel3D(p_loc);};
	virtual Vector3D GetVelD(const Vector3D& p_loc) const {return GetVel3D(p_loc, 1);};
	virtual Vector3D GetDisplacementD(const Vector3D& p_loc) const;
	virtual Vector3D GetDisplacement(const Vector3D& p_loc) const;
	virtual void GetIntDuDq(Matrix& dudq);
	//virtual Vector3D GetDisplacement3DD(const Vector3D& p_loc) const;
	//virtual Vector3D GetDisplacement3D(const Vector3D& p_loc, const Vector& xg) const;
	//virtual void GetIntDuDq(Matrix& dudq);

	//for drawing of fores and vectors (see MBSload!!)
	//virtual Vector3D GetPosD(const Vector3D& p_loc) const;

	//r,rx,rxx,rxy,ry (flagD for Visualization)
	virtual Vector3D GetPos3D(const Vector3D& p_loc, const Vector& xg) const;
	virtual Vector3D GetPos3D(const Vector3D& p_loc, int flagD) const;
	virtual Vector3D GetPosD(const Vector3D& p_loc) const;
	virtual Vector3D GetPos(const Vector3D& p_loc) const { return GetPos3D(p_loc, 0); }
	virtual Vector3D GetPosx3D(const Vector3D& p_loc, const Vector& xg) const;

	//virtual Vector3D GetdPosdx(const Vector3D& p_loc, Vector3D& dpdx) const { return GetPosx3D(p_loc, 0); }

	virtual Vector3D GetPosx3D(const Vector3D& p_loc, int flagD=0) const;
	virtual Vector3D GetPosy3D(const Vector3D& p_loc, const Vector& xg) const;
	virtual Vector3D GetPosy3D(const Vector3D& p_loc, int flagD=0) const;
	virtual Vector3D GetPosz3D(const Vector3D& p_loc, const Vector& xg) const;
	virtual Vector3D GetPosz3D(const Vector3D& p_loc, int flagD=0) const;
	virtual Vector3D GetPosxy3D(const Vector3D& p_loc, const Vector& xg) const;
	virtual Vector3D GetPosxy3D(const Vector3D& p_loc, int flagD=0) const;
	virtual Vector3D GetPosxz3D(const Vector3D& p_loc, const Vector& xg) const;
	virtual Vector3D GetPosxz3D(const Vector3D& p_loc, int flagD=0) const;

	virtual Vector3D GetPosyP3D(const Vector3D& p_loc, int flagD=0) const;
	virtual Vector3D GetPoszP3D(const Vector3D& p_loc, int flagD=0) const;

	//for sliding joint; $JG2012-02-21
	virtual void GetdPosdx(const Vector3D& ploc, Vector3D& dpdx);
	
	//r0 in reference element
	//sollte Ref-Conf heißen
	virtual Vector3D GetInitPos3D(const Vector3D& p_loc) const;
	virtual Vector3D GetInitPosx3D(const Vector3D& p_loc) const;
	virtual Vector3D GetInitPosy3D(const Vector3D& p_loc) const;
	virtual Vector3D GetInitPosz3D(const Vector3D& p_loc) const;
	virtual Vector3D GetInitPosxy3D(const Vector3D& p_loc) const;
	virtual Vector3D GetInitPosxz3D(const Vector3D& p_loc) const;

	//Plotscaling: l/2, h/2
	virtual Vector3D GetPos3D0D(const Vector3D& p_loc) const;
	virtual Vector3D GetPos3D0D(const Vector3D& p_loc, double defscale) const;

	//----------------
	//GetRotMatrix & derivative 
	//----------------

	virtual Matrix3D GetRotMatrix(const Vector3D& p_loc) const { return GetRotMatrix(p_loc, 0); }
	virtual Matrix3D GetRotMatrixD(const Vector3D& p_loc) const { return GetRotMatrix(p_loc, 0, 1); }
	virtual Matrix3D GetRotMatrix(const Vector3D& p_loc, int on_reference_element, int flagD=0) const;
	
	virtual Matrix3D GetRotMatrixP(const Vector3D& p_loc) const {return GetRotMatrixP(p_loc, 0);};
	virtual Matrix3D GetRotMatrixPD(const Vector3D& p_loc) const {return GetRotMatrixP(p_loc, 1);};
	virtual Matrix3D GetRotMatrixP(const Vector3D& p_loc, int flagD) const; // { mbs->UO() << "Error: ANCFBeamShear3D::GetRotMatrixP not implemented yet!\n"; assert(0); return Matrix3D(0.);}

	virtual Vector3D GetTwistAndCurvature(const Vector3D& p_loc) const { return GetTwistAndCurvature(p_loc, 0); }
	virtual Vector3D GetTwistAndCurvature(const Vector3D& p_loc, int on_reference_element, int flagD=0) const;
	virtual Vector3D GetAngularVel(const Vector3D& p_loc) const;

	bool IsStraightBeamInReferenceConfiguration() const
	{
		return is_straight_beam_in_reference_configuration;
	}

	virtual void GetdRotvdqT(const Vector3D& vloc, const Vector3D& ploc, Matrix& drotvdqT) { return GetdRotvdqT(vloc, ploc, drotvdqT, 0); }
	void GetdRotvdqT(const Vector3D& vloc, const Vector3D& ploc, Matrix& drotvdqT, int use_transposed_rotmatrix) const;
	void GetdRotdqT(const Vector3D& ploc, Matrix& d);

	//----------------
	//rot matrix columns A = (e_1 e_2 e_3) and their derivatives with respect to xi and diffvars, where
	//----------------


	// caution: ei denotes the frame vectors, which are not normalized!!!!!!!!!!!!!!
#include"ANCFBeamShear3DGetddeidxidDiffVar.h"
#include"ANCFBeamShear3DGetdeidDiffVar.h"
#include"ANCFBeamShear3DGetdeidxi.h"
#include"ANCFBeamShear3DGetei.h"


	double BeamGAkyz(const double beamGAky, const double beamGAkz)
	// ACHTUNG: Literaturstudie notwendig - dieser Faktor wird für den letzten Term bzgl. Dickendeformation verwendet (GAkyz*2E_yz*2deltaE_yz).
	{
		return 0.5*(beamGAky + beamGAkz);
	}

	double Power(double x, double y) const
	{
		if (y == 2) return x*x;
		if (y == 3) return x*x*x;
		else return pow(x,y);
	}
	double Sqrt(double x) const
	{
		if (x < 0) mbs->UO() << "Sqrt(" << x << ")\n";
		return sqrt(x);
	}

	virtual void GetJacobi(Matrix3D& jac, const Vector3D& ploc, const Vector& xg0) const;


	virtual void GetCoordinates(Vector& dc) const
	{
		for (int i = 1; i <= SOS(); i++)
			dc(i) = XG(i);
	}

	virtual void GetCoordinatesP(Vector& dc) const
	{
		for (int i = 1; i <= SOS(); i++)
			dc(i) = XGP(i);
	}

	virtual void GetDrawCoordinates(Vector& dc) const
	{
		for (int i = 1; i <= SOS(); i++)
			dc(i) = XGD(i);
	}

	virtual void GetDrawCoordinatesP(Vector& dc) const
	{
		for (int i = 1; i <= SOS(); i++)
			dc(i) = XGPD(i);
	}

	virtual Vector3D GetRefPosD() const {return GetPos3D(Vector3D(0.,0.,0.),1);};//

	virtual Vector3D GetNodeLocPos(int i) const; //returns LOCAL position of i-th node

	virtual Vector3D GetNodePos(int i) const; //returns position of i-th node
	virtual Vector3D GetNodePosD(int i) const; //returns position of i-th node (draw mode)

	virtual Vector3D GetDOFPosD(int idof) const; //returns position of i-th DOF
	virtual Vector3D GetDOFDirD(int idof) const; //returns direction of action of i-th DOF
	
	virtual void SetComputeCoordinates()
	{
		for (int i = 1; i <= SOS(); i++)
			xg(i) = XG(i);
	}
	virtual void SetComputeCoordinatesP()
	{
		for (int i = 1; i <= SOS(); i++)
			xg(i) = XGP(i);
	}

	virtual double GetBeamEIy() const {return GetMaterial().BeamEIy(); }
	virtual double GetBeamEA() const {return GetMaterial().BeamEA(); }
	virtual double GetRhoA() const {return GetMaterial().BeamRhoA(); }

	virtual void GetdPosdqT(const Vector3D& ploc, Matrix& d)
	{
		d.SetSize(NS()*Dim(),Dim());
		d.FillWithZeros();

		for (int i = 1; i <= NS(); i++)
		{
			d((i-1)*Dim()+1,1)=GetSF(i,ploc);
			d((i-1)*Dim()+2,2)=GetSF(i,ploc);
			d((i-1)*Dim()+3,3)=GetSF(i,ploc);
		}
	}

	virtual void Gradu(const Vector3D& ploc, const Vector& u, Matrix3D& gradu) const;

	virtual void EvalM(Matrix& m, double t);
	
	virtual void GetGamma(const double& xi, Vector3D& gamma, int flagD=0);
	virtual void GetDeltaGamma(const double& xi, Matrix& dgammadq, Vector3D& gamma);	
	virtual void GetKappa(const double& xi, Vector3D& kappa, int flagD=0);
	virtual void GetDeltaKappa(const double& xi, Matrix& dkappa, Vector3D& kappa);
	virtual void GetThicknessStrain(const double& xi, Vector3D& thicknessstrain, int flagD=0);
	virtual void GetDeltaThicknessStrain(const double& xi, Matrix& dthicknessstraindq, Vector3D& thicknessstrain);

	virtual void EvalF2(Vector& f, double t);
	//residual of variational formulation in continuum mechanics formulation
	virtual void EvalF2CM(Vector& f, double t);
	//residual of variational formulation in structural mechanics formulation
	virtual void EvalF2SM(Vector& f, double t);

	virtual double PostNewtonStep(double t);		// changes plastic strains, returns yieldfunction(sigma)/Em

	virtual void GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables);
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & local_position, bool flagD);
	virtual void DrawElement();


protected:
	
	// deprecated set function!!!
	void SetANCFBeamShear3DGeneric(int materialnumi, const Vector3D& si, const Vector3D& coli);
	
	// default set function
	// gets cross section size via beam3dproperties, and calculates length via CalculateElementLength()
	void SetANCFBeamShear3DGeneric(int materialnumi, const Vector3D& coli);
	double CalculateElementLength() const;

	//temporary storage for acceleration:
	Vector xg, xgd; // [e]; e,x
	Matrix massmatrix, Hmatrix; //M = int(rho*((S)^T).S, dV,V); H = int(S,dV,V)
	int is_straight_beam_in_reference_configuration; //$EDC$[varaccess,EDCvarname="straight_beam",int_bool,EDCfolder="ShearBeam",tooltiptext="is straight beam in reference configuration"]

	//element data variables
	// old:
	//double lx;		//EDC$[varaccess,EDCvarname="lx",EDCfolder="Geometry",tooltiptext="beam dimension in x-direction"]
	//double ly;		//EDC$[varaccess,EDCvarname="ly",EDCfolder="Geometry",tooltiptext="beam dimension in y-direction"]
	//double lz;		//EDC$[varaccess,EDCvarname="lz",EDCfolder="Geometry",tooltiptext="beam dimension in z-direction"]
	// new:
	//EDC Vector3D size;	//$EDC$[varaccess,minval=0,EDCvarname="body_dimensions",EDCfolder="Geometry",tooltiptext="dimensions of the beam. [L_x (length), L_y (width), L_z (height)]"]
	//EDC int materialnum; //$EDC$[varaccess,EDCvarname="material_number",EDCfolder="Physics",tooltiptext="material number which contains the main material properties of the beam"] //material number in MBS which contains material information

	int perform_reduced_integration;	//$EDC$[varaccess,EDCvarname="reduced_integration",EDCfolder="ShearBeam",int_bool,tooltiptext="reduced integration in cont. mech. formulation (CMF)"]
	int BeamFormulation;		//$EDC$[varaccess,EDCvarname="beamformulation",EDCfolder="ShearBeam",tooltiptext="2 = CMF, 4 = SMF"]
	//tooltiptext="1=CMF standard, 2=CMF locking free, 3=CMF nu=0, 4=SMF Simo Reissner"
	int calclinear;		//$EDC$[varaccess,EDCvarname="calc_linear",EDCfolder="ShearBeam",int_bool,tooltiptext="linearized strain computation in cont. mech. formulation (CMF)"]

	int n1; //$EDC$[varaccess,EDCvarname="node_number1",EDCfolder="Geometry",tooltiptext="global number of node 1 (left), node must already exist"] 
	int n2; //$EDC$[varaccess,EDCvarname="node_number2",EDCfolder="Geometry",tooltiptext="global number of node 2 (right), node must already exist"] 

	int n3; // just for quadratic

  int nnodes;				//$EDC$[varaccess,readonly,EDCvarname="nnodes",EDCfolder="ShearBeam",tooltiptext="number of nodes"]

	//integration points
	//Vector x1,w1,x2,w2,x3,w3; //memory for numerical integration
	//order of num. Int. in x-direction:
	int order_axial;						//$EDC$[varaccess,EDCvarname="integration_order_axial",EDCfolder="ShearBeam",tooltiptext="axial integration order"]
	int order_crosssectional;		//$EDC$[varaccess,EDCvarname="integration_order_cross_section",EDCfolder="ShearBeam",tooltiptext="cross section integration order, taks effect only in cont. mech. formulation (CMF)"]

	Vector q0; //initial vector for positions
	//EDC Vector q0;		//$EDC$[varaccess,EDCvarname="node1_reference_position",EDCfolder="Initialization",vecstart=1,vecend=3,tooltiptext="position of node 1 (left) in reference configuration."]
	//EDC Vector q0;		//$EDC$[varaccess,EDCvarname="node1_reference_slope_2",EDCfolder="Initialization",vecstart=4,vecend=6,tooltiptext="slope vector 2 of node 1 (left) in reference configuration."]
	//EDC Vector q0;		//$EDC$[varaccess,EDCvarname="node1_reference_slope_3",EDCfolder="Initialization",vecstart=7,vecend=9,tooltiptext="slope vector 3 of node 1 (left) in reference configuration."]
	//EDC Vector q0;		//$EDC$[varaccess,EDCvarname="node2_reference_position",EDCfolder="Initialization",vecstart=10,vecend=12,tooltiptext="position of node 2 (right) in reference configuration."]
	//EDC Vector q0;		//$EDC$[varaccess,EDCvarname="node2_reference_slope_2",EDCfolder="Initialization",vecstart=13,vecend=15,tooltiptext="slope vector 2 of node 2 (right) in reference configuration."]
	//EDC Vector q0;		//$EDC$[varaccess,EDCvarname="node2_reference_slope_3",EDCfolder="Initialization",vecstart=16,vecend=18,tooltiptext="slope vector 3 of node 2 (right) in reference configuration."]


	Vector penalty;    //penalty for specific terms in SMF
	//EDC Vector penalty;		//$EDC$[varaccess,EDCvarname="penalty_kappa",EDCfolder="ShearBeam",vecstart=1,vecend=3,tooltiptext="penalty term for kappa [kappa1,kappa2,kappa3]"]
	//EDC Vector penalty;		//$EDC$[varaccess,EDCvarname="penalty_gamma",EDCfolder="ShearBeam",vecstart=4,vecend=6,tooltiptext="penalty term for gamma [gamma1,gamma2,gamma3]"]
	//EDC Vector penalty;		//$EDC$[varaccess,EDCvarname="penalty_E",EDCfolder="ShearBeam",vecstart=1,vecend=3,tooltiptext="penalty term for green lagrange strains (E) [Eyy,Ezz,Eyz]"]

	//EDC double density;				//$EDC$[varaccess,remove,EDCvarname="density",EDCfolder="Physics"]
	//EDC double mass;				//$EDC$[varaccess,remove,EDCvarname="mass",EDCfolder="Physics"]
	//EDC double damping_m;		//$EDC$[varaccess,remove,EDCvarname="mass_prop_damping",EDCfolder="Physics"]
};//$EDC$[endclass,ANCFBeamShear3DGeneric]


#pragma endregion


#pragma region ANCFBeamShear3DLinear
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ANCFBeamShear3DLinear ANCFBeamShear3DLinear ANCFBeamShear3DLinear ANCFBeamShear3DLinear
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ANCFBeamShear3DLinear: public ANCFBeamShear3DGeneric //$EDC$[beginclass,classname=ANCFBeamShear3DLinear,parentclassname=ANCFBeamShear3DGeneric,addelementtype=TAEBody,addelementtypename=ANCFBeamShear3DLinear,
//texdescription="ANCFBeamShear3DLinear is an ANCF beam finite element for multibody dynamics systems which is capable of large deformations and can be used for static as well as dynamic investigations. The beam finite element can reproduce axial, bending, shear and torsional deformation. A linear interpolation for the geometry and the displacement along the beam axis is chosen.\\
//The definition of the beam finite element is based on the absolute nodal coordinate formulation (ANCF), which uses slope vectors for the parameterization of the orientation of the cross section instead of rotational parameters. Two different formulations for the elastic forces of the beam elements are presented:\\
//(1) A structural mechanics based formulation of the elastic forces based on Reissner's nonlinear rod theory including generalized strain measures. A term accounting for thickness and cross section deformation is included and shear locking is prevented.\\
//(2) A continuum mechanics based formulation of the elastic forces for a St.Venant Kirchhoff material which avoids the Poisson and shear locking phenomenon.",
//figure="Geometry_linear, The geometric description of the elements is based on a position vector $\mathbf{r}^{(i)}$ and two slope vectors $\mathbf{r}_{,\eta}^{(i)}$ and $\mathbf{r}_{,\zeta}^{(i)}$ in the $i$-th node. These vectors are defined on a scaled and straight reference element, given in coordinates $(\xi,\eta,\zeta)$.",
//texdescriptionNode="The element needs 2 nodes of type 'Node3DS2S3'. The element is described by two nodes at the end points of the beam (node 1 = left node, node 2 = right node). See Fig. \ref{ANCFBeamShear3DLinearfigure1} for a sketch of the two-noded linear beam element and the degrees of freedom per node.",
//texdescriptionDOF="The degrees of freedom of the $i$-th node are the nodal displacements and change of slope vectors and read as follows
//\begin{equation}
//\mathbf{q}^{(i)} = [\mathbf{u}^{(i)\,T}\quad \mathbf{u}^{(i)\,T}_{,\eta}\quad \mathbf{u}^{(i)\,T}_{,\zeta}]^T.
//\end{equation}
//Hence, nine degrees of freedom are specified in each node, therefore the two-noded linear beam element has 18 degrees of freedom.",
//texdescriptionGeometry="The deformed geometry of the ANCF beam finite elements is defined by position and two slope vectors in each node, see Fig. \ref{ANCFBeamShear3DLinearfigure1}. The slope vectors $\mathbf{r}_{,\eta}^{(i)}$ and $\mathbf{r}_{,\zeta}^{(i)}$ are no unit vectors, therefor a cross section deformation is not prohibited. The displacement along the beam axis is interpolated with linear shape functions, while the orientation of the cross section is interpolated linearly. The slope vectors are the derivative vectors with respect to the coordinate system of the scaled straight reference element, see Fig. \ref{fig:ANCFBeamShear3DLinearfigure2}.",
//figure="DiffConfig_linear, Different configurations of the finite beam element: (a) scaled straight reference element and (b) the reference element depicted in the global coordinate system.",
//modus="{CMF}{The definition of the elastic forces is based on a continuum mechanics based formulation for a St.Venant Kirchhoff material using the relation between the nonlinear Green-Lagrange strain  tensor and the second Piola-Kirchhoff stress tensor. The beam is defined as any other solid finite element and the volume integration can be chosen via the variables order\_axial and order\_crosssectional in this modus. Using the parameter perform\_reduced\_integration, the standard integration order is reduced, in order to reduce stiffening effects or computation time.}",
//modus="{SMF}{The definition of the elastic forces is based on a structural mechanics based formulation based on Reissner's nonlinear rod theory including generalized strain measures, namely the axial strain, the shear strains, the torsional strain, and the bending strains. The integration along the beam axis is performed as follows: two Lobatto integration points are used for the integration of the elastic forces covering cross section deformation and one Gauss point is used for the integration of the terms accounting for axial deformation, bending, shear and torsion.}",
//example="ANCFBeamShear3DLinear.txt",
//reference="K. Nachbagauer, P. Gruber, J. Gerstmayr. Structural and Continuum Mechanics Approaches for a 3D Shear Deformable ANCF Beam Finite Element: Application to static and linearized dynamic examples. Journal for Computational and Nonlinear Dynamics, 8, 021004, DOI:10.1115/1.4006787, 2012.",
//reference="K. Nachbagauer. Development of shear and cross section deformable beam finite elements applied to large deformation and dynamics problems, Johannes Kepler University Linz, 2012.",
//texdescriptionComments="In general: For further details on the definition of the elastic forces, the strain measures or the cross section deformation see reference \cite{ANCFBeamShear3DLinearreference1}.\\
//Cross section deformation: In order to penalize a possible cross section deformation of the beam, an additional term is added to the classical strain energy and can be varied by the penalization factors named penalty. See reference \cite{ANCFBeamShear3DLinearreference1} for more details.
//Examples: Find static and linearized dynamic applications of the beam element as well as nonlinear dynamic examples and buckling tests in reference \cite{ANCFBeamShear3DLinearreference2}."]
{
public:

	virtual void ElementDefaultConstructorInitialization();

	ANCFBeamShear3DLinear(MBS* mbsi) : ANCFBeamShear3DGeneric(mbsi)
	{
		ElementDefaultConstructorInitialization();
		nnodes = 2;
	}

	ANCFBeamShear3DLinear(const ANCFBeamShear3DLinear& e) : ANCFBeamShear3DGeneric(e.mbs)
	{
		CopyFrom(e);
	};

	virtual int SetElementData(ElementDataContainer& edc)
	{
		int rv = ANCFBeamShear3DGeneric::SetElementData(edc);
		SetANCFBeamShear3DLinear(n1, n2, materialnum, col);
		return rv;
	}

	// default set function
	void SetANCFBeamShear3DLinear(int n1i, int n2i, int materialnumi, const Vector3D& coli);

	// deprecated set function!!!
	void SetANCFBeamShear3DLinear(const Vector& xc1, const Vector& xc2, int n1i, int n2i, int materialnumi, const Vector3D& si, const Vector3D& coli);

	virtual Element* GetCopy()
	{
		Element* ec = new ANCFBeamShear3DLinear(*this);
		return ec;
	}

	virtual const char* GetElementSpec() const {return "ANCFBeamShear3DLinear";}

	virtual const int& NodeNum(int i) const
	{
		if (i == 1) return n1;
		else if(i == 2) return n2;
		else 
		{
			mbs->UO() << "Error in ANCFBeamShear3DLinear::NodeNum: node " << i << " does not exist!";
			return n1;
		}
	}

	virtual int& NodeNum(int i)
	{
		if (i == 1) return n1;
		else if(i == 2) return n2;
		else 
		{
			mbs->UO() << "Error in ANCFBeamShear3DLinear::NodeNum: node " << i << " does not exist!";
			return n1;
		}
	}

	//Link local with global node numbers
	virtual void LinkToElements();

	//Shapefunctions
	virtual void GetShapes(Vector& sf, const Vector3D& ploc) const;
	virtual double GetSF(int i, const Vector3D& ploc) const;
	virtual void GetShapesx(Vector& sf, const Vector3D& ploc) const;
	virtual double GetSFx(int i, const Vector3D& ploc) const;
	virtual void GetShapesy(Vector& sf, const Vector3D& ploc) const;
	virtual double GetSFy(int i, const Vector3D& ploc) const;
	virtual void GetShapesxy(Vector& sf, const Vector3D& ploc) const;
	virtual double GetSFxy(int i, const Vector3D& ploc) const;
	virtual void GetShapesz(Vector& sf, const Vector3D& ploc) const;
	virtual double GetSFz(int i, const Vector3D& ploc) const;
	virtual void GetShapesxz(Vector& sf, const Vector3D& ploc) const;
	virtual double GetSFxz(int i, const Vector3D& ploc) const;

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter
}; //$EDC$[endclass,ANCFBeamShear3DLinear]

#pragma endregion


#pragma region ANCFBeamShear3DQuadratic
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ANCFBeamShear3DQuadratic ANCFBeamShear3DQuadratic ANCFBeamShear3DQuadratic ANCFBeamShear3DQuadratic
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ANCFBeamShear3DQuadratic: public ANCFBeamShear3DGeneric //$EDC$[beginclass,classname=ANCFBeamShear3DQuadratic,parentclassname=ANCFBeamShear3DGeneric,addelementtype=TAEBody,addelementtypename=ANCFBeamShear3DQuadratic,
//texdescription="ANCFBeamShear3DQuadratic is an ANCF beam finite element for multibody dynamics systems which is capable of large deformations and can be used for static as well as dynamic investigations. The beam finite element can reproduce axial, bending, shear and torsional deformation. A quadratic interpolation for the geometry and the displacement along the beam axis is chosen.\\
//The definition of the beam finite element is based on the absolute nodal coordinate formulation (ANCF), which uses slope vectors for the parameterization of the orientation of the cross section instead of rotational parameters. Two different formulations for the elastic forces of the beam elements are presented:\\
//(1) A structural mechanics based formulation of the elastic forces based on Reissner's nonlinear rod theory including generalized strain measures. A term accounting for thickness and cross section deformation is included and shear locking is prevented.\\
//(2) A continuum mechanics based formulation of the elastic forces for a St.Venant Kirchhoff material which avoids the Poisson and shear locking phenomenon.",
//figure="Geometry_quadratic, The geometric description of the elements is based on a position vector $\mathbf{r}^{(i)}$ and two slope vectors $\mathbf{r}_{,\eta}^{(i)}$ and $\mathbf{r}_{,\zeta}^{(i)}$ in the $i$-th node. These vectors are defined on a scaled and straight reference element, given in coordinates $(\xi,\eta,\zeta)$.",
//texdescriptionNode="The element needs 3 nodes of type 'Node3DS2S3'. The element is described by three nodes: at the end points and the mid point of the beam (node 1 = left node, node 2 = right node, node 3 = mid point). See Fig. \ref{ANCFBeamShear3DLinearfigure1} for a sketch of the three-noded quadratic beam element and the degrees of freedom per node.",
//texdescriptionDOF="The degrees of freedom of the $i$-th node are the nodal displacements and change of slope vectors and read as follows
//\begin{equation}
//\mathbf{q}^{(i)} = [\mathbf{u}^{(i)\,T}\quad \mathbf{u}^{(i)\,T}_{,\eta}\quad \mathbf{u}^{(i)\,T}_{,\zeta}]^T.
//\end{equation}
//Hence, nine degrees of freedom are specified in each node, therefore the three-noded quadratic beam element has 27 degrees of freedom.",
//texdescriptionGeometry="The deformed geometry of the ANCF beam finite elements is defined by position and two slope vectors in each node, see Fig. \ref{ANCFBeamShear3DLinearfigure1}. The slope vectors $\mathbf{r}_{,\eta}^{(i)}$ and $\mathbf{r}_{,\zeta}^{(i)}$ are no unit vectors, therefor a cross section deformation is not prohibited. The displacement along the beam axis is interpolated with quadratic shape functions, while the orientation of the cross section is interpolated linearly. The slope vectors are the derivative vectors with respect to the coordinate system of the scaled straight reference element, see Fig. \ref{fig:ANCFBeamShear3DLinearfigure2}.",
//figure="DiffConfig_quadratic, Different configurations of the finite beam element: (a) scaled straight reference element and (b) the reference element depicted in the global coordinate system.",
//modus="{CMF}{The definition of the elastic forces is based on a continuum mechanics based formulation for a St.Venant Kirchhoff material using the relation between the nonlinear Green-Lagrange strain  tensor and the second Piola-Kirchhoff stress tensor. The beam is defined as any other solid finite element and the volume integration can be chosen via the variables order\_axial and order\_crosssectional in this modus. Using the parameter perform\_reduced\_integration, the standard integration order is reduced, in order to reduce stiffening effects or computation time.}",
//modus="{SMF}{The definition of the elastic forces is based on a structural mechanics based formulation based on Reissner's nonlinear rod theory including generalized strain measures, namely the axial strain, the shear strains, the torsional strain, and the bending strains. The integration along the beam axis is performed as follows: two Lobatto integration points are used for the integration of the elastic forces covering cross section deformation and one Gauss point is used for the integration of the terms accounting for axial deformation, bending, shear and torsion.}",
//example="ANCFBeamShear3DQuadratic.txt",
//texdescriptionComments="In general: For further details on the definition of the elastic forces, the strain measures or the cross section deformation see reference \cite{ANCFBeamShear3DLinearreference1}.\\
//Cross section deformation: In order to penalize a possible cross section deformation of the beam, an additional term is added to the classical strain energy and can be varied by the penalization factors named penalty. See reference \cite{ANCFBeamShear3DLinearreference1} for more details.\\
//Examples: Find static and linearized dynamic applications of the beam element as well as nonlinear dynamic examples and buckling tests in reference \cite{ANCFBeamShear3DLinearreference2}."]
{
public:

	virtual void ElementDefaultConstructorInitialization();

	ANCFBeamShear3DQuadratic(MBS* mbsi) : ANCFBeamShear3DGeneric(mbsi)
	{
		ElementDefaultConstructorInitialization();
		nnodes = 3;
	}

	ANCFBeamShear3DQuadratic(const ANCFBeamShear3DQuadratic& e) : ANCFBeamShear3DGeneric(e.mbs)
	{
		CopyFrom(e);
	};

	virtual int SetElementData(ElementDataContainer& edc)
	{
		int rv = ANCFBeamShear3DGeneric::SetElementData(edc);
		SetANCFBeamShear3DQuadratic(n1, n2, n3, materialnum, col);
		return rv;
	}

	// default set function
	void SetANCFBeamShear3DQuadratic(int n1i, int n2i, int n3i, int materialnumi, const Vector3D& coli);

	// deprecated set function!!!
	void SetANCFBeamShear3DQuadratic(const Vector& xc1, const Vector& xc2, const Vector& xc3, int n1i, int n2i, int n3i, int materialnumi, const Vector3D& si, const Vector3D& coli);

	virtual Element* GetCopy()
	{
		Element* ec = new ANCFBeamShear3DQuadratic(*this);
		return ec;
	}

	virtual const char* GetElementSpec() const {return "ANCFBeamShear3DQuadratic";}

	virtual const int& NodeNum(int i) const
	{
		if (i == 1) return n1;
		else if(i == 2) return n2;
		else if(i == 3) return n3;
		else 
		{
			mbs->UO() << "Error in ANCFBeamShear3DQuadratic::NodeNum: node " << i << " does not exist!";
			return n1;
		}
	}

	virtual int& NodeNum(int i)
	{
		if (i == 1) return n1;
		else if(i == 2) return n2;
		else if(i == 3) return n3;
		else 
		{
			mbs->UO() << "Error in ANCFBeamShear3DQuadratic::NodeNum: node " << i << " does not exist!";
			return n1;
		}
	}

	//Link local with global node numbers
	virtual void LinkToElements();

	//Shapefunctions
	virtual void GetShapes(Vector& sf, const Vector3D& ploc) const;
	virtual double GetSF(int i, const Vector3D& ploc) const;
	virtual void GetShapesx(Vector& sf, const Vector3D& ploc) const;
	virtual double GetSFx(int i, const Vector3D& ploc) const;
	virtual void GetShapesy(Vector& sf, const Vector3D& ploc) const;
	virtual double GetSFy(int i, const Vector3D& ploc) const;
	virtual void GetShapesxy(Vector& sf, const Vector3D& ploc) const;
	virtual double GetSFxy(int i, const Vector3D& ploc) const;
	virtual void GetShapesz(Vector& sf, const Vector3D& ploc) const;
	virtual double GetSFz(int i, const Vector3D& ploc) const;
	virtual void GetShapesxz(Vector& sf, const Vector3D& ploc) const;
	virtual double GetSFxz(int i, const Vector3D& ploc) const;

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	//EDC Vector q0;		//$EDC$[varaccess,EDCvarname="node3_reference_position",EDCfolder="Initialization",vecstart=19,vecend=21,tooltiptext="position of node 3 (middle) in reference configuration."]
	//EDC Vector q0;		//$EDC$[varaccess,EDCvarname="node3_reference_slope_2",EDCfolder="Initialization",vecstart=22,vecend=24,tooltiptext="slope vector 2 of node 3 (middle) in reference configuration."]
	//EDC Vector q0;		//$EDC$[varaccess,EDCvarname="node3_reference_slope_3",EDCfolder="Initialization",vecstart=25,vecend=27,tooltiptext="slope vector 3 of node 3 (middle) in reference configuration."]
	//EDC int n3;				//$EDC$[varaccess,EDCvarname="node_number3",EDCfolder="Geometry",tooltiptext="global number of node 3 (middle), node must already exist"] 

};//$EDC$[endclass,ANCFBeamShear3DQuadratic]

#pragma endregion