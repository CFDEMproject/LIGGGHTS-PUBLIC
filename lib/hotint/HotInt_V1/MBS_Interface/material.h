//#**************************************************************
//#
//# filename:             material.h           
//#
//# author:               Gerstmayr Johannes, Aigner Larissa
//#
//# generated:						July 2004
//# updated:							September 2010
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
//#**************************************************************

#ifndef MATERIAL__H
#define MATERIAL__H

class FiniteElement3D;
class FiniteElement2D;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                                  MATERIAL
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
typedef enum
{
	TMatEmpty = 0, 
	TMatRigidBody = 1, 
	TMatFlexGeneral = 2, 
	TMat2D = 4, 
	TMatBeam = 8, 
	TMatPlate = 16, 
	TMatShell = 32, 
	TMatBeam2D = 64, 
	//TMatPlast = 128,			//$ PG 2013-7-31: deleted
	//TMatSpecial = 256,		//$ PG 2013-7-31: deleted
	TMatOrthotropic = 128,
	TMatAnisotropic = 256, 
	TMatPlaneStrain = 512, 
	//TMatNonlinElast = 1024	//$ PG 2013-7-31: deleted
} TMBSMaterial;   // these flags are not used/set when material is instantiated by script language 


const int MAXInelasticVariablesCount = 20;


//$EDC$[beginclass,classname=Material,addelementtypename=Material,
//texdescription="Material is the basic Object for defining material properties for standard finite elements (in contrast to structural finite elements such as beams and plates).",
//texdescriptionComments="For static problems define the elastic properties $\mathtt{Solid.youngs\_modulus}$ and $\mathtt{Solid.poisson\_ratio}$, whereas for dynamic problems also $\mathtt{Solid.density}$ is required. If the problem is planar ($\mathtt{Solid.plane}$ is set to $1$), then the plane strain case is assumed unless $\mathtt{Solid.plane\_stress}$ is set to $1$. If the material is inelastic, then also the properties in the subtree $\mathtt{Inelasticity}$ have to be set.",
//example="Material.txt"]
class Material
{
public:
	
	struct OrthotropicConstants
	{
	public:
		double& E1() {return e1;} double& E2() {return e2;} double& E3() {return e3;}
		double& NU12() {return nu12;}	double& NU13() {return nu13;}	double& NU23() {return nu23;}
		double& G12() {return g12;}	double& G13() {return g13;}	double& G23() {return g23;}
		
		const double E1() const {return e1;} const double& E2() const {return e2;} const double& E3() const {return e3;}
		const double NU12() const {return nu12;}	const double& NU13() const {return nu13;}	const double& NU23() const {return nu23;}
		const double G12() const {return g12;}	const double& G13() const {return g13;}	const double& G23() const {return g23;}

	private:
		double e1, e2, e3;
		double nu12, nu13, nu23;
		double g12, g13, g23;
	};

	enum InelasticitySolutionMethod { ISM_Default = 0, ISM_FixedPoint = 1, ISM_ReturnMapping = 2, ISM_ConsistentTangentStiffness = 3, ISM_NotSpecified = 2359819237};
	enum InelasticityType {IT_LinearElastic = 0, IT_ElastoPlastic = 1, IT_NonlinElasticSimoHughes = 2, IT_NotSpecified = 2359819237};

public:
	Material()
	{
		mbs = 0;
		InitConstructor();
	}

	//recommended
	Material(MBS* mbsi)
	{
		mbs = mbsi;
		InitConstructor();
	}

	//deprecated - use Material(MBS* mbsi) instead
	Material(MBS* mbsi, double densityI)
	{
		mbs = mbsi;
		InitConstructor();

		type = TMatRigidBody;
		density = densityI;

		mbs->UO() << "WARNING: constructor deprecated - use __thiscall Material::Material(class MBS *) instead of " << __FUNCSIG__  << ".\n";
	}

	//deprecated - use Material(MBS* mbsi) instead
	Material(MBS* mbsi, double densityI, double youngs_modulusI, double poisson_ratioI, int plane = 0, int planestress = 1)
	{
		mbs = mbsi;
		InitConstructor();

		if (!plane)
		{
			type = TMatFlexGeneral;
		}
		else
		{
			type = TMBSMaterial(TMatFlexGeneral + TMat2D);
		}

		if (!planestress)
		{
			//type = TMatPlaneStrain;
			type = TMBSMaterial(type | TMatPlaneStrain);
		}

		density = densityI;
		youngs_modulus = youngs_modulusI;
		poisson_ratio = poisson_ratioI;
		
		mbs->UO() << "WARNING: constructor deprecated - use __thiscall Material::Material(class MBS *) instead of " << __FUNCSIG__  << ".\n";
	}

	//deprecated - use Material(MBS* mbsi) instead
	Material(MBS* mbsi, double densityI, double youngs_modulusI, double poisson_ratioI, Vector3D color, int plane = 0, int planestress = 1)
	{
		Material(mbsi, densityI, youngs_modulusI, poisson_ratioI, plane, planestress);		
		materialcolor = color;
		
		mbs->UO() << "WARNING: constructor deprecated - use __thiscall Material::Material(class MBS *) instead of " << __FUNCSIG__  << ".\n";
	}

	//depricated - old Set-Function
	virtual void SetMaterial(int nonlinmattypeI, double densityI, double youngs_modulusI, double poisson_ratioI);

	//recommended Set-Functions
	//rigid
	virtual void SetMaterialRigid(double densityI);

	//CM: continuum mechanics
	virtual void SetMaterialCM(double densityI, double youngs_modulusI, double poisson_ratioI, Vector3D color, int plane = 0, int planestress = 1);
	virtual void SetMaterialCMInelastic(double densityI, double youngs_modulusI, double poisson_ratioI, InelasticityType inelasticity_typeI, 
		InelasticitySolutionMethod inelasticity_solution_methodI, double yield_stressI, double tangent_moduleI, Vector3D color, int plane = 0, int planestress = 1);
	virtual void SetMaterialCMOrthotropic(double densityI, Vector3D color, const Material::OrthotropicConstants& oc);
	
	//$ PG 2013-6-6: [ these old functions
	//virtual void SetMaterialSMBeam(double _rho, double _EA, double _EIy, double _EIz=0, double _GJkx=0, double _GAky=0, double GAkz=0, double _rhoA=0, double _rhoIx=0, double _rhoIy=0, double _rhoIz=0);
	//virtual void SetMaterialSMBeam2D(double _rho, double _EA, double _EIy, double _GJkx=0, double _GAky=0, double _rhoA=0, double _rhoIx=0, double _rhoIy=0);
	//$ PG 2013-6-6: have been replaced by Beam3DProperties::SetBeam3DProperties and Beam2DProperties::SetBeam2DProperties ]

	//$ PG 2012-12-19: Set material type and parameters according to user specified edc data
	virtual void SetMaterialFromEDCData();    
		
	//$ PG 2012-12-19: CheckConsistency added
	virtual int CheckConsistency(mystr& errorstr);   //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute

	Material(const Material& s)
	{
		CopyFrom(s);
	};

	Material& operator=(const Material& s) 
	{
		if (this == &s) {return *this;}
		CopyFrom(s);
		return *this;
	}

	virtual Material* GetCopy() const
	{
		Material* ec = new Material();
		ec->CopyFrom(*this);
		return ec;
	}

	virtual void CopyFrom(const Material& s);

	virtual ~Material() {};

	virtual void InitConstructor();

	virtual mystr GetTypeName() const {	return "Material";} //$EDC$[funcaccess,EDCvarname="material_type",tooltiptext="specification of material type. Once the material is added to the mbs, you MUST NOT change this type anymore!"]
	virtual mystr GetTypeSpec() const 
	{

		mystr addon = "";
		//if (type&TMatSpecial) addon += "_special";				//$ PG 2013-7-31: deleted
		//if (type&TMatPlast) addon += "_plast";						//$ PG 2013-7-31: deleted
		//if (type&TMatNonlinElast) addon += "_nlelast";		//$ PG 2013-7-31: deleted
		if (GetInelasticityType() == IT_ElastoPlastic) addon += "_plast";
		if (GetInelasticityType() == IT_NonlinElasticSimoHughes) addon += "_nlelast";

		if (type&TMatRigidBody) return mystr("rigid")+addon;
		if (type&TMatFlexGeneral) return mystr("general3D")+addon;
		if (type&TMat2D) return mystr("2D")+addon;
		if (type&TMatBeam) return mystr("beam3D")+addon;
		//if (type&TMatBeam2D) return mystr("beam2D")+addon;
		if (type&TMatPlate) return mystr("plate")+addon;

		return addon;
	}

	virtual const mystr& GetMaterialName() const {return materialname;}
	virtual mystr& GetMaterialName() {return materialname;}
	virtual Vector3D GetMaterialColor() const {return materialcolor;}
	virtual void SetMaterialColor(Vector3D col) {materialcolor=col;}
	virtual TMBSMaterial GetType() const {return type;}
	virtual void SetType(TMBSMaterial t) {type = t;}
	virtual void AddType(TMBSMaterial t) {type = TMBSMaterial(type|t);}
	virtual const MBS* GetMBS() const {return mbs;}
	virtual MBS* GetMBS() {return mbs;}

	// introduced by YV to avoid frequent initialization of Matrix3D
	static const Matrix3D dummyMatrix3D;
	static const Matrix2D dummyMatrix2D;
	static const Vector dummyVector;

	//compute stress (piola2) from strain
	//* initial strains added by YV
	//* changed Matrix3D inel_strain to Vector inel_variables by PG
	virtual void ComputeStressFromStrain(
		const Matrix3D& strain, Matrix3D& stress,
		const Vector& inel_variables = dummyVector,
		const Matrix3D& initial_strain = dummyMatrix3D
		) const;
	//void Material::ComputeStressFromStrainD(
	//	const Matrix3D& strain, Matrix3D& stress,
	//	const Vector& inel_variables = dummyVector,
	//	const Matrix3D& initial_strain = dummyMatrix3D
	//	);
	virtual void ComputeStressFromStrain2D(
		const Matrix2D& strain, Matrix2D& stress,
		const Vector& inel_variables = dummyVector,
		const Matrix2D& initial_strain = dummyMatrix2D
		) const;
	//virtual double ComputeAveragedContactStress(FiniteElement3D* const fe);	
	//$ PG 2013-8-2: ]
	virtual double ComputeInelasticVariablesUpdateFromStrain(Vector& inel_variables_update, const Vector& inel_variables, const Matrix3D& strain) const;
	virtual double ComputeYieldFunction(const Matrix3D& strain, const Vector& inelastic_variables) const;
	virtual double ComputeYieldFunction2D(const Matrix2D& strain, const Vector& inelastic_variables) const;
	virtual double ComputeInelasticVariablesUpdateFromStrain2D(Vector& inel_variables_update, const Vector& inel_variables, const Matrix2D& strain) const;
	// AP Dez. 2012
	// Material law due to return mapping algorithm by Simo and Hughes, plane stress and strain case
	// plane strain case differs from ComputeInelasticVariablesFromStrain2D by definition of TangentModule
	// plane stress case differs from ComputeInelasticVariablesFromStrain2D substancially
	void ComputeInelasticVariablesFromStrainPlaneStress2D(Vector& inel_variables, const Matrix2D& strain);
	void ComputeInelasticVariablesFromStrainPlaneStrain2D(Vector& inel_variables, const Matrix2D& strain);
	void ComputeInelasticVariablesFromStrain1D(Vector& inel_variables, double strain);


	//PG: This routine addresses geometrical properties only.
	//*   Hence, it should be provided by the element - not by the material.
	virtual void ComputeStrain(Matrix3D& strain, const Matrix3D& grad_u, int isgeomnonlin);
	virtual void ComputeStrain2D(Matrix2D& strain, const Matrix2D& grad_u, int isgeomnonlin);

	
	//6x6 in 3D, or 3x3 in 2D
	virtual void ComputeElasticityMatrix(Matrix& C) const;
	virtual void ComputeConsistentTangentStiffnessMatrix(Matrix& C, const Vector& strain, const Vector& inelastic_variables) const;
	virtual void ComputeConsistentTangentStiffnessMatrix2D(Matrix& C, const Matrix2D& strain, const Vector& inelastic_variables) const;

	//$ PG 2013-8-1: flag TMatNonlinElast and member nonlinmattype have been deleted. nonlinear elasticity is now specified via setting inelasticity_type = IT_NonlinElasticSimoHughes.
	//$ PG 2013-8-1: therefore, the following function may be deleted
	//virtual int IsNonlinElasticMaterial() const {if (type&TMatNonlinElast) return 1; else return 0;}

	virtual int IsOrthotropicMaterial() const {if (type&TMatOrthotropic) return 1; else return 0;}
	virtual void SetOrthotropicConstants(const OrthotropicConstants oc) {orthotropic_constants = oc;}
	virtual const OrthotropicConstants& GetOrthotropicConstants() const {return orthotropic_constants;}


	virtual int IsInelasticMaterial() const { return inelasticity_type != IT_LinearElastic;}
	virtual int GetInelasticVariablesCount() const
	{
		if (!IsInelasticMaterial()) 
			return 0;

		switch (GetInelasticityType())
		{
		case IT_ElastoPlastic:
			if (!IsPlanarMaterial())
				return 7;     // 6*strain+1*hardening_parameter
			else
				return 4;     // 3*strain+1*hardening_parameter
			break;

		default:
			return 0;
		}
	}

	virtual int IsPlaneStrain() const;
	virtual int IsPlanarMaterial() const;

	virtual int GetMaterialNumber() const {return materialnumber;}
	virtual void SetMaterialNumber(int mnum) {materialnumber = mnum;}

	//general:
	virtual const double& Density() const {return density;}
	virtual double& Density() {return density;}
	virtual const double& DampingK() const {return dampingK;}
	virtual double& DampingK() {return dampingK;}
	virtual const double& DampingM() const {return dampingM;}
	virtual double& DampingM() {return dampingM;}

	//flexible:
	virtual const double& YoungsModulus() const {return youngs_modulus;}
	virtual double& YoungsModulus() {return youngs_modulus;}
	virtual const double& PoissonRatio() const {return poisson_ratio;}
	virtual double& PoissonRatio() {return poisson_ratio;}

	virtual double GetMu() const { return youngs_modulus / (2. + 2.*poisson_ratio); }
	virtual double GetLambda() const { return poisson_ratio*youngs_modulus / (1. + poisson_ratio) / (1. - 2*poisson_ratio); }

	//beam:
	virtual bool IsMaterialOfBeam() const;
	virtual bool IsMaterialOfBeamWithRectangularCrossSection() const;
	virtual bool IsMaterialOfBeamWithCircularCrossSection() const;
	virtual double GetBeamThicknessY() const { return IsMaterialOfBeamWithRectangularCrossSection() ? beamCrossSectionSize(1) : 0; }
	virtual double GetBeamThicknessZ() const {return IsMaterialOfBeamWithRectangularCrossSection() ? beamCrossSectionSize(2) : 0; }
	virtual double GetBeamRadius() const {return IsMaterialOfBeamWithCircularCrossSection() ? beamCrossSectionSize(1) : 0; }
//$ PG 2012-12-6:[ BeamA() depricated, but used in old models and elements
	virtual const double& BeamA() const {return beamA;   }
	virtual double& BeamA() {return beamA;   }
//$ PG 2012-12-6:]
	virtual const double& BeamRhoA() const {return beamRhoA;   }
	virtual double& BeamRhoA() {return beamRhoA;   }
	virtual const double& BeamRhoIx() const {return beamRhoIx;   }
	virtual double& BeamRhoIx() {return beamRhoIx;   }
	virtual const double& BeamRhoIy() const {return beamRhoIy;   }
	virtual double& BeamRhoIy() {return beamRhoIy;   }
	virtual const double& BeamRhoIz() const {return beamRhoIz;   }
	virtual double& BeamRhoIz() {return beamRhoIz;   }
	virtual const double& BeamEA() const {return beamEA;   }
	virtual double& BeamEA() {return beamEA;   }
	virtual const double& BeamEIy() const {return beamEIy; }
	virtual double& BeamEIy() {return beamEIy; }
	virtual const double& BeamEIz() const {return beamEIz; }
	virtual double& BeamEIz() {return beamEIz; }
	virtual const double& BeamGAky() const {return beamGAky;}
	virtual double& BeamGAky() {return beamGAky;}
	virtual const double& BeamGAkz() const {return beamGAkz;}
	virtual double& BeamGAkz() {return beamGAkz;}
	virtual const double& BeamGJkx() const {return beamGJkx;}
	virtual double& BeamGJkx() {return beamGJkx;}

	// plastic:
	virtual int IsHardeningMaterial() const {if (tangentmodule>0) return 1; else return 0;}
	virtual double& YieldStress() {return yieldstress;}
	virtual const double& YieldStress() const {return yieldstress;}
	virtual double& TangentModule() {return tangentmodule;}
	virtual const double& TangentModule() const {return tangentmodule;}
	//virtual double DoNonlinStepBeam2D(double& stress, double eps, double& epsplast, double epsplast_old) const;
	virtual double DoNonlinStepBeam2D(double& stress, double eps, double& epsplast, double& internalvar, 
		double epsplast_old, double internalvar_old, int use_tangent_stiffness=0) const;

	//$ PG 2013-8-7: replaced by ComputeInelasticVariablesFromStrain(...)
	//virtual double DoNonlinStep(FiniteElement3D* const fe, const Matrix3D& F, Vector& inel_strain, Vector& inel_strain_old) const 
	//{
	//	//------------------------------
	//	// equations used:             |
	//	//                             |
	//	// int_A (P : d F / d q) dA    |
	//	// P = F * S                   |
	//	// S = D : E                   |
	//	// E = 1/2*(F^T F-I)           |
	//	//------------------------------

	//	return 0;
	//}

	//virtual double DoNonlinStep(FiniteElement2D* const fe, const Matrix2D& F, Vector& inel_strain, Vector& inel_strain_old) const 
	//{
	//	//------------------------------
	//	// equations used:             |
	//	//                             |
	//	// int_A (P : d F / d q) dA    |
	//	// P = F * S                   |
	//	// S = D : E                   |
	//	// E = 1/2*(F^T F-I)           |
	//	//------------------------------

	//	return 0;
	//}

	//$ PG 2013-8-1: deleted nonlinmattype, nonlinear elastic material is now specified via setting inelasticity_type = IT_NonlinElasticSimoHughes
	//$ PG 2013-8-1: the use of the following two functions has been replaced in all models and elements, they are not needed anymore and shall be erased from this file
	//virtual const int& NonlinMatType() const {return nonlinmattype;}
	//virtual int& NonlinMatType() {return nonlinmattype;}

	virtual void GetStress2D(const Vector3D& strain, Vector3D& stress, const Vector3D& inel_strain, const double contactstress = -1) const;
	virtual void GetLinearizedStress2D(const Vector3D& strain, Vector3D& stress, const Vector3D& inel_strain, const Vector3D& stress_lin) const;

	//General access functions:
	virtual int NParameters() const {return 21;} //this should be the number of possible parameters 

	virtual void GetParam(int i, double& value, mystr& name, int& group, int& flag);

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all material data

	virtual int SetElementData(ElementDataContainer& edc); //set material data
	virtual void GetElementDataAuto(ElementDataContainer& edc);
	virtual int SetElementDataAuto(ElementDataContainer& edc);
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); // Read access to a single element variable 
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); // Write access to a single element variable
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); // returns all available directly (ReadWrite-) accessable variables

	virtual mystr GetInelasticityTypeString(InelasticityType it) const;
	virtual InelasticityType GetInelasticityType(const mystr& it_str) const;
	virtual mystr GetInelasticitySolutionMethodString(InelasticitySolutionMethod ism) const;
	virtual InelasticitySolutionMethod GetInelasticitySolutionMethod(const mystr& ism_str) const;
	
	virtual InelasticitySolutionMethod GetInelasticitySolutionMethod() const { return inelasticity_solution_method; }
	virtual mystr GetInelasticitySolutionMethodString() const { return inelasticity_solution_method_str; }
	virtual void SetInelasticitySolutionMethod(const InelasticitySolutionMethod ism)
	{
		inelasticity_solution_method = ism;
		inelasticity_solution_method_str = GetInelasticitySolutionMethodString(ism);
	}	

	virtual InelasticityType GetInelasticityType() const { return inelasticity_type;	}
	virtual mystr GetInelasticityTypeString() const { return inelasticity_type_str; }
	virtual void SetInelasticityType(const InelasticityType it)
	{
		inelasticity_type = it;
		inelasticity_type_str = GetInelasticityTypeString(it); 
	}

protected:

	//general:
	Vector3D materialcolor;				//$EDC$[varaccess,EDCvarname="color",EDCfolder="Graphics",tooltiptext="material color (as yet used with FEMesh, only)"]
	

	MBS* mbs;
	//used for distinguishing base and derived classes: See GetTypeName()
	mystr materialname;						//$EDC$[varaccess,EDCvarname="name",tooltiptext="name of the material"]
	int materialnumber; //material number in MBS system
	TMBSMaterial type;  //fast internal use 

	//solid:
	double density;								//$EDC$[varaccess,EDCvarname="density",EDCfolder="Solid",tooltiptext="density (rho) for gravitational force"]
	double youngs_modulus;				//$EDC$[varaccess,EDCvarname="youngs_modulus",EDCfolder="Solid",tooltiptext="Youngs modulus"]
	double poisson_ratio;					//$EDC$[varaccess,EDCvarname="poisson_ratio",EDCfolder="Solid",tooltiptext="Poisson ratio"]
	int plane;										//$EDC$[varaccess,EDCvarname="plane",EDCfolder="Solid",int_bool,tooltiptext="true: 2D, false: 3D"]
	int plane_stress;							//$EDC$[varaccess,EDCvarname="plane_stress",EDCfolder="Solid",int_bool,tooltiptext="for 2D-Elements only; 1: plane stress, 0: plane strain"]
	double dampingM;
	double dampingK;

	//beam:    //DEPRECATED, should be moved to class Beam3DProperties : public Material and variables corresponding to y-dir should be also copied to Beam2DProperties : public Material
	int beamCrossSectionType;			
	Vector beamCrossSectionSize;
  double beamEA;								
	double beamEIy;								
	double beamEIz;								
	double beamGAky;							
	double beamGAkz;							
	double beamGJkx;							
	double beamRhoA;							
	double beamRhoIx;							
	double beamRhoIy;							
	double beamRhoIz;							

//$ PG 2012-12-6:[ beamA is DEPRECATED, but used in old models and elements
	double beamA;    //cross section area
//$ PG 2012-12-6:]

	//inelasticity:
	double yieldstress;						//$EDC$[varaccess,EDCvarname="yield_stress",EDCfolder="Inelasticity",tooltiptext="Yield Stress s_y, e.g., |dev s| <= s_y"]
	double tangentmodule;					//$EDC$[varaccess,EDCvarname="tangent_module",EDCfolder="Inelasticity",tooltiptext="Modulus of hardening H"]
	
	//$ PG 2013-7-31: merged former 'nonlinmattype' with inelasticity_type
	InelasticityType inelasticity_type;				//$!EDC$[varaccess,EDCvarname="inelasticity_type",EDCfolder="Inelasticity",tooltiptext="0 = not specified, 1 = Prandtl Reuss plasticity + isotropic hardening"]
	mystr inelasticity_type_str;	//$EDC$[varaccess,EDCvarname="inelasticity_type",EDCfolder="Inelasticity",tooltiptext="linear_elastic, elasto_plastic (= Prandtl Reuss plasticity + isotropic hardening), nonlinear_elastic_Simo_Hughes (see Simo and Hughes, Computational Inelasticity 1998: S=lambda/2*(J*J-1)/C + mu*(1-1/C))"] // $ MSax 2013-08-06 : bugfix for docu creation: changed S=lambda/2*(J²-1)*C^(-1) + mu*(1-C^(-1)) --> S=lambda/2*(J*J-1)/C + mu*(1-1/C))

	InelasticitySolutionMethod inelasticity_solution_method;
	mystr inelasticity_solution_method_str; //$EDC$[varaccess,EDCvarname="inelasticity_solution_method",EDCfolder="Inelasticity",tooltiptext="fixed_point, return_mapping, consistent_tangent_stiffness (see Simo and Hughes, Computational Inelasticity 1998)"]

	// nonlinear elasticity 
	//$ PG 2013-7-31: Removed 'nonlinmattype', former nonlinmattype = 1 is now met by inelasticity_type = IT_SimonHughes (see enum InelasticityType)
	//int nonlinmattype;	//type of nonlinear elastic material: 1...Simo and Hughes, Computational Inelasticity 1998: S=lambda/2*(J²-1)*C^(-1) + mu*(1-C^(-1))

	//anisotropy: use these engineering constants:
	OrthotropicConstants orthotropic_constants;

};
//$EDC$[endclass,Material]










////////////////////////////////////////////
// class Beam3DProperties : public Material
////////////////////////////////////////////


//$EDC$[beginclass,classname=Beam3DProperties,parentclassname=Material,addelementtypename=Beam3DProperties,
//texdescription="Beam3DProperties defines material and geometric properties for beam structural finite elements.",
//texdescriptionComments="First, specify the $\mathtt{cross\_section\_type}$ of the beam, which may be either rectangular (if set to $1$), circular (if set to $2$) or polygonal (if set to $3$). In either case the $\mathtt{cross\_section\_size}$ is a vector of $2$, $1$, or $2\,n$ entries, where $n$ confers to the number of vertices of a closed polygon. Then specify the stiffnesses and moments of inertias, as they are neede by your beam and problem.",
//example="Beam3DProperties.txt"]
class Beam3DProperties: public Material
{
public:
	Beam3DProperties() : Material () { InitConstructor(); }
	Beam3DProperties(MBS* mbs) : Material (mbs) { InitConstructor(); }
	Beam3DProperties(const Beam3DProperties& bp) : Material (bp) { CopyFrom(bp); }
	//virtual ~Beam3DProperties() { Material::~Material(); }
	virtual void CopyFrom(const Beam3DProperties& s) 
	{
		Material::CopyFrom(s);
	}
	virtual Beam3DProperties* GetCopy() const
	{
		Beam3DProperties* ec = new Beam3DProperties(mbs);
		ec->CopyFrom(*this);
		return ec;
	}
	virtual mystr GetTypeName() const;
	virtual void InitConstructor();
	virtual void SetMaterialFromEDCData(); 

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual void SetBeam3DProperties(int _beamCrossSectionType /*1: rectangular, 2: circular, 3: polygonal*/, Vector _beamCrossSectionSize, double _rho, double _EA, double _EIy, double _EIz=0, double _GJkx=0, double _GAky=0, double _GAkz=0, double _rhoA=0, double _rhoIx=0, double _rhoIy=0, double _rhoIz=0);
	virtual void SetBeam3DProperties(double _rho, double _EA, double _EIy, double _EIz=0, double _GJkx=0, double _GAky=0, double GAkz=0, double _rhoA=0, double _rhoIx=0, double _rhoIy=0, double _rhoIz=0);

protected:
	//EDC int beamCrossSectionType;			//$EDC$[varaccess,EDCvarname="cross_section_type",EDCfolder="",tooltiptext="1: rectangular, 2: circular, 3: polygonal"]
	//EDC Vector beamCrossSectionSize;	//$EDC$[varaccess,EDCvarname="cross_section_size",EDCfolder="",variable_length_vector,tooltiptext="vector length of cross_section_size depends on cross_section_type: length 1 for circular cross section, length 2 for rectangular cross section (y and z extension), and length 2*n for polygonal cross section (p1y,p1z,p2y,p2z,...,pny,pnz)"];
  //EDC double beamEA;								//$EDC$[varaccess,EDCvarname="EA",EDCfolder="",tooltiptext="youngs modulus * area"]
	//EDC double beamEIy;								//$EDC$[varaccess,EDCvarname="EIy",EDCfolder="",tooltiptext="bending stiffness w.r.t. y-axis (2D-beam)"]
	//EDC double beamEIz;								//$EDC$[varaccess,EDCvarname="EIz",EDCfolder="",tooltiptext="bending stiffness w.r.t. z-axis"]
	//EDC double beamGAky;							//$EDC$[varaccess,EDCvarname="GAky",EDCfolder="",tooltiptext="shear stiffness including shear correction factor ky (2D-beam)"]
	//EDC double beamGAkz;							//$EDC$[varaccess,EDCvarname="GAkz",EDCfolder="",tooltiptext="shear stiffness including shear correction factor kz"]
	//EDC double beamGJkx;							//$EDC$[varaccess,EDCvarname="GJkx",EDCfolder="",tooltiptext="torsional stiffness including shear correction factor kx"]
	//EDC double beamRhoA;							//$EDC$[varaccess,EDCvarname="RhoA",EDCfolder="",tooltiptext="density * area"]
	//EDC double beamRhoIx;							//$EDC$[varaccess,EDCvarname="RhoIx",EDCfolder="",tooltiptext="density * second area of moment w.r.t. x-axis"]
	//EDC double beamRhoIy;							//$EDC$[varaccess,EDCvarname="RhoIy",EDCfolder="",tooltiptext="density * second area of moment w.r.t. y-axis (2D-beam)"]
	//EDC double beamRhoIz;							//$EDC$[varaccess,EDCvarname="RhoIz",EDCfolder="",tooltiptext="density * second area of moment w.r.t. z-axis"]

	// to be removed from edc of base class
	//EDC double density;								//$EDC$[varaccess,EDCvarname="density",remove,EDCfolder="Solid",tooltiptext="density (rho) for gravitational force"]
	//EDC double youngs_modulus;				//$EDC$[varaccess,EDCvarname="youngs_modulus",remove,EDCfolder="Solid",tooltiptext="Youngs modulus"]
	//EDC double poisson_ratio;					//$EDC$[varaccess,EDCvarname="poisson_ratio",remove,EDCfolder="Solid",tooltiptext="Poisson ratio"]
	//EDC int plane_stress;							//$EDC$[varaccess,EDCvarname="plane_stress",remove,EDCfolder="Solid",int_bool,tooltiptext="for 2D-Elements only; 1: plane stress, 0: plane strain"]
	//EDC double yieldstress;						//$EDC$[varaccess,EDCvarname="yield_stress",remove,EDCfolder="Inelasticity",tooltiptext="Yield Stress s_y, e.g., |dev s| <= s_y"]
	//EDC double tangentmodule;					//$EDC$[varaccess,EDCvarname="tangent_module",remove,EDCfolder="Inelasticity",tooltiptext="Modulus of hardening H"]
	//EDC int plane;										//$EDC$[varaccess,EDCvarname="plane",EDCfolder="Solid",remove]

	//EDC double density;								//$EDC$[varaccess,EDCvarname="density",EDCfolder="",tooltiptext="density (rho) for gravitational force"] // changed MSax 3013-02-21 because needed for gravity load for Beam3D
};
//$EDC$[endclass,Beam3DProperties]

//$EDC$[beginclass,classname=Beam2DProperties,parentclassname=Material,addelementtypename=Beam2DProperties]
class Beam2DProperties: public Material
{
public:
	Beam2DProperties() : Material () { InitConstructor(); }
	Beam2DProperties(MBS* mbs) : Material (mbs) { InitConstructor(); }
	Beam2DProperties(const Beam2DProperties& bp) : Material (bp) { CopyFrom(bp); }
	//virtual ~Beam2DProperties() { Material::~Material(); }
	virtual void CopyFrom(const Beam2DProperties& s) {}
	virtual Beam2DProperties* GetCopy() const
	{
		Beam2DProperties* ec = new Beam2DProperties(mbs);
		ec->CopyFrom(*this);
		return ec;
	}
	virtual mystr GetTypeName() const;
	virtual void InitConstructor();
	virtual void SetMaterialFromEDCData(); 

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual void SetBeam2DProperties(int _beamCrossSectionType, Vector _beamCrossSectionSize, double _rho, double _EA, double _EIy, double _GAky, double _rhoA, double _rhoIx, double _rhoIy);
	virtual void SetBeam2DProperties(double _rho, double _EA, double _EIy, double _GAky=0, double _rhoA=0, double _rhoIx=0, double _rhoIy=0);

protected:
	//EDC int beamCrossSectionType;			//$EDC$[varaccess,EDCvarname="cross_section_type",EDCfolder="",tooltiptext="1: rectangular, 2: circular, 3: polygonal"]
	//EDC Vector beamCrossSectionSize;	//$EDC$[varaccess,EDCvarname="cross_section_size",EDCfolder="",variable_length_vector,tooltiptext="vector length of cross_section_size depends on cross_section_type: length 1 for circular cross section, length 2 for rectangular cross section (y and z extension), and length 2*n for polygonal cross section (p1y,p1z,p2y,p2z,...,pny,pnz)"];
  //EDC double beamEA;								//$EDC$[varaccess,EDCvarname="EA",EDCfolder="",tooltiptext="youngs modulus * area"]
	//EDC double beamEIy;								//$EDC$[varaccess,EDCvarname="EIy",EDCfolder="",tooltiptext="bending stiffness w.r.t. y-axis (2D-beam)"]
	//EDC double beamGAky;							//$EDC$[varaccess,EDCvarname="GAky",EDCfolder="",tooltiptext="shear stiffness including shear correction factor ky (2D-beam)"]
	//EDC double beamGJkx;							//$EDC$[varaccess,EDCvarname="GJkx",EDCfolder="",tooltiptext="torsional stiffness including shear correction factor kx"]
	//EDC double beamRhoA;							//$EDC$[varaccess,EDCvarname="RhoA",EDCfolder="",tooltiptext="density * area"]
	//EDC double beamRhoIx;							//$EDC$[varaccess,EDCvarname="RhoIx",EDCfolder="",tooltiptext="density * second area of moment w.r.t. x-axis"]
	//EDC double beamRhoIy;							//$EDC$[varaccess,EDCvarname="RhoIy",EDCfolder="",tooltiptext="density * second area of moment w.r.t. y-axis (2D-beam)"]

	// to be removed from edc of base class
	//EDC double density;								//$EDC$[varaccess,EDCvarname="density",remove,EDCfolder="Solid",tooltiptext="density (rho) for gravitational force"]
	//EDC double youngs_modulus;				//$EDC$[varaccess,EDCvarname="youngs_modulus",remove,EDCfolder="Solid",tooltiptext="Youngs modulus"]
	//EDC double poisson_ratio;					//$EDC$[varaccess,EDCvarname="poisson_ratio",remove,EDCfolder="Solid",tooltiptext="Poisson ratio"]
	//EDC int plane_stress;							//$EDC$[varaccess,EDCvarname="plane_stress",remove,EDCfolder="Solid",int_bool,tooltiptext="for 2D-Elements only; 1: plane stress, 0: plane strain"]
	//EDC double yieldstress;						//$EDC$[varaccess,EDCvarname="yield_stress",remove,EDCfolder="Inelasticity",tooltiptext="Yield Stress s_y, e.g., |dev s| <= s_y"]
	//EDC double tangentmodule;					//$EDC$[varaccess,EDCvarname="tangent_module",remove,EDCfolder="Inelasticity",tooltiptext="Modulus of hardening H"]
	//EDC int inelasticity_type;				//$EDC$[varaccess,EDCvarname="inelasticity_type",remove,EDCfolder="Inelasticity",tooltiptext="0 = not specified, 1 = Prandtl Reuss plasticity + isotropic hardening"]
	//EDC int plane;										//$EDC$[varaccess,EDCvarname="plane",EDCfolder="Solid",remove]
};
//$EDC$[endclass,Beam2DProperties]


#endif