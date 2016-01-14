//#**************************************************************
//#
//# filename:             material.cpp
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

#include "mbs_interface.h"
#include "element.h"
#include "material.h"
#include "solversettings_auto.h"

const Matrix3D Material::dummyMatrix3D = Matrix3D(0);
const Matrix2D Material::dummyMatrix2D = Matrix2D(0);
const Vector Material::dummyVector = Vector(0);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Material   Material   Material   Material   Material   Material   
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Material::InitConstructor()
{
	//fill in default values
	materialname = GetTypeName();
	type = TMatEmpty;
	plane = 0;
	plane_stress = 0;

	materialcolor.Set(0, 0, 1); //default color - blue

	density = 0;
	dampingM = 0;
	dampingK = 0;

	youngs_modulus = 0;
	poisson_ratio = 0;

	beamCrossSectionType = 1;
	beamCrossSectionSize.SetLen(2);
	beamCrossSectionSize.SetAll(0);
	beamA = 0;
	beamRhoA = 0;
	beamRhoIx = 0;
	beamRhoIy = 0;
	beamRhoIz = 0;
	beamEA = 0;
	beamEIy = 0;
	beamEIz = 0;
	beamGAky = 0;
	beamGAkz = 0;
	beamGJkx = 0;

	yieldstress = 0;
	tangentmodule = 0;

	inelasticity_type = IT_LinearElastic;
	inelasticity_solution_method = ISM_Default;
}

void Material::CopyFrom(const Material& s)
{
	mbs = s.mbs;

	Material::InitConstructor(); // (AD) explicit call of class Init() function

	plane = s.plane;  
	plane_stress = s.plane_stress;

	materialname = s.materialname;
	materialnumber = s.materialnumber;
	type = s.type;

	materialcolor = s.materialcolor;

	density = s.density;
	dampingM = s.dampingM;
	dampingK = s.dampingK;

	youngs_modulus = s.youngs_modulus;
	poisson_ratio = s.poisson_ratio;

	beamCrossSectionType = s.beamCrossSectionType;
	beamCrossSectionSize = s.beamCrossSectionSize;
	beamA = s.beamA;
	beamRhoA = s.beamRhoA;
	beamRhoIx = s.beamRhoIx;
	beamRhoIy = s.beamRhoIy;
	beamRhoIz = s.beamRhoIz;
	beamEA = s.beamEA;
	beamEIy = s.beamEIy;
	beamEIz = s.beamEIz;
	beamGAky = s.beamGAky;
	beamGAkz = s.beamGAkz;
	beamGJkx = s.beamGJkx;

	yieldstress = s.yieldstress;
	tangentmodule = s.tangentmodule;

	//PG: merged inelasticity_type and nonlinmattype
	inelasticity_type = s.inelasticity_type;
	inelasticity_type_str = s.inelasticity_type_str;
	//nonlinmattype = s.nonlinmattype;  //$ PG 2013-7-31: deleted

	inelasticity_solution_method = s.inelasticity_solution_method;
	inelasticity_solution_method_str = s.inelasticity_solution_method_str;

	orthotropic_constants = s.orthotropic_constants;
}

void Material::GetElementData(ElementDataContainer& edc) 		//fill in all material data
{
	ElementData ed;
	ed.SetInt(GetMaterialNumber(), "Material_number"); ed.SetLocked(1); edc.Add(ed);

	//int nparam = NParameters();
	//for (int i=1; i <= nparam; i++)
	//{
	//	mystr name;
	//	double val;
	//	int flag;
	//	int group;
	//	GetParam(i, val, name, group, flag);

	//	int useit = 1; //check if this item should be included in saved file or in dialog
	//	if (IsPlanarMaterial() && (flag&64) || (!IsPlanarMaterial() && (flag&8))) useit = 0;
	//	if (!(type&TMatSpecial) && (flag&32)) useit = 0;
	//	if ((!(type&TMatBeam)) && (flag&16)) useit = 0;

	//	if (useit)
	//	{
	//		if (flag&1) //integer
	//		{
	//			ed.SetInt((int)val, name.c_str()); edc.Add(ed);
	//		}
	//		else if (flag&1) //check
	//		{
	//			ed.SetBool((int)val, name.c_str()); edc.Add(ed);
	//		}
	//		else //double value
	//		{
	//			ed.SetDouble(val, name.c_str()); edc.Add(ed);
	//		}
	//	}
	//}
	GetElementDataAuto(edc);
}

int Material::SetElementData(ElementDataContainer& edc) //set material data according to ElementDataContainer
{
	int rv = 	SetElementDataAuto(edc);
	SetMaterialFromEDCData();
	return rv;
}

void Material::SetMaterialFromEDCData()    //set material type and material parameters by user specified EDC variables
{
	if (plane)
	{
		type = TMBSMaterial(TMat2D+TMatFlexGeneral);
		if (!plane_stress)
		{
			type = (TMBSMaterial)(type + TMatPlaneStrain);
		}
	}
	else
	{
		type = TMatFlexGeneral;
	}
	inelasticity_type = GetInelasticityType(inelasticity_type_str);
	inelasticity_solution_method = GetInelasticitySolutionMethod(inelasticity_solution_method_str);
}

//$ PG 2013-1-21: Function Deprecated, Rigid should not use class Material!
void Material::SetMaterialRigid(double densityI)
{
	//$ PG 2012-12-19: Call of Init() should be removed: 
	//Init();

	type = TMatRigidBody;
	density = densityI;
}

void Material::SetMaterialCM(double densityI, double youngs_modulusI, double poisson_ratioI, Vector3D color, int planeI, int planestressI)
{
	//$ JG: changed because old code was confusing!
	
	//$ PG 2012-12-19: Call of Init() should be removed:
	//Init();
	
	density = densityI;
	youngs_modulus = youngs_modulusI;
	poisson_ratio = poisson_ratioI;

	inelasticity_type = IT_LinearElastic;
	inelasticity_type_str = "linear_elastic";
	
	materialcolor = color;
		
	if (planeI)
	{
		type = TMBSMaterial(TMat2D+TMatFlexGeneral);
		if (!planestressI)
		{
			type = (TMBSMaterial)(type + TMatPlaneStrain);
		}
	}
	else
	{
		type = TMatFlexGeneral;
	}
}

void Material::SetMaterialCMOrthotropic(double densityI, Vector3D color, const Material::OrthotropicConstants& oc)
{
	//$JG
	//$ PG 2012-12-19: Call of Init() should be removed:
	//Init();
	
	type = (TMBSMaterial)(TMatFlexGeneral + TMatOrthotropic);

	density = densityI;
	materialcolor = color;

	orthotropic_constants = oc;

	//GetMBS()->UO(UO_LVL_all) << "ortho.E1=" << orthotropic_constants.E1 << "\n";
	//GetMBS()->UO(UO_LVL_all) << "ortho.E2=" << orthotropic_constants.E2 << "\n";
	//GetMBS()->UO(UO_LVL_all) << "ortho.E3=" << orthotropic_constants.E3 << "\n";
	//GetMBS()->UO(UO_LVL_all) << "ortho.G12=" << orthotropic_constants.G12 << "\n";

	//GetMBS()->UO(UO_LVL_all) << "ortho.nu12=" << orthotropic_constants.nu12 << "\n";
	//GetMBS()->UO(UO_LVL_all) << "ortho.nu13=" << orthotropic_constants.nu13 << "\n";
	//GetMBS()->UO(UO_LVL_all) << "ortho.nu23=" << orthotropic_constants.nu23 << "\n";


	youngs_modulus = 2.1e11;		//$ PG 2012-12-19: Why???
	poisson_ratio = 0.3;

}

void Material::SetMaterialCMInelastic(double densityI, double youngs_modulusI, double poisson_ratioI,
																			Material::InelasticityType inelasticity_typeI, Material::InelasticitySolutionMethod inelasticity_solution_methodI, double yield_stressI, double tangent_moduleI,
																			Vector3D color, int planeI, int planestressI)
//$ PG 2013-7-31: standard set function for inelastic material of any kind (nonlinear elastic, elastoplastic, etc.)
{
	SetMaterialCM(densityI, youngs_modulusI, poisson_ratioI, color, planeI, planestressI);

	//$ PG 2013-7-31: TMatPlast deleted
	//SetType((TMBSMaterial)(GetType() + TMatPlast));

	SetInelasticityType(inelasticity_typeI);
			
	if (inelasticity_typeI == IT_ElastoPlastic)
	{
		yieldstress = yield_stressI;
		tangentmodule = tangent_moduleI;
	}

	SetInelasticitySolutionMethod(inelasticity_solution_methodI);
}

mystr Material::GetInelasticityTypeString(Material::InelasticityType it) const
{
	switch (it)
	{
	case IT_LinearElastic:
		return mystr("linear_elastic");
	case IT_ElastoPlastic:
		return mystr("elasto_plastic");
	case IT_NonlinElasticSimoHughes:
		return mystr("nonlinear_elastic_Simo_Hughes");
	default:
		return mystr("");
	}
}

Material::InelasticityType Material::GetInelasticityType(const mystr& it_str) const
{
	if (it_str.Compare(mystr("linear_elastic")))
		return IT_LinearElastic;
	else if (it_str.Compare(mystr("elasto_plastic")))
		return IT_ElastoPlastic;
	else if (it_str.Compare(mystr("nonlinear_elastic_Simo_Hughes")))
		return IT_NonlinElasticSimoHughes;
	else
		return IT_NotSpecified;
}

mystr Material::GetInelasticitySolutionMethodString(Material::InelasticitySolutionMethod ism) const
{
	switch (ism)
	{
	case ISM_Default:
		return mystr("default");
	case ISM_FixedPoint:
		return mystr("fixed_point");
	case ISM_ReturnMapping:
		return mystr("return_mapping");
	case ISM_ConsistentTangentStiffness:
		return mystr("consistent_tangent_stiffness");
	default:
		return mystr("");
	}
}

Material::InelasticitySolutionMethod Material::GetInelasticitySolutionMethod(const mystr& ism_str) const
{
	if (ism_str.Compare(mystr("default")))
		return ISM_Default;
	else if (ism_str.Compare(mystr("fixed_point")))
		return ISM_FixedPoint;
	else if (ism_str.Compare(mystr("return_mapping")))
		return ISM_ReturnMapping;
	else if (ism_str.Compare(mystr("consistent_tangent_stiffness")))
		return ISM_ConsistentTangentStiffness;
	else
		return ISM_NotSpecified;
}

//$ PG 2013-7-31: set function depricated!
void Material::SetMaterial(int nonlinmattypeI, double densityI, double youngs_modulusI, double poisson_ratioI)
{
	InitConstructor();
	type = (TMBSMaterial)(TMatFlexGeneral); //$ JG: changed, erased TMatNonlinElast

	//$ PG 2013-7-31: nonlinmattype deleted, is now part of inelasticity_type, where nonlinmattype == 1 matches inelasticity_type = IT_NonlinElasticSimoHughes
	if (nonlinmattypeI == 1)
	{
		SetInelasticityType(IT_NonlinElasticSimoHughes);
	}

	density = densityI;
	youngs_modulus = youngs_modulusI;
	poisson_ratio = poisson_ratioI;

	mbs->UO() << "WARNING in " << __FILE__ << ": void __thiscall Material::SetMaterial(...) deprecated - use SetMaterialCM or SetMaterialSM instead.\n";
}

int Material::CheckConsistency(mystr& errorstr)   //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = 0;

	return rv;
}







// ++++++++++++++++++++++++++++++++++++++
// ++++++++++  Beam3DProperties +++++++++
// ++++++++++++++++++++++++++++++++++++++

//PG, new version
void Beam3DProperties::SetBeam3DProperties(int _beamCrossSectionType /*1: rectangular, 2: circular, 3: polygonal*/, Vector _beamCrossSectionSize, double _rho, double _EA, double _EIy, double _EIz, double _GJkx, double _GAky, double _GAkz, double _rhoA, double _rhoIx, double _rhoIy, double _rhoIz)
{
	SetBeam3DProperties(_rho, _EA, _EIy, _EIz, _GJkx, _GAky, _GAkz, _rhoA, _rhoIx, _rhoIy, _rhoIz);

	beamCrossSectionType = _beamCrossSectionType; /*1: rectangular, 2: circular, 3: polygonal*/
	beamCrossSectionSize = _beamCrossSectionSize; //Vector of size 1 if circular, 2 if rectangular, and 2*n if polyognal (n... number of nodes: [p1y,p1z,p2y,p2z,...,pny,pnz])
}

//JG, depricated, because it does not initialize size
void Beam3DProperties::SetBeam3DProperties(double _rho, double _EA, double _EIy, double _EIz, double _GJkx, double _GAky, double _GAkz, double _rhoA, double _rhoIx, double _rhoIy, double _rhoIz)
{
	SetType(TMBSMaterial(TMatBeam));

	Density() = _rho;
	BeamEA() = _EA;
	BeamEIy() = _EIy;
	BeamEIz() = _EIz;
	BeamGJkx()= _GJkx;
	BeamGAky()= _GAky;
	BeamGAkz()= _GAkz;
	BeamRhoA() = _rhoA;
	BeamRhoIx() = _rhoIx;
	BeamRhoIy() = _rhoIy;
	BeamRhoIz() = _rhoIz;
}


// ++++++++++++++++++++++++++++++++++++++
// ++++++++++  Beam2DProperties +++++++++
// ++++++++++++++++++++++++++++++++++++++

void Beam2DProperties::SetBeam2DProperties(int _beamCrossSectionType, Vector _beamCrossSectionSize, double _rho, double _EA, double _EIy, double _GAky, double _rhoA, double _rhoIx, double _rhoIy)
{
	SetBeam2DProperties(_rho, _EA, _EIy, _GAky, _rhoA, _rhoIx, _rhoIy);

	beamCrossSectionType = _beamCrossSectionType; //1: rectangular, 2: circular, 3: polygonal
	beamCrossSectionSize = _beamCrossSectionSize; //Vector of size 1 if circular, 2 if rectangular, and 2*n if polyognal (n... number of nodes: [p1y,p1z,p2y,p2z,...,pny,pnz])
}

void Beam2DProperties::SetBeam2DProperties(double _rho, double _EA, double _EIy, double _GAky, double _rhoA, double _rhoIx, double _rhoIy)
{
	//$ PG 2012-12-19: Should be removed: Init();
	SetType(TMBSMaterial(TMatBeam+TMat2D));
	plane = 1;

	Density() = _rho;
	BeamEA() = _EA;
	BeamEIy() = _EIy;
	BeamGAky()= _GAky;
	BeamRhoA() = _rhoA;
	BeamRhoIx() = _rhoIx;
	BeamRhoIy() = _rhoIy;
}




























//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/// Computation

void Material::ComputeStressFromStrain(
									const Matrix3D& strain, Matrix3D& stress,
									const Vector& inel_variables,	const Matrix3D& initial_strain
									) const
{
	double Em = YoungsModulus();
	double nu = PoissonRatio();
	double la = Em * nu / ((1.+nu)*(1.-2.*nu));
	double mu = Em / 2. / (1.+nu);

	if (GetInelasticityType() == IT_LinearElastic)
	{ //material is linear elastic
		//orthotropic material just works with linear elasticity!!!!
		Matrix3D effective_strain = strain - initial_strain;

		if (!IsOrthotropicMaterial())
		{
			stress = ((2.*mu) * effective_strain + Matrix3D(la * effective_strain.Trace()));

		}
		else
		{


			ConstVector<6> stressv(6);
			ConstMatrix<36> C(6,6);

			ComputeElasticityMatrix(C);

			ConstVector<6> strainv(effective_strain(1,1),effective_strain(2,2),effective_strain(3,3),  
														 2*effective_strain(2,3), 2*effective_strain(1,3), 2*effective_strain(1,2));

			Mult(C, strainv, stressv);

			stress(1,1) = stressv(1); stress(1,2) = stressv(6); stress(1,3) = stressv(5);
			stress(2,1) = stressv(6); stress(2,2) = stressv(2); stress(2,3) = stressv(4);
			stress(3,1) = stressv(5); stress(3,2) = stressv(4); stress(3,3) = stressv(3);

			//if (once < 100)
			//{
			//	once++;
			//	Matrix3D stressold=stress;
			//  stress = ((2.*mu) * effective_strain + Matrix3D(la * effective_strain.Trace()));
			//	GetMBS()->UO(UO_LVL_all) << "stressdiff=" << stress-stressold << "\n";

			//	GetMBS()->UO(UO_LVL_all) << "C=" << C << "\n";
			//	GetMBS()->UO(UO_LVL_all) << "mu=" << mu << "\n";
			//	GetMBS()->UO(UO_LVL_all) << "la=" << la << "\n";
			//	GetMBS()->UO(UO_LVL_all) << "strain=" << strain << "\n";
			//	GetMBS()->UO(UO_LVL_all) << "effective strain=" << effective_strain << "\n";
			//	//GetMBS()->UO(UO_LVL_all) << stress << "\n";
			//}

		}

	}
	else if  (GetInelasticityType() == IT_NonlinElasticSimoHughes)
	{ //material is nonlinear elastic (something different from Hookes law), but not inelastic (eg., no strain splitting)

		Matrix3D C = 2.*strain; //C = 2E+I = F^T F
		C(1,1)+=1.; C(2,2)+=1.; C(3,3)+=1.;
		if (!C.Invert()) {GetMBS()->UO() << "ERROR: cannot invert strain tensor in Material\n";}

		//PG: F not needed in the interface
		// since J^2 = det(F)^2 = det(F^T F) = det(C)
		double J_square = C.Det();
		//double J = gradient_F.Det();

		stress = 0.5*la*(J_square - 1.)*C + mu*(Matrix3D(1.)-C); //C is inverse of C !!!!

	}
	else if (GetInelasticityType() == IT_ElastoPlastic)
	{ //material is inelastic (or plastic), but elastic law is linear (Hookes law)
		
		assert (inel_variables.Length() >= 6);

		Matrix3D inel_strain;
		StrainVectorToMatrix3D(inel_strain, inel_variables);

		if (GetInelasticitySolutionMethod() == ISM_ConsistentTangentStiffness || GetInelasticitySolutionMethod() == ISM_ReturnMapping)   //use the updated instead of the old inelastic strain
		{
			double hardening_param = inel_variables(7);

			Matrix3D trial_strain = strain - initial_strain - inel_strain;
			Matrix3D trial_stress = 2.*mu*trial_strain + la*trial_strain.Trace()*Matrix3D(1.);
			trial_stress += - trial_stress.Trace()/3.*Matrix3D(1.);    //deviator of trial stress
			double norm_trial_stress = sqrt(trial_stress.InnerProduct(trial_stress));
			double scaled_yield_stress = sqrt(2./3.)*YieldStress();

			double beta = norm_trial_stress - scaled_yield_stress*(1. + TangentModule()*hardening_param);
			if ( beta > 0 )   // new plasticity occurs
			{
				double xi = 2.*mu + pow(scaled_yield_stress*TangentModule(), 2);
				inel_strain = inel_strain + (beta/xi/norm_trial_stress)*trial_stress;
			}
		}

		//-----------------------------------------------------
		//EK 2012-06-25 if we have orthotropic + inel material: computation via matrix C
		Matrix3D s = strain - inel_strain - initial_strain;
		
		if (!IsOrthotropicMaterial())
			//old
			stress = ((2*mu) * s + Matrix3D(la * s.Trace()));
		else
		{
			ConstMatrix<36> C;
			this->ComputeElasticityMatrix(C);
			
			ConstVector<6> elastic_strain(s(1,1),s(2,2),s(3,3), 2*s(2,3), 2*s(1,3), 2*s(1,2));
			ConstVector<6> stress_vec(6);
			stress_vec = C*elastic_strain;

			stress(1,1) = stress_vec(1); stress(1,2) = stress_vec(6); stress(1,3) = stress_vec(5); 
			stress(2,1) = stress_vec(6); stress(2,2) = stress_vec(2); stress(2,3) = stress_vec(4); 
			stress(3,1) = stress_vec(5); stress(3,2) = stress_vec(4); stress(3,3) = stress_vec(3); 
		}
		//-----------------------------------------------------
	}
	else
	{
		assert(0 && "This inelasticity type has not been implemented yet.");
	}
}

double Material::ComputeYieldFunction(const Matrix3D& strain, const Vector& inelastic_variables) const
{
	double mu = GetMu();
	double lambda = GetLambda();
	Matrix3D inel_strain;
	StrainVectorToMatrix3D(inel_strain, inelastic_variables);
	Matrix3D trial_strain = strain - inel_strain;
	Matrix3D trial_stress = 2*mu*trial_strain + lambda*trial_strain.Trace()*Matrix3D(1.);
	trial_stress += (-trial_stress.Trace()/3.)*Matrix3D(1.);     //deviator of trial stress
	double norm_trial_stress = sqrt(trial_stress.InnerProduct(trial_stress));
	double scaled_yield_stress = sqrt(2./3.)*YieldStress();
	double hardening_param = inelastic_variables(7);
	double yield_function = sqrt(3./2.)*(norm_trial_stress - scaled_yield_stress*(1 + TangentModule()*hardening_param));
	
	return yield_function;
}

double Material::ComputeInelasticVariablesUpdateFromStrain(Vector& inel_variables_update, const Vector& inel_variables, const Matrix3D& strain) const
{ // calculate the actual state of inelastic variables as a function of the total strain
	if (GetInelasticityType() == IT_ElastoPlastic)
	{
		double inel_strain_update_norm = 0.;   // return value in case of ISM_ReturnMapping and ISM_FixedPoint
		double mu = GetMu();
		double lambda = GetLambda();
		Matrix3D trial_stress;

		Matrix3D inel_strain;
		StrainVectorToMatrix3D(inel_strain, inel_variables);
		double hardening_param = inel_variables(7);

		Matrix3D trial_strain = strain - inel_strain;
		trial_stress = 2*mu*trial_strain + lambda*trial_strain.Trace()*Matrix3D(1.);
		trial_stress += - (trial_stress.Trace()/3.)*Matrix3D(1.);    //deviator of trial stress
		double norm_trial_stress = sqrt(trial_stress.InnerProduct(trial_stress));
		double scaled_yield_stress = sqrt(2./3.)*YieldStress();

		double beta = norm_trial_stress - scaled_yield_stress*(1 + TangentModule()*hardening_param);
		if ( beta > 0. )   // new plasticity occurs
		{
			double xi = 2.*mu + pow(scaled_yield_stress*TangentModule(),2);
			Matrix3D inel_strain_update = (beta/xi/norm_trial_stress)*trial_stress;
			Matrix3DToStrainVector(inel_variables_update, inel_strain_update);

			inel_strain_update_norm = sqrt(inel_strain_update.InnerProduct(inel_strain_update));
			inel_variables_update(7) = scaled_yield_stress*TangentModule()*inel_strain_update_norm;
		}
		else
		{
			inel_variables_update.SetAll(0.);
		}

		switch (GetInelasticitySolutionMethod())
		{
		case ISM_FixedPoint:
			//return norm of inelastic strain update
			return inel_strain_update_norm;
		case ISM_ConsistentTangentStiffness:
		case ISM_ReturnMapping:
			//in this case return 0, since the update formula is substituted in the equations of motion (i.e., in EvalF2()), PostNewtonStep is just used for the final update of the old inelastic variables
			return 0.;
		default:
			return 1e99;
		}
	}
	else
	{
		assert (0 && "This type of inelasticity has not yet been implemented!\n");
	}
	return 0.;
}


void Material::ComputeStressFromStrain2D(
		const Matrix2D& strain, Matrix2D& stress,
		const Vector& inel_variables,	const Matrix2D& initial_strain
		) const
{
	double Em = YoungsModulus();
	double nu = PoissonRatio();
	double la = GetLambda();
	double mu = GetMu();
	if (IsPlanarMaterial() && !IsPlaneStrain()) //plane stress
	{
		//plane stress case can be implemented identical to the plane strain case, up to a reparameterization of the lame parameter lambda (la)
		la *= mu/(mu + .5*la);
		// hence, also Young's modulus 'Em' and Poisson's ratio 'nu' have to be adapted (however, Youngs modulus 'Em' ain't gonna be used hereafter)
		nu = .5*la/(la+mu);
	}


	if (GetInelasticityType() == IT_LinearElastic)
	{ //material is linear elastic
		Matrix2D effective_strain = strain - initial_strain;
	//EK 2012-04-25 if we have orthotropic with plane strain/stress -> computation via matrix C
		if (!IsOrthotropicMaterial())
			stress = ((2.*mu) * effective_strain + Matrix2D(la * effective_strain.Trace()));
		else
		{
			ConstMatrix<9> C;
			this->ComputeElasticityMatrix(C);
			Vector3D elastic_strain(effective_strain(1,1), effective_strain(2,2), 2*effective_strain(1,2));
			Vector3D stress_vec;
			stress_vec = C*elastic_strain;
			stress(1,1) = stress_vec(1); stress(1,2) = stress_vec(3); 
			stress(2,1) = stress_vec(3); stress(2,2) = stress_vec(2); 
		}
	}
	else if  (GetInelasticityType() == IT_NonlinElasticSimoHughes)
	{ //material is nonlinear elastic (something different from Hookes law), but not inelastic (eg., no strain splitting)

		Matrix2D C = 2.*strain; //C = 2E+I = F^T F
		C(1,1)+=1.; C(2,2)+=1.;
		if (!C.Invert()) {GetMBS()->UO() << "ERROR: cannot invert strain tensor in Material\n";}

		//PG: F not needed in the interface
		// since J^2 = det(F)^2 = det(F^T F) = det(C)
		double J_square = C.Det();
		//double J = gradient_F.Det();

		stress = 0.5*la*(J_square - 1.)*C + mu*(Matrix2D(1.)-C); //C is inverse of C !!!!
	}
	else if (GetInelasticityType() == IT_ElastoPlastic)
	{ 
		assert (inel_variables.Length() >= GetInelasticVariablesCount());

		Matrix2D inel_strain;
		StrainVectorToMatrix2D(inel_strain, inel_variables);
		
		Matrix2D inel_strain_update;
		if (GetInelasticitySolutionMethod() == ISM_ConsistentTangentStiffness || GetInelasticitySolutionMethod() == ISM_ReturnMapping)   //use the updated instead of the old inelastic strain
		{
			Vector3D el_strain_vec;      //  ( eps_11, eps_22, 2*eps_12 )     (engineering strain notation)  //$ PG 2013-11-25: make sure that only elastic and plastic stresses (and not strains) may be added
			el_strain_vec(1) = strain(1,1);
			el_strain_vec(2) = strain(2,2);
			el_strain_vec(3) = 2*strain(1,2);

			el_strain_vec(1) -= initial_strain(1,1);
			el_strain_vec(2) -= initial_strain(2,2);
			el_strain_vec(3) -= 2*initial_strain(1,2);

			Vector3D inel_strain_vec;    //  ( eps^p_11, eps^p_22, eps^p_12 )  (engineering stress notation), s.t. sigma^p = 2*mu*p
			inel_strain_vec(1) = inel_variables(1);
			inel_strain_vec(2) = inel_variables(2);
			inel_strain_vec(3) = 0.5*inel_variables(3);	
			
			double hardening_param = inel_variables(4);

			double la2mu = 2.*mu + la;
			Matrix3D elasticity_tensor(
				la2mu, la, 0,
				la, la2mu, 0,
				0,  0,    mu);

			// generate deviator, 
			// in case of plane strain, sigma_zz = nu(sigma_xx+sigma_yy), tr(sigma) = (1+nu)*(sigma_xx+sigma_yy), dev(sigma) = sigma-tr(sigma)/3 I computes as
			double dval = (nu+1.)/3.;
			Matrix3D deviator(
				1.-dval, -dval,   0,
				-dval,   1.-dval, 0,
				0,       0,       1);

			Vector3D trial_stress = deviator*elasticity_tensor*el_strain_vec - 2.*mu*inel_strain_vec;

			// compute norm of deviatoric matrix, where only x and y-components are provided
			// zz-component = -(xx+yy)
			Matrix3D norm_mat(2.);
			norm_mat(1,2) = 1.;
			norm_mat(2,1) = 1.;
			double norm_trial_stress = sqrt(trial_stress*norm_mat*trial_stress);
			double scaled_yield_stress = sqrt(2./3.)*YieldStress();
			double beta = norm_trial_stress - scaled_yield_stress*(1 + TangentModule()*hardening_param);

			if ( beta > 0 )   // new plasticity occurs
			{
				double xi = 2.*mu + pow(scaled_yield_stress*TangentModule(), 2);
				Vector3D inel_strain_update_vec = (sqrt(3./2.)*beta/xi/norm_trial_stress)*trial_stress;

				inel_strain_update(1,1) = inel_strain_update_vec(1);
				inel_strain_update(2,2) = inel_strain_update_vec(2);
				inel_strain_update(1,2) = inel_strain_update_vec(3);
				inel_strain_update(2,1) = inel_strain_update_vec(3);
			}
		}

		Matrix2D s = strain - initial_strain;
		inel_strain += inel_strain_update;

		if (!IsOrthotropicMaterial())
		{
			stress = ((2*mu) * (s - inel_strain) + Matrix2D(la * s.Trace()));
		}
		else
		{
			ConstMatrix<9> C;
			this->ComputeElasticityMatrix(C);
			
			s -= inel_strain;
			
			Vector3D elastic_strain(s(1,1),s(2,2), 2*s(1,2));
			Vector3D stress_vec;
			stress_vec = C*elastic_strain;
			stress(1,1) = stress_vec(1); stress(1,2) = stress_vec(3); 
			stress(2,1) = stress_vec(3); stress(2,2) = stress_vec(2); 
		}
	}
	else
	{
		//material is inelastic (or plastic) and elastic law is nonlinear (something different from Hookes law)
		assert(0 && "error: This inelasticity type has not yet been implemented.");
	}
}

double Material::ComputeYieldFunction2D(const Matrix2D& strain, const Vector& inelastic_variables) const
{
		double Em = YoungsModulus();
		double nu = PoissonRatio();
		double la = GetLambda();
		double mu = GetMu();
		if (!IsPlaneStrain()) //plane stress
		{
			//plane stress case can be implemented identical to the plane strain case, up to a reparameterization of the lame parameter lambda (la)
			la *= mu/(mu + .5*la);
			nu = .5*la/(la+mu);
		}

		Vector3D inel_strain;    //  ( comp_11, comp_22, comp_12 )
		inel_strain(1) = inelastic_variables(1);
		inel_strain(2) = inelastic_variables(2);
		inel_strain(3) = 0.5*inelastic_variables(3);

		Vector3D el_strain;      //  ( comp_11, comp_22, comp_12 )
		el_strain(1) = strain(1,1);
		el_strain(2) = strain(2,2);
		el_strain(3) = 2*strain(1,2);

		double la2mu = 2.*mu + la;
		Matrix3D elasticity_tensor(
			la2mu, la, 0,
			la, la2mu, 0,
			0,  0,    mu);

		// generate deviator, 
		// in case of plane strain, sigma_zz = nu(sigma_xx+sigma_yy), sigma_xx + sigma_yy + sigma_zz = (1+nu)*(sigma_xx+sigma_yy)
		double dval = (nu+1.)/3.;
		Matrix3D deviator(1.);
		deviator(1,1) -= dval;
		deviator(1,2) -= dval;
		deviator(2,1) -= dval;
		deviator(2,2) -= dval;

		Vector3D trial_stress = deviator*elasticity_tensor*el_strain - 2.*mu*inel_strain;

		// to compute norm of deviatoric matrix, where only x and y-components are provided
		// zz-component = -(xx+yy)
		Matrix3D norm_mat(2.);
		norm_mat(1,2) = 1.;
		norm_mat(2,1) = 1.;
		double norm_trial_stress = sqrt(trial_stress*norm_mat*trial_stress);
		double scaled_yield_stress = sqrt(2./3.)*YieldStress();
		double hardening_param = inelastic_variables(4);
		return sqrt(3./2.)*(norm_trial_stress - scaled_yield_stress*(1 + TangentModule()*hardening_param));
}

double Material::ComputeInelasticVariablesUpdateFromStrain2D(Vector& inel_variables_update, const Vector& inel_variables, const Matrix2D& strain) const
{ // calculate the actual state of inelastic variables as a function of the total strain in the planar case
	
	assert(IsPlanarMaterial());

	if (GetInelasticityType() == IT_ElastoPlastic)
	{ 
		double Em = YoungsModulus();
		double nu = PoissonRatio();
		double la = GetLambda(); 
		double mu = GetMu();
		if (!IsPlaneStrain()) //plane stress
		{
			//plane stress case can be implemented identical to the plane strain case, up to a reparameterization of the lame parameter lambda (la)
			la *= mu/(mu + .5*la);
			// hence, also Young's modulus 'Em' and Poisson's ratio 'nu' have to be adapted (however, Youngs modulus 'Em' ain't gonna be used hereafter)
			nu = .5*la/(la+mu);
		}

		Vector3D inel_strain;    //  ( eps^p_11, eps^p_22, eps^p_12 )  (engineering stress notation), s.t. sigma^p = 2*mu*p
		inel_strain(1) = inel_variables(1);
		inel_strain(2) = inel_variables(2);
		inel_strain(3) = 0.5*inel_variables(3);					
		double hardening_param = inel_variables(4);

		Vector3D el_strain;      //  ( eps_11, eps_22, 2*eps_12 )     (engineering strain notation)  //$ PG 2013-11-25: make sure that only elastic and plastic stresses (and not strains) may be added
		el_strain(1) = strain(1,1);
		el_strain(2) = strain(2,2);
		el_strain(3) = 2*strain(1,2);

		double la2mu = 2.*mu + la;
		Matrix3D elasticity_tensor(
			la2mu, la, 0,
			la, la2mu, 0,
			0,  0,    mu);

		// generate deviator, 
		// in case of plane strain, sigma_zz = nu(sigma_xx+sigma_yy), tr(sigma) = (1+nu)*(sigma_xx+sigma_yy), dev(sigma) = sigma-tr(sigma)/3 I computes as
		double dval = (nu+1.)/3.;
		Matrix3D deviator(1.);
		deviator(1,1) -= dval;
		deviator(1,2) -= dval;
		deviator(2,1) -= dval;
		deviator(2,2) -= dval;

		Vector3D trial_stress = deviator*elasticity_tensor*el_strain - 2.*mu*inel_strain;

		// compute norm of deviatoric matrix, where only x and y-components are provided
		// zz-component = -(xx+yy)
		Matrix3D norm_mat(2.);
		norm_mat(1,2) = 1.;
		norm_mat(2,1) = 1.;
		double norm_trial_stress = sqrt(trial_stress*norm_mat*trial_stress);

		double scaled_yield_stress = sqrt(2./3.)*YieldStress();
		double beta = norm_trial_stress - scaled_yield_stress*(1 + TangentModule()*hardening_param);
		if ( beta > 0 )   // new plasticity occurs
		{
			double xi = 2.*mu + pow(scaled_yield_stress*TangentModule(), 2);
			Vector3D inel_strain_update = (sqrt(3./2.)*beta/xi/norm_trial_stress)*trial_stress;
			inel_variables_update(1) = inel_strain_update(1);
			inel_variables_update(2) = inel_strain_update(2);
			inel_variables_update(3) = 2*inel_strain_update(3);
			
			// update hardening parameter
			double inel_strain_update_norm = sqrt(inel_strain_update*norm_mat*inel_strain_update);
			inel_variables_update(4) = scaled_yield_stress*TangentModule()*inel_strain_update_norm;
		}
		else
		{
			inel_variables_update.SetAll(0.);
		}
	}
	else
	{
		assert (0 && "This type of inelasticity has not yet been implemented!\n");
	}
	return 0.;
}

// AP Dez. 2012
// Material law due to return mapping algorithm by Simo and Hughes, 
// plane stress case differs from ComputeInelasticVariablesFromStrain2D substancially
void Material::ComputeInelasticVariablesFromStrainPlaneStress2D(Vector& inel_variables, const Matrix2D& strain)
{ // calculate the actual state of inelastic variables as a function of the total strain in the planar case
	
	assert(IsPlanarMaterial());

	if (GetInelasticityType() == IT_ElastoPlastic)
	{ 
		double Em = YoungsModulus();
		double nu = PoissonRatio();
		// standard Lame parameters
		double mu = Em / 2. / (1.+nu);
		double la = Em * nu / ((1.+nu)*(1.-2.*nu));
		// modified parameters for plane stress case
		la *= mu/(mu + .5*la);
		nu = .5*la/(la+mu);

		Vector3D inel_strain;    //  ( comp_11, comp_22, comp_12 )
		inel_strain(1) = inel_variables(1);
		inel_strain(2) = inel_variables(2);
		inel_strain(3) = inel_variables(3);					
		double hardening_param = inel_variables(4);
		double yield_function = inel_variables(5);

		Vector3D total_strain;      //  ( comp_11, comp_22, comp_12 )
		total_strain(1) = strain(1,1);
		total_strain(2) = strain(2,2);
		//AP Nov 2012: do not use engineering strains, but strain vector with eps_xy in third component
		//             now total_strain and inel_strain have the same form, and can be added, subtracted...
		total_strain(3) = strain(1,2);

		Matrix3D elasticity_tensor(2.*mu);
		elasticity_tensor(1,1) += la;
		elasticity_tensor(1,2) += la;
		elasticity_tensor(2,1) += la;
		elasticity_tensor(2,2) += la;
		//AP Nov 2012: do not use engineering strains, but strain vector with eps_xy in third component
		//             elasticity tensor has 2*mu on (3,3) position instead of standard 1*mu

		// Matrix P according to SimoHughes page 93, generate deviator, 
		// in case of plane stress, sigma_zz = 0
		Matrix3D P(2.,-1.,0.,
			-1.,2.,0.,
			0.,0.,6.);
		P *= 1./3.;

		// trial stress
		Vector3D trial_stress = elasticity_tensor*(total_strain - inel_strain);
		// deviator of trial stress
		Vector3D dev_trial_stress = P*trial_stress;

		// elementary computation:
		// ||dev sigma_0||^2 = (trial_stress^T*P*trial_stress)
		double sqrnorm_dev_trial_stress = dev_trial_stress*trial_stress;
		double norm_dev_trial_stress = sqrt(sqrnorm_dev_trial_stress);

		// value of (sqared) yield function:   ||dev sigma_0||^2 - 2/3(sigma_y + K hardening_param)^2
		// K .. plastic modulus
		double K = (YoungsModulus()*TangentModule())/(YoungsModulus()-TangentModule());
		double beta = sqrnorm_dev_trial_stress - 2./3.*sqr(YieldStress()+K*hardening_param);
		if ( beta > 0 )   // new plasticity occurs
		{
			// solve quadratic equation for gamma
			// update = gamma * dev sigma_0 / ||dev sigma_0||
			Vector3D C_dev_trial_stress = elasticity_tensor*dev_trial_stress;

			double A = (C_dev_trial_stress*P*C_dev_trial_stress)/sqrnorm_dev_trial_stress - 4./9.*K*K;
			double B = -2/norm_dev_trial_stress*(dev_trial_stress*C_dev_trial_stress) - 4./3.*(YieldStress()+K*hardening_param)*sqrt(2./3.)*K;
			double C = beta;
			double gamma = (-B - sqrt(B*B-4*A*C))/(2*A);
			Vector3D inel_strain_update = (gamma/norm_dev_trial_stress)*dev_trial_stress;
			inel_strain += inel_strain_update;


			inel_variables(1) = inel_strain(1);
			inel_variables(2) = inel_strain(2);
			inel_variables(3) = inel_strain(3);

			hardening_param += sqrt(2./3.)*gamma;
			inel_variables(4) = hardening_param;
			// update the trial stress
			trial_stress = elasticity_tensor*(total_strain - inel_strain);
			dev_trial_stress = P*trial_stress;
		}
		//yield_function (as additional 5th inelastic variable)
		norm_dev_trial_stress = sqrt(dev_trial_stress*trial_stress);
		inel_variables(5) = norm_dev_trial_stress - sqrt(2./3.)*(YieldStress()+K*hardening_param);
	}
	else
	{
		assert (0 && "This type of inelasticity has not yet been implemented!\n");
	}
}

// AP Dez. 2012
// Material law due to return mapping algorithm by Simo and Hughes, plane stress and strain case
// plane strain case differs from ComputeInelasticVariablesFromStrain2D by definition of TangentModule
// here, TangentModule is elastoplastic tangent module
void Material::ComputeInelasticVariablesFromStrainPlaneStrain2D(Vector& inel_variables, const Matrix2D& strain)
{ // calculate the actual state of inelastic variables as a function of the total strain in the planar case
	
	assert(IsPlanarMaterial());
	
	if (GetInelasticityType() == IT_LinearElastic || GetInelasticityType() == IT_NonlinElasticSimoHughes)
	{
		return;
	}
	else if (GetInelasticityType() == IT_ElastoPlastic)
	{ 
		double Em = YoungsModulus();
		double nu = PoissonRatio();
		// standard Lame parameters
		double mu = Em / 2. / (1.+nu);
		double la = Em * nu / ((1.+nu)*(1.-2.*nu));

		Vector3D inel_strain;    //  ( comp_11, comp_22, comp_12 )
		inel_strain(1) = inel_variables(1);
		inel_strain(2) = inel_variables(2);
		inel_strain(3) = inel_variables(3);					
		double hardening_param = inel_variables(4);
		double yield_function = inel_variables(5);

		Vector3D total_strain;      //  ( comp_11, comp_22, comp_12 )
		total_strain(1) = strain(1,1);
		total_strain(2) = strain(2,2);
		//AP Nov 2012: do not use engineering strains, but strain vector with eps_xy in third component
		//             now total_strain and inel_strain have the same form, and can be added, subtracted...
		total_strain(3) = strain(1,2);

		Matrix3D elasticity_tensor(2.*mu);
		elasticity_tensor(1,1) += la;
		elasticity_tensor(1,2) += la;
		elasticity_tensor(2,1) += la;
		elasticity_tensor(2,2) += la;
		//AP Nov 2012: do not use engineering strains, but strain vector with eps_xy in third component
		//             modify elasticity tensor accordingly, 2*mu on (3,3) position instead of 1*mu
		//elasticity_tensor(3,3) -= mu;

		// Matrix P which generates deviator of stress sigma = C epsilon in case of plane strain, 
		// in case of plane stress, sigma_zz = 0
		Matrix3D P(1.-(1.+nu)/3.,-(1.+nu)/3.,0.,
			-(1.+nu)/3.,1.-(1.+nu)/3.,0.,
			0.,0.,1);

		// deviator of trial stress
		Vector3D dev_trial_stress = P*elasticity_tensor*total_strain - 2*mu*inel_strain;

		// matrix which represents norm for deviatoric stress vector
		Matrix3D norm(2.,1.,0.,
			1.,2.,0.,
			0.,0.,2);
		double sqrnorm_dev_trial_stress = dev_trial_stress*norm*dev_trial_stress;
		double norm_dev_trial_stress = sqrt(sqrnorm_dev_trial_stress);

		// value of yield function:   ||dev sigma_0|| - sqrt(2/3)(sigma_y + K hardening_param)
		// K .. plastic modulus
		double K = (YoungsModulus()*TangentModule())/(YoungsModulus()-TangentModule());
		double beta = norm_dev_trial_stress - sqrt(2./3)*(YieldStress()+K*hardening_param);
		if ( beta > 0 )   // new plasticity occurs
		{
			// solve linear update equation
			double gamma = beta/(2*mu+2./3.*K);
			Vector3D inel_strain_update = (gamma/norm_dev_trial_stress)*dev_trial_stress;
			inel_strain += inel_strain_update;

			inel_variables(1) = inel_strain(1);
			inel_variables(2) = inel_strain(2);
			inel_variables(3) = inel_strain(3);

			hardening_param += sqrt(2./3.)*gamma;
			inel_variables(4) = hardening_param;
			// different stress updates for plane strain/plane stress
			dev_trial_stress -= 2*mu*inel_strain_update;
		}
		//yield_function (as additional 5th inelastic variable)
		norm_dev_trial_stress = sqrt(dev_trial_stress*norm*dev_trial_stress);
		inel_variables(5) = norm_dev_trial_stress - sqrt(2./3.)*(YieldStress()+K*hardening_param);
	}
	else
	{
		assert (0 && "This type of inelasticity has not yet been implemented!\n");
	}
}



void Material::ComputeInelasticVariablesFromStrain1D(Vector& inel_variables, const double strain)
{ // calculate the actual state of inelastic variables as a function of the total strain in the planar case

	assert(IsPlanarMaterial());

	if (GetInelasticityType() == IT_LinearElastic || GetInelasticityType() == IT_NonlinElasticSimoHughes)
	{
		return;
	}
	else if (GetInelasticityType() == IT_ElastoPlastic)
	{ 
		double Em = YoungsModulus();
		double inel_strain;   
		inel_strain = inel_variables(1);
		double hardening_param = inel_variables(2);
		double yield_function = inel_variables(3);

		double el_strain;    
		el_strain = strain;

		double trial_stress = Em*(el_strain - inel_strain);

		double norm_trial_stress = fabs(trial_stress);

		double scaled_yield_stress = YieldStress();
		double beta = norm_trial_stress - scaled_yield_stress*(1 + TangentModule()*hardening_param);
		if ( beta > 0 )   // new plasticity occurs
		{
			double xi = Em + pow(scaled_yield_stress*TangentModule(), 2);
			double inel_strain_update = (beta/xi/norm_trial_stress)*trial_stress;
			inel_strain += inel_strain_update;

			inel_variables(1) = inel_strain;

			double inel_strain_update_norm = fabs(inel_strain_update);
			hardening_param += scaled_yield_stress*TangentModule()*inel_strain_update_norm;
			inel_variables(2) = hardening_param;
			trial_stress += - Em * inel_strain_update;
		}
		//yield_function (as additional 5th inelastic variable)
		norm_trial_stress = fabs(trial_stress);
		inel_variables(3) = norm_trial_stress - scaled_yield_stress*(1 + TangentModule()*hardening_param);
	}
	else
	{
		assert (0 && "This type of inelasticity has not yet been implemented!\n");
	}
}


//PG: This routine is deprecated, since it addresses geometrical properties only.
//*   It should be provided by the element - not by the material.
void Material::ComputeStrain(Matrix3D& strain, const Matrix3D& grad_u, int isgeomnonlin)
{
	if (isgeomnonlin)
	{
		Matrix3D Fh = grad_u + Matrix3D(1);
		strain = 0.5*(Fh.GetTp()*Fh-Matrix3D(1));
	}
	else
	{
		strain = 0.5*(grad_u.GetTp() + grad_u);
	}
}

void Material::ComputeStrain2D(Matrix2D& strain, const Matrix2D& grad_u, int isgeomnonlin)
{
	if (isgeomnonlin)
	{
		Matrix2D Fh = grad_u + Matrix2D(1);
		strain = 0.5*(Fh.GetTp()*Fh-Matrix2D(1));
	}
	else
	{
		strain = 0.5*(grad_u.GetTp() + grad_u);
	}
}

//$ PG 2013-8-2: replaced by DoAveragingStep/DoAveragingStep2D
//// Computes value of stress component zz averaged over element
//double Material::ComputeAveragedContactStress(FiniteElement3D* const fe)
//{
//	return 0;
//}

void Material::GetParam(int i, double& value, mystr& name, int& group, int& flag)
{
	//flag&1  ==integer
	//flag&2  ==positive with zero
	//flag&4  ==check (0 or 1)
	//flag&8  ==for plane problems only
	//flag&16 ==for beam problems only
	//flag&32 ==for special materials only
	//flag&64 ==for 3D problems only

	switch(i)
	{
	case 1: {value = density; name = "density"; group = 1; flag = 2; break;}
	case 2: {value = dampingM; name = "mass_prop_damping"; group = 1; flag = 2; break;}
	case 3: {value = dampingK; name = "stiffness_prop_damping"; group = 1; flag = 2; break;}

	case 4: {value = youngs_modulus; name = "youngs_modulus"; group = 2; flag = 2; break;}
	case 5: {value = poisson_ratio; name = "poisson_ratio"; group = 2; flag = 0; break;}
	case 6: {value = !IsPlaneStrain(); name = "planestress"; group = 3; flag = 8+4; break;}

	case 7: {value = beamA   ; name = "beamA";    group = 4; flag = 2+16; break;}
	case 8: {value = beamEIy ; name = "beamEIy";  group = 4; flag = 2+16; break;}
	case 9: {value = beamEIz ; name = "beamEIz";  group = 4; flag = 2+16+64; break;}
	case 10:{value = beamGAky; name = "beamGAky"; group = 4; flag = 2+16; break;}
	case 11:{value = beamGAkz; name = "beamGAkz"; group = 4; flag = 2+16+64; break;}
	case 12:{value = beamGJkx; name = "beamGJkx"; group = 4; flag = 2+16+64; break;}
	case 13:{value = beamEA  ; name = "beamEA";   group = 4; flag = 2+16; break;}

	case 14:{value = beamRhoA; name = "beamRhoA";   group = 4; flag = 2+16; break;}
	case 15:{value = beamRhoIx; name = "beamRhoIx"; group = 4; flag = 2+16+64; break;}
	case 16:{value = beamRhoIy; name = "beamRhoIy"; group = 4; flag = 2+16; break;}
	case 17:{value = beamRhoIz; name = "beamRhoIz"; group = 4; flag = 2+16+64; break;}



	case 21: {value = (int)inelasticity_type;     name = "inelasticity_type";     group = 3; flag = 1; break;}
	default:
		value=0; name=""; group=0; flag=0;
		GetMBS()->UO() << "Error in Material, parameter not found!\n";
	}
}

void Material::ComputeElasticityMatrix(Matrix& C) const 
{
	//assert ( (C.Getcols() == 6 && !IsPanarMaterial())||(C.Getcols() == 3 && IsPanarMaterial()) );

	double Em = YoungsModulus();
	double nu = PoissonRatio();

	if (!IsPlanarMaterial())
	{	//3D
		if (!IsOrthotropicMaterial())
		{
			C.SetSize(6,6);
			C.SetAll(0.);
			C(1,1) = 1.-nu;	C(1,2) = nu;	  C(1,3) = nu;
			C(2,1) = nu;	  C(2,2) = 1.-nu;	C(2,3) = nu;
			C(3,1) = nu;	  C(3,2) = nu;	  C(3,3) = 1.-nu;
			C(4,4) = (1.-2.*nu)/2.;
			C(5,5) = (1.-2.*nu)/2.;
			C(6,6) =(1.-2.*nu)/2.;
			C *= Em/((1.+nu)*(1.-2.*nu));
		}
		else
		{
			C.SetSize(6,6);
			C.SetAll(0.);
			const OrthotropicConstants& oc = orthotropic_constants;


			C(1,1) = (-oc.E2()+Sqr(oc.NU23())*oc.E3())*Sqr(oc.E1());    C(1,2) = -(oc.NU12()*oc.E2()+oc.NU13()*oc.NU23()*oc.E3())*oc.E1()*oc.E2();	  C(1,3) = -(oc.NU12()*oc.NU23()+oc.NU13())*oc.E1()*oc.E2()*oc.E3();
			C(2,1) = -(oc.NU12()*oc.E2()+oc.NU13()*oc.NU23()*oc.E3())*oc.E1()*oc.E2();	  C(2,2) = -(oc.E1()-Sqr(oc.NU13())*oc.E3())*Sqr(oc.E2());    C(2,3) = -(oc.NU23()*oc.E1()+oc.NU12()*oc.NU13()*oc.E2())*oc.E2()*oc.E3();
			C(3,1) = -(oc.NU12()*oc.NU23()+oc.NU13())*oc.E1()*oc.E2()*oc.E3();	  C(3,2) = -(oc.NU23()*oc.E1()+oc.NU12()*oc.NU13()*oc.E2())*oc.E2()*oc.E3();	  C(3,3) = -(oc.E1()-Sqr(oc.NU12())*oc.E2())*oc.E2()*oc.E3();
			C *= 1./(-oc.E1()*oc.E2()+oc.E1()*Sqr(oc.NU23())*oc.E3()+Sqr(oc.NU12())*Sqr(oc.E2())+2*oc.NU12()*oc.E2()*oc.NU13()*oc.NU23()*oc.E3()+Sqr(oc.NU13())*oc.E2()*oc.E3());


			C(4,4) = oc.G23();
			C(5,5) = oc.G13();
			C(6,6) = oc.G12();

		}
	}
	else
	{
		C.SetSize(3,3);
		if (!IsOrthotropicMaterial())
		{	
			if (!IsPlaneStrain())
			{ //plane stress
				double k = Em/(1.-nu*nu);
				C(1,1) = k;    C(1,2) = k*nu; C(1,3) = 0;
				C(2,1) = k*nu; C(2,2) = k; 		C(2,3) = 0;
				C(3,1) = 0; 	 C(3,2) = 0;  	C(3,3) = k * (1.-nu)/2.;
			}
			else
			{ //plane strain
				double k = Em/((1.+nu)*(1.-2*nu));
				C(1,1) = k*(1-nu); C(1,2) = k*nu;	    C(1,3) = 0;
				C(2,1) = k*nu;		 C(2,2) = k*(1-nu);	C(2,3) = 0;
				C(3,1) = 0; 			 C(3,2) = 0;				C(3,3) = k * (1.-2*nu)/2.;
			}
		}
		else
		{ 
			if (IsPlaneStrain())
			{
				//EK 2012-01-24 : plane strain case of orthotropic material 
				const OrthotropicConstants& oc = orthotropic_constants;
				C.SetAll(0.);
				double oce2 = oc.E2();
				double nu12 = oc.NU12();

				C(1,1) = oc.E1() * oce2 * (oce2 * Sqr(oc.NU23())-oc.E3());    C(1,2) = -oc.E1() * oce2 * (oc.E3() * nu12 + oce2 * oc.NU13() * oc.NU23());	 
				C(2,1) = -oc.E1() * oce2 *(nu12 * oc.E3() + oce2 * oc.NU13() * oc.NU23());	  C(2,2) = Sqr(oce2) * ( oc.E1() * Sqr(oc.NU13())-oc.E3());    
				C *= 1./(-oce2 * oc.E3()+Sqr(oc.NU23()) * Sqr(oce2)+oce2*Sqr(oc.NU13())*oc.E1()+oc.E1()*Sqr(nu12)*oc.E3()+2*oc.E1()*nu12*oc.NU13()*oc.NU23()*oce2);
				
				C(3,3) = oc.G12();
			}
		}
	}
}

void Material::ComputeConsistentTangentStiffnessMatrix2D(Matrix& C, const Matrix2D& strain, const Vector& inelastic_variables) const
{
	assert(0.); //$ PG 2013-8-14: to be implemented
}

// computes the consistent tangent stiffness matrix C (I - d actual_plastic_strain/d strain)
// C [at input]: elasticity tensor C,
// C [at output]: consistent tangent stiffness matrix C (I - d actual_plastic_strain/d strain)
// strain: (eps_xx, eps_yy, eps_zz, 2 eps_yz, 2 eps_zx, 2 eps_xy)
// inelastic_variables: (eps^p_xx, eps^p_yy, eps^p_zz, 2 eps^p_yz, 2 eps^p_zx, 2 eps^p_xy, hardening_parameter)
void Material::ComputeConsistentTangentStiffnessMatrix(Matrix& C, const Vector& strain, const Vector& inelastic_variables) const
{
	// the consisten tangent stiffness matrix is returned in C, depending on strain and inelastic varibles (in integration point)

	assert(!IsPlanarMaterial() && "use ComputeTangentStiffnessMatrix2D for 2D case");

	if (GetInelasticityType() == IT_ElastoPlastic && GetInelasticitySolutionMethod() == ISM_ConsistentTangentStiffness)
	{
		ConstMatrix<6*6> K(6,6), N(6,6);   //initialized with 0
		for (int i=1; i<=3; i++)
		{
			int i3 = i+3;
			N(i,i) = 1;
			N(i3,i3) = 2;
			K(i,i) = 1;
			K(i3,i3) = 1;
			for (int j=1; j<=3; j++)
			{
				K(i,j) -= 1./3.;
			}
		}
				
		ConstVector<6> inelastic_strain(6,1);
		for (int i=1; i<=3; i++)
		{
			int i3 = i+3;
			inelastic_strain(i) = inelastic_variables(i);
			inelastic_strain(i3) = 0.5*inelastic_variables(i3);
		}
		double hardening_parameter = inelastic_variables(7);

		double mu2 = 2*GetMu();
		ConstVector<6> dev_trial_stress = (ConstVector<6>&) (K*(C*strain) - mu2*inelastic_strain);

		double norm_dev_trial_stress_sq = dev_trial_stress*(N*dev_trial_stress);
		double norm_dev_trial_stress = sqrt(norm_dev_trial_stress_sq);
		double scaled_yield_stress = sqrt(2./3.)*YieldStress();
		double yield_function_of_trial_stress = norm_dev_trial_stress - scaled_yield_stress*(1 + hardening_parameter*TangentModule());

		if (yield_function_of_trial_stress > 0)
		{
			dev_trial_stress *= 1./norm_dev_trial_stress;
			double xi =  sqrt(3./2.)*mu2/(mu2 + pow(scaled_yield_stress*TangentModule(),2));
			double beta = yield_function_of_trial_stress/norm_dev_trial_stress;

			ConstVector<6> N_times_dev_trial_stress((ConstVector<6>&)(N*dev_trial_stress));
			dev_trial_stress *= xi*(1-beta);
			ConstMatrix<6*6> dPdE(6,6);
			for (int i=1; i<=6; i++)
			{
				for (int j=1; j<=6; j++)
				{
					dPdE(i,j) = dev_trial_stress(i) * N_times_dev_trial_stress(j);
				}
			}
			double beta_xi = beta*xi;
			for (int i=1; i<=6; i++)
			{
				dPdE(i,i) += beta_xi;
			}
			C -= dPdE*K*C; //C = Celast - 2*mu*dPdE   (2*mu was already muptiplied in xi)
		}
	}
}

void Material::GetLinearizedStress2D(const Vector3D& strain, Vector3D& stress, const Vector3D& inel_strain, const Vector3D& stress_lin) const
{
	ConstMatrix<36> C(6,6);
	ComputeElasticityMatrix(C);

	stress = C*(strain);
}

void Material::GetStress2D(const Vector3D& strain, Vector3D& stress, const Vector3D& inel_strain, const double contactstress) const
{
	ConstMatrix<36> C(6,6);
	ComputeElasticityMatrix(C);

	stress = C*(strain);
}

// new plastic strain, internal variable are computed
// new stress is computed
// input: from actual iteration step: eps, epsplast, internalvar, 
// input: from last time step: epsplast_old, internalvar_old
double Material::DoNonlinStepBeam2D(double& stress, double eps, double& epsplast, double& internalvar, double epsplast_old, double internalvar_old, int use_tangent_stiffness) const
{
	// discontinuous accuracy
	double da = GetMBS()->GetSolSet().discontinuousaccuracy; //JG: old DOption8

	if (yieldstress < 0 ) {GetMBS()->UO() << "Material:: Warning -- yieldstress " << yieldstress << " < 0 for inelastic material!\n"; return 0;}
	if (TangentModule() < 0) {GetMBS()->UO() << "Material:: Warning -- Tangent module " << TangentModule() << " < 0 for inelastic material, softening?!\n"; return 0;}

	stress = YoungsModulus()*(eps - epsplast);
	double depsplast = 0;//epsplast - epsplast_old;
	double dinternalvar = 0;//internalvar - internalvar_old;
	double depsplastold = depsplast;
	double smin, smax, help, yieldfun;
	double error = 0;

	if (!IsHardeningMaterial()) // perfect plasticity
	{
		smax = yieldstress;
		smin = -yieldstress;
		help = 0;
		yieldfun = fabs(stress) - yieldstress;
		if (stress > smax ) 
		{
			depsplast = yieldfun / YoungsModulus();
			error = yieldfun / YoungsModulus();
		}
		else if (stress < smin ) 
		{
			depsplast = -yieldfun / YoungsModulus();
			error = yieldfun / YoungsModulus();
		}
		else if (yieldfun < -da/**YoungsModulus()*/ && fabs(depsplast) > da /YoungsModulus() && !use_tangent_stiffness) 
		{
			error = fabs(depsplast);
			depsplast = epsplast_old - epsplast;
		}
		dinternalvar = 0;
	}
	else // isotropic hardening
	{
		help = YoungsModulus()*TangentModule()/(YoungsModulus()-TangentModule());
		smax = yieldstress + help*internalvar;
		smin = -yieldstress - help*internalvar;
		yieldfun = fabs(stress) - smax;
		if (stress > smax ) 
		{
			dinternalvar = yieldfun/(YoungsModulus()+help);
			depsplast = dinternalvar;
			error = yieldfun / YoungsModulus();
		}
		else if (stress < smin ) 
		{
			dinternalvar = yieldfun/(YoungsModulus()+help);
			depsplast = -dinternalvar;
			error = yieldfun / YoungsModulus();
		}
		else if (yieldfun < -da/*YoungsModulus()*/ && fabs(depsplast) > da/YoungsModulus() && !use_tangent_stiffness ) 
		{
			error = fabs(depsplast);
			depsplast = 1*(epsplast_old - epsplast);
			dinternalvar = 1*(internalvar_old - internalvar);
		}
		internalvar = internalvar + dinternalvar;
	}
	epsplast = epsplast + depsplast;

	return error;
}

int Material::IsPlanarMaterial() const
{
	if ((type&TMat2D) /*|| (type&TMatBeam2D)*/)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

int Material::IsPlaneStrain() const
{
	if (type&TMatPlaneStrain)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}




bool Material::IsMaterialOfBeamWithRectangularCrossSection() const 
{
	return IsMaterialOfBeam() && beamCrossSectionType == 1;
}

bool Material::IsMaterialOfBeamWithCircularCrossSection() const 
{
	return IsMaterialOfBeam() && beamCrossSectionType == 2;
}

bool Material::IsMaterialOfBeam() const 
{
	return type & TMatBeam; 
}


////////////////////////////////////////////
// class Beam3DProperties : public Material
////////////////////////////////////////////

mystr Beam3DProperties::GetTypeName() const 
{
	return "Beam3DProperties";
}

void Beam3DProperties::InitConstructor()
{
	materialname = GetTypeName();
	type = TMatBeam;

	beamCrossSectionType = 1;
	beamCrossSectionSize.SetLen(2);
	beamCrossSectionSize.SetAll(0);
	beamA = 0;
	beamRhoA = 0;
	beamRhoIx = 0;
	beamRhoIy = 0;
	beamRhoIz = 0;
	beamEA = 0;
	beamEIy = 0;
	beamEIz = 0;
	beamGAky = 0;
	beamGAkz = 0;
	beamGJkx = 0;
}

void Beam3DProperties::GetElementData(ElementDataContainer& edc)
{
	Material::GetElementData(edc);
}

int Beam3DProperties::SetElementData(ElementDataContainer& edc)
{
	int rv = Material::SetElementData(edc);
	return rv;
}

void Beam3DProperties::SetMaterialFromEDCData()    //set material type and material parameters by user specified EDC variables
{
	// do not change anything
	SetType((TMBSMaterial)(TMatBeam));
	inelasticity_type = GetInelasticityType(inelasticity_type_str);
	inelasticity_solution_method = GetInelasticitySolutionMethod(inelasticity_solution_method_str);
}

////////////////////////////////////////////
// class Beam2DProperties : public Material
////////////////////////////////////////////

mystr Beam2DProperties::GetTypeName() const 
{
	return "Beam2DProperties";
}

void Beam2DProperties::InitConstructor()
{
	materialname = GetTypeName();
	type = TMatBeam;

	beamCrossSectionType = 1;
	beamCrossSectionSize.SetLen(2);
	beamCrossSectionSize.SetAll(0);
	beamA = 0;
	beamRhoA = 0;
	beamRhoIx = 0;
	beamRhoIy = 0;
	beamRhoIz = 0;
	beamEA = 0;
	beamEIy = 0;
	beamEIz = 0;
	beamGAky = 0;
	beamGAkz = 0;
	beamGJkx = 0;
}

void Beam2DProperties::GetElementData(ElementDataContainer& edc)
{
	Material::GetElementData(edc);
}

int Beam2DProperties::SetElementData(ElementDataContainer& edc)
{
	int rv = Material::SetElementData(edc);
	return rv;
}

void Beam2DProperties::SetMaterialFromEDCData()    //set material type and material parameters by user specified EDC variables
{
	// do not change anything
	SetType((TMBSMaterial)(TMatBeam+TMat2D));
	inelasticity_type = GetInelasticityType(inelasticity_type_str);
	inelasticity_solution_method = GetInelasticitySolutionMethod(inelasticity_solution_method_str);
}

