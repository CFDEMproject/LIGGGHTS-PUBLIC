//#**************************************************************
//#
//# filename:             femath.h
//#
//# author:               Gerstmayr Johannes
//#
//# generated:            15.05.97
//# description:          Classes for linear and nonlinear algebra which is
//#												thought to be used in finite element or similar applications
//#												There are 2D and 3D and arbitrary size Vectors (Vector2D, Vector3D, Vector),
//#												arbitrary size matrices (Matrix), and a nonlinear system (NumNLSys)
//#												and a nonlinear solver (NumNLSolver)
//# remarks:							Indizes run from 1 to n except in Vector3D/2D
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
 
#ifndef MFEMATH__H
#define MFEMATH__H

#include <math.h>
#include "release_assert.h"
#include "solversettings_auto.h"

/////////////////////////
//turns off critical asserts (faster):
//PG: Here we define the strategy for asserting femath code.
//We still have to discuss if this is the optimal setting.
//Meanwhile, I suggest: it's best to assert in debug-mode
//and if the preprocessor flag "__ASSERT_IN_RELEASE_MODE__"
//is set (e.g., in "MBSKernelLib\preprocessor_includes.h").
//This is due to the fact, that some of the models may
//be executed in release-mode only.
#ifndef _DEBUG
#ifndef __ASSERT_IN_RELEASE_MODE__
#define __QUICKMATH
#endif
#else
struct UserOutputInterface;
extern UserOutputInterface * global_uo;
#endif
/////////////////////////


#include "mbs_interface.h"

#include "femathHelperFunctions.h"

typedef TArray <DVector*> DArray;
typedef TArray <IVector*> IArray;
typedef TArray <Vector*> VArray;
typedef TArray <Vector3D> AVector3D;

struct SuperMatrix;
//#include "..\SuperLU\supermatrix.h"
#include "..\SuperLU\SuperLUmain.h" //for external solver

#include "linalg.h"
#include "linalgeig.h"
#include "elementdata.h"      //$AD 2011-03-24: new entry for use of parser in mathfunc
#include "parser.h"                           //$AD 2011-03-24: new entry for use of parser in mathfunc
#include "mathfunc.h"
#include <sstream>

#include "..\workingmodule\WorkingModuleBaseClass.h"

class NumNLSys;


class SparseJacMat
{
public:
	SparseJacMat(): J_vv_band_pivot(), J_vv(), J_zv(), J_vz(), J_zz()
		, sparseinv(0)
	{
		LUcomputed = 0;
		LUcomputed_bw = 0;
		Init();
	};

	void Init()
	{
		J_vv.Init();
		J_vv_band_pivot.Init();
		LUcomputed = 0;
		LUcomputed_bw = 0;
		useband = 0;
		n_vv_written = 0;

		sparseinv = 0;

	}

	~SparseJacMat()
	{
		delete sparseinv;
	}

	void Apply(const Vector& R, Vector& d, NumNLSys* nls);

	int Factorize(NumNLSys* nls);

	void ApplyTransform(Vector& v, int mode, NumNLSys* nls);

	SparseMatrix J_vv;
	SparseMatrix J_zvS, J_vzS;

	Matrix J_zv, J_vz, J_zz;
	Matrix J_sch;
	Matrix J_vv_band;
	Matrix J_vv_band_aux;
	IVector J_vv_band_pivot;
	Matrix minv_test;

	IVector LUindx; //temporary storage full LU decomposition
	Vector LUvv;    //temporary storage full LU decomposition

	//++++++++++++++++++++++++++++
	// SUPERLU; PARDISO

	SparseInverse *sparseinv;
	//++++++++++++++++++++++++++++

	int lu;
	int LUcomputed;
	int LUcomputed_bw;
	int n_vv_written; //set to 1 if state has been already written!

	//for renumbering of band-diagonal part (rigid bodies!!!)
	IVector resortvector;
	int resortsize;

	Vector temp;   //temporary variables for Apply
	Vector temp2;  //temporary variables for Apply
	Vector Rsort;  //temporary variables for Apply
	Vector vtrans; //temporary variable for ApplyTransform

	Matrix mtemp1, mtemp2;  //temporary variables for Factorize
	Vector tempcol;         //temporary variables for Factorize

	int useband;
};

class SaveJac
{
public:
	SaveJac(): oldjacmat(), oldsjacmat(), sparsejacmat()
	{
		oldjacmat.Init();
		oldsjacmat.Init();
		sparsejacmat.Init();
		Init();
	}

	void Init()
	{
		oldjac=0; //jacobimatrix existiert
		oldjacsize=0; //größe
		oldjacage=0; //alter
		maxjacage=10000000;//10000000; //maximales alter
		nlsinfo = 0; //nonlinear solve info (z.b. zeitschrittweite von timeint.cpp)
		jacfailedcnt = 0;
		lastbuilt = 0;

		oldjacmat.SetSize(0,0);
		oldsjacmat.SetSize(0,0);
		sparsejacmat.Init();
	}

	int oldjac;
	int oldjacsize;
	int oldjacage;
	int maxjacage;
	Matrix oldjacmat;
	SparseMatrix oldsjacmat;
	SparseJacMat sparsejacmat;
	double nlsinfo;
	int jacfailedcnt;
	int lastbuilt;
};

#include "options_class_auto.h"

class NumNLSolver;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Numerical nonlinear solver
//Needs a function which is solved for x f(x)=0
//Jacobion matrix is calculated by numerical differentiation
//$!YV 2012-11-29: inherits MBS interface, which is then implemented in the class hierarchy upwards
class NumNLSys : public WorkingModuleBaseClass, public MBS
{

public:
	NumNLSys() {jaccol = 0; TIstages=1; jacfullnewton=0; fullnewtoncnt = 0; doresort = 0; hotint_options = new HOTINTOptions(this);};
	~NumNLSys() { DeleteOptions(); }
	virtual void NLF(const Vector& x, Vector& f)=0; //nonlinear function F(x) = 0
	virtual void SetSolver(NumNLSolver* s);
	virtual void Jacobian(Matrix& m, Vector& x);
	virtual void Jacobian(SparseMatrix& m, Vector& x) {};

	int NLS_GetJacCol() const;
	void NLS_SetJacCol(int i) {jaccol = i;}
	int NLS_GetTIstages() const {return TIstages;}
	void NLS_SetTIstages(int i) {TIstages=i;};
	int NLS_IsJacFullNewton() const {return jacfullnewton;}
	void NLS_SetJacFullNewton(int i) {jacfullnewton=i;};
	int FullNewtonCnt() const {return fullnewtoncnt;}
	int& FullNewtonCnt() {return fullnewtoncnt;}
	void ResetNLSys() {fullnewtoncnt = 0;}

	virtual const IVector& GetResortVector() const {return resortvector;}
	virtual IVector& GetResortVector() {return resortvector;}
	virtual int GetResortSize() const {return resortsize;};
	virtual void SetResortSize(int s) {resortsize = s;};

	virtual int UseSparseSolver() const;
	virtual int SolveUndeterminedSystem() const;
	virtual int& SolveUndeterminedSystem();
	virtual double EstimatedConditionNumber() const;
	virtual double& EstimatedConditionNumber();
	virtual int UseSparseJac() const;
	virtual int UseSuperLU() const;

	virtual int TransformJacApply() const {return 0;}
	virtual void ApplyTransform(const Vector& v, Vector& Av, int mode) {}; //compute Av=A^T*v in mode==0 and Av=A*v in mode==1

	//[I|D|T]Options
	virtual void SetIOption(int index, int data); 
	virtual const int& GetIOption(int index) const; 
	virtual int& GetIOption(int index); 
	virtual void SetDOption(int index, double data); 
	virtual const double& GetDOption(int index) const; 
	virtual double& GetDOption(int index); 

	virtual void SetTOption(int index, const char* data); 
	virtual const char* GetTOption(int index) const; 
	virtual void InitializeOptions();
	//+++++++++++++++++++++++++++++++++++++++++++++++++
	// implemented in auto - generated file "options_auto.h"
	virtual void InitializeOptionsAuto();  
	//+++++++++++++++++++++++++++++++++++++++++++++++++
	virtual void DeleteOptions();

	virtual HOTINTOptions* GetOptions() { return hotint_options; }
	virtual const HOTINTOptions* GetOptions() const { return hotint_options; }

	/////////////////////////////
	//#output details of computation:              (changes in the 'Edit All Options' window take effect immediately, since the solset.* values are returned here!)
	
	// definition of user output level for detailed solver logs
	UO_MSGLVL UO_LVL_SolverDetails() const { return UO_LVL_all; }

	// solver details are only printed to logfile, if file_output_level is set greater or equal to UO_LVL_SolverDetails()
	// solver details are only printed to window, if output_level is set greater or equal to UO_LVL_SolverDetails()
	bool SolverPrintsDetailedOutput()
	{ 
		return 
			(
			SolverUO().GetGlobalFileMessageLevel() >= UO_LVL_SolverDetails() || 
			SolverUO().GetGlobalMessageLevel() >= UO_LVL_SolverDetails()
			)
			&&
			AnySolverOutputChecked();
	}

	bool AnySolverOutputChecked()
	{ 
		return
			GetOptions()->LoggingOptions()->SolverGeneralInformation() ||
			GetOptions()->LoggingOptions()->SolverNewtonIterationJacobiMatrix() ||
			GetOptions()->LoggingOptions()->SolverNewtonIterationResidualVector() ||
			GetOptions()->LoggingOptions()->SolverNewtonIterationSolutionVector() ||
			GetOptions()->LoggingOptions()->SolverPostNewtonIterationDataVector() ||
			GetOptions()->LoggingOptions()->SolverStepSolutionVector() ||
			GetOptions()->LoggingOptions()->SolverStepSolutionVectorIncrement();
	}

	UserOutput& SolverUO() {	return UOfull(UO_LVL_SolverDetails()); }
	UserOutput& SolverUO() const	{	return UOfull(UO_LVL_SolverDetails()); }
	
	//virtual UserOutput& UO(int message_level = UO_LVL_all) {uo.SetLocalMessageLevel(message_level); return uo;}; //(AD) 
	virtual UserOutput& UOfull(int message_level = UO_LVL_all, int output_prec = -1) const; //$ AD 2011-02 output_prec
	virtual UserOutputInterface& UO(int message_level = UO_LVL_all, int output_prec = -1) const { return UOfull(message_level,output_prec); }

	//option members:
	int ioptions[301];
	double doptions[301];
	char* toptions[301];


private:
	int jaccol; 
	int jacfullnewton;
	int TIstages;
	NumNLSolver* solver;
	int fullnewtoncnt;

	Vector f1; //temporary variables for Jacobian
	Vector f2; //temporary variables for Jacobian

protected:
	int resortsize;
	int doresort;
	IVector resortvector;
	HOTINTOptions* hotint_options;
};


class NumNLSolver : public NumSolverInterface
{
	SolverSettings* solset; //this part of the solversettings are referenced to the settings, which can be online modified by user; do only use solset for options which should show
													//immediate effect of change, e.g. the logging of newton iterations or errors; DO NOT USE E.G. FOR NEWTON PARAMETERS (e.g. TOLERANCES)==>THIS LEADS TO ARBITRARY RESULTS
public:
	NumNLSolver():iv() {};
	NumNLSolver(NumNLSys* nlsi, SolverSettings* solsetI):iv()
	{
		nls = nlsi;
		solset = solsetI;
		nls->SetSolver(this);

		relativeaccuracy = 1e-8;
		maxmodnewtonsteps=50;
		maxrestartnewtonsteps=25;
		maxfullnewtonsteps=20;
		numdiffepsi = 1e-8;
		newtonits = 0;
		trustregion = 0;
		//trustregionitmax = 6;
		trustregiondiv = 1./10.;
		output=0;
		modifiednewton = 1; // Sets solving-method to modified Newton
		symmetricjacobian = 1;
		nonlinsolveinfo = 0;
		stopmnr = 0;
		usesparsesolver = 0;
		solveundeterminedsystem = 0;
		estimatedconditionnumber = 1e12;

		contractivity = 0;

		jaccount=0;
		fulljaccnt = 0;

		jac_condnum=-1;
		error_msg = "";

		bandsize = 0;

		//DestroyOldJacs();
		//sjac[0].oldjacmat = Matrix();
		//sjac[1].oldjacmat = Matrix();

		//log.SetAllInactive();
	}

	void ResetSolver()
	{
		relativeaccuracy = 1e-8;
		maxmodnewtonsteps=50;
		maxrestartnewtonsteps=25;
		maxfullnewtonsteps=20;
		numdiffepsi = 1e-8;
		newtonits = 0;
		trustregion = 0;
		//trustregionitmax = 6;
		trustregiondiv = 1./10.;
		output=0;
		modifiednewton = 1; // Sets solving-method to modified Newton
		symmetricjacobian = 1;
		nonlinsolveinfo = 0;
		stopmnr = 0;

		contractivity = 0;

		jaccount=0;
		fulljaccnt = 0;

		jac_condnum=-1;
		error_msg = "";

		bandsize = 0;

		DestroyOldJacs();
		//sjac[0].oldjacmat = Matrix();
		//sjac[1].oldjacmat = Matrix();
		sjac[0].Init();
		sjac[1].Init();

	}

	//Newton Solver
	int NLSolve(Vector& x0);

	int Factorize(Matrix& minv, SparseMatrix& sminv, SparseJacMat& sparseminv);

	//set res=Jac^(-1)*f if jac exists, else return 0
	int ApplyJac(const Vector& f, Vector& res);

	int ChooseJac()
	{
		if (sjac[0].oldjac && RelApproxi(sjac[0].nlsinfo,nonlinsolveinfo,2e-2)) 
		{
			return 1;
		}
		else if (sjac[1].oldjac && RelApproxi(sjac[1].nlsinfo,nonlinsolveinfo,2e-2)) 
		{
			return 2;
		}
		else return 0; 
	}

	//virtual void Jacobian(Matrix& m, Vector& x);

	int GetNewtonIts() const {return newtonits;};

	int ModifiedNewton() const {return modifiednewton;};
	int& ModifiedNewton() {return modifiednewton;};

	double AbsoluteAccuracy() const {return absoluteaccuracy;};
	double& AbsoluteAccuracy() {return absoluteaccuracy;};

	double RelativeAccuracy() const {return relativeaccuracy;};
	double& RelativeAccuracy() {return relativeaccuracy;};

	//int& MaxNewtonSteps() {return maxnewtonsteps;};   //$ PG 2013-9-19: who changed that???
	//int MaxNewtonSteps() const {return MaxModNewtonSteps()+MaxRestartNewtonSteps()+MaxFullNewtonSteps();};    //$ PG 2013-9-19: does not make any sense in case of full newton!
	int MaxNewtonSteps() const 
	{
		if (!ModifiedNewton())
		{
			return MaxFullNewtonSteps();
		}

		return MaxModNewtonSteps()+MaxRestartNewtonSteps()+MaxFullNewtonSteps();
	}

	int MaxModNewtonSteps() const {return maxmodnewtonsteps;};
	int& MaxModNewtonSteps() {return maxmodnewtonsteps;};

	int MaxRestartNewtonSteps() const {return maxrestartnewtonsteps;};
	int& MaxRestartNewtonSteps() {return maxrestartnewtonsteps;};

	int MaxFullNewtonSteps() const {return maxfullnewtonsteps;};
	int& MaxFullNewtonSteps() {return maxfullnewtonsteps;};

	double NumDiffepsi() const {return numdiffepsi;};
	double& NumDiffepsi() {return numdiffepsi;};

	int TrustRegion() const {return trustregion;};
	int& TrustRegion() {return trustregion;};

	//int TrustRegionItMax() const {return trustregionitmax;};
	//int& TrustRegionItMax() {return trustregionitmax;};

	double TrustRegionDiv() const {return trustregiondiv;};
	double& TrustRegionDiv() {return trustregiondiv;};

	int Output() const {return output;};
	int& Output() {return output;};

	int SymmetricJacobian() const {return symmetricjacobian;};
	int& SymmetricJacobian() {return symmetricjacobian;};

	double NLSolveInfo() const {return nonlinsolveinfo;};
	double& NLSolveInfo() {return nonlinsolveinfo;}; 

	double Contractivity() const {return contractivity;}

	const Matrix& GetJac() const 
	{
		if (sjac[0].lastbuilt) return sjac[0].oldjacmat;
		else return sjac[1].oldjacmat;
	}
	void DestroyOldJac(double info) 
	{
		if (info == sjac[0].nlsinfo)
		{
			sjac[0].oldjac = 0; sjac[0].nlsinfo = -1; sjac[0].lastbuilt = 0; sjac[1].lastbuilt = 1;
		}
		else if (info == sjac[1].nlsinfo)
		{
			sjac[1].oldjac = 0; sjac[1].nlsinfo = -1; sjac[0].lastbuilt = 1; sjac[1].lastbuilt = 0;
		}

	}
	void DestroyOldJacs() 
	{
		sjac[0].oldjac = 0; sjac[0].nlsinfo = -1; sjac[0].lastbuilt = 0;
		sjac[1].oldjac = 0; sjac[1].nlsinfo = -1; sjac[1].lastbuilt = 0;
	}

	int ModifiedNewtonActive() const {return ModifiedNewton() && (!stopmnr);}

	int GetJacCount() const {return jaccount;};
	int GetFullJacCnt() const {return fulljaccnt;};
	double GetJacCondnum() const {return jac_condnum;};
	const mystr& GetErrorMsg() const {return error_msg;};
	mystr& GetErrorMsg() {return error_msg;};

	int GetBandSize() const {return bandsize;}
	void SetBandSize(int bs) {bandsize = bs;}

	int UseSparseSolver() const { return usesparsesolver; }
	int& UseSparseSolver() { return usesparsesolver; }
	int SolveUndeterminedSystem() const { return solveundeterminedsystem;	}
	int& SolveUndeterminedSystem() { return solveundeterminedsystem; }
	double EstimatedConditionNumber() const	{	return estimatedconditionnumber; }
	double& EstimatedConditionNumber() { return estimatedconditionnumber; }
	

	//print detailed output concerning NLSolve computations? (time critical!)
	bool SolverPrintsDetailedOutput() const { return nls->SolverPrintsDetailedOutput(); }
	UserOutput& SolverUO() { return nls->SolverUO(); }
	UserOutput& SolverUO() const { return nls->SolverUO(); }
	HOTINTOptions* GetOptions() { return nls->GetOptions(); }


private:

	//$ PG 2013-11-6: assemble Jacobian (and print some information), set jaccount++
	void AssembleJacobian(Matrix* minv, SparseMatrix* sminv, Vector& x0);

	//$ PG 2013-10-17: calculate condition number of jacobi matrix, IS SLOW, should only be used for designing elements, constraints, for adjusting models, or for general debugging purpose
	void OutputConditionNumber(const SparseMatrix& M) const;
	void OutputConditionNumber(const Matrix& M) const;
	void OutputJacobian(const SparseMatrix& M) const;
	void OutputJacobian(const Matrix& M) const;

	int newtonits;
	NumNLSys* nls;

	SaveJac sjac[2];
	
	//solver settings
	int modifiednewton; //indicates modified Newton algorithm is used
	double absoluteaccuracy;
	double relativeaccuracy;
	double numdiffepsi;
	int maxmodnewtonsteps, maxrestartnewtonsteps, maxfullnewtonsteps;
	int stopmnr;

	int trustregion; //turn on or off
	//int trustregionitmax;
	double trustregiondiv;

	int output;

	int symmetricjacobian;

	double nonlinsolveinfo;
	int usesparsesolver;
	int solveundeterminedsystem;
	double estimatedconditionnumber;

	int jaccount;
	double contractivity;
	int fulljaccnt;
	double jac_condnum;
  mystr error_msg;

	//sparse:
	int bandsize; //size of leftupper system (first n unknowns) which has small bandwidth, esp. for MBS


	Vector x0start, f, xd, hv; //temporary variables in NLsolve
	IVector iv;                //temporary variables in NLsolve

};



#endif
