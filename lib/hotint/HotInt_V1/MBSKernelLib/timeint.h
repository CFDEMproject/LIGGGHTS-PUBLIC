//#**************************************************************
//#
//# filename:             timeint.h
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						09.06.99
//# description:          Class for implicit and explicit Runge Kutta Time integration, 
//#												variable stepsize and arbitrary order
//# remarks:						  The file tableaus.txt must be provided in Project/Release or Debug !!!
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
 

#ifndef TIMEINT__H
#define TIMEINT__H

#include "mystring.h"
#include "timeintlog.h"
#include "femath.h"

//# maximum stages for static variables in K-Version timeint
#define global_maxstages 21



class IRK_Tableau
{
public:
	IRK_Tableau() {};

	mystr name;
	mystr info;
	int nstages;
	int ODE_order;
	int DAE_order_y[10]; //DAE_order[1]=index1 order, DAE_order[2]=index2 order, DAE_order[3]=index3 order, ...
	int DAE_order_z[10]; //DAE_order[1]=index1 order, DAE_order[2]=index2 order, DAE_order[3]=index3 order, ...

	int Ainvertable;
	int implicit;

	Vector b;
	Vector c;
	Matrix A;
	Matrix A2;
	Matrix Ainv;
	Vector bAinv;  ///b^T * A^(-1)
};





//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class TimeIntVars
{
public:
	TimeIntVars() {Init();}
	virtual void Init()
	{
		start_clock_time = 0;
		timetogo = 0;
		lastdraw = 0;
		lastprintres = 0;
		last_showstatustext = 0;
		timetogo2 = 0;
		outputlevel = 3; //2
		maxnewtonit=0;
	}

	//time, drawing, log messages:
	double start_clock_time; //stores when clock was started
	double timetogo; //estimated time to go
	double lastdraw; //time when last was drawn
	double lastprintres; //time when last was printed results
	double last_showstatustext; //time when last was printed results
	double timetogo2; //damping term for time to go estimation
	int outputlevel; //output level, how much is printed per step
	
	//temporary solver variables:
	int maxnewtonit; //maximum number of newton iterations since last printing
};



//structure of data and state in one time step
struct TIstepsave
{
	Vector x0; //state variables
	Vector data;  //data variables
	Vector k0; //state variables in differentiated form
	double time;
};


class TimeInt : public NumNLSys
{
public:
	TimeInt();

	virtual ~TimeInt()
	{
		for (int i=1; i <= tableaus.Length(); i++)
		{
			delete tableaus.Elem(i);
		}
	}
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//#main functions for time integration
	//#x = [u v q z] // \dot u = v, M (TIx0) \dot v = F2(TIx0,t), \dot q = F2(TIx0,t), G(TIx0) = 0
	//#len(u) = len(v), [len(u)=0 | len(q) = 0], len(z) arbitrary
	virtual void PrecomputeEvalFunctions() {}; //$ PG 2012-1-15: precomputation for EvalF(), EvalF2(), EvalG(), EvalM(), and also for CqTLambda() in case of constraints
	virtual void EvalF(const Vector& x, Vector& f, double t) {};  //*first order equations: \dot q=f(q,t), return vector:  len(q)
	virtual void EvalM(const Vector& x, Matrix& m, double t) {};  //*evaluate mass matrix M(TIx0), return matrix: len(u)*len(u)
	virtual void EvalM(const Vector& x, SparseMatrix& m, double t) {};  //*evaluate mass matrix M(TIx0), return matrix: len(u)*len(u)
	virtual void SetBWM(int i) {bandwidthm = i;}
	virtual void EvalF2(const Vector& x, Vector& f, double t) {}; //*second order equations: M \ddot u = F2(u,\dot u,t), len(u)
	virtual void EvalMinvF2(const Vector& x, Vector& f, double t) {}; //*second order equations: M \ddot u = F2(u,\dot u,t), len(u)
	//virtual void EvalF2(const Vector& x, Vector& f, double t, int& minind, int& maxind) {}; //only for jacobian, provide also min and max entry
	virtual void EvalF2(const Vector& x, SparseVector& f, double t, 
		IVector& rowind, IVector& clearind, IVector& elems) {}; //only for jacobian, rowind and clearind for speedup, rowind must be initialized with zeros and l=f.Length()
	virtual void EvalFelastic(const Vector& x, Vector& f, double t) {};

	virtual void EvalG(const Vector& x, Vector& f, double t) {}; //*evaluate constraints: len(z)
	virtual void EvalG(const Vector& x, SparseVector& f, double t, IVector& rowind, IVector& clearind) {}; //only for jacobian, rowind and clearind for speedup, rowind must be initialized with zeros and l=f.Length()

	virtual int GetSecondOrderSize() const {return 0;};  //*size of second order equations, len(u)
	virtual int GetSecondOrderSize_RS() const {return GetSecondOrderSize();};  //*size of second order equations, len(u), after resorting!!!
	virtual int GetFirstOrderSize() const {return 0;};  //*size of first order equations
	virtual int GetImplicitSize() const {return 0;};  //* len(z)
	virtual int GetSystemSize() const {return GetSecondOrderSize()*2
		+GetFirstOrderSize()+GetImplicitSize();}; //*len(z)+len(u)+len(v)+len(q)
	virtual int GetDataSize() const {return 0;} //size of data variables for each time step of system


	virtual void NLF(const Vector& x, Vector& f); //*
	virtual int NonlinSubStep(); //*
	virtual void FinishStep(); //*
	virtual void TIDrawAndStore(); //*
	virtual double GetEnergy() const {return 0;};
	virtual void Initialize() {}; //initialize MBS system
	virtual void SetInitialConditions() {}; //initial conditions and some initializations of MBS
	virtual void InitFirst();			//initialize at program start when UserInterface is linked, before config file
	virtual void InitializeAfterConfigLoaded() {};
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//#internal functions for simulation
	void LoadTableaus(const mystr& tableauname);  // nur filename mit tableaus
	int GetTableau(mystr tabname, int stages); //tableau-nummer auswählen
	int LoadTableauList(mystr tabname, 
		TArray<int>& list, 
		int minstage, int maxstage);
	IRK_Tableau* GetActTableau() const {return tableaus[act_tableau];}

	virtual void TIInit();  //*
	virtual int FullStep(); //*
	virtual int FirstStep();
	virtual int GeneralIStep();

	virtual double GetCharacteristicSol() const {return TIx0(1);}

	virtual void PrintTimingList();
	virtual void ResetComputation();
	virtual int Integrate();//JG  Integrate(const SolverSettings& solversettings);

	virtual int IntegrateStep();  //integrate one step (adaptive or constant stepsize)
	virtual int StaticComputation();  //try to solve static problem

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	virtual int ExplicitIntegration();  //do explicit integration, return 1 if success
	virtual int AdaptiveExplicitStep(); //perform one adaptive explicit step, return 1 if success
	virtual int ExplicitStep();					//do one explicit nonlinear step, return 1 if success
	virtual int GeneralEStep();					//Do one general explicit step
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //# Computation Steps with changing SolverOptions...
	//Time mapping functions
	virtual int ParseStepEndTimes() {return 0;};										// parses all EDC for ComputationSteps entries - implemented in MBS
	virtual int DoubleCheckWithLoadSettings() {return 0;};

	TArray<double> CSEndTimes;																			// array the holds all endtimes of the computation steps
	virtual double GetTMaxCompSteps() {return CSEndTimes.Last();};
	virtual double GetEndTimeCompStep(int stepnumber) 
	{
		if (stepnumber <= CSEndTimes.Length()) return CSEndTimes(stepnumber);
		else return solset.endtime; // in case last computation step finishes BEFORE end of simulation
	};

	int computationstepnumber;																			// number of the computation step ("current")
	int& ComputationStepNumber() {return computationstepnumber;};   // quick access to the computationstepnumber (variable, no recalculaiton)
	double computationsteptime;                                     // relative time within the current computation step ("percentage")
	double& ComputationStepTime() {return computationsteptime;};    // quick access to the computaionsteptime (variable, no recalculation)

	virtual int GetCSNumber(double globaltime);											// calculates number of the Computation Step for given global time // returns -1 it globaltime < 0 or globaltime > last step_end_time
  virtual double GetCSTime(double globaltime);										// calculates relative time in corresponding step   rv(0..1] // returns 0. only for globaltime == 0 or globaltime > last step_end_time
	virtual int IsComputationStepFinished(double globaltime);       // returns 1 if globaltime is any step_end_time

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //# Computation Steps with changing SolverOptions...
  //replace solver options
	virtual int ApplyComputationStepSettings(int stepnumber) {return 0;}  // changes the solver options - implemented in MBS
	virtual void ComputationStepFinished() {};                            //function is called when a computation step is finished

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	virtual void StepPrintAndDraw(int rv); //print step information and draw results each step

	virtual void SetComputationSolverOptions() {};
	virtual int GetTMaxLoadStepsLength() const {return 0;} //overloaded in mbs.h
	//Compute weights for Lagrangian interpolation, tau=0..1
	double LagrangeWeight(int l, double tau, const Vector& c) const;

	virtual int FullAdaptiveStep();

  virtual void StartTimeStep() {}; //function is called when computation of time step is started
  virtual void EndTimeStep() {}; //function is called when computation of time step is started

	virtual double PostNewtonStep(double ) {return 0;};
	virtual void PostprocessingStep() {};

	//#save and restore inner variables for calculating the error, restarting, ...
	virtual void SaveState() {};
	virtual void RestoreState() {};
	virtual void TISaveState(TIstepsave& stepsave);
	virtual void TIRestoreState(const TIstepsave& stepsave);
	virtual void LocalJacobianF(Matrix& m, Vector& x);  //Compute (df/dx)
	virtual void LocalJacobianF2(Matrix& m, Vector& x);  //Compute (df2/dx)
	virtual void LocalJacobianF2(SparseMatrix& m, Vector& x);  //Compute (df2/dx)
	virtual void LocalJacobianG(Matrix& m, Vector& x);  //Compute (dg/dx)
	virtual void LocalJacobianG(SparseMatrix& m, Vector& x);  //Compute (df2/dx)
	virtual void LocalJacobianM(Matrix& m, Vector& x);  //Compute (d(M*TIkv)/dx)
	virtual	void Jacobian(Matrix& m, Vector& x);  //Compute approximate Jacobian (if modified newton) or call base function
	virtual	void Jacobian(SparseMatrix& m, Vector& x);  //Compute approximate Jacobian (if modified newton) or call base function
	virtual	void StaticJacobian(SparseMatrix& m, Vector& x);  //Compute approximate Jacobian for static case (if modified newton) or call base function
	virtual	void FullJacobian(Matrix& m, Vector& x);  //Compute approximate Jacobian (if modified newton) or call base function
	virtual void ForceJacobianRecomputation() {numsol.DestroyOldJacs();} //Destroy actual jacobians and compute new jacobian, e.g. if stiffness or configuration of system changes has been changed significantly


	virtual void ComputeStiffnessAndDampingMatrix(SparseMatrix& k, SparseMatrix& d, Vector& x);  //Compute (df2/dx),(df2/dxp), directly fill into final stiffness and damping matrix
	virtual void ComputeGyroscopicMatrix(SparseMatrix& gy); // $ MSax 2013-07-25 : added

	virtual void StaticJacobianF2(SparseMatrix& m, Vector& x) {};  //Compute (df2/dx), directly fill into final stiffness matrix

	virtual void ComputeSparseMatrixUsage(IVector& usage_per_dof) {}; //compute entries of each row in sparse matrix

	virtual int UseSparseMass() const {return UseSparseJac();}
	virtual int MaxSparseBandwidth() const {return 32;} //default value
	virtual void SetReducedBandsize(int reducedbandsizeI) {reducedbandsize = reducedbandsizeI;}
	virtual int DoStaticComputation() const {return solset.dostaticcomputation;}
	virtual int& DoStaticComputation() {return solset.dostaticcomputation;}
	virtual int DoImplicitIntegration() const {return solset.doimplicitintegration;}
	virtual int& DoImplicitIntegration() {return solset.doimplicitintegration;}
	virtual double LoadFact() {return loadfact;}
	

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//#Access to solution and computation

	//#starting the simulation:
	virtual void SetStartVector(const Vector& x0);
	//#get characteristic value for error; this error is attempted to be lower than the desired accuracy
	virtual double GetError() {return TIx0.GetNorm();};

	virtual void TIFinished() {TIfinished = 1;}
	virtual void TIFinishedWithError() {TIfinished = -1;}
	virtual const int& GetTIFinished() const {return TIfinished;}
	virtual int& GetTIFinished() {return TIfinished;}
	virtual int GetTIit() const {return TIit;}

	virtual int GetNewtonItSum() const {return log.TInewtonitsum;}
	virtual int GetNewtonIt() const {return log.TInewtonit;}
	virtual void SetNewtonItSum(int ii) {log.TInewtonitsum = ii;}
	virtual void SetNewtonIt(int ii) {log.TInewtonit = ii;}
	//virtual Vector GetNonlinSol() {return Vector(0);} //erase
	virtual int GetNonlinIts() const {return log.TInonlinit;}
	
	virtual const Vector& GetSolVector() const {return TIx0;} //return actual states
	virtual Vector& GetSolVector() {return TIx0;}
	virtual void SetSolVector(const Vector& sol) {TIx0 = sol;}
	virtual const double& GetSol(int i) const {return TIx0(i);}

	virtual const Vector& GetDataVector() const {return TIdata;} //return non-state data vector (plastic strains, contact conditions, discrete states, etc.)
	virtual Vector& GetDataVector() {return TIdata;}   
	virtual void SetDataVector(const Vector& data) {TIdata = data; TIdrawdata = data;} //set non-state data vector (plastic strains, contact conditions, discrete states, etc.)

	//$EK 2012-09-17: return non-state data draw vector (plastic strains, contact conditions, discrete states, etc.)
	virtual const Vector& GetDataDrawVector() const {return TIdrawdata;} 
	virtual Vector& GetDataDrawVector() {return TIdrawdata;}   

	virtual const double& GetDataAct(int i) const {return TIdata(i);};
	virtual double& GetDataAct(int i) {return TIdata(i);};

	virtual const double& GetDataDraw(int i) const {return TIdrawdata(i);};
	virtual double& GetDataDraw(int i) {return TIdrawdata(i);};

	virtual const Vector& GetLastSolVector() const {return TIlaststep_state.x0;} //return solution after end of last step == beginning of this step
	virtual const Vector& GetLastNLItSolVector() const {return TIlastnonlinit_state.x0;} //return solution after last nonlinear iteration

	virtual Vector& GetLastDataVector() {return TIlaststep_state.data;} //return data vector after end of last step == beginning of this step
	virtual const Vector& GetLastDataVector() const {return TIlaststep_state.data;} //return data vector after end of last step == beginning of this step
	virtual Vector& GetTempDataVector() {return TItemp_state.data;} //return data vector after StartTimeStep()
	virtual const Vector& GetTempDataVector() const {return TItemp_state.data;} //return data vector after StartTimeStep()
	virtual Vector& GetLastNLItDataVector() {return TIlastnonlinit_state.data;} //return data vector after last nonlinear iteration
	virtual const Vector& GetLastNLItDataVector() const {return TIlastnonlinit_state.data;} //return data vector after last nonlinear iteration
	virtual const int& GetNumberOfSolutionDataSteps() const { return TInumber_of_solution_data_steps; } //return number of stored soltuion data steps (which is particularly used for communication with data manager of WCDriver)
	virtual void SetNumberOfSolutionDataSteps(const int n_steps) { TInumber_of_solution_data_steps = n_steps; } //set number of stored soltuion data steps (which is particularly used for communication with data manager of WCDriver)

	virtual const Vector& GetVelocityAndAccelerationVector() const {return TIk0;}   //for 2nd order components, 1..ss // $ MSax 2013-07-16 : renamed from GetAccelerationVector to GetVelocityAndAccelerationVector
	virtual const double& GetVelocityAndAcceleration(int i) const {return TIk0(i);} //1..ss: 1..sos:velocities, sos+1 .. 2*sos:accelerations, 2*sos+1 .. 2*sos+es:first order velocities, 2*sos+es+1..:algebraic variables // $ MSax 2013-07-16 : renamed from GetAcceleration to GetVelocityAndAcceleration

	virtual double GetStepSize() const {return TIstep;}
	virtual void SetStepSize(double s) {TIstepnew = s;}
	virtual double GetStepSizeNew() {return TIstepnew;} // necessary for time discrete input output elements

	virtual double GetMinStepSize() const {return TIminstep;}


	virtual double GetStepRecommendation() const {return TIsteprecommendation;}
	virtual void SetStepRecommendation(double s) {TIsteprecommendation = s;}

	virtual void SetTime(double tt) {TItime = tt; TIdrawtime = tt;}
	virtual double GetTime() const {return TItime;}

	virtual double GetStepEndTime() const {return TItime+TIstep;}
	virtual double GetDrawTime() const {return TIdrawtime;}		// for MBS interface
	virtual double GetActualDrawTime() const {return TIdrawtime;}	// for WinCompDriverInterface
	virtual void SetDrawTime(double t) {TIdrawtime = t;}
	virtual const Vector& GetDrawVector() const {return TIx0draw;}
	virtual Vector& GetDrawVector() {return TIx0draw;}
	virtual const double& GetDrawValue(int i) const {return TIx0draw(i);}
	virtual double& GetDrawValue(int i) {return TIx0draw(i);}

	virtual double GetCompTime() const {return TIcomptime;}
	virtual double GetClockTime() const;
	virtual int GetStepChanges() const {return log.changestep;}
	virtual int GetJacCount() const {return log.jaccount;}

	//#write problem dependent solution data
	virtual void WriteSol() {}

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//main functions for data handling

	// these functions are due to the data storage
	// after each call of ComputationFeedBack::ResultsUpdated()
	// the function StoreResultsIsOn() is called to check if the results should be stored
	virtual int StoreResultsIsOn();

	// here the working object should remove all stored states
	virtual int RemoveResults();

	// here the working object should save its actual state
	virtual void StoreResults(DataSaver & storage, double & m_TimePoint);

	// here the working object should replace its actual state with the one from the DataLoader
	// the function is called from the driving module with the user request
	virtual void LoadResults(DataLoader & loader, int m_TimePointNumber);

	// this line identifies the current data storage structure
	// when the data set or the way of storage is changed,
	// new version identifier should be introduced
	//virtual const char * GetDataFormatVersionIdentifier() {return "Version0.1";};

	//TIME-INT VARIABLES
	virtual const SolverSettings& GetSolSet() const {return solset;}		//$ YV 2012-11-29: MBSSolSet does not exist any longer
	virtual SolverSettings& GetSolSet() {return solset;}
	virtual void SetSolSet(const SolverSettings& solsetI) {solset = solsetI;}

	const double& WriteSolTau() const {return TIwritesolevery;} 
	double& WriteSolTau() {return TIwritesolevery;} 

	int MaxDiscontinuousIt() const {return solset.maxdiscontinuousit;};
	int& MaxDiscontinuousIt() {return solset.maxdiscontinuousit;};

	double DiscontinuousAccuracy() const {return solset.discontinuousaccuracy;};
	double& DiscontinuousAccuracy() {return solset.discontinuousaccuracy;};

	//virtual int StopCalculation() const {return 1; };
	virtual int StopCalculation() const {return bStopWhenPossible;};
	virtual int PauseCalculation() const {return bPause;};

	double AbsAccuracy() const {return solset.absaccuracy;};
	double& AbsAccuracy() {return solset.absaccuracy;};

	int FullAdaptive() const {return solset.fulladaptive;};
	int& FullAdaptive() {return solset.fulladaptive;};

	int KVersion() const {return TIkv;};
	int& KVersion() {return TIkv;};

	void SetStepDiscontinuous() {TIdiscontstep++;}
	int IsStepDiscontinuous() const {return (TIdiscontstep!=0);}

	virtual int IsJacobianComputation() const {return computejacflag;}
	virtual void SetJacobianComputationFlag(int i) {computejacflag = i;}

	virtual int NLS_UseSparseSolver() const {return solset.nls_usesparsesolver;};
	virtual int& NLS_UseSparseSolver() {return solset.nls_usesparsesolver;};

	virtual int NLS_SolveUndeterminedSystem() const {return solset.nls_solve_undetermined_system;};
	virtual double NLS_EstimatedConditionNumber() const	{return solset.nls_estimated_condition_number;}

	//NEWTON-METHOD:
	int NLS_ModifiedNewton() const {return solset.nls_modifiednewton;};
	int& NLS_ModifiedNewton() {return solset.nls_modifiednewton;};

	double NLS_AbsoluteAccuracy() const {return solset.nls_absoluteaccuracy;};
	double& NLS_AbsoluteAccuracy() {return solset.nls_absoluteaccuracy;};

	double NLS_RelativeAccuracy() const {return solset.nls_relativeaccuracy;};
	double& NLS_RelativeAccuracy() {return solset.nls_relativeaccuracy;};

	double NLS_NumDiffepsi() const {return solset.nls_numdiffepsi;};
	double& NLS_NumDiffepsi() {return solset.nls_numdiffepsi;};

	int NLS_MaxNewtonSteps() const {return solset.nls_maxmodnewtonsteps + solset.nls_maxrestartnewtonsteps + solset.nls_maxfullnewtonsteps;};

	int NLS_MaxModNewtonSteps() const {return solset.nls_maxmodnewtonsteps;};
	int& NLS_MaxModNewtonSteps() {return solset.nls_maxmodnewtonsteps;};

	int NLS_MaxRestartNewtonSteps() const {return solset.nls_maxrestartnewtonsteps;};
	int& NLS_MaxRestartNewtonSteps() {return solset.nls_maxrestartnewtonsteps;};

	int NLS_MaxFullNewtonSteps() const {return solset.nls_maxfullnewtonsteps;};
	int& NLS_MaxFullNewtonSteps() {return solset.nls_maxfullnewtonsteps;};

	int NLS_TrustRegion() const {return solset.nls_trustregion;};
	int& NLS_TrustRegion() {return solset.nls_trustregion;};

	double NLS_TrustRegionDiv() const {return solset.nls_trustregiondiv;};
	double& NLS_TrustRegionDiv() {return solset.nls_trustregiondiv;};

	int NLS_SymmetricJacobian() const {return solset.nls_symmetricjacobian;};
	int& NLS_SymmetricJacobian() {return solset.nls_symmetricjacobian;};

	const NumNLSolver& NumSolver() const {return numsol;}
	NumNLSolver& NumSolver() {return numsol;}
	//void RestartNumSolver() {numsol = NumNLSolver(this);}

	//void SetJacCount() {jaccount++;}
	void SetJacCount(int ii) {log.jaccount=ii;}

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//#graphics:
	virtual void DrawSystem() {};
	virtual RenderContext * GetRC() { return pCurrentRC; }
	virtual ControlWindowContext * GetCWC() { return p2DDrawWindow; }
	virtual void MyDrawLineH(const Vector3D& p1, const Vector3D&p2, const Vector3D& vy2, 
		double t, double h, int drawouterface=1) const;

	void DrawLegend(double yoff=0) const;

	void DrawColorQuads(const TArray<Vector3D>& p, const TArray<double>& v, int n1, int n2, 
		int colormode, int drawlines=0, int vres=1); //draw n1*n2 quads with FEcolors, interpolate values v with higher resultion vres

	void DrawQuad(const Vector3D& p1,const Vector3D& p2,const Vector3D& p3,const Vector3D& p4) const;

	void DrawTrig(const Vector3D& p1,const Vector3D& p2,const Vector3D& p3) const;

	void DrawHex(const Vector3D& p1,const Vector3D& p2,const Vector3D& p3,const Vector3D& p4,
		const Vector3D& p5,const Vector3D& p6,const Vector3D& p7,const Vector3D& p8, int drawouterfaces=1) const;
	
	void DrawCube(const Vector3D& p0, const Vector3D& v1, const Vector3D& v2, const Vector3D& v3) const; //draws cube with reference point and 3 axes
	

	virtual void MyDrawRectangle(const Vector3D& p1, const Vector3D& p2, const Vector3D& p3, const Vector3D& p4, double thickness, const Vector3D* colline, const Vector3D* colfill=0) const;

	virtual void MyDrawLine(const Vector3D& p1, const Vector3D&p2, double thickness) const;
	virtual void MyDrawLine(const Vector3D& p1, const Vector3D&p2, double thickness, const Vector3D& col) const; 
	virtual void MyDrawLine(const Vector3D& p1, const Vector3D&p2, double t, double h) const
	{
		MyDrawLineH(p1,p2,p2-p1,t,h);
	}
	virtual void MyDrawCircleXY(const Vector3D p, double r, const Vector3D& col, int res=12, double thickness=1) const;
	virtual void MyDrawArrow(const Vector3D& p1, const Vector3D&p2, const Vector3D& col, 
		double linethickness = -1, double headsize = -1, int resolution=8) const;
	void DrawArrow(const Vector3D& p1, const Vector3D&p2, double linethickness = -1, double headsize = -1, int resolution = 8) const;
	void DrawColorArrow(double v, const Vector3D& p1, const Vector3D&p2, double linethickness = -1, double headsize = -1, int resolution = 8);

	
	virtual void DrawColorZyl(double v, const Vector3D& pz1,const Vector3D& pz2, double rz, int tile);     //draw a zylinder in 3D, with color specified by value v, endpoints p1 and p2, Radius r, number of surfaces (discretized): tile
	virtual void DrawZyl(const Vector3D& pz1, const Vector3D& pz2, double rz, int tile=8) const;
	virtual void DrawCone(const Vector3D& pz1, const Vector3D& pz2, double rz, int tile=8, int drawconelines = 0) const;
	virtual void DrawZyl(const Vector3D& pz1, const Vector3D& pz2, const Vector3D& pz1dir,const Vector3D& pz2dir, double rz, int leftend = 1, int rightend = 1, int tile=8) const;

	//tile ... number of tiles per circumference, fill:1 ...draw whole sphere, fill=0.5 ... draw half sphere
	virtual void DrawSphere(const Vector3D& p1, double r, int tile=8, double fill=1) const;
	virtual void DrawColorSphere(double value, const Vector3D& p, double r, int tile=8, double fill=1);


	virtual void DrawPolygon(const TArray<Vector3D>& p, int drawlines=0, double linewidth=1) const;
	virtual void DrawPolygonOutline(const TArray<Vector3D>& p, double linewidth=1) const;

	virtual void SetColor(const Vector3D& col); //set new color (set actcolor)
	virtual void SetLineColor(const Vector3D& col) {colline = col;}; //set new color (set actcolor)
	virtual const Vector3D& GetLineColor() {return colline;}; //set new color (set actcolor)
	virtual void SetTransparency(int transp_on) {transparency_on = transp_on;};
	virtual int GetTransparency() const {return transparency_on;};
	virtual void ChooseColor(float R, float G, float B) const; //choose color for painting
	virtual void SetDrawlines(int i) {drawlines = i;}
	virtual void SetDrawlinesH(int i) {drawlinesh = i;}
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//!AD: 2012-12-13: functions for 2D Draw window
	virtual void DrawSystem2D() {};
	//virtual void AddDrawComponent(ControlWindowContext::DrawComponent dcomp);
	virtual void AddDrawComponent_Line(int mbs_elnr, TMBSElementType mbs_type, int sub_elnr, TDrawElementType sub_type, Vector2D p1, Vector2D p2, Vector3D col);
	virtual void AddDrawComponent_Rect(int mbs_elnr, TMBSElementType mbs_type, int sub_elnr, TDrawElementType sub_type, Vector2D center, Vector2D size, Vector3D col_border, Vector3D col_background);
	virtual void AddDrawComponent_Ellipse(int mbs_elnr, TMBSElementType mbs_type, int sub_elnr, TDrawElementType sub_type, Vector2D center, Vector2D size, Vector3D col_border, Vector3D col_background);
	virtual void AddDrawComponent_Text(int mbs_elnr, TMBSElementType mbs_type, int sub_elnr, TDrawElementType sub_type, Vector2D center, Vector2D size, Vector3D col_text, mystr& text, TTextAllign positioning);
// functions that allow the 2D Draw Window to change properties of MBS Elements
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	virtual int GetDrawResolution() const {return GetIOption(139);} //depreciated --> use options!

	//to be overwritten in MBS:
	//get values in MBS_EDC:
	virtual void MBS_EDC_TreeGetDouble(double& data, const char* name) const {};
	virtual void MBS_EDC_TreeGetInt(int& data, const char* name) const {};
	virtual const char* MBS_EDC_TreeGetString(const char* name) const {return 0;};
	//set values in MBS_EDC:
	virtual void MBS_EDC_TreeSetDouble(double data, const char* name) {};
	virtual void MBS_EDC_TreeSetInt(int data, const char* name) {};
	virtual void MBS_EDC_TreeSetString(char* data, const char* name) {};


	virtual const double& GetTImincol() const {return TImincol;}
	virtual double& GetFEmincol() {return TImincol;}
	virtual const double& GetTImaxcol() const {return TImaxcol;}
	virtual double& GetFEmaxcol() {return TImaxcol;}
	virtual int IsFlagComputeMinMaxFECol() const {return flag_compute_minmax_FEcol;}
	virtual void SetFlagComputeMinMaxFECol(int flag) {flag_compute_minmax_FEcol = flag;}

	//+++++++++++++++++++++++++++++++++++++++++++++++++
	//element data access:
	virtual int GetNElements() const {return 0;} //overwritten in mbs.h
	virtual void GetElementData(int i, int type, int value, ElementDataContainer& edc) {};
	virtual int SetElementData(int i, int type, int value, ElementDataContainer& edc) {return 0;};

	//$AD 2013-07-08: Manipulate Content of arrays in IOElements from 2D Draw Window
	virtual void InsertIOElemConNode(int elemnr, int list_idx, int input_nr, Vector2D pos) = 0;
	virtual void DeleteIOElemConNode(int elemnr, int list_idx) = 0;

	//function call:
	//virtual int CallCompFunction(int action, int option, int value, ElementDataContainer* edc) {return 0;};
	virtual int CallCompFunction(int action, int option = 0, int value = 0, ElementDataContainer* edc = NULL) {return 0;};

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // access functions for Sensor / PlotTool (Plottool has no classes MBS, Sensor, etc included)
	virtual TArray<double>* GetSensorValuesArrayPtr(int i) =0;        // get the i-th sensors (as registered in mbs) own stored values from internal array
	virtual TArray<double>* GetSensorTimesArrayPtr(int i) =0;         // get the i-th sensors (as registered in mbs) own stored times from internal array
	virtual mystr GetSensorName(int i) {return 0;};                   // get the i-th sensors (as registered in mbs) name
	virtual int GetSensorNumber(mystr& name) {return 0;};             // get the number of the (first) sensor that matches the sensorname

	//+++++++++++++++++++++++++++++++++++++++++++++++++
	virtual Vector3D FEColor(double val) const;
	virtual void SetFEColor(double val) const;
	virtual void DrawIsoQuad(Vector3D* p, const Vector3D& n, double* v); //draw quad with iso-lines
	virtual void DrawIsoQuadTex(Vector3D* p, const Vector3D& n, double* v, int res, char* teximage);

	// update minimum and maximum values by val
	virtual void UpdateFEMinMaxCol(double val);

	virtual void SetTexStoreResolution(int n);
	virtual int GetTexStoreResolution() const {return texstorenn;}
	virtual int GetTexStoreMaxResolution();

	//drawing:
	// int TIdrawmode;				//*YV: commented out as it is not needed any longer
	Vector3D actcolor, colline;
	virtual Vector3D & ColLine() { return colline; }	// for interface access
	int drawlines;
	int drawlinesh;
	int transparency_on; //activate / deactivate transpareny mode, factor is an option

	//texture:
	int TIdrawtexres;
	int texstorenn;  //size of texture, minimum requirement: 64x64 for opengl
	int reduceimage; //reduce size of image, such that only a part of it is used
	char* texImage;  //pointer to texture image
	unsigned int texName;

	//public Time-Integration variables:
	double ada_err_sum;
	double acterror, acterror_damp;
	double comperr;
	double TItimeperstep;
	int TIissubstep; //in fulladaptive, currently computing the 
	double oldenergy;
	double TIcondnum; //rough estimate condition number of jacobian

	//drawing/storing/writing results:
	double laststoredata;
	double lastloadresults;
	int drawnow;

	int system_matrices_written; //flag, system matrices are written at first time step

//graphics:
	RenderContext * pCurrentRC;
	ControlWindowContext* p2DDrawWindow;               //!AD: new 2012-12-13

	double tidraw_offx;
	double tidraw_offy;
	double tidraw_offz;


	//PG:
	//// output in log file
	//virtual ofstream& LogOut() { return logout; }
	ofstream logout;	// log-file
	mystr logout_name;	// log-file name

protected:
	TimeIntLog log;

protected:
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//solver settings, can be accessed in derived class:
	SolverSettings solset;

private:

	virtual void RenderScene(RenderContext * pRC) 
	{	
		pCurrentRC = pRC;
		if (!solset.withgraphics){return;};

		DrawSystem();
	};

	virtual void RenderControlWindow(ControlWindowContext* pCWC)
	{
		p2DDrawWindow = pCWC ;               //!AD: new 2012-12-13
		if (!solset.withgraphics){return;};

		DrawSystem2D();
	}


  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//new structure vor TimeInt variables during computation
	TimeIntVars TIV;
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	
	protected:
	//erase:
	//Vector TIx0last, TIk0last;
	//Vector TINLsubstep_k0save, TINLsubstep_x0save;
	//TIx0old, TIx0save, TIk0old, TIk0save;
	//double TItimesave; //actual simulation time

	//general time ingetration:
	//actual state, old state (NonlinSubStep), saved state (Save/RestoreState), state for drawing
	Vector TIx0;   //actual states for computation: [q v x z], lengths [sos sos es is]
	private:
	Vector TIk0;     //actual states for k-version of time integration
	Vector TIdata;   //data variables for actual step
	double TItime;   //actual time for computation

	Vector TIx0draw;   //actual states for drawing;
	Vector TIdrawdata; //data variables for drawing;
	double TIdrawtime; //actual time for drawing result;

	Vector TIx0start;  //starting vector for time integration and static solver
	Vector TIk0start;  //k-version of starting vector for time integration and static solver

	TIstepsave TIlaststep_state;		//state at end of last successful time step
	TIstepsave TIlastnonlinit_state;//state after last nonlinear iteration
	TIstepsave TItemp_state;				//temporary state (for postprocessing) // AP: AND uses this state for the plastic strains after transport, before nonlinear iteration

	double TIstep; //actual timestepsize
	double TIstepnew; //suggested timestepsize for next timestep
	double TImaxstep; //max timestepsize
	double TIminstep; //min timestepsize


	double TIsteprecommendation; //step suggested by mbs system

	protected:
	double TIcomptime; //elapsed computational time
	double TIactlocerr, TInonls;
	int TIit;
	int TIrejectedsteps;
	int TIdiscontstep;

	private:
	int TIfinished;		//can be set 1 if integration shall be stopped; set -1 if integration stopped with error

	double TImincol;
	double TImaxcol;
	int flag_compute_minmax_FEcol;


	double loadfact; //for static computation, do not apply all load at once ..., for dynamic computation == 1 !!!!
	int computejacflag; //is 1 during computation of jacobian

	//Tableau data, variable stepsize
	TArray<int> stage_tableaus;
	TArray<IRK_Tableau*> tableaus;
	int act_tableau; //actual used table
	int act_stage;

	//variables to control order selection
	int steps_since_orderchange;
	int min_newtonits;
	int newjacobians_since_orderchange;
	int last_jaccount;
	int reduce_order; //reduce as soon as possible

	double TIlastwritesol;
	double TIwritesolevery;
	int TInumber_of_solution_data_steps;   //$ PG 2012-4-13: store number of stored solution data steps

	//old:
	int TInoincs;
	double TIbasestep;
	int impl_equ_flag; //for implicit differential equations
	int nls_output;

protected:
	NumNLSolver numsol;

private:
	int TIkv; 	//K-version integration
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//local vectors/matrices for Jacobian/LocalJacobian
	SparseMatrix mlocgs;
	SparseMatrix mlocf2s;
	SparseMatrix mlocms;

	Matrix mlocf;
	Matrix mlocg;
	Matrix mlocf2;
	Matrix mlocm;
	Matrix mlocm2;
	Matrix diag;
	Matrix TIm; //for K-version NLF
	SparseMatrix TIm_sparse;

	Vector locjacf0,locjacf20,locjacg0; //for non-symmetric jacobian
	Vector locjacf1,locjacf21,locjacg1;
	Vector locjacf2,locjacf22,locjacg2, xx;
	SparseVector slocjacf20, slocjacf21, slocjacf22;
	SparseVector slocjacg0, slocjacg1, slocjacg2;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	int bandwidthm;
	int reducedbandsize; //reduce bandsize-part of the matrix for certain flexible ACRS models connected with rigid bodies
	int TIm_initialized; //initialize only at beginning

	IVector temprowind, tempclearind;
	//local vectors/matrices for GeneralIstep

	//local vectors/matrices for NLF
	Vector fg;
	Vector gv[global_maxstages];
	Vector xv[global_maxstages];
	Vector iv[global_maxstages];
	Vector gu[global_maxstages];
	Vector gue[global_maxstages];
	Vector xi[global_maxstages];
  Vector kue[global_maxstages]; //first order
  Vector ku[global_maxstages];  //second order u
  Vector kv[global_maxstages];  //second order u
	Vector xh;
	Vector xh2;
	Vector evalfs[global_maxstages];
	Vector evalf;
	Vector evalfe;
	Vector evalg;
	Vector TIu0;
	Vector TIv0;
	Vector TIu0e;
	Vector TIi0;
	Vector TIku0;
	Vector TIku0e;
	Vector TIkv0;
	Vector minvf[global_maxstages];
	//Matrix m[global_maxstages];
};

#endif
