//#**************************************************************
//#
//# filename:             mbs.h
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						July 2004
//# description:          Driver and model for timeintegration
//#                       Model of a rigid arm and a hydraulic zylinder
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
 



#ifndef MBS__H
#define MBS__H

#include "timeint.h"
#include "geomelements.h"
#include "ElementsAndModelsLibraryInterface.h"

class CEDCParser;

void TIMBSWarningHandle(const char* warn, int use_instant_message_text=0);

//$!YV 2012-11-28:	MultiBodySystem is now the class itself, and MBS is an interface for its functionality,
//									which may be required by elements and model functions
class MultiBodySystem : public TimeInt
{
public: 
	MultiBodySystem();

	virtual ~MultiBodySystem();

	virtual void Destroy();

	virtual void ClearSystem() 
	{
		Destroy();
	};
	virtual void InitFirst();
	virtual void InitGlobalEDCVariables();

	//# Timeint specific derived functions (must)
	virtual void PerformComputation();
	virtual int ComputeSystem(); //do integration/static computation or eigenmode analysis or other

	virtual void PerformSingleStaticDynamicComputation(); // perform single computation
	virtual void PerformParameterVariation(); // perform parameter variation
	//$!PG 2011-2-22:[DoParameterVariationComputation is needed in PerformParameterVariation()
	virtual bool DoParameterVariationComputation(bool const isgeometric, const double par, const double varstep, const double varstepfact, const double varend) const;
	//$!PG 2011-2-22:]
	//$ YV 2012-11-28[
	//virtual int PerformOptimization(MBS* mbs); // perform optimization of parameters using sensor values as cost function (e.g.: Newton, Genetic)
	//virtual int PerformSensitivityAnalysis(MBS* mbs); // analyze sensitivity of sensor values with respect to parameters
	//$ YV ]
	virtual double EvalSensorCostFunctionVal(); // evaluate cost functions of sensors
	virtual int ComputeEigenmodes(); // compute eigenmodes of system, rv==0 --> ok, rv==1 --> not ok

	virtual int RepeatedlyPerformComputation(int flags=0); //repeatedly call computation (integration) and do pre/postcomputation, sensor evaluation, etc....

	virtual int PreComputationOperations(int flags=0); //flag&1==>always set initial conditions;return 1==OK, return 0==FAILED ==> ABORT simulation
	
	virtual void PostComputationOperations(int integrate_rv); //write log messages, so some sensor computations, etc.
	virtual void DoFinalSensorEvaluations(); //do final sensor evaluations (min/max/diff, etc.) for sensors after computation
	virtual void WriteSolutionVector(); //write solution vector into file if option is set

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// initialization

	virtual void SetInitialConditions();
	virtual void Initialize();
	virtual void InitializeAfterConfigLoaded();

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// load save files
	virtual void OpenFiles(int flag); //flag&1==append
	virtual void OpenLogFile(int flag); //flag&1==append
	virtual bool LogFileNameChanged();
	virtual void CloseLogFile();
	virtual void CloseFiles();
	virtual int SaveSolutionVector(const mystr& filename); //rv=1 --> OK
	virtual int LoadInitialVector(const mystr& filename);  //rv=1 --> OK

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//MBS_communication
	virtual void SetOptions2EDCOptions(); //copy numbered i/d/t options to EDC Tree options 
	virtual void SetEDCOptions2Options(); //copy EDC Tree options to numbered i/d/t options

	virtual void SetSolverDialogOptions();
	virtual void SetComputationSolverOptions();

	virtual void SolverOptions2EDC(ElementDataContainer* edc);
	virtual void EDC2SolverOptions(const ElementDataContainer* edc);
	virtual void EDC2Options(const ElementDataContainer* edc); // implemented in auto - generated file "options_mbs_auto.h"
	virtual void Options2EDC(ElementDataContainer* edc); // implemented in  auto - generated file "options_mbs_auto.h"


	virtual void GetElementData(int i, int type, int value, ElementDataContainer& edc);
	virtual int SetElementData(int i, int type, int value, ElementDataContainer& edc);

	virtual TimeIntLog GetTimeIntLog() const {return log; }

	//$AD 2013-07-08: Manipulate Content of arrays in IOElements from 2D Draw Window
	virtual void InsertIOElemConNode(int elemnr, int list_idx, int input_nr, Vector2D pos);
	virtual void DeleteIOElemConNode(int elemnr, int list_idx);
	//$AD 2013-07-10: Change Position of a single element (MBS element, conNode, ...) 
	virtual void MoveConNode2D(int elemnr, int list_idx, double delta_x, double delta_y);
	virtual void MoveElement(int elemnr, double delta_x, double delta_y, double delta_z = 0.);


	//call function/perform action in MBS-System:
	virtual int CallCompFunction(int action, int option=0, int value=0, ElementDataContainer* edc=NULL);

	virtual int CheckSystemConsistency();//check system; rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	virtual int CheckConsistencyLTG();
	virtual int CheckConsistencyElements();
	virtual int CheckConsistencySensors();
	virtual int CheckConsistencyNodes();
	virtual int CheckConsistencyMaterials();


	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//element and file access:
	virtual void SaveToFile(mystr filename, int save_options=1);
	virtual int LoadFromFile(mystr filename); //return 1 if success, 0 otherwise
	virtual int IsLoadSaveMode() const {return isloadsavemode;}
	//$ YV 2012-12-30: this function was used only in script_parser - moved to script parser
	//virtual void LoadError(const mystr& str) {UO()<< "Error in load file: " << str.c_str() << "\n";}
	virtual void EDCError(const mystr& str) {UO()<< "Error in element data: " << str.c_str() << "\n";}


	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//these functions are called from WCDriver via WinCompInterface:
	virtual MBSObjectFactoryInterface * GetObjectFactory();
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//mbs_edc_options
	virtual ElementDataContainer* GetMBS_EDC_Options();
	virtual const ElementDataContainer* GetMBS_EDC_Options() const;
	virtual void SetMBS_EDC_Options(const ElementDataContainer& edc);

	virtual ElementDataContainer* GetMBS_EDC_Variables();
	virtual const ElementDataContainer* GetMBS_EDC_Variables() const;
	virtual void SetMBS_EDC_Variables(const ElementDataContainer& edc);

	//get values in MBS_EDC:
	virtual void MBS_EDC_TreeGetDouble(double& data, const char* name) const;
	virtual void MBS_EDC_TreeGetInt(int& data, const char* name) const;
	virtual const char* MBS_EDC_TreeGetString(const char* name) const;
	//set values in MBS_EDC:
	virtual void MBS_EDC_TreeSetDouble(double data, const char* name);
	virtual void MBS_EDC_TreeSetInt(int data, const char* name);
	virtual void MBS_EDC_TreeSetString(const char* data, const char* name);
	
	//$ YV 2012-11-28: commented out
	//virtual const MBS* GetMBS() const {return this;}

	virtual const CMBSParser& MBSParser() const;
	virtual CMBSParser& MBSParser();
	virtual const CEDCParser& EDCParser() const;
	virtual CEDCParser& EDCParser();

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//model data container
	virtual ElementDataContainer* GetModelDataContainer() {return edc_modeldata;}
	virtual ElementDataContainer* GetModelDataContainer_args() {return &modeldata_edc_args;}
	virtual void SetModelDataContainer(const ElementDataContainer& edc)
	{
		if (edc_modeldata != 0)
		{ 
			delete edc_modeldata;
			edc_modeldata = 0;
		}

		edc_modeldata = new ElementDataContainer;

		edc_modeldata->CopyFrom(edc);
	}
	virtual void AddReplaceModelDataEDC(/*const */ElementDataContainer& edc);//$ RL 2011-6-16: function for adding or replacing model data edc with data from file
	virtual int ReadModelData(mystr filename); //read model data from file; return 1 in case of success, 0 else
	virtual int File2Str(const char* filename, mystr& str); // this function reads a file into a mystr
	virtual int File2EDC(const char* filename, ElementDataContainer* edc_file); //this function reads a file and stores it in an edc

	virtual int LoadVectorFromFile(const char* filename, int col, ElementDataContainer* return_value);	//$ DR 2013-07-04: reads a column of a file and stores the vector in an edc

	// this function reads a stl-file and writes the data in the edc
	virtual int STLFile2EDC(char* file, int binary, ElementDataContainer* return_value); //DR 2013-01-14 

	//this function computes mass, moment of inertia, volume and center of gravity based on the data about the geometry and the material
	virtual int ComputeInertia(ElementDataContainer* data, ElementDataContainer* return_value); //$ DR 2013-01-30

	//$ YV 2012-12-30: moved this function to femathhelperfunctions.h
	//$ RL 2012-7-25:[
	/*
	virtual int DoesFileExist(mystr filename, int warn = 1);// this function tries to open a file and returns 1 if the file was OK, otherwise 0; if the flag 'warn' is set to 1, an warning text message appears if file was not found
	virtual const int DoesFileExist(mystr filename, int warn = 1) const;// this function tries to open a file and returns 1 if the file was OK, otherwise 0; if the flag 'warn' is set to 1, an warning text message appears if file was not found
	*/
//$ RL 2012-7-25:] 
	virtual int AddModelData(mystr filename); // add/replace model data to edc_model_data from a file; return 1 in case of success, 0 else
	virtual int RemoveVariableFromEDCTree(ElementDataContainer& edc, mystr& varname, int warning); // remove a variable or branch from the edc by name - this is used to prevent ComputationSteps to be written to file and to be remembered from a previous model		
	virtual int IsModelData_Initialized() const {return modeldata_initialized;}
	virtual void SetModelData_Initialized(int flag) {modeldata_initialized = flag;}
	virtual int AddFileSensor(const char* filename, const mystr sensor_name, const int ncolumns, const int colTime, const int colSignal, const int interpolation=0, const int nrOfHeaderLines=0); // 3 elements are created to "measure" data from columns of a file (time->IOMathFunction->MBSSensor), accurate during whole computation, e.g. sensor data is used during computation not only for viewing
	virtual int AddFileSensor2(const char* filename, const mystr sensor_name, const int colTime, const int colSignal, const int interpolation=0, mystr comment = mystr("%"));                     // just one MBSSensor is created to "measure" data from columns of a file, accurate at begin/end of solver time steps, recommended for use e.g. to plot signals from files to compare the solution with the file values

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//access to model functions

	//$ YV 2013-01-02: some of the functions below are no longer necessary - they are replaced
	// by the models library in the elements and models module

	//virtual int GetNumberOfModelFunctions() const ;					     //return number of available model functions
	//virtual char* GetModelFunctionName(int index0) const;        //return model name string
	//virtual char* GetModelFunctionDescription(int index0) const; //return model description string
	//virtual int GetModelFunctionOption(int index0) const;				 //return option of model function
	//virtual int CallModelFunction(int index0);										//Call model function with index "index0"
	virtual int GetModelFunctionsIndex0();

	//virtual int CallModelFunction_InitModelData(int index0);			//Call initialization function for model data
	//virtual int ModelFunctionHasInitModelData(int index0);

	virtual void ModelFunctionChanged();

	// this function returns the current hotint version, as defined in "WorkingModule\hotint_version.h", 
	virtual const HotintVersionInfo& GetHotintVersion() const {return WCDInterface::GetHotintVersion();}

	//$ YV 2013-01-02: models library from the models and elements module
private:
	MBSModelsLibraryInterface * pModelsLibrary;

public:
	MBSModelsLibraryInterface * GetModelsLibrary() { return pModelsLibrary; }
	void SetModelsLibrary(MBSModelsLibraryInterface * pModelsLibrary) { this->pModelsLibrary = pModelsLibrary; }

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//# output and drawing
	virtual void WriteSol();
	virtual void DrawSystem();
	virtual void DrawSystem2D();


	virtual void WriteSolDataInfo();
	virtual int ReadSolDataInfo(); // returns number of stored solution data files

	virtual void DrawBodies();
	virtual void DrawConstraints();
	virtual void DrawSensors();
	virtual void DrawLoads();
	virtual void DrawNodes();

	virtual void PrecomputeEvalFunctions(); //$ PG 2012-1-15: precomputation for EvalF(), EvalF2(), EvalG(), EvalM(), and also for CqTLambda() in case of constraints

	virtual void EvalF(const Vector& x, Vector& f, double t);  //*first order equations: \dot q=f(q,t), return vector:  len(q)
	virtual void EvalG(const Vector& x, Vector& f, double t); //*evaluate constraints: len(z)
	virtual void EvalG(const Vector& x, SparseVector& f, double t, IVector& rowind, IVector& clearind); //only for jacobian, rowind and clearind for speedup, rowind must be initialized with zeros and l=f.Length()

	virtual void EvalM(const Vector& x, Matrix& m, double t);
	virtual void EvalM(const Vector& x, SparseMatrix& m, double t);

	virtual void EvalF2(const Vector& x, Vector& f, double t); //*second order equations: M \ddot u = F2(u,\dot u,t), len(u)
	virtual void EvalMinvF2(const Vector& x, Vector& f, double t); //*second order equations: M \ddot u = F2(u,\dot u,t), len(u)
	virtual void EvalF2(const Vector& x, SparseVector& f, double t, 
		IVector& rowind, IVector& clearind, IVector& elems); //only for jacobian, rowind and clearind for speedup, rowind must be initialized with zeros and l=f.Length()
	virtual void LocalJacobianF2(SparseMatrix& m, Vector& x);
	virtual void StaticJacobianF2(SparseMatrix& m, Vector& x); //new version, which replaces local jacobian f2 for static comptutation
	virtual void LocalJacobianF2(Matrix& m, Vector& x);  //Compute (df2/dx)
	virtual void LocalJacobianG(SparseMatrix& m, Vector& x);
	virtual void ComputeSparseMatrixUsage(IVector& usage_per_dof); //compute entries of each row in sparse matrix

	virtual double GetKineticEnergy();  //compute kinetic energy of system
	virtual double GetPotentialEnergy();//compute potential (strain) energy of system

	virtual int GetImplicitSize() const;
	virtual int GetSecondOrderSize() const;  //*size of second order equations, len(u)
	virtual int GetSecondOrderSize_RS() const;  //*size of second order equations, len(u), after resorting!!!!
	virtual int GetFirstOrderSize() const;  //*size of first order equations
	virtual int GetResortSize() const;  //*size of first order equations
	virtual int GetDataSize() const; //size of data variables for each time step of system

	//# Timeint specific derived functions: for discontinuities
	virtual void StartTimeStep();
	virtual void EndTimeStep();
	//function is called when computation of time step is started
	virtual double PostNewtonStep(double t);
	virtual void PostprocessingStep();
	virtual double GetError();
	
	virtual void ComputationFinished(); //$ MaSch 2013-08-08
	
	//virtual int UseSparse() const {return 0;}
	virtual int ResortActive() const {return resort_active;}
	virtual void SetResortActive(int flag) {resort_active = flag;}
	virtual int ResortMode2() const {return 1;}

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //# Computation Steps with changing SolverOptions...
	//Time mapping
	virtual int ParseStepEndTimes();       // parses all EDC for ComputationSteps entries - implemented in MBS
	virtual int DoubleCheckWithLoadSettings();
	virtual ElementDataContainer* GetComputationStepEDC(int computation_step_number);
	virtual int GetMaxLoadSteps() const;   // maximal number of loadsteps from all loads - there MUST be at least that many computation steps...
  virtual int AddDummyComputationStep(); // additional computationsteps if loadsteps supercede computation steps, interval always 1 sec*
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //# Computation Steps with changing SolverOptions...
  //replace solver options
	virtual int ApplyComputationStepSettings(int stepnumber);  // changes the solver options - implemented in MBS


	virtual int GetJacCol() 
	{
		int jc = NLS_GetJacCol();
		int ns = GetActTableau()->nstages;
		if (ns > 1 && NLS_IsJacFullNewton())
		{
			return 0;
			if (jc > (MBSsos+MBSes)*ns && MBSis != 0)
			{
				jc = (jc-(MBSsos+MBSes)*ns-1)%MBSis+1+MBSsos+MBSes;
			} 
			else if (jc > (MBSsos)*ns && MBSes!=0)
			{
				jc = (jc-(MBSsos)*ns-1)%MBSes+1+MBSsos;
			}
			else if (MBSsos!=0)
			{
				jc = (jc-1)%MBSsos+1;
			}
		}
		return jc;
	}

	//add additional drawing objects to scene
	virtual int Add(const GeomElement& de)
	{
		GeomElement* den = de.GetCopy();
		return drawelements.Add(den);
	}
	//add additional drawing objects to scene
	virtual int Add(const GeomElement& de, int elnum)
	{
		GeomElement* den = de.GetCopy();
		den->SetElnum(elnum);
		return drawelements.Add(den);
	}
	virtual GeomElement* GetDrawElement(int i) const {return drawelements(i);}
	virtual GeomElement* GetDrawElement(int i) {return drawelements(i);}
	virtual int NDrawElements() const {return drawelements.Length();}
	virtual void DeleteDrawElement(int i);

	virtual void DrawElementAdd();

	virtual int GetProhibitRedraw() const {return prohibit_redraw;}
	virtual void SetProhibitRedraw(int i) {prohibit_redraw = i;}

	virtual int MaxIndex() const {return maxindex;};
	virtual void SetMaxIndex(int i) {maxindex=i;};
	virtual void AddControllerSensors(); // add sensors for all controller objects
	//assemble system
	//assign degrees of freedom
	virtual void Assemble();
	virtual void BuildLTGLists();
	virtual void LabelIOElementOutputs(); //$ AD: new 2013-02-01 see Log 388
	virtual int FilterIOElements(); //$ AD: new 2013-02-25 see log 4??
	virtual void ComputeResortVector();
	virtual int UseDependencies() const {return use_dependencies;};
	virtual void SetUseDependencies(int flag) {use_dependencies = flag;};
	virtual int UseJacobianElementWise() const {return GetSolSet().element_wise_jacobian;}

	//element access via graphical interface or for storage:
	virtual int GetNElements() const {return NE();}

	//add a new element to the system
	virtual int AddElement(Element* e);
	//add a new body to the system, bodies are sorted such that they are before the joints
	//virtual int AddBody(Element* e); //$ DR 2013-05-23 removed deprecated function
	//insert an element at position i -> move elements backwards
	//virtual void InsertElement(int i, Element* e);	//$ DR 2013-05-23 removed deprecated function
	//add a element to the system which is not computed, only for drawing and solution
	virtual int AddAuxElement(Element* e);

	virtual void DeleteNode(int i);
	//add a new node to the system
	virtual int AddNode(Node* n); 
	//add new node to system if it does not exist in body n.body
	virtual int AddBodyNode(Node* n); 
	//add new node to system if it does not exist in searchtree
	virtual int AddNode(Node* n, SearchTree& tree);
	//add new node to system if it does not exist in searchtree and it does not exist in body n.body
	virtual int AddBodyNode(Node* n, SearchTree& tree);
	//
	virtual Node& GetNode(int i) 
	{
		if (i > 0 && i <= NNodes()) return *nodes(i); //$ DR 2013-02-22 added if condition and error handling
		else 
		{
			mystr error_msg = mystr("FATAL ERROR: You tried to access the node with number ") +mystr(i) + mystr(", which does not exist! HOTINT will close now.");			
			UO().InstantMessageText(error_msg);
			assert(0);
			exit(0);
		}
	}
	//
	virtual const Node& GetNode(int i) const 
	{
		if (i > 0 && i <= NNodes()) return *nodes(i); //$ DR 2013-02-22 added if condition and error handling
		else 
		{
			mystr error_msg = mystr("FATAL ERROR: You tried to access the node with number ") +mystr(i) + mystr(", which does not exist! HOTINT will close now.");
			UO().InstantMessageText(error_msg);
			assert(0);
			exit(0);
		}
	}
	//
	virtual int NNodes() const {return nodes.Length();}
	//$ DR added the following 2 functions:
	virtual Node* GetNodePtr(int i) 
	{
		if (i > 0 && i <= NNodes()) return nodes(i);	//$ DR 2013-02-22 update of the error handling
		else 
		{
			mystr error_msg = mystr("FATAL ERROR: You tried to access the node with number ") +mystr(i) + mystr(", which does not exist! HOTINT will close now.");
			UO().InstantMessageText(error_msg);
			assert(0);
			exit(0);
		}
	}
	//
	virtual const Node* GetNodePtr(int i) const
	{
		if (i > 0 && i <= NNodes()) return nodes(i);	//$ DR 2013-02-22 update of the error handling
		else 
		{
			mystr error_msg = mystr("FATAL ERROR: You tried to access the node with number ") +mystr(i) + mystr(", which does not exist! HOTINT will close now.");
			UO().InstantMessageText(error_msg);
			assert(0);
			exit(0);
		}
	}
	//
	virtual void DeleteElement(int i);
	//
	virtual Element& GetElement(int i) //$ DR 2013-02-22 added if condition and error handling
	{	
		if (i > 0 && i <= NE())	return *elements(i);
		else 
		{
			mystr error_msg = mystr("FATAL ERROR: You tried to access the element with number ") +mystr(i) + mystr(", which does not exist! HOTINT will close now.");
			UO().InstantMessageText(error_msg);
			assert(0);
			exit(0);
		}
	}
	//
	virtual const Element& GetElement(int i) const
	{	
		if (i > 0 && i <= NE())	return *elements(i); //$ DR 2013-02-22 added if condition and error handling
		else 
		{
			mystr error_msg = mystr("FATAL ERROR: You tried to access the element with number ") +mystr(i) + mystr(", which does not exist! HOTINT will close now.");
			UO().InstantMessageText(error_msg);
			assert(0);
			exit(0);
		}
	}
	//
	virtual Element* GetElementPtr(int i) 
	{
		if (i > 0 && i <= NE()) return elements(i); //$ DR 2013-02-22 update of the error handling
		else 
		{
			mystr error_msg = mystr("FATAL ERROR: You tried to access the element with number ") +mystr(i) + mystr(", which does not exist! HOTINT will close now.");
			UO().InstantMessageText(error_msg);
			assert(0);
			exit(0);
			//return 0;
		}
	}
	//
	virtual const Element* GetElementPtr(int i) const 
	{
		if (i > 0 && i <= NE()) return elements(i); //$ DR 2013-02-22 update of the error handling
		else
		{
			mystr error_msg = mystr("FATAL ERROR: You tried to access the element with number ") +mystr(i) + mystr(", which does not exist! HOTINT will close now.");
			UO().InstantMessageText(error_msg);
			assert(0);
			exit(0);
			//return 0;
		}
	}
	//
	virtual Element& GetAuxElement(int i) 
	{
		if (i > 0 && i <= NAuxE()) return *auxelements(i); //$ DR 2013-02-22 added if condition and error handling
		else 
		{
			mystr error_msg = mystr("FATAL ERROR: You tried to access the aux-element with number ") +mystr(i) + mystr(", which does not exist! HOTINT will close now.");
			UO().InstantMessageText(error_msg);
			assert(0);
			exit(0);
		}
	}
	//
	virtual const Element& GetAuxElement(int i) const 
	{
		if (i > 0 && i <= NAuxE()) return *auxelements(i); //$ DR 2013-02-22 added if condition and error handling
		else 
		{
			mystr error_msg = mystr("FATAL ERROR: You tried to access the aux-element with number ") +mystr(i) + mystr(", which does not exist! HOTINT will close now.");
			UO().InstantMessageText(error_msg);
			assert(0);
			exit(0);
		}
	}
	//
	virtual Element* GetAuxElementPtr(int i) 
	{
		if (i > 0 && i <= NAuxE()) return auxelements(i);  //$ DR 2013-02-22 update of the error handling
		else 
		{
			mystr error_msg = mystr("FATAL ERROR: You tried to access the aux-element with number ") +mystr(i) + mystr(", which does not exist! HOTINT will close now.");
			UO().InstantMessageText(error_msg);
			assert(0);
			exit(0);
			//return 0;
		}
	}
	//
	virtual const Element* GetAuxElementPtr(int i) const 
	{
		if (i > 0 && i <= NAuxE()) return auxelements(i); //$ DR 2013-02-22 update of the error handling
		else 
		{
			mystr error_msg = mystr("FATAL ERROR: You tried to access the aux-element with number ") +mystr(i) + mystr(", which does not exist! HOTINT will close now.");
			UO().InstantMessageText(error_msg);
			assert(0);
			exit(0);
			//return 0;
		}
	}
	//

	virtual int NE() const {return elements.Length();}
	virtual int NAuxE() const {return auxelements.Length();}
	virtual int NConstraints() const;


	virtual const double& GetXact(int i) const
	{
		return xact(i);
	}
	virtual double& GetXact(int i)
	{
		return xact(i);
	}

	virtual const Vector& GetXact() const {return xact;};
	//virtual Vector& GetXact() {return xact;};
	//get actual state from timeint
	virtual void SetActState(const Vector& x) 
	{
		xact.LinkWith(x);
	}
	virtual void SetActState(double* ptr, int len) 
	{
		xact.LinkWith(ptr, len);
	}
	virtual void SetGlobalInitConditions(Vector& x0);

	virtual void SetGlobalInitData(Vector& xdata);

	virtual double GetMagnifyYZ() const {return magnify_yz;};
	virtual void SetMagnifyYZ(double amagyz) {magnify_yz = amagyz;}

	virtual const IVector& GetResortVector() const {return resortvector;}
	virtual IVector& GetResortVector() {return resortvector;}


	virtual void ComputeMaxSparseBandwidth(); 
	virtual int MaxSparseBandwidth() const {return elementbandwidth;} //default value

	virtual void ComputeGyroscopicMatrix(SparseMatrix& gy); // $ MSax 2013-07-25 : added

	virtual void SetTransformJacApply(int i) {transformJacApply = i;};
	virtual int TransformJacApply() const;
	virtual void ApplyTransform(const Vector& v, Vector& Av, int mode); //compute Av=A^T*v in mode==0 and Av=A*v in mode==1

	//useroutput
	mutable UserOutput uout;

	virtual Box3D GetBoundingBoxD() const;
	virtual float GetSceneMaxAbsCoordinate() { return (float)(GetBoundingBoxD().Radius() + GetBoundingBoxD().Center().Norm());} //for zoom all

	virtual void SetCenterObject(int co, const Vector3D& offset) {centerobject = co; centerobject_offset = offset;} //camera moves with object
	virtual int GetCenterObject() const {return centerobject;}
	virtual Vector3D GetCenterObjectOffset() const {return centerobject_offset;}
	virtual double CharacteristicLength() const {return 1;} //scale dimensions, if atomistic (1e-10) or galactic dimensions (1e10)

	virtual int TestCnt() const {return log.testcnt;}

	//solver settings for 
	//$ YV 2012-11-29
	//virtual const SolverSettings& MBSSolSet() const {return solset;}
	//virtual SolverSettings& MBSSolSet() {return solset;}

	virtual int AddSensor(Sensor * s);
	//
	virtual Sensor & GetSensor(int i)
	{
		if (i > 0 && i <= NSensors())	return *sensors(i);	//$ DR 2013-02-22 added if condition and error handling
		else 
		{
			mystr error_msg = mystr("FATAL ERROR: You tried to access the sensor with number ") +mystr(i) + mystr(", which does not exist! HOTINT will close now.");
			UO().InstantMessageText(error_msg);
			assert(0);
			exit(0);
		}
	}
	//
	virtual const Sensor & GetSensor(int i) const 
	{
		if (i > 0 && i <= NSensors())	return *sensors(i);	//$ DR 2013-02-22 added if condition and error handling
		else 
		{
			mystr error_msg = mystr("FATAL ERROR: You tried to access the sensor with number ") +mystr(i) + mystr(", which does not exist! HOTINT will close now.");
			UO().InstantMessageText(error_msg);
			assert(0);
			exit(0);
		}
	}
	//
	virtual Sensor* GetSensorPtr(int i) 
	{
		if (i > 0 && i <= NSensors()) return sensors(i);	//$ DR 2013-02-22 update of the error handling
		else 
		{
			mystr error_msg = mystr("FATAL ERROR: You tried to access the sensor with number ") +mystr(i) + mystr(", which does not exist! HOTINT will close now.");
			UO().InstantMessageText(error_msg);
			assert(0);
			exit(0);
			//return 0;
		}
	}

	virtual const Sensor* GetSensorPtr(int i) const 
	{
		if (i > 0 && i <= NSensors()) return sensors(i);  //$ DR 2013-02-22 update of the error handling
		else 
		{
			mystr error_msg = mystr("FATAL ERROR: You tried to access the sensor with number ") +mystr(i) + mystr(", which does not exist! HOTINT will close now.");
			UO().InstantMessageText(error_msg);
			assert(0);
			exit(0);
			//return 0;
		}
	}


	virtual int NSensors() const {return sensors.Length();}
	//
	virtual void DeleteSensor(int i);

	//$ DR 2012-10:[ loads moved from element to mbs
	virtual int AddLoad(const MBSLoad& li);
	virtual MBSLoad& GetLoad(int i) 
	{
		if (i > 0 && i <= NLoads())	return *loads(i);	//$ DR 2013-02-22 added if condition and error handling
		else 
		{
			mystr error_msg = mystr("FATAL ERROR: You tried to access the load with number ") +mystr(i) + mystr(", which does not exist! HOTINT will close now.");
			UO().InstantMessageText(error_msg);
			assert(0);
			exit(0);
		}
	}
	virtual const MBSLoad& GetLoad(int i) const 
	{
		if (i > 0 && i <= NLoads())	return *loads(i);	//$ DR 2013-02-22 added if condition and error handling
		else 
		{
			mystr error_msg = mystr("FATAL ERROR: You tried to access the load with number ") +mystr(i) + mystr(", which does not exist! HOTINT will close now.");
			UO().InstantMessageText(error_msg);
			assert(0);
			exit(0);
		}
	}
	virtual void DeleteLoad(int i);
	virtual int NLoads() const {return loads.Length();}
	virtual MBSLoad* GetLoadPtr(int i) 
	{
		if (i > 0 && i <= NLoads()) return loads(i); //$ DR 2013-02-22 update of the error handling
		else 
		{
			mystr error_msg = mystr("FATAL ERROR: You tried to access the load with number ") +mystr(i) + mystr(", which does not exist! HOTINT will close now.");
			UO().InstantMessageText(error_msg);
			assert(0);
			exit(0);
			//return 0;
		}
	}
	virtual const MBSLoad* GetLoadPtr(int i) const 
	{
		if (i > 0 && i <= NLoads()) return loads(i); //$ DR 2013-02-22 update of the error handling
		else 
		{
			mystr error_msg = mystr("FATAL ERROR: You tried to access the load with number ") +mystr(i) + mystr(", which does not exist! HOTINT will close now.");
			UO().InstantMessageText(error_msg);
			assert(0);
			exit(0);
			//return 0;
		}
	}
	//$ DR 2012-10:] loads moved from element to mbs


	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // access functions for Sensor / PlotTool (Plottool has no classes MBS, Sensor, etc included)
	virtual TArray<double>* GetSensorValuesArrayPtr(int i);  // get the i-th sensors (as registered in mbs) own stored values from internal array
	virtual TArray<double>* GetSensorTimesArrayPtr(int i);   // get the i-th sensors (as registered in mbs) own stored times from internal array
	virtual mystr GetSensorName(int i);                      // get the i-th sensors (as registered in mbs) name
	virtual int GetSensorNumber(mystr& name);                // get the number of the (first) sensor that matches the sensorname

	//add sensor for the 6 rigid body degrees of freedom:
	virtual void AddRigidBodyDOFSensor(int bodynum, const Vector3D& pos_rel, mystr general_str);

	//add sensor for the 3 rigid 2d body degrees of freedom:
	//virtual void AddRigid2DBodyDOFSensor(int bodynum, mystr general_str, int usepos = 1, int usevel = 0, int displacement = 0); // old version
	virtual void AddRigid2DBodyDOFSensor(int bodynum, Vector2D& pos_rel, mystr general_str, int usepos = 1, int usevel = 0, int displacement = 0);

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//
	virtual int AddMaterial(const Material& m);
	//
	virtual Material& GetMaterial(int i) 
	{
		if (i > 0 && i <= NMaterials()) return *materials(i); //$ DR 2013-02-22 added if condition and error handling
		else 
		{
			mystr error_msg = mystr("FATAL ERROR: You tried to access a material or beam property, which does not exist! HOTINT will close now.");
			UO().InstantMessageText(error_msg);
			assert(0);
			exit(0);
		}
	}
	//
	virtual const Material& GetMaterial(int i) const 
	{
		if (i > 0 && i <= NMaterials()) return *materials(i); //$ DR 2013-02-22 added if condition and error handling
		else 
		{
			mystr error_msg = mystr("FATAL ERROR: You tried to access a material or beam property, which does not exist! HOTINT will close now.");
			UO().InstantMessageText(error_msg);
			assert(0);
			exit(0);
		}
	}
	//
	virtual int NMaterials() const {return materials.Length();}
	//
	virtual int DeleteMaterial(int i); //returns 1 if successful
	//
	virtual Material* GetMaterialPtr(int i) 
	{
		if (i > 0 && i <= NMaterials()) return materials(i);  //$ DR 2013-02-22 update of the error handling
		else 
		{
			mystr error_msg = mystr("FATAL ERROR: You tried to access a material or beam property, which does not exist! HOTINT will close now.");
			UO().InstantMessageText(error_msg);
			assert(0);
			exit(0);
			//return 0;
		}
	}
	//
	virtual const Material* GetMaterialPtr(int i) const 
	{
		if (i > 0 && i <= NMaterials()) return materials(i);  //$ DR 2013-02-22 update of the error handling
		else 
		{
			mystr error_msg = mystr("FATAL ERROR: You tried to access a material or beam property, which does not exist! HOTINT will close now.");
			UO().InstantMessageText(error_msg);
			assert(0);
			exit(0);
			//return 0;
		}
	}

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	virtual const double& GetVarparameter() const {return varparameter;};
	virtual double& GetVarparameter() {return varparameter;};
	virtual const double& GetVarparameter2() const {return varparameter2;};
	virtual double& GetVarparameter2() {return varparameter2;};

	virtual void SetPtr2PerformComputation_Function(int(*ptr2PerformComputation_FunctionI)(MBS*)) {ptr2PerformComputation_Function = ptr2PerformComputation_FunctionI;}
	virtual int CallPtr2PerformComputation_Function(MBS* mbs) {return ptr2PerformComputation_Function(mbs);}

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//Methods used by some old-style models (don't use anymore!)
	virtual void SetShowelem(const int __se) { showelem  = __se; }  // depreciated
	virtual void SetShowelem2(const int __se) { showelem2 = __se; }  // depreciated
	virtual void SetShowelem3(const int __se) { showelem3 = __se; }  // depreciated
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//for animation driver:
	virtual double GetObjData(int row, int col) const 
	{
		int r = row+1;
		if (r < 1) r = 1;
		if (r > objdata.Getrows()) r = objdata.Getrows();

		return objdata(r,col);
	}
	virtual double GetObjDataStepSize() const {return objdata_stepsize;}

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//general solution file: header functions
	virtual int WriteSolutionFileHeader() const {return GetSolSet().sol_write_sol_file_header && !isSolutionFileHeaderWritten;}
	virtual int WriteSolParFileHeader() const {return GetSolSet().sol_write_sol_file_header && !isSolParFileHeaderWritten;}
	virtual mystr GetSolutionFileHeaderString(); // returns the header string
  virtual mystr GetSolParFileHeaderString(); // returns the header string for SolParFile
	virtual mystr GetSolutionDataInfoFileHeaderString(); // returns the header string for the solution-data info file

	virtual ofstream& SolFile() {return sol;}
	virtual ofstream& SolParFile() {return solpar;}
	virtual ofstream& LogFile() {return logout;}

	//* YV added
	//{
	// actually selected field variable for contour plotting (NULL if body color is to be painted)
	FieldVariableDescriptor * GetActualPostProcessingFieldVariable();
	virtual double GetEigenValue(int index) const // eigenvalues are stored in this vector
	{
		if(eigval.Length() == 0)
		{
			return 0.; // if not computed, set to zero
		}
		return eigval(index);
	} 
	virtual int DoFinalEigenValueComputation()
	{
		return mbs_edc_options->TreeGetInt("SolverOptions.Solution.Sensor.postproc_compute_eigenvalues",0);
	}
	virtual void StopByElement(){ simulationStatus.SetStatusFlag(TSimulationStoppedByElement); TIFinished(); } //$ MS + RL 2012-2-29: //$ MsSch 2013-08-19 
	virtual SimulationStatus & GetSimulationStatus() { return simulationStatus; } //$ MaSch 2013-08-19
private:
	int indexOfActualPostProcessingFieldVariable;		// is synchronized with TOption(107)

public:
	virtual int GetIndexOfActualPostProcessingFieldVariable() { return indexOfActualPostProcessingFieldVariable; }
	virtual void SetIndexOfActualPostProcessingFieldVariable(int index);

	//$ MS + RL 2012-2-29: virtual int IscoreComputationFinished() { return coreComputationFinished; } //$ MS 2011-3-16: 
	//}

	//$ AD: new 2013-02-25: functions that let MBS respond to a key or mouse input ( triggered by IOBlocks or GLDrawWindow )
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//functions for response on Input by key or mouse          
	virtual int RespondToKey(int key);
	//virtual void RespondToMouse(double x, double y);

private:
	//actual solution vector (mirror from timeint)
	Vector xact;

	//storage for elements, nodes, sensors, etc.:
	TArray<Element*> elements;
	TArray<Element*> auxelements;
	TArray<Node*> nodes;
	TArray<Sensor*> sensors;		//$ YV 2012-06: Sensor base class is used instead of MBSSensor
	TArray<GeomElement*> drawelements;
	TArray<Material*> materials;
	TArray<MBSLoad*> loads;			//$ DR 2012-10: loads moved from element to mbs
	TArray<int> ioelements;    //$ AD: shortlist containing the elementnumbers of IOElements ( function FilterIOElements() )

	// YV, 17.11.10
	// list of field variables, which are available for post-processing (contour plotting) and for the sensors
	TArray<FieldVariableDescriptor> availableFieldVariables;
	// the list is built by asking all the existing elements at the assembly stage (is called from Assemble())
	void CollectAvailableFieldVariables();
	// the list should be available to the user interface
	virtual const TArray<FieldVariableDescriptor> & GetAvailableFieldVariables() { return availableFieldVariables; }

	//temporary vectors:
	Vector tempf;   //for LocalJacobianF2
	Matrix tempm;   //for LocalJacobianF2
	IVector colref; //for LocalJacobianF2 and G

	//+++++++++++++++++
	int maxindex;  // maximum index of a set of DAE == maximum number of necessary differentiations such that the resulting system contains no algebraic equations
	int elementbandwidth;

	int MBSsos, MBSes, MBSis, MBSss;
	int resortsize_add;
	int resort_active;
	int transformJacApply; //transformation of Jacobian according to ACRS formulation
	int isfirstcomputation;
	int isinconsistent; //0 --> OK, 1 --> can not compute, 2 --> can not draw and not compute
	int use_dependencies;
	int isloadsavemode;

	//time integration parameters:
	//SolverSettings mbssolset;

	//sensor:
	int showelem;  //depreciated
	int showelem2; //depreciated
	int showelem3; //depreciated

	//drawing
	double magnify_yz; // in order to make thin lines visible!!!
	int prohibit_redraw; //do not redraw if this flag is set

	int centerobject;							//use this element to center the computation
	Vector3D centerobject_offset; //use this offset relative to object center

	//cutting planes  //$ PG 2012-4-8:[ Cutting plane may also be handled by OpenGL directly (if ViewingOptions.CuttingPlane.use_open_gl is set to true) 
	TArray<Vector3D> ncut;    //normal vector of cutting plane i (i=1,..,n_cutting_planes   -   cut objects which are in front of plane w.r.t. ncut)
	TArray<Vector3D> pcut;    //distance of cutting plane i from origin multiplied with normal vector
	IVector use_cutting_plane;   //switch on(1)/off(0) cutting plane i
public:
	int UseCuttingPlanes();   //sets ncut, pcut, use_cutting_plane according to EDC-options, and returns if any of the n_cutting_planes is indeed used
	int CuttingPlanesAllow(const Vector3D& pos);   // returns, if object (with golbal reference position pos) is allowed for drawing due to cutting planes

private:
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	double varparameter;
	double varparameter2;
	int (*ptr2PerformComputation_Function)(MBS*); //set this pointer, in order to make a user/customized computation function
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	Vector eigval; // eigenvalues are stored in this vector [rad/s]
	//+++++++++++++++++++++++++++++++++++
	//file in/output:
public:
	ofstream sol; //general solution file
	ofstream solpar; //general solution file for parametric investigations
private:
	ofstream stateout; //file to save last step of computation
	ifstream statein;	 //file to open first step

	Vector stored_initialconditions; //if these initial conditions are set, 
	//																 at every run of the computation, the initial conditions are set


	//animation driver: (special viewer mode, not generally used in HOTINT)
	Matrix objdata;						//AnimationFile, contains objdata for n objects
	double objdata_stepsize;	//AnimationFile stepsize

	ElementDataContainer* edc_modeldata; //use this container to modify/load/save your model data!
	ElementDataContainer* mbs_edc_options; //these options are a mirror of the old options in ti_misc.cpp; only edit these options; they are copied to the timeint i/d/toptions before simulation start

	ElementDataContainer* mbs_edc_global_variables; //this container contains global variables: e.g. "pi", physical constants; evtl. material constants

	ElementDataContainer modeldata_edc_args; //arguments for modeldata from hotint startup
	ElementDataContainer mbs_edc_args; //arguments for mbsedc from hotint startup

	//++++++++++++++++++++++++++++++++++++++++
	//parsing:
	CMBSParser* mbs_parser; //this is a mathematical object parser, which parses mathematical objects (e.g. "1+sin(2*pi*x)"); the evaluation of parsed objects is efficient
	
	CEDCParser* edc_parser; //this is the parser for EDC-input files for HOTINT, which contains hierarchical data structures, commands and can contain parsable mathematical objects

	int modeldata_initialized; //if model is called first time-->initialize model data before, afterwards reused and can be changed in windows interface
	int isSolutionFileHeaderWritten; // is set to 1 when header is already written ==> transfer this variable into mbsedc
	int isSolParFileHeaderWritten; // is set to 1 when header is already written ==> transfer this variable into mbsedc

	//$ MS + RL 2012-2-29: int coreComputationFinished; //$ MS 2011-3-16: is set to 1 when Time-Integration is finished and set to 0 when Time-Integration is started
	//TSimulationStatus simulationStatus; //$ YV 2012-11-28
	SimulationStatus simulationStatus; //$ MaSch 2013-08-19

	//$ YV 2012-12-12: added
	virtual double EvaluateParsedFunction1D(const mystr & parsedFunctionExpression, const mystr & parsedFunctionVariables, double t);
	virtual double ExpressionToDouble(mystr & expression);

	//$ DR 2013-01-14: added as new interface to script parser
	virtual int ExecuteParserCommand(mystr & command, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);

//$ YV 2013-01-02: integration rules live now with the elements; the library is a static member of FiniteElementGeneric
/*
private:
	IntegrationRulesLibrary integrationRulesLibrary;
public:
	IntegrationRulesLibrary * GetIntegrationRulesLibrary() { return &integrationRulesLibrary; }
*/

	//$ MaSch 2013-08-26
	//++++++++++++++++++++++++++++++++++++++++
	//TCP/IP:
	private:
		TCPIPHotInt* serversocket;
	public:
		virtual TCPIPHotInt* & GetServerSocket() {return serversocket;}

};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//$JG2012-01-27: put the following into CEDCParser class!
//transform edc into string, works already for TreeEDCs!
//void EDC2String(const ElementDataContainer& edc, mystr& str, const mystr& indent = "");
//transforms string into edc until a keyword is read, returns remaining string:
//int String2EDC(MBS* mbs, mystr& str, ElementDataContainer& edc, ElementDataContainer& testedc); 
//
////transforms a string automatically into a tree edc; only accepts: double, string (""), vector, matrix
//int String2TreeEDC(MBS* mbs, mystr str, ElementDataContainer& edc); 
//int String2TreeEDC_SIM(MBS* mbs, mystr& str, ElementDataContainer& edc); 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#endif
