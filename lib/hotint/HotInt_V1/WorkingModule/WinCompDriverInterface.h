//#**************************************************************
//# filename:             WinCompDriverInterface.h
//#
//# author:               Yury Vetyukov, Johannes Gerstmayr
//#
//# generated:						2003
//# description:          
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
 
/*
This is an interface class between an OS-independent working (computational) module
of an application and an OS-dependent driving module of an application.

Architecture:
The working module of the program contains a class, derived from WCDInterface.
A function creating dynamically an object of this kind (working object)
WCDInterface * CreateWCDObject();
must be available to the driving module. A pointer to the created working object
is returned. In Windows the working module can be compiled into a standalone
DLL, which exports the function CreateWCDObject().

The driving module is linked (statically or dynamically) to the working module.
After the application is initialized, the function CreateWCDObject() is called
and the main object of the working module is created.

The interface WCDInterface allows direct calls, when something is needed from
the working module. The working module performes OS-dependent actions via
the feed-back interfaces, provided by the driving module.

The multithreading of the application is handled in the driving module.
Initially only the user interface thread exists. As the computation starts,
the working thread is initiated.

---------------         ----------------
|   working   |         | ------------ |
|   object    |-------->| |feed-back | |
|             |         | |interfaces| |
| ----------- |  calls  | ------------ |
| |interface| |<--------|   driving    |
| ----------- |         |   module     |
---------------         ----------------
*/

//#include "tarray.h"
//#include "mystring.h"
//#include "elementdata.h"
//#include "FieldVariableDescriptor.h"
#include "mbs_interface.h"
#include "ElementsAndModelsLibraryInterface.h"
#include "hotint_version_info.h"


#ifndef __WCD_INTERFACE_H_INCLUDED__
#define __WCD_INTERFACE_H_INCLUDED__


#define WORKING_OBJECT_CREATION_FUNCTION_NAME "CreateWCDObject"
#ifdef COMPILE_AND
#define WORKING_DLL_NAME "WorkingModuleHOTINT.dll"
#else
#define WORKING_DLL_NAME "WorkingModule.dll"
#endif
//#define HOTINT_PUBLIC_DOMAIN




//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// This version of the interface is designed for the 3D scenes

#include "renderContext.h"

struct WCDInterface
{
	virtual ~WCDInterface()
	{
	}

	// FEED-BACK INTERFACES
	// these interfaces are provided by the driving module in order to perform
	// the OS-dependent actions. The working module operates with pointers to these objects,
	// which are implemented in the driving module.

	// bringing INFORMATION to the user
	struct UserInterface
	{
		virtual ~UserInterface() {}
		// message box; does not return until the user responds
		virtual void InstantMessageText(const char * ) = 0;
		// set text on something like status bar
		virtual void StatusText(const char * ) = 0;
		// log-file or text running on the screen like in the console window; '\n' means end-of-line
		virtual void AddText(const char * ) = 0;
		// if the text-running window has been switch off, this command will switch it on
		virtual void AssureTextIsVisible() = 0;
		
		// create directory (only, if it does not exist).
		// returns 0 in case of error (wrong path), 1 in case it already existed, and 2 if created
		virtual int CreateMissingDirectory(const char * ) = 0;

		//call a function in WCDriver:
		virtual int CallWCDriverFunction(int action, int option = 0, int value = 0, ElementDataContainer* edc = NULL) = 0;

		//go to sleep for x milliseconds
		virtual void SleepX(int xMilliseconds) = 0;

		//is used in MBS::PreComputationOperations
		virtual int GetModelModified() const = 0;
	};

	// these methods are available as the COMPUTATION procedure starts
	struct ComputationFeedBack
	{
		// informs the driving module that the scene must be repainted
		virtual void ResultsUpdated(int flag=0) = 0;
		// informs the driving module that the computation is paused
		virtual void PausedComputation() = 0;
		// informs the driving module that the computation is finished
		virtual void FinishedComputation() = 0;
		// AH: informs the driving module that new data is sent to plottool
	};

	// DIRECT calls from the driving module to the working object
	// these methods must be implemented in the working module

	// the driving module will give the working module a pointer to the
	// UserInterface interface immediately after the working object is created;
	// and the text output can be done starting from this point
	virtual void SetUserInterface(UserInterface* ) = 0;
	virtual UserInterface* GetUserInterface() const = 0;

	// the maximal absolute value of a coordinate in the scene
	// the function is called once when the graphic system is initialized
	virtual float GetSceneMaxAbsCoordinate() = 0;
	// if the scene is in fact plane, it might be useful
	// to prevent the user from rotating it
	virtual int AllowRotation() = 0;
	// tells the working module to execute all the painting instructions
	// for the currently provided RenderContext
	virtual void RenderScene(RenderContext * pRC) = 0;

  // Sets the Pointer to the Dialog in the 
	virtual void RenderControlWindow(ControlWindowContext* pCWC) = 0;    //!AD: 2012-12-13

	// prints data to the UserInterface context on the user request
	virtual void PrintData(const char * pDataSpecification) = 0;

	// start the computation from the current position
	// this function is run in a separate (working) thread
	// and should return only when the computation is stopped
	// the ComputationFeedBack interface should be used to inform the driving module
	// about the progress in the computation; the UserInterface can be also used
	virtual void Go(ComputationFeedBack * pCFB) = 0;
	virtual void SetPCFB(ComputationFeedBack * p_cfb) = 0;


	//call initialize function after config file loaded:
	virtual void InitializeAfterConfigLoaded() = 0;

	// check if the computation is still in progress
	virtual int IsComputationInProgress() = 0;

	// tells to stop the computation when possible
	virtual void StopWhenPossible() = 0;

	// reset computation, new computation is started with new initial conditions at starttime
	virtual void ResetComputation() = 0;

	// returns if calculation is paused
	virtual int IsPaused() = 0;

	// tells to pause the computation
	virtual void Pause() = 0;
	
	// tells to resume the computation
	virtual void Resume() = 0;

	// print timing statistics
	virtual void PrintTimingList() = 0;

	//old options:
	virtual void SetIOption(int index, int data) = 0; 
	virtual const int& GetIOption(int index) const = 0; 
	virtual int& GetIOption(int index) = 0; 
	virtual void SetDOption(int index, double data) = 0; 
	virtual const double& GetDOption(int index) const = 0; 
	virtual double& GetDOption(int index) = 0; 
	virtual void SetTOption(int index, const char* data) = 0; 
	virtual const char* GetTOption(int index) const = 0; 

	//new options:
	//get values in MBS_EDC:
	virtual void MBS_EDC_TreeGetDouble(double& data, const char* name) const = 0;
	virtual void MBS_EDC_TreeGetInt(int& data, const char* name) const = 0;
	virtual const char* MBS_EDC_TreeGetString(const char* name) const = 0;
	//set values in MBS_EDC:
	virtual void MBS_EDC_TreeSetDouble(double data, const char* name) = 0;
	virtual void MBS_EDC_TreeSetInt(int data, const char* name) = 0;
	virtual void MBS_EDC_TreeSetString(const char* data, const char* name) = 0;

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//exchange data: elements, MBS-System
	virtual int GetNElements() const = 0;
	virtual void GetElementData(int i, int type, int value, ElementDataContainer& edc) = 0;
	virtual int SetElementData(int i, int type, int value, ElementDataContainer& edc) = 0;

	//$AD 2013-07-08: Manipulate Content of arrays in IOElements from 2D Draw Window
	virtual void InsertIOElemConNode(int elemnr, int list_idx, int input_nr, Vector2D pos) = 0;
	virtual void DeleteIOElemConNode(int elemnr, int list_idx) = 0;
	//$AD 2013-07-10: Change Position of a single element (MBS element, conNode, ...) 
	virtual void MoveConNode2D(int elemnr, int list_idx, double delta_x, double delta_y) = 0;
	virtual void MoveElement(int elemnr, double delta_x, double delta_y, double delta_z = 0.) = 0;
	//$ YV 2012-11-29: commented out, these functions go into the general MBS interface
	//virtual const class CMBSParser& MBSParser() const = 0;
  //virtual class CMBSParser& MBSParser() = 0;

	virtual int CallCompFunction(int action, int option = 0, int value = 0, ElementDataContainer* edc = NULL) = 0;
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // access functions for Sensor / PlotTool (Plottool has no classes MBS, Sensor, etc included)
	virtual TArray<double>* GetSensorValuesArrayPtr(int i) = 0;       // get the i-th sensors (as registered in mbs) own stored values from internal array
	virtual TArray<double>* GetSensorTimesArrayPtr(int i) = 0;        // get the i-th sensors (as registered in mbs) own stored times from internal array
	virtual mystr GetSensorName(int i) {return 0;};                   // get the i-th sensors (as registered in mbs) name
	virtual int GetSensorNumber(mystr& name) {return 0;};             // get the number of the (first) sensor that matches the sensorname

	// user interface part of the software needs access to the object factory
	virtual MBSObjectFactoryInterface * GetObjectFactory() = 0;
	
	// handling of variable types for post-processing (contour plotting)
	virtual const TArray<FieldVariableDescriptor> & GetAvailableFieldVariables() = 0;
	virtual int GetIndexOfActualPostProcessingFieldVariable() = 0;
	virtual void SetIndexOfActualPostProcessingFieldVariable(int index) = 0;

	// these two structures are responsible for the streamed data storage/loading
	struct DataSaver
	{
		virtual void SetTime(double time) = 0;
		virtual DataSaver & operator << (int) = 0;
		virtual DataSaver & operator << (double) = 0;
		// this function accepts a pointer to a C-style zero terminated string
		virtual DataSaver & operator << (const char *) = 0;
	};

	struct DataLoader
	{
		virtual double GetTime() = 0;		// should be called before other
		// data reading operations
		virtual DataLoader & operator >> (int &) = 0;
		virtual DataLoader & operator >> (double &) = 0;
		// here the working module obtains a new pointer  to
		// the the string is C-style zero terminated string char * p,
		// so that delete[] p should be executed afterwards
		virtual DataLoader & operator >> (char * & p) = 0;
	};

	// these functions are due to the data storage
	// after each call of ComputationFeedBack::ResultsUpdated()
	// the function StoreResultsIsOn() is called to check if the results should be stored
	virtual int StoreResultsIsOn() = 0;
	// here the working object should remove all saved states
	virtual int RemoveResults() = 0;
	// here the working object should save its actual state
	virtual void StoreResults(DataSaver & storage, double& m_TimePoint) = 0;
	// here the working object should replace its actual state with the one from the DataLoader
	// the function is called from the driving module with the user request
	virtual void LoadResults(DataLoader & loader, int m_TimePointNumber) = 0;
	// this function returns the actual drawing time
	virtual double GetActualDrawTime() const = 0;
	// this function returns the number of stored solution data files
	virtual int ReadSolDataInfo() = 0;
	// here the current version of Hotint is returned
	// this is particularly needed for validation of config files or data sets
	virtual const HotintVersionInfo& GetHotintVersion() const {return hotint_version;};

	virtual void SetModelData_Initialized(int flag) {};

};

#endif // __WCD_INTERFACE_H_INCLUDED__