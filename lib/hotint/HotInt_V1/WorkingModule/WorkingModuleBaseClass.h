//#**************************************************************
//# filename:							WorkingModuleBaseClass.h
//#
//# author:               Yury Vetyukov
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
 
#ifndef __SAMPLE_WORKING_MODULE_H_INCLUDED__
#define __SAMPLE_WORKING_MODULE_H_INCLUDED__

#include "WinCompDriverInterface.h"
#include "UserOutput.h"

// Default implementation of the working module
// the functionality should be added in the derived class

class WorkingModuleBaseClass: public WCDInterface
{

	// these are the functions from the interface
	void SetUserInterface(UserInterface * pui);
	virtual void InitializeAfterConfigLoaded();
	virtual UserInterface* GetUserInterface() const
	{
		return uo.pUI;
	}


	virtual void ButtonDrawMode(char*& text) {};
public:
	virtual int IsComputationInProgress() { return bComputationIsInProgress; }
private:
	virtual void StopWhenPossible() { bStopWhenPossible = 1; }
	virtual int IsPaused() { return bPause; }
	virtual void Pause() { bPause = 1; }
	virtual void Resume() { bPause = 0; }
	virtual void Go(ComputationFeedBack * p_cfb);
	virtual void SetPCFB(ComputationFeedBack * p_cfb) {pCFB = p_cfb;}

	// TODO:
	// these functions should be redefined in the derived class
	virtual void RenderScene(RenderContext * pRC) {}
	virtual void RenderControlWindow(ControlWindowContext* pCWC) {}    //!AD: 2012-12-13
	virtual float GetSceneMaxAbsCoordinate() { return 1;}
	virtual int AllowRotation() { return 1; }
	virtual void PrintData(const char * pDataSpecification) {}
	virtual int StoreResultsIsOn() { return 0; }
	virtual void StoreResults(DataSaver & storage, double& m_TimePoint) {}
	virtual void LoadResults(DataLoader & loader, int m_TimePointNumber) {}
	//virtual double GetDrawTime() const { return 0.; }
	virtual void ResetComputation() {};
	virtual void PrintTimingList() {};

protected:
	mutable UserOutput uo;
	const char  * ENDL;

	int bComputationIsInProgress;
	int bComputeEigenmodes;
	int bStopWhenPossible;
	int bPause;    // 0: don't pause, 1: pause calculation when possible, 2: calculation paused

	virtual void InitFirst() {};
	ComputationFeedBack * pCFB;

	// TODO:
	// this function should be redefined with the derived class
	// with the implementation of the actual computation code
	virtual void PerformComputation()
	{
		while(!bStopWhenPossible && !bPause)
			;
	}

	virtual void PerformComputation2()
	{
	//	while(!bStopWhenPossible && !bPause)
	//		;
	}

public:
	WorkingModuleBaseClass(): ENDL("\n"), bComputationIsInProgress(0), bComputeEigenmodes(0) {}
	ComputationFeedBack* Get_pCFB() {return pCFB;}

};


#endif	//__SAMPLE_WORKING_MODULE_H_INCLUDED__