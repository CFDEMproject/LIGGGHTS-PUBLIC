//#**************************************************************
//# filename:							WorkingModuleBaseClass.cpp
//#
//# author:               Johannes Gerstmayr
//#
//# generated:			  2007
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
 
#include "StdAfx.h"
#include "mbs.h"

// Default implementation of the working module
// the functionality should be added in the derived class

void WorkingModuleBaseClass::SetUserInterface(UserInterface * pui)
{
	uo.pUI = pui;
	uo << "Working module initialized\n";
	InitFirst();
}
void WorkingModuleBaseClass::InitializeAfterConfigLoaded()
{
	bPause = 0;
}

void WorkingModuleBaseClass::Go(ComputationFeedBack * p_cfb) 
{
	pCFB = p_cfb;

	// initializing the computation
	bComputationIsInProgress = 1;
	bStopWhenPossible = 0;
	bPause = 0;
	uo.pUI->StatusText("Running computation");

	PerformComputation();

	uo.pUI->StatusText("Finished computation");
	// finalizing the computation
	bComputationIsInProgress = 0;
	pCFB->FinishedComputation();
}

