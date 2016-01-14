//#***************************************************************************************
//# filename:     serviceFunctions.cpp
//#
//# author:				Johannes Gerstmayr, Yuri Vetyukov
//# 
//# generated:      
//# description:  
//#                       
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


// here some global functions are defined, which replace the same functions in the kernel

#include "mbs_interface.h"
#include "femathhelperfunctions.h"

// functions for starting and stopping timers need to be implemented in the client dll differently:
// here we cannot access the timing arrays directly

TArray<double> * ptrTMtspent;
TArray<double> * ptrTMtstart;

#ifdef timeron
void TMStartTimer(const int& i)
{
	(*ptrTMtstart)(i) = GetClockTime();
}

void TMStopTimer(const int& i)
{
	(*ptrTMtspent)(i) += GetClockTime() - (*ptrTMtstart)(i);
}
#endif

extern UserOutputInterface * global_uo;

void TIMBSWarningHandle(const char* warn, int use_instant_message_text)
{
//$ AD 2011-11: flag may have negative values, then no output
	if(use_instant_message_text < 0) return;

	global_uo->SetLocalMessageLevel(UO_LVL_warn);

	if(!use_instant_message_text)
		(*global_uo) << "WARNING:" << warn << "\n";
	else
		global_uo->InstantMessageText(mystr("WARNING:")+mystr(warn));
}
