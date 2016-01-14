//#**************************************************************
//#
//# filename:             timeintlog.h
//#
//# author:               Peter Gruber
//#
//# generated:						26.11.2013
//# description:          log information concerning time integration & multibody system
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
//#***************************************************************************************
 


#pragma once

class TimeIntLog
{
public:
	TimeIntLog() {Init();}
	virtual void Init()
	{
		writesolcnt = 0;
		evalfcnt = 0;
		evalf_jaccnt = 0;
		testcnt = 0;
		changestep = 0;
		jaccount = 0;
		TInewtonit = 0;
		TInewtonitsum = 0;
		TInonlinit = 0;
	}

public:
	int writesolcnt;
	int evalfcnt; //number of right-hand-side evaluations
	int evalf_jaccnt; //number of right-hand-side evaluations for jacobian
	int testcnt; //for test purposes
	int changestep; //counter for number of changed step sizes
	int jaccount; //counter for jacobians
	int TInewtonit;
	int TInewtonitsum;
	int TInonlinit;
};