//#***************************************************************************************
//# filename:     stepsettings.h
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

#pragma once

typedef enum 
{ 
	TCSRStep = 0,          // step function, loadfactor is set immediately
	TCSRLinear = 1,        // linear interpolation, loadfactor is reached at end of interval
	TCSRExponential = 2    // exponential interpolation
} TCSRamp;

// struct for StepSettings in MBSLoad, Element, Element-derived,
class StepSettings // valid for one interval - interval times are defined in TimeInt::CSEndTimes
{
public:	//lifecycle
	StepSettings(): loadfactor(0.),rampmode(TCSRStep) {}
	StepSettings(StepSettings& other) { CopyFrom(other); }
	~StepSettings() {}
public: //lifecycle II
	StepSettings(double loadfactori, TCSRamp rampmodei) {	SetStepSettings(loadfactori, rampmodei); }
	StepSettings(double loadfactori) { SetStepSettings(loadfactori, TCSRLinear); }
	void CopyFrom(StepSettings& other) 
	{ 
		loadfactor = other.loadfactor;
		rampmode = other.rampmode;
	}
	StepSettings* GetCopy() { return new StepSettings(*this); } 
	void SetStepSettings(double loadfactori, TCSRamp rampmodei) { loadfactor = loadfactori; rampmode = rampmodei; }

public:
	void SetLoadFactor(double lf) { loadfactor = lf; }  
	double const GetLoadFactor() const { return loadfactor; }
	double& LoadFactor() { return loadfactor; }
	
	void SetRampMode(TCSRamp rm) { rampmode = rm; }  
	TCSRamp const GetRampMode() const { return rampmode; }
	TCSRamp& RampMode() { return rampmode; }

private:
// for calculation of LoadStepFactor
	double loadfactor;   // value at end of interval    
	TCSRamp rampmode;        // interpolation mode (see above) - weighting factor for 
};