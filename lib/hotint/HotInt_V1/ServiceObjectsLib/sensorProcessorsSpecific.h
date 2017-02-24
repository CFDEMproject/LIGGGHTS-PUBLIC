//#**************************************************************
//#
//# filename:             sensorProcessorsSpecific.h
//#
//# author:               Gerstmayr Johannes
//#												Vetyukov Yury
//#
//# generated:						February-June 2012
//# description:          
//#                       
//# remarks:						  HotInt sensors - specific sensor processors
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

#pragma once

#include "sensors.h"

class OffsetScaleSensorProcessor : public SensorProcessor
{
	double scale;
	double offset;

public:
	OffsetScaleSensorProcessor(double scale, double offset)
	{
		this->scale = scale;
		this->offset = offset;
	}
	virtual double ProcessCurrentValue(double time, double value)
	{
		return value * scale + offset;
	}
	virtual SensorProcessor * GetCopy()
	{
		OffsetScaleSensorProcessor * sp = new OffsetScaleSensorProcessor(scale, offset);
		return sp;
	}
};

class MaxAbsValueSensorProcessor : public SensorProcessor
{
	double currentValue;

public:
	MaxAbsValueSensorProcessor()
	{
		currentValue = 0;
	}
	virtual double ProcessCurrentValue(double time, double value)
	{
		currentValue = max(currentValue, fabs(value));
		return currentValue;
	}
	virtual SensorProcessor * GetCopy()
	{
		MaxAbsValueSensorProcessor * sp = new MaxAbsValueSensorProcessor();
		return sp;
	}
};

class ReferenceComparerSensorProcessor : public SensorProcessor
{
	MathFunction referenceSignal;
	MBS * mbs;

public:
	ReferenceComparerSensorProcessor(MBS * mbs, MathFunction & referenceSignal)	// the provided object should be deleted by the client code
	{
		this->referenceSignal = referenceSignal;
		this->mbs = mbs;
	}
	virtual double ProcessCurrentValue(double time, double value)
	{
		return value - referenceSignal.Evaluate(time);
	}
	virtual SensorProcessor * GetCopy()
	{
		ReferenceComparerSensorProcessor * sp = new ReferenceComparerSensorProcessor(mbs, referenceSignal);
		return sp;
	}
};

//$ DR 2012-08-24 ReferenceSensorComparerSensorProcessor added
// computes the difference of the actual value of the sensor and the actual value of the reference sensor
// output = sensor - reference
class ReferenceSensorComparerSensorProcessor : public SensorProcessor
{
	Sensor *  sensP;

public:
	ReferenceSensorComparerSensorProcessor(Sensor * referenceSensor)	
	{
		this->sensP = referenceSensor;
	}
	virtual double ProcessCurrentValue(double time, double value)
	{
		return value - this->sensP->GetCurrentValueWithSensorProcessing(time);
	}
	virtual SensorProcessor * GetCopy()
	{
		ReferenceSensorComparerSensorProcessor * sp = new ReferenceSensorComparerSensorProcessor(sensP);
		return sp;
	}
};


class FFTSensorProcessor : public SensorProcessor
{
public:
	virtual bool NeedsSignalTimeHistoryForPostComputationProcessing() { return true; }
	virtual bool NeedsFileForPostComputationProcessing() { return true; }
	virtual TArray<double> DoPostComputationProcessing(TArray<double> & sensorHistoryTimes, TArray<double> & sensorHistoryValues, ofstream & outputFile)
	{
		// TODO
		assert(0);
		return TArray<double>();
	}
	virtual FFTSensorProcessor * GetCopy()
	{
		FFTSensorProcessor * sp = new FFTSensorProcessor();
		return sp;
	}
};