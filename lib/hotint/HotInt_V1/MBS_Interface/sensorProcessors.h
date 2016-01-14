//#**************************************************************
//#
//# filename:             sensorProcessors.h
//#
//# author:               Gerstmayr Johannes
//#												Vetyukov Yury
//#
//# generated:						February-June 2012
//# description:          
//#                       
//# remarks:						  HotInt sensors - sensor processor class
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

// - base class for a processor of a time signal measured by a sensor
// - both instant processing of the present sensor values (min-max, scale-shift, etc.)
//	 and final sensor computation (fft, reference data, etc.) are possible
class SensorProcessor
{
public:
	// a virtual destructor will be needed when a derived class needs some clean up in its destructor
	virtual ~SensorProcessor() {}
	// processing of a single physically measured value - one after another, returns the processed value
	// may also store data in the local variables of the sensor processor and give it out at the postcomputation stage
	virtual double ProcessCurrentValue(double time, double value) { return value; }
	// tells whether this processor performs operation on the whole time history of the signal;
	// if yes, then the sensor must store data in memory
	virtual bool NeedsSignalTimeHistoryForPostComputationProcessing() { return false; }
	// tells whether this processor needs a file opened for writing
	// to store the results of the postcomputation processing
	// if yes, then the sensor must store data in memory
	virtual bool NeedsFileForPostComputationProcessing() { return false; }
	// - performes postcomputation signal processing
	// - may process the whole time signal of the sensor or give out data, accumulated during the computation
	//   (the first two arguments are empty if none of the sensor processors has requested time history)
	// - output is the vector of values,
	//   which are then assembled into "SignalProcessingEvaluationData" and may be used for optimization or sensitivity analysis
	// - may write the results to the provided file stream
	//   (the last argument is NULL if none of the sensor processors has requested a file stream)
	virtual TArray<double> DoPostComputationProcessing(TArray<double> & sensorHistoryTimes, TArray<double> & sensorHistoryValues, ofstream * outputFile) { return TArray<double>(); }
	// a sensor processor must be able to clone itself
	virtual SensorProcessor * GetCopy() = 0;
};