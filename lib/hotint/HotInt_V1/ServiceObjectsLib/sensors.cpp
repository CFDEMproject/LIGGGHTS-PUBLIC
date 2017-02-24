//#**************************************************************
//#
//# filename:             sensors.cpp
//#
//# author:               Gerstmayr Johannes
//#												Vetyukov Yury
//#
//# generated:						February-June 2012
//# description:          
//#                       
//# remarks:						  HotInt sensors - 2nd generation
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


#include "mbs_interface.h"
#include "element.h"
#include "sensors.h"
#include "rendercontext.h"

// adds a copy of the sensor processor to the processing chain;
void Sensor::AddSensorProcessor(SensorProcessor & sp)
{
	sensorProcessors.Add(sp.GetCopy());
	// post computation sensor processors require the time signal to be saved in the memory
	if(sp.NeedsSignalTimeHistoryForPostComputationProcessing())
		signalStorageMode = (SensorSignalStorageMode)(signalStorageMode | SSM_InternalArray);
}

void Sensor::RemoveSensorProcessors()
{
	for(int i = 1; i <= sensorProcessors.Length(); i++)
		delete sensorProcessors(i);
	sensorProcessors.SetLen(0);
}

// these functions need to be normally overridden in the derived class (the base implementation must be called)
void Sensor::CopyFrom(const Sensor& s)
{
	// here we copy just the settings - the sensor does not have any data yet
	signalStorageMode = s.signalStorageMode;
	for(int i = 1; i <= s.sensorProcessors.Length(); i++)
		AddSensorProcessor(*s.sensorProcessors(i));
	visibleFlag = s.visibleFlag;
	drawColor = s.drawColor;
	drawDimension = s.drawDimension;
	mbs = s.mbs;
	name = s.name;
	precision = s.precision;
	openSensorWatchPlotToolAtStartUpFlag = s.openSensorWatchPlotToolAtStartUpFlag;
	ownOutputFile = s.ownOutputFile;
}

// default constructor
Sensor::Sensor(MBS * mbs)
{
	ownOutputFile = NULL;
	this->mbs = mbs;
	name = "sensor";
	precision = 17;
	signalStorageMode = SensorSignalStorageMode(SSM_CommonFile | SSM_InternalArray);		// YV->AD: why SSM_InternalArray per default?
	drawDimension = Vector3D(0.001,6,0);
	visibleFlag = true;
	openSensorWatchPlotToolAtStartUpFlag = false;
}

Sensor::Sensor(const Sensor & s)
{
	CopyFrom(s);
	ownOutputFile = NULL;
}

// performs an evaluation of the present value
void Sensor::Evaluate(double time)
{
	lastSensorSignalTime = time;
	lastSensorSignalValue = GetCurrentValueWithSensorProcessing(time);
	// the newly computed time point may be added to the internal storage
	if(GetSignalStorageMode() & SSM_InternalArray)
		AddLastSignalToHistory();
}

// the following function evaluates the present value with sensor processing but does not store the results anywhere;
// should be used during iterations within a time step
double Sensor::GetCurrentValueWithSensorProcessing(double time)
{
	double value = GetCurrentValue(time);
	// signal processors are invoked
	for(int i = 1; i <= sensorProcessors.Length(); i++)
		value = sensorProcessors(i)->ProcessCurrentValue(time, value);
	return value;
}

// this function is called by mbs after the computation is finished
// the post-processing procedures of the attached signal processors are invoked one after another;
// these may write the results into the provided file stream (managed by mbs) and provide some sensor evaluation vectors,
// which are then combined into a common vector signalProcessingEvaluationData
void Sensor::ApplyPostComputationSensorProcessing(ofstream * outputFile)
{
	sensorProcessingEvaluationData.SetLen(0);
	for(int i = 1; i <= sensorProcessors.Length(); i++)
	{
		TArray<double> newEvaluationData = sensorProcessors(i)->DoPostComputationProcessing(signalHistoryTimes, signalHistoryValues, outputFile);
		// and the newly obtained data is merged into the global data vector
		int oldLength = sensorProcessingEvaluationData.Length();
		sensorProcessingEvaluationData.SetLen(oldLength + newEvaluationData.Length());
		for(int j = 1; j <= newEvaluationData.Length(); j++)
			sensorProcessingEvaluationData(oldLength + j) = newEvaluationData(j);
	}
}

// derived sensor classes may redefine the default drawing behavior
// label (text caption) for the sensor will be provided by mbs
void Sensor::Draw(const char * textLabel)
{
	if (!GetVisible())
		return;

	double rad = GetDrawDimension().X();
	Vector3D col(0.7,0.5,0.1);

	for (int i = 1; i <= GetNumberOfDrawingPositions(); i++)
	{
		Vector3D v0 = GetDrawPosition(1);

		double s = mbs->GetDOption(118);
		mbs->SetColor(Vector3D(0.9,0.5,0.45)); //steel blue

		Vector3D v1(-s, 0,0);
		Vector3D v2( s, 0,0);
		Vector3D v3( 0,-s,0);
		Vector3D v4( 0, s,0);
		Vector3D v5( 0,0,-s);
		Vector3D v6( 0,0, s);

		double d = mbs->GetDOption(114); //global line thickness
		mbs->MyDrawLine(v1+v0,v2+v0,d);
		mbs->MyDrawLine(v3+v0,v4+v0,d);
		mbs->MyDrawLine(v5+v0,v6+v0,d);

		mbs->GetRC()->PrintText3D((float)v0.X(), (float)v0.Y(), (float)v0.Z(), textLabel);

		//draw sphere after sensor frame, because of transparency
		mbs->SetColor(col);
		mbs->DrawSphere(v0, rad, (int)GetDrawDimension().Y());
	}
}

//$ DR 2012-10
void Sensor::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	GetElementDataAuto(edc);
}

//$ DR 2012-10
int Sensor::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = 1;
	mystr type_old = GetTypeName();

	rv = SetElementDataAuto(edc);

	mystr type_new = edc.TreeGetString("sensor_type");
	if(!type_new.Compare(type_old))
	{
		mbs->UO().InstantMessageText("ERROR: You MUST NOT change the type of the sensor!");
		return 0;
	}

	return rv;
}