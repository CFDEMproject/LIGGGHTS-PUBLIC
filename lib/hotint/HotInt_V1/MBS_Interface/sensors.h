//#**************************************************************
//#
//# filename:             sensors.h
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

#pragma once

#include "sensorProcessors.h"

// bit-wise flags for possible data storage options in a sensor - can be combined
enum SensorSignalStorageMode
	{
		SSM_None = 0,					// no data storage
		SSM_CommonFile = 1,		// general solution file of the simulation - prescribed output frequency is governed by MBS
		SSM_OwnFile = 2,			// separate file with the results of this particular sensor - prescribed output frequency is governed by MBS
		SSM_InternalArray = 4	// internal array in memory - all time signal points; will be set automatically if a signal processor with post computation is applied
	};

// - this is a base class for actually used sensor classes
// - these sensors can have only one scalar value at a time (for multiple-value sensors other mechanisms will be involved)
// - the base class includes data storing functions
// - the base class is responsible for the interaction with possible signal processors
// - the base class provides interaction with mbs
// - there are no copy features for the sensors: the instances are created with "new" and added to mbs per pointer
// - tha base class is abstract and cannot be instantiated
class Sensor //$EDC$[beginclass,classname=Sensor]
{
private:
	// data of the sensor
	SensorSignalStorageMode signalStorageMode;
	// these are the processors for the signals attached to the present sensor;
	// the processors are applied one after another for the actual value, and are invoked at the post computation stage
	TArray<SensorProcessor*> sensorProcessors;

	// computation results
	// this value and time point correspond to the last call to Evaluate() - instant value of the sensor between time steps;
	// will be used by MBS::WriteSol()
	double lastSensorSignalTime;
	double lastSensorSignalValue;
	// the arrays below contain measured signal time history if signalStorageMode & SSM_InternalArray != 0
	// the two arrays are made protected to be available in MBSSensor - as soon as it is removed from the system, they can be made private again
protected:
	TArray<double> signalHistoryValues;		// stored measured values
	TArray<double> signalHistoryTimes;			// times of the measurements
	int sensnum; //$EDC$[varaccess,EDCvarname="sensor_number",EDCfolder="",readonly,tooltiptext="number of the sensor in the mbs"] //number in MBS-system

private:
	// at the post-computation stage the sensor processors may generate a set of values,
	// which characterize the measured time signal and which might then be used for optimization, sensitivity analysis, etc.
	TArray<double> sensorProcessingEvaluationData;
	
	// draw settings
	bool visibleFlag;
	Vector3D drawColor;
	Vector3D drawDimension;

	// system variables
	MBS * mbs;
	mystr name;		//$EDC$[varaccess,EDCvarname="name",EDCfolder="",tooltiptext="name of the sensor for the output files and for the plot tool"]
	bool openSensorWatchPlotToolAtStartUpFlag;		// if this flag is set, then a window with the time history of the sensor signal is opened by mbs at the initialization stage after assembling the model
	int precision;	// precision of output of the sensor signal
	ofstream * ownOutputFile;		// stream to which the sensor signal history is written if SSM_OwnFile flag is set - is created, opened, written to, closed and deleted by mbs

	// adds a measured value for the present time instant
	void AddLastSignalToHistory()
	{
		signalHistoryValues.Add(lastSensorSignalValue);
		signalHistoryTimes.Add(lastSensorSignalTime);
	}

public:
	// access to the above variables
	MBS * GetMBS() { return mbs; }
	const MBS * GetMBS() const { return mbs; }
	SensorSignalStorageMode GetSignalStorageMode() const { return signalStorageMode; }
	void SetSignalStorageMode(SensorSignalStorageMode mode) { signalStorageMode = mode; }
	void ModifySignalStorageMode(SensorSignalStorageMode flagsSet, SensorSignalStorageMode flagsRemove = SSM_None)
	{
		signalStorageMode = SensorSignalStorageMode( (signalStorageMode | flagsSet) & (~flagsRemove) );
	}
	TArray<double>* GetSignalHistoryValuesPtr() { return &signalHistoryValues; }
	TArray<double>* GetSignalHistoryTimesPtr() { return &signalHistoryTimes; }

	virtual int GetOwnNum() const {return sensnum;}
	virtual void SetOwnNum(int i) {sensnum = i;}

	// draw settings
	bool GetVisible() const { return visibleFlag; }
	void SetVisible(bool visibleFlag) { this->visibleFlag = visibleFlag; }
	Vector3D GetDrawColor() const { return drawColor; }
	void SetDrawColor(Vector3D drawColor) { this->drawColor = drawColor; }
	Vector3D GetDrawDimension() const { return drawDimension; }
	void SetDrawDimension(Vector3D drawDimension) { this->drawDimension = drawDimension; }
	// name of the sensor used in the user interface (sensors list, plot tool, etc.)
	// and in the text files with the values; should be set in the constructor of a derived class
	mystr GetSensorName() const { return name; }
	//$ PG 2013-1-16:[
	// name of the sensor MAY be automatically set in SensorName::AfterSetFunction,
	// which is overwritten in the derived classes, and called at the end of the specific set functions of the derived classes
	// (such as FieldVariableElementSensor::SetFVESPos3D, SingleElementDataSensor::SetSingleElementDataSensor, LoadSensor::SetLoadSensor),
	// i.e., SetSensorName has to be called AFTER those set-functions in the models-cpp-files to take effect
	void SetSensorName(mystr name) { this->name = name; }
	//$ PG 2013-1-16:]
	virtual bool IsSensorOnEigensystem() { return 0; } //$ PG 2013-11-27: potentially overridden in derived class
	bool GetOpenSensorWatchPlotToolAtStartUpFlag() const { return openSensorWatchPlotToolAtStartUpFlag; }
	void SetOpenSensorWatchPlotToolAtStartUpFlag(bool openSensorWatchPlotToolAtStartUpFlag) { this->openSensorWatchPlotToolAtStartUpFlag = openSensorWatchPlotToolAtStartUpFlag; }
	// precision of output to the text file
	int GetPrecision() const { return precision; }
	void SetPrecision(int precision) { this->precision = precision; }
	// output file
	ofstream * GetOwnOutputFile() const { return ownOutputFile; }
	void SetOwnOutputFile(ofstream * ownOutputFilei) { ownOutputFile = ownOutputFilei; }

	// adds a copy of a sensor processor to the processing chain;
	void AddSensorProcessor(SensorProcessor & sp);
	// if some of the signal processors require a file stream to write the results of the postcomputation evaluation,
	// the sensor needs to report it;
	// the function is virtual for compatibility with the obsolete MBSSensor
	virtual bool NeedsFileForPostComputationProcessing()
	{
		for(int i = 1; i <= sensorProcessors.Length(); i++)
			if(sensorProcessors(i)->NeedsFileForPostComputationProcessing())
				return true;
		return false;
	}
	// clears the list of attached sensor processors
	void RemoveSensorProcessors();

	// access to computation results; some functions below are virtual for compatibility with the obsolete MBSSensor
	virtual bool HasSensorProcessingEvaluationData() { return sensorProcessingEvaluationData.Length() != 0; }
	virtual TArray<double> & GetSignalProcessingEvaluationData() { return sensorProcessingEvaluationData; }
	double GetLastValue() { return lastSensorSignalValue; }

	// building dependencies requires that the elements, which affect the signal of this sensor, are known
	// these numbers may be modified by mbs when elements are inserted/deleted
	virtual int GetNumberOfRelatedElements() { return 0; }
	virtual int& GetRelatedElementNumber(int nElement) { assert(0); static int dummy = 0; return dummy; }
	virtual int GetNumberOfRelatedNodes() { return 0; }			//$ DR 2013-05-21
	virtual int& GetRelatedNodeNumber(int nNode) { assert(0); static int dummy = 0; return dummy; } //$ DR 2013-05-21
	virtual int GetNumberOfRelatedLoads() { return 0; }			//$ DR 2013-05-21
	virtual int& GetRelatedLoadNumber(int nLoad) { assert(0); static int dummy = 0; return dummy; } //$ DR 2013-05-21
	virtual int GetNumberOfRelatedSensors() { return 0; }			//$ DR 2013-05-21
	virtual int& GetRelatedSensorNumber(int nSensor) { assert(0); static int dummy = 0; return dummy; } //$ DR 2013-05-21


protected:
	// the functions below are assumed to be overridden
	// the following function computes the present value of the sensor
	// and it must be implemented in the derived sensor class;
	// this function is unavailable outside of the class as this value is not yet processed by sensor processors
	virtual double GetCurrentValue(double time) = 0;
public:
	// literal identifier of the type of the sensor
	virtual mystr GetTypeName() {return "Sensor";};	//$EDC$[funcaccess,EDCvarname="sensor_type",tooltiptext="specification of sensor type. Once the sensor is added to the mbs, you MUST NOT change this type anymore!"]
	// the derived sensor classes may provide particular positions for drawing (if it makes sense);
	// then this function tells how many points should be plotted for this sensor in the 3D scene
	virtual int GetNumberOfDrawingPositions() { return 0; }
	// the drawing positions are indexed as there might be several ones for a given sensor
	virtual Vector3D GetDrawPosition(int i) { assert(0); return Vector3D(); }
	// for the auto generated documentation //$ DR 2012-10 added
	virtual mystr GetSensorTypeTexDescription() {return mystr("");};

protected:
	// default constructor is just for copy making
	Sensor() {}
	// this function should be called after setting the data of a sensor by
	// - set functions
	// - setting the data from an EDC
	// - modifying or setting up a sensor via user interface
	// in this function additional initialization is performed, e.g. a name is generated
	virtual void AfterSetFunction() {}

public:
	// creating / copying - must be overridden in the derived class
	virtual Sensor * GetCopy() = 0;
	// these functions need to be normally overridden in the derived class (the base implementation must be called)
	virtual void CopyFrom(const Sensor& s);

	// constructor / destructor
	// default constructor with member initialization
	Sensor(MBS * mbs);
	Sensor(const Sensor & s);
	// default destructor - deletes signal processors
	virtual ~Sensor() { RemoveSensorProcessors(); }
	// performs an evaluation of the present value, then sensor processing,
	// and saves to time history if the corresponding flag is set;
	// should be used between time steps
	// the function is virtual for compatibility with the obsolete MBSSensor
	virtual void Evaluate(double time);
	// the following function evaluates the present value with sensor processing but does not store the results anywhere;
	// should be used during iterations within a time step
	// the function is virtual for compatibility with the obsolete MBSSensor
	double GetCurrentValueWithSensorProcessing(double time);
	// derived sensor classes may redefine the default drawing behavior
	// label (text caption) for the sensor will be provided by mbs
	virtual void Draw(const char * textLabel);
	// testing the data integrity
	virtual bool IsConsistent(mystr & errorStr) { return true; }
	// this function is called by mbs after the computation is finished
	// the post-processing procedures the attached signal processors are invoked one after another;
	// these may write the results into the provided file stream (managed by mbs) and provide some sensor evaluation vectors,
	// which are then combined into a common vector signalProcessingEvaluationData
	// the function is virtual for compatibility with the obsolete MBSSensor
	virtual void ApplyPostComputationSensorProcessing(ofstream * outputFile);

	// sensor data setting/retrieving; these functions need to be normally inherited (the base class implementation is to be called in the derived class)
	virtual void GetElementData(ElementDataContainer& edc);
	virtual int SetElementData(ElementDataContainer& edc);
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

};//$EDC$[endclass,Sensor]