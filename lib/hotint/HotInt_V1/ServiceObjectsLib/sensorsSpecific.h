//#**************************************************************
//#
//# filename:             sensorsSpecific.h
//#
//# author:               Gerstmayr Johannes
//#												Vetyukov Yury
//#
//# generated:						February-June 2012
//# description:          
//#                       
//# remarks:						  HotInt sensors - particular sensor types
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

// added by YV
// this is an abstract base class for sensors, which are attached to a particular element in MBS
class ElementSensor : public Sensor //$EDC$[beginclass,classname=ElementSensor,parentclassname=Sensor]
{
	int elementNumber; //$EDC$[varaccess,EDCvarname="element_number",EDCfolder="",tooltiptext="number of the element, to which the sensor is applied"]

protected:
	int& GetElementNumber() { return elementNumber; }
	Element & GetElement() { return GetMBS()->GetElement(GetElementNumber()); }		//$ DR 2012-09-12 changed GetElement(elementNumber) to GetElement(GetElementNumber())

	// default constructor is just for copy making
	ElementSensor() {}
	
public:
	ElementSensor(MBS * mbs) : Sensor(mbs)
	{
	}
	void SetElementNumber(int elementNumber)
	{
		this->elementNumber = elementNumber;
	}
	virtual void CopyFrom(const Sensor& s)
	{
		Sensor::CopyFrom(s);
		ElementSensor & S = (ElementSensor&)s;
		elementNumber = S.elementNumber;
	}
	virtual bool IsConsistent(mystr & errorStr)
	{
		if(elementNumber <= 0 || elementNumber > GetMBS()->NE())
		{
			errorStr = mystr("Sensor ") +mystr(GetOwnNum()) + mystr(" has an invalid element number: ") + mystr(elementNumber) + "\n";
			return false;
		}
		return true;
	}
	virtual int GetNumberOfRelatedElements() { return 1; }
	virtual int& GetRelatedElementNumber(int nElement) { return elementNumber; }
	virtual void GetElementData(ElementDataContainer & edc)
	{
		Sensor::GetElementData(edc);
	}
	virtual int SetElementData(ElementDataContainer& edc) //$ DR 2012-10 return value changed from bool to int
	{
		int rv = Sensor::SetElementData(edc);
		AfterSetFunction();	// DR 2012-06-25 added
		return rv;
	}
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

};//$EDC$[endclass,ElementSensor]

// added by YV
// these sensors measure a particular field variable for a given element
// either a local node number or a local position may be specified
// the local position may be either 2D or 3D
class FieldVariableElementSensor : public ElementSensor //$EDC$[beginclass,classname=FieldVariableElementSensor,parentclassname=ElementSensor]
{
protected:
	FieldVariableDescriptor fieldVariableDescriptor;
	int localNodeNumber;			//$EDC$[varaccess,EDCvarname="node_number",EDCfolder="",tooltiptext="local node number. If > 0, then the position of this node is used."]
	Vector3D localPosition;		//$EDC$[varaccess,EDCvarname="local_position",EDCfolder="",tooltiptext="local position at which the field variable is evaluated."]
	//bool localPosition2D;			//[varaccess,EDCvarname="local_position2D",EDCfolder="",tooltiptext="Only the two first components in local_position are used!"]
	int dimension;	//$ DR 2013-06-19 replaces bool localPosition2D

	virtual double GetCurrentValue(double time);

	virtual int GetNumberOfDrawingPositions() { return 1; }
	virtual Vector3D GetDrawPosition(int i)
	{
		return GetElement().GetPosD(localPosition);
	}

	virtual void AfterSetFunction()
	{
		//SetSensorName(GetTypeName());	//$ DR 2012-06-25 commented out
		if (localNodeNumber == 0 && (GetElementNumber() <= GetMBS()->NE()))
		{ 
			dimension = GetElement().Dim();
		}
	}

	// default constructor is just for copy making
	FieldVariableElementSensor() {}

public:
	FieldVariableElementSensor(MBS * mbs) : ElementSensor(mbs)
	{
	}
	// three set-functions for three variants of specification of the local position
	void SetFVESPos3D(int elementNumber, FieldVariableDescriptor fieldVariableDescriptor_, Vector3D localPosition);
	void SetFVESPos2D(int elementNumber, FieldVariableDescriptor fieldVariableDescriptor_, Vector2D localPosition);
	void SetFVESNode(int elementNumber, FieldVariableDescriptor fieldVariableDescriptor_, int localNodeNumber);
	virtual void CopyFrom(const Sensor& s);
	virtual Sensor * GetCopy()
	{
		FieldVariableElementSensor * S = new FieldVariableElementSensor();
		S->CopyFrom(*this);
		return S;
	}
	virtual mystr GetTypeName();
	virtual mystr GetSensorTypeTexDescription(); 	// for the auto generated documentation //$ DR 2012-10 added
	virtual bool IsConsistent(mystr & errorStr);
	virtual void GetElementData(ElementDataContainer& edc);
	virtual int SetElementData(ElementDataContainer& edc); //$ DR 2012-10 return value changed from bool to int
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	Vector3D GetLocalPosition()	{ return localPosition; }
};//$EDC$[endclass,FieldVariableElementSensor]

// added by YV
// this sensor reads values from a given element using Element::ReadSingleElementData access function.
class SingleElementDataSensor : public ElementSensor 
{
protected:
	ReadWriteElementDataVariableType RWdata;		// corresponding identifier in the internal presentation
	int successfully_converted;
	//virtual void AfterSetFunction();
	virtual double GetCurrentValue(double time);

	// default constructor is just for copy making
	SingleElementDataSensor() {}

public:
	SingleElementDataSensor(MBS * mbs) : ElementSensor(mbs)
	{
	}
	void SetSingleElementDataSensor(int elementNumber, const char * variable_name);
	virtual bool IsConsistent(mystr & errorStr);
	virtual void GetElementData(ElementDataContainer& edc);
	virtual int SetElementData(ElementDataContainer& edc);
	virtual mystr GetTypeName(){return "ElementSensor"; }
	virtual mystr GetSensorTypeTexDescription(); 	// for the auto generated documentation //$ DR 2012-11 added

	virtual void CopyFrom(const Sensor& s);
	virtual Sensor * GetCopy()
	{
		SingleElementDataSensor * S = new SingleElementDataSensor();
		S->CopyFrom(*this);
		return S;
	}
};

// added by YV; to be implemented
// prospective description:
// abstract base class for sensors, which are attached to a particular node in MBS
class NodeSensor// : public Sensor
{
};

// added by YV; to be implemented
// prospective description:
// computes an average for all elements, which meet at this node;
// the field variable can be computed either directly in the node or with some offset into the elements
class FieldVariableNodeSensor : public NodeSensor
{
};

// added by YV; to be implemented
// prospective description:
// nodal DOF, etc.
class SpecialValueNodeSensor : public NodeSensor
{
};




// added by DR
// sensor that consists of a list of sensors
// this replaces the old MBSMultipleSensor
class MultipleSensor : public Sensor //$EDC$[beginclass,classname=MultipleSensor,parentclassname=Sensor]
{
	TArray<int> sensorNrs; //$EDC$[varaccess,variable_length_vector,EDCvarname="sensor_numbers",EDCfolder="",tooltiptext="number of the sensors, that are used for computation"]
	int operation;	
	Vector weights; //$EDC$[varaccess,variable_length_vector,EDCvarname="weights",EDCfolder="",tooltiptext="weights for e.g. a weighted sum. This vector must have the same length as sensor_numbers or must be empty!"]

protected:
	int& GetSensorNumber(int i) { return sensorNrs(i); }
	Sensor& GetSensor(int i) { return GetMBS()->GetSensor(GetSensorNumber(i)); }		

	// default constructor is just for copy making
	MultipleSensor() {}
	
public:
	MultipleSensor(MBS * mbs) : Sensor(mbs)
	{
	}

	void SetMultipleSensor(TArray<int> sensorNumbers, mystr operationMinMaxAverage, Vector weights = Vector(0,1,0.));
	virtual mystr GetTypeName() {return "MultipleSensor";};	
	virtual mystr GetSensorTypeTexDescription(); 	// for the auto generated documentation 

	virtual void CopyFrom(const Sensor& s)
	{
		Sensor::CopyFrom(s);
		MultipleSensor & S = (MultipleSensor&)s;
		sensorNrs = S.sensorNrs;
		operation = S.operation;
	}
	
	virtual Sensor * GetCopy()
	{
		MultipleSensor * S = new MultipleSensor();
		S->CopyFrom(*this);
		return S;
	}

	virtual bool IsConsistent(mystr & errorStr);
	
	virtual void GetElementData(ElementDataContainer& edc);
	virtual int SetElementData(ElementDataContainer& edc); 
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter


	virtual int GetOperationFlagFromTextualIdentifier(mystr operationMinMaxAverage);
	virtual mystr TextualIdentifierOfOperation();

	virtual double GetCurrentValue(double time);

	virtual int GetNumberOfDrawingPositions() { return 1; }
	virtual Vector3D GetDrawPosition(int i)
	{
		return GetSensor(1).GetDrawPosition(i);
	}
	virtual int GetNumberOfRelatedSensors() { return sensorNrs.Length(); }			//$ DR 2013-05-21
	virtual int& GetRelatedSensorNumber(int nSensor) { return sensorNrs(nSensor); } //$ DR 2013-05-21


};//$EDC$[endclass,MultipleSensor]


// added by DR
class LoadSensor : public Sensor //$EDC$[beginclass,classname=LoadSensor,parentclassname=Sensor]
{
	int loadNumber; //$EDC$[varaccess,EDCvarname="load_number",EDCfolder="",tooltiptext="number of the load, to which the sensor is applied"]

protected:
	int& GetLoadNumber() { return loadNumber; }
	MBSLoad& GetLoad() { return GetMBS()->GetLoad(GetLoadNumber()); }		//$ DR 2012-09-12 changed GetElement(elementNumber) to GetElement(GetElementNumber())

	// default constructor is just for copy making
	LoadSensor() {}
	
public:
	LoadSensor(MBS * mbs) : Sensor(mbs)
	{
	}

	virtual void SetLoadSensor(int load_nr)
	{
		loadNumber = load_nr;
		AfterSetFunction(); //$ PG 2013-1-16: call added, in order to be consistent with FieldVariableElementSensor (see , e.g., FieldVariableElementSensor::SetFVESPos3D)
	}

	virtual void CopyFrom(const Sensor& s)
	{
		Sensor::CopyFrom(s);
		LoadSensor & S = (LoadSensor&)s;
		loadNumber = S.loadNumber;
	}

	virtual Sensor * GetCopy()
	{
		LoadSensor * S = new LoadSensor();
		S->CopyFrom(*this);
		return S;
	}

	virtual bool IsConsistent(mystr & errorStr)
	{
		if(loadNumber <= 0 || loadNumber > GetMBS()->NLoads())
		{
			errorStr = mystr("Sensor ") +mystr(GetOwnNum()) +mystr(" has an invalid load number: ") + mystr(loadNumber) + "\n";
			return false;
		}
		return true;
	}

	virtual double GetCurrentValue(double time);

	//virtual int GetNumberOfRelatedElements() { return 1; }
	//virtual int& GetRelatedElementNumber(int nElement) { return elementNumber; }
	virtual void GetElementData(ElementDataContainer & edc)
	{
		Sensor::GetElementData(edc);
	}
	virtual int SetElementData(ElementDataContainer& edc) 
	{
		int rv = Sensor::SetElementData(edc);
		AfterSetFunction();	// DR 2012-06-25 added
		return rv;
	}
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual int GetNumberOfDrawingPositions() { return 0; } // if just the load is known, but no element, it is not possible to draw the load
	virtual Vector3D GetDrawPosition(int i)		{	return Vector3D(0);	}

	virtual mystr GetTypeName(){return "LoadSensor"; }
	virtual mystr GetSensorTypeTexDescription(); 	// for the auto generated documentation 

	virtual int GetNumberOfRelatedLoads() { return 1; }			//$ DR 2013-05-21
	virtual int& GetRelatedLoadNumber(int nLoad) {return loadNumber;} //$ DR 2013-05-21

}; //$EDC$[endclass,LoadSensor]


// added by PG
// this is class for sensors, which are attached to a particular MBS quantity such as a global DOF, global Eigenvalue, iteration numbers, or scalar performance values.
// it should be constantly extended by all hotint developers by modifying
// 1) typedef enum {...} Object,
// 2) mystr SystemSensor::GetObjectString(Object object), and
// 3) SystemSensor::Object SystemSensor::GetObject(const mystr& object_str)
// 4) double SystemSensor::GetCurrentValue(double time)
// 5) EDCAuto-comment regarding the class member mystr object_str
class SystemSensor : public Sensor //$EDC$[beginclass,classname=SystemSensor,parentclassname=Sensor]
{
public:

// the following typedef defines all possible objects being tracked by the system sensor
// it should be constantly extended by all hotint users in particular accordance with
// 1) mystr SystemSensor::GetObjectString(Object object), and
// 2) SystemSensor::Object SystemSensor::GetObject(const mystr& object_str)
// 3) double SystemSensor::GetCurrentValue(double time)
	typedef enum
	{
		None,
		RHSEvaluations,
		RHSEvaluationsJacobian,
		Jacobians,
		NewtonIterations,
		DiscontinuousIterations,
		DOF,
		EV,
		Unknown
	} Object;

	SystemSensor(MBS * mbs) : Sensor(mbs)
	{
		SetObject(None);
		global_index = 0;
		SetSensorName(mystr("Systemsensor"));
	}
	virtual void CopyFrom(const Sensor& s)
	{
		Sensor::CopyFrom(s);
		SystemSensor & S = (SystemSensor&)s;
		global_index = S.global_index;
		object = S.object;
		object_str = S.object_str;
	}
	virtual Sensor * GetCopy()
	{
		SystemSensor * S = new SystemSensor();
		S->CopyFrom(*this);
		return S;
	}
	virtual bool IsConsistent(mystr & errorStr);
	virtual mystr GetTypeName(){return "SystemSensor"; }
	virtual mystr GetSensorTypeTexDescription();
	virtual void GetElementData(ElementDataContainer & edc)
	{
		Sensor::GetElementData(edc);
	}
	virtual int SetElementData(ElementDataContainer& edc) //$ DR 2012-10 return value changed from bool to int
	{
		int rv = Sensor::SetElementData(edc);
		SetObject(GetObject(object_str));
		AfterSetFunction();
		return rv;
	}
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual bool IsSensorOnEigensystem() { return (object == EV); }

	// main operation functions: setting and evaluating
	virtual void SetSystemSensor(Object object, int global_index=0);
	virtual double GetCurrentValue(double time);

protected:
	// internal access functions
	int GetGlobalIndex() { return global_index; }
	void SetGlobalIndex(int global_index);
	Object GetObject() { return object; }
	void SetObject(Object object); // sets object and object_str at once
	bool IsGlobalIndexObject() { return (object==EV || object==DOF); }

	// conversion between object and object_str
	mystr GetObjectString(Object object);
	Object GetObject(const mystr& object_str);
	
	SystemSensor() {}  // just for copy making
	
private:
	mystr object_str; //$EDC$[varaccess,EDCvarname="object",EDCfolder="",tooltiptext="Object tracked by systemsensor. Is either 'DOF' (global degree of freedom), 'EV' (global eigenvalue), 'jacobians', 'newton_iterations', 'discontinuous_iterations', 'rhs_evaluations', or 'rhs_evaluations_jacobian'"]
	Object object;
	int global_index; //$EDC$[varaccess,EDCvarname="global_index",EDCfolder="",tooltiptext="Number of the global index. Has to be set if (and only if) object is 'DOF' or 'EV'."]

};//$EDC$[endclass,SystemSensor]