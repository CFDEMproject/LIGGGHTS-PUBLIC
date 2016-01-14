//#**************************************************************
//#
//# filename:             sensorsSpecific.cpp
//#
//# author:               Gerstmayr Johannes
//#												Vetyukov Yury
//#
//# generated:						February-June 2012
//# description:          
//#                       
//# remarks:						  HotInt sensors - specific sensors
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
#include "sensorsSpecific.h"
#include "elementdataaccess.h"
#include "timeintlog.h"

double FieldVariableElementSensor::GetCurrentValue(double time)
{
	if(localNodeNumber == 0)
	{
		//if(!localPosition2D)
		//	return GetElement().GetFieldVariableValue(fieldVariableDescriptor, localPosition, false);
		//else
		//	return GetElement().GetFieldVariableValue(fieldVariableDescriptor, (const Vector2D &)localPosition, false);
		if(dimension==3)
			return GetElement().GetFieldVariableValue(fieldVariableDescriptor, localPosition, false);
		else if(dimension==2)
			return GetElement().GetFieldVariableValue(fieldVariableDescriptor, (const Vector2D &)localPosition, false);
		else
			return GetElement().GetFieldVariableValue(fieldVariableDescriptor, (const double &)localPosition.X(), false);
	}
	return GetElement().GetFieldVariableValue(fieldVariableDescriptor, localNodeNumber, false);
}

void FieldVariableElementSensor::SetFVESPos3D(int elementNumber, FieldVariableDescriptor fieldVariableDescriptor_, Vector3D localPosition)
{
	SetElementNumber(elementNumber);
	fieldVariableDescriptor = fieldVariableDescriptor_;
	localNodeNumber = 0;
	this->localPosition = localPosition;
	//localPosition2D = false;
	dimension = 3;
	AfterSetFunction();
}

void FieldVariableElementSensor::SetFVESPos2D(int elementNumber, FieldVariableDescriptor fieldVariableDescriptor_, Vector2D localPosition)
{
	SetElementNumber(elementNumber);
	fieldVariableDescriptor = fieldVariableDescriptor_;
	localNodeNumber = 0;
	this->localPosition.X() = localPosition.X();
	this->localPosition.Y() = localPosition.Y();
	this->localPosition.Z() = 0;
	//localPosition2D = true;
	dimension = 2;
	AfterSetFunction();
}

void FieldVariableElementSensor::SetFVESNode(int elementNumber, FieldVariableDescriptor fieldVariableDescriptor_, int localNodeNumber)
{
	//$ YV 2012-10-22 with a node number we do not need to differ between 2D and 3D, as we use the same function from the element to retrieve the value
	SetElementNumber(elementNumber);
	fieldVariableDescriptor = fieldVariableDescriptor_;
	this->localNodeNumber = localNodeNumber;
	// AD: include dimension of node here ?
	AfterSetFunction();
}

void FieldVariableElementSensor::CopyFrom(const Sensor& s)
{
	ElementSensor::CopyFrom(s);
	FieldVariableElementSensor & S = (FieldVariableElementSensor&)s;
	fieldVariableDescriptor = S.fieldVariableDescriptor;
	localPosition = S.localPosition;
	//localPosition2D = S.localPosition2D;
	localNodeNumber = S.localNodeNumber;
	dimension = S.dimension;
}

mystr FieldVariableElementSensor::GetTypeName()
{
	return "FVElementSensor"; //$ DR 2012-10 changed in order to be consistent with CEDCParser
	//if(localNodeNumber > 0)
	//	return "Element_" + mystr(GetElementNumber()) + "-" + "Node_" + mystr(localNodeNumber) + "-" + fieldVariableDescriptor.GetTextualIdentifier();
	//return "Element_" + mystr(GetElementNumber()) + "-" + fieldVariableDescriptor.GetTextualIdentifier();
}

mystr FieldVariableElementSensor::GetSensorTypeTexDescription()	// for the auto generated documentation //$ DR 2012-10 added
{
	mystr descr = mystr("The FieldVariableElementSensor evaluates the value of a field variable at a specified position. There are two possibilities to define this position: \n")
		+ mystr("\\begin{itemize} \n")
		+ mystr("   \\item element number + local position \n")
		+ mystr("   \\item element number + local node number \n")
		+ mystr("\\end{itemize} \n")
		+ mystr("The descriptions of the elements above include a list of available field variables for each element. ")
		+ mystr("Possible field variables are e.g.")
		+ mystr("\\begin{itemize} \n")
		+ mystr("   \\item position, velocity and displacement \n")
		+ mystr("   \\item bryant\\_angle and angular\\_velocity \n")
		+ mystr("   \\item beam\\_axial\\_extension, beam\\_torsion, beam\\_curvature \n")		
		+ mystr("   \\item many more \n")		
		+ mystr("\\end{itemize} \n");
	return descr;
}

bool FieldVariableElementSensor::IsConsistent(mystr & errorStr)
{
	if(!ElementSensor::IsConsistent(errorStr))
		return false;

	// AD: 2013-02-19: included checks for SetFVESNode
	if (localNodeNumber == 0)
	{ // Element & posiiton 
		//if( !( (GetElement().Dim() == 2 && localPosition2D) || (GetElement().Dim() == 3 && !localPosition2D)) )
		//{
		//	errorStr = "Body and sensor dimension are different (2D <==> 3D)\n";
		//	return false;
		//}
	}
	else
	{ // element & node
		return true; // dimension of node is not checked against dimension of element...
	}

	TArray<FieldVariableDescriptor> availableFieldVariables;
	GetElement().GetAvailableFieldVariables(availableFieldVariables);
	if(! (availableFieldVariables.Find(fieldVariableDescriptor)))
	{
		errorStr = "The FieldVariable '" + mystr(fieldVariableDescriptor.GetTextualIdentifier()) + "' is not available in Element " + mystr(GetElementNumber()) +"!\n";
		return false;
	}

	return true;
}

void FieldVariableElementSensor::GetElementData(ElementDataContainer& edc)
{
	ElementSensor::GetElementData(edc);

	ElementData ed;
	ed.SetText(fieldVariableDescriptor.GetTextualIdentifierWithoutComponents(), "field_variable"); ed.SetToolTipText("name of the field variable, e.g. 'position', see the documentation of the elements for the available field variables"); edc.Add(ed);
	ed.SetText(fieldVariableDescriptor.GetTextualIdentifierComponentsOnly(), "component");  ed.SetToolTipText("component of the field variable, e.g. 'x'"); edc.Add(ed);
}

int FieldVariableElementSensor::SetElementData(ElementDataContainer& edc) //$ DR 2012-10 return value changed from bool to int
{
	int rv = 1;
	rv = ElementSensor::SetElementData(edc);

	mystr FV_temp;
	GetElemDataText(GetMBS(),edc,"field_variable",FV_temp,1);

	mystr FVcomp_temp = mystr("");
	GetElemDataText(GetMBS(),edc,"component",FVcomp_temp,1);

	fieldVariableDescriptor.InitializeFromTextualIdentifier(FV_temp + mystr(FieldVariableDescriptor::GetComponentsDelimiter()) + FVcomp_temp);
	
	AfterSetFunction();	// DR 2012-06-25 added
	return rv;
}




// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// SingleElementDataSensor
// 2012-11 added by YV
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


double SingleElementDataSensor::GetCurrentValue(double time)
{
	RWdata.time = time;
	if (successfully_converted) {GetElement().ReadSingleElementData(RWdata);}
	else {RWdata.value = 0;}
	return RWdata.value;
}

void SingleElementDataSensor::SetSingleElementDataSensor(int elementNumber, const char * variable_name)
{
	SetElementNumber(elementNumber);

	successfully_converted = RWdata.GetVariableNameAndComponents(variable_name); //convert variable name into name and components in RWdata structure
	RWdata.RWaccess = TRWElementDataRead; //JG 2013-01-11: we only want to read values in sensors

	//if (!successfully_converted) 
	//{
	//	UO() << "ERROR: SingleElementDataSensor: variable name='" << variable_name << "' is not legal. Sensor returns 0.\n";
	//}

	AfterSetFunction();   //$ PG 2013-1-16: call added, in order to be consistent with FieldVariableElementSensor (see , e.g., FieldVariableElementSensor::SetFVESPos3D)
}

bool SingleElementDataSensor::IsConsistent(mystr & errorStr)
{
	if(!ElementSensor::IsConsistent(errorStr))
	{
		return false;
	}
	if(!successfully_converted)
	{
		errorStr = mystr("Variable of Sensor ") +mystr(GetOwnNum()) +mystr(" could not be converted to a ReadWriteElementDataVariableType data structure");
		return false;
	}
	TArrayDynamic<ReadWriteElementDataVariableType> available_variables;
	GetElement().GetAvailableSpecialValues(available_variables);

	for(int i = 1; i <= available_variables.Length(); i++)
	{		
		if(available_variables(i).IsAvailable(RWdata))
		{
			return true;
		}
	}
	errorStr = "Variable " + RWdata.variable_name + " is not supported by the specified element";
	return false;
}

void SingleElementDataSensor::CopyFrom(const Sensor& s)
{
	ElementSensor::CopyFrom(s);
	SingleElementDataSensor & S = (SingleElementDataSensor&)s;
	successfully_converted = S.successfully_converted;
	RWdata = S.RWdata;
}

//$ DR 2012-11-06 added
void SingleElementDataSensor::GetElementData(ElementDataContainer& edc)	//fill in all element data 
{
	ElementSensor::GetElementData(edc);

	mystr text = RWdata.variable_name;
	if(RWdata.comp1)	
	{
		text+=mystr("[")+mystr(RWdata.comp1);
		if(RWdata.comp2)
		{
			text+= mystr(",")+mystr(RWdata.comp2);
		}
		text+=mystr("]");
	}
	ElementData ed;
	ed.SetText(text,"value"); 
	ed.SetToolTipText("special value of the element, use ''[]'' to access vector or matrix values, e.g. force[1] or stress[2,3]"); 
	edc.TreeAdd("",ed);
}

//$ DR 2012-11-06 added
int SingleElementDataSensor::SetElementData(ElementDataContainer& edc) //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
{
	int rv = 1;
	rv = ElementSensor::SetElementData(edc);

	mystr full_string="";
	rv = GetElemDataText(GetMBS(), edc, "value",full_string, 1);
	if(rv)
	{
		rv=RWdata.GetVariableNameAndComponents(full_string);
		successfully_converted = rv;
	}

	AfterSetFunction();	// DR 2012-06-25 added
	return rv;
}

//$ DR 2012-11-06 added
mystr SingleElementDataSensor::GetSensorTypeTexDescription()	// for the auto generated documentation //$ DR 2012-10 added
{
	mystr descr = mystr("The ElementSensor returns special values evaluated in the element. It can be used e.g. for measuring a specific degree of freedom of an element. \n")
		+ mystr("The descriptions of the elements above include a list of available special values for each element.");
	return descr;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// MultipleSensor
// 2012-11 added by DR
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// sensor that consists of a list of sensors
// this replaces the old MBSMultipleSensor
void MultipleSensor::SetMultipleSensor(TArray<int> sensorNumbers, mystr operationMinMaxAverage, Vector weights_i)
{
	sensorNrs = sensorNumbers;
	operation = GetOperationFlagFromTextualIdentifier(operationMinMaxAverage);
	weights = weights_i;
	AfterSetFunction(); //$ PG 2013-1-16: call added, in order to be consistent with FieldVariableElementSensor (see , e.g., FieldVariableElementSensor::SetFVESPos3D)
}

int MultipleSensor::GetOperationFlagFromTextualIdentifier(mystr operation_str)
{
	// this list must coincide with the function GetCurrentValue
	if(operation_str.Compare("average")) { return 1;}
	if(operation_str.Compare("minimum")) { return 2;}
	if(operation_str.Compare("maximum")) { return 3;}
	if(operation_str.Compare("sum"))		 { return 4;}
	else return 0;
}

mystr MultipleSensor::TextualIdentifierOfOperation()
{
	// this list must coincide with the function GetOperationFlagFromTextualIdentifier
	switch(operation)
	{
	case 1:	// average
		return "average";
		break;
	case 2:	// minimum
		return "minimum";
		break;
	case 3:	// maximum
		return "maximum";
		break;
	case 4:	// sum
		return "sum";
		break;
	}
	return "";
}

mystr MultipleSensor::GetSensorTypeTexDescription()	// for the auto generated documentation 
{
	mystr descr = mystr("The MultipleSensor applies mathematical operations to a list of sensors. \n");
	descr+= mystr("The sensor can be used, e.g. to get the maximum or average value of a list of sensors. \n");
	descr+= mystr("The following mathematical operations are possible (use these words for 'operation'): \\\\");
	descr+= mystr("\n \\begin{itemize} \n");

	mystr op = "start";
	int i=1;
	while(!op.Compare(""))
	{
		operation = i;
		op = TextualIdentifierOfOperation();
		i++;
		if(!op.Compare(""))
		{
			descr+= mystr("  \\item ") +op + mystr("\n");
		}
	}
	descr+= mystr("\\end{itemize} \n");
	descr+= mystr("If weights are used, than the value of each sensor is multiplied with the weight before the mathematical operation is performed.\n");
	descr+= mystr("To compute a weighted sum of the first 4 sensors, the entries would be e.g. sensor\\_numbers = \[1,2,3,4\] and weights = \[0.125,0.125,0.25,0.5\].\n");
	return descr;
}

double MultipleSensor::GetCurrentValue(double time)
{
	bool use_weights = weights.Length();
	double w=1.;

	// initialize with first value
	if(use_weights){w=weights(1);}
	double value = w*this->GetSensor(1).GetCurrentValueWithSensorProcessing(time);

	switch(operation)
	{
	case 1:	// average
		for(int i=2; i<=sensorNrs.Length(); i++)
		{
			if(use_weights){w=weights(i);}
			value += w*this->GetSensor(i).GetCurrentValueWithSensorProcessing(time);
		}
		value /= (sensorNrs.Length());	
		break;
	case 2:	// minimum
		for(int i=2; i<=sensorNrs.Length(); i++)
		{
			if(use_weights){w=weights(i);}
			value = min(value,w*this->GetSensor(i).GetCurrentValueWithSensorProcessing(time));
		}
		break;
	case 3:	// maximum
		for(int i=2; i<=sensorNrs.Length(); i++)
		{
			if(use_weights){w=weights(i);}
			value = max(value,w*this->GetSensor(i).GetCurrentValueWithSensorProcessing(time));
		}
		break;
	case 4:	// sum
		for(int i=2; i<=sensorNrs.Length(); i++)
		{
			if(use_weights){w=weights(i);}
			value += w*this->GetSensor(i).GetCurrentValueWithSensorProcessing(time);
		}
		break;
	}

	return value;
}

bool MultipleSensor::IsConsistent(mystr & errorStr)
{
	if(sensorNrs.Length()<1)
	{
		errorStr = "MultipleSensor: There has to be at least one sensor in the list!\n";
		return false;
	}

	for(int i=1; i<=sensorNrs.Length(); i++)
	{
		if((sensorNrs(i)<1) || (sensorNrs(i)> GetMBS()->NSensors()) || (sensorNrs(i)==GetOwnNum()))
		{
			errorStr = mystr("MultipleSensor: It is not possible to use sensor ") +mystr(sensorNrs(i)) +mystr(" for the MultipleSensor!\n");
			return false;
		}
	}

	if(!((weights.Length()==0) || (weights.Length() == sensorNrs.Length())))
	{
		errorStr = mystr("MultipleSensor: The length of weights must be equal to the length of the sensor_numbers or zero!\n");
		return false;	
	}


	return true;
}

void MultipleSensor::GetElementData(ElementDataContainer & edc)
{
	Sensor::GetElementData(edc);

	ElementData ed;
	ed.SetText(TextualIdentifierOfOperation() ,"operation"); ed.SetToolTipText("mathematical operation that is applied to the sensor values, e.g. 'maximum','average',..."); edc.Add(ed);

}
int MultipleSensor::SetElementData(ElementDataContainer& edc) 
{
	int rv = 1;
	rv = Sensor::SetElementData(edc);

	mystr temp;
	GetElemDataText(GetMBS(),edc,"operation",temp,1);
	operation = GetOperationFlagFromTextualIdentifier(temp);
	
	AfterSetFunction();	// DR 2012-06-25 added
	return rv;

}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// LoadSensor
// 2012-11 added by DR
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double LoadSensor::GetCurrentValue(double time)
{
	if(GetLoadNumber())
	{
		return GetLoad().Evaluate(time);
	}
	return -1;
}

mystr LoadSensor::GetSensorTypeTexDescription()	// for the auto generated documentation 
{
	mystr descr = mystr("The LoadSensor can be applied to loads in order to measure their time dependency. \n");
	descr+= mystr("The value $F(t)$ of a load at time $t$ is computed (see the description of the loads for more details) as: \n");
	descr+= mystr("	\\begin{equation} \n");
	descr+= mystr("	   F(t) = f(t) \\vec{F} \n");
	descr+= mystr("	\\end{equation} \n");
	descr+= mystr("The LoadSensor returns the value of the factor $f(t)$ and not the value $F(t)$. ");
	descr+= mystr("If the LoadSensor is used for a scalar load (e.g. GCLoad), than $f(t)$ and $F(t)$ are equal. If the LoadSensor is used for a load vector (e.g. ForceVector3D) than $f(t)$ and $F(t)$ may not be equal.  \\\\");
	descr+= mystr("The LoadSensor can not be shown in the graphical output, because the load does not have a position by itself and may be applied to several elements or nodes.");

	return descr;
}



// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// SystemSensor
// 2013-11 added by PG
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void SystemSensor::SetSystemSensor(Object object, int global_index)
{
	SetObject(object);
	SetGlobalIndex(global_index);
	AfterSetFunction();

	mystr name("Systemsensor");
	if (object != None && object != Unknown)
	{
		name += mystr("-")+object_str;
		if (IsGlobalIndexObject())
		{
			name += mystr(global_index);
		}
	}
	SetSensorName(name);
}

void SystemSensor::SetObject(Object object)
{
	this->object = object;
	if (object != Unknown)
	{
		this->object_str = GetObjectString(object);
	}
}

void SystemSensor::SetGlobalIndex(int global_index)
{
	this->global_index = global_index;
}

mystr SystemSensor::GetSensorTypeTexDescription()	// for the auto generated documentation 
{
	mystr descr = mystr("The SystemSensor can be applied to global degrees of freedom, eigenvalues, several iteration numbers or performance indicators. \n");
	descr+= mystr("It returns the value of the specified quantity at time $t$, ");
	descr+= mystr("and can not be shown in the graphical output, because a system quantity does in general not have a position by itself.");

	return descr;
}

bool SystemSensor::IsConsistent(mystr & errorStr)
{
	if(object == Unknown)
	{
		errorStr = mystr("Sensor ") +mystr(GetOwnNum()) + mystr(" has an unknown object '") + object_str + "'\n";
		return false;
	}

	if(IsGlobalIndexObject())
	{
		if ((global_index <= 0)
			//|| (object==EV && ??? size of eigenvalue vector not known at the time this function is called)
			|| (object==DOF && global_index > GetMBS()->GetXact().Length()))
		{
			errorStr = mystr("Sensor ") +mystr(GetOwnNum()) + mystr(" has an invalid global index number: ") + mystr(global_index) + "\n";
			return false;
		}
	}
	else
	{
		if (global_index != 0)
		{
			GetMBS()->UO(UO_LVL_warn) << "WARNING: SystemSensor: global index was set to " << global_index << ", but object '" << object_str << "' does not need a global index!\n";
		}
	}
	return true;
}

// the following two functions are converting object to object_str and vice versa.
// they should be constantly extended by all hotint users in particular accordance with
// 1) typedef enum {...} Object (see class definition in header), and
// 2) double SystemSensor::GetCurrentValue(double time)
mystr SystemSensor::GetObjectString(Object object)
{
	switch(object)
	{
	case None:
		return mystr("none");
	case NewtonIterations:
		return mystr("newton_iterations");
	case DiscontinuousIterations:
		return mystr("discontinuous_iterations");
	case Jacobians:
		return mystr("jacobians");
	case RHSEvaluationsJacobian:
			return mystr("rhs_evaluations_jacobian");
	case RHSEvaluations:
		return mystr("rhs_evaluations");
	case DOF:
		return mystr("DOF");
	case EV:
		return mystr("EV");		
	}

	return mystr("unknown");
}
SystemSensor::Object SystemSensor::GetObject(const mystr& object_str)
{
	if (object_str.Compare("none", 1))
	{
		return None;
	}
	else if (object_str.Compare("newton_iterations", 1))
	{
		return NewtonIterations;
	}
	else if (object_str.Compare("discontinuous_iterations", 1))
	{
		return DiscontinuousIterations;
	}
	else if (object_str.Compare("jacobians", 1))
	{
		return Jacobians;
	}
	else if (object_str.Compare("rhs_evaluations_jacobian", 1))
	{
		return RHSEvaluationsJacobian;
	}
	else if (object_str.Compare("rhs_evaluations", 1))
	{
		return RHSEvaluations;
	}
	else if (object_str.Compare("DOF", 1))
	{
		return DOF;
	}
	else if (object_str.Compare("EV", 1))
	{
		return EV;
	}
	else
	{
		GetMBS()->UO(UO_LVL_warn) << "WARNING: SystemSensor: object '" << object_str << "' unknown!\n";
	}

	return Unknown;
}

// the following function is called via Sensor::Evaluate(..) or Sensor::GetCurrentValueWithSensorProcessing(..)
// at the end of each time/load step or after an eigenvalue computation, i.e. whenever sensor values are written to memory or solution files.
// the function should be constantly extended by all hotint users in particular accordance with
// 1) typedef enum {...} Object (see class definition in header),
// 2) mystr SystemSensor::GetObjectString(Object object), and
// 3) SystemSensor::Object SystemSensor::GetObject(const mystr& object_str)
double SystemSensor::GetCurrentValue(double time)   
{
	switch(object)
	{
	case NewtonIterations:
		return GetMBS()->GetTimeIntLog().TInewtonit;
	case DiscontinuousIterations:
		return GetMBS()->GetTimeIntLog().TInonlinit;
	case RHSEvaluations:
		return GetMBS()->GetTimeIntLog().evalfcnt;
	case RHSEvaluationsJacobian:
		return GetMBS()->GetTimeIntLog().evalf_jaccnt;
	case Jacobians:
		return GetMBS()->GetTimeIntLog().jaccount;
	case DOF:
		return GetMBS()->GetXact(global_index);
	case EV:
		return GetMBS()->GetEigenValue(global_index);
	}

	return 0.;
}

