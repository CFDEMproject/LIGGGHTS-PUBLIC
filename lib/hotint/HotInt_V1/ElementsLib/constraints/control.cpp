//#**************************************************************
//#
//# filename:             control.cpp
//#
//# author:               Gerstmayr Johannes, Rafael Ludwig
//#
//# generated:						17.October 2004
//# description:          Driver and model for timeintegration
//#                       Model of a rigid arm and a hydraulic zylinder
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
 
#include "element.h"
#include "constraint.h"
#include "sensors.h"
#include "control.h"
#include "elementdataaccess.h"
#include "graphicsconstants.h"
#include "rendercontext.h"

const double DEFAULT_IOE_SIZE = 20;

//! AD: these strings are the default SYMBOL TEXTs, used as label in the 2D IO Blocks Window
//! AD: shortened names to fit nicer into the element rectangle
const char* CSText_InputOutputElement = "IOElement";
const char* CSText_LinearTransformation = "y=Ax+b"; //"IOLinearTransformation";
const char* CSText_STransferFunction = "S-Fcn"; //"IOContinuousTransferFunction";
const char* CSText_LinearODE = "lin ODE"; //"IOLinearODE";
const char* CSText_IOTime = "Time"; //"IOTime";
const char* CSText_IODisplay = "..."; //"IODisplay";
const char* CSText_IOResponseElement = "Response";
const char* CSText_IOKeyResponseElement = "Key";
const char* CSText_IOMouseResponseElement = "Mouse";
const char* CSText_InputOutputElementDiscrete = "Discr."; //"IODiscrElement";
const char* CSText_ZTransferFunction = "Z-Fcn"; //"IODiscreteTransferFunction";
const char* CSText_RandomSource = "Random"; //"IORandomSource";
const char* CSText_Quantizer = "Quant."; //"IOQuantizer";
const char* CSText_SDActor = "SDActor";
const char* CSText_ControllerInterface = "Contr."; //"ControllerInterface"; //$ RL 2010-04
const char* CSText_IOSaturate = "Saturate"; //"IOSaturate"; //$ RL 2011-01
const char* CSText_IODeadZone = "DeadZone"; //"IODeadZone"; //$ RL 2011-01
const char* CSText_IOProduct = "X"; //"IOProduct"; //$ RL 2011-05
const char* CSText_IOPulseGenerator = "Pulse"; //"IOPulseGenerator"; //$ RL 2011-07
const char* CSText_IOStopComputation = "Stop"; //"IOStopComputation";//$ RL 2012-2-1: //$ MS 2012-2-1: 
const char* CSText_IOElementDataModifier = "Modify"; 
const char* CSText_IOTimeWindow = "TWindow"; //"IOTimeWindow";
const char* CSText_IOMathFunction = "f(x)";
const char* CSText_IOMinMax = "MinMax";
const char* CSText_IOTCPIPBlock = "TCPIP"; 

const char* Quantizer::SymbolText() const 
{
	return CSText_Quantizer;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// DISCRETE SYSTEMS DISCRETE SYSTEMS DISCRETE SYSTEMS DISCRETE SYSTEMS DISCRETE SYSTEMS DISCRETE SYSTEMS DISCRETE SYSTEMS DISCRETE SYSTEMS

const char* InputOutputElementDiscrete::SymbolText() const 
{
	return CSText_InputOutputElementDiscrete;
}

double InputOutputElementDiscrete::Roundval(double x) const
{
	double d = 0.;
	if (modf(x,&d)>=.5)
		return x>=0?ceil(x):floor(x);
	else
		return x<0?ceil(x):floor(x);
}

double InputOutputElementDiscrete::GetDiscreteTime() const 
{
	double tol = discrete_time_tol_inv;
	double val = (double)GetK()*deltaT + toff;
	val = Roundval(val*tol)/tol;
	return val;
}
double InputOutputElementDiscrete::GetNextDiscreteTime() const 
{
	double tol = discrete_time_tol_inv;
	double val = (double)(GetK()+1.)*deltaT + toff;
	val = Roundval(val*tol)/tol;
	return val;
}

void InputOutputElementDiscrete::StartTimeStep()	
{
	double time_to_next_event = GetNextDiscreteTime() - GetMBS()->GetTime();
	if(GetMBS()->GetStepSizeNew() >= time_to_next_event)
	{
		if (time_to_next_event < GetMBS()->GetStepRecommendation())
		{
			GetMBS()->SetStepRecommendation(time_to_next_event);
		}
	}	
}
void InputOutputElementDiscrete::EndTimeStep()
{
	// call this function after TItime += TIstep and before the sensor data is written!
	if(isDiscreteEvent(GetNextDiscreteTime())) //t == (k+1)* deltaT?
	{
		UpdateDiscreteElementState();  // update state variables xk
		SetK(GetK()+1); // k --> k+1
		UpdateDiscreteElementInput(GetDiscreteTime());  // update uk	
	}
}

int InputOutputElementDiscrete::isDiscreteEvent(double discrete_time) const
{
	//if(GetMBS()->GetStepEndTime() >= (discrete_time-discrete_time_tol))return 1; // values of end of time step are computed --> update of discrete element
	if(GetMBS()->GetTime() >= (discrete_time-discrete_time_tol))return 1;	         // values of actual time step must be computed --> update of discrete element
	else return 0;	                                                               // solver endtime point before than next disctrete time point --> no update of discrete element
};

void InputOutputElementDiscrete::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	InputOutputElement::GetElementData(edc);
	ElementData ed; 
	ed.SetDouble(deltaT, "discrete_time_step"); ed.SetToolTipText("Sample time \"dT\""); /*edc.Add(ed);*/ edc.TreeAdd("IOBlock",ed); // change: MSax 28-02-2013
	ed.SetDouble(toff, "discrete_time_offset"); ed.SetToolTipText("Sample offset \"off\": Tk = k*dT + off"); /*edc.Add(ed);*/ edc.TreeAdd("IOBlock",ed); // change: MSax 28-02-2013
}

int InputOutputElementDiscrete::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = InputOutputElement::SetElementData(edc);
  GetElemDataDouble(mbs, edc, "IOBlock.discrete_time_step", deltaT, 0);
	GetElemDataDouble(mbs, edc, "IOBlock.discrete_time_offset", toff, 0);
	return rv;
}

void InputOutputElementDiscrete::DrawBlockSymbol() 
{ 
	// TODO: add specific symbol 
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
const char* ZTransferFunction::SymbolText() const 
{
	return CSText_ZTransferFunction;
}

int ZTransferFunction::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = InputOutputElementDiscrete::CheckConsistency(errorstr);
	if(rv){	return rv;}

	return rv;
}

void ZTransferFunction::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	InputOutputElementDiscrete::GetElementData(edc);

	ElementData ed;
	ed.SetVector(num.GetVecPtr(), num.Length(), "num_coeffs"); ed.SetVariableLength(); ed.SetToolTipText("Coefficients of numerator polynomial of z-function"); /*, b0+b1*z+b2*z\\^2+...+bn*z\\^n*/  edc.TreeAdd("IOBlock",ed); // change: MSax 28-02-2013
	ed.SetVector(den.GetVecPtr(), den.Length(), "den_coeffs"); ed.SetVariableLength(); ed.SetToolTipText("Coefficients of denominator polynomial of z-function"); /*, a0+a1*z+a2*z\\^2+...+am*z\\^m*/  edc.TreeAdd("IOBlock",ed); // change: MSax 28-02-2013

	//ed.SetVector(init_vec.GetVecPtr(), init_vec.Length(), "init_vec"); ed.SetVariableLength(); ed.SetToolTipText("Initial state vector.");  edc.TreeAdd("IOBlock",ed); // added: MSax 28-02-2013
	
	//SetElemDataVector(edc, num, "Num_coeffs"); edc.Last().SetVariableLength(); edc.Last().SetToolTipText("Coefficients of numerator polynomial of z-function, b0+b1*z+b2*z^2+...+bn*z^n");
	//SetElemDataVector(edc, den, "Den_coeffs"); edc.Last().SetVariableLength(); edc.Last().SetToolTipText("Coefficients of denominator polynomial of z-function, a0+a1*z+a2*z^2+...+am*z^m");

	//XXX: ElementData ed; 
	//XXX: ed.SetDouble(yk, "Initial_output"); ed.SetToolTipText("Initial output value."); edc.Add(ed);
}

int ZTransferFunction::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = InputOutputElementDiscrete::SetElementData(edc);

	Vector numI, denI;
	//double y0;
	//GetElemDataVector(mbs, edc, "Num_coeffs", numI, 1);
	//GetElemDataVector(mbs, edc, "Den_coeffs", denI, 1);

	GetElemDataVector(mbs,edc,"IOBlock.num_coeffs",numI,1);
	GetElemDataVector(mbs,edc,"IOBlock.den_coeffs",denI,1);
	//XXX: GetElemDataDouble(mbs, edc, "Initial_output", y0, 1);
	//GetElemDataVector(mbs,edc,"IOBlock.init_vec",init_vec,1); // $MSax 2013-03-01: added
	
	if (num.Length() != den.Length())
	{
		rv = 0; //not successful
	}
	else
	{
		SetZTransferFunction(numI, denI); //, init_vec); // $MSax 2013-03-01: changed --> added 3rd. input argument init_vec
		//XXX: SetZTransferFunction(numI, denI, y0);  // not implemented yet!
	}

	return rv;
}

void ZTransferFunction::DrawBlockSymbol() 
{ 
	return InputOutputElement::DrawBlockSymbol(); // $ MSax 2013-03-01: added
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void RandomSource::SetRandomSource(double amplitudeI, double offsetI, TRandomType methodI, double seedI, double init_val, int Nbits)
{
	InitRandomSource(amplitudeI, offsetI, methodI, seedI, init_val, Roundval(pow(2., Nbits)-1.));
}
//                                                                                                                     rand_max_i ... 2^N-1
void RandomSource::InitRandomSource(double amplitudeI, double offsetI, TRandomType methodI, double seedI, double init_val, int rand_maxi)
{
	amplitude = amplitudeI;
	offset = offsetI;
	seed = seedI;
	method = methodI;

	SetNOutputs(1);
	SetNStates(0);                	
	int dimIODiscrete = DataS(); // dimension of root element

	if(seed > 1.0)
	{
		seed = 1.0; // max. value of seed
		GetMBS()->UO(UO_LVL_warn) << "WARNING: InitRandomSource: \"seed\" bigger than 1.0 --> automatically set to 1.0\n";
	}
	if(seed < 0.0)
	{
		seed = 0.0; // min. value of seed
		GetMBS()->UO(UO_LVL_warn) << "WARNING: InitRandomSource: \"seed\" negative --> automatically set to 0.0\n";
	}	
	
	
	Vector init;
	init.SetLen(0); // new init vector inclusive base element variables
	init = init.Append(GetDataInit());
	init = init.Append(init_val); // x(t=0) = y(t=0)
  
	// initialize random signal generator	
	if(method == (TRandomType)Trand)
	{
		if(seed)srand((int)seed); //initialize random generator
		SetNDiscreteStates(1); // yk = xk
	}
	else //if(method == (TRandomType)TLFSR
	{				
		SetRandMax(rand_maxi);
		int random_value = (int)(seed * rand_maxi); 
		if(!random_value)
		{
			random_value++;                 // intialize state with value unequal 0
		}
		init = init.Append((double)random_value); // internal state of LSFR random generator, unequal 0.
		SetNDiscreteStates(2);            // yk = xk, state
	}
	SetDataInit(init);
}

// random integer number € [0, RAND_MAX]; method TLFSR returns numbers unequal zero
int RandomSource::GetRandomInt()
{
	if(method == (TRandomType)Trand)
	{
		return rand();
	}
	else //if method == (TRandomType)TLFSR
	{		
		int shift_left = (int)(log10((double)abs(~GetRandMax()))/log10(2.));			// NOT(7FFF) = 32768 = 2^15 = 1000 0000 0000 0000 ==> 16 = log2(2^15)+1 = log10(2^15)/log10(2) + 1
		unsigned int rand_val = (unsigned int)XDiscrete(2);												// get current internal state
		unsigned int newBit = 0;
		
		switch(GetRandMax())
		{
		case 1:  newBit = rand_val&1;		break;//Fibonacci LFSR, Bit1 ... 1 Bits
		case 3:  newBit = rand_val & (rand_val&1)^((rand_val&2) >> 1);	break;	//Fibonacci LFSR, (Bit1) XOR (Bit2) ... 2 Bits
		case 7:  newBit = rand_val & (rand_val&1)^((rand_val&2) >> 1);	break;	//Fibonacci LFSR, (Bit1) XOR (Bit2) ... 3 Bits
		case 15:  newBit = rand_val & (rand_val&1)^((rand_val&2)>> 1);	break;	//Fibonacci LFSR, (Bit1) XOR (Bit2) ... 4 Bits
		case 31:  newBit = rand_val & (rand_val&1)^((rand_val&4) >> 2);	break;	//Fibonacci LFSR, (Bit1) XOR (Bit3) ... 5 Bits
		case 63:  newBit = rand_val & (rand_val&1)^((rand_val&2) >> 1);	break;	//Fibonacci LFSR, (Bit1) XOR (Bit2) ... 6 Bits
		case 127:  newBit = rand_val & (rand_val&1)^((rand_val&2) >> 1);	break;	//Fibonacci LFSR, (Bit1) XOR (Bit2) ... 7 Bits
		case 255:  newBit = rand_val & (rand_val&1)^((rand_val&4) >> 2)^((rand_val&8) >> 3)^((rand_val&16) >> 4);	break;	//Fibonacci LFSR, (Bit1) XOR (Bit3) XOR (Bit4) XOR (Bit5)... 8 Bits
		case 511:  newBit = rand_val & (rand_val&1)^((rand_val&16) >> 4);	break;	//Fibonacci LFSR, (Bit1) XOR (Bit5) ... 9 Bits
		case 1023:  newBit =  rand_val & (rand_val&1)^((rand_val&8) >> 3);	break;	//Fibonacci LFSR, (Bit1) XOR (Bit4) ... 10 Bits
		case 2047:  newBit = rand_val & (rand_val&1)^((rand_val&4) >> 2);	break;	//Fibonacci LFSR, (Bit1) XOR (Bit3) ... 11 Bits
		case 4095:  newBit = rand_val & (rand_val&1)^((rand_val&2) >> 1)^((rand_val&16) >> 4)^((rand_val&64) >> 6);	break;	//Fibonacci LFSR, (Bit1) XOR (Bit2) XOR (Bit5) XOR (Bit7) ... 12 Bits
		case 8191:  newBit = rand_val & (rand_val&1)^((rand_val&2) >> 1)^((rand_val&8) >> 3)^((rand_val&128) >> 7);	break;	//Fibonacci LFSR, (Bit1) XOR (Bit2) XOR (Bit4) XOR (Bit8) ... 13 Bits
		case 16383:  newBit = rand_val & (rand_val&1)^((rand_val&2) >> 1)^((rand_val&8) >> 3)^((rand_val&32) >> 5);	break;	//Fibonacci LFSR, (Bit1) XOR (Bit2) XOR (Bit4) XOR (Bit6) ... 14 Bits
		case 32767:  newBit = rand_val & (rand_val&1)^((rand_val&2) >> 1);	break;	//Fibonacci LFSR, (Bit1) XOR (Bit2) ... 15 Bits
		case 65535:  newBit = rand_val & (rand_val&1)^((rand_val&4) >> 2)^((rand_val&8) >> 3)^((rand_val&32) >> 5);	break;	//Fibonacci LFSR, (Bit1) XOR (Bit3) XOR (Bit4) XOR (Bit6) ... 16 Bits
		case 131071:  newBit = rand_val & (rand_val&1)^((rand_val&8) >> 3);	break;	//Fibonacci LFSR, (Bit1) XOR (Bit4) ... 17 Bits
		case 262143:  newBit = rand_val & (rand_val&1)^((rand_val&128) >> 7);	break;	//Fibonacci LFSR, (Bit1) XOR (Bit8) ... 18 Bits
		case 524287:  newBit = rand_val & (rand_val&1)^((rand_val&2) >> 1)^((rand_val&4) >> 2)^((rand_val&32) >> 5);	break;	//Fibonacci LFSR, (Bit1) XOR (Bit2) XOR (Bit3) XOR (Bit6) ... 19 Bits
		case 1048575:  newBit = rand_val & (rand_val&1)^((rand_val&8) >> 3);	break;	//Fibonacci LFSR, (Bit1) XOR (Bit4) ... 20 Bits
		case 2097151:  newBit = rand_val & (rand_val&1)^((rand_val&4) >> 2);	break;	//Fibonacci LFSR, (Bit1) XOR (Bit3) ... 21 Bits
		case 4194303:  newBit = rand_val & (rand_val&1)^((rand_val&2) >> 1);	break;	//Fibonacci LFSR, (Bit1) XOR (Bit2) ... 22 Bits
		case 8388607:  newBit = rand_val & (rand_val&1)^((rand_val&32) >> 5);	break;	//Fibonacci LFSR, (Bit1) XOR (Bit6) ... 23 Bits
		case 16777215:  newBit = rand_val & (rand_val&1)^((rand_val&2) >> 1)^((rand_val&8) >> 3)^((rand_val&16) >> 4);	break;	//Fibonacci LFSR, (Bit1) XOR (Bit2) XOR (Bit4) XOR (Bit5) ... 24 Bits
		default: mbs->UO(UO_LVL_err) << "XOR-operation only implemented for maximal 24 Bits.";
		}

		newBit = newBit << shift_left;                                       // newBit € {1000 0000 0000 0000, 0000 0000 0000 0000}
		rand_val = (rand_val + newBit) >> 1;                                 // one shift left (doesn't work with rand_val = 0)
		XDiscrete(2) = (double)rand_val;                                     // store new internal state
		return rand_val;
	}
}

//function is called every discrete event
// xk+1 = f(xk,uk)
void RandomSource::UpdateDiscreteElementState()
{		
	int randval = GetRandomInt();
	assert(randval >= 0 && randval <= GetRandMax());
	XDiscrete(1) = amplitude*((double)randval/(double)GetRandMax() - 0.5) + offset;
}

const char* RandomSource::SymbolText() const 
{
	return CSText_RandomSource;
}

int RandomSource::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = InputOutputElementDiscrete::CheckConsistency(errorstr);
	if(rv){	return rv;}

	return rv;
}

void RandomSource::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	InputOutputElementDiscrete::GetElementData(edc);

	ElementData ed; 
	ed.SetDouble(amplitude, "max_amplitude"); ed.SetToolTipText("Max. amplitude of random value."); edc.TreeAdd("IOBlock",ed);
	ed.SetDouble(offset, "mean_value"); ed.SetToolTipText("Offset of random signal."); edc.TreeAdd("IOBlock",ed);
	
	ed.SetBool((TRandomType)method, "method"); ed.SetToolTipText("Random generator method."); edc.TreeAdd("IOBlock",ed);

	int bits = (int)Roundval(log10(rand_max+1.)/log10(2.)); // Bits = log2(max_rand+1)
	ed.SetInt(bits, "bits"); ed.SetToolTipText("Number of bits for random signal."); edc.TreeAdd("IOBlock",ed);
	
	ed.SetBool(isConstAmplitude, "constant_amplitude"); ed.SetToolTipText("Output values are +amplitude or -amplitude if flag is activate."); edc.TreeAdd("IOBlock",ed);

	ed.SetDouble(seed, "seed"); ed.SetToolTipText("seed € [0.,1.]... initialization of random generator"); edc.TreeAdd("IOBlock",ed); //$ MSax 2013-02-28: added

	ed.SetDouble(init_val, "init_val"); ed.SetToolTipText("initial value of the generator x(t=0) = y(t=0)"); edc.TreeAdd("IOBlock",ed); //$ MSax 2013-02-28: added
}

int RandomSource::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = InputOutputElementDiscrete::SetElementData(edc);


	GetElemDataDouble(mbs, edc, "IOBlock.max_amplitude", amplitude, 1);
	GetElemDataDouble(mbs, edc, "IOBlock.mean_value", offset, 1);
	
	int methodI;
	GetElemDataBool(mbs, edc, "IOBlock.method", methodI, 1);

	int bits;
	GetElemDataInt(mbs, edc, "IOBlock.bits", bits, 1);
		
	GetElemDataBool(mbs, edc, "IOBlock.constant_amplitude", isConstAmplitude, 1);
	
	GetElemDataDouble(mbs, edc, "IOBlock.seed", seed, 1); //$ MSax 2013-02-28: added

	GetElemDataDouble(mbs, edc, "IOBlock.init_val", init_val, 1); //$ MSax 2013-02-28: added


	// TODO: add init_data

	SetRandomSource(amplitude, offset, (TRandomType)methodI, seed, init_val, bits);
	SetConstantAmplitude(isConstAmplitude);

	inputs.SetLen(0);  //$ MSax 2013-04-03: added
	input_types.SetLen(0); //$ MSax 2013-04-03: added
	input_localnum.SetLen(0); //$ MSax 2013-04-03: added

	return rv;
}

void RandomSource::DrawBlockSymbol() 
{ 
	return InputOutputElement::DrawBlockSymbol();
}


// DISCRETE SYSTEMS DISCRETE SYSTEMS DISCRETE SYSTEMS DISCRETE SYSTEMS DISCRETE SYSTEMS DISCRETE SYSTEMS DISCRETE SYSTEMS DISCRETE SYSTEMS
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void InputOutputElement::CopyFrom(const Element& e)
{
	Constraint::CopyFrom(e);
	const InputOutputElement& ce = (const InputOutputElement&)e;

	inputs = ce.inputs;
	input_types = ce.input_types;
	input_localnum = ce.input_localnum;

	n_output = ce.n_output;
	n_state = ce.n_state;
	ref_pos=ce.ref_pos;
	draw_dim=ce.draw_dim;

	rotation = ce.rotation;
	input_nodes.SetLen(ce.input_nodes.Length());
	input_nodes_num.SetLen(ce.input_nodes_num.Length());
	for (int i=1; i <= ce.input_nodes.Length(); i++)
	{
		input_nodes(i) = ce.input_nodes(i);
		input_nodes_num(i) = ce.input_nodes_num(i);
	}
	colbackground = ce.colbackground;
	colforeground = ce.colforeground;
}

void InputOutputElement::InitConstructor()
{
	draw_element = 1;
	type = TMBSElement(TController+TConstraint);
	elementname = GetElementSpec();
	elements.SetLen(0);

	inputs.SetLen(0);
	input_localnum.SetLen(0);
	input_types.SetLen(0);
	n_output = 0;
	n_state = 0;

	x_init = Vector(SS());
	rotation = 0;
	input_nodes.SetLen(0);
	input_nodes_num.SetLen(0);

	//AddInput(1, 1); // $ MSax 2013-02-28: added

	//drawing properties:
	ref_pos = Vector2D(0.,0.);
	// default size for elements made much larger to represent pixels in 2D View
	draw_dim = Vector3D(DEFAULT_IOE_SIZE,DEFAULT_IOE_SIZE,0.);
	// draw_dim = Vector3D(1.,1.,0.);
	rotation = 0;
	
	colbackground = Vector3D(-1.,-1.,-1.); //-1. ... no background color
	colforeground = colblack;
}


void InputOutputElement::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	//if(1)
	//{
		Constraint::GetElementData(edc);
	//}
	//else 
	//{
	//  old:
	//	Constraint::GetElementData(edc);

	//	ElementData ed;

	//	ed.SetInt(GetNInputs(), "Number_of_inputs"); ed.SetLocked(1); edc.Add(ed);
	//	ed.SetInt(GetNOutputs(), "Number_of_outputs"); ed.SetLocked(1); edc.Add(ed);
	//	ed.SetInt(GetNStates(), "Number_of_states"); ed.SetLocked(1); edc.Add(ed);
	//	
	//	ed.SetVector2D(ref_pos.X(), ref_pos.Y(), "Drawing_ref_pos"); ed.SetToolTipText("reference position for drawing"); edc.Add(ed);
	//	ed.SetVector3D(draw_dim.X(), draw_dim.Y(), draw_dim.Z(), "Draw_size"); ed.SetToolTipText("drawing parameters: x=width, y=height, z=unused"); edc.Add(ed);

	//	SetElemDataIVector(edc, inputs, "Input_element_numbers"); edc.Last().SetVariableLength(); edc.Last().SetToolTipText("Only valid element numbers permitted!");
	//	SetElemDataIVector(edc, input_types, "Input_element_types"); edc.Last().SetVariableLength(); edc.Last().SetToolTipText("1=InputOutputElement, 2=Sensor");
	//	SetElemDataIVector(edc, input_localnum, "Input_local_number"); edc.Last().SetVariableLength(); edc.Last().SetToolTipText("Local number i = output number i of referred element");
	//}

	Matrix mat;
	mat.SetMatrix2n(input_nodes);

	ElementData ed;
	ed.SetVector2DList(mat.GetMatPtr(), input_nodes.Length(), "input_nodes");
	edc.TreeAdd("Graphics",ed);
}

int InputOutputElement::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	//InitConstructor();
	int rv = Constraint::SetElementData(edc);
	//int old = 0;
	//if(old)
	//{
	//	GetElemDataVector2D(mbs, edc, "Drawing_ref_pos", ref_pos, 0); 
	//	GetElemDataVector3D(mbs, edc, "Draw_size", draw_dim, 0); 

	//	GetElemDataIVector(mbs, edc, "Input_element_numbers", inputs, 1);
	//	GetElemDataIVector(mbs, edc, "Input_element_types", input_types, 1);
	//	GetElemDataIVector(mbs, edc, "Input_local_number", input_localnum, 1);
	//}

	double *mp;
	int n;
	edc.TreeGetVector2DList("Graphics.input_nodes",&mp,n);
	Matrix mat(n,2,mp);
	mat.GetVector2DList(input_nodes);

	return rv;
}

int InputOutputElement::GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables)
{
	// call base class routines
	Constraint::GetAvailableSpecialValues(available_variables);

	// Automatic entries for this class 
	InputOutputElement::GetAvailableSpecialValuesAuto(available_variables);

	// Manual entries for this class
	for(int i=1; i<=GetNOutputs();i++)
	{
		available_variables.Add(ReadWriteElementDataVariableType("IOBlock.output",i,0,0.,mystr("IOBlock.output[i] ... measures the i-th output of this IOBlock"))) ;
	}
	for(int i=1; i<=GetNInputs();i++)
	{
		available_variables.Add(ReadWriteElementDataVariableType("IOBlock.input",i,0,0.,mystr("IOBlock.input[i] ... measures the i-th input of this IOBlock"))) ;
	}
	return 0;
}

int InputOutputElement::ReadSingleElementData(ReadWriteElementDataVariableType& RWdata) 		
{
	// call base class routine 	
	int rv = Constraint::ReadSingleElementData(RWdata);
	if (rv == 1) return 1;

	// manual things to read  
	if( RWdata.variable_name == mystr("IOBlock.output") )
	{
		if (RWdata.comp1 > 0 && RWdata.comp1 <= GetNOutputs()) //range check
		{
			RWdata.value = GetOutput(RWdata.time, RWdata.comp1); 
			return 1; 
		}
		else return -2; 
	}
	else if(RWdata.variable_name == mystr("IOBlock.input"))
	{
		if (RWdata.comp1 > 0 && RWdata.comp1 <= GetNInputs()) //range check
		{
			RWdata.value = GetInput(RWdata.time, RWdata.comp1); 
			return 1; 
		}
		else return -2; 
	}

	return ReadSingleElementDataAuto(RWdata);
}

const char* InputOutputElement::SymbolText() const 
{
	return CSText_InputOutputElement;
}

int InputOutputElement::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = Constraint::CheckConsistency(errorstr);
	if(rv==2){	return rv;}
	
	if ((GetNInputs() != input_types.Length()) || (inputs.Length() != input_localnum.Length())) 
	{
		errorstr = "Inconsistent system: lengths of vectors for definition of input don't match in element " + elementname + ". ";
		return 2;
	}
	if(GetExpectedNumberOfInputs()!= -1 && GetNInputs() != GetExpectedNumberOfInputs())
	{
		errorstr = "Inconsistent system: please connect " + mystr(GetExpectedNumberOfInputs()) + " inputs to element " + elementname + ". ";
		return 2;
	}

	for(int i=1; i<=input_types.Length();i++)
	{
		if(input_types(i)!=IOInputTypeSensor && input_types(i)!=IOInputTypeElement)
		{
			errorstr = "Inconsistent system: unexpected type of " + mystr(i) + "-th input. ";
			return 2;			
		}
		else if(input_types(i)==IOInputTypeSensor)
		{			
			if(inputs(i) < 1 /*||inputs(i) > NSensors()*/)
			{
				errorstr = "Inconsistent system: invalid sensor number found in element " + elementname + ". ";
				return 2;				
			}
		}
		else if(input_types(i)==IOInputTypeElement)
		{
			if(inputs(i)>mbs->GetNElements())
			{
				errorstr = "Inconsistent system: input element number not well-defined in element " + elementname + ". ";
				return 2;
			}
			else if(!mbs->GetElement(inputs(i)).IsType(TController))
			{
				errorstr = "Inconsistent system: can't connect the element with number " + mystr(inputs(i)) + " to element " + elementname + "(no output found). ";
				return 2;	
			}
			else if(input_localnum(i) < 1 || mbs->GetElement(inputs(i)).IsType(TController) && ((InputOutputElement*)mbs->GetElementPtr(inputs(i)))->GetNOutputs() < input_localnum(i))
			{
				errorstr = "Inconsistent system: can't connect the element with number " + mystr(inputs(i)) + " to element " + elementname + "(unexpected output number " + mystr(input_localnum(i)) + "). ";
				return 2;				
			}
		}
	}
	return rv;
}

Vector2D InputOutputElement::GetInputPosD(int i) const //return absolute position of input #i
{
	double phi = rotation*MY_PI*0.5;
	Matrix3D rot; rot.Set22(cos(phi), -sin(phi),sin(phi), cos(phi));

	double y = 0;
	if (GetNInputs() > 1)
	{
	//	y = draw_dim.Y() - (draw_dim.Y()/((double)GetNInputs()+0.5)) * (double)i;
			y = draw_dim.Y()/2. - (draw_dim.Y()/((double)GetNInputs())) * ((double)i - 0.5);

	}
	double x =draw_dim.X()*0.5;
	Vector2D dpos(-x,y);
	if(rotation == 1)//rotation, 1==90°, 2==180°, 3==270°, 4=360°)
	{
		dpos = Vector2D(-y,-x);
	}
	else if(rotation == 2)
	{
		dpos = Vector2D(x,y);
	}
	else if(rotation == 3)
	{
		dpos = Vector2D(-y,x);
	}
	return GetRefPos2DD() + dpos; //rot * Vector2D(-draw_dim.X()*0.5,y);
} 

Vector2D InputOutputElement::GetOutputPosD(int i) const //return absolute position of input #i
{
	double phi = rotation*MY_PI*0.5;
	Matrix3D rot; rot.Set22(cos(phi), -sin(phi),sin(phi), cos(phi));

	double y = 0;
	if (GetNOutputs() > 1)
	{
		y = draw_dim.Y()/2. - (draw_dim.Y()/((double)GetNOutputs())) * ((double)i - 0.5);
	}
	return GetRefPos2DD() + rot * Vector2D(draw_dim.X()*0.5,y);
}

void InputOutputElement::DrawElement() 
{
	//DrawElement2D();
	if (!GetMBS()->GetIOption(144)) return; //do not draw control object

	//draw text in center of symbol:
	double phi = rotation*MY_PI*0.5;
	Matrix3D rot = RotMatrix3(phi);
	//Matrix rot(cos(phi), -sin(phi), sin(phi), cos(phi));
	
	//draw rectangle:
	double b = draw_dim.X();
	double h = draw_dim.Y();
	Vector2D rp = GetRefPos2DD();	
	Vector3D rp3 = ToP3D(rp);

	Vector3D p1(-0.5*b,-0.5*h,0.); //ll 2
	Vector3D p2( 0.5*b,-0.5*h,0.); //lr 3
	Vector3D p3( 0.5*b, 0.5*h,0.); //ur 4
	Vector3D p4(-0.5*b, 0.5*h,0.); //ul 1
	
	if(colbackground(1) < 0)
	{
		GetMBS()->MyDrawRectangle(rp3+rot*p4, rp3+rot*p1, rp3+rot*p2, rp3+rot*p3, 1, &colforeground);	//$ RL 2011-11-17: draw only lines
	}
	else
	{
		GetMBS()->MyDrawRectangle(rp3+rot*p4, rp3+rot*p1, rp3+rot*p2, rp3+rot*p3, 1, &colforeground, &colbackground);	//$ RL 2011-11-17: draw color rectangle
	}
	Vector2D tp = rp - Vector2D(draw_dim.X()*0.3,-draw_dim.Y()*0.3);
	GetMBS()->GetRC()->PrintText3D((float)tp.X(), (float)tp.Y(), 0., GetElementName().c_str());
	
	if (GetMBS()->GetIOption(123))
	{
		char text[40];
		sprintf(text,"element %d", GetOwnNum());

		GetMBS()->GetRC()->PrintText3D((float)tp.X(),(float)(tp.Y()-draw_dim.Y()*0.3),0.,text);
	}

	//draw inputs and outputs:
	Vector3D v1=rot*Vector3D(-0.1*b, 0.1*h, 0.);
	Vector3D v2=rot*Vector3D(-0.1*b,-0.1*h, 0.);

	if(GetNInputs() == input_types.Length())
	{
		for (int i=1; i <= GetNInputs(); i++)
		{
			//GetMBS()->MyDrawCircleXY(ToP3D(GetInputPosD(i)), draw_dim.Y()*0.1, color);
			Vector3D pi(0.);
			if(Dim() == 2)
			{
				pi = ToP3D(GetInputPosD(i));
			}
			else
			{
				pi = GetInputPos3DD(i);
			}

			GetMBS()->MyDrawLine(pi + v1, pi, 2., colforeground);
			GetMBS()->MyDrawLine(pi + v2, pi, 2., colforeground);

			Vector3D po(0.,0.,0.);
			int show = 0;
			if (GetInputType(i) == 1) //if InputOutputElement, draw line to last IO-Element
			{
				Element *elp = GetMBS()->GetElementPtr(inputs(i));
				if(! elp->IsType(TController))	//$ DR 2013-03-28 added this check
				{
					GetMBS()->UO(UO_LVL_warn) << "WARNING: input " << i << " of element " << elp->GetOwnNum() << " should be an other IO-Element, but is not. Drawing of controll elements stopped!\n";
					return;
				}
				InputOutputElement* ioe = (InputOutputElement*)GetMBS()->GetElementPtr(inputs(i));	
				if (ioe->Dim()==2)
				{
					po = ToP3D(ioe->GetOutputPosD(GetInputLocalNum(i)));
				}
				else
				{
					po = ioe->GetOutputPos3DD(GetInputLocalNum(i));
				}
				show = 1;
			}
			else if (GetMBS()->GetIOption(127) && GetInputType(i) == 2)
			{
				// sensor als eingang
				if(GetMBS()->GetSensor(inputs(i)).GetNumberOfDrawingPositions() > 0)
				{
					po = GetMBS()->GetSensor(inputs(i)).GetDrawPosition(1);
					show = 1;
				}
			}
			if (show)
			{ // zwischenpunkte
				Vector3D last = pi;
				for (int j=1; j <= input_nodes.Length(); j++)
				{
					if (input_nodes_num(j) == i)
					{
						GetMBS()->MyDrawLine(ToP3D(input_nodes(j)), last, 1., colforeground);
						last = ToP3D(input_nodes(j));
					}
				}
				GetMBS()->MyDrawLine(last, po, 1., colforeground);
			}
		}
	}

	for (int i=1; i <= GetNOutputs(); i++)
	{
		Vector2D v2d = GetOutputPosD(i);
		Vector3D po = ToP3D(GetOutputPosD(i));

		GetMBS()->MyDrawLine(po + Vector3D(v1(1),2.0/GetNOutputs()*v1(2),v1(3)), po, 1., colforeground);
		GetMBS()->MyDrawLine(po + Vector3D(v2(1),2.0/GetNOutputs()*v2(2),v2(3)), po, 1., colforeground);
	}
};



void InputOutputElement::SetOutputName(int output_nr, mystr& name) 
{
	if(output_names.Length() != GetNOutputs())      // conditional reset ( is called when a name is added for the first time )
	{
		output_names.Reset();
		for (int i=1; i<=GetNOutputs(); i++)
		{
			output_names.Add("");
		}
	}
	if (name.Length()>=80)
	{
		name.Left(77);
		name += "...";
	}
	output_names.Set(output_nr, name);
}

char* InputOutputElement::GetOutputName(int output_nr)
{
	if( output_names.Length() > 0 && output_nr <= output_names.Length() ) 
		return output_names.Get(output_nr);
	else
		return NULL;
}

void InputOutputElement::DrawElement2D()
{
// ASSUMPTION: the default element size is 24x24 (const double DEFAULT_IOE_SIZE = 24;)
	const int textboxheight = 8;
	const int nodediameter = 1;

	
// Basic Size, rotation
	Vector3D refpos = GetRefPosD();
	double phi = rotation*MY_PI*0.5;
	Matrix2D rot;	rot.SetSkew(phi);

	Vector2D center, size;            // used for frame, labels 
	Vector2D nodecenter, nodesize;    // individual nodes / lines


// Drawing order for nice pictures
// 1: connection lines
// 2: Element Frame
// 3: Nodes
// 4: Text

// ****************
// (1) CONNECTIONS
// ****************
// loop over Input Nodes
	center = Vector2D(refpos.X(),refpos.Y()); 
	size = Vector2D(draw_dim.X(), draw_dim.Y());

	int last_array_pos_to_insert = 1;      // array pos to insert the node at

	for (int i=1; i <= GetNInputs(); i++)
	{
		nodecenter = Vector2D(-0.5*draw_dim.X(), 0.5*draw_dim.Y() - (draw_dim.Y()/((double)GetNInputs())) * ((double)i - 0.5)); // relative position
		nodecenter = rot * nodecenter + center;
		nodesize = Vector2D(nodediameter,nodediameter);         
	//	nodesize = Vector2D(nodediameter/2.,nodediameter/2.);             // construction nodes smaller then IO Nodes

		if (GetInputType(i) == 1) // connected to an other IO-Element
		{
			//InputOutputElement* ioe = (InputOutputElement*)GetMBS()->GetElementPtr(inputs(i));	
			//$ DR 2013-03-28 added the following check
			Element *elp = GetMBS()->GetElementPtr(inputs(i));
			if(! elp->IsType(TController))
			{
				GetMBS()->UO(UO_LVL_warn) << "WARNING: input " << i << " of element " << elp->GetOwnNum() << " should be an other IO-Element, but is not.";
			}
			else
			{
				InputOutputElement* ioe = (InputOutputElement*)elp;

				Vector2D final = ioe->GetOutputPosD(GetInputLocalNum(i));
				Vector2D start = nodecenter;
				Vector2D end;

				for (int j=1; j <= input_nodes.Length(); j++)
				{
					if (input_nodes_num(j) == i)
					{
						end = input_nodes(j);

						last_array_pos_to_insert = j;                                                   
						int identifier_subelnr = i * 65536 + last_array_pos_to_insert; // combines Input Number and Array Position to insert a split node

// AD: NOTE - I dont want to change the data struct and access functions right now, but to identify the connection line unambiguous there are two integer numbers required. 
//            As these integers are short I store them both in the variable sub_elnr.

						GetMBS()->AddDrawComponent_Line( GetOwnNum(), TIOBlock,								// MBS identifier
							identifier_subelnr, TConnectionLine,																// DRAW identifier
							start, end,																													// extent
							Vector3D(0.,0.,0.));																								// color

						GetMBS()->AddDrawComponent_Ellipse( GetOwnNum(), TIOBlock,						// MBS identifier
							identifier_subelnr, TConstructionNode,															// DRAW identifier
							end, nodesize,																											// extent
							Vector3D(0.,0.,0.), Vector3D(0.,0.,0.));														// color
						start = end;
						
						last_array_pos_to_insert++;																						// after a node is found increase the array position to insert new nodes
					}
				}

				int identifier_subelnr = i * 65536 + last_array_pos_to_insert;						// combines Input Number and Array Position to insert a split node
				GetMBS()->AddDrawComponent_Line( GetOwnNum(), TIOBlock,										// MBS identifier
					identifier_subelnr, TConnectionLine,																		// DRAW identifier
					start, final,																														// extent
					Vector3D(0.,0.,0.));																										// color
			}
		}
	}

// ******************
// (2) ELEMENT FRAME
// ******************
	center = Vector2D(refpos.X(),refpos.Y()); 
	size = Vector2D(draw_dim.X(), draw_dim.Y());
	GetMBS()->AddDrawComponent_Rect( GetOwnNum(), TIOBlock,												// MBS identifier
																	 1, TElementFrame,														// DRAW identifier
																	 center, size,																// extent
																	 Vector3D(0.,0.,0.), colbackground);          // colors : border always black, area as set ...
	// specific Symbol
	if(size.X() > 2)
	DrawBlockSymbol();

// ***************************
// (3) INPUT AND OUTPUT NODES
// ***************************
	center = Vector2D(refpos.X(),refpos.Y());
	size = Vector2D(draw_dim.X(), draw_dim.Y());
	nodesize = Vector2D(nodediameter,nodediameter);
	// Input Nodes
	for (int i=1; i <= GetNInputs(); i++)
	{
		nodecenter = Vector2D(-0.5*draw_dim.X(), 0.5*draw_dim.Y() - (draw_dim.Y()/((double)GetNInputs())) * ((double)i - 0.5)); // relative position
		nodecenter = rot * nodecenter + center;
		
		if (GetInputType(i) == 1) // connected to an other IO-Element
		{
			GetMBS()->AddDrawComponent_Ellipse( GetOwnNum(), TIOBlock,								// MBS identifier
																					i, TInputNode,												// DRAW identifier
																					nodecenter, nodesize,									// extent
																					Vector3D(0.,0.,0.), Vector3D(0.,0.,0.));		// color
		}
		else if (GetInputType(i) == 2) // connected to a sensor
		{
			GetMBS()->AddDrawComponent_Ellipse( GetOwnNum(), TIOBlock,								// MBS identifier
																					i, TInputNode,												// DRAW identifier
																					nodecenter, nodesize,									// extent
																					Vector3D(0.,0.,0.), Vector3D(0.,0.,0.));		// color
		}
		else 
		{
			GetMBS()->AddDrawComponent_Ellipse( GetOwnNum(), TIOBlock,								// MBS identifier
																					i, TInputNode,												// DRAW identifier
																					nodecenter, nodesize,									// extent
																					Vector3D(0.,0.,0.), Vector3D(1.,0.,0.));		// color
		}		
	}
	// Output Nodes
	for (int i=1; i <= GetNOutputs(); i++)
	{
		nodecenter = Vector2D( 0.5*draw_dim.X(), 0.5*draw_dim.Y() - (draw_dim.Y()/((double)GetNOutputs())) * ((double)i - 0.5)); // relative position
		nodecenter = rot * nodecenter + center;

		GetMBS()->AddDrawComponent_Ellipse( GetOwnNum(), TIOBlock,									// MBS identifier
																				i, TOutputNode,													// DRAW identifier
																				nodecenter, nodesize,										// extent
																				Vector3D(0.,0.,0.), Vector3D(0.,0.,0.));			// color
	}

// ***********
// (4) LABELS
// ***********
	// Element Name below the element // "about a 3rd to 4th of the default height"
	center = Vector2D(refpos.X(),refpos.Y()-0.5*draw_dim.Y() - textboxheight/2);	
	size = Vector2D(draw_dim.X(), textboxheight);                                 
	GetMBS()->AddDrawComponent_Text( GetOwnNum(), TIOBlock,												// MBS identifier
																	 1, TElementName,															// DRAW identifier
																	 center, size,																// extent
															     Vector3D(0.,0.,0.),													// color of this label is alway black
																	 GetElementName(),                            // label element with element name
																	 TTextAllign (HCenter+VTop) );								// other properties
	// MBS-element number 
	//b	123	PostProcOptions.Bodies.show_element_numbers	0	"1|(0) ... (Don't) show element body numbers."
	//b	124	PostProcOptions.Connectors.show_constraint_numbers	0	"1|(0) ... (Don't) show constraint number."
	if (GetMBS()->GetIOption(124))
	{
		size = Vector2D( 6, 5); // rectangle of constant absolute size in top right corner ( "for about 4 digits" )
		center = Vector2D( refpos.X()+0.5*draw_dim.X()-0.5*size.X(), refpos.Y()+0.5*draw_dim.Y()-0.5*size.Y());
		if(draw_dim.X() > 5)
		{
		GetMBS()->AddDrawComponent_Rect( GetOwnNum(), TIOBlock,												// MBS identifier
																		 2, TElementName,															// DRAW identifier
																		 center, size,																// extent
																		 colforeground, colbackground);								// colors : border always black, area as set ...
		GetMBS()->AddDrawComponent_Text( GetOwnNum(), TIOBlock,												// MBS identifier
																		 2, TElementName,															// DRAW identifier
																		 center, size,																// extent
																		 colforeground,																// color of this label is alway black
																		 mystr(GetOwnNum()),	  											// label element with element name
																		 TTextAllign (HCenter+VTop) );								// other properties
		}
	}

	center = Vector2D(refpos.X(),refpos.Y()); 
	size = Vector2D(draw_dim.X(), draw_dim.Y());
	// Input Nodes
	for (int i=1; i <= GetNInputs(); i++)
	{
		nodecenter = Vector2D(-0.5*draw_dim.X(), 0.5*draw_dim.Y() - (draw_dim.Y()/((double)GetNInputs())) * ((double)i - 0.5)); // relative position
		nodecenter = rot * nodecenter + center;
		Vector2D shift(-DEFAULT_IOE_SIZE/2-nodediameter , textboxheight/2 + nodediameter);
		Vector2D textbox(DEFAULT_IOE_SIZE, textboxheight);

		if (GetInputType(i) == 2) // connected to a sensor
		{
			if (inputs(i) > 0 && inputs(i) <= GetMBS()->NSensors())
			{
				mystr name = GetMBS()->GetSensor(inputs(i)).GetSensorName();
				mystr label = mystr("S")+mystr(inputs(i))+mystr(":")+name;
				GetMBS()->AddDrawComponent_Text( GetOwnNum(), TIOBlock,									// MBS identifier
																				 i, TInputNode,													// DRAW identifier
																				 nodecenter+shift, textbox,							// extent
																				 Vector3D(0.,0.,0.),										// color of this label is alway black
																				 label,																	// label input node with sensor name
			  																 TTextAllign (HRight+VBottom) );				// other properties 
			}
		}
	}
	// Output Nodes
	for (int i=1; i <= GetNOutputs(); i++)
	{
		nodecenter = Vector2D( 0.5*draw_dim.X(), 0.5*draw_dim.Y() - (draw_dim.Y()/((double)GetNOutputs())) * ((double)i - 0.5)); // relative position
		nodecenter = rot * nodecenter + center;
		Vector2D shift(+DEFAULT_IOE_SIZE/2+nodediameter , textboxheight/2 + nodediameter);
		Vector2D textbox(DEFAULT_IOE_SIZE, textboxheight);
		if (GetOutputName(i) != NULL)
		{
		GetMBS()->AddDrawComponent_Text( GetOwnNum(), TIOBlock,											// MBS identifier
																		 i, TOutputNode,														// DRAW identifier
																		 nodecenter+shift, textbox,									// extent
																		 Vector3D(0.,0.,0.),												// color of this label is alway black
																		 mystr(GetOutputName(i)),										// label output with load identifier
			  														 TTextAllign (HLeft+VBottom) );							// other properties 
		}
	}

	return;
};

void InputOutputElement::DrawBlockSymbol()
{
	Vector3D refpos = GetRefPosD();
	Vector2D center(refpos.X(),refpos.Y());
	Vector2D size(draw_dim.X(), draw_dim.Y());

	GetMBS()->AddDrawComponent_Text( GetOwnNum(), TIOBlock,												// MBS identifier
																	 1, TSymbol,																	// DRAW identifier
																	 center, size,																// extent
																	 colforeground,																// color
																	 mystr(SymbolText()),	  											// label element with element name
																	 TTextAllign (HCenter+VCenter) );							// other properties
}

void InputOutputElement::GetDirectFeedThroughElements(TArray<int>& elnums) const //return all elements which are depending on this element by direct feed through
{
	if (IsDirectFeedThrough())
	{
		for (int i=1; i <= GetNInputs(); i++)
		{
			if (GetInputType(i) == IOInputTypeElement)
			{
				elnums.Add(GetInputNum(i));
				GetMBS()->GetElement(GetInputNum(i)).GetDirectFeedThroughElements(elnums);
			}
		}
	}
} 


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
const char* LinearTransformation::SymbolText() const 
{
	return CSText_LinearTransformation;
}

void LinearTransformation::SetLinearTransformation(const Matrix& A, const Vector& b)
{
	//inputs must be added separately
	A_coeff.Resize(A.Getrows(), A.Getcols(),1);
	A_coeff = A;
	b_coeff = b;
	SetNOutputs(b.Length());

	//check if sizes are consistent!
	//if (A.Getrows()*A.Getcols() != 0) assert(A.Getrows() == b.Length()); //-->moved to CheckConsistency
}

void LinearTransformation::SetConstant(double val)
{
	//inputs must be added separately
	A_coeff = Matrix(1,0);		b_coeff = Vector(val);
	SetNOutputs(1);
}

void LinearTransformation::SetGain(double val)
{
	//inputs must be added separately
	A_coeff = Matrix(val);
	b_coeff = Vector(1);
	b_coeff(1) = 0;
	SetNOutputs(1);
}

void LinearTransformation::SetAdder(const Vector& signs)
{
	//inputs must be added separately
	A_coeff = Matrix(1,signs.Length());
	for (int i=1; i <= signs.Length(); i++) {A_coeff(1,i) = signs(i);}
	b_coeff = Vector(1);
	b_coeff(1) = 0.;
	SetNOutputs(1);
}


//To be overwritten in derived class:
void LinearTransformation::CopyFrom(const Element& e)
{
	InputOutputElement::CopyFrom(e);
	const LinearTransformation& ce = (const LinearTransformation&)e;

	A_coeff.Resize(ce.A_coeff.Getrows(), ce.A_coeff.Getcols(),1);
	A_coeff = ce.A_coeff;
	b_coeff = ce.b_coeff;
}

void LinearTransformation::InitConstructor()
{
	InputOutputElement::InitConstructor();
	elementname = GetElementSpec();
	SetNStates(0);
	Matrix A_init(4,4);
	SetLinearTransformation(A_init,Vector(4));
}

int LinearTransformation::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = InputOutputElement::CheckConsistency(errorstr);
	if(rv==2){return rv;}
	if(A_coeff.Getrows()*A_coeff.Getcols() != 0)
	{
		if(A_coeff.Getcols() != GetNInputs()) 
		{
			errorstr = "Inconsistent system: columns of A_matrix and number of inputs do not match " + mystr("(check element with number ") + mystr(GetOwnNum()) + mystr(")");
			return 2;
		}
		else if(A_coeff.Getrows() != b_coeff.Length())
		{	
			errorstr = "Inconsistent system: rows of A_matrix and b_vector do not match " + mystr("(check element with number ") + mystr(GetOwnNum()) + mystr(")");
			return 2;
		}
		else if(A_coeff.Getrows() != GetNOutputs())
		{	
			errorstr = "Inconsistent system: rows of A_matrix and number of output do not match " + mystr("(check element with number ") + mystr(GetOwnNum()) + mystr(")");
			return 2;		
		}
	}
	else if(A_coeff.Getrows()*A_coeff.Getcols() == 0 && GetNInputs() != 0)
	{
		errorstr = "Inconsistent system: no rows or columns of A_matrix defined, but number of inputs is not zero " + mystr("(check element with number ") + mystr(GetOwnNum()) + mystr(")");
		return 2;		
	}
	return rv;
}

void LinearTransformation::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	InputOutputElement::GetElementData(edc);

	ElementData ed;
}

int LinearTransformation::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	//TODO: resize matrix, if user change the size!!!
	int rv = InputOutputElement::SetElementData(edc);

	//// possible extensions e.g. with radio buttons
	//SetConstant(val);
	//SetGain(val);
	//SetAdder(signs);


	//if(A_coeff.Getrows() == 1 && A_coeff.Getcols() == 1)
	//{
	//	if( A_coeff(1,1) == 0.0)
	//	{
	//		SetConstant(b_coeff(1));
	//	}
	//	else if(b_coeff.Length()==1 && b_coeff(1) == 0.0)
	//	{
	//		SetGain(A_coeff(1,1));
	//	}
	//}
	//else
	//{
		SetLinearTransformation(A_coeff, b_coeff);
	//}
	return rv;
}

double LinearTransformation::GetOutput(double t, int i) const 
{
	//to be overwritten in specific class!
	int elnum = GetOwnNum();
	double yi = b_coeff(i);
	for (int j=1; j <= GetNInputs(); j++)
	{
		//y_i = sum_j=1..nI {A_ij*u_j}  + b_i
		yi += A_coeff(i, j)*GetInput(t, j);
	}

	/*if (GetOwnNum() == 25 && GetMBS()->GetTime()>0.005)
	{
		GetMBS()->UO() << "Gain: t=" << GetMBS()->GetTime() << "\n";
		GetMBS()->UO() << "u=" << GetInput(t, 1) << "\n";
		GetMBS()->UO() << "y=" << yi << "\n";
	}*/

	return yi;
}

void LinearTransformation::DrawBlockSymbol() 
{ 
	return InputOutputElement::DrawBlockSymbol();

	Vector3D refpos = GetRefPosD();
	Vector2D center(refpos.X(),refpos.Y());
	Vector2D size(draw_dim.X(), draw_dim.Y());

	GetMBS()->AddDrawComponent_Line( GetOwnNum(), TIOBlock,												// MBS identifier
			  													 1, TSymbol,																	// DRAW identifier
																	 center + Vector2D( -0.35*size.X(), -0.4*size.Y() ), center + Vector2D( -0.35*size.X(), 0.4*size.Y() ),	// extent
																	 Vector3D(0.,0.,0.) );												// other properties

	GetMBS()->AddDrawComponent_Line( GetOwnNum(), TIOBlock,												// MBS identifier
																	 1, TSymbol,																	// DRAW identifier
																	 center + Vector2D( -0.4*size.X(), -0.35*size.Y() ), center + Vector2D( 0.4*size.X(), -0.35*size.Y() ),	// extent
																	 Vector3D(0.,0.,0.) );												// other properties

	GetMBS()->AddDrawComponent_Line( GetOwnNum(), TIOBlock,												// MBS identifier
																	 1, TSymbol,																	// DRAW identifier
																	 center + Vector2D( -0.4*size.X(), -0.2*size.Y() ), center + Vector2D( 0.4*size.X(), 0.2*size.Y() ),	// extent
																	 Vector3D(1.,0.,0.) );												// other properties
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Quantizer::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	InputOutputElement::GetElementData(edc);

	ElementData ed;
	ed.SetDouble(roundval, "rounding_value"); ed.SetToolTipText("Max. amplitude of random value."); edc.TreeAdd("IOBlock",ed);
	//edc.Add(ed);
}

int Quantizer::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = InputOutputElement::SetElementData(edc);

	GetElemDataDouble(mbs, edc, "IOBlock.rounding_value", roundval, 1);

	SetNOutputs(1);
	SetNStates(0);

	return rv;
}

void Quantizer::DrawBlockSymbol() 
{ 
	InputOutputElement::DrawBlockSymbol();
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
const char* STransferFunction::SymbolText() const 
{
	return CSText_STransferFunction;
}

int STransferFunction::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = InputOutputElement::CheckConsistency(errorstr);
	if(rv){	return rv;}

	if(num.Length() != den.Length())
	{
		errorstr = "Error: Numerator and denominator vectors must have same size. Transfer function is set automatically to \"1/s\"."; 
		SetSTransferFunction(Vector(1.,0.), Vector(0.,1.)); // set to valid value
		return 1; 
	}

	if(den(den.Length()) == 0)
	{
		errorstr = "Error: Denominator coefficient corresponding to Laplace variable with highest exponent is zero. Transfer function is set automatically to \"1/s\".\n"; 
		SetSTransferFunction(Vector(1.,0.), Vector(0.,1.)); // set to valid value
		return 1;
	}

	if (GetNStates() == 0)
	{
		errorstr = "Error: Transfer function is fully algebraic. Transfer function is set automatically to \"1/s\".";
		SetSTransferFunction(Vector(1.,0.), Vector(0.,1.)); // set to valid value
		return 1;
	}

	//check if sizes are consistent!
	if(GetNStates() != x_init.Length())
	{
		errorstr = "Error: Wrong size of intitialization vector. Transfer function is set automatically to \"1/s\".";
		SetSTransferFunction(Vector(1.,0.), Vector(0.,1.)); // set to valid value
		return 1;
	}
	
	return rv;
}

void STransferFunction::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	InputOutputElement::GetElementData(edc);

	ElementData ed;
	//SetElemDataVector(edc, x_init, "Initial_vector"); edc.Last().SetVariableLength(); edc.Last().SetToolTipText("Initial values of time-domain variables");
}

int STransferFunction::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = InputOutputElement::SetElementData(edc);
	SetSTransferFunction(num, den);
	// alternative with inital state "initVector"
	//GetElemDataVector(mbs, edc, "Initial_vector", initVector, 1);
	//SetSTransferFunction(numI, denI, initVector);
	return rv;
}

void STransferFunction::EvalF(Vector& f, double t) 
{
	//compute f(i) = sum_j (A_ij * XG(j)) + b_i
	/*if (GetOwnNum() == 26 && GetMBS()->GetTime()>0.005)
	{
		int test=1;
	}*/
	if(GetNInputs() >=1) //$ RL 2011-01
	{
		int n = ES();
		double u = GetInput(t, 1);

		//f = A*x+u*b;
		//A=[0 ..... 0 -den1    ]
		//  [1 0 ... 0 -den2    ]
		//  [0 0 ... 1 -den(n)]
		//bb(i) = num(i) - num(n+1)*den(i)

		f(1) += -den(1)*XG(n) + (num(1)-num(n+1)*den(1))*u;

		/*if (GetOwnNum() == 26 && GetMBS()->GetTime()>0.005)
		{
			GetMBS()->UO() << "integrator: t=" << GetMBS()->GetTime() << "\n";
			GetMBS()->UO() << "u=" << u << "\n";
			GetMBS()->UO() << "xg=" << XG(1) << "\n";
			GetMBS()->UO() << "f=" << f(1) << "\n";
		}*/

		for(int i = 2; i<=n; i++ )
		{
			f(i) += XG(i-1) - den(i)*XG(n) + (num(i)-num(n+1)*den(i))*u; //A(i,i-1) * x(i-1) + b(i)*u;  // ones in A-matrix
		}
	}
};

double STransferFunction::GetOutput(double t, int i) const 
{
	//to be overwritten in specific class!
	if (i != 1) 
	{
		assert(0);
		return 0;
	}

	//y = [0 0 .... 1][x_1 x_2 ... x_n] + b_n*u;
	int n = ES();
	double u = 0;
	if (num(n+1) != 0.) u = GetInput(t, 1);

	if (n==0) return num(n+1)*u;
	else return XG(n) + num(n+1)*u;
}

void STransferFunction::DrawBlockSymbol() 
{ 
	InputOutputElement::DrawBlockSymbol();
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
const char* LinearODE::SymbolText() const 
{
	return CSText_LinearODE;
}

int LinearODE::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = InputOutputElement::CheckConsistency(errorstr);
	if(rv){	return rv;}

	return rv;
}

void LinearODE::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	InputOutputElement::GetElementData(edc);

	ElementData ed;
	ed.SetMatrix(A_coeff.GetMatPtr(),A_coeff.Getrows(),A_coeff.Getcols(),"A_coeffs"); ed.SetVariableLength(); ed.SetToolTipText("Coefficients of state matrix A, x_dot = A*x + B*u"); edc.TreeAdd("IOBlock",ed); // change: MSax 28-02-2013
	ed.SetMatrix(B_coeff.GetMatPtr(),B_coeff.Getrows(),B_coeff.Getcols(),"B_coeffs"); ed.SetVariableLength(); ed.SetToolTipText("Coefficients of input matrix B, x_dot = A*x + B*u"); edc.TreeAdd("IOBlock",ed); // change: MSax 28-02-2013
	ed.SetMatrix(C_coeff.GetMatPtr(),C_coeff.Getrows(),C_coeff.Getcols(),"C_coeffs"); ed.SetVariableLength(); ed.SetToolTipText("Coefficients of output matrix C, y = C*x + D*u"); edc.TreeAdd("IOBlock",ed); // change: MSax 28-02-2013
	ed.SetMatrix(D_coeff.GetMatPtr(),D_coeff.Getrows(),D_coeff.Getcols(),"D_coeffs"); ed.SetVariableLength(); ed.SetToolTipText("Coefficients of output matrix D, y = C*x + D*u"); edc.TreeAdd("IOBlock",ed); // change: MSax 28-02-2013

	ed.SetVector(x_init.GetVecPtr(),x_init.GetLen(),"initital_vector"); ed.SetVariableLength(); ed.SetToolTipText("Initial values of time-domain variables"); edc.TreeAdd("IOBlock",ed); // change: MSax 28-02-2013

	//ElementData ed;

	//SetElemDataMatrix(edc, A_coeff, "A_coeffs"); edc.Last().SetVariableLength(); edc.Last().SetToolTipText("Coefficients of state matrix A, x_dot = A*x + B*u");
	//SetElemDataMatrix(edc, B_coeff, "B_coeffs"); edc.Last().SetVariableLength(); edc.Last().SetToolTipText("Coefficients of input matrix B, x_dot = A*x + B*u");
	//SetElemDataMatrix(edc, C_coeff, "C_coeffs"); edc.Last().SetVariableLength(); edc.Last().SetToolTipText("Coefficients of output matrix C, y = C*x + D*u");
	//SetElemDataMatrix(edc, D_coeff, "D_coeffs"); edc.Last().SetVariableLength(); edc.Last().SetToolTipText("Coefficients of output matrix D, y = C*x + D*u");

	//SetElemDataVector(edc, x_init, "initital_vector"); edc.Last().SetVariableLength(); edc.Last().SetToolTipText("Initial values of time-domain variables");
}

int LinearODE::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = InputOutputElement::SetElementData(edc);

	Vector numI, denI, initVector;
	GetElemDataMatrix(mbs, edc, "IOBlock.A_coeffs", A_coeff, 1); // change: MSax 28-02-2013
	GetElemDataMatrix(mbs, edc, "IOBlock.B_coeffs", B_coeff, 1); // change: MSax 28-02-2013
	GetElemDataMatrix(mbs, edc, "IOBlock.C_coeffs", C_coeff, 1); // change: MSax 28-02-2013
	GetElemDataMatrix(mbs, edc, "IOBlock.D_coeffs", D_coeff, 1); // change: MSax 28-02-2013

	GetElemDataVector(mbs, edc, "IOBlock.initital_vector", initVector, 1); // change: MSax 28-02-2013

	int ni = B_coeff.Getcols(); //number of inputs
	int ns = A_coeff.Getcols(); //number of states
	int no = C_coeff.Getrows(); //number of outputs

	if (!(A_coeff.Getrows() == ns) || !(B_coeff.Getrows() == ns) ||
		!(C_coeff.Getcols() == ns) || !(D_coeff.Getrows() == no) ||
		!(D_coeff.Getcols() == ni) || !(initVector.Length() == ns)) rv = 0;


	if (rv) SetLinearODE(A_coeff, B_coeff, C_coeff, D_coeff, initVector);

	return rv;
}

void LinearODE::DrawBlockSymbol() 
{ 
	InputOutputElement::DrawBlockSymbol();
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
const char* IOMathFunction::SymbolText() const 
{
	return mathfunc.GetTypeName();
}

void IOMathFunction::SetIOMathFunction(const MathFunction& mathfuncI)
{
	mathfunc = mathfuncI;

	//elementspec = mystr("IOMathFunction");// + mystr(mathfunc.GetTypeName());

	//SetNOutputs(1);
	//SetNStates(0);
}

Element* IOMathFunction::GetCopy()
{
	Element* ec = new IOMathFunction(mbs);
	ec->CopyFrom(*this);

	return ec;
}
//To be overwritten in derived class:
void IOMathFunction::CopyFrom(const Element& e)
{
	InputOutputElement::CopyFrom(e);
	const IOMathFunction& ce = (const IOMathFunction&)e;

	mathfunc = ce.mathfunc;	
	pieceWiseSwitchOnlyInPostNewton = ce.pieceWiseSwitchOnlyInPostNewton;
	//pieceWiseIndex = ce.pieceWiseIndex;
	pieceWiseIteration = ce.pieceWiseIteration;
}

void IOMathFunction::InitConstructor()
{
	InputOutputElement::InitConstructor();
	elementname = GetElementSpec();

	//Vector datainit(DataS()); 
	//datainit.SetAll(-1);
	//SetDataInit(datainit);
	pieceWiseSwitchOnlyInPostNewton = 0;
	pieceWiseIndex = -1;
	pieceWiseIteration = 0;
	SetNOutputs(1);
	SetNStates(0);
}

int IOMathFunction::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = InputOutputElement::CheckConsistency(errorstr);
	if(rv){	return rv;}

	return rv;
}

void IOMathFunction::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	InputOutputElement::GetElementData(edc);
	//mathfunc.GetElementData(edc); 		//fill in all element data

	ElementDataContainer edc_mf;
	mathfunc.GetElementData(edc_mf);
	ElementData ed;
	ed.SetEDC(&edc_mf,"MathFunction"); ed.SetToolTipText("mathematical function"); edc.TreeAdd("IOBlock",ed);	
	//ElementData ed;
	//SetElemDataMathFunc(edc, mathfunc, "Math_function"); edc.Last().SetToolTipText("Set coefficients of math function");
}

int IOMathFunction::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	
	int rv = InputOutputElement::SetElementData(edc);

	//const ElementData* ed = edc.TreeFind("IOBlock.MathFunction");
	ElementData* ed = edc.TreeFind("IOBlock.MathFunction");
	if(ed && ed->IsEDC())
	{
		//mathfunc.SetElementData(mbs, *ed->GetEDC());
		mathfunc.SetElementData(mbs,*ed->GetEDC());	//$ DR 2012-12-12
		SetIOMathFunction(mathfunc);
	}
	return rv;
}

void IOMathFunction::DrawBlockSymbol() 
{ 
	return InputOutputElement::DrawBlockSymbol();
}

double IOMathFunction::GetOutput(double t, int i) const
{
	//to be overwritten in specific class!
	double u = GetInput(t, i);

	//------------------------------------
  //begin: special case: input is IOTime
	if(mathfunc.GetFuncMode() == TMFpiecewiseconst  && u == GetMBS()->GetStepEndTime() && u > 0.0)
	{
		if(GetInputType(i) == IOInputTypeElement)
		{
			mystr inputname(((const InputOutputElement&)(GetMBS()->GetElement(GetInputNum(i)))).SymbolText());
			if(inputname.Compare(mystr(CSText_IOTime)))
			{
				// use time tolerance ONLY if input element is IOTime and time is equal to solver - StepEndTime and piecewise constant Mathfunction is used!!!
				u = u - discrete_time_tol; // if solver needs value at next time point for integration; value changes AFTERWARDS integration of time step is finished (then TItime is increased by FinishStep).
			}
		}
	}
  //end: special case: input is IOTime
	//----------------------------------
	if(pieceWiseSwitchOnlyInPostNewton && pieceWiseIndex>0)
	{
		return mathfunc.InterpolatePiecewise(u, pieceWiseIndex);
	}
	else
	{
		return mathfunc.Evaluate(u);
	}
}
//$ RL 2012-7-25:[ 
double IOMathFunction::PostNewtonStep(double t)
{
	double error=0.;
	if(pieceWiseSwitchOnlyInPostNewton && mathfunc.GetFuncMode() == TMFpiecewiselinear &&  !pieceWiseIteration == 1 )
	{	
		double u = GetInput(t, 1);
		int pieceWiseIndexOld = pieceWiseIndex;
		pieceWiseIndex = mathfunc.FindIndexPiecewise(u);
		if(pieceWiseIndexOld != -1 && pieceWiseIndex != pieceWiseIndexOld)
		{
			mbs->ForceJacobianRecomputation();

			pieceWiseIteration++;
			double fnew = mathfunc.InterpolatePiecewise(u, pieceWiseIndex);
			double fold = mathfunc.InterpolatePiecewise(u, pieceWiseIndexOld);
						
			//error = fabs(fnew-fold); // difference of new and old value
			double f_nominal = mathfunc.GetXVector().MaxNorm();
			if(!f_nominal){	f_nominal = 1;} // in case of zero - vector: use 1.0 as nominal value to omit division by zero
			error = fabs((fnew-fold)/f_nominal); // difference of new and old value
		}
		else
		{
			pieceWiseIteration = 0;
		}
	}
	return error;
}

void IOMathFunction::PostprocessingStep() 
{
	pieceWiseIteration = 0;
};
//$ RL 2012-7-25:] 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void IOSaturate::GetElementData(ElementDataContainer& edc)
{
	InputOutputElement::GetElementData(edc);
	ElementData ed;
	ed.SetDouble(upperLimit, "upper_limit"); ed.SetToolTipText("Upper limit of saturate."); edc.TreeAdd("IOBlock",ed);
	ed.SetDouble(lowerLimit, "lower_limit"); ed.SetToolTipText("Lower limit of saturate."); edc.TreeAdd("IOBlock",ed);

}
int IOSaturate::SetElementData(ElementDataContainer& edc)
{
	int rv = InputOutputElement::SetElementData(edc);
  GetElemDataDouble(mbs, edc, "IOBlock.upper_limit", upperLimit, 0);
	GetElemDataDouble(mbs, edc, "IOBlock.lower_limit", lowerLimit, 0);

	SetNOutputs(1);
	SetNStates(0);
	return rv;
}


const char* IOSaturate::SymbolText() const 
{
	return CSText_IOSaturate;
}

void IOSaturate::DrawBlockSymbol() 
{ 
	InputOutputElement::DrawBlockSymbol();
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
const char* IODeadZone::SymbolText() const 
{
	return CSText_IODeadZone;
}

void IODeadZone::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	InputOutputElement::GetElementData(edc);
	ElementData ed; 
	ed.SetDouble(start_deadzone, "start_deadzone"); ed.SetToolTipText("Start of dead zone."); edc.TreeAdd("IOBlock",ed);
	ed.SetDouble(end_deadzone, "end_deadzone"); ed.SetToolTipText("End of dead zone."); edc.TreeAdd("IOBlock",ed);
}

int IODeadZone::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = InputOutputElement::SetElementData(edc);
  GetElemDataDouble(mbs, edc, "IOBlock.start_deadzone", start_deadzone, 0); // lower limit of deadzone
	GetElemDataDouble(mbs, edc, "IOBlock.end_deadzone", end_deadzone, 0);     // upper limit of deadzone
	return rv;
}

void IODeadZone::DrawBlockSymbol() 
{ 
	InputOutputElement::DrawBlockSymbol();
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
const char* IOTime::SymbolText() const 
{
	return CSText_IOTime;
}

void IOTime::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	InputOutputElement::GetElementData(edc);
	ElementData ed;
}

int IOTime::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = InputOutputElement::SetElementData(edc);

	inputs.SetLen(0);  //$ MSax 2013-08-28: added
	input_types.SetLen(0); //$ MSax 2013-08-28: added
	input_localnum.SetLen(0); //$ MSax 2013-08-28: added

	return rv;
}

void IOTime::DrawBlockSymbol()
{
	Vector3D refpos = GetRefPosD();
	Vector2D center(refpos.X(),refpos.Y());
	Vector2D size(draw_dim.X(), draw_dim.Y());

	GetMBS()->AddDrawComponent_Ellipse( GetOwnNum(), TIOBlock,										// MBS identifier
																			1, TSymbol,																// DRAW identifier
																			center, 0.7* size,												// extent
																			colforeground, colbackground );						// color
//#define animate_clock
#ifdef animate_clock
	double angle = GetMBS()->GetDrawTime()*2.*MY_PI;                              // 1 full revelation per second
	GetMBS()->AddDrawComponent_Line( GetOwnNum(), TIOBlock,												// MBS identifier
																	 1, TSymbol,																	// DRAW identifier
																	 center, center + Vector2D(0.30*size.Y()*sin(angle),0.30*size.Y()*cos(angle)),	// extent
																	 colforeground );															// color
#else
	GetMBS()->AddDrawComponent_Line( GetOwnNum(), TIOBlock,												// MBS identifier
																	 1, TSymbol,																	// DRAW identifier
																	 center, center + Vector2D(0.,0.30*size.Y()),	// extent
																	 colforeground );															// color
	
	GetMBS()->AddDrawComponent_Line( GetOwnNum(), TIOBlock,												// MBS identifier
																	 1, TSymbol,																	// DRAW identifier
																	 center, center + Vector2D(0.20*size.X(),0.),	// extent
																	 colforeground );															// color
#endif
}
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
const char* IOPulseGenerator::SymbolText() const 
{
	return CSText_IOPulseGenerator;
}

void IOPulseGenerator::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	InputOutputElement::GetElementData(edc);
	
	ElementData ed;
	ed.SetDouble(amplitude, "amplitude"); ed.SetToolTipText("Amplitude of rectangle pulse generator."); edc.TreeAdd("IOBlock",ed);
	ed.SetDouble(toffs, "offset"); ed.SetToolTipText("Time offset (s)."); edc.TreeAdd("IOBlock",ed);
	ed.SetDouble(period, "period"); ed.SetToolTipText("Period of signal (s)."); edc.TreeAdd("IOBlock",ed);
	ed.SetDouble(pulseWidth, "pulse_width"); ed.SetToolTipText("Pulse width (s)."); edc.TreeAdd("IOBlock",ed);
	ed.SetInt(useExternalTime, "use_external_time_source"); ed.SetLocked(1); ed.SetToolTipText("1|(0) ... (Don't) use external input as time source."); edc.TreeAdd("IOBlock",ed);

}

int IOPulseGenerator::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = InputOutputElement::SetElementData(edc);

	GetElemDataDouble(mbs, edc, "IOBlock.amplitude", amplitude, 0);
	GetElemDataDouble(mbs, edc, "IOBlock.offset", toffs, 0);
	GetElemDataDouble(mbs, edc, "IOBlock.period", period, 0);
	GetElemDataDouble(mbs, edc, "IOBlock.pulse_width", pulseWidth, 0);
	GetElemDataInt(mbs, edc, "IOBlock.use_external_time_source", useExternalTime, 0);

	assert(period > 0. && pulseWidth >= 0.); // $ MSax 2013-02-28: added from set function

	inputs.SetLen(0);  //$ MSax 2013-04-03: added
	input_types.SetLen(0); //$ MSax 2013-04-03: added
	input_localnum.SetLen(0); //$ MSax 2013-04-03: added

	return rv;
}

void IOPulseGenerator::DrawBlockSymbol() 
{ 
	return InputOutputElement::DrawBlockSymbol();

	Vector3D refpos = GetRefPosD();
	Vector2D center(refpos.X(),refpos.Y());
	Vector2D size(draw_dim.X(), draw_dim.Y());

	Vector2D p1, p2;
	size *= 0.16; // scale overall size of the symbol
	for (int i=1; i<=7; i++)
	{
		switch (i)
		{
		case 1: p1 = center + Vector2D(-2.*size.X(),-2.*size.Y()); p2 = center + Vector2D(   -size.X(),-2.*size.Y()); break;
		case 2: p1 = center + Vector2D(   -size.X(),-2.*size.Y()); p2 = center + Vector2D(   -size.X(), 2.*size.Y()); break;
		case 3: p1 = center + Vector2D(   -size.X(), 2.*size.Y()); p2 = center + Vector2D(         0.0, 2.*size.Y()); break;
		case 4: p1 = center + Vector2D(         0.0, 2.*size.Y()); p2 = center + Vector2D(         0.0,-2.*size.Y()); break;
		case 5: p1 = center + Vector2D(         0.0,-2.*size.Y()); p2 = center + Vector2D(    size.X(),-2.*size.Y()); break;
		case 6: p1 = center + Vector2D(    size.X(),-2.*size.Y()); p2 = center + Vector2D(    size.X(), 2.*size.Y()); break;
		case 7: p1 = center + Vector2D(    size.X(), 2.*size.Y()); p2 = center + Vector2D( 2.*size.X(), 2.*size.Y()); break;
		default: break;
		}
		GetMBS()->AddDrawComponent_Line( GetOwnNum(), TIOBlock,											// MBS identifier
																		 1, TSymbol,																// DRAW identifier
																	   p1, p2,																		// extent
																		 colforeground );														// color
	}
}

//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
const char* IOProduct::SymbolText() const 
{
	return CSText_IOProduct;
}

void IOProduct::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	InputOutputElement::GetElementData(edc);
	ElementData ed; 
	ed.SetVector(exp.GetVecPtr(), exp.Length(), "exponents"); ed.SetVariableLength(); ed.SetToolTipText("Exponent of inputs. y=u1\\^exp1*u2\\^exp2*...*un\\^expn+offset.");  edc.TreeAdd("IOBlock",ed); //$ MSax2013-02-28
	ed.SetDouble(offset, "offset"); ed.SetToolTipText("Output offset.");  edc.TreeAdd("IOBlock",ed);
}

int IOProduct::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = InputOutputElement::SetElementData(edc);
	GetElemDataVector(mbs, edc, "IOBlock.exponents", exp, 0);
	GetElemDataDouble(mbs, edc, "IOBlock.offset", offset, 0);

	if(exp.Length() <= 0) //$ MSax 2013-02-28 added from set function
	{
		exp = Vector(1.,1.);
		mbs->UO(UO_LVL_err).InstantMessageText("Error in SetIOProduct: exponent vector is empty.");
	}
	return rv;
}

double IOProduct::GetOutput(double t, int i) const 
{
	assert(GetNInputs() == exp.Length());

	double outval = 1.;
	for(int i = 1; i<=GetNInputs(); i++)
	{
		// multiplication
		double u = GetInput(t, i);
		if(exp(i)<0)
		{
			assert(u != 0.);
		}
		////TODO: if necessary, add strategy if input is zero (e.g. like code below)
		//if(u == 0.)
		//{
		//	// division	
		//	u = divide_by_zero_tol; // omit division by zero
		//	mbs->UO(UO_LVL_warn) << "Warning in IOProduct: Division by zero. Divisor is set to tolerance\n";
		//}
		//else if(fabs(u) < divide_by_zero_tol)
		//{				
		//	u = u/fabs(u)*divide_by_zero_tol; // omit division by zero
		//	mbs->UO(UO_LVL_warn) << "Warning in IOProduct: Division by value close to zero. Divisor is set to tolerance\n";
		//}
		outval *= pow(u, exp(i));
		// just for debugging
		//if(1)
		//{
		//	mbs->UO() << "Product: u = " << u  << ", exp = " << exp(i) << "\n";
		//	if(i == GetNInputs())
		//	{					
		//		mbs->UO() << "Product: y = " << outval  << "\n\n";
		//	}
		//}
	}
	outval += offset;

	return outval;
}

void IOProduct::DrawBlockSymbol() 
{ 
	return InputOutputElement::DrawBlockSymbol();

	Vector3D refpos = GetRefPosD();
	Vector2D center(refpos.X(),refpos.Y());
	Vector2D size(draw_dim.X(), draw_dim.Y());


	size *= 0.16; // scale overall size of the symbol

	GetMBS()->AddDrawComponent_Line( GetOwnNum(), TIOBlock,												// MBS identifier
																	 1, TSymbol,																	// DRAW identifier
																	 center + Vector2D( size.X(), size.Y() ), center + Vector2D( -size.X(), -size.Y() ),	// extent
																	 colforeground );															// color
	GetMBS()->AddDrawComponent_Line( GetOwnNum(), TIOBlock,												// MBS identifier
																	 1, TSymbol,																	// DRAW identifier
																	 center + Vector2D( -size.X(), size.Y() ), center + Vector2D(  size.X(), -size.Y() ),	// extent
																	 colforeground );															// color
}
//
//+++++++Interface Data element: ControllerInterfaceData ++++++++++++++++++++++++++++++++++++++++
void ControllerInterfaceData::GetElementData(ElementDataContainer& edc)
{
	ElementData ed;
	//adresses to functions are not in user menu.
	//void (*Reg_Init)();	  // pointer on initialize function, not necessary to show these pointer in data container
	//void (*Reg_Update)();	// pointer on update function, not necessary to show these pointer in data container
	//ed.SetDouble(deltaT, "discrete_update_time");	edc.Add(ed); // already defined in InputOutputElementDiscrete::GetElementData(edc);
	//ed.SetInt(nOut, "number_of_outputs");	edc.Add(ed); // not necessary, already written into menu dialog from InputOutputElement::GetElementData(edc);
	//ed.SetInt(nInp, "number_of_inputs");	edc.Add(ed);// not necessary, already written into menu dialog from InputOutputElement::GetElementData(edc);	
}

int ControllerInterfaceData::SetElementData(ElementDataContainer& edc)
{
	int rv = 0;
	//adresses to functions are not in user menu.
	//void (*Reg_Init)();	  // pointer on initialize function
	//void (*Reg_Update)();	// pointer on update function
	//GetElemDataDouble(mbs, edc, "discrete_update_time", deltaT, 0);	 // already defined in InputOutputElementDiscrete::SetElementData(edc)
	//GetElemDataInt(mbs, edc, "number_of_outputs", nOut, 0);		 // not necessary, already written into menu dialog from InputOutputElement::GetElementData(edc);
	//GetElemDataInt(mbs, edc, "number_of_outputs", nInp, 0);		 // not necessary, already written into menu dialog from InputOutputElement::GetElementData(edc);
	return rv;
}

//extern char* CSText_ControllerInterface;
extern const char* CSText_ControllerInterface;
const char* ControllerInterface::SymbolText() const 
{
	return CSText_ControllerInterface;
}

int ControllerInterface::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = InputOutputElementDiscrete::CheckConsistency(errorstr);
	if(rv){	return rv;}	

	return rv;
}

void ControllerInterface::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	InputOutputElementDiscrete::GetElementData(edc);
	data.GetElementData(edc);
}

int ControllerInterface::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = InputOutputElementDiscrete::SetElementData(edc);
	rv += data.SetElementData(edc);
	return rv;
}

void ControllerInterface::DrawBlockSymbol() 
{ 
	InputOutputElement::DrawBlockSymbol();
}

// INTERFACE INTERFACE INTERFACE INTERFACE INTERFACE INTERFACE INTERFACE INTERFACE INTERFACE
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// IOTimeWindow
//$ DR 2011-12:[ IOTimeWindow added
//const char* CSText_IOTimeWindow = "Window"; //"IOTimeWindow"; // AD: moved to top of file ...

const char* IOTimeWindow::SymbolText() const 
{
	return CSText_IOTimeWindow;
}

int IOTimeWindow::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = InputOutputElement::CheckConsistency(errorstr);
	if(rv){	return rv;}

	return rv;
}

void IOTimeWindow::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	InputOutputElement::GetElementData(edc);
	ElementData ed;
	ed.SetDouble(t_start, "t_start"); ed.SetToolTipText("Start time (s)."); edc.TreeAdd("IOBlock",ed);
	ed.SetDouble(t_end, "t_end"); ed.SetToolTipText("End time (s)."); edc.TreeAdd("IOBlock",ed);
}

int IOTimeWindow::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = InputOutputElement::SetElementData(edc);

	GetElemDataDouble(mbs, edc, "IOBlock.t_start", t_start, 0);
	GetElemDataDouble(mbs, edc, "IOBlock.t_end", t_end, 0);

	if (t_end > t_start)
	{
		SetIOTimeWindow(t_start, t_end);
	}
	else
	{
		SetIOTimeWindow(t_start);
	}

	return rv;
}

double IOTimeWindow::GetOutput(double t, int i) const 
{
	// y = u(2)		if t_start <= u(1) <= t_start + delta_t
	// y = u(2)		if t_start <= u(1) and delta_t < 0
	// y = 0	else
	
	if(reached_end)
	{
		return 0;	
	}
	else
	{
		if (GetInput(t, 1) < t_start )
		{
			return 0;
		}
		else 
		{
			if(delta_t < 0.)
			{
				return GetInput(t, 2);
			}
			else if (t_start + delta_t < GetInput(t, 1))
			{
				reached_end = 1;
				return 0;
			}
			else
			{
				return GetInput(t, 2);
			}
		}
	}
}

void IOTimeWindow::DrawBlockSymbol() 
{ 
	InputOutputElement::DrawBlockSymbol();
}

//$ DR 2011-12:] IOTimeWindow added
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
const char* IOStopComputation::SymbolText() const 
{
	return CSText_IOStopComputation;
}

void IOStopComputation::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	InputOutputElement::GetElementData(edc);

	ElementData ed;
}

int IOStopComputation::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = InputOutputElement::SetElementData(edc);

	return rv;
}

void IOStopComputation::DrawBlockSymbol() 
{ 
	InputOutputElement::DrawBlockSymbol();
}

//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// functions for class IOElementDataModifier

Element* IOElementDataModifier::GetCopy()
{
	Element* ec = new IOElementDataModifier(mbs);
	ec->CopyFrom(*this);

	return ec;
}

void IOElementDataModifier::CopyFrom(const Element& e)
{
	InputOutputElement::CopyFrom(e);
	const IOElementDataModifier& ce = (const IOElementDataModifier&)e;

	RWdata = ce.RWdata;
	element_number = ce.element_number;
	variable_name = ce.variable_name;
	modify_at_start_time_step_only = ce.modify_at_start_time_step_only;
}

void IOElementDataModifier::InitConstructor()
{
	InputOutputElement::InitConstructor();

	type = TMBSElement(type + TIOElementDataModifier);

	elementname = GetElementSpec();
	
	RWdata.comp1 = 0;
	RWdata.comp2 = 0;
	RWdata.value = 0.;

	RWdata.RWaccess = TRWElementDataNoAccess;
	element_number = 1;
	variable_name = mystr("");

	modify_at_start_time_step_only = 0;
}

const char* IOElementDataModifier::SymbolText() const 
{
	return CSText_IOElementDataModifier;
}

void IOElementDataModifier::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	InputOutputElement::GetElementData(edc);
  ElementData ed;
  ed.SetText(variable_name.c_str(),"mod_variable_name"); ed.SetToolTipText("variable name of the modified element data"); edc.TreeAdd("IOBlock",ed);
  ed.SetInt(element_number,"mod_element_number"); ed.SetToolTipText("element number of the modified element or constraint"); edc.TreeAdd("IOBlock",ed);

	//ElementData ed;
}

int IOElementDataModifier::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	//TODO: resize matrix, if user change the size!!!
	int rv = InputOutputElement::SetElementData(edc);

	ElementData ed;
  GetElemDataText(GetMBS(), edc, "IOBlock.mod_variable_name",variable_name, 1);
	GetElemDataInt(GetMBS(), edc, "IOBlock.mod_element_number",element_number, 1);

	SetIOElementDataModifier(element_number, variable_name.c_str());
	
	return rv;
}

void IOElementDataModifier::DrawBlockSymbol() 
{ 
	return InputOutputElement::DrawBlockSymbol();

	Vector3D refpos = GetRefPosD();
	Vector2D center(refpos.X(),refpos.Y());
	Vector2D size(draw_dim.X(), draw_dim.Y());


	//dial
	Vector2D dial = center + Vector2D( 0., -0.1*draw_dim.Y());
	GetMBS()->AddDrawComponent_Ellipse( GetOwnNum(), TIOBlock,										// MBS identifier
																			1, TSymbol,																// DRAW identifier
																			dial , 0.5*size,													// extent
																			Vector3D(0.,0.,0.), Vector3D(.5,.5,.5));	// color
	
	GetMBS()->AddDrawComponent_Line( GetOwnNum(), TIOBlock,												// MBS identifier
																	 1, TSymbol,																	// DRAW identifier
																	 dial+Vector2D(0., 0.15*size.Y()) , dial+Vector2D(0., 0.25*size.Y()),							// extent
																	 Vector3D(0.,0.,0.) );												// color

	int nr_of_ticks = 7;
	double angle_first=0.;
	double angle_last=180.;
	for (int i=1; i<=7; i++)
	{
		double angle = (MY_PI / 180.) * ((angle_last-angle_first) / (double) (nr_of_ticks-1) ) * (i-1.);
		Vector2D dir( size.X()*cos(angle), size.Y()*sin(angle) );

		GetMBS()->AddDrawComponent_Line( GetOwnNum(), TIOBlock,											// MBS identifier
																		 1, TConnectionLine,												// DRAW identifier
																		 dial + 0.3*dir , dial + 0.4*dir,						// extent
																		 Vector3D(0.,0.,0.) );											// color

	}
}


// SM: 10012013
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
const char* IODisplay::SymbolText() const 
{
	return CSText_IODisplay;
}

void IODisplay::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	InputOutputElement::GetElementData(edc);
	ElementData ed;
}

int IODisplay::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = InputOutputElement::SetElementData(edc);
	return rv;
}

void IODisplay::DrawBlockSymbol()
{
	Vector3D refpos = GetRefPosD();
	Vector2D center(refpos.X(),refpos.Y());
	Vector2D size(draw_dim.X(), draw_dim.Y());

	double t = GetMBS()->GetDrawTime();
	double value = XDataD(1); // $ MSax 2013-03-01: added
	mystr text(value,ndigits);

	GetMBS()->AddDrawComponent_Text( GetOwnNum(), TIOBlock,												// MBS identifier
																	 1, TSymbol,																	// DRAW identifier
																	 center, 0.8 * size,													// extent
															     Vector3D(0.,0.,0.),													// color of this label is alway black
																	 text,															          // display value
																	 TTextAllign (HCenter+VCenter) );							// other properties
};


// ELEMENTS for Response to input Key or Mouse during computation
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
const char* IOResponseElement::SymbolText() const 
{
	return CSText_IOResponseElement;
}

void IOResponseElement::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	InputOutputElement::GetElementData(edc);
	ElementData ed;
}

int IOResponseElement::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = InputOutputElement::SetElementData(edc);
	return rv;
}

void IOResponseElement::DrawBlockSymbol() 
{
	Vector3D refpos = GetRefPosD();
	Vector2D center(refpos.X(), refpos.Y());
	Vector2D size(draw_dim.X(), draw_dim.Y());

	int drawn_value = XData(1);

// display current value
	if (IsToggle())
	{
		Vector3D color_onoff;
		if (drawn_value == 0) color_onoff = Vector3D(1.,0.,0.);
		if (drawn_value == 1) color_onoff = Vector3D(0.,1.,0.);

		GetMBS()->AddDrawComponent_Rect( GetOwnNum(), TIOBlock,												// MBS identifier
																		 1, TSymbol,																	// DRAW identifier
																		 center, 0.5*size,														// extent
																		 Vector3D(0.,0.,0.), color_onoff );           // colors
	}
	else
	{
		mystr text_curr(drawn_value);
		GetMBS()->AddDrawComponent_Text( GetOwnNum(), TIOBlock,												// MBS identifier
																		 1, TSymbol,																	// DRAW identifier
																		 center, 0.8 * size,													// extent
																		 Vector3D(0.,0.,0.),													// color of this label is alway black
																		 text_curr,													          // display value
																		 TTextAllign (HCenter+VCenter) );							// other properties
	}

// marker mode in upper left corner
	mystr text_mode;
	Vector2D mode_m(refpos.X()-0.4*draw_dim.X(), refpos.Y()+0.4*draw_dim.Y());
	Vector2D mode_s(0.2*draw_dim.X(), 0.2*draw_dim.Y());
	if (GetInputMode() == TMouseResponse) text_mode = mystr("M");
	if (GetInputMode() == TKeyResponse) text_mode = mystr("K");
	GetMBS()->AddDrawComponent_Text( GetOwnNum(), TIOBlock,												// MBS identifier
																	 2, TSymbol,																	// DRAW identifier
																	 mode_m, mode_s,															// extent
															     Vector3D(0.,0.,0.),													// color of this label is alway black
																	 text_mode,													          // display value
																	 TTextAllign (HCenter+VTop) );								// other properties
// marker range at bottom
	mystr text_range;
	Vector2D range_c(refpos.X(), refpos.Y()-0.35*draw_dim.Y());
	Vector2D range_s(draw_dim.X(), 0.2*draw_dim.Y());
	if (IsToggle()) text_range = mystr("[on/off]");
	else text_range = mystr("[") + mystr(lower_bound) + mystr(",") + mystr(upper_bound) + mystr("]");
	GetMBS()->AddDrawComponent_Text( GetOwnNum(), TIOBlock,												// MBS identifier
																	 3, TSymbol,																	// DRAW identifier
																	 range_c,	range_s,														// extent
															     Vector3D(0.,0.,0.),													// color of this label is alway black
																	 text_range,												          // display value
																	 TTextAllign (HCenter+VCenter) );							// other properties
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
const char* IOMinMax::SymbolText() const 
{
	return CSText_IOMinMax;
}

int IOMinMax::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = InputOutputElement::CheckConsistency(errorstr);
	if(rv){	return rv;}

	return rv;
}

void IOMinMax::StartTimeStep()
{
	double t = mbs->GetTime();
	double u = GetInput(t, 1);

	if(t<start_time)
	{
		current_value = u;
	}
	else
	{
		if(first_time_step)
		{
			current_value = u;
			if(mode >= 4) {current_value = fabs(current_value);}
			first_time_step = 0;
			WriteToXData();
		}
		else
		{
			switch(mode)
			{
			case 1:	{current_value = min(XData(1),u); break;}
			case 2:	{current_value = max(XData(1),u); break;}
			case 3:	{current_value = XData(1)+u; break;}
			case 4:	{current_value = min(XData(1),fabs(u)); break;}
			case 5:	{current_value = max(XData(1),fabs(u)); break;}
			case 6:	{current_value = XData(1)+fabs(u); break;}
			}
			WriteToXData();
		}
	}
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// functions for class IOTCPIPBlock

Element* IOTCPIPBlock::GetCopy()
{
	Element* ec = new IOTCPIPBlock(mbs);
	ec->CopyFrom(*this);

	return ec;
}

void IOTCPIPBlock::CopyFrom(const Element& e)
{
	InputOutputElement::CopyFrom(e);
	const IOTCPIPBlock& ce = (const IOTCPIPBlock&)e;

	port_number = ce.port_number;
	ip_address = ce.ip_address;
	timeout = ce.timeout;
	serversocket = ce.serversocket;
	tcpdata = ce.tcpdata;
	addtcpdata = ce.addtcpdata;
	use_default_protocol = ce.use_default_protocol;
	use_additional_communication = ce.use_additional_communication;
	separate_additional_communication = ce.separate_additional_communication;
	additional_communication_length = ce.additional_communication_length;
	additional_communication_unit_size = ce.additional_communication_unit_size;
	isinitialized = ce.isinitialized;

}

void IOTCPIPBlock::InitConstructor()
{
	InputOutputElement::InitConstructor();

	//type = TMBSElement(type + TIOTCPIPBlock);

	elementname = GetElementSpec();
	port_number = 50000;
	ip_address = mystr("127.0.0.1");
	use_default_protocol = 1;
	use_additional_communication = 1;
	separate_additional_communication = 0;
	additional_communication_length = 1;
	additional_communication_unit_size = 8;
	timeout = 10000;
	isinitialized = 0;
}

const char* IOTCPIPBlock::SymbolText() const 
{
	return CSText_IOTCPIPBlock;
}

int IOTCPIPBlock::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = InputOutputElement::SetElementData(edc);

	Vector datainit(DataS());
	datainit.SetAll(0.);
	SetDataInit(datainit);

	return rv;
}

void IOTCPIPBlock::Initialize() // allocate memory and initialize TCP/IP connection
{ 
	double tnew;
	static double told;
	if(!isinitialized)
	{
		told=-1;
		//ensure consistency
		if(use_default_protocol)
		{
			use_additional_communication = 1;
			separate_additional_communication = 0;
			additional_communication_length = 1;
			additional_communication_unit_size = 8;
		}
		else if(!use_additional_communication)
		{
			separate_additional_communication = 0;
			additional_communication_length = 0;
			additional_communication_unit_size = 8;
		}
		if(additional_communication_unit_size!=1 && additional_communication_unit_size!=2 && additional_communication_unit_size!=4 && additional_communication_unit_size!=8)
		{
			additional_communication_unit_size = 8;
			mbs->UO(UO_LVL_warn) << "Warning in IOTCPIPBlock: additional_communication_unit_size was reset to default (8 bytes)\n";
		}

		tcpdata = new double[(GetNInputs()>GetNOutputs() ? GetNInputs() : GetNOutputs()) + additional_communication_length + 1 ];

		if(use_additional_communication)
		{
			addtcpdata = new char[additional_communication_length*additional_communication_unit_size];
		}			

		if(!serversocket)
		{
			serversocket = new TCPIPHotInt;
			serversocket->Set(mbs,SocketType::ST_server,ip_address.c_str(),port_number,1,1);
			serversocket->SetAutoReconnect(0);
			serversocket->SetAcceptTimeOut(timeout);
			serversocket->SetUpConnection();
			serversocket->SetRecvTimeOut(timeout);
			serversocket->SetSendTimeOut(timeout);
		}

		isinitialized = 1;
	}

	if(use_default_protocol)
	{
		//set reset-flag in additional communication per default
		SetCommunicationFlag(0,Freset,1);
		//communication of the reset flag and initial values
		tnew = mbs->GetTime(); //ensure that this is only done once per simulation run (single computation or one computation of a parameter variation)
		if(abs(tnew-told)>1E-14)
		{
			Communicate(tnew);
			told=tnew;
		}
		SetCommunicationFlag(0,Fneutral,1);
	}
	else
	{
		//to be implemented as required
	}
}

void IOTCPIPBlock::CloseAndCleanUp()
{
	if(isinitialized && mbs->GetSimulationStatus().GetStatusFlag(TSimulationProcessFinished))
	{
		SetCommunicationFlag(0,Fclose,1);
		OutgoingDataCommunication(GetMBS()->GetTime());
		delete [] tcpdata;
		if(use_additional_communication)
			delete [] addtcpdata;
		
		if(serversocket)
		{
			serversocket->CloseConnection();
			delete [] serversocket;
			serversocket = NULL;
		}

		isinitialized = 0;
	}
}

IOTCPIPBlock::~IOTCPIPBlock()
{
	CloseAndCleanUp();
}

void IOTCPIPBlock::ComputationFinished()
{
	CloseAndCleanUp();
	//if(isinitialized)
	//{
	//	if(GetMBS()->GetMBS_EDC_Options()->TreeGetInt("SolverOptions.ParameterVariation.activate"))
	//	{
	//		//without the following, when the parameter variation is finished, the client will just run into a recv timeout

	//		//if(parameter variation is finished)
	//		//{
	//		//	SetCommunicationFlag(0,Fclose,1);
	//		//	OutgoingDataCommunication(GetMBS()->GetStepEndTime());
	//		//}
	//	}
	//	else
	//	{
	//		//SetCommunicationFlag(0,Fclose,1);
	//		//OutgoingDataCommunication(GetMBS()->GetStepEndTime());
	//		//TSimulationStatus stat = mbs->GetSimulationStatus();
	//		CloseAndCleanUp();
	//	}
	//}
}

void IOTCPIPBlock::OutgoingDataCommunication(double t)
{
	//outgoing data transfer
	*(tcpdata)=t;
	int nin = GetNInputs();
	for(int i=0; i<nin; ++i)
		*(tcpdata+i+1)=GetInput(t,i+1);
	if(use_additional_communication && !separate_additional_communication)
	{
		double* temp = (double*) addtcpdata;
		for(int i=0; i<additional_communication_length; ++i)
			*(tcpdata+nin+i+1) = *(temp+i);
	}
	serversocket->SendData((char*)tcpdata,(1+nin+additional_communication_length)*8);
}

void IOTCPIPBlock::IncomingDataCommunication()
{
	//incoming data transfer
	int nout = GetNOutputs();
	serversocket->RecvData((char*)tcpdata,(nout+additional_communication_length)*8);
	for(int i=0; i<nout; ++i)
		SetOutput(*(tcpdata+i),i+1);
	if(use_additional_communication && !separate_additional_communication)
	{
		double* temp = (double*) addtcpdata;
		for(int i=0; i<additional_communication_length; ++i)
			*(temp+i) = *(tcpdata+nout+i);
	}	
}

void IOTCPIPBlock::Communicate(double t)
{
	if(use_additional_communication && !separate_additional_communication)
	{
		OutgoingDataCommunication(t);
		IncomingDataCommunication();
		ReactToAdditionalCommunication(t);	
	}
	else if(use_additional_communication && separate_additional_communication)
	{
		serversocket->SendData(addtcpdata,additional_communication_length*additional_communication_unit_size);
		serversocket->RecvData(addtcpdata,additional_communication_length*additional_communication_unit_size);
		int commdata = ReactToAdditionalCommunication(t);
		if(commdata)
		{
			OutgoingDataCommunication(t);
			IncomingDataCommunication();
		}
	}
	else
	{
		OutgoingDataCommunication(t);
		IncomingDataCommunication();
	}
}

void IOTCPIPBlock::SetCommunicationFlag(int pos, CommunicationFlag flag, int convertdouble)
{
	if(convertdouble)
	{
		double temp = (double)(int)flag;
		*(double*)(addtcpdata+pos) = temp;
	}
	else
		*(int*)(addtcpdata+pos) = (int)flag;
}

CommunicationFlag IOTCPIPBlock::GetCommunicationFlag(int pos, int convertdouble)
{
	if(convertdouble)
	{
		return (CommunicationFlag)(int)*(double*)(addtcpdata+pos);
	}
	else
		return (CommunicationFlag)(*(int*)(addtcpdata+pos));
}

//called by Communicate(double t)
//read received data and respond somehow
//possibly rewrite addtcpdata and call Communicate(t) recursively for complex communication schemes
int IOTCPIPBlock::ReactToAdditionalCommunication(double t)
{
	//default
	if(use_default_protocol)
	{
		switch(GetCommunicationFlag(0,1))
		{
			case Fneutral:
				break;
			case Ferror:
				mbs->UO().InstantMessageText("An error occurred on the client side!\n");
				assert(0);
				exit(0);
				break;
			default:
				mbs->UO().InstantMessageText("Received an unknown command from the client!\n");
				assert(0);
				exit(0);
				break;
		}
	}
	else
	{
		//to be implemented as required
	}

	return 1; //returns a flag deciding if DataCommunication(t) should be called afterwards in case of separate addíonal communication
}
