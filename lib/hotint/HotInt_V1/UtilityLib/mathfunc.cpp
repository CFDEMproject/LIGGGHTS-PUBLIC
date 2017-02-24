//#**************************************************************
//#
//# filename:             mathfunc.cpp
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						October 2006
//# description:          general mathematical function
//#                       mostly for time or space-curves of arbitrary shape
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
//#**************************************************************

//$JG2012-01: because of elementdata access, mbs is included now:
//#include "../WorkingModule/stdafx.h"
//#include "ioincludes.h"
//
//#include <assert.h>
//#include <memory.h>
//
//#include <math.h>
//
//#include "tarray.h"    
//#include "mystring.h"  
//#include "elementdata.h"
//#include "femath.h"    
//#include "myfile.h"
//
////for output possibility within femath
//#include "..\workingmodule\WorkingModuleBaseClass.h"

#include "mbs_interface.h"
#include "../Parser/parser.h"
#include "mathfunc.h"
#include "elementdataaccess.h"
#include "myfile.h"

extern UserOutputInterface * global_uo;


char* MathFuncStringList[] = {"NoFunction", "PieceConstant", "PieceLinear", "PieceQuadratic", "UserDefined", "Expression", "Polynomial", "Harmonic", "Sin", "Cos", "StaticDynamicFricion"};
int MathFuncStringListLength = 6; //*** number of elements in MathFuncStringList, 

void MathFunction::CopyFrom(const MathFunction& e)
{
	funcmode = e.funcmode;
	dim = e.dim;
	vectime = e.vectime;
	valX = e.valX;
	valY = e.valY;
	valZ = e.valZ;
	coeff = e.coeff;
	ufunc = e.ufunc;
	parsedFunctionExpression = e.parsedFunctionExpression;
	parsedFunctionVariables = e.parsedFunctionVariables;
	mbsPI = e.mbsPI;

	//$JG2013-4-29: copy parsed function
	pf = e.pf;
}

const char* MathFunction::GetTypeName() const
{
	if (funcmode < MathFuncStringListLength)
	{
		return MathFuncStringList[(int)funcmode];
	}
	else
	{
		return "Unknown";
	}
}

TMathFuncType MathFunction::StringToMathFuncType(const mystr& funcname)
{
	for (int i=0; i < MathFuncStringListLength; i++)
	{
		if (mystr(MathFuncStringList[i]) == funcname)
		{
			return (TMathFuncType)i;
		}
	}
	return TMFempty;
}

void MathFunction::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	ElementData ed;

	//structure of ElementData:
	//Mathfunction
	//{
	//  piecewise_mode = 0/1/2     %modus for piecewise interpolation: -1=not piecewise==>use parsed function, 0=constant, 1=linear, 2=quadratic
	//  piecewise_points = [ ... ] %supporting points (e.g. time or place) for piecewise interpolation
	//  piecewise_values = [ ... ] %values at supporting points
	//  piecewise_diff_values = [ ... ] %differential values at supporting points - for quadratic interpolation
	//  parsed_function = " ... " %string representing parsed function, e.g. "A*sin(omega*t)"
	//  parsed_function_parameter = " ... " %string representing parameter of parsed function, e.g. "t"
	//  user_defined_function = yes %not editable, hard-coded userdefined function!
  //}
	int piecewise_mode = -1;

	//if (funcmode == TMFpiecewiseconst || funcmode == TMFpiecewiselinear || funcmode == TMFpiecewisequad) 
	//{
	//}
	if (funcmode == TMFpiecewiseconst) 
	{
		piecewise_mode = 0;
	}
	else if (funcmode == TMFpiecewiselinear) 
	{
		piecewise_mode = 1;
	}
	else if (funcmode == TMFpiecewisequad) 
	{
		piecewise_mode = 2;
	}
	ed.SetInt(piecewise_mode, "piecewise_mode", -1, 2); ed.SetToolTipText("modus for piecewise interpolation: -1=not piecewise, 0=constant, 1=linear, 2=quadratic"); edc.Add(ed);
	ed.SetVector(vectime.GetVecPtr(), vectime.Length(), "piecewise_points"); ed.SetVariableLength(); ed.SetToolTipText("supporting points (e.g. time or place) for piecewise interpolation"); edc.Add(ed);
	ed.SetVector(valX.GetVecPtr(), valX.Length(), "piecewise_values"); ed.SetVariableLength(); ed.SetToolTipText("values at supporting points"); edc.Add(ed);
	ed.SetVector(valY.GetVecPtr(), valY.Length(), "piecewise_diff_values"); ed.SetVariableLength(); ed.SetToolTipText("differential values at supporting points - for quadratic interpolation"); edc.Add(ed);


	if (funcmode == TMFexpression)
	{
		ed.SetText(parsedFunctionExpression.c_str(), "parsed_function"); ed.SetToolTipText("string representing parsed function, e.g. 'A*sin(omega*t)'"); edc.Add(ed);
		ed.SetText(parsedFunctionVariables.c_str(), "parsed_function_parameter"); ed.SetToolTipText("string representing parameter of parsed function, e.g. 't'"); edc.Add(ed);
	}
	else
	{
		//initialize with zero strings:
		ed.SetText("", "parsed_function"); ed.SetToolTipText("string representing parsed function, e.g. 'A*sin(omega*t)'"); edc.Add(ed);
		ed.SetText("", "parsed_function_parameter"); ed.SetToolTipText("string representing parameter of parsed function, e.g. 't'"); edc.Add(ed);
	}

	if (funcmode == TMFuserdefined)
	{
		ed.SetBool(1,"user_defined_function"); ed.SetLocked(1); ed.SetToolTipText("not editable, hard-coded userdefined function!"); edc.Add(ed);
	}
}

int MathFunction::SetElementData(MBS* mbs, ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = 1;

	int udf = 0; //user defined function
	GetElemDataBool(mbs, edc, "user_defined_function", udf, 0); //do not warn, if not exists
	if (udf) return rv; //this case is not treated!

	int piecewise_mode;
	Vector t, u, v;
	GetElemDataInt(mbs, edc, "piecewise_mode", piecewise_mode, 1);
	GetElemDataVector(mbs, edc, "piecewise_points", t, 1);
	GetElemDataVector(mbs, edc, "piecewise_values", u, 1);
	GetElemDataVector(mbs, edc, "piecewise_diff_values", v, 1);

	mystr pf, pfp;
	GetElemDataText(mbs, edc, "parsed_function", pf, 1);
	GetElemDataText(mbs, edc, "parsed_function_parameter", pfp, 1);
	
	if (piecewise_mode == 0) {SetPiecewise(t, u, 0);}
	else if (piecewise_mode == 1) {SetPiecewise(t, u, 1);}
	else if (piecewise_mode == 2) {SetPiecewiseQuadratic(t, u, v);}
	else
	{
		SetExpression(pf, pfp,mbs);
	}

	return rv;
}

//the following two functions is out-dated, use Get/SetElementData instead!
int MathFunction::SetData(int mode, const Matrix& data)
{
	int rv = 1;
	funcmode = mode;
	if (funcmode == TMFpolynomial)
	{
		if (data.Getcols() == 1 && data.Getrows() > 0)
			data.GetColVec(1, coeff);
		else
		{
			funcmode = TMFempty;
			rv = 0;
		}
	}
	else if (funcmode == TMFpiecewiseconst || funcmode == TMFpiecewiselinear) //piecewise constant/linear
	{
		if (data.Getcols() == 2 && data.Getrows() > 0)
		{
			data.GetColVec(1, vectime);
			data.GetColVec(2, valX);
		}
		else
		{
			funcmode = TMFempty;
			rv = 0;
		}
	}
	else if (funcmode == TMFpiecewisequad) //piecewise quadratic
	{
		if (data.Getcols() == 3 && data.Getrows() > 0)
		{
			data.GetColVec(1, vectime);
			data.GetColVec(2, valX); //pos
			data.GetColVec(3, valY); //vel
		}
		else
		{
			funcmode = TMFempty;
			rv = 0;
		}
	}
	else if (funcmode == TMFharmonic) //harmonic
	{
		if (data.Getcols() == 3 && data.Getrows() > 0)
		{
			data.GetColVec(1, valX); //frequency
			data.GetColVec(2, valY); //phase
			data.GetColVec(3, valZ); //amplitude
		}
		else
		{
			funcmode = TMFempty;
			rv = 0;
		}
	}
	return rv;
}

//the following two functions is out-dated, use Get/SetElementData instead!
void MathFunction::GetData(int& mode, Matrix& data) const
{
	mode = funcmode;
	if (funcmode == TMFempty)
	{
		data.SetSize(0,0);
	}
	else if (funcmode == TMFpolynomial)
	{
		if (coeff.Length() != 0)
		{
			data.SetSize(coeff.Length(),1);
			data.SetColVec(coeff, 1);
		}
		else data.SetSize(0,0);
	}
	else if (funcmode == TMFpiecewiseconst || funcmode == TMFpiecewiselinear) //piecewise constant/linear
	{
		if (vectime.Length() != 0)
		{
			data.SetSize(vectime.Length(),2);
			data.SetColVec(vectime, 1);
			data.SetColVec(valX, 2);
		}
		else data.SetSize(0,0);
	}
	else if (funcmode == TMFpiecewisequad) //piecewise quadratic
	{
		if (vectime.Length() != 0)
		{
			data.SetSize(vectime.Length(),3);
			data.SetColVec(vectime, 1);
			data.SetColVec(valX, 2); //pos
			data.SetColVec(valY, 3); //vel
		}
		else data.SetSize(0,0);
	}
	else if (funcmode == TMFharmonic) //harmonic
	{
		if (valX.Length() != 0)
		{
			data.SetSize(valX.Length(),3);
			data.SetColVec(valX, 1);
			data.SetColVec(valY, 2);
			data.SetColVec(valZ, 3);
		}
		else data.SetSize(0,0);
	}
}

int MathFunction::SetData(const mystr& funcname, const Vector& coefficients)
{
	int rv = 1; //1=no error
	coeff = coefficients;
	funcmode = StringToMathFuncType(funcname);
	
	return rv;
}

void MathFunction::GetData(mystr& funcname, Vector& coefficients) const//return functionname
{
	funcname = GetTypeName();
	coefficients = coeff;
}



void MathFunction::SetPiecewise(const Vector& times, const Vector& coeffs, int interp)
{
	if (interp == 0) funcmode = TMFpiecewiseconst;					//piecewise constant, with jumps
	else if (interp == 1) funcmode = TMFpiecewiselinear;		//piecewise linear

	vectime = times;
	valX = coeffs;
	if (valX.Length() != vectime.Length()) valX.SetLen(vectime.Length());
}

void MathFunction::SetPiecewiseQuadratic(const Vector& times, const Vector& coeffs_p, const Vector& coeffs_v)
{
	funcmode = TMFpiecewisequad;			//piecewise quadratic

	vectime = times;
	valX = coeffs_p;
	valY = coeffs_v;

	if (valX.Length() != vectime.Length()) valX.SetLen(vectime.Length());
	if (valY.Length() != vectime.Length()) valY.SetLen(vectime.Length());
}

int MathFunction::SetPiecewiseFromFile(const char* filename, int ncolumns, int column1, int column2, int interp, int nrOfHeaderLines, double offset1, double offset2, int offset1_start_index, int offset2_start_index)
{
	// *** retvalues: -1: not found | 1: #error | 0: OK
	int rv = 0; //OK
	CMFile ifile(filename, TFMread);
	TArray<double> data_col1, data_col2;
	rv=ifile.ReadTwoColumnsFromFile(data_col1, data_col2,ncolumns, column1, column2, nrOfHeaderLines, offset1, offset2, offset1_start_index, offset2_start_index);
	if(rv)
	{
		return rv; //error
	}
	else
	{
		SetPiecewise(Vector(data_col1), Vector(data_col2), interp);
	}
	return rv; //OK
}	

int MathFunction::SetPiecewiseFromFile2(const char* filename, int column1, int column2, int interp, mystr comment)
{
	int rv = 0; // OK
	CMatrixFile cf(filename, (TFileMode)TFMread);
	if(!cf.IsGood()) 
	{
		return 1; // error
	}
	cf.ReadSpecifiedColumns(IntVec2(column1,column2));
	SetPiecewise(Vector(cf.Column(column1)), Vector(cf.Column(column2)), interp);
	return rv;
}


//f = c1*t^0 + c2*t^1 + c3*t^2 + ....
void MathFunction::SetPolynomial(const Vector& coeffs)
{
	coeff = coeffs;
	funcmode = TMFpolynomial;
}

void MathFunction::SetConstant(double x)
{
	coeff = Vector(x);
	funcmode = TMFpolynomial;
}

void MathFunction::SetVecTimes(const Vector& times) //for polynomial turn on/off times
{
	vectime = times;
}

void MathFunction::SetHarmonic(const Vector& frequency, const Vector& phase, const Vector& amplitude)
{
	valX = frequency;
	valY = phase;
	valZ = amplitude;

	if (valY.Length() != valX.Length()) valY.SetLen(valX.Length());
	if (valZ.Length() != valX.Length()) valZ.SetLen(valX.Length());

	funcmode = TMFharmonic; //harmonic function
}

//void MathFunction::SetSinFunction(double amplitude, double freq_rad, double phase) 
//{
//	funcmode = TMFsin;
//	coeff.SetLen(3);
//	coeff(1) = amplitude;
//	coeff(2) = freq_rad;
//	coeff(3) = phase;
//}
//
//void MathFunction::SetCosFunction(double amplitude, double freq_rad, double phase)  //***
//{
//	funcmode = TMFcos;
//	coeff.SetLen(3);
//	coeff(1) = amplitude;
//	coeff(2) = freq_rad;
//	coeff(3) = phase;
//}

//void MathFunction::SetStaticDynamicFricionFunction(double staticFriction, double viscousFriction, double zeroZone) //***
//{
//	funcmode = TMFstaticDynamicFricion;
//	coeff.SetLen(4);
//	coeff(1) = staticFriction;
//	coeff(2) = viscousFriction;
//	coeff(3) = zeroZone;   // to avoid problems around speed 0 (no step function around zero velocity ==> ramp)
//	coeff(4) = staticFriction + viscousFriction * zeroZone; // "corner"
//}

// sets a userdefined function, additional parameters cam be passed in Vector coeffs
void MathFunction::SetUserDefined(double (*func)(double,const Vector&),const Vector& coeffs)
{
	funcmode = TMFuserdefined;
	ufunc = func;
	coeff = coeffs;
}

// sets a userdefined expression with one variable - uses parseobect
void MathFunction::SetExpression(mystr expression, mystr variable, MBS * mbs)
{
	funcmode = TMFexpression;
 
	//$!AD: Elements with forcemode may end up here with unused mathfunctions ( expression "", variable "" ), so do not set the parsedfunction then...
	if( expression.Length() > 0 && variable.Length() > 0)
	{
		pf.SetParsedFunction1D(&(mbs->MBSParser()), expression, variable); //$JG2013-4-29: 
	}

	//erase two lines: //$JG2013-04-29:
	parsedFunctionExpression = expression;
	parsedFunctionVariables = variable; 

	mbsPI = mbs;
}

int MathFunction::FindIndexPiecewise(double t) const //$ MS 2011-5-23
{
	switch (funcmode)
	{
	case TMFpiecewiselinear:
		{
			int i = 1;

			int inc = vectime.Length();
			double linc = floor(log((double)inc)/log(2.));
			inc = (int)pow(2,linc);
			i = inc;
			if (i < 1) i=1;

			while (inc >= 1) 
			{
				inc /= 2;
				if (t > vectime(i))
				{
					if((i+inc) <= vectime.Length()) 
					{
						i+=inc;
					}
				}
				else if ((i >= 2) && (t <= vectime(i-1))) 
				{
					i-=inc;
				}

				if (i<1) i=1;
				if (i>vectime.Length()) i=vectime.Length();
			}

			//$ LA 2011-06-01 is this really necessary?
			//$ RL 2012-7-24: yes, the returned index 'i' is element of the set {2,3,...,vectime.Length()}
			if (i > vectime.Length()) 
			{
				i = vectime.Length();
			}

			if(i < 2){i = 2;}//$ RL 2012-7-24: for extrapolation of values left of vectime(1)
			
			return i;
			break;
		}

	case TMFpiecewiseconst:
		{
			int i = 1;

			while (i <= vectime.Length() && t >= vectime(i)) {i++;}
			i--;
			if (i < 1) {i = 1;}
			return i;
			break;
		}
	case TMFstaticDynamicFricion:
		{
			//input variable: velocity
			assert(0 && "ERROR in mathfunc.cpp: SetStaticDynamicFricionFunction not anymore supported by MathFunction, use TMFExpression instead!"); 

			if(coeff(3) > 0.0 && fabs(t) <= coeff(3))//$ RL 2011-5-23:  code from GetSegmentNumber(double t, int diff=0) moved to FindIndexPiecewise
			{
				return 0;  // zero zone
			}else
			{				
				return 1;  // static + viscous friction
			}
			break;
		}
	default:
		{
			(*global_uo) << "ERROR: MathFunction::FindIndexPiecewise not defined for this Function-Mode! \n";
			return 0;
			break;
		}
	}
}

double MathFunction::InterpolatePiecewise(double t, int index) const //$ MS 2011-5-23
{
	switch (funcmode)
	{
	case TMFpiecewiselinear:
		{
			//if (index == vectime.Length()) 
			//{
			//	return valX(index);
			//}

			//$ RL 2012-7-24:[ commented out
			//if (index > vectime.Length()) //$ LA 2011-06-01 last interval also piecewise linear!
			//{
			//	index = vectime.Length();
			//	return valX(index);
			//}
			//else
			//{
			//
			//	double t1 = 0;
			//	double v1 = 0;
			//	if (index >= 2) 
			//	{
			//		t1 = vectime(index-1);
			//		v1 = valX(index-1);
			//	}
			//$ RL 2012-7-24:]
				
			
				//$ RL 2012-7-24: index is element of the set {2,3,...,vectime.Length()}
				double t1 = vectime(index-1);
				double v1 = valX(index-1);
				double dt = vectime(index)-t1;
				if (dt <= 0) dt = 1;
				double val = (1.-(t-t1)/dt)*v1 + (1.+(t-vectime(index))/dt)*valX(index);
				return val;
			//$ RL 2012-7-24: commented out: }
			break;
		}
	default:
		{
			(*global_uo) << "MathFunction::InterpolatePiecewise not defined for this Function-Mode! \n";
			return 0;
			break;
		}
	}
}

double MathFunction::Evaluate(double t, int diff) const //if diff=1, then differentiate w.r.t. time
{
	//t .. time or input parameter
	double val = 0;
	switch (funcmode)
	{
	case TMFempty: //no function
		{
			//val = 0;
			break;
		}
	case TMFpolynomial: //polynomial, coeffs contain polynomial coeffs a_1 * 1 + a_2 * t + a_3 * t^2 + ...
		{
			//if (t > 1) t = 1; //hack!!!!

			if (vectime.Length() == 2)
			{
				if (t < vectime(1)) t = vectime(1);
				else if (t > vectime(2)) t = vectime(2);
			}

			double pt = 1;
			if (!diff)
			{
				for (int i=1; i <= coeff.Length(); i++)
				{
					val += coeff(i) * pt;
					pt *= t;
				}
			}
			else if (diff==1) //a_2 * 1 + a_3 * 2* t + a_4 * 3* t^2 + ...
			{
				for (int i=2; i <= coeff.Length(); i++)
				{
					val += coeff(i) * pt * ((double)i-1); 
					pt *= t;
				}
			}
			break;
		}
	case TMFpiecewiseconst: //time points, not interpolated, differentiation not possible
		{
			if (vectime.Length() != 0 && !diff)
			{
				int index = 1;
				index = FindIndexPiecewise(t);
				val = valX(index);
			}
			break;
		}
	case TMFpiecewiselinear: //time points linearly interpolated
		{
			if (vectime.Length() != 0 && !diff)
			{	
				int index = 1;
				index = FindIndexPiecewise(t);
				val = InterpolatePiecewise(t, index);
			}
			break;
		}
	case TMFpiecewisequad: //time points quadratic interpolated, not yet implemented!
		{
			if (vectime.Length() != 0)
			{
				int i = 1;

				while (i <= vectime.Length() && t > vectime(i)) {i++;}

				if (i > vectime.Length()) 
				{
					i = vectime.Length();
					
					if (diff)
						val = valY(i); //vel
					else
						val = valX(i); //pos
				}
				else
				{
					double t1 = 0, t2;
					double x1 = 0;
					double x2 = 0;
					double v1 = 0;
					double v2 = 0;

					if (i >= 2) 
					{
						t1 = vectime(i-1);
						x1 = valX(i-1); //pos
						v1 = valY(i-1); //vel
					}
					double dt = vectime(i)-t1;
					if (dt <= 0) dt = 1;

					x2 = valX(i); //pos
					v2 = valY(i); //vel
					t2 = vectime(i);

					if (diff)
						val = (1.-(t-t1)/dt)*v1 + (1.+(t-t2)/dt)*v2; //linearly interpolate velocity
					else
					{
						val = x1 + (t-t1)*v1 + 0.5*(v2-v1)*Sqr(t-t1)/dt; //linearly interpolate velocity
					}
				}
			}
			break;
		}
	case TMFharmonic: //harmonic function
		{
			if (valX.Length() != 0 && valX.Length() == valY.Length() && valX.Length() == valZ.Length())
			{
				val = 0;
				for (int i=1; i <= valX.Length(); i++)
				{
					if (!diff)
						val += valZ(i)*sin(valX(i)*t+valY(i));
					else
						val += valZ(i)*valX(i)*cos(valX(i)*t+valY(i));
				}
			}
			break;
		}
	case TMFsin: 
		{
			//input variable:t
			//Evaluate val = coeff(1) * Sin(coeff(2)*t + coeff(3))

			val = coeff(1) * sin(coeff(2)*t + coeff(3));

			break;
		}
	case TMFcos: //***
		{
			//input variable:t
			val = coeff(1) * cos(coeff(2)*t + coeff(3));

			break;
		}
	case TMFstaticDynamicFricion: //***
		{
			assert(0 && "ERROR in mathfunc.cpp: SetStaticDynamicFricionFunction not anymore supported by MathFunction, use TMFExpression instead!"); 

			////input variable: velocity
			//if(coeff(3) > 0.0 && fabs(t) <= coeff(3))
			//{
			//	val = coeff(4)/(coeff(3))*t;
			//}else
			//{
			//	//if(-1e-10<t && t<1e-10)t=0.0;
			//	val = coeff(1)*Sgn(t) + coeff(2)*t;
			//}
			//break;
		}
	case TMFuserdefined: //***
		{
			val = ufunc(t,coeff);
			break;
		}
	case TMFexpression:
		{
			//val = mbsPI->EvaluateParsedFunction1D(parsedFunctionExpression, parsedFunctionVariables, t);               //$JG2013-4-29: old slow code of YV
			val = pf.Evaluate(t);                        //$JG2013-4-29: new code, faster: Parsed function is only evaluated, not interpreted!
			break;
		}
	default: ;
	}
	return val;
}

// =====================================================================================================
// PieceWiseMathFunction
// multiple serial MathFunctions define one (discontinous) function
//
//$ DR 2012-06: new class PieceWiseMathFunction added

void PieceWiseMathFunction::CopyFrom(const PieceWiseMathFunction& e)
{
	start_values = e.start_values;
	MathFuncPointArray = e.MathFuncPointArray;
}

void PieceWiseMathFunction::AddMathFunction(MathFunction* mathFpointer, double start_value)
{
	if(start_values.Length()!=0)	// otherwise it is the first on, can be added always
	{
		double max = start_values(start_values.Length());
		if(start_value <= max) 
		{
			(*global_uo) << "Error in PieceWiseMathFunction:: The start_value of the added MathFunction* is smaller than or equal to the maximum already added start_value!\n";
			return;
		}
	}

	MathFunction* mathfun= new MathFunction;
	mathfun = mathFpointer;
	MathFuncPointArray.Add(mathfun);

	start_values.Add(start_value);
}

double PieceWiseMathFunction::Evaluate(double x, int left_limit) const
{
	int i = GetIndexOfMathFunctionPointer(x);
	//(*global_uo) << "PieceWiseMathFunction::GetIndexOfMathFunctionPointer: i = " << i <<" \n";
	if(i>0)
	{
		if(left_limit && i>1)
		{
			return MathFuncPointArray(i-1)->Evaluate(x);
		}
		else
		{
			return MathFuncPointArray(i)->Evaluate(x);
		}
	}
	else
	{
		return -1;
	}
}

int PieceWiseMathFunction::GetIndexOfMathFunctionPointer(double x) const
{
	if(x < start_values(1))
	{
		(*global_uo) << "Error in PieceWiseMathFunction::GetIndexOfMathFunctionPointer(x) x is smaller than the smallest start_value!\n";
		return -1;	//x is smaller than lowest start_value
	}
	else
	{
		for(int i=1; i<= start_values.Length()-1; i++)
		{
			if(x < start_values(i+1))
			{
				return i;
			}
		}
	}
	return start_values.Length();	// the last MathFunction is valid up to infinity
}