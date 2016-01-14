//#**************************************************************
//#
//# filename:             mathfunc.h
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						October 2006
//# description:          general mathematical function
//#                       mostly for time or space-curves of arbitrary shape
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

#ifndef MATHFUNC__H 
#define MATHFUNC__H

typedef enum {TMFempty=0, TMFpiecewiseconst=1, TMFpiecewiselinear=2, 
TMFpiecewisequad=3, TMFuserdefined=4, TMFexpression=5, 
//the following modes are depreciated! do not use, use TMFexpression instead!!!
TMFpolynomial=6, TMFharmonic=7, TMFsin=8, TMFcos=9, TMFstaticDynamicFricion=10} TMathFuncType; //***

#include "../Parser/parser.h"
//class CParsedFunction;

//general function that evaluates for certain parameter x or t
class MathFunction
{
public:
	MathFunction() {Init();}

	MathFunction(const MathFunction& e)
	{
		Init();
		CopyFrom(e);
	}
	MathFunction& operator=(const MathFunction& e) 
	{
		if (this == &e) {return *this;}
		CopyFrom(e);
		return *this;
	}
	//To be overwritten in derived class:
	virtual MathFunction* GetCopy()
	{
		MathFunction* ec = new MathFunction();
		ec->CopyFrom(*this);
		return ec;
	}
	virtual ~MathFunction() { };

	//To be overwritten in derived class:
	virtual void CopyFrom(const MathFunction& e);

	virtual void Init() {dim = 1; funcmode = TMFempty; vectime.SetLen(0); valX.SetLen(0); valY.SetLen(0); valZ.SetLen(0); coeff.SetLen(0); ufunc=NULL;}

	virtual const char* GetTypeName() const;

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(MBS* mbs, ElementDataContainer& edc); //set element data according to ElementDataContainer	

	virtual int SetData(int mode, const Matrix& data);

	virtual void GetData(int& mode, Matrix& data) const;

	virtual int SetData(const mystr& funcname, const Vector& coefficients);

	virtual void GetData(mystr& funcname, Vector& coefficients) const;

	virtual TMathFuncType StringToMathFuncType(const mystr& funcname);

	virtual void SetPiecewise(const Vector& times, const Vector& coeffs, int interp = 0);

	virtual void SetPiecewiseQuadratic(const Vector& times, const Vector& coeffs_p, const Vector& coeffs_v);

	virtual int SetPiecewiseFromFile(const char* filename, int ncolumns, int column1, int column2, int interp = 0,  int nrOfHeaderLines = 0, double offset1 = 0., double offset2 = 0., int offset1_start_index=1, int offset2_start_index=1);	// read x- and y- data from columns col1 and col2 of file, error if return value = 1, specified number of header lines can be cutted off, additional time or signal offset is considered during reading from file 
	
	virtual int SetPiecewiseFromFile2(const char* filename, int column1, int column2, int interp = 0, mystr comment = mystr("%"));	// read x- and y- data from columns col1 and col2 of file, error if return value = 1, number of columns is not needed, comment char for cutting off header available


	//f = c1*t^0 + c2*t^1 + c3*t^2 + ....
	virtual void SetPolynomial(const Vector& coeffs);
	virtual void SetVecTimes(const Vector& times); //for polynomial turn on/off times

	virtual void SetConstant(double x);
  
	virtual void SetHarmonic(const Vector& frequency, const Vector& phase, const Vector& amplitude);

	virtual int GetFuncMode() const {return funcmode;}
	virtual int GetMaxFuncMode() const {return 5;} //maximum number of usable functions!

	virtual int Dim() const {return dim;}
	//return 0 or 1 depending if condition f >= 1 is fulfilled:
	virtual int CheckCondition(double t) const {return (Evaluate(t) >= 1); }
	//evaluate given function
	virtual double Evaluate(double t, int diff=0) const;

	//virtual void SetSinFunction(double amplitude, double freq_rad, double phase); 
	//virtual void SetCosFunction(double amplitude, double freq_rad, double phase); //***
	//virtual void SetStaticDynamicFricionFunction(double staticFriction, double viscousFriction, double zeroZone = 0); //***
  virtual void SetUserDefined(double (*func)(double,const Vector&),const Vector& coeffs = Vector(0)); // sets a userdefined function, additional parameters cam be passed in Vector coeffs
  virtual void SetExpression(mystr expression, mystr variable, MBS * mbs); // sets a userdefined expression - uses parseobect for evaluation
               
	virtual const Vector& GetCoeffVector() const {return coeff;}
	virtual Vector& GetCoeffVector() {return coeff;}
	
	virtual const Vector& GetTimeVector() const {return vectime;}
	virtual Vector& GetTimeVector() {return vectime;}

	virtual const Vector& GetXVector() const {return valX;}
	virtual Vector& GetXVector() {return valX;}

	//$ YV 2012-12-10
	/*
	virtual const CParsedFunction& PF() const {return pf;}
	virtual CParsedFunction& PF(){return pf;}
	*/
	virtual CParsedFunction* PF() {return &pf;} //$JG2013-04-29

	virtual int FindIndexPiecewise(double t) const;  //$ MS 2011-5-23
	virtual double InterpolatePiecewise(double t, int index) const; //$ MS 2011-5-23: 

private:
	int funcmode;					//see Evaluate function
	int dim;							//dimension of return value
	Vector vectime;				//vector of time points
	Vector valX;					//value X for time point vectime
	Vector valY;					//value Y for time point vectime
	Vector valZ;					//value Z for time point vectime

	Vector coeff;					//coefficients for mathematical operations - (AD) additional parameters for userdefined function
	double (*ufunc)(double,const Vector&);			// pointer to a User-defined function double = f(double,Vector), eg. f( time, other_coefficients)
	
	CParsedFunction pf; //$JG2013-4-29: for acceleration of parsed function!!!

	//$JG2013-4-29: shoud be erased soon:
	//++++++++++++++++
	//$ YV 2012-12-12: we store definition of a parsed function instead of CParsedFunction itself
	mystr parsedFunctionExpression;
	mystr parsedFunctionVariables;
	//++++++++++++++++
	MBSParserInterface * mbsPI; //$JG2013-4-29: used to communicate with MBSParser (for translation)
};

// =====================================================================================================
// PieceWiseMathFunction
// multiple serial MathFunctions define one (discontinous) function
//
//$ DR 2012-06: new class PieceWiseMathFunction added
class PieceWiseMathFunction
{
public:
		PieceWiseMathFunction() {Init();}

	PieceWiseMathFunction(const PieceWiseMathFunction& e)
	{
		Init();
		CopyFrom(e);
	}
	PieceWiseMathFunction& operator=(const PieceWiseMathFunction& e) 
	{
		if (this == &e) {return *this;}
		CopyFrom(e);
		return *this;
	}
	//To be overwritten in derived class:
	virtual PieceWiseMathFunction* GetCopy()
	{
		PieceWiseMathFunction* ec = new PieceWiseMathFunction();
		ec->CopyFrom(*this);
		return ec;
	}
	virtual ~PieceWiseMathFunction() { };

	//To be overwritten in derived class:
	virtual void CopyFrom(const PieceWiseMathFunction& e);

	virtual void Init() {};

	virtual void AddMathFunction(MathFunction* mathFpointer, double start_value);
	virtual double Evaluate(double x, int left_limit=0) const;

	virtual int GetIndexOfMathFunctionPointer(double x) const;

private:
	TArray<double> start_values;				//MathFunction i is valid for points with x >= start_value(i)
	TArray<MathFunction*> MathFuncPointArray;
};


#endif