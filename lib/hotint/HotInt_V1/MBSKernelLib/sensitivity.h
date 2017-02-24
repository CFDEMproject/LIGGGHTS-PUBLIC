//#**************************************************************
//#
//# filename:             sensitivity.h
//#
//# author:               Rafael Ludwig
//#
//# generated:						February 2011
//# description:          Sensitivity analysis of sensor values from certain model parameters.
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
//#***************************************************************************************
 
#ifndef SENSITIVITY__H
#define SENSITIVITY__H
#include "MBS.h"
#include "performcomputation.h"
class Sensitivity: public PerformComputation
{
public:
	//constructor
	Sensitivity(MultiBodySystem* mbsi):PerformComputation(mbsi), dx(), f1(), f2(), S()
	{
		abs_diff_val = 0.;
		rel_diff_val = 0.;
		type_output = 0;
		type_diff = 0;
		number_of_params = 0;
		sensor_names(0);   // names of parameters for sensitivity analysis
	}
	//destructor
	~Sensitivity()
	{
		for (int i=1; i <= sensor_names.Length(); i++)
		{
			delete sensor_names(i);
			sensor_names(i) = 0;
		}
	}
  virtual mystr GetHeaderString() const; //returns header for output file
	virtual void Init();//initialize sensitivity analysis

	virtual void GetFunctionValues(Vector& f);                 // get function values from sensors
	virtual int PerformSensitivityAnalysis();          // analyze sensitivity of sensor values with respect to parameters
private:
	//MBS* mbs;	// multi body system
	double abs_diff_val; // absolute value "D" for computation of dx=K*x + D (==>f'(x) = df/dx). 
	double rel_diff_val; // relative tolerance for step size
	int type_diff; // -1 ... df/dx =~ (f(x)-f(x-dx))/dx        backward differentiation
	               //  0 ... df/dx =~ (f(x+dx/2)-f(x-dx/2))/dx central differentiation
	               //  1 ... df/dx =~ (f(x+dx)-f(x))/dx        forward differentiation
	int number_of_params; // number of parameters
	// RL 2011-02: moved to class PerformComputation: TArray<mystr*> param_names;   // names of parameters
	TArray<mystr*> sensor_names;  // names of sensors
	int type_output; // 0 ... use final sensor values, 1 ... sensor computation values
	Vector dx;  // step size for differentiation ... dx = x *tol_rel + tol_abs
  //df = f2-f1
	Vector f1;   // function value 1 for difference computation
	Vector f2;   // function value 2 for difference computation
	Matrix S;   // S = dfi/dxj of all computation
	//test: Matrix Smax;   // Smax = dfi/dxj ... maximal sensitivity
	int NCompSensors; // number of sensors with active sensor computation
};
#endif
