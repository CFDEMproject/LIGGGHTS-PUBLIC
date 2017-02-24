//#**************************************************************
//#
//# filename:             makefft.cpp
//#
//# author:               Rafael Ludwig
//#
//# generated:						April 2011
//# description:          
//#                       
//# remarks:						  This file constains the algorithms defined in makefft.h. 
//#                       
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

#include "element.h"


//please move this function to some appropriate place in math/utilities !!!! //JG2013-01-07
// generate vectors with constant step-size dtIn; only time points are considered which are within the interval [tStart,tEnd]
int constTimeStep(Vector& times, Vector& values, double dtIn, double tStart, double tEnd, int interpolation_order)
{
	if(tStart >= tEnd || times.Length()!=values.Length() || times.Length() < 2){assert(0 && "Problems during linear interpolation of t-y data to constant time step.");}
	MathFunction mf;
	mf.SetPiecewise(times, values,interpolation_order);

	// start time point
	double t = tStart;
	if(t<times(1))
	{
		t = times(1);
	}	
	TArray<double> tt(1000); 
	TArray<double> tval(1000); 
	while(t <= tEnd && t <= times.Get(times.GetLen()))
	{
		tt.Add(t);
		tval.Add(mf.Evaluate(t));
		t += dtIn;
	}
	times = Vector(tt);
	values = Vector(tval);
	return 1;
}


void makefft(const double sampletime, const Vector& samples, Vector& frequencies, Vector& amplitudes, Vector& phase)
{
  assert(0 && "FFT is not available at the moment!");
}

