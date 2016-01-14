//#**************************************************************
//#
//# filename:             fft_utilities.h
//#
//# author:               Rafael Ludwig
//#
//# generated:						April 2011
//# description:          
//#                       
//# remarks:						  Include this file for computations of FFT in HOTINT.
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

#ifndef FFT_UTILITIES__H
#define FFT_UTILITIES__H

// generate vectors with constant step-size dtIn; only time points are considered which are within the interval [tStart,tEnd], interpolation_order = 0|(1) ... constant/linear interpolation
int constTimeStep(Vector& times, Vector& values, double dtIn, double tStart =-1.e30, double tEnd = 1.e30, int interpolation_order = 1);

// input:   sampletime  ... constant sample time of vector sampletime
//          samples     ... vector with constant sample time
// output:  frequencies ... vector with frequencies
//          amplitudes  ... vector with amplitudes corresponding to frequencies
//          phase (rad) ... vector with phase corresponding to frequencies
void makefft(double sampletime, const Vector& samples, Vector& frequencies, Vector& amplitudes, Vector& phase);


 

#endif

//define useful FFT-functions here
