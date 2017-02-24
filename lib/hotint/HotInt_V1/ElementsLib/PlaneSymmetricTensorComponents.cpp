//#**************************************************************
//# filename:             PlaneSymmetricTensorComponents.cpp
//#
//# author:               Vetyukov Yury
//#
//# generated:						February 2011
//# description:          Utility class
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
 
#include "PlaneSymmetricTensorComponents.h"

void PlaneSymmetricTensorComponents::SetComponents(double t11, double t12, double t22)
{
	this->t11 = t11;
	this->t12 = t12;
	this->t22 = t22;
}

void PlaneSymmetricTensorComponents::SetComponent(int alphaBeta, double value)
{
	switch(alphaBeta)
	{
	case 1: t11 = value; break;
	case 2: t12 = value; break;
	case 3: t22 = value; break;
	}
}

PlaneSymmetricTensorComponents & PlaneSymmetricTensorComponents::operator*=(double k)
{
	t11 *= k;
	t12 *= k;
	t22 *= k;
	return *this;
}

PlaneSymmetricTensorComponents & PlaneSymmetricTensorComponents::operator+=(const PlaneSymmetricTensorComponents & T)
{
	t11 += T.t11;
	t12 += T.t12;
	t22 += T.t22;
	return *this;
}

PlaneSymmetricTensorComponents & PlaneSymmetricTensorComponents::operator-=(const PlaneSymmetricTensorComponents & T)
{
	t11 -= T.t11;
	t12 -= T.t12;
	t22 -= T.t22;
	return *this;
}

double PlaneSymmetricTensorComponents::Convolute(
			const PlaneSymmetricTensorComponents & T,double C1,double C2,
			const PlaneSymmetricTensorComponents & A)
{
	return A.t11*A.t11*(C1 + C2)*t11*T.t11 + A.t22*A.t22*(C1 + C2)*t22*T.t22 + 
		A.t12*A.t12*(C2*t22*T.t11 + 4*C1*t12*T.t12 + 2*C2*t12*T.t12 + C2*t11*T.t22) + 
		2*A.t12*A.t22*(C1 + C2)*(t22*T.t12 + t12*T.t22) + 
		A.t11*(2*A.t12*(C1 + C2)*(t12*T.t11 + t11*T.t12) + 
		A.t22*(C1*t22*T.t11 + 2*C2*t12*T.t12 + C1*t11*T.t22));
}

double PlaneSymmetricTensorComponents::SelfConvolute(
			double C1,double C2,
			const PlaneSymmetricTensorComponents & A) const
{
	return A.t11*A.t11*(C1 + C2)*t11*t11 + 4*A.t12*A.t22*(C1 + C2)*t12*t22 + 
		A.t22*A.t22*(C1 + C2)*t22*t22 + 
		2*A.t11*(t12*(2*A.t12*(C1 + C2)*t11 + A.t22*C2*t12) + A.t22*C1*t11*t22) + 
		2*A.t12*A.t12*((2*C1 + C2)*t12*t12 + C2*t11*t22);
}

double PlaneSymmetricTensorComponents::Trace(const PlaneSymmetricTensorComponents & A)
{
    return t11 * A.t11 + 2 * t12 * A.t12 + t22 * A.t22;
}

void PlaneSymmetricTensorComponents::BuildDeviator(const PlaneSymmetricTensorComponents & Acovariant, const PlaneSymmetricTensorComponents & Acontravariant)
{
	PlaneSymmetricTensorComponents A = Acovariant;
	A *= Trace(Acontravariant) / 3.;
	*this -= A;
}

double PlaneSymmetricTensorComponents::Det() const
{
    return t11 * t22 - t12 * t12;
}

PlaneSymmetricTensorComponents PlaneSymmetricTensorComponents::Inverse() const
{
    double det = Det();
    return PlaneSymmetricTensorComponents(t22 / det, -t12 / det, t11 / det);
}