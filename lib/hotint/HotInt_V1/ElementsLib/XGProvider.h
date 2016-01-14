//#**************************************************************
//#
//# filename:             XGProvider.h
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

#pragma once

// this structure and its derived versions switch the behavior of various functions,
// which are implemented by finite elements, between "compute", "draw", displacement, total and initial coordinates;
// the architecture is specialized for the scenario, used by ANCFThinPlate3D element:
// the initial values of the degrees of freedom for the reference configuration are stored separately,
// and XG() of the element operate with the increments of the degrees of freedom (displacements)
// and their time derivatives (velocities);
// such behavior is more consistent for the geometrically linear solutions.

struct XGProvider
{
	void SetXGProvider(Element * element_, const Vector * xgReferenceConfiguration_)
	{
		element = element_;
		xgReferenceConfiguration = xgReferenceConfiguration_;
		assert(xgReferenceConfiguration->Length() == element->FlexDOF());
	}

	// "displacement" of the degree of freedom from the reference configuration
	virtual double XGdispl(int iloc) const = 0;
	// total value of the degree of freedom
	virtual double XGcoord(int iloc) const = 0;
	// time derivative of the degree of freedom
	virtual double XGP(int iloc) const = 0;
	// sometimes it is necessary to know whether this is a "draw" configuration
	virtual bool IsDrawConfiguration() const = 0;

protected:
	Element * element;
	const Vector * xgReferenceConfiguration;			// degrees of freedom in the reference configuration (without time derivatives)
};

// initial state of the element
struct XGProviderInit : public XGProvider
{
	virtual double XGdispl(int iloc) const { return 0; }
	virtual double XGcoord(int iloc) const { return (*xgReferenceConfiguration)(iloc); }
	virtual double XGP(int iloc) const { return element->GetXInit()(iloc + element->SOS()); }
	virtual bool IsDrawConfiguration() const { return false; }
};

// accessing "compute" displacements and coordinates
struct XGProviderCompute : public XGProvider
{
	virtual double XGdispl(int iloc) const { return element->XG(iloc); }
	virtual double XGcoord(int iloc) const { return (*xgReferenceConfiguration)(iloc) + element->XG(iloc); }
	virtual double XGP(int iloc) const { return element->XGP(iloc); }
	virtual bool IsDrawConfiguration() const { return false; }
};

// accessing "draw" displacements and coordinates
struct XGProviderDraw : public XGProvider
{
	virtual double XGdispl(int iloc) const { return element->XGD(iloc); }
	virtual double XGcoord(int iloc) const { return (*xgReferenceConfiguration)(iloc) + element->XGD(iloc); }
	virtual double XGP(int iloc) const { return element->XGPD(iloc); }
	virtual bool IsDrawConfiguration() const { return true; }
};

// "caching" the current state of another XGProvider for the accelerated access;
// only the coordinate values are cached;
// for optimal performance the object should be created on the stack similar to ConstVector;
// the performance should be the same or even better as with the ConstVector,
// as double& Vector::operator() is also virtual
template <int data_size>
struct XGProviderCoordCached : public XGProvider
{
	XGProviderCoordCached(const XGProvider & xg, int length)
	{
		draw = xg.IsDrawConfiguration();
		for(int i = 1; i <= length; i++)
			data[i-1] = xg.XGcoord(i);
	}

	virtual double XGdispl(int iloc) const { assert(0); return 0; }
	virtual double XGcoord(int iloc) const { return data[iloc-1]; }
	virtual double XGP(int iloc) const { assert(0); return 0; }
	virtual bool IsDrawConfiguration() const { return draw; }

protected:
	double data[data_size];
	bool draw;
};

// the same as before, but all data are cached
template <int data_size>
struct XGProviderCached : public XGProvider
{
	XGProviderCached(const XGProvider & xg, int length)
	{
		draw = xg.IsDrawConfiguration();
		for(int i = 1; i <= length; i++)
		{
			dataDispl[i-1] = xg.XGdispl(i);
			dataCoord[i-1] = xg.XGcoord(i);
			dataVel[i-1] = xg.XGP(i);
		}
	}

	virtual double XGdispl(int iloc) const { return dataDispl[iloc-1]; }
	virtual double XGcoord(int iloc) const { return dataCoord[iloc-1]; }
	virtual double XGP(int iloc) const { return dataVel[iloc-1]; }
	virtual bool IsDrawConfiguration() const { return draw; }

protected:
	double dataDispl[data_size];
	double dataCoord[data_size];
	double dataVel[data_size];
	bool draw;
};