//#***************************************************************************************
//# filename:     useroutputinterface.h
//#
//# author:				Johannes Gerstmayr, Yuri Vetyukov
//# 
//# generated:      
//# description:  
//#                       
//# comments:      
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

#include "elementdata.h"

typedef enum { UO_LVL_0 = 0,    // no Output
							 UO_LVL_err = 1,  // necessary output (Errors, start/end simulation)
							 UO_LVL_warn = 2, // almost necessary output (Warnings)
							 UO_LVL_multsim = 3, // multiple simulation (parameter variation/optimization)
							 UO_LVL_sim = 4,  // simulation output (solver)
							 UO_LVL_ext = 5,  // extended output (useful information)
							 UO_LVL_all = 6,  // complete information
							 UO_LVL_dbg1 = 7, // debug level 1
							 UO_LVL_dbg2 = 8,  // debug level 2
							 UO_LVL_max = 1000 //this should always be the maximum ==> used to have output in any case!
						 } UO_MSGLVL;

class Vector3D;
class Vector2D;
class Vector;
class SparseVector;
class Matrix3D;
class Matrix;
class SparseMatrix;

struct UserOutputInterface
{
	virtual UserOutputInterface & operator <<(const char * pStr) = 0;
	virtual UserOutputInterface & operator <<(int x) = 0;
	virtual UserOutputInterface & operator <<(double x) = 0;
	virtual UserOutputInterface & operator <<(const Vector3D & v) = 0;
	virtual UserOutputInterface & operator <<(const Vector2D & v) = 0;
	virtual UserOutputInterface & operator <<(const Vector & v) = 0;
	virtual UserOutputInterface & operator <<(const IVector & v) = 0;
	virtual UserOutputInterface & operator <<(const TArray<double> & v) = 0;
	virtual UserOutputInterface & operator <<(const SparseVector & v) = 0;
	virtual UserOutputInterface & operator <<(const Matrix3D & m) = 0;
	virtual UserOutputInterface & operator <<(const Matrix & m) = 0;
	virtual UserOutputInterface & operator <<(const SparseMatrix & m) = 0;

	virtual void InstantMessageText(const char* pStr) = 0;
	virtual int GetGlobalMessageLevel() = 0;
	// cms elements wish to use this
	virtual int CallWCDriverFunction(int action, int option = 0, int value = 0, ElementDataContainer* edc = NULL) = 0;
	virtual void SaveLocalMessageLevel() = 0;
	virtual void SetLocalMessageLevel(UO_MSGLVL message_level) = 0;
	virtual void ResetLocalMessageLevel() = 0;
	virtual int PrintMsg() = 0;
};