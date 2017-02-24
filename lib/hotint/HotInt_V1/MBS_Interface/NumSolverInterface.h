//#***************************************************************************************
//# filename:     NumSolverInterface.h
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

struct NumSolverInterface
{
	virtual int GetNewtonIts() const = 0;
	virtual double NumDiffepsi() const = 0;
	virtual double& NumDiffepsi() = 0;
	virtual int UseSparseSolver() const = 0;
	virtual int& UseSparseSolver() = 0;
	virtual int SymmetricJacobian() const = 0;
};