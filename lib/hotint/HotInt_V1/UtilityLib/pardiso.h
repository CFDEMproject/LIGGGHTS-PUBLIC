//#**************************************************************
//#
//# filename:             pardiso.h
//#
//# author:               Astrid Pechstein
//#
//# generated:          
//# description:          
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

#ifndef PARDISO_INCLUDES_H
#define PARIDSO_INCLUDES_H

#include "mkl_includes.h"

class PardisoInverse
{
public:
	PardisoInverse() : sparseA(0), size(0), mtype(11), solver(0), error(0), perm(0), isFirstStep(1), isFactorized(0)
	{
	}

	PardisoInverse:: PardisoInverse(const SparseMatrix & matA) : sparseA(&matA), size(0), mtype(11), solver(0), error(0), perm(0), isFirstStep(1), isFactorized(0)
	{
		sparseA =&matA;
		size = matA.Getcols();
		SetDefaultIParm();

		perm = new _INTEGER_t[size];
	}


	~PardisoInverse();

	void SetDefaultIParm();

	void SetSparseMatrix(const SparseMatrix& matA)
	{
		sparseA =&matA;
		size = matA.Getcols();

		SetDefaultIParm();
		delete [] perm;

		perm = new _INTEGER_t[size];
	}

	int Solve(Vector& q);
	int Factorize();
	int Apply(Vector& q);

	int& IsFirstStep() {return isFirstStep; }
	const int& IsFirstStep() const {return isFirstStep; }


private:
	// pointer to original sparsematrix
	const SparseMatrix* sparseA;

	// size, has to be quadratic
	_INTEGER_t size;
	// Pardiso-Data
	_MKL_DSS_HANDLE_t pt[64];

	_INTEGER_t mtype;
	int solver;
	_INTEGER_t iparm[64];

	double dparm[64];

	_INTEGER_t error;

	_INTEGER_t* perm;

	int isFirstStep;
	int isFactorized;

};


#endif