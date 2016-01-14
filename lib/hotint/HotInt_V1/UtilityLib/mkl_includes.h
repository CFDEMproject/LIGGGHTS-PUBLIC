//#**************************************************************
//#
//# filename:             mkl_includes.cpp
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

#ifndef LAPACK_INCLUDES_H
#define LAPACK_INCLUDES_H

#ifdef USE_MKL
#include <mkl_lapack.h>
#include <mkl_pardiso.h>
#else
typedef int MKL_INT;
typedef int _INTEGER_t;
typedef int _MKL_DSS_HANDLE_t;
#endif

#endif // LAPACK_INCLUDES_H