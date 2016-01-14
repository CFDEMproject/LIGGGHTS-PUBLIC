//#**************************************************************
//#
//# filename:             lapack_routines.h
//#
//# author:               Astrid Pechstein
//#
//# generated:						August 2010
//# description:          Interface to LAPACK routines
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
//#**************************************************************

#ifndef	LAPACK_ROUTINES_HPP
#define LAPACK_ROUTINES_HPP



// solve equation A x = q
// solution vector x is stored in q
// ATTENTION: Lapack matrices are stored columnwise, TRANSPOSED Hotint matrix has to be provided!!!!
int LapackSolve(int n, int nrhs, double*A, double*q);

// estimate the reciprocal of the condition number of a square matrix A
int LapackEstimateReciprocalConditionNumber(int n, double* A, double* rcond);

// invert Matrix A
// inverse is stored in invA
// solution vector x is stored in q
// Matrix needs not be transposed as in LapackSolve routine
int LapackInvert(int n, double*A, double*invA);

// factorization of Matrix A
// The factorization has the form
//    A = P * L * U
// where P is a permutation matrix, L is lower triangular with unit
// diagonal elements (lower trapezoidal if m > n), and U is upper
// triangular (upper trapezoidal if m < n).
// L and U are stored in A
// Matrix needs not be transposed as in LapackSolve routine
int LapackLU(int m, int n, double*A, int* permutation);

// solve over- or underdetermined system A x = q with m times n matrix A
// solution vector x is stored in q.
// DGELSY computes the minimum-norm solution to a real linear least
// squares problem:
//     minimize || A * X - B ||
// using a complete orthogonal factorization of A.  A is an M-by-N
// matrix which may be rank-deficient.
//
// Several right hand side vectors b and solution vectors x can be
// handled in a single call; they are stored as the columns of the
// M-by-NRHS right hand side matrix B and the N-by-NRHS solution
// matrix X.
//
// The routine first computes a QR factorization with column pivoting:
//     A * P = Q * [ R11 R12 ]
//                 [  0  R22 ]
// with R11 defined as the largest leading submatrix whose estimated
// condition number is less than 1/RCOND.  The order of R11, RANK,
// is the effective rank of A.
//
// Then, R22 is considered to be negligible, and R12 is annihilated
// by orthogonal transformations from the right, arriving at the
// complete orthogonal factorization:
//    A * P = Q * [ T11 0 ] * Z
//                [  0  0 ]
// The minimum-norm solution is then
//    X = P * Z**T [ inv(T11)*Q1**T*B ]
//                 [        0         ]
// where Q1 consists of the first RANK columns of Q.
int LapackUndeterminedSystemSolve(int m, int n, int nrhs, double* A, double* q, double rcond, int& rank);

// solves the eigenvalue problem A u = lambda u
// eigenvectors are stored ROWWISE in matrix A
// lwork has to be at least 4*n
// A has to be symmetric, upper right triangle has to be provided (rest optional)
int LapackEVPSPD(int n, double* A, double* lami, double* work, int lwork);

// solves the generalized eigenvalue problem A u = lambda B u
// eigenvectors are stored ROWWISE in matrix A
// lwork has to be at least 4*n
// A, B have to be symmetric, upper right triangle has to be provided (rest optional)
int LapackGenEVPSPD(int n, double* A, double* B,  double* lami, double* work, int lwork)  ;

// solves the generalized eigenvalue problem A u = lambda B u
// lwork has to be at least 8*n+16
// A, B are overwritten
// ev contain right eigenvectors stored row-wise
// eigenvalues are given by (alphar + i alphai) / beta
int LapackGenEVP(int n, double* A, double* B,  double* alphar, double* alphai, double* beta, double* ev, double* work, int lwork)  ;

// $ MSax 2013-07-09 : added
// solves the generalized eigenvalue problem A u = lambda u
// lwork has to be at least 8*n+16
// A is overwritten
// ev contain right eigenvectors stored row-wise
// eigenvalues are given by wr + i wi
int LapackGenEVP_CGEEV(int n, double* A, double* wr, double* wi, double* ev, double* work, int lwork)  ; 	// $ MSax 2013-07-09 : added solver for A*x=lambda*x; nonsymmetric matrix


#endif // LAPACK_ROUTINES_HPP
