//#**************************************************************
//#
//# filename:             lapack_routines.cpp
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

#include <assert.h>
#include <memory.h>

#include <math.h>
#include <iostream>

#include "mkl_includes.h"

#include "lapack_routines.h"


#ifndef USE_MKL
extern "C" {
	void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
	
	void dgels_(char *trans, int *m, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, double *work, int *lwork, int *info);

	void dgelsy_(int *m, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, int *jpvt, double *rcond, int *rank, double *work, int *lwork, int *info);

	void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
	
	void dsygv_(int *itype, char *jobz, char* uplo, int *n, double *a, int* lda, double* b, int *ldb, 
		double *w, double *work, int *lwork, int *info);

	void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *info);

	void dggev_(char *jobvl, char* jobvr, int *n, double *a, int *lda, double *b, int *ldb, 
		double *alphar, double *alphai, double* beta, double *vl, int *ldvl, double* vr, int *ldvr, 
		double *work, int *lwork, int *info);

	void dgecon_(char *norm, int *n, double *a, int *lda, double *anorm, double *rcond, double *work, int *iwork, int *info);

	double dlange_(char *norm, int *m, int *n, double *a, int *lda, double *work);

	// $ MSax 2013-07-09 : added solver for A*x=lambda*x; nonsymmetric matrix
	int dgeev_(char *jobvl, char *jobvr, int *n, double *	a, int *lda, double *wr, double *wi, double *vl, 
		int *ldvl, double *vr, int *ldvr, double *work, int *lwork, int *info);
}
#endif

using namespace std;

// solve equation A x = q
// solution vector x is stored in q
// ATTENTION: Lapack matrices are stored columnwise, TRANSPOSED Hotint matrix has to be provided!!!!
int LapackSolve(int n, int nrhs, double*A, double*q)
{
	MKL_INT mkl_n(n);
	MKL_INT mkl_nrhs(nrhs);
	MKL_INT* permutation;
	permutation = new MKL_INT[n];
	MKL_INT info;
	MKL_INT lda = n; // dimesion of q
	MKL_INT ldq = n; // dimension of invA

	dgesv_( &mkl_n, &mkl_nrhs, A, &lda, permutation, q, &ldq, &info );
	delete permutation;
	return info;
}

// estimate the reciprocal of the condition number of a square matrix A
int LapackEstimateReciprocalConditionNumber(int n, double* A, double* rcond)
{
	char cnorm = '1';   // this norm is chosen for calculation of min-max eigenvalues of A^T A

	MKL_INT info;
	MKL_INT mkl_n = n;
	MKL_INT mkl_lda = n;
	
	double* work = new double[4*n];     // (double) work array - for dlange_
	MKL_INT* mkl_ipv = new MKL_INT[n];  // for dgetrf_: (integer) permutation vector, for dgecon_: integer work array
	
	// Computes the norm of x
	double anorm = dlange_(&cnorm, &mkl_n, &mkl_n, A, &mkl_lda, work);

	// Modifies x in place with a LU decomposition
	dgetrf_(&mkl_n, &mkl_n, A, &mkl_lda, mkl_ipv, &info);

	// Computes the reciprocal condition number
	if (info == 0) 
	{
		MKL_INT dummy_info;  //info of dgecon_ is subset of info in dgetrf_, thus redundant.
		dgecon_(&cnorm, &mkl_n, A, &mkl_lda, &anorm, rcond, work, mkl_ipv, &dummy_info);
	}
	else
	{
		rcond = 0;
	}

	delete[] mkl_ipv;
	delete[] work;
	
	return info;
}

// Solves (factorizes) undetermined system (overdetermined: least squares, underdetermined: minimum norm solution) using LAPACK, returns 1 if successfull
// Solution vector x is stored in q, which at input represents the right hand side vector.
// The routine first computes a QR factorization with column pivoting: */
//     A * P = Q * [ R11 R12 ] */
//                 [  0  R22 ] */
// with R11 defined as the largest leading submatrix whose estimated */
// condition number is less than 1/RCOND.  The order of R11, RANK, */
// is the effective rank of A. */
//
// Then, R22 is considered to be negligible, and R12 is annihilated */
// by orthogonal transformations from the right, arriving at the */
// complete orthogonal factorization: */
//    A * P = Q * [ T11 0 ] * Z */
//                [  0  0 ] */
// The minimum-norm solution is then */
//    X = P * Z' [ inv(T11)*Q1'*B ] */
//               [        0       ] */
// where Q1 consists of the first RANK columns of Q. */
int LapackUndeterminedSystemSolve(int m, int n, int nrhs, double* A, double* q, double rcond, int& rank)
{
	// matrix A is m x n is already inserted in column wise ordering, such as is standard in LAPACK (in contrast to Hotint)
	MKL_INT mkl_m(m); //number of rows
	MKL_INT mkl_n(n); //number of columns
	MKL_INT mkl_nrhs(nrhs);  //number of right hand sides (stored in the array q)
	MKL_INT mkl_lda = m; //leading dimesion of A (number of indices to jump from one to the other column of the matrix)
	MKL_INT mkl_ldq = max(m,n); //leading dimension of q (offset, by which one jumps from one to the other solution/right_hand_side)
	
	MKL_INT mkl_rank; //effective rank of A (output: int& rank)
	MKL_INT info;     // return value
	
	MKL_INT *mkl_jpvt = new MKL_INT[n]; //permutation array; initially jpvt[i] must be 0 for all i 
	for(int i=0; i<n; ++i)
	{
		mkl_jpvt[i] = 0;
	}

	//workspace query: determine optimal size of work array in dependence of m,n,nrhs
	//result is written to the first element of the work-array
	//nothing else done here 
	double optimal_work_size;
	MKL_INT mkl_lwork = -1;	
	dgelsy_(&mkl_m, &mkl_n, &mkl_nrhs, A, &mkl_lda, q, &mkl_ldq, mkl_jpvt, &rcond, &mkl_rank, &optimal_work_size, &mkl_lwork, &info);
	mkl_lwork = MKL_INT(optimal_work_size);

	//compute solution
	double* work = new double[mkl_lwork];
	dgelsy_(&mkl_m, &mkl_n, &mkl_nrhs, A, &mkl_lda, q, &mkl_ldq, mkl_jpvt, &rcond, &mkl_rank, work, &mkl_lwork, &info);
	
	delete[] work;
	delete[] mkl_jpvt;

	rank = int(mkl_rank);   //output rank of resulting factorization (rank of R11), which might not be full, i.e., rank < min(m,n).

	return info;
}

// invert Matrix A
// inverse is stored in invA
// solution vector x is stored in q
// Matrix needs not be transposed as in LapackSolve routine
int LapackInvert(int n, double*A, double*invA)
{
	MKL_INT mkl_n(n);
	MKL_INT* permutation;
	permutation = new MKL_INT[n];
	MKL_INT info;
	MKL_INT lda = n; // dimesion of A
	MKL_INT ldinva = n; // dimension of invA
	MKL_INT nrhs = n;  // number of right hand sides

	for (int i=0; i<n*n; i++)
		invA[i] = 0.;
	for (int i=0; i<n; i++)
		invA[i*n+i] = 1.;

	dgesv_( &mkl_n, &nrhs, A, &lda, permutation, invA, &ldinva, &info );
	delete permutation;
	return info;
}

// factorization of Matrix A
// The factorization has the form
//    A = P * L * U
// where P is a permutation matrix, L is lower triangular with unit
// diagonal elements (lower trapezoidal if m > n), and U is upper
// triangular (upper trapezoidal if m < n).
// L and U are stored in A
// Matrix needs not be transposed as in LapackSolve routine
int LapackLU(int m, int n, double*A, int* permutation)
{
	MKL_INT mkl_m(m);
	MKL_INT mkl_n(n);
	MKL_INT info;
	MKL_INT lda(m); // leading dimesion of A
	MKL_INT* mkl_permutation;
	int min_m_n = min(m,n);

	mkl_permutation = new MKL_INT[min_m_n];
	dgetrf_( &mkl_m, &mkl_n, A, &lda, mkl_permutation, &info);

	for (int i=0; i<min_m_n; i++)
	{
		permutation[i] = mkl_permutation[i];
	}

	return info;
}

// solves the generalized eigenvalue problem A u = lambda B u
// eigenvectors are stored ROWWISE in matrix A
// lwork has to be at least 4*n
int LapackGenEVPSPD(int n, double* A, double* B,  double* lami, double* work, int lwork)  
{
	char jobzm = 'V'; // compute eigenmodes also
	char uplo = 'L'; // use upper right triangular matrix (LAPACK matrix is stored column-wise, thus LAPACK 'L' is upper right, 'U' is lower left in Hotint

	if (lwork < 4*n)
	{
		cerr << "length lwork " << lwork << " of array work is to small, needs to be at least 4*n = " << 4*n << "!\n" << flush;
		return(1000);
	}
	MKL_INT mkl_n(n);
	MKL_INT mkl_lwork(lwork);
	MKL_INT info; 
	MKL_INT itype =1; // problem is of form A u = lambda B u
	MKL_INT lda = n;  // size of A
	MKL_INT ldb = n;  // size of B


	dsygv_(&itype,&jobzm,&uplo , &mkl_n , A , &lda, B, &ldb, lami, work, &mkl_lwork, &info); 

	return(info); 
}


// solves the eigenvalue problem A u = lambda u
// eigenvectors are stored ROWWISE in matrix A
// lwork has to be at least 4*n
int LapackEVPSPD(int n, double* A, double* lami, double* work, int lwork)  
{
	char jobz = 'V'; // compute eigenmodes also
	char uplo = 'L'; // use upper right triangular matrix (LAPACK matrix is stored column-wise, thus LAPACK 'L' is upper right, 'U' is lower left in Hotint

	if (lwork < 4*n)
	{
		cerr << "length lwork " << lwork << " of array work is to small, needs to be at least 4*n = " << 4*n << "!\n" << flush;
		return(1000);
	}
	MKL_INT mkl_n(n);
	MKL_INT mkl_lwork(lwork);
	MKL_INT info; 
	MKL_INT lda = n;  // size of A
	MKL_INT ldb = n;  // size of B

	dsyev_(&jobz, &uplo, &mkl_n, A, &lda, lami, work, &mkl_lwork, &info);

	return(info); 
}

// solves the generalized eigenvalue problem A u = lambda B u
// computes right eigen vectors
// lwork has to be at least 8*n+16
// A, B are overwritten
// ev contain right eigenvectors stored row-wise
// eigenvalues are given by (alphar + i alphai) / beta
int LapackGenEVP(int n, double* A, double* B,  double* alphar, double* alphai, double* beta, double* ev, double* work, int lwork)  
{
	char jobvl = 'V'; // right generalized eigenvectors are computed -> Lapack-matrices are stored columnwise, thus lapack l is Hotint right
	char jobvr = 'N'; // left generalized eigenvectors are not computed -> Lapack-matrices are stored columnwise, thus lapack l is Hotint right
	
	if (lwork < 8*n+16)
	{
		cerr << "length lwork " << lwork << " of array work is to small, needs to be at least 8*n+16 = " << 8*n+16 << "!\n" << flush;
		return(1000);
	}

	MKL_INT mkl_n(n);
	MKL_INT mkl_lwork(lwork);
	MKL_INT info; 
	MKL_INT lda = n; // size of A
	MKL_INT ldb = n; // size of B
	MKL_INT ldvl = n; // size of vl, where Hotint-right eigenvectors are stored
	MKL_INT ldvr = 1; // size of vr, where Hotint-left eigenvectors are stored
	double* vl = ev; // eigenvectors from left are not computed
	double* vr = 0;

	dggev_(&jobvl, &jobvr, &mkl_n, A, &lda, B, &ldb, alphar, alphai, beta, vl, &ldvl, vr, &ldvr, work, &mkl_lwork, &info);

	return(info); 
}


int LapackGenEVP_CGEEV(int n, double* A, double* wr, double* wi, double* ev, double* work, int lwork)  
{
	char jobvl = 'V'; // right generalized eigenvectors are computed -> Lapack-matrices are stored columnwise, thus lapack l is Hotint right
	char jobvr = 'N'; // left generalized eigenvectors are not computed -> Lapack-matrices are stored columnwise, thus lapack l is Hotint right
	
	if (lwork < 8*n+16)
	{
		cerr << "length lwork " << lwork << " of array work is to small, needs to be at least 8*n+16 = " << 8*n+16 << "!\n" << flush;
		return(1000);
	}

	MKL_INT mkl_n(n);
	MKL_INT mkl_lwork(lwork);
	MKL_INT info; 
	MKL_INT lda = n; // size of A
	MKL_INT ldvl = n; // size of vl, where Hotint-right eigenvectors are stored
	MKL_INT ldvr = 1; // size of vr, where Hotint-left eigenvectors are stored
	double* vl = ev; // eigenvectors from left are not computed
	double* vr = 0;

	dgeev_(&jobvl, &jobvr, &mkl_n, A, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &mkl_lwork, &info);

	return(info); 
}


