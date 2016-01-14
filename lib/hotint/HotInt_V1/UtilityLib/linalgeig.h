//#**************************************************************
//#
//# filename:             linalgeig.h
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						April 2006
//# description:          Eigenvalue solver from TNT library
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

#ifndef LINALGEIG__H
#define LINALGEIG__H


class MultiBodySystem;
/** 

Computes eigenvalues and eigenvectors of a real (non-complex)
matrix. 
<P>
If A is symmetric, then A = V*D*V' where the eigenvalue matrix D is
diagonal and the eigenvector matrix V is orthogonal. That is,
the diagonal values of D are the eigenvalues, and
V*V' = I, where I is the identity matrix.  The columns of V 
represent the eigenvectors in the sense that A*V = V*D.

<P>
If A is not symmetric, then the eigenvalue matrix D is block diagonal
with the real eigenvalues in 1-by-1 blocks and any complex eigenvalues,
a + i*b, in 2-by-2 blocks, [a, b; -b, a].  That is, if the complex
eigenvalues look like
<pre>

u + iv     .        .          .      .    .
.      u - iv     .          .      .    .
.        .      a + ib       .      .    .
.        .        .        a - ib   .    .
.        .        .          .      x    .
.        .        .          .      .    y
</pre>
then D looks like
<pre>

u        v        .          .      .    .
-v        u        .          .      .    . 
.        .        a          b      .    .
.        .       -b          a      .    .
.        .        .          .      x    .
.        .        .          .      .    y
</pre>
This keeps V a real matrix in both symmetric and non-symmetric
cases, and A*V = V*D.



<p>
The matrix V may be badly
conditioned, or even singular, so the validity of the equation
A = V*D*inverse(V) depends upon the condition number of V.

<p>
(Adapted from JAMA, a Java Matrix Library, developed by jointly 
by the Mathworks and NIST; see  http://math.nist.gov/javanumerics/jama).
**/


class Eigenvalue
{
public:
	// Check for symmetry, then construct the eigenvalue decomposition
	// @param A    Square double (non-complex) matrix
	Eigenvalue(const Matrix &A);

	// Return the eigenvector matrix
	// @return     V
	const Matrix& GetV() 
	{
		return V;
	}

	// Return the double parts of the eigenvalues
	// @return     double(diag(D))
	const Vector& GetRealEigenvalues() 
	{
		return d;
	}

	// Return the imaginary parts of the eigenvalues in parameter e_.
	// @pararm e_: new matrix with imaginary parts of the eigenvalues.
	const Vector& GetImagEigenvalues() 
	{
		return e;
	}


	/** 
	Computes the block diagonal eigenvalue matrix.
	If the original matrix A is not symmetric, then the eigenvalue 
	matrix D is block diagonal with the real eigenvalues in 1-by-1 
	blocks and any complex eigenvalues,
	a + i*b, in 2-by-2 blocks, [a, b; -b, a].  That is, if the complex
	eigenvalues look like
	<pre>

	u + iv     .        .          .      .    .
	.      u - iv     .          .      .    .
	.        .      a + ib       .      .    .
	.        .        .        a - ib   .    .
	.        .        .          .      x    .
	.        .        .          .      .    y
	</pre>
	then D looks like
	<pre>

	u        v        .          .      .    .
	-v        u        .          .      .    . 
	.        .        a          b      .    .
	.        .       -b          a      .    .
	.        .        .          .      x    .
	.        .        .          .      .    y
	</pre>
	This keeps V a double matrix in both symmetric and non-symmetric
	cases, and A*V = V*D.

	@param D: upon return, the matrix is filled with the block diagonal 
	eigenvalue matrix.

	*/

	void GetD(Matrix &D);

private:

	/** Row and column dimension (square matrix).  */
	int n;

	int issymmetric; /* boolean*/

	/** Arrays for internal storage of eigenvalues. */

	Vector d;         /* double part */
	Vector e;         /* img part */

	/** Array for internal storage of eigenvectors. */
	Matrix V;

	/** Array for internal storage of nonsymmetric Hessenberg form.
	@serial internal storage of nonsymmetric Hessenberg form.
	*/
	Matrix H;


	/** Working storage for nonsymmetric algorithm.
	@serial working storage for nonsymmetric algorithm.
	*/
	Vector ort;


	// Symmetric Householder reduction to tridiagonal form.

	void Tred2();

	// Symmetric tridiagonal QL algorithm.
	void TQL2();

	// Nonsymmetric reduction to Hessenberg form.
	void OrtHes ();

	// Complex scalar division.

	void Cdiv(double xr, double xi, double yr, double yi, double& cdivr, double& cdivi);

	// Nonsymmetric reduction from Hessenberg to double Schur form.
	void HQR2 ();



};

struct MBS;

// Solves the generalized Eigenvalue problem
//
//    A u = lambda M u
//
// using EITHER 
// the LOBPCG I method by Knyazev 
// (see paper Toward the optimal preconditioned Eigensolver, Algorithm 5.1)
// OR
// Matlabs eigs for sparse Matrices, Arnoldi iteration
// Matrix A positive semidefinite
// Matrix M positive definite, 
//   if M is Null-pointer, identity matrix is used
//   A has to be provided
// computed are the smallest eigenvalues and corresponding modes
class SparseEigenvalueSolver
{
public:
	SparseEigenvalueSolver() : A(0), M(0), maxsteps(0), precision(0), mbs(0), use_preconditioner(0), lambda_precond(0), can_be_stopped(0)
	{ ; }

	SparseEigenvalueSolver(MBS* mbsi) : A(0), M(0), maxsteps(0), precision(0), mbs(mbsi), use_preconditioner(0), lambda_precond(0), can_be_stopped(0)
	{ ; }

	void Set(SparseMatrix* matA, SparseMatrix* matM, int maxstepsI, double precisionI, int init_random)
	{
		A = matA; M = matM;
		maxsteps = maxstepsI;
		precision = precisionI;
		initialize_random = init_random;
	}

	void SetPreconditioner(double lambda)
	{
		use_preconditioner = 1;
		lambda_precond = lambda;
	}

	void SetMatrices(SparseMatrix* matA, SparseMatrix* matM)
	{
		A = matA; M = matM;
	}
	int& MaxSteps() {return maxsteps;}
	double& Precision() {return precision;}
	int& InitializeRandom() {return initialize_random;}


	void ApplyMat(SparseMatrix* mat, Vector& u, Vector& Mu)
	{
		if (mat)
			Mult(*mat, u, Mu);
		else
			Mu = u;
	}

	// compute nev smallest eigenvalues of the generalized eigenvalue system 
	// A u = lambda M u
	// by LOBPCG method (Knyazev)
	// Vector eigenvalues of length nev
	//   and
	// Array<Vector*> of length nev, vectors of matrix dimensions
	// have to be provided!!
	// returns 1 if converged correctly, 0 else
	int ComputeEigenModes(int nev, Vector& eigenvalues, Matrix& eigmodes, TArray<int>& unconstraineddofs, int nzeromodes = 0);

	// compute generalized eigenmodes in Matlab
	// neig smallest eigenvalues are computed iteratively
	// eigenmodes are stored ROWwise in Matrix eigmodes
	// eigval and eigmodes have to be provided with appropriate size!
	// returns 1 if converged correctly, 0 else
	int ComputeEigenModesMatlab(const mystr& pathMatlab, int neig, 
		Vector& eigval, Matrix& eigmodes, TArray<int>& unconstraineddofs,
		int nzeromodes=0);

	int& CanBeStopped() {return can_be_stopped;}
	const int& CanBeStopped() const {return can_be_stopped;}

private:
	MBS* mbs;							// pointer to Multibody system, for output..
	SparseMatrix* A;     // operator A, usually stiffness matrix, positive definite?
	SparseMatrix* M;     // operator M, usually mass matrix, positive semidefinite

	int maxsteps;
	double precision;
	int initialize_random;

	int use_preconditioner;  // use preconditioner: (A + lambda M)^-1
	double lambda_precond;   // lambda for preconditioner

	int can_be_stopped;  // eigenmode computation can be stopped via MBS

};

#endif
