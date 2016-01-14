
/*
* -- SuperLU routine (version 2.0) --
* Univ. of California Berkeley, Xerox Palo Alto Research Center,
* and Lawrence Berkeley National Lab.
* November 15, 1997
*
*/

//#include "stdafx.h"

#include "slu_ddefs.h"
#include "SuperLUmain.h"


#include "ioincludes.h"


ofstream slu_outfile("superLUout.txt");

void WriteSLUErrorMessage(char* str)
{
	slu_outfile << "Super LU message='" << str << "'\n" << flush;
}

mystr systemstr;


mystr& GetSuperLUMessage()
{
	return systemstr;
}

void PrintSuperMatrix(const SuperMatrix *A)
{
	if (A->Stype == SLU_SC)
	{
		SCformat* sc = (SCformat*)(A->Store);
		systemstr += "nnz = ";
		systemstr += sc->nnz;
		systemstr += "\nVals = (";
		for ( int i = 0; i < sc->nnz; i++)
		{
			double* vals = (double*)sc->nzval;
			systemstr += vals[i];
			systemstr += ", ";
		}
		systemstr += ")\n";
		//for(int i = 0; i < 5; i++)
		//{
		//	systemstr += sc->rowind[i];
		//	systemstr += ", ";
		//}
	}
	else if (A->Stype == SLU_NR)
	{
		NRformat* sc = (NRformat*)(A->Store);
		systemstr += "nnz = ";
		systemstr += sc->nnz;
		systemstr += "\nVals = (";
		for ( int i = 0; i < sc->nnz; i++)
		{
			double* vals = (double*)sc->nzval;
			systemstr += vals[i];
			systemstr += ", ";
		}
		systemstr += ")\n";
	}
	else if (A->Stype == SLU_NC)
	{
		NCformat* sc = (NCformat*)(A->Store);
		systemstr += "nnz = ";
		systemstr += sc->nnz;
		systemstr += "\nVals = (";
		for ( int i = 0; i < sc->nnz; i++)
		{
			float* vals = (float*)sc->nzval;
			double a = vals[i];
			char str[50];
			sprintf_s(str, 50, "%g", a);
			systemstr += str;
			systemstr += ", ";
		}
		systemstr += ")\n";
		//for(int i = 0; i < 5; i++)
		//{
		//	systemstr += sc->rowind[i];
		//	systemstr += ", ";
		//}
	}


}

void PrintIntVec(int* vec, int n)
{
	systemstr += "(";
	for ( int i=0; i<n; i++)
	{
		systemstr += vec[i];
		systemstr += ", ";
	}
	systemstr += ")\n";
}

int CallSuperLUSparseCol(int nrows, int ncols, int nvals, int* col_ptr, int* row_ind, double* vals, double* rhs)
{
	SuperMatrix A, L, U, B;
	//int      *asub, *xa;
	int      *perm_r; /* row permutations from partial pivoting */
	int      *perm_c; /* column permutation vector */
	int      nrhs, info;
	superlu_options_t options;
	SuperLUStat_t stat;

	nrhs = 1; //number of right-hand sides

	/* Create matrix A in the format expected by SuperLU. */
	dCreate_CompCol_Matrix(&A, nrows, ncols, nvals, vals, row_ind, col_ptr, SLU_NC, SLU_D, SLU_GE);

	/* Create right-hand side matrix B. */
	//if ( !(rhs = doubleMalloc(m * nrhs)) ) ABORT("Malloc fails for rhs[].");
	//for (i = 0; i < m; ++i) rhs[i] = 1.0;

	dCreate_Dense_Matrix(&B, ncols, nrhs, rhs, ncols, SLU_DN, SLU_D, SLU_GE);

	if ( !(perm_r = intMalloc(ncols)) ) ABORT("Malloc fails for perm_r[].");
	if ( !(perm_c = intMalloc(nrows)) ) ABORT("Malloc fails for perm_c[].");

	/* Set the default input options. */
	set_default_options(&options);
	options.ColPerm = NATURAL;
	//options.Trans = TRANS;

	/* Initialize the statistics variables. */
	StatInit(&stat);

	dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

	DNformat     *Bstore;
	double       *dp;

	Bstore = (DNformat *) B.Store;
	dp = (double *) Bstore->nzval;

	for (int i = 0; i < nrows; i++) 
	{
		rhs[i] = dp[i];
	}


	//dPrint_CompCol_Matrix("A", &A);
	//dPrint_CompCol_Matrix("U", &U);
	//dPrint_SuperNode_Matrix("L", &L);
	//print_int_vec("\nperm_r", m, perm_r);

	/* De-allocate storage */
	//SUPERLU_FREE (rhs);
	SUPERLU_FREE (perm_r);
	SUPERLU_FREE (perm_c);

	//Destroy_CompCol_Matrix(&A);
	Destroy_SuperMatrix_Store(&B);


	Destroy_SuperNode_Matrix(&L);
	Destroy_CompCol_Matrix(&U);

	StatFree(&stat);
	return 1;
}

int CallSuperLUSparseRow(int nrows, int ncols, int nvals, int* row_ptr, int* col_ind, double* vals, double* rhs)
{
	//note that SuperLU internally solves A^T q = f;
	SuperMatrix A, L, U, B;
	//int      *asub, *xa;
	int      *perm_r; /* row permutations from partial pivoting */
	int      *perm_c; /* column permutation vector */
	int      nrhs, info;
	superlu_options_t options;
	SuperLUStat_t stat;

	nrhs = 1; //number of right-hand sides

	/* Create matrix A in the format expected by SuperLU. */
	dCreate_CompRow_Matrix(&A, nrows, ncols, nvals, vals, col_ind, row_ptr, SLU_NR, SLU_D, SLU_GE); //HOTINT-sparse format in sparse-row format

	//create right-hand-side from array given in rhs
	dCreate_Dense_Matrix(&B, ncols, nrhs, rhs, ncols, SLU_DN, SLU_D, SLU_GE);

	if ( !(perm_r = intMalloc(ncols)) ) ABORT("Malloc fails for perm_r[].");
	if ( !(perm_c = intMalloc(nrows)) ) ABORT("Malloc fails for perm_c[].");

	/* Set the default input options. */
	set_default_options(&options);
	
	options.ColPerm = MMD_ATA; //works good for FE-problem
	//options.ColPerm = MMD_AT_PLUS_A / NATURAL==> does not work!

	//options.Fact = DOFACT;
	//options.Fact = SamePattern;

	/* Initialize the statistics variables. */
	StatInit(&stat);
	dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

	//read out solution from dense matrix B:
	DNformat     *Bstore;
	double       *dp;
	Bstore = (DNformat *) B.Store;
	dp = (double *) Bstore->nzval;
	for (int i = 0; i < nrows; i++) 
	{
		rhs[i] = dp[i];
	}



	/* De-allocate storage */
	//SUPERLU_FREE (rhs);
	SUPERLU_FREE (perm_r);
	SUPERLU_FREE (perm_c);

	//Destroy_CompCol_Matrix(&A);
	Destroy_SuperMatrix_Store(&B);


	Destroy_SuperNode_Matrix(&L);
	
	Destroy_CompCol_Matrix(&U);

	StatFree(&stat);
	return 1;
}

int CallSuperLUFactorize(int nrows, int ncols, int nvals, int* row_ptr, int* col_ind, double* vals, 
												 SuperMatrix*& A, SuperMatrix*& L, SuperMatrix*& U, int*& perm_r, int*& perm_c, int isFirstStep)
{
	if (nrows > 20000)
	{
		slu_outfile << "Super LU Factorize with dimension" << nrows << "..." << flush;
	}

	if (!isFirstStep)
	{	
		SUPERLU_FREE (perm_r);
		SUPERLU_FREE (perm_c);
		Destroy_SuperNode_Matrix(L);
		Destroy_CompCol_Matrix(U);
	}
	else
	{
		A = new SuperMatrix;
		L = new SuperMatrix;
		U = new SuperMatrix;
	}

	SuperMatrix B;
	int      nrhs, info;
	superlu_options_t options;
	SuperLUStat_t stat;

	nrhs = 1; //number of right-hand sides
	double* rhs = new double[ncols];
	for (int i = 0; i < ncols; i++) {rhs[i] = 0;}

	/* Create matrix A in the format expected by SuperLU. */
	dCreate_CompRow_Matrix(A, nrows, ncols, nvals, vals, col_ind, row_ptr, SLU_NR, SLU_D, SLU_GE); //HOTINT-sparse format in sparse-row format

	//create right-hand-side from array given in rhs
	dCreate_Dense_Matrix(&B, ncols, nrhs, rhs, ncols, SLU_DN, SLU_D, SLU_GE);

	if ( !(perm_r = intMalloc(ncols)) ) ABORT("Malloc fails for perm_r[].");
	if ( !(perm_c = intMalloc(nrows)) ) ABORT("Malloc fails for perm_c[].");

	/* Set the default input options. */
	set_default_options(&options);
	options.ColPerm = MMD_ATA; //options.ColPerm = NATURAL;
	//options.SymmetricMode = YES;
	
	
	// options.Fact = DOFACT;

	/* Initialize the statistics variables. */
	StatInit(&stat);
	dgssv(&options, A, perm_c, perm_r, L, U, &B, &stat, &info);

	/* De-allocate storage */
	Destroy_SuperMatrix_Store(&B);

	delete [] rhs;

	StatFree(&stat);

	if (nrows > 20000)
	{
		slu_outfile << "finished\n" << flush;
	}

	return info;
}



int CallSuperLUApply(int nrows, int ncols, double* rhs,
										 SuperMatrix*& A, SuperMatrix*& L, SuperMatrix*& U, int*& perm_r, int*& perm_c)
{
	if (nrows > 20000)
	{
		slu_outfile << "Super LU Apply with dimension" << nrows << "..." << flush;
	}

	SuperMatrix B;
	int      nrhs, info;
	superlu_options_t options;
	SuperLUStat_t stat;

	nrhs = 1; //number of right-hand sides

	//create right-hand-side from array given in rhs
	dCreate_Dense_Matrix(&B, ncols, nrhs, rhs, ncols, SLU_DN, SLU_D, SLU_GE);

	/* Set the default input options. */
	set_default_options(&options);
	options.ColPerm = MMD_ATA; //options.ColPerm = NATURAL;
	//options.Fact = FACTORED;
	options.SymmetricMode = YES;

	/* Initialize the statistics variables. */
	StatInit(&stat);
	// dgssv(&options, A, perm_c, perm_r, L, U, &B, &stat, &info);
	// AUS NETGEN  -- DANKE JOACHIM; GEHT EH
	dgstrs( TRANS, L, U, perm_c, perm_r, &B, &stat, &info );

	//read out solution from dense matrix B:
	DNformat     *Bstore;
	double       *dp;
	Bstore = (DNformat *) B.Store;
	dp = (double *) Bstore->nzval;
	for (int i = 0; i < ncols; i++) 
	{
		rhs[i] = dp[i];
	}


	/* De-allocate storage */
	Destroy_SuperMatrix_Store(&B);

	StatFree(&stat);

	if (nrows > 20000)
	{
		slu_outfile << "finished\n" << flush;
	}

	return info;
}

void DestroySuperLUMatrices(SuperMatrix*& A, SuperMatrix*& L, SuperMatrix*& U)
{
	Destroy_SuperNode_Matrix(L);
	Destroy_CompCol_Matrix(A);
	Destroy_CompCol_Matrix(U);

}
void DestroySuperLUFactorization(SuperMatrix*& L, SuperMatrix*& U)
{
	Destroy_SuperNode_Matrix(L);
	Destroy_CompCol_Matrix(U);
}


