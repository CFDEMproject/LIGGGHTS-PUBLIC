
//#include "stdafx.h"
#include "ioincludes.h"

#include <assert.h>
#include <string.h>
#include <time.h>
#include <sys/timeb.h>
#include <math.h>

#include "mystring.h" 

int CallSuperLUSparseCol(int nrows, int ncols, int nvals, int* col_ptr, int* row_ind, double* vals, double* rhs);

int CallSuperLUSparseRow(int nrows, int ncols, int nvals, int* row_ptr, int* col_ind, double* vals, double* rhs);

int CallSuperLUFactorize(int nrows, int ncols, int nvals, int* row_ptr, int* col_ind, double* vals, 
												 SuperMatrix*& A, SuperMatrix*& L, SuperMatrix*& U, int*& perm_r, int*& perm_c, int isFirstStep);

int CallSuperLUApply(int nrows, int ncols, double* rhs,
										 SuperMatrix*& A, SuperMatrix*& L, SuperMatrix*& U, int*& perm_r, int*& perm_c);

void DestroySuperLUMatrices(SuperMatrix*& A, SuperMatrix*& L, SuperMatrix*& U);
void DestroySuperLUFactorization(SuperMatrix*& L, SuperMatrix*& U);
mystr& GetSuperLUMessage();


//namespace SuperLU
//{
//	double TestFunction1(double x);
//}