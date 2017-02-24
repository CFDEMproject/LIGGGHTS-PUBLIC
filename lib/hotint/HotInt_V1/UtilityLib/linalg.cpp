//#**************************************************************
//#
//# filename:             linalg.cpp
//#
//# author:               Gerstmayr Johannes
//#
//# generated:            15.05.97
//# description:          Classes for linear and nonlinear algebra which is
//#												thought to be used in finite element or similar applications
//#												There are 2D and 3D and arbitrary size Vectors (Vector2D, Vector3D, Vector),
//#												arbitrary size matrices (Matrix), and a nonlinear system (NumNLSys)
//#												and a nonlinear solver (NumNLSolver)
//# remarks:							Indizes run from 1 to n except in Vector3D/2D
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
#include "mbs_interface.h"

#include "../WorkingModule/stdafx.h"
#include "ioincludes.h"

#include "release_assert.h"
#include <memory.h>

#include <math.h>

#include "tarray.h"    
#include "mystring.h"  
#include "femathhelperfunctions.h"

#include "lapack_routines.h"
#include "pardiso.h"

#include "..\superlu\superlumain.h"

//for output possibility within femath
//#include "..\workingmodule\WorkingModuleBaseClass.h"
extern UserOutputInterface * global_uo;

//count news in sparse matrices!!!
int gensparsemat = 0;
void GenSparsemat() {gensparsemat++;}
int GetGenSparsemat() {return gensparsemat;}

#ifdef gencnt
extern int *genvec;
extern int *genmat;
#endif



//-----------------------------------------------------------------
//--------------    CLASS VECTOR    -------------------------------
//-----------------------------------------------------------------
//Generation of user-defiened Vector
//Operations: +,-,*,<< (also Matrices)


Vector::Vector(const Vector& veci, int noownmem)
{
	if (noownmem)
	{
		l = veci.l;
		lalloc = 0;
		ownmemory = 0;
		if (l==0) vec = NULL;
		else vec = veci.vec;
	}
	else
	{
		Vector::Vector(veci);
	}
}

Vector::Vector(int leni, int nofill, double fill)
{
	vec=NULL;
	GenVector(leni);
	if (!nofill && l!=0)
	{
		for (int i=0; i<l; i++) {vec[i]=fill; }
	}
};

Vector::Vector(const Vector& veci)
{
	l=veci.l;
	lalloc = l;
	if (l==0)
	{
		vec=NULL; 
	}
	else
	{
#ifdef gencnt
		(*genvec)++;
#endif
		ownmemory = 1;
		vec=new double[l];
		memcpy(vec,veci.vec,l*sizeof(double));
	}
};

Vector::Vector(const char* vector, char osbc, char csbc, char commac) // (RL)
{
	//mystr must be formatted like following examples:  [1 2 4] or [1,2,3]
	mystr el = "\n"; //end of line
	char elc = 'n';
	mystr osb = osbc;  //open square bracket, for vector
	mystr csb = csbc;  //closing square bracket
	mystr semicolon = ";"; //double stop
	mystr comma = commac; //comma
	int endstr = 0;
	int pos = 0;
	mystr dataname;
	char ch;

	mystr dummycomment; //dummy comment string
	mystr str(vector);
	mystr text = str.GetStringInBrackets(pos, osbc, csbc);					

	TArray<double> ta_vec(0);
	if(text.Find(el)==-1  && text.Find(semicolon)==-1)
	{
		// vector has no semicolon and is written in one single line
		mystr s_value;						
		int posvec = 0;
		int possval = 0;
		int comma_found = 0;
		while(posvec != -1)
		{
			comma_found = text.GetUntil(posvec, commac, s_value); // 'value,'
			s_value.ReadLeadingSpaces(possval);
			ta_vec.Add(s_value.MakeDouble());

			if(comma_found)posvec++;
		}											
	}			
	SetVector(ta_vec);
}

int Vector :: IsValid(double maxval)
{
	for (int i=0; i<l; i++)
	{
		if (IsNaN(vec[i]) || fabs(vec[i])>=maxval) {return 0;}
	}
	return 1;
}


void Vector::GenVector(int leni)
{
	l=leni;
	lalloc = leni;
	if (l!=0) 
	{
#ifdef gencnt
		(*genvec)++;
#endif
		ownmemory = 1;
		vec=new double[l];
	}
	else
	{
		vec = NULL;
	}
};

Vector& Vector::operator= (const Vector& veci)
{
	if (this == &veci) {return *this;}
	if (veci.l!=0)
	{
		if ((((lalloc < veci.l && ownmemory) || (l < veci.l && !ownmemory)) && lalloc_act) || (l != veci.l && !lalloc_act)) 
		{
#ifdef gencnt
			(*genvec)++;
#endif
			l=veci.l;
			if (vec!=NULL && ownmemory && lalloc) {delete [] vec; }
			ownmemory = 1;
			vec=new double[l];
			lalloc = l;
		} 
		else 
		{
			l=veci.l; 
		}
		memcpy(vec,veci.vec,l*sizeof(double));
	}
	else 
	{
		l=0; 
		//lalloc = 0; vec = NULL;
	}
	return *this;
};

Vector& Vector::operator= (const Vector3D& veci)
{
	if ((lalloc < 3 && lalloc_act) || (l != 3 && !lalloc_act)) 
	{
#ifdef gencnt
		(*genvec)++;
#endif
		l=3;
		if (vec!=NULL && ownmemory && lalloc) {delete [] vec; }
		ownmemory = 1;
		vec=new double[l];
		lalloc = l;
	} 
	else {l=3; }
	vec[0] = veci.X(); vec[1] = veci.Y(); vec[2] = veci.Z();
	return *this;
};

int Vector :: operator == (const Vector& vec1)
{
	if (vec1.GetLen() != GetLen()) 
	{cout << "Warning: vectors compared with different length!" << endl;}

	int i;
	for (i = 1; i <= vec1.GetLen(); i++)
	{
		if (vec1.Get(i) != Get(i)) return 0;
	}

	return 1;
}

void Vector::SetLen(int leni)
{
	if ((((lalloc < leni && ownmemory) || (!ownmemory && leni > l)) && lalloc_act) || (l != leni && !lalloc_act)) 
	{
		if (vec!=NULL && ownmemory && lalloc) {delete [] vec;}
		GenVector(leni);
	}
	l = leni;
}

double Vector::GetNorm() const
{
	return sqrt((*this)*(*this));
}

double Vector::MaxNorm() const
{
	double res=-1e100;
	for(int i=0; i<l; i++)
	{
		res=Maximum(fabs(vec[i]),res);
	}
	if (res==-1E100) {res=0;}
	return res;
}

double Vector::MinNorm() const
{
	double res=1E100;
	for(int i=0; i<l; i++)
	{
		if (fabs(vec[i])<res && vec[i]!=0) {res=fabs(vec[i]);}
	}
	if (res==1E100) {res=0;}
	return res;
}
Vector Vector::GetNormed2DNormalVec() const
{
	Vector v(-(*this)(2),(*this)(1));
	double len=v.GetNorm();
	if (len!=0) {v*=1/len;}
	return v;
}
Vector Vector::GetNormed2DVec() const
{
	Vector v((*this));
	double len=v.GetNorm();
	if (len!=0) {v*=1/len;}
	return v;
}

double Vector::Sum() const
{
	int i;
	double s = 0;
	for (i = 0; i < l; i++)
	{
		s += vec[i];
	}
	return s;
}

double Vector::Sum(int n) const
{
	int i;
	double s = 0;
	for (i = 0; i < n && i < l; i++)
	{
		s += vec[i];
	}
	return s;
}


void Vector::FillWithZeros()
{
	for (int i=0; i<l; i++) {vec[i]=0; }
}
/*
void Vector::FillInVector(const Vector& v, const IVector& ref)
{
for(int i=1; i<=ref.Length(); i++)
{
(*this)[ref[i]]+=v[i];
}
}
*/
Vector Vector::Append(const Vector& v) const
{
	if (v.l == 0) {return *this;}
	int i;
	Vector res(l+v.l);
	for (i=1;i<=l;i++)
	{
		res(i)=(*this)(i);
	}
	for (i=1;i<=v.l;i++)
	{
		res(i+l)=v(i);
	}
	return res;
}

void Vector::Insert(const Vector& v, int pos)
{
	if (v.l == 0) {return;}
	int i;
	for (i=1;i<=v.l;i++)
	{
		Elem(i+pos-1) = v(i);
	}
}

Vector Vector::SubVector(int from, int to) const
{
	int i;
	Vector res(to-from+1);
	for (i=from;i<=to;i++)
	{
		res(i-from+1)=(*this)(i);
	}
	return res;
}

void Vector::MultAdd(const double& k, const Vector& v)
{
#ifndef __QUICKMATH
	release_assert(v.l==l);
#endif
	for (int i=0; i<l; i++)
	{
		vec[i]+=k*v.vec[i];
	}
}

Vector& Vector::operator+= (const Vector& v1)
{
#ifndef __QUICKMATH
	release_assert(v1.l==l);
#endif
	for (int i=0;i<l;i++)
	{
		vec[i]+=v1.vec[i];
	}
	return *this;
}

Vector& Vector::operator-= (const Vector& v1)
{
#ifndef __QUICKMATH
	release_assert(v1.l==l);
#endif
	for (int i=0;i<l;i++)
	{
		vec[i]-=v1.vec[i];
	}
	return *this;
}

Vector& Vector::operator*= (const double& v)
{
	for (int i=0;i<l;i++)
	{
		vec[i]*=v;
	}
	return *this;
}

ostream& operator<<(ostream& os, const Vector& v)
{
	if (v.GetLen()<=70)
	{
		os << "[";
		for (int i=1; i<v.GetLen(); i++)
		{
			os << v(i) << ", ";
		}
		if (v.GetLen()!=0)
		{
			os << v(v.GetLen());
		}
		os << "]";
	} else
	{
		os << "[\n";
		for (int i=1; i<=v.GetLen(); i++)
		{
			os << v(i);
			if (i!=v.GetLen()) {os << ", ";}
			if (i%8==0 && i!=v.GetLen()) {os << "\n  col(" << i << "): ";}
		}
		os << "]\n";
	}
	return os;
}

Vector operator+ (const Vector& v1, const Vector& v2)
{
	Vector result=v1;
	result+=v2;
	return result;
}

Vector operator- (const Vector& v1, const Vector& v2)
{
	Vector result=v1;
	result-=v2;
	return result;
}

Vector operator* (const Vector& vec, const double& val)
{
	Vector result=vec;
	result*=val;
	return result;
}

Vector operator* (const double& val, const Vector& vec)
{
	Vector result=vec;
	result*=val;
	return result;
}

double operator* (const Vector& vec1, const Vector& vec2)
{
#ifndef __QUICKMATH
	release_assert(vec1.l==vec2.l);
#endif
	double result=0;
	for (int i=0;i<vec1.l;i++)
	{
		result+=vec1.vec[i]*vec2.vec[i];
	}
	return result;
}

// Vector entries are equally distributed values of piecewise linear function on [-1, 1]
// Interpolate returns value of function in point ploc € [-1, 1]
// Vector entries specify function values, ordered as
//   [ 1          2         3   ... width ]       -->      [x=-1         x=1]
double Vector::Interpolate(const double ploc) const
{
	int imin;
	double hx = 2./(l-1.);
	// lower indices of cell containing ploc
	imin = int((1.-ploc) / hx ) +1;
	if (imin >= l) imin = l-1;
	// "mesh sizes"
	double xmin = hx * (imin-1)-1;
	double val = ( (*this)(imin+1)*(hx-ploc+xmin) + 
								(*this)(imin)*(ploc-xmin) ) / hx;
	return val;

}

mystr& Vector::MakeString(mystr intro)
{
	mystr* buffer = new mystr(intro);
	*buffer +=("["); *buffer += mystr(l); *buffer +=("]");
	*buffer +=(" = { ");
	for (int i=0; i < this->l; i++)
	{
		if( i ) *buffer += ", ";
		*buffer += mystr(Get0(i));
	}
	*buffer += mystr(" };");
	return *buffer;
}



//-----------------------------------------------------------------
//--------------    CLASS MATRIX    -------------------------------
//-----------------------------------------------------------------
//Generation of Matrices with user-definied size
//Operations: +,-,*,<<

Matrix::Matrix(double x)
{
	rows=1; cols=1;
#ifdef gencnt
	(*genmat)++;
#endif
	mat=new double[1];
	lalloc = 1;
	mat[0]=x;
}
Matrix::Matrix(double x11, double x12, double x21, double x22)
{
	rows=2; cols=2;
#ifdef gencnt
	(*genmat)++;
#endif
	mat=new double[4];
	lalloc = 4;
	mat[0]=x11; mat[1]=x12;
	mat[2]=x21; mat[3]=x22;
}
Matrix::Matrix(double x11, double x12, double x13,
							 double x21, double x22, double x23,
							 double x31, double x32, double x33)
{
	rows=3; cols=3;
#ifdef gencnt
	(*genmat)++;
#endif
	mat=new double[9];
	lalloc = 9;
	mat[0]=x11; mat[1]=x12; mat[2]=x13;
	mat[3]=x21; mat[4]=x22; mat[5]=x23;
	mat[6]=x31; mat[7]=x32; mat[8]=x33;
}

Matrix::Matrix(const Matrix3D mat3d)
{
	rows=3; cols=3;
#ifdef gencnt
	(*genmat)++;
#endif
	mat=new double[9];
	lalloc = 9;
	mat[0]=mat3d(1,1); mat[1]=mat3d(1,2); mat[2]=mat3d(1,3);
	mat[3]=mat3d(2,1); mat[4]=mat3d(2,2); mat[5]=mat3d(2,3);
	mat[6]=mat3d(3,1); mat[7]=mat3d(3,2); mat[8]=mat3d(3,3);
}

int bigmats = 0;
int bigmatval = 10000;
Matrix::Matrix(int rowsi, int colsi, int nofill, double fill)
{
	rows=rowsi;
	cols=colsi;
	if (rows*cols!=0)
	{
#ifdef gencnt
		(*genmat)++;
#endif
		mat=new double[rows*cols];
		lalloc = rows*cols;
		/*		if (lalloc > bigmatval) 
		{
		bigmats++;
		(*global_uo) << "BigMat=" << lalloc*8/1000 << "K\n";
		}*/
		if (!nofill)
		{
			int size=rows*cols;
			for (int i=0; i<size; i++) {mat[i]=fill; }
		}
	}
	else
	{
		mat = NULL;
		lalloc=0;
	}
};


Matrix::Matrix(const Matrix& mati)
{
	rows=mati.rows;
	cols=mati.cols;
	if (mati.mat==NULL)
	{
		mat = NULL;
		lalloc = 0;
	}
	else
	{
#ifdef gencnt
		(*genmat)++;
#endif
		mat=new double[rows*cols];
		lalloc = rows*cols;
		/*		if (lalloc > bigmatval) 
		{
		bigmats++;
		(*global_uo) << "BigMat=" << lalloc*8/1000 << "K\n";
		}*/
		memcpy(mat,mati.mat,rows*cols*sizeof(double));
	}
};

void Matrix::CopyFrom(const Matrix& m, int row1, int col1, int row2, int col2)
{
	Resize(row2-row1+1, col2-col1+1, 1);
	int i=0;
	for (int row = row1; row <= row2; row++)
	{
		for (int col = col1; col <= col2; col++)
		{
			mat[i++] = m(row,col);
		}
	}
}

void Matrix::CopyFrom(const Matrix& m, int row1, int col1, int row2, int col2, const IVector& r)
{
	Resize(row2-row1+1, col2-col1+1, 1);
	int i=0;
	for (int row = row1; row <= row2; row++)
	{
		for (int col = col1; col <= col2; col++)
		{
			mat[i++] = m(r(row),r(col));
		}
	}
}

void Matrix::CopyFrom(const SparseMatrix& m, int row1, int col1, int row2, int col2)
{
	Resize(row2-row1+1, col2-col1+1, 1);

	for (int i = row1; i <= row2; i++)
	{
		int lastsparsecol = 1;
		for (int j = col1; j <= col2; j++)
		{
			Elem0(i-row1,j-col1) = m.Get(i,j,lastsparsecol);
		}
	}
}

void Matrix::CopyFrom(const SparseMatrix& m, int row1, int col1, int row2, int col2, const IVector& r)
{
	Resize(row2-row1+1, col2-col1+1, 1);
	for (int i = row1; i <= row2; i++)
	{
		for (int j = col1; j <= col2; j++)
		{
			Elem0(i-row1,j-col1) = m(r(i),r(j));
		}
	}
}

Matrix& Matrix::operator= (const Matrix& mati)
{
	if (this == &mati) {return *this;}
	if (rows!=mati.rows || cols!=mati.cols)
	{
		if ((mati.rows*mati.cols > lalloc && lalloc_act) || (mati.rows*mati.cols != rows*cols && !lalloc_act))
		{
			if (mat!=NULL && lalloc != 0) {delete [] mat; mat=NULL; }
			rows=mati.rows;
			cols=mati.cols;

			if (rows*cols!=0)
			{
#ifdef gencnt
				(*genmat)++;
#endif
				mat=new double[rows*cols];
				lalloc = rows*cols;
				/*				if (lalloc > bigmatval) 
				{
				bigmats++;
				(*global_uo) << "BigMat=" << lalloc*8/1000 << "K\n";
				}*/
			}
		}
		else
		{
			rows=mati.rows;
			cols=mati.cols;
		}
	}
	if (mati.mat==NULL || rows*cols==0)
	{
		mat = NULL;
	}
	else
	{
		memcpy(mat,mati.mat,rows*cols*sizeof(double));
	}
	return *this;
};

void Matrix::Resize(int rowsi, int colsi, int nofill, double fill)
{
	if ((rowsi*colsi > MaxAllocatedSize() && lalloc_act) || (rowsi*colsi != rows*cols && !lalloc_act) || mat==NULL)
	{
#ifndef __QUICKMATH
		release_assert(IsConstSizeMatrix() == 0);
#endif
		if (mat!=NULL && lalloc != 0) {delete [] mat; mat=NULL; }
		rows=rowsi;
		cols=colsi;
		if (rows*cols!=0)
		{
#ifdef gencnt
			(*genmat)++;
#endif
			lalloc = rows*cols;
			mat = new double[rows*cols];
			if (mat == NULL) 
			{
				(*global_uo) << "ERROR: memory allocation failed: Matrix::Resize(...) !!!!\n";
				(*global_uo) << "MEM=" << (8*rows*cols)/(1024*1024) << "MB, mat=[" << rows << "x" << cols << "]\n";
				release_assert(0);
			}
		}
	}
	else
	{
		rows=rowsi;
		cols=colsi;
	}
	if (!nofill)
	{
		int size=rows*cols;
		for (int i=0; i<size; i++) {mat[i]=fill; }
	}
}

double Matrix::MaxNorm() const
{
	double res=0;
	for(int i=0; i<rows*cols; i++)
	{
		res=Maximum(fabs(mat[i]),res);
	}
	return res;
}

double Matrix::MinNorm() const
{
	double res=10E100;
	for(int i=0; i<rows*cols; i++)
	{
		if (fabs(mat[i])<res && mat[i]!=0) {res=fabs(mat[i]);}
	}
	if (res==10E100) {res=1;}
	return res;
}

double Matrix::Norm2() const
{
	double res=0;
	for(int i=0; i<rows*cols; i++)
	{
		res+=mat[i]*mat[i];
	}
	return sqrt(res);
}

void Matrix::FillWithZeros()
{
	for(int i=0; i<rows*cols; i++)
	{
		mat[i]=0;
	}
}

void Matrix::SetAll(double x)
{
	for(int i=0; i<rows*cols; i++)
	{
		mat[i]=x;
	}
}

void Matrix::FillWithZeros(int row, int col, int nrows, int ncols)
{
	for (int i=row; i<row+nrows; i++)
	{
		for (int j=col; j<col+ncols; j++)
		{
			Elem(i,j) = 0;
		}
	}
}

/*
void Matrix::FillInMatrix(const Matrix& K, const IVector& ref1, const IVector& ref2)
{
for(int i=1; i<=K.Getrows(); i++)
{
for(int j=1; j<=K.Getcols(); j++)
{
(*this)(ref1[i],ref2[j])+=K(i,j);
}
}
}

void Matrix::FillInBMatrix(const Matrix& K, const IVector& ref1, const IVector& ref2)
{
for(int i=1; i<=K.Getrows(); i++)
{
for(int j=1; j<=K.Getcols(); j++)
{
if (ref1[i]<=ref2[j]) {Bset(ref1[i],ref2[j],BGet(ref1[i],ref2[j])+K(i,j));}
}
}
}
*/
Matrix Matrix::BToMatrix()
{
	Matrix m(rows,rows);
	for (int i=1; i<=rows; i++)
	{
		for (int j=1; j<=rows; j++)
		{
			m(i,j)=BGet(i,j);
		}
	}
	return m;
}

Matrix Matrix::MatrixToB()
{
	int bw = GetBandWidth()+1;
	Matrix m(rows,bw);

	for (int i=1; i<=rows; i++)
	{
		for (int j=1; j<=bw; j++)
		{
			if (j+i-1 <= cols)
				m.Bset(i,j+i-1,(*this)(i,j+i-1));
		}
	}
	return m;
}

Vector Matrix::GetColVec(int col) const
{
#ifndef __QUICKMATH
	release_assert((col>0) && (col<=cols));
#endif
	Vector v(rows);
	int end=(rows-1)*(cols)+col-1;
	int j=1;
	for (int i=col-1;i<=end;i+=cols)
	{
		v(j++)=mat[i];
	}
	return v;
}

Vector Matrix::GetRowVec(int row) const
{
#ifndef __QUICKMATH
	release_assert((row>0) && (row<=rows));
#endif
	Vector v(cols);
	int offset=(row-1)*cols-1;
	for (int i=1;i<=cols;i++)
	{
		v(i)=mat[offset+i];
	}
	return v;
}

void Matrix::GetColVec(int col, Vector& v) const
{
#ifndef __QUICKMATH
	release_assert((col>0) && (col<=cols));
#endif
	v.SetLen(rows);
	int end=(rows-1)*(cols)+col-1;
	int j=1;
	for (int i=col-1;i<=end;i+=cols)
	{
		v(j++)=mat[i];
	}
}

void Matrix::GetRowVec(int row, Vector& v) const
{
#ifndef __QUICKMATH
	release_assert((row>0) && (row<=rows));
#endif
	v.SetLen(cols);
	int offset=(row-1)*cols-1;
	for (int i=1;i<=cols;i++)
	{
		v(i)=mat[offset+i];
	}
}

double Matrix::Det() const
{
#ifndef __QUICKMATH
	release_assert(IsSquare());
#endif
	if (cols==1) {return mat[0];}
	if (cols==2) {return mat[0]*mat[3]-mat[1]*mat[2];}
	if (cols==3) {return mat[0]*mat[4]*mat[8]-mat[0]*mat[5]*mat[7]-mat[3]*mat[1]*mat[8]
	+mat[3]*mat[2]*mat[7]+mat[6]*mat[1]*mat[5]-mat[6]*mat[2]*mat[4];}
#ifndef __QUICKMATH
	release_assert(cols>3);
#endif
	return 0;
}

//Returns the diagonal values
Vector Matrix::GetDiag() const
{
#ifndef __QUICKMATH
	release_assert(IsSquare());
#endif

	Vector ev(cols);

	for (int i=1; i<= cols; i++)
	{
		ev(i) = Get(i,i);  
	}
	return ev;
}

double Matrix::Trace() const
{
#ifndef __QUICKMATH
	release_assert(IsSquare());
#endif

	double sum=0;

	for (int i=1; i<= cols; i++)
	{
		sum += Get(i,i);  
	}
	return sum;
}

void Matrix::SetColVec(const Vector& x, int col)
{
#ifndef __QUICKMATH
	release_assert((col>0) && (col<=cols));
#endif
	int end=(rows-1)*(cols)+col-1;
	int j=1;
	for (int i = col-1; i <= end; i += cols)
	{
		mat[i]=x(j++);
	}
}

//Sets the column-vector at clomn col, taken from Matrix m, column mcol
void Matrix::SetColVec(const Matrix& m, int mcol, int col)
{
#ifndef __QUICKMATH
	release_assert((col>0) && (col<=cols) && (mcol <= m.cols) && (cols == m.cols));
#endif

	int end=(rows-1)*(cols)+col-1;
	int diff = mcol-col;
	for (int i = col-1; i <= end; i += cols)
	{
		mat[i]=m.mat[i+diff];
	}
}

void Matrix::SetRowVec(const Vector& x, int row)
{
#ifndef __QUICKMATH
	release_assert((row>0) && (row <= rows));
#endif
	int offset = (row-1)*cols-1;
	for (int i = 1;i <= cols; i++)
	{
		mat[offset+i] = x(i);
	}
}

double Matrix::MultCol(const Matrix& m, int mcol, int col)
{
#ifndef __QUICKMATH
	release_assert((col>0) && (col<=cols) && (mcol <= m.cols) && (cols == m.cols));
#endif

	double val = 0;
	int end=(rows-1)*(cols)+col-1;
	int diff = mcol-col;

	for (int i = col-1; i <= end; i += cols)
	{
		val += mat[i]*m.mat[i+diff];
	}
	return val;
}

int Matrix::GetBandWidth() const
{
	if (!IsSquare()) {return rows;}
	int bw=0; //minimal, wenn Nullmatrix
	for (int i=1;i<=rows;i++)
	{
		for (int j=0;j<=rows-i;j++)
		{
			if ((*this)(j+1,i+j)!=0 ||
				(*this)(i+j,j+1)!=0) {bw=i;}
		}
	}
	return bw;
}

int Matrix::IsSymmetric()
{
#ifndef __QUICKMATH
	release_assert(IsSquare());
#endif

	for (int i=0; i<rows; i++)
	{
		for (int j=0; j<i; j++)
		{
			if (mat[i*cols+j]!=mat[j*cols+i]) {return 0;}
		}
	}
	return 1;
}

void Matrix::MakeSymmetric()
{
#ifndef __QUICKMATH
	release_assert(IsSquare());
#endif

	double x;
	for (int i=0; i<rows; i++)
	{
		for (int j=0; j<i; j++)
		{
			x = 0.5*(mat[i*cols+j] + mat[j*cols+i]);
			mat[i*cols+j] = x;
			mat[j*cols+i] = x;
		}
	}
}

int Matrix::IsSymmetric(double eps)
{
#ifndef __QUICKMATH
	release_assert(IsSquare());
#endif

	for (int i=0; i<rows; i++)
	{
		for (int j=0; j<i; j++)
		{
			if (fabs(mat[i*cols+j] - mat[j*cols+i]) > eps) {return 0;}
		}
	}
	return 1;
}


void Matrix::TpYs()
{
	if (IsSquare())
	{
		for (int i=0; i<rows; i++)
		{
			for (int j=0; j<i; j++)
			{
				Swap(mat[i*cols+j],mat[j*cols+i]);
			}
		}
	} 
	else
	{
		Matrix mt(cols,rows);

		for (int i=0; i<rows; i++)
		{
			for (int j=0; j<cols; j++)
			{
				mt.Elem0(j,i) = Get0(i,j);
				//mt.mat[j*rows+i]=mat[i*cols+j];
			}
		}
		//CopyFrom(mt,1,1,cols,rows);
		*this = mt;
	}
}

Matrix Matrix::GetTp() const 
{
	Matrix mt(cols,rows);

	for (int i=0; i<rows; i++)
	{
		for (int j=0; j<cols; j++)
		{
			mt.mat[j*rows+i]=mat[i*cols+j];
		}
	}
	return mt;
}

//$ SW 2013-10-9: added
//calculates *= A where A is a square matrix (e.g. a rotation matrix) of dimension 3x3, 2x2 or 1x1
void Matrix::ApplySqrMat(const MatrixXD &A)
{
#ifndef __QUICKMATH
	release_assert("ERROR: ApplySqrMat dimensions wrong" && (cols==A.Getrows() && A.Getrows()==A.Getcols()));
#endif
	if (A.Getrows()==3)
	{
		double tmp1,tmp2,tmp3;
		for (int i=0; i<rows; i++)
		{
			tmp1= mat[i*3]*A(1,1)+mat[i*3+1]*A(2,1)+mat[i*3+2]*A(3,1);
			tmp2= mat[i*3]*A(1,2)+mat[i*3+1]*A(2,2)+mat[i*3+2]*A(3,2);
			tmp3= mat[i*3]*A(1,3)+mat[i*3+1]*A(2,3)+mat[i*3+2]*A(3,3);
			mat[i*3]= tmp1;
			mat[i*3+1]= tmp2;
			mat[i*3+2]= tmp3;
		}
		return;
	}
	if (A.Getrows()==2)
	{
		double tmp1,tmp2;
		for (int i=0; i<rows; i++)
		{
			tmp1= mat[i*2]*A(1,1)+mat[i*2+1]*A(2,1);
			tmp2= mat[i*2]*A(1,2)+mat[i*2+1]*A(2,2);
			mat[i*2]= tmp1;
			mat[i*2+1]= tmp2;
		}
		return;
	}
	if (A.Getrows()==1)
	{
		operator*=(A(1,1));
		return;
	}
	#ifndef __QUICKMATH
	release_assert("ERROR: ApplySqrMat works only for matrices of dimension 1, 2 or 3" && false);
	#endif
}

Matrix& Matrix::operator+= (const Matrix& m)
{
#ifndef __QUICKMATH
	release_assert(m.rows==rows && m.cols==cols);
#endif
	int size=cols*rows;
	for (int i=0;i<size;i++)
	{
		mat[i]+=m.mat[i];
	}
	return *this;
}

Matrix& Matrix::operator-= (const Matrix& m)
{
#ifndef __QUICKMATH
	release_assert(m.rows==rows && m.cols==cols);
#endif
	int size=cols*rows;
	for (int i=0;i<size;i++)
	{
		mat[i]-=m.mat[i];
	}
	return *this;
}

Matrix& Matrix::operator*= (const double& val)
{
	int size=cols*rows;
	for (int i=0;i<size;i++)
	{
		mat[i]*=val;
	}
	return *this;
}

void Matrix::AddSubmatrix(const Matrix& sm, int sr, int sc, int r, int c, int nr, int nc, double x)
{
	int off = -cols+c-1;
	for(int i=0; i<nr; i++) 
	{
		for(int j=0; j<nc; j++) 
		{
			mat[(i+r)*cols+j+off]+=x*sm(i+sr,j+sc);
			//Elem(i+r,j+c)+=x*sm(i+sr,j+sc);
		}
	}
}

void Matrix::SetSubmatrix(const Matrix& sm, int r, int c)
{
	int off = -cols+c-1;
	for(int i=0; i<sm.rows; i++) 
	{
		for(int j=0; j<sm.cols; j++) 
		{
			mat[(i+r)*cols+j+off] = sm.Get0(i,j);
		}
	}
}

void Matrix::SetSubmatrix(const Matrix& sm, int r, int c, double x)
{
	int off = -cols+c-1;
	for(int i=0; i<sm.rows; i++) 
	{
		for(int j=0; j<sm.cols; j++) 
		{
			mat[(i+r)*cols+j+off] = x*sm.Get0(i,j);
		}
	}
}

void Matrix::SetSubmatrix(const Matrix3D& sm, int r, int c)
{
	int off = -cols+c-1;
	for(int i=0; i<sm.Getrows(); i++) 
	{
		for(int j=0; j<sm.Getcols(); j++) 
		{
			mat[(i+r)*cols+j+off] = sm.Get0(i,j);
			//=>mat[(i+r-1)*cols+c+j-1]+=x*sm(i,j);
			//=>Elem(i+r,j+c)+=x*sm(i+sr,j+sc);
		}
	}
}

void Matrix::MulCol(int col, double val)
{
#ifndef __QUICKMATH
	release_assert((col>0) && (col<=cols));
#endif

	int end=(rows-1)*(cols)+col-1;
	for (int i=col-1;i<=end;i+=cols) {mat[i]*=val; }
}

void Matrix::MulRow(int row, double val)
{
#ifndef __QUICKMATH
	release_assert((row>0) && (row<=rows));
#endif
	for (int i=(row-1)*cols;i<=row*cols-1;i++) {mat[i]*=val; }
}

void Matrix::InsertMatrix(int row, int col, const Matrix& m)
{
#ifndef __QUICKMATH
	release_assert((row>0) && (col>0) && (row+m.rows-1<=rows) && (col+m.cols-1<=cols));
#endif
	int i,j;
	for (i = 1; i <= m.rows; i++)
	{
		for (j = 1; j <= m.cols; j++)
		{
			Elem(row+i-1,col+j-1) = m(i,j);
		}
	}
}

void Matrix::InsertMatrix(int row, int col, const Matrix& m, int fromrow, int fromcol, int mrows, int mcols)
{
#ifndef __QUICKMATH
	release_assert((row>0) && (col>0) && (row+mrows-1<=rows) && (col+mcols-1<=cols) && (mrows<=m.rows) && (mcols<=m.cols));
#endif
	int i,j;
	for (i = 1; i <= mrows; i++)
	{
		for (j = 1; j <= mcols; j++)
		{
			Elem(row+i-1,col+j-1) = m(fromrow-1+i,fromcol-1+j);
		}
	}
}


void Matrix::AddMatrix(int row, int col, const Matrix& m)
{
#ifndef __QUICKMATH
	release_assert((row>0) && (col>0) && (row+m.rows-1<=rows) && (col+m.cols-1<=cols));
#endif
	int i,j;
	for (i = 1; i <= m.rows; i++)
	{
		for (j = 1; j <= m.cols; j++)
		{
			Elem(row+i-1,col+j-1) += m(i,j);
		}
	}
}

void Matrix::AddMatrix(const TArray<int>& rowref, const TArray<int>& colref, const Matrix& m)
{
	int i,j;
	int k;
	for (i = 0; i < m.rows; i++)
	{
		k = (rowref.Get0(i)-1)*cols;
		for (j = 0; j < m.cols; j++)
		{
			mat[k+colref.Get0(j)-1] += m.Get0(i,j);
		}
	}
}

void Matrix::AddRowVec(int row, const Vector& vec)
{
#ifndef __QUICKMATH
	release_assert((row>0) && (row<=rows) && vec.GetLen()==cols);
#endif
	int j=1;
	for (int i=(row-1)*cols;i<=row*cols-1;i++)
	{
		mat[i]+=vec(j++);
	}
}

void Matrix::AddColVec(int col, const Vector& vec)
{
#ifndef __QUICKMATH
	release_assert((col>0) && (col<=cols) && vec.GetLen()==rows);
#endif

	int end=(rows-1)*(cols)+col-1;
	int j=1;
	for (int i = col-1; i <= end; i += cols)
	{
		mat[i] += vec(j++);
	}
}

//Adds fromcol multiplied with fact to column toCol
void Matrix::AddColVec(int fromCol, int toCol, double fact)
{
#ifndef __QUICKMATH
	release_assert((fromCol>0) && (fromCol<=cols) && (toCol > 0) && (toCol <= cols));
#endif

	int end=(rows-1)*(cols)+toCol-1;
	int diff = fromCol-toCol;

	for (int i = toCol-1; i <= end; i += cols)
	{
		mat[i] += fact*mat[i+diff];
	}
}

void Matrix::AddRowVec(int fromRow, int toRow, double Fact)
{
#ifndef __QUICKMATH
	release_assert((toRow>0) && (toRow<=rows) &&
		(fromRow>0) && (fromRow<=rows));
#endif
	int j=(fromRow-1)*cols-(toRow-1)*cols;
	int end=toRow*cols-1;
	for (int i=(toRow-1)*cols;i<=end;i++)
	{
		mat[i]+=Fact*mat[j+i];
	}
}

void Matrix::AddRowVec(int fromRow, int toRow, double Fact, int fromCol, int toCol)
{
#ifndef __QUICKMATH
	release_assert((toRow>0) && (toRow<=rows) &&
		(fromRow>0) && (fromRow<=rows));
#endif

	int j=(fromRow-1)*cols-(toRow-1)*cols;
	int end=(toRow-1)*cols+toCol-1;
	for (int i=(toRow-1)*cols+fromCol-1;i<=end;i++)
	{
		mat[i]+=Fact*mat[j+i];
	}
}

void Matrix::AddRowVecB(int fromRow, int toRow, const double& Fact, int bw)
{
#ifndef __QUICKMATH
	release_assert((toRow>0) && (toRow<=rows) &&
		(fromRow>0) && (fromRow<=rows));
#endif
	int j=(fromRow-1)*cols-(toRow-1)*cols;
	int minRow = Minimum(fromRow, toRow);
	int maxRow = Maximum(fromRow, toRow);
	int start = (toRow-1)*cols + Maximum(minRow-bw,0);
	int end = toRow*cols-1 - Maximum(cols-maxRow-bw,0);
	//cout << "save=" << Maximum(minRow-bw,0)+Maximum(cols-maxRow-bw,0) << endl;
	for (int i = start; i <= end; i++)
	{
		mat[i]+=Fact*mat[j+i];
	}
}

void Matrix::SwapRows(int r1, int r2)
{
	if (r1==r2) {return;}
#ifndef __QUICKMATH
	release_assert((r1>0) && (r1<=rows) && (r2>0) && (r2<=rows));
#endif
	for (int i=1; i<=cols; i++) {Swap((*this)(r1,i),(*this)(r2,i));}
}

void Matrix::SwapCols(int c1, int c2)
{
	if (c1==c2) {return;}
#ifndef __QUICKMATH
	release_assert((c1>0) && (c1<=cols) && (c2>0) && (c2<=cols));
#endif
	for (int i=1; i<=rows; i++) {Swap((*this)(i,c1),(*this)(i,c2));}
}

//slow convert doubles matrix ....
//returns 1, if sucsess
int Matrix::Invert()
{
	//für 2*2 Matrix optimierter Lösungsalgorithmus
	if (cols==1 && rows==1)
	{
		if (mat[0] == 0) 
		{
			return 0;
		}
		mat[0] = 1/mat[0];
		return 1;
	}

	if (cols==2 && rows==2)
	{
		double det=mat[0]*mat[3]-mat[1]*mat[2];
		if (det==0) {return 0;}
		det=1/det;
		Swap(mat[0],mat[3]);
		mat[0]*=det;
		mat[1]*=-det;
		mat[2]*=-det;
		mat[3]*=det;
		return 1;
	}
	//TMStartTimer(9);

#ifndef __QUICKMATH
	release_assert(cols==rows && cols*rows!=0 && mat!=NULL);
#endif

	//Insert identity-matrix on left-hand-side
	mystatic Matrix m;
	m.SetSize(rows,cols*2);

	int i,j,k;
	for (j=1; j<=rows; j++)
	{
		for (i=1; i<=cols; i++) {m(j,i)=(*this)(j,i);}
		for (i=cols+1; i<=2*cols; i++)
		{
			if (i-cols==j) {m(j,i)=1;}
			else {m(j,i)=0;};
		}
	}

	// Solve lower triangular Matrix

	for (j=1; j<=cols; j++)
	{
		double pivot=fabs(m(j,j));
		int pivotpos=j;
		for (k=j+1; k<=rows; k++)
		{
			if (fabs(m(k,j)) > pivot) {pivotpos=k; pivot=fabs(m(k,j));}
		}
		//cout << "Pivot = " << pivot << endl;
		if (pivot == 0) 
		{
			//(*global_uo) << "Invert: problems with column " << j << "\n"; 
			//TMStopTimer(9);
			return 0;
		}

		m.SwapRows(pivotpos,j);
		m.MulRow(j,1/m(j,j));
		double mij;
		for (i=j+1; i<=rows; i++)
		{
			mij = m(i,j);
			if (mij !=0)
			{
				//				m.AddRowVec(j,i,-mij);
				m.AddRowVec(j,i,-mij);
			}
		}
		//#ifdef __sgi
		if (cols >= 500 && j%50==0)
		{
			if (j%500==0) cout << j << flush;
			else cout << "." << flush;
		}
		//#endif

	}
	for (j=cols-1; j>=1; j--)
	{
		for (i=1; i<=j; i++)
		{
			if (m(i,j+1)!=0)
			{
				m.AddRowVec(j+1,i,-m(i,j+1));
			}
		}
	}

	// remove unity-matrix on the right-hand-side
	for (j=1; j<=rows; j++)
	{
		for (i=1; i<=cols; i++)
		{
			(*this)(j,i)=m(j,i+cols);
		}
	}

	//cout << "Invertierte-Matrix: " << m << endl;

	//TMStopTimer(9);

	return 1;
}

//returns 1, if sucsess
int Matrix::Invert2()
{
	if (rows*cols == 0) return 1;
	//	TMStartTimer(9);
#ifndef __QUICKMATH
	release_assert(cols==rows && mat!=NULL);
#endif

	//Insert identity-matrix on left-hand-side
	mystatic Matrix m;
	m.SetSize(rows,cols);
	m.FillWithZeros();

	double mij;
	int i,j,k;
	int n=rows;
	int maxj=0;


	for (j=1; j<=n; j++)
	{
		m(j,j) = 1;
	}

	// Solve lower triangular Matrix

	for (j=1; j<=n; j++)
	{
		double pivot=fabs(Get(j,j));
		int pivotpos=j;
		for (k=j+1; k<=n; k++)
		{
			if (fabs(Get(k,j)) > pivot) {pivotpos=k; pivot=fabs(Get(k,j));}
		}
		if (pivot == 0) 
		{
			//(*global_uo) << "invert:problems with pivot in col " << j << "\n";
			//(*global_uo) << "m=\n" << *this << "\n";
			//TMStopTimer(9);
			return 0;
		}

		maxj=Maximum(pivotpos,maxj);

		m.SwapRows(pivotpos,j);
		SwapRows(pivotpos,j);
		m.MulRow(j,1/Get(j,j));
		MulRow(j,1/Get(j,j));

		for (i=j+1; i<=n; i++)
		{
			mij = Get(i,j);
			if (mij !=0)
			{
				AddRowVec(j,i,-mij,j,n); //j..n
				m.AddRowVec(j,i,-mij,1,maxj); //1..j
			}
		}

	}
	//(*global_uo) << "m1=\n"; PrintMatrix01(*this);

	//backsubstitution, unfortunately this takes most of the time!!!
	for (j=n-1; j>=1; j--)
	{
		for (i=1; i<=j; i++)
		{
			mij = Get(i,j+1);
			if (mij!=0)
			{
				m.AddRowVec(j+1,i,-mij); //1..n
			}
		}
	}
	//(*global_uo) << "m2=\n"; PrintMatrix01(*this);

	*this = m;

	//(*global_uo) << "Jac_inv=\n"; PrintMatrix01(m);

	//TMStopTimer(9);
	return 1;
}

int Matrix::Solve(const Vector& fv, Vector& q)
{
#ifndef __QUICKMATH
	release_assert(fv.GetLen()==rows && cols==rows && cols*rows!=0 && mat!=NULL);
#endif
	mystatic Matrix m;
	mystatic Vector f;

	m=*this;
	//Matrix mstore = *this;
	f=fv;

	if (f.GetLen() == 1) 
	{
		if (mat[0] == 0) {return 0;}
		q(1) = f(1)/mat[0];
		return 1;
	}
	TMStartTimer(11);

	int i,j,k;
	// Solve lower triangular matrix
	for (j=1; j<=cols; j++)
	{
		double pivot=fabs(m(j,j));
		int pivotpos=j;
		for (k=j+1; k<=rows; k++)
		{
			if (fabs(m(k,j)) > pivot) {pivotpos=k; pivot=fabs(m(k,j));}
		}
		if (pivot == 0) 
		{
			(*global_uo) << "Invert: problems with column " << j << "\n"; 
			//(*global_uo) << *this << "\n";
			TMStopTimer(11);
			return 0;
		}

		m.SwapRows(pivotpos,j); Swap(f(pivotpos),f(j));

		f(j)*=1/m(j,j); //diese Zeile muß vor MulRow(j,1/m(j,j)) stehen!!!
		m.MulRow(j,1/m(j,j));

		for (i=j+1; i<=rows; i++)
		{
			if (m(i,j)!=0)
			{
				f(i)+=-m(i,j)*f(j);
				m.AddRowVec(j,i,-m(i,j));
			}
		}
		//    if (j>50 && j%50==1) {lout << "Gauss-Pivot-Solve: " << j << " Zeilen von " << rows << "\n" << flush;}
	}

	//q.SetAll(0); //initialized with zeros
	for (j=rows; j>=1; j--)
	{
		q(j)=f(j);
		double sum = 0; 
		for (int i=j+1; i<=rows; i++)
		{
			sum += m(j,i)*q(i);
		}
		q(j) -= sum;
	}
	//  if (rows>50) {lout << "Gauss-Pivot-Solve: ready\n" << flush;}
	TMStopTimer(11);
	return 1;
}


int Matrix::SolveLapack(Vector& q)
{
	mystatic Matrix m;
	m=*this;
	m.TpYs();
	int info = LapackSolve(rows, 1, &m(1,1), &q(1));
	if (!info) return 1;
	else return 0;
}

int Matrix::InvertLapack()
{
	mystatic Matrix m;
	m=*this;
	int info = LapackInvert(rows, &m(1,1), mat);
	if (!info) return 1;
	else return 0;
}

//$ PG 2013-10-30: function estimates reciprocal of the condition number via LAPACK, returns 0 on success!!!
int Matrix::EstimateReciprocalConditionNumber(double& rcond) const
{
#ifndef __QUICKMATH
	release_assert(cols==rows && cols*rows!=0); //only for sqare matrices
#endif
	
	int info = LapackEstimateReciprocalConditionNumber(rows, GetTp().GetMatPtr(), &rcond);	
	return info;
}

//$ PG 2013-10-30: used for solving an undetermined system (underdetermined: minimum norm solution |q|->min subject to Aq=f; overdetermined: least squares solution |Aq-f|^2->min)
// f (input) right hand side
// q (output) solution
// rcond (input) reciprocal of estimated condition number of A (this matrix), used for recognizing linearly dependendt equations numerically
// rank (output) computed rank of matrix A, depending on rcond.
int Matrix::UndeterminedSystemSolve(const Vector& f, Vector& q, double rcond, int& rank)
{
#ifndef __QUICKMATH
	release_assert(f.GetLen()==rows && cols*rows!=0);
#endif

	double* m = new double[rows*cols];
	int offset = 0;
	for (int i=0; i<cols; i++)
	{
		for (int j=0; j<rows; j++)
		{
			m[offset + j] = this->mat[j*cols + i];
		}
		offset += rows;
	}
	
	int N = max(rows,cols);
	
	q.SetLen(N);
	q.Copy(f,1,1,rows);
	q.SetAll(0.,rows+1,N);   // if rows == N, nothing is done here
	
	int info = LapackUndeterminedSystemSolve(rows, cols, 1, m, &(q(1)), rcond, rank);

	delete[] m;

	if (!info) return 1;
	else return 0;
}

//$ PG 2013-10-2: [ modified function Matrix::LU in order not to use source code, which is not compatible with the HOTINT licence.
//$ PG 2013-10-2: browse the HOTINT CVS repository for older versions of linalg.cpp and linalg.h prior to 2013-10-2 in order to recover former functionality. ]
int Matrix::LU(IVector& indx)
{
	mystatic Matrix m;
	m=(*this).GetTp();
	
	indx.SetLen(min(cols,rows));
	int info = LapackLU(rows, cols, &m(1,1), indx.GetDataPtr());
	*this=m.GetTp();
	
	if (!info)
		return 1;
	
	return 0;
}
//$ PG 2013-10-2: ]


int Matrix::LUBCKSUB(IVector& indx, Vector& b)
{
	//void lubksb(float **a, int n, int *indx, float b[])
	//Solves the set of n linear equations A·X = B. Here a[1..n][1..n] is input, not as the matrix
	//A but rather as its LU decomposition, determined by the routine ludcmp. indx[1..n] is input
	//as the permutation vector returned by ludcmp. b[1..n] is input as the right-hand side vector
	//B, and returns with the solution vector X. a, n, and indx are not modified by this routine
	//and can be left in place for successive calls with different right-hand sides b. This routine takes
	//into account the possibility that b will begin with many zero elements, so it is efficient for use
	//in matrix inversion.

	//(*global_uo) << "indx=" << indx << "\n";
	//(*global_uo) << "b1=" << b << "\n";

	int i,ii=0,ip,j;
	int n=rows;
	double sum;
	for (i=1;i<=n;i++) 
	{ 
		/*
		if (i != indx(i))
			Swap(b(i),b(indx(i)));

		if (ii)
			for (j=ii;j<=i-1;j++) b(i) -= Get(i,j)*b(j);
		else if (b(i)) ii=i; */
		
		ip=indx(i);
		sum=b(ip);
		b(ip)=b(i);
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= Get(i,j)*b(j);
		else if (sum) ii=i; 
		b(i)=sum;	
	}
	//(*global_uo) << "b2=" << b << "\n";

	if (ii) //if whole vector b=0 --> solution: b=0 ...
	{
		for (i=n;i>=1;i--) 
		{ 
			sum=b(i);
			for (j=i+1;j<=n;j++) 
				sum -= Get(i,j)*b(j);
			if (Get(i,i)==0) return 0;
			b(i)=sum/Get(i,i);
		} 
	}
	//(*global_uo) << "b3=" << b << "\n";
	return 1;
}

int Matrix::Bsolve(Vector f, Vector& q)
{
#ifndef __QUICKMATH
	release_assert(f.GetLen()==rows && cols==rows && cols*rows!=0 && mat!=NULL);
#endif
	Matrix m=*this;

	int bw = GetBandWidth();

	int i,j,k;
	//Solve lower triangular matrix
	for (j=1; j<=cols; j++)
	{
		double pivot=fabs(m(j,j));
		int pivotpos=j;
		for (k=j+1; k<=rows; k++)
		{
			if (fabs(m(k,j)) > pivot) {pivotpos=k; pivot=fabs(m(k,j));}
		}
		if (pivot == 0) 
		{
			(*global_uo) << "Invert: problems with column " << j << "\n"; return 0;
		}

		m.SwapRows(pivotpos,j); Swap(f(pivotpos),f(j));

		f(j)*=1/m(j,j); //diese Zeile muß vor MulRow(j,1/m(j,j)) stehen!!!
		m.MulRow(j,1/m(j,j));

		for (i=j+1; i<=rows; i++)
		{
			if (m(i,j)!=0)
			{
				f(i)+=-m(i,j)*f(j);
				m.AddRowVecB(j,i,-m(i,j),bw);
			}
		}
		//    if (j>50 && j%50==1) {lout << "Gauss-Pivot-Solve: " << j << " Zeilen von " << rows << "\n" << flush;}
	}

	q=Vector(rows); //initialized with zeros
	for (j=rows; j>=1; j--)
	{
		q(j)=f(j)-m.GetRowVec(j)*q;
	}
	//  if (rows>50) {lout << "Gauss-Pivot-Solve: ready\n" << flush;}
	return 1;
}

int Matrix::CholeskyBsolve(Vector f, Vector& q)
{
	//(*fout) << "size of B-Matrix = " << rows*cols*8 << " bytes \n";

#ifndef __QUICKMATH
	release_assert(f.GetLen()==rows && cols*rows!=0 && mat!=NULL);
#endif
	Matrix m=*this;

	int i,j,k;

	// generate upper right triangulat matrix (LT)
	for (k=1; k<=rows; k++)
	{
		if (m(k,1)<=0) {(*global_uo) << "Invert: problems with row " << k << "\n"; return 0;}
		double lkkinv=1/sqrt(m(k,1));
		for (i=1; i<=cols; i++)
		{
			m(k,i)*=lkkinv;
		}
		for (i=1; (i<=cols && i+k<=rows); i++)
		{
			for (j=1; j<=cols-i; j++)
			{
				m(i+k,j)-=m(k,i+1)*m(k,j+i);
			}
		}
		//    if (k>200 && k%200==1) {lout << "CholeskySolve: " << k << " Zeilen von " << rows << "\n" << flush;}
	}

	Vector c(rows);  //mot 0n init.

	//L*c=f für c lösen:  (m=LT !!!)
	for (j=1; j<=rows; j++)
	{
		int stoff=0;
		if (j>cols) {stoff=j-cols;}

		c(j)=f(j);
		for (i=1+stoff; i<j; i++)
		{
			c(j)-=c(i)*m(i,j-i+1);
		}
		c(j)*=1/m(j,1);
	}

	q=Vector(rows);  //with 0n init.

	//LT*q=c für q lösen:
	for (j=rows; j>=1; j--)
	{
		int stop=cols;
		if (rows-j+1<=cols) {stop=rows-j+1;}

		q(j)=c(j);
		for (i=2; i<=stop; i++)
		{
			q(j)-=q(j+i-1)*m(j,i);
		}
		q(j)*=1/m(j,1);
	}

	return 1;
}


void Matrix::PrintToMatlab(ostream& os) const
{
	os << "[";
	for (int i=1; i<=Getrows(); i++)
	{
		for (int j=1; j<=Getcols(); j++)
		{
			char str[32];
			sprintf(str,"%1.24g",(*this)(i,j));
			os << str;
			if (j!=Getcols()) {os << " ";}
		}
		if (i!=Getrows()) {os << ";\n";}
	}
	os << "]";
}


void Matrix::PrintToMatlabSparse(ostream& os, double tol)
{
	int mode = 2;

	//count non-zero elements
	int n = 0;
	for (int i=1; i<=Getrows(); i++)
	{
		for (int j=1; j<=Getcols(); j++)
		{
			if (fabs((*this)(i,j)) > tol) n++;
		}
	}

	IVector iv(n);
	IVector jv(n);
	Vector sv(n);

	int cnt = 0;
	for (int i=1; i<=Getrows(); i++)
	{
		for (int j=1; j<=Getcols(); j++)
		{
			if (fabs((*this)(i,j)) > tol)
			{
				cnt++;

				iv(cnt) = i;
				jv(cnt) = j;
				sv(cnt) = (*this)(i,j);
			}
		}
	}

	if (mode == 1)
	{

		os << "sparse([";
		for (int i=1; i<=n; i++)
		{
			char str[32];
			sprintf(str,"%d",iv(i));
			os << str;
			if (i!=n) {os << ";";}
		}
		os << "],[";
		for (int i=1; i<=n; i++)
		{
			char str[32];
			sprintf(str,"%d",jv(i));
			os << str;
			if (i!=n) {os << ";";}
		}
		os << "],[";
		for (int i=1; i<=n; i++)
		{
			char str[32];
			sprintf(str,"%1.16g",sv(i));
			os << str;
			if (i!=n) {os << ";";}
		}
		os << "]";
		os << "," << Getrows() << "," << Getcols();
		os << ")";
	}
	else
	{
		for (int i=1; i<=n; i++)
		{
			char str[32];
			sprintf(str,"%d ",iv(i));
			os << str;
			sprintf(str,"%d ",jv(i));
			os << str;
			sprintf(str,"%1.16g\n",sv(i));
			os << str;
		}
		os << Getrows() << " " << Getcols() << " 0\n"; //highest entry specifies size of matrix ... 0 is added to remaining entries
	}

}

void Matrix::PrintToMaple(ostream& os)
{
	os << "matrix(" << Getrows()
		<< ", " << Getcols() << ", [";

	for (int i=1; i<=Getrows(); i++)
	{
		os << "[";
		for (int j=1; j<=Getcols(); j++)
		{
			char str[32];
			sprintf(str,"%1.16g",(*this)(i,j));
			os << str;
			if (j!=Getcols()) {os << ",";}
		}
		os << "]";
		if (i!=Getrows()) {os << ",\n";}
	}
	os << "]);\n" << endl;
}

mystr Matrix::PrintToMathematica()
{
	mystr os;

	int rows = Getrows();
	int cols = Getcols();

	os += "{";
	for (int i=1; i<=rows; i++)
	{
		os += "{";
		for (int j=1; j<=cols; j++)
		{
			//char str[32];
			//sprintf(str,"%1.16g",(*this)(i,j));
			//os += str;
			os += (*this)(i,j);
			if (j < cols)	os += ",";
		}
		os += mystr("}");
		if (i < rows)	os += mystr(",");
	}
	os += mystr("}");

	os.Replace("e","*^");   // for copy&pasting the string into mathematica

	return os;
}

	// Write Matrix to string to export to C  "double intro[rows][cols] = {.,.,.};"
mystr& Matrix::MakeString(mystr intro)
{
	mystr* buffer = new mystr(intro);
	*buffer +=("["); *buffer += mystr(Getrows()*Getcols()); *buffer +=("]"); 
	*buffer +=(" = \n{ ");
	for (int i=0; i < this->Getrows(); i++)
	{
		for (int j=0; j < this->Getcols(); j++)
		{
			if( j ) *buffer += ", ";
			*buffer += mystr(Get0(i,j));
		}
		if (i<this->Getrows()-1) *buffer += ",\n";
	}
	*buffer += mystr(" };");
	return *buffer;
}

ostream& operator<<(ostream& os, const Matrix& m)
{
	double max=0;
	int i;
	for (i=1; i<=m.Getrows(); i++)
	{
		for (int j=1; j<=m.Getcols(); j++)
		{
			if (fabs(m(i,j))>max) {max=fabs(m(i,j));}
		}
	}
	if (max==0) {max=1;}
	max=(int)log10(max);
	max=pow(10,max);

	os << max << " *\n";

	for (i=1; i<=m.Getrows(); i++)
	{
		os << "[";
		for (int j=1; j<=m.Getcols(); j++)
		{
			char str[32];
			sprintf(str,"% 1.9f",m(i,j)/max);
			os << str;
			if (j!=m.Getcols()) {os << ",";}
		}
		os << "]\n";
	}

	return os;
}

Matrix operator+ (const Matrix& m1, const Matrix& m2)
{
#ifndef __QUICKMATH
	release_assert(m1.cols==m2.cols && m1.rows==m2.rows);
#endif
	Matrix result=m1;
	result+=m2;
	return result;
}

Matrix operator- (const Matrix& m1, const Matrix& m2)
{
#ifndef __QUICKMATH
	release_assert(m1.cols==m2.cols && m1.rows==m2.rows);
#endif
	Matrix result=m1;
	result-=m2;
	return result;
}

Matrix operator* (const Matrix& mat, const double& val)
{
	Matrix result=mat;
	result*=val;
	return result;
}

Matrix operator* (const double& val, const Matrix& mat)
{
	Matrix result=mat;
	result*=val;
	return result;
}

Matrix operator* (const Matrix& m1, const Matrix& m2)
{
#ifndef __QUICKMATH
	release_assert(m1.cols==m2.rows);
#endif

	Matrix res(m1.rows,m2.cols);

	double *rp=res.mat;
	double *mp1=m1.mat;
	double *mp2=m2.mat;
	int rows=res.rows;
	int cols=res.cols;

	for (int i=0;i<rows;i++)
	{
		for (int j=0;j<cols;j++)
		{
			for (int k=0;k<m1.cols;k++)
			{
				rp[i*cols+j]+=mp1[i*m1.cols+k]*mp2[k*cols+j];
			}
		}
	}
	return res;
}

void Mult(const Matrix& m1, const Matrix& m2, Matrix& res)
{
#ifndef __QUICKMATH
	release_assert(m1.cols==m2.rows);
#endif

	res.SetSize(m1.rows,m2.cols);

	double *rp=res.mat;
	double *mp1=m1.mat;
	double *mp2=m2.mat;
	int rows=res.rows;
	int cols=res.cols;

	for (int i=0;i<rows;i++)
	{
		for (int j=0;j<cols;j++)
		{
			rp[i*cols+j] = 0;
			for (int k=0;k<m1.cols;k++)
			{
				rp[i*cols+j] += mp1[i*m1.cols+k] * mp2[k*cols+j];
			}
		}
	}
}

void MultTp(const Matrix& m1, const Matrix& m2, Matrix& res)
{
#ifndef __QUICKMATH
	release_assert(m1.rows==m2.rows);
#endif

	res.SetSize(m1.cols,m2.cols);

	double *rp=res.mat;
	double *mp1=m1.mat;
	double *mp2=m2.mat;
	int rows=res.rows;
	int cols=res.cols;

	for (int i=0;i<rows;i++)
	{
		for (int j=0;j<cols;j++)
		{
			rp[i*cols+j] = 0;
			for (int k=0;k<m1.rows;k++)
			{
				rp[i*cols+j] += mp1[k*m1.cols+i] * mp2[k*cols+j];
			}
		}
	}
}

void MultSym(const Matrix& m1, const Matrix& m2, Matrix& res) //computes res=m1*m2, result is assumed to be symmetric
{
#ifndef __QUICKMATH
	release_assert(m1.cols==m2.rows);
#endif

	res.SetSize(m1.rows,m2.cols);

	double *rp=res.mat;
	double *mp1=m1.mat;
	double *mp2=m2.mat;
	int rows=res.rows;
	int cols=res.cols;

	for (int i=0;i<rows;i++)
	{
		for (int j=0;j<=i;j++)
		{
			rp[i*cols+j] = 0;
			for (int k=0;k<m1.cols;k++)
			{
				rp[i*cols+j] += mp1[i*m1.cols+k] * mp2[k*cols+j];
			}
		}
	}
	for (int i=0;i<rows;i++)
	{
		for (int j=0;j<i;j++)
		{
			rp[j*cols+i] = rp[i*cols+j];
		}
	}
}

void MultSymTp(const Matrix& m1, const Matrix& m2, Matrix& res) //computes res=m1*m2, result is assumed to be symmetric
{
#ifndef __QUICKMATH
	release_assert(m1.rows==m2.rows);
#endif

	res.SetSize(m1.cols,m2.cols);

	double *rp=res.mat;
	double *mp1=m1.mat;
	double *mp2=m2.mat;
	int rows=res.rows;
	int cols=res.cols;

	if (1)
	{
		//+++++++++++++++++++++++++++++++++++++++++
		//new, better caching:
		Matrix m1t = m1; m1t.TpYs();
		Matrix m2t = m2; m2t.TpYs();
		mp1=m1t.mat;
		mp2=m2t.mat;

		int m1tcolsblocks = 8*((int)(m1t.cols/8));

		for (int i=0;i<rows;i++)
		{
			for (int j=0;j<=i;j++)
			{
				rp[i*cols+j] = 0;
				for (int k=0;k<m1tcolsblocks;k+=8)
				{
					rp[i*cols+j] += 
						mp1[i*m1t.cols+k  ] * mp2[j*m2t.cols+k  ]+
            mp1[i*m1t.cols+k+1] * mp2[j*m2t.cols+k+1]+
            mp1[i*m1t.cols+k+2] * mp2[j*m2t.cols+k+2]+
            mp1[i*m1t.cols+k+3] * mp2[j*m2t.cols+k+3]+
            mp1[i*m1t.cols+k+4] * mp2[j*m2t.cols+k+4]+
            mp1[i*m1t.cols+k+5] * mp2[j*m2t.cols+k+5]+
            mp1[i*m1t.cols+k+6] * mp2[j*m2t.cols+k+6]+
            mp1[i*m1t.cols+k+7] * mp2[j*m2t.cols+k+7];
				}
				for (int k=m1tcolsblocks;k<m1t.cols;k++)
				{
					rp[i*cols+j] += mp1[i*m1t.cols+k] * mp2[j*m2t.cols+k];
				}
			}
		}
		//+++++++++++++++++++++++++++++++++++++++++
	}
	else
	{
		for (int i=0;i<rows;i++)
		{
			for (int j=0;j<=i;j++)
			{
				rp[i*cols+j] = 0;
				for (int k=0;k<m1.rows;k++)
				{
					rp[i*cols+j] += mp1[k*m1.cols+i] * mp2[k*cols+j];
				}
			}
		}
	}
	for (int i=0;i<rows;i++)
	{
		for (int j=0;j<i;j++)
		{
			rp[j*cols+i] = rp[i*cols+j];
		}
	}
}

void Mult(const MatrixXD& m1, const Matrix& m2, Matrix& res)
{
#ifndef __QUICKMATH
	release_assert(m1.cols==m2.rows);
#endif

	res.SetSize(m1.rows,m2.cols);

	double *rp=res.mat;
	const double *mp1=m1.mat;
	double *mp2=m2.mat;
	int rows=res.rows;
	int cols=res.cols;

	for (int i=0;i<rows;i++)
	{
		for (int j=0;j<cols;j++)
		{
			rp[i*cols+j] = 0;
			for (int k=0;k<m1.cols;k++)
			{
				rp[i*cols+j] += mp1[i*m1.cols+k] * mp2[k*cols+j];
			}
		}
	}
}

void Mult(const Matrix& m1, const MatrixXD& m2, Matrix& res)
{
#ifndef __QUICKMATH
	release_assert(m1.cols==m2.rows);
#endif

	res.SetSize(m1.rows,m2.cols);

	double *rp=res.mat;
	double *mp1=m1.mat;
	const double *mp2=m2.mat;
	int rows=res.rows;
	int cols=res.cols;

	for (int i=0;i<rows;i++)
	{
		for (int j=0;j<cols;j++)
		{
			rp[i*cols+j] = 0;
			for (int k=0;k<m1.cols;k++)
			{
				rp[i*cols+j] += mp1[i*m1.cols+k] * mp2[k*cols+j];
			}
		}
	}
}

Vector operator* (const Vector& v, const Matrix& m)
{
#ifndef __QUICKMATH
	release_assert(m.rows==v.l);
#endif
	Vector res(m.cols,NOFILL);

	for (int i=1;i<=res.l;i++)
	{
		res(i)=0;
		for (int j=1;j<=v.l;j++)
		{
			res(i)+=m(j,i)*v(j);
		}
		//res(i)=(m.GetColVec(i))*v;
	}
	return res;
}

Vector operator* (const Matrix& m, const Vector& v)
{
#ifndef __QUICKMATH
	release_assert(m.cols==v.l);
#endif

	Vector res(m.rows,NOFILL);

	for (int i=1;i<=res.l;i++)
	{
		res(i)=0;
		for (int j=1;j<=v.l;j++)
		{
			res(i)+=m(i,j)*v(j);
		}
	}
	return res;
}

void Mult(const Matrix& m, const Vector& v, Vector& res) //computes res=m*v
{
#ifndef __QUICKMATH
	release_assert(m.Getcols()==v.GetLen());
#endif
	res.SetLen(m.Getrows());

	double* mm = m.mat; 
	double* vv = v.vec;
	int rl=res.GetLen();
	int vl=v.GetLen();
	for (int i=0; i<rl; i++)
	{
		double val=0;
		double* mr = &mm[i*m.cols];
		for (int j=0; j<vl; j++)
		{
			val+=mr[j]*vv[j];
			//val+=m(i+1,j+1)*v(j+1);
		}
		res.vec[i] = val;
	}
}
void MultTp(const Matrix& m, const Vector& v, Vector& res) //computes res=m.GetTp()*v
{
#ifndef __QUICKMATH
	release_assert(m.Getrows()==v.GetLen());
#endif
	res.SetLen(m.Getcols());

	int rl=res.GetLen();
	int vl=v.GetLen();
	for (int i=0; i<rl; i++)
	{
		double val=0;
		for (int j=0; j<vl; j++)
		{
			val+=m(j+1,i+1)*v(j+1);
		}
		res.vec[i] = val;
	}
}

void Mult(const Matrix& m, const Vector& v, Vector3D& res)
{
	int cols = m.Getcols();
	int rows = m.Getrows();

	#ifndef __QUICKMATH
	release_assert(cols == 3 && rows == v.Length());
	#endif

	for (int i=0; i < rows; i++)
	{
		res.vec[i] = 0;
		for (int j=0; j < cols; j++)
		{
			res.vec[i] += m.mat[i*cols+j]*v.vec[j];
		}
	}
}

void MultTp(const Matrix& m, const Vector& v, Vector3D& res)
{
	int cols = m.Getcols();
	int rows = m.Getrows();

	#ifndef __QUICKMATH
	release_assert(cols == 3 && rows == v.Length());
	#endif

	for (int i=0; i < cols; i++)
	{
		res.vec[i] = 0;
		for (int j=0; j < rows; j++)
		{
			res.vec[i] += m.mat[i+j*cols]*v.vec[j];
		}
	}
}

void Mult(const Vector& v, const Matrix& m, Vector& res) //computes res=v*m
{
#ifndef __QUICKMATH
	release_assert(m.Getrows()==v.GetLen());
#endif
	res.SetLen(m.Getcols());

	int rl=res.GetLen();
	int vl=v.GetLen();
	for (int i=0; i<rl; i++)
	{
		double val=0;
		for (int j=0; j<vl; j++)
		{
			val+=m(j+1,i+1)*v(j+1);
		}
		res.vec[i] = val;
	}
}


void MultBW(const Matrix& m, const Vector& v, Vector& res, int bw) //computes res=m*v
{
	TMStartTimer(12);
#ifndef __QUICKMATH
	release_assert(m.Getcols()==v.GetLen());
#endif
	res.SetLen(m.Getrows());

	if (bw==0) bw = m.cols+m.rows;
	double* mm = m.mat; 
	double* vv = v.vec;
	int rl=res.GetLen();
	int vl=v.GetLen();
	for (int i=0; i<rl; i++)
	{
		double val=0;
		int bwoff = 0;
		int jend = vl;
		if (i > bw) bwoff = i-bw; 
		if (i+bw < vl) jend = i+bw;
		double* mr = &mm[i*m.cols];
		for (int j=bwoff; j<jend; j++)
		{
			val+=mr[j]*vv[j];
			//val+=m(i+1,j+1)*v(j+1);
		}
		res.vec[i] = val;
	}
	TMStopTimer(12);
}

void Mult(const Vector& v1, const Vector& v2, Matrix& res) //computes res=v1*v2^T, where v1 and v2 are interpreted as a column vectors
{
	double rows = v1.GetLen();
	double cols = v2.GetLen();
	res.SetSize(rows, cols);

	double* mm = res.mat;
	for (int i=0; i<rows; i++)
	{
		double rowfactor = v1.vec[i];
		int off = i*cols;
		for (int j=0; j<cols; j++)
		{
			mm[off+j] = rowfactor*v2.vec[j];
		}
	}
}

void Mult(const Matrix& m, const Vector3D& v, Vector& res) //computes res=m*v
{
#ifndef __QUICKMATH
	release_assert(m.Getcols()==v.GetLen());
#endif
	res.SetLen(m.Getrows());

	double* mm = m.mat;
	for (int i=0; i<res.GetLen(); i++)
	{
		int off = i*m.cols;
		res.vec[i] = mm[off]*v.vec[0]+mm[off+1]*v.vec[1]+mm[off+2]*v.vec[2];
	}
}

void Mult(const Matrix& m, const Vector2D& v, Vector& res) //computes res=m*v
{
#ifndef __QUICKMATH
	release_assert(m.Getcols()==v.GetLen());
#endif
	res.SetLen(m.Getrows());

	double* mm = m.mat;
	for (int i=0; i<res.GetLen(); i++)
	{
		int off = i*m.cols;
		res.vec[i] = mm[off]*v.vec[0]+mm[off+1]*v.vec[1];
	}
}

Vector Matrix::Bmult(Vector& v) const
{
#ifndef __QUICKMATH
	release_assert(v.GetLen()==rows);
#endif
	Vector res(rows); //initialized with zeros

	for (int i=1;i<=rows;i++)
	{
		int start=i-cols;
		int end=i+cols-1;
		if (start < 1) {start = 1;}
		if (end > rows) {end = rows;}
		for (int j=start;j<=end;j++)
		{
			res(i)+=v(j)*BGet(i,j);
		}
	}
	return res;
}

int Matrix::BCG(Vector f, Vector& q, const Vector& qs)
{
	Vector x = qs;
	Vector r(f.GetLen());
	Vector rq(f.GetLen());
	Vector p(f.GetLen());
	Vector pq(f.GetLen());
	Vector Ap(f.GetLen());
	Vector App(f.GetLen());

	double alpha, beta, rqr;

	//initial conditions:
	r = f-(*this)*x;
	rq = r;
	p = r;
	pq = rq;

	//error condition
	Vector err(f.GetLen());
	int end = 0;
	int it = 0;

	int maxit = f.GetLen();

	while(!end && it < maxit)
	{
		Ap = ((*this)*p);

		alpha = (rq*r)/(pq*Ap);
		err = alpha * p;
		x = x + err;

		rqr = rq*r;
		r = r - alpha * Ap;
		rq = rq - alpha * GetTp() * p;

		beta = (rq*r)/rqr;
		p = r + beta * p;
		pq = rq + beta * pq;

		it++;
		if (it % 100 == 0) {(*global_uo) << "    CG:" << it << " its, err=" << fabs(err.MaxNorm()) << "\n";}
		if (fabs(err.MaxNorm()) < 1E-10) {end = 1;}
	}
	(*global_uo) << "CG: " << it << " its, err=" << fabs(err.MaxNorm()) << "\n";

	q = x;
	if (it < maxit) {return 1;}
	else {return 0;}
}


double eigenvalue_err;
int Matrix::Eigenvalues(Vector& ev)
{
	eigenvalue_err = 0;
	//evout << "Compute Eigenvalue:" << endl;
	int end = 0;
	double tol = 1e-5;
	const int maxit = 500;
	int i = 0;

	Matrix M(*this);
	Matrix Q,R;
	Vector evstart = M.GetDiag();
	double err;
	double relval = evstart.GetNorm();
	if (relval == 0) {relval = 1;}

	while(i < maxit && !end)
	{
		i++;
		M.QRDecomposition(Q,R);
		M = R*Q;
		ev = M.GetDiag();
		err = (ev-evstart).GetNorm()/relval;
		if (err < tol) {end = 1;}
		evstart = ev;

		eigenvalue_err = err;
		//evout << "err=" << err << endl;
	}

	//::cout << "EV-iterations=" << i << endl;
	int retv=1;
	if (i >= maxit) 
	{
		retv=0; 
		//evout << "ERROR: Eigenvalues not converged!!!" << endl;
	}


	return retv;

};

int Matrix::EigenvaluesLapack(Vector& lami) const
{
	mystatic Matrix m;
	m=*this;
	lami.SetLen(Getrows());
	mystatic Vector work;
	work.SetLen(4*Getrows());
	int info = LapackEVPSPD(Getrows(), &(m(1,1)), &(lami(1)), &(work(1)), work.Length());
	if (info) { (*global_uo) << "Matrix::EigenvaluesLapack: Eigenvalues could not be computed, Lapack error message " << info << "\n"; return 0;}
	else return 1;
}

int Matrix::QRDecomposition(Matrix& Q, Matrix& R)
{
#ifndef __QUICKMATH
	release_assert(IsSquare());
#endif
	if (Q.Getcols()!=cols)
	{
		Q=Matrix(cols,cols);
	}
	else
	{
		Q.FillWithZeros();
	}
	if (R.Getcols()!=cols)
	{
		R = Matrix(cols,cols);
	}
	else
	{
		R.FillWithZeros();
	}
	mystatic Matrix Mu;

	if (Mu.Getcols()!=cols)
	{
		Mu = Matrix(cols,cols);
	}
	else
	{
		Mu.FillWithZeros();
	}

	//Vector v;

	int i,j;

	//Mu.SetColVec(this->GetColVec(1),1);
	Mu.SetColVec(*this,1,1);

	//Q.SetColVec(Mu.GetColVec(1),1);
	Q.SetColVec(Mu,1,1);

	for(i=1; i<= cols; i++)
	{
		R(i,i) = 1;
	}

	for (j=2; j<=cols; j++)
	{
		//Vector sumv(cols);

		//Mu.SetColVec(this->GetColVec(j),j);
		Mu.SetColVec(*this,j,j);

		for (i=1; i<=j-1; i++)
		{
			//double u2=Mu.GetColVec(i)*Mu.GetColVec(i);
			double u2=Mu.MultCol(Mu,i,i);

			if (u2 == 0) 
			{
				//cout << "ERROR: QR-Decomposition, Matrix singular!!!" << endl; 
				return 0;
			}
			R(j,i)=0;
			//R(i,j)=(1./u2)*Mu.GetColVec(i)*this->GetColVec(j);
			R(i,j)=(1./u2)*Mu.MultCol(*this,j,i);

			//Mu.SetColVec(Mu.GetColVec(j)-R(i,j)*Mu.GetColVec(i),j);
			Mu.AddColVec(i,j,-R(i,j));

		}
		for (i=1; i <= cols; i++)
		{
			Q(i,j) = Mu(i,j);
		}
	}
	Vector Dm(cols);
	for (i=1; i <= cols; i++)
	{
		//Dm(i,i) = sqrt(Mu.GetColVec(i)*Mu.GetColVec(i));
		Dm(i) = sqrt(Mu.MultCol(Mu,i,i));
		R.MulRow(i,Dm(i));
		Q.MulCol(i,1./Dm(i));
	}
	return 1;

}


// Evaluate piecewise bilinear function on unit square [-1,1]^2 at point p,
// Matrix entries specify function values
//   [ 1          2         3   ... width ]             [        y=1         ]
//   [ w+1        w+2       w+3 ... 2w    ]   ----->    [x=-1             x=1]
//
//   [ (h-1)w+1  (h-1)w+1  (h-1)w+2  ...  hw]           [        y=-1        ]
double Matrix::Interpolate(const Vector2D& ploc) const
{
	int imin, jmin;
	double hy = 2./(Getrows()-1), hx = 2./(Getcols()-1.);
	// lower indices of cell containing ploc
	imin = int((1.-ploc.Y()) / hy ) +1; // int((ploc.Y()+1.) / hy ) +1;
	jmin = int((ploc.X()+1.) / hx ) +1;
	if (imin >= Getrows()) imin = Getrows()-1;
	if (jmin >= Getcols()) jmin = Getcols()-1;
	// "mesh sizes"
	double xmin = hx * (jmin-1)-1, ymin = 1-hy*(imin);
	double val = (*this)(imin, jmin)*(hx-ploc.X()+xmin)*(ploc.Y()-ymin)/(hx*hy) + 
								(*this)(imin, jmin+1)*(ploc.X()-xmin)*(ploc.Y()-ymin)/(hx*hy) +
								(*this)(imin+1, jmin)*(hx-ploc.X()+xmin)*(hy-ploc.Y()+ymin)/(hx*hy) +
								(*this)(imin+1, jmin+1)*(ploc.X()-xmin)*(hy-ploc.Y()+ymin)/(hx*hy);
	return val;
}

// Evaluate piecewise bilinear function on unit square [-1,1]^2 at x-coord xloc, 
// Vector entries(i) correspond to y-coordinates
// Matrix entries specify function values
//   [ 1          2         3   ... width ]             [        y=1         ]
//   [ w+1        w+2       w+3 ... 2w    ]   ----->    [x=-1             x=1]
//
//   [ (h-1)w+1  (h-1)w+1  (h-1)w+2  ...  hw]           [        y=-1        ]
void Matrix::InterpolateX(const double xloc, Vector& interpol_vector) const
{
	interpol_vector.SetLen(Getrows());
	int  jmin;
	double hy = 2./(Getrows()-1), hx = 2./(Getcols()-1.);
	// lower indices of cell containing ploc
	jmin = int((xloc+1.) / hx ) +1;
	if (jmin >= Getcols()) jmin = Getcols()-1;
	// "mesh sizes"
	double xmin = hx * (jmin-1)-1;
	for (int i=1; i<=Getrows(); i++)
	{
	interpol_vector(i) = (*this)(i, jmin)*(hx-xloc+xmin)/(hx) + 
								(*this)(i, jmin+1)*(xloc-xmin)/(hx);
	}
}

void PrintMatrix01(const Matrix& m)
{
	for (int i=1; i<=m.Getrows(); i++)
	{
		mystr str="";
		for (int j=1; j<=m.Getcols(); j++)
		{
			if (fabs(m.Get(i,j)) <= 1e-14) str+= "0"; 
			else str += "1";
		}
		str+="\n";
		(*global_uo) << str;
	}
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++              SPARSE MATRIX        ++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void SparseMatrix::GenMat(int lalloc, int rowsi)
{
#ifdef gen_sparsemat
	GenSparsemat();
#endif
	dnull = 0;
	mat = new double[lalloc];
	colind = new int[lalloc];
	rowind = new int[rowsi+1]; //for convention with Harwell Boeing format
	rowlen = new int[rowsi];

	//(*global_uo) << "genmat=" << lalloc*12./1.e6 << "MB\n";
}

int reallocate_sparsematrix_cnt = 0;
void SparseMatrix::ReAllocate()
{
	if (reallocate_sparsematrix_cnt%100 == 0)	(*global_uo) << "Reallocate sparse matrix " << reallocate_sparsematrix_cnt++ << "times ==> switch to SolverOptions.Static.experimental_sparse_jacobian\n";

	//calculate space needed
	int incelems = 0;
	const int incfact = 2;
	for (int row=0; row < rows; row++)
	{
		int rlmax = RowLenMax0(row);
		//(*global_uo) << rlmax << ",";
		if (rlmax/incfact < rowlen[row])
		{
			int addtorow = rowlen[row]*incfact-rlmax;
			incelems += addtorow;
		}
	}
	//(*global_uo) << "]\n";

	int newsize = lalloc+incelems;
	//(*global_uo) << "Reallocate=" << newsize << "\n";

	if (newsize==0) {return; }

	double* nmat = new double[newsize];
	int* ncolind = new int[newsize];
	int* nrowind = new int[rows+1]; //for convention with Harwell Boeing format
	int* nrowlen = new int[rows];
	//(*global_uo) << "rnewmat=" << newsize*12./1.e6 << "MB\n";

	if (mat!=NULL && lalloc!=0)
	{
		//generate new data:
		nrowind[0] = 0;
		for (int row=0; row < rows; row++)
		{
			int addtorow = 0;
			int rlmax = RowLenMax0(row); //accesses lalloc!!!
			if (rlmax/incfact < rowlen[row])
			{
				addtorow = rowlen[row]*incfact-rlmax;
			}

			int newmaxlen = rlmax+addtorow;
			if (row < rows-1) nrowind[row+1] = nrowind[row]+newmaxlen;
			nrowlen[row] = rowlen[row];

			for (int i = 0; i < rowlen[row]; i++)
			{
				nmat[nrowind[row]+i] = mat[rowind[row]+i];
				ncolind[nrowind[row]+i] = colind[rowind[row]+i];
			}

			int nrlmax; if (row < rows-1) nrlmax = nrowind[row+1]-nrowind[row];	else nrlmax = lalloc-nrowind[row];
			for (int i = rowlen[row]; i < nrlmax; i++)
			{
				nmat[nrowind[row]+i] = 0; //not necessary
				ncolind[nrowind[row]+i] = sparseMatrixVoid; //not necessary
			}
		}
	}

	Destroy();

	lalloc = newsize;
	mat = nmat;
	colind = ncolind;
	rowind = nrowind;
	rowlen = nrowlen;
}

void SparseMatrix::GetMatrix(Matrix& m) const
{
	m.SetSize(rows,cols);
	m.FillWithZeros();
	for (int i = 0; i < rows; i++)
	{
		//(*global_uo) << "row=" << i << ": ";
		for (int j = 0; j < rowlen[i]; j++)
		{
			//(*global_uo) << "(" << colind[rowind[i]+j] << "," << mat[rowind[i]+j] << "),";

			m.Elem0(i,colind[rowind[i]+j]) = mat[rowind[i]+j];
		}
		//(*global_uo) << "\n";
	}
}


void SparseMatrix::CopyFrom(const Matrix& mati)
{
	int nrows = mati.Getrows();
	int ncols = mati.Getcols();

	int neededmem = 0;
	for (int i = 0; i < nrows*ncols; i++)
	{
		if (fabs(mati.GetMatPtr()[i]) > sparsetolzero)
		{
			neededmem++;
		}
	}

	if (neededmem > lalloc || nrows > rows)
	{
		Destroy();

		lalloc = neededmem;
		GenMat(lalloc, nrows);
	}
	rows = mati.Getrows();
	cols = mati.Getcols();

	int cnt = 0;
	for (int i = 0; i < rows; i++)
	{
		rowind[i] = cnt;
		int colcnt = 0;
		for (int j = 0; j < cols; j++)
		{
			if (fabs(mati.GetMatPtr()[i*cols+j]) > sparsetolzero)
			{
				colind[cnt] = j;
				colcnt++;
				mat[cnt++] = mati.GetMatPtr()[i*cols+j];
			}
		}
		rowlen[i] = colcnt;
	}
}

void SparseMatrix::CopyFrom(const Matrix& mati, int row1, int col1, int row2, int col2)
{
	int nrows = row2-row1+1;
	int ncols = col2-col1+1;

	int neededmem = 0;
	for (int i = row1; i <= row2; i++)
	{
		for (int j = col1; j <= col2; j++)
		{
			if (fabs(mati(i,j)) > sparsetolzero)
			{
				neededmem++;
			}
		}
	}

	if (neededmem > lalloc || nrows > rows)
	{
		Destroy();

		lalloc = neededmem;
		GenMat(lalloc, nrows);
	}
	rows = nrows;
	cols = ncols;

	int cnt = 0;
	for (int i = 0; i < rows; i++)
	{
		rowind[i] = cnt;
		int colcnt = 0;
		for (int j = 0; j < cols; j++)
		{
			if (fabs(mati(i+row1,j+col1)) > sparsetolzero)
			{
				colind[cnt] = j;
				colcnt++;
				mat[cnt++] = mati(i+row1,j+col1);
			}
		}
		rowlen[i] = colcnt;
	}
}

//row1, col1, row2, col2 are 1-based!!!
void SparseMatrix::CopyFrom(const SparseMatrix& mati, int row1, int col1, int row2, int col2)
{

	//(*global_uo) << "sparsematrix::copyfrom(sparsematrix) does not work yet!\n";

	int nrows = row2-row1+1;
	int ncols = col2-col1+1;

	int neededmem = 0;
	for (int i = row1-1; i <= row2-1; i++)
	{
		for (int j = 0; j < mati.rowlen[i]; j++)
		{
			int c = mati.colind[mati.rowind[i]+j];
			if (c >= col1-1 && c <= col2-1) neededmem++;
		}
	}

	if (neededmem > lalloc || nrows > rows)
	{
		Destroy();

		lalloc = neededmem;
		GenMat(lalloc, nrows);
	}

	rows = nrows;
	cols = ncols;

	int cnt = 0;
	for (int i = 0; i < rows; i++)
	{
		rowind[i] = cnt;
		int matirow = i+row1-1;
		int rind = mati.rowind[matirow];
		int k = 0;
		while (k < mati.rowlen[matirow] && mati.colind[k+rind] < col1-1) k++;

		int j = 0;
		while (mati.colind[k+j+rind] <= col2-1 && (k+j) < mati.rowlen[matirow])
		{
			colind[j+rowind[i]] = mati.colind[rind+j+k]-col1+1;
			mat[j+rowind[i]] = mati.mat[rind+j+k];
			cnt++; j++;
		}
		rowlen[i] = j;
	}
}

//generate new matrix from sparse matrix, only from row1 to row2, col1 to col2 (incl.), resort with vector r
void SparseMatrix::CopyFrom(const SparseMatrix& mati, int row1, int col1, int row2, int col2, const IVector& r)
{
	int nrows = row2-row1+1;
	int ncols = col2-col1+1;

	mystatic IVector invr;
	invr.SetLen(r.Length());
	invr.SetAll(-1);
	for (int i=1; i <= r.Length(); i++) invr(r(i)) = i;

	int neededmem = 0;


	for (int i = 0; i < nrows; i++)
	{
		int matirow = r(i+row1)-1;
		int rind = mati.rowind[matirow];

		//read out necessary entries and transform by r:
		for (int j = 0; j < mati.rowlen[matirow]; j++)
		{
			int c = mati.colind[rind+j]; //0-based
			int cr = invr(c+1); //1-based
			if (cr >= col1 && cr <= col2) neededmem++;
		}
	}


	if (neededmem > lalloc || nrows > rows)
	{
		Destroy();

		lalloc = neededmem;
		GenMat(lalloc, nrows);
	}

	rows = nrows;
	cols = ncols;



	mystatic TArray<double> rsort; //resorted values
	mystatic TArray<int> rsortind; //resorted column indices 0-based for destination


	//sm    1  2  3  4  5  6  7  8  9
	//rind  1     3        6     8
	//mat   11    12       13    14
	//r     1  2  8  3  4  5  9  6  7
	//invr  1  2  4  5  6  8  9  3  7
	//m     11    14 12          13    
	//

	int c, cr;
	int cnt = 0;

	//(*global_uo) << "lalloc=" << lalloc << ", nrows=" << nrows << "\n";

	//mat[i++] = m(r(row),r(col));


	//TMStartTimer(24);
	for (int i = 0; i < rows; i++)
	{
		rsortind.SetLen(0);
		rsort.SetLen(0); 

		int matirow = r(i+row1)-1;
		int rind = mati.rowind[matirow];

		//read out necessary entries and transform by r:
		for (int j = 0; j < mati.rowlen[matirow]; j++)
		{
			c = mati.colind[rind+j]; //0-based
			cr = invr(c+1); //1-based
			if (cr >= col1 && cr <= col2)
			{
				rsortind.Add(invr(c+1)-col1); //0-based index for destination
				//runsortind.Add(c-col1+1); //0-based index for destination
				rsort.Add(mati.mat[rind+j]); //entry
			}
		}
		//sort entries
		Quicksort(rsortind, rsort);

		//fill in entries in i-th row
		rowind[i] = cnt;
		int cnt2 = 0;

		/*
		(*global_uo) << "row=" << i << ", rowind=" << rowind[i] << "\n";
		//(*global_uo) << "rsortind_len=(" << rsortind.Length() << "=" << rsort.Length() << ")\n";
		(*global_uo) << "rsortind=" << rsortind << "\n";
		(*global_uo) << "rsort   =" << rsort << "\n";
		*/


		for (int j=0; j < rsortind.Length(); j++)
		{
			if (rsort.Get0(j) != 0)
			{
				colind[cnt2+rowind[i]] = rsortind.Get0(j);
				//(*global_uo) << "(" << cnt2+rowind[i] << "," << colind[cnt2+rowind[i]] << ") ";

				mat[cnt2+rowind[i]] = rsort.Get0(j);
				cnt2++;
			}
		}
		//(*global_uo) << "\n";

		rowlen[i] = cnt2;
		cnt+=cnt2;
	}
	//TMStopTimer(24);
}

void SparseMatrix::GetBandWidths(int& lb, int& ub) const 
{
	lb = 0;
	ub = 0;
	for (int i = 0; i < rows; i++)
	{
		//(*global_uo) << "\nrow " << i << ": ";
		//(*global_uo) << "(len=" << rowlen[i] << ") ";

		for (int j = 0; j < rowlen[i]; j++)
		{

			//(*global_uo) << colind[rowind[i]+j] << ", ";
			if (mat[rowind[i]+j] != 0)
			{
				if (colind[rowind[i]+j]>i) {ub = Maximum(ub,colind[rowind[i]+j]-i);}
				if (colind[rowind[i]+j]<i) {lb = Maximum(lb,i-colind[rowind[i]+j]);}
			}
		}

		//(*global_uo) << "\n";
	}
}

void SparseMatrix::ComputeTranspose(SparseMatrix& res) const //computes res=m1^T
{
	int n = CountEntries();
	int rows = Getrows();
	int cols = Getcols();

	int rowlen = n/cols + 1;
	res.SetSize(cols, rows, rowlen);
	if (rowlen*cols < n) (*global_uo) << "Warning: SparseMatrix::ComputeTranspose allocated insufficient memory!\n";

	//IVector ncolind(cols);
	TArray<int> ncolind;
	ncolind.SetLen(cols);
	int i;

	for (i=1; i<=cols; i++)
	{
		ncolind(i) = 0;
	}

	for (i=1; i<=rows; i++)
	{
		for (int j=1; j <= RowLen(i); j++)
		{
			int ci = GetRowColind(i,j);
			ncolind(ci+1)++;
		}
	}

	
	int actrowind = 0;
	for (i=1; i<=cols; i++)
	{
		res.rowlen[i-1] = 0;
		res.rowind[i-1] = actrowind;
		actrowind += ncolind(i);
	}

	for (i=1; i<=rows; i++) //run over rows of m1 ==> columns of new matrix
	{
		for (int j=1; j <= RowLen(i); j++)
		{
			int resrow = GetRowColind(i,j);
			res.colind[res.rowind[resrow] + res.rowlen[resrow]] = i-1;
			res.mat[res.rowind[resrow] + res.rowlen[resrow]] = GetRowEntry(i,j);
			res.rowlen[resrow]++;
		}
	}
	

}

double MultRow(const SparseMatrix& m1, const SparseMatrix& m2, int row1, int row2) //mult row1 of m1 times row2 of m2
{
	double val = 0;
	int m2j = 1;
	int m1rowlen = m1.RowLen(row1);
	int m2rowlen = m2.RowLen(row2);

	for (int m1j=1; m1j <= m1rowlen; m1j++)
	{
		int m1ci = m1.GetRowColind(row1, m1j);
		int m2ci = m2.GetRowColind(row2, m2j);
		int end = 0;
		while (m2j < m2rowlen && m2ci < m1ci && !end)
		{
			m2j++;
			m2ci = m2.GetRowColind(row2, m2j);
			if (m2ci > m1ci) 
			{
				m2j--;
				m2ci = m2.GetRowColind(row2, m2j);
				end = 1;
			}
		}
		if (m1ci == m2ci) val += m1.GetRowEntry(row1,m1j) * m2.GetRowEntry(row2, m2j);
	}
	return val;
}

void Mult(const SparseMatrix& m1, const SparseMatrix& m2, SparseMatrix& res) //computes res=m1*m2
{
	release_assert(m2.Getrows() == m1.Getcols());

	SparseMatrix m2T;
	m2.ComputeTranspose(m2T);

	int minrowsize = m1.GetLAlloc()/m1.Getrows() + m2T.GetLAlloc()/m2.Getrows();
	res.SetSize(m1.Getrows(), m2T.Getcols(), minrowsize);

	int n = 0;
	int i, j;
	for (i=1; i<=m1.Getrows(); i++)
	{
		for (j=1; j<=m2T.Getrows(); j++)
		{
			double val = 0;

			val = MultRow(m1, m2T, i, j);

			/*
			int m2j = 1;
			int m2rowlen = m2T.RowLen(i);
			for (int m1j=1; m1j <= m1.RowLen(i); m1j++)
			{
				int m1ci = m1.GetRowColind(i, m1j);
				int m2ci = m2T.GetRowColind(j, m2j);
				int end = 0;
				while (m2j < m2rowlen && m2ci < m1ci && !end)
				{
					m2j++;
					m2ci = m2T.GetRowColind(j, m2j);
					if (m2ci > m1ci) 
					{
						m2j--;
						m2ci = m2T.GetRowColind(j, m2j);
						end = 1;
					}
				}
				if (m1ci == m2ci) val += m1.GetRowEntry(i,m1j) * m2T.GetRowEntry(j, m2j);
			}
			*/
			if (val != 0)
			{
				res(i,j) = val;
			}
		}
	}

}

void Mult(const SparseMatrix& m, const Vector& v, Vector& res)
{
#ifndef __QUICKMATH
	release_assert(m.Getcols()==v.GetLen());
#endif
	TMStartTimer(12);

	//if (res.GetLen() != m.rows) 
	res.SetLen(m.rows);

	for(int row = 0; row < m.rows; row++)
	{
		res.vec[row] = 0;
		for(int i = 0; i < m.rowlen[row]; i++)
		{
			res.vec[row] += m.mat[m.rowind[row]+i]*v.vec[m.colind[m.rowind[row]+i]];
		}
	}
	TMStopTimer(12);
}

Matrix operator* (const SparseMatrix& m1, const Matrix& m2)
{
#ifndef __QUICKMATH
	release_assert(m1.cols==m2.rows);
#endif

	Matrix res(m1.rows,m2.cols);


	for(int row = 0; row < res.rows; row++)
	{
		for(int col = 0; col < res.cols; col++)
		{
			//res.mat[row*res.cols+col] = 0; //matrix is already filled with zeros
			for(int i = 0; i < m1.rowlen[row]; i++)
			{
				res.mat[row*res.cols+col] += m1.mat[m1.rowind[row]+i]*m2.mat[m1.colind[m1.rowind[row]+i]*m2.cols+col];
			}
		}
	}
	return res;
}

void Mult(const SparseMatrix& m1, const Matrix& m2, Matrix& res) //computes res=m1*m2
{
#ifndef __QUICKMATH
	release_assert(m1.cols==m2.rows);
#endif

	res.SetSize(m1.rows,m2.cols);

	for(int row = 0; row < res.rows; row++)
	{
		for(int col = 0; col < res.cols; col++)
		{
			res.mat[row*res.cols+col] = 0;
			for(int i = 0; i < m1.rowlen[row]; i++)
			{
				res.mat[row*res.cols+col] += m1.mat[m1.rowind[row]+i]*m2.mat[m1.colind[m1.rowind[row]+i]*m2.cols+col];
			}
		}
	}
}

//alex for double indices/not sorted indices rowref/colref
void SparseMatrix::AddMatrix(const TArray<int>& rowref, const TArray<int>& colref, int lrow, int lcol, const Matrix& m)
{
	int k;
	if (colref.Length() == 0) return;
	int rlmax, rind;
	//(*global_uo) << "colref=" << colref << "\n";
	for (int i=1; i <= lrow; i++)
	{
		int row = rowref(i)-1;
		k = 0; //temporary variable for searching index of actual row
		int j = 1;
#if _DEBUG
		int j_old = j-1;
#endif
		while(j <= lcol)
		{
#if _DEBUG
			assert (j > j_old  || k == 0);
			j_old = j;
#endif
			if (m(i,j) != 0)
			{
				rlmax = RowLenMax0(row);
				rind = rowind[row];

				if (k!=0 && k < rowlen[row] && rowlen[row] != 0 && colind[k+rind] > colref(j)-1) k = 0;

				while (k < rowlen[row] && colind[k+rind] < colref(j)-1) k++;

				if ((k < rowlen[row]) && (colref(j)-1 == colind[k+rind]))
				{
					mat[k+rind] += m(i,j);
					j++; 
					k++; //$ AD: NOTE: this line may cause infinite loops when the DOF of a new entry occurs multiple times in the matrix m resp. colref array ("double indices")
				}
				else
				{
					//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					if (k == rowlen[row] && k < rlmax) 
					{
						//add new entry, because index not found

						//int sk = 0;
						int lastcol = -1;
						if (k > 0) lastcol = colind[k-1+rind];

						while(j <= lcol && k < rlmax && lastcol < (colref(j)-1)) 
//$ AD: NOTE: further entries are likely to be missing as well
//            bailout 1: k < rlmax --> array full
//						bailout 2: lastcol < (colref(j)-1) --> in case colref is not ascending order
						{
							if (m(i,j) != 0)
							{
								//sk = 1;
								rowlen[row]++;
								mat[k+rind] = m(i,j);
								colind[k+rind] = colref(j)-1;
								lastcol = colref(j)-1;
								k++;
							}
							j++;
						}
						//if (sk) k--;
						//if (j <= lcol && lastcol > colref(j)-1) k = 0; //old
						if (j <= lcol && lastcol >= colref(j)-1) k = 0; //new proposal JG 4.4.2013
					}
					//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					else //other case, e.g. no space for entry
					{
						(*this)(row+1,colref(j)) += m(i,j); // += also possible ...
						j++;
					}
				}
			}
			else
			{
				j++;
			}
		}
	}
}

//Take Submatrix sm at (sr,sc), multiply with x and add it to local matrix at (r,c), size: (nr x nc)
void SparseMatrix::AddSubmatrix(const Matrix& sm, int sr, int sc, int r, int c, int nr, int nc, double x)
{
	if (nr==0 || nc ==0) {return;}
	int k;
	int rlmax, rind;
	int offr = r-sr;
	int offc = c-sc;
	for (int i=sr; i < sr+nr; i++)
	{
		int row = i-1+offr;
		k = 0;
		int j = sc;

		while (j < sc+nc)
		{
			if (sm(i,j) != 0)
			{
				//(*this)(row+1,j+offc) += x*sm(i,j);
				//j++;
				rlmax = RowLenMax0(row);
				rind = rowind[row];

				while (k < rowlen[row] && colind[k+rind] < j-1+offc) k++;

				if ((k < rowlen[row]) && (j-1+offc == colind[k+rind]))
				{
					mat[k+rind] += x*sm(i,j);
					j++; k++;
				}
				else
				{
					if (k == rowlen[row] && k < rlmax) 
					{
						while(j < sc+nc && k < rlmax)
						{
							rowlen[row]++;
							mat[k+rind] = x*sm(i,j);
							colind[k+rind] = j-1+offc;
							j++; k++;
						}
					}
					else
					{
						(*this)(row+1,j+offc) += x*sm(i,j);
						j++;
					}
				}

			}
			else j++;
		}
	}
}

//Take Submatrix sm at (sr,sc), multiply with x and add it to local matrix at (r,c), size: (nr x nc)
void SparseMatrix::AddSubmatrix(const SparseMatrix& sm, int sr, int sc, int r, int c, int nr, int nc, double x)
{
	if (nr==0 || nc ==0) {return;}
	int k;
	int rlmax, rind;

	int offr = r-sr;
	int offc = c-sc;

	//int failed = 0;

	for (int i=sr; i < sr+nr; i++)
	{
		int row = i-1+offr;
		k = 0;
		int scol = 0; //sparse column
		int ksm = 0;
		int rindsm = sm.rowind[i-1];
		int rowlensm = sm.rowlen[i-1];

		//find starting index in sm-matrix
		while(sm.colind[ksm+rindsm] < sc-1 && ksm < rowlensm)
		{
			ksm++;
		}
		//do until end index sc+nc reached
		while(sm.colind[ksm+rindsm] < sc+nc-1 && ksm < rowlensm)
		{
			//(*this)(row+1,sm.colind[ksm+rindsm]+1+offc) += x*sm.mat[ksm+rindsm]; 
			//ksm++;

			rlmax = RowLenMax0(row);
			rind = rowind[row];

			int actcol = sm.colind[ksm+rindsm]+offc;
			//find starting index in this matrix
			while (k < rowlen[row] && colind[k+rind] < actcol) k++;

			if ((k < rowlen[row]) && (colind[k+rind] == actcol))
			{
				mat[k+rind] += x*sm.mat[ksm+rindsm];
				k++; ksm++;
			}
			else
			{
				if (k == rowlen[row] && k < rlmax && sm.colind[ksm+rindsm] < sc+nc-1 && ksm < rowlensm) 
				{
					while(sm.colind[ksm+rindsm] < sc+nc-1 && k < rlmax && ksm < rowlensm)
					{
						rowlen[row]++;
						mat[k+rind] = x*sm.mat[ksm+rindsm];
						colind[k+rind] = sm.colind[ksm+rindsm]+offc;
						k++; ksm++;
					}
				}
				else
				{
					//failed++;
					//rescue strategy, should not happen; is very inefficient!
					(*this)(row+1,sm.colind[ksm+rindsm]+1+offc) += x*sm.mat[ksm+rindsm]; 
					k++; ksm++;
				}
			}

		}
	}
	//(*global_uo) << "AddSubmatrix rescue count=" << failed << "\n";
}


//set column vector, only entries of v from vr1 to vr2, insert at column col
void SparseMatrix::SetColVector(const Vector& v, int vr1, int vr2, int col)
{
	int row;
	for (int i = vr1; i <= vr2; i++)
	{
		if (fabs(v(i)) > sparsetolzero)
		{
			row = i-1;
			if ((rowlen[row] == 0 || col-1 > colind[rowind[row]+rowlen[row]-1]) && rowlen[row] < RowLenMax0(row))
			{
				mat[rowind[row]+rowlen[row]] = v(i);
				colind[rowind[row]+rowlen[row]] = col-1;
				rowlen[row]++;
			}
			else
			{
				(*this)(i,col) = v(i);
			}
		}
	}
}

void SparseMatrix::SetColVector(const SparseVector& v, int col)
{
	int row;
	double val;
	for (int i = 1; i <= v.NEntries(); i++)
	{
		v.GetEntry(i, row, val);
		if (fabs(val) > sparsetolzero)
		{
			row--;
			if ((rowlen[row] == 0 || col-1 > colind[rowind[row]+rowlen[row]-1]) && rowlen[row] < RowLenMax0(row))
			{
				mat[rowind[row]+rowlen[row]] = val;
				colind[rowind[row]+rowlen[row]] = col-1;
				rowlen[row]++;
			}
			else
			{
				(*this)(row+1,col) = val;
			}
		}
	}
}

//eliminate zero entries in sparse matrix (caused by zero entries in stiffness matrices or by adding of pieces of matrices!)
//otherwise, SuperLU does not work!!!
int SparseMatrix::EliminateZeroEntries(double mindouble)
{
	int n = 0; //global count
	
	if (rows != 0) release_assert(rowind[0] == 0);
	
	for (int i = 0; i < rows; i++)
	{
		int rz = 0; //zeros in row
		int oldrowind = rowind[i];
		if (i > 0) rowind[i] = rowind[i-1]+rowlen[i-1];  //there might be emtpy entries in the SparseMatrix structure!

		for (int j = 0; j < rowlen[i]; j++)
		{
			mat[rowind[i]+j-rz] = mat[oldrowind+j];
			colind[rowind[i]+j-rz] = colind[oldrowind+j];

			//if (fabs(mat[rowind[i]+j]) <= mindouble) 
			if (mat[oldrowind+j] == 0.) 
			{
				rz++;
			}
		}
		rowlen[i] -= rz; //reduce actual length
		n += rz;
	}
	return n;
}

void SparseMatrix::PrintData() const
{
	for (int i=0; i < rows; i++)
	{
		(*global_uo) << "row " << i << ", rowptr=" << rowind[i] << ", rowlen=" << rowlen[i] << ": ";
		for (int j = 0; j < rowlen[i]; j++)
		{
			(*global_uo) << "(" << colind[rowind[i]+j] << "," << mat[rowind[i]+j] << "), ";
		}
		(*global_uo) << "\n";
	}
}


ostream& operator<<(ostream& os, const SparseMatrix& m)
{
	double max=0;
	int i;
	max = m.MaxNorm();
	if (max==0) {max=1;}
	max=(int)log10(max);
	max=pow(10,max);

	os << max << " *\n";

	for (i=1; i<=m.Getrows(); i++)
	{
		os << "[";
		for (int j=1; j<=m.Getcols(); j++)
		{
			char str[32];
			sprintf(str,"% 1.4f",m(i,j)/max);
			os << str;
			if (j!=m.Getcols()) {os << ",";}
		}
		os << "]";
		os << "max=" << m.RowLenMax0(i) << "\n";
	}

	return os;
}


void SparseMatrix::PrintToMatlabSparse(ostream& os)
{
	int mode = 2;


	if (mode == 1)
	{
		os << "sparse([";
		//write rowindex in matrix (deliminated by space and ";")
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < rowlen[i]; j++)
			{
				char str[32];
				sprintf(str,"%d",i+1);
				os << str;
				if (!(i == rows-1 && j == rowlen[i]-1)) os << ";";
			}
		}
		os << "],[";
		//write colindex in matrix (deliminated by space and ";")
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < rowlen[i]; j++)
			{
				char str[32];
				sprintf(str,"%d",colind[rowind[i]+j]+1);
				os << str;
				if (!(i == rows-1 && j == rowlen[i]-1)) os << ";";
			}
		}
		os << "],[";
		//write matrix-entry in matrix (deliminated by space and ";")
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < rowlen[i]; j++)
			{
				char str[32];
				sprintf(str,"%1.16g",mat[rowind[i]+j]);
				os << str;
				if (!(i == rows-1 && j == rowlen[i]-1)) os << ";";
			}
		}
		os << "]";
		os << "," << Getrows() << "," << Getcols();
		os << ")";
	}
	else
	{
		//write rowindex, colindex, mat (deliminated by space)
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < rowlen[i]; j++)
			{
				char str[32];
				sprintf(str,"%d ",i+1);
				os << str;
				sprintf(str,"%d ",colind[rowind[i]+j]+1);
				os << str;
				sprintf(str,"%1.16g\n",mat[rowind[i]+j]);
				os << str;
			}
		}
		os << Getrows() << " " << Getcols() << " 0\n"; //highest entry specifies size of matrix ... 0 is added to remaining entries
	}
}





// AP -- SuperLU

SuperLUInverse:: ~SuperLUInverse()
{
	DestroySuperLUFactorization(SuperL, SuperU);
	delete[] perm_r;
	delete[] perm_c;
}

int SuperLUInverse:: Solve(Vector& q)
{
	rows = sparseA->Getrows();
	cols = sparseA->Getcols();

	int* colind = NULL;
	double* mat = NULL;
	int* rowptr = NULL;

	const_cast<SparseMatrix*>(sparseA)->EliminateZeroEntries();
	sparseA->GetMatPtrs(colind, rowptr, mat);

	int nnonzero = sparseA->CountEntries();
	rowptr[rows] = nnonzero;

	int output = 0;
	if (rows > 50000) output = 1;
	if (output) (*global_uo) << "SuperLU Solve ...";
	CallSuperLUSparseRow(rows, cols, nnonzero, rowptr, colind, mat, q.GetVecPtr());
	if (output) (*global_uo) << "done!\n";

	return 1;
}

int SuperLUInverse:: Factorize()
{

	int* colind = NULL;
	double* mat = NULL;
	int* rowptr = NULL;

	const_cast<SparseMatrix*>(sparseA)->EliminateZeroEntries();
	const_cast<SparseMatrix*>(sparseA)->EliminateZeroEntries();
	sparseA->GetMatPtrs(colind, rowptr, mat);
	int nnonzero = sparseA->CountEntries();
	rowptr[rows] = nnonzero;
	rows = sparseA->Getrows();
	cols = sparseA->Getcols();

	delete[] perm_r;
	delete[] perm_c;
	perm_r = new int(rows);
	perm_c = new int(cols);

	//int* perm_r_ptr = perm_r.GetDataPtr();
	//int* perm_c_ptr = perm_c.GetDataPtr();
	int rv = CallSuperLUFactorize(rows, cols, nnonzero, rowptr, colind, mat, 
											 SuperA, SuperL, SuperU, perm_r, perm_c, isFirstStep);

	isFactorized = 1;
	if (rv == 0) return 1;
	else return 0;
}

int SuperLUInverse:: SuperLUApply(Vector& q)
{
	if (!isFactorized) Factorize();

	rows = sparseA->Getrows();
	cols = sparseA->Getcols();

	//int* perm_r_ptr = perm_r.GetDataPtr();
	//int* perm_c_ptr = perm_c.GetDataPtr();
	int rv = CallSuperLUApply(rows, cols, q.GetVecPtr(1), SuperA, SuperL, SuperU, perm_r, perm_c);

	//(*global_uo) << "SuperLUApply returns " << rv << "\n";

	if (rv == 0) return 1;
	else return 0;
}



#ifdef USE_MKL
#include <mkl_pardiso.h> 

PardisoInverse::~PardisoInverse()
{
	//int maxfct=1, mnum=1;
	//int phase=-1;
	//int* ja=0, *ia=0;
	double* a=0;
	//int nrhs=1;
	//int msglevel=1;
	double*x=0;
	double*b=0;
	_INTEGER_t maxfct=1, mnum=1, phase=-1;
	_INTEGER_t* ja=0, *ia=0;
	_INTEGER_t nrhs=1, msglevel=1;


	pardiso (pt, &maxfct, &mnum, &mtype, &phase, &size, a, ia, ja, perm, &nrhs, iparm, &msglevel, b, x, &error);

	delete[] perm;
}

void PardisoInverse::SetDefaultIParm()
{
	iparm[0] = 0;  // use default values
	iparm[1] = 0;  // 0..minimum degree ordering
	iparm[2] = 1;  // number of processors used
	iparm[3] = 0; // direct..0, or PCG..1/2
	iparm[4] = 0; // 0..dont use user-permutation, 1..user-defined perm, 2..perm returned
	iparm[5] = 0; // 1 .. store solution in b, 0 .. store solution in x
	iparm[7] = 1; // iterative refinement steps
	iparm[8] = 0; // not used, has to be zero
	iparm[9] = 8; // accuracy
	if (mtype == 11 || mtype == 13) // unsymmetric matrices
	{
		iparm[10] = 1; // scaling
		iparm[12] = 1; // improved accuracy
	}
	else //  mtype = -2, -4, 6, symmetric matrices
	{
		iparm[10] = 0; // no scaling
		iparm[12] = 0; // no improved accuracy
	}
	iparm[11] = 0; // not used, has to be zero
	iparm[17] = -1; // report number of non-zero entries if < 0
	iparm[18] = -1; // report number of MFlops
	iparm[20] = 0; // 1x1 pivoting
	iparm[26] = 0; // check arrays ia, ja
	iparm[27] = 0; // 0..double precision, 1..single
	iparm[59] = 0; // in-core..0, out of core=1
}


int PardisoInverse::Solve(Vector& q)
{
	// initialize internal data pointer (necessary for MKL-pardiso)
	for (int i=0; i<64; i++)
		pt[i] = 0;

	//int maxfct=1, mnum=1;
	//int phase=13;
	//int* ja, *ia;
	double* a;
	//int nrhs=1;
	//int msglevel=1;
	_INTEGER_t maxfct=1, mnum=1, phase=13;
	int* ja,* ia;
	_INTEGER_t nrhs=1, msglevel=1;


	Vector sol(q.Length());

	iparm[5] = 1; // replace right hand side by solution
	double*b=q.GetVecPtr();
	double*x=sol.GetVecPtr();

	const_cast<SparseMatrix*>(sparseA)->EliminateZeroEntries();
	sparseA->GetMatPtrs(ja, ia, a);
	int nnz = sparseA->CountEntries();
	_INTEGER_t* ia2 = new _INTEGER_t[size+1];
	for (int i=0; i<size; i++)
		ia2[i] = ia[i]+1;
	ia2[size] = nnz+1;
	_INTEGER_t* ja2 = new _INTEGER_t[nnz];
	for (int i=0; i<nnz; i++)
		ja2[i] = ja[i]+1;


	PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &size, a, ia2, ja2, perm, &nrhs, iparm, &msglevel, b, x, &error);

	isFactorized = 1;
	isFirstStep = 0;
	//q = sol;

	delete [] ia2;
	delete [] ja2;
	if (!error) return 1;
	else return 0;
}

int PardisoInverse::Factorize()
{
		// initialize internal data pointer (necessary for MKL-pardiso)
	for (int i=0; i<64; i++)
		pt[i] = 0;

	//int maxfct=1, mnum=1;
	//int phase=13;
	//int* ja, *ia;
	double* a;
	//int nrhs=1;
	//int msglevel=1;
	_INTEGER_t maxfct=1, mnum=1, phase=12;
	int* ja,* ia;
	_INTEGER_t nrhs=1, msglevel=1;

	// initialize internal data pointer (necessary for MKL-pardiso)
	for (int i=0; i<64; i++)
		pt[i] = 0;

	double*b=0;
	double*x=0;

	const_cast<SparseMatrix*>(sparseA)->EliminateZeroEntries();
	sparseA->GetMatPtrs(ja, ia, a);
	int nnz = sparseA->CountEntries();
	_INTEGER_t* ia2 = new _INTEGER_t[size+1];
	for (int i=0; i<size; i++)
		ia2[i] = ia[i]+1;
	ia2[size] = nnz+1;
	_INTEGER_t* ja2 = new _INTEGER_t[nnz];
	for (int i=0; i<nnz; i++)
		ja2[i] = ja[i]+1;

	PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &size, a, ia2, ja2, perm, &nrhs, iparm, &msglevel, b, x, &error);
	//PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &size, a, ia2, ja2, perm, &nrhs, iparm, &msglevel, b, x, &error);

	isFactorized = 1;
	isFirstStep = 0;
	delete [] ia2;
	delete [] ja2;
	if (!error) return 1;
	else return 0;
}
int PardisoInverse::Apply(Vector& q)
{
	//int maxfct=1, mnum=1;
	//int phase=33;
	//int *ja, *ia;
	double* a;
	//int nrhs=1;
	//int msglevel=1;
	_INTEGER_t maxfct=1, mnum=1, phase=33;
	int* ja,* ia;
	_INTEGER_t nrhs=1, msglevel=1;

	iparm[5] = 1; // replace right hand side by solution

	double*b = q.GetVecPtr();
	double*x = new double[q.Length()];
	sparseA->GetMatPtrs(ja, ia, a);
	int nnz = sparseA->CountEntries();
	_INTEGER_t* ia2 = new _INTEGER_t[size+1];
	for (int i=0; i<size; i++)
		ia2[i] = ia[i]+1;
	ia2[size] = nnz+1;
	_INTEGER_t* ja2 = new _INTEGER_t[nnz];
	for (int i=0; i<nnz; i++)
		ja2[i] = ja[i]+1;

	PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &size, a, ia2, ja2, perm, &nrhs, iparm, &msglevel, b, x, &error);
	//PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &size, a, ia2, ja2, perm, &nrhs, iparm, &msglevel, b, x, &error);

	delete [] ia2;
	delete [] ja2;
	delete [] x;

	if (!error) return 1;
	else return 0;
}

#else

void PardisoInverse::SetDefaultIParm()
{
}

	PardisoInverse::~PardisoInverse()
{
	delete[] perm;
}

int PardisoInverse::Solve(Vector& q)
{
	(*global_uo) << "Error in PardisoInverse::Solve: Pardiso Solver not available!\n";
	return 0;
}

int PardisoInverse::Factorize()
{
	(*global_uo) << "Error in PardisoInverse::Factorize: Pardiso Solver not available!\n";
	return 0;
}
int PardisoInverse::Apply(Vector& q)
{
	(*global_uo) << "Error in PardisoInverse::Apply: Pardiso Solver not available!\n";
	return 0;
}

#endif //USE_MKL for Pardiso


SparseInverse::SparseInverse(const SparseMatrix & matA) 
{
#ifdef USE_MKL
	pardisoinv = new PardisoInverse(matA);
	superluinv = NULL;
#else
	superluinv = new SuperLUInverse(matA);
	pardisoinv = NULL;
#endif
}


///--------------------
/// SparseInverse wrapper - choose pardiso or superlu
SparseInverse::~SparseInverse()
{
	delete pardisoinv;
	delete superluinv;
}

void SparseInverse::SetSparseMatrix(const SparseMatrix& matA)
{
#ifdef USE_MKL
	pardisoinv->SetSparseMatrix(matA);
#else
	superluinv->SetSparseMatrix(matA);
#endif
}

int SparseInverse::Solve(Vector& q)
{
#ifdef USE_MKL
	return pardisoinv->Solve(q);
#else
	return superluinv->Solve(q);
#endif

}
int SparseInverse::Factorize()
{
#ifdef USE_MKL
	return pardisoinv->Factorize();
#else
	return superluinv->Factorize();
#endif
}
int SparseInverse::Apply(Vector& q)
{
#ifdef USE_MKL
	return pardisoinv->Apply(q);
#else
	return superluinv->Apply(q);
#endif

}

int& SparseInverse::IsFirstStep() 
{
#ifdef USE_MKL
	return pardisoinv->IsFirstStep();
#else
	return superluinv->IsFirstStep();
#endif
}
const int& SparseInverse::IsFirstStep() const 
{
#ifdef USE_MKL
	return pardisoinv->IsFirstStep();
#else
	return superluinv->IsFirstStep();
#endif
}


int SuperLUSolve(SparseMatrix& A, Vector& q)
{
	int* colind = NULL;
	double* mat = NULL;
	int* rowptr = NULL;


/*
	Matrix AF = A.GetMatrix();
	(*global_uo) << "AF=" << AF << "\n";
	A.PrintData();
	*/

	A.EliminateZeroEntries();
	A.GetMatPtrs(colind, rowptr, mat);

	//A.CopyFrom(AF);
	//A.PrintData();
	//A.PrintToMatlabSparse(fout);

	int nnonzero = A.CountEntries();
	int rows = A.Getrows();
	int cols = A.Getcols();
	rowptr[rows] = nnonzero;

	int output = 0;
	if (rows > 50000) output = 1;
	if (output) (*global_uo) << "SuperLU Solve ...";
	CallSuperLUSparseRow(rows, cols, nnonzero, rowptr, colind, mat, q.GetVecPtr());
	if (output) (*global_uo) << "done!\n";

	return 1;
}

int SuperLUFactorize(SparseMatrix& A, SuperMatrix*& SuperA, SuperMatrix*& L, SuperMatrix*& U, int*& perm_r, int*& perm_c, int isFirstStep)
{
	int* colind = NULL;
	double* mat = NULL;
	int* rowptr = NULL;

	A.EliminateZeroEntries();
	A.GetMatPtrs(colind, rowptr, mat);
	int nnonzero = A.CountEntries();
	int rows = A.Getrows();
	int cols = A.Getcols();
	rowptr[rows] = nnonzero;

	int rv = CallSuperLUFactorize(rows, cols, nnonzero, rowptr, colind, mat, 
											 SuperA, L, U, perm_r, perm_c, isFirstStep);

	// (*global_uo) << "SuperLUFactorize returns " << rv << "\n";

	if (rv == 0) return 1;
	else return 0;
}

int SuperLUApply(SparseMatrix& A, SuperMatrix*& SuperA, SuperMatrix*& L, SuperMatrix*& U, int*& perm_r, int*& perm_c, Vector& q)
{
	int rows = A.Getrows();
	int cols = A.Getcols();

	int rv = CallSuperLUApply(rows, cols, q.GetVecPtr(), SuperA, L, U, perm_r, perm_c);

	// (*global_uo) << "SuperLUApply returns " << rv << "\n";

	if (rv == 0) return 1;
	else return 0;
	//return 1;
}



void GetExtremeValues(const Vector& data, TArray<int>& max_indices,TArray<double>& max_values, TArray<int>& min_indices,TArray<double>& min_values, int compute_minima)
{
	max_indices.SetLen(0);
	max_values.SetLen(0);
	
	if(compute_minima)
	{
		min_indices.SetLen(0);
		min_values.SetLen(0);
	}

	if(data.GetLen()==1)
	{
		return;
	}
	int sign = 1; //search maxima
	double dataoldold = data.Get(1)*sign;
	double dataold    = data.Get(2)*sign;
  
	// check if first value is maxima or minima
	if(dataoldold > dataold)
	{
		max_indices.Add(1);
		max_values.Add(dataoldold);
	}
	if(compute_minima && dataoldold < dataold)
	{
		min_indices.Add(1);
		min_values.Add(dataoldold);
	}
	
	// check local minima inside
	for(int j=1;j<=2;j++)
	{
		if(j==2)sign=-1; //search minima
		
		for(int i=3;i<=data.Length();i++)
		{
			if(dataoldold < dataold && dataold > data.Get(i)*sign)
			{
				// data(i-1)*factor is already a maximum 
				if(sign>0)
				{
					max_indices.Add(i-1);
					max_values.Add(dataold);
				}
				else if(compute_minima)
				{
					min_indices.Add(i-1);
					min_values.Add(dataold);
				}				
			}
			dataoldold = dataold;
			dataold = data.Get(i)*sign;
		}
	}

	// check if first value is maxima or minima
	if(data.Get(data.Length()) > data.Get(data.Length()-1))
	{
		max_indices.Add(data.Length());
		max_values.Add(data.Get(data.Length()));
	}
	if(compute_minima && data.Get(data.Length()) < data.Get(data.Length()-1))
	{
		min_indices.Add(data.Length());
		min_values.Add(data.Get(data.Length()));
	}	
}

//symmetrical smoothing
void SmoothenDoubleArray(Vector& values, Vector& result, int smoothing_length)
{
	if (smoothing_length == 0)
	{
		result = values;
		return;
	}
	for(int i=1; i<=values.Length(); i++)
	{
		double avg = 0;
		int cnt = 0;
		for (int j = -smoothing_length; j<=smoothing_length; j++)
		{
			if (i+j >= 1 && i+j<=values.Length())
			{
				avg += values(i+j);
				cnt++;
			}
		}
		if (cnt)
			result(i) = avg/(double)cnt;
		else
			result(i) = 0;
	}
}

/*
//does not work: old version:

//generate new matrix from sparse matrix, only from row1 to row2, col1 to col2 (incl.), resort with vector r
void SparseMatrix::CopyFrom(const SparseMatrix& mati, int row1, int col1, int row2, int col2, const IVector& r)
{
int nrows = row2-row1+1;
int ncols = col2-col1+1;

mystatic IVector invr;
invr.SetLen(r.Length());
for (int i=1; i <= r.Length(); i++) invr(r(i)) = i;

int neededmem = 0;


for (int i = 0; i < nrows; i++)
{
int matirow = r(i+row1)-1;
int rind = mati.rowind[matirow];

//read out necessary entries and transform by r:
for (int j = 0; j < mati.rowlen[matirow]; j++)
{
int c = mati.colind[rind+j]; //0-based
int cr = r(c+1); //1-based
if (cr >= col1 && cr <= col2) neededmem++;
}
}


if (neededmem > lalloc || nrows > rows)
{
Destroy();

lalloc = neededmem;
GenMat(lalloc, nrows);
}

rows = nrows;
cols = ncols;



mystatic TArray<double> rsort; //resorted values
mystatic TArray<int> rsortind; //resorted column indices 0-based for destination


//sm    1  2  3  4  5  6  7  8  9
//rind  1     3        6     8
//mat   11    12       13    14
//r     1  2  8  3  4  5  9  6  7
//invr  1  2  4  5  6  8  9  3  7
//m     11    14 12          13    
//

int c, cr;
int cnt = 0;

//(*global_uo) << "lalloc=" << lalloc << ", nrows=" << nrows << "\n";

//TMStartTimer(24);
for (int i = 0; i < rows; i++)
{
rsortind.SetLen(0);
rsort.SetLen(0); 

int matirow = r(i+row1)-1;
int rind = mati.rowind[matirow];

//read out necessary entries and transform by r:
for (int j = 0; j < mati.rowlen[matirow]; j++)
{
c = mati.colind[rind+j]; //0-based
cr = r(c+1); //1-based
if (cr >= col1 && cr <= col2)
{
rsortind.Add(invr(c+1)-col1); //0-based index for destination
//runsortind.Add(c-col1+1); //0-based index for destination
rsort.Add(mati.mat[rind+j]); //entry
}
}
//sort entries
//TMStartTimer(25);
Quicksort(rsortind, rsort);
//TMStopTimer(25);

//fill in entries in i-th row
rowind[i] = cnt;
int cnt2 = 0;



for (j=0; j < rsortind.Length(); j++)
{
if (rsort.Get0(j) != 0)
{
colind[cnt2+rowind[i]] = rsortind.Get0(j);
//(*global_uo) << "(" << cnt2+rowind[i] << "," << colind[cnt2+rowind[i]] << ") ";

mat[cnt2+rowind[i]] = rsort.Get0(j);
cnt2++;
}
}
//(*global_uo) << "\n";
*/