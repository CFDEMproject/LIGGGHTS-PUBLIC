//#**************************************************************
//#
//# filename:             linalg.h
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
#pragma once


#define gen_sparsemat_off

int GenVecCnt();
int GenMatCnt();

void GenSparsemat();

int GetGenSparsemat();


//-----------------------------------------------------------------
//--------------    CLASS VECTOR    -------------------------------
//-----------------------------------------------------------------
//Generation of a vector with user-defined length
//Operations: +,-,* (also matrices), <<
class Vector
{
public:

	//Constructors, Destructor
	Vector(): l(0),vec(NULL),lalloc(0),ownmemory(0) {};
	Vector(int leni, int nofill=0, double fill=0);
	Vector(const Vector& veci);
	Vector(const Vector& veci, int noownmem);
	Vector(double x0, double x1, double x2, double x3, double x4, double x5, double x6)
	{
		GenVector(7); vec[0]=x0; vec[1]=x1; vec[2]=x2; vec[3]=x3; vec[4]=x4; vec[5]=x5; vec[6]=x6; 
	}

	Vector(double x0, double x1, double x2, double x3, double x4, double x5)
	{
		GenVector(6); vec[0]=x0; vec[1]=x1; vec[2]=x2; vec[3]=x3; vec[4]=x4; vec[5]=x5; 
	}
	Vector(double x0, double x1, double x2, double x3, double x4)
	{
		GenVector(5); vec[0]=x0; vec[1]=x1; vec[2]=x2; vec[3]=x3; vec[4]=x4;
	}
	Vector(double w, double x, double y, double z)
	{
		GenVector(4);vec[0]=w; vec[1]=x; vec[2]=y; vec[3]=z;
	}
	Vector(double x, double y, double z)
	{
		GenVector(3);vec[0]=x; vec[1]=y; vec[2]=z;
	}
	Vector(double x, double y)
	{
		GenVector(2);vec[0]=x; vec[1]=y;
	}
	Vector(double x)
	{
		GenVector(1);vec[0]=x;
	}
	Vector(const TArray<double>& ta){SetVector(ta);} //(RL)
	
	
	void SetVector(const TArray<double>& ta) //(RL)
	{
		GenVector(ta.Length());
		for (int i=1; i <= ta.Length(); i++)
		{
			vec[i-1]=ta(i);
		}
	}

	Vector(const char* vector, char osbc = '[', char csbc = ']', char commac = ',');

	Vector(Vector3D& v1, Vector3D& v2)
	{
		GenVector(6); vec[0]=v1.X(); vec[1]=v1.Y(); vec[2]=v1.Z(); vec[3]=v2.X(); vec[4]=v2.Y(); vec[5]=v2.Z(); 
	}
	virtual ~Vector() 
	{
		if (vec!=NULL && ownmemory && lalloc) 
		{delete [] vec; vec=NULL;} 
	};

	virtual void Init()
	{
		l = 0;
		vec = NULL;
		lalloc = 0;
		ownmemory = 0;
	};
	//Generates a vector with length "leni" (leni=0 possible)
	virtual void GenVector(int leni);
	//Link with other vector, no own memory!!!
	virtual void LinkWith(const Vector& v)
	{
		if (vec!=NULL && ownmemory && lalloc) 
		{delete [] vec; vec=NULL;} 
		ownmemory = 0;
		vec = v.vec;
		l = v.l;
		lalloc = 0;
	}
	virtual void LinkWith(double* ptr, int len)
	{
		if (vec!=NULL && ownmemory && lalloc) 
		{delete [] vec; vec=NULL;} 
		ownmemory = 0;
		vec = ptr;
		l = len;
		lalloc = 0;
	}
	//Link with other vector, no own memory!!!
	virtual void LinkWith(const Vector& v, int pos, int newlen)
	{
		if (vec!=NULL && ownmemory && lalloc) {delete [] vec; vec=NULL;} 
		ownmemory = 0;
		vec = v.vec+(pos-1); //no sizeof(double) !!!!!!!!!!!!
		l = newlen;
		lalloc = 0;
		if (newlen == 0) vec = NULL;
	}

	virtual int IsValid(double maxval);

	virtual void SetAll(double x)
	{
		for(int i=0; i< l; i++) vec[i]=x;
	}

	virtual void SetAll(double x, int i1, int i2)
	{
		for(int i = i1-1; i <= i2-1; i++) vec[i]=x;
	}

	virtual void SetRandom()
	{
		for (int i=0; i<l; i++) vec[i] = (double)rand()/(double)RAND_MAX;
	}

	virtual int IsZero(double eps) const //check if it is smaller than tolerance eps
	{
		int found = 0;
		int i = 0;
		if (eps == 0)
		{
			while (!found && i < l)
			{
				if (vec[i++] != 0) found = 1;
			}
		}
		else
		{
			while (!found && i < l)
			{
				if (fabs(vec[i++]) > eps) found = 1;
			}
		}
		return !found;
	}

	virtual void Mult(double x)
	{
		for(int i=0; i< l; i++) vec[i]*=x;
	}

	virtual void SetLen(int leni);

	//Generate vector with 3 components
	virtual void Set3D(double x,double y,double z)
	{
		if (l!=3)
		{
			if (vec!=NULL && ownmemory && lalloc) {delete[] vec;}
			GenVector(3);
		}
		vec[0]=x; vec[1]=y; vec[2]=z;
	};

	//Generate vector with 2 components
	virtual void Set2D(double x,double y)
	{
		if (l!=2)
		{
			if (vec!=NULL && ownmemory && lalloc) {delete[] vec;}
			GenVector(2);
		}
		vec[0]=x; vec[1]=y;
	};

	//Generate vector with 6 components
	virtual void Set6D(double v1, double v2, double v3, double v4, double v5, double v6)
	{
		if (l!=6)
		{
			if (vec!=NULL && ownmemory && lalloc) {delete[] vec;}
			GenVector(6);
		}
		vec[0]=v1;
		vec[1]=v2;
		vec[2]=v3;
		vec[3]=v4;
		vec[4]=v5;
		vec[5]=v6;
	};

	//Assignment-operator
	virtual Vector& operator= (const Vector& veci);

	//Special Assignment-operator
	virtual Vector& operator= (const Vector3D& veci);

	//Logical operator
	virtual int operator== (const Vector& vec1);



	//Referencing access-operator, ZERO-based
	virtual double& operator[](int elem)
	{
#ifndef __QUICKMATH
		release_assert((elem>=0) && (elem<l));
#endif
		return vec[elem];
	};


	//Referencing ccess-operator for constant access, ZERO-based
	virtual const double& operator[](int elem) const
	{
#ifndef __QUICKMATH
		release_assert((elem>=0) && (elem<l));
#endif
		return vec[elem];
	};


	//Referencing access-operator
	virtual double& operator()(int elem)
	{
#ifndef __QUICKMATH
		release_assert((elem>0) && (elem<=l));
#endif
		return vec[elem-1];
	};


	//Referencing access-operator for constant access
	virtual const double& operator()(int elem) const
	{
#ifndef __QUICKMATH
		release_assert((elem>0) && (elem<=l));
#endif
		return vec[elem-1];
	};


	//Returns the elem-component of a vector
	virtual double& Elem(int elem)
	{
#ifndef __QUICKMATH
		release_assert((elem>0) && (elem<=l));
#endif
		return vec[elem-1];
	};

	//Returns the constant elem-component of a vector
	virtual const double& Get(int elem) const
	{
#ifndef __QUICKMATH
		release_assert((elem>0) && (elem<=l));
#endif
		return vec[elem-1];
	};

	//Returns the elem-component of a vector, 0-based
	virtual double& Elem0(int elem)
	{
#ifndef __QUICKMATH
		release_assert((elem>=0) && (elem<l));
#endif
		return vec[elem];
	};

	//Returns the constant elem-component of a vector, 0-based
	virtual const double& Get0(int elem) const
	{
#ifndef __QUICKMATH
		release_assert((elem>=0) && (elem<l));
#endif
		return vec[elem];
	};

	virtual double* GetVecPtr() const {return vec;}
	virtual double* GetVecPtr(int i) const {return &(vec[i-1]);}

	//Returns the length of a vector
	virtual int GetLen() const {return l; };
	virtual int Length() const {return l; };

	//Returns the quadratic norm of a vector
	virtual double GetNorm() const;

	//Returns the maximum-norm of a vector
	virtual double MaxNorm() const;

	//Returns the minimum-norm of a vector
	virtual double MinNorm() const;

	//Returns the normalized planar perpendicular vector of a vector
	virtual Vector GetNormed2DNormalVec() const;

	//Returns the normalized planar vector of a vector
	virtual Vector GetNormed2DVec() const;

	//Returns the cross sum of the components of a vector
	virtual double Sum() const;

	//Returns the cross sum of the first n components of a vector
	virtual double Sum(int n) const;

	//Fills a vector with zeros
	virtual void FillWithZeros();

	//Inserts a vector v at reference point ref
	//    void FillInVector(const Vector& v, const IVector& ref);

	//Inserts a vector v at reference point ref and returns it
	virtual Vector Append(const Vector& v) const; 

	virtual void Insert(const Vector& v, int pos); 
	virtual void Copy(const Vector& v, int vpos, int thispos, int len)
	{
		int i;
		for (i = 0; i < len; i++)
		{
			Elem(i+thispos) = v(i+vpos);
		}
	}

	//Returns a fraction of a vector
	virtual Vector SubVector(int from, int to) const;

	//add k*v to vector
	virtual void MultAdd(const double& k, const Vector& v); 

	//Arithmetic operations with one parameter
	virtual Vector& operator+= (const Vector& v1);
	virtual Vector& operator-= (const Vector& v1);
	virtual Vector& operator*= (const double& v);

	//Arithmetic operations with two parameters
	friend Vector operator+ (const Vector& v1, const Vector& v2);
	friend Vector operator- (const Vector& v1, const Vector& v2);
	friend Vector operator* (const Vector& vec, const double& val);
	friend Vector operator* (const double& val, const Vector& vec);
	friend Vector operator* (const Vector& v, const Matrix& m);
	friend Vector operator* (const Matrix& m, const Vector& v);
	friend double operator* (const Vector& vec1, const Vector& vec2);
	friend Vector3D operator* (const Matrix3D& m, const Vector& v);
	friend Vector2D operator* (const Matrix2D& m, const Vector& v);
	friend void Mult(const Matrix2D& m, const Vector2D& v, Vector& res);

	friend void Mult(const Matrix& m, const Vector& v, Vector& res); //computes res=m*v
	friend void MultTp(const Matrix& m, const Vector& v, Vector& res); //computes res=m.GetTp()*v
	friend void Mult(const Matrix& m, const Vector& v, Vector3D& res); //computes res=m*v
	friend void MultTp(const Matrix& m, const Vector& v, Vector3D& res); //computes res=m*v
	friend void Mult(const Vector& v, const Matrix& m, Vector& res); //computes res=m*v
	friend void Mult(const Vector& v1, const Vector& v2, Matrix& res); //computes res=v1*v2^T, where v1 and v2 are interpreted as a column vectors
	friend void MultBW(const Matrix& m, const Vector& v, Vector& res, int bw); //computes res=m*v with bandwidth bw
	friend void Mult(const Matrix3D& m, const Vector3D& v, Vector& res);
	friend void Mult(const MatrixXD& m, const Vector& v, Vector& res);
	friend void Mult(const Matrix& m, const Vector3D& v, Vector& res); //computes res=m*v
	friend void Mult(const Matrix& m, const Vector2D& v, Vector& res); //computes res=m*v

	friend void Mult(const SparseMatrix& m, const Vector& v, Vector& res); //computes res=m*v

	//transform strain and stress vectors to matrices and vice versa (by PG)
	friend void StrainVectorToMatrix2D(Matrix2D& m, const Vector& v);
	friend void StrainVectorToMatrix2D(Matrix2D& m, const Vector3D& v);
	friend void StressVectorToMatrix2D(Matrix2D& m, const Vector& v);
	friend void Matrix2DToStrainVector(Vector& v, const Matrix2D& m);
	friend void Matrix2DToStressVector(Vector& v, const Matrix2D& m);

	friend void StrainVectorToMatrix3D(Matrix3D& m, const Vector& v);
	friend void StressVectorToMatrix3D(Matrix3D& m, const Vector& v);
	friend void Matrix3DToStrainVector(Vector& v, const Matrix3D& m);
	friend void Matrix3DToStressVector(Vector& v, const Matrix3D& m);


	// Vector entries are equally distributed values of piecewise linear function on [-1, 1]
	// Interpolate returns value of function in point ploc € [-1, 1]
	// Vector entries specify function values, ordered as
	//   [ 1          2         3   ... width ]       -->      [x=-1         x=1]
	virtual double Interpolate(const double ploc) const;

	//Output parameter
	mystr& MakeString(mystr intro = mystr(""));
	friend ostream& operator<<(ostream& os, const Vector& v);

protected:

	double* vec;
	int l;
	int lalloc;
	int ownmemory;
};

template <int data_size>
class ConstVector: public Vector
{
public:

	//Constructors, Destructor
	ConstVector()
	{
		lalloc = 0;
		ownmemory = 0;
		l = data_size;
		vec = &constdata[0];
		SetAll(0.);
	}
	ConstVector(const ConstVector& veci)
	{
		lalloc = 0;
		ownmemory = 0;
		vec = &constdata[0];

#ifndef __QUICKMATH
		release_assert(data_size >= veci.l);
#endif
		l = veci.l;
		for (int i=0; i < l; i++)
		{
			constdata[i] = veci.vec[i];
		}
	}

	ConstVector(int leni, int nofill=0, double fill=0)
	{
		lalloc = 0;
		ownmemory = 0;
		l = leni;
		vec = &constdata[0];
		if (!nofill) SetAll(fill);
	}

	ConstVector(double x0)
	{
		release_assert(data_size >= 1);
		InitConstVector(); vec[0]=x0;
	}
	ConstVector(double x0, double x1)
	{
		release_assert(data_size >= 2);
		InitConstVector(); vec[0]=x0; vec[1]=x1;
	}
	ConstVector(double x0, double x1, double x2)
	{
		release_assert(data_size >= 3);
		InitConstVector(); vec[0]=x0; vec[1]=x1; vec[2]=x2; 
	}
	ConstVector(double x0, double x1, double x2, double x3)
	{
		release_assert(data_size >= 4);
		InitConstVector(); vec[0]=x0; vec[1]=x1; vec[2]=x2; vec[3]=x3;
	}
	ConstVector(double x0, double x1, double x2, double x3, double x4)
	{
		release_assert(data_size >= 5);
		InitConstVector(); vec[0]=x0; vec[1]=x1; vec[2]=x2; vec[3]=x3; vec[4]=x4; 
	}
	ConstVector(double x0, double x1, double x2, double x3, double x4, double x5)
	{
		release_assert(data_size >= 6);
		InitConstVector(); vec[0]=x0; vec[1]=x1; vec[2]=x2; vec[3]=x3; vec[4]=x4; vec[5]=x5; 
	}
	ConstVector(double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7, double x8, double x9)
	{
		release_assert(data_size >= 10);
		InitConstVector(); vec[0]=x0; vec[1]=x1; vec[2]=x2; vec[3]=x3; vec[4]=x4; vec[5]=x5; vec[6]=x6; vec[7]=x7; vec[8]=x8; vec[9]=x9;
	}
	virtual void InitConstVector()
	{
		lalloc = 0;
		ownmemory = 0;
		l = data_size;
		vec = &constdata[0];
	}
	virtual int GetDataSize() const {return data_size;}


	virtual ConstVector<data_size> & operator= (const Vector& veci)
	{
		if (this == &veci) {return *this;}
		if (veci.Length()!=0)
		{
#ifndef __QUICKMATH
			release_assert(data_size>=veci.Length());
#endif
			l = veci.Length();
			memcpy(vec,&veci[0],l*sizeof(double));
		}
		else 
		{
			this->l=0;
			//lalloc = 0; vec = NULL;
		}
		return *this;
	}

private:
	double constdata[data_size];
};




//-----------------------------------------------------------------
//--------------    CLASS MATRIX    -------------------------------
//-----------------------------------------------------------------
//Generation of a Matrix with user-defined size
//Operations: +,-,*,<<
class Matrix
{
public:

	//Constructors,destructor
	Matrix()
	{
		rows=0;
		cols=0;
		mat=NULL;
		lalloc = 0;
	};

	Matrix(int rowsi, int colsi, int nofill=0, double fill=0);
	Matrix(double x);
	Matrix(double x11, double x12, double x21, double x22);
	Matrix(double x11, double x12, double x13,
		double x21, double x22, double x23,
		double x31, double x32, double x33);

	Matrix(const Matrix& mati);
	Matrix(const Matrix3D mat3d);

	virtual ~Matrix() 
	{
		if (mat!=NULL && lalloc != 0) {delete [] mat; mat=NULL;} 
	};


	Matrix(int r, int c, double * ptr)
	{
		lalloc = 0;
		mat = ptr;
		rows = r;
		cols = c;
	}

	virtual void Init()
	{
		rows=0;
		cols=0;
		mat=NULL;
		lalloc = 0;
	};
	virtual int IsConstSizeMatrix() const {return 0;}

	//Setsize of matrix, don't change if ok
	virtual void SetSize(int rowsi, int colsi)
	{
		if (rows!=rowsi || cols!=colsi || mat == NULL) {Resize(rowsi,colsi,1);}
	}
	virtual double* GetMatPtr() const {return mat;}
	virtual int MaxAllocatedSize() const {return lalloc;}

	//Set new size of matrix, delete old one
	virtual void Resize(int rowsi, int colsi, int nofill=0, double fill=0);

	virtual void CopyFrom(const Matrix& m, int row1, int col1, int row2, int col2);

	virtual void CopyFrom(const Matrix& m, int row1, int col1, int row2, int col2, const IVector& r);

	virtual void CopyFrom(const SparseMatrix& m, int row1, int col1, int row2, int col2);

	virtual void CopyFrom(const SparseMatrix& m, int row1, int col1, int row2, int col2, const IVector& r);

	//Access-functions for private data
	virtual int Getrows() const {return rows; };
	virtual int Getcols() const {return cols; };

	//Assignment-operator
	virtual Matrix& operator= (const Matrix& mati);

	//Assignment-operator
	virtual Matrix& operator= (const double& val)
	{
		for (int i=0; i < rows*cols; i++)
		{
			mat[i] = val;
		}
		return *this;
	}

	//Returns Determinate of matrix
	virtual double Det() const;

	//Returns the transposed matrix
	virtual Matrix GetTp() const;


	virtual double* operator[](int row0) //zero based for [][] access
	{
#ifndef __QUICKMATH
		release_assert((row0>=0) && (row0<rows));
#endif
		return &(mat[row0*cols]);
	};

	virtual const double* operator[](int row0) const //zero based for [][] access
	{
#ifndef __QUICKMATH
		release_assert((row0>=0) && (row0<rows));
#endif
		return &(mat[row0*cols]);
	};

	/*
	//Referencing access-operator on element using absolute indexing
	double& operator[](int elem)
	{
	#ifndef __QUICKMATH
	release_assert((elem>0) && (elem<=rows*cols));
	#endif
	return mat[elem-1];
	};


	//Referencing constant access-operator on element using absolute indexing
	const double& operator[](int elem) const
	{
	#ifndef __QUICKMATH
	release_assert((elem>0) && (elem<=rows*cols));
	#endif
	return mat[elem-1];
	};
	*/


	//Referencing access-operator on element using row- and column-values
	virtual double& operator()(int row, int col)
	{
#ifndef __QUICKMATH
		release_assert((row>0) && (row<=rows) && (col>0) && (col<=cols));
#endif
		return mat[(row-1)*cols+col-1];
	};


	//Referencing constant access-operator on element using row- and column-values
	virtual const double& operator()(int row, int col) const
	{
#ifndef __QUICKMATH
		release_assert((row>0) && (row<=rows) && (col>0) && (col<=cols));
#endif
		return mat[(row-1)*cols+col-1];
	};


	//Referencing access-operator on element using row- and column-values
	virtual double& Elem(int row, int col)
	{
#ifndef __QUICKMATH
		release_assert((row>0) && (row<=rows) && (col>0) && (col<=cols));
#endif
		return mat[(row-1)*cols+col-1];
	};


	//Referencing constant access-operator on element using row- and column-values
	virtual const double& Get(int row, int col) const
	{
#ifndef __QUICKMATH
		release_assert((row>0) && (row<=rows) && (col>0) && (col<=cols));
#endif
		return mat[(row-1)*cols+col-1];
	};

	//Referencing access-operator on element using row- and column-values, 0-based
	virtual const double& Get0(int row, int col) const
	{
#ifndef __QUICKMATH
		release_assert((row>=0) && (row<rows) && (col>=0) && (col<cols));
#endif
		return mat[row*cols+col];
	};

	//Referencing access-operator on element using row- and column-values, 0-based
	virtual double& Elem0(int row, int col)
	{
#ifndef __QUICKMATH
		release_assert((row>=0) && (row<rows) && (col>=0) && (col<cols));
#endif
		return mat[row*cols+col];
	};


	virtual void LinkWith(int r, int c, double * ptr)
	{
		if (mat!=NULL && lalloc != 0) 
		{delete [] mat; mat=NULL;} 
		lalloc = 0;
		mat = ptr;
		rows = r;
		cols = c;
	}

	//Access-operator for symmetric band matrix, returns value of element
	virtual double BGet(int row, int col) const
	{
		if (row > col) {Swap (row,col);}
		release_assert((row>0) && (row<=rows) && (col>0) && (col<=rows));
		if ((col>row+cols-1)) {return 0;}
		return mat[(row-1)*cols+col-row];
	};


	//Access-operator for symmetric band matrix, sets value of element
	virtual void Bset(int row, int col, const double& v)
	{
		if (row > col) {Swap (row,col);}
		release_assert((row>0) && (row<=rows) && (col<row+cols));
		if ((col>row+cols-1)) {return;}
		mat[(row-1)*cols+col-row]=v;
	};


	//Converts symmetric band-matrix to matrix
	virtual Matrix BToMatrix();

	//Converts symmetric matrix to symmetric band-matrix
	virtual Matrix MatrixToB();

	//Fills matrix with zeros
	virtual void FillWithZeros();

	virtual void SetAll(double x);

	//Fills matrix with zeros
	virtual void FillWithZeros(int row, int col, int nrows, int ncols);

	//Returns the maximum-norm (largest value in matrix)
	virtual double MaxNorm() const;


	//Returns the minimum-norm (smallest value in matrix)
	virtual double MinNorm() const;

	//Returns the quadratic-norm
	virtual double Norm2() const;

	//Returns true if matrix is symmetric
	virtual int IsSymmetric();
	virtual int IsSymmetric(double eps); //version with small tolerance epsilon

	virtual void MakeSymmetric(); //makes a matrix really symmetric (if there are small errors)

	//Returns true if matrix is quadratic
	virtual int IsSquare() const {return rows==cols;}


	//Returns the bandwidth of a matrix
	virtual int GetBandWidth() const;

	virtual void GetBandMatrix(Matrix& bm, int lu)
	{
		int col_band = lu*2+1;
		int n = rows;

		bm.SetSize(n,col_band);

		int ind;
		for (int i = 1; i <= n; i++)
		{
			for (int j = 1; j <= col_band; j++)
			{
				ind  = i-(lu+1)+j;
				if (ind >= 1 && ind <= n)
				{
					bm(i,j) = Get(i, ind);
				}
				else bm(i,j) = 0;
			}
		}
	}

	//Returns the column-vector at clomn col
	virtual Vector GetColVec(int col) const;
	virtual void GetColVec(int col, Vector& v) const;

	//Returns the row-vector at row row
	virtual Vector GetRowVec(int row) const;
	virtual void GetRowVec(int row, Vector& v) const;


	//Returns the diagonal values
	virtual Vector GetDiag() const;

	//Set Matrix to diagonal matrix with value v
	virtual void SetDiagMatrix(double v)
	{
#ifndef __QUICKMATH
		release_assert((rows==cols));
#endif
		for (int i=1; i <= rows; i++)
		{
			for (int j=1; j <= rows; j++)
			{
				if (i!=j) Elem(i,j) = 0;
				else Elem(i,j) = v;
			}
		}
	}

	//Returns the Trace
	virtual double Trace() const;

	//Sets the column-vector at clomn col
	virtual void SetColVec(const Vector& x, int col);

	//Sets the row-vector at row row
	virtual void SetRowVec(const Vector& x, int row);

	//Transposes a matrix
	virtual void TpYs();

	virtual void Mult(double x)
	{
		for(int i=0; i< rows*cols; i++) mat[i]*=x;
	}

	//Take Submatrix sm at (sr,sr), multiply with x and add it to local matrix at (r,c), size: (nr x nc)
	virtual void AddSubmatrix(const Matrix& sm, int sr, int sc, int r, int c, int nr, int nc, double x);
	//Set whole Submatrix sm at *this matrix at (r,c)
	virtual void SetSubmatrix(const Matrix& sm, int r, int c);
	virtual void SetSubmatrix(const Matrix& sm, int r, int c, double x);
	virtual void SetSubmatrix(const Matrix3D& sm, int r, int c);

	//Multiplies the cloumn col with val
	virtual void MulCol(int col, double val);

	//Multiplies the row row with val
	virtual void MulRow(int row, double val);


	//Insert (smaller) Matrix at position (row,col), delete entries
	virtual void InsertMatrix(int row, int col, const Matrix& m);

	//Insert (smaller) Matrix at position (row,col), mrows and mcols of matrix m, delete entries
	virtual void InsertMatrix(int row, int col, const Matrix& m, int fromrow, int fromcol, int mrows, int mcols);

	//Add (possibly smaller) Matrix at position (row,col)
	virtual void AddMatrix(int row, int col, const Matrix& m);
	//Add (possibly smaller) Matrix at position given by references
	virtual void AddMatrix(const TArray<int>& rowref, const TArray<int>& colref, const Matrix& m);

	//Adds the vector vec to row row
	//Adds the vector vec to row row
	virtual void AddRowVec(int row, const Vector& vec);

	virtual void AddColVec(int col, const Vector& vec);

	//Adds fromrow multiplied with fact to row row
	virtual void AddRowVec(int fromRow, int toRow, double fact);
	//Adds fromrow multiplied with fact to row row
	virtual void AddRowVec(int fromRow, int toRow, double fact, int fromCol, int toCol);

	//Multiply column col with column colm of matrix m
	virtual double MultCol(const Matrix& m, int mcol, int col);

	//Adds fromrow multiplied with fact to row row
	virtual void AddColVec(int fromCol, int toCol, double fact);

	//Sets the column-vector at clomn col, taken from Matrix m, column mcol
	virtual void SetColVec(const Matrix& m, int mcol, int col);


	//Adds fromrow multiplied with fact to row row, including bandwidth
	virtual void AddRowVecB(int fromRow, int toRow, const double& fact, int bw);

	//Swaps the rows r1 and r2
	virtual void SwapRows(int r1, int r2);

	//Swaps the columns c1 and c2
	virtual void SwapCols(int c1, int c2);

	//Multiplies band-matrix with vector vec
	virtual Vector Bmult(Vector& v) const;


	//Inverts matrix, returns 1 if successfull
	virtual int Invert();

	//Inverts matrix, returns 1 if successfull, faster than invert()
	virtual int Invert2();

	//Inverts matrix, returns 1 if successfull, using LAPACK
	virtual int InvertLapack();

	//Computes an estimate of the reciprocal of the condition number of this matrix, returns 1 if successfull, using LAPACK
	virtual int EstimateReciprocalConditionNumber(double& rcond) const;

	//$ PG 2013-10-30: used for solving an undetermined system (underdetermined: minimum norm solution |q|->min subject to Aq=f; overdetermined: least squares solution |Aq-f|^2->min)
	// f (input) right hand side
	// q (output) solution
	// rcond (input) reciprocal of estimated condition number of A (this matrix), used for recognizing linearly dependendt equations numerically
	// rank (output) computed rank of matrix A, depending on rcond.
	virtual int UndeterminedSystemSolve(const Vector& f, Vector& q, double rcond, int& rank);

	//Solving equation mat*q=f using gaussian elimination method, returns 1 if successfull
	virtual int Solve(const Vector& fv, Vector& q);  //K*q=f-->q

	//Solving equation mat*q=f using LAPACK, returns 1 if successfull
	virtual int SolveLapack(Vector& q);  

	//LU-Decomposition of matrix, returns 1 if successfull, indx is vector of row interchanges
	virtual int LU(IVector& indx);

	//Solve A*x=b, A is decomposed with LU, result is x=b
	virtual int LUBCKSUB(IVector& indx, Vector& b);

	//Solve A*x=b, A is decomposed with LU, result is x=b
	//virtual int LUBCKSUB(IVector& indx, SparseVector& b, const Vector& bdense);

	//Solving equation mat*q=f using gaussian elimination method
	//and considering band-structure, returns 1 if successfull
	virtual int Bsolve(Vector f, Vector& q);  //K*q=f-->q


	//Solving the equation mat*q=f using a band-matrix based on Cholesky-Method, returns 1 if successfull
	virtual int CholeskyBsolve(Vector f, Vector& q);

	virtual int QRDecomposition(Matrix& Q, Matrix& R);

	virtual int Eigenvalues(Vector& ev);
	virtual int EigenvaluesLapack(Vector& lami) const;

	//writes out matrix using MatLab format
	virtual void PrintToMatlab(ostream& os) const;
	virtual void PrintToMatlabSparse(ostream& os, double tol=1e-16);

	//writes out matrix using Maple format
	virtual void PrintToMaple(ostream& os);

	//writes out matrix using Mathematica (Table) format
	virtual mystr PrintToMathematica();




	//Biconjugate gradient method
	virtual int BCG(Vector f, Vector& q, const Vector& qs);

	// Matrix entries are equally distributed values of piecewise bilinear function on the unit square
	// Interpolate returns value of function in point ploc € [-1, 1] x [-1, 1]
	// Matrix entries specify function values, ordered as
	//   [ 1          2         3   ... width ]             [        y=1         ]
	//   [ w+1        w+2       w+3 ... 2w    ]   ----->    [x=-1             x=1]
	//
	//   [ (h-1)w+1  (h-1)w+1  (h-1)w+2  ...  hw]           [        y=-1        ]
	virtual double Interpolate(const Vector2D& ploc) const;
	// Interpolate at given x-coord xloc, 
	// y-coord corresponds to vector entries (at exit, vector length = number of rows of matrix
	void InterpolateX(const double xloc, Vector& interpol_vector) const;

	//Arithmetic operations with 1 parameter
	virtual Matrix& operator+= (const Matrix& m);
	virtual Matrix& operator-= (const Matrix& m);
	virtual Matrix& operator*= (const double& val);

	//$ SW 2013-10-9: added
	//calculates *= A where A is a square matrix (e.g. a rotation matrix) of dimension 3x3, 2x2 or 1x1
	virtual void ApplySqrMat(const MatrixXD &A);

	// Write Matrix to string to export to C  "double intro[rows][cols] = {.,.,.};"
	mystr& MakeString(mystr intro = mystr(""));

	//Arithmetic operations with 2 parameters
	friend Matrix operator+ (const Matrix& m1, const Matrix& m2);
	friend Matrix operator- (const Matrix& m1, const Matrix& m2);
	friend Matrix operator* (const Matrix& mat, const double& val);
	friend Matrix operator* (const double& val, const Matrix& mat);
	friend Matrix operator* (const Matrix& m1, const Matrix& m2);
	friend Vector operator* (const Vector& v, const Matrix& m);
	friend Vector operator* (const Matrix& m, const Vector& v);
	friend void Mult(const Matrix& m, const Vector& v, Vector& res); //computes res=m*v
	friend void Mult(const Matrix& m, const Vector& v, Vector3D& res); //computes res=m*v
	friend void MultTp(const Matrix& m, const Vector& v, Vector3D& res); //computes res=m*v
	friend void Mult(const Vector& v, const Matrix& m, Vector& res); //computes res=m*v
	friend void Mult(const Vector& v1, const Vector& v2, Matrix& res); //computes res=v1*v2^T, where v1 and v2 are interpreted as a column vectors
	friend void MultBW(const Matrix& m, const Vector& v, Vector& res, int bw); //computes res=m*v with bandwidth bw
	friend void Mult(const Matrix& m, const Vector3D& v, Vector& res); //computes res=m*v
	friend void Mult(const Matrix& m, const Vector2D& v, Vector& res); //computes res=m*v
	friend void Mult(const Matrix& m1, const Matrix& m2, Matrix& res); //computes res=m1*m2
	friend void MultTp(const Matrix& m1, const Matrix& m2, Matrix& res); //computes res=m1^T*m2
	friend void MultSym(const Matrix& m1, const Matrix& m2, Matrix& res); //computes res=m1*m2, result is assumed to be symmetric
	friend void MultSymTp(const Matrix& m1, const Matrix& m2, Matrix& res); //computes res=m1^T*m2, result is assumed to be symmetric
	friend void Mult(const MatrixXD& m1, const Matrix& m2, Matrix& res); //computes res=m1*m2
	friend void Mult(const Matrix& m1, const MatrixXD& m2, Matrix& res); //computes res=m1*m2
	friend Matrix operator* (const SparseMatrix& m1, const Matrix& m2);
	friend void Mult(const SparseMatrix& m1, const Matrix& m2, Matrix& res); //computes res=m1*m2

	//writes out matrix with constandt width and factor normalized
	friend ostream& operator<<(ostream& os, const Matrix& m);

	virtual void SetMatrix2n(TArray<Vector2D>& list) // $ MSax 2013-07-09 : added
	{
		SetSize(list.Length(),2);
		Vector2D tmp;
		for (int i=1; i <= list.Length(); i++)
		{
			tmp = list.Get(i);
			mat[2*i-2]=tmp.X();
			mat[2*i-1]=tmp.Y();
		}
	}

	virtual void GetVector2DList(TArray<Vector2D>& list) // $ MSax 2013-07-09 : added
	{
		list.Init();
		for (int i=1; i <= rows; i++)
		{
			list.Add(Vector2D(mat[2*i-2],mat[2*i-1]));
		}
	}

protected:

	double* mat;
	int rows, cols;
	int lalloc;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


template <int data_size>
class ConstMatrix: public Matrix
{
public:

	//Constructors, Destructor
	ConstMatrix()
	{
		InitConstMatrix();
	}

	ConstMatrix(const ConstMatrix& mati)
	{
		InitConstMatrix();

#ifndef __QUICKMATH
		release_assert(data_size >= mati.rows*mati.cols);
#endif
		rows = mati.rows;
		cols = mati.cols;
		for (int i=0; i < rows*cols; i++)
		{
			constdata[i] = mati.mat[i];
		}
	}

	virtual int MaxAllocatedSize() const {return data_size;}
	virtual int IsConstSizeMatrix() const {return 1;}

	ConstMatrix(int rowsi, int colsi, int nofill=0, double fill=0)
	{
#ifndef __QUICKMATH
		release_assert(data_size >= rowsi*colsi);
#endif
		InitConstMatrix(); 
		Resize(rowsi, colsi, nofill, fill);
	}

	ConstMatrix(double x0)
	{
#ifndef __QUICKMATH
		release_assert(data_size >= 1);
#endif
		InitConstMatrix(); 

		rows=1; cols=1;
		mat[0]=x0;
	}

	ConstMatrix(double x11, double x12, double x21, double x22)
	{
#ifndef __QUICKMATH
		release_assert(data_size >= 4);
#endif
		InitConstMatrix(); 

		rows=2; cols=2;
		mat[0]=x11;
		mat[1]=x12;
		mat[2]=x21;
		mat[3]=x22;
	}

	ConstMatrix(double x11, double x12, double x13,
		double x21, double x22, double x23,
		double x31, double x32, double x33)
	{
#ifndef __QUICKMATH
		release_assert(data_size >= 9);
#endif
		InitConstMatrix(); 

		rows=3; cols=3;
		mat[0]=x11;
		mat[1]=x12;
		mat[2]=x13;
		mat[3]=x21;
		mat[4]=x22;
		mat[5]=x23;
		mat[6]=x31;
		mat[7]=x32;
		mat[8]=x33;
	}

	virtual void InitConstMatrix()
	{
		rows=0;
		cols=0;
		mat = &constdata[0];
		lalloc = 0;
	}
	virtual int GetDataSize() const {return data_size;}

private:
	double constdata[data_size];
};




Matrix GetDiag(int n, double val = 1);
double Mises(const Matrix& m);
void PrintMatrix01(const Matrix& m);



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//-----------------------------------------------------------------
//--------------    CLASS SPARSEVECTOR    -------------------------------
//-----------------------------------------------------------------
//Generation of a vector with user-defined length
//Operations: +,-,* (also matrices), <<
class SparseVector
{
public:

	//Constructors, Destructor
	SparseVector(): l(0),vec(NULL),lalloc(0), rowind(NULL), nelem(0), dnull(0.) {};

	SparseVector(int leni, int ilalloc) 
	{
		GenVector(leni, ilalloc);
	};

	SparseVector(const SparseVector& veci)
	{
		GenVector(veci.l, veci.lalloc);
		nelem = veci.nelem;
		if (nelem!=0)
		{
			memcpy(vec,veci.vec,nelem*sizeof(double));
			memcpy(rowind,veci.rowind,nelem*sizeof(int));
		}
	}

	virtual ~SparseVector() 
	{
		Destroy();
	};

	//Generates a vector with length "leni" (leni=0 possible)
	void GenVector(int leni, int ilalloc = 0)
	{
		dnull = 0.;

		l=leni;
		nelem = 0;
		lalloc = ilalloc;

		if (lalloc)
		{
			vec = new double[lalloc];
			rowind = new int[lalloc];
		}
		else
		{
			vec = NULL;
			rowind = NULL;
		}
	}
	void Destroy()
	{
		if (vec!=NULL && lalloc) {delete [] vec; vec=NULL;} 
		if (rowind!=NULL && lalloc) {delete [] rowind; rowind=NULL;} 
	}

	//deletes data!!!
	void SetLen(int len, int ilalloc = 0) 
	{
		if (ilalloc > lalloc)
		{
			Destroy();
			GenVector(len, ilalloc);
		}
		else
		{
			l = len;
		}		
	}

	void CopyFrom(const Vector& v)
	{
		int nonzero = 0;
		for (int i = 0; i < v.Length(); i++)
			if (v[i] != 0.) nonzero++;

		SetLen(v.Length(), nonzero);
		FillWithZeros();

		for (int i = 0; i < v.Length(); i++)
		{
			if (v[i] != 0.) AddEntry(i+1, v[i]);
		}

	}

	int Length() const {return l;}
	int GetLen() const {return l;}

	int NEntries() const {return nelem;}
	int NAlloc() const {return lalloc;}

	//const Referencing access-operator
	const double& operator()(int row) const 
	{
		for (int i=0; i < nelem; i++)
		{
			if (row == rowind[i]+1) return vec[i];
		}
		//static double dnull=0;
		return dnull;
	};

	//Referencing access-operator
	double& operator()(int row)
	{
		int i=0;
		while (row < rowind[i]+1 && i < nelem)
		{
			i++;
		}
		if (row == rowind[i]+1) return vec[i];

		//element does not exist:
		if (nelem == lalloc) //resize?
		{
			int oldlalloc = lalloc;
			if (lalloc == 0) lalloc = 1;
			else lalloc *= 2;

			double* new_vec = new double[lalloc];
			int* new_rowind = new int[lalloc];

			if (oldlalloc != 0)
			{
				memcpy(new_vec,vec,oldlalloc*sizeof(double));
				memcpy(new_rowind,rowind,oldlalloc*sizeof(int));

				delete [] vec;
				delete [] rowind;
			}
			vec = new_vec;
			rowind = new_rowind;
		}

		if (i < nelem)
		{
			for (int j=nelem; j > i; j--)
			{
				rowind[j] = rowind[j-1];
				vec[j] = vec[j-1];
			}
		}
		nelem++;
		vec[i] = 0;
		rowind[i] = row-1;
		return vec[i];
	};



	//read Entry i:
	void GetEntry(int i, int& row, double& d) const 
	{
#ifndef __QUICKMATH
		release_assert((i>0) && (i<=nelem));
#endif
		row = rowind[i-1]+1;
		d = vec[i-1];
	}

	void SetEntry(int i, int row, double d) const 
	{
#ifndef __QUICKMATH
		release_assert((i>0) && (i<=nelem));
#endif
		rowind[i-1] = row-1;
		vec[i-1] = d;
	}

	double& Entry(int i)
	{
#ifndef __QUICKMATH
		release_assert((i>0) && (i<=nelem));
#endif
		return vec[i-1];
	}

	const double& Entry(int i) const 
	{
#ifndef __QUICKMATH
		release_assert((i>0) && (i<=nelem));
#endif
		return vec[i-1];
	}

	//read Entry i:
	void AddEntry(int row, double d)
	{
		if (nelem == lalloc)
		{
			int oldlalloc = lalloc;
			if (lalloc == 0) lalloc = 1;
			else lalloc *= 2;

			double* new_vec = new double[lalloc];
			int* new_rowind = new int[lalloc];

			if (oldlalloc != 0)
			{
				memcpy(new_vec,vec,oldlalloc*sizeof(double));
				memcpy(new_rowind,rowind,oldlalloc*sizeof(int));

				delete [] vec;
				delete [] rowind;
			}
			vec = new_vec;
			rowind = new_rowind;
		}
		rowind[nelem] = row-1;
		vec[nelem] = d;
		nelem++;
	}

	int GetRow(int i) const 
	{
#ifndef __QUICKMATH
		release_assert((i>0) && (i<=nelem));
#endif
		return rowind[i-1]+1;
	}

	void FillWithZeros() 
	{
		nelem = 0;
	}

private:

	//as vector:
	int l;
	double* vec;

	//dynamic:
	int lalloc;
	int nelem;
	int* rowind;

	double dnull;
};












//-----------------------------------------------------------------
//--------------    CLASS SPARSEMATRIX    -------------------------------
//-----------------------------------------------------------------
//Generation of a Matrix with user-defined size
//Operations: +,-,*,<<

//column based sparse matrix!
//element i, j: find x such that colind[x]==j, mat[rowind[i]+x]
//allocated length of row: through rowind[i+1]-rowind;
//used length of row: rowlen
//column index is stored in colind for every element mat

//undefined element:
#define sparseMatrixVoid -1
//increase factor for reallocate:
#define sparseMatrixInc 2

const double sparsetolzero = 1e-20;

class SparseMatrix
{
public:

	//Constructors,destructor
	SparseMatrix()
	{
		Init();
	};

	void GenMat(int lalloc, int rowsi);

	void Destroy()
	{
		if (mat!=NULL) {delete [] mat; mat=NULL;} 
		if (rowind!=NULL) {delete [] rowind; rowind=NULL;} 
		if (rowlen!=NULL) {delete [] rowlen; rowlen=NULL;} 
		if (colind!=NULL) {delete [] colind; colind=NULL;} 
	}

	SparseMatrix(int rowsi, int colsi, int initcols=1)
	{
		if (rowsi*colsi == 0) {Init(); return;}

		if (initcols < 1) {initcols = 1;}

		rows = rowsi;
		cols = colsi;

		lalloc = initcols*rowsi;
		GenMat(lalloc, rowsi);

		for (int i=0; i < rows; i++)
		{
			rowind[i] = i*initcols;
			rowlen[i] = 0;
			for (int j=0; j < initcols; j++)
			{
				colind[rowind[i]+j] = sparseMatrixVoid;
			}
		}
	}

	SparseMatrix(const SparseMatrix& mati)
	{
		lalloc = mati.lalloc;
		rows = mati.rows;
		cols = mati.cols;

		GenMat(lalloc, rows);

		for (int i = 0; i < rows; i++)
		{
			rowind[i] = mati.rowind[i];
			rowlen[i] = mati.rowlen[i];
		}
		for (int i = 0; i < lalloc; i++)
		{
			colind[i] = mati.colind[i];
			mat[i] = mati.mat[i];
		}
	}

	virtual ~SparseMatrix() 
	{
		Destroy();
	};

	SparseMatrix(const Matrix& mati)
	{
		Init();
		CopyFrom(mati);
	}

	void CopyFrom(const Matrix& mati);

	void CopyFrom(const Matrix& mati, int row1, int col1, int row2, int col2);

	void CopyFrom(const SparseMatrix& mati, int row1, int col1, int row2, int col2);

	//generate new matrix from sparse matrix, only from row1 to row2, col1 to col2 (incl.), resort with vector r
	void CopyFrom(const SparseMatrix& mati, int row1, int col1, int row2, int col2, const IVector& r);

	void Init()
	{
		rows=0;
		cols=0;
		mat=NULL;
		rowind = NULL;
		rowlen = NULL;
		colind = NULL;
		lalloc = 0;
		dnull = 0.;
	};

	void GetMatrix(Matrix& m) const;

	Matrix GetMatrix() const
	{
		//very slow, only for printing!!!!!
		Matrix m;
		GetMatrix(m);
		return m;
	}

	void SetSize(int rowsi, int colsi, int minrowsize = 1)
	{
		int resize = 0;
		if (rows!=rowsi || cols!=colsi) resize = 1;
		else
		{
			int i=0;
			while (!resize && i<rows) 
			{
				if (RowLenMax0(i) < minrowsize) resize = 1;
				i++;
			}
		}

		if (resize)
		{
			if (rowsi == 0 || colsi ==0)
			{
				rows = rowsi;
				cols = colsi;
				return;
			}

			int rowsize = lalloc/rowsi; //round off!
			if (rowsize < minrowsize) rowsize = minrowsize;

			int neededmem = rowsize*rowsi;

			if (neededmem > lalloc || rowsi > rows)
			{
				Destroy();

				lalloc = neededmem;
				GenMat(lalloc, rowsi);
			}
			rows = rowsi;
			cols = colsi;

			int cnt = 0;
			for (int i = 0; i < rows; i++)
			{
				rowind[i] = cnt;
				cnt += rowsize;
				rowlen[i] = 0;
			}
		}
	}

	void SetSizePerColumn(int rowsi, int colsi, const TArray<int>& size_per_column)
	{
		//first it could be checked if resize is necessary!
		if (rowsi == 0 || colsi ==0)
		{
			rows = rowsi;
			cols = colsi;
			return;
		}

		int neededmem = 0;
		for (int i=1; i<=size_per_column.Length(); i++)
		{
			neededmem += size_per_column(i);
		}

		if (neededmem > lalloc || rowsi > rows)
		{
			Destroy();

			lalloc = neededmem;
			GenMat(lalloc, rowsi);
		}
		rows = rowsi;
		cols = colsi;

		int cnt = 0;
		for (int i = 0; i < rows; i++)
		{
			rowind[i] = cnt;
			cnt += size_per_column(i+1);
			rowlen[i] = 0;
		}
	}

	Vector GetMaxlenVector() const
	{
		Vector v(rows);
		for (int i=0; i < rows; i++)
			v.Elem0(i) = RowLenMax0(i);
		return v;
	}

	int CheckMatrix() const
	{
		int n = 0;
		for (int i = 0; i < rows; i++)
		{
			for (int j = 1; j < rowlen[i]; j++)
			{
				if (colind[rowind[i]+j] <= colind[rowind[i]+j-1] ||
					(colind[rowind[i]+j] > cols-1)) n++;
			}
		}
		return n;
	}

	int CountZeroEntries() const
	{
		int n = 0;
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < rowlen[i]; j++)
			{
				if (mat[rowind[i]+j] == 0) n++;
			}
		}
		return n;
	}

	int CountEntries() const
	{
		int n=0;
		for (int i = 0; i < rows; i++)
		{
			n += rowlen[i];
		}
		return n;
	}

	int EliminateZeroEntries(double mindouble=0);

	void PrintData() const;

	int GetMaxRowlen() const
	{
		int n=0;
		for (int i = 0; i < rows; i++)
		{
			n = Maximum(rowlen[i],n);
		}
		return n;
	}

	//give each row twice space if more than 50% is used
	void ReAllocate();

	//Access-functions for private data
	int Getrows() const {return rows; };
	int Getcols() const {return cols; };
	int GetLAlloc() const {return lalloc;};

	//Assignment-operator
	SparseMatrix& operator= (const SparseMatrix& mati)
	{
		if (this == &mati) {return *this;}

		if (mati.lalloc > lalloc || mati.rows > rows)
		{
			Destroy();

			if (rows!=0 && cols!=0)
			{
				lalloc = mati.lalloc;

				GenMat(lalloc, mati.rows);
			}
		}

		rows=mati.rows;
		cols=mati.cols;

		if (mati.mat==NULL || rows*cols==0)
		{
			mat = NULL;
		}
		else
		{
			memcpy(mat,mati.mat,rows*cols*sizeof(double));
			memcpy(colind,mati.colind,rows*cols*sizeof(int));
			memcpy(rowind,mati.rowind,rows*sizeof(int));
			memcpy(rowlen,mati.rowlen,rows*sizeof(int));
		}
		return *this;
	}

	//fill with zeros, keep structure
	void FillWithZeros()
	{
		for (int i=0; i < rows; i++)
		{
			rowlen[i] = 0;
		}
	}

	//Referencing access-operator on element using row- and column-values
	double& operator()(int row, int col)
	{
#ifndef __QUICKMATH
		release_assert((row>0) && (row<=rows) && (col>0) && (col<=cols));
#endif
		//increase when necessary???
		int i = 0;
		while (i < rowlen[row-1] && colind[rowind[row-1]+i] < col-1) {i++;}

		if (i < rowlen[row-1] && colind[rowind[row-1]+i] == col-1)
		{
			return mat[rowind[row-1]+i];
		}
		else
		{
			int ind = i;

			if (rowlen[row-1] == RowLenMax0(row-1)) ReAllocate();

			for (i = rowlen[row-1]-1; i >= ind; i--)
			{
				colind[rowind[row-1]+i+1] = colind[rowind[row-1]+i];
				mat[rowind[row-1]+i+1] = mat[rowind[row-1]+i];
			}
			rowlen[row-1]++;
			mat[rowind[row-1]+ind] = 0;
			colind[rowind[row-1]+ind] = col-1;
			return mat[rowind[row-1]+ind];
		}
	};


	//Referencing constant access-operator on element using row- and column-values
	const double& operator()(int row, int col) const
	{
#ifndef __QUICKMATH
		release_assert((row>0) && (row<=rows) && (col>0) && (col<=cols));
#endif
		for (int i=0; i < rowlen[row-1]; i++)
		{
			if (colind[rowind[row-1]+i] == col-1)
				return mat[rowind[row-1]+i];
		}
		//static double dnull=0;
		return dnull;
	};

	//Referencing constant access-operator on element using row- and column-values
	const double& Get(int row, int col) const
	{
#ifndef __QUICKMATH
		release_assert((row>0) && (row<=rows) && (col>0) && (col<=cols));
#endif
		for (int i=0; i < rowlen[row-1]; i++)
		{
			if (colind[rowind[row-1]+i] == col-1)
				return mat[rowind[row-1]+i];
		}
		//static double dnull=0;
		return dnull;
	};

	//Referencing constant access-operator on element using row- and column-values
	//uses last sparsecolindex of same row for faster access
	const double& Get(int row, int col, int& lastsparsecol) const
	{
#ifndef __QUICKMATH
		release_assert((row>0) && (row<=rows) && (col>0) && (col<=cols));
#endif
		int i=lastsparsecol-1;
		while (i < rowlen[row-1] && colind[rowind[row-1]+i] < col-1) i++;

		lastsparsecol = i+1; //i+1 for faster incremental access;
		//if (lastsparsecol >1 ) lastsparsecol--;

		if (i < rowlen[row-1] && colind[rowind[row-1]+i] == col-1) return mat[rowind[row-1]+i];

		//static double dnull=0;
		return dnull;
	};

	//maximum number of elements in row-1
	int RowLenMax0(int row) const
	{
		if (row < rows-1) 
		{
			return rowind[row+1]-rowind[row];
		}
		else 
		{
			return lalloc-rowind[row];
		}	
	}

	//Returns the maximum-norm (largest value in matrix)
	double MaxNorm() const
	{
		double res=0;
		for(int row = 0; row < rows; row++)
		{
			for(int i = 0; i < rowlen[row]; i++)
			{
				res = Maximum(fabs(mat[rowind[row]+i]),res);
			}
		}
		return res;
	}

	//Returns the quadratic-norm
	double Norm2() const
	{
		double res=0;
		for(int row = 0; row < rows; row++)
		{
			for(int i = 0; i < rowlen[row]; i++)
			{
				res += Sqr(mat[rowind[row]+i]);
			}
		}
		return sqrt(res);
	}

	double FillInFact() const
	{
		if (rows*cols == 0) return 0;
		int used = 0;

		for (int i = 0; i < rows; i++)
		{
			used += rowlen[i];
		}
		return (double)used/(double)(rows*cols);
	}

	//Returns the lower and upper bandwith of matrix
	void GetBandWidths(int& lb, int& ub) const;

	int GetBandWidth() const 
	{
		int a,b;
		GetBandWidths(a,b);
		return 1+Maximum(a,b);
	}

	void GetBandMatrix(Matrix& bm, int lu) const
	{
		if (0) //frage: warum funktioniert die erste version nicht???
		{
			int col_band = lu*2+1;
			int n = rows;

			bm.SetSize(n,col_band);

			int ind;
			int lastsparseind;
			for (int i = 1; i <= n; i++)
			{
				lastsparseind = 1;
				for (int j = 1; j <= col_band; j++)
				{
					ind  = i-(lu+1)+j;
					if (ind >= 1 && ind <= n)
					{
						bm(i,j) = Get(i, ind);
						//bm(i,j) = Get(i, ind, lastsparseind);
					}
					else bm(i,j) = 0;
				}
			}
		}
		else
		{
			int col_band = lu*2+1;
			int n = rows;

			bm.SetSize(n,col_band);
			bm.FillWithZeros();

			int ind, k, row;
			for (int i = 1; i <= n; i++)
			{
				row = i-1;
				int rind = rowind[row];
				for (int j = 0; j < rowlen[row]; j++)
				{
					k = colind[rind+j]+1;
					ind  = k-i+(lu+1);
					if (k >= 1 && k <= n && ind >= 1 && ind <= col_band)
					{
						bm(i,ind) = mat[rind+j];
					}
				}
			}
		}
	}

	void Mult(double x)
	{
		for(int row = 0; row < rows; row++)
		{
			for(int i = 0; i < rowlen[row]; i++)
			{
				mat[rowind[row]+i] *= x;
			}
		}
	}

	//Add (possibly smaller) Matrix m of size (lrow x lcol) at reference position
	void AddMatrix(const TArray<int>& rowref, const TArray<int>& colref, int lrow, int lcol, const Matrix& m);

	//Take Submatrix sm at (sr,sr), multiply with x and add it to local matrix at (r,c), size: (nr x nc)
	void AddSubmatrix(const Matrix& sm, int sr, int sc, int r, int c, int nr, int nc, double x);

	//Take Submatrix sm at (sr,sc), multiply with x and add it to local matrix at (r,c), size: (nr x nc)
	void AddSubmatrix(const SparseMatrix& sm, int sr, int sc, int r, int c, int nr, int nc, double x);
	//set column vector, only entries of v from vr1 to vr2, insert at column col
	void SetColVector(const Vector& v, int vr1, int vr2, int col);

	void SetColVector(const SparseVector& v, int col);

	SparseMatrix& operator*= (const double& val)
	{
		//global_uo << "operator *= for sparse matrices not yet tested\n";
		for(int row = 0; row < rows; row++)
		{
			for(int i = 0; i < rowlen[row]; i++)
			{
				mat[rowind[row]+i] *= val;
			}
		}
		return *this;
	}

	SparseMatrix& operator+= (const SparseMatrix& m)
	{
		//global_uo << "operator += for sparse matrices not yet tested\n";
		AddSubmatrix(m,1,1,1,1,rows,cols,1);
		return *this;
	}

	SparseMatrix& operator-= (const SparseMatrix& m)
	{
		AddSubmatrix(m,1,1,1,1,rows,cols,-1);
		return *this;
	}

	/*
	//Arithmetic operations with 1 parameter
	Matrix& operator-= (const Matrix& m);

	friend void Mult(const SparseMatrix& m1, const Matrix& m2, Matrix& res); //computes res=m1*m2
	friend void Mult(const Matrix3D& m1, const SparseMatrix& m2, Matrix& res); //computes res=m1*m2
	*/

	friend Matrix operator* (const SparseMatrix& m1, const Matrix& m2);
	friend void Mult(const SparseMatrix& m, const Vector& v, Vector& res); //computes res=m*v
	friend void Mult(const SparseMatrix& m1, const Matrix& m2, Matrix& res); //computes res=m1*m2

	friend void Mult(const SparseMatrix& m1, const SparseMatrix& m2, SparseMatrix& res); //computes res=m1*m2
	virtual void ComputeTranspose(SparseMatrix& res) const; //computes res=m1^T
	friend double MultRow(const SparseMatrix& m1, const SparseMatrix& m2, int row1, int row2); //mult row1 of m1 times row2 of m2


	void PrintToMatlabSparse(ostream& os);

	//writes out matrix with constandt width and factor normalized
	friend ostream& operator<<(ostream& os, const SparseMatrix& m);

	void GetMatPtrs(int*& colindI, int*& rowptrI, double*& matI) const
	{
		colindI = colind;
		rowptrI = rowind;
		matI = mat;
	}
	int* GetRowLenArray() const {return rowlen;}

	virtual int RowLen(int i) const {return rowlen[i-1];}
	virtual double GetRowEntry(int row, int i) const {return mat[rowind[row-1]+i-1];}
	virtual int GetRowColind(int row, int i) const {return colind[rowind[row-1]+i-1];}

private:

	int rows, cols;
	int lalloc;

	double* mat;
	int* colind;
	int* rowind;
	int* rowlen;

	double dnull;
};

//$ PG 2013-10-2: deleted class HarwellBoeingMatrix
class HarwellBoeingMatrix { public:	HarwellBoeingMatrix() {assert(0 && "ERROR: HarwellBoeingMatrix removed");} };

//struct SuperMatrix;


// AP
class SuperLUInverse
{
public:
	SuperLUInverse() : sparseA(0), rows(0), cols(0), perm_r(0), perm_c(0), isFirstStep(1), isFactorized(0),
		SuperA(0), SuperL(0), SuperU(0)
	{
		;
	}

SuperLUInverse:: SuperLUInverse(const SparseMatrix & matA) : sparseA(&matA), isFactorized(0), isFirstStep(1),
		SuperA(0), SuperL(0), SuperU(0)
	{
		sparseA =&matA;
		rows = matA.Getrows();
		cols = matA.Getcols();
		
		isFactorized = 0;
		isFirstStep = 1;
		perm_r = new int(rows);
		perm_c = new int(cols);
	}


	~SuperLUInverse();

	void SetSparseMatrix(const SparseMatrix& matA)
	{
		sparseA =&matA;
		rows = matA.Getrows();
		cols = matA.Getcols();
		
		isFactorized = 0;
		isFirstStep = 1;
		delete[] perm_r;
		delete[] perm_c;
		perm_r = new int(rows);
		perm_c = new int(cols);
	}

	int Solve(Vector& q);
	int Factorize();
	int SuperLUApply(Vector& q);
	int Apply(Vector& q) { return SuperLUApply(q);}

	int& IsFirstStep() {return isFirstStep; }
	const int& IsFirstStep() const {return isFirstStep; }


private:
	// pointer to original sparsematrix
	const SparseMatrix* sparseA;
	// size, has to be quadratic
	int rows, cols;
	// SuperLU-Data
	SuperMatrix* SuperA;
	SuperMatrix* SuperU;
	SuperMatrix* SuperL;
	int* perm_r;
	int* perm_c;
	int isFirstStep;
	// is factorization done?
	int isFactorized;

};
// END AP


//solve A*q = rhs; q contains right-hand-side on input and the solution on output
int SuperLUSolve(SparseMatrix& A, Vector& q);
int SuperLUFactorize(SparseMatrix& A, SuperMatrix*& SuperA, SuperMatrix*& L, SuperMatrix*& U, int*& perm_r, int*& perm_c, int isFirstStep = 1);
int SuperLUApply(SparseMatrix& A,     SuperMatrix*& SuperA, SuperMatrix*& L, SuperMatrix*& U, int*& perm_r, int*& perm_c, Vector& q);


//$ DR 2013-10-04:[ the following functions are not available anymore
//int banddec(int n, int ld, int ud, Matrix& pmat, IVector& perm, int&  signd);
//int bandsol(int n, int ld, int ud, Matrix& pmat, Vector& b, const IVector& perm);
//int Band_LU_Dec(Matrix& a, int n, int m1, int m2, Matrix& al, IVector& indx, double& d);
//void Band_LU_Bks(const Matrix& a, int n, int m1, int m2, const Matrix& al, const IVector& indx, Vector& b);
//void Band_LU_Bks0(const Matrix& a, int n, int m1, int m2, const Matrix& al, const IVector& indx, Vector& b);
//$ DR 2013-10-04:]

class PardisoInverse;

// The SparseInverse uses Pardiso or SuperLU internally
class SparseInverse
{
public:
	SparseInverse() : pardisoinv(NULL), superluinv(NULL)
	{
	}

	SparseInverse(const SparseMatrix & matA) ;

	~SparseInverse();
	void SetSparseMatrix(const SparseMatrix& matA);
	int Solve(Vector& q);
	int Factorize();
	int Apply(Vector& q);
	int& IsFirstStep() ;
	const int& IsFirstStep() const ;

private:
	SuperLUInverse* superluinv;
	PardisoInverse* pardisoinv;
};


//symmetrical smoothing
void SmoothenDoubleArray(Vector& values, Vector& result, int smoothing_length);


// this function computes the local maxima and minima values if a vector 'data' from left to right (in case of similar points, the first point from left to right is used - see '°' below)
// if only max-values are computed, set flag compute_minima to zero and use dummys for TArrays 'min_indices', and 'min_values'
//
//max1
//o    max2
// \    °--             maxN 
//  \  /   \    max3    °
//   °      °   °   .../
//           \ / \ /
//            °   °
//        min1  min2
void GetExtremeValues(const Vector& data, TArray<int>& max_indices,TArray<double>& max_values, TArray<int>& min_indices,TArray<double>& min_values, int compute_minima=1);