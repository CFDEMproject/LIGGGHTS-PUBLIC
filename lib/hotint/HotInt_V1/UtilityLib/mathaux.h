//#**************************************************************
//#
//# filename:             mathaux.h
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

#include <math.h>

const int MYMAXINT = 2147483647; //change for 64 bit ...
const double MY_PI = 3.1415926535897932384;

inline int IsNaN(double x)
{
	return !(x==x);
}


inline double Norm3D(double x, double y, double z);
inline double Maximum(double a, double b, double c)
{
	if (a > b && a > c) return a;
	else if (b > c) return b;
	else return c;
}
inline double Minimum(double a, double b, double c)
{
	if (a < b && a < c) return a;
	else if (b < c) return b;
	else return c;
}

template <class T>
inline T Maximum(const T& a, const T& b)
//berechnet das Maximum zweier Objekte
{
	return (a > b) ? a : b;
}

template <class T>
inline T Minimum(const T& a, const T& b)
//berechnet das Minimum zweier Objekte
{
	return (a < b) ? a : b;
}
int IsNaN(double x);

template <class T>
inline void Swap(T & a, T & b)
//vertauscht zwei Objekte
{
	T temp = a;
	a = b;
	b = temp;
}

template <class T>
inline int Sgn(T a)
//wendet die Signum-Funktion auf ein Objekt an
{
	if (a > 0) return 1;
	if (a < 0) return -1;
	return 0;
}

inline double Sqr(double a)
//berechnet das Quadrat eines Objektes
{ return a * a; }

inline double sqr(double a)
//berechnet das Quadrat eines Objektes
{ return a * a; }

template <class T>
inline T Cub(T a)
//berechnet das a^3
{ return a * a * a; }

inline double To3(double a)
//berechnet das a^3
{ 
	return a * a * a; 
}

inline double To4(double a)
//berechnet das a^4
{ 
	double a2 = a*a;
	return a2 * a2; 
}

inline double To5(double a)
//berechnet das a^4
{ 
	double a2 = a*a;
	return a2 * a2 * a; 
}

inline double To6(double a)
//berechnet das a^4
{ 
	double a2 = a*a;
	return a2 * a2 * a2; 
}


inline double Fact(int x)
{
	double f=1;
	for (int i=2; i<=x; i++)
	{
		f*=i;
	}
	return f;
}

inline double MOtoX(double x)
//berechnet (-1)^x, x wird gerundet!!! (x ganzzahlig)
{ 
	if ((int)x % 2 == 0) {return 1;}
	else {return -1;}
}

inline double MyPow(double x,double y)
//Berechnet x hoch y
{ 
	if (y == 2.) {return x*x;}
	if (y == -1.) {return 1/x;}
	else {return pow(x,y);}
}
inline double RelApproxi(double d1, double d2, double acc = 1e-8)
{
	double val = fabs(d1)+fabs(d2);
	if (val <= acc) val=1;
	if (fabs(d1-d2)/val < acc) return 1;
	return 0;
}

int GenVecCnt();
int GenMatCnt();
int NewtonIts();

Matrix GetRot(double phi);
void GetRot(double phi, Matrix& m);
void GetRotT(double phi, Matrix& m);

void Rotate(Vector2D& v, double phi); 

Matrix GetRot1(double p);
Matrix GetRot2(double p);
Matrix GetRot3(double p);
Matrix GetRot(double p1,double p2,double p3);

const int int_to_bit_list[18] = {0,1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536};

inline int GetIntBit(int i) {return int_to_bit_list[i];}
inline char GetCharBit(int i) {return (char)int_to_bit_list[i];}

const int NOFILL = 0;


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ abstract base class for integerX, derive integer arrays of fixed length from this class +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! no object of this class my exist, only pointer to this class on derived object is valid !

class TIntX 
{
	static const int length = 0; // length defined in baseclass & derived class, GetLen() operator takes values from derived class

public: // constructor & destructor
	TIntX() {;}
	virtual ~TIntX() {;}

public: // functions valid for all derived classes, 

//Sets all entries to zero
	virtual void Zeros() { for (int i=1; i <= GetLen(); i++) Get(i) = 0; }

//Inverts the order if the IntX
//keeps first elements fixed
	virtual void Invert() { for (int i=2; i <= ((GetLen()+1)/2); i++) Swap(i,GetLen()+2-i); }

//Set function which takes any number of arguments 
//surplus entries are truncated, missing remain uninitialized
	int Set(int i1,...)
	{
		va_list ilist;
		va_start(ilist, i1);
		Get(1)=i1;
		for (int i=2; i<=GetLen(); i++) Get(i) = va_arg(ilist,int);	
		va_end(ilist);
		return GetLen();
	}

//Set function which copies from an existing array 
//surplus entries are truncated, missing remain uninitialized
	virtual int SetV(const int* ii) { for (int i=1; i<=GetLen(); i++) Get(i) = ii[i-1]; return GetLen(); } 

//copies data from an other derived class object
//discards tail of sorce if source is longer, 
//does NOT assign values to remaining entries if source is shorter
	virtual int CopyFrom(const TIntX* source)
	{
		int min = source->GetLen();
		if(GetLen() < min) min = GetLen();
		for(int i=1; i <= min; i++) Get(i) = source->Get(i);
		return min;
	}

	virtual int CopyFrom(const TIntX& source)
	{
		int min = source.GetLen();
		if(GetLen() < min) min = GetLen();
		for(int i=1; i <= min; i++) Get(i) = source.Get(i);
		return min;
	}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ some conversion functions between derived of TIntX, derived classes must be declared above that
//+	-> cpp file																																
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*
// copy from one derived to another
	virtual TIntX* GetCopy() const
	{ // can not use new TIntX ...
		TIntX* copy;
		switch(GetLen())
		{
			case 2: copy = new int2(data[0],data[1]); break;
			case 3: copy = new int3(data[0],data[1],data[2]); break;
			case 4: copy = new int4(data[0],data[1],data[2],data[3]); break;
			default: break;
		}
		return copy;
	}

// eliminate double entries
	virtual TInt* Shorten()
	{
		TIntX* short;
		switch(length)
		{
			case 3: 
			case 4:
			default: break;
				}
		return short;
	}
	
	*/

// assignment operator
	virtual TIntX& operator= (const TIntX& other)
	{
		if (&other == this) return *this;
		this->CopyFrom((&other));
		return *this;
	}

// access operators
	virtual int& operator[](int i) { return Get(i); }
	virtual const int& operator[](int i) const { return Get(i); }
	virtual int& operator()(int i) { return Get(i); }
	virtual const int& operator()(int i) const { return Get(i); }
	virtual int& GetMod(int i) { return Get(((i-1)%GetLen())+1); }
	virtual const int& GetMod(int i) const { return Get(((i-1)%GetLen())+1); }

	virtual int& Get(int i) = 0; 
	virtual const int& Get(int i) const = 0; 
	
	virtual int GetLen() {return length;}
	virtual const int GetLen() const {return length;} 

// identity operator
	virtual int operator== (const TIntX& other) 
	{
		for (int i=1; i<= other.GetLen(); i++)
		{
			if (Get(i) != other.Get(i)) return 0;
		}
		return 1;
	}

// sort function - own QS(TIntX*)...
	virtual void Sort()
	{
//Sorts an array a(1..a.GetLen()) into ascending numerical order by Shell’s method (diminishing increment
//sort). a is replaced on output by its sorted rearrangement. Normally, the argument n should
//be set to the size of array a, but if n is smaller than this, then only the first n elements of a
//are sorted. This feature is used in selip.
		int n = GetLen();
		int i,j,inc;
		int v;
		inc=1; //Determine the starting increment.
		do 
		{
			inc *= 3;
			inc++;
		} while (inc <= n);

		do 
		{ //Loop over the partial sorts.
			inc /= 3;
			for (i=inc+1;i<=n;i++) 
			{ //Outer loop of straight insertion.
				v=Get(i);
				j=i;
				while (Get(j-inc) > v) 
				{ //Inner loop of straight insertion.
					Get(j)=Get(j-inc);
					j -= inc;
					if (j <= inc) break;
				}
				Get(j)=v;
			}
		} while (inc > 1);
	}

// find (first) occurance of number "f" in list
// returns position or zero
	virtual int Find(int f) const 
	{
		for(int i=1; i <= GetLen(); i++)
		{
			if (Get(i) == f) return (i);
		}
		return 0;
	}
	
// DONT USE FOR LARGE LENGTH
// returns number of redundant entries 
// flag not set: counts number of redundant entries -> tripple entry == 2 redundant
// flag set: counts pairs -> tripple entry == 3 pairs
	virtual int CountRedundantEntries(int flag_pairs = 0) const
	{
		int counter = 0;
		for (int i=2; i <= GetLen(); i++)
		{
			for (int j=1; j <= i-1; j++)
			{
				if ( Get(i) == Get(j) ) 
				{
					counter ++;
					if (flag_pairs == 0) 
						break; // break loop, only raise count once, not for every pair... 
				}
			}
		}
		return counter;
	}

// counts entries of the "tail" number at the end of the sequence
	virtual int CountTailEntries(int tail = 0)
	{
		int counter = 0;
		for (int i = GetLen(); i>1 ; i++)
		{
			if (Get(i) == tail) 
				counter ++;
			else 
				break;
		}
		return counter;
	}

// identify double entries
// returns position of n-th redundant entry or 0 if not found (not that many double entries)
	virtual int IdentifyRedundantEntry(int n)
	{
		int counter = 0;
		for (int i=2; i <= GetLen(); i++)
		{
			for (int j=1; j <= i-1; j++)
			{
				if ( Get(i) == Get(j) ) 
				{
					counter ++;
					if (counter == n) return i; 
					if (counter > n) return 0;
					break;
				}
			}
		}
	return 0;
	}

// remove double entry
// kicks out i-th entry, moves following elements forward, fill tail with chosen number
	virtual void RemoveEntry(int i, int tail = 0)
	{
		if (i <= 0) return;
		if (i > GetLen()) return;
		if (i == GetLen()) { Get(GetLen()) = tail; return; }
		MoveBlock(i+1,GetLen(),-1,tail);
	}

// move block of entries definde by (first,last) n steps to the right (positive x), fill with chosen number
// (AD) needs to be revised for several non standard cases where array limits are exceeded...
	virtual void MoveBlock(int first, int last, int n, int tail = 0)
	{
		if (first > last) return; 
		if (first < 1) return;
		if (last > GetLen()) return;
		if (n == 0) return;

		if (n > 0) // >>
		{
			for (int i=last; i >= first+n; i--) // highest allowed destination is last, lowest allowed source is first
				Get(i) = Get(i-n);
			for (int i=Minimum(last,first+n-1); i >= first; i--) // highest allowed destination is Min(first+n,last), lowest allowed destination is first
				Get(i) = tail;
		}
		if (n < 0) // <<  
		{
			for (int i=Maximum(1,first+n); i <= last+n; i++) // lowest allowed destination is 1, highest allowed source is last 
				Get(i) = Get(i-n);
			for (int i=Maximum(first,last+n+1); i <= last; i++) // lowest allowed destination is max(last+n+1,first), hightest allowed destination is last
				Get(i) = tail;
		}
	}

// remove all double entries
// conserves class, eliminates all multiple entries, sets end to value tail
// returns number of unique entries
	virtual int RemoveAllRedundantEntries(int tail = 0)
	{
		int len = GetLen();
		int count = CountRedundantEntries();
		int tailc = CountTailEntries(tail);
		for (int i=1; i<= count; i++)
			this->RemoveEntry(this->IdentifyRedundantEntry(1),tail); // always first double entry where there are still "numbers", not "tail"...
		return GetLen()-count-tailc;
	}

	virtual int IsCyclicEqual(const TIntX& other) const
	{
		if(GetLen() == other.GetLen()) // same length
		{
			int i,j;
			for(j=0; j<GetLen(); j++) // offset
			{
				for(i=1; i<=GetLen(); i++) // loop elem
				{
					if(Get(i) != other.GetMod(i+j)) break;
				}
				if(i==(GetLen()+1)) return 1;
			}
		}
		else return 0;// not same length
		return 0;
	}

//Swaps entries 
	virtual void Swap(int i, int j) { int dummy = Get(i); Get(i) = Get(j); Get(j)=dummy; }
	
private: // pure virtual function Abstractor: 
				 // MUST override in derived class (or only derived* allowed)
				 // use     private: void Abstractor() {;};           in derived class
				 // then you can make an actual object of that derived class

//this function makes TIntX an Abstract Base Class (ABC)
	virtual void Abstractor() = 0;

protected:
//	int length; 
//	int* data; 
};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ classes derived from TintX, these classes will replace current int2, int3 and int4 soon +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// new derived class MUST have:
// overloaded constructors: (own class&)==CC, TIntX*, TIntX&
// overloaded Abstractor()
// own accessfunctions to be TArray-save (TArray::Resize (memcpy) & TIntX::data (int*) don't mix)
// new derived class SHOULD have:
// overloaded Initialization: Init(), Zeros()
// overloaded Inversion: Invert()
// overloaded copy assignment operator (CAO):   operator =

class int2 : public TIntX
{
	static const int length = 2;
public:
	int2() { Zeros(); }
	int2(int i1,int i2) { Set(i1,i2); }
	int2(const int2& source) { idata[0]=source.idata[0]; idata[1]=source.idata[1]; } //CC(1)
	int2(const TIntX* source) { Zeros(); CopyFrom(source); }
	int2(const TIntX& source) { Zeros(); CopyFrom(&source); }
	~int2() {;} //D(2)

	virtual int& Get(int i) { return idata[i-1]; }
	virtual const int& Get(int i) const { return idata[i-1]; }
	virtual int GetLen() {return length;}
	virtual const int GetLen() const {return length;} 

	void Zeros() {Set(0,0);};
//	virtual int2& operator= (const int2& other) //CAO(3)
//	{
//		if (&other == this) return *this;
//		idata[0]=other.idata[0]; idata[1]=other.idata[1]; return *this;
//	}
	void Invert() {;};
private:
	void Abstractor() {;};
	int idata[2];

public: // from old class int2
	void Swap() { TIntX::Swap(1,2); }
	void Sort() { if (idata[1]<idata[0]) Swap(); }
};

class int3 : public TIntX
{
	static const int length = 3;
public:
	int3() { Zeros(); }
	int3(int i1,int i2, int i3) { Set(i1,i2,i3); }
	int3(const int3& source) { idata[0]=source.idata[0]; idata[1]=source.idata[1]; idata[2]=source.idata[2]; }
	int3(const TIntX* source) { CopyFrom(source); }
	int3(const TIntX& source) { CopyFrom(&source); }
	~int3() {;}

	virtual int& Get(int i) { return idata[i-1]; }
	virtual const int& Get(int i) const { return idata[i-1]; }
	virtual int GetLen() {return length;}
	virtual const int GetLen() const {return length;} 

	void Zeros() {Set(0,0,0);};
	void Invert() { Swap(2,3); }
private:
	void Abstractor() {;};
	int idata[3];

public: // from old class int3
	virtual void ChangeOrder() { Swap(1,3); }
};

class int4 : public TIntX
{
	static const int length = 4;
public:
	int4() { Zeros();}
	int4(int i1,int i2, int i3, int i4) { Set(i1,i2,i3,i4);}
	int4(const int4& source) { idata[0]=source.idata[0]; idata[1]=source.idata[1]; idata[2]=source.idata[2]; idata[3]=source.idata[3]; }
	int4(const TIntX* source) { CopyFrom(source);}
	int4(const TIntX& source) { CopyFrom(&source);}
	~int4() {;}

	virtual int& Get(int i) { return idata[i-1]; }
	virtual const int& Get(int i) const { return idata[i-1]; }
	virtual int GetLen() {return length;}
	virtual const int GetLen() const {return length;} 

	void Zeros() {Set(0,0,0,0);};
	void Invert() {Swap(2,4);};
private:
	void Abstractor() {;};
	int idata[4];

public: // from old class int4
	virtual void ChangeOrder() {Swap(2,4);}
	virtual int IsEqual(int4 other) { return (*this == other); }
};


// copy of old classes int2, int3, int4
// data already declared private...
/*
class int2
{
public:
	int2() {ii[0]=0; ii[1]=0;}
	int2(const int2& i2) {ii[0] = i2.Get(1); ii[1]=i2.Get(2);}
	int2(int a1, int a2) {ii[0]=a1; ii[1]=a2;}
	int2& operator= (const int2& i2)
	{
		if (&i2 == this) return *this;
		
		ii[0] = i2.Get(1); 
		ii[1] = i2.Get(2);

		return *this;
	}

	int operator== (const int2& i2) const
	{
		return (i2.Get(1) == ii[0] && i2.Get(2) == ii[1]);
	}
	void Sort() 
	{
		if (ii[1] < ii[0]) Swap();
	}

	void Swap() 
	{
		int i = ii[0];
		ii[0] = ii[1];
		ii[1] = i;
	}

	virtual int& operator[](int elem) { return ii[elem]; }
	virtual const int& operator[](int elem) const { return ii[elem]; }
	virtual int& operator()(int elem) { return ii[elem-1]; }
	virtual const int& operator()(int elem) const { return ii[elem-1]; }
	virtual int& Get(int elem) { return ii[elem-1]; }
	virtual const int& Get(int elem) const { return ii[elem-1]; }
	virtual int& GetMod(int elem) { return ii[(elem-1)%2]; }
	virtual const int& GetMod(int elem) const { return ii[(elem-1)%2]; }

//private:
	int ii[2];
//	int i1,i2;
};
*/
/*
class int3
{
public:
	int3() {ii[0]=0; ii[1]=0; ii[2]=0;}
	int3(int a1, int a2, int a3) {ii[0]=a1; ii[1]=a2; ii[2]=a3;}

	int3(const int3& i3) {ii[0] = i3.ii[0]; ii[1]=i3.ii[1]; ii[2]=i3.ii[2];}
	int3& operator= (const int3& i3)
	{
		if (&i3 == this) return *this;
		
		ii[0] = i3.ii[0]; 
		ii[1] = i3.ii[1];
		ii[2] = i3.ii[2];

		return *this;
	}

	virtual int Find(int i) const 
	{
		if (i == ii[0]) return 1;
		else if (i == ii[1]) return 2;
		else return 3;
	}

	virtual void ChangeOrder()
	{
		Swap(ii[0],ii[2]);
	}
	
	virtual int& operator[](int elem)
	{
#ifndef __QUICKMATH
		assert((elem>=0) && (elem<=2));
#endif
		return ii[elem];
	};
	virtual const int& operator[](int elem) const
	{
#ifndef __QUICKMATH
		assert((elem>=0) && (elem<=2));
#endif
		return ii[elem];
	};

	virtual int& operator()(int elem)
	{
#ifndef __QUICKMATH
		assert((elem>=1) && (elem<=3));
#endif
		return ii[elem-1];
	};
	virtual const int& operator()(int elem) const
	{
#ifndef __QUICKMATH
		assert((elem>=1) && (elem<=3));
#endif
		return ii[elem-1];
	};

	virtual int& Get(int elem)
	{
#ifndef __QUICKMATH
		assert((elem>=1) && (elem<=3));
#endif
		return ii[elem-1];
	};
	virtual const int& Get(int elem) const
	{
#ifndef __QUICKMATH
		assert((elem>=1) && (elem<=3));
#endif
		return ii[elem-1];
	};

	virtual int& GetMod(int elem)
	{
		return ii[(elem-1)%3];
	}

	virtual const int& GetMod(int elem) const
	{
		return ii[(elem-1)%3];
	}
	virtual int IsCyclicEqual(int3 other) const
	{
		if (ii[0]==other.ii[0] &&
				ii[1]==other.ii[1] &&
				ii[2]==other.ii[2]) return 1;
		else if (ii[0]==other.ii[1] &&
				ii[1]==other.ii[2] &&
				ii[2]==other.ii[0]) return 1;
		else if (ii[0]==other.ii[2] &&
				ii[1]==other.ii[0] &&
				ii[2]==other.ii[1]) return 1;

		return 0;
	}

private:
	int ii[3];
};
*/
/*
class int4
{
public:
	int4() {};
	int4(int v1, int v2, int v3, int v4)
	{ii[0] = v1; ii[1] = v2; ii[2] = v3; ii[3] = v4; }

	virtual int& operator()(int elem)
	{
#ifndef __QUICKMATH
		assert((elem>=1) && (elem<=4));
#endif
		switch (elem)
		{
		case 1: return ii[0]; break;
		case 2: return ii[1]; break;
		case 3: return ii[2]; break;
		case 4: return ii[3]; break;
		default: ;
		}
		return ii[0]; //code should not be reached!!!
	};
	virtual const int& operator()(int elem) const
	{
#ifndef __QUICKMATH
		assert((elem>=1) && (elem<=4));
#endif
		switch (elem)
		{
		case 1: return ii[0]; break;
		case 2: return ii[1]; break;
		case 3: return ii[2]; break;
		case 4: return ii[3]; break;
		default: ;
		}
		return ii[0]; //code should not be reached!!!
	};
	
	virtual int& Get(int elem)
	{
#ifndef __QUICKMATH
		assert((elem>=1) && (elem<=4));
#endif
		switch (elem)
		{
		case 1: return ii[0]; break;
		case 2: return ii[1]; break;
		case 3: return ii[2]; break;
		case 4: return ii[3]; break;
		default: ;
		}
		return ii[0]; //code should not be reached!!!
	};
	virtual const int& Get(int elem) const
	{
#ifndef __QUICKMATH
		assert((elem>=1) && (elem<=4));
#endif
		switch (elem)
		{
		case 1: return ii[0]; break;
		case 2: return ii[1]; break;
		case 3: return ii[2]; break;
		case 4: return ii[3]; break;
		default: ;
		}
		return ii[0]; //code should not be reached!!!
	};

	virtual void Invert()
	{
		int i = ii[1];
		ii[1] = ii[3];
		ii[3] = i;
	}
	virtual void ChangeOrder() //==Invert()
	{
		int i = ii[1];
		ii[1] = ii[3];
		ii[3] = i;
	}
	//return 1 if both int4 are exactly the same
	virtual int IsEqual(int4 other)
	{
		if (ii[0]==other.ii[0] &&
				ii[1]==other.ii[1] &&
				ii[2]==other.ii[2] &&
				ii[3]==other.ii[3]) return 1;
		return 0;
	}
	//return 1 if both int4 are the same or if any cyclic change is the same
	virtual int IsCyclicEqual(int4 other)
	{
		if (ii[0]==other.ii[0] &&
				ii[1]==other.ii[1] &&
				ii[2]==other.ii[2] &&
				ii[3]==other.ii[3]) return 1;
		else if (ii[0]==other.ii[1] &&
						 ii[1]==other.ii[2] &&
						 ii[2]==other.ii[3] &&
						 ii[3]==other.ii[0]) return 1;
		else if (ii[0]==other.ii[2] &&
						 ii[1]==other.ii[3] &&
						 ii[2]==other.ii[0] &&
						 ii[3]==other.ii[1]) return 1;
		else if (ii[0]==other.ii[3] &&
						 ii[1]==other.ii[0] &&
						 ii[2]==other.ii[1] &&
						 ii[3]==other.ii[2]) return 1;

		return 0;
	}

private:
	int ii[4];
//	int i1,i2,i3,i4;
};
*/

class double2
{
public:
	double2() {d1=0;d2=0;}
	double2(double d1i, double d2i) {d1 = d1i; d2 = d2i;}

	void Swap() {double d3 = d1; d1 = d2; d2 = d3;}

	double d1;
	double d2;
};

int IsEqual(const IVector& v1, const IVector& v2);

//find entry in vector v; if found, return pos, else return 0
int Find(int i, const IVector& v);

//interpolate several values by means of cosine functions
double InterpolateSmooth(double t, double t0, double t1, double A);

//smooth interpolation over time t of function values f0,f1,f2,f3 at times t0,t1,t2,t3
double InterpolateFnt2(double t, double t0, double t1, double t2, double t3, double f0, double f1, double f2, double f3);

//smooth interpolation over time t of function values f0,f1,f2,f3 at times t0,t1,t2,t3
double InterpolateFntArray(double t, TArray<double> t_vec, TArray<double> f_vec);

inline double Heaviside(const double& t)
//applies the Heaviside function to variable t
{
	if (t < 0) return 0.;
	return 1.;
}

