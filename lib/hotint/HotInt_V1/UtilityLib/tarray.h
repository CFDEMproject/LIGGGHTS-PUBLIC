//#**************************************************************
//#
//# filename:             tarray.h
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						19.06.97
//# description:          Dynamical array class
//# remarks:						  Array enlarges dynamically its size, trying to avoid many new commands
//#												index runs from 1 to n, the used data size is length (Length()), the available
//#												space for data is size (Size()).
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


#ifndef TARRAY__H
#define TARRAY__H

#pragma warning(disable: 4996)		//$ YV 2012-12-27: to avoid complaints concerning depricated functions




#include "ioincludes.h"

#include "release_assert.h"
#include <string.h>
#include  <stdio.h>
#include  <stdarg.h> 


/////////////////////////
//PG: Here we define the strategy for asserting tarray code.
//We still have to discuss if this is the optimal setting.
//Meanwhile, I suggest: it's best to assert in debug-mode
//and if the preprocessor flag "__ASSERT_IN_RELEASE_MODE__"
//is set (e.g., in "MBSKernelLib\preprocessor_includes.h").
//This is due to the fact, that some of the models may
//only be executed in release-mode.
#ifdef  _DEBUG

#define __ASSERT_TARRAY

#else  /*_DEBUG*/

#ifdef  __ASSERT_IN_RELEASE_MODE__
#define __ASSERT_TARRAY
#endif  /*__ASSERT_IN_RELEASE_MODE__*/

#endif  /*_DEBUG*/
/////////////////////////



template <class T>
class TArray
{

public:

	//Konstruktoren, Destruktor:
	TArray() { data = NULL; size = 0; length = 0; inc = 0; }

	TArray(int sizei, int inci = 0)
	{
		data = NULL;
		length = 0;
		size = sizei;
		inc = inci;
		if (size>0) {data = new T[size];}
	}

	TArray(const TArray<T> &t)
	{
		data = NULL;
		size = 0;
		length=0;
		inc=t.inc;
		CopyFrom(t);
	}

	virtual ~TArray () 
	{
		if (data) {delete[] data;}
	}

	TArray & operator = (const TArray &t)
		//Zuweisungsoperator
	{
		if (this==&t) {return *this;}
		inc=t.inc;
		CopyFrom(t);
		return *this;
	}

	void CopyFrom(const TArray<T>& t, int from=0, int to=0);
	//kopiert einen Teil des Objektes vom Startindex from bis zum Index to
	//bei form=to=0 wird das gesamte Objekt kopiert
	void ReSize (int minsize);
	void Init() {data = NULL; size = 0; length = 0; inc = 0;}
	//das Array wird auf eine Mindestgröße vergrößert
	T& operator[] (int i);  //(PG): deprecated -- use operator() instead
	T& operator() (int i);
	T& Elem(int i);
	const T& Get(int i) const;
	T& Elem0(int i);
	const T& Get0(int i) const;
	//referenzierter Zugriffsoperator
	const T& operator[] (int i) const;	//(PG): deprecated -- use operator() instead
	const T& operator() (int i) const;

	T& Last();
	const T& Last() const;

	//referenzierter konstanter Zugriffsoperator
	int Length() const { return length; }
	//gibt die Länge des Arrays zurück
	void SetLen(int nl) { length = nl; if (length > size) {ReSize(length);};}
	//setzt die Länge des Arrays

	// set to length 1 and set value
	void Set1(T t1) { SetLen(1); data[0] = t1; }
	// set to length 2 and set values
	void Set2(T t1, T t2) { SetLen(2); data[0] = t1; data[1]=t2; }
	// set to length 3 and set values
	void Set3(T t1, T t2, T t3) { SetLen(3); data[0] = t1; data[1]=t2; data[2] = t3; }
	// set to length 4 and set values
	void Set4(T t1, T t2, T t3, T t4) { SetLen(4); data[0] = t1; data[1]=t2; data[2] = t3; data[3] = t4; }
	// set to length 6 and set values
	void Set6(T t1, T t2, T t3, T t4, T t5, T t6) { SetLen(6); data[0] = t1; data[1]=t2; data[2] = t3; data[3] = t4; data[4] = t5; data[5] = t6; }

	// set function with variable number of arguments
	// set to length n and then set n values
	// !not tested! for non-fundamental types (AD)
	void SetXN(int n = 0, ...)
	{
		if (n == 0) { return Init(); }

		va_list ilist;
		va_start(ilist, n);

		SetLen(n);
		for (int i=1; i<=n; i++)
		{
			data[i-1] = va_arg(ilist,T);	
		}
		va_end(ilist);
	}

	int GetItemsInContainer() const { return length; }
	//wie len, kompatibel zum Windows-Array TArrayAsVector
	int Size() { return size; }
	//ermittelt die tatsächlich im Speicher vorhandene Länge des Arrays
	int Add(const T& t)
		//fügt ein Element am Schluß der Liste hinzu und erhöht die Länge
	{
		length++;
		if (length > size) {ReSize(length);}
		data[length-1] = t;
		return length;
	}
	void SetAll(const T& t)
	{
		for (int i=0; i < length; i++)
		{
			data[i] = t;
		}
	}

	void CopyFrom(const T* t, int len)
		//copy data from C-Array; NOT TESTED!!!
	{
		SetLen(0);
		for (int i=0; i < len; i++)
		{
			Add(t[i]);
		}
	}
	void Insert(int i, const T& t)
		//fügt ein Element an der Stelle i ein, verschiebt die hinteren Elemente
	{
		length++;
		if (length > size) {ReSize(length);}
		for (int j=length; j > i; j--)
		{
			data[j-1] = data[j-2];
		}
		data[i-1] = t;
	}
	void Erase(int i)
		//löscht das Element an der Stelle i, verschiebt die hinteren Elemente; wenn i==0, tue nichts->mit find kombinierbar!
	{
		if (i == 0) return;
		for (int j=i; j < length; j++)
		{
			data[j-1] = data[j];
		}
		length--;
	}

// deletes all elements marked in array flags (1 -> delete, 0 -> keep)
// array flags is forced to same size as this
	void EraseMany(TArray<int>& flags)
	{
		flags.SetLen(Length()); // force same length...
		
		int newnum = 0;
		for(int i=1; i<= Length(); i++)
		{
			if(flags(i) == 0) // dont delete
			{			
				newnum++;
				if(newnum < i)
					data[newnum-1] = data[i-1];
			}			
		}
		SetLen(newnum);
	}

// Rearrange the TArray as specified in mask 
// mask holds the new index numbers of the items
// mask is expected to: contain all numbers 1 .. m or 0 if element should be deleted
// mask must be same length as this TArray !! mask must be ascending (no repetitions allowed)
	int RearrangeArray(TArray<int>& mask) 
	{
// compute length of rearranged TArray
		int highest = 0;
		for(int i=1; i<= length; i++)
		{
			if(mask(i) != 0)
			{
				if(mask(i) > highest)
				{
					highest = mask(i); // this may cause empty entries in new array: sequence must be complete! easiest way: ascending (or 0)
				}
				else 
					mask(i) = 0; // this sets all repetitions to zero 
			}
		}
// make copy & delete original
		TArray copy_of_this(*this);
		this->Flush();
		this->SetLen(highest);
// reassemble in new order
		for(int i=1; i <= copy_of_this.Length(); i++)
		{
			if(mask(i) != 0)
			{
			
				(*this)(mask(i)) = copy_of_this(i);
				//Add(copy_of_this(i));
			}
		}
		return highest;
	}

	int Find(const T& elem) const
		//find elem in TArray
	{
		for (int j=0; j < length; j++)
		{
			if (data[j] == elem) return j+1;
		}
		return 0;
	}

	//insert an element, if it does not already exist in array
	int AddIfNotExists(const T& t)
	{
		int f = Find(t);
		if (f) return f;

		length++;
		if (length > size) {ReSize(length);}
		data[length-1] = t;
		return length;
	}


	void Merge(const TArray<T>& t, int from=0, int to=0);
	//fügt ein anderes Array zum aktuellen Array hinzu

	void Flush();
	//löscht die gesamte Liste. Speicher wird nur von Liste freigegeben, 
	//    nicht aber daten der elemente wenn die daten pointer waren

	void Append(const TArray<T>& t, int from=0, int to=0);
	//fügt ein anderes Array zum aktuellen Array hinzu

// set operations - methods & operators, 
	TArray<T>& Union(const TArray<T>& other);															// attention: this routine uses operator new 
  TArray<T>& Intersection(const TArray<T>& other);											// attention: this routine uses operator new
	TArray<T>& Difference(const TArray<T>& other);												// attention: this routine uses operator new

	void operator+= (const TArray<T>& other); 														// attention: this operator alters the original array
	void operator&= (const TArray<T>& other);													    // attention: this operator alters the original array
	void operator-= (const TArray<T>& other);															// attention: this operator alters the original array

	template<class T>
	friend TArray<T> operator+ (const TArray<T>& a, const TArray<T>& b);  // attention: this operator copies array on the stack
	template<class T>
	friend TArray<T> operator& (const TArray<T>& a, const TArray<T>& b);  // attention: this operator copies array on the stack
	template<class T>
	friend TArray<T> operator- (const TArray<T>& a, const TArray<T>& b);  // attention: this operator copies array on the stack


	T* GetDataPtr() {return data;}

protected:
	T * data;
	int size, inc, length;
};

////partial specialization of template class: Flush for T pointer
//template <class T>
//class TArray<T *>
//{
//public:
//	void Flush();
//		
//protected:
//	T ** data;
//	int size, inc, length;
//};

template <class T>
ostream& operator<<(ostream& os, const TArray<T>& t)
{
	os << "[";
	for (int i=1; i<=t.Length(); i++)
	{
		os << t[i] << " ";
	}
	os << "]" << flush;
	return os;
}

template <class T>
void TArray<T> :: CopyFrom(const TArray<T>& t, int from, int to)
{
	if (from==0) {from=1;}
	if (to==0) {to=t.length;}
	//if (from==0 && to==0) {from=1; to=t.length;}
	if (t.length == 0) {length=0; return;}
	if ((to-from+1) > size) ReSize(to-from+1);
	memcpy (data, &t.data[from-1], (to-from+1) * sizeof (T));
	length=to-from+1;
}

template <class T>
T & TArray<T> :: operator[] (int i)
{
	if (i > size) {ReSize(i);}
	if (i > length) {length=i;}
	return data[i-1];
}

template <class T>
const T & TArray<T> :: operator[] (int i) const
{
#ifdef __ASSERT_TARRAY
	release_assert(i>0 && i<=size);
#endif
	return data[i-1];
}

template <class T>
T & TArray<T> :: operator() (int i)
{
	if (i > size) 
	{ReSize(i);}
	if (i > length) {length=i;}
	return data[i-1];
}

template <class T>
const T & TArray<T> :: operator() (int i) const
{
#ifdef __ASSERT_TARRAY
//#ifdef _DEBUG
	release_assert(i>0 && i<=size);
#endif
	return data[i-1];
}

template <class T>
T & TArray<T> :: Elem(int i)
{
	if (i > size) {ReSize(i);}
	if (i > length) {length=i;}
	return data[i-1];
}

template <class T>
const T & TArray<T> :: Get(int i) const
{
#ifdef __ASSERT_TARRAY
	release_assert(i>0 && i<=size);
#endif
	return data[i-1];
}

template <class T>
T & TArray<T> :: Elem0(int i)
{
	if (i >= size) {ReSize(i+1);}
	if (i >= length) {length=i+1;}
	return data[i];
}

template <class T>
const T & TArray<T> :: Get0(int i) const
{
#ifdef __ASSERT_TARRAY
	release_assert(i>=0 && i<size);
#endif
	return data[i];
}

template <class T>
T & TArray<T> :: Last()
{
	return data[length-1];
}

template <class T>
const T & TArray<T> :: Last() const
{
	return data[length-1];
}

template <class T>
void TArray<T> :: Merge(const TArray<T>& t, int from, int to)
{
	if (from==0 && to==0) {from=1; to=t.length;}
	//ReSize(length+to-from+1);
	//int pos=length+1;
	//for(int i=0; i<=to-from; i++)
	//{
	//	(*this)[pos+i]= t[i+from];
	//}
//$ AD 2011-04: changed: for many consecutive merges the ReSize did always double size :)
	if (from < 1) from = 1;
	if (to > t.length) to = length;
	for(int i=from; i <= to; i++)
	{
		Add(t(i));
	}
}

template <class T>
void TArray<T> :: Append(const TArray<T>& t, int from, int to)
{
	if (from==0 && to==0) {from=1; to=t.length;}
	if (from < 1) from = 1;
	if (to > t.length) to = length;
	for(int i=from; i <= to; i++)
	{
		Add(t(i));
	}
}

template <class T>
TArray<T>& TArray<T> :: Union(const TArray<T>& other) 
{
	TArray<T>* setunion = new TArray<T>(*this);
	*setunion += other;
	return *setunion;
}

template <class T>
TArray<T>& TArray<T> :: Intersection(const TArray<T>& other) 
{
	TArray<T>* setinter = new TArray<T>(*this);
  *setinter &= other;
	return *setinter;
}

template <class T>
TArray<T>& TArray<T> :: Difference(const TArray<T>& other)
{
  TArray<T>* setdiff = new TArray<T>(*this);
	*setdiff -= other;
	return *setdiff;
}

template <class T>
void TArray<T>::operator+= (const TArray<T>& other)               // attention: this operator alters the original array
{
	for(int j=1; j<=other.Length(); j++)
	{
		if(!Find(other(j)))
			Add(other(j));
	}
}

template <class T>
void TArray<T>::operator&= (const TArray<T>& other)               // attention: this operator alters the original array
{
	for(int i=length; i>0; i--)
	{
		TArray<T>& nco = (TArray<T>&) other;
		const T& elem = this->Get(i);
		if(!nco.Find(elem))
			Erase(i);
	}
}

template <class T>
void TArray<T>::operator-= (const TArray<T>& other)               // attention: this operator alters the original array
{
	for(int i=length; i>0; i--)
	{
		TArray<T>& nco = (TArray<T>&) other;
		const T& elem = this->Get(i);
		if( nco.Find(elem))
			Erase(i);
	}
}

template<class T>
TArray<T> operator+ (const TArray<T>& a, const TArray<T>& b)

{
	TArray<T> c(a);
	c += b;
	return c;
}

template <class T>
TArray<T> operator & (const TArray<T>& a, const TArray<T>& b)
{
	TArray<T> c(a);
	c &= b;
	return c;
}

template <class T>
TArray<T> operator - (const TArray<T>& a, const TArray<T>& b)
{
	TArray<T> c(a);
	c -= b;
	return c;
}


template <class T>
void TArray<T> :: Flush()
{
	if (data) delete [] data;
	data = NULL;
	length = 0;
	size = 0;
	inc = 0;
}

template <class T>
void TArray<T> :: ReSize(int minsize)
{
	T* ndata;
	if (minsize==0) {return; }
	int oldsize = size;

	size = (inc) ? size + inc : 2 * size;
	if (size < minsize) {size = minsize;}

	ndata = new T[size];

	if (data!=NULL && oldsize!=0)
	{
		memcpy (ndata, data, oldsize * sizeof (T));
		//    memcpy (ndata, data, length * sizeof (T));
	}

	if (data!=NULL)
		delete [] data;

	data = ndata;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Integer Arrays that are initialized with a constant size (with constructors)
template<int n>
class IntVectorWithLengthN : public TArray<int> // for all possible lengths - TArray MUST contain Integers
{
public:
	IntVectorWithLengthN():TArray(n) { SetAll(0); }
	void SetArithmetic(int a1, int d, int n = 0) 
	{ 
	  if(n>0) SetLen(n); // new length
		int add = 0;
		for(int i=1; i<=Length(); i++)
		{		
			Elem(i) = a1 + add;
			add += d;
		}
	}
	void SetGeometric(int b1, int q, int n = 0)
	{
	  if(n>0) SetLen(n); // new length
		int mult = 1;
		for(int i=1; i<=Length(); i++)
		{		
			Elem(i) = b1 * mult;
			mult *= q;
		}
	}
};

class IntVec1 : public IntVectorWithLengthN<1> // length 1, data type integer, length can be changed as in normal TArray
{
public:
	IntVec1():IntVectorWithLengthN() {;}
	IntVec1(int i1) { Set1(i1); }
};

class IntVec2 : public IntVectorWithLengthN<2> // length 2, data type integer, length can be changed as in normal TArray
{
public:
	IntVec2():IntVectorWithLengthN() {;}
	IntVec2(int i1, int i2) { Set2(i1,i2); }
};

class IntVec3 : public IntVectorWithLengthN<3> // length 3, data type integer, length can be changed as in normal TArray
{
public:
	IntVec3():IntVectorWithLengthN() {;}
	IntVec3(int i1, int i2, int i3) { Set3(i1,i2,i3); }
};

// special case: sequence of natural numbers
class NaturalNumbers : public TArray<int> 
{
public:
	NaturalNumbers():TArray(0) {}
	NaturalNumbers(int n):TArray(0) 
	{
		for(int i=1; i<=n; i++)
		{	
			Add(i);
		}
	}
};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Double Arrays that are initialized with a constant size (with constructors)
template<int n>
class DblVectorWithLengthN : public TArray<double> // for all possible lengths - TArray MUST contain doubles
{
public:
	DblVectorWithLengthN():TArray(n) { SetAll(0.); }
	void SetArithmetic(double a1, double d, int n = 0) 
	{ 
	  if(n>0) SetLen(n); // new length
		double add = 0.;
		for(int i; i<=Length(); i++)
		{		
			Elem(i) = a1 + add;
			add += d;
		}
	}
	void SetGeometric(double b1, double q, int n = 0)
	{
	  if(n>0) SetLen(n); // new length
		double mult = 1.;
		for(int i; i<=Length(); i++)
		{		
			Elem(i) = b1 * mult;
			mult *= q;
		}
	}
};

class DblVec1 : public DblVectorWithLengthN<1> // length 1, data type double, length can be changed as in normal TArray
{
public:
	DblVec1():DblVectorWithLengthN() {;}
	DblVec1(double d1) { Set1(d1); }
};

class DblVec2 : public DblVectorWithLengthN<2> // length 2, data type double, length can be changed as in normal TArray
{
public:
	DblVec2():DblVectorWithLengthN() {;}
	DblVec2(double d1, double d2) { Set2(d1,d2); }
};

class DblVec3 : public DblVectorWithLengthN<1> // length 3, data type double, length can be changed as in normal TArray
{
public:
	DblVec3():DblVectorWithLengthN() {;}
	DblVec3(double d1, double d2, double d3) { Set3(d1,d2,d3); }
};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Array for Templates that dynamically copy memory by calling constructors, derived from TArray<T>
//requires implemented operator = in template class!

// tested (models_alex.cpp) for mystr
// vector, matrix not working 100%: problem with operator delete data[i]; 

template<class T>
class TArrayDynamic : public TArray<T>
{
public:
	TArrayDynamic():TArray() {};
	TArrayDynamic(int sizei, int inci = 0)
	{
		data = NULL;
		length = 0;
		size = sizei;
		inc = inci;
		if (size>0) {data = new T[size];}
	}
	virtual ~TArrayDynamic ()
	{
		if (data) delete []data;
		//if (data) 
		//{
		//	for(int i=0; i<length; i++)
		//	delete data[i];
		//}
		data = NULL; // important to set to NULL otherwise deleted again by base class destructor
	}

	TArrayDynamic & operator = (const TArrayDynamic &t)
	{
		if (this==&t) {return *this;}
		inc=t.inc;
		CopyFrom(t);
		return *this;
	}

	void CopyFrom(const TArrayDynamic<T>& t, int from=0, int to=0);
	void ReSize (int minsize);
	int Add (const T& t);
	int Insert (int i, const T& t);
};

template <class T>
void TArrayDynamic<T> :: CopyFrom(const TArrayDynamic<T>& t, int from, int to)
{
	if (from==0 && to==0) {from=1; to=t.length;}
	if (t.length == 0) {length=0; return;}
	if ((to-from+1) > size) ReSize(to-from+1);

	for (int i=from; i <= to; i++)
	{
		data[i-from] = t.data[i-1];
	}
	length=to-from+1;
}

template <class T>
int TArrayDynamic<T> :: Add(const T& t)
		//fügt ein Element am Schluß der Liste hinzu und erhöht die Länge
{
	length++;
	if (length > size) {ReSize(length);}
	data[length-1] = t;
	return length;
}

template <class T>
int TArrayDynamic<T> :: Insert(int i, const T& t)
		//fügt ein Element am Schluß der Liste hinzu und erhöht die Länge
{
length++;
if (length > size) {ReSize(length);}
	for (int j=length; j > i; j--)
	{
		data[j-1] = data[j-2];
	}
	data[i-1] = t;
return length;
}

template <class T>
void TArrayDynamic<T> :: ReSize(int minsize)
{
	T* ndata;
	if (minsize==0) {return; }
	int oldsize = size;

	size = (inc) ? size + inc : 2 * size;
	if (size < minsize) {size = minsize;}

	ndata = new T[size];

	if (data!=NULL && oldsize!=0)
	{
		for (int i=0; i < oldsize; i++)
		{
			ndata[i] = data[i];
		}
	}

	if (data!=NULL)
		delete [] data;

	data = ndata;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Array of variable size, lines and columns
//Get(i,j) or Elem(i,j) accesses the i-th line, j-th column
template <class T>
class TMatrix
{
public:

	//Konstruktoren, Destruktor:
	TMatrix(): data() {nlines=1; ncols=1; };


	TMatrix(int li, int ci): data(li, 0)
	{
		nlines = li;
		ncols = ci;
	};

	virtual ~TMatrix () {
		data.Flush();
	};

	int SetDim(int li, int ci)
	{
		data.Flush();
		nlines = li;
		ncols = ci;
		for (int i=0; i<li; i++)
		{
			TArray<T>* newline = new TArray<T>(ci,0);
			newline->SetLen(ci);
			data.Add(newline);	
		}
		return data.Length();
	}

	//number of lines
	int Length() const {return data.Length();}
	int NLines() const {return data.Length();}
	//number of cols in i-th line
	int EntryLength(int i) const 
	{
#ifdef __ASSERT_TARRAY
		release_assert(i>0 && i<=data.Length());
#endif
		return data.Get(i)->Length();
	}
	int NCols(int i) const 
	{
#ifdef __ASSERT_TARRAY
		release_assert(i>0 && i<=data.Length());
#endif
		return data.Get(i)->Length();
	}

	void GenLine(int i)
	{
#ifdef __ASSERT_TARRAY
		release_assert(i>0);
#endif
		if (i>data.Length())
		{
			for (int k=data.Length(); k <= i; k++)
			{
				TArray<T>* newline = new TArray<T>(ncols,0);
				data.Add(newline);
			}
		}
	}

	//Get read-access to the element of the i-th line, j-th column
	const T& Get(int i, int j) const 
	{
#ifdef __ASSERT_TARRAY
		release_assert(i>0 && i<=data.Length());
		release_assert(j>0 && j<=data.Get(i)->Length());
#endif
		return data.Get(i)->Get(j);
	}

	//Get write-access the element of i-th line, j-th column
	T& Elem(int i, int j)
	{
		GenLine(i);

#ifdef __ASSERT_TARRAY
		release_assert(j>0);
#endif
		return data.Elem(i)->Elem(j);
	}

	//Get read-access to the element of the i-th line, j-th column
	const T& operator()(int i, int j) const 
	{
#ifdef __ASSERT_TARRAY
		release_assert(i>0 && i<=data.Length());
		release_assert(j>0 && j<=data.Get(i)->Length());
#endif
		return data.Get(i)->Get(j);
	}

	//Get write-access the element of i-th line, j-th column
	T& operator()(int i, int j)
	{
#ifdef __ASSERT_TARRAY
		release_assert(i>0);
#endif
		GenLine(i);

#ifdef __ASSERT_TARRAY
		release_assert(j>0 && j<=data.Get(i)->Length());
#endif
		return data.Elem(i)->Elem(j);
	}

	// MaSch 08/2012:
	//Get const reference (read access) to i-th line (i.e. a TArray<T>)
	const TArray<T>& GetLine(int i) const 
	{
#ifdef __ASSERT_TARRAY
		release_assert(i>0 && i<=data.Length());
#endif
		return *data.Get(i);
	}

	// MaSch 08/2012:
	//Get const reference (read access) to i-th line (i.e. a TArray<T>)
	const TArray<T>& operator()(int i) const 
	{
#ifdef __ASSERT_TARRAY
		release_assert(i>0 && i<=data.Length());
#endif
		return *data.Get(i);
	}


	//add element to i-th line:
	void Add(int i, const T& t) 
	{
#ifdef __ASSERT_TARRAY
		release_assert(i>0);
#endif
		GenLine(i);
		data.Elem(i)->Add(t);
	}

	void Flush()
	{
		for (int i=1; i <= data.Length(); i++)
		{
			data.Elem(i)->Flush();
			delete data.Elem(i);
		}
		data.Flush();
	}

	void print()
	{
		for (int i=1; i <= data.Length(); i++)
		{
			for (int j=1; j <= data.Get(i)->Length(); j++)
			{
				::cout << (*this)(i,j) << "  ";
			}
			::cout << "\n";
		}

	}

private:
	TArray<TArray<T>*> data;
	int nlines;
	int ncols;
};





//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//cstr list
class MyStrList
{
public:
	MyStrList(int initlen=0): strlist(initlen)
	{
		for (int j=1; j <= initlen; j++) strlist(j) = 0;
	}

	MyStrList(const MyStrList& ed):strlist()
	{
		CopyFrom(ed);
	};
	MyStrList& operator=(const MyStrList& ed) 
	{
		if (this == &ed) {return *this;}
		CopyFrom(ed);
		return *this;
	}
	virtual MyStrList* GetCopy() const
	{
		MyStrList* ec = new MyStrList();
		ec->CopyFrom(*this);
		return ec;
	}

	virtual void CopyFrom(const MyStrList& e)
	{
		strlist.SetLen(0);		
		for (int i=1; i <= e.strlist.Length(); i++)
		{
			size_t length = strlen(e.strlist(i));
			char* str = new char[length + 1];
			strcpy(str, e.strlist(i));

			strlist.Add(str);
		}
	}

	virtual ~MyStrList()
	{
		Reset();
	}

	virtual void Reset()
	{
		for (int i=1; i <= strlist.Length(); i++)
		{
			if (strlist(i) != 0) delete [] strlist(i);
			strlist(i) = 0;
		}
		strlist.SetLen(0);
	}

	virtual int Find(const char* name) const 
	{
		for (int i=1; i <= strlist.Length(); i++)
		{
			if (strcmp(name, strlist(i)) == 0) return i;
		}
		return 0;
	}
	virtual int Length() const {return strlist.Length();}
	virtual int Add(const char* s) 
	{
		size_t length = strlen(s);
		char* str = new char[length + 1];
		strcpy(str, s);

		return strlist.Add(str);
	}
	virtual const char* Get(int i) const 
	{
#ifdef __ASSERT_TARRAY
		release_assert(i>0 && i<=Length());
#endif
		return strlist(i);
	}
	virtual char* Get(int i) 
	{
#ifdef __ASSERT_TARRAY
		release_assert(i>0 && i<=Length());
#endif
		return strlist(i);
	}
	const char* operator() (int i) const
	{
#ifdef __ASSERT_TARRAY
		release_assert(i>0 && i<=Length());
#endif
		return strlist(i);
	}
	char* operator() (int i)
	{
#ifdef __ASSERT_TARRAY
		release_assert(i>0 && i<=Length());
#endif
		return strlist(i);
	}

	virtual void Set(int i, const char* s)
	{
		size_t length = strlen(s);
		char* str = new char[length + 1];
		strcpy(str, s);

		if (i > Length())
		{
			for (int j=Length(); j <= i; j++) strlist(j) = 0;
		}
		strlist(i) = str;
	}
	virtual const char* Last() const {return strlist(strlist.Length());}
	virtual void Delete(int i) 
	{
		if (strlist(i) != 0) delete [] strlist(i);
		for (int j=i; j < strlist.Length(); j++)
		{
			strlist(j) = strlist(j+1);
		}
		strlist.SetLen(strlist.Length()-1);
	}

private:
	TArray<char*> strlist; 
};




#endif
