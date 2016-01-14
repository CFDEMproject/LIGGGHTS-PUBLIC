//#***************************************************************************************
//# filename:     femathHelperFunctions.h
//#
//# author:				Yury Vetyukov, Johannes Gerstmayr
//# 
//# generated:      
//# description:  
//#                       
//# comments:      
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
//#***************************************************************************************
#pragma once

#include <time.h>
#include <sys/timeb.h>


// timers

const int TMlen = 30;

#define timeron

double GetClockTime();
void TMResetTimer();
void TMPrintTimer();
#ifdef timeron
void TMStartTimer(const int& i);
void TMStopTimer(const int& i);
#else
inline void TMStartTimer(const int& i) {};
inline void TMStopTimer(const int& i) {};
#endif
double TMGetTimer(int i);

struct UserOutputInterface;

// this function tries to open a file and returns 1 if the file was OK, otherwise 0; 
// if pUO != NULL, warning text message appears if file was not found, the same with the log file
int DoesFileExist(mystr filename, UserOutputInterface * pUO = NULL, ofstream * logFile = NULL);

const int lalloc_act = 1; //assures that lalloc is utilized (for decreasing the size of a Matrix without new memory allocation)



//Sorts an array a(1..a.Length()) into ascending numerical order by Shell’s method (diminishing increment
//sort). a is replaced on output by its sorted rearrangement. b is sorted according to the order of a.

void QuicksortDouble(TArray<double>& a, TArray<int>& b);
//Sorts an array a(1..a.Length()) into ascending numerical order by Shell’s method (diminishing increment
//sort). a is replaced on output by its sorted rearrangement. b is sorted according to the order of a.

void Quicksort(TArray<int>& a, TArray<int>& b, TArray<double>& c);
//Sorts an array a(1..a.Length()) into ascending numerical order by Shell’s method (diminishing increment
//sort). a is replaced on output by its sorted rearrangement. b and c are sorted according to the order of a.

void QuicksortDouble(TArray<double>& a, TArray<double>& b, TArray<int>& c);
//Sorts an array a(1..a.Length()) into ascending numerical order by Shell’s method (diminishing increment
//sort). a is replaced on output by its sorted rearrangement. b and c are sorted according to the order of a.

int GetQSortIndices(TArray<int>& arr, TArray<int>& indices);
int FindRedundantEntries_QSort(TArray<int>& arr, TArray<int>& flags_redundant,int flag_sortascending = 0);
int RemoveRedundantEntries(TArray<int>& arr, int flag_sortascending = 0); // removes redundant entries from a TArray<int>, optional ascending sort 
int FindRedundantEntries_QSort(TArray<double>& arr, TArray<int>& flags_redundant,int flag_sortascending = 0);
int RemoveRedundantEntries(TArray<double>& arr, int flag_sortascending = 0); // removes redundant entries from a TArray<double>, optional ascending sort 
int RemoveRedundantEntries_ErrorTolerance(TArray<double>& arr, double tol = 1e-10, int flag_sortascending = 0); // removes redundant entries (with error tolerance) from a TArray<double>, optional ascending sort 

void GetIntegrationRule(Vector& x, Vector& w, int order);
void GetIntegrationRuleLobatto(Vector& x, Vector& w, int order);
void GetIntegrationRuleTrig(Vector& x1, Vector& x2, Vector& w, int order); 	//int F(r,s) dr ds = sum(F(ri, si)*w(i))
void GetIntegrationRuleTet(Vector& x1, Vector& x2, Vector& x3, Vector& w, int order); 	//int F(r,s,t) dr ds dt = sum(F(ri, si, ti)*w(i))

void Quicksort(TArray<int>& x);
//Sorts an array x(1..x.Length()) into ascending numerical order by Shell’s method. x is replaced on output by its sorted rearrangement.

template <typename T>
void Quicksort(TArray<int>& a, TArray<T>& b)
{
	int n = a.Length();
	int i,j,inc;
	int v;
  T vb;
	inc=1; //Determine the starting increment.
	do 
	{
		inc *= 3;
		inc++;
	} while (inc <= n);
	do 
	{ //Loop over the partial sorts.
		inc /= 3;
		for (i=inc;i<n;i++) 
		{ //Outer loop of straight insertion.
			v=a.Get0(i);
			vb=b.Get0(i);
			j=i;
			while (a.Get0(j-inc) > v) 
			{ //Inner loop of straight insertion.
				a.Elem0(j)=a.Get0(j-inc);
				b.Elem0(j)=b.Get0(j-inc);
				j -= inc;
				if (j < inc) break;
			}
			a.Elem0(j)=v;
			b.Elem0(j)=vb;
		}
	} while (inc > 1);

}