//#**************************************************************
//# filename:             DataStorage.h
//#
//# author:               Gerstmayr, Vetyukov
//#
//# generated:						
//# description:					interface for the DataStorage class.     
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
 
#if !defined(AFX_DATASTORAGE_H__399321A2_93FB_4A66_AC5B_B62B95851816__INCLUDED_)
#define AFX_DATASTORAGE_H__399321A2_93FB_4A66_AC5B_B62B95851816__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "..\WorkingModule\WinCompDriverInterface.h"

// here we implement a data storage unit for a single time point
// for the usage conventions see the interface file included above

class DataStorage : 
	public WCDInterface::DataLoader, 
	public WCDInterface::DataSaver  
{
	double PointTime;	// time value to which this data corresponds
	char * pData;
	int nSize;
	int nReservedMemorySize;
	char * pCurrentReadPosition;

	// this function allocates new memory if needed
	// returns the position to write the new data
	char * NewDataWritePos(int nAdditionalSize);

	// maximal size between all the DataStorage objectes
	static int nMaxSize;
	static long nMemoryUsed;
	
	// implementation for the interface DataSaver
	virtual void SetTime(double time_);
	virtual WCDInterface::DataSaver & operator << (int);
	virtual WCDInterface::DataSaver & operator << (double);
	virtual WCDInterface::DataSaver & operator << (const char *);

	// implementation for the interface DataLoader
	virtual double GetTime();
	virtual WCDInterface::DataLoader & operator >> (int &);
	virtual WCDInterface::DataLoader & operator >> (double &);
	virtual WCDInterface::DataLoader & operator >> (char * &);

public:
	DataStorage();
	DataStorage(const DataStorage & ds) { *this = ds; }
	virtual ~DataStorage();

	static long GetMemoryUsed() { return nMemoryUsed; }

	virtual void AllocateNewData(int sizebytes) 
	{
		delete[] pData; 
		nReservedMemorySize = sizebytes;
		pData = new char[sizebytes];
	}

	DataStorage & operator=(const DataStorage & ds);

	void Serialize(CArchive & ar);
};

#endif // !defined(AFX_DATASTORAGE_H__399321A2_93FB_4A66_AC5B_B62B95851816__INCLUDED_)
