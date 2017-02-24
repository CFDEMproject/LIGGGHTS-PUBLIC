//#**************************************************************
//# filename:             DataStorage.cpp
//#
//# author:               Gerstmayr, Vetyukov
//#
//# generated:						
//# description:          
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
 
#include "stdafx.h"
#include "WCDriver3DDlg.h"
#include "WCDriver3D.h"
#include "DataStorage.h"
#include "memory.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

// from the very beginning the data storage units
//  allocate this quantity of bytes
int DataStorage::nMaxSize = 50;

// counter of the memory usage
long DataStorage::nMemoryUsed = 0;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


DataStorage::DataStorage() :
PointTime(-1),
nSize(0),
nReservedMemorySize(nMaxSize),
pCurrentReadPosition(NULL)
{
	pData = new char[nMaxSize];
	nMemoryUsed += nMaxSize;
}

DataStorage::~DataStorage()
{
	delete[] pData;
	nMemoryUsed -= nReservedMemorySize;
}

DataStorage & DataStorage::operator=(const DataStorage & ds)
{
	PointTime = ds.PointTime;
	nSize = ds.nSize;
	nReservedMemorySize = ds.nReservedMemorySize;
	delete[] pData;
	pData = new char[nReservedMemorySize];
	memcpy(pData,ds.pData,nSize);
	pCurrentReadPosition = NULL;

	return *this;
}


// implementation

char * DataStorage::NewDataWritePos(int nAdditionalSize)
{
	if(nSize + nAdditionalSize > nReservedMemorySize)
	{
		nMemoryUsed += nReservedMemorySize;
		char * pNewData =
			new char[nReservedMemorySize *= 2];
		memcpy(pNewData,pData,nSize);
		delete[] pData;
		pData = pNewData;
	}
	char * pWritePos = pData + nSize;
	nSize += nAdditionalSize;
	nMaxSize = max(nMaxSize,nSize);

	pCurrentReadPosition = pData;

	return pWritePos;
}

void DataStorage::SetTime(double time_)
{
	PointTime = time_;
}

WCDInterface::DataSaver & DataStorage::operator << (int i)
{
	*((int*)NewDataWritePos(sizeof(int))) = i;
	return *this;
}

WCDInterface::DataSaver & DataStorage::operator << (double d)
{
	*((double*)NewDataWritePos(sizeof(double))) = d;
	return *this;
}

WCDInterface::DataSaver & DataStorage::operator << (const char * p)
{
	int nBytes = strlen(p) + 1;		// 1 byte for the terminating zero
	memcpy(NewDataWritePos(nBytes),p,nBytes);
	return *this;
}

double DataStorage::GetTime()
{
	pCurrentReadPosition = pData;
	return PointTime;
}

WCDInterface::DataLoader & DataStorage::operator >> (int & i)
{
	ASSERT(pCurrentReadPosition);
	i = *((int*)pCurrentReadPosition);
	pCurrentReadPosition += sizeof(int);
	ASSERT(pCurrentReadPosition - pData <= nSize);
	return *this;
}

WCDInterface::DataLoader & DataStorage::operator >> (double & d)
{
	ASSERT(pCurrentReadPosition);
	d = *((double*)pCurrentReadPosition);
	pCurrentReadPosition += sizeof(double);
	ASSERT(pCurrentReadPosition - pData <= nSize);
	return *this;
}

WCDInterface::DataLoader & DataStorage::operator >> (char * & p)
{
	ASSERT(pCurrentReadPosition);
	const char * pStr = (char*)pCurrentReadPosition;
	int nBytes = strlen(pStr) + 1;		// 1 byte for the terminating zero
	pCurrentReadPosition += nBytes;
	p = new char[nBytes];
	memcpy(p,pStr,nBytes);
	ASSERT(pCurrentReadPosition - pData <= nSize);
	return *this;
}

void DataStorage::Serialize(CArchive & ar)
{
	if(ar.IsStoring())
	{
		ar << PointTime;
		ar << nSize;
		ar.Write(pData,nSize);
	}
	else
	{
		ar >> PointTime;
		delete[] pData;
		ar >> nSize;
		nReservedMemorySize = nSize;
		pData = new char[nSize];
		ar.Read(pData,nSize);
		nMaxSize = max(nMaxSize,nSize);
	}
}