
//#**************************************************************
//#
//# filename:             parseobjlist.cpp
//#
//# autor:                Gerstmayr Johannes
//#
//# generated:            22.11.98
//# last change:          22.11.98
//# description:          list of parseobjects
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

#include "mbs_interface.h"
#include "parser.h"
#include "parseobjlist.h"


CParseObjList::CParseObjList()
{
	InitConstructor();
}

//$ RL 2011-5-30:[ //$ JG 2011-5-31:[
CParseObjList& CParseObjList::operator=(const CParseObjList& e) 
{
	if (this == &e) {return *this;}
	CopyFrom(e);
	return *this;
}

CParseObjList::CParseObjList(const CParseObjList& e) 
{
	InitConstructor();
	CopyFrom(e);
}

void CParseObjList::CopyFrom(const CParseObjList& e) 
{
	//RL: no delete function --> parsed objects already existing
	assert(0); //RL: copyFrom not needed yet (variableList stored in mbs parser)
}

void CParseObjList::InitConstructor()
{
	parseObjList = 0;
  parseObjList = new AParseObjList();
}

int CParseObjList::IsObject(mystr str)
{
  for (int i = 1; i <= parseObjList->Length(); i++)
  {
    if((*parseObjList)[i]->Name() == str) {return i;}
  }
  return 0;
}

int CParseObjList::Length() const
{
  return parseObjList->Length();
}

CParseObj* CParseObjList::Get(int i)
{
  if (i > Length() || i <= 0) 
  {
    SysError("CParseObjList_Get");
    return (*parseObjList)[1];
  }
  return (*parseObjList)[i];
}

void CParseObjList::Add(CParseObj* obj)
{
  parseObjList->Add(obj);
}

void CParseObjList::Flush()
{
  for (int i = 1; i <= parseObjList->Length(); i++)
  {
    delete (*parseObjList)[i];
  }
  parseObjList->Flush();
}
