
//#**************************************************************
//#
//# filename:             parseobjlist.h
//#
//# autor:                Gerstmayr Johannes
//#
//# generated:            22.11.98
//# last change:          20.12.98
//# description:          list of parseobjects
//# remarks:              index from 1 to obj n
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


#ifndef PARSEOBJLIST__H
#define PARSEOBJLIST__H


typedef TArray<CParseObj*> AParseObjList;

//the CParseObjList contains a list of objects which have been parsed
//this can be 
//  - a list of commands, 
//  - a list of parameters or 
//  - a list of variables ==> see CVariableList
class CParseObjList
{
public:
  CParseObjList();
	
  virtual ~CParseObjList()
  {
	  if (this->parseObjList != 0)
	  {
		  delete parseObjList;
	  }
  }

	CParseObjList& operator=(const CParseObjList& e);

	CParseObjList(const CParseObjList& e);

	virtual void CopyFrom(const CParseObjList& e);

	virtual void InitConstructor();

  virtual int IsObject(mystr str); //check if string 'str' occurs in one of the objects of parseObjList

  virtual int Length() const; 

  virtual CParseObj* Get(int i); //get i-th entry of parseObjList

  virtual void Add(CParseObj* obj);

  virtual void Flush(); //clear list

private:
 AParseObjList* parseObjList;
};



#endif
