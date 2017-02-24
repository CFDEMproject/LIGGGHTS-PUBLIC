
//#**************************************************************
//#
//# filename:             variable.cpp
//#
//# autor:                Gerstmayr Johannes
//#
//# generated:            29.11.98
//# last change:          29.11.98
//# description:          implementation class for variables
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


//+++++++++++++++++++++  CVariable +++++++++++++++++++++

CVariable::CVariable() 
{
  mathObj = NULL;
};

CVariable::CVariable(const mystr& namei)
{
  mathObj = NULL;
  name = namei;
}

mystr CVariable::Name() 
{
  return name;
}

void CVariable::SetName(mystr namei)
{
	name = namei;
}

void CVariable::Assign(CMathObj* mathObji)
{
  if (mathObj != NULL) 
  {
    //mathObj->DeleteMO(); 
    //delete mathObj;
  }
  mathObj = mathObji;
}

CMathObj* CVariable::MathObj() 
{
  return mathObj;
}

void CVariable::SetMathObj(CMathObj* mathObj) 
{
  this->mathObj = mathObj;
}

//+++++++++++++++++++++  CVariableList +++++++++++++++++++++

CVariableList::CVariableList() : CParseObjList() 
{
}

CVariable* CVariableList::Object(mystr str)
{
  int i = IsObject(str);
  if (i)
  {
    return (CVariable*)Get(i);
  } else
  {
    CVariable* var = new CVariable(str);
    Add(var);
    return var;
  }
}

CVariable* CVariableList::GetVariable(int i)
{
	return (CVariable*)Get(i);
}


