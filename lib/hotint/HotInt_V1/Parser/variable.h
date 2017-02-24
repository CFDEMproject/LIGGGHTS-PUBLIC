
//#**************************************************************
//#
//# filename:             variable.h
//#
//# autor:                Gerstmayr Johannes
//#
//# generated:            22.11.98
//# last change:          22.11.98
//# description:          base class for variables
//# remarks:       variablelist: index from 1 to obj n
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
//**************************************************************


#ifndef VARIABLE__H
#define VARIABLE__H

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//CVariable
//the CVariable is a ParseObj, which contains:
// * a variable name
// * a pointer to a CMathObj which contains the variable information (either CMNumber or a symbolic math CMathObj structure)
// * the CVariable is the interface of the CMVariable to the stored CMathObj
//
// * during assignment, the variable is either directly evaluated to double, or the symbolic structure is pertained,
//   however, if the variable itself is contained, it is eliminated first 
//   ==>this can produce an error, e.g. if a=a and a is not initialized
//   the code a=1; a=a+b+1 works!
class CVariable : public CParseObj
{
public:
  CVariable();
  CVariable(const mystr& namei);

  virtual TParseObj ParseObjType() {return TPOVar;}
  virtual CParseObj* GetCopy() 
  {
    CVariable* o = new CVariable();
    o->mathObj = mathObj;
    o->name = name;
    return o;
  }
  virtual mystr Name();
	virtual void SetName(mystr namei);
  virtual void Assign(CMathObj* mathObji);

	//add: 
	virtual double GetDouble() {return MathObj()->Double();}
	virtual int EvaluableDouble() {return MathObj()->EvaluableDouble();}

//remove private
//private:
  virtual CMathObj* MathObj();
	virtual void SetMathObj(CMathObj* mathObj);

protected:
  CMathObj* mathObj;
  mystr name;
};

//the CVariableList contains a list of CParseObj, which are CVariables
class CVariableList : public CParseObjList
{
public:
  CVariableList();

  virtual CVariable* Object(mystr str);

  virtual CVariable* GetVariable(int i);
};


#endif
