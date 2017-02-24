//#**************************************************************
//#
//# filename:             parseobj.h
//#
//# autor:                Gerstmayr Johannes
//#
//# generated:             1.11.98
//# last change:          20.12.98
//# description:          base class for parseobjects
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

#ifndef PARSEOBJ__H
#define PARSEOBJ__H

/*
GetMathObj
  function, bracket, unoperator
GenParseObjectList
klassenaufteilung

GenObjList von auﬂen erweiterbar um neue Befehle!!!
TCATest: only for 'else', is not necessary needed!
*/
//**************************************************************

typedef enum {TPOStatement=1, TPOMathObj=2, TPOVar=4, TPOOpenBr=8, TPOCloseBr=16,
      TPOComma=32, TPOSemicolon=64, TPOOpenState=128, TPOCloseState=256,
      TPOString=512, TPOVoid=1024, TPOTest=2048, TPOAssign=4096, TPOEndObj=8192,
      TPOExMathObj=2+4+4096, TPOCommand=16384, TPOElse=32768, 
      TPOOpenArgBr=0x10000, TPOCloseArgBr=0x20000 /*, TPOElementData=0x40000*/}
  TParseObj;

typedef enum {TMOMathObj = 1, TMONumber = 2, TMOVariable = 4, TMOFunction = 8,
      TMOOperator = 16, TMOUnOperator = 32, TMOMatrix = 64, TMOOpenBr = 128,
      TMOCloseBr = 256, TMOEnd = 512, TMOBase = 1024, TMOBracket = 2048,
      TMOVecFunction = 4096, TMOMatElem = 8192,
      TMOOperand = 2+4+8+32+64+128}
  TMathObj;

//typedef enum {TPCAArgs=1, TPCAMathObj=2, TPCANoArg=3}
//  TPCommandArg;

//this is the most general object in the parser which can be parsed, stored and/or computed
class CParseObj //base class
{
public:
  CParseObj() {};
  virtual mystr Name() {return mystr("@CParseObj");}
  virtual TParseObj ParseObjType() {return TPOVoid;}
  virtual CParseObj* GetCopy() {return new CParseObj();}
};

class CPString : public CParseObj
{
public:
  CPString() : CParseObj() {str = "";}
  CPString(const mystr& s) : CParseObj() {str = s;}
  virtual mystr Name() {return str;}
  virtual TParseObj ParseObjType() {return TPOString;}
  virtual CParseObj* GetCopy() {return new CPString();}

  virtual const mystr& String() {return str;}

private:
  mystr str;
};

//parse object, needed if "=" is expected
class CAssign : public CParseObj
{
public:
  virtual mystr Name() {return mystr("=");}
  virtual TParseObj ParseObjType() {return TPOAssign;}
  virtual CParseObj* GetCopy() {return new CAssign();}
};

//+++++++++++++++++++++++++++++++++++++++++++++
//the following objects are mostly for definition of a sequence of objects defining a command/function/etc.
//e.g. a function consists of "identifier" + "(" + CMathObj + "," + CMathObj + ")"
//     for example: f(x,y)
// in this way, the parser can have a list of expected objects
// some of these objects are used to tell the parser that a bracket is closed or that the end of the string has been reached



class CPVoid : public CParseObj
{
public:
  virtual mystr Name() {return mystr("@Void");}
  virtual TParseObj ParseObjType() {return TPOVoid;}
  virtual CParseObj* GetCopy() {return new CPVoid();}
};

class CEndObj : public CParseObj
{
public:
  virtual mystr Name() {return mystr("@EndObj");}
  virtual TParseObj ParseObjType() {return TPOEndObj;}
  virtual CParseObj* GetCopy() {return new CEndObj();}
};

class COpenState : public CParseObj
{
public:
  virtual mystr Name() {return mystr("{");}
  virtual TParseObj ParseObjType() {return TPOOpenState;}
  virtual CParseObj* GetCopy() {return new COpenState();}
};

class CCloseState : public CParseObj
{
public:
  virtual mystr Name() {return mystr("}");}
  virtual TParseObj ParseObjType() {return TPOCloseState;}
  virtual CParseObj* GetCopy() {return new CCloseState();}
};

class COpenArgBr : public CParseObj
{
public:
  virtual mystr Name() {return mystr("[");}
  virtual TParseObj ParseObjType() {return TPOOpenArgBr;}
  virtual CParseObj* GetCopy() {return new COpenArgBr();}
};

class CCloseArgBr : public CParseObj
{
public:
  virtual mystr Name() {return mystr("]");}
  virtual TParseObj ParseObjType() {return TPOCloseArgBr;}
  virtual CParseObj* GetCopy() {return new CCloseArgBr();}
};

class COpenBr : public CParseObj
{
public:
  virtual mystr Name() {return mystr("(");}
  virtual TParseObj ParseObjType() {return TPOOpenBr;}
  virtual CParseObj* GetCopy() {return new COpenBr();}
};

class CCloseBr : public CParseObj
{
public:
  virtual mystr Name() {return mystr(")");}
  virtual TParseObj ParseObjType() {return TPOCloseBr;}
  virtual CParseObj* GetCopy() {return new CCloseBr();}
};

class CComma : public CParseObj
{
public:
  virtual mystr Name() {return mystr(",");}
  virtual TParseObj ParseObjType() {return TPOComma;}
  virtual CParseObj* GetCopy() {return new CComma();}
};

class CSemicolon : public CParseObj
{
public:
  virtual mystr Name() {return mystr(";");}
  virtual TParseObj ParseObjType() {return TPOSemicolon;}
  virtual CParseObj* GetCopy() {return new CSemicolon();}
};

class CPElse : public CParseObj
{
public:
  virtual mystr Name() {return mystr("else");}
  virtual TParseObj ParseObjType() {return TPOElse;}
  virtual CParseObj* GetCopy() {return new CPElse();}
};


#endif
