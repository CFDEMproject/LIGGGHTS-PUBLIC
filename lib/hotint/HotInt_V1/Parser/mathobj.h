
//#**************************************************************
//#
//# filename:             mathobj.h
//#
//# autor:                Gerstmayr Johannes
//#
//# generated:            23.11.98
//# last change:          20.12.98
//# description:          base class for mathematical objects
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

#ifndef MATHOBJ__H
#define MATHOBJ__H

typedef enum {TMPLow, TMPBase, TMPBracket, TMPAssign, TMPCondOr,
              TMPCondAnd, TMPEqual, TMPOr, TMPAnd, 
              TMPPlus, TMPMul, TMPPow, TMPFunction, TMPHigh}
  TMathPriority;

class CVariable;
  
double fact(double n);

typedef TArray<CVariable*> ACVariable;

//CMathObj: a mathematical object, which can be parsed
//in contrast to other CParseObj, the CMathObj can be computed (double, Vector/Matrix, string)
class CMathObj : public CParseObj
{
public:
  CMathObj();
    //name of object eg. '+' or 'sin'
  virtual mystr Name();
    //get copy of this object (parseobj-function)
  virtual CParseObj* GetCopy();
    //get copy of this object and of all arguments
  virtual CMathObj* CopyMO();
    //delete this object and all arguments
  virtual void DeleteMO();
    //type of mathobj, eg. TMOOperator
  virtual TMathObj MOType();
  //if mathobj hass no base obj, then add one
  virtual CMathObj* AddBase();

    //is always TPOMathObj
  virtual TParseObj ParseObjType();
  //get Argument x
  virtual CMathObj* Arg(int arg = 0);
  //get number of arguments, -1 == undefined
  virtual int GetNOArg();
  //set argument x
  virtual void SetArg(CMathObj* arg, int argnum = 0);
  //while parsing, add an object at the next position (->operators at arg1!)
  virtual void AddArg(CMathObj* arg);

	//get parent object of this object:
  virtual CMathObj* Parent();
	//set parent object of this object:
  virtual void SetParent(CMathObj* par);

	//set in value of variable var, for assign
  virtual void Eliminate(CVariable* var);
  //replace oldarg with replacearg, set parent in replacearg
  virtual void ReplaceArg(CMathObj* oldarg, CMathObj* replacearg);
  //Evaluate expression to a number, error, if not possible
  virtual double Double();
  //check, if evaluation to double is possible
	//evaluable returns 1, if : 1) all objects of operator are evaluable to double 2) operand is evaluable to double (must be a number) 3) variable is evaluable to double
  virtual int EvaluableDouble();

    //get mystr of mathobj
  virtual mystr String() {return Name();};
    //get mystr of mathobj with expressed variables ==> with this function, the MathObj can be printed!
  virtual mystr StringAll() {return String();};

    //return priority of obj in sense of operators
    //it defines the affinity to other objects
  virtual TMathPriority Priority() {return TMPLow;}
    //executes generation of matrix etc., e.g. when assigning

  virtual CMathObj* Differentiate(CMathObj* var);
  //StdSimplify is for simplification of eg. 1*x=x, 2+2=4, ...
  virtual CMathObj* StdSimplify();
  
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //special matrix funcitons:
  //virtual void Elem(int row, int col, CMathObj* replace); nur für variable
	//get CMathObj at specific row and column
  virtual CMathObj* Get(int row, int col); //for mateval
	//get double entry of matrix at specific row and column
  virtual double MatDouble(int row, int col);
	//get dimensionality (rc='row' and 'column') of MathObj: rc=[1,1] ==> double value
  virtual void Dim(int2& rc);

private:
  CMathObj* parent;
};

//typedef CMathObj* CMsthObjP;
typedef TArray<CMathObj*> ACMathObj;
//+++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++   OPERATORS     ++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++
//an operator links two CMathObjs
//e.g. Mathematical operations '+', '*', '^', '-'
//Compare operators: '==', '<', '>=', "!=",  etc.
//conditional operators: '&&', '||'
//binary operators: '&', '|'
class COperator : public CMathObj
{
public:
  COperator();

  virtual mystr Name();
  virtual CParseObj* GetCopy() {return new COperator();}

  virtual void SetArg(CMathObj* arg, int argnum = 0);

  virtual CMathObj* Arg(int arg = 0);

  //virtual double MatDouble(int row, int col);

  virtual double Double();

  virtual mystr String();

  virtual mystr StringAll();

  virtual void AddArg(CMathObj* arg);

  virtual int GetNOArg();

  virtual TMathObj MOType();

  virtual TMathPriority Priority();

  virtual int Commutative() {return 1;}

  virtual CMathObj* StdSimplify();

protected:
  CMathObj* arg0; //left
  CMathObj* arg1; //right
};

class CMEqual : public COperator
{
public:
  CMEqual() : COperator() {};
  virtual CParseObj* GetCopy() {return new CMEqual();}
  virtual mystr Name() {return mystr("==");}
  virtual double Double() {return Arg(0)->Double() == Arg(1)->Double();}
  virtual TMathObj MOType() {return TMOOperator;}
  virtual TMathPriority Priority() {return TMPEqual;}
};

class CMUnEqual : public COperator
{
public:
  CMUnEqual() : COperator() {};
  virtual CParseObj* GetCopy() {return new CMUnEqual();}
  virtual mystr Name() {return mystr("!=");}
  virtual double Double() {return Arg(0)->Double() != Arg(1)->Double();}
  virtual TMathObj MOType() {return TMOOperator;}
  virtual TMathPriority Priority() {return TMPEqual;}
};

class CMGrTh : public COperator
{
public:
  CMGrTh() : COperator() {};
  virtual CParseObj* GetCopy() {return new CMGrTh();}
  virtual mystr Name() {return mystr(">");}
  virtual double Double() {return Arg(0)->Double() > Arg(1)->Double();}
  virtual TMathObj MOType() {return TMOOperator;}
  virtual TMathPriority Priority() {return TMPEqual;}
  virtual int Commutative() {return 0;}
};

class CMLeTh : public COperator
{
public:
  CMLeTh() : COperator() {};
  virtual CParseObj* GetCopy() {return new CMLeTh();}
  virtual mystr Name() {return mystr("<");}
  virtual double Double() {return Arg(0)->Double() < Arg(1)->Double();}
  virtual TMathObj MOType() {return TMOOperator;}
  virtual TMathPriority Priority() {return TMPEqual;}
  virtual int Commutative() {return 0;}
};

class CMGrEq : public COperator
{
public:
  CMGrEq() : COperator() {};
  virtual CParseObj* GetCopy() {return new CMGrEq();}
  virtual mystr Name() {return mystr(">=");}
  virtual double Double() {return Arg(0)->Double() >= Arg(1)->Double();}
  virtual TMathObj MOType() {return TMOOperator;}
  virtual TMathPriority Priority() {return TMPEqual;}
  virtual int Commutative() {return 0;}
};

class CMLeEq : public COperator
{
public:
  CMLeEq() : COperator() {};
  virtual CParseObj* GetCopy() {return new CMLeEq();}
  virtual mystr Name() {return mystr("<=");}
  virtual double Double() {return Arg(0)->Double() <= Arg(1)->Double();}
  virtual TMathObj MOType() {return TMOOperator;}
  virtual TMathPriority Priority() {return TMPEqual;}
  virtual int Commutative() {return 0;}
};

class CMCondAnd : public COperator
{
public:
  CMCondAnd() : COperator() {};
  virtual CParseObj* GetCopy() {return new CMCondAnd();}
  virtual mystr Name() {return mystr("&&");}
  virtual double Double() {return Arg(0)->Double() && Arg(1)->Double();}
  virtual TMathObj MOType() {return TMOOperator;}
  virtual TMathPriority Priority() {return TMPCondAnd;}
};

class CMCondOr : public COperator
{
public:
  CMCondOr() : COperator() {};
  virtual CParseObj* GetCopy() {return new CMCondOr();}
  virtual mystr Name() {return mystr("||");}
  virtual double Double() {return Arg(0)->Double() || Arg(1)->Double();}
  virtual TMathObj MOType() {return TMOOperator;}
  virtual TMathPriority Priority() {return TMPCondOr;}
};

class CMAnd : public COperator
{
public:
  CMAnd() : COperator() {};
  virtual CParseObj* GetCopy() {return new CMAnd();}
  virtual mystr Name() {return mystr("&");}
  virtual double Double() {return (int)Arg(0)->Double() & (int)Arg(1)->Double();}
  virtual TMathObj MOType() {return TMOOperator;}
  virtual TMathPriority Priority() {return TMPAnd;}
};

class CMOr : public COperator
{
public:
  CMOr() : COperator() {};
  virtual CParseObj* GetCopy() {return new CMOr();}
  virtual mystr Name() {return mystr("|");}
  virtual double Double() {return (int)Arg(0)->Double() | (int)Arg(1)->Double();}
  virtual TMathObj MOType() {return TMOOperator;}
  virtual TMathPriority Priority() {return TMPOr;}
};

class CMPlus : public COperator
{
public:
  CMPlus() : COperator() {};
  virtual CParseObj* GetCopy() {return new CMPlus();}

  virtual mystr Name() {return mystr("+");}

  virtual double Double()
  {
    return Arg(0)->Double() + Arg(1)->Double();
  }

  virtual CMathObj* Differentiate(CMathObj* var);
  virtual CMathObj* StdSimplify();

  virtual double MatDouble(int row, int col)
  {
    return Arg(0)->MatDouble(row,col) + Arg(1)->MatDouble(row,col);
  }
  virtual void Dim(int2& rc);

  virtual TMathObj MOType() {return TMOOperator;}
  virtual TMathPriority Priority() {return TMPPlus;}
};

class CMMinus : public COperator
{
public:
  CMMinus() : COperator() {};

  virtual CParseObj* GetCopy() {return new CMMinus();}

  virtual mystr Name() {return mystr("-");}

  virtual double Double()
  {
    return Arg(0)->Double() - Arg(1)->Double();
  }

  virtual double MatDouble(int row, int col)
  {
    return Arg(0)->MatDouble(row,col) - Arg(1)->MatDouble(row,col);
  }

  virtual CMathObj* Differentiate(CMathObj* var);
  virtual CMathObj* StdSimplify();

  virtual void Dim(int2& rc);

  virtual TMathObj MOType() {return TMOOperator;}
  virtual TMathPriority Priority() {return TMPPlus;}

  //virtual mystr String();
  virtual int Commutative() {return 0;}
};

class CMMul : public COperator
{
public:
  CMMul() : COperator() {};

  virtual CParseObj* GetCopy() {return new CMMul();}

  virtual mystr Name() {return mystr("*");}

  virtual double Double()
  {
    return Arg(0)->Double() * Arg(1)->Double();
  }

  virtual CMathObj* Differentiate(CMathObj* var);
  virtual CMathObj* StdSimplify();

  virtual double MatDouble(int row, int col);

  virtual void Dim(int2& rc);

  virtual TMathObj MOType() {return TMOOperator;}
  virtual TMathPriority Priority() {return TMPMul;}
};

class CMDiv : public COperator
{
public:
  CMDiv() : COperator() {};

  virtual CParseObj* GetCopy() {return new CMDiv();}

  virtual mystr Name() {return mystr("/");}

  virtual double Double()
  {
    return Arg(0)->Double() / Arg(1)->Double();
  }

  virtual TMathObj MOType() {return TMOOperator;}
  virtual TMathPriority Priority() {return TMPMul;}

  virtual CMathObj* Differentiate(CMathObj* var);
  virtual CMathObj* StdSimplify();

  //virtual mystr String();
  virtual int Commutative() {return 0;}
};

class CMPow : public COperator
{
public:
  CMPow() : COperator() {};

  virtual CParseObj* GetCopy() {return new CMPow();}

  virtual mystr Name() {return mystr("^");}

  virtual double Double();

  virtual TMathObj MOType() {return TMOOperator;}
  virtual TMathPriority Priority() {return TMPPow;}

  virtual CMathObj* Differentiate(CMathObj* var);
  virtual CMathObj* StdSimplify();

  virtual int Commutative() {return 0;}
};

//$JG2013-5-3: removed modulo because conflict to comment sign '%' ==> should be replaced by command mod(x,y) (replaces x % y)
//class CMMod : public COperator
//{
//public:
//  CMMod() : COperator() {};
//
//  virtual CParseObj* GetCopy() {return new CMMod();}
//
//  virtual mystr Name() {return mystr("%");}
//
//  virtual double Double();
//
//  virtual TMathObj MOType() {return TMOOperator;}
//  virtual TMathPriority Priority() {return TMPMul;}
//  virtual int Commutative() {return 0;}
//};

//+++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++   OPERANDS      ++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++

class COperand : public CMathObj
{
public:
  COperand() : CMathObj() {};
  virtual CParseObj* GetCopy() {return new COperand();}
  virtual mystr Name() {return mystr("@Operand");}
  virtual TMathObj MOType() {return TMOOperand;}
};

//the CMNumber is the container for numeric values and usually stays at the end of a symbolic object tree
class CMNumber : public COperand
{
public:
	CMNumber(): COperand() {n = 0;}
  CMNumber(double ni) : COperand() {n = ni;}
  virtual CParseObj* GetCopy() {return new CMNumber(n);}
  virtual mystr Name();
  virtual TMathObj MOType() {return TMONumber;}

  virtual double Double() {return n;}
  virtual int EvaluableDouble() {return 1;}

	virtual CMathObj* Differentiate(CMathObj* var);
  virtual CMathObj* StdSimplify();

  virtual void Set(double ni) {n = ni;}

private:
  double n;
};

//the MathObj CMVariable contains a pointer to a CVariable, which is stored in a CVariableList
//the name of the MathObj is identical to the CVariable
//this object is necessary, because a CVariable cannot be computed (CVariable is the 'container' for the variable)
class CMVariable : public COperand
{
public:
  CMVariable();

  CMVariable(CVariable* vari);
  
  virtual CParseObj* GetCopy();
  
  virtual mystr Name();
  
  virtual TMathObj MOType();

  virtual double Double();
  
  virtual int EvaluableDouble();

  virtual void Eliminate(CVariable* vari);
  
  virtual CVariable* Variable();

  virtual mystr StringAll();

  virtual CMathObj* Differentiate(CMathObj* var);
  virtual CMathObj* StdSimplify();

  virtual void Elem(int row, int col, CMathObj* replace);

  virtual CMathObj* Get(int row, int col);
  virtual double MatDouble(int row, int col);
  virtual void Dim(int2& rc);

private:
  CVariable* var;
};


//+++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++   FUNCTIONS     ++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++

class CFunction : public CMathObj
{
public:
  CFunction();

  virtual CParseObj* GetCopy() 
  {
    CFunction* o = new CFunction();
    o->arg0 = arg0;
    return o;
  }

  virtual mystr Name();

  virtual CMathObj* Arg(int arg = 0);

  virtual void SetArg(CMathObj* arg, int argnum = 0);

  virtual int GetNOArg();

  virtual TMathObj MOType();

  virtual double Double();

  virtual mystr String();

  virtual mystr StringAll();

  virtual TMathPriority Priority() {return TMPFunction;}

  virtual CMathObj* StdSimplify();

protected:
  CMathObj* arg0;
};

//is elimated during parsing
class CMUnPlus : public CFunction
{
public:
  CMUnPlus() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMUnPlus();}
  virtual mystr Name() {return mystr("+");}
  virtual double Double() {return Arg(0)->Double();}
  virtual TMathPriority Priority() {return TMPPlus;}
  virtual mystr String();
  virtual mystr StringAll();
  virtual TMathObj MOType() {return TMOUnOperator;}
  virtual void Dim(int2& rc);
  virtual double MatDouble(int row, int col)
  {
    return Arg(0)->MatDouble(row,col);
  }
};

class CMUnMinus : public CFunction
{
public:
  CMUnMinus() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMUnMinus();}
  virtual mystr Name() {return mystr("-");}
  virtual double Double() {return -Arg(0)->Double();}
  virtual TMathPriority Priority() {return TMPPlus;}
  virtual mystr String();
  virtual mystr StringAll();
  virtual TMathObj MOType() {return TMOUnOperator;}
  virtual CMathObj* Differentiate(CMathObj* var);
  virtual void Dim(int2& rc);
  virtual double MatDouble(int row, int col)
  {
    return -Arg(0)->MatDouble(row,col);
  }
};

class CMSqrt : public CFunction
{
public:
  CMSqrt() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMSqrt();}
  virtual mystr Name() {return mystr("sqrt");}
  virtual CMathObj* Differentiate(CMathObj* var);
  virtual double Double();
};

class CMSqr : public CFunction
{
public:
  CMSqr() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMSqr();}
  virtual mystr Name() {return mystr("sqr");}
  virtual CMathObj* Differentiate(CMathObj* var);
  virtual double Double();
};

class CMSin : public CFunction
{
public:
  CMSin() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMSin();}
  virtual mystr Name() {return mystr("sin");}
  virtual CMathObj* Differentiate(CMathObj* var);
  virtual double Double();
};

class CMCos : public CFunction
{
public:
  CMCos() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMCos();}
  virtual mystr Name() {return mystr("cos");}
  virtual CMathObj* Differentiate(CMathObj* var);
  virtual double Double();
};

class CMTan : public CFunction
{
public:
  CMTan() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMTan();}
  virtual mystr Name() {return mystr("tan");}
  virtual CMathObj* Differentiate(CMathObj* var);
  virtual double Double();
};

class CMASin : public CFunction
{
public:
  CMASin() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMASin();}
  virtual mystr Name() {return mystr("asin");}
  virtual double Double();
};

class CMACos : public CFunction
{
public:
  CMACos() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMACos();}
  virtual mystr Name() {return mystr("acos");}
  virtual double Double();
};

class CMATan : public CFunction
{
public:
  CMATan() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMATan();}
  virtual mystr Name() {return mystr("atan");}
  virtual double Double();
};

class CMSinh : public CFunction
{
public:
  CMSinh() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMSinh();}
  virtual mystr Name() {return mystr("sinh");}
  virtual double Double();
};

class CMCosh : public CFunction
{
public:
  CMCosh() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMCosh();}
  virtual mystr Name() {return mystr("cosh");}
  virtual double Double();
};

class CMTanh : public CFunction
{
public:
  CMTanh() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMTanh();}
  virtual mystr Name() {return mystr("tanh");}
  virtual double Double();
};

class CMExp : public CFunction
{
public:
  CMExp() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMExp();}
  virtual mystr Name() {return mystr("exp");}
  virtual CMathObj* Differentiate(CMathObj* var);
  virtual double Double();
};

  //natural logarithm  ln(e) == 1
class CMLn : public CFunction
{
public:
  CMLn() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMLn();}
  virtual mystr Name() {return mystr("ln");}
  virtual CMathObj* Differentiate(CMathObj* var);
  virtual double Double();
};

 //natural logarithm  log(e) == 1//$ RL 2011-6-10:[ base of log is e, log10 with base 10 added.
class CMLog : public CFunction
{
public:
  CMLog() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMLog();}
  virtual mystr Name() {return mystr("log");}
  virtual CMathObj* Differentiate(CMathObj* var);
  virtual double Double();
};

  //log10: log(10) == 1
class CMLog10 : public CFunction
{
public:
  CMLog10() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMLog10();}
  virtual mystr Name() {return mystr("log10");}
  virtual CMathObj* Differentiate(CMathObj* var);
  virtual double Double();
};
//$ RL 2011-6-10:] base of log is e, log10 with base 10 added.

class CMFact : public CFunction
{
public:
  CMFact() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMFact();}
  virtual mystr Name() {return mystr("fact");}
  virtual double Double();
};

class CMAbs : public CFunction
{
public:
  CMAbs() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMAbs();}
  virtual mystr Name() {return mystr("abs");}
  virtual double Double();
};

//$ RL 2011-6-10: [ fabs added (same as abs)
class CMFAbs : public CFunction 
{
public:
  CMFAbs() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMFAbs();}
  virtual mystr Name() {return mystr("fabs");}
  virtual double Double();
};
//$ RL 2011-6-10: ] fabs added (same as abs)

class CMRound : public CFunction
{
public:
  CMRound() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMRound();}
  virtual mystr Name() {return mystr("round");}
  virtual double Double();
};

class CMFloor : public CFunction
{
public:
  CMFloor() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMFloor();}
  virtual mystr Name() {return mystr("floor");}
  virtual double Double();
};

class CMCeil : public CFunction
{
public:
  CMCeil() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMCeil();}
  virtual mystr Name() {return mystr("ceil");}
  virtual double Double();
};

class CMHeaviside : public CFunction
{
public:
  CMHeaviside() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMHeaviside();}
  virtual mystr Name() {return mystr("heaviside");}
  virtual double Double();
};

class CMSgn : public CFunction
{
public:
  CMSgn() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMSgn();}
  virtual mystr Name() {return mystr("sgn");}
  virtual double Double();
};

  //get digits after Comma
class CMPPD : public CFunction
{
public:
  CMPPD() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMPPD();}
  virtual mystr Name() {return mystr("ppd");}
  virtual double Double();
};

/*
class CM$ : public CFunction
{
public:
  CM$() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CM$();}
  virtual mystr Name() {return mystr("$");}
  virtual double Double();
};
*/
class CMNot : public CFunction
{
public:
  CMNot() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMNot();}
  virtual mystr Name() {return mystr("!");}
  virtual double Double() {return !Arg(0)->Double();};
};

class CMTranspose : public CFunction
{
public:
  CMTranspose() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMTranspose();}
  virtual mystr Name() {return mystr("transpose");}
//  virtual double Double();

  virtual CMathObj* Get(int row, int col);
  virtual void Dim(int2& rc);
  virtual double MatDouble(int row, int col);  
};

//$ AD 2013-10-20: vector length
class CMVAbs : public CFunction
{
public:
	CMVAbs() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMVAbs();}
  virtual mystr Name() {return mystr("vabs");}
  virtual double Double();
};

//$ AD 2013-10-20: rows of a matrix
class CMMatRows : public CFunction
{
public:
	CMMatRows() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMMatRows();}
  virtual mystr Name() {return mystr("rows");}
  virtual double Double();
};

//$ AD 2013-10-20: columns of a matrix
class CMMatCols : public CFunction
{
public:
	CMMatCols() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMMatCols();}
  virtual mystr Name() {return mystr("cols");}
  virtual double Double();
};


//$ AD 2013-10-20: random number
class CMRandom : public CFunction
{
public:
	CMRandom() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMRandom();}
  virtual mystr Name() {return mystr("random");}
  virtual double Double();
};

class CMRandomSeed: public CFunction
{
public:
	CMRandomSeed() : CFunction() {};
  virtual CParseObj* GetCopy() {return new CMRandomSeed();}
  virtual mystr Name() {return mystr("random_seed");}
  virtual double Double();
};

class CBaseObj : public CFunction
{
public:
  CBaseObj();
  virtual mystr Name();
  virtual CParseObj* GetCopy() {return new CBaseObj();}
  virtual TMathObj MOType();
  virtual TMathPriority Priority() {return TMPBase;}
};

//+++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++      ELSE       ++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++

class CMEnd : public CMathObj
{
public:
  CMEnd();
  virtual CParseObj* GetCopy() {return new CMEnd();}
  virtual mystr Name();
  virtual TMathObj MOType();
};

class CMOpenBr : public CMathObj
{
public:
  CMOpenBr() : CMathObj() {};
  virtual mystr Name() {return mystr("(");}
  virtual CParseObj* GetCopy() {return new CMOpenBr();}
  virtual TMathObj MOType() {return TMOOpenBr;};
};

class CMCloseBr : public CMathObj
{
public:
  CMCloseBr() : CMathObj() {};
  virtual mystr Name() {return mystr(")");}
  virtual CParseObj* GetCopy() {return new CMCloseBr();}
  virtual TMathObj MOType() {return TMOCloseBr;};
};

class CMBracket : public CFunction
{
public:
  CMBracket() : CFunction() {};
  virtual mystr Name() {return mystr("()");};
  virtual CParseObj* GetCopy() {return new CMBracket();}
  virtual TMathObj MOType() {return TMOBracket;};
  virtual TMathPriority Priority() {return TMPBracket;}

  virtual mystr String();
  virtual mystr StringAll();
};

//matrix is read while parsing. if no elements are 
//defined (e.g. x=matrix(a,b) ) the size of the matrix is defined
//by two mathobjs. while executing, the matrix x is generated
//when first accessing to an element; when generated once, generated = 1,
//and the size is unchangeable! 
//if the matrix is generated through e.g. matrix(2,2, 1,0;1,0)
//then the matrix is generated instantly and generated = 1 before 
//executing starts
class CMatrix : public CMathObj
{
public:
  CMatrix();

  CMatrix(CMathObj* cmorowsi, CMathObj* cmocolsi);

  virtual CParseObj* GetCopy(); 

  virtual mystr Name() {return mystr("matrix");}

  virtual CMathObj* Arg(int arg = 0);

  virtual void SetArg(CMathObj* arg, int argnum = 0);

  virtual int GetNOArg();

  virtual TMathObj MOType() {return TMOMatrix;};

  virtual double Double();

  virtual mystr String();

  virtual mystr StringAll();

  virtual TMathPriority Priority() {return TMPFunction;}

    //important, for assigning not only the values of a matrix
  virtual int EvaluableDouble() {return 0;}
  //virtual void MatExecute();

//special matrix funcitons:
  virtual void Dim(int2& rc);
  virtual double MatDouble(int row, int col);  

//for symbolic computations
  virtual CMathObj* Get(int row, int col);
//only for assigning a variable
  virtual void Elem(int row, int col, CMathObj* replace);

    //generate matrix dimensions, if not yet generated
  virtual void Generate();
  virtual void SetDimensions(CMathObj* mor, CMathObj* moc);
  virtual void AddElem(CMathObj* cm);
  
private:

  CMathObj* cmorows;
  CMathObj* cmocols;
  int generated;
  int rows;
  int cols;
  
  ACMathObj* mathobjs;
};

  //function with an arbitrary number of arguments;
class CVecFunction : public CMathObj
{
public:
  CVecFunction();

  virtual CParseObj* GetCopy();

  virtual mystr Name() {return mystr("@VecFunction");};

  virtual CMathObj* Arg(int arg = 0);

  virtual void SetArg(CMathObj* arg, int argnum = 0);

  virtual int GetNOArg() {return 0;};

  virtual TMathObj MOType() {return TMOFunction;};

  virtual double Double() {return 0;};

  virtual mystr StringAll();

  virtual mystr String();

  virtual TMathPriority Priority() {return TMPFunction;}

private:
  ACMathObj* args;
};


// matrix element access function
//arg(1) and arg(2) are matrix rows and cols, arg(0) is the
//mathobj of the matrix
class CMMatElem : public CVecFunction
{
public:
  CMMatElem() : CVecFunction() {};
  virtual mystr Name() {return mystr("[]");};
  virtual CParseObj* GetCopy();
  virtual TMathObj MOType() {return TMOMatElem;};
  virtual int GetNOArg() {return 3;};

  virtual void SetVar(CMathObj* cm);
  virtual double Double();
  virtual mystr String();
  virtual mystr StringAll();
};


#endif

