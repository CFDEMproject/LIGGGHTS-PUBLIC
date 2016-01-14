
//#**************************************************************
//#
//# filename:             command.cpp
//#
//# autor:                Gerstmayr Johannes
//#
//# generated:            23.11.98
//# last change:          20.12.98
//# description:          implementation for commands
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
//**************************************************************

#include "mbs_interface.h"
#include "parser.h"

extern int ParserErrorFlag();


//++++++++++++++++++++ CCASSIGN ++++++++++++++++++++++

CCAssign :: CCAssign(CParseObj* vari)
{
  var = vari;
}

TParseObj CCAssign :: CommandoArg(int i) 
{
  if (i == 0) {return TPOAssign;}
  else {return TPOMathObj;} // for i==1
}
  
void CCAssign :: Add(CParseObj* parseObj, int i) 
{
  if (i == 1) {mathObj = (CMathObj*)parseObj;}
  //else *delete parseObj
};

void CCAssign :: Execute(CParser* /*parser*/)
{

  CMathObj* mo2;
  if (mathObj->EvaluableDouble())
  {
    mo2 = new CMNumber(mathObj->Double());
  } else
  {
      //copy obj included baseobj
    mathObj->AddBase();
    mo2 = mathObj->Parent()->CopyMO();

    if (var->ParseObjType() == TPOVar)
    {
      mo2->Eliminate((CVariable*)var);
    } else
    {
      CVariable* cv = ((CMVariable*)((CMMatElem*)var)->Arg(0))->Variable();
      mo2->Eliminate(cv);
    }
    mo2=mo2->Arg(0);
  }
  if (var->ParseObjType() == TPOVar)
  {
    ((CVariable*)var)->Assign(mo2);
  } else  //it can only be a MatElem function
  {
    ((CMMatElem*)var)->SetVar(mo2);
  }
}


//++++++++++++++++++++ CCFOR ++++++++++++++++++++++
CCFor::CCFor()
{
}

TParseObj CCFor::CommandoArg(int i)
{
  TParseObj a = TPOMathObj;
  if (i == 0) {a = TPOOpenBr;}
  if (i == 2 || i == 4 || i == 6) {a = TPOSemicolon;}
  if (i == 8) {a = TPOCloseBr;}
  if (i == 9) {a = TPOStatement;}
  return a;
}

void CCFor::Add(CParseObj* parseObj, int i)
{
  if (i == 1) 
  {
    if (parseObj->ParseObjType() == TPOMathObj)
    {
      CMathObj* mo = (CMathObj*)parseObj;
      if (mo->MOType() == TMOVariable)
      {
        var = ((CMVariable*)mo)->Variable();
      }
      else {MyError("Syntax Error, awaited countervariable");}
    }
    else {MyError("Syntax Error, awaited countervariable");}
  }
  else if (i == 3) {mathObj1 = (CMathObj*)parseObj;}
  else if (i == 5) {mathObj2 = (CMathObj*)parseObj;}
  else if (i == 7) {mathObj3 = (CMathObj*)parseObj;}
  else if (i == 9) {statement = parseObj;}
  //*delete parseObj;
}

void CCFor::Execute(CParser* parser)
{
  double d1 = mathObj1->Double();
  double d2 = mathObj2->Double();
  double d3 = mathObj3->Double();

  CMNumber* num = new CMNumber(d1);
  var->Assign(num);
  while (num->Double() <= d2 && !ParserErrorFlag())
  {
    ((CCommand*)statement)->Execute(parser);
    num->Set(num->Double()+d3);
  }
}

//++++++++++++++++++++ CCStatement ++++++++++++++++++++++
CCStatement::CCStatement() 
{
  statements = new AParseObjList();
}

CParseObj* CCStatement::GetCopy() 
{
  CCStatement* ccs = new CCStatement();

  //statements are not copied, this makes no sense!

  return ccs;
}


mystr CCStatement::Name() 
{
  return mystr("@Statement");
}

int CCStatement::GetNOArg(int actNOArg) 
{
  return actNOArg+1;
}

TParseObj CCStatement::EndArg() 
{
  return TPOCloseState;
}
  
TParseObj CCStatement::CommandoArg(int /*i*/)
{
  return TPOStatement;
}

void CCStatement::Add(CParseObj* parseObj, int /*i*/)
{
  if (parseObj->ParseObjType() != TPOSemicolon)
  {
    statements->Add(parseObj);
  }
}

void CCStatement::Execute(CParser* parser)
{
  int i = 1;
  while (i <= statements->Length() && !ParserErrorFlag())
  {
    if ((*statements)[i]->ParseObjType() == TPOCommand)
    {
      CCommand* state=(CCommand*)((*statements)[i]);
      state->Execute(parser);
    } else
    {
      MyError((mystr)"expected command instead of '"+
              (*statements)[i]->Name()+(mystr)"'");
    }
    i++;
  }
}

//++++++++++++++++++++ CCIf ++++++++++++++++++++++
CCIf::CCIf()
{
  statement1 = NULL;
  statement2 = NULL;
  withelse = 1;
}

int CCIf::GetNOArg(int /*actNOArg*/) 
{
  if (withelse) {return 6;}
  else {return 4;}
};

TParseObj CCIf::CommandoArg(int i)
{

  if (i == 0) {return TPOOpenBr;}
  else if (i == 1) {return TPOMathObj;}
  else if (i == 2) {return TPOCloseBr;}
  else if (i == 4) {return TPOTest;}
  else {return TPOStatement;}
}

void CCIf::Add(CParseObj* parseObj, int i)
{
  if (i == 1) {mathObj = (CMathObj*)parseObj;}
  else if (i == 3) {statement1 = parseObj;}
  else if (i == 4) 
  {
    if (parseObj->ParseObjType() != TPOElse)
    {
      withelse = 0;
    }
  }
  else if (i == 5) {statement2 = parseObj;}
}

void CCIf::Execute(CParser* parser)
{
  if (mathObj->Double() == 1)
  {
    ((CCommand*)statement1)->Execute(parser);
  } else if (statement2 != NULL)
  {
    ((CCommand*)statement2)->Execute(parser);
  }
}

//++++++++++++++++++++ CCWhile ++++++++++++++++++++++
CCWhile::CCWhile()
{
}

TParseObj CCWhile::CommandoArg(int i)
{

  if (i == 0) {return TPOOpenBr;}
  else if (i == 1) {return TPOMathObj;}
  else if (i == 2) {return TPOCloseBr;}
  else {return TPOStatement;}
}

void CCWhile::Add(CParseObj* parseObj, int i)
{
  if (i == 1) {mathObj = (CMathObj*)parseObj;}
  else if (i == 3) {statement = parseObj;}
}

void CCWhile::Execute(CParser* parser)
{
  const int maxit = 65536;
  int i=0;
  while (mathObj->Double() == 1 && i < maxit && !ParserErrorFlag())
  {
    ((CCommand*)statement)->Execute(parser);
    i++;
  }
  if (i == maxit)
  {MyError("Execution Error: to many iterations in while loop (65536) !");}
}

//++++++++++++++++++++ PRINT ++++++++++++++++++++++

TParseObj CCPrint::CommandoArg(int /*i*/) 
{
  return TPOString;
}
  
void CCPrint::Add(CParseObj* parseObj, int /*i*/) 
{
  str = (CPString*) parseObj;
}

void CCPrint::Execute(CParser* parser)
{
  parser->Write(str->String());
}
    
//++++++++++++++++++++ PRINT ++++++++++++++++++++++

TParseObj CCInclude::CommandoArg(int /*i*/) 
{
  return TPOString;
}
  
void CCInclude::Add(CParseObj* parseObj, int /*i*/) 
{
  str = (CPString*) parseObj;
}

void CCInclude::Execute(CParser* parser)
{
  parser->ParseFile(str->String());
}
  
//++++++++++++++++++++ LF ++++++++++++++++++++++

void CCLF::Execute(CParser* parser)
{
  parser->Write(mystr("\n"));
}
  
//++++++++++++++++++++ INFO ++++++++++++++++++++++

void CCInfo::Execute(CParser* parser)
{
  parser->WriteInfo();
}
  
//++++++++++++++++++++ DIFFERENTIATE ++++++++++++++++++++++

CCDifferentiate::CCDifferentiate()
{
  expr = NULL;
  result = NULL;
  var = NULL;
}

CParseObj* CCDifferentiate::GetCopy()
{
  CCDifferentiate* ccr= new CCDifferentiate();
  ccr->expr = expr;
  ccr->result = result;
  ccr->var = var;
  return ccr;
}

TParseObj CCDifferentiate::CommandoArg(int i)
{
  if (i == 0) {return TPOOpenBr;}
  else if (i == 1) {return TPOMathObj;}
  else if (i == 2) {return TPOComma;}
  else if (i == 3) {return TPOMathObj;}
  else if (i == 4) {return TPOComma;}
  else if (i == 5) {return TPOMathObj;}
  else {return TPOCloseBr;}
}

void CCDifferentiate::Add(CParseObj* parseObj, int i)
{
  if (i == 1) {expr = (CMathObj*)parseObj;}
  else if (i == 3) 
    {
      var = (CMathObj*)parseObj;
      if (var->MOType() != TMOVariable && 
          var->MOType() != TMOMatElem)
	{
	  MyError("Syntax Error: Differentiate awaits a variable as second parameter!");
	}
    }
  else if (i == 5) 
    {
      if (parseObj->ParseObjType() == TPOMathObj &&
	  ((CMathObj*)parseObj)->MOType() == TMOVariable)
	{
	  result = ((CMVariable*)parseObj)->Variable();
	}
    }
 }

void CCDifferentiate::Execute(CParser* parser)
{
  result->Assign(expr->Differentiate(var));
}


//$JG2013-5-10: calculate is not needed: use Print (evaluates and prints) and PrintSymbolic (prints formula/object) instead
/*
//++++++++++++++++++++ CCCalculate ++++++++++++++++++++++
CCCalculate::CCCalculate()
{
}
  
mystr CCCalculate::Name() 
{
  return mystr("calc");
}

int CCCalculate::GetNOArg(int actNOArg) 
{
  return 1;
}

TParseObj CCCalculate::CommandoArg(int i) 
{
  return TPOMathObj;
}

void CCCalculate::Add(CParseObj* parseObj, int i) 
{
  mathObj = (CMathObj*) parseObj;
}

void CCCalculate::Execute(CParser* parser)
{
  int2 rc;
  mathObj->Dim(rc);
  if (rc.Get(1) == 1 && rc.Get(2) == 1)
  {
    parser->Write(mathObj->String()+(mystr)" = "+mystr(mathObj->Double())+mystr("\n"));
  } else
  {
    mystr str = mathObj->String()+mystr(" = \n  (");
    for (int j = 1; j <= rc.Get(1); j++)
    {
      for (int i = 1; i <= rc.Get(2); i++)
      {
        str += mathObj->MatDouble(j,i);
        if (!(i == rc.Get(2))) {str += mystr(",");}
      } 
      if (j == rc.Get(1)) {str += mystr(")\n");}
      else {str += ";\n   ";}
    }
    parser->Write(str);
    
  }
  
}
*/
//$JG2013-5-10: not needed any more, use EDC commands in future
/*
//++++++++++++++++++++ CCPrintm ++++++++++++++++++++++
CCPrintm::CCPrintm()
{
}
  
mystr CCPrintm::Name() 
{
  return mystr("printm");
}

int CCPrintm::GetNOArg(int actNOArg) 
{
  return 1;
}

TParseObj CCPrintm::CommandoArg(int i) 
{
  return TPOMathObj;
}

void CCPrintm::Add(CParseObj* parseObj, int i) 
{
  mathObj = (CMathObj*) parseObj;
}

void CCPrintm::Execute(CParser* parser)
{
  //CMathObj* smo = mathObj->StdSimplify();
  parser->Write(mathObj->StringAll()+mystr("\n"));
  //parser->Write(smo->StringAll()+mystr("\n"));
}
*/


