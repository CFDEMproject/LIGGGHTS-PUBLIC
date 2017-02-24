
//#**************************************************************
//#
//# filename:             command.h
//#
//# autor:                Gerstmayr Johannes
//#
//# generated:            22.11.98
//# last change:          20.12.98
//# description:          base class for commands
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

#ifndef COMMAND__H
#define COMMAND__H

class CParser;

//a command can be executed by parser
//a code block usually consists of a list of commands, definitions (functions, maybe variables as well in future) 
//and variable assignments
class CCommand : public CParseObj
{
public:
  CCommand() {};
  virtual mystr Name() {return mystr("@CCommand");}
  virtual CParseObj* GetCopy() {return new CCommand();}
  virtual TParseObj ParseObjType() {return TPOCommand;}
  virtual int GetNOArg(int actNOArg) {return 0;}; //number of arguments needed for command

  virtual TParseObj EndArg() {return TPOVoid;} //only used in Statement: EndArg="}" ==> maybe erase this feature in future?
  virtual TParseObj CommandoArg(int i) {return TPOVoid;} //returns a specific argument type for each i-th argument: "(", ",", MathObj, ")", ";", etc. 
  virtual void Add(CParseObj* parseObj, int i) {}; //process i-th argument (e.g. store MathObj) or do nothing with argument (e.g. with ",")

  virtual void Execute(CParser* /*parser*/) {}; //execute command: this is done during execution of parsed objects
};


class CCAssign : public CCommand
{
public:
  CCAssign(CParseObj* vari);

  virtual CParseObj* GetCopy() 
  {
    CCAssign* o = new CCAssign(var);
    o->mathObj = mathObj;
    return o;
  }

  virtual mystr Name() {return mystr("@CCAssign");}

  virtual int GetNOArg(int actNOArg) {return 2;};

  virtual TParseObj CommandoArg(int i);
  
  virtual void Add(CParseObj* parseObj, int i) ;

	//assign variable: if variable content can be evaluated ==> evaluate to double
	//  otherwise, check if recursive assignment (e.g. a=a and a is previously undefined)
  virtual void Execute(CParser* parser);
  
private:
  CMathObj* mathObj;
  CParseObj* var;
};

//has structure: for '(' countervariable (eg. i) ';' 
//  startvalue for which statement is first executed ';' 
//  lastvalue for which statement is executed ';' 
//  increasement after one step ')' statement
class CCFor : public CCommand
{
public:
  CCFor();

  virtual CParseObj* GetCopy() 
  {
    CCFor* o = new CCFor();
    o->var = var;
    o->mathObj1 = mathObj1;
    o->mathObj2 = mathObj2;
    o->mathObj3 = mathObj3;
    o->statement = statement;
    return o;
  }

  virtual mystr Name() {return mystr("for");}

  virtual int GetNOArg(int actNOArg) {return 10;};

  virtual TParseObj CommandoArg(int i);

  virtual void Add(CParseObj* parseObj, int i);

  virtual void Execute(CParser* parser);

private:
  CVariable* var;
  CMathObj* mathObj1;
  CMathObj* mathObj2;
  CMathObj* mathObj3;
  CParseObj* statement;
};

class  CCStatement : public CCommand
{
public:
  CCStatement();

  virtual CParseObj* GetCopy();

  virtual mystr Name();
  //first arg is statement after "{" !
  virtual int GetNOArg(int actNOArg);

  virtual TParseObj EndArg();

  virtual TParseObj CommandoArg(int i);

  //End-statement is not added
  virtual void Add(CParseObj* parseObj, int i);

  virtual void Execute(CParser* parser);

private:
  AParseObjList* statements;
};


class CCIf : public CCommand
{
public:
  CCIf();

  virtual CParseObj* GetCopy() 
  {
    CCIf* o = new CCIf();
    o->mathObj = mathObj;
    o->statement1 = statement1;
    o->statement2 = statement2;
    return o;
  }

  virtual mystr Name() {return mystr("if");}

  virtual int GetNOArg(int actNOArg);

  virtual TParseObj CommandoArg(int i);

  virtual void Add(CParseObj* parseObj, int i);

  virtual void Execute(CParser* parser);

private:
  CMathObj* mathObj;
  CParseObj* statement1;
  CParseObj* statement2;
  int withelse;
};

class CCWhile : public CCommand
{
public:
  CCWhile();

  virtual CParseObj* GetCopy() 
  {
    CCWhile* o = new CCWhile();
    o->mathObj = mathObj;
    o->statement = statement;
    return o;
  }

  virtual mystr Name() {return mystr("while");}

  virtual int GetNOArg(int actNOArg) {return 4;};

  virtual TParseObj CommandoArg(int i);

  virtual void Add(CParseObj* parseObj, int i);

  virtual void Execute(CParser* parser);

private:
  CMathObj* mathObj;
  CParseObj* statement;
};

  //write mystr without carrige return
class CCPrint : public CCommand
{
public:
  CCPrint() : CCommand() {};
  virtual mystr Name() {return mystr("print");};
  virtual CParseObj* GetCopy() 
  {
    CCPrint* o = new CCPrint();
    o->str = str;
    return o;
  }
  virtual int GetNOArg(int actNOArg) {return 1;};

  virtual TParseObj CommandoArg(int i);
  virtual void Add(CParseObj* parseObj, int i);

  virtual void Execute(CParser* parser);
  
private:
  CPString* str;
};
  //include a file and execute it 
class CCInclude : public CCommand
{
public:
  CCInclude() : CCommand() {};
  virtual mystr Name() {return mystr("include");};
  virtual CParseObj* GetCopy() 
  {
    CCInclude* o = new CCInclude();
    o->str = str;
    return o;
  }
  virtual int GetNOArg(int actNOArg) {return 1;};

  virtual TParseObj CommandoArg(int i);
  virtual void Add(CParseObj* parseObj, int i);

  virtual void Execute(CParser* parser);
  
private:
  CPString* str;
};

  //print line feed to output file
class CCLF : public CCommand
{
public:
  CCLF() : CCommand() {};
  virtual mystr Name() {return mystr("lf");};
  virtual CParseObj* GetCopy() {return new CCLF();}
  virtual int GetNOArg(int actNOArg) {return 0;}
  virtual TParseObj CommandoArg(int i) {return TPOVoid;}
  virtual void Add(CParseObj* parseObj, int i) {}
  virtual void Execute(CParser* parser);
};

  //print line feed to output file
class CCInfo : public CCommand
{
public:
  CCInfo() : CCommand() {};
  virtual mystr Name() {return mystr("info");};
  virtual CParseObj* GetCopy() {return new CCInfo();}
  virtual int GetNOArg(int actNOArg) {return 0;}
  virtual TParseObj CommandoArg(int i) {return TPOVoid;}
  virtual void Add(CParseObj* parseObj, int i) {}
  virtual void Execute(CParser* parser);
};

//Differentiate and put into variable
class CCDifferentiate : public CCommand
{
public:
  CCDifferentiate();

  virtual CParseObj* GetCopy();

  virtual mystr Name() {return mystr("diff");}

  virtual int GetNOArg(int actNOArg) {return 7;}

  virtual TParseObj CommandoArg(int i);

  virtual void Add(CParseObj* parseObj, int i);

  virtual void Execute(CParser* parser);

private:
  CMathObj* expr;
  CMathObj* var;
  CVariable* result;
};

/*class CCCalculate : public CCommand
{
public:
  CCCalculate();
  virtual mystr Name();
  virtual CParseObj* GetCopy() 
  {
    CCCalculate* o = new CCCalculate();
    o->mathObj = mathObj;
    return o;
  }
  virtual int GetNOArg(int actNOArg);

  virtual TParseObj CommandoArg(int i);
  virtual void Add(CParseObj* parseObj, int i);

  virtual void Execute(CParser* parser);
  
private:
  CMathObj* mathObj;
};

class CCPrintm : public CCommand
{
public:
  CCPrintm();
  virtual mystr Name();
  virtual CParseObj* GetCopy() 
  {
    CCPrintm* o = new CCPrintm();
    o->mathObj = mathObj;
    return o;
  }
  virtual int GetNOArg(int actNOArg);

  virtual TParseObj CommandoArg(int i);
  virtual void Add(CParseObj* parseObj, int i);

  virtual void Execute(CParser* parser);
  
private:
  CMathObj* mathObj;
};

*/


#endif
