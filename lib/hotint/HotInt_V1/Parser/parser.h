
//#**************************************************************
//#
//# filename:             parser.h
//#
//# autor:                Gerstmayr Johannes
//#
//# generated:            23.11.98
//# last change:          20.12.98
//# description:          base class for parser
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

#ifndef PARSER__H
#define PARSER__H

#include "parseobj.h"
#include "mathobj.h"
#include "parseobjlist.h"
#include "command.h"
#include "variable.h"


void MyError();

void SysError(mystr string1);

void MyError(mystr string1);

void MyWarning(mystr string1);

//the CParser has the following functionality
//  - convert string into list of statements (assignments, commands, etc.): ParseString(...)
//    and executes the object tree if commands are contained ==> CHANGE?
//  - convert list of statements into CParseObjList
//  - convert string into MathObj Tree: GetMathObjTree()
class CParser
{
public:
  CParser();

	virtual int ParserErrorFlag();
	virtual void SetParserErrorFlag(int flag);


  virtual void ParseFile(const mystr& file);

	//initialization: register all available functions and commands
  virtual void GenParseObjList();
	//register a single parsable object
  virtual void AddParseObj(CParseObj* cpo);

  virtual void Write(mystr str);
  //initialise all for new parsing
  virtual void Init();

	//starting routine for parsing, generates tree, executes, etc.
	//parse string (multi-line) into commands, definitions, assignments, etc.
  virtual void ParseString(mystr parsestring);
	//parse a string into a CMathObj (String2MathObj) and compute to double
	virtual double CalculateString(mystr parsestring);
	//convert string into a CMathObj
	virtual CMathObj* String2MathObj(mystr parsestring, int parselinepos = 0);
 
	//parses a statement until ';', EOF, '}'
  //parses recursively parseobjects
  //returns NULL, if end.
  virtual CParseObj* ParseObjTree();
	//read a CMathObj-structure (error, if other than CMathObj; terminated by other objects (e.g. else))
  virtual CMathObj* GetMathObjTree();
	//read a single mathematical object 
  virtual CMathObj* GetMObject();
  //get single object of type getobject, error if object type does not match
  virtual CParseObj* GetObject(TParseObj getobject);
  //get a single object (non-CMathObj) until separator or get separator itself
  virtual CParseObj* GetObject();

    //speichert letzte parseposition im mystr
  virtual void StorePos();
  virtual void RestorePos();

  virtual void ReadSpaces();

    // read comment beginning with // ending at end of line
  virtual void ReadLineComment();
    // read comment of form /**/
  virtual void ReadComment();

  virtual int ReadChar(char& ch);

  virtual int PeekChar(char& ch);
    //read mystr until '"', return mystr without '"'
  virtual mystr ReadString();

  virtual int IsNumber(char ch);

  virtual int IsNumberChar(char ch);

  virtual int IsCharacter(char ch);

  virtual int IsSeparator(char ch);

  virtual int IsEOL(char ch);

  virtual int IsEOF(char ch);

  virtual int IsNumber(const mystr& str);

  virtual int IsVariable(const mystr& str);

  virtual double Number(const mystr& str);

	virtual CVariableList* GetVariableList() {return variableList;}

	virtual CParseObj* GetVariable(mystr& str)
	{
		return GetVariableList()->Object(str);
	}

  virtual void WriteInfo();
  virtual void Error(const mystr& string1);

  virtual void Error(const mystr& string1, const mystr& string2);

private:
  mystr actstring;
  int storeactpos, storelines, storeactlinepos;
  int actpos,      lines,      actlinepos;
  int endpos;
  CParseObjList* parseObjList; //this is the list of available objects of the parser
  CVariableList* variableList; //this is the global list of variables in the parser

  //ofstream* outFile;  //JG2013-04-29 erased, not necessary
};


//parsed function is represented by CMathObj-Tree, which contains mathematical structure and links to variables
class CParsedFunction
{
public:
	CParsedFunction() {InitConstructor();};

	void InitConstructor()
	{
		mathobj = 0;
		actual_parser = 0;
		parsed_function = "";
		parsed_function_parameter = "";
	}
	//$ RL 2011-5-31:[
	CParsedFunction& operator=(const CParsedFunction& e) 
	{
		if (this == &e) {return *this;}
		CopyFrom(e);
		return *this;
	}
	virtual void CopyFrom(const CParsedFunction& e);
	//$ RL 2011-5-31:]

	//$JG2013-4-29: Destructor missing!

	//set a one-dimensional function with one parameter; the string is parsed and stored in "mathobj"
	void SetParsedFunction1D(CParser* parser, const mystr& mathstring, const mystr& parameter_name);
	
	//evaluate the parsed mathobj for the given parameter
	double Evaluate(const double& parameter) const; //$ RL 2011-5-31: const added.

	virtual const CVariableList& GetVariableList() const {return variableList;}
	virtual CVariableList& GetVariableList() {return variableList;}

	virtual const mystr& GetParsedFunction() const {return parsed_function;}
	virtual const mystr& GetParsedFunctionParameter() const {return parsed_function_parameter;}

private:
	TArray<CMNumber*> parameters;//$ RL 2011-5-31: pointers to parameters are used now (they keep constant in Evaluate function)
  CVariableList variableList;	
	CMathObj* mathobj;

	mystr parsed_function;  //storage of parsed function
	mystr parsed_function_parameter; //storage of parsed parameter

	CParser* actual_parser;   //$ PG 2013-8-13: Bugfix @{//$JG2013-4-29: possible workaround, copy pointers:} ... but then the variablesList has to be adapted also, otherwise modifications in parameters will not take effect
};

//parser linked to mbs:
//in the .h header file, no MBS* operations may be done! otherwise, WorkingModule.dll can not be loaded!!!
class CMBSParser: public CParser
{
public:
	CMBSParser();

	void SetMBS(MBS * mbsi);

	virtual void Init();

	virtual CParseObj* GetVariable(mystr& str);

	virtual int NElementDataContainers() const
	{
		return 5; 
	}

	virtual void SetLocalEDC(ElementDataContainer* local_edci);
	virtual void SetLocalEDC2(ElementDataContainer* local_edci);
	virtual ElementDataContainer* GetLocalEDC() const;
	virtual ElementDataContainer* GetLocalEDC2() const;

	virtual ElementDataContainer* GetElementDataContainers(int i);

	double ExpressionToDouble(mystr & data, int& error_flag, int line = 0);
	int2 EvaluateExpression(mystr & data, int& error_flag, double& d_var, Vector& v_var, Matrix& m_var, int line = 0);

private:
	//this list contains all EDCs linked to MBS
	//sorted regarding priority
	//local_edc (linked only if model is parsed)
	//global_model_edc
	//mbs_options_edc
	//global_variables_edc (pi, physical constants, materials): pi, mat.steel.Em, mat.alu.Em, etc.
	
	//TArray<ElementDataContainer*> edcs; 

	MBS * mbs;
	ElementDataContainer* local_edc; //this is the local (sub-) data container, which is currently read
	ElementDataContainer* local_edc2; //this is the root of the currently read data container

};

/*
//not used up to now
class CEDVariable : public CVariable
{
//CMVariable--> override: mathobj->Double() bzw. mathobj->Evaluate 
//derive ElementData from CVariable, not CMVariable, use with EDC functions

public:
	CEDVariable(): CVariable() {elementdata = 0;};
	CEDVariable(const mystr& namei):CVariable(namei)
	{
		elementdata = 0;
	}

  virtual TParseObj ParseObjType() {return TPOVar;}
  virtual CParseObj* GetCopy() 
  {
    CEDVariable* o = new CEDVariable();
    o->mathObj = mathObj;
    o->name = name;
		o->elementdata = elementdata;
    return o;
  }

	virtual double GetDouble()
	{
		if (elementdata->IsDouble()) {return elementdata->GetDouble();}
		else if (elementdata->IsInt()) {return elementdata->GetInt();}
		else if (elementdata->IsBool()) {return elementdata->GetBool();}
		
		assert(0);
		return 0;
	}
	virtual int EvaluableDouble() {return elementdata != 0;}

	virtual void Assign(ElementData* ed)
	{
		elementdata = ed;
	}

private:
	ElementData* elementdata;
};
*/


#endif
