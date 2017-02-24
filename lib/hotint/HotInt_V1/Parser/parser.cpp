
//#**************************************************************
//#
//# filename:             parser.cpp
//#
//# autor:                Gerstmayr Johannes
//#
//# generated:            23.11.98
//# last change:          23.11.98
//# description:          implementation for parser
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
#include "myfile.h"

extern UserOutputInterface * global_uo;

//ofstream pout("aoutput.txt");

int errorflag=0;

int glob_objcnt = 0;

//const int& ParserErrorFlag() {return errorflag;}
int ParserErrorFlag() {return errorflag;}
void SetParserErrorFlag(int flag) {errorflag = flag;}

int CParser::ParserErrorFlag() 
{
	return ::ParserErrorFlag();
}

void CParser::SetParserErrorFlag(int flag) 
{
	::SetParserErrorFlag(flag);
}


void MyError()
{
  SetParserErrorFlag(1);
  //(*global_uo) << "test";
  //(*global_uo) << "Error occured" << endl << flush;
}

void SysError(mystr string1)
{
  if (!ParserErrorFlag())
  {
    (*global_uo) << "Parser error: " << string1 << " \n";
  }
  MyError();
}

void MyError(mystr string1)
{
  if (!ParserErrorFlag())
  {
    (*global_uo) << string1 << "\n";
  }
  MyError();
}

void MyWarning(mystr string1)
{
  (*global_uo) << string1 << "\n";
}

CParser::CParser()
{
  GenParseObjList();
  
  variableList = new CVariableList();
  Init();

  //outFile = &pout;  //JG2013-04-29
  //(*outFile) << "solutionfile for parser:" << flush << endl;
}

  //starting routine for parsing, generates tree, executes, etc.
void CParser::ParseString(mystr parsestring)
{
  //if (ParserErrorFlag()) {return;}
    
  //  (*global_uo) << "parsestring=" << parsestring << flush;
  actstring = mystr("{\n")+parsestring+mystr("}");
  actpos = 0;
  actlinepos = 0;
  lines = 0;
  
  endpos = actstring.Length(); // one pos after end!

    //recursively generate parseobjtree
  CParseObj* parseTree = ParseObjTree();

	
    //execute parsetree
  if (!ParserErrorFlag() && parseTree->ParseObjType() == TPOCommand)
  {
    ((CCommand*)parseTree)->Execute(this);
  }
	
}

CMathObj* CParser::String2MathObj(mystr parsestring, int parselinepos)
{
  actstring = parsestring;
  actpos = 0;
  actlinepos = parselinepos;
  lines = 0;
  
  endpos = actstring.Length(); // one pos after end!

	return GetMathObjTree();
}

double CParser::CalculateString(mystr parsestring)
{
 // actstring = parsestring;
 // actpos = 0;
 // actlinepos = 0;
 // lines = 0;
 // 
 // endpos = actstring.Length(); // one pos after end!


	//CMathObj* mo = GetMathObjTree();

	CMathObj* mo = String2MathObj(parsestring);
	
    //execute parsetree
  if (mo != 0 && !ParserErrorFlag()) // && parseTree->ParseObjType() == TPOMathObj)
  {
		return mo->Double();
  }
	return 0;	
}

  //parses a statement until ';', EOF, '}'
  //parses recursively parseobjects
  //returns NULL, if end.
CParseObj* CParser::ParseObjTree()
{
  //if (ParserErrorFlag()) {return new CPVoid();}
    
  if (actpos < endpos)
  {
    TParseObj comname = (TParseObj)(TPOCommand + TPOVar + TPOCloseState + TPOOpenState
                        + TPOVoid + TPOSemicolon);

    CParseObj* obj = GetObject(comname);
//pout << "  objname=" << obj->Name() << endl;
    CCommand* com = NULL;

      //read commando start word
    if (obj->ParseObjType() == TPOVar)
    {
      //check if MatElem occurs ([i,j])
      StorePos();
      CParseObj* obj2 = GetObject();
      if (obj2->ParseObjType() == TPOOpenArgBr)
      {
        CMathObj* cmr = GetMathObjTree();
        GetObject(TPOComma);
        CMathObj* cmc = GetMathObjTree();
        GetObject(TPOCloseArgBr);
        CMVariable* cmvar = new CMVariable((CVariable*)obj);
        
        CMMatElem* me = new CMMatElem();
        me->SetArg(cmvar,0);
        me->SetArg(cmr,1);
        me->SetArg(cmc,2);
        com = new CCAssign(me);
      } else 
      {
        RestorePos();
        com = new CCAssign((CVariable*)obj);
      }
      
      
    } else if (obj->ParseObjType() == TPOCommand)
    {
      com = (CCommand*)obj;
    } else if (obj->ParseObjType() == TPOCloseState)
    {
      return new CCloseState();
    } else if (obj->ParseObjType() == TPOOpenState)
    {
      com = new CCStatement();
    } else if (obj->ParseObjType() == TPOSemicolon)
    {
      return obj;
    }

      //read commando parameters
    if (com != NULL)
    {
      int i = 0;
      int end = 0;
      while (!end && i < com->GetNOArg(i) && !ParserErrorFlag())
      {
        TParseObj tcomarg = com->CommandoArg(i);
        CParseObj* comarg;

        if (tcomarg == TPOStatement)
        {
          comarg = ParseObjTree();
        } else // math-object
        if (tcomarg == TPOMathObj)
        {
          comarg = GetMathObjTree();
        } else // test-object
        if (tcomarg == TPOTest)
        {
          StorePos();
          comarg = GetObject();
        } else // other objects like ',', '(', ...
        {
          comarg = GetObject(tcomarg);
        }

        if (comarg != NULL)
        {
          if (comarg->ParseObjType() == com->EndArg()) 
          {
            end = 1;
          } else
          {
            com->Add(comarg, i);
              //if test-object failed, NOArg is <= than i
            if (com->GetNOArg(i) <= i && tcomarg == TPOTest) 
            {
              RestorePos();
              end = 1;
            }
          }
        } else
        {
          Error("Syntax Error, no arguments");
        }
        i++;
      } 
    } else
    {      
      Error("Syntax Error, unknown statement");
      return NULL;
    }
    
    return com; 

  } else 
  {
    return NULL;
  }    
}
CMathObj* CParser::GetMathObjTree()
{
  CMathObj* base = new CBaseObj();
  CMathObj* actobj;
  CMathObj* obj;
  TMathObj await = TMOOperand;
  int end = 0;
  int newloop;
  int brackets = 0;

  actobj = base;

  while (!end && !ParserErrorFlag())
  {
//(*global_uo) << "   actobj=" << actobj->Name() << " par=";
//(*global_uo) << actobj->Parent() << endl << flush;
//(*global_uo) << " s=" << base->String() << endl;
    
    newloop = 1;
    StorePos();
    obj = GetMObject();
//(*global_uo) << "obj=" << obj << endl;

    if (await == TMOOperand && newloop)
    {
      if (obj->MOType() & TMOOperand)
      {
        if (obj->MOType() == TMONumber ||
            obj->MOType() == TMOVariable)
        {
          actobj->AddArg(obj);
          actobj = obj;
          await = TMOOperator;
        } else
        if (obj->MOType() == TMOFunction)
        {
          actobj->AddArg(obj);
          actobj = obj;
        } else
        if (obj->MOType() == TMOOpenBr)
        {
          CMathObj* cm = new CMBracket();
          actobj->AddArg(cm);
          actobj = cm;
          brackets++;
          //*delete obj;
        } else
        if (obj->MOType() == TMOMatrix)
        {
          actobj->AddArg(obj);
          actobj = obj;

          GetObject(TPOOpenBr);
          CMathObj* rows = GetMathObjTree();
          GetObject(TPOComma);
          CMathObj* cols = GetMathObjTree();

          CMatrix* cm = (CMatrix*)obj;
          cm->SetDimensions(rows, cols);

          int mend = 0;
          
          while (!mend)
          {
            CParseObj* cp = GetObject();
            if (cp->ParseObjType() != TPOCloseBr)
            {
              if (cp->ParseObjType() == TPOComma ||
                  cp->ParseObjType() == TPOSemicolon)
              {
                cm->AddElem(GetMathObjTree());
              } else
              {
                Error("Syntax Error: closing bracket or Comma awaited!");
                mend = 1;
              }
            } else {mend = 1;}
          }
          
          await = TMOOperator;          
        }         
        newloop=0;
      } else if (obj->MOType() == TMOOperator &&
                 obj->Name()==mystr("+"))
      {
        //can be left away, its dummy!!!
        /*CMathObj* cm = new CMUnPlus();
        actobj->AddArg(cm);
        actobj = cm;*/
        //*delete obj;      
      } else if (obj->MOType() == TMOOperator &&
                 obj->Name()==mystr("-"))
      {
        CMathObj* cm = new CMUnMinus();
        actobj->AddArg(cm);
        actobj = cm;
        //*delete obj;      
      } else
      {
        Error((mystr)"Syntax Error, awaited operand, not "+obj->Name());
      }
    }

    if (newloop)
    {
      if (await == TMOOperator)
      {
          //else UnOperator!!!!
        if (obj->MOType() == TMOOperator)
        {
          //find right inputposition ("point before line"):
          COperator* op = (COperator*)obj;
          TMathPriority p = op->Priority();

          int priend = 0;
          while (!priend)
          {
//pout << "actobj=" << actobj->Name() << ", pri=" << actobj->Priority() << "\n";
            if (actobj->Parent()->Priority() >= p)
            {
              actobj = actobj->Parent();
            } else {priend = 1;}
            
          }

          //insert operator: (comp with (*)!!!)
          CMathObj* parent = actobj->Parent();
          obj->SetArg(actobj,0);
          parent->ReplaceArg(actobj,obj);
          actobj = obj;
          await = TMOOperand;

          newloop = 0;
        } else
        if (obj->MOType() == TMOMatElem)
        {
          obj->SetArg(GetMathObjTree(),1);
          GetObject(TPOComma);
          obj->SetArg(GetMathObjTree(),2);
          GetObject(TPOCloseArgBr);

            //is like in operator:  (*)        
          CMathObj* parent = actobj->Parent();
          obj->SetArg(actobj,0);
          parent->ReplaceArg(actobj,obj);
          actobj = obj;
          await = TMOOperator;

          newloop = 0;
        } else
        if (obj->MOType() == TMOCloseBr && brackets > 0)
        {
          brackets--;
          int end = 0;
          while (!end)
          {
            if (actobj->Parent() != NULL)
            {
              actobj = actobj->Parent();
              if (actobj->MOType() == TMOBracket)
              {
                CMathObj* lasto=actobj->Arg(0);
                CMathObj* br = actobj;
                actobj = actobj->Parent();
                actobj->ReplaceArg(br, lasto);
                actobj = lasto;
                //*delete br;
                end = 1;
              }
            } else 
            {
              Error("Syntax Error, too many Closing brackets");
              end = 1;
            }
          }
          
          
        } else
        if (obj->MOType() == TMOEnd ||
            obj->MOType() == TMOVariable ||
            obj->MOType() == TMOCloseBr)
        {
          RestorePos();
          newloop = 0;
          end = 1;
        } else
        {
          Error((mystr)"Syntax Error, awaited operator, not "+obj->Name());
        }
      }
    }
    if (newloop) {} //dummy

  }

  if (!ParserErrorFlag())
  {
    actobj = base->Arg();
    //actobj->SetParent(NULL);
    //delete base;
    return actobj;
  }
	else
	{
		global_uo->SaveLocalMessageLevel();
		global_uo->SetLocalMessageLevel(UO_LVL_err);
		(*global_uo) << "ERROR: while parsing line " << actlinepos << ".\n";
		global_uo->ResetLocalMessageLevel();
		
		return new CMEnd();
	}

}


CMathObj* CParser::GetMObject()
{
  CParseObj* o = GetObject();

//  pout << "mo: " << o->Name() << endl << flush;

  switch (o->ParseObjType())
  {
    case TPOMathObj: 
      return (CMathObj*)o; 
      //break;
    case TPOOpenBr:
      delete o;
      return new CMOpenBr();
      //break;
    case TPOCloseBr:
      delete o;
      return new CMCloseBr();
      //break;
    case TPOVar:
      return new CMVariable((CVariable*)o);
      //break;
    case TPOOpenArgBr:
      delete o;
      return new CMMatElem();
      //break;
    default: 
      if (o->ParseObjType() == TPOEndObj ||
          o->ParseObjType() == TPOCommand ||
          o->ParseObjType() == TPOComma ||
          o->ParseObjType() == TPOCloseState ||
          o->ParseObjType() == TPOOpenState ||
          o->ParseObjType() == TPOElse ||
          o->ParseObjType() == TPOCloseArgBr ||
          o->ParseObjType() == TPOSemicolon)
      {
        return new CMEnd();
      } else
      {
        Error("Syntax Error, not a valid mathematical expression", o->Name());
        return new CMathObj();
      }
  }

}

  //get object of type getobject
CParseObj* CParser::GetObject(TParseObj getobject)
{
  CParseObj* obj;

  obj = GetObject();

//  pout << "o: " << obj->Name() << endl << flush;

  if (obj->ParseObjType() & getobject) {return obj;}
  else
  {
    mystr str;
    Error((mystr)"Syntax Error, unexpected '" + obj->Name()+(mystr)"' ");
    return obj;
  }

}

  //get an object until separator or separator itself
CParseObj* CParser::GetObject()
{
  //object bis zum nächsten Trennsymbol lesen!!!!

  glob_objcnt++;
  
  mystr readstring = "";
  ReadSpaces();

  char actchar, nextchar;
  int end;

  if (!ReadChar(actchar)) {return new CEndObj();}
  readstring+=actchar;
  end = !PeekChar(nextchar);

    //read comment
  if (actchar == '/')
  {
    if (nextchar == '/')
    {
      ReadChar(actchar);
      ReadLineComment();
      return GetObject();
    } else
    if (nextchar == '*')
    {
      ReadChar(actchar);
      ReadComment();
      return GetObject();
    }
  }
  
  if (actchar == '"')
  {
    return new CPString(ReadString());
  }
    //read objects which consist of one or two separators
  if (IsSeparator(actchar))
  {
    //test for double-separator-char object (==,:= , <=, >=, etc.)
    mystr teststring = readstring + (mystr)nextchar;
    int iso_test = parseObjList->IsObject(teststring);
    int iso_read = parseObjList->IsObject(readstring);
    
    if (iso_test)
    {
      ReadChar(nextchar); //was previously only peeked!!
      return parseObjList->Get(iso_test)->GetCopy();
    }
    else if (iso_read)
    {
      return parseObjList->Get(iso_read)->GetCopy();
    }
    else {Error("Syntax Error",readstring); return new CPVoid();}
  }

  if (IsNumberChar(actchar))
  {
    while (!end && !IsSeparator(nextchar) && IsNumberChar(nextchar))
    {
      if (!ReadChar(actchar)) 
      {
        MyError("Error: unexpected end of parsing");
        return new CEndObj();
      }
      readstring+=actchar;
      end = !PeekChar(nextchar);
    }
    ReadSpaces();
    StorePos();
    mystr storestring = readstring;

      //test, if exponent follows
    PeekChar(nextchar);
    if (nextchar == 'e' || nextchar == 'E')
    {
      ReadChar(actchar);
      readstring+=actchar;
      ReadSpaces();
      PeekChar(nextchar);
      if (nextchar == '-' || nextchar == '+')
      {
        ReadChar(actchar);
        readstring+=actchar;
        ReadSpaces();
        PeekChar(nextchar);
      }
        //if no number after exponent->e was form another object (eg. else)
      if (!IsNumber(nextchar)) 
      {
        RestorePos();
        readstring = storestring;
        end = 1;
        //MyError("Error: expected number after exponent sign");
        //return new CEndObj();
      }
      while (!end && !IsSeparator(nextchar) && IsNumber(nextchar))
      {
        if (!ReadChar(actchar)) 
        {
          MyError("Error: unexpected end of parsing");
          return new CEndObj();
        }
        readstring+=actchar;
        end = !PeekChar(nextchar);
      }
    }

    if (IsNumber(readstring))
    {
      return new CMNumber(Number(readstring));
    } else
    {
      return new CEndObj();
    }
  }

  while (!end && !IsSeparator(nextchar))
  {
      if (!ReadChar(actchar)) 
      {
        MyError("Error: unexpected end of parsing");
        return new CEndObj();
      }
    readstring+=actchar;
    end = !PeekChar(nextchar);
  }
  int iso_read = parseObjList->IsObject(readstring);

  if (iso_read)
  {
    return parseObjList->Get(iso_read)->GetCopy();
  } 
	else if (IsVariable(readstring))
  {
//!AD: 2012-11-29 for right hand side access to vector or matrix component: [
		// if variable name is followed by a '[', continue to read, string will be parsed to vector access in GetVariable...
		if (!end && nextchar=='[')
		{  
			int brackets_open = 0;
			while (!end && !(actchar==']' && brackets_open==0))
			{
				if (!ReadChar(actchar)) 
				{
					MyError("Error: unexpected end of parsing");
					return new CEndObj();
				}
				readstring+=actchar;
				if (actchar=='[') brackets_open++;
				if (actchar==']') brackets_open--;
				end = !PeekChar(nextchar);
			}
		//	if(!end) ReadChar(actchar);
			return GetVariable(readstring); // use sabe GetVariable routine - could also implement a GetVariableComponent-like function
//!AD ]
		}
		else
		{
		  return GetVariable(readstring);
    //return GetVariableList()->Object(readstring);
		}
  }
  else {Error("Syntax Error", readstring); return new CPVoid();}

}

void CParser::ParseFile(const mystr& file)
{
  //  (*global_uo) << "parsefile\n" << flush;
  CMFile fin(file, TFMread);
  mystr parsestr;

  fin.RWF(parsestr);

//  (*outFile) << "string = '" << parsestr << "'" << endl;
  
  ParseString(parsestr);
}

void CParser::GenParseObjList()
{
  parseObjList = new CParseObjList();
  CParseObj* po;

  po = new CAssign(); parseObjList->Add(po);
  po = new CComma(); parseObjList->Add(po);
  po = new CSemicolon(); parseObjList->Add(po);

  po = new COpenBr(); parseObjList->Add(po);
  po = new CCloseBr(); parseObjList->Add(po);
  po = new COpenState(); parseObjList->Add(po);
  po = new CCloseState(); parseObjList->Add(po);
  po = new COpenArgBr(); parseObjList->Add(po);
  po = new CCloseArgBr(); parseObjList->Add(po);

  //po = new CCCalculate(); parseObjList->Add(po);
  po = new CCFor(); parseObjList->Add(po);
  po = new CCIf(); parseObjList->Add(po);
  po = new CPElse(); parseObjList->Add(po);
  po = new CCWhile(); parseObjList->Add(po);
  po = new CCPrint(); parseObjList->Add(po);
  //po = new CCPrintm(); parseObjList->Add(po);
  po = new CCLF(); parseObjList->Add(po);
  
  po = new CMEqual(); parseObjList->Add(po);
  po = new CMGrTh(); parseObjList->Add(po);
  po = new CMLeTh(); parseObjList->Add(po);
  po = new CMGrEq(); parseObjList->Add(po);
  po = new CMLeEq(); parseObjList->Add(po);

  po = new CMAnd(); parseObjList->Add(po);
  po = new CMOr(); parseObjList->Add(po);
  po = new CMCondAnd(); parseObjList->Add(po);
  po = new CMCondOr(); parseObjList->Add(po);
  
  po = new CMPlus(); parseObjList->Add(po);
  po = new CMMinus(); parseObjList->Add(po);
  po = new CMMul(); parseObjList->Add(po);
  po = new CMDiv(); parseObjList->Add(po);
  po = new CMPow(); parseObjList->Add(po);
  //po = new CMMod(); parseObjList->Add(po);

  po = new CMSqrt(); parseObjList->Add(po);
  po = new CMSqr(); parseObjList->Add(po);
  po = new CMSin(); parseObjList->Add(po);
  po = new CMCos(); parseObjList->Add(po);
  po = new CMTan(); parseObjList->Add(po);
  po = new CMSinh(); parseObjList->Add(po);
  po = new CMCosh(); parseObjList->Add(po);
  po = new CMTanh(); parseObjList->Add(po);
  po = new CMASin(); parseObjList->Add(po);
  po = new CMACos(); parseObjList->Add(po);
  po = new CMATan(); parseObjList->Add(po);
  po = new CMExp(); parseObjList->Add(po);
  po = new CMLn(); parseObjList->Add(po);
  po = new CMLog(); parseObjList->Add(po);
  po = new CMLog10(); parseObjList->Add(po);//$ RL 2011-6-10: log10(10.) = 1, log(e) = 1

  po = new CMNot(); parseObjList->Add(po);
  po = new CMFact(); parseObjList->Add(po);
  po = new CMAbs(); parseObjList->Add(po);
  po = new CMFAbs(); parseObjList->Add(po);//$ RL 2011-6-10: fabs added (same as abs)

  po = new CMRound(); parseObjList->Add(po);
  po = new CMFloor(); parseObjList->Add(po);
  po = new CMCeil(); parseObjList->Add(po);
  po = new CMHeaviside(); parseObjList->Add(po);
  po = new CMSgn(); parseObjList->Add(po);
  po = new CMPPD(); parseObjList->Add(po);

	po = new CMVAbs(); parseObjList->Add(po); //$ AD 2013-10-20: vector length
	po = new CMMatRows(); parseObjList->Add(po); //$ AD 2013-10-20: rows of a matrix
	po = new CMMatCols(); parseObjList->Add(po); //$ AD 2013-10-20: columns of a matrix
	po = new CMRandom(); parseObjList->Add(po); //$ AD 2013-10-20: random number
  po = new CMatrix(); parseObjList->Add(po);  
  po = new CMTranspose(); parseObjList->Add(po);  

  po = new CCDifferentiate(); parseObjList->Add(po);  
  //po = new CCInclude(); parseObjList->Add(po); //$JG2013-5-10: is available in EDCParser
  po = new CCInfo(); parseObjList->Add(po);

  //for (int i = 1; i <= parseObjList->Length(); i++)
  //{
  //  (*global_uo) << "Obj " << i << ": " << parseObjList->Get(i)->Name() << "\n";
  //}

}

void CParser::AddParseObj(CParseObj* cpo)
{
  parseObjList->Add(cpo);
}


void CParser::Write(mystr str)
{
  //(*outFile) << str.c_str(); //JG2013-04-29
  (*global_uo) << str.c_str();
}

void CParser::Init()
{
  GetVariableList()->Flush();
  SetParserErrorFlag(0);
}


void CParser::StorePos()
{
  storeactpos = actpos;
  storelines = lines;
  storeactlinepos = actlinepos;
}

void CParser::RestorePos()
{
  actpos = storeactpos;
  lines = storelines;
  actlinepos = storeactlinepos;
}

int CParser::ReadChar(char& ch)
{
  if (actpos < endpos)
  {
    ch = actstring[actpos++];

    actlinepos++;
    if (IsEOL(ch)) {lines++; actlinepos = 0;}
    
    return 1;
  } else {ch = 0; return 0;}
}

int CParser::PeekChar(char& ch)
{
  if (actpos < endpos)
  {
    ch = actstring[actpos];
    return 1;
  } else {ch = 0; return 0;}
}

void CParser::ReadSpaces()
{
  char ch;
  PeekChar(ch);
  while ((ch == ' ' || IsEOL(ch)) &&  !IsEOF(ch)) 
  {
    ReadChar(ch);
    PeekChar(ch);
  }
}

void CParser::ReadLineComment()
{
  char ch;
  ReadChar(ch);
  while (!IsEOL(ch) && !IsEOF(ch)) 
  {
    ReadChar(ch);
  }  
}

void CParser::ReadComment()
{
  char a = ' ';
  char b;
  ReadChar(b);
  while (!(a == '*' && b == '/') && !IsEOF(b)) 
  {
    a = b;
    ReadChar(b);
  }
  if (IsEOF(b))
  {
    Error("Syntax Error: you probably forgot to close a comment '/*' with '*/'");
  }
}

mystr CParser::ReadString()
{
  mystr retstr;
  char ch;
  ReadChar(ch);
  while (ch != '"' && !IsEOF(ch)) 
  {
    retstr += ch;
    ReadChar(ch);
  }
  if (IsEOF(ch))
  {
    Error("Syntax Error: you probably forgot to end a string");
  }
  return retstr;
}

int CParser::IsNumber(char ch)
{
  return ((int)ch <= (int)'9') && ((int)ch >= (int)'0');
}

int CParser::IsNumberChar(char ch)
{
  return ((int)ch <= (int)'9') && ((int)ch >= (int)'0') ||
          (ch == '.');
}

int CParser::IsCharacter(char ch)
{
  int ci  = (int) ch;
  return ((ci >= (int)'A') && (ci <= (int)'Z')) ||
         ((ci >= (int)'a') && (ci <= (int)'z')) ||
         (ch == 'Ä') || (ch == 'ä') || (ch == 'Ö') || (ch == 'ö') ||
         (ch == 'Ü') || (ch == 'ü') || (ch == 'é') || (ch == 'É') ||
         (ch == 'á') || (ch == 'Á') || (ch == '_') || (ch == '§') ||
         (ch == '$') || (ch == 'ß') || (ch == '.') || (ch == (char)39);
}

int CParser::IsSeparator(char ch)
{
  if ( ((ch >= 'A') && (ch <= 'Z')) ||
       ((ch >= 'a') && (ch <= 'z')) || 
       ((ch <= '9') && (ch >= '0'))  ) {return 0;}

  return (ch == '+') || (ch == '-') || (ch == '*') || (ch == '/') || (ch == '=') ||
         (ch == ':') || (ch == '(') || (ch == ')') || (ch == '[') || (ch == ']') ||
         (ch == '{') || (ch == '}') || (ch == '<') || (ch == '>') || (ch == '^') ||
         (ch == ' ') || (ch == ',') || (ch == ';') || (ch == '!') || (int)ch == 0 || 
         (ch == '%') || (ch == '&') || (ch == '|') || IsEOL(ch) || IsEOF(ch);
}

int CParser::IsEOL(char ch)
{
  return ch == (char)10 || ch == '\n';
}

int CParser::IsEOF(char ch)
{
  return ch == (char)EOF || ch == (char)0;
}

int CParser::IsNumber(const mystr& str)
{
  int pos = 0;
  int end = str.Length();
  while (pos < end && IsNumber(str[pos])) {pos++;}

  if (pos == end) {return 1;}
  else if (str[pos] == '.')
  {
    pos++;
    while (pos < end && IsNumber(str[pos])) {pos++;}

    if (pos == end) {return 1;}
  }
  if (str[pos] == 'e' || str[pos] == 'E')
  {
    pos++;
    if (str[pos] == '+' || str[pos] == '-')
    {
      pos++;
    }
    while (pos < end && IsNumber(str[pos])) {pos++;}

    if (pos == end) {return 1;}
  }
  
  Error("Syntax Error, \"" + str + "\" is not a valid Number");

  return 0;
}

int CParser::IsVariable(const mystr& str)
{
  if (str.Length() == 0 || !IsCharacter(str[0])) {return 0;}

  for (int i=1; i < str.Length(); i++)
  {
    char c = str[i];
    if (!IsCharacter(c) && !IsNumber(c)) {return 0;}
  }
  return 1;
}

double CParser::Number(const mystr& str)
{
  if (IsNumber(str))
  {
    mystr str1 = str;
    return atof(str1.c_str());
  }
  return 0.;
}

void CParser::WriteInfo()
{
  Write(mystr("Smart parser version 0.1 by Dipl.-Ing. Johannes Gerstmayr\n"));
  Write(mystr("NO. of variables = ")+mystr(GetVariableList()->Length())+mystr("\n"));
  Write(mystr("NO. of objects   = ")+mystr(parseObjList->Length())+mystr("\n"));
}

void CParser::Error(const mystr& string1)
{
  if (!ParserErrorFlag())
  {
    int p2 = actlinepos+1;
    Write(string1+mystr(" in line ")+mystr(lines)+mystr(", pos ")+mystr(p2)+mystr("\n"));
    //pout << flush; //JG2013-04-29: not necessary
  }
  MyError();
}

void CParser::Error(const mystr& string1, const mystr& string2)
{
  if (!ParserErrorFlag())
  {
    int p2 = actlinepos+1;
    Write(string1+mystr(" in line ")+mystr(lines)+mystr(", pos ")+mystr(p2)+
	  mystr(", unknown '")+string2+mystr("' !\n"));
    //pout << flush;//JG2013-04-29: not necessary
  }
  MyError();
}

//$ RL 2011-5-31:[
void CParsedFunction::CopyFrom(const CParsedFunction& e)
{
	//$JG2013-4-29: crashes, because of different DLLs
	//parameters = e.parameters; //$ RL 2011-5-31: only pointers to parameters are copied by TArray.

	////$JG2013-4-29: possible workaround, copy pointers:
	int parlen = e.parameters.Length();
	parameters.SetLen(parlen);
	for (int i=1; i<= parlen; i++)
	{
		parameters(i) = new CMNumber(e.parameters(i)->Double());
	}

	//$ PG 2013-8-13:[ Bugfix @{//$JG2013-4-29: possible workaround, copy pointers:}
	//$ ... but then the variablesList of the parser has to be adapted also,
	//$ otherwise modifications in parameters will not take effect!
	//$ the member actual_parser has been added to this class just for this purpose.
	//$ here we copy it, and redirect its variableList to the new 'TArray<CMNumbers*> parameters'.
	//$ also the method CVariable::SetMathObj(..) has been added to the class CVariable just for this purpose.
	actual_parser = 0;
	if (e.actual_parser)
	{
		actual_parser = e.actual_parser;

		for (int i=1; i<=actual_parser->GetVariableList()->Length(); i++)
		{
			CVariable* var = actual_parser->GetVariableList()->GetVariable(i);
			int idx = e.parameters.Find(static_cast <CMNumber*>(var->MathObj()));
			if (idx)
			{
				var->SetMathObj(parameters(idx));
			}
		}
	}
	//$ PG 2013-8-13:]

	int len = e.variableList.Length();
	if(len)
	{
		assert(0);                             //RL: not needed yet (variableList stored in mbs parser)
		variableList = e.variableList;         //RL: should be tested in case of variableList has nonzero length
	}
	
	mathobj = 0;
	if(e.mathobj)
	{
		mathobj = e.mathobj;		
	}

	parsed_function = e.parsed_function;
	parsed_function_parameter = e.parsed_function_parameter;

}
//$ RL 2011-5-31:]

void CParsedFunction::SetParsedFunction1D(CParser* parser, const mystr& mathstring, const mystr& parameter_name)
{
	//store strings for later use:
	parsed_function = mathstring;
	parsed_function_parameter = parameter_name;

	parameters.SetLen(1);   //$ JG 2011-5-31: //$ RL 2011-5-31: default value
	CMNumber* num = new CMNumber;
	parameters(1) = num;    //$ JG 2011-5-31: //$ RL 2011-5-31: pointer
	parameters(1)->Set(0.); //$ JG 2011-5-31: //$ RL 2011-5-31: initial value

	CVariable* var = new CVariable(parameter_name); 
	var->Assign(parameters(1));
	parser->GetVariableList()->Add(var);        // new variable for each parsed function
	mathobj = parser->String2MathObj(mathstring);
	var->SetName(mystr("@")+var->Name());       // rename variable to avoid conflicts, "@variable" is for use in functions only 

	actual_parser = parser;   //$ PG 2013-8-13: Bugfix @{//$JG2013-4-29: possible workaround, copy pointers:} ... but then the variablesList has to be adapted also, otherwise modifications in parameters will not take effect
}

//evaluate the parsed mathobj for the given parameter
double CParsedFunction::Evaluate(const double& parameter) const//$ RL 2011-5-31: const added.
{
	parameters(1)->Set(parameter);
	return mathobj->Double();
}
CMBSParser::CMBSParser():CParser()
{
	mbs = 0;
};

void CMBSParser::SetMBS(MBS * mbsi)
{
	mbs = mbsi;
}

void CMBSParser::Init()
{
	CParser::Init(); //variables list is flushed!
	local_edc = 0;
	local_edc2 = 0;
	//clear local variable lists / tree
}

void CMBSParser::SetLocalEDC(ElementDataContainer* local_edci)
{
	local_edc = local_edci;
}
void CMBSParser::SetLocalEDC2(ElementDataContainer* local_edci)
{
	local_edc2 = local_edci;
}

ElementDataContainer* CMBSParser::GetLocalEDC() const
{
	return local_edc;
}
ElementDataContainer* CMBSParser::GetLocalEDC2() const
{
	return local_edc2;
}

ElementDataContainer* CMBSParser::GetElementDataContainers(int i)
{
	if (i == 1) {return GetLocalEDC();}
	if (i == 2) {return GetLocalEDC2();}
	else if (i == 3) {return (mbs)->GetModelDataContainer();}
	else if (i == 4) {return (mbs)->GetMBS_EDC_Options();}
	else if (i == 5) {return (mbs)->GetMBS_EDC_Variables();}
	else {assert(0); return 0;}
}

CParseObj* CMBSParser::GetVariable(mystr& str)
{
	//TODO AD://check if variable exists in mathparser: e.g. pi or other
	//return GetVariableList()->Object(readstring);

	// extract the indixes from the full string, trim the textual identifier
	mystr comp1_str, comp2_str; 
//	int n_comp = str.Trim_And_Get_Indices(comp1,comp2);
	int n_comp = str.Trim_And_Get_Indices(comp1_str,comp2_str);
	
	// evaluate the indices strings to integers
	int comp1=0, comp2=0;
	if (n_comp==1 || n_comp==2)
	{
		CMathObj* comp1_MO = String2MathObj(comp1_str,0);
		comp1 = (int) comp1_MO->Double();
	}
	if (n_comp==2)
	{
		CMathObj* comp2_MO = String2MathObj(comp2_str,0);
		comp2 = (int) comp2_MO->Double();
	}


	if (GetVariableList()->IsObject(str))
	{
		return GetVariableList()->Object(str); //$AD variable in variable list has priority (use in Mathfunc)
	}

//	double val = 0;
	int found = 0;
	CVariable* var = new CVariable(str);

	//search variable in data containers, assign value to "val"
	int i=1;
	while (i <= NElementDataContainers() && ! found)
	{
		ElementDataContainer* edc = GetElementDataContainers(i);
		if (edc && edc->TreeFind(str))
		{
			found = 1;
//!AD: 2012-12-07: check if the variable is scalar, vector or matrix
			ElementData* ed = edc->TreeFind(str);
			if (ed->IsBool() || ed->IsInt() || ed->IsDouble())
			{
				if (n_comp==0) // scalar value, return a double
				{
					double val = edc->TreeGetDouble(str);
					CMNumber* num = new CMNumber(val);
					var->Assign(num);
				}
				else
					MyError(mystr("ERROR: access to component ")+mystr(comp1)+mystr(" of scalar object '")+str+mystr("' with [...] operator: can not access component of a scalar object!\n"));
			}
			else if (ed->IsVector() || ed->IsVectorXYZ())
			{
				if (n_comp==0) // entire Vector
				{
					double * vec; int len=0; int rv;
					rv = edc->TreeGetVector(str,&vec,len);

					CMatrix* vector = new CMatrix();
					vector->SetDimensions(new CMNumber(1),new CMNumber(len));  // 1 row, len columns
					vector->Generate();                                        // must generate to evaluate dimensions !!!
					for(int i=1; i<=len; i++)
					{
						vector->AddElem(new CMNumber(vec[i-1]));
					}
					var->Assign(vector);
				}
				else if (n_comp==1) // single vector component
				{
					double * vec; int len=0; int rv;
					rv = edc->TreeGetVector(str,&vec,len);
					if (rv)
					{
						if (comp1 > 0 && comp1 <= len)
						{
							CMNumber* component = new CMNumber(vec[comp1-1]);
							var->Assign(component);
						}
						else
						{
							MyError(mystr("ERROR: access to component ")+mystr(comp1)+mystr(" of object '")+str+mystr("' with [...] operator: out of range!\n"));
						}
					}
					else
					{
						MyError(mystr("ERROR: access to component of object '")+str+mystr("' with [...] operator not possible, because it is no vector!\n"));
					}
				}
				else
					MyError(mystr("ERROR: access to component ")+mystr(comp1)+mystr(',')+mystr(comp2)+mystr(" of vector object '")+str+mystr("' with [.,.] operator: can not access component of a vector object!\n"));
			}
			else if (ed->IsMatrix())
			{
				if (n_comp==0) // entire Matrix
				{
					double * vec; int rows=0; int cols=0; int rv;
					rv = edc->TreeGetMatrix(str,&vec,rows,cols);

					CMatrix* matrix = new CMatrix();
					matrix->SetDimensions(new CMNumber(rows),new CMNumber(cols));  // rows, cols
					matrix->Generate();                                            // must generate to evaluate dimensions !!!
					for(int i=1; i<=rows*cols; i++)
					{
						matrix->AddElem(new CMNumber(vec[i-1]));
					}
					var->Assign(matrix);
				}
				else if (n_comp==1) // single row/ single column ????, or set one component to zero
				{
					MyError(mystr("ERROR: access to components of object '")+str+mystr("' with [...] operator not possible: USE [] or [..,..] for Matrix access!\n"));
				}
				else if (n_comp==2) // single entry in the Matrix
				{
					double * vec; int rows=0; int cols=0; int rv;
					rv = edc->TreeGetMatrix(str,&vec,rows,cols);
					if (rv)
					{
						if (comp1 > 0 && comp1 <= rows && comp2 > 0 && comp2 <= cols)
						{
							CMNumber* component = new CMNumber(vec[(comp1-1)*cols+(comp2-1)]);
							var->Assign(component);
						}
						else
						{
							MyError(mystr("ERROR: access to component (")+mystr(comp1)+mystr(',')+mystr(comp2)+mystr(") of object '")+str+mystr("' with [..,..] operator: out of range!\n"));
						}
					}
					else
					{
						MyError(mystr("ERROR: access to component of object '")+str+mystr("' with [..,..] operator not possible, because it is no vector!\n"));
					}
				}
				else
					MyError(mystr("ERROR: access to components of object '")+str+mystr("' with [...,...] operator not possible: NOT IMPLEMENTED!\n"));
			}
		}



//!AD,JG: 2012-11-29 for right hand side access to vector or matrix component: [
		////	if (n_comp == 0)
		////	{
		////		val = edc->TreeGetDouble(str);
		////	}
		////	else if (n_comp == 1)
		////	{
		////		double * v;
		////		int len;
		////		int rv;
		////		rv = edc->TreeGetVector(str, &v, len);
		////		if (rv)
		////		{
		////			if (comp1 > 0 && comp1 <= len)
		////			{
		////				val = v[comp1-1];
		////			}
		////			else
		////			{
		////				MyError(mystr("ERROR: access to component ")+mystr(comp1)+mystr(" of object '")+str+mystr("' with [...] operator: out of range!\n"));
		////			}
		////		}
		////		else
		////		{
		////			MyError(mystr("ERROR: access to component of object '")+str+mystr("' with [...] operator not possible, because it is no vector!\n"));
		////		}
		////	}
		////	else if (n_comp == 2)
		////	{
		////		MyError(mystr("ERROR: access to components of object '")+str+mystr("' with [...,...] operator not possible: NOT IMPLEMENTED!\n"));
		////	}
		////}
//!AD,JG: ]
	i++;
	}
	if (!found)
	{
		//mbs->UO(UO_LVL_err) << "ERROR: did not find variable '" << str << "' in expression parser, assume '" << str << "' = 0 !!!!\n";

		//$ DR 2013-01-23 error handling added, in order to stop the parser as soon as an error is detected
		mbs->UO(UO_LVL_err) << "ERROR: did not find variable '" << str << "' in expression parser. Stopped parsing!\n";
		SetParserErrorFlag(1);
		return new CPVoid();
	}

	//now generate new number object for further processing with mathparser
	//CVariable* var = new CVariable(str);
	//CMNumber* num = new CMNumber(val);
	//var->Assign(num);
	//==> CMNumber can only be right-hand-side object (e.g.: a = 3) and not left-hand-side!
	//to change this, a new CMObject must be created, which has a number that points to the EDC Data Structure (e.g. CMNumberPtr)
	//WARNING: this could lead to problems in expression based MathFunctions!!!

	return var;
	//return CEDVariable(str);
}


double CMBSParser::ExpressionToDouble(mystr & data, int& error_flag, int line)
{
	error_flag = 0;
	CMathObj* mathobj;
	double val;
	SetParserErrorFlag(0);

	mathobj = String2MathObj(data, line);
	//check mathobj->EvaluableDouble() in future
	val = mathobj->Double();
	mbs->UO(UO_LVL_dbg1) << "Expression parser: '" << data << "' = " << val << "\n";

	if (ParserErrorFlag())
	{
		error_flag = 1;
	}

	//ed.SetDouble(val, elemname);
	return val;
}

//!AD to replace ExpressionToDouble
// returns Rows/Cols of result, fills corresponding *_var
int2 CMBSParser::EvaluateExpression(mystr & data, int& error_flag, double& d_var, Vector& v_var, Matrix& m_var, int line)
{
	error_flag = 0;
	CMathObj* mathobj;
	SetParserErrorFlag(0);

	mathobj = String2MathObj(data, line);

	int t1 = mathobj->MOType();  // can not use the Type to determine MOType Number/Vector/Scalar
	int2 rc;	mathobj->Dim(rc);  // possible to use Dimensions
	int t2 = mathobj->ParseObjType(); // whatever this is...

	//$ DR 2013-01-23 error handling added, in order to stop the parser as soon as an error is detected
	if(ParserErrorFlag())
	{
		error_flag = 1;
		return rc;
	}

  if(rc(1)==1 && rc(2)==1) //mathobj->MOType() == TMONumber)
	{
		// scalar
		d_var = mathobj->Double();
		//mbs->UO(UO_LVL_dbg1) << "Expression parser: '" << data << "' = " << d_var << " (scalar)\n"; //JG2012-04-29 : this line will slow down HOTINT for parsed functions extremely!
	}
	else if(rc(1)==1 || rc(2)==1) //mathobj->MOType() == TMOMatrix)
	{
		// vector
		v_var.SetLen(rc(1)*rc(2));
		v_var.SetAll(0.);
		for (int r=1; r <=rc(1); r++)
		{
			for (int c=1; c<=rc(2); c++) 
			{
				v_var.Elem( (r-1)*rc(2)+c) = mathobj->MatDouble(r,c);
			}
		}
	}
	else  
	{
		//true Matrix
		m_var.SetSize(rc(1),rc(2));
		m_var.SetAll(0.);
		for (int r=1; r <=rc(1); r++)
		{
			for (int c=1; c<=rc(2); c++) 
			{
				m_var.Elem(r,c) = mathobj->MatDouble(r,c);
			}
		}
	}
	return rc;
}