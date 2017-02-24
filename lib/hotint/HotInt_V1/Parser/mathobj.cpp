
//#**************************************************************
//#
//# filename:             mathobj.cpp
//#
//# autor:                Gerstmayr Johannes
//#
//# generated:            23.11.98
//# last change:          23.11.98
//# description:          implementation class for mathematical objects
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
#include "parseobj.h"
#include "mathobj.h"
#include "parser.h"

extern int ParserErrorFlag();
extern UserOutputInterface * global_uo;

double fact(double n)
{
  double x=1;
  for (int i = 2; i <= n; i++)
  {
    x*=(double)i;
  }
  return x;
}

CMathObj::CMathObj()
{
  parent = NULL;
}

CParseObj* CMathObj::GetCopy() 
{
  return new CMathObj();
}

mystr CMathObj::Name() 
{
  return mystr("@MathObj");
}

CMathObj* CMathObj::CopyMO()
{
  CMathObj* co= (CMathObj*)GetCopy();
  co->parent = NULL;
  
  int i;
  for (i = 0; i < GetNOArg(); i++)
  {
    co->SetArg(Arg(i)->CopyMO(),i);    
  }
  return co;
}

void CMathObj::DeleteMO()
{
  int i;
  for (i = 0; i < GetNOArg(); i++)
  {
    Arg(i)->DeleteMO();
    delete Arg(i);
  }
}

CMathObj* CMathObj::AddBase()
{
  if (Parent() == NULL || Parent()->MOType() != TMOBase)
    {
      (*global_uo) << "addbase done\n";
      SetParent(new CBaseObj());
    }
  return this;
}

TParseObj CMathObj::ParseObjType() 
{
  return TPOMathObj;
}

CMathObj* CMathObj::Arg(int /*arg*/) 
{
  return NULL;
}

int CMathObj::GetNOArg() 
{
  return 0;
}

void CMathObj::SetArg(CMathObj* /*arg*/, int /*argnum*/) 
{
};

void CMathObj::AddArg(CMathObj* arg) 
{
  SetArg(arg);
};


void CMathObj::ReplaceArg(CMathObj* oldarg, CMathObj* replacearg)
{
  for (int i = 0; i < GetNOArg(); i++)
  {
    if (Arg(i) == oldarg) {SetArg(replacearg, i);}
  }
} 

CMathObj* CMathObj::Parent() 
{
  if (parent == NULL) 
  {
    MyError(mystr("parent of object '")+Name()+mystr("' is NULL !"));
  }
  return parent;
}

void CMathObj::SetParent(CMathObj* par) 
{
  parent = par;
}

TMathObj CMathObj::MOType() 
{
  return TMOMathObj;
}

void CMathObj::Eliminate(CVariable* var) 
{
  for (int i = 0; i < GetNOArg(); i++)
  {
    Arg(i)->Eliminate(var);
  }
};

double CMathObj::Double()
{
  return 0;
}

int CMathObj::EvaluableDouble()
{
  for (int i = 0; i < GetNOArg(); i++)
  {
    if (!Arg(i)->EvaluableDouble()) {return 0;}
  }
  return 1;
}

/*
void CMathObj::MatExecute()
{
  for (int i = 0; i < GetNOArg(); i++)
  {
    Arg(i)->MatExecute();
  }
}
*/

CMathObj* CMathObj::Differentiate(CMathObj* var)
{
  MyWarning(mystr("Warning: differentiation for mathematical object '")+Name()+mystr("' is not defined!"));
  return CopyMO();
}

CMathObj* CMathObj::StdSimplify()
{
  MyWarning(mystr("Warning: stdsimplify for mathematical object '")+Name()+mystr("' is not defined!"));
  return CopyMO();
}

//void CMathObj::Elem(int /*row*/, int /*col*/, CMathObj* /*replace*/)
/*
{
  MyError(mystr("Execution Error: cannot modify matrix-element of '")+
          String()+mystr("'"));
}
*/
CMathObj* CMathObj::Get(int /*row*/, int /*col*/)
{
  MyError(mystr("Execution Error: cannot get matrix-element of '")+
          String()+mystr("'"));
  return new CMEnd();
}

double CMathObj::MatDouble(int /*row*/, int /*col*/)
{
  MyError(mystr("Execution Error: cannot calculate matrix-element of '")+
          String()+mystr("'"));
  return 0;
}

void CMathObj::Dim(int2& rc)
{
  rc.Get(1) =1;
  rc.Get(2) =1;
/*
  MyError(mystr("Execution Error: dim is not defined for '")+
          Name()+mystr("'"));*/
}
//++++++++++++++++++++++Operator++++++++++++++++++++++

COperator::COperator() : CMathObj()
{
  arg0 = NULL; arg1 = NULL;
}

mystr COperator::Name() 
{
  return mystr("@Operator");
}

void COperator::SetArg(CMathObj* arg, int argnum)
{
  arg->SetParent(this);
  if (argnum == 0) {arg0 = arg;}
  else if (argnum == 1) {arg1 = arg;}
  else {SysError("COperator_SetArg");}
}

CMathObj* COperator::Arg(int arg)
{
  if (arg == 0) {return arg0;}
  else if (arg == 1) {return arg1;}
  else {SysError("COperator_Arg"); return arg0;}
}

double COperator::Double() 
{
  return 0;
}

mystr COperator::String()
{
  mystr retstr;

  TMathPriority p = Priority();
  
  if (Arg(0)->MOType() == TMOOperator &&
      (Arg(0)->Priority() < p ||
       (Arg(0)->Priority() == p && p == TMPEqual)))
  {
    retstr += mystr("(") + Arg(0)->String() + mystr(")");
  } else
  {
    retstr += Arg(0)->String();
  }  

  retstr += Name();
  
  if ((Arg(1)->MOType() == TMOOperator &&
      (Arg(1)->Priority() < p || 
       (Arg(1)->Priority() == p && 
        (!Commutative() || p == TMPEqual) )) ||
      (Arg(1)->MOType() == TMONumber && Arg(1)->Double()<0 && p>=TMPPlus)) )
      // || (Arg(1)->MOType() == TMOUnOperator && Arg(1)->Priority() <= p))
  {
    retstr += mystr("(") + Arg(1)->String() + mystr(")");
  } else
  {
    retstr += Arg(1)->String();
  }  

  return retstr;
}

mystr COperator::StringAll()
{
  mystr retstr;

  TMathPriority p = Priority();
  
  if (Arg(0)->MOType() == TMOOperator &&
      (Arg(0)->Priority() < p ||
       (Arg(0)->Priority() == p && p == TMPEqual)))
  {
    retstr += mystr("(") + Arg(0)->StringAll() + mystr(")");
  } else
  {
    retstr += Arg(0)->StringAll();
  }  

  retstr += Name();
  
  if ((Arg(1)->MOType() == TMOOperator &&
      (Arg(1)->Priority() < p || 
       (Arg(1)->Priority() == p && 
        (!Commutative() || p == TMPEqual) )) ||
      (Arg(1)->MOType() == TMONumber && Arg(1)->Double()<0 && p>=TMPPlus)) )
      //|| (Arg(1)->MOType() == TMOUnOperator && Arg(1)->Priority() <= p))
  {
    retstr += mystr("(") + Arg(1)->StringAll() + mystr(")");
  } else
  {
    retstr += Arg(1)->StringAll();
  }  

  return retstr;
}

void COperator::AddArg(CMathObj* arg) 
{
  SetArg(arg,1);
}

int COperator::GetNOArg() {return 2;}

TMathObj COperator::MOType() 
{
  return TMOOperator;
}

TMathPriority COperator::Priority() 
{
  return TMPLow;
}

CMathObj* COperator::StdSimplify()
{
  CMathObj* co = (CMathObj*)GetCopy();

  co->SetArg(arg0->StdSimplify(),0);
  co->SetArg(arg1->StdSimplify(),1);

  return co;
}

//+++++++++++++++  CMPlus ++++++++++++++++++++
void CMPlus::Dim(int2& rc)
{
  int2 r0,r1;
  Arg(0)->Dim(r0);
  Arg(1)->Dim(r1);

  rc.Get(1)=r0.Get(1);
  rc.Get(2)=r0.Get(2);
  
  if (r0.Get(1) != r1.Get(1) || r0.Get(2) != r1.Get(2)) 
  {
    MyError(mystr("Execution Error: incompatible dim in matrix addition '")+
            String()+mystr("'"));
  }
}

CMathObj* CMPlus::Differentiate(CMathObj* var)
{
  CMPlus* cp = new CMPlus();
  
  cp->SetArg(arg0->Differentiate(var),0);
  cp->SetArg(arg1->Differentiate(var),1);
  return cp;  
}

CMathObj* CMPlus::StdSimplify()
{
  CMathObj* a0ss = arg0->StdSimplify();
  CMathObj* a1ss = arg1->StdSimplify();

  int a0e = a0ss->EvaluableDouble();
  int a1e = a1ss->EvaluableDouble();

  if (a0e && a1e) {return new CMNumber(a0ss->Double()+a1ss->Double());}
  else if (a0e)
    {
      double d = a0ss->Double();
      if (d == 0)
	{
	  return a1ss;
	}
      else
	{
	  CMPlus* cp = new CMPlus();
	  cp->SetArg(new CMNumber(d),0);
	  cp->SetArg(a1ss,1);
	  return cp;  
	}
    }
  else if (a1e)
    {
      double d = a1ss->Double();
      if (d == 0)
	{
	  return a0ss;
	}
      else
	{
	  CMPlus* cp = new CMPlus();
	  cp->SetArg(a0ss,0);
	  cp->SetArg(new CMNumber(d),1);
	  return cp;  
	}
    }
  else
    {  
      CMPlus* cp = new CMPlus();
      
      cp->SetArg(a0ss,0);
      cp->SetArg(a1ss,1);
      return cp;  
    }
}

//+++++++++++++++  CMMinus ++++++++++++++++++++
void CMMinus::Dim(int2& rc)
{
  int2 r0,r1;
  Arg(0)->Dim(r0);
  Arg(1)->Dim(r1);

  rc.Get(1)=r0.Get(1);
  rc.Get(2)=r0.Get(2);
  
  if (r0.Get(1) != r1.Get(1) || r0.Get(2) != r1.Get(2)) 
  {
    MyError(mystr("Execution Error: incompatible dim in matrix subtraction '")+
            String()+mystr("'"));
  }
}

CMathObj* CMMinus::Differentiate(CMathObj* var)
{
  CMMinus* cp = new CMMinus();
  
  cp->SetArg(arg0->Differentiate(var),0);
  cp->SetArg(arg1->Differentiate(var),1);
  return cp;  
}

CMathObj* CMMinus::StdSimplify()
{
  CMathObj* a0ss = arg0->StdSimplify();
  CMathObj* a1ss = arg1->StdSimplify();

  int a0e = a0ss->EvaluableDouble();
  int a1e = a1ss->EvaluableDouble();

  if (a0e && a1e) {return new CMNumber(a0ss->Double()-a1ss->Double());}
  else if (a0e)
    {
      double d = a0ss->Double();
      if (d == 0)
	{
	  CMUnMinus* cmu = new CMUnMinus();
	  cmu->SetArg(a1ss);
	  return cmu;
	}
      else
	{
	  CMMinus* cp = new CMMinus();
	  cp->SetArg(new CMNumber(d),0);
	  cp->SetArg(a1ss,1);
	  return cp;  
	}
    }
  else if (a1e)
    {
      double d = a1ss->Double();
      if (d == 0)
	{
	  return a0ss;
	}
      else
	{
	  CMMinus* cp = new CMMinus();
	  cp->SetArg(a0ss,0);
	  cp->SetArg(new CMNumber(d),1);
	  return cp;  
	}
    }
  else
    {  
      CMMinus* cp = new CMMinus();
      
      cp->SetArg(a0ss,0);
      cp->SetArg(a1ss,1);
      return cp;  
    }
}
//+++++++++++++++  CMMUL ++++++++++++++++++++
//return matrix element which is the result of a multiplication of two matrices
double CMMul::MatDouble(int row, int col)
{
  int2 rc, r0;
  Dim(rc);
  double val = 0;
  Arg(0)->Dim(r0);
  
  if (row < 1 || row > rc.Get(1) || col < 1 || col > rc.Get(2) || ParserErrorFlag())
  {
    MyError(mystr("Execution Error: tried to get a matrix-element out of dimension'")+
            String()+mystr("'"));
    return 0;
  }
  
  for (int i = 1; i <= r0.Get(2); i++)
  {
    val += Arg(0)->MatDouble(row,i) * Arg(1)->MatDouble(i,col);
  }
  return val;
}

void CMMul::Dim(int2& rc)
{
  int2 r0,r1;
  Arg(0)->Dim(r0);
  Arg(1)->Dim(r1);

  rc.Get(1)=r0.Get(1);
  rc.Get(2)=r1.Get(2);
  
  if (r0.Get(2) != r1.Get(1)) 
  {
    MyError(mystr("Execution Error: incompatible dimensions in matrix multiplication '")+
            String()+mystr("'"));
  }
}


CMathObj* CMMul::Differentiate(CMathObj* var)
{
  //product rule:
  CMMul* cm0 = new CMMul();
  CMMul* cm1 = new CMMul();
  CMPlus* cp = new CMPlus();
  cp->SetArg(cm0,0);
  cp->SetArg(cm1,1);

  cm0->SetArg(arg0->Differentiate(var),0);
  cm0->SetArg(arg1->CopyMO(),1);
  cm1->SetArg(arg0->CopyMO(),0);
  cm1->SetArg(arg1->Differentiate(var),1);

  return cp;  
}

CMathObj* CMMul::StdSimplify()
{
  CMathObj* a0ss = arg0->StdSimplify();
  CMathObj* a1ss = arg1->StdSimplify();

  int a0e = a0ss->EvaluableDouble();
  int a1e = a1ss->EvaluableDouble();

  if (a0e && a1e) {return new CMNumber(a0ss->Double()*a1ss->Double());}
  else if (a0e)
    {
      double d = a0ss->Double();
      if (d == 0) {return new CMNumber(0);}
      else if (d == 1) {return a1ss;}
      else
	{
	  CMMul* cp = new CMMul();
	  cp->SetArg(new CMNumber(d),0);
	  cp->SetArg(a1ss,1);
	  return cp;  
	}
    }
  else if (a1e)
    {
      double d = a1ss->Double();
      if (d == 0) {return new CMNumber(0);}
      else if (d == 1) {return a0ss;}
      else
	{
	  CMMul* cp = new CMMul();
	  cp->SetArg(a0ss,0);
	  cp->SetArg(new CMNumber(d),1);
	  return cp;  
	}
    }
  else
    {  
      CMMul* cp = new CMMul();
      
      cp->SetArg(a0ss,0);
      cp->SetArg(a1ss,1);
      return cp;  
    }
}
CMathObj* CMDiv::Differentiate(CMathObj* var)
{
  //product rule:
  CMMul* cm0 = new CMMul();
  CMMul* cm1 = new CMMul();
  CMPlus* cp = new CMPlus();
  CMDiv* cd = new CMDiv();
  CMPow* arg0sqr = new CMPow();
  

  cp->SetArg(cm0,0);
  cp->SetArg(cm1,1);
  arg0sqr->SetArg(arg1->CopyMO(),0);
  arg0sqr->SetArg(new CMNumber(2),1);
  cd->SetArg(cp,0);
  cd->SetArg(arg0sqr,1);

  cm0->SetArg(arg0->Differentiate(var),0);
  cm0->SetArg(arg1->CopyMO(),1);
  cm1->SetArg(arg0->CopyMO(),0);
  cm1->SetArg(arg1->Differentiate(var),1);

  return cd;  
}

CMathObj* CMDiv::StdSimplify()
{

  int a0e = arg0->EvaluableDouble();
  int a1e = arg1->EvaluableDouble();

  if (a0e && a1e) {return new CMNumber(Double());}
  else if (a0e)
    {
      double d = arg0->Double();
      if (d == 0) {return new CMNumber(0);}
      else
	{
	  CMDiv* cp = new CMDiv();
	  cp->SetArg(new CMNumber(d),0);
	  cp->SetArg(arg1->StdSimplify(),1);
	  return cp;  
	}
    }
  else if (a1e)
    {
      double d = arg1->Double();
      if (d == 1) {return arg0->StdSimplify();}
      else
	{
	  CMDiv* cp = new CMDiv();
	  cp->SetArg(arg0->StdSimplify(),0);
	  cp->SetArg(new CMNumber(d),1);
	  return cp;  
	}
    }
  else
    {  
      CMDiv* cp = new CMDiv();
      
      cp->SetArg(arg0->StdSimplify(),0);
      cp->SetArg(arg1->StdSimplify(),1);
      return cp;  
    }
}
//+++++++++++++++  CMPOW ++++++++++++++++++++

double CMPow::Double()
{
  return pow(Arg(0)->Double(), Arg(1)->Double());
}

CMathObj* CMPow::Differentiate(CMathObj* var)
{
  if (arg1->MOType() == TMONumber)
    {
      CMMul* cm0 = new CMMul();
      CMMul* cm1 = new CMMul();
      CMPow* cp = new CMPow();
      double n = arg1->Double();

      cp->SetArg(arg0->CopyMO(),0);
      cp->SetArg(new CMNumber(n-1.),1);
      cm0->SetArg(cp,0);
      cm0->SetArg(new CMNumber(n),1);
      cm1->SetArg(cm0,0);
      cm1->SetArg(arg0->Differentiate(var),1);
      
      return cm1;
    } else
      {
	CMMul* cm0 = new CMMul();
	CMMul* cm1 = new CMMul();
	CMPow* cp1 = new CMPow();
	CMMinus* cmi = new CMMinus();

	CMMul* cm2 = new CMMul();
	CMMul* cm3 = new CMMul();
	CMPow* cp2 = new CMPow();
	CMLn* ln1 = new CMLn();

	CMPlus* cpl = new CMPlus();
	
	cmi->SetArg(arg1->CopyMO(),0);
	cmi->SetArg(new CMNumber(1),1);
	cp1->SetArg(arg0->CopyMO(),0);
	cp1->SetArg(cmi,1);
	cm0->SetArg(cp1,0);
	cm0->SetArg(arg1->CopyMO(),1);
	cm1->SetArg(cm0,0);
	cm1->SetArg(arg0->Differentiate(var),1);

	cp2->SetArg(arg0->CopyMO(),0);
	cp2->SetArg(arg1->CopyMO(),1);
	cm2->SetArg(cp2,0);
	cm2->SetArg(arg1->Differentiate(var),1);
	ln1->SetArg(arg0->CopyMO());
	cm3->SetArg(cm2,0);
	cm3->SetArg(ln1,1);

	cpl->SetArg(cm3,0);
	cpl->SetArg(cm1,1);
      
	return cpl;
	
      }
}

CMathObj* CMPow::StdSimplify()
{
  int a0e = arg0->EvaluableDouble();
  int a1e = arg1->EvaluableDouble();

  if (a0e && a1e) {return new CMNumber(Double());}
  else if (a0e && arg0->Double() == 0) {return new CMNumber(0);}
  else if (a1e)
    {
      double d = arg1->Double();
      if (d == 0) {return new CMNumber(1);}
      else if (d == 1) {return arg0->StdSimplify();}
      else
	{
	  CMPow* cp = new CMPow();
	  cp->SetArg(arg0->StdSimplify(),0);
	  cp->SetArg(new CMNumber(d),1);
	  return cp;  
	}
    }
  else
    {  
      CMPow* cp = new CMPow();
      
      cp->SetArg(arg0->StdSimplify(),0);
      cp->SetArg(arg1->StdSimplify(),1);
      return cp;  
    }
}
//$JG2013-5-3: removed modulo because conflict to comment sign '%'
////+++++++++++++++  CMMod ++++++++++++++++++++
//
//double CMMod::Double()
//{
//  double a = Arg(0)->Double();
//  double b = Arg(1)->Double();
//  if ((int)a == a && (int)b == b) {return (int)a % (int)b;}
//  else 
//  {
//    MyError(mystr("Calculation Error: Mod accepts only whole numbers '")+
//            String()+mystr("' !"));
//    return 0;
//  }
//}

//+++++++++++++++  CMNUMBER  ++++++++++++++++++++

mystr CMNumber::Name() 
{
  return mystr(n);
}

CMathObj* CMNumber::Differentiate(CMathObj* var)
{
  return new CMNumber(0);
}

CMathObj* CMNumber::StdSimplify()
{
  return new CMNumber(n);
}

//+++++++++++++++CMVARIABLE++++++++++++++++++++

CMVariable::CMVariable() : COperand()
{
  var = NULL;
};

CMVariable::CMVariable(CVariable* vari) : COperand()
{
  var = vari;
}

CParseObj* CMVariable::GetCopy() 
{
  return new CMVariable(var);
}

mystr CMVariable::Name() 
{
  return var->Name();
}

TMathObj CMVariable::MOType() 
{
  return TMOVariable;
}

double CMVariable::Double() 
{
  if (var->MathObj() != NULL)
  {
    return var->MathObj()->Double();
  } else
  {
    MyError(mystr("Evaluation Error: cannot calculate '")+Name()+mystr("'"));
    return 0;
  }
}

int CMVariable::EvaluableDouble() 
{
  if (var->MathObj() != NULL)
  {
    return var->MathObj()->EvaluableDouble();
  } else {return 0;}
  
}

void CMVariable::Eliminate(CVariable* vari) 
{
  if (var == vari)
  {                     
    Parent()->ReplaceArg(this, var->MathObj());
  }
};

CVariable* CMVariable::Variable() 
{
  return var;
}

mystr CMVariable::StringAll()
{
  if (var->MathObj() == NULL) {return Name();}
  else {return var->MathObj()->StringAll();}
}

CMathObj* CMVariable::Differentiate(CMathObj* var)
{
  if (var->MOType() == TMOVariable)
    {
      if (((CMVariable*)var)->Variable() == Variable())
	{
	  return new CMNumber(1);
	}
    }
  if (Variable()->MathObj() == NULL) {return new CMNumber(0);}
  else
    {
      return Variable()->MathObj()->Differentiate(var);
    }
}

CMathObj* CMVariable::StdSimplify()
{
  if (var->MathObj() == NULL) {return CopyMO();}
  else
    {
      return var->MathObj()->StdSimplify();
    }
}

void CMVariable::Elem(int row, int col, CMathObj* replace)
{
  if (var->MathObj()->MOType() == TMOMatrix)
  {
    ((CMatrix*)var->MathObj())->Elem(row,col,replace);
  } else
  {
    MyError(mystr("Execution Error: tried to set a matrix-element in '")+
            var->MathObj()->String()+mystr("'"));    
  }
}

CMathObj* CMVariable::Get(int row, int col)
{
  return var->MathObj()->Get(row,col);
}

double CMVariable::MatDouble(int row, int col)
{
  return var->MathObj()->MatDouble(row,col);
}

void CMVariable::Dim(int2& rc)
{
  var->MathObj()->Dim(rc);
}

//+++++++++++++++++++++++  CFunction ++++++++++++++++++

CFunction::CFunction() : CMathObj() 
{
  arg0 = NULL;
};

mystr CFunction::Name() 
{
  return mystr("@Function");
}

CMathObj* CFunction::Arg(int /*arg*/) 
{
  return arg0;
};

void CFunction::SetArg(CMathObj* arg, int /*argnum*/)
{
  arg->SetParent(this);
  arg0 = arg;
};

int CFunction::GetNOArg() 
{
  return 1;
}

TMathObj CFunction::MOType() 
{
  return TMOFunction;
}

/*
void CFunction::ReplaceArg(CMathObj* oldarg, CMathObj* replacearg)
{
  if (arg0 == oldarg) {SetArg(replacearg);}
  // *delete oldarg;
}*/

double CFunction::Double() 
{
  return arg0->Double();
}

mystr CFunction::String()
{
  if (arg0 != NULL)
  {
    return Name()+(mystr)"("+arg0->String()+(mystr)")";
  } else
  {
    return Name()+(mystr)"("+(mystr)")";
  }
  
}

mystr CFunction::StringAll()
{
  if (arg0 != NULL)
  {
    return Name()+(mystr)"("+arg0->StringAll()+(mystr)")";
  } else
  {
    return Name()+(mystr)"("+(mystr)")";
  }
  
}

CMathObj* CFunction::StdSimplify()
{
  CMathObj* cf = (CMathObj*)GetCopy();

  cf->SetArg(arg0->StdSimplify(),0);

  return cf;
}


//++++++++++++++++++  CMUnOperators +++++++++++++++++++
mystr CMUnPlus::String()
{
    //plus is not needed
  return mystr("+") + Arg(0)->String();
}

mystr CMUnPlus::StringAll()
{
    //plus is not needed
  return mystr("+") + Arg(0)->StringAll();
}

void CMUnPlus::Dim(int2& rc)
{
  Arg(0)->Dim(rc);
}

mystr CMUnMinus::String()
{
  mystr str;
  if (Arg(0)->MOType() == TMOOperator &&
      Arg(0)->Priority() <= TMPPlus)
  {
    str = Name() + mystr("(") + Arg(0)->String() + mystr(")");
  } else
  {
    str = Name() + Arg(0)->String();
  }  
  
  if (Parent()->MOType() == TMOOperator &&
      Parent()->Priority() >= TMPPlus) 
    {
     return mystr("(")+str+mystr(")");
     }
  return str;

}

mystr CMUnMinus::StringAll()
{
  mystr str;
  if (Arg(0)->MOType() == TMOOperator &&
      Arg(0)->Priority() <= TMPPlus)
  {
    str = Name() + mystr("(") + Arg(0)->StringAll() + mystr(")");
  } else
  {
    str = Name() + Arg(0)->StringAll();
  }  
   
  if (Parent()->MOType() == TMOOperator &&
      Parent()->Priority() >= TMPPlus) 
    {
      str = mystr("(")+str+mystr(")");
      }
  return str;
}

CMathObj* CMUnMinus::Differentiate(CMathObj* var)
{
  CMUnMinus* cp = new CMUnMinus();
  
  cp->SetArg(arg0->Differentiate(var),0);
  return cp;  
}


void CMUnMinus::Dim(int2& rc)
{
  Arg(0)->Dim(rc);
}

//++++++++++++++++++  CMFUNCTIONS +++++++++++++++++++

double CMSqrt::Double()
{
  double a = Arg(0)->Double();
  if (a > 0) {return sqrt(a);}
	else if(a == 0){return 0.;}//$ RL 2011-6-1: line added.
  else 
  {
    MyError(mystr("Calculation Error: negative argument of sqrt is impossible: '")+
            String()+mystr("' !"));
    return 0;
  }
}

CMathObj* CMSqrt::Differentiate(CMathObj* var)
{
  CMPow* po = new CMPow();
  CMMul* mu1 = new CMMul();
  CMMul* mu2 = new CMMul();

  
  po->SetArg(arg0->CopyMO(),0);
  po->SetArg(new CMNumber(-0.5),1);
  mu1->SetArg(new CMNumber(0.5),0);
  mu1->SetArg(po,1);
  mu2->SetArg(mu1,0);
  mu2->SetArg(arg0->Differentiate(var),1);

  return mu2;  
}


double CMSqr::Double()
{
  double a = Arg(0)->Double();
  return a*a;
}

CMathObj* CMSqr::Differentiate(CMathObj* var)
{
  CMMul* cm1 = new CMMul();
  CMMul* cm2 = new CMMul();
  
  cm1->SetArg(new CMNumber(2),0);
  cm1->SetArg(arg0->CopyMO(),1);
  cm2->SetArg(cm1,0);
  cm2->SetArg(arg0->Differentiate(var),1);

  return cm2;  
}

double CMSin::Double()
{
  return sin(Arg(0)->Double());
}

CMathObj* CMSin::Differentiate(CMathObj* var)
{
  CMCos* co = new CMCos();
  CMMul* cm = new CMMul();
  

  co->SetArg(arg0->CopyMO(),0);
  cm->SetArg(co,0);
  cm->SetArg(arg0->Differentiate(var),1);

  return cm;  
}

double CMCos::Double()
{
  return cos(Arg(0)->Double());
}

CMathObj* CMCos::Differentiate(CMathObj* var)
{
  CMUnMinus* umi = new CMUnMinus();
  CMSin* si = new CMSin();
  CMMul* cm = new CMMul();
  

  si->SetArg(arg0->CopyMO(),0);
  umi->SetArg(si,0);
  cm->SetArg(si,0);
  cm->SetArg(arg0->Differentiate(var),1);

  return cm;  
}

double CMTan::Double()
{
  return tan(Arg(0)->Double());
}

CMathObj* CMTan::Differentiate(CMathObj* var)
{
  CMPow* po = new CMPow();
  CMCos* co = new CMCos();
  CMMul* cm = new CMMul();
  

  co->SetArg(arg0->CopyMO(),0);
  po->SetArg(co,0);
  po->SetArg(new CMNumber(-2),1);
  cm->SetArg(po,0);
  cm->SetArg(arg0->Differentiate(var),1);

  return cm;  
}

double CMSinh::Double()
{
  return sinh(Arg(0)->Double());
}

double CMCosh::Double()
{
  return cosh(Arg(0)->Double());
}

double CMTanh::Double()
{
  return tanh(Arg(0)->Double());
}

double CMASin::Double()
{
  double a = Arg(0)->Double();
  if (fabs(a) <= 1) {return asin(a);}
  else 
  {
    MyError(mystr("Calculation Error: asin accepts values from -1 to +1 '")+
            String()+mystr("' !"));
    return 0;
  }
}

double CMACos::Double()
{
  double a = Arg(0)->Double();
  if (fabs(a) <= 1) {return acos(a);}
  else 
  {
    MyError(mystr("Calculation Error: acos accepts values from -1 to +1 '")+
            String()+mystr("' !"));
    return 0;
  }
}

double CMATan::Double()
{
  double a = Arg(0)->Double();
	return atan(a);
	//$ AD,CZ: atan does not have a argument restriction
  ////if (fabs(a) <= MY_PI/2) {return atan(a);}
  ////else 
  ////{
  ////  MyError(mystr("Calculation Error: atan accepts values from -pi/2 to +pi/2 '")+
  ////          String()+mystr("' !"));
  ////  return 0;
  ////}
}

double CMExp::Double()
{
  return exp(Arg(0)->Double());
}

CMathObj* CMExp::Differentiate(CMathObj* var)
{
  CMExp* ex = new CMExp();
  CMMul* cm = new CMMul();
  

  ex->SetArg(arg0->CopyMO(),0);
  cm->SetArg(ex,0);
  cm->SetArg(arg0->Differentiate(var),1);

  return cm;  
}

double CMLn::Double()
{
  double a = Arg(0)->Double();
  if (a > 0) {return log(a);}
  else 
  {
    MyError(mystr("Calculation Error: ln accepts values from >0 to +infinity '")+
            String()+mystr("' !"));
    return 0;
  }
}

CMathObj* CMLn::Differentiate(CMathObj* var)
{
  CMDiv* di = new CMDiv();
  
  di->SetArg(arg0->Differentiate(var),0);
  di->SetArg(arg0->CopyMO(),1);

  return di;  
}

double CMLog::Double()
{
  double a = Arg(0)->Double();
  if (a > 0) {return log(a);}
  else 
  {
    MyError(mystr("Calculation Error: ln accepts values from >0 to +infinity '")+
            String()+mystr("' !"));
    return 0;
  }
}

CMathObj* CMLog::Differentiate(CMathObj* var)
{
  CMDiv* di = new CMDiv();
  
  di->SetArg(arg0->Differentiate(var),0);
  di->SetArg(arg0->CopyMO(),1);

  return di;  
}

double CMLog10::Double()
{
  double a = Arg(0)->Double();
  if (a > 0) {return log10(a);}
  else 
  {
    MyError(mystr("Calculation Error: log10 accepts values from >0 to +infinity '")+
            String()+mystr("' !"));
    return 0;
  }
}

CMathObj* CMLog10::Differentiate(CMathObj* var)
{
  CMDiv* di = new CMDiv();
  CMMul* mu = new CMMul();
  
  mu->SetArg(new CMNumber(log(10.)),0);
  mu->SetArg(arg0->CopyMO(),1);
  di->SetArg(arg0->Differentiate(var),0);
  di->SetArg(mu,1);

  return di;  
}

double CMFact::Double()
{
  double a = Arg(0)->Double();
  if (a >= 0 && a < 170 && ((int)a == a)) {return fact(a);}
  else 
  {
    MyError(mystr("Calculation Error: fact accepts whole numbers from 0 to 170 '")+
            String()+mystr("' !"));
    return 0;
  }
}

double CMAbs::Double()
{
  return fabs(Arg(0)->Double());
}

//$ RL 2011-6-10: [ fabs added (same as abs)
double CMFAbs::Double()
{
  return fabs(Arg(0)->Double());
}
//$ RL 2011-6-10: ] fabs added (same as abs)
double CMRound::Double()
{
  double x = Arg(0)->Double();
  if (x >= 0) {return floor(x);}
  return ceil(x);
}

double CMFloor::Double()
{
  return floor(Arg(0)->Double());
}

double CMCeil::Double()
{
  return ceil(Arg(0)->Double());
}

double CMHeaviside::Double()
{
  double x = Arg(0)->Double();
  if (x >= 0) {return 1;}
  return 0;
}

double CMSgn::Double()
{
  double a = Arg(0)->Double();
  if (a > 0) {return 1;}
  if (a < 0) {return -1;}
  return 0;
}

double CMPPD::Double()
{
  double a = Arg(0)->Double();
  if (a > 0) {return a-floor(a);}
  if (a < 0) {return a-ceil(a);}
  return 0;
}

double CMVAbs::Double()
{
	int2 rc;
	Arg()->Dim(rc);

	double len = 0;
	for(int i=1; i<=rc(2); i++)
	{
		double comp = Arg()->Get(1,i)->Double();
		len += pow(comp,2.);
	}
	return pow(len,.5);
}

double CMMatRows::Double()
{
	int2 rc;
	Arg()->Dim(rc);
	return rc(1);
}

double CMMatCols::Double()
{
	int2 rc;
	Arg()->Dim(rc);
	return rc(2);
}


double CMRandom::Double()
{
	double a = Arg(0)->Double();
	double r = rand()/(double)RAND_MAX;
	return r*a;
}

double CMRandomSeed::Double()
{
	double s = Arg(0)->Double();
	srand(s);
	return s;
}

//++++++++++++++++++  TRANSPOSE  +++++++++++++++++

CMathObj* CMTranspose::Get(int row, int col)
{
  return Arg()->Get(col,row);
}

void CMTranspose::Dim(int2& rc)
{
  Arg()->Dim(rc);
  rc.Swap();
}

double CMTranspose::MatDouble(int row, int col)
{
  return Arg()->MatDouble(col, row);
}


//++++++++++++++++++  CBaseObj +++++++++++++++++++

CBaseObj::CBaseObj() : CFunction()
{
};

mystr CBaseObj::Name() 
{
  return mystr("@BaseObj");
}

TMathObj CBaseObj::MOType() 
{
  return TMOBase;
}

//++++++++++++++++++  CMEnd +++++++++++++++++++

CMEnd::CMEnd() : CMathObj() 
{
};

mystr CMEnd::Name() 
{
  return mystr("@CMEnd");
}

TMathObj CMEnd::MOType() 
{
  return TMOEnd;
}

//++++++++++++++++++  CMBracket +++++++++++++++++++

mystr CMBracket::String()
{
  return mystr("_(")+Arg(0)->String()+mystr(")_");
}

mystr CMBracket::StringAll()
{
  return mystr("_(")+Arg(0)->StringAll()+mystr(")_");
}


//++++++++++++++++++  CMatrix +++++++++++++++++++

CMatrix::CMatrix(): CMathObj()
{
  generated = 0;
  rows = 1;
  cols = 0;
  mathobjs = new ACMathObj();
}

CMatrix::CMatrix(CMathObj* cmorowsi, CMathObj* cmocolsi)
{
  CMatrix();
  cmorows = cmorowsi; cmocols = cmocolsi;
}

CParseObj* CMatrix::GetCopy() 
{
  CMatrix* o = new CMatrix();

  for (int i = 1; i <= mathobjs->Length(); i++)
  {
    o->mathobjs->Elem(i) = mathobjs->Get(i);
  }
  o->rows = rows;
  o->cols = cols;
  o->cmorows = cmorows;
  o->cmocols = cmocols;
  o->generated = generated;
  return o;
}

CMathObj* CMatrix::Arg(int arg)
{
  if (arg < mathobjs->Length())
  {
    return mathobjs->Get(arg+1);
  } else 
  {
    MyError("Execution Error: range fault in matrix (Arg)");
    return new CMEnd();
  }
}

void CMatrix::SetArg(CMathObj* arg, int argnum)
{
  if (argnum < mathobjs->Length())
  {
    mathobjs->Elem(argnum + 1) = arg;
  } else {MyError("Execution Error: range fault in matrix (SetArg)");}
}

int CMatrix::GetNOArg()
{
  return mathobjs->Length();
}

double CMatrix::Double()
{
  Generate();
  int2 rc;
  Dim(rc);
  if (rc.Get(1) == 1 && rc.Get(2) == 1) {return mathobjs->Get(1)->Double();}
  else {MyError("Execution Error: expected value instead of matrix");}
  return 0;
}
/*
double CMatrix::Double()
{
  MyError(mystr("Execution Error: cannot get R1-value of '")+
          String()+mystr("'"));
}
*/
mystr CMatrix::String()
{  
  Generate();
  int i,j;
  mystr str = mystr("matrix(r=")+mystr(rows)+mystr(",c=")+
               mystr(cols)+mystr(",\n  (");
  for (i = 1; i <= rows; i++)
  {
    for (j = 1; j <= cols; j++)
    {
      str += mathobjs->Get(j+(i-1)*cols)->String();
      if (!(j == cols)) {str += mystr(",");}
    } 
    if (i == rows) {str += mystr(")");}
    else {str += ";\n   ";}
  } 
  return str;
}

mystr CMatrix::StringAll()
{  
  Generate();
  int i,j;
  mystr str = mystr("matrix(r=")+mystr(rows)+mystr(",c=")+
               mystr(cols)+mystr(",\n  (");
  for (i = 1; i <= rows; i++)
  {
    for (j = 1; j <= cols; j++)
    {
      str += mathobjs->Get(j+(i-1)*cols)->StringAll();
      if (!(j == cols)) {str += mystr(",");}
    } 
    if (i == rows) {str += mystr(")");}
    else {str += ";\n   ";}
  } 
  return str;
}

void CMatrix::Dim(int2& rc)
{
  Generate();
  rc.Get(1) = (int)rows;
  rc.Get(2) = (int)cols; 
}

void CMatrix::Elem(int row, int col, CMathObj* replace)
{
  Generate();
  mathobjs->Elem(col+(row-1)*cols) = replace;  
}

CMathObj* CMatrix::Get(int row, int col)
{
  Generate();
  if (col+(row-1)*cols <= mathobjs->Length())
  {
    return mathobjs->Get(col+(row-1)*cols);
  } else
  {
    MyError("Execution Error: tried to access a non-initialized matrix element!");
    return new CMEnd();
  }
}

double CMatrix::MatDouble(int row, int col)
{
  Generate();
  return mathobjs->Get(col+(row-1)*cols)->Double();
}


void CMatrix::Generate()
{
  if (!generated) 
  {
    rows = (int)cmorows->Double();
    cols = (int)cmocols->Double();
    generated = 1;
  }
}

void CMatrix::SetDimensions(CMathObj* mor, CMathObj* moc)
{
  cmorows = mor; cmocols = moc;
}

void CMatrix::AddElem(CMathObj* cm)
{
  mathobjs->Add(cm);
}
  
//++++++++++++++++++  CVecFunction +++++++++++++++++++
CVecFunction::CVecFunction(): CMathObj()
{
  args = new ACMathObj(GetNOArg());
  for(int i = 1; i <= GetNOArg(); i++)
  {
    args->Elem(i) = NULL;
  }
}

CParseObj* CVecFunction::GetCopy() 
{
  CVecFunction* o = new CVecFunction();
  
  o->args = new ACMathObj(GetNOArg());

  for(int i = 1; i <= GetNOArg(); i++)
  {
    args->Elem(i) = o->args->Get(i);
  }
  return o;
}

CMathObj* CVecFunction::Arg(int arg)
{
  if (arg >= 0 && arg < GetNOArg())
  {
    return args->Get(arg+1);
  } else {SysError("range fault in CVecFunction_Arg"); return new CMEnd();}
}

void CVecFunction::SetArg(CMathObj* arg, int argnum)
{
  if (argnum >= 0 && argnum < GetNOArg())
  {
    args->Elem(argnum+1) = arg;
    arg->SetParent(this);
  } else {SysError("range fault in CVecFunction_Arg");}  
}

mystr CVecFunction::String()
{
  mystr rs=Name()+mystr("(");
  for (int i=1; i <= GetNOArg(); i++)
  {
    rs += args->Get(i)->String();
    if (i < GetNOArg()) {rs += mystr(",");} 
  }
  rs += mystr(")");
  return rs;
}

mystr CVecFunction::StringAll()
{
  mystr rs=Name()+mystr("(");
  for (int i=1; i <= GetNOArg(); i++)
  {
    rs += args->Get(i)->StringAll();
    if (i < GetNOArg()) {rs += mystr(",");} 
  }
  rs += mystr(")");
  return rs;
}

//++++++++++++++++++  CMMatElem +++++++++++++++++++

CParseObj* CMMatElem::GetCopy() 
{
  CMMatElem* o = new CMMatElem();
  
  for(int i = 0; i < GetNOArg(); i++)
  {
    SetArg(o->Arg(i), i);
  }
  return o;
}

mystr CMMatElem::String()
{
  mystr rs;

  if (Arg(0)->MOType() == TMOVariable)
  {
    rs = Arg(0)->String();
  } else
  {
    rs = mystr("(")+Arg(0)->String()+mystr(")");
  }

  rs += mystr("[")+Arg(1)->String()+mystr(",")+Arg(2)->String()+mystr("]");

  return rs;
}

mystr CMMatElem::StringAll()
{
  mystr rs;

  if (Arg(0)->MOType() == TMOVariable)
  {
    rs = Arg(0)->StringAll();
  } else
  {
    rs = mystr("(")+Arg(0)->StringAll()+mystr(")");
  }

  rs += mystr("[")+Arg(1)->StringAll()+mystr(",")+Arg(2)->StringAll()+mystr("]");

  return rs;
}

double CMMatElem::Double()
{
  return Arg(0)->MatDouble(Arg(1)->Double(),Arg(2)->Double());
}

void CMMatElem::SetVar(CMathObj* cm)
{
  if (Arg(0)->MOType() == TMOVariable)
  {
    ((CMVariable*)Arg(0))->Elem(Arg(1)->Double(), Arg(2)->Double(), cm);
  } else {MyError("System-Error: CMMatElem::SetVar");}
}

