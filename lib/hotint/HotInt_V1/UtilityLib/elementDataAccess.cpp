//#***************************************************************************************
//# filename:     elementDataAccess.cpp
//#
//# author:				Johannes Gerstmayr, Yuri Vetyukov
//# 
//# generated:      
//# description:  
//#                       
//# comments:      
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
//#***************************************************************************************

//#include "stdafx.h"
#include "mbs_interface.h"
#include "elementData.h"
#include "mathfunc.h"
#include "elementDataAccess.h"

void WarnElementData(MBS* mbs, const ElementDataContainer& edc, const char* name)
{
	if (mbs == 0) return;

	int en = 0;
	int pos = edc.Find("Element_number");
	if (pos)
	{
		const ElementData& ed = edc.Get(pos);
		if (ed.IsInt())
		{
			en = ed.GetInt();
		}
	}
	mbs->UO() << "Error in element";
	if (en) mbs->UO() << " " << en;
	mbs->UO() << ": Data entry '" << name << "' not found!\n";
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//the following functions are used for simplified EDC access, especially for access by Set/GetElementData functions of Elements
//flag&1 --> needed value
//$JG2012-01-19: because element data EDCs are now hierarchical, these functions have been changed from Find(...) to TreeFind(...) needed to change, 
//                because from now on this function searches in the whole tree, not only in the actual list
//                WARNING: note that the integer return value now is changed from the local position in the EDC to 0/1 only!!!
int GetElemDataBool(MBS* mbs, const ElementDataContainer& edc, const char* name, int& v, int flag)
{
	//old: int pos = edc.Find(name); if (pos) ...

	const ElementData* ed = edc.TreeFind(name); //$JG2012-01-19; because element data EDCs are now hierarchical, these functions have been changed from Find(...) to TreeFind(...) needed to change, because from now on this function searches in the whole tree, not only in the actual list
	int found = (ed != 0);
	if (ed)
	{
		if (ed->IsBool())
		{
			v = ed->GetBool();
		}
	}
	else if (flag&1)
		WarnElementData(mbs, edc, name);

	return found; //old: return pos
}

int GetElemDataInt(MBS* mbs, const ElementDataContainer& edc, const char* name, int& v, int flag)
{
	const ElementData* ed = edc.TreeFind(name); //$JG2012-01-19; needed to change, because from now on this function searches in the whole tree, not only in the actual list
	int found = (ed != 0);
	if (ed)
	{
		if (ed->IsInt())
		{
			v = ed->GetInt();
		}
	}
	else if (flag&1)
		WarnElementData(mbs, edc, name);

	return found;
}

int GetElemDataDouble(MBS* mbs, const ElementDataContainer& edc, const char* name, double& v, int flag)
{
	const ElementData* ed = edc.TreeFind(name); //$JG2012-01-19; needed to change, because from now on this function searches in the whole tree, not only in the actual list
	int found = (ed != 0);
	if (ed)
	{
		if (ed->IsDouble())
		{
			v = ed->GetDouble();
		}
	}
	else if (flag&1)
		WarnElementData(mbs, edc, name);

	return found;
}

int GetElemDataText(MBS* mbs, const ElementDataContainer& edc, const char* name, mystr& str, int flag)
{
	const ElementData* ed = edc.TreeFind(name); //$JG2012-01-19; needed to change, because from now on this function searches in the whole tree, not only in the actual list
	int found = (ed != 0);
	if (ed)
	{
		if (ed->IsText())
		{
			str = ed->GetText();
		}
	}
	else if (flag&1)
		WarnElementData(mbs, edc, name);

	return found;
}

int GetElemDataVector2D(MBS* mbs, const ElementDataContainer& edc, const char* name, Vector2D& v, int flag)
{
	const ElementData* ed = edc.TreeFind(name); //$JG2012-01-19; needed to change, because from now on this function searches in the whole tree, not only in the actual list
	int found = (ed != 0);
	if (ed)
	{
		if (ed->IsVector() && ed->GetVectorLen() == 2)
		{
			ed->GetVector(v.X(), v.Y());
		}
	}
	else if (flag&1)
		WarnElementData(mbs, edc, name);

	return found;
}

int GetElemDataVector3D(MBS* mbs, const ElementDataContainer& edc, const char* name, Vector3D& v, int flag)
{
	const ElementData* ed = edc.TreeFind(name); //$JG2012-01-19; needed to change, because from now on this function searches in the whole tree, not only in the actual list
	int found = (ed != 0);
	if (ed)
	{
		if (ed->IsVector() && ed->GetVectorLen() == 3)
		{
			ed->GetVector(v.X(), v.Y(), v.Z());
		}
	}
	else if (flag&1)
		WarnElementData(mbs, edc, name);

	return found;
}

int GetElemDataVector2D(MBS* mbs, const ElementDataContainer& edc, const char* name, double& v1, double& v2, int flag)
{
	Vector2D v(v1, v2);
	int found = GetElemDataVector2D(mbs, edc, name, v, flag);

	v1 = v.X();
	v2 = v.Y();

	return found;
}

int GetElemDataVector3D(MBS* mbs, const ElementDataContainer& edc, const char* name, double& v1, double& v2, double& v3, int flag)
{
	Vector3D v(v1, v2, v3);
	int found = GetElemDataVector3D(mbs, edc, name, v, flag);

	v1 = v.X();
	v2 = v.Y();
	v3 = v.Z();

	return found;
}

int GetElemDataVector(MBS* mbs, const ElementDataContainer& edc, const char* name, Vector& v, int flag)
{
	const ElementData* ed = edc.TreeFind(name); //$JG2012-01-19; needed to change, because from now on this function searches in the whole tree, not only in the actual list
	int found = (ed != 0);
	if (ed)
	{
		if (ed->IsVector())
		{
			v.SetLen(ed->GetVectorLen());
			for (int i=1; i <= ed->GetVectorLen(); i++)
			{
				v(i) = ed->GetVectorVal(i);
			}
		}
	}
	else if (flag&1)
		WarnElementData(mbs, edc, name);

	return found;
}

int GetElemDataIVector(MBS* mbs, const ElementDataContainer& edc, const char* name, IVector& v, int flag)
{
	const ElementData* ed = edc.TreeFind(name); //$JG2012-01-19; needed to change, because from now on this function searches in the whole tree, not only in the actual list
	int found = (ed != 0);
	if (ed)
	{
		if (ed->IsVector()) 
		{
			v.SetLen(ed->GetVectorLen());
			for (int i=1; i <= ed->GetVectorLen(); i++)
			{
				v(i) = (int)ed->GetVectorVal(i);
			}
		}
	}
	else if (flag&1)
		WarnElementData(mbs, edc, name);

	return found;
}

int GetElemDataMatrix(MBS* mbs, const ElementDataContainer& edc, const char* name, Matrix& v, int flag)
{
	const ElementData* ed = edc.TreeFind(name); //$JG2012-01-19; needed to change, because from now on this function searches in the whole tree, not only in the actual list
	int found = (ed != 0);
	if (ed)
	{
		if (ed->IsMatrix()) 
		{
			v.SetSize(ed->GetMatrixRows(), ed->GetMatrixCols());

			for (int i=1; i <= ed->GetMatrixRows(); i++)
			{
				for (int j=1; j <= ed->GetMatrixCols(); j++)
				{
					v(i,j) = ed->GetMatrixVal(i,j);
				}
			}
		}
	}
	else if (flag&1)
		WarnElementData(mbs, edc, name);

	return found;
}

int GetElemDataMathFunc(MBS* mbs, const ElementDataContainer& edc, const mystr& funcname, MathFunction& mathfunc, int flag) //depreciated!!!!, do not use anymore

{
	int rv;

	//mathfunc
	int mode;
	Matrix data;

	const ElementData* ed = edc.TreeFind(mystr(funcname+mystr("_type")).c_str()); //$JG2012-01-19; needed to change, because from now on this function searches in the whole tree, not only in the actual list
	int found = (ed != 0);
	//int pos = edc.Find(mystr(funcname+mystr("_type")).c_str());
	if (found)
	{
		mystr fntype = funcname+mystr("_type");
		mystr fndata = funcname+mystr("_data");

		GetElemDataInt(mbs, edc, fntype.c_str(), mode, 1);
		GetElemDataMatrix(mbs, edc, fndata.c_str(), data, 1);

		rv = mathfunc.SetData(mode, data); //set data only accepts valid data!
	}
	else
	{
		mystr funcstr;
		Vector coeffs;

		GetElemDataText(mbs, edc, "Math_function_name", funcstr, 1);
		GetElemDataVector(mbs, edc, "Mathfunc_coefficients", coeffs, 1);

		rv = mathfunc.SetData(funcstr, coeffs);
	}





	return rv;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//ElementData access

void SetElemDataVector2D(ElementData& ed, const Vector2D& v, const char* name)
{
	ed.SetVector2D(v.X(), v.Y(), name); 
}

void SetElemDataVector3D(ElementData& ed, const Vector3D& v, const char* name)
{
	ed.SetVector3D(v.X(), v.Y(), v.Z(), name); 
}

void SetElemDataIVector(ElementData& ed, const IVector& v, const char* name)
{
	Vector vv(v.Length());
	for (int i=1; i <= v.Length(); i++) vv(i) = v(i);

	ed.SetVector(vv.GetVecPtr(), vv.Length(), name); 
	ed.SetValuesInt();
}

void SetElemDataVector(ElementData& ed, const Vector& v, const char* name)
{
	ed.SetVector(v.GetVecPtr(), v.Length(), name); 
}

void SetElemDataMatrix(ElementData& ed, const Matrix& v, const char* name)
{
	ed.SetMatrix(v.GetMatPtr(), v.Getrows(), v.Getcols(), name); 
}


//ElementDataContainer access

void SetElemDataVector2D(ElementDataContainer& edc, const Vector2D& v, const char* name, const char* tooltiptext)
{
	ElementData ed;
	ed.SetVector2D(v.X(), v.Y(), name); 
	ed.SetToolTipText(tooltiptext);
	edc.Add(ed);
}

void SetElemDataVector3D(ElementDataContainer& edc, const Vector3D& v, const char* name, const char* tooltiptext)
{
	ElementData ed;
	ed.SetVector3D(v.X(), v.Y(), v.Z(), name); 
	ed.SetToolTipText(tooltiptext);
	edc.Add(ed);
}

void SetElemDataIVector(ElementDataContainer& edc, const IVector& v, const char* name, const char* tooltiptext)
{
	ElementData ed;
	Vector vv(v.Length());
	for (int i=1; i <= v.Length(); i++) vv(i) = v(i);

	ed.SetVector(vv.GetVecPtr(), vv.Length(), name); 
	ed.SetValuesInt();
	ed.SetToolTipText(tooltiptext);
	edc.Add(ed);
}

void SetElemDataVector(ElementDataContainer& edc, const Vector& v, const char* name, const char* tooltiptext)
{
	ElementData ed;

	ed.SetVector(v.GetVecPtr(), v.Length(), name); 
	ed.SetToolTipText(tooltiptext);
	edc.Add(ed);
}

void SetElemDataMatrix(ElementDataContainer& edc, const Matrix& v, const char* name, const char* tooltiptext)
{
	ElementData ed;

	ed.SetMatrix(v.GetMatPtr(), v.Getrows(), v.Getcols(), name); 
	ed.SetToolTipText(tooltiptext);
	edc.Add(ed);
}

void SetElemDataMathFunc(ElementDataContainer& edc, MathFunction& mathfunc, const mystr& funcname) //depreciated!!!!, do not use anymore
{
	ElementData ed;

	//mathfunc
	int mode;
	Matrix data; 
	mathfunc.GetData(mode, data);

	if (mathfunc.GetFuncMode() <= mathfunc.GetMaxFuncMode())
	{
		mystr fntype = funcname+mystr("_type");
		mystr fndata = funcname+mystr("_data");
		ed.SetInt(mode, fntype.c_str(), 0, mathfunc.GetMaxFuncMode()); ed.SetToolTipText("0=no func., 1=polynomial, 2/3/4=piecewise const./linear/quad., 5=harmonic"); edc.Add(ed);
		SetElemDataMatrix(edc, data, fndata.c_str()); 
		edc.Last().SetVariableLength();
		edc.Last().SetToolTipText("columns: 1=coeff | 2/3=time,value | 4=time,pos,vel | 5=freq.,phase,amplitude");
	}
	else
	{
		mystr funcstr;
		Vector coeffs;

		switch (mathfunc.GetFuncMode())
		{
		case TMFsin: //***
			{
				mathfunc.GetData(funcstr, coeffs);
				ed.SetText("Math_function_name", funcstr); ed.SetLocked(1); edc.Add(ed);
				SetElemDataVector(edc, coeffs, "Mathfunc_coefficients"); edc.Last().SetVariableLength(); edc.Last().SetToolTipText("coefficients = [A omega phi], y = A*Sin(omega*x + phi)");

				break;
			}
		case TMFcos: //***
			{
				mathfunc.GetData(funcstr, coeffs);
				ed.SetText("Math_function_name", funcstr); ed.SetLocked(1); edc.Add(ed);
				SetElemDataVector(edc, coeffs, "Mathfunc_coefficients"); edc.Last().SetVariableLength(); edc.Last().SetToolTipText("coefficients = [A omega phi], y = A*Cos(omega*x + phi)");

				break;
			}
		case TMFstaticDynamicFricion: //***
			{
				assert(0 && "ERROR in mbs_communication.cpp: SetStaticDynamicFricionFunction not anymore supported by MathFunction, use TMFExpression instead!"); 

				mathfunc.GetData(funcstr, coeffs);
				ed.SetText("Math_function_name", funcstr); ed.SetLocked(1); edc.Add(ed);
				SetElemDataVector(edc, coeffs, "Mathfunc_coefficients"); edc.Last().SetVariableLength(); edc.Last().SetToolTipText("coefficients = [staticFrictionCoeff dynamicFrictionCoeff zeroZone]");

				break;
			}
		default: ;
		}
	}
	//+++++++++++++++++++++++++++
}

Vector3D EDCTreeGetVector3D(ElementDataContainer& edc, const char* name, Vector3D default_val)
{
	double vx = 0;double vy = 0;double vz = 0;
	if (edc.TreeGetVector3D(name, vx, vy, vz))
	{
		return Vector3D(vx,vy,vz);		
	}
	return default_val;
}

Vector EDCTreeGetVector(ElementDataContainer& edc, const char* name, Vector default_val)
{
	double* v = 0;
	int len;
	if (edc.TreeGetVector(name, &v, len))
	{
		Vector ta(len);
		for(int i=1;i<=len; i++)
		{
			ta(i) = v[i-1];
		}
		return ta;
	}
	return default_val;
}

mystr GetRotUnitStr(int type) //0=rad, 1=degree
{
	if (type == 0)
	{
		return mystr("(rad)");
	}
	return mystr("(°)");
}