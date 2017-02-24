//#**************************************************************
//#
//# filename:             elementdata.cpp
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						March 2010
//# description:          Utilities for general data structure ElementData and ElementDataContainer
//#                       
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

#include "tarray.h"
#include "mystring.h"
#include "elementdata.h"

ElementData::ElementData()
{
	dataname = 0;
	vdata = 0;
	vdatalen = 0;
	tdata = 0;
	locked = 0;
	hastooltip = 0;
	tooltiptext = 0;
	idata = 0;
	idata2 = 0;
	idata3 = 0;
	groupnum = 0;
	minval = 1;
	maxval = 0;
	oneLimit = 0;
	edc = 0;
	type = 0;
};

ElementData::ElementData(const ElementData& e)
{
	dataname = 0;
	vdata = 0;
	vdatalen = 0;
	tdata = 0;
	idata = 0;
	idata2 = 0;
	idata3 = 0;
	groupnum = 0;
	minval = 1;
	maxval = 0;
	oneLimit = 0;
	locked = 0;
	hastooltip = 0;
	tooltiptext = 0;
	edc = 0;
	type = 0;

	CopyFrom(e);
};


void ElementData::CopyFrom(const ElementData& e)
{
	Reset(); //free data

	vdatalen = e.vdatalen;
	type = e.type;
	idata = e.idata;
	idata2 = e.idata2;
	idata3 = e.idata3;
	groupnum = e.groupnum;
	locked = e.locked;
	minval = e.minval;
	maxval = e.maxval;
	oneLimit = e.oneLimit;
	hastooltip = e.hastooltip;

	if (e.edc != 0)
	{
		edc = e.edc->GetCopy();
	} 
	else 
	{
		edc = 0;
	}
	if (e.dataname != 0)
	{
		int length = strlen(e.dataname);
		dataname = new char[length + 1];
		strcpy(dataname, e.dataname);
	} else 
	{
		dataname = 0;
	}
	if (e.vdata != 0)
	{
		vdata = new double[vdatalen];
		for (int i=0; i < vdatalen; i++)
		{
			vdata[i] = e.vdata[i];
		}
	} else 
	{
		vdata = 0;
	}
	if (e.tdata != 0)
	{
		int length = strlen(e.tdata);
		tdata = new char[length + 1];
		strcpy(tdata, e.tdata);
	} 
	else 
	{
		tdata = 0;
	}
	if (e.tooltiptext != 0)
	{
		int length = strlen(e.tooltiptext);
		tooltiptext = new char[length + 1];
		strcpy(tooltiptext, e.tooltiptext);
	} 
	else 
	{
		tooltiptext = 0;
	}
}

void ElementData::Reset()
{
	locked = 0;
	minval = 1;
	maxval = 0;
	hastooltip = 0;

	if (vdatalen != 0 && vdata != 0)
	{
		delete[] vdata;
		vdata = 0;
	}
	if (dataname != 0)
	{
		delete[] dataname;
		dataname = 0;
	}
	if (tdata != 0)
	{
		delete[] tdata;
		tdata = 0;
	}
	if (tooltiptext != 0)
	{
		delete[] tooltiptext;
		tooltiptext = 0;
	}
	if (edc != 0)
	{
		delete edc;
		edc = 0;
	}
}


void ElementData::SetVector2D(double vx, double vy, const char* name,double minvalI, double maxvalI, bool onlyLowerBorder, bool onlyUpperBorder)
{
	Reset();
	SetDataName(name);
	type = 4+32;
	vdatalen = 2;
	vdata = new double[vdatalen];
	vdata[0] = vx;
	vdata[1] = vy;
	SetMinMaxVal(minvalI, maxvalI, onlyLowerBorder, onlyUpperBorder);
}
void ElementData::SetVector3D(double vx, double vy, double vz, const char* name, double minvalI, double maxvalI, bool onlyLowerBorder, bool onlyUpperBorder)
{
	Reset();
	SetDataName(name);
	type = 4+32;
	vdatalen = 3;
	vdata = new double[vdatalen];
	vdata[0] = vx;
	vdata[1] = vy;
	vdata[2] = vz;
	SetMinMaxVal(minvalI, maxvalI, onlyLowerBorder, onlyUpperBorder);
}

void ElementData::SetVector(double* v, int len, const char* name, double minvalI, double maxvalI, bool onlyLowerBorder, bool onlyUpperBorder)
{
	Reset();
	SetDataName(name);
	type = 4;
	vdatalen = len;

	if (len > 0)
	{
		vdata = new double[vdatalen];
		for (int i = 0; i < vdatalen; i++)
		{
			vdata[i] = v[i];
		}
	}
	SetMinMaxVal(minvalI, maxvalI, onlyLowerBorder, onlyUpperBorder);
}

void ElementData::SetMatrix(double* v, int rows, int cols, const char* name)
{
	Reset();

	SetDataName(name);
	type = 1024;
	vdatalen = rows*cols;

	idata = rows;
	idata2 = cols;

	if (vdatalen > 0)
	{
		vdata = new double[vdatalen];
		for (int i = 0; i < vdatalen; i++)
		{
			vdata[i] = v[i];
		}
	}
}

void ElementData::SetVector2DList(double* v, int n, const char* name) // $ MSax 2013-07-09 : added
{
	Reset();

	SetDataName(name);
	type = 1024+4096+16384; //matrix & variable length & fixed number of columns
	vdatalen = 2*n;

	idata = n;
	idata2 = 2; //number of columns = 2

	if (vdatalen > 0)
	{
		vdata = new double[vdatalen];
		for (int i = 0; i < vdatalen; i++)
		{
			vdata[i] = v[i];
		}
	}
}

void ElementData::SetMinMaxVal(double minvalI, double maxvalI, bool onlyLowerBorder, bool onlyUpperBorder)
{
	// 4 cases for boundaries are possible:
	// oneLimit = 0: 
	//			minval<=maxval		upper and lower boundary are active:	[...]
	//			minval >maxval		no boundary is active:								 ...
	// onLimit = 1:
	//			minval<=maxval		only upper boundary is active:				 ...]
	//			minval >maxval		only lower boundary is active:				[...

	if(onlyLowerBorder)
	{
		minval = minvalI;
		maxval = minval-1;
	}
	else if(onlyUpperBorder)
	{
		maxval = maxvalI;
		minval = maxval-1;
	}
	else
	{
		maxval = maxvalI;
		minval = minvalI;
	}

	oneLimit = onlyUpperBorder || onlyLowerBorder;
}

void ElementData::SetInt(int num, const char* name, int minvalI, int maxvalI, bool onlyLowerBorder, bool onlyUpperBorder)
{
	Reset();
	SetDataName(name);
	type = 1;
	idata = num;
	double min = (double)minvalI;
	double max = (double)maxvalI;
	SetMinMaxVal(min, max, onlyLowerBorder, onlyUpperBorder);
}


void ElementData::SetEDC(ElementDataContainer* edcI, const char* name)
{
	Reset();
	SetDataName(name);
	type = 8192;

	edc = edcI->GetCopy();
}

void ElementData::SetEDCno_copy(ElementDataContainer* edcI, const char* name)
{
	Reset();
	SetDataName(name);
	type = 8192;

	edc = edcI;
}

ElementDataContainer* ElementData::GetEDC()
{
	return edc;
}

const ElementDataContainer* ElementData::GetEDC() const
{
	return edc;
}

int ElementDataContainer::Find(const char* name) const
{
	for (int i=1; i <= elemdata.Length(); i++)
	{
		if (_stricmp(name, elemdata(i)->GetDataName()) == 0) return i; //stricmp ignores case!!!!
	}
	return 0;
}
//$ RL 2011-6-23:[find all indices in an EDC of entries with same name.
void ElementDataContainer::FindAll(const char* name, TArray<int>& indices) const
{
	indices.SetLen(0);
	for (int i=1; i <= elemdata.Length(); i++)
	{
		if (_stricmp(name, elemdata(i)->GetDataName()) == 0)
		{
			indices.Add(i); //stricmp ignores case!!!!
		}
	}
}
//$ RL 2011-6-23:]find all indices in an EDC of entries with same name.


const ElementData* ElementDataContainer::TreeFind(const char* name) const
{
	mystr str(name);
	if (str.Find('.') != -1)
	{
		int pos = 0;
		mystr dataname;
		str.GetUntil(pos, '.', dataname);
		mystr str_rest = str.SubString(pos+1, str.Length()-1);

	  int i = Find(dataname);
		if (i && Get(i).IsEDC())
		{
			return Get(i).GetEDC()->TreeFind(str_rest.c_str());
		}
		return 0;

	}
	else //no tree structure searched for
	{
	  int i = Find(name);
		if (i) return GetPtr(i);
		return 0;
	}
}

ElementData* ElementDataContainer::TreeFind(const char* name)
{
	mystr str(name);
	if (str.Find('.') != -1)
	{
		int pos = 0;
		mystr dataname;
		str.GetUntil(pos, '.', dataname);
		mystr str_rest = str.SubString(pos+1, str.Length()-1);

	  int i = Find(dataname);
		if (i && Get(i).IsEDC())
		{
			return Get(i).GetEDC()->TreeFind(str_rest.c_str());
		}
		return 0;

	}
	else //no tree structure searched for
	{
	  int i = Find(name);
		if (i) return GetPtr(i);
		return 0;
	}
}

//delete an element in tree, name is e.g. "root.child1.child2.leave"; the name of the deleted object is "leave"
void ElementDataContainer::TreeDelete(const char* name)
{
	mystr str(name);
	if (str.Find('.') != -1)
	{
		int pos = 0;
		mystr dataname;
		str.GetUntil(pos, '.', dataname);
		mystr str_rest = str.SubString(pos+1, str.Length()-1);

	  int i = Find(dataname);
		if (i && Get(i).IsEDC())
		{
			Get(i).GetEDC()->TreeDelete(str_rest.c_str());
			if(Get(i).GetEDC()->Length()==0) {Delete(i);}		// if it was the last leave, delete the child
		}
	}
	else
	{
	  int i = Find(name);
		if(i)	{	Delete(i);}		// delete the leave
	}
}

extern void TIMBSWarningHandle(const char* warn, int use_instant_message_text);

//find element, assume that it is double; if not exists use default value
//accepts integer and double!
double ElementDataContainer::TreeGetDouble(const char* name, double default_val) const
{
	const ElementData* ed = TreeFind(name);
	if (ed)
	{
		if (ed->IsBool()) return (double)ed->GetInt();
		if (ed->IsInt()) return (double)ed->GetInt();
		if (ed->IsDouble()) return ed->GetDouble();
	}

	TIMBSWarningHandle(mystr("'")+mystr(name)+mystr("' not found in TreeGetDouble"),EDCWarningLevel());

	return default_val;
}

//find element, assume that it is int; if not exists use default value
//accepts integer and double!
int ElementDataContainer::TreeGetInt(const char* name, int default_val) const
{
	const ElementData* ed = TreeFind(name);
	if (ed)
	{
		if (ed->IsBool()) return ed->GetInt();
		if (ed->IsInt()) return ed->GetInt();
		if (ed->IsDouble()) return (int)ed->GetDouble();
	}

	TIMBSWarningHandle(mystr("'")+mystr(name)+mystr("' not found in TreeGetInt"),EDCWarningLevel());
	
	return default_val;
}

//find element, assume that it is int; if not exists use default value
//accepts integer and double!
int ElementDataContainer::TreeGetBool(const char* name, int default_val) const
{
	const ElementData* ed = TreeFind(name);
	if (ed)
	{
		if (ed->IsInt()) return ed->GetInt();
		if (ed->IsBool()) return ed->GetInt();
		if (ed->IsDouble()) return (int)ed->GetDouble();
	}

	TIMBSWarningHandle(mystr("'")+mystr(name)+mystr("' not found in TreeGetBool"),EDCWarningLevel());
	
	return default_val;
}

//find element, assume that it is char*; if not exists use default value
const char* ElementDataContainer::TreeGetString(const char* name, char* default_val) const
{
	const ElementData* ed = TreeFind(name);
	if (ed)
	{
		if(ed->GetText())
		{
			return ed->GetText();
		}
	}

	TIMBSWarningHandle(mystr("'")+mystr(name)+mystr("' not found in TreeGetString"),EDCWarningLevel());

	return default_val;
}


int ElementDataContainer::TreeGetVector3D(const char* name, double& vx,  double& vy,  double& vz) const
{
	const ElementData* ed = TreeFind(name);
	if (ed)
	{
		if (ed->IsVector() && ed->GetVectorLen() == 3)
		{
			ed->GetVector(vx, vy, vz);
			return 1;
		}
	}

	TIMBSWarningHandle(mystr("'")+mystr(name)+mystr("' not found in TreeGetVector3D"),EDCWarningLevel());

	return 0;
}

int ElementDataContainer::TreeGetVector2D(const char* name, double& vx,  double& vy) const
{
	const ElementData* ed = TreeFind(name);
	if (ed)
	{
		if (ed->IsVector() && ed->GetVectorLen() == 2)
		{
			ed->GetVector(vx, vy);
			return 1;
		}
	}

	TIMBSWarningHandle(mystr("'")+mystr(name)+mystr("' not found in TreeGetVector2D"),EDCWarningLevel());

	return 0;
}

//find element, assume that it is Vector; if not exists use default value
int ElementDataContainer::TreeGetVector(const char* name, double** v, int& len) const
{
	len = 0;
	*v = 0;
	const ElementData* ed = TreeFind(name);
	if (ed)
	{
		if (ed->IsVector())
		{
			*v = ed->GetVector();
			len = ed->GetVectorLen();
			return 1;
		}
	}

	TIMBSWarningHandle(mystr("'")+mystr(name)+mystr("' not found in TreeGetVector"),EDCWarningLevel());

	return 0;
}

//find element, assume that it is MAtrix; if not exists use default value
int ElementDataContainer::TreeGetMatrix(const char* name, double** v, int& rows, int& cols) const
{
	rows = 0;
	cols = 0;
	*v = 0;
	const ElementData* ed = TreeFind(name);
	if (ed)
	{
		if (ed->IsMatrix())
		{
			*v = ed->GetMatrix();
			rows = ed->GetMatrixRows();
			cols = ed->GetMatrixCols();
			return 1;
		}
	}

	TIMBSWarningHandle(mystr("'")+mystr(name)+mystr("' not found in TreeGetMatrix"),EDCWarningLevel());

	return 0;
}

//find element, assume that it is MAtrix; if not exists use default value
int ElementDataContainer::TreeGetVector2DList(const char* name, double** v, int& n) const // $ MSax 2013-07-09 : added
{
	n = 0;
	*v = 0;
	const ElementData* ed = TreeFind(name);
	if (ed)
	{
		if (ed->IsMatrix())
		{
			//if(ed->GetMatrixCols() != 2)
			//{
			//	ed->SetMatrixSize(n, 2);
			//	for (int i=1; i<=n; i++)
			//	{
			//		ed->SetMatrixVal(i,1,0);
			//		ed->SetMatrixVal(i,2,0);
			//	}
			//}

			if(ed->GetMatrixCols() != 2 && ed->GetMatrixCols() !=0)
			{
				TIMBSWarningHandle(mystr("'")+mystr(name)+mystr("' has wrong dimensions"),EDCWarningLevel());
				return 0;
			}

			*v = ed->GetMatrix();
			n = ed->GetMatrixRows();

			return 1;
		}
	}

	TIMBSWarningHandle(mystr("'")+mystr(name)+mystr("' not found in TreeGetVector2DList"),EDCWarningLevel());

	return 0;
}


int ElementDataContainer::Add(const ElementData& ed) 
{
	return elemdata.Add(ed.GetCopy());
}


void ElementDataContainer::SplitIntoTreeAndElementName(const char* name, mystr& tree, mystr& elementname)
{
	mystr str(name);
	int pos = 0;
	int posold;

	//either finds '.' and gets position of last '.', or returns -1
	while (pos != -1)
	{
		posold = pos;
		pos = str.Find(pos+1, '.');
	}

	//in case of 0, the '.' is either at first position (invalid), or the '.' has not been found
	if (posold > 0) //split tree string into tree and elementname: str="part1.part2.part3" ==> tree = "part1.part2", elementname="part3" 
	{
		tree = str.SubString(0, posold-1);
		elementname = str.SubString(posold+1, str.Length()-1);
	}
	else
	{
		tree = "";
		elementname = str;
	}
}

	//add or replace element in EDC: double; e.g. name="tree1.leave2.data"
void ElementDataContainer::TreeSetDouble(const char* name, double val)
{
	mystr tree, elementname;
	SplitIntoTreeAndElementName(name, tree, elementname);

	ElementData ed;
	ed.SetDouble(val, elementname);
	TreeSet(tree, ed);
}

//add or replace element in EDC: int; e.g. name="tree1.leave2.data"
void ElementDataContainer::TreeSetInt(const char* name, int val)
{
	mystr tree, elementname;
	SplitIntoTreeAndElementName(name, tree, elementname);

	ElementData ed;
	ed.SetInt(val, elementname);
	TreeSet(tree, ed);
}

//add or replace element in EDC: bool; e.g. name="tree1.leave2.data"
void ElementDataContainer::TreeSetBool(const char* name, int val)
{
	mystr tree, elementname;
	SplitIntoTreeAndElementName(name, tree, elementname);

	ElementData ed;
	ed.SetBool(val, elementname);
	TreeSet(tree, ed);
}

//add or replace element in EDC: string; e.g. name="tree1.leave2.data"
void ElementDataContainer::TreeSetString(const char* name, const char* val)
{
	mystr tree, elementname;
	SplitIntoTreeAndElementName(name, tree, elementname);

	ElementData ed;
	ed.SetText(val, elementname);
	TreeSet(tree, ed);
}


	//add or replace element in EDC with COMMENT: double; e.g. name="tree1.leave2.data"
void ElementDataContainer::TreeSetDoubleC(const char* name, double val, const char* comment)
{
	mystr tree, elementname;
	SplitIntoTreeAndElementName(name, tree, elementname);

	ElementData ed;
	ed.SetDouble(val, elementname);
	ed.SetToolTipText(comment);
	TreeSet(tree, ed);
}

//add or replace element in EDC with COMMENT: int; e.g. name="tree1.leave2.data"
void ElementDataContainer::TreeSetIntC(const char* name, int val, const char* comment)
{
	mystr tree, elementname;
	SplitIntoTreeAndElementName(name, tree, elementname);

	ElementData ed;
	ed.SetInt(val, elementname);
	ed.SetToolTipText(comment);
	TreeSet(tree, ed);
}

//add or replace element in EDC with COMMENT: bool; e.g. name="tree1.leave2.data"
void ElementDataContainer::TreeSetBoolC(const char* name, int val, const char* comment)
{
	mystr tree, elementname;
	SplitIntoTreeAndElementName(name, tree, elementname);

	ElementData ed;
	ed.SetBool(val, elementname);
	ed.SetToolTipText(comment);
	TreeSet(tree, ed);
}

//add or replace element in EDC with COMMENT: string; e.g. name="tree1.leave2.data"
void ElementDataContainer::TreeSetStringC(const char* name, const char* val, const char* comment)
{
	mystr tree, elementname;
	SplitIntoTreeAndElementName(name, tree, elementname);

	ElementData ed;
	ed.SetText(val, elementname);
	ed.SetToolTipText(comment);
	TreeSet(tree, ed);
}

//find element, assume that it is Vector; if not exists use default value
void ElementDataContainer::TreeSetVectorC(const char* name, double* v, int len, const char* comment)
{
	mystr tree, elementname;
	SplitIntoTreeAndElementName(name, tree, elementname);

	ElementData ed;
	ed.SetVector(v,len,elementname);
	ed.SetToolTipText(comment);
	TreeSet(tree, ed);
}

void ElementDataContainer::TreeAdd(const char* tree, const ElementData& ed)
//void ElementDataContainer::TreeAdd(mystr tree, const ElementData& ed)
{
	mystr str(tree);
	if (str.Length() != 0)
	{
		int pos = 0;
		mystr str_rest;
		if (str.Find('.') != -1) //split tree string into part1 and rest: str="part1.part2.part3" ==> str_new = "part1", str_rest="part2.part3" 
		{
			mystr str_new;
			str.GetUntil(pos, '.', str_new);
			str_rest = str.SubString(pos+1, str.Length()-1);
			str = str_new;
		}

		//search if node exists
		int i = Find(str);
		if (i && Get(i).IsEDC()) 
		{
			Get(i).GetEDC()->TreeAdd(str_rest.c_str(), ed);
		}
		else
		{
			//create a new element data container at node
			ElementData ed_edc;
			ElementDataContainer edc_new;

			edc_new.TreeAdd(str_rest.c_str(), ed);
			ed_edc.SetEDC(&edc_new, str);

			Add(ed_edc);
		}

	}
	else //no tree structure searched for
	{
		Add(ed);
	}
}
//
//void ElementDataContainer::TreeAdd1(mystr tree, const ElementData& ed)
//{
//	mystr str(tree);
//	if (str.Length() != 0)
//	{
//		int pos = 0;
//		mystr str_rest;
//		if (str.Find('.') != -1) //split tree string into part1 and rest: str="part1.part2.part3" ==> str_new = "part1", str_rest="part2.part3" 
//		{
//			mystr str_new;
//			str.GetUntil(pos, '.', str_new);
//			str_rest = str.SubString(pos+1, str.Length()-1);
//			str = str_new;
//		}
//
//		//search if node exists
//		int i = Find(str);
//		if (i && Get(i).IsEDC()) 
//		{
//			Get(i).GetEDC()->TreeAdd2(str_rest.c_str(), ed);
//		}
//		else
//		{
//			//create a new element data container at node
//			ElementData ed_edc;
//			ElementDataContainer edc_new;
//
//			edc_new.TreeAdd2(str_rest.c_str(), ed);
//			ed_edc.SetEDC(&edc_new, str);
//
//			Add(ed_edc);
//		}
//
//	}
//	else //no tree structure searched for
//	{
//		Add(ed);
//	}
//}
//
//void ElementDataContainer::TreeAdd2(mystr tree, const ElementData& ed)
//{
//	mystr str(tree);
//	if (str.Length() != 0)
//	{
//		int pos = 0;
//		mystr str_rest;
//		if (str.Find('.') != -1) //split tree string into part1 and rest: str="part1.part2.part3" ==> str_new = "part1", str_rest="part2.part3" 
//		{
//			mystr str_new;
//			str.GetUntil(pos, '.', str_new);
//			str_rest = str.SubString(pos+1, str.Length()-1);
//			str = str_new;
//		}
//
//		//search if node exists
//		int i = Find(str);
//		if (i && Get(i).IsEDC()) 
//		{
//			Get(i).GetEDC()->TreeAdd3(str_rest.c_str(), ed);
//		}
//		else
//		{
//			//create a new element data container at node
//			ElementData ed_edc;
//			ElementDataContainer edc_new;
//
//			edc_new.TreeAdd3(str_rest.c_str(), ed);
//			ed_edc.SetEDC(&edc_new, str);
//
//			Add(ed_edc);
//		}
//
//	}
//	else //no tree structure searched for
//	{
//		Add(ed);
//	}
//}
//
//void ElementDataContainer::TreeAdd3(mystr tree, const ElementData& ed)
//{
//	mystr str(tree);
//	if (str.Length() != 0)
//	{
//		int pos = 0;
//		mystr str_rest;
//		if (str.Find('.') != -1) //split tree string into part1 and rest: str="part1.part2.part3" ==> str_new = "part1", str_rest="part2.part3" 
//		{
//			mystr str_new;
//			str.GetUntil(pos, '.', str_new);
//			str_rest = str.SubString(pos+1, str.Length()-1);
//			str = str_new;
//		}
//
//		//search if node exists
//		int i = Find(str);
//		if (i && Get(i).IsEDC()) 
//		{
//			Get(i).GetEDC()->TreeAdd3(str_rest.c_str(), ed);
//		}
//		else
//		{
//			//create a new element data container at node
//			ElementData ed_edc;
//			ElementDataContainer edc_new;
//
//			edc_new.TreeAdd3(str_rest.c_str(), ed);
//			ed_edc.SetEDC(&edc_new, str);
//
//			Add(ed_edc);
//		}
//
//	}
//	else //no tree structure searched for
//	{
//		Add(ed);
//	}
//}



//set an element in the tree of edc; if the tree does not exist, build the tree and insert data; if the data exists: replace the data
void ElementDataContainer::TreeSet(const char* tree, const ElementData& ed, int warn_if_items_dont_exist, int copy_ToolTipText)
{
	mystr str = mystr(tree);
	if (str.Length() != 0) {str += ".";}

	str += ed.GetDataName();
  
	ElementData* edfound = TreeFind(str);
	if (edfound != 0)
	{
		if(edfound->IsLocked())	//$ DR 2013-04-09
		{
			TIMBSWarningHandle(mystr("The entry '") + local_warn_tree + str + mystr("' is marked as readonly in the original ElementDataContainer. The value is therefore not changed! If the entry is the variable of the parameter variation, ignore this warning.\n"),0);	
		}
		else	// replace the data
		{
			if (edfound->GetType() == ed.GetType())
			{
				//copy data, but keep tooltip text
				mystr old_tooltip("");
				if (edfound->HasToolTip()) old_tooltip = edfound->GetToolTipText();
				edfound->CopyFrom(ed);
				if(!copy_ToolTipText)
				{
					edfound->SetToolTipText(old_tooltip);
				}
			}
			//treat int/bool/double the same; this means that the type of the old edfound is kept:
			else if ( (edfound->IsInt() && ed.IsBool()) || (edfound->IsInt() && ed.IsDouble()) || (edfound->IsBool() && ed.IsInt()) || (edfound->IsBool() && ed.IsDouble())
				|| (edfound->IsDouble() && ed.IsInt())  || (edfound->IsDouble() && ed.IsBool()) )
			{
				double val;
				if (ed.IsBool()) val = ed.GetBool();
				if (ed.IsInt()) val = ed.GetInt();
				if (ed.IsDouble()) val = ed.GetDouble();


				if (edfound->IsBool())   {edfound->SetBool((int)val);} //JG, do not replace data name because of case sensitivity
				if (edfound->IsInt())    {edfound->SetInt((int)val);} //JG, do not replace data name because of case sensitivity
				if (edfound->IsDouble()) {edfound->SetDouble(val);} //JG, do not replace data name because of case sensitivity
				//if (edfound->IsBool())   {edfound->SetBool(val, ed.GetDataName());}
				//if (edfound->IsInt())    {edfound->SetInt(val, ed.GetDataName());}
				//if (edfound->IsDouble()) {edfound->SetDouble(val, ed.GetDataName());}

			}
			// write vector2d/3d into vector-ed
			else if ( edfound->IsVector() && ed.IsVectorXYZ() )
			{
				if (ed.GetVectorLen() == 2)
				{
					edfound->SetVector2D(ed.GetVectorVal(1), ed.GetVectorVal(2), ed.GetDataName());
				}
				if (ed.GetVectorLen() == 3)
				{
					edfound->SetVector3D(ed.GetVectorVal(1), ed.GetVectorVal(2), ed.GetVectorVal(3), ed.GetDataName());
				}
			}
			// write first 2 or 3 values of vector into vector2d/3d-ed
			else if ( edfound->IsVectorXYZ() && ed.IsVector() )
			{
				double* dummy = ed.GetVector();
				if (edfound->GetVectorLen() == 2 && ed.GetVectorLen() >= 2 )
				{
					edfound->SetVector2D(dummy[0], dummy[1], ed.GetDataName());
				}
				if (edfound->GetVectorLen() == 3 && ed.GetVectorLen() >= 3 )
				{
					edfound->SetVector3D(dummy[0], dummy[1], dummy[2], ed.GetDataName());
				}
				if (edfound->GetVectorLen() == 3 && ed.GetVectorLen() == 2 )
				{
					edfound->SetVector3D(dummy[0], dummy[1], 0., ed.GetDataName());
				}
			}
			//different types of vectors; length must be same, if vector does not have variable length
			else if ( edfound->IsVector() && ed.IsVector() )
			{
				//independently of exact types (intvalues, etc.), replace vector
				double* dummy = ed.GetVector();
				if ((edfound->GetVectorLen() == ed.GetVectorLen())||(edfound->IsVariableLength()))	//$ DR 2012-1: variable length added
				{
					int oldtype = edfound->GetType();
					edfound->SetVector(dummy, ed.GetVectorLen(), ed.GetDataName());
					edfound->SetType(oldtype);
				}
				else
				{
					TIMBSWarningHandle(mystr("'") + local_warn_tree+mystr(ed.GetDataName())+mystr("': incompatible vector '")+mystr(edfound->GetDataName())+mystr("', check length and type of default entry!!"),EDCWarningLevel());
				}
			}
			else if(edfound->IsMatrix() && ed.IsVector() && (edfound->IsVariableLength() || edfound->GetMatrixCols() == ed.GetVectorLen() && edfound->GetMatrixRows() == 1))
			{	
				//$ SW 2013-11-5: added the case that the matrix has fixed column size and the vector is empty.
				// matrix = vector & matrix has fixed column size => vector = [] or length(vector) == columns(matrix)
				if (edfound->IsVariableLength() && edfound->IsColumnNumberFixed() && ed.GetVectorLen() != 0 && ed.GetVectorLen() != edfound->GetMatrixCols())
				{
					TIMBSWarningHandle(mystr("'") + local_warn_tree + mystr(ed.GetDataName()) + mystr("': incompatible vector '") + mystr(edfound->GetDataName()) + mystr("', check dimensions!"),EDCWarningLevel());
				}
				//matrix = vector & matrix has fixed column size && vector = []
				else if (edfound->IsVariableLength() && edfound->IsColumnNumberFixed() && ed.GetVectorLen() == 0)
				{
					double* dummy = ed.GetVector();
					int oldtype = edfound->GetType();
					edfound->SetMatrix(dummy, 0, edfound->GetMatrixCols(), ed.GetDataName());
					edfound->SetType(oldtype);
				}
				else 
				{
					double* dummy = ed.GetVector();
					int oldtype = edfound->GetType();
					edfound->SetMatrix(dummy, 1, ed.GetVectorLen(), ed.GetDataName());
					edfound->SetType(oldtype);
				}
			}
			else if(edfound->IsMatrix() && ed.IsMatrix() && (edfound->IsVariableLength() || (edfound->GetMatrixCols() == ed.GetMatrixCols() && edfound->GetMatrixRows() == ed.GetMatrixRows())))
			{
				double* dummy = ed.GetMatrix();
				int oldtype = edfound->GetType();
				edfound->SetMatrix(dummy, ed.GetMatrixRows(), ed.GetMatrixCols(), ed.GetDataName());
				edfound->SetType(oldtype);
			}
			else if(edfound->IsMatrix() && ed.IsVector() && ed.GetVectorLen() == 0)
			{
				// $ MSax 2013-07-11: do nothing, because ed is an empty matrix
			}
			else if(edfound->IsVector() && edfound->GetVectorLen() == 1 && (ed.IsDouble() || ed.IsInt()) )
			{
				//$ AD 2013-10-10: write a double in a length 1 vector 
				double dummy;
				if(ed.IsDouble()) dummy = ed.GetDouble();
				else if(ed.IsInt()) dummy = (double) ed.GetInt();
				int oldtype = edfound->GetType();
				edfound->SetVector(&dummy, 1, ed.GetDataName());
				edfound->SetType(oldtype);
			}
			else if(edfound->IsMatrix() && edfound->GetMatrixRows() == 1 && edfound->GetMatrixCols() == 1 && (ed.IsDouble() || ed.IsInt()) )
			{
				//$ AD 2013-10-10: write a double in a length 1 vector 
				double dummy;
				if(ed.IsDouble()) dummy = ed.GetDouble();
				else if(ed.IsInt()) dummy = (double) ed.GetInt();
				int oldtype = edfound->GetType();
				edfound->SetMatrix(&dummy, 1, 1, ed.GetDataName());
				edfound->SetType(oldtype);
			}

			else
			{
				//here we should put a warning, because user could set a string options to a number value
				//TIMBSWarningHandle(mystr("'")+mystr(ed.GetDataName())+mystr("' is incompatible with data type of '")+mystr(edfound->GetDataName())+mystr("'!!"),EDCWarningLevel());
				TIMBSWarningHandle(mystr("'")+local_warn_tree+mystr("' is incompatible with data type of '")+mystr(edfound->GetDataName())+mystr("'!!"),EDCWarningLevel());
				//assert(0);
			}
		}
	}
	else
	{
		if (warn_if_items_dont_exist) 
		{
			//TIMBSWarningHandle(mystr("'") + local_warn_tree + str + mystr("' added (check spelling)."), EDCWarningLevel());
			TIMBSWarningHandle(mystr(" The entry '") + local_warn_tree + str + mystr("' does not exist in the original ElementDataContainer (check spelling)."), EDCWarningLevel());	//$ DR 2013-02-20 new warning, also used for skript language
		}
		TreeAdd(tree, ed);
	}
}

void ElementDataContainer::TreeReplaceEDCDataWith(const ElementDataContainer* edc_new, int warn_if_items_dont_exist, mystr warn_str)
{
	//SetLocalWarnTree(warn_str);
	TreeReplaceEDCDataWithRecursive(edc_new, warn_if_items_dont_exist);
}

void ElementDataContainer::TreeReplaceEDCDataWithRecursive(const ElementDataContainer* edc_new, int warn_if_items_dont_exist)
{
	for(int i = 1; i<=edc_new->Length(); i++)
	{
		ElementData ed_new = edc_new->Get(i);
		mystr name = ed_new.GetDataName();
		int j = Find(name);
		if(j && Get(j).IsEDC()) // same variable names exist
		{			
			Get(j).GetEDC()->TreeReplaceEDCDataWithRecursive(ed_new.GetEDC(), warn_if_items_dont_exist);// replace sub-tree
		}
		else
		{
			TreeSet("",ed_new,warn_if_items_dont_exist, 0);
		}
	}
}  

//$ DR 2012-06-27
//void ElementDataContainer::PrintEDCRecursive(mystr filename)
void ElementDataContainer::PrintEDCRecursive(ofstream &edcout)
{
	//ofstream edcout(filename);
	for(int i = 1; i<=Length(); i++)
	{
		ElementData ed_tmp = Get(i);
		if(Get(i).IsEDC())
		{
			ElementDataContainer *edcp = Get(i).GetEDC();
			edcout <<  "================== \n"; 
			edcout << ed_tmp.GetDataName() << ":\n"; 
			edcp->PrintEDCRecursive(edcout);
		}
		else
		{
			edcout <<  "------------------- \n"; 
			edcout << ed_tmp.GetDataName() << ":\n"; 
			if(ed_tmp.GetType() ==1) 
			{
				edcout << "int = " << ed_tmp.GetInt() << "\n";
			}
			if(ed_tmp.GetType() ==2) 
			{
				edcout << "double = " << ed_tmp.GetDouble() << "\n";
			}
			if(ed_tmp.GetType() ==4) 
			{
				edcout << "vec = " << ed_tmp.GetVector() << "\n";
			}
			if(ed_tmp.GetType() == 8) 
			{
				edcout << "txt = " << ed_tmp.GetText() << "\n";
			}
			if(ed_tmp.GetType() == 16) 
			{
				edcout << "bool = " << ed_tmp.GetBool() << "\n";
			}
		}
	}
}




