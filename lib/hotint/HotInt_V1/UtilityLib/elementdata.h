//#**************************************************************
//#
//# filename:             elementdata.h
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

#ifndef ELEMENTDATA__H 
#define ELEMENTDATA__H

class mystr;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++       data exchance interface for element data      ++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class ElementDataContainer;

class ElementData
{
public:
	ElementData();

	ElementData(const ElementData& e);

	ElementData& operator=(const ElementData& ed) 
	{
		if (this == &ed) {return *this;}
		CopyFrom(ed);
		return *this;
	}

	virtual ElementData* GetCopy() const
	{
		ElementData* ed = new ElementData();
		ed->CopyFrom(*this);
		return ed;
	}

	virtual void CopyFrom(const ElementData& e);

	virtual void Reset();

	virtual ~ElementData()
	{
		Reset();
	}

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//functions to define the ElementData for specific data type:
	virtual void SetBool(int num, const char* name)
	{
		Reset();
		SetDataName(name);
		type = 16;
		idata = num;
	}
	virtual void SetBool(int num) //only change value!
	{
		idata = num;
	}
	virtual void SetBoolGroup(int num, int group, const char* name)
	{
		Reset();
		SetDataName(name);
		type = 16+512;
		idata = num;
		groupnum = group;
	}
	virtual void SetGroup(int groupnumI) {groupnum = groupnumI; type = type|512;} //only add group number

	virtual void SetVariableLength() {type = type|4096;} //vectors of variable length (editable)

	virtual void SetInt(int num) //only change value!
	{
		idata = num;
	}
	virtual void SetInt(int num, const char* name, int minvalI=1, int maxvalI=0, bool onlyLowerBorder=0, bool onlyUpperBorder=0);

	virtual void SetDouble(double num) //only change value!
	{
		vdata[0] = num;
	}
	virtual void SetDouble(double num, const char* name, double minvalI=1, double maxvalI=0, bool onlyLowerBorder=0, bool onlyUpperBorder=0)
	{
		Reset();
		SetDataName(name);
		type = 2;
		vdatalen = 1;
		vdata = new double[vdatalen];
		vdata[0] = num;
		SetMinMaxVal(minvalI, maxvalI, onlyLowerBorder, onlyUpperBorder);
	}

	virtual void SetVectorLen(int i) //only change length!
	{
		if (vdatalen != 0 && vdata != 0)
		{
			delete[] vdata;
			vdata = 0;
		}

		vdatalen = i;
		vdata = new double[vdatalen];
	}
	virtual void SetVectorVal(int i, double v) //only change value!
	{
		vdata[i-1] = v;
	}

	virtual void SetVector2D(double vx, double vy, const char* name, double minvalI=1, double maxvalI=0, bool onlyLowerBorder=0, bool onlyUpperBorder=0);
	virtual void SetVector3D(double vx, double vy, double vz, const char* name, double minvalI=1, double maxvalI=0, bool onlyLowerBorder=0, bool onlyUpperBorder=0);
	virtual void SetVector(double* v, int len, const char* name, double minvalI=1, double maxvalI=0, bool onlyLowerBorder=0, bool onlyUpperBorder=0);
	virtual void SetMatrix(double* v, int rows, int cols, const char* name);
	virtual void SetVector2DList(double* v, int n, const char* name);

	virtual void SetText(const char* textdata, const char* name)
	{
		Reset();
		SetDataName(name);
		type = 8;

		size_t length = strlen(textdata);
		tdata = new char[length + 1];
		strcpy(tdata, textdata);
	}

	//virtual void SetEDC(ElementDataContainer* edcI, mystr name);
	virtual void SetEDC(ElementDataContainer* edcI, const char* name);

	//same as SetEDC, but do not copy edcI (just set reference); be careful with usage!!!
	virtual void SetEDCno_copy(ElementDataContainer* edcI, const char* name);

	virtual void SetWCDriverAction(int action, int option, int value, const char* name) 
	{
		Reset();
		SetDataName(name);
		type = 128;
		idata = action;
		idata2 = option;
		idata3 = value;
	}

	virtual void SetCompAction(int action, int option, int value, const char* name)
	{
		Reset();
		SetDataName(name);
		type = 256;
		idata = action;
		idata2 = option;
		idata3 = value;
	}

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	virtual int GetType() const {return type;}
	virtual void SetType(int typeI) {type = typeI;}

	virtual int IsLocked() const {return locked;}
	virtual void SetLocked(int flag) {locked = flag;}
	virtual int HasToolTip() const {return hastooltip;}
	virtual void SetToolTipText(const char* text) 
	{
		if (tooltiptext != 0) delete[] tooltiptext;

		hastooltip = 1;
		size_t length = strlen(text);
		tooltiptext = new char[length + 1];
		strcpy(tooltiptext, text);
	}
	virtual const char* GetToolTipText() const {return tooltiptext;}

	virtual int IsElementInfo() const {return (type&64) != 0;}
	virtual int IsBool() const {return (type&16) != 0;}
	virtual int IsInt() const {return (type&1) != 0;}
	virtual int IsDouble() const {return (type&2) != 0;}
	virtual int IsVector() const {return (type&4) != 0;}
	virtual int IsVectorXYZ() const {return (type&32) != 0;}
	virtual int IsText() const {return (type&8) != 0;}
	virtual int IsWCDriverAction() const {return (type&128) != 0;}
	virtual int IsCompAction() const {return (type&256) != 0;}
	virtual int IsGroup() const {return (type&512) != 0;}
	virtual int IsMatrix() const {return (type&1024) != 0;}
	virtual int IsValuesInt() const {return (type&2048) != 0;} //vector/matrix values are integer!
	virtual int IsVariableLength() const {return (type&4096) != 0;} //vectors of variable length (editable)
	virtual int IsEDC() const {return (type&8192) != 0;} //pointer to other ElementData
	virtual int IsColumnNumberFixed() const {return (type&16384) != 0;}

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//functions that get element data:
	virtual char* GetDataName() const {return dataname;}
	virtual void SetDataName(const char* name) 
	{
		size_t length = strlen(name);
		dataname = new char[length + 1];
		strcpy(dataname, name);
	}

	virtual int GetInt() const {return idata;}
	virtual int GetMinMaxVal(double& min, double& max, bool& justoneLimit) {min = minval; max = maxval; justoneLimit = oneLimit; return (min <= max);}
	virtual void SetMinMaxVal(double min, double max, bool onlyLowerBorder, bool onlyUpperBorder); //$ DR 2012-11 added

	virtual int GetBool() const {return idata;}
	virtual int GetGroup() const {return groupnum;}
	virtual double GetDouble() const {return vdata[0];}
	virtual double* GetVector() const {return vdata;}	
	virtual void GetVector(double& vx, double& vy) const 
	{
		vx = vdata[0];
		vy = vdata[1];
	}	
	virtual void GetVector(double& vx, double& vy, double& vz) const 
	{
		vx = vdata[0];
		vy = vdata[1];
		vz = vdata[2];
	}	
	virtual int GetVectorLen() const {return vdatalen;}
	virtual double GetVectorVal(int i) const {return vdata[i-1];}

	virtual double* GetMatrix() const {return vdata;}
	virtual int GetMatrixRows() const {return idata;}
	virtual int GetMatrixCols() const {return idata2;}

	virtual double GetMatrixVal(int row, int col) const //row and col are 1-based!!!
	{
		return vdata[(row-1)*idata2 + col-1];
	}

	virtual char* GetText() const {return tdata;}

	virtual void GetWCDriverAction(int& action, int& option, int& value) const
	{
		action = idata ;
		option = idata2;
		value  = idata3;
	}

	virtual void GetCompAction(int& action, int& option, int& value) const 
	{
		action = idata ;
		option = idata2;
		value  = idata3;
	}

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//set some data, but do not change type of ElementData
	virtual void SetText(const char* textdata) //only change value!
	{
		if (tdata != 0)
		{
			delete[] tdata;
			tdata = 0;
		}

		size_t length = strlen(textdata);
		tdata = new char[length + 1];
		strcpy(tdata, textdata);
	}

	virtual ElementDataContainer* GetEDC();
	virtual const ElementDataContainer* GetEDC() const;

	//virtual const ElementDataContainer* GetEDC() const;

	virtual void SetValuesInt() {type = type|2048;} //vector/matrix values are integer!

	virtual void SetMatrixSize(int rows, int cols)
	{
		if (vdatalen != 0 && vdata != 0)
		{
			delete[] vdata;
			vdata = 0;
		}

		idata = rows;
		idata2 = cols;
		vdatalen = rows*cols;

		if (vdatalen > 0)
			vdata = new double[vdatalen];
	}
	virtual void SetMatrixVal(int row, int col, double value) //row and col are 1-based!!!
	{
		vdata[(row-1)*idata2 + col-1] = value;
	}

private:
	char* dataname;	   //name of data unit: name that appears in dialog window
	int type;				   //type 1==int, 2==double, 4==vector, 8==text, 16==BOOL, 32==vector2D/3D,
	;									 //    64==elementtype, 128==actionDialogFunctionWCdriver, 256=actionCompFunction
	;									 //    512==checkgroup=idata2, 1024==matrix, 2048==vdata is integer, 4096==variable length vector, 
	;									 //    8192==ElementDataContainer, 16384==matrix with fixed number of columns
	int hastooltip;		 //data has additional tooltiptext
	char* tooltiptext; //tooltip text to be displayed in dialogs
	int locked;				 //data can not be changed
	int groupnum;			 //number of group where data belongs

	int idata;			   //integer data  (for action: action_number)
	int idata2;			   //integer data2 (for action: option)
	int idata3;			   //integer data3 (for action: value)
	double minval, maxval;
	bool oneLimit;
	// 4 cases for boundaries are possible:
	// oneLimit = 0: 
	//			minval<=maxval		upper and lower boundary are active:	[...]
	//			minval >maxval		no boundary is active:								 ...
	// onLimit = 1:
	//			minval<=maxval		only upper boundary is active:				 ...]
	//			minval >maxval		only lower boundary is active:				[...

	double* vdata;     //double/vector3D/vector data
	int vdatalen;			 //vector len
	char* tdata;		   //text data
	ElementDataContainer* edc;
};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ElementDataContainer
{
public:
	ElementDataContainer():elemdata()
	{
		type = 1; //standard==element data
		SetEDCWarningLevel(1);
		//defaultint = 0;
		//defaultdouble = 0.;
		//defaultchar = 0;
	}
	ElementDataContainer(const ElementDataContainer& ed):elemdata()
	{
		CopyFrom(ed);
	};
	ElementDataContainer& operator=(const ElementDataContainer& ed) 
	{
		if (this == &ed) {return *this;}
		CopyFrom(ed);
		return *this;
	}
	//To be overwritten in derived class:
	virtual ElementDataContainer* GetCopy() const
	{
		ElementDataContainer* ec = new ElementDataContainer();
		ec->CopyFrom(*this);
		return ec;
	}

	virtual void CopyFrom(const ElementDataContainer& e)
	{
		elemdata.SetLen(0);		
		for (int i=1; i <= e.elemdata.Length(); i++)
		{
			elemdata.Add(e.elemdata(i)->GetCopy());
		}
		type = e.type;
		flag_warning_level = e.flag_warning_level;
		//defaultint = e.defaultint;
		//defaultdouble = e.defaultdouble;
		//defaultchar = e.defaultchar;
	}

	virtual ~ElementDataContainer()
	{
		Reset();
	}

	virtual void Reset()
	{
		for (int i=1; i <= elemdata.Length(); i++)
		{
			if (elemdata(i) != 0) delete elemdata(i);
			elemdata(i) = 0;
		}
		elemdata.SetLen(0);
		type = 1;
		flag_warning_level = 1;
	}

	virtual int GetType() const {return type;}

	virtual int Length() const {return elemdata.Length();}

	//find in list of edc, not hierarchically; return index of ElementData
	virtual int Find(const char* name) const;
	virtual void FindAll(const char* name, TArray<int>& indices) const;// find all indices in an EDC of entries with same name.

	//find an element in tree, name is e.g. "root.child1.child2.leave"; the name of the found object is "leave"
	//return pointer to ElementData
	virtual const ElementData* TreeFind(const char* name) const;
	virtual ElementData* TreeFind(const char* name);

	//delete an element in tree, name is e.g. "root.child1.child2.leave"; the name of the deleted object is "leave"
	virtual void TreeDelete(const char* name);

	//find element, assume that it is double; if not exists use default value
	virtual double TreeGetDouble(const char* name, double default_val = 0.) const;
	//find element, assume that it is int; if not exists use default value
	virtual int TreeGetInt(const char* name, int default_val = 0) const;
	//find element, assume that it is bool; if not exists use default value (int and bool are internally the same!)
	virtual int TreeGetBool(const char* name, int default_val = 0) const;
	//find element, assume that it is char*; if not exists use default value
	virtual const char* TreeGetString(const char* name, char* default_val = "") const;
	//find element, assume that it is Vector3D; if not exists use 
	virtual int TreeGetVector3D(const char* name, double& vx,  double& vy,  double& vz) const;
	//find element, assume that it is Vector2D; if not exists use 
	virtual int TreeGetVector2D(const char* name, double& vx,  double& vy) const;
	//find element, assume that it is Vector; if not exists use 
	virtual int TreeGetVector(const char* name, double** v, int& len) const;
	//find element, assume that it is Matrix; if not exists use 
	virtual int TreeGetMatrix(const char* name, double** v, int& rows, int& cols) const;
	//find element, assume that it is Matrix; if not exists use 
	virtual int TreeGetVector2DList(const char* name, double** v, int& n) const;

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//add or replace element in EDC: double; e.g. name="tree1.leave2.data"
	virtual void TreeSetDouble(const char* name, double val);
	//add or replace element in EDC: int; e.g. name="tree1.leave2.data"
	virtual void TreeSetInt(const char* name, int val);
	//add or replace element in EDC: bool; e.g. name="tree1.leave2.data" (int and bool are internally the same!)
	virtual void TreeSetBool(const char* name, int val);
	//add or replace element in EDC: string; e.g. name="tree1.leave2.data"
	virtual void TreeSetString(const char* name, const char* val);

	//add or replace element in EDC: vector
	virtual void TreeSetVectorC(const char* name, double* v, int len, const char* comment);
	//add or replace element in EDC with COMMENT: double; e.g. name="tree1.leave2.data"
	virtual void TreeSetDoubleC(const char* name, double val, const char* comment);
	//add or replace element in EDC with COMMENT: double; e.g. name="tree1.leave2.data"
	virtual void TreeSetIntC(const char* name, int val, const char* comment);
	//add or replace element in EDC with COMMENT: double; e.g. name="tree1.leave2.data"
	virtual void TreeSetBoolC(const char* name, int val, const char* comment);
	//add or replace element in EDC with COMMENT: double; e.g. name="tree1.leave2.data"
	virtual void TreeSetStringC(const char* name, const char* val, const char* comment);
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//split hierarchical name into tree and elementname; e.g. name="tree1.leave2.data"; tree may also be empty ("")
	virtual void SplitIntoTreeAndElementName(const char* name, mystr& tree, mystr& elementname);

	//add an element to list of edc, not hierarchically!
	virtual int Add(const ElementData& ed);

	//replace edc-data (i.e., replace common, add new, and keep old edc-data) 
	void TreeReplaceEDCDataWith(const ElementDataContainer* edc_new, int warn_if_items_dont_exist=0, mystr warn_str="");
		
private:
	//full pathname (for warnings)
	mystr local_warn_tree;
	
	//called by public method TreeReplaceEDCDataWith(ElementDataContainer* edc_new, int warn_if_items_dont_exist)
	void TreeReplaceEDCDataWithRecursive(const ElementDataContainer* edc_new, int warn_if_items_dont_exist);

public:
	//set the value of private member local_warn_tree;
	void SetLocalWarnTree(mystr tree_str) {local_warn_tree = tree_str;}

	//add an element into the tree of edc; the tree is given e.g. as "root.child1"
	//virtual void TreeAdd(mystr tree, const ElementData& ed);
//	virtual void TreeAdd(mystr tree, const ElementData& ed);
	virtual void TreeAdd(const char* tree, const ElementData& ed);
	
	//set an element in the tree of edc; if the tree does not exist, build the tree and insert data; if the data exists: replace the data
	void TreeSet(const char* tree, const ElementData& ed, int warn_if_items_dont_exist=0,int copy_ToolTipText=1);

	virtual const ElementData& Get(int i) const {return *(elemdata(i));} //get element number i
	virtual ElementData& Get(int i) {return *(elemdata(i));}//get element number i
	virtual ElementData* GetPtr(int i) {return elemdata(i);}//get element* number i
	virtual const ElementData* GetPtr(int i) const {return elemdata(i);}//get element* number i

	virtual const ElementData& Last() const {return *(elemdata(Length()));} //get last element
	virtual ElementData& Last() {return *(elemdata(Length()));}//get last element

	virtual void Delete(int i) 
	{
		delete elemdata(i);
		for (int j=i; j < elemdata.Length(); j++)
		{
			elemdata(j) = elemdata(j+1);
		}
		elemdata.SetLen(elemdata.Length()-1);
	}

	int EDCWarningLevel() const
	{
		if (flag_warning_level > 0) 
		{
			flag_warning_level--;
			return 1; // popup
		}
		else if (flag_warning_level == 0)
		{
			return 0; // text output
		}
		else //if (flag_warning_level < 0)
		{
			return -1; // no warning
		}
	}
	int GetEDCWarningLevel() { return flag_warning_level; }
	void SetEDCWarningLevel(int flag) { flag_warning_level = flag; }

	//void PrintEDCRecursive(mystr filename);			//$ DR 2012-06-27
	void PrintEDCRecursive(ofstream &outstr);			//$ DR 2012-06-27

private:
	TArray<ElementData*> elemdata;
	int type; //1==element data, 2==node data, 4==MBS data, 8==element initvector, 16==element actual DOF state
	//          32==GeomObject, 64==force

	mutable int flag_warning_level; // can be changed via EDCWarningLevel() - !const!
	//int defaultint;				//only default return value, don't use!!!!!
	//double defaultdouble;	//only default return value, don't use!!!!!
	//char* defaultchar;		//only default return value, don't use!!!!!

};







#endif 
//ELEMENTDATA__H

