//#**************************************************************
//# filename:             CustomEditDialog.h
//#
//# author:               Gerstmayr, Vetyukov
//#
//# generated:						
//# description:          
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
 

#pragma once

class CustomEditDialog;
class CTreeViewCustomEditDialog;

const int CEedit_height				= 18;
const int CEedit_width				= 80;
const int CEstatictext_width	= 200;
const int CEbutton_width			= 130; //3*CEbutton_width < 3*CEedit_width + CEstatictext_width
const int CEedit_vspace				= 2;
const int CEedit_hspace				= 5; 
const int CEedit_xoff					= 8;
const int CEedit_yoff					= 8;
const int CEbutton_offy				= 30;

const int CEid_offset1				= 10000;
const int CEid_offset2				= 11000;
const int CEid_offset_group		= 12000;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// MyEdit

class MyEdit : public CEdit
{
	DECLARE_DYNAMIC(MyEdit)

public:
	MyEdit();
	virtual ~MyEdit();
	virtual void SetMinMax(double min, double max, bool justOneLimit=0) {minval = min; maxval = max; oneLimit=justOneLimit;};

	double minval, maxval;
	bool oneLimit;
	// 4 cases for boundaries are possible:
	// oneLimit = 0: 
	//			minval<=maxval		upper and lower boundary are active:	[...]
	//			minval >maxval		no boundary is active:								 ...
	// onLimit = 1:
	//			minval<=maxval		only upper boundary is active:				 ...]
	//			minval >maxval		only lower boundary is active:				[...

protected:
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnKillFocus(CWnd* pNewWnd);
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// MyButton

class MyCButton : public CButton
{
	DECLARE_DYNAMIC(MyCButton)

public:
	MyCButton();
	virtual ~MyCButton();
	virtual void SetCustomEditDialog(CustomEditDialog* cedI, int edcnumI);
	virtual void SetCustomEditDialog(CTreeViewCustomEditDialog* cedI, int edcnumI);
	
	CustomEditDialog* ced;
	CTreeViewCustomEditDialog* tvced;
	int edcnum;

protected:
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClicked();
};




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// CustomEditDialog-Dialogfeld

class CustomEditDialog : public CDialog
{
	DECLARE_DYNAMIC(CustomEditDialog)

public:
	CustomEditDialog(CWnd* pParent = NULL);   // Standardkonstruktor
	virtual ~CustomEditDialog();

// Dialogfelddaten
	enum { IDD = IDD_CUSTOMEDITDIALOG };

	//functions not added by windows:
	virtual void Create(CWnd * pParent);
	virtual void SetWCDI(WCDInterface * pWCDI_) { pWCDI = pWCDI_; }
	virtual void SetGLDrawWnd(CGLDrawWnd * pGLDrawWnd_) { pGLDrawWnd = pGLDrawWnd_; }
	virtual void WriteData(); //write data back to ElementDataContainer* edc
	virtual void SetElementDataContainer(ElementDataContainer* edcI) {edc = edcI;}
	virtual void SetDialogName(const CString& name) {dialogname = name;}
	virtual void ButtonClicked(int edcnum);
	virtual void SetDeleteButton(int flag) {add_deletebutton = flag;};
	virtual int GetDeleteFlag() const {return deleteflag;}; //show parent that element/force shall be deleted!

	virtual void AddToolTipText(const char* str, int itemID);
	virtual int FindToolTipID(int itemID)
	{
		for (int i=1; i <= tooltiptextID.Length(); i++)
		{
			if (tooltiptextID(i) == itemID) return i;
		}
		return 0;
	}

	virtual void CreateDialogItems();
	virtual void ClearDialogItems();
	virtual void Delete();

	int add_deletebutton;
	int deleteflag;

//tooltips:
	
private:
	WCDInterface* pWCDI;
	CGLDrawWnd* pGLDrawWnd;
	ElementDataContainer* edc;
	ElementDataContainer edc_original;
	CString dialogname;

	//Test:
	TArray<MyEdit*> cedit;
	TArray<CStatic*> cstatic;
	TArray<CButton*> cbutton;
	TArray<int> ceditID;
	TArray<int> cstaticID;
	TArray<int> cbuttonID;
	TArray<int> groupbutton; //button number for group

	//add tooltip texts for IDs, only works for CButton and CEdit fields!
	TArray<char*> tooltiptext;
	TArray<int> tooltiptextID;

	

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV-Unterstützung
	BOOL OnToolTipNotify(UINT id, NMHDR *pNMHDR, LRESULT *pResult);
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedOk();
	afx_msg void OnBnClickedCancel();
	virtual BOOL OnInitDialog();
};


#pragma once


