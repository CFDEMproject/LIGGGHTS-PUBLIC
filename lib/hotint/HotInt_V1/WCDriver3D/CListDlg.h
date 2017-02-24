//#**************************************************************
//# filename:             ClistDlg.h
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
#include "afxwin.h"


// CustomListDialog-Dialogfeld

class CustomListDialog : public CDialog
{
	DECLARE_DYNAMIC(CustomListDialog)

public:
	CustomListDialog(CWnd* pParent = NULL);   // Standardkonstruktor
	virtual ~CustomListDialog();
	void SetUseEditItem(int flag) {use_edit_item = flag;}
	void InitList(const char* dialognameI, const char* selectinfotextI, int use_custom_check1I=0, const CString& checkstring1I="");
	void AddString(const CString& str, int value = -1);
	void SetPixelSize(int width, int height) {wwindow = width; hwindow=height;}

	CString dialogname;
	CString selectinfotext;
	CString checkstring1;
	int item_selected;
	int custom_check1;
	int use_custom_check1;
	int use_edit_item;
	TArray<CString*> items;
	TArray<int> values;

	int hwindow; //increase window height
	int wwindow; //increase window width
	int iscancel;

// Dialogfelddaten
	enum { IDD = IDD_CUSTOMLISTDIALOG };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV-Unterstützung

	DECLARE_MESSAGE_MAP()
public:
	BOOL m_check_customlist;
	CListBox listctrl_custom;
	virtual BOOL OnInitDialog();
	afx_msg void OnBnClickedCancel();
	afx_msg void OnBnClickedOk();
	afx_msg void OnLbnDblclkCustomList();
	afx_msg void OnEnKillfocusEdit2();
	afx_msg void OnLbnSelchangeCustomList();
	CString m_edit_text;
};
