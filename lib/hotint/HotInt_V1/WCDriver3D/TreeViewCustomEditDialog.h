//#**************************************************************
//# filename:             TreeViewCustomEditDialog.h
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
#include "afxcmn.h"

// this is a helper class: a background for the dynamically created controls should support tooltips
class CStaticBackground : public CStatic
{
	DECLARE_DYNAMIC(CStaticBackground)

protected:
	DECLARE_MESSAGE_MAP()
	BOOL OnToolTipNotify(UINT id, NMHDR *pNMHDR, LRESULT *pResult);
};

// CTreeViewCustomEditDialog

class CTreeViewCustomEditDialog : public CDialog
{
	DECLARE_DYNAMIC(CTreeViewCustomEditDialog)
	int nVerticalScroll;
	int nVerticalSpaceForTheControls;
	int nCurrentBackgroundHeight;
	bool flagDialogInitialized;
	CStaticBackground wndStaticBackground;

	// matching between items in the tree control, edc's and their paths
	struct TreeItemData
	{
		HTREEITEM hti;
		ElementDataContainer * edc;
		char * path;
	};
	TArray<TreeItemData> treeItemsData;
	TreeItemData currentTreeItemData;

	// data for the dynamically created controls
	struct CreatedControlData
	{
		int id;
		char * toolTipText;
		char * path;
	};
	TArray<CreatedControlData> createdControlsData;

public:
	CTreeViewCustomEditDialog(CWnd* pParent = NULL);   // Standardkonstruktor
	virtual ~CTreeViewCustomEditDialog();

	enum { IDD = IDD_TREEVIEWCUSTOMEDITDIALOG };

	virtual void Create(CWnd * pParent);
	virtual void SetWCDI(WCDInterface * pWCDI_) { pWCDI = pWCDI_; }
	virtual void SetGLDrawWnd(CGLDrawWnd * pGLDrawWnd_) { pGLDrawWnd = pGLDrawWnd_; }
	virtual void SetParser(class CEDCParser* pParser_) { pParser = pParser_; }         //AD: Added for ClipBoardOperations 2013-02-21

	virtual void WriteData(); //write data back to ElementDataContainer* edc
	virtual void SetElementDataContainer(ElementDataContainer* edcI) {edc = edcI;}
	virtual void SetUseApply(int useapplyI, int callcompfn_actionI, int callcompfn_optionI, 
		int callcompfn_valueI) 
	{
		useapply = useapplyI;
		callcompfn_action = callcompfn_actionI;
		callcompfn_option = callcompfn_optionI;
		callcompfn_value = callcompfn_valueI;
	}
	virtual void SetDialogName(const CString& name) {dialogname = name;}
	virtual void ButtonClicked(int edcnum);
	virtual void SetDeleteButton(int flag) {add_deletebutton = flag;};
	virtual int GetDeleteFlag() const {return deleteflag;}; //show parent that element/force shall be deleted!

	void AddCreatedControlData(const char* toolTipText, const char * path, int itemID);
	char * FindToolTipText(int itemID);

	virtual ElementDataContainer* GetRootEDC() {return edc;}
	virtual ElementDataContainer* GetCurrentEDC() {return currentTreeItemData.edc;}

	virtual void InsertTreeItem(TreeItemData tid);

	virtual void UpdateDialogItems();
	virtual void ClearDialogItems();
	virtual void Delete();

	void PositionControls();		// places the control to appropriate positions in the dialog
	void UpdateScrollBarInfo();

	int add_deletebutton;
	int deleteflag;
	
private:
	WCDInterface* pWCDI;
	CGLDrawWnd* pGLDrawWnd;
	class CEDCParser* pParser;         //AD: Added for ClipBoardOperations 2013-02-21
	ElementDataContainer* edc;

	int useapply;
	int	callcompfn_action;
	int	callcompfn_option;
	int	callcompfn_value;

	CString dialogname;

	//Test:
	TArray<MyEdit*> cedit;
	TArray<CStatic*> cstatic;
	TArray<CButton*> cbutton;
	TArray<int> groupbutton; //button number for group

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV-Unterstützung
	BOOL OnToolTipNotify(UINT id, NMHDR *pNMHDR, LRESULT *pResult);
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedOk();
	afx_msg void OnBnClickedCancel();
	virtual BOOL OnInitDialog();
public:
	CTreeCtrl m_treectrl_control;
public:
	afx_msg void OnTvnSelchangedTreecustomedit(NMHDR *pNMHDR, LRESULT *pResult);
public:
	afx_msg void OnBnClickedApply();
	afx_msg void OnVScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar);
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg void OnContextMenu(CWnd* /*pWnd*/, CPoint /*point*/);
	afx_msg void OnCopyToClipboard();
};


#pragma once
