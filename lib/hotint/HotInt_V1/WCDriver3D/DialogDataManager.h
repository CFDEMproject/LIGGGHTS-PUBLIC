//#**************************************************************
//# filename:             DialogDataManager.h
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
 

#if !defined(AFX_DIALOGDATAMANAGER_H__B6E6D790_FAD8_43A9_9639_A4017D526EB7__INCLUDED_)
#define AFX_DIALOGDATAMANAGER_H__B6E6D790_FAD8_43A9_9639_A4017D526EB7__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// DialogDataManager.h : header file
//

#include <afxtempl.h>

#include "resource.h"
#include "DataStorage.h"
#include "gldrawwnd.h"

/////////////////////////////////////////////////////////////////////////////
// CDialogDataManager dialog

#include "..\WorkingModule\WinCompDriverInterface.h"

#include "dialogsavespecial.h"

// this dialog operates with the solution data

class CDialogDataManager : public CDialog
{
	CArray<DataStorage,DataStorage&> DataStorageArray;
	CArray<double> TimePointArray;
	WCDInterface * pWCDI;
	CWnd * pGLDrawWnd;

	CString ProgramDirectory;	// directory where the program is running

	void ActualizeTimePoint();

	bool bAnimating;

	bool bRetrievingData;		// a flag, meaning that the status text should not be updated yet

	void SaveFile(CDialogSaveSpecial & info);

// Construction
public:
	CDialogDataManager();   // standard constructor
	void Create(CWnd * pParent);
	void SetWCDI(WCDInterface * pWCDI_) { pWCDI = pWCDI_; }
	void SetGLDrawWnd(CWnd * pGLDrawWnd_) { pGLDrawWnd = pGLDrawWnd_; }

	void AddEntry(DataStorage &, double);
	void RemoveAll();

	bool RetrievingData() { return bRetrievingData; }

	// the actual implementation of	WCDInterface::ComputationFeedBack::RedrawNewResults()
	bool RedrawNewResults();
	void CheckScreenLocation();
	void SetScrollPosToLastTimePoint();

	CString GetHotintDataVersionString() const;

// Dialog Data
	//{{AFX_DATA(CDialogDataManager)
	enum { IDD = IDD_DIALOG_DATA_MANAGER };
	CScrollBar	m_ScrollTimePoint;
	int		m_TimePointNumber;
	CString	m_strDataRequest;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CDialogDataManager)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	void OnOK();
	void OnCancel();

	// Generated message map functions
	//{{AFX_MSG(CDialogDataManager)
	afx_msg void OnButtonAnimate();
	afx_msg void OnButtonLoadFile();
	afx_msg void OnButtonSaveFile();
	afx_msg void OnKillfocusEditTimePointNumber();
	afx_msg void OnClose();
	afx_msg void OnHScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar);
	virtual BOOL OnInitDialog();
	afx_msg void OnTimer(UINT_PTR nIDEvent);
	afx_msg void OnButtonSubmitPrintDataCommand();
	afx_msg void OnButtonRetrieveData();
	afx_msg void OnButtonSaveFileSpecial();
	//}}AFX_MSG
	afx_msg LRESULT OnUpdate(WPARAM, LPARAM);
	DECLARE_MESSAGE_MAP()
public:
	double m_animation_delay;
	afx_msg void OnEnKillfocusEditAnimationDelay();
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_DIALOGDATAMANAGER_H__B6E6D790_FAD8_43A9_9639_A4017D526EB7__INCLUDED_)
