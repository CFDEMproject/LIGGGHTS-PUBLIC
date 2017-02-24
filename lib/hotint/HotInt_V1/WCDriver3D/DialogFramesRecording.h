//#**************************************************************
//# filename:             DialogFramesRecording.h
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


#if !defined(AFX_DIALOGFRAMESRECORDING_H__830C73BC_4C26_41D3_822F_EB2EC1BFC852__INCLUDED_)
#define AFX_DIALOGFRAMESRECORDING_H__830C73BC_4C26_41D3_822F_EB2EC1BFC852__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// DialogFramesRecording.h : header file
//

#include "resource.h"
#include "afxwin.h"

/////////////////////////////////////////////////////////////////////////////
// CDialogFramesRecording dialog

class CDialogFramesRecording : public CDialog
{
	BOOL bCurrentCheckRecordFrames;
	void UpdateEnableControls();

// Construction
public:
	CDialogFramesRecording(CWnd* pParent = NULL);   // standard constructor
	void SetWindowDimensions(const CRect & r)
	{
		m_nWindowSizeX = r.Width();
		m_nWindowSizeY = r.Height();
	}

	//void Serialize(CArchive & ar); //AD removed 2012-11-09

	//convert some internal settings (window position, open dialog windows, etc.) to EDC
	virtual void Configuration2EDC(ElementDataContainer& edc);
	//convert EDC to some internal settings (window position, open dialog windows, etc.)
	virtual void EDC2Configuration(const ElementDataContainer& edc);

private:  
	WCDInterface* pWCDI;

public:
	void SetWCDI(WCDInterface * pWCDI_) { pWCDI = pWCDI_; }
	WCDInterface* GetWCDI() { return pWCDI; }
	virtual void LoadData();  //Get data from WCDinterface
	virtual void WriteData(); //Put data to WCDinterface

// Dialog Data
	//{{AFX_DATA(CDialogFramesRecording)
	enum { IDD = IDD_DIALOG_FRAMES_RECORDING };
	BOOL m_bCheckRecordFrames;                //B163
	//CString	m_strPathToImageFiles;            //T102 // AD: old code - remove this to find all access to the variable
	BOOL m_bCheckIncludeOutputWindow;         //B168

	CString	m_strPathToSingleImageFiles;      //T111
  CString m_strSingleImageFileName;					//T112
	CString	m_strPathToVideoImageFiles;				//T113
  CString m_strVideoImageFileName;					//T114
	int m_radio_fileformat;                   //format of the saved files

	int		m_nRecordEachFrameOf;               //I166
	int		m_nWindowSizeX;
	int		m_nWindowSizeY;
	int		m_nFrameCounter;                    
	BOOL	m_bShowFrameNumbers;                //B165
	//BOOL	m_bProcessImage;                    //B164 // AD: no longer needed ? file format is picked via radio controöl ?
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CDialogFramesRecording)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CDialogFramesRecording)
	afx_msg void OnCheckRecordFrames();
	virtual BOOL OnInitDialog();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedOk();
	afx_msg void OnBnClickedButtonBrowsesingle();
	afx_msg void OnBnClickedButtonBrowsevideo();
private:
	void OnBrowse(CString& pathstring, const char* title); // returns the full path from the BrowseInfoWindow. CString& version possible...

public:
	afx_msg void OnBnClickedRadioFormat();
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_DIALOGFRAMESRECORDING_H__830C73BC_4C26_41D3_822F_EB2EC1BFC852__INCLUDED_)
