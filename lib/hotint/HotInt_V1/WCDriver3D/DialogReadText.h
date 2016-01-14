//#**************************************************************
//# filename:             DialogReadText.h
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
 

#if !defined(AFX_DIALOGREADTEXT_H__947C227F_731A_4438_8B49_82AC09435E8B__INCLUDED_)
#define AFX_DIALOGREADTEXT_H__947C227F_731A_4438_8B49_82AC09435E8B__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// DialogReadText.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CDialogReadText dialog

class CDialogReadText : public CDialog
{
// Construction
public:
	CDialogReadText(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(CDialogReadText)
	enum { IDD = IDD_DIALOG_READ_TEXT };
	CString	m_Text;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CDialogReadText)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CDialogReadText)
		// NOTE: the ClassWizard will add member functions here
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_DIALOGREADTEXT_H__947C227F_731A_4438_8B49_82AC09435E8B__INCLUDED_)
