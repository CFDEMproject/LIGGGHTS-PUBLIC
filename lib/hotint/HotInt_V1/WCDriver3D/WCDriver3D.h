//#**************************************************************
//# filename:             WCDriver3D.h
//#
//# author:               Gerstmayr, Vetyukov
//#
//# generated:						
//# description:        main header file for the WCDriver3D application  
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
 


#if !defined(AFX_WCDriver3D_H__CD254043_2EFD_4173_A955_7E28AAE1FF5E__INCLUDED_)
#define AFX_WCDriver3D_H__CD254043_2EFD_4173_A955_7E28AAE1FF5E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#ifndef __AFXWIN_H__
	#error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"		// main symbols

/////////////////////////////////////////////////////////////////////////////
// CWCDriver3DApp:
// See WCDriver3D.cpp for the implementation of this class
//
void TIMBSWarningHandle(const char* warn, int use_instant_message_text=0);

class CWCDriver3DApp : public CWinApp
{
public:
	CWCDriver3DApp();
	CWCDriver3DDlg* dlg;

	//!AD: 03-01-2013: test to add hotkeys [
	HACCEL m_haccel;
	virtual BOOL ProcessMessageFilter(int code, LPMSG lpMsg);
	//!AD: ]

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CWCDriver3DApp)
	public:
	virtual BOOL InitInstance();
	virtual int ExitInstance();
	//}}AFX_VIRTUAL

// Implementation

	//{{AFX_MSG(CWCDriver3DApp)
		// NOTE - the ClassWizard will add and remove member functions here.
		//    DO NOT EDIT what you see in these blocks of generated code !
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};


/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_WCDriver3D_H__CD254043_2EFD_4173_A955_7E28AAE1FF5E__INCLUDED_)
