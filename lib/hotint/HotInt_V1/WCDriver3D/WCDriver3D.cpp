//#**************************************************************
//# filename:             WCDriver3D.cpp
//#
//# author:               Gerstmayr, Vetyukov
//#
//# generated:						
//# description:          Defines the class behaviors for the application.
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
 


#include "stdafx.h"
#include "mbs_interface.h"
#include "WCDriver3DDlg.h"
#include "WCDriver3D.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CWCDriver3DApp

BEGIN_MESSAGE_MAP(CWCDriver3DApp, CWinApp)
	//{{AFX_MSG_MAP(CWCDriver3DApp)
		// NOTE - the ClassWizard will add and remove mapping macros here.
		//    DO NOT EDIT what you see in these blocks of generated code!
	//}}AFX_MSG
	ON_COMMAND(ID_HELP, CWinApp::OnHelp)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CWCDriver3DApp construction

void TIMBSWarningHandle(const char* warn, int use_instant_message_text)
{
	AfxMessageBox(warn);
}


CWCDriver3DApp::CWCDriver3DApp()
{
	// TODO: add construction code here,
	// Place all significant initialization in InitInstance
}

/////////////////////////////////////////////////////////////////////////////
// The one and only CWCDriver3DApp object

CWCDriver3DApp theApp;

/////////////////////////////////////////////////////////////////////////////
// CWCDriver3DApp initialization

BOOL CWCDriver3DApp::InitInstance()
{
	AfxOleInit();
	// Standard initialization
	// If you are not using these features and wish to reduce the size
	//  of your final executable, you should remove from the following
	//  the specific initialization routines you do not need.

	//!AD: 03-01-2013: test to add hotkeys
	m_haccel = LoadAccelerators(AfxGetInstanceHandle(), MAKEINTRESOURCE(IDR_ACCELERATOR1));

	//CWCDriver3DDlg dlg;
	dlg = new CWCDriver3DDlg();
	m_pMainWnd = dlg;
	int nResponse = dlg->DoModal();
	if (nResponse == IDOK)
	{
		// TODO: Place code here to handle when the dialog is
		//  dismissed with OK
	}
	else if (nResponse == IDCANCEL)
	{
		// TODO: Place code here to handle when the dialog is
		//  dismissed with Cancel
	}

	// Since the dialog has been closed, return FALSE so that we exit the
	//  application, rather than start the application's message pump.
	return FALSE;
}

int CWCDriver3DApp::ExitInstance()
{
	delete dlg;
	return 0;
}

//!AD: 03-01-2013: test to add hotkeys
BOOL CWCDriver3DApp::ProcessMessageFilter(int code, LPMSG lpMsg) 
{
    if(m_haccel)
    {
        if (::TranslateAccelerator(m_pMainWnd->m_hWnd, m_haccel, lpMsg)) 
            return(TRUE);
    }
	
    return CWinApp::ProcessMessageFilter(code, lpMsg);
}