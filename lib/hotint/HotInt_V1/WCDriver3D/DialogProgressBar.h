//#**************************************************************
//# filename:             DialogProgressBar.h
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
#include "resource.h"
#include "afxcmn.h"
#include "afxwin.h"
#include "mystring.h"

// CDialogProgressBar-Dialogfeld

class CDialogProgressBar : public CDialog
{
	DECLARE_DYNAMIC(CDialogProgressBar)

public:
	CDialogProgressBar(CWnd* pParent = NULL);   // Standardkonstruktor
	virtual ~CDialogProgressBar();

public:
	void Reset();
	void OnClose();
	void SetTotalTicks(int totalticks);
	void SetActualTicks(int actticke);
	void SetCaptionText(mystr& caption);
	void SetTextOver(mystr& over);
	void SetTextUnder(mystr& under);

	void SetWCDI(WCDInterface * pWCDI_) { pWCDI = pWCDI_; }
	void SetGLDrawWnd(CGLDrawWnd * pGLDrawWnd_) { pGLDrawWnd = pGLDrawWnd_; }
	void SetWCDDlg(CWCDriver3DDlg* adial) {WCDDlg = adial;}

// Dialogfelddaten
	enum { IDD = IDD_DIALOG_PROGRESSBAR };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV-Unterstützung

	DECLARE_MESSAGE_MAP()

public:
	CProgressCtrl m_progress_control;
	CEdit m_textover_control;
	CEdit m_textunder_control;

private:  
	WCDInterface* pWCDI;
	CGLDrawWnd* pGLDrawWnd;
	CWCDriver3DDlg* WCDDlg;
};
