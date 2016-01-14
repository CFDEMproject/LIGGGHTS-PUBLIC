//#**************************************************************
//# filename:             OutputDialog.h
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

// COutputDialog-Dialogfeld

class COutputDialog : public CDialog
{
	DECLARE_DYNAMIC(COutputDialog)

public:
	COutputDialog(CWnd* pParent = NULL);   // Standardkonstruktor
	virtual ~COutputDialog();

	CRect DialogRect;
	int storedwidth;

// Dialogfelddaten
	enum { IDD = IDD_OUTPUTDIALOG };

	void Create(CWnd * pParent);
	void OnCancel();
	void OnClose();

	void UpdateText(const CString& str_full);
	void AppendText(const CString& str_append);
	void ReplaceLastLine(const CString& str_replace);
	int ReachedLimitText(const CString& str_append);
	void DoubleLimitText();

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV-Unterstützung

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnSize(UINT nType, int cx, int cy);

public:
	CEdit m_control_outputcedit;
};
