//#**************************************************************
//# filename:             DialogReadText.cpp 
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
 

#include "stdafx.h"
#include "WCDriver3DDlg.h"
#include "WCDriver3D.h"
#include "DialogReadText.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CDialogReadText dialog


CDialogReadText::CDialogReadText(CWnd* pParent /*=NULL*/)
	: CDialog(CDialogReadText::IDD, pParent)
{
	//{{AFX_DATA_INIT(CDialogReadText)
	m_Text = _T("");
	//}}AFX_DATA_INIT
}


void CDialogReadText::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CDialogReadText)
	DDX_Text(pDX, IDC_EDIT_TEXT, m_Text);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CDialogReadText, CDialog)
	//{{AFX_MSG_MAP(CDialogReadText)
		// NOTE: the ClassWizard will add message map macros here
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CDialogReadText message handlers
