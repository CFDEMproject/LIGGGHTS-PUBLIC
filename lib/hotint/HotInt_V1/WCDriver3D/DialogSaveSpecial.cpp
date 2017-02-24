//#**************************************************************
//# filename:             DialogSaveSpecial.cpp
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
#include "DialogSaveSpecial.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CDialogSaveSpecial dialog


CDialogSaveSpecial::CDialogSaveSpecial(int NumberOfDataUnits_,CWnd* pParent /*=NULL*/) :
CDialog(CDialogSaveSpecial::IDD, pParent),
NumberOfDataUnits(NumberOfDataUnits_)
{
	//{{AFX_DATA_INIT(CDialogSaveSpecial)
	m_FirstDataUnit = 1;
	m_LastDataUnit = NumberOfDataUnits;
	m_SaveEachOf = 1;
	//}}AFX_DATA_INIT
}


void CDialogSaveSpecial::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_EDIT_TOTAL_NUMBER, NumberOfDataUnits);
	//{{AFX_DATA_MAP(CDialogSaveSpecial)
	DDX_Text(pDX, IDC_EDIT_FIRST_UNIT, m_FirstDataUnit);
	DDV_MinMaxInt(pDX, m_FirstDataUnit, 1, NumberOfDataUnits);
	DDX_Text(pDX, IDC_EDIT_LAST_UNIT, m_LastDataUnit);
	DDV_MinMaxInt(pDX, m_LastDataUnit, m_FirstDataUnit, NumberOfDataUnits);
	DDX_Text(pDX, IDC_EDIT_SAVE_EACH_OF, m_SaveEachOf);
	DDV_MinMaxInt(pDX, m_SaveEachOf, 1, NumberOfDataUnits);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CDialogSaveSpecial, CDialog)
	//{{AFX_MSG_MAP(CDialogSaveSpecial)
		// NOTE: the ClassWizard will add message map macros here
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CDialogSaveSpecial message handlers
