//#**************************************************************
//# filename:             DialogProgressBar.cpp
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
#include "DialogProgressBar.h"

// CDialogProgressBar-Dialogfeld

IMPLEMENT_DYNAMIC(CDialogProgressBar, CDialog)

CDialogProgressBar::CDialogProgressBar(CWnd* pParent /*=NULL*/)
	: CDialog(CDialogProgressBar::IDD, pParent)
{
  
}

CDialogProgressBar::~CDialogProgressBar()
{

}

void CDialogProgressBar::Reset()
{
	CProgressCtrl* pPC = (CProgressCtrl*) GetDlgItem(IDC_PROGRESS);
	pPC->SetRange(0,100);
	pPC->SetPos(50);

  CEdit* pE = (CEdit*) GetDlgItem(IDC_EDIT_OVER);
	pE->SetWindowTextA("not linked to a process");
	pE = (CEdit*) GetDlgItem(IDC_EDIT_UNDER);
	pE->SetWindowTextA("...");
}
void CDialogProgressBar::SetTotalTicks(int totalticks)
{
	m_progress_control.SetRange(0,totalticks);
}

void CDialogProgressBar::SetActualTicks(int actticks)
{
	m_progress_control.SetPos(actticks);
}

void CDialogProgressBar::SetCaptionText(mystr& caption)
{
	SetWindowTextA(caption.c_str());
}

void CDialogProgressBar::SetTextOver(mystr& over)
{
	CEdit* pEdit = (CEdit*) GetDlgItem(IDC_EDIT_OVER);
	pEdit->SetWindowTextA(over.c_str());
//	m_textover_control.SetWindowText(over.c_str());
}
void CDialogProgressBar::SetTextUnder(mystr& under)
{
	m_textunder_control.SetWindowText(under.c_str());
}

void CDialogProgressBar::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_PROGRESS, m_progress_control);
	DDX_Control(pDX, IDC_EDIT_OVER, m_textover_control);
	DDX_Control(pDX, IDC_EDIT_UNDER, m_textunder_control);
}


BEGIN_MESSAGE_MAP(CDialogProgressBar, CDialog)
	ON_WM_CLOSE()
END_MESSAGE_MAP()


// CDialogProgressBar-Meldungshandler

void CDialogProgressBar::OnClose()
{
	DestroyWindow();
}