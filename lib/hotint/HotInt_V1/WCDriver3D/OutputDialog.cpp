//#**************************************************************
//# filename:             OutputDialog.cpp
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
#include "OutputDialog.h"

// COutputDialog-Dialogfeld

IMPLEMENT_DYNAMIC(COutputDialog, CDialog)
COutputDialog::COutputDialog(CWnd* pParent /*=NULL*/)
	: CDialog(COutputDialog::IDD, pParent)
{
	storedwidth = 200;
}

COutputDialog::~COutputDialog()
{
}

void COutputDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_EDIT_COMPUTATION_OUTPUT, m_control_outputcedit);
}

void COutputDialog::Create(CWnd * pParent)
{
	CDialog::Create(IDD,pParent);

	CRect r;
	pParent->GetWindowRect(&r);
	CRect rw;
	GetWindowRect(&rw);

	CPoint p(r.left,r.bottom);
	if (p.x < rw.Width()) p.x = rw.Width();

//	ClientToScreen(&r);
	rw.top = p.y - r.Height();
	rw.left = p.x - storedwidth;
	//rw.right = p.x;
	rw.right = rw.left+storedwidth;
	rw.bottom = p.y;
	MoveWindow(rw,FALSE);
	ShowWindow(SW_SHOW);

	GetClientRect(&DialogRect);
}

void COutputDialog::OnClose() 
{
	DestroyWindow();
}

void COutputDialog::OnCancel() 
{
	DestroyWindow();
}

// checks if textlimit is reached 
int COutputDialog::ReachedLimitText(const CString& str_append)
{
	int limittext = m_control_outputcedit.GetLimitText();
	int cedittext = m_control_outputcedit.GetWindowTextLengthA();
	int appendlen = str_append.GetLength();
	if ((cedittext+appendlen) > limittext ) 
		return 1;
	else return 0;

}
void COutputDialog::DoubleLimitText()
{
  m_control_outputcedit.LimitText(m_control_outputcedit.GetLimitText()*2);
}

// OLD VERSION
/*
void COutputDialog::UpdateText(const CString& str_full)
{
	m_edit_computationoutput = str_full;
	UpdateData(FALSE);
	CEdit * pEdit = (CEdit*)GetDlgItem(IDC_EDIT_COMPUTATION_OUTPUT);
	pEdit->SetSel(100000000,100000000,FALSE);
}
*/

// NEW VERSION (AD): this is a full update!
void COutputDialog::UpdateText(const CString& str_full)
{
	m_control_outputcedit.SetWindowTextA(str_full);	// replaces entire text
	int pos = m_control_outputcedit.LineIndex(m_control_outputcedit.GetLineCount()-1);					// index number of character at beginning of last line
	m_control_outputcedit.SetSel(pos,pos,FALSE); //this command is quite slow!!!

	//Mögliche einfache Verbesserung:
	//mitzählen der Zeilen auf WCDriver-Seite (wird ja ohnehin abgefangen), merken des Index auf die 100-letzte Zeile
	//Übergeben der 100-letzten Zeilen an COutputDialog
	//eventuell könnte der gesamte String von WCDriver im Falle des Aktivieren des Scrollbars gesetzt werden

	//int scroll_limit = m_control_outputcedit.GetScrollLimit(1);
	//m_control_outputcedit.SetScrollPos(1,scroll_limit); // does not work!!!
	//char intstr[50]; sprintf_s(intstr,"%d",scroll_limit);
	//m_control_outputcedit.SetWindowTextA(str_full+CString(", scroll limit=")+CString(intstr)+CString("\n"));	// replaces entire text
}
// NEW VERSION (AD): appends text
void COutputDialog::AppendText(const CString& str_append)
{
	if(ReachedLimitText(str_append)) DoubleLimitText();
	m_control_outputcedit.SetSel(-1,-1,FALSE);  
	m_control_outputcedit.ReplaceSel(str_append);
}
// NEW VERSION (AD): replaces text of last line
void COutputDialog::ReplaceLastLine(const CString& str_replace)
{
	int pos = m_control_outputcedit.LineIndex(m_control_outputcedit.GetLineCount()-1);
	int endpos = m_control_outputcedit.GetWindowTextLengthA();	
	m_control_outputcedit.SetSel(pos,endpos,FALSE); 
	m_control_outputcedit.ReplaceSel(str_replace);
}

BEGIN_MESSAGE_MAP(COutputDialog, CDialog)
	ON_WM_SIZE()
END_MESSAGE_MAP()


// COutputDialog-Meldungshandler

void COutputDialog::OnSize(UINT nType, int cx, int cy)
{
	if(cx > 50 && cy > 50)
	{
		int vertical_shift = DialogRect.bottom-cy;
		int horizontal_shift = DialogRect.right-cx;

		DialogRect.right = cx;
		DialogRect.bottom = cy;

		CRect r;

		CWnd * pEditConsole = GetDlgItem(IDC_EDIT_COMPUTATION_OUTPUT);
		if(pEditConsole)
		{
			pEditConsole->GetWindowRect(&r);
			ScreenToClient(&r);
			pEditConsole->MoveWindow(r.left,r.top,r.Width()-horizontal_shift,r.Height()-vertical_shift,FALSE);
		}

		/*
		CWnd * pB = GetDlgItem(IDC_BUTTONOK_OUTPUT);
		if (pB)
		{
			pB->GetWindowRect(&r);
			ScreenToClient(r);
			int h = r.Height();
			int w = r.Width();
			r.left = (DialogRect.left + DialogRect.right - w)/2;
			r.right = (DialogRect.left + DialogRect.right + w)/2;
			r.bottom -= vertical_shift;
			r.top -= vertical_shift;
			
			pB->MoveWindow(r,FALSE);
		}*/

		RedrawWindow();

	} 

	CDialog::OnSize(nType, cx, cy);
}
