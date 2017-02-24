//#**************************************************************
//# filename:             DlgOneEditCtrl.cpp
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
#include "DlgOneEditCtrl.h"


// CDialogOneEditControl-Dialogfeld

IMPLEMENT_DYNAMIC(CDialogOneEditControl, CDialog)
CDialogOneEditControl::CDialogOneEditControl(CWnd* pParent /*=NULL*/)
	: CDialog(CDialogOneEditControl::IDD, pParent)
{
	statictext = "";
	dialogname = "";
	edittext = "";
	value = 0;
	intvalue = 0;
	type = 1; //1==int, 2==double, 3=text, 4=watch
	minval = -1e100;
	maxval = 1e100;
	idok = 0;
}

void CDialogOneEditControl::InitDouble(const char* dialognameI, const char* statictextI, double* outputval, 
																			 double min, double max)
{
	statictext = CString(statictextI);
	dialogname = CString(dialognameI);
	minval = min;
	maxval = max;
	type = 2;
	value = outputval;
	idok = 0;
}

void CDialogOneEditControl::InitInt(const char* dialognameI, const char* statictextI, int* outputval, 
																			 double min, double max)
{
	statictext = CString(statictextI);
	dialogname = CString(dialognameI);
	minval = min;
	maxval = max;
	type = 1;
	intvalue = outputval;
	idok = 0;
}

void CDialogOneEditControl::InitText(const char* dialognameI, const char* statictextI)
{
	statictext = CString(statictextI);
	dialogname = CString(dialognameI);
	type = 3;
	idok = 0;
}

void CDialogOneEditControl::InitWatch(const char* dialognameI, int sensornumberI, int sensorwatchdialognumberI)
{
	dialogname = CString(dialognameI);
	sensornumber = sensornumberI;
	sensorwatchdialognumber = sensorwatchdialognumberI;
	type = 4;
	idok = 0;
}

CDialogOneEditControl::~CDialogOneEditControl()
{
}

void CDialogOneEditControl::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}


BEGIN_MESSAGE_MAP(CDialogOneEditControl, CDialog)
	ON_BN_CLICKED(IDOK, OnBnClickedOk)
	ON_BN_CLICKED(IDCANCEL, OnBnClickedCancel)
END_MESSAGE_MAP()


// CDialogOneEditControl-Meldungshandler

void CDialogOneEditControl::OnBnClickedOk()
{
	if (type == 1 && intvalue != 0)
	{
		double v = GetDlgItemInt(IDC_EDIT_TEXT);
		if (v >= minval && v <= maxval)
		{
			idok = 1;
		}
		else
		{
			idok = 0;
			v = minval;
		}
		*intvalue = (int)v;
	}
	else if (type == 2 && value != 0)
	{
		CString str;
		GetDlgItemText(IDC_EDIT_TEXT, str);
		double v = atof((LPCTSTR)str);
		if (v >= minval && v <= maxval)
		{
			idok = 1;
		}
		else
		{
			idok = 0;
			v = minval;
		}
		*value = v;
	}
	else if (type == 3)
	{
		idok = 1;
		GetDlgItemText(IDC_EDIT_TEXT, edittext);
	}
	else
	{
		idok = 0;
	}

	if (idok) OnOK(); //only accept if values are ok
	else
	{
		char str[200];
		if (type == 1) sprintf(str, "The input value must be between %d and %d !!!", (int)minval, (int)maxval);
		else if (type == 2) sprintf(str, "The input value must be between %g and %g !!!", minval, maxval);
		else sprintf(str, "An error occured during the input, check the values !!!");

		AfxMessageBox(str);
	}
}

void CDialogOneEditControl::OnBnClickedCancel()
{
	idok = 0;
	OnCancel();
}

BOOL CDialogOneEditControl::OnInitDialog()
{
	CDialog::OnInitDialog();

	SetDlgItemText(IDC_STATIC_TEXT, statictext);
	SetWindowText(dialogname);

	if (type == 1 || type == 2)
	{
		double vinit = 0;
		if (minval != -1e100) vinit = minval;

		char str[64];
		if (type == 1) sprintf(str, "%d", (int)vinit);
		if (type == 2) sprintf(str, "%.16g", vinit);

		SetDlgItemText(IDC_EDIT_TEXT, str);
	}
	
	return TRUE;  // return TRUE unless you set the focus to a control
}

void CDialogOneEditControl::SetSensorText(const char* str) 
{
	SetDlgItemText(IDC_EDIT_TEXT, str);
	//RedrawWindow();	
}


void CDialogOneEditControl::Create(CWnd * pParent)
{
	CDialog::Create(IDD,pParent);

	UpdateData(FALSE);


	CRect r, cr, cer;
	int min_width = 200;

	if (type == 4) //watch
	{
		SetWindowText(dialogname);

		CWnd* w = GetDlgItem(IDC_STATIC_TEXT);
		w->ModifyStyle(WS_VISIBLE, 0, 0);
		w = GetDlgItem(IDOK);
		w->ModifyStyle(WS_VISIBLE, 0, 0);
		w = GetDlgItem(IDCANCEL);
		w->ModifyStyle(WS_VISIBLE, 0, 0);
		
		CEdit* ce = (CEdit*)GetDlgItem(IDC_EDIT_TEXT);
		//DWORD sty = ce->GetStyle();
		//ce->ModifyStyle(0, sty|ES_READONLY);
		ce->SetReadOnly(TRUE);

		ce->GetWindowRect(cer);
		ScreenToClient(cer);
		int dx = cer.left - 1;
		int dy = cer.top - 1;
		cer.left -= dx;
		cer.right -= dx;
		if (cer.Width() < min_width) cer.right += min_width - cer.Width();
		cer.top -= dy;
		cer.bottom -= dy;
		ce->MoveWindow(cer, FALSE);
		ce->ShowWindow(SW_SHOW);


		pParent->GetWindowRect(&r);
		CPoint p(r.left,r.bottom - 5);
		GetWindowRect(&r);
		p.x += r.Width();

		GetClientRect(&cr);
		int wi = cer.Width();
		if (wi < min_width) wi = min_width;

		int decx = cr.Width() - wi - 2;
		if (decx < 0 ) decx = 0;
		int decy = cr.Height() - cer.Height() - 2;
		if (decy < 0 ) decy = 0;

		r.left = p.x - r.Width() + decx;
		r.top = p.y - r.Height() + decy;
		r.right = p.x;
		r.bottom = p.y;
		MoveWindow(r,FALSE);
		ShowWindow(SW_SHOW);

	}
	else
	{
		pParent->GetWindowRect(&r);
		CPoint p((r.left+r.right)/2,(r.bottom+r.top)/2);
		GetWindowRect(&r);
		p.x += r.Width();

		r.top = p.y - r.Height();
		r.left = p.x - r.Width();
		r.right = p.x;
		r.bottom = p.y;
		MoveWindow(r,FALSE);
		ShowWindow(SW_SHOW);
	}

}

void CDialogOneEditControl::OnOK()
{
	if (type == 4)
		DestroyWindow();
	else
		CDialog::OnOK();
}

void CDialogOneEditControl::OnClose() 
{
	if (type == 4)
		DestroyWindow();
	else
		CDialog::OnClose();
}

void CDialogOneEditControl::OnCancel() 
{
	if (type == 4)
		DestroyWindow();
	else
		CDialog::OnCancel();
}
 



