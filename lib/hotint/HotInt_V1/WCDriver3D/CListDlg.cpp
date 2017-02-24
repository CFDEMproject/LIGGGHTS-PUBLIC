//#**************************************************************
//# filename:             CListDlg.cpp
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
#include "CListDlg.h"


// CustomListDialog-Dialogfeld

IMPLEMENT_DYNAMIC(CustomListDialog, CDialog)
CustomListDialog::CustomListDialog(CWnd* pParent /*=NULL*/)
	: CDialog(CustomListDialog::IDD, pParent)
	, m_check_customlist(FALSE)
	, m_edit_text(_T(""))
{
	iscancel = 0; //tell: oneditkillfocus that it is called from cancel
	use_edit_item = 0;
	hwindow = 200;
	wwindow = 180; //180
}

CustomListDialog::~CustomListDialog()
{
	for (int i=1; i <= items.Length(); i++)
	{
		delete items(i);
	}
}

void CustomListDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Check(pDX, IDC_CUSTOMLIST_CHECK1, m_check_customlist);
	DDX_Control(pDX, IDC_CUSTOM_LIST, listctrl_custom);
	DDX_Text(pDX, IDC_EDIT2, m_edit_text);
}


BEGIN_MESSAGE_MAP(CustomListDialog, CDialog)
	ON_BN_CLICKED(IDCANCEL, OnBnClickedCancel)
	ON_BN_CLICKED(IDOK, OnBnClickedOk)
	ON_LBN_DBLCLK(IDC_CUSTOM_LIST, OnLbnDblclkCustomList)
	ON_EN_KILLFOCUS(IDC_EDIT2, OnEnKillfocusEdit2)
	ON_LBN_SELCHANGE(IDC_CUSTOM_LIST, OnLbnSelchangeCustomList)
END_MESSAGE_MAP()


void CustomListDialog::InitList(const char* dialognameI, const char* selectinfotextI, 
																int use_custom_check1I, const CString& checkstring1I)
{
	selectinfotext = CString(selectinfotextI);
	dialogname = CString(dialognameI);
	use_custom_check1 = use_custom_check1I;
	checkstring1 = checkstring1I;
	custom_check1 = 0;
	item_selected = 0;
}

// CustomListDialog-Meldungshandler

BOOL CustomListDialog::OnInitDialog()
{
	CDialog::OnInitDialog();


	//set new size of window:

	//position and size window:
	CRect rw, rwold, cwr, itemr;
	
	GetWindowRect(&rw);
	rwold = rw;
	GetClientRect(&cwr);

	hwindow += rw.Height() - cwr.Height();
	wwindow += rw.Width() - cwr.Width();

	if (hwindow < rw.Height()) hwindow = rw.Height();
	if (wwindow < rw.Width()) wwindow = rw.Width();

	CPoint p((rw.left+rw.right)/2,(rw.top+rw.bottom)/2);

	rw.left = p.x - (long)(0.5*wwindow);
	rw.top = p.y - (long)(0.*hwindow);
	rw.right = p.x + (long)(0.5*wwindow);
	rw.bottom = p.y + (long)(1.*hwindow);
	
	int offx=0, offy=0;
	if (rw.top < 1) offy = 1-rw.top;
	rw.top += offy;
	rw.bottom += offy;
	if (rw.left < 1) offx = 1-rw.left;
	rw.left += offx;
	rw.right += offx;

	MoveWindow(rw,FALSE);
	ShowWindow(SW_SHOW);

	int dy = rw.Height() - rwold.Height();
	int dx = rw.Width() - rwold.Width();
	listctrl_custom.GetWindowRect(&itemr);
	ScreenToClient(&itemr);
	itemr.bottom += dy;
	itemr.right += dx;
	int listwidth = itemr.Width();
	listctrl_custom.MoveWindow(itemr);


	listctrl_custom.ResetContent();
	for (int i=1; i <= items.Length(); i++)
	{
		listctrl_custom.AddString(*(items(i)));
	}
	listctrl_custom.SetSel(0);

	SetDlgItemText(IDC_STATIC_LISTINFO, selectinfotext);
	SetWindowText(dialogname);
	CButton* b = (CButton*)GetDlgItem(IDC_CUSTOMLIST_CHECK1);

	if (use_custom_check1)
	{
		b->SetWindowText(checkstring1);
		b->SetCheck(0);
		b->ModifyStyle(0,BS_AUTOCHECKBOX|BS_LEFT|WS_CHILD|WS_VISIBLE);

		b->GetWindowRect(&itemr);
		ScreenToClient(&itemr);
		itemr.top += dy;
		itemr.bottom += dy;
		b->MoveWindow(itemr);
	}
	else
	{
		b->SetButtonStyle(0);
	}

	CEdit* ce = (CEdit*)GetDlgItem(IDC_EDIT2);
	if (use_edit_item)
	{
		ce->GetWindowRect(&itemr);
		ScreenToClient(&itemr);
		int ewidth = itemr.Width();
		itemr.left += dx;
		itemr.right += dx;
		ce->MoveWindow(itemr);
		ce->ModifyStyle(ES_LEFT,ES_NUMBER|ES_AUTOHSCROLL|WS_TABSTOP|ES_RIGHT|WS_VISIBLE|WS_BORDER);		

		GetDlgItem(IDC_STATIC_LISTINFO)->GetWindowRect(&itemr);
		ScreenToClient(&itemr);
		itemr.right = itemr.left+listwidth-ewidth;
		GetDlgItem(IDC_STATIC_LISTINFO)->MoveWindow(itemr);
		GetDlgItem(IDC_STATIC_LISTINFO)->RedrawWindow();
	}
	else
	{
		//ce->ModifyStyle(0,0);
	}

	b = (CButton*)GetDlgItem(IDCANCEL);
	b->GetWindowRect(&itemr);
	ScreenToClient(&itemr);
	itemr.top += dy;
	itemr.bottom += dy;
	b->MoveWindow(itemr);

	b = (CButton*)GetDlgItem(IDOK);
	b->GetWindowRect(&itemr);
	ScreenToClient(&itemr);
	itemr.top += dy;
	itemr.bottom += dy;
	b->MoveWindow(itemr);

	UpdateData(FALSE);
	RedrawWindow();

	return TRUE;  // return TRUE unless you set the focus to a control
}

void CustomListDialog::AddString(const CString& str, int value) 
{
	CString* str2 = new CString(str);
	items.Add(str2);

	if (value != -1)
		values.Add(value);
	else
		values.Add(items.Length());
};


void CustomListDialog::OnBnClickedCancel()
{
	iscancel = 1;
	item_selected = 0;
	OnCancel();
}

void CustomListDialog::OnBnClickedOk()
{
	UpdateData(TRUE);

	custom_check1 = m_check_customlist;

	int num = atoi(m_edit_text);
	if (num > 0 && num <= items.Length()) 
	{
		item_selected = num;
	}
	else
	{
		int sel = listctrl_custom.GetCurSel();
		if (sel == LB_ERR)
		{
			item_selected = 0;
		}
		else
		{
			item_selected = values(sel+1); //GetCurSel() is 0-based!!!
		}
	}

	if (!item_selected)
	{
		AfxMessageBox("Invalid item selected!");
	}
	else 
	{
		OnOK();
	}
}

void CustomListDialog::OnLbnDblclkCustomList()
{
	OnBnClickedOk();
}

void CustomListDialog::OnEnKillfocusEdit2()
{
	if (!use_edit_item /*|| iscancel*/) return;

	UpdateData(TRUE);

	int num = atoi(m_edit_text);
	if (num > 0 && num <= items.Length())
	{
		listctrl_custom.SetCurSel(num-1); //SetCurSel() is 0-based!!!
	}
	else
	{
		//AfxMessageBox("An invalid item number has been entered!");
		CEdit* ce = (CEdit*)GetDlgItem(IDC_EDIT2);
		ce->SetWindowText("");
	}
}

void CustomListDialog::OnLbnSelchangeCustomList()
{
	if (!use_edit_item) return;

	UpdateData(TRUE);

	int sel = listctrl_custom.GetCurSel();
	if (sel == LB_ERR)
	{
		item_selected = 0;
	}
	else
	{
		item_selected = sel+1; //GetCurSel() is 0-based!!!
		char str[32];
		sprintf(str, "%d", item_selected);
		m_edit_text = str;
		UpdateData(FALSE);
	}

}
